"""Stage-resolved partition coefficient using RYAN'S HAND-TRACED pachytene zones.
Usage: python pc_zone_worker.py <image_id>

This is the successor to the lost-axis fig70 analysis, now on defensible staging: Ryan drew the
pachytene region on each gonad (trace_pachytene.py, 2026-08-23); zones = equal-length thirds of his
polyline (early / mid / late), assigned per nucleus in staging/zones/<iid>_zones.csv.

Per gonad:
  nuclei     = lamin envelopes from the staging cache (the PI-requested masking)
  granules   = Imaris-calibrated recipe on raw PGL inside the 2.5 um perinuclear shell
  zone of a VOXEL = zone of its nearest germline nucleus (EDT indices), so granules inherit the zone
                    of the nucleus they dock on - same convention as the original fig70
  per zone   = distance-matched PC (raw SYP) + z-shift control + PC_specific = PC / zsh (per-gonad,
               per-zone paired floor)
Output: one JSON row -> coloc_analysis/pc_zone_rows/<iid>.json"""
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import chload
import crescent_axis as ca

SP = chload.SP
SPACING = tuple(SP)
VVOL = float(SP.prod())
CYTO_UM = 2.5
SMOOTH_UM, BG_UM, K_THRESH, VMIN, VMAX = 0.15, 0.513, 4.48, 0.003, 15.0
DZ = 6
DBINS = np.arange(0.0, CYTO_UM + 1e-6, 0.25)
ZONES = ["early", "mid", "late"]
OUTDIR = r"C:/Users/ryane/coloc_analysis/pc_zone_rows"
os.makedirs(OUTDIR, exist_ok=True)


def segment_granules(pgl, cyto):
    sm = ndi.gaussian_filter(pgl, sigma=SMOOTH_UM / SP)
    th = sm - ndi.gaussian_filter(sm, sigma=BG_UM / SP)
    rr = th[cyto]
    neg = rr[rr < 0]
    noise = float(np.sqrt(np.mean(neg ** 2))) if neg.size > 50 else float(np.std(rr))
    lab, _ = ndi.label((th > K_THRESH * noise) & cyto)
    if lab.max() == 0:
        return np.zeros(pgl.shape, bool), 0
    vol = np.bincount(lab.ravel()) * VVOL
    keep = np.where((vol >= VMIN) & (vol <= VMAX))[0]
    keep = keep[keep != 0]
    return np.isin(lab, keep), len(keep)


def pc_dm(syp, gran, outside, dt, bg, region):
    g0 = gran & region
    o0 = outside & region
    if int(g0.sum()) < 30 or int(o0.sum()) < 200:
        return None
    num = wsum = 0.0
    for lo, hi in zip(DBINS[:-1], DBINS[1:]):
        b = (dt > lo) & (dt <= hi)
        g = g0 & b
        o = o0 & b
        if int(g.sum()) < 15 or int(o.sum()) < 40:
            continue
        den = float(syp[o].mean()) - bg
        if den <= 0:
            continue
        num += int(g.sum()) * ((float(syp[g].mean()) - bg) / den)
        wsum += int(g.sum())
    return round(num / wsum, 4) if wsum > 0 else None


def main(iid):
    tr = json.load(open(r"C:/Users/ryane/coloc_analysis/staging/pachytene_traces.json", encoding="utf-8"))
    if tr.get(iid, {}).get("status") != "traced":
        raise SystemExit(f"{iid}: trace status is {tr.get(iid, {}).get('status')!r}, not 'traced' - refusing")
    if chload.is_excluded(iid):
        raise SystemExit(f"{iid}: excluded by the mask-contamination audit (exclusions.json) - refusing")
    zones = pd.read_csv(rf"C:/Users/ryane/coloc_analysis/staging/zones/{iid}_zones.csv")
    dapi, syp_c, lam, lab_c, fb = ca.load_crops(iid)          # cached crop, same frame as the trace
    # need the PGL channel, not in the cache: reload the nd2 crop the same way
    rd = ca.find_run(iid)
    csv = os.path.join(rd, f"{iid}__nuclei.csv")
    nc = pd.read_csv(csv)
    germ = [int(x) for x in nc[nc.in_germline.astype(bool)].nucleus_id]
    import tifffile
    full = tifffile.imread(os.path.join(rd, f"{iid}__nuclei_labels.tif"))
    gn = np.isin(full, germ)
    obj = ndi.find_objects(gn.astype(np.uint8))[0]
    sl = tuple(slice(max(0, o.start - ca.PAD), min(dim, o.stop + ca.PAD)) for o, dim in zip(obj, gn.shape))
    ch = chload.load(ca.find_nd2(iid, csv))
    pgl = ch["pgl"][sl].astype(np.float32)
    syp = ch["syp"][sl].astype(np.float32)
    del ch
    assert pgl.shape == lam.shape, f"crop mismatch {pgl.shape} vs {lam.shape}"

    env = lam > 0
    dt, inds = ndi.distance_transform_edt(~env, sampling=SPACING, return_indices=True)
    cyto = ndi.binary_fill_holes(dt <= CYTO_UM) & ~env
    gran, ngran = segment_granules(pgl, cyto)
    outside = cyto & ~gran
    # per-voxel zone via nearest envelope voxel's nucleus id
    nn = lam[tuple(inds)]
    zmap = dict(zip(zones.nucleus_id.astype(int), zones.zone))
    maxl = int(lam.max())
    lut = np.zeros(maxl + 1, np.int8)          # 0=none 1=early 2=mid 3=late
    code = {"early": 1, "mid": 2, "late": 3}
    for nid, z in zmap.items():
        if 0 <= nid <= maxl and z in code:
            lut[nid] = code[z]
    zone_vox = lut[nn]
    bg = float(np.percentile(syp, 3))
    gz = np.zeros_like(gran)
    gz[DZ:] = gran[:-DZ]
    gz &= cyto
    ozr = cyto & ~gz
    out = {"image_id": iid,
           "sex": chload.parse_iid(iid)["sex"], "treat": chload.parse_iid(iid)["treat"],
           "batch": chload.parse_iid(iid)["batch"],
           "n_gran": ngran,
           "off_axis_cut_um": float(zones.off_axis_cut_um.iloc[0]) if "off_axis_cut_um" in zones else None,
           "n_zoned_nuclei": int(zones.zone.isin(list(code)).sum())}   # provenance: which zoning produced this row
    # "pach" = the three zones pooled: the whole hand-traced pachytene region as ONE region. This is the
    # contamination-free replacement for the whole-gonad PC (2026-08-24 mask audit: whole-gonad masks
    # carry sperm and somatic nuclei; the traced region does not).
    for z, c in list(code.items()) + [("pach", None)]:
        reg = (zone_vox > 0) if c is None else (zone_vox == c)
        dm = pc_dm(syp, gran, outside, dt, bg, reg)
        zs = pc_dm(syp, gz, ozr, dt, bg, reg)
        out[f"{z}_PC"] = dm
        out[f"{z}_zsh"] = zs
        out[f"{z}_PCspec"] = round(dm / zs, 4) if (dm is not None and zs is not None and zs > 0) else None
        out[f"{z}_vox"] = int(reg.sum())
    with open(os.path.join(OUTDIR, iid + ".json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out))


if __name__ == "__main__":
    main(sys.argv[1])
