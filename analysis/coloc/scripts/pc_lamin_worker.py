"""Per-gonad worker: LAMIN-masked partition coefficient (PI-requested masking change), computed
alongside the DAPI-masked version for a paired comparison. Usage: python pc_lamin_worker.py <image_id>

Masking by lamin: the validated DAPI/Cellpose nuclei are used as SEEDS for a 3D watershed whose
elevation map is the smoothed LMN-1 channel, restricted to a 2 um band around the DAPI nuclei. The
watershed boundary snaps to the lamin ridge = the true nuclear envelope. Everything inside the envelope
is excluded; the cytoplasm shell, granule segmentation (Imaris-calibrated recipe, unchanged), and the
distance-matched PC are all anchored to the lamin surface instead of the chromatin edge.
Per-nucleus safety gate: a lamin region must be 1.0-3.0x its DAPI volume, else that nucleus falls back
to its DAPI mask (count reported). Outputs one JSON row to coloc_analysis/pc_lamin_rows/<iid>.json with
BOTH metrics (lamin-masked and DAPI-masked) + rotation/z-shift controls for each.
QC: use scripts/qc_mask_audit.py (all-depth overlay + per-label envelope test) on any gonad."""
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage.segmentation import watershed

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
import chload

SP = chload.SP
SPACING = tuple(SP)
VVOL = float(SP.prod())
CYTO_UM = 2.5
LAM_BAND_UM = 2.0            # watershed search band beyond the DAPI nucleus
LAM_SMOOTH_UM = 0.15
SMOOTH_UM, BG_UM, K_THRESH, VMIN, VMAX = 0.15, 0.513, 4.48, 0.003, 15.0
DZ = 6
DBINS = np.arange(0.0, CYTO_UM + 1e-6, 0.25)
PAD = 30
OUTDIR = r"C:/Users/ryane/coloc_analysis/pc_lamin_rows"
os.makedirs(OUTDIR, exist_ok=True)


sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import crescent_axis as _ca

find_run = _ca.find_run      # shared data-location resolution (one list of run/nd2 dirs, incl. N2 dry-ice)
find_nd2 = _ca.find_nd2


def lamin_nuclei(lab_crop, germ_ids, lamin):
    """expand each germline DAPI nucleus to its lamin envelope via seeded watershed."""
    gn_dapi = np.isin(lab_crop, germ_ids)
    dt = ndi.distance_transform_edt(~gn_dapi, sampling=SPACING)
    band = dt <= LAM_BAND_UM
    elev = ndi.gaussian_filter(lamin.astype(np.float32), sigma=LAM_SMOOTH_UM / SP)
    markers = np.where(gn_dapi, lab_crop, 0).astype(np.int32)
    BGL = int(lab_crop.max()) + 1
    markers[dt > LAM_BAND_UM] = BGL
    ws = watershed(elev, markers=markers, mask=band | gn_dapi | (dt > LAM_BAND_UM))
    lam_lab = np.where(ws == BGL, 0, ws)
    # per-nucleus sanity gate vs DAPI volume
    dapi_vol = np.bincount(np.where(gn_dapi, lab_crop, 0).ravel())
    lam_vol = np.bincount(lam_lab.ravel(), minlength=len(dapi_vol))
    fell_back = 0
    out = np.zeros_like(gn_dapi)
    for nid in germ_ids:
        dv = dapi_vol[nid] if nid < len(dapi_vol) else 0
        lv = lam_vol[nid] if nid < len(lam_vol) else 0
        if dv > 0 and 1.0 <= lv / dv <= 3.0:
            out |= lam_lab == nid
        else:
            out |= (lab_crop == nid) & gn_dapi
            fell_back += 1
    return out, gn_dapi, fell_back


def segment_granules(pgl, cyto):
    sm = ndi.gaussian_filter(pgl, sigma=SMOOTH_UM / SP)
    th = sm - ndi.gaussian_filter(sm, sigma=BG_UM / SP)
    rr = th[cyto]; neg = rr[rr < 0]
    noise = float(np.sqrt(np.mean(neg ** 2))) if neg.size > 50 else float(np.std(rr))
    lab, _ = ndi.label((th > K_THRESH * noise) & cyto)
    if lab.max() == 0:
        return np.zeros(pgl.shape, bool), 0
    vol = np.bincount(lab.ravel()) * VVOL
    keep = np.where((vol >= VMIN) & (vol <= VMAX))[0]; keep = keep[keep != 0]
    return np.isin(lab, keep), len(keep)


def pc_dm(syp, gran, outside, dt, bg):
    if int(gran.sum()) < 50 or int(outside.sum()) < 500:
        return None
    num = wsum = 0.0
    for lo, hi in zip(DBINS[:-1], DBINS[1:]):
        b = (dt > lo) & (dt <= hi); g = gran & b; o = outside & b
        if int(g.sum()) < 20 or int(o.sum()) < 50:
            continue
        den = float(syp[o].mean()) - bg
        if den <= 0:
            continue
        num += int(g.sum()) * ((float(syp[g].mean()) - bg) / den); wsum += int(g.sum())
    return round(num / wsum, 4) if wsum > 0 else None


def trio(syp, gran, cyto, dt, bg):
    dm = pc_dm(syp, gran, cyto & ~gran, dt, bg)
    rg = np.rot90(gran, 2, axes=(1, 2)) & cyto        # rotation null: same mask, elsewhere in the cytoplasm
    rot = pc_dm(syp, rg, cyto & ~rg, dt, bg)          # its own outside set, like the z-shift control
    gz = np.zeros_like(gran); gz[DZ:] = gran[:-DZ]; gz &= cyto
    zsh = pc_dm(syp, gz, cyto & ~gz, dt, bg)
    return dm, rot, zsh


def metric_set(nuc_mask, pgl, syp, bg):
    dt = ndi.distance_transform_edt(~nuc_mask, sampling=SPACING)
    cyto = ndi.binary_fill_holes(dt <= CYTO_UM) & ~nuc_mask
    if cyto.sum() < 5000:
        return None
    gran, ng = segment_granules(pgl, cyto)
    if ng < 40:
        return None
    dm, rot, zsh = trio(syp, gran, cyto, dt, bg)
    return {"n_gran": ng, "PC_dm": dm, "rot": rot, "zsh": zsh,
            "PC_specific": round(dm / zsh, 4) if (dm is not None and zsh is not None and zsh > 0) else None}


def main(iid):
    rd = find_run(iid)
    lab = tifffile.imread(os.path.join(rd, f"{iid}__nuclei_labels.tif"))
    csv = os.path.join(rd, f"{iid}__nuclei.csv")
    nc = pd.read_csv(csv)
    germ = [int(x) for x in nc[nc.in_germline.astype(bool)].nucleus_id]
    gn = np.isin(lab, germ)
    obj = ndi.find_objects(gn.astype(np.uint8))[0]
    sl = tuple(slice(max(0, o.start - PAD), min(dim, o.stop + PAD)) for o, dim in zip(obj, gn.shape))
    ch = chload.load(find_nd2(iid, csv))
    pgl = ch["pgl"][sl].astype(np.float32); syp = ch["syp"][sl].astype(np.float32)
    lamin = ch["lamin"][sl]; lab_c = lab[sl]
    del ch
    bg = float(np.percentile(syp, 3))
    nuc_lam, nuc_dapi, fell_back = lamin_nuclei(lab_c, germ, lamin)
    res = {"image_id": iid,
           "sex": chload.parse_iid(iid)["sex"], "treat": chload.parse_iid(iid)["treat"],
           "batch": chload.parse_iid(iid)["batch"],
           "n_nuclei": len(germ), "lamin_fallback_n": fell_back,
           "lam_vol_ratio": round(float(nuc_lam.sum()) / max(1, float(nuc_dapi.sum())), 3),
           "lamin": metric_set(nuc_lam, pgl, syp, bg),
           "dapi": metric_set(nuc_dapi, pgl, syp, bg)}
    with open(os.path.join(OUTDIR, iid + ".json"), "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res))


if __name__ == "__main__":
    main(sys.argv[1])
