"""AMOUNT of SYP-3 in P granules, not just its partitioning ratio (Ryan, 2026-08-24: "I can see overlap in
males but not herms, the graph says they're the same"). The partition coefficient divides granule SYP by
neighbouring cytoplasmic SYP; in HS males both rise, in herms both are near zero, so the ratio hides the
difference the eye sees. Per gonad, background-subtracted:
    G = mean SYP in granule voxels,  N = mean SYP inside nuclear envelopes,  C = mean SYP in the
    cytoplasm shell excluding granules,  bg = camera/dark floor (0.5th pct of extra-worm voxels)
Reported:  gran_vs_nuc = (G-bg)/(N-bg)   fraction of the nuclear SYP level found in granules
           cyto_vs_nuc = (C-bg)/(N-bg)   how much SYP has left the nucleus into the cytoplasm
           gran_abs    = (G-bg)/exposure_ms  absolute granule SYP in counts per ms (laser ~constant)
           PC          = gran_vs_nuc / cyto_vs_nuc (sanity: reproduces the whole-gonad, non-distance-matched PC)"""
import json
import os
import sys

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import chload
import crescent_axis as ca

SP = chload.SP
SPACING = tuple(SP)
VVOL = float(SP.prod())
CYTO = 2.5
SMOOTH_UM, BG_UM, K, VMIN, VMAX = 0.15, 0.513, 4.48, 0.003, 15.0
OUT = r"C:/Users/ryane/coloc_analysis/amount_metrics.csv"
META = pd.read_csv(r"C:/Users/ryane/coloc_analysis/acquisition_metadata.csv").set_index("iid")
IIDS = [i for i in META.index if not chload.is_excluded(i)]

rows = []
if os.path.exists(OUT):
    rows = pd.read_csv(OUT).to_dict("records")
done = {r["image_id"] for r in rows}
for iid in IIDS:
    if iid in done:
        continue
    dapi, syp_c, lam, lab_c, fb = ca.load_crops(iid)
    rd = ca.find_run(iid)
    csv = os.path.join(rd, f"{iid}__nuclei.csv")
    nc = pd.read_csv(csv)
    germ = [int(x) for x in nc[nc.in_germline.astype(bool)].nucleus_id]
    full = tifffile.imread(os.path.join(rd, f"{iid}__nuclei_labels.tif"))
    gn = np.isin(full, germ)
    obj = ndi.find_objects(gn.astype(np.uint8))[0]
    sl = tuple(slice(max(0, o.start - ca.PAD), min(dim, o.stop + ca.PAD)) for o, dim in zip(obj, gn.shape))
    ch = chload.load(ca.find_nd2(iid, csv))
    pgl = ch["pgl"][sl].astype(np.float32)
    syp = ch["syp"][sl].astype(np.float32)
    del ch
    assert pgl.shape == lam.shape
    env = lam > 0
    dt = ndi.distance_transform_edt(~env, sampling=SPACING)
    cyto = ndi.binary_fill_holes(dt <= CYTO) & ~env
    sm = ndi.gaussian_filter(pgl, sigma=SMOOTH_UM / SP)
    th = sm - ndi.gaussian_filter(sm, sigma=BG_UM / SP)
    rr = th[cyto]
    neg = rr[rr < 0]
    noise = float(np.sqrt(np.mean(neg ** 2))) if neg.size > 50 else float(np.std(rr))
    lab2, _ = ndi.label((th > K * noise) & cyto)
    vol = np.bincount(lab2.ravel()) * VVOL
    keep = np.where((vol >= VMIN) & (vol <= VMAX))[0]
    keep = keep[keep != 0]
    gran = np.isin(lab2, keep)
    far = dt > 10.0
    bg = float(np.percentile(syp[far], 0.5)) if far.sum() > 1e4 else float(np.percentile(syp, 0.5))
    G = float(syp[gran].mean()) - bg
    N = float(syp[env].mean()) - bg
    C = float(syp[cyto & ~gran].mean()) - bg
    exp = float(META.loc[iid, "exp_545"])
    p = chload.parse_iid(iid)
    row = {"image_id": iid, "sex": p["sex"], "treat": p["treat"], "batch": p["batch"], "exp_545": exp,
           "bg": round(bg, 1), "G": round(G, 1), "N": round(N, 1), "C": round(C, 1), "n_gran": len(keep),
           "gran_vs_nuc": round(G / N, 4), "cyto_vs_nuc": round(C / N, 4),
           "gran_abs_per_ms": round(G / exp, 4), "nuc_abs_per_ms": round(N / exp, 4),
           "PC_plain": round(G / C, 4)}
    rows.append(row)
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"{iid[-18:]:20s} {p['sex']:4s} {p['treat']:4s} G={G:6.1f} N={N:6.1f} C={C:6.1f} | gran/nuc={G/N:.3f} cyto/nuc={C/N:.3f} PC={G/C:.2f}", flush=True)
print("DONE", flush=True)
