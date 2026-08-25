"""Per-granule SYP-3 excess: which fraction of granules "light up"? (Ryan, 2026-08-24)
The mean over all granules is dominated by blur from the adjacent nucleus and does not separate sexes.
The eye counts granules that are much brighter than their local surroundings. Per granule:
    excess = (mean SYP in granule) - (mean SYP in non-granule cytoplasm at the SAME distance from the
             envelope, 0.25 um bins)         ... removes the distance-dependent nuclear blur
    excess_n = excess / (nuclear SYP - bg)   ... in units of the nuclear level, exposure-free
Per gonad we report the distribution: fraction of granules with excess_n > 0.25 and > 0.5, and the
90th percentile of excess_n. Output: granule_tail_v3.csv (per gonad) + granule_tail_pergranule_v3.parquet
(per granule, with distance to the envelope and PGL brightness for the radius and bleed-through checks).
v4 (2026-08-24): labels with no lamin envelope (sperm / debris, nucleus_filter.no_ring_ids) are dropped from
the envelope set before the shell is built; v3 files are kept for the record."""
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
import nucleus_filter as NF

SP = chload.SP
SPACING = tuple(SP)
VVOL = float(SP.prod())
CYTO = 2.5
SMOOTH_UM, BG_UM, K, VMIN, VMAX = 0.15, 0.513, 4.48, 0.003, 15.0
DB = np.arange(0.0, CYTO + 1e-6, 0.25)
OUT = r"C:/Users/ryane/coloc_analysis/granule_tail_v4.csv"
META = pd.read_csv(r"C:/Users/ryane/coloc_analysis/acquisition_metadata.csv").set_index("iid")
IIDS = [i for i in META.index if not chload.is_excluded(i)]

PARQ = r"C:/Users/ryane/coloc_analysis/granule_tail_pergranule_v4.parquet"
rows, per = [], []
if os.path.exists(OUT):
    rows = [r for r in pd.read_csv(OUT).to_dict("records") if r["image_id"] in IIDS]   # drop rows excluded since
    if os.path.exists(PARQ):                 # resume must keep the per-granule rows of finished gonads
        old = pd.read_parquet(PARQ)
        per = [old[old.image_id.isin(IIDS)]]
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
    lamin_raw = ch["lamin"][sl]
    del ch
    junk = NF.no_ring_ids(NF.ring_scores(lab_c, germ, lamin_raw))      # v4: no-envelope labels out
    lam = np.where(np.isin(lam, junk), 0, lam)
    env = lam > 0
    dt = ndi.distance_transform_edt(~env, sampling=SPACING)
    cyto = ndi.binary_fill_holes(dt <= CYTO) & ~env
    sm = ndi.gaussian_filter(pgl, sigma=SMOOTH_UM / SP)
    th = sm - ndi.gaussian_filter(sm, sigma=BG_UM / SP)
    rr = th[cyto]
    neg = rr[rr < 0]
    noise = float(np.sqrt(np.mean(neg ** 2))) if neg.size > 50 else float(np.std(rr))
    lab2, nlab = ndi.label((th > K * noise) & cyto)
    if nlab == 0:
        print(f"{iid}: no granules segmented, skipping", flush=True)
        continue
    vol = np.bincount(lab2.ravel()) * VVOL
    keep = np.where((vol >= VMIN) & (vol <= VMAX))[0]
    keep = keep[keep != 0]
    gran = np.isin(lab2, keep)
    far = dt > 10.0
    bg = float(np.percentile(syp[far], 0.5)) if far.sum() > 1e4 else float(np.percentile(syp, 0.5))
    N = float(syp[env].mean()) - bg
    # local (distance-matched) cytoplasm SYP per distance bin
    outside = cyto & ~gran
    binidx = np.clip(np.digitize(dt, DB) - 1, 0, len(DB) - 2)
    cyto_by_bin = np.array([float(syp[outside & (binidx == b)].mean()) - bg if (outside & (binidx == b)).sum() > 50 else np.nan
                            for b in range(len(DB) - 1)])
    # per-granule mean SYP and mean distance bin
    idx = np.arange(1, nlab + 1)
    g_syp = ndi.mean(syp, lab2, idx) - bg
    pbg = float(np.percentile(pgl[far], 0.5)) if far.sum() > 1e4 else float(np.percentile(pgl, 0.5))
    g_pgl = ndi.mean(pgl, lab2, idx) - pbg          # per-granule PGL brightness, for bleed-through estimation
    g_bin = ndi.mean(binidx.astype(np.float32), lab2, idx)
    g_dist = ndi.mean(dt, lab2, idx)                # mean distance of the granule's voxels from the envelope, um
    g_dmin = ndi.minimum(dt, lab2, idx)             # closest approach to the envelope, um
    g_vol = vol[1:nlab + 1]
    ok = np.isin(idx, keep)
    g_bin_i = np.clip(np.round(g_bin[ok]).astype(int), 0, len(DB) - 2)
    local = cyto_by_bin[g_bin_i]
    excess_n = (g_syp[ok] - local) / N
    fin = np.isfinite(excess_n)
    g_pgl_n = (g_pgl[ok] / max(float(pgl[env].mean()) - pbg, 1.0))[fin]   # PGL in units of nuclear-region PGL (haze-normalised)
    g_pgl_raw = g_pgl[ok][fin]
    excess_n = excess_n[fin]
    p = chload.parse_iid(iid)
    row = {"image_id": iid, "sex": p["sex"], "treat": p["treat"], "batch": p["batch"], "n_gran": int(ok.sum()), "n_no_ring_dropped": len(junk),
           "frac_excess_gt_0.25": round(float((excess_n > 0.25).mean()), 4),
           "frac_excess_gt_0.5": round(float((excess_n > 0.5).mean()), 4),
           "frac_excess_gt_1.0": round(float((excess_n > 1.0).mean()), 4),
           "excess_p50": round(float(np.percentile(excess_n, 50)), 4),
           "excess_p90": round(float(np.percentile(excess_n, 90)), 4),
           "excess_p99": round(float(np.percentile(excess_n, 99)), 4)}
    rows.append(row)
    per.append(pd.DataFrame({"image_id": iid, "sex": p["sex"], "treat": p["treat"], "excess_n": excess_n,
                             "pgl_raw": g_pgl_raw, "pgl_n": g_pgl_n, "N_syp": N, "exp_477": float(META.loc[iid, "exp_477"]),
                             "exp_545": float(META.loc[iid, "exp_545"]), "vol_um3": g_vol[ok][fin],
                             "dist_um": g_dist[ok][fin], "dmin_um": g_dmin[ok][fin]}))
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"{iid[-18:]:20s} {p['sex']:4s} {p['treat']:4s} n={ok.sum():5d} | >0.25: {row['frac_excess_gt_0.25']:.3f}  >0.5: {row['frac_excess_gt_0.5']:.3f}  >1.0: {row['frac_excess_gt_1.0']:.3f} | p90={row['excess_p90']:.3f}", flush=True)
if per:
    new = [p_ for p_ in per if p_ is not per[0] or not os.path.exists(PARQ)]   # frames produced in THIS run
    fresh = {p_.image_id.iloc[0] for p_ in new if len(p_)}
    kept = per[0][~per[0].image_id.isin(fresh)] if (per and os.path.exists(PARQ)) else per[0].iloc[0:0]
    pd.concat([kept] + new).to_parquet(PARQ, index=False)          # reprocessed gonads replace their old rows
print("DONE", flush=True)
