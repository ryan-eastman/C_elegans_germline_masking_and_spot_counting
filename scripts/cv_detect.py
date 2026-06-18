"""Cross-gonad robust-detection experiment. Reuse each gonad's saved cellpose labels (skip the 18-min
segmentation) and re-run only SpotMAX detection on the RAD-51 channel with otsu / li / triangle, then
sweep effect_size_min. For each (method, effect_size_min) report the per-gonad matched per-nucleus
mean vs Imaris coloc and the pooled cross-gonad agreement (CCC/bias). Finds the combo that is robust
across all 6 gonads, not overfit to one. Writes results_cv_detect/grid.csv + pooled.csv."""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd
import tifffile

import germquant.io.nd2_reader as R
from germquant.spots import detect_spots
from germquant.validate.imaris_ims import read_imaris_ims
from germquant.validate.imaris_xlsx import per_nucleus_rad51
from scripts.cv_analyze import _ccc, _greedy

_ap = argparse.ArgumentParser()
_ap.add_argument("--manifest", default="cv_manifest.json")
_ap.add_argument("--ndir", default=r"E:\Madeleine\N2\20251105_N2_noHS", help="folder with the .nd2 files")
_ap.add_argument("--out", default="results_cv_detect")
_args = _ap.parse_args()

SP = (0.2, 0.108333, 0.108333)
NDIR = _args.ndir
OUT = _args.out
TOL = 2.5
METHODS = ["threshold_otsu", "threshold_li", "threshold_triangle"]
ESMINS = [0.0, 2.0, 3.0, 4.0, 5.0]
os.makedirs(OUT, exist_ok=True)


def matched_for(g):
    """matched (our_nucleus_id, coloc_count) for a gonad, from saved nuclei + .ims + xlsx."""
    nuc = pd.read_csv(glob.glob(os.path.join(g["out_dir"], "*__nuclei.csv"))[0])
    germ = nuc[nuc["in_germline"] == True]  # noqa: E712
    _, im_surf, geom = read_imaris_ims(g["ims"])
    xsurf, _ = per_nucleus_rad51(g["xlsx"])
    emin = geom["ext_min_xyz"]
    x_img = np.c_[xsurf["X"] - emin[0], xsurf["Y"] - emin[1], xsurf["Z"] - emin[2]]
    surf = im_surf[["x_um", "y_um", "z_um"]].to_numpy()
    coloc = np.full(len(im_surf), np.nan)
    for xi, sj in _greedy(x_img, surf, TOL):
        coloc[sj] = xsurf["RAD51"].to_numpy()[xi]
    oxyz = germ[["centroid_x_um", "centroid_y_um", "centroid_z_um"]].to_numpy()
    oid = germ["nucleus_id"].to_numpy()
    return [(int(oid[oj]), coloc[sj]) for sj, oj in _greedy(surf, oxyz, TOL) if not np.isnan(coloc[sj])]


gonads = json.load(open(_args.manifest))
matched = {g["name"]: matched_for(g) for g in gonads}
nd2_of = {g["name"]: os.path.join(NDIR, g["name"] + ".nd2") for g in gonads}
lab_of = {g["name"]: glob.glob(os.path.join(g["out_dir"], "*__nuclei_labels.tif"))[0] for g in gonads}

grid, pooled = [], {}
for g in gonads:
    name = g["name"]
    print(f"\n=== {name}: reading channel + labels ...", flush=True)
    rad = np.asarray(R.read_stack(nd2_of[name]).data[2]).astype("float32")
    lab = tifffile.imread(lab_of[name]).astype("int32")
    cmean = float(np.mean([c for _, c in matched[name]]))
    for method in METHODS:
        ps, _ = detect_spots(rad, lab, SP, thresholding_method=method, effect_size_min=0.0,
                             merge_z_columns=True)
        for es in ESMINS:
            f = ps[ps["effect_size"].fillna(np.inf) >= es]  # match production detect.py (keep NaN-ES)
            per = f.groupby("nucleus_id").size()
            oc = [per.get(i, 0) for i, _ in matched[name]]
            cc = [c for _, c in matched[name]]
            grid.append({"gonad": name, "method": method, "esmin": es,
                         "total": int(len(f)), "matched_mean": float(np.mean(oc)), "coloc_mean": cmean})
            pooled.setdefault((method, es), ([], []))
            pooled[(method, es)][0].extend(oc)
            pooled[(method, es)][1].extend(cc)
        print(f"  {method:20s} detected={len(ps):6d}  "
              f"matched_mean@es0={np.mean([ps.groupby('nucleus_id').size().get(i,0) for i,_ in matched[name]]):.2f} "
              f"(coloc {cmean:.2f})", flush=True)
    del rad, lab

pd.DataFrame(grid).to_csv(os.path.join(OUT, "grid.csv"), index=False)
rows = []
for (method, es), (o, c) in pooled.items():
    o, c = np.array(o, float), np.array(c, float)
    rows.append({"method": method, "esmin": es, "n": len(o), "ccc": _ccc(o, c),
                 "pearson": float(np.corrcoef(o, c)[0, 1]) if np.std(o) > 0 else np.nan,
                 "bias": float(o.mean() - c.mean())})
pl = pd.DataFrame(rows).sort_values("ccc", ascending=False)
pl.to_csv(os.path.join(OUT, "pooled.csv"), index=False)
print("\n=== pooled cross-gonad agreement (sorted by CCC) ===")
print(pl.to_string(index=False))
best = pl.iloc[0]
print(f"\nBEST: {best['method']} esmin={best['esmin']}  CCC={best['ccc']:.3f} bias={best['bias']:+.2f}")
print("\n=== per-gonad matched_mean at best config ===")
gd = pd.DataFrame(grid)
print(gd[(gd["method"] == best["method"]) & (gd["esmin"] == best["esmin"])]
      [["gonad", "total", "matched_mean", "coloc_mean"]].to_string(index=False))
