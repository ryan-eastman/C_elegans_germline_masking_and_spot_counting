"""Diagnose the cross-gonad spot-count instability. For each gonad, match traced Imaris surfaces to
our nuclei (and attach xlsx coloc counts), then sweep effect_size_min on OUR spots and recompute the
matched per-nucleus mean. Asks: can a single effect_size_min threshold bring all 6 gonads into
agreement with Imaris coloc, or is the per-image otsu detection too unstable for one threshold?
Also reports each gonad's effect-size distribution (to see noise-flooded vs under-detected images)."""
import json

import numpy as np
import pandas as pd

from germquant.validate.imaris_ims import read_imaris_ims
from germquant.validate.imaris_xlsx import per_nucleus_rad51
from scripts.cv_analyze import _ccc, _greedy, _load  # reuse

TOL = 2.5
THRS = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]


def gonad_data(g):
    nuc = _load(g["out_dir"], "*__nuclei.csv")
    spots = _load(g["out_dir"], "*__spots.csv")
    germ = nuc[nuc["in_germline"] == True]  # noqa: E712
    im_spots, im_surf, geom = read_imaris_ims(g["ims"])
    xsurf, _ = per_nucleus_rad51(g["xlsx"])
    emin = geom["ext_min_xyz"]
    x_img = np.c_[xsurf["X"] - emin[0], xsurf["Y"] - emin[1], xsurf["Z"] - emin[2]]
    surf_xyz = im_surf[["x_um", "y_um", "z_um"]].to_numpy()
    coloc = np.full(len(im_surf), np.nan)
    for xi, sj in _greedy(x_img, surf_xyz, TOL):
        coloc[sj] = xsurf["RAD51"].to_numpy()[xi]
    our_xyz = germ[["centroid_x_um", "centroid_y_um", "centroid_z_um"]].to_numpy()
    our_ids = germ["nucleus_id"].to_numpy()
    matched = []  # (our_nucleus_id, coloc_count)
    for sj, oj in _greedy(surf_xyz, our_xyz, TOL):
        if not np.isnan(coloc[sj]):
            matched.append((int(our_ids[oj]), coloc[sj]))
    return spots, matched


import argparse

_ap = argparse.ArgumentParser()
_ap.add_argument("--manifest", default="cv_manifest.json")
gonads = json.load(open(_ap.parse_args().manifest))
data = {g["name"]: gonad_data(g) for g in gonads}

# effect-size distribution per gonad
print("=== effect-size distribution per gonad (on our spots) ===")
for name, (spots, matched) in data.items():
    es = spots["effect_size"].dropna()
    short = name.replace("20251105_N2_nohs_", "").strip()
    print(f"  {short:8s} n_spots={len(spots):6d} es_pop={len(es):6d} "
          f"median={es.median():.2f} p90={es.quantile(.9):.2f} max={es.max():.1f}"
          if len(es) else f"  {short}: no effect_size")

print("\n=== effect_size_min sweep: matched per-nucleus mean per gonad (coloc in []) ===")
hdr = "  thr  " + "".join(f"{n.replace('20251105_N2_nohs_','').strip():>10}" for n in data)
print(hdr)
colocs = {n: np.mean([c for _, c in m]) for n, (_, m) in data.items()}
print("  coloc " + "".join(f"{colocs[n]:>10.2f}" for n in data))
agg = []
for thr in THRS:
    line, all_o, all_c = [], [], []
    for name, (spots, matched) in data.items():
        f = spots[spots["effect_size"].fillna(np.inf) >= thr]  # match production (keep NaN-ES)
        per = f.groupby("nucleus_id").size()
        om = np.mean([per.get(i, 0) for i, _ in matched]) if matched else np.nan
        line.append(om)
        for i, c in matched:
            all_o.append(per.get(i, 0))
            all_c.append(c)
    ccc = _ccc(all_o, all_c)
    bias = np.mean(all_o) - np.mean(all_c)
    agg.append((thr, ccc, bias))
    print(f"  {thr:>4.0f}  " + "".join(f"{v:>10.2f}" for v in line))

print("\n=== aggregate paired agreement vs coloc at each threshold ===")
print("  thr    CCC    bias(ours-coloc)")
for thr, ccc, bias in agg:
    print(f"  {thr:>4.0f}  {ccc:>6.3f}  {bias:>+8.2f}")
best = max(agg, key=lambda t: t[1])
print(f"\nbest CCC at effect_size_min={best[0]} (CCC={best[1]:.3f}, bias={best[2]:+.2f})")
