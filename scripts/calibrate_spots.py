"""Calibration analysis for the SpotMAX over-detection vs Imaris on 20251105_N2_nohs_HERM_001.

Investigates four things on the SAME image:
  (A) Imaris per-nucleus count under better assignment (sphere-inside vs nearest-centroid).
  (B) Canonical Imaris coloc number for the matching gonad (xlsx HERM_1/2/3 — naming ambiguous).
  (C) effect_size_min sweep on OUR spots: where does per-nucleus count match Imaris?
  (D) Spot-splitting check: are we detecting multiple peaks per real focus?
"""
import glob
import os

import numpy as np
import pandas as pd

from germquant.validate.imaris_ims import read_imaris_ims

OUT = r"C:\Users\ryane\C_elegans_germline_masking_and_spot_counting\results_spotmax_test"
IMS = r"E:\Madeleine\N2\20251105_N2_noHS\20251105_N2_nohs_HERM _001.ims"
TOL = 3.0


def _load(p):
    h = glob.glob(os.path.join(OUT, p))
    return pd.read_csv(h[0]) if h else None


nuc = _load("*__nuclei.csv")
spots = _load("*__spots.csv")
germ = nuc[nuc["in_germline"] == True].copy()  # noqa: E712
im_spots, im_surf, geom = read_imaris_ims(IMS)
surf_xyz = im_surf[["x_um", "y_um", "z_um"]].to_numpy()
surf_r = (3 * im_surf["volume_um3"].to_numpy() / (4 * np.pi)) ** (1 / 3)  # sphere-equiv radius
print(f"surfaces n={len(im_surf)}  sphere-equiv radius mean={surf_r.mean():.2f} um")

# (A) Imaris per-surface counts: nearest-centroid<=TOL vs inside-sphere(radius)
sp_xyz = im_spots[["x_um", "y_um", "z_um"]].to_numpy()
cnt_near = np.zeros(len(im_surf), int)
cnt_inside = np.zeros(len(im_surf), int)
for p in sp_xyz:
    d = np.linalg.norm(surf_xyz - p, axis=1)
    j = int(d.argmin())
    if d[j] <= TOL:
        cnt_near[j] += 1
    inside = np.where(d <= surf_r)[0]
    if len(inside):
        cnt_inside[inside[np.argmin(d[inside])]] += 1
print("\n(A) Imaris RAD-51/nucleus on THIS gonad (146 surfaces):")
print(f"   nearest-centroid<= {TOL}um : mean={cnt_near.mean():.2f}  assigned={cnt_near.sum()}/{len(im_spots)}")
print(f"   inside-sphere(r_vol)       : mean={cnt_inside.mean():.2f}  assigned={cnt_inside.sum()}/{len(im_spots)}")

# (B) canonical xlsx for the matching gonad
from germquant.validate.imaris_xlsx import summarize
xdir = r"E:\Madeleine\202511_imaris_export_het_mutants_Crest"
print("\n(B) Imaris xlsx coloc (canonical) for 20251105 N2 noHS HERM:")
for n in (1, 2, 3):
    f = os.path.join(xdir, f"20251105_N2_nohs_HERM_{n}.xlsx")
    if os.path.exists(f):
        s = summarize(f)
        print(f"   HERM_{n}: nuclei={s['n_nuclei']}  RAD51/nuc mean={s['rad51_per_nucleus_mean']:.2f}")

# match each Imaris surface -> our nearest germline nucleus (for paired sweep)
pred_xyz = germ[["centroid_x_um", "centroid_y_um", "centroid_z_um"]].to_numpy()
match_our_idx = []
used = np.zeros(len(pred_xyz), bool)
for g in surf_xyz:
    d = np.linalg.norm(pred_xyz - g, axis=1)
    d[used] = np.inf
    j = int(d.argmin())
    if d[j] <= TOL:
        used[j] = True
        match_our_idx.append(germ.iloc[j]["nucleus_id"])
matched_ids = match_our_idx  # list, NOT set: 1:1 here, but a set would silently dedupe collisions

# (C) effect_size_min sweep on OUR spots
print("\n(C) effect_size_min sweep (our spots; paired = matched-to-Imaris nuclei):")
print(f"   {'thr':>5} {'n_spots':>8} {'all_germ_mean':>13} {'paired_mean':>11}  (Imaris inside-sphere ~"
      f"{cnt_inside.mean():.2f}, HERM_1 coloc 5.81)")
for thr in (0.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0):
    s = spots[spots["effect_size"] >= thr]
    per = s.groupby("nucleus_id").size()
    all_mean = germ["nucleus_id"].map(per).fillna(0).mean()
    paired_mean = np.mean([per.get(i, 0) for i in matched_ids]) if matched_ids else float("nan")
    print(f"   {thr:>5} {len(s):>8} {all_mean:>13.2f} {paired_mean:>11.2f}")

# (D) Z-AXIS spot-splitting: do NOT use a 3D NN <0.3um test — it is blind to z-splitting because the
# z-voxel is 0.2um (a 1-2 voxel z-pair is 0.2-0.4um, above 0.3um). Test per-(x,y)-column stacking.
DX = DY = 0.108333
print("\n(D) z-axis spot-splitting check (per-(x,y)-voxel column):")
inn = spots[spots["nucleus_id"] > 0].copy()
inn["xp"] = np.round(inn["x_um"].to_numpy() / DX).astype(int)
inn["yp"] = np.round(inn["y_um"].to_numpy() / DY).astype(int)
col_sizes = inn.groupby(["xp", "yp"]).size()
multi = col_sizes[col_sizes > 1]
frac_stacked = int(multi.sum() - len(multi)) / max(len(inn), 1)  # spots beyond 1-per-column
print(f"   in-nucleus spots={len(inn)}  unique (x,y) columns={len(col_sizes)}  "
      f"multi-spot columns={len(multi)}")
print(f"   spots that are z-stacked duplicates (col size>1, beyond first): {frac_stacked:.1%}")
# collapse each (x,y) column to one spot -> compare to Imaris total
collapsed = len(col_sizes) + int((spots['nucleus_id'] <= 0).sum())  # +any non-nucleus spots kept
print(f"   total spots: raw ours={len(spots)}  column-collapsed~={collapsed}  "
      f"Imaris={len(im_spots)}  (raw ratio={len(spots)/len(im_spots):.2f})")
# per-axis NN within nuclei (z separated from xy) to show the splitting is along z
znn, xynn = [], []
for nid, g in inn.groupby("nucleus_id"):
    if len(g) < 2:
        continue
    xy = g[["x_um", "y_um"]].to_numpy()
    z = g["z_um"].to_numpy()
    for k in range(len(g)):
        dxy = np.hypot(xy[:, 0] - xy[k, 0], xy[:, 1] - xy[k, 1])
        dxy[k] = np.inf
        j = int(dxy.argmin())
        xynn.append(dxy[j])
        znn.append(abs(z[j] - z[k]))
if znn:
    print(f"   nearest-neighbour: median xy-dist={np.median(xynn):.3f}um  median |z-dist|={np.median(znn):.3f}um  "
          f"(splitting is along z if z-dist << xy-dist)")
