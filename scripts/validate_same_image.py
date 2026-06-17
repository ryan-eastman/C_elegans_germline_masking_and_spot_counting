"""Same-image validation: our pipeline for 20251105_N2_nohs_HERM_001 vs the Imaris .ims ground
truth for the EXACT same image.

IMPORTANT: Imaris masked only a SUBSET of nuclei (146) — "enough to get usable data", NOT every
nucleus. So nucleus COUNT and precision/F1 are meaningless here (our unmatched nuclei are real
nuclei Imaris simply didn't trace, not false positives). The valid comparisons are:
  * RECALL of the Imaris subset: for each Imaris surface, did our segmentation find a nucleus there?
  * PAIRED per-nucleus RAD-51: on matched nuclei, our n_spots vs Imaris spots-per-surface.
  * Per-nucleus RAD-51 DISTRIBUTION: ours vs Imaris (and the canonical xlsx coloc number, 5.84).
"""
import glob
import os

import numpy as np
import pandas as pd

from germquant.validate.imaris_ims import read_imaris_ims, spots_per_surface

OUT = r"C:\Users\ryane\C_elegans_germline_masking_and_spot_counting\results_spotmax_test"
IMS = r"E:\Madeleine\N2\20251105_N2_noHS\20251105_N2_nohs_HERM _001.ims"
TOL_UM = 3.0


def _load(pattern):
    hits = glob.glob(os.path.join(OUT, pattern))
    return pd.read_csv(hits[0]) if hits else None


nuc = _load("*__nuclei.csv")
spots = _load("*__spots.csv")
if nuc is None:
    raise SystemExit(f"no pipeline output in {OUT} yet")

germ = nuc[nuc["in_germline"] == True] if "in_germline" in nuc.columns else nuc  # noqa: E712
im_spots, im_surf, geom = read_imaris_ims(IMS)
im_counts = spots_per_surface(im_spots, im_surf, max_dist_um=TOL_UM)  # per Imaris surface

print("=== SAME-IMAGE VALIDATION: 20251105_N2_nohs_HERM_001 ===")
print(f"voxel (x,y,z) um {geom['voxel_xyz_um'].round(4)}  shape {geom['shape_xyz']}")
print("NOTE: Imaris masked only a SUBSET of nuclei -> counts are NOT expected to match.\n")

print("--- nucleus counts (context only, NOT a match metric) ---")
print(f"  ours all / germline   : {len(nuc)} / {len(germ)}")
print(f"  Imaris (partial subset): {len(im_surf)}")

# Match each Imaris surface to the nearest of OUR germline nuclei (greedy, unique).
pred_xyz = germ[["centroid_x_um", "centroid_y_um", "centroid_z_um"]].to_numpy()
gt_xyz = im_surf[["x_um", "y_um", "z_um"]].to_numpy()
pred_n_spots = germ["n_spots"].to_numpy() if "n_spots" in germ.columns else np.full(len(germ), np.nan)

used = np.zeros(len(pred_xyz), dtype=bool)
pairs = []  # (imaris_count, our_n_spots)
for i, g in enumerate(gt_xyz):
    if len(pred_xyz) == 0:
        break
    d = np.linalg.norm(pred_xyz - g, axis=1)
    d[used] = np.inf
    j = int(d.argmin())
    if d[j] <= TOL_UM:
        used[j] = True
        pairs.append((im_counts[i], pred_n_spots[j]))

n_match = len(pairs)
recall = n_match / max(len(im_surf), 1)
print("\n--- RECALL of Imaris masked subset (did we find their nuclei?) ---")
print(f"  matched {n_match}/{len(im_surf)}  recall={recall:.3f}  (within {TOL_UM} um)")

print("\n--- per-nucleus RAD-51 DISTRIBUTION ---")
if "n_spots" in germ.columns:
    ns = germ["n_spots"].dropna()
    print(f"  ours (all germline)    : mean={ns.mean():.2f} median={ns.median():.0f} (>=1: {(ns > 0).mean():.0%})")
print(f"  Imaris .ims (centroid) : mean={im_counts.mean():.2f} median={np.median(im_counts):.0f} "
      f"(>=1: {(im_counts > 0).mean():.0%})")
print(f"  Imaris coloc (canonical lab method): this gonad HERM_1 = 5.81/nuc; N2 pooled = 5.84/nuc")

if n_match >= 3:
    a = np.array([p[0] for p in pairs], float)  # imaris
    b = np.array([p[1] for p in pairs], float)  # ours
    r = np.corrcoef(a, b)[0, 1] if np.std(a) > 0 and np.std(b) > 0 else float("nan")
    print("\n--- PAIRED per-nucleus RAD-51 (matched nuclei only) ---")
    print(f"  n_pairs={n_match}  ours mean={b.mean():.2f}  Imaris mean={a.mean():.2f}  "
          f"bias(ours-Imaris)={b.mean() - a.mean():+.2f}  Pearson r={r:.3f}")
