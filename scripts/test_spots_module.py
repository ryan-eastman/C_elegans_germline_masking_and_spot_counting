"""Exercise germquant.spots.detect_spots on the cached crop and sanity-check the tidy tables."""
import os

import numpy as np

from germquant.spots import detect_spots

SP = (0.2, 0.108333, 0.108333)
cache = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "smoke_crop.npz")
d = np.load(cache)
rad, lab = d["rad"], d["lab"].astype("int32")

per_spot, per_nuc = detect_spots(rad, lab, SP, marker="RAD-51", effect_size_min=0.0)
n_nuc = len(per_nuc)
print(f"per_spot: {len(per_spot)} rows, cols={list(per_spot.columns)}")
print(f"per_nuc : {n_nuc} nuclei (one row each), cols={list(per_nuc.columns)}")
print(f"  total spots = {int(per_nuc['n_spots'].sum())}  mean/nucleus = {per_nuc['n_spots'].mean():.2f}")
print(f"  nuclei with >=1 spot = {(per_nuc['n_spots'] > 0).sum()}/{n_nuc}")
es = per_spot["effect_size"]
print(f"  effect_size: {es.notna().sum()}/{len(es)} populated; "
      f"range [{np.nanmin(es):.2f}, {np.nanmax(es):.2f}]" if es.notna().any() else "  effect_size: none")
print("\nper_nuc head:")
print(per_nuc.head().to_string(index=False))
print("\nper_spot head:")
print(per_spot.head().to_string(index=False))
