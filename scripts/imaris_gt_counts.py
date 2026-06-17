"""Compute the Imaris ground-truth RAD-51-per-nucleus distribution for the N2 control gonads,
using the lab's exact coloc-matching. This is the number to calibrate our pipeline's count against.
"""
import glob
import os

import pandas as pd

from germquant.validate.imaris_xlsx import summarize

EXPORT_DIR = r"E:\Madeleine\202511_imaris_export_het_mutants_Crest"
files = sorted(glob.glob(os.path.join(EXPORT_DIR, "*n2_nohs*.xlsx")))
print(f"N2 no-HS Imaris exports found: {len(files)}")

rows = []
for f in files:
    try:
        rows.append(summarize(f))
        print(f"  OK {os.path.basename(f)}")
    except Exception as e:  # noqa: BLE001
        print(f"  FAIL {os.path.basename(f)}: {type(e).__name__}: {e}")

if rows:
    df = pd.DataFrame(rows)
    print("\n=== per-gonad ===")
    print(df.to_string(index=False))
    print("\n=== N2 no-HS RAD-51 per nucleus (ground truth) ===")
    print(f"  gonads             : {len(df)}")
    print(f"  mean of per-gonad means : {df['rad51_per_nucleus_mean'].mean():.2f}")
    print(f"  median of per-gonad     : {df['rad51_per_nucleus_median'].median():.2f}")
    print(f"  total spots / total nuclei: "
          f"{df['n_spots_assigned'].sum() / max(df['n_nuclei'].sum(), 1):.2f}")
    print(f"  mean frac nuclei w/ >=1 : {df['frac_nuclei_with_rad51'].mean():.2f}")
