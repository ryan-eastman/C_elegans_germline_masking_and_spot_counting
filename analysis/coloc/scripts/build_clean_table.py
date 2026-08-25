"""Build pc_clean13.csv from pc_lamin_vs_dapi.csv using exclusions.json (single source of truth).
Run after any pc_lamin_worker rerun and before any publication figure (fig_pub_partition reads it)."""
import sys

import pandas as pd

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
import chload

CA = r"C:/Users/ryane/coloc_analysis"
d = pd.read_csv(f"{CA}/pc_lamin_vs_dapi.csv")
d["excl"] = [chload.is_excluded(i) for i in d.image_id]
clean = d[~d.excl].drop(columns="excl")
clean.to_csv(f"{CA}/pc_clean13.csv", index=False)
print(f"clean table: {len(clean)} of {len(d)} gonads -> {CA}/pc_clean13.csv")
print("excluded:", sorted(chload.parse_iid(i)["short"] for i in d.image_id[d.excl]))
