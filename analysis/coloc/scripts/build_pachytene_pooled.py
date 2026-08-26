"""pc_pachytene_pooled.csv: the pooled hand-traced pachytene PC per clean-13 gonad (pach_* columns of
pc_zone_all.csv, pc_zone_worker.py) next to the whole-gonad lamin PC (pc_clean13.csv), one row per gonad.
Read by build_filter_compare.py. Run after collect_tables.py and build_clean_table.py."""
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
import chload

CA = r"C:/Users/ryane/coloc_analysis"
vox_um3 = float(np.prod(chload.SP))
z = pd.read_csv(f"{CA}/pc_zone_all.csv")
c = pd.read_csv(f"{CA}/pc_clean13.csv")[["image_id", "lamin_PC_dm", "lamin_zsh", "lamin_PC_specific"]]
d = z.merge(c, on="image_id", how="inner")
d["short"] = d.image_id.str.replace(r"^.*_(LMN1|lmn1)_", "", regex=True)
d["pach_vox_um3"] = (d.pach_vox * vox_um3).round().astype(int)
d = d[["short", "sex", "treat", "batch", "pach_PC", "pach_zsh", "pach_PCspec", "pach_vox_um3",
       "lamin_PC_dm", "lamin_zsh", "lamin_PC_specific"]].sort_values(["sex", "treat", "short"])
d.to_csv(f"{CA}/pc_pachytene_pooled.csv", index=False)
print(f"wrote pc_pachytene_pooled.csv ({len(d)} gonads)")
