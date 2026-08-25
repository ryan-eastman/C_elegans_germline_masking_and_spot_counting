"""Build session_covariates.csv (per-gonad SYP-3 image contrast + acquisition + PC values), the table
behind the session-confound figure. Run after collect_tables.py so the PC columns are current."""
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
import chload

CA = r"C:/Users/ryane/coloc_analysis"
rows = []
for f in sorted(glob.glob(f"{CA}/staging/cache/*ccw77*.npz")):
    iid = os.path.basename(f)[:-4]
    z = np.load(f)
    syp, lam = z["syp"], z["lam"] > 0
    nuc, cyt, floor = np.percentile(syp[lam], 50), np.percentile(syp[~lam], 50), np.percentile(syp, 1)
    p = chload.parse_iid(iid)
    rows.append({"image_id": iid, "batch": p["batch"], "sex": p["sex"], "treat": p["treat"], "short": p["short"],
                 "excl": chload.is_excluded(iid), "syp_nuc_p50": float(nuc), "syp_cyto_p50": float(cyt),
                 "syp_contrast": float((nuc - floor) / max(cyt - floor, 1))})
d = pd.DataFrame(rows)
m = pd.read_csv(f"{CA}/acquisition_metadata.csv")[["iid", "exp_545", "pow_545", "exp_477"]]
pc = pd.read_csv(f"{CA}/pc_lamin_vs_dapi.csv")[["image_id", "lamin_PC_dm", "lamin_zsh", "lamin_PC_specific"]]
d = d.merge(m, left_on="image_id", right_on="iid").merge(pc, on="image_id").drop(columns="iid")
d.to_csv(f"{CA}/session_covariates.csv", index=False)
print(f"session_covariates.csv: {len(d)} gonads")
