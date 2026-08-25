"""Build acquisition_metadata.csv from the nd2 text metadata: per-channel exposure (ms) and laser power
for every ccw77 gonad on the local drive, plus the audit exclusion flag."""
import glob
import os
import re
import sys

import nd2
import pandas as pd

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
import chload

CA = r"C:/Users/ryane/coloc_analysis"
FOLDERS = [r"E:/20260708_ccw77_IF_pgl1_syp3_LMN1", r"E:/20260622_ccw77_IF_syp3_pgl1_LMN1_dapi"]
rows = []
for folder in FOLDERS:
    for f in sorted(glob.glob(folder + "/*.nd2")):
        iid = os.path.basename(f)[:-4]
        p = chload.parse_iid(iid)
        with nd2.ND2File(f) as h:
            desc = h.text_info.get("description", "")
        r = {"iid": iid, "batch": p["batch"], "sex": p["sex"], "treat": p["treat"], "excl": chload.is_excluded(iid)}
        for pl in re.split(r"Plane #\d+:", desc)[1:]:
            name = re.search(r"Name:\s*(\S+)", pl).group(1)
            exp = re.search(r"Exposure:\s*([\d.]+)\s*ms", pl)
            pw = re.search(r"ExW:%s; Power:\s*([\d.]+)" % name, pl)
            r[f"exp_{name}"] = float(exp.group(1)) if exp else None
            r[f"pow_{name}"] = float(pw.group(1)) if pw else None
        rows.append(r)
pd.DataFrame(rows).to_csv(f"{CA}/acquisition_metadata.csv", index=False)
print(f"acquisition_metadata.csv: {len(rows)} gonads")
