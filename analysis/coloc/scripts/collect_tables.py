"""Rebuild the aggregate tables from the per-gonad JSON rows (the workers write one JSON per gonad;
this is the ONLY place that assembles them, so a rerun of any worker is followed by this script).
  pc_lamin_rows/*.json  ->  pc_lamin_vs_dapi.csv   (then build_clean_table.py -> pc_clean13.csv)
  pc_zone_rows/*.json   ->  pc_zone_all.csv
Rows whose gonad is excluded by exclusions.json, or whose trace is no longer 'traced', are dropped
here rather than silently carried forward."""
import glob
import json
import os
import sys

import pandas as pd

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
import chload

CA = r"C:/Users/ryane/coloc_analysis"

# whole-gonad lamin/DAPI rows
rows = []
for f in sorted(glob.glob(f"{CA}/pc_lamin_rows/*.json")):
    r = json.load(open(f))
    if "batch" not in r:      # rows without the batch key predate the rotation-null fix: never mix versions
        raise SystemExit(f"{os.path.basename(f)} is from an older worker version; rerun the lane before collecting")
    p = chload.parse_iid(r["image_id"])   # all 22 kept here; build_clean_table applies exclusions.json
    row = {"image_id": r["image_id"], "sex": p["sex"], "treat": p["treat"], "batch": p["batch"],
           "n_nuclei": r["n_nuclei"], "fallback_n": r["lamin_fallback_n"], "vol_ratio": r["lam_vol_ratio"]}
    for mk in ["lamin", "dapi"]:
        m = r.get(mk) or {}
        for k in ["n_gran", "PC_dm", "rot", "zsh", "PC_specific"]:
            row[f"{mk}_{k}"] = m.get(k)
    rows.append(row)
lam = pd.DataFrame(rows).sort_values(["sex", "treat", "image_id"])
lam.to_csv(f"{CA}/pc_lamin_vs_dapi.csv", index=False)
print(f"pc_lamin_vs_dapi.csv: {len(lam)} gonads")

# zone rows (traced + not excluded + ccw77 only)
traces = json.load(open(f"{CA}/staging/pachytene_traces.json", encoding="utf-8"))
rows = []
for f in sorted(glob.glob(f"{CA}/pc_zone_rows/*.json")):
    r = json.load(open(f))
    iid = r["image_id"]
    if traces.get(iid, {}).get("status") != "traced" or chload.is_excluded(iid) or "ccw77" not in iid:
        print(f"  dropping zone row {iid[-22:]} (trace status / exclusion / strain)")
        continue
    p = chload.parse_iid(iid)
    r.update(sex=p["sex"], treat=p["treat"], batch=p["batch"])
    rows.append(r)
zone = pd.DataFrame(rows).sort_values(["sex", "treat", "image_id"])
zone.to_csv(f"{CA}/pc_zone_all.csv", index=False)
print(f"pc_zone_all.csv: {len(zone)} gonads")
