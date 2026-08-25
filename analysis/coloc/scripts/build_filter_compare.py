"""Assemble pc_filter_compare.csv: whole-gonad lamin-masked PC under three nucleus sets (all in_germline
labels / no-envelope objects removed / plus off-trace territories removed; pc_filter_compare.py rows)
next to the pooled hand-traced pachytene PC (pc_zone_worker.py, pach_* columns), clean-13 gonads, and
print the male / hermaphrodite HS-vs-noHS comparison for each version."""
import glob
import json

import numpy as np
import pandas as pd
from scipy import stats

CA = r"C:/Users/ryane/coloc_analysis"
rows = []
for f in sorted(glob.glob(f"{CA}/pc_lamin_rows_filtered/*.json")):
    r = json.load(open(f))
    d = {"image_id": r["image_id"], "short": r["image_id"].split("LMN1_")[-1].split("lmn1_")[-1], "sex": r["sex"],
         "treat": r["treat"], "batch": r["batch"], "n_all": r["n_all"], "n_no_ring": r["n_no_ring"],
         "n_territories": r["n_territories"], "n_kept": r["n_ring_terr"], "territory_note": r["territory_note"]}
    for k in ["all", "ring", "ring_terr"]:
        m = r[k] or {}
        d[f"{k}_PC"] = m.get("PC_dm")
        d[f"{k}_zsh"] = m.get("zsh")
        d[f"{k}_PCspec"] = m.get("PC_specific")
        d[f"{k}_ngran"] = m.get("n_gran")
    rows.append(d)
d = pd.DataFrame(rows).sort_values(["sex", "treat", "short"])
pool = pd.read_csv(f"{CA}/pc_pachytene_pooled.csv")[["short", "pach_PC", "pach_zsh", "pach_PCspec", "pach_vox_um3"]]
d = d.merge(pool, on="short", how="left")
d.to_csv(f"{CA}/pc_filter_compare.csv", index=False)
pd.set_option("display.width", 250)
print(d[["short", "n_all", "n_no_ring", "n_territories", "n_kept", "all_PCspec", "ring_PCspec", "ring_terr_PCspec", "pach_PCspec",
         "all_PC", "ring_PC", "ring_terr_PC", "pach_PC"]].to_string(index=False))
print(f"\n{len(d)} gonads")
for col in ["all_PCspec", "ring_PCspec", "ring_terr_PCspec", "pach_PCspec", "all_PC", "ring_PC", "ring_terr_PC", "pach_PC"]:
    out = []
    for sex in ["male", "herm"]:
        a = d[(d.sex == sex) & (d.treat == "noHS")][col].dropna()
        b = d[(d.sex == sex) & (d.treat == "HS")][col].dropna()
        p = stats.mannwhitneyu(a, b)[1] if len(a) > 1 and len(b) > 1 else np.nan
        out.append(f"{sex} noHS {a.mean():.3f} (n{len(a)}) HS {b.mean():.3f} (n{len(b)}) d {b.mean() - a.mean():+.3f} P {p:.3f}")
    print(f"{col:18s} | " + " | ".join(out))
