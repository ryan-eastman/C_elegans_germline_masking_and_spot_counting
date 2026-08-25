"""Re-assign zones for every 'traced' gonad from the saved polylines (no GUI). Run after any change to
assign_zones (e.g. the adaptive off-axis cutoff) so zones/*.csv reflect the current rule."""
import json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import trace_pachytene as tp
tr = tp.load_traces()
for iid, t in sorted(tr.items()):
    if t["status"] != "traced" or len(t["points_um"]) < 2:
        continue
    df = pd.read_csv(os.path.join(tp.STAGING, f"{iid}_nuclei.csv"))
    z, L = tp.assign_zones(df, np.asarray(t["points_um"]))
    z.to_csv(os.path.join(tp.ZONE_DIR, f"{iid}_zones.csv"), index=False)
    n = z.zone.value_counts()
    print(f"{iid[-20:]:22s} cut={z.off_axis_cut_um.iloc[0]:5.1f}  early/mid/late={n.get('early',0)}/{n.get('mid',0)}/{n.get('late',0)}  off_axis={n.get('off_axis',0)}")
