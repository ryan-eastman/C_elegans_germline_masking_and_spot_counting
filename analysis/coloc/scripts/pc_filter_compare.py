"""Whole-gonad lamin-masked PC under three nucleus sets, side by side, WITHOUT touching the staged rows:
  all        = every in_germline label (what pc_lamin_worker.py uses today)
  ring       = minus labels with no lamin envelope           (nucleus_filter.no_ring_ids)
  ring+terr  = ring, minus labels outside the traced territory (nucleus_filter.territory_ids)
Output: pc_lamin_rows_filtered/<iid>.json. Usage: python pc_filter_compare.py <image_id>"""
import json
import os
import sys

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import chload
import nucleus_filter as NF
import pc_lamin_worker as W

OUTDIR = r"C:/Users/ryane/coloc_analysis/pc_lamin_rows_filtered"
os.makedirs(OUTDIR, exist_ok=True)


def main(iid):
    rd = W.find_run(iid)
    lab = tifffile.imread(os.path.join(rd, f"{iid}__nuclei_labels.tif"))
    csv = os.path.join(rd, f"{iid}__nuclei.csv")
    nc = pd.read_csv(csv)
    germ = [int(x) for x in nc[nc.in_germline.astype(bool)].nucleus_id]
    gn = np.isin(lab, germ)
    obj = ndi.find_objects(gn.astype(np.uint8))[0]
    sl = tuple(slice(max(0, o.start - W.PAD), min(dim, o.stop + W.PAD)) for o, dim in zip(obj, gn.shape))
    ch = chload.load(W.find_nd2(iid, csv))
    pgl = ch["pgl"][sl].astype(np.float32)
    syp = ch["syp"][sl].astype(np.float32)
    lamin = ch["lamin"][sl]
    del ch
    lab_c = lab[sl]
    bg = float(np.percentile(syp, 3))          # same crop, same background as the worker

    scores = NF.ring_scores(lab_c, germ, lamin)
    junk = set(NF.no_ring_ids(scores))
    ring = [g for g in germ if g not in junk]
    tr = json.load(open(r"C:/Users/ryane/coloc_analysis/staging/pachytene_traces.json", encoding="utf-8")).get(iid, {})
    if tr.get("status") == "traced":
        kept, n_terr, note = NF.territory_ids(lab_c, ring, tr["points_um"])
    else:
        kept, n_terr, note = ring, None, "no trace: territory filter not applied"
    p = chload.parse_iid(iid)
    res = {"image_id": iid, "sex": p["sex"], "treat": p["treat"], "batch": p["batch"],
           "n_all": len(germ), "n_no_ring": len(junk), "n_ring": len(ring), "n_territories": n_terr,
           "n_ring_terr": len(kept), "territory_note": note, "ring_thr": round(scores.attrs["ring_thr"], 1)}
    for name, ids in [("all", germ), ("ring", ring), ("ring_terr", kept)]:
        nuc_lam, _, fb = W.lamin_nuclei(lab_c, ids, lamin)
        m = W.metric_set(nuc_lam, pgl, syp, bg)
        res[name] = m
        res[f"{name}_fallback_n"] = fb
        print(f"  {name:10s} n_nuc={len(ids):4d} fallback={fb:3d} {m}", flush=True)
    with open(os.path.join(OUTDIR, iid + ".json"), "w") as f:
        json.dump(res, f, indent=1)
    print(f"{p['short']}: no-ring {len(junk)}, territories {n_terr}, kept {len(kept)} of {len(germ)} {note}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
