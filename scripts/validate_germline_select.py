"""Sanity-check germline isolation on real per-nucleus tables. A good method DROPS the off-gonad junk
while KEEPING the spot-bearing nuclei (RAD-51 spots mark germline cells — an independent, non-circular
check, since selection never reads the spot count).

    python scripts/validate_germline_select.py --glob "results_cv_hs/*/*__nuclei.csv"
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

from germquant.germline.select import select_germline


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", default="results_cv_hs/*/*__nuclei.csv",
                    help="glob for per-nucleus CSVs to evaluate")
    files = sorted(glob.glob(ap.parse_args().glob))
    print(f"{len(files)} gonad tables\n")
    for method in ("syp_seeded_cc", "multi_cc", "largest_cc"):
        recalls, drops = [], []
        print(f"=== {method} ===")
        print(f"{'gonad':<14} {'n':>4} {'kept':>4} {'drop%':>5} {'spotRecall':>10}")
        for f in files:
            df = pd.read_csv(f)
            if "n_spots" not in df.columns:
                continue
            out, _ = select_germline(df.drop(columns=["in_germline"], errors="ignore"), method=method)
            m = out["in_germline"].to_numpy(bool)
            spots = df["n_spots"].fillna(0).to_numpy(float)
            sp = spots > 0
            recall = (sp & m).sum() / max(sp.sum(), 1)          # spot-bearing nuclei kept
            drop_pct = 100 * (~m).mean()
            stem = os.path.basename(f).replace("__nuclei.csv", "")[-14:]
            recalls.append(recall); drops.append(drop_pct)
            print(f"{stem:<14} {len(df):>4} {int(m.sum()):>4} {drop_pct:>4.0f}% {recall:>10.3f}")
        if recalls:
            print(f"  --> mean spot-recall={np.mean(recalls):.3f}  min={np.min(recalls):.3f}  "
                  f"mean drop={np.mean(drops):.0f}%\n")


if __name__ == "__main__":
    main()
