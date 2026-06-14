"""Pick the germline-isolation default empirically on the 14 real per-nucleus tables. A good method
DROPS a meaningful fraction (the off-gonad junk) while KEEPING the foci-bearing nuclei (foci mark
germline cells — an independent, non-circular check since selection never reads n_foci).

    python scripts/validate_germline_select.py
"""
import glob
import os

import numpy as np
import pandas as pd

from germquant.germline.select import select_germline

GLOB = "results_nas/*/*__nuclei.csv"


def main():
    files = sorted(glob.glob(GLOB))
    print(f"{len(files)} real gonad tables\n")
    for method in ("syp_seeded_cc", "multi_cc", "largest_cc"):
        recalls, drops, kepts = [], [], []
        print(f"=== {method} ===")
        print(f"{'gonad':<10} {'n':>4} {'kept':>4} {'drop%':>5} {'fociRecall':>10} {'droppedFoci':>11}")
        for f in files:
            df = pd.read_csv(f)
            out, _ = select_germline(df, method=method)
            m = out["in_germline"].to_numpy(bool)
            foci = df["n_foci"].to_numpy(float)
            fp = foci > 0
            recall = (fp & m).sum() / max(fp.sum(), 1)       # foci-bearing nuclei kept
            dropf = foci[~m].sum()                            # total foci thrown away
            drop_pct = 100 * (~m).mean()
            stem = os.path.basename(f).split("syp1_")[-1].replace("__nuclei.csv", "")
            recalls.append(recall); drops.append(drop_pct); kepts.append(int(m.sum()))
            print(f"{stem:<10} {len(df):>4} {int(m.sum()):>4} {drop_pct:>4.0f}% {recall:>10.3f} {dropf:>11.0f}")
        print(f"  --> mean foci-recall={np.mean(recalls):.3f}  min={np.min(recalls):.3f}  "
              f"mean drop={np.mean(drops):.0f}%\n")


if __name__ == "__main__":
    main()
