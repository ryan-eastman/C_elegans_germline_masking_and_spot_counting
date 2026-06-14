"""Does the chromatin-polarization signal exist on REAL data, and where? Bin each real gonad's
nuclei by axis position and show the mean polarization profile + zone calls. If polarization is
elevated at one end (the transition zone), the per-nucleus metric works and the failure is the
threshold/aggregation; if it's flat, the metric itself is too weak on real data."""
import glob
import os

import numpy as np
import pandas as pd

GLOB = "results_nas/*/*__nuclei.csv"
NB = 10


def main():
    files = sorted(glob.glob(GLOB))
    print(f"{len(files)} real gonads; polarization profile along axis (10 distal->proximal bins)\n")
    for f in files:
        df = pd.read_csv(f)
        if "chromatin_polarization" not in df or "axis_position_norm" not in df:
            continue
        d = df.dropna(subset=["axis_position_norm", "chromatin_polarization"])
        if len(d) < 30:
            continue
        b = np.clip((d["axis_position_norm"] * NB).astype(int), 0, NB - 1)
        stem = os.path.basename(f).split("syp1_")[-1].replace("__nuclei.csv", "")
        pol = d.groupby(b)["chromatin_polarization"].mean()
        vol = d.groupby(b)["volume_um3"].median() if "volume_um3" in d else None
        cnt = d.groupby(b).size()
        polline = " ".join(f"{pol.get(i, float('nan')):.2f}" for i in range(NB))
        print(f"{stem:10s} pol =[{polline}]")
        if vol is not None:
            volline = " ".join(f"{vol.get(i, float('nan')):.0f}" for i in range(NB))
            print(f"{'':10s} vol =[{volline}]  (um3 median; smaller distal = mitotic/TZ?)")
        cntline = " ".join(f"{cnt.get(i, 0):3d}" for i in range(NB))
        print(f"{'':10s} n   =[{cntline}]  (count per bin; denser distal = mitotic/TZ?)")
    print("\n(distal = bin 0. Looking for any feature elevated/depressed at one end = the TZ/mitotic band)")


if __name__ == "__main__":
    main()
