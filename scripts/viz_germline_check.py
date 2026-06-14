"""High-res check: are the RED (dropped) nuclei genuine off-gonad junk, or did SYP-based isolation
wrongly drop the SYP-negative DISTAL tip (mitotic/TZ) of the gonad? Render the DNA MIP with green
(kept) / red (dropped) nucleus outlines, large, so it can be judged by eye."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from skimage.segmentation import find_boundaries

import germquant.io.nd2_reader as R
from germquant.germline import select_germline

REAL = r"data/raw_examples/madeleline images/20251105_N2_noHS/20251105_N2_nohs_HERM _001.nd2"
IMG = "20251105_N2_nohs_HERM _001"
LABELS = f"results_germline/{IMG}__nuclei_labels.tif"
NUC = f"results_germline/{IMG}__nuclei.csv"


def norm(a):
    lo, hi = np.percentile(a, 1), np.percentile(a, 99.5)
    return np.clip((a - lo) / (hi - lo + 1e-9), 0, 1)


def main():
    labels = tifffile.imread(LABELS)
    df = pd.read_csv(NUC)
    # recompute in_germline with the CURRENT default method (no re-segmentation needed)
    df, _ = select_germline(df.drop(columns=["in_germline"], errors="ignore"))
    keep = df["in_germline"].astype(bool)
    excl = set(df.loc[~keep, "nucleus_id"].astype(int))
    st = R.read_stack(REAL)
    dna = norm(st.data[0].max(0))
    syp = norm(st.data[1].max(0))
    lab_mip = labels.max(0)
    bnd = find_boundaries(lab_mip, mode="inner")
    is_excl = np.isin(lab_mip, list(excl))

    fig, ax = plt.subplots(1, 2, figsize=(26, 13))
    for a, base, t in ((ax[0], dna, "DNA + green=germline / red=off-gonad"),
                       (ax[1], syp, "SYP (central element) — note SYP-negative distal tip")):
        a.imshow(base, cmap="gray")
        ov = np.zeros((*bnd.shape, 4))
        ov[bnd & ~is_excl] = (0, 1, 0, 1)
        ov[bnd & is_excl] = (1, 0, 0, 1)
        a.imshow(ov); a.set_title(t, fontsize=13); a.axis("off")
    fig.tight_layout(); fig.savefig("results_fullres/germline_check.png", dpi=130, bbox_inches="tight")
    print(f"kept={int(keep.sum())} dropped={len(excl)}  wrote results_fullres/germline_check.png")


if __name__ == "__main__":
    main()
