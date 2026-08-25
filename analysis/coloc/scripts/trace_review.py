"""Review montage: every gonad's mid-z MIP with the drawn pachytene trace overlaid, 5x5 grid.
Green line = trace, dot = START. Panels flagged [EXCLUDED] are not in the clean-13 PC set anyway."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STAGING = r"C:/Users/ryane/coloc_analysis/staging"
PXY, DS = 0.1083, 4
import sys
sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
import chload   # exclusions come from coloc_analysis/exclusions.json via chload.is_excluded
tr = json.load(open(os.path.join(STAGING, "pachytene_traces.json"), encoding="utf-8"))
iids = sorted(tr.keys())
import math
nc_ = 5
nr_ = max(1, math.ceil(len(iids) / nc_))
fig, axes = plt.subplots(nr_, nc_, figsize=(22, 4.4 * nr_), facecolor="black", squeeze=False)
for ax, iid in zip(axes.ravel(), iids):
    z = np.load(os.path.join(STAGING, "cache", f"{iid}.npz"))
    D = z["dapi"][:, ::DS, ::DS]
    nz = D.shape[0]
    d = D[nz // 4:max(nz // 4 + 1, 3 * nz // 4)].max(0).astype(np.float32)
    p0, p1 = np.percentile(d, [1, 99.6])
    dn = np.clip((d - p0) / (p1 - p0 + 1e-9), 0, 1)
    h, w = dn.shape
    ax.imshow(dn, cmap="gray", extent=[0, w * PXY * DS, h * PXY * DS, 0])
    P = tr[iid].get("points_um", [])
    if len(P) >= 2:
        xs, ys = zip(*P)
        ax.plot(xs, ys, "-", color="lime", lw=1.6)
        ax.plot(xs[0], ys[0], "o", color="lime", ms=6)
        ax.plot(xs[-1], ys[-1], "s", color="magenta", ms=6)
    short = iid.split("LMN1_")[-1] if "LMN1_" in iid else iid.split("lmn1_")[-1] if "lmn1_" in iid else iid.split("lmn_")[-1]
    is_n2 = "N2" in iid
    tag = " [N2]" if is_n2 else ""
    excl = chload.is_excluded(iid)
    ax.set_title(f"{short}{tag}" + (" [EXCLUDED from PC]" if excl else ""), fontsize=9,
                 color="#f66" if excl else "w")
    ax.axis("off")
for ax in axes.ravel()[len(iids):]:
    ax.axis("off")
fig.suptitle("Pachytene traces: review grid (green=trace, circle=START, square=END; red titles = already excluded from the PC analysis)",
             color="w", fontsize=13)
fig.tight_layout()
out = os.path.join(STAGING, "trace_review_grid.png")
fig.savefig(out, dpi=110, facecolor="black")
print("WROTE", out)
