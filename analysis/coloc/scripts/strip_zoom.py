"""Zoomed verification strip for a row range: python strip_zoom.py <iid> <row0> <row1> [per_row]
Reads staging/<iid>_nuclei.csv, reloads DAPI+SYP from E:/NAS, renders big tiles per row."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import chload
import crescent_axis as ca

iid, r0, r1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
PER = int(sys.argv[4]) if len(sys.argv) > 4 else 7
SP = chload.SP
df = pd.read_csv(os.path.join(ca.OUTDIR, f"{iid}_nuclei.csv"))
rd = ca.find_run(iid)
lab = tifffile.imread(os.path.join(rd, f"{iid}__nuclei_labels.tif"))
nc = pd.read_csv(os.path.join(rd, f"{iid}__nuclei.csv"))
germ = [int(x) for x in nc[nc.in_germline.astype(bool)].nucleus_id]
gn = np.isin(lab, germ)
obj = ndi.find_objects(gn.astype(np.uint8))[0]
sl = tuple(slice(max(0, o.start - ca.PAD), min(dim, o.stop + ca.PAD)) for o, dim in zip(obj, gn.shape))
ch = chload.load(ca.find_nd2(iid, os.path.join(rd, f"{iid}__nuclei.csv")))
dapi = ch["dapi"][sl].astype(np.float32); syp = ch["syp"][sl].astype(np.float32)
del ch
R = int(np.ceil(3.4 / SP[1]))
rows = list(range(r0, r1 + 1))
strip = np.zeros((len(rows) * 2 * R, PER * 2 * R, 3), np.float32)
for i, ri in enumerate(rows):
    g = df[df.row == ri].sort_values("s_um").head(PER)
    for k, r in enumerate(g.itertuples()):
        zc, yc, xc = int(r.cz / SP[0]), int(r.cy / SP[1]), int(r.cx / SP[2])
        dd = dapi[zc, max(0, yc - R):yc + R, max(0, xc - R):xc + R]
        ss = syp[zc, max(0, yc - R):yc + R, max(0, xc - R):xc + R]
        if dd.shape != (2 * R, 2 * R):
            continue
        dn = np.clip((dd - np.percentile(dd, 5)) / (np.percentile(dd, 99.5) - np.percentile(dd, 5) + 1e-9), 0, 1)
        sn = np.clip((ss - np.percentile(ss, 30)) / (np.percentile(ss, 99.5) - np.percentile(ss, 30) + 1e-9), 0, 1)
        tile = np.stack([np.maximum(dn, sn * 0.8), dn, dn], -1)
        colr = ca.CMAP.get(r.cls, (1, 1, 1))
        tile[:3, :] = colr; tile[-3:, :] = colr; tile[:, :3] = colr; tile[:, -3:] = colr
        strip[i * 2 * R:(i + 1) * 2 * R, k * 2 * R:(k + 1) * 2 * R] = tile
fig, ax = plt.subplots(figsize=(PER * 1.6, len(rows) * 1.6), facecolor="black")
ax.imshow(strip); ax.axis("off")
for i, ri in enumerate(rows):
    ax.text(-6, i * 2 * R + R, f"row {ri}", color="w", fontsize=10, ha="right", va="center")
ax.set_title(f"{iid[-22:]}  rows {r0}-{r1}  (DAPI gray, SYP-3 red; border orange=auto TZ, blue=auto pachytene)", color="w", fontsize=10)
fig.tight_layout()
out = os.path.join(ca.OUTDIR, f"{iid}_zoom_rows{r0}-{r1}.png")
fig.savefig(out, dpi=100, facecolor="black")
print("WROTE", out)
