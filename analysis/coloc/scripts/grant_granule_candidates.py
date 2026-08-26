"""Contact sheet of candidate P-granule panels for the grant figure: several windows along the
pachytene arm of HS_male_008 x several planes, PGL-1 in green with granule outlines (yellow) and the
envelope (red). Pick one by (window, plane) and pass it to fig_grant_assets.py."""
import os
import sys

import numpy as np
import pandas as pd
import tifffile
from PIL import Image
from scipy import ndimage as ndi
from skimage.segmentation import find_boundaries
from skimage.transform import resize

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import chload
import pc_lamin_worker as W

IID = "20260708_ccw77_IF_pgl1_syp3_lmn1_HS_male_008"
OUT = r"C:/Users/ryane/coloc_analysis/grant/assets"
SP = chload.SP
WIN_UM = 24.0
# candidate centres (y, x) in crop um along the traced arm (early -> mid pachytene)
CENTRES = [(150.0, 146.0), (120.0, 143.0), (95.0, 132.0), (85.0, 122.0), (76.0, 100.0), (73.0, 72.0)]
PLANES_UM = [6.0, 8.0, 10.0]

rd = W.find_run(IID)
lab = tifffile.imread(os.path.join(rd, f"{IID}__nuclei_labels.tif"))
csv = os.path.join(rd, f"{IID}__nuclei.csv")
nc = pd.read_csv(csv)
germ = [int(x) for x in nc[nc.in_germline.astype(bool)].nucleus_id]
gn = np.isin(lab, germ)
obj = ndi.find_objects(gn.astype(np.uint8))[0]
sl = tuple(slice(max(0, o.start - W.PAD), min(dim, o.stop + W.PAD)) for o, dim in zip(obj, gn.shape))
ch = chload.load(W.find_nd2(IID, csv))
pgl = ch["pgl"][sl].astype(np.float32)
lamin = ch["lamin"][sl]
del ch
lab_c = lab[sl]
nuc_lam, _, _ = W.lamin_nuclei(lab_c, germ, lamin)
dt = ndi.distance_transform_edt(~nuc_lam, sampling=tuple(SP))
cyto = ndi.binary_fill_holes(dt <= W.CYTO_UM) & ~nuc_lam
gran, ng = W.segment_granules(pgl, cyto)
half = int(round(WIN_UM / 2 / SP[1]))
tiles = []
for (cy_um, cx_um) in CENTRES:
    row = []
    cy, cx = int(round(cy_um / SP[1])), int(round(cx_um / SP[2]))
    for zu in PLANES_UM:
        z = int(round(zu / SP[0]))
        w = (z, slice(cy - half, cy + half), slice(cx - half, cx + half))
        P = pgl[w]
        lo, hi = np.percentile(P, 1), np.percentile(P, 99.8)
        g = np.clip((P - lo) / (hi - lo + 1e-9), 0, 1)
        rgb = np.stack([g * 0.15, g, g * 0.15], -1)
        rgb = resize(rgb, (rgb.shape[0] * 2, rgb.shape[1] * 2, 3), order=1, preserve_range=True)
        E = resize(nuc_lam[w].astype(float), rgb.shape[:2], order=0) > 0.5
        G = resize(gran[w].astype(float), rgb.shape[:2], order=0) > 0.5
        rgb[ndi.binary_dilation(find_boundaries(E, mode="inner"))] = (1.0, 0.2, 0.2)
        rgb[ndi.binary_dilation(find_boundaries(G, mode="inner"))] = (1.0, 0.9, 0.1)
        n_gr = int(ndi.label(gran[w])[1])
        row.append((rgb, n_gr))
    tiles.append(row)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(len(CENTRES), len(PLANES_UM), figsize=(3.2 * len(PLANES_UM), 3.2 * len(CENTRES)), facecolor="black")
for i, row in enumerate(tiles):
    for j, (rgb, n_gr) in enumerate(row):
        a = ax[i, j]
        a.imshow(np.clip(rgb, 0, 1))
        a.axis("off")
        a.set_title(f"win {i} (y{CENTRES[i][0]:.0f},x{CENTRES[i][1]:.0f})  z {PLANES_UM[j]:.0f} um  granules {n_gr}", color="w", fontsize=8)
plt.tight_layout()
fig.savefig(os.path.join(OUT, "granule_candidates.png"), dpi=100, facecolor="black")
print("wrote granule_candidates.png")
