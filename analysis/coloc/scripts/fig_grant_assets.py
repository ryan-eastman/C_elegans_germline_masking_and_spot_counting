"""Raster assets for the grant pipeline figure (SVG built by fig_grant_pipeline.py).
All crops from the calibration gonad HS_male_008, mid-pachytene window centred on Ryan's trace.
Outputs PNGs to coloc_analysis/grant/assets/ (2x nearest/bicubic upscaled so they print crisply at ~20 mm)."""
import json
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
os.makedirs(OUT, exist_ok=True)
SP = chload.SP
WIN_UM = 24.0
CENTRE_UM = (85.0, 122.0)        # (y, x) crop-um, mid-pachytene point of the trace
ZONE_COL = {"early": (86, 180, 233), "mid": (0, 158, 115), "late": (204, 121, 167)}   # Okabe-Ito
UP = 2


def stretch(a, lo=1, hi=99.7):
    a0, a1 = np.percentile(a, lo), np.percentile(a, hi)
    return np.clip((a - a0) / (a1 - a0 + 1e-9), 0, 1)


def up(a, order=1):
    return resize(a, (a.shape[0] * UP, a.shape[1] * UP), order=order, preserve_range=True, anti_aliasing=False)


def save(name, rgb):
    Image.fromarray((np.clip(rgb, 0, 1) * 255).astype(np.uint8)).save(os.path.join(OUT, name))
    print("wrote", name, rgb.shape[1], "x", rgb.shape[0])


def paint_edges(rgb, mask, color, thick=1):
    e = find_boundaries(mask, mode="inner")
    if thick > 1:
        e = ndi.binary_dilation(e, iterations=thick - 1)
    rgb[e] = color
    return rgb


rd = W.find_run(IID)
lab = tifffile.imread(os.path.join(rd, f"{IID}__nuclei_labels.tif"))
csv = os.path.join(rd, f"{IID}__nuclei.csv")
nc = pd.read_csv(csv)
germ = [int(x) for x in nc[nc.in_germline.astype(bool)].nucleus_id]
gn = np.isin(lab, germ)
obj = ndi.find_objects(gn.astype(np.uint8))[0]
sl = tuple(slice(max(0, o.start - W.PAD), min(dim, o.stop + W.PAD)) for o, dim in zip(obj, gn.shape))
ch = chload.load(W.find_nd2(IID, csv))
dapi = ch["dapi"][sl].astype(np.float32)
pgl = ch["pgl"][sl].astype(np.float32)
syp = ch["syp"][sl].astype(np.float32)
lamin = ch["lamin"][sl].astype(np.float32)
del ch
lab_c = lab[sl]
nuc_lam, nuc_dapi, _ = W.lamin_nuclei(lab_c, germ, lamin)

# ---- window -------------------------------------------------------------------------------------
half = int(round(WIN_UM / 2 / SP[1]))
cy, cx = int(round(CENTRE_UM[0] / SP[1])), int(round(CENTRE_UM[1] / SP[2]))
win = (slice(None), slice(cy - half, cy + half), slice(cx - half, cx + half))
env_w = nuc_lam[win]
z = int(np.argmax(env_w.reshape(env_w.shape[0], -1).sum(1)))      # plane with most envelope in window
print("window plane z =", z, f"({z * SP[0]:.1f} um)")
D, P, S, L = dapi[win][z], pgl[win][z], syp[win][z], lamin[win][z]
E = env_w[z]
Ed = nuc_dapi[win][z]

# granules + shell in this window (same recipe as the worker, on the window's cyto)
dt = ndi.distance_transform_edt(~nuc_lam, sampling=tuple(SP))
cyto = ndi.binary_fill_holes(dt <= W.CYTO_UM) & ~nuc_lam
gran, ng = W.segment_granules(pgl, cyto)
G = gran[win][z]
SH = cyto[win][z]

# 1. merged 3-channel (R = SYP-3 mCherry, G = PGL-1 GFP, B = DAPI), single plane
merge = np.stack([up(stretch(S)), up(stretch(P)), up(stretch(D, 1, 99.9))], -1)
save("merge_24um.png", merge)
# 1b. the z-stack as a stack of single-plane slices (offset cards, deepest at the back), for step 1
STACK_UM = [4.0, 6.0, 8.0, 10.0]
cards = []
for zu in STACK_UM:
    zz = int(round(zu / SP[0]))
    Dz, Pz, Sz = dapi[win][zz], pgl[win][zz], syp[win][zz]
    cards.append(np.stack([up(stretch(Sz)), up(stretch(Pz)), up(stretch(Dz, 1, 99.9))], -1))
side = cards[0].shape[0]
off = int(side * 0.16)
canvas = np.full((side + off * (len(cards) - 1), side + off * (len(cards) - 1), 3), 1.0, np.float32)   # white
n = len(cards)
for k, card in enumerate(cards[::-1]):          # deepest plane drawn first (back, top-right)
    y0 = off * k
    x0 = off * (n - 1 - k)
    canvas[y0:y0 + side, x0:x0 + side] = card
    bw = max(2, side // 150)                     # thin white border so the slices read as separate sheets
    canvas[y0:y0 + bw, x0:x0 + side] = 1.0
    canvas[y0 + side - bw:y0 + side, x0:x0 + side] = 1.0
    canvas[y0:y0 + side, x0:x0 + bw] = 1.0
    canvas[y0:y0 + side, x0 + side - bw:x0 + side] = 1.0
save("stack_slices.png", canvas)
# 2. lamin gray + envelope edges (red)
lam_rgb = np.stack([up(stretch(L, 5, 99.8))] * 3, -1)
save("lamin_raw_24um.png", lam_rgb.copy())
save("lamin_masks_24um.png", paint_edges(lam_rgb.copy(), up(E, 0) > 0.5, (1.0, 0.2, 0.2), 2))
# 2b. DAPI label (blue) vs lamin envelope (red) on lamin gray
both = paint_edges(lam_rgb.copy(), up(Ed, 0) > 0.5, (0.25, 0.45, 1.0), 2)
save("masks_dapi_vs_lamin_24um.png", paint_edges(both, up(E, 0) > 0.5, (1.0, 0.2, 0.2), 2))
# 3. P-granule panel: its OWN window and plane (chosen from grant_granule_candidates.py: early-pachytene
#    arm, many crisp perinuclear granules). PGL green + granule outlines (yellow) + envelope (red) +
#    shell boundary (cyan).
GRAN_CENTRE_UM, GRAN_Z_UM = (150.0, 146.0), 10.0
gcy, gcx = int(round(GRAN_CENTRE_UM[0] / SP[1])), int(round(GRAN_CENTRE_UM[1] / SP[2]))
gz = int(round(GRAN_Z_UM / SP[0]))
gwin = (gz, slice(gcy - half, gcy + half), slice(gcx - half, gcx + half))
Pg, Dg, Sg = pgl[gwin], dapi[gwin], syp[gwin]
Eg, Gg, SHg = nuc_lam[gwin], gran[gwin], cyto[gwin]
pg = np.stack([up(stretch(Pg, 1, 99.8)) * 0.15, up(stretch(Pg, 1, 99.8)), up(stretch(Pg, 1, 99.8)) * 0.15], -1)
pg = paint_edges(pg, up(SHg, 0) > 0.5, (0.0, 0.85, 0.95), 1)
pg = paint_edges(pg, up(Eg, 0) > 0.5, (1.0, 0.2, 0.2), 2)
pg = paint_edges(pg, up(Gg, 0) > 0.5, (1.0, 0.9, 0.1), 2)
save("granules_24um.png", pg)
# 3b. granule outlines on the merged image of the same window
mgi = np.stack([up(stretch(Sg)), up(stretch(Pg)), up(stretch(Dg, 1, 99.9))], -1)
mg = paint_edges(mgi, up(Eg, 0) > 0.5, (1.0, 1.0, 1.0), 1)
save("granules_on_merge_24um.png", paint_edges(mg, up(Gg, 0) > 0.5, (1.0, 0.9, 0.1), 2))
print("granule panel window", GRAN_CENTRE_UM, "plane", GRAN_Z_UM, "um; granules in window:", int(ndi.label(Gg)[1]))
# 4. single channels (spares)
save("dapi_24um.png", np.stack([up(stretch(D, 1, 99.9))] * 3, -1))
save("pgl_24um.png", np.stack([np.zeros_like(up(P)), up(stretch(P)), np.zeros_like(up(P))], -1))
save("syp_24um.png", np.stack([up(stretch(S)), np.zeros_like(up(S)), np.zeros_like(up(S))], -1))
# 5. wider 40 um window merge + masks (spare)
half2 = int(round(20.0 / SP[1]))
w2 = (slice(cy - half2, cy + half2), slice(cx - half2, cx + half2))
m2 = np.stack([stretch(syp[z][w2]), stretch(pgl[z][w2]), stretch(dapi[z][w2], 1, 99.9)], -1)
save("merge_40um.png", m2)
l2 = np.stack([stretch(lamin[z][w2], 5, 99.8)] * 3, -1)
save("lamin_masks_40um.png", paint_edges(l2, nuc_lam[z][w2], (1.0, 0.2, 0.2), 1))

# 6. whole gonad: DAPI max projection (downsampled) + zone-coloured nuclei + trace polyline
ds = 4
mp = stretch(dapi.max(0), 5, 99.8)
mp = resize(mp, (mp.shape[0] // ds, mp.shape[1] // ds), order=1, anti_aliasing=True)
whole = np.stack([mp] * 3, -1) * 0.75
zones = pd.read_csv(rf"C:/Users/ryane/coloc_analysis/staging/zones/{IID}_zones.csv")
zoned = zones[zones.zone.isin(ZONE_COL)]
rad = max(2, int(round(2.0 / SP[1] / ds)))
yy, xx = np.ogrid[-rad:rad + 1, -rad:rad + 1]
disk = (yy ** 2 + xx ** 2) <= rad * rad
for _, r in zoned.iterrows():
    yi, xi = int(round(r.cy / SP[1] / ds)), int(round(r.cx / SP[2] / ds))
    y0, y1 = max(0, yi - rad), min(whole.shape[0], yi + rad + 1)
    x0, x1 = max(0, xi - rad), min(whole.shape[1], xi + rad + 1)
    d = disk[(y0 - (yi - rad)):(y1 - (yi - rad)), (x0 - (xi - rad)):(x1 - (xi - rad))]
    col = np.array(ZONE_COL[r.zone]) / 255
    sub = whole[y0:y1, x0:x1]
    sub[d] = sub[d] * 0.3 + col * 0.7
tr = json.load(open(r"C:/Users/ryane/coloc_analysis/staging/pachytene_traces.json", encoding="utf-8"))[IID]["points_um"]
pts = np.array(tr)
from skimage.draw import line as sk_line
for a, b in zip(pts[:-1], pts[1:]):
    r0, c0 = int(round(a[1] / SP[1] / ds)), int(round(a[0] / SP[2] / ds))
    r1, c1 = int(round(b[1] / SP[1] / ds)), int(round(b[0] / SP[2] / ds))
    rr, cc = sk_line(r0, c0, r1, c1)
    ok = (rr >= 0) & (rr < whole.shape[0]) & (cc >= 0) & (cc < whole.shape[1])
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            r2, c2 = np.clip(rr[ok] + dr, 0, whole.shape[0] - 1), np.clip(cc[ok] + dc, 0, whole.shape[1] - 1)
            whole[r2, c2] = (1.0, 1.0, 1.0)
save("gonad_zones.png", whole)
# whole gonad merged max projection (spare)
wm = np.stack([stretch(syp.max(0)), stretch(pgl.max(0)), stretch(dapi.max(0), 5, 99.8)], -1)
wm = resize(wm, (wm.shape[0] // ds, wm.shape[1] // ds, 3), order=1, anti_aliasing=True)
save("gonad_merge.png", wm)
# whole gonad LMN-1 + envelope edges (spare). ONE plane, not a projection: projecting the 3D masks
# merges nuclei stacked in z into blobs. Plane = the one with the most envelope voxels; 2x downsample.
ds2 = 2
zb = int(np.argmax(nuc_lam.reshape(nuc_lam.shape[0], -1).sum(1)))
wl = np.stack([stretch(lamin[zb], 5, 99.8)] * 3, -1)
wl = resize(wl, (wl.shape[0] // ds2, wl.shape[1] // ds2, 3), order=1, anti_aliasing=True)
em = resize(nuc_lam[zb].astype(float), (wl.shape[0], wl.shape[1]), order=0) > 0.5
save("gonad_lamin_masks.png", paint_edges(wl, em, (1.0, 0.2, 0.2), 1))
print("whole-gonad lamin plane:", round(zb * SP[0], 1), "um")
json.dump({"window_um": WIN_UM, "plane_um": round(z * SP[0], 1), "centre_um": CENTRE_UM, "n_granules_gonad": int(ng),
           "gonad_px_per_um": 1 / (SP[1] * ds), "whole_shape": [int(v) for v in whole.shape[:2]]},
          open(os.path.join(OUT, "meta.json"), "w"), indent=1)
print("DONE")
