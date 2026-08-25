"""Recover the per-stage z-shift (artifact floor) group means from fig70's pixels.
The fig70 stage figure plotted z-shift means as alpha-0.7 'x' markers on dotted lines (noHS #2c6fbb,
HS #c0392b). The male floors are known EXACTLY from the analysis log (noHS 1.23/1.18/1.07, HS
1.24/1.29/1.17), so the male panel serves as calibration: fit value = a*y_pixel + b on those 6 known
points, check residuals, then apply the same mapping (sharey axes) to the herm panel markers.
Outputs zsh_floors.csv with all 12 floors + calibration residuals."""
import numpy as np
from PIL import Image

FIG = r"C:/Users/ryane/coloc_analysis/fig70_stage_resolved_pachytene.png"
OUT = r"C:/Users/ryane/coloc_analysis/zsh_floors.csv"
# alpha-0.7 blend over white: c*0.7 + 255*0.3
BLEND = {
    "noHS": np.array([44, 111, 187]) * 0.7 + 255 * 0.3,   # (107.3, 154.2, 207.4)
    "HS":   np.array([192, 57, 43]) * 0.7 + 255 * 0.3,    # (210.9, 116.4, 106.6)
}
KNOWN_MALE = {"noHS": [1.23, 1.18, 1.07], "HS": [1.24, 1.29, 1.17]}
TOL = 14  # per-channel color tolerance

img = np.asarray(Image.open(FIG).convert("RGB")).astype(float)
H, W, _ = img.shape
print(f"figure {W}x{H}")

# panel split: male panel = left half, herm = right half (two subplots)
panels = {"male": (0, W // 2), "herm": (W // 2, W)}


def marker_ys(panel_x0, panel_x1, color):
    """find 3 stage x-clusters of matching pixels, return their median y (pixel)."""
    sub = img[:, panel_x0:panel_x1]
    m = np.all(np.abs(sub - color) < TOL, axis=2)
    ys, xs = np.where(m)
    if len(xs) < 30:
        return None
    # cluster columns: find 3 densest x-groups (markers add pixel mass at the 3 stage positions)
    hist, edges = np.histogram(xs, bins=120)
    # take the 3 highest well-separated peaks
    order = np.argsort(hist)[::-1]
    peaks = []
    for b in order:
        c = 0.5 * (edges[b] + edges[b + 1])
        if all(abs(c - p) > (panel_x1 - panel_x0) * 0.15 for p in peaks):
            peaks.append(c)
        if len(peaks) == 3:
            break
    peaks = sorted(peaks)
    out = []
    for p in peaks:
        sel = np.abs(xs - p) < 8
        out.append(float(np.median(ys[sel])))
    return out  # y pixels for stages [early, mid, late]


# 1) male panel: get pixel ys for both colors, fit calibration on the 6 known values
ys_m = {t: marker_ys(*panels["male"], BLEND[t]) for t in ["noHS", "HS"]}
ypix = np.array(ys_m["noHS"] + ys_m["HS"])
vals = np.array(KNOWN_MALE["noHS"] + KNOWN_MALE["HS"])
A = np.vstack([ypix, np.ones_like(ypix)]).T
(a, b), res, *_ = np.linalg.lstsq(A, vals, rcond=None)
fit = a * ypix + b
resid = vals - fit
print("male calibration: recovered vs known:")
for i, (v, f) in enumerate(zip(vals, fit)):
    print(f"  known {v:.3f}  from-pixels {f:.3f}  resid {v - f:+.4f}")
print(f"max |residual| = {np.max(np.abs(resid)):.4f} PC units")

# 2) herm panel: apply mapping
ys_h = {t: marker_ys(*panels["herm"], BLEND[t]) for t in ["noHS", "HS"]}
rows = ["sex,treat,stage,zsh_floor,source"]
for t in ["noHS", "HS"]:
    for s, v in zip(["early", "mid", "late"], KNOWN_MALE[t]):
        rows.append(f"male,{t},{s},{v},log_exact")
    hv = [a * y + b for y in ys_h[t]]
    for s, v in zip(["early", "mid", "late"], hv):
        rows.append(f"herm,{t},{s},{v:.3f},pixel_recovered")
    print(f"herm {t}: early/mid/late = {[round(x,3) for x in hv]}")
open(OUT, "w").write("\n".join(rows) + "\n")
print("WROTE", OUT)
