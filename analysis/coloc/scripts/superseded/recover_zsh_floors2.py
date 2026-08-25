"""Recover per-stage z-shift floors from fig70 pixels, v2: analytic layout (I authored fig_stage.py:
figsize 12.5x6.2 dpi 200, subplots(1,2, sharey), subplots_adjust(bottom=0.2, top=0.9, wspace=0.08),
default left 0.125 right 0.9, ylim (0.6, 2.7), xlim (-0.4, 2.4), markers at data x = 0,1,2).
So marker pixel columns are computable exactly; detection = alpha-0.7 blended color within +/-6 px of
those columns; value from the analytic y-mapping. Male knowns validate with NO free parameters."""
import numpy as np
from PIL import Image

FIG = r"C:/Users/ryane/coloc_analysis/fig70_stage_resolved_pachytene.png"
OUT = r"C:/Users/ryane/coloc_analysis/zsh_floors.csv"
BLEND = {
    "noHS": np.array([44, 111, 187]) * 0.7 + 255 * 0.3,
    "HS":   np.array([192, 57, 43]) * 0.7 + 255 * 0.3,
}
KNOWN_MALE = {"noHS": [1.23, 1.18, 1.07], "HS": [1.24, 1.29, 1.17]}
TOL = 14

img = np.asarray(Image.open(FIG).convert("RGB")).astype(float)
H, W, _ = img.shape  # 1240, 2500

L, R, B, T, WS = 0.125, 0.9, 0.2, 0.9, 0.08
w = (R - L) / (2 + WS)                     # axes width (figure fraction)
ax_x = [(L, L + w), (L + w * (1 + WS), R)]  # [male, herm]
YLIM = (0.6, 2.7); XLIM = (-0.4, 2.4)


def col_px(panel, xdata):
    x0, x1 = ax_x[panel]
    fx = (xdata - XLIM[0]) / (XLIM[1] - XLIM[0])
    return (x0 + fx * (x1 - x0)) * W


def val_from_ypx(ypx):
    yfig = 1.0 - ypx / H
    fy = (yfig - B) / (T - B)
    return YLIM[0] + fy * (YLIM[1] - YLIM[0])


def detect(panel, color):
    vals = []
    for s in range(3):
        c = col_px(panel, s)
        band = img[:, int(c - 6): int(c + 7)]
        m = np.all(np.abs(band - color) < TOL, axis=2)
        ys = np.where(m.any(axis=1))[0]
        if len(ys) == 0:
            vals.append(np.nan); continue
        # z-shift floors are all < 1.6: restrict to that y-window so the solid PC line
        # (higher on the plot) cannot capture the blob
        ymin_px = int(H * (1.0 - (B + (T - B) * (1.6 - YLIM[0]) / (YLIM[1] - YLIM[0]))))
        m[:ymin_px] = False
        # marker is a contiguous ~12px blob; take the densest contiguous run's center
        rows = np.where(m.sum(axis=1) > 0)[0]
        if len(rows) == 0:  # marker overprinted by another element (e.g. a crossing solid line)
            vals.append(np.nan); continue
        splits = np.split(rows, np.where(np.diff(rows) > 3)[0] + 1)
        blob = max(splits, key=lambda r: m[r].sum())
        yc = float(np.average(blob, weights=m[blob].sum(axis=1)))
        vals.append(val_from_ypx(yc))
    return vals


print("male validation (no free parameters):")
maxres = 0.0
for t in ["noHS", "HS"]:
    got = detect(0, BLEND[t])
    for s, k, g in zip(["early", "mid", "late"], KNOWN_MALE[t], got):
        r = k - g; maxres = max(maxres, abs(r))
        print(f"  male {t:4s} {s:5s}: known {k:.3f}  recovered {g:.3f}  resid {r:+.4f}")
print(f"max |residual| = {maxres:.4f}")

rows = ["sex,treat,stage,zsh_floor,source"]
for t in ["noHS", "HS"]:
    for s, v in zip(["early", "mid", "late"], KNOWN_MALE[t]):
        rows.append(f"male,{t},{s},{v},log_exact")
    hv = detect(1, BLEND[t])
    print(f"herm {t}: {[round(x, 3) for x in hv]}")
    for s, v in zip(["early", "mid", "late"], hv):
        rows.append(f"herm,{t},{s},{v:.3f},pixel_recovered")
open(OUT, "w").write("\n".join(rows) + "\n")
print("WROTE", OUT)
