"""Shared publication figure style: Arial, journal sizing (180 mm double column), no titles or
footnote paragraphs, panel letters, Okabe-Ito colorblind-safe colors, PDF (vector) + 600 dpi PNG."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COL_NOHS = "#0072B2"     # Okabe-Ito blue
COL_HS = "#D55E00"       # Okabe-Ito vermillion

RC = {
    "font.family": "Arial",
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.linewidth": 0.6,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "legend.fontsize": 7,
    "pdf.fonttype": 42,          # editable text in Illustrator
    "ps.fonttype": 42,
    "svg.fonttype": "none",
}


def apply():
    plt.rcParams.update(RC)


def clean_axes(ax):
    ax.spines[["top", "right"]].set_visible(False)


def panel_letter(ax, letter, dx=-0.28, dy=1.05):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="bottom", ha="left")


def dots_with_mean(ax, x, values, color, jitter=0.09, seed=0):
    """open circles per gonad + black mean bar with SEM whiskers."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    rng = np.random.RandomState(abs(seed + int(abs(x) * 10)) % (2**31))
    ax.scatter(np.full(len(v), x) + rng.uniform(-jitter, jitter, len(v)), v,
               s=14, facecolors="none", edgecolors=color, linewidths=0.9, zorder=3)
    m = v.mean()
    ax.hlines(m, x - 0.18, x + 0.18, color="black", lw=1.2, zorder=4)
    if len(v) > 1:
        sem = v.std(ddof=1) / np.sqrt(len(v))
        ax.vlines(x, m - sem, m + sem, color="black", lw=0.8, zorder=4)
    return m


def p_bracket(ax, x0, x1, y, p, h=0.02):
    ax.plot([x0, x0, x1, x1], [y, y + h, y + h, y], c="black", lw=0.7, clip_on=False)
    if p is None or not np.isfinite(p):
        txt = "P = n/a"
    elif p >= 0.001:
        txt = f"P = {p:.3f}"
    else:
        txt = "P < 0.001"
    ax.text((x0 + x1) / 2, y + h * 1.4, txt, ha="center", va="bottom", fontsize=6.5)


def save(fig, stem):
    import os
    import shutil
    for ext, dpi in [("pdf", None), ("png", 600)]:
        out = f"{stem}.{ext}"
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
    desk = r"C:/Users/ryane/OneDrive/Desktop/germquant_figures"
    if os.path.isdir(desk):
        shutil.copy(f"{stem}.png", desk)
        shutil.copy(f"{stem}.pdf", desk)
    print("WROTE", stem, "(.pdf + .png)")
