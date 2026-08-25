"""Publication figure 2: SYP-3 partitioning by position within pachytene (hand-traced thirds).
Journal style; y axis states the exact quantity; per-gonad points + mean +/- SEM; exact P values."""
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import pubstyle as ps

ps.apply()
import matplotlib.pyplot as plt

CA = r"C:/Users/ryane/coloc_analysis"
sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
import chload
d = pd.read_csv(f"{CA}/pc_zone_all.csv")
d = d[[("ccw77" in i) and not chload.is_excluded(i) for i in d.image_id]]   # clean ccw77 gonads only
ZONES = ["early", "mid", "late"]
XT = [0, 1, 2]

fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.5), sharey=True)
for k, (sex, ttl) in enumerate([("male", "Male"), ("herm", "Hermaphrodite")]):
    ax = axes[k]
    ax.axhline(1.0, ls=(0, (4, 3)), c="0.6", lw=0.7, zorder=1)
    for treat, color, off in [("noHS", ps.COL_NOHS, -0.13), ("HS", ps.COL_HS, 0.13)]:
        sub = d[(d.sex == sex) & (d.treat == treat)]
        means = []
        for j, z in enumerate(ZONES):
            v = sub[f"{z}_PCspec"].dropna()
            m = ps.dots_with_mean(ax, XT[j] + off, v, color, jitter=0.05, seed=j)
            means.append(m)
        ax.plot(np.array(XT) + off, means, color=color, lw=1.0, zorder=2,
                label="No HS" if treat == "noHS" else "HS")
    for j, z in enumerate(ZONES):
        a = d[(d.sex == sex) & (d.treat == "noHS")][f"{z}_PCspec"].dropna()
        b = d[(d.sex == sex) & (d.treat == "HS")][f"{z}_PCspec"].dropna()
        if len(a) >= 2 and len(b) >= 2:
            p = stats.mannwhitneyu(a, b)[1]
            ytop = max(a.max(), b.max()) + 0.035
            ps.p_bracket(ax, XT[j] - 0.13, XT[j] + 0.13, ytop, p, h=0.012)
    n_no = len(d[(d.sex == sex) & (d.treat == "noHS")])
    n_hs = len(d[(d.sex == sex) & (d.treat == "HS")])
    ax.set_title(f"{ttl} (No HS n = {n_no}, HS n = {n_hs})", fontsize=7.5)
    ax.set_xticks(XT)
    ax.set_xticklabels(["Early", "Mid", "Late"])
    ax.set_xlabel("Position within pachytene\n(equal thirds of traced region)")
    ax.set_xlim(-0.55, 2.55)
    ps.clean_axes(ax)
    ps.panel_letter(ax, "AB"[k], dx=-0.14 if k else -0.26)
axes[0].set_ylabel("Granule-specific SYP-3 enrichment\n(partition coefficient / z-shift control)")
axes[0].set_ylim(0.95, 1.62)
axes[0].legend(frameon=False, loc="upper center", ncol=2, handlelength=1.4, bbox_to_anchor=(0.62, 1.0), columnspacing=1.0)
fig.subplots_adjust(wspace=0.12, bottom=0.26)
ps.save(fig, f"{CA}/figpub2_stage")
for sex in ["male", "herm"]:
    for z in ZONES:
        a = d[(d.sex == sex) & (d.treat == "noHS")][f"{z}_PCspec"].dropna()
        b = d[(d.sex == sex) & (d.treat == "HS")][f"{z}_PCspec"].dropna()
        if len(a) >= 2 and len(b) >= 2:
            print(f"{sex} {z}: noHS {a.mean():.3f} HS {b.mean():.3f} p={stats.mannwhitneyu(a,b)[1]:.4f}")
