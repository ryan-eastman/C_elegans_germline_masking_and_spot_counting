"""Supplementary figure 1S: WHOLE-GONAD SYP-3 partitioning into P granules, clean-13 dataset (all
in_germline labels; carries sperm / somatic contamination in several males, see README "Mask audit";
the primary figure 1 is the pooled hand-traced pachytene version, fig_pub_pachytene_pooled.py).
Two panels, journal style (Arial, no titles/footers, explicit axis formulas, exact P values).
Data: pc_clean13.csv (lamin-envelope masking)."""
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import pubstyle as ps

ps.apply()
import matplotlib.pyplot as plt

CA = r"C:/Users/ryane/coloc_analysis"
d = pd.read_csv(f"{CA}/pc_clean13.csv")

fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.5))
POS = {("male", "noHS"): 0.0, ("male", "HS"): 0.75, ("herm", "noHS"): 2.0, ("herm", "HS"): 2.75}

for k, (metric, ylab, ylim) in enumerate([
        ("lamin_PC_dm",
         "SYP-3 partition coefficient\n(SYP$_{granule}$ - bkg) / (SYP$_{cytoplasm}$ - bkg)",
         (0.9, 2.35)),
        ("lamin_PC_specific",
         "Granule-specific SYP-3 enrichment\n(partition coefficient / z-shift control)",
         (0.9, 1.62))]):
    ax = axes[k]
    ax.axhline(1.0, ls=(0, (4, 3)), c="0.6", lw=0.7, zorder=1)
    ns = {}
    for (sex, treat), x in POS.items():
        v = d[(d.sex == sex) & (d.treat == treat)][metric].dropna()
        ns[(sex, treat)] = len(v)
        ps.dots_with_mean(ax, x, v, ps.COL_NOHS if treat == "noHS" else ps.COL_HS)
    for sex, x0, x1 in [("male", 0.0, 0.75), ("herm", 2.0, 2.75)]:
        a = d[(d.sex == sex) & (d.treat == "noHS")][metric].dropna()
        b = d[(d.sex == sex) & (d.treat == "HS")][metric].dropna()
        p = stats.mannwhitneyu(a, b)[1]
        ytop = max(a.max(), b.max()) + (ylim[1] - ylim[0]) * 0.05
        ps.p_bracket(ax, x0, x1, ytop, p, h=(ylim[1] - ylim[0]) * 0.02)
    ax.set_xticks(list(POS.values()))
    ax.set_xticklabels([f"No HS\n(n = {ns[('male','noHS')]})", f"HS\n(n = {ns[('male','HS')]})",
                        f"No HS\n(n = {ns[('herm','noHS')]})", f"HS\n(n = {ns[('herm','HS')]})"])
    for x, lab in [(0.375, "Male"), (2.375, "Hermaphrodite")]:
        ax.text(x, -0.24, lab, transform=ax.get_xaxis_transform(), ha="center", fontsize=7.5,
                fontweight="bold")
    ax.set_xlim(-0.5, 3.25)
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylab)
    ps.clean_axes(ax)
    ps.panel_letter(ax, "AB"[k])

fig.subplots_adjust(wspace=0.42, bottom=0.22)
ps.save(fig, f"{CA}/figpub1S_whole_gonad")
