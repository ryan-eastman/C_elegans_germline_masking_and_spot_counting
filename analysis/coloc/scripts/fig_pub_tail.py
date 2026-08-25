"""Publication figure 6: fraction of P granules that visibly "light up" with SYP-3.
Per granule, SYP-3 excess over the distance-matched local cytoplasm, in units of the nuclear SYP-3
level; a granule with excess > 0.5 holds at least half the nuclear level above its surroundings.
Panel A: fraction of granules with excess > 0.5, per gonad, by sex and treatment.
Panel B: 90th-percentile granule excess per gonad (the brightest granules)."""
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import pubstyle as ps

ps.apply()
import matplotlib.pyplot as plt

CA = r"C:/Users/ryane/coloc_analysis"
d = pd.read_csv(f"{CA}/granule_tail_v3.csv")   # v3 = current per-granule pass (adds distance + PGL columns)
POS = {("male", "noHS"): 0.0, ("male", "HS"): 0.75, ("herm", "noHS"): 2.0, ("herm", "HS"): 2.75}

fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.5))
for k, (metric, ylab) in enumerate([
        ("frac_excess_gt_0.5", "Fraction of P granules with SYP-3\n> 0.5 x nuclear level above local cytoplasm"),
        ("excess_p90", "90th-percentile granule SYP-3 excess\n(units of nuclear SYP-3 level)")]):
    ax = axes[k]
    ns = {}
    for (sex, treat), x in POS.items():
        v = d[(d.sex == sex) & (d.treat == treat)][metric].dropna()
        ns[(sex, treat)] = len(v)
        ps.dots_with_mean(ax, x, v, ps.COL_NOHS if treat == "noHS" else ps.COL_HS)
    ymax = d[metric].max()
    for sex, x0, x1 in [("male", 0.0, 0.75), ("herm", 2.0, 2.75)]:
        a = d[(d.sex == sex) & (d.treat == "noHS")][metric].dropna()
        b = d[(d.sex == sex) & (d.treat == "HS")][metric].dropna()
        if len(a) >= 2 and len(b) >= 2:
            p = stats.mannwhitneyu(a, b)[1]
            ps.p_bracket(ax, x0, x1, max(a.max(), b.max()) + ymax * 0.06, p, h=ymax * 0.025)
    # HS male vs HS herm
    a = d[(d.sex == "male") & (d.treat == "HS")][metric].dropna()
    b = d[(d.sex == "herm") & (d.treat == "HS")][metric].dropna()
    if len(a) >= 2 and len(b) >= 2:
        ps.p_bracket(ax, 0.75, 2.75, ymax * 1.22, stats.mannwhitneyu(a, b)[1], h=ymax * 0.025)
    ax.set_xticks(list(POS.values()))
    ax.set_xticklabels([f"No HS\n(n = {ns[('male','noHS')]})", f"HS\n(n = {ns[('male','HS')]})",
                        f"No HS\n(n = {ns[('herm','noHS')]})", f"HS\n(n = {ns[('herm','HS')]})"])
    for x, lab in [(0.375, "Male"), (2.375, "Hermaphrodite")]:
        ax.text(x, -0.24, lab, transform=ax.get_xaxis_transform(), ha="center", fontsize=7.5, fontweight="bold")
    ax.set_xlim(-0.5, 3.25)
    ax.set_ylim(0, ymax * 1.4)
    ax.set_ylabel(ylab)
    ps.clean_axes(ax)
    ps.panel_letter(ax, "AB"[k])
fig.subplots_adjust(wspace=0.42, bottom=0.22)
ps.save(fig, f"{CA}/figpub6_granule_tail")
print(d.groupby(["sex", "treat"])[["frac_excess_gt_0.25", "frac_excess_gt_0.5", "frac_excess_gt_1.0", "excess_p90"]].mean().round(3).to_string())
