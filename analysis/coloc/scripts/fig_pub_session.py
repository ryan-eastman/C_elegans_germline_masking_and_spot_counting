"""Publication figure 3: heat-shock effect on granule-specific SYP-3 enrichment, blocked by imaging session.
Uses ALL 22 gonads for the within-session view (session is the dominant covariate; excluded gonads are
shown as hollow grey so their contribution is visible, and the clean-13 P values are reported)."""
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import chload
import pubstyle as ps

ps.apply()
import matplotlib.pyplot as plt

CA = r"C:/Users/ryane/coloc_analysis"
d = pd.read_csv(f"{CA}/pc_lamin_vs_dapi.csv")
d["batch"] = [chload.parse_iid(i)["batch"] for i in d.image_id]
d["excl"] = [chload.is_excluded(i) for i in d.image_id]
LAB = {"20260622": "Session 1 (22 Jun)", "20260708": "Session 2 (8 Jul)"}

fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.6), sharey=True)
for k, batch in enumerate(["20260622", "20260708"]):
    ax = axes[k]
    ax.axhline(1.0, ls=(0, (4, 3)), c="0.6", lw=0.7, zorder=1)
    pos = {("male", "noHS"): 0.0, ("male", "HS"): 0.75, ("herm", "noHS"): 2.0, ("herm", "HS"): 2.75}
    ns = {}
    for (sex, treat), x in pos.items():
        s = d[(d.batch == batch) & (d.sex == sex) & (d.treat == treat)]
        keep = s[~s.excl].lamin_PC_specific.dropna()
        drop = s[s.excl].lamin_PC_specific.dropna()
        color = ps.COL_NOHS if treat == "noHS" else ps.COL_HS
        if len(keep):
            ps.dots_with_mean(ax, x, keep, color)
        if len(drop):
            rng = np.random.RandomState(7)
            ax.scatter(np.full(len(drop), x) + rng.uniform(-0.09, 0.09, len(drop)), drop, s=14,
                       facecolors="none", edgecolors="0.6", linewidths=0.8, zorder=2)
        ns[(sex, treat)] = len(keep)
    for sex, x0, x1 in [("male", 0.0, 0.75), ("herm", 2.0, 2.75)]:
        s = d[(d.batch == batch) & (d.sex == sex) & (~d.excl)]
        a = s[s.treat == "noHS"].lamin_PC_specific.dropna()
        b = s[s.treat == "HS"].lamin_PC_specific.dropna()
        if len(a) >= 2 and len(b) >= 2:
            p = stats.mannwhitneyu(a, b)[1]
            y = max(a.max(), b.max()) + 0.04
            ps.p_bracket(ax, x0, x1, y, p, h=0.012)
    ax.set_xticks(list(pos.values()))
    ax.set_xticklabels([f"No HS\n(n = {ns[('male','noHS')]})", f"HS\n(n = {ns[('male','HS')]})",
                        f"No HS\n(n = {ns[('herm','noHS')]})", f"HS\n(n = {ns[('herm','HS')]})"])
    for x, lab in [(0.375, "Male"), (2.375, "Hermaphrodite")]:
        ax.text(x, -0.26, lab, transform=ax.get_xaxis_transform(), ha="center", fontsize=7.5, fontweight="bold")
    ax.set_title(LAB[batch], fontsize=8)
    ax.set_xlim(-0.5, 3.25)
    ps.clean_axes(ax)
    ps.panel_letter(ax, "AB"[k], dx=-0.14 if k else -0.26)
axes[0].set_ylabel("Granule-specific SYP-3 enrichment\n(partition coefficient / z-shift control)")
axes[0].set_ylim(0.95, 1.55)
axes[0].scatter([], [], facecolors="none", edgecolors="0.6", label="excluded (mask contamination)")
axes[0].legend(frameon=False, loc="upper left", fontsize=6.5)
fig.subplots_adjust(wspace=0.12, bottom=0.25)
ps.save(fig, f"{CA}/figpub3_session")
