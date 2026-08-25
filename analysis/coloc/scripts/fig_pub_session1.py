"""Publication figure 4: heat-shock effect using ONLY the 22-Jun imaging session (all 12 gonads).
Gonads flagged in the mask-contamination audit are drawn hollow; audited-clean gonads filled.
Means/SEM and P values use all 12 (the session has too few clean gonads to test otherwise)."""
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
d = d[d.batch == "20260622"]
POS = {("male", "noHS"): 0.0, ("male", "HS"): 0.75, ("herm", "noHS"): 2.0, ("herm", "HS"): 2.75}

fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.5))
for k, (metric, ylab, ylim) in enumerate([
        ("lamin_PC_dm", "SYP-3 partition coefficient\n(SYP$_{granule}$ - bkg) / (SYP$_{cytoplasm}$ - bkg)", (0.9, 2.3)),
        ("lamin_PC_specific", "Granule-specific SYP-3 enrichment\n(partition coefficient / z-shift control)", (0.9, 1.55))]):
    ax = axes[k]
    ax.axhline(1.0, ls=(0, (4, 3)), c="0.6", lw=0.7, zorder=1)
    ns = {}
    for (sex, treat), x in POS.items():
        s = d[(d.sex == sex) & (d.treat == treat)]
        v = s[metric].to_numpy(float)
        color = ps.COL_NOHS if treat == "noHS" else ps.COL_HS
        rng = np.random.RandomState(int(x * 10) + 3)
        xs = x + rng.uniform(-0.09, 0.09, len(v))
        clean = ~s.excl.to_numpy()
        ax.scatter(xs[clean], v[clean], s=14, color=color, zorder=3)
        ax.scatter(xs[~clean], v[~clean], s=14, facecolors="none", edgecolors=color, linewidths=0.9, zorder=3)
        m = v.mean()
        ax.hlines(m, x - 0.18, x + 0.18, color="black", lw=1.2, zorder=4)
        if len(v) > 1:
            ax.vlines(x, m - v.std(ddof=1) / np.sqrt(len(v)), m + v.std(ddof=1) / np.sqrt(len(v)), color="black", lw=0.8, zorder=4)
        ns[(sex, treat)] = len(v)
    for sex, x0, x1 in [("male", 0.0, 0.75), ("herm", 2.0, 2.75)]:
        a = d[(d.sex == sex) & (d.treat == "noHS")][metric].dropna()
        b = d[(d.sex == sex) & (d.treat == "HS")][metric].dropna()
        if len(a) >= 2 and len(b) >= 2:
            p = stats.mannwhitneyu(a, b)[1]
            y = max(a.max(), b.max()) + (ylim[1] - ylim[0]) * 0.05
            ps.p_bracket(ax, x0, x1, y, p, h=(ylim[1] - ylim[0]) * 0.02)
    ax.set_xticks(list(POS.values()))
    ax.set_xticklabels([f"No HS\n(n = {ns[('male','noHS')]})", f"HS\n(n = {ns[('male','HS')]})",
                        f"No HS\n(n = {ns[('herm','noHS')]})", f"HS\n(n = {ns[('herm','HS')]})"])
    for x, lab in [(0.375, "Male"), (2.375, "Hermaphrodite")]:
        ax.text(x, -0.24, lab, transform=ax.get_xaxis_transform(), ha="center", fontsize=7.5, fontweight="bold")
    ax.set_xlim(-0.5, 3.25)
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylab)
    ps.clean_axes(ax)
    ps.panel_letter(ax, "AB"[k])
axes[1].scatter([], [], s=14, color="0.3", label="audit-clean mask")
axes[1].scatter([], [], s=14, facecolors="none", edgecolors="0.3", label="mask flagged")
axes[1].legend(frameon=False, loc="upper left", fontsize=6.5, handletextpad=0.3)
fig.subplots_adjust(wspace=0.42, bottom=0.22)
ps.save(fig, f"{CA}/figpub4_session1_only")
for sex in ["male", "herm"]:
    a = d[(d.sex == sex) & (d.treat == "noHS")].lamin_PC_specific
    b = d[(d.sex == sex) & (d.treat == "HS")].lamin_PC_specific
    p = stats.mannwhitneyu(a, b)[1] if len(a) > 1 and len(b) > 1 else float("nan")
    print(f"{sex}: noHS {a.mean():.3f} (n{len(a)})  HS {b.mean():.3f} (n{len(b)})  delta {b.mean()-a.mean():+.3f}  p={p:.3f}")
