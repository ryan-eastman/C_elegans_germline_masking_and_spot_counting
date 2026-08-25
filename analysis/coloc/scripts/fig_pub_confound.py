"""Publication figure 5: the imaging session as a confound.
A: granule-specific enrichment by session, split by treatment (marker = sex). The HS group shifts up in
   session 2 while No HS does not: the treatment effect differs by session.
B: SYP-3 image contrast (nuclear / cytoplasmic median, background-floored) by session: session 2 was
   imaged dimmer.
C: enrichment against that contrast, per gonad: no relationship (Spearman shown), so the session offset
   is not an image-brightness artifact in the metric itself."""
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import pubstyle as ps

ps.apply()
import matplotlib.pyplot as plt

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
import chload

CA = r"C:/Users/ryane/coloc_analysis"
d = pd.read_csv(f"{CA}/session_covariates.csv")
d["excl"] = [chload.is_excluded(i) for i in d.image_id]      # live exclusion, not the frozen column
SES = {"20260622": 0, "20260708": 1}
SLAB = ["Session 1\n(22 Jun)", "Session 2\n(8 Jul)"]
MK = {"male": "o", "herm": "s"}
d["sx"] = d.batch.astype(str).map(SES)
clean = d[~d.excl]                                            # all statistics use the audit-clean 13

fig, ax = plt.subplots(1, 3, figsize=(7.09, 2.5))

# A: PC_specific by session and treatment
a = ax[0]
a.axhline(1.0, ls=(0, (4, 3)), c="0.6", lw=0.7, zorder=1)
for treat, color, off in [("noHS", ps.COL_NOHS, -0.16), ("HS", ps.COL_HS, 0.16)]:
    means = []
    for s in (0, 1):
        sub = d[(d.sx == s) & (d.treat == treat)]
        rng = np.random.RandomState(s * 3 + int(off > 0))
        for sex in ["male", "herm"]:
            ss = sub[sub.sex == sex]
            for excl_flag, ec in [(False, color), (True, "0.65")]:      # excluded gonads: grey, not in stats
                v = ss[ss.excl == excl_flag].lamin_PC_specific.to_numpy()
                a.scatter(np.full(len(v), s + off) + rng.uniform(-0.06, 0.06, len(v)), v, s=14, marker=MK[sex],
                          facecolors="none", edgecolors=ec, linewidths=0.9, zorder=3)
        kept = sub[~sub.excl].lamin_PC_specific
        m = kept.mean()
        means.append(m)
        a.hlines(m, s + off - 0.12, s + off + 0.12, color="black", lw=1.2, zorder=4)
    a.plot([off, 1 + off], means, color=color, lw=1.0, zorder=2, label="No HS" if treat == "noHS" else "HS")
a.set_xticks([0, 1]); a.set_xticklabels(SLAB)
a.set_xlim(-0.55, 1.55); a.set_ylim(0.95, 1.5)
a.set_ylabel("Granule-specific SYP-3 enrichment\n(partition coefficient / z-shift control)")
a.legend(frameon=False, loc="upper left", fontsize=6.5, handlelength=1.2)
a.scatter([], [], marker="o", facecolors="none", edgecolors="0.3", label="male")
a.scatter([], [], marker="s", facecolors="none", edgecolors="0.3", label="hermaphrodite")
a.legend(frameon=False, loc="upper left", fontsize=6.5, handlelength=1.2, ncol=1)
ps.clean_axes(a); ps.panel_letter(a, "A")

# B: image contrast by session
b = ax[1]
for s in (0, 1):
    sub = d[d.sx == s]
    for treat, color in [("noHS", ps.COL_NOHS), ("HS", ps.COL_HS)]:
        rng = np.random.RandomState(s + 11)
        for excl_flag, ec in [(False, color), (True, "0.65")]:
            v = sub[(sub.treat == treat) & (sub.excl == excl_flag)].syp_contrast.to_numpy()
            b.scatter(np.full(len(v), s) + rng.uniform(-0.12, 0.12, len(v)), v, s=14, facecolors="none",
                      edgecolors=ec, linewidths=0.9, zorder=3)
    kept = sub[~sub.excl].syp_contrast
    m = kept.mean()
    b.hlines(m, s - 0.2, s + 0.2, color="black", lw=1.2, zorder=4)
    if len(kept) > 1:
        b.vlines(s, m - kept.std(ddof=1) / np.sqrt(len(kept)), m + kept.std(ddof=1) / np.sqrt(len(kept)), color="black", lw=0.8)
p_b = stats.mannwhitneyu(clean[clean.sx == 0].syp_contrast, clean[clean.sx == 1].syp_contrast)[1]
ps.p_bracket(b, 0, 1, d.syp_contrast.max() + 0.25, p_b, h=0.12)
b.set_xticks([0, 1]); b.set_xticklabels(SLAB)
b.set_xlim(-0.55, 1.55); b.set_ylim(1.0, 7.2)
b.set_ylabel("SYP-3 image contrast\n(nuclear median / cytoplasmic median)")
ps.clean_axes(b); ps.panel_letter(b, "B")

# C: enrichment vs contrast
c = ax[2]
for s, ec in [(0, "0.55"), (1, "black")]:
    sub = d[d.sx == s]
    for treat, color in [("noHS", ps.COL_NOHS), ("HS", ps.COL_HS)]:
        for excl_flag, col_e in [(False, color), (True, "0.65")]:
            v = sub[(sub.treat == treat) & (sub.excl == excl_flag)]
            c.scatter(v.syp_contrast, v.lamin_PC_specific, s=16, marker="o" if s == 0 else "D",
                      facecolors=col_e if (s == 1 and not excl_flag) else "none", edgecolors=col_e, linewidths=0.9, zorder=3)
rho, p_c = stats.spearmanr(clean.syp_contrast, clean.lamin_PC_specific)
c.text(0.97, 0.05, f"Spearman $\\rho$ = {rho:.2f}, P = {p_c:.2f}", transform=c.transAxes, ha="right", fontsize=6.5)
c.scatter([], [], marker="o", facecolors="none", edgecolors="0.3", label="Session 1")
c.scatter([], [], marker="D", facecolors="0.3", edgecolors="0.3", label="Session 2")
c.scatter([], [], marker="o", facecolors="none", edgecolors="0.65", label="excluded (not in statistics)")
c.legend(frameon=False, loc="upper left", fontsize=6.5, handletextpad=0.3)
c.set_xlabel("SYP-3 image contrast")
c.set_ylabel("Granule-specific SYP-3 enrichment")
c.set_ylim(0.95, 1.5)
ps.clean_axes(c); ps.panel_letter(c, "C")

fig.subplots_adjust(wspace=0.55, bottom=0.24, left=0.09, right=0.99)
ps.save(fig, f"{CA}/figpub5_session_confound")
print(f"contrast by session P={p_b:.3f}; PC vs contrast rho={rho:.2f} P={p_c:.2f}")
