"""fig77: lamin-envelope masking result (PI-requested), 3 panels.
A: raw distance-matched PC by group (lamin masking). B: granule-specific enrichment = PC / per-gonad
z-shift floor (noHS anchors near 1; per-gonad pairing makes this the valid corrected-scale test).
C: agreement between lamin- and DAPI-anchored metrics per gonad (r=0.997) = masking choice does not
drive the result. Whole-gonad metrics; n=22."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

CA = r"C:/Users/ryane/coloc_analysis"
OUT = f"{CA}/fig77_lamin_masked_result.png"
DESK = r"C:/Users/ryane/OneDrive/Desktop/germquant_figures"
d = pd.read_csv(f"{CA}/pc_lamin_vs_dapi.csv")
order = [("male", "noHS"), ("male", "HS"), ("herm", "noHS"), ("herm", "HS")]
labels = ["male\nnoHS", "male\nHS", "herm\nnoHS", "herm\nHS"]
col = {"noHS": "#2c6fbb", "HS": "#c0392b"}
rng = np.random.RandomState(3)

fig, ax = plt.subplots(1, 3, figsize=(15, 5.6))

for k, (metric, ttl, ylab) in enumerate([
        ("lamin_PC_dm", "A  raw partition coefficient", "SYP-3 distance-matched PC"),
        ("lamin_PC_specific", "B  granule-specific enrichment", "PC / per-gonad z-shift floor")]):
    a = ax[k]
    a.axhline(1.0, ls="--", c="#888", lw=1.2, zorder=0)
    for i, (sex, treat) in enumerate(order):
        v = d[(d.sex == sex) & (d.treat == treat)][metric].dropna()
        a.bar(i, v.mean(), width=0.62, color=col[treat], alpha=0.75, zorder=1)
        a.errorbar(i, v.mean(), yerr=v.std() / np.sqrt(len(v)), c="k", capsize=3, lw=1.1, zorder=3)
        a.scatter(np.full(len(v), i) + rng.uniform(-0.13, 0.13, len(v)), v, s=22, color="#222",
                  alpha=0.55, zorder=4)
    for sex, i0, i1, y in [("male", 0, 1, None), ("herm", 2, 3, None)]:
        va = d[(d.sex == sex) & (d.treat == "noHS")][metric].dropna()
        vb = d[(d.sex == sex) & (d.treat == "HS")][metric].dropna()
        p = stats.mannwhitneyu(va, vb)[1]
        y = max(va.max(), vb.max()) + 0.09
        a.plot([i0, i0, i1, i1], [y, y + 0.03, y + 0.03, y], c="k", lw=1.1)
        star = "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        a.text((i0 + i1) / 2, y + 0.045, f"{star} p={p:.3f}", ha="center", fontsize=9.5)
    a.set_xticks(range(4)); a.set_xticklabels(labels)
    a.set_ylabel(ylab); a.set_title(ttl, loc="left", fontweight="bold", fontsize=11)
    a.spines[["top", "right"]].set_visible(False)
ax[0].set_ylim(0.85, 2.35); ax[1].set_ylim(0.85, 1.62)

c = ax[2]
r = np.corrcoef(d.lamin_PC_specific, d.dapi_PC_specific)[0, 1]
lims = [1.0, 1.5]
c.plot(lims, lims, ls="--", c="#c0392b", lw=1.1, label="identity")
for treat, m in [("noHS", "o"), ("HS", "s")]:
    s = d[d.treat == treat]
    c.scatter(s.dapi_PC_specific, s.lamin_PC_specific, s=42, marker=m, color=col[treat],
              alpha=0.8, label=treat)
c.set_xlim(*lims); c.set_ylim(*lims); c.set_aspect("equal")
c.set_xlabel("DAPI-anchored (chromatin mask)"); c.set_ylabel("lamin-anchored (envelope mask)")
c.set_title(f"C  masking agreement  (r = {r:.3f})", loc="left", fontweight="bold", fontsize=11)
c.legend(frameon=False, fontsize=9, loc="upper left")
c.spines[["top", "right"]].set_visible(False)

fig.suptitle("Lamin-envelope masking confirms the male-specific heat-shock SYP-3 enrichment in p-granules",
             fontsize=13.5, y=0.995)
foot = ("Nuclei masked to the LMN-1 envelope (watershed on lamin, seeded by validated DAPI nuclei; median envelope/chromatin volume 1.09, "
        "3.5% weak-lamin fallbacks).\nGranule segmentation unchanged (Imaris-calibrated). B uses each gonad's own z-shift floor "
        "(paired normalization): noHS sits near 1 (no enrichment), male HS is elevated (p=0.004).\n"
        "C: per-gonad values are nearly identical under either nuclear mask, the result does not depend on the masking choice. n=22 gonads.")
fig.text(0.5, 0.012, foot, ha="center", va="bottom", fontsize=8.4, color="#333")
fig.subplots_adjust(bottom=0.24, top=0.87, wspace=0.3, left=0.06, right=0.98)
fig.savefig(OUT, dpi=200)
import os, shutil
if os.path.isdir(DESK):
    shutil.copy(OUT, os.path.join(DESK, os.path.basename(OUT)))
print("WROTE", OUT)
