"""fig78: whole-gonad SYP-3 partition coefficient, CONTAMINATION-FILTERED (supersedes fig77 reporting).

Change vs fig77: a per-gonad visual audit (22 independent reviews, 2026-08-22) found fused somatic /
carcass nuclei, spermatid fields and second germline limbs inside the germline mask of 9 of the 12 0622
gonads. That material sits in the cytoplasm compartment that the partition coefficient divides by, so it
biases the denominator. Those 9 are excluded here. The male effect roughly DOUBLES, i.e. contamination
was diluting it toward the null rather than creating it.

Also corrected: fig77 reported hermaphrodites as a null result. After exclusions herms are n=2 vs 2,
where an exact Mann-Whitney cannot return anything below p=0.333. Herms are UNDERPOWERED, not negative,
and the panel now says so explicitly.

Panels: A raw PC by group; B granule-specific enrichment (PC / per-gonad z-shift floor, 1.0 = no
enrichment); C the exclusion sensitivity ladder (all-22 -> clean-13 -> within-batch 0708) showing the
effect grows as contamination is removed. Batch offset is annotated, not hidden."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

CA = r"C:/Users/ryane/coloc_analysis"
OUT = f"{CA}/fig78_clean13_partition.png"
DESK = r"C:/Users/ryane/OneDrive/Desktop/germquant_figures"

DROP = {"HS_herm_14", "HS_herm_16", "HS_male_07", "HS_male_11", "HS_male_13",
        "HS_male_15", "noHS_herm_01", "noHS_herm_02", "noHS_male_05"}

d = pd.read_csv(f"{CA}/pc_lamin_vs_dapi.csv")
d["short"] = [i.split("_LMN1_")[-1] if "_LMN1_" in i else i.split("_lmn1_")[-1] for i in d.image_id]
d["batch"] = np.where(d.image_id.str.contains("20260622"), "0622", "0708")
d["excl"] = (d.batch == "0622") & d.short.isin(DROP)
clean = d[d.excl == False].copy()
clean.to_csv(f"{CA}/pc_clean13.csv", index=False)

order = [("male", "noHS"), ("male", "HS"), ("herm", "noHS"), ("herm", "HS")]
labels = ["male\nnoHS", "male\nHS", "herm\nnoHS", "herm\nHS"]
col = {"noHS": "#2c6fbb", "HS": "#c0392b"}
rng = np.random.RandomState(4)

fig, ax = plt.subplots(1, 3, figsize=(15.5, 5.8))

for k, (metric, ttl, ylab, ylim) in enumerate([
        ("lamin_PC_dm", "A  raw partition coefficient", "SYP-3 distance-matched PC", (0.85, 2.35)),
        ("lamin_PC_specific", "B  granule-specific enrichment", "PC / per-gonad z-shift floor", (0.9, 1.62))]):
    a = ax[k]
    a.axhline(1.0, ls="--", c="#888", lw=1.2, zorder=0)
    for i, (sex, treat) in enumerate(order):
        v = clean[(clean.sex == sex) & (clean.treat == treat)][metric].dropna()
        a.bar(i, v.mean(), width=0.62, color=col[treat], alpha=0.75, zorder=1)
        if len(v) > 1:
            a.errorbar(i, v.mean(), yerr=v.std() / np.sqrt(len(v)), c="k", capsize=3, lw=1.1, zorder=3)
        a.scatter(np.full(len(v), i) + rng.uniform(-0.13, 0.13, len(v)), v, s=26, color="#222",
                  alpha=0.6, zorder=4)
        a.text(i, ylim[0] + 0.02, f"n={len(v)}", ha="center", fontsize=8.5, color="#444")
    for sex, i0, i1 in [("male", 0, 1), ("herm", 2, 3)]:
        va = clean[(clean.sex == sex) & (clean.treat == "noHS")][metric].dropna()
        vb = clean[(clean.sex == sex) & (clean.treat == "HS")][metric].dropna()
        p = stats.mannwhitneyu(va, vb)[1]
        y = max(va.max(), vb.max()) + 0.07
        a.plot([i0, i0, i1, i1], [y, y + 0.025, y + 0.025, y], c="k", lw=1.1)
        if sex == "male":
            txt = f"* p={p:.3f}"
        else:
            txt = f"p={p:.2f} (floor {1/3:.2f})\nUNDERPOWERED n=2v2"
        a.text((i0 + i1) / 2, y + 0.04, txt, ha="center", fontsize=8.6,
               color="#111" if sex == "male" else "#a33")
    a.set_xticks(range(4)); a.set_xticklabels(labels)
    a.set_ylabel(ylab); a.set_ylim(*ylim)
    a.set_title(ttl, loc="left", fontweight="bold", fontsize=11)
    a.spines[["top", "right"]].set_visible(False)

# panel C: sensitivity ladder
c = ax[2]
sets = [("all 22\n(fig77)", d), ("clean 13\n(9 excluded)", clean), ("0708 only\n(4v4)", d[d.batch == "0708"])]
xs = np.arange(len(sets))
deltas, ps, ns = [], [], []
for _, sub in sets:
    s = sub[sub.sex == "male"]
    a_ = s[s.treat == "noHS"].lamin_PC_specific.dropna()
    b_ = s[s.treat == "HS"].lamin_PC_specific.dropna()
    deltas.append(b_.median() - a_.median())
    ps.append(stats.mannwhitneyu(a_, b_)[1])
    ns.append((len(b_), len(a_)))
c.bar(xs, deltas, width=0.55, color=["#b8c9d9", "#5590c0", "#0b6ba8"])
for x, dv, p, n in zip(xs, deltas, ps, ns):
    c.text(x, dv + 0.006, f"+{dv:.3f}\np={p:.3f}\n{n[0]}v{n[1]}", ha="center", fontsize=9)
c.set_xticks(xs); c.set_xticklabels([s[0] for s in sets])
c.set_ylabel("male HS - noHS  (median PC_specific)")
c.set_ylim(0, max(deltas) * 1.45)
c.set_title("C  effect grows as contamination is removed", loc="left", fontweight="bold", fontsize=11)
c.spines[["top", "right"]].set_visible(False)

fig.suptitle("Heat shock enriches SYP-3 in p-granules in males: effect roughly doubles once contaminated gonads are excluded",
             fontsize=13, y=0.995)
foot = ("Excluded: 9 of 12 gonads from the 2026-06-22 batch whose germline mask contained fused somatic/carcass nuclei, spermatid fields\n"
        "or a second germline limb (per-gonad visual audit); that material lands in the cytoplasm compartment the PC divides by.\n"
        "Male HS-vs-noHS: +0.118 (all 22, p=0.004) -> +0.184 (clean 13, p=0.032) -> +0.226 (0708 only, p=0.029); HS>noHS in 4/4 sex x batch cells.\n"
        "Batch offset: male-HS median PC_specific 1.227 (0622) vs 1.384 (0708), larger than the whole herm difference, so the 0708 test is the cleanest read.\n"
        "Hermaphrodites are NOT null here: n=2v2 cannot return p below 0.333, and the herm effect (+0.176) is comparable to the male one.")
fig.text(0.5, 0.012, foot, ha="center", va="bottom", fontsize=8.4, color="#333")
fig.subplots_adjust(bottom=0.27, top=0.87, wspace=0.28, left=0.06, right=0.98)
fig.savefig(OUT, dpi=200)
import os, shutil
if os.path.isdir(DESK):
    shutil.copy(OUT, os.path.join(DESK, os.path.basename(OUT)))
print("clean n =", len(clean), "| male", ns[1], "delta", round(deltas[1], 3), "p", round(ps[1], 4))
print("WROTE", OUT)
