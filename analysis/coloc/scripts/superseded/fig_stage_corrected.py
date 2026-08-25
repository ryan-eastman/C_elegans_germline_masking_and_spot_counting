"""Stage-resolved SYP-3 partition figure, ARTIFACT-FLOOR-CORRECTED (v2 of fig70).
Metric change (flagged): PC_specific = PC_dm / z-shift floor, where the z-shift control (granule mask
slid +1.2 um in z) measures the enrichment produced by axial bleed + PGL-correlated background alone.
Dividing by it anchors "no granule-specific enrichment" at 1.0, so the noHS baseline reads near 1.
Floors are per (sex, treat, stage) group: male floors exact from the analysis log; herm floors recovered
from fig70's plotted markers (validated on male knowns, max residual 0.011). Per-gonad values from
pc_stage_recovered.csv (recovered verbatim from the batch log). Stats recomputed on corrected values."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

CA = r"C:/Users/ryane/coloc_analysis"
OUT = f"{CA}/fig76_stage_zsh_normalized.png"
DESK = r"C:/Users/ryane/OneDrive/Desktop/germquant_figures"
STAGES = ["early", "mid", "late"]
XLAB = ["early\npachytene", "mid\npachytene", "late\npachytene"]

d = pd.read_csv(f"{CA}/pc_stage_recovered.csv")
fl = pd.read_csv(f"{CA}/zsh_floors.csv")
floors = {(r.sex, r.treat, r.stage): r.zsh_floor for r in fl.itertuples()}
d["floor"] = [floors[(r.sex, r.treat, r.stage)] for r in d.itertuples()]
d["PC_corr"] = d.PC_dm_raw / d.floor

fig, axes = plt.subplots(1, 2, figsize=(12.5, 6.2), sharey=True)
col = {"HS": "#c0392b", "noHS": "#2c6fbb"}
xs = np.arange(3)

for ax, sex, title in [(axes[0], "male", "MALES"), (axes[1], "herm", "hermaphrodites")]:
    ax.axhline(1.0, ls="--", c="#888", lw=1.2, zorder=0)
    for treat in ["noHS", "HS"]:
        sub = d[(d.sex == sex) & (d.treat == treat)]
        m = [sub[sub.stage == s].PC_corr.mean() for s in STAGES]
        se = [sub[sub.stage == s].PC_corr.std() / max(1, np.sqrt(sub[sub.stage == s].PC_corr.notna().sum()))
              for s in STAGES]
        ax.errorbar(xs, m, yerr=se, color=col[treat], lw=2.4, marker="o", ms=8, capsize=4,
                    label=treat, zorder=4)
        for k, s in enumerate(STAGES):
            v = sub[sub.stage == s].PC_corr.dropna().values
            ax.scatter(np.full(len(v), xs[k]) + np.random.RandomState(2).uniform(-0.07, 0.07, len(v)),
                       v, s=18, color=col[treat], alpha=0.35, zorder=2)
    for k, s in enumerate(STAGES):
        # primary difference test = HS vs noHS on RAW PC (unaffected by the floor constants);
        # corrected-scale MW shown beneath in gray for transparency
        ar = d[(d.sex == sex) & (d.treat == "noHS") & (d.stage == s)].PC_dm_raw.dropna()
        br = d[(d.sex == sex) & (d.treat == "HS") & (d.stage == s)].PC_dm_raw.dropna()
        ac = d[(d.sex == sex) & (d.treat == "noHS") & (d.stage == s)].PC_corr.dropna()
        bc = d[(d.sex == sex) & (d.treat == "HS") & (d.stage == s)].PC_corr.dropna()
        if len(ar) >= 2 and len(br) >= 2:
            p = stats.mannwhitneyu(ar, br)[1]
            pc = stats.mannwhitneyu(ac, bc)[1]
            star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            ytop = max(np.nanmean(ac), np.nanmean(bc)) + 0.16
            ax.text(xs[k], ytop, f"{star} p={p:.3f}", ha="center", va="bottom", fontsize=8.5,
                    color="#111" if p < 0.05 else "#888")
            ax.text(xs[k], ytop - 0.008, f"\n(corr {pc:.3f})", ha="center", va="top", fontsize=7,
                    color="#999")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xticks(xs); ax.set_xticklabels(XLAB)
    ax.set_xlim(-0.4, 2.4)
    ax.spines[["top", "right"]].set_visible(False)

axes[0].set_ylabel("granule-specific SYP-3 enrichment\n(PC / z-shift artifact floor)")
axes[0].set_ylim(0.7, 1.85)
axes[0].legend(loc="upper right", frameon=False, fontsize=9)
axes[1].text(0.03, 0.955, "dashed line = artifact floor (no granule-specific enrichment)",
             transform=axes[1].transAxes, fontsize=8, color="#666", va="top")
fig.suptitle("Heat shock drives granule-specific SYP-3 enrichment in males, peaking in mid-pachytene",
             fontsize=14, y=0.99)
foot = ("Each PC divided by its condition's z-shift control (mask slid +1.2 µm in z = axial-bleed + background floor), "
        "so 1.0 = no granule-specific enrichment; noHS sits near 1 as expected for baseline.\n"
        "p = Mann-Whitney HS vs noHS on raw PC (the valid difference test; floor constants don't affect it); "
        "gray = same test on corrected values (conservative: the HS floor contains real HS signal).\n"
        "Floors are group-level (male: exact from logs; herm: recovered from fig70, validated ±0.011). "
        "Per-gonad paired floors (more power) pending E:-drive/VPN reconnection. n=22 gonads.")
fig.text(0.5, 0.015, foot, ha="center", va="bottom", fontsize=8.3, color="#333")
fig.subplots_adjust(bottom=0.2, top=0.9, wspace=0.08)
fig.savefig(OUT, dpi=200)
import os, shutil
if os.path.isdir(DESK):
    shutil.copy(OUT, os.path.join(DESK, os.path.basename(OUT)))

print("corrected group means (PC / floor):")
for sex in ["male", "herm"]:
    for treat in ["noHS", "HS"]:
        sub = d[(d.sex == sex) & (d.treat == treat)]
        vals = [round(sub[sub.stage == s].PC_corr.mean(), 3) for s in STAGES]
        print(f"  {sex:5s} {treat:4s}: early={vals[0]} mid={vals[1]} late={vals[2]}")
print("WROTE", OUT)
