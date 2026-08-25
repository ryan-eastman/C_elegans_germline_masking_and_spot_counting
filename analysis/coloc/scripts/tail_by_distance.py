"""Two questions from the per-granule v3 table (granule_tail_pergranule_v3.parquet):
 1. RADIUS: fraction of granules "lit" (SYP-3 excess > 0.5 nuclear) as a function of distance from the
    nuclear envelope, per group. The right detection window keeps noHS herms at ~0 (no false positives)
    and HS males high.
 2. BLEED-THROUGH: in noHS herms granules carry no real SYP-3, so any dependence of excess on PGL
    brightness there is GFP->mCherry bleed. Fit excess = k * pgl (robust, Theil-Sen) on noHS herms, then
    subtract k * pgl from EVERY granule and re-tabulate the lit fractions."""
import numpy as np
import pandas as pd
from scipy import stats

CA = r"C:/Users/ryane/coloc_analysis"
g = pd.read_parquet(f"{CA}/granule_tail_pergranule_v3.parquet")
g["grp"] = g.sex + " " + g.treat
SHELLS = [(0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.5)]

print("=== fraction of granules with SYP-3 excess > 0.5 nuclear, by distance of closest approach to the envelope ===")
print(f"{'group':12s} " + " ".join(f"{a:.1f}-{b:.1f}um" for a, b in SHELLS) + "   (n granules per shell)")
for grp, s in g.groupby("grp"):
    fr, ns = [], []
    for a, b in SHELLS:
        m = (s.dmin_um >= a) & (s.dmin_um < b)
        fr.append((s.excess_n[m] > 0.5).mean() if m.sum() > 50 else np.nan)
        ns.append(int(m.sum()))
    print(f"{grp:12s} " + " ".join(f"{f:9.3f}" for f in fr) + "   " + str(ns))
print("\nwhere do the granules sit? median closest-approach distance (um) per group:",
      g.groupby("grp").dmin_um.median().round(2).to_dict())

# ---- bleed-through calibration on noHS herms ----
cal = g[(g.sex == "herm") & (g.treat == "noHS")]
# PGL brightness normalised by that gonad's SYP exposure ratio is not needed: excess_n is already in nuclear-SYP
# units and pgl_raw in counts; bleed adds counts to the 545 channel proportional to 477 counts x (exp545/exp477)
cal = cal.assign(pgl_scaled=cal.pgl_raw * cal.exp_545 / cal.exp_477 / cal.N_syp)   # bleed predictor in nuclear-SYP units
slope, intercept, lo, hi = stats.theilslopes(cal.excess_n, cal.pgl_scaled)
print(f"\n=== bleed-through calibration (noHS herms, n={len(cal)} granules) ===")
print(f"excess_n = {intercept:+.4f} + {slope:.4f} * PGL(scaled)   [95% CI slope {lo:.4f}..{hi:.4f}]")
print(f"Spearman excess vs PGL in noHS herms: rho={stats.spearmanr(cal.excess_n, cal.pgl_scaled)[0]:.3f}")
g["pgl_scaled"] = g.pgl_raw * g.exp_545 / g.exp_477 / g.N_syp
g["excess_corr"] = g.excess_n - slope * g.pgl_scaled

print("\n=== per-gonad lit fraction (>0.5), raw vs bleed-corrected ===")
rows = []
for iid, s in g.groupby("image_id"):
    rows.append({"short": iid[-14:], "grp": s.grp.iloc[0], "raw": (s.excess_n > 0.5).mean(),
                 "corrected": (s.excess_corr > 0.5).mean(), "p90_raw": s.excess_n.quantile(.9), "p90_corr": s.excess_corr.quantile(.9)})
t = pd.DataFrame(rows).sort_values("grp")
print(t.round(3).to_string(index=False))
print("\n=== group means, corrected ===")
print(t.groupby("grp")[["raw", "corrected", "p90_raw", "p90_corr"]].mean().round(3).to_string())
for col in ["corrected", "p90_corr"]:
    a = t[t.grp == "male noHS"][col]; b = t[t.grp == "male HS"][col]; h = t[t.grp == "herm HS"][col]
    print(f"{col}: male HS vs noHS p={stats.mannwhitneyu(a,b)[1]:.3f}; HS male vs HS herm p={stats.mannwhitneyu(b,h)[1]:.3f}")
t.to_csv(f"{CA}/granule_tail_bleedcorrected.csv", index=False)
