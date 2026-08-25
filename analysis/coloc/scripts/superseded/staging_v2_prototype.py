"""PROTOTYPE replacement staging decision layer (NOT production-ready - see limits below).

Replaces the old first-crossing count rule ("first row with <=2 TZ nuclei"), which had no hysteresis,
used an absolute count instead of a fraction, and relied on gap_int/thick (a feature that runs BACKWARDS
on some gonads). Instead: build a smoothed per-row profile of SC loading (fraction of nuclei whose
syp_io exceeds the midpoint between the gonad's own low and high syp modes) and take pachytene start as
the logistic midpoint - the row where the profile first crosses halfway from its distal baseline to its
pachytene plateau, and stays there.

MEASURED PERFORMANCE against 7 human-verified pachytene starts (5 ccw77 + 2 N2 dry-ice):
    old count rule : leave-one-out MAE 12.6 rows (~52 um)
    this prototype : leave-one-out MAE  6.9 rows, all-data MAE 4.9 rows
    6 of 7 gonads land within 3 rows; HS_male_008 fails by +20.
KNOWN LIMIT (do not paper over): HS_male_008 has weak SYP-3 staining, so SC loading does not rise until
row 35 while chromatin morphology says 15. All four candidate features AGREE on 35, so the failure
cannot be caught by a feature-disagreement abstention rule (spread-vs-error correlation is -0.20, i.e.
none). Any gonad with weak/delayed SC staining will be staged late and will not announce itself.
=> Not accurate enough for early/mid/late thirds (that needs <=2-3 rows). Use only with the strip-image
   check, or fix the staining first.
STILL REQUIRED before this is trustworthy: the axis bug (geodesic_axis returns an iso-distance FIELD,
not a traced path, so branched/fused gonads collapse two arms into one row - 53% foreign nuclei on N2
noHS_male_001). Fixing the decision layer alone is not enough."""
import sys

import numpy as np
import pandas as pd

WIN = 5
FRAC = 0.5


def sc_profile(d, win=WIN):
    """per-row fraction of nuclei with SC loaded, smoothed."""
    s = d.syp_io.dropna()
    lo, hi = np.percentile(s, [15, 85])
    thr = 0.5 * (lo + hi)                      # midpoint between the gonad's own low/high syp modes
    p = d.assign(x=d.syp_io > thr).groupby("row").x.mean()
    p = p.reindex(range(int(d.row.max()) + 1)).interpolate().bfill().ffill().fillna(0)
    return p.rolling(win, center=True, min_periods=1).mean()


def pachytene_start(sm, frac=FRAC):
    """logistic midpoint: first sustained crossing of baseline + frac*(plateau - baseline)."""
    v = sm.values
    n = len(v)
    base = np.nanpercentile(v[:max(3, n // 8)], 50)
    plat = np.nanpercentile(v[n // 3:], 75)
    if not np.isfinite(plat) or plat <= base:
        return None, dict(base=base, plateau=plat, reason="no rise detected")
    tgt = base + frac * (plat - base)
    for i in range(n - 3):
        if v[i] < tgt <= v[i + 1] and np.nanmean(v[i + 1:i + 4]) >= tgt:
            return i + 1, dict(base=round(float(base), 3), plateau=round(float(plat), 3),
                               target=round(float(tgt), 3))
    return None, dict(base=base, plateau=plat, reason="never crossed")


if __name__ == "__main__":
    iid = sys.argv[1]
    d = pd.read_csv(rf"C:/Users/ryane/coloc_analysis/staging/{iid}_nuclei.csv")
    sm = sc_profile(d)
    start, info = pachytene_start(sm)
    print(f"{iid}: pachytene_start_row={start}  {info}")
