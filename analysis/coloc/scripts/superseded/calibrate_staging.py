"""Calibrate the TZ band + landmark rule against the 5 human-verified pachytene starts, jointly.
Tuning on one gonad is how the earlier version fooled itself (male_008 was fine, the other four were off
by 8-23 rows), so this sweeps the parameters once and scores them on ALL five simultaneously.

Verified truth (independent visual review, 2026-08-22):
  HS_herm_011 16 | HS_male_008 15 | HS_male_013 14 | noHS_herm_004 25 | noHS_male_002 35
Uses the bug-fixed axis (component filter + volume-trend polarity) and the cached crops, so it is fast.
Reports the best (lo, hi, TZ_ESTABLISH, TZ_MAX) by mean absolute row error, plus a leave-one-out check
so we can see whether the fit generalises or is just memorising."""
import itertools
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import crescent_axis as ca

TRUTH = {
    "20260708_ccw77_IF_pgl1_syp3_lmn1_HS_herm_011": 16,
    "20260708_ccw77_IF_pgl1_syp3_lmn1_HS_male_008": 15,
    "20260708_ccw77_IF_pgl1_syp3_lmn1_HS_male_013": 14,
    "20260708_ccw77_IF_pgl1_syp3_lmn1_noHS_herm_004": 25,
    "20260708_ccw77_IF_pgl1_syp3_lmn1_noHS_male_002": 35,
}


def prep(iid):
    """features -> component filter -> axis -> orientation -> rows (all bug-fixed), band-independent."""
    dapi, syp, lam, lab_c, fb = ca.load_crops(iid)
    df = ca.nucleus_features(lam, dapi, syp)
    df, comp = ca.largest_component(df)
    s, L = ca.geodesic_axis(df)
    df["s_um"] = s
    df = df.sort_values("s_um").reset_index(drop=True)
    flip, ori = ca.orient(df, df.s_um.to_numpy(), L)
    if flip:
        df["s_um"] = L - df.s_um
        df = df.sort_values("s_um").reset_index(drop=True)
    df["s_norm"] = df.s_um / L
    df["score"] = df.gap_int / df.thick.clip(lower=1e-6)
    core = df[(df.s_norm > 0.35) & (df.s_norm < 0.65)]
    df["p_ref"] = float(core.score.median())
    row_w = float(2 * df.r_eq.median())
    df["row"] = (df.s_um // row_w).astype(int)
    return df, dict(flip=flip, kept=comp["kept_frac"], slope=ori["vol_slope_um3_per_um"], row_w=row_w)


def detect(df, lo, hi, establish, tzmax, persist=2):
    p = df.p_ref.iloc[0]
    tz = (df.score >= lo * p) & (df.score < hi * p) & (df.s_norm < 0.5) & (df.vol_env < 2.2 * df.vol_env.median())
    rt = df.assign(tz=tz).groupby("row").tz.agg(["sum", "size"]).reset_index()
    rt.columns = ["row", "n_tz", "n"]
    start = None
    for i in range(len(rt) - persist + 1):
        if all(rt.n_tz.iloc[i + k] >= establish for k in range(persist)):
            start = i
            break
    if start is None:
        hits = [i for i in range(len(rt)) if rt.n_tz.iloc[i] >= max(1, establish - 1)]
        start = hits[0] if hits else None
    if start is None:
        return None
    for i in range(start, len(rt) - persist + 1):
        if all(rt.n_tz.iloc[i + k] <= tzmax for k in range(persist)):
            return int(rt.row.iloc[i])
    return None


print("preparing 5 verified gonads (bug-fixed axis)...", flush=True)
DATA = {}
for iid in TRUTH:
    df, info = prep(iid)
    DATA[iid] = df
    print(f"  {iid[-16:]:18s} n={len(df):4d} kept={info['kept']} flip={info['flip']} "
          f"slope={info['slope']} row_w={info['row_w']:.2f}", flush=True)

grid = list(itertools.product([0.20, 0.25, 0.30, 0.35, 0.40], [0.55, 0.62, 0.72, 0.82],
                              [2, 3, 4], [1, 2, 3]))
res = []
for lo, hi, est, tzmax in grid:
    if hi <= lo:
        continue
    errs, miss = [], 0
    for iid, truth in TRUTH.items():
        got = detect(DATA[iid], lo, hi, est, tzmax)
        if got is None:
            miss += 1; errs.append(99)
        else:
            errs.append(abs(got - truth))
    res.append((np.mean(errs), max(errs), miss, lo, hi, est, tzmax))
res.sort()
print(f"\ntop 8 of {len(res)} parameter sets (mean |error| in rows, over all 5 gonads):")
print(f"{'mae':>5s} {'max':>4s} {'miss':>4s}  lo    hi   establish tzmax")
for r in res[:8]:
    print(f"{r[0]:5.2f} {r[1]:4d} {r[2]:4d}  {r[3]:.2f}  {r[4]:.2f}   {r[5]}        {r[6]}")

best = res[0]
_, _, _, lo, hi, est, tzmax = best
print(f"\nBEST: lo={lo} hi={hi} establish={est} tzmax={tzmax}")
print(f"{'gonad':18s} {'truth':>5s} {'got':>5s} {'err':>5s}")
for iid, truth in TRUTH.items():
    got = detect(DATA[iid], lo, hi, est, tzmax)
    print(f"{iid[-16:]:18s} {truth:5d} {str(got):>5s} {'' if got is None else got - truth:>5}")

# leave-one-out: fit on 4, test on the held-out one
print("\nleave-one-out generalisation:")
loo = []
for held in TRUTH:
    sub = {k: v for k, v in TRUTH.items() if k != held}
    scored = []
    for lo2, hi2, est2, tz2 in grid:
        if hi2 <= lo2:
            continue
        e = []
        for iid, truth in sub.items():
            g = detect(DATA[iid], lo2, hi2, est2, tz2)
            e.append(99 if g is None else abs(g - truth))
        scored.append((np.mean(e), lo2, hi2, est2, tz2))
    scored.sort()
    _, lo2, hi2, est2, tz2 = scored[0]
    g = detect(DATA[held], lo2, hi2, est2, tz2)
    err = None if g is None else abs(g - TRUTH[held])
    loo.append(99 if err is None else err)
    print(f"  held out {held[-16:]:18s} fit=({lo2},{hi2},{est2},{tz2}) -> got {g} truth {TRUTH[held]} err {err}")
print(f"LOO mean |error| = {np.mean(loo):.2f} rows  ({np.mean(loo) * 4.1:.0f} um)")
