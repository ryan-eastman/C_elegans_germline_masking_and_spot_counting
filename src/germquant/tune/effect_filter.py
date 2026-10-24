"""The SpotMAX effect-size filter, applied after a run made with ``effect_size_min`` 0 (design 6.8).

``spots.detect.detect_spots_spotmax`` applies ``effect_size_min`` as its very last step, after the
detection, the features and the z-split merge, as a plain row filter: a spot whose effect size is below
the minimum is dropped, a missing effect size counts as passing, and nothing is filtered when no spot has
an effect size or the minimum is 0 (the chain with the optional pack keys, spots/spotmax_options.py, ends
the same way). Parameter sets that differ only in ``effect_size_min`` therefore share everything before
that filter, so the spot tuner runs SpotMAX once per group of such sets, with ``effect_size_min`` 0, and
applies each set's minimum afterwards with ``filter_effect_size`` (``germquant.tune.spots``).

This module re-implements the filter of the RAD-51 harness (benchmark/tuning.py), which is sealed (design
rule R8): it does not import it. tests/test_tune_effect_filter.py imports both and proves them equal on
synthetic data, and tests/test_tune_spots.py proves that a filtered run equals a direct SpotMAX run at
each effect size.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def filter_effect_size(per_spot: pd.DataFrame, per_nuc: pd.DataFrame, effect_size_min: float):
    """Apply SpotMAX's effect-size filter to the tables of a run made with ``effect_size_min`` 0, exactly
    as ``detect_spots_spotmax`` applies it before building its tables.

    Rows whose effect size (a missing one counting as infinite) is below the minimum are dropped, in order,
    when the minimum is above 0 and at least one spot has an effect size. Then the spot ids are numbered
    again from 0, and the per-nucleus counts and mean spot intensities are computed again over every
    nucleus of ``per_nuc``, with the expressions of ``detect._finish_tables`` (one row per nucleus, a zero
    count being a real zero). Returns (per_spot, per_nucleus); the inputs are not changed."""
    from germquant.spots.backends import PER_NUC_COLS

    e = float(effect_size_min)
    ps = per_spot
    if e > 0 and len(ps) and ps["effect_size"].notna().any():
        keep = ps["effect_size"].astype(float).fillna(np.inf).to_numpy() >= e
        ps = ps.loc[keep]
    ps = ps.reset_index(drop=True)
    if "spot_id" in ps.columns:
        ps = ps.assign(spot_id=np.arange(len(ps)))
    if not len(per_nuc):
        return ps, per_nuc.copy()
    marker = per_nuc["marker"].iloc[0]
    detector = per_nuc["detector"].iloc[0]
    nuc_ids = [int(v) for v in per_nuc["nucleus_id"]]
    counts = ps.groupby("nucleus_id").size()
    means = ps.groupby("nucleus_id")["spot_mean_intensity"].mean()
    pn = pd.DataFrame({
        "nucleus_id": nuc_ids,
        "marker": marker,
        "n_spots": [int(counts.get(n, 0)) for n in nuc_ids],
        "mean_spot_intensity": [float(means.get(n, np.nan)) for n in nuc_ids],
        "detector": detector,
    }, columns=PER_NUC_COLS)
    return ps, pn
