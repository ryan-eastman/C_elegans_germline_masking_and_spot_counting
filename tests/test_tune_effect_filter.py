"""tune.effect_filter re-implements the SpotMAX effect-size filter of the RAD-51 harness (design 6.8; P5
acceptance 11): on synthetic tables it equals benchmark.tuning.filter_effect_size, which this test imports
next to it. germquant.tune itself never imports the harness (the seal: tests/test_tune_seal.py)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from germquant.benchmark import tuning as harness
from germquant.spots.backends import PER_NUC_COLS, PER_SPOT_COLS
from germquant.tune import effect_filter


def _tables(seed: int, n_spots: int, n_nuclei: int, nan_frac: float = 0.15):
    rng = np.random.default_rng(seed)
    nuc = rng.integers(1, n_nuclei + 1, size=n_spots)
    es = rng.gamma(2.0, 2.0, size=n_spots)
    es[rng.random(n_spots) < nan_frac] = np.nan
    per_spot = pd.DataFrame({
        "spot_id": np.arange(n_spots), "nucleus_id": nuc, "marker": "foci",
        "z_um": rng.random(n_spots) * 5, "y_um": rng.random(n_spots) * 20, "x_um": rng.random(n_spots) * 20,
        "effect_size": es, "spot_mean_intensity": rng.random(n_spots) * 100 + 50, "detector": "spotmax",
    }, columns=PER_SPOT_COLS)
    ids = list(range(1, n_nuclei + 1))
    counts = per_spot.groupby("nucleus_id").size()
    per_nuc = pd.DataFrame({"nucleus_id": ids, "marker": "foci", "n_spots": [int(counts.get(i, 0)) for i in ids],
                            "mean_spot_intensity": np.nan, "detector": "spotmax"}, columns=PER_NUC_COLS)
    return per_spot, per_nuc


@pytest.mark.parametrize("seed, n_spots, n_nuclei", [(1, 60, 8), (2, 5, 9), (3, 0, 4), (4, 40, 1)])
def test_equals_benchmark_filter(seed, n_spots, n_nuclei):
    per_spot, per_nuc = _tables(seed, n_spots, n_nuclei)
    es = per_spot["effect_size"].dropna()
    minima = [0.0, -1.0, 0.5, 2.0, 7.5, (float(es.max()) + 1.0) if len(es) else 3.0]
    if len(es):
        minima.append(float(es.median()))
    before = (per_spot.copy(), per_nuc.copy())
    for e in minima:
        got = effect_filter.filter_effect_size(per_spot, per_nuc, e)
        want = harness.filter_effect_size(per_spot, per_nuc, e)
        pd.testing.assert_frame_equal(got[0], want[0])
        pd.testing.assert_frame_equal(got[1], want[1])
    pd.testing.assert_frame_equal(per_spot, before[0])          # the inputs are never changed
    pd.testing.assert_frame_equal(per_nuc, before[1])


def test_missing_effect_sizes_and_empty_tables():
    per_spot, per_nuc = _tables(5, 30, 6, nan_frac=0.0)
    all_nan = per_spot.assign(effect_size=np.nan)
    for e in (0.0, 3.0, 100.0):                    # no spot has an effect size: nothing is filtered
        got = effect_filter.filter_effect_size(all_nan, per_nuc, e)
        assert len(got[0]) == len(per_spot)
        pd.testing.assert_frame_equal(got[1], harness.filter_effect_size(all_nan, per_nuc, e)[1])
    one_nan = per_spot.copy()
    one_nan.loc[3, "effect_size"] = np.nan          # a missing effect size counts as passing
    got = effect_filter.filter_effect_size(one_nan, per_nuc, float(per_spot["effect_size"].max()) + 1)
    assert got[0]["spot_id"].tolist() == [0] and int(got[1]["n_spots"].sum()) == 1
    assert len(got[1]) == len(per_nuc)
    empty_nuc = per_nuc.iloc[0:0]
    got = effect_filter.filter_effect_size(per_spot, empty_nuc, 2.0)
    want = harness.filter_effect_size(per_spot, empty_nuc, 2.0)
    pd.testing.assert_frame_equal(got[0], want[0])
    pd.testing.assert_frame_equal(got[1], want[1])
