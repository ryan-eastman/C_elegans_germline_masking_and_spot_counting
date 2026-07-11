"""Construction-truth recovery: plant a KNOWN SYP<->PGL-1 overlap gradient and confirm the coloc
metrics recover it. This is ground truth *by construction* — stronger evidence than the tiny unit
cases that the headline numbers (frac_granules_overlapping, Manders, the overlap null) are correct,
without needing any hand-annotated data.
"""
import numpy as np
import pytest

from germquant.coloc import colocalize

SP = (0.4, 0.2, 0.2)


def _scene(n_granules: int, n_overlap: int):
    """n_granules equal-size PGL-1 blobs along x; the first `n_overlap` are fully covered by the SC
    mask (planted coincidence), the rest are SC-free. Returns masks + intensity images + region."""
    Z, Y = 6, 20
    X = 4 * n_granules + 8
    gr = np.zeros((Z, Y, X), bool)
    labels = np.zeros((Z, Y, X), np.int32)
    sc = np.zeros((Z, Y, X), bool)
    for i in range(n_granules):
        x0 = 4 + i * 4
        sl = (slice(2, 4), slice(8, 12), slice(x0, x0 + 3))
        labels[sl] = i + 1
        gr[sl] = True
        if i < n_overlap:
            sc[sl] = True
    syp = sc.astype(np.float32) * 100.0
    pgl = gr.astype(np.float32) * 100.0
    region = np.ones((Z, Y, X), bool)
    return sc, gr, labels, syp, pgl, region


@pytest.mark.parametrize("n,k", [(10, 0), (10, 3), (10, 7), (10, 10)])
def test_overlap_fraction_tracks_planted(n, k):
    sc, gr, labels, syp, pgl, region = _scene(n, k)
    row, per = colocalize(sc, gr, labels, syp, pgl, region, SP, n_random=50, rng_seed=0)

    assert row["n_granules"] == n
    # the headline object-overlap fraction equals the planted fraction exactly
    assert abs(row["frac_granules_overlapping_sc"] - k / n) < 1e-9
    assert int(per["overlaps"].sum()) == k
    # Manders M2 (fraction of PGL-1 intensity inside the SC mask) = k/n for equal-size blobs
    assert abs(row["manders_m2"] - k / n) < 1e-6
    if k > 0:
        # SC is fully inside granules here -> all SYP intensity coincides with the PGL-1 mask
        assert abs(row["manders_m1"] - 1.0) < 1e-6


def test_overlap_null_is_significant_only_with_real_overlap():
    # no planted overlap -> SC mask empty -> null is undefined (NaN), overlap is a real zero
    row0, _ = colocalize(*_scene(10, 0), spacing=SP, n_random=100, rng_seed=1)
    assert row0["overlap_volume_um3"] == 0.0
    assert np.isnan(row0["overlap_pvalue"])

    # strong planted overlap -> observed far exceeds the translation null -> significant
    row1, _ = colocalize(*_scene(10, 8), spacing=SP, n_random=100, rng_seed=1)
    assert row1["overlap_pvalue"] <= 0.05
    assert row1["overlap_zscore"] > 2
