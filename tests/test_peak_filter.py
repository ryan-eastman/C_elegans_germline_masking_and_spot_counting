"""The exact replica of SpotMAX's minimum-distance peak filter (germquant.spots.peak_filter):
- it returns what spotmax.filters.filter_valid_points_min_distance returns (same rows, same order,
  same dtype, same keep mask) on random inputs with tied intensities, duplicate points, 2D and 3D
  points, scalar / list / float32 radii, and on pairs that sit exactly on the ellipsoid boundary;
- inputs outside the verified domain go to the fallback (SpotMAX's own function) untouched;
- fast_min_distance_filter installs it in spotmax.filters only inside the block, counts nested
  blocks, and restores SpotMAX's function on exit and on error;
- the SpotMAX spot stage gives the same tables with and without it on a flooded synthetic image;
- the flood guard's warning carries the candidate count and adds no column.
"""
from __future__ import annotations

import contextlib
import logging

import numpy as np
import pandas as pd
import pytest

from germquant.spots import peak_filter
from germquant.spots.detect import _cap_candidates
from germquant.spots.peak_filter import fast_min_distance_filter, filter_valid_points_min_distance

RADII = [np.array([1.5, 2.76923077, 2.76923077]), np.array([1.0, 1.0, 1.0]), np.array([2.0, 2.0, 2.0]),
         np.array([1.5, 1.5, 1.5]), np.array([np.sqrt(2), np.sqrt(2), 1.0]), np.array([0.7, 3.3, 2.2]),
         np.array([1.5, 2.76923077, 2.76923077], dtype=np.float32), 2.5, [1, 2, 2]]


def _spotmax_filter():
    smf = pytest.importorskip("spotmax.filters")
    fn = smf.filter_valid_points_min_distance
    return getattr(fn, "__wrapped__", fn)          # SpotMAX's own function even inside a block


def _same(a, b):
    assert a.dtype == b.dtype and a.shape == b.shape
    assert np.array_equal(a, b)


def test_replica_matches_spotmax_on_random_inputs():
    original = _spotmax_filter()
    rng = np.random.default_rng(20251105)
    for trial in range(120):
        ndim = 2 if trial % 5 == 0 else 3
        n = int(rng.integers(0, 300))
        box = int(rng.integers(3, 30))
        pts = rng.integers(0, box, size=(n, ndim)).astype(np.int64)
        if trial % 7 == 0 and n:
            pts = np.concatenate([pts, pts[: n // 3]])          # duplicate points
        kind = trial % 3
        inten = (rng.integers(0, 4, size=len(pts)).astype(np.float32) if kind == 0   # heavy ties
                 else rng.random(len(pts)).astype(np.float32) if kind == 1 else None)
        radii = RADII[trial % len(RADII)]
        if ndim == 2 and np.ndim(radii) and len(radii) == 3:
            radii = np.asarray(radii)[1:]
        _same(filter_valid_points_min_distance(pts, radii, intensities=inten),
              original(pts, radii, intensities=inten))
        got_pts, got_mask = filter_valid_points_min_distance(pts, radii, intensities=inten,
                                                             return_valid_points_mask=True)
        ref_pts, ref_mask = original(pts, radii, intensities=inten, return_valid_points_mask=True)
        _same(got_pts, ref_pts)
        _same(got_mask, ref_mask)


def test_replica_matches_spotmax_on_boundary_pairs():
    # offsets that land exactly on the ellipsoid (sum == 1, not < 1) or a hair inside it
    original = _spotmax_filter()
    base = np.array([[10, 10, 10]], np.int64)
    offsets = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [2, 0, 0], [1, 1, 0], [1, 1, 1], [0, 2, 1]], np.int64)
    pts = np.concatenate([base, base + offsets, base - offsets])
    inten = np.linspace(10, 1, len(pts)).astype(np.float32)
    for radii in ([1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [np.sqrt(2)] * 3, [np.sqrt(3)] * 3, [1.5, 2.0, 2.0],
                  [2.0000000001, 1.0, 1.0], [1.9999999999, 1.0, 1.0]):
        radii = np.asarray(radii)
        _same(filter_valid_points_min_distance(pts, radii, intensities=inten),
              original(pts, radii, intensities=inten))


def test_replica_falls_back_outside_the_verified_domain():
    calls = []

    def fallback(points, min_distance, intensities=None, return_valid_points_mask=False, debug=False):
        calls.append(np.asarray(points).dtype.kind)
        return "fallback"

    radii = np.array([1.5, 2.0, 2.0])
    ints = np.zeros((4, 3), np.int64)
    for pts, rad, inten in (
        (ints.astype(float), radii, None),                    # float coordinates
        (ints.astype(np.uint16), radii, None),                # unsigned coordinates (wrap on subtraction)
        (ints, radii, np.ones(3, np.float32)),                # intensity array of another length
        (ints, np.array([0.0, 2.0, 2.0]), None),              # a zero radius
        (ints, np.array([np.nan, 2.0, 2.0]), None),           # a non-finite radius
        (ints, np.array([1.0, 2.0]), None),                   # radii that do not match the points
        (ints.tolist(), radii, None),                         # not an array
    ):
        assert filter_valid_points_min_distance(pts, rad, intensities=inten, _fallback=fallback) == "fallback"
    assert len(calls) == 7
    # inside the domain the fallback is never consulted
    out = filter_valid_points_min_distance(ints, radii, intensities=np.ones(4, np.float32), _fallback=fallback)
    assert len(calls) == 7 and out.shape == (1, 3)


def test_fast_filter_block_installs_counts_and_restores():
    smf = pytest.importorskip("spotmax.filters")
    original = smf.filter_valid_points_min_distance
    with fast_min_distance_filter():
        inside = smf.filter_valid_points_min_distance
        assert inside is not original and inside.__wrapped__ is original
        with fast_min_distance_filter():                       # nested: still the same replica
            assert smf.filter_valid_points_min_distance is inside
        assert smf.filter_valid_points_min_distance is inside  # inner exit does not restore
        pts = np.array([[0, 0, 0], [0, 1, 0], [5, 5, 5]], np.int64)
        _same(inside(pts, np.array([1.5, 2.0, 2.0]), intensities=np.array([3, 2, 1], np.float32)),
              original(pts, np.array([1.5, 2.0, 2.0]), intensities=np.array([3, 2, 1], np.float32)))
    assert smf.filter_valid_points_min_distance is original
    with pytest.raises(RuntimeError):
        with fast_min_distance_filter():
            raise RuntimeError("boom")
    assert smf.filter_valid_points_min_distance is original and peak_filter._depth == 0


def _flooded_volume(seed=3):
    """Three box nuclei filled with a bright smooth texture: with the mean threshold SpotMAX proposes
    hundreds of peak candidates, its minimum-distance filter drops some of them, and the flood guard
    caps the rest (max_spot_candidates=200 in the test)."""
    from scipy import ndimage as ndi

    rng = np.random.default_rng(seed)
    shape = (10, 64, 96)
    labels = np.zeros(shape, np.int32)
    labels[1:9, 4:60, 4:30] = 1
    labels[1:9, 4:60, 34:62] = 2
    labels[1:9, 4:60, 66:92] = 3
    img = rng.normal(100.0, 5.0, shape)
    texture = ndi.gaussian_filter(rng.normal(0.0, 1.0, shape), (0.7, 1.2, 1.2))
    img[labels > 0] += 300.0 + 2000.0 * texture[labels > 0]
    return np.clip(img, 0, 65535).astype(np.uint16), labels


def test_spotmax_stage_identical_with_and_without_the_replica(monkeypatch, caplog):
    pytest.importorskip("spotmax")
    from germquant.spots.detect import detect_spots_spotmax

    img, labels = _flooded_volume()
    params = dict(thresholding_method="threshold_mean", effect_size_min=3.0, max_spot_candidates=200)
    spacing = (0.2, 0.108333333333333, 0.108333333333333)
    seen = []
    real_greedy = peak_filter._greedy_keep_mask

    def spy(pts, md):
        keep = real_greedy(pts, md)
        seen.append((len(pts), int(keep.sum())))
        return keep

    monkeypatch.setattr(peak_filter, "_greedy_keep_mask", spy)
    with caplog.at_level(logging.WARNING, logger="germquant.spots.detect"):
        fast = detect_spots_spotmax(img, labels, spacing, **params)
    assert len(seen) == 1 and seen[0][0] > 300 and seen[0][1] < seen[0][0]   # the replica ran and dropped peaks
    assert "FLOODED" in caplog.text                                         # and the cap was hit
    monkeypatch.setattr(peak_filter, "fast_min_distance_filter", contextlib.nullcontext)
    slow = detect_spots_spotmax(img, labels, spacing, **params)
    assert len(seen) == 1                                       # SpotMAX's own filter ran the second time
    pd.testing.assert_frame_equal(fast[0], slow[0], check_exact=True)
    pd.testing.assert_frame_equal(fast[1], slow[1], check_exact=True)
    assert len(fast[0]) > 0


def test_flood_guard_warns_with_the_candidate_count(caplog):
    img = np.zeros((1, 1, 10), np.float32)
    img[0, 0, :] = np.arange(10)
    df = pd.DataFrame({"z": [0] * 10, "y": [0] * 10, "x": list(range(10)), "nucleus_id": [1] * 10})
    with caplog.at_level(logging.WARNING, logger="germquant.spots.detect"):
        out = _cap_candidates(df, img, cap=3)
    assert list(out.columns) == list(df.columns) and len(out) == 3
    msg = " ".join(r.getMessage() for r in caplog.records)
    assert "FLOODED: 10 spot candidates" in msg and "max_spot_candidates=3" in msg and "dropping 7" in msg
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="germquant.spots.detect"):
        _cap_candidates(df, img, cap=10)                        # at the cap: no warning
    assert not caplog.records
