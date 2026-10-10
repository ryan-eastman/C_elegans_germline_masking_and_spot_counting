"""An exact, fast replica of SpotMAX's minimum-distance peak filter, used by the SpotMAX spot stage.

Why this exists. SpotMAX (1.3.1, the pinned version) finds spot candidates with
``skimage.feature.peak_local_max`` and then thins them with
``spotmax.filters.filter_valid_points_min_distance``: sort the peaks by descending intensity, walk
them in that order, keep a peak unless an earlier kept peak lies inside its ellipsoid
(``sum(((p_j - p_i) / min_distance) ** 2) < 1``) and drop the later peaks inside the ellipsoid of every
kept one. SpotMAX does this with a masked numpy array and a full pass over ALL peaks for every kept
peak, so the cost grows with the square of the number of peaks. A normal RAD-51 gonad has tens of
thousands of peaks and the pass takes seconds to minutes; a flooded one (n2rad51_noHS_male_002:
479,901 peaks) needs hours, before the candidate cap of the spot stage is ever reached.
docs/PERFORMANCE.md has the profile.

What this module does. ``filter_valid_points_min_distance`` below returns exactly what SpotMAX's
function returns (same rows, same order, same dtype): the same intensity sort call (so ties are
broken the same way), the same greedy walk, and the same floating point expression for the
ellipsoid test, evaluated on the same integer difference vectors. Only the search for the peaks near
a kept peak changes: a KD-tree on the coordinates divided by ``min_distance`` proposes every pair
closer than ``1 + _KD_MARGIN`` in those units (a superset of the pairs the test can accept, because
the rounding error of the scaled coordinates is many orders of magnitude below the margin), and the
SpotMAX test then decides each proposed pair. Inputs outside the domain this was verified on (points
that are not a 2D numpy array of signed integers, an intensity array of another length, a radius
that is not a positive finite number, coordinates too large for the margin argument) go to SpotMAX's
own function unchanged.

``fast_min_distance_filter()`` installs the replica in ``spotmax.filters`` for the duration of a
``with`` block (the SpotMAX spot stage wraps its ``spot_detection`` call in it) and puts SpotMAX's
function back afterwards. Nested and concurrent blocks are counted, so the original is restored only
when the last one exits. tests/test_peak_filter.py checks the replica against the installed SpotMAX
function on random inputs with ties, duplicates and boundary distances.
"""
from __future__ import annotations

import contextlib
import logging
import threading
import time

import numpy as np

log = logging.getLogger(__name__)

# KD-tree query radius margin, in units of min_distance. The rounding error of a scaled coordinate
# (|q| * 2.2e-16, with |q| capped by _MAX_SCALED below) is below 1e-7, far inside this margin.
_KD_MARGIN = 1e-6
_MAX_SCALED = 1e8

_lock = threading.Lock()
_depth = 0
_saved = None


def filter_valid_points_min_distance(points, min_distance, intensities=None,
                                     return_valid_points_mask=False, debug=False, *, _fallback=None):
    """Drop-in replica of ``spotmax.filters.filter_valid_points_min_distance`` (SpotMAX 1.3.1).

    ``points`` is an (N, D) array of peak coordinates, ``min_distance`` the D ellipsoid radii (or one
    number), ``intensities`` the N peak intensities (the points are walked brightest first when
    given, in input order otherwise). Returns the kept points in walk order, and with
    ``return_valid_points_mask`` also the boolean keep mask in walk order, exactly as SpotMAX does.
    ``debug`` is accepted for signature compatibility and ignored, as in SpotMAX. ``_fallback`` is
    the function used for inputs outside the verified domain (default: SpotMAX's own)."""
    if not _fast_path_ok(points, min_distance, intensities):
        if _fallback is None:
            from spotmax.filters import filter_valid_points_min_distance as _fallback
        return _fallback(points, min_distance, intensities=intensities,
                         return_valid_points_mask=return_valid_points_mask, debug=debug)
    t0 = time.perf_counter()
    md = np.asarray(min_distance)
    if intensities is not None:
        # the SpotMAX sort, call for call: argsort (default kind) reversed, so equal intensities
        # come out in the same order as in SpotMAX
        sorting_idxs_descending = np.flip(intensities.argsort())
        points = points[sorting_idxs_descending]
    valid = _greedy_keep_mask(points, md)
    valid_points = points[valid]
    log.debug("minimum-distance peak filter: %d candidates in, %d kept, %.2f s",
              len(points), len(valid_points), time.perf_counter() - t0)
    if return_valid_points_mask:
        return valid_points, valid
    return valid_points


def _fast_path_ok(points, min_distance, intensities) -> bool:
    """True when the inputs are inside the domain the replica was verified on."""
    if not isinstance(points, np.ndarray) or points.ndim != 2 or points.dtype.kind != "i":
        return False
    if intensities is not None and (not isinstance(intensities, np.ndarray) or intensities.ndim != 1
                                    or len(intensities) != len(points)):
        return False
    md = np.asarray(min_distance)
    if md.ndim > 1 or md.size not in (1, points.shape[1]) or md.dtype.kind not in "iuf":
        return False
    if not (np.all(np.isfinite(md)) and np.all(md > 0)):
        return False
    # SpotMAX divides through a masked array, which masks any entry with |a| * tiny >= |b|; with
    # integer coordinates that can only happen for a vanishing radius, so require a sane one
    if float(np.min(md)) < 1e-6:
        return False
    if len(points) and float(np.abs(points).max()) / float(np.min(md)) > _MAX_SCALED:
        return False
    return True


def _greedy_keep_mask(pts: np.ndarray, md: np.ndarray) -> np.ndarray:
    """The keep mask of SpotMAX's greedy walk over ``pts`` (already in walk order)."""
    from scipy.spatial import cKDTree

    n = len(pts)
    valid = np.ones(n, dtype=bool)
    if n < 2:
        return valid
    tree = cKDTree(pts / md)
    pairs = tree.query_pairs(r=1.0 + _KD_MARGIN, output_type="ndarray")   # rows (i, j) with i < j
    if len(pairs) == 0:
        return valid
    first, second = pairs[:, 0], pairs[:, 1]
    # SpotMAX's ellipsoid test on the difference vectors (np.square of the division, summed over
    # the coordinate axis of a C-contiguous (M, D) array, strict < 1): the same expression it
    # evaluates on its (N, D) array, so a boundary pair is decided the same way
    diff = pts[second] - pts[first]
    close = np.square(diff / md).sum(axis=1) < 1
    first, second = first[close], second[close]
    if len(first) == 0:
        return valid
    order = np.argsort(first, kind="stable")
    first, second = first[order], second[order]
    starts = np.searchsorted(first, np.arange(n + 1)).tolist()
    later = second.tolist()
    keep = bytearray(b"\x01") * n
    # walk in order: a kept point drops every later point inside its ellipsoid; a dropped point
    # never drops anything (in SpotMAX it is masked out of every later pass)
    for i in np.unique(first).tolist():
        if keep[i]:
            for j in later[starts[i]:starts[i + 1]]:
                keep[j] = 0
    return np.frombuffer(bytes(keep), dtype=np.uint8).astype(bool)


@contextlib.contextmanager
def fast_min_distance_filter():
    """Inside the block, ``spotmax.filters.filter_valid_points_min_distance`` is the exact replica
    above (with SpotMAX's own function as the fallback); afterwards SpotMAX's function is back.
    Nested and concurrent blocks are reference counted. Without SpotMAX, or with a SpotMAX that has
    no such function, the block changes nothing."""
    global _depth, _saved
    try:
        import spotmax.filters as smf
    except Exception:  # noqa: BLE001 - no SpotMAX: nothing to replace; the caller's own import reports it
        smf = None
    installed = False
    if smf is not None:
        with _lock:
            if _depth == 0:
                original = getattr(smf, "filter_valid_points_min_distance", None)
                if original is not None:
                    _saved = original
                    smf.filter_valid_points_min_distance = _make_replica(original)
                    installed = True
            else:
                installed = True
            if installed:
                _depth += 1
    if not installed:
        if smf is not None:
            log.debug("spotmax.filters has no filter_valid_points_min_distance; not replaced")
        yield
        return
    try:
        yield
    finally:
        with _lock:
            _depth -= 1
            if _depth == 0:
                smf.filter_valid_points_min_distance = _saved
                _saved = None


def _make_replica(original):
    """The function installed in spotmax.filters: the replica, falling back to ``original``."""

    def filter_valid_points_min_distance_replica(points, min_distance, intensities=None,
                                                 return_valid_points_mask=False, debug=False):
        return filter_valid_points_min_distance(
            points, min_distance, intensities=intensities,
            return_valid_points_mask=return_valid_points_mask, debug=debug, _fallback=original)

    filter_valid_points_min_distance_replica.__wrapped__ = original
    return filter_valid_points_min_distance_replica
