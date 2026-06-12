"""Germline axis linearization: assign each nucleus a distal->proximal position.

Method (Libuda-lab Gonad Linearization, reimplemented): fit the germline's central axis as
a smooth **principal curve**, project each nucleus centroid onto it, and report position as
**arc length** along that curve (0 = distal tip, 1 = proximal end after normalization).

Why a curve and not a line: the *C. elegans* gonad arm is U-shaped. A single straight PCA
axis (the v0 here) cuts the corner — the two arms collapse onto each other, positions become
a chord instead of arc length, and a healthy curved gonad is spuriously flagged low-confidence.
We instead:
  1. order the nuclei along the gonad with the diameter path of their minimum spanning tree
     (robust to curvature — a U-shape's backbone is a single long path, not two overlapping
     PCA ranks);
  2. fit a smoothing spline through that backbone (falls back to the smoothed polyline, then
     to straight PCA, if the spline can't be fit);
  3. assign every nucleus the cumulative **arc length** of its nearest point (foot-point) on
     the curve.

Orientation (which end is distal) is a heuristic: the distal tip is densely packed with small
mitotic nuclei, so we orient so the small-nucleus end = position 0. FLAG: confirm orientation
against a known landmark before trusting absolute direction.
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def linearize_germline(
    nuclei: pd.DataFrame,
    *,
    confidence_min: float = 0.60,
) -> tuple[pd.DataFrame, float, list[str]]:
    """Add axis_position_norm (0..1) + axis_position_um. Return (df, confidence, flags).

    confidence is an R²-like goodness-of-fit of the principal curve (1 = nuclei lie tightly on
    the fitted axis, 0 = no coherent axis); below ``confidence_min`` a flag routes the gonad to
    a manual axis-draw. Unlike the old PCA variance-ratio, this does NOT penalise curvature.
    """
    flags: list[str] = []
    df = nuclei.copy()
    if len(df) < 5:
        df["axis_position_norm"] = float("nan")
        df["axis_position_um"] = float("nan")
        return df, 0.0, ["axis:too_few_nuclei"]

    pts = df[["centroid_z_um", "centroid_y_um", "centroid_x_um"]].to_numpy(dtype=float)
    center = pts.mean(axis=0)
    scale = float(np.std(pts - center)) or 1.0          # single isotropic scale (keeps shape)
    sp = (pts - center) / scale                         # work in scaled space; ×scale back to µm

    curve, method = _fit_principal_curve(sp)            # (S,3) dense samples or None
    if curve is None:
        proj, conf, flag = _pca_fallback(sp)
        if flag:
            flags.append(flag)
        arc_um = proj * scale
        method = "pca_line"
    else:
        # arc length along the dense curve, and each nucleus's foot-point arc length
        seg = np.linalg.norm(np.diff(curve, axis=0), axis=1)
        s_curve = np.concatenate([[0.0], np.cumsum(seg)])         # (S,) scaled arc length
        # nearest curve sample per nucleus (foot-point projection)
        d2 = ((sp[:, None, :] - curve[None, :, :]) ** 2).sum(axis=2)   # (N,S)
        nearest = d2.argmin(axis=1)
        resid2 = d2[np.arange(len(sp)), nearest]
        arc = s_curve[nearest]
        arc_um = arc * scale
        # R²-like confidence: residual scatter off the curve vs total scatter
        total = float((sp ** 2).sum(axis=1).mean())
        conf = float(np.clip(1.0 - resid2.mean() / total, 0.0, 1.0)) if total > 0 else 0.0

    # orient distal(0) -> proximal(1): distal tip = small mitotic nuclei
    if "volume_um3" in df and np.ptp(arc_um) > 0:
        lo_end = arc_um < np.percentile(arc_um, 20)
        hi_end = arc_um > np.percentile(arc_um, 80)
        vol = df["volume_um3"].to_numpy()
        if lo_end.any() and hi_end.any() and vol[lo_end].mean() > vol[hi_end].mean():
            arc_um = arc_um.max() - arc_um
            flags.append("axis:orientation_flipped_by_nucleus_size_heuristic")

    a0, a1 = float(arc_um.min()), float(arc_um.max())
    df["axis_position_um"] = arc_um - a0
    df["axis_position_norm"] = (arc_um - a0) / (a1 - a0) if a1 > a0 else float("nan")

    if conf < confidence_min:
        flags.append(f"axis:low_confidence_{conf:.2f}_needs_manual_axis")
    log.debug("axis: method=%s conf=%.2f length=%.1fµm", method, conf, a1 - a0)
    return df, float(conf), flags


def _mst_backbone_order(sp: np.ndarray) -> np.ndarray | None:
    """Order points along the diameter path of their minimum spanning tree.

    Returns ordered indices into ``sp`` for the nuclei lying on the backbone path, or None if
    the graph machinery is unavailable / degenerate. The path runs between the two graph-most-
    distant nuclei — for a U-shaped gonad this traces the arm rather than cutting the corner.
    """
    try:
        from scipy.sparse.csgraph import minimum_spanning_tree, shortest_path
        from scipy.spatial.distance import squareform, pdist
    except Exception:  # pragma: no cover - scipy always present in this project
        return None
    n = len(sp)
    if n < 3:
        return np.arange(n)
    d = squareform(pdist(sp))
    mst = minimum_spanning_tree(d)
    tree = mst + mst.T                                   # undirected
    # graph diameter: farthest node from 0 -> A, farthest from A -> B, path A..B
    dist0, _ = shortest_path(tree, method="D", indices=0, return_predecessors=True)
    a = int(np.argmax(np.where(np.isfinite(dist0), dist0, -1)))
    distA, predA = shortest_path(tree, method="D", indices=a, return_predecessors=True)
    b = int(np.argmax(np.where(np.isfinite(distA), distA, -1)))
    path = []
    cur = b
    guard = 0
    while cur != -9999 and cur != a and guard <= n:
        path.append(cur)
        cur = int(predA[cur])
        guard += 1
    path.append(a)
    path.reverse()
    return np.asarray(path) if len(path) >= 2 else None


def _fit_principal_curve(sp: np.ndarray, n_samples: int = 400) -> tuple[np.ndarray | None, str]:
    """Fit a smooth principal curve through the MST backbone; return dense samples (S,3).

    Tries a smoothing B-spline; falls back to the (smoothed) backbone polyline. Returns
    (samples, method) or (None, 'none') so the caller can fall back to straight PCA.
    """
    order = _mst_backbone_order(sp)
    if order is None or len(order) < 2:
        return None, "none"
    bb = sp[order]
    # collapse consecutive duplicates (splprep needs strictly increasing parameter)
    keep = np.concatenate([[True], (np.linalg.norm(np.diff(bb, axis=0), axis=1) > 1e-9)])
    bb = bb[keep]
    if len(bb) < 2:
        return None, "none"

    try:
        from scipy.interpolate import splev, splprep

        k = min(3, len(bb) - 1)
        u = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(bb, axis=0), axis=1))])
        u = u / u[-1]
        # s≈m smoothing is sane because coordinates are pre-scaled to unit std
        tck, _ = splprep(bb.T, u=u, k=k, s=len(bb))
        uu = np.linspace(0, 1, n_samples)
        samples = np.vstack(splev(uu, tck)).T
        return samples, "spline"
    except Exception as e:  # noqa: BLE001 - spline can fail on pathological backbones
        log.debug("spline fit failed (%s); using smoothed polyline backbone", e)

    # fallback: moving-average-smoothed backbone, densely re-sampled by arc length
    if len(bb) >= 5:
        w = 3
        kernel = np.ones(w) / w
        bb = np.vstack([np.convolve(bb[:, j], kernel, mode="same") for j in range(3)]).T
    seg = np.linalg.norm(np.diff(bb, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] <= 0:
        return None, "none"
    uu = np.linspace(0, s[-1], n_samples)
    samples = np.vstack([np.interp(uu, s, bb[:, j]) for j in range(3)]).T
    return samples, "polyline"


def _pca_fallback(sp: np.ndarray) -> tuple[np.ndarray, float, str | None]:
    """Straight principal-axis projection (the v0 behaviour) when no curve could be fit."""
    u, s, vt = np.linalg.svd(sp, full_matrices=False)
    proj = sp @ vt[0]
    var_ratio = float((s[0] ** 2) / (s ** 2).sum())
    flag = None if var_ratio >= 0.6 else "axis:pca_fallback_straight_line"
    return proj, var_ratio, flag
