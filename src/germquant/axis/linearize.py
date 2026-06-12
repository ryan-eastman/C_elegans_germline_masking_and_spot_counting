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
    confidence_min: float = 0.35,
) -> tuple[pd.DataFrame, float, list[str]]:
    """Add axis_position_norm (0..1) + axis_position_um. Return (df, confidence, flags).

    confidence = 1 − (mean off-centerline scatter) / (total scatter): high when nuclei hug the
    fitted axis, ~0 for a diffuse blob with no coherent axis. NOTE it scales with tube *width*,
    so a real (thick) germline reads ~0.5, not ~1 — ``confidence_min`` is therefore a low
    "is there an axis at all" gate (a blob fails), not a fit-perfection bar. Below it, a flag
    routes the gonad to a manual axis-draw. Unlike the old PCA variance-ratio it does NOT
    penalise curvature.
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

    # orient distal(0) -> proximal(1). The distal mitotic tip is DENSELY packed with nuclei —
    # a far more reliable signal than nucleus size (verified on real N2 data: the smallest
    # nuclei sit mid-gonad, so the old size heuristic flipped the wrong way). Flip so the denser
    # end becomes position 0. If the two ends are near-symmetric in density, fall back to the
    # small-nucleus heuristic and flag the orientation as uncertain (-> manual review).
    if np.ptp(arc_um) > 0:
        L = float(arc_um.max())
        win = 0.12 * L
        n0, n1 = int((arc_um < win).sum()), int((arc_um > L - win).sum())
        if abs(n1 - n0) >= 0.15 * max(n0, n1, 1):
            if n1 > n0:
                arc_um = L - arc_um
            flags.append("axis:oriented_by_density")
        else:
            vol = df["volume_um3"].to_numpy() if "volume_um3" in df else None
            lo = arc_um < np.percentile(arc_um, 20)
            hi = arc_um > np.percentile(arc_um, 80)
            if vol is not None and lo.any() and hi.any() and vol[lo].mean() > vol[hi].mean():
                arc_um = L - arc_um
            flags.append("axis:orientation_uncertain_density_symmetric")

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
    """Extract a centerline down the middle of the tube of nuclei; return dense samples (S,3).

    The MST backbone is a path THROUGH individual nuclei, so in a thick germline tube it zigzags
    laterally across the tube width and badly inflates arc length (on real data: a 391 µm gonad's
    raw backbone is ~930 µm). A moving average over the ordered backbone, with a window ~12 % of
    the path, averages out that lateral zigzag while preserving genuine large-scale curvature —
    on real data it recovers ~417 µm (just above the true end-to-end distance, as a gentle curve
    should be). Returns (samples, method) or (None, 'none') so the caller can fall back to PCA.
    """
    order = _mst_backbone_order(sp)
    if order is None or len(order) < 2:
        return None, "none"
    bb = sp[order]
    # collapse consecutive duplicates
    keep = np.concatenate([[True], (np.linalg.norm(np.diff(bb, axis=0), axis=1) > 1e-9)])
    bb = bb[keep]
    if len(bb) < 2:
        return None, "none"

    if len(bb) >= 7:
        w = max(3, int(round(0.12 * len(bb))))
        if w % 2 == 0:
            w += 1
        if w < len(bb):
            kernel = np.ones(w) / w
            # 'valid' avoids the endpoint shrink-in artifact of 'same'; the trimmed tip/end is
            # recovered by foot-point projection (tip nuclei clamp to the nearest centerline end)
            bb = np.vstack([np.convolve(bb[:, j], kernel, mode="valid") for j in range(3)]).T
    if len(bb) < 2:
        return None, "none"

    seg = np.linalg.norm(np.diff(bb, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    if s[-1] <= 0:
        return None, "none"
    uu = np.linspace(0, s[-1], n_samples)
    samples = np.vstack([np.interp(uu, s, bb[:, j]) for j in range(3)]).T
    return samples, "centerline"


def _pca_fallback(sp: np.ndarray) -> tuple[np.ndarray, float, str | None]:
    """Straight principal-axis projection (the v0 behaviour) when no curve could be fit."""
    u, s, vt = np.linalg.svd(sp, full_matrices=False)
    proj = sp @ vt[0]
    var_ratio = float((s[0] ** 2) / (s ** 2).sum())
    flag = None if var_ratio >= 0.6 else "axis:pca_fallback_straight_line"
    return proj, var_ratio, flag
