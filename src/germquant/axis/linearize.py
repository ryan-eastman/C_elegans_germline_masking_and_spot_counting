"""Germline axis linearization: assign each nucleus a distal->proximal position.

Method (Libuda-lab Gonad Linearization, reimplemented): fit the germline's central axis,
project each nucleus centroid onto it, normalize 0 (distal tip) -> 1 (proximal end). The
published method uses a hand-drawn axis; v1 here fits the axis automatically by principal
curve (PCA seed). When the fit is low-confidence it raises a QC flag so that gonad can be
routed to a 30-second manual axis-draw (semi-automated v1 — see ARCHITECTURE.md Risk 3).

Orientation (which end is distal) is a heuristic: the distal tip is densely packed with
small mitotic nuclei, so we orient so the small-nucleus end = position 0. FLAG: confirm
orientation against a known landmark before trusting absolute direction.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def linearize_germline(
    nuclei: pd.DataFrame,
    *,
    confidence_min: float = 0.60,
) -> tuple[pd.DataFrame, float, list[str]]:
    """Add axis_position_norm (0..1) + axis_position_um. Return (df, confidence, flags)."""
    flags: list[str] = []
    df = nuclei.copy()
    if len(df) < 5:
        df["axis_position_norm"] = float("nan")
        df["axis_position_um"] = float("nan")
        return df, 0.0, ["axis:too_few_nuclei"]

    pts = df[["centroid_z_um", "centroid_y_um", "centroid_x_um"]].to_numpy(dtype=float)
    center = pts.mean(axis=0)
    centered = pts - center

    # PCA: principal axis = direction of greatest variance (germline is elongated)
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    axis_dir = vt[0]
    var_ratio = float((s[0] ** 2) / (s**2).sum())  # how 1D the gonad is = confidence

    proj = centered @ axis_dir  # signed distance along axis (µm)

    # orient distal(0) -> proximal(1): distal tip = small mitotic nuclei.
    if "volume_um3" in df:
        lo_end = proj < np.percentile(proj, 20)
        hi_end = proj > np.percentile(proj, 80)
        if df["volume_um3"].to_numpy()[lo_end].mean() > df["volume_um3"].to_numpy()[hi_end].mean():
            proj = -proj
            flags.append("axis:orientation_flipped_by_nucleus_size_heuristic")

    p0, p1 = proj.min(), proj.max()
    df["axis_position_um"] = proj - p0
    df["axis_position_norm"] = (proj - p0) / (p1 - p0) if p1 > p0 else float("nan")

    if var_ratio < confidence_min:
        flags.append(f"axis:low_confidence_{var_ratio:.2f}_needs_manual_axis")
    return df, var_ratio, flags
