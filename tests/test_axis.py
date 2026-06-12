"""Axis linearization on a curved (hairpin) gonad — the case a straight PCA axis gets wrong.

A real C. elegans gonad arm is U-shaped: the two arms run close together in space but are far
apart *along the germline*. The principal-curve linearizer must (a) measure arc length, not the
chord, and (b) keep the two arms separated in position even though they nearly touch.
"""
import numpy as np
import pandas as pd

from germquant.axis import linearize_germline


def _hairpin(n=120, L=40.0, gap=6.0, jitter=0.6, seed=0):
    """Centroids along a hairpin centerline: up arm1, semicircular cap, down arm2.

    Returns (df, t) where t in [0,1] is the true arc parameter for each nucleus (ground truth
    ordering). All coordinates in microns.
    """
    rng = np.random.default_rng(seed)
    r = gap / 2.0
    # piecewise lengths: arm1 (L) -> cap (pi*r) -> arm2 (L)
    cap = np.pi * r
    total = 2 * L + cap
    s = np.linspace(0, total, n)
    xs, ys = [], []
    for si in s:
        if si <= L:                                   # arm 1: x=-r, y up
            x, y = -r, si
        elif si <= L + cap:                            # semicircular cap over the top
            a = (si - L) / r                           # angle 0..pi
            x, y = -r * np.cos(a), L + r * np.sin(a)
        else:                                          # arm 2: x=+r, y down
            x, y = r, L - (si - L - cap)
        xs.append(x)
        ys.append(y)
    x = np.array(xs) + rng.normal(0, jitter, n)
    y = np.array(ys) + rng.normal(0, jitter, n)
    z = rng.normal(0, jitter, n)                        # thin in z
    df = pd.DataFrame({
        "nucleus_id": np.arange(1, n + 1),
        "centroid_z_um": z, "centroid_y_um": y, "centroid_x_um": x,
        "volume_um3": np.full(n, 5.0),
    })
    return df, s / total


def test_arc_length_not_chord():
    df, t = _hairpin()
    out, conf, flags = linearize_germline(df)
    length = float(out["axis_position_um"].max())
    # the two endpoints are only ~gap apart in space (chord ~6µm); the true path is ~2*L
    assert length > 60, f"germline length {length:.1f}µm collapsed toward the chord (straight-PCA bug)"
    assert conf > 0.6, f"curved gonad spuriously low-confidence: {conf:.2f}"


def test_ordering_follows_the_curve():
    df, t = _hairpin()
    out, conf, flags = linearize_germline(df)
    pos = out["axis_position_norm"].to_numpy()
    # position must be monotonic in the true arc parameter (up to global flip)
    rho = np.corrcoef(np.argsort(np.argsort(pos)), np.argsort(np.argsort(t)))[0, 1]
    assert abs(rho) > 0.95, f"axis ordering does not follow the germline curve (|rho|={abs(rho):.2f})"


def test_opposite_arms_stay_separated():
    df, t = _hairpin(jitter=0.3)
    out, _, _ = linearize_germline(df)
    pos = out["axis_position_norm"].to_numpy()
    x = df["centroid_x_um"].to_numpy()
    y = df["centroid_y_um"].to_numpy()
    # pick a nucleus low on each arm (same y, opposite x) — Euclidean-close, should be arc-far
    lowL = np.argmin(np.where(x < 0, y, np.inf))       # arm1, smallest y
    lowR = np.argmin(np.where(x > 0, y, np.inf))       # arm2, smallest y
    assert abs(pos[lowL] - pos[lowR]) > 0.5, (
        "opposite arms collapsed to the same position — the U-shape was cut by a straight axis"
    )
