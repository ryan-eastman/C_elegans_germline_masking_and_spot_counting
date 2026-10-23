"""Length accounting: every stretch of SC is counted once.

Within a nucleus, a stretch where two strands run within ``dedupe.lateral_um`` laterally (xy) and
``dedupe.axial_um`` axially (z) of each other, parallel (|cos| > ``parallel_cos``), for at least ``min_run_um``,
stays on the stronger strand (integrated ridge strength) and is cut from the other; remainders shorter than
``min_remainder_um`` are dropped and the cut ends become ``contact`` ends. Crossings (not parallel) stay on
both strands. Across nuclei, the same rule runs between neighbouring nuclei with the priority given to the
strand whose own label holds more of it, so a ribbon traced twice (once from each side of a shared shell) keeps
one copy. A nucleus total is the sum of its strands after this step.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .follow import RawStrand, arclen


def resample(points: np.ndarray, ds: float, values: np.ndarray | None = None):
    s = arclen(points)
    if s[-1] <= 0:
        return points[:1].copy(), (None if values is None else values[:1].copy()), np.zeros(1)
    n = max(2, int(np.ceil(s[-1] / ds)) + 1)
    q = np.linspace(0.0, s[-1], n)
    P = np.stack([np.interp(q, s, points[:, k]) for k in range(3)], 1)
    V = None if values is None else np.interp(q, s, values)
    return P, V, q


def tangents_of(P: np.ndarray) -> np.ndarray:
    if len(P) < 2:
        return np.tile([0.0, 0.0, 1.0], (len(P), 1))
    T = np.gradient(P, axis=0)
    n = np.linalg.norm(T, axis=1, keepdims=True)
    return T / np.maximum(n, 1e-12)


def strength_sum(st: RawStrand) -> float:
    if st.strength is None or len(st.strength) == 0:
        return st.length_um
    return float(np.mean(st.strength) * st.length_um)


def _scaled(P, lat, ax):
    return np.stack([P[:, 0] / ax, P[:, 1] / lat, P[:, 2] / lat], 1)


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Inclusive index runs of True."""
    out = []
    i, n = 0, len(mask)
    while i < n:
        if mask[i]:
            j = i
            while j + 1 < n and mask[j + 1]:
                j += 1
            out.append((i, j))
            i = j + 1
        else:
            i += 1
    return out


def cut_duplicates(st: RawStrand, kept_P: np.ndarray, kept_T: np.ndarray, p, ds: float = 0.05) -> list[RawStrand]:
    """Remove from ``st`` the parallel stretches that duplicate the kept samples; return the pieces."""
    lat = float(p["dedupe.lateral_um"])
    ax = float(p["dedupe.axial_um"])
    pc = float(p["dedupe.parallel_cos"])
    min_run = float(p["dedupe.min_run_um"])
    min_rem = float(p["dedupe.min_remainder_um"])
    if len(kept_P) == 0 or len(st.points) < 2:
        return [st]
    P, V, q = resample(st.points, ds, st.strength)
    T = tangents_of(P)
    tree = cKDTree(_scaled(kept_P, lat, ax))
    d, k = tree.query(_scaled(P, lat, ax), distance_upper_bound=1.0)
    near = np.isfinite(d)
    par = np.zeros(len(P), bool)
    if near.any():
        kk = k[near]
        par[near] = np.abs((T[near] * kept_T[kk]).sum(1)) > pc
    dup = near & par
    cut = np.zeros(len(P), bool)
    for a, b in _runs(dup):
        if q[b] - q[a] >= min_run - 1e-9:
            cut[a:b + 1] = True
    if not cut.any():
        return [st]
    pieces = []
    keep_runs = _runs(~cut)
    for a, b in keep_runs:
        if q[b] - q[a] < min_rem:
            continue
        e0 = st.end_types[0] if a == 0 else "contact"
        e1 = st.end_types[1] if b == len(P) - 1 else "contact"
        pieces.append(RawStrand(points=P[a:b + 1].copy(), end_types=[e0, e1],
                                strength=None if V is None else V[a:b + 1].copy(), links=list(st.links),
                                members=st.members, nucleus_id=st.nucleus_id, source=st.source,
                                verified=st.verified, import_id=st.import_id))
    return pieces


def dedupe_within(strands: list[RawStrand], p) -> list[RawStrand]:
    """Greedy over strands in descending integrated strength (ties: longer, then input order)."""
    if len(strands) < 2:
        return list(strands)
    order = sorted(range(len(strands)), key=lambda i: (-strength_sum(strands[i]), -strands[i].length_um, i))
    kept: list[RawStrand] = []
    kP, kT = [], []
    for i in order:
        st = strands[i]
        if st.source != "auto":
            kept.append(st)
            continue
        allP = np.vstack(kP) if kP else np.zeros((0, 3))
        allT = np.vstack(kT) if kT else np.zeros((0, 3))
        for piece in cut_duplicates(st, allP, allT, p):
            kept.append(piece)
            P, _, _ = resample(piece.points, 0.05)
            kP.append(P)
            kT.append(tangents_of(P))
    return kept


def drop_axial_ghosts(strands: list[RawStrand], p) -> tuple[list[RawStrand], int]:
    """Remove axial shadows: a short (``dedupe.ghost_max_um``), steep (mean |t_z| at least ``dedupe.ghost_tz_min``)
    traced strand with an end inside the axial PSF reach (``dedupe.ghost_reach_um``, axial and lateral) of a longer
    strand is the elongated PSF beyond that strand's end, or at a crossing, not a ribbon."""
    max_len = float(p["dedupe.ghost_max_um"])
    tz_min = float(p["dedupe.ghost_tz_min"])
    ax, lat = (float(v) for v in p["dedupe.ghost_reach_um"])
    if len(strands) < 2 or max_len <= 0:
        return list(strands), 0
    keep = []
    dropped = 0
    for i, st in enumerate(strands):
        L = st.length_um
        if st.source != "auto" or L > max_len or len(st.points) < 3:
            keep.append(st)
            continue
        T = tangents_of(st.points)
        if float(np.mean(np.abs(T[:, 0]))) < tz_min:
            keep.append(st)
            continue
        ghost = False
        mine = np.stack([st.points[0], st.points[-1]])
        for j, ot in enumerate(strands):
            if j == i or ot.length_um <= L:
                continue
            d = mine[:, None, :] - ot.points[None, :, :]
            if np.any((d[..., 0] / ax) ** 2 + (d[..., 1] ** 2 + d[..., 2] ** 2) / lat ** 2 <= 1.0):
                ghost = True
                break
        if ghost:
            dropped += 1
        else:
            keep.append(st)
    return keep, dropped


def dedupe_between(group_a: list[RawStrand], group_b: list[RawStrand], frac_a: list[float], frac_b: list[float],
                   p) -> tuple[list[RawStrand], list[RawStrand]]:
    """Two neighbouring nuclei: priority (fraction of the strand inside its own label, then strength); a strand
    is compared only with the kept strands of the OTHER nucleus. Returns the new groups."""
    items = [("a", i, st, frac_a[i]) for i, st in enumerate(group_a)] + \
            [("b", i, st, frac_b[i]) for i, st in enumerate(group_b)]
    items.sort(key=lambda t: (-t[3], -strength_sum(t[2]), -t[2].length_um, t[0], t[1]))
    kept = {"a": [], "b": []}
    samples = {"a": ([], []), "b": ([], [])}
    for g, _i, st, _f in items:
        other = "b" if g == "a" else "a"
        oP, oT = samples[other]
        allP = np.vstack(oP) if oP else np.zeros((0, 3))
        allT = np.vstack(oT) if oT else np.zeros((0, 3))
        pieces = [st] if st.source != "auto" else cut_duplicates(st, allP, allT, p)
        for piece in pieces:
            kept[g].append(piece)
            P, _, _ = resample(piece.points, 0.05)
            samples[g][0].append(P)
            samples[g][1].append(tangents_of(P))
    # keep the input order of each group as far as possible (pieces in priority order)
    return kept["a"], kept["b"]
