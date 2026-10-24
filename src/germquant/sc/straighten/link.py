"""One global end matching per nucleus: join strands end to end across gaps and crossings the tracer broke.

Candidate ends: types ``free``, ``contact`` and ``loop``; the outward direction of an end is taken from its last
``dir_window_um`` (0.15 to 0.9 um). A pair is admissible when the gap is at most ``max_gap_um``, both outward
directions lie within ``max_angle_deg`` of the link, and the closest approach of the two end tangent lines is at
most ``max_skew_um`` (the skew criterion; ``null`` switches it off). Cost ``gap_cost * gap + angle_a + angle_b +
skew_cost * skew`` (degrees, um); an unlinked end costs ``unlinked_cost``; ``networkx.max_weight_matching`` with
weight ``2 * unlinked_cost - cost`` on a graph built in sorted end order. A strand's two ends never pair with
each other; chains become strands and a closed chain is cut at its most expensive link (``cycle_cut`` ends).

v2 (protocol changelog entry 21): with ``link.junction_dirs`` an end whose last ``JUNCTION_REACH_UM`` holds points
within one PSF FWHM (axial and lateral) of another strand takes its anchor and outward direction from the first
point beyond them (inside a junction the traced end hooks toward the neighbour), the gap still counted from the
tip; with ``link.max_chain_um`` a matched chain longer than the bound is cut at its most expensive link, again
until none is longer, and the cut ends are free.
"""
from __future__ import annotations

import math

import numpy as np

from .follow import LINKABLE_ENDS, RawStrand, arclen

JUNCTION_REACH_UM = 1.0          # entry 21 (a): how far from the tip junction points are looked for
FWHM = 2.0 * math.sqrt(2.0 * math.log(2.0))


def end_direction(points: np.ndarray, side: int, window=(0.3, 1.2), tip_um: float = 0.3):
    """(anchor position, outward unit direction, anchor offset) of end ``side`` (0 = first point, 1 = last).

    The traced tip of a PSF-blurred end drifts along the optical axis, so the end is represented by an anchor
    ``tip_um`` inside the tip (at most a quarter of the strand) and a direction from the points between
    ``window`` (um from the tip) to the anchor."""
    P = points if side == 0 else points[::-1]
    s = arclen(P)
    tip = min(float(tip_um), s[-1] / 4.0)
    k = int(np.searchsorted(s, tip))
    k = min(max(k, 0), len(P) - 1)
    anchor = P[k]
    lo, hi = max(window[0], tip), max(window[1], tip + 0.3)
    sel = (s >= lo) & (s <= hi)
    if sel.sum() < 2:
        sel = s >= tip
    if sel.sum() < 2:
        sel = np.ones(len(P), bool)
    d = anchor - P[sel].mean(0)
    n = np.linalg.norm(d)
    if n < 1e-9 and len(P) > 1:
        d = P[0] - P[1]
        n = np.linalg.norm(d)
    return anchor.copy(), (d / n if n > 1e-12 else np.array([0.0, 0.0, 1.0])), float(s[k])


def junction_shift(points: np.ndarray, side: int, others: np.ndarray | None, fwhm_um) -> float:
    """Arc length from end ``side`` to the first point beyond the points of the strand's last
    ``JUNCTION_REACH_UM`` that lie within one PSF FWHM (``fwhm_um`` = (axial, lateral)) of ``others`` (the other
    strands' points); 0 when there are none."""
    if others is None or not len(others) or len(points) < 3:
        return 0.0
    from scipy.spatial import cKDTree

    P = points if side == 0 else points[::-1]
    s = arclen(P)
    near_end = s <= JUNCTION_REACH_UM
    if not near_end.any():
        return 0.0
    fz, fxy = (float(v) for v in fwhm_um)
    sc = np.array([1.0 / fz, 1.0 / fxy, 1.0 / fxy])
    tree = cKDTree(np.asarray(others, float) * sc)
    d, _ = tree.query(P[near_end] * sc, distance_upper_bound=1.0)
    hit = np.nonzero(np.isfinite(d))[0]
    if not len(hit):
        return 0.0
    k = int(hit.max()) + 1                       # the first point beyond the junction points
    if k >= len(P) or s[-1] - s[k] < 0.3:        # nothing usable beyond them: keep the v1 end
        return 0.0
    return float(s[k])


def end_geometry(points: np.ndarray, side: int, window, tip_um: float, shift_um: float = 0.0):
    """end_direction measured from ``shift_um`` inside the end: (anchor, direction, anchor arc from the tip)."""
    if shift_um <= 0:
        return end_direction(points, side, window, tip_um)
    P = points if side == 0 else points[::-1]
    s = arclen(P)
    k = int(np.searchsorted(s, shift_um))
    Q = P[k:] if side == 0 else P[k:]
    pos, d, tip = end_direction(Q, 0, window, tip_um)
    return pos, d, float(s[k]) + tip


def _angle(a, b) -> float:
    return math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(a, b))))))


def skew_distance(pa, da, pb, db) -> float:
    """Closest approach of the lines pa + l da and pb + m db (um)."""
    w = pb - pa
    c = np.cross(da, db)
    nc = np.linalg.norm(c)
    if nc < 1e-6:
        return float(np.linalg.norm(np.cross(w, da)))
    return float(abs(np.dot(w, c)) / nc)


def link_candidates(strands: list[RawStrand], p) -> tuple[list, list]:
    win = tuple(float(v) for v in p["link.dir_window_um"])
    tip_um = float(p["link.tip_um"])
    max_gap = float(p["link.max_gap_um"])
    max_gap_cc = float(p["link.max_gap_contact_um"])
    max_ang = float(p["link.max_angle_deg"])
    max_skew = p["link.max_skew_um"]
    gc = float(p["link.gap_cost_per_um"])
    sc = float(p["link.skew_cost_per_um"])
    unl = float(p["link.unlinked_cost"])
    junction = bool(p.get("link.junction_dirs", False))
    fwhm = tuple(FWHM * float(v) for v in p["psf.fallback_sigma_um"])
    ends = []
    for i, st in enumerate(strands):
        others = None
        if junction:
            rest = [o.points for j, o in enumerate(strands) if j != i and len(o.points)]
            others = np.vstack(rest) if rest else None
        for side in (0, 1):
            if st.end_types[side] in LINKABLE_ENDS:
                shift = junction_shift(st.points, side, others, fwhm) if junction else 0.0
                pos, d, tip = end_geometry(st.points, side, win, tip_um, shift)
                ends.append((i, side, pos, d, tip, st.end_types[side]))
    cands = []
    for a in range(len(ends)):
        ia, sa, pa, da, ta, ea = ends[a]
        for b in range(a + 1, len(ends)):
            ib, sb, pb, db, tb, eb = ends[b]
            if ia == ib:
                continue
            w = pb - pa
            span = float(np.linalg.norm(w))
            gap = max(span - (ta + tb), 0.0)          # the tip-to-tip gap estimated from the anchors
            limit = max_gap_cc if (ea == "contact" and eb == "contact") else max_gap
            if gap > limit:
                continue
            if span < 0.05 or span <= ta + tb:
                # tips that meet or overlap (PSF-blurred tips drift past each other): judge the collinearity of
                # the two ends instead of the direction of a vanishing link
                ang = _angle(da, -db)
                aa = ab = 0.5 * ang
            else:
                u = w / span
                aa = _angle(da, u)
                ab = _angle(db, -u)
            if aa > max_ang or ab > max_ang:
                continue
            skew = skew_distance(pa, da, pb, db)
            if max_skew is not None:
                if skew > float(max_skew):
                    continue
                cost = gc * gap + aa + ab + sc * skew
            else:
                cost = gc * gap + aa + ab
            if cost >= 2 * unl:
                continue
            cands.append({"a": a, "b": b, "gap_um": gap, "angle_a_deg": aa, "angle_b_deg": ab,
                          "skew_um": skew, "cost": cost, "mid_um": [float(v) for v in 0.5 * (pa + pb)]})
    return ends, cands


def link_strands(strands: list[RawStrand], p) -> list[RawStrand]:
    """The matched chains as new strands (unlinked strands pass through unchanged)."""
    if len(strands) < 2:
        return list(strands)
    import networkx as nx

    ends, cands = link_candidates(strands, p)
    if not cands:
        return list(strands)
    unl = float(p["link.unlinked_cost"])
    G = nx.Graph()
    G.add_nodes_from(range(len(ends)))
    for c in cands:
        G.add_edge(c["a"], c["b"], weight=2 * unl - c["cost"], rec=c)
    match = nx.max_weight_matching(G, maxcardinality=False)
    tip_um = float(p["link.tip_um"])
    partner: dict[tuple, tuple] = {}
    for a, b in sorted(tuple(sorted(e)) for e in match):
        rec = G[a][b]["rec"]
        ea, eb = (ends[a][0], ends[a][1]), (ends[b][0], ends[b][1])
        partner[ea] = (eb, rec)
        partner[eb] = (ea, rec)
    # components of the strand graph: paths or cycles
    seen: set[int] = set()
    paths: list[tuple[list, list, str | None]] = []

    def walk(start: int, enter_side: int):
        seq, recs = [], []
        cur, side = start, enter_side
        while True:
            seq.append((cur, side))
            seen.add(cur)
            exit_end = (cur, 1 - side)
            nxt = partner.get(exit_end)
            if nxt is None:
                return seq, recs, False
            (j, sj), rec = nxt
            if j in seen:
                return seq, recs + [rec], True
            recs.append(rec)
            cur, side = j, sj

    order = list(range(len(strands)))
    for i in order:
        if i in seen:
            continue
        heads = [s for s in (0, 1) if (i, s) not in partner]
        if not heads:
            continue
        seq, recs, _ = walk(i, heads[0])
        paths.append((seq, recs, None))
    for i in order:
        if i in seen:
            continue
        seq, recs, closed = walk(i, 0)
        # a closed chain: cut it at its most expensive link
        if closed and recs:
            k = int(np.argmax([r["cost"] for r in recs]))
            # the chain links are recs[0] .. recs[-1] where recs[-1] closes the loop back to seq[0];
            # cutting link k means starting the path at seq[k + 1]
            n = len(seq)
            start = (k + 1) % n
            seq2 = seq[start:] + seq[:start]
            recs2 = [recs[(start + j) % n] for j in range(n - 1)]
            paths.append((seq2, recs2, "cycle_cut"))
        else:
            paths.append((seq, recs, None))
    max_chain = p.get("link.max_chain_um", None)
    if max_chain is None:
        return [_assemble(strands, seq, recs, cut, tip_um) for seq, recs, cut in paths]
    # the chain bound (entry 21 (b)): a chain longer than max_chain_um is cut at its most expensive link, again
    # until none is longer or no link is left; the cut ends are free
    out: list[RawStrand] = []
    todo = [(seq, recs, (cut, cut) if cut else (None, None)) for seq, recs, cut in paths]
    while todo:
        seq, recs, (c0, c1) = todo.pop(0)
        st = _assemble(strands, seq, recs, None, tip_um, ends_override=(c0, c1))
        if not recs or st.length_um <= float(max_chain):
            out.append(st)
            continue
        k = int(np.argmax([r["cost"] for r in recs]))
        todo.insert(0, (seq[k + 1:], recs[k + 1:], ("free", c1)))
        todo.insert(0, (seq[:k + 1], recs[:k], (c0, "free")))
    return out


def _trim_ends(P: np.ndarray, S, front: bool, back: bool, tip_um: float):
    """Drop the drifting tip (``tip_um``, at most a quarter of the piece) at linked ends."""
    if not (front or back) or len(P) < 4:
        return P, S
    s = arclen(P)
    tip = min(tip_um, s[-1] / 4.0)
    keep = np.ones(len(P), bool)
    if front:
        keep &= s >= tip
    if back:
        keep &= s <= s[-1] - tip
    if keep.sum() < 2:
        return P, S
    return P[keep], (None if S is None else S[keep])


def _past_previous(prev: np.ndarray, P: np.ndarray, n_dir: int = 6) -> int:
    """Index of the first point of ``P`` that lies beyond the last point of ``prev`` along prev's end direction
    (taken over its last ``n_dir`` points); len(P) when none does. Tips that met or overlapped are joined without
    a backward step, so the assembled polyline (and its spline) never folds back."""
    d = prev[-1] - prev[max(len(prev) - n_dir, 0)]
    nd = float(np.linalg.norm(d))
    if nd < 1e-9:
        return 0
    ahead = (P - prev[-1]) @ (d / nd) > 0
    return int(np.argmax(ahead)) if ahead.any() else len(P)


def _assemble(strands, seq, recs, cut_type, tip_um: float = 0.0, ends_override=(None, None)) -> RawStrand:
    parts, strg, links = [], [], []
    members = 0
    n = len(seq)
    all_strength = True
    for k, (j, side) in enumerate(seq):
        st = strands[j]
        P = st.points if side == 0 else st.points[::-1]
        S = None if st.strength is None else (st.strength if side == 0 else st.strength[::-1])
        P, S = _trim_ends(P, S, front=k > 0, back=k < n - 1, tip_um=tip_um)
        links.extend(st.links)
        members += st.members
        if parts:
            i0 = _past_previous(parts[-1], P)
            P = P[i0:]
            S = None if S is None else S[i0:]
            if len(P) < 2:
                continue                  # the piece lies entirely within the previous one's tip
        parts.append(P)
        if S is None:
            all_strength = False
        else:
            strg.append(S)
    links.extend({k: v for k, v in r.items() if k not in ("a", "b")} for r in recs)
    pts = np.vstack(parts)
    first_j, first_side = seq[0]
    last_j, last_side = seq[-1]
    e0 = strands[first_j].end_types[first_side]
    e1 = strands[last_j].end_types[1 - last_side]
    if cut_type is not None:
        e0 = e1 = cut_type
    if ends_override[0] is not None:
        e0 = ends_override[0]
    if ends_override[1] is not None:
        e1 = ends_override[1]
    base = strands[first_j]
    return RawStrand(points=pts, end_types=[e0, e1],
                     strength=np.concatenate(strg) if all_strength and strg else None,
                     links=links, members=members, nucleus_id=base.nucleus_id, source=base.source,
                     verified=base.verified, import_id=base.import_id)
