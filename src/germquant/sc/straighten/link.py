"""One global end matching per nucleus: join strands end to end across gaps and crossings the tracer broke.

Candidate ends: types ``free``, ``contact`` and ``loop``; the outward direction of an end is taken from its last
``dir_window_um`` (0.15 to 0.9 um). A pair is admissible when the gap is at most ``max_gap_um``, both outward
directions lie within ``max_angle_deg`` of the link, and the closest approach of the two end tangent lines is at
most ``max_skew_um`` (the skew criterion; ``null`` switches it off). Cost ``gap_cost * gap + angle_a + angle_b +
skew_cost * skew`` (degrees, um); an unlinked end costs ``unlinked_cost``; ``networkx.max_weight_matching`` with
weight ``2 * unlinked_cost - cost`` on a graph built in sorted end order. A strand's two ends never pair with
each other; chains become strands and a closed chain is cut at its most expensive link (``cycle_cut`` ends).
"""
from __future__ import annotations

import math

import numpy as np

from .follow import LINKABLE_ENDS, RawStrand, arclen


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
    ends = []
    for i, st in enumerate(strands):
        for side in (0, 1):
            if st.end_types[side] in LINKABLE_ENDS:
                pos, d, tip = end_direction(st.points, side, win, tip_um)
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
            if span < 0.05:
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
    out: list[RawStrand] = []

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
        out.append(_assemble(strands, seq, recs, None, tip_um))
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
            out.append(_assemble(strands, seq2, recs2, "cycle_cut", tip_um))
        else:
            out.append(_assemble(strands, seq, recs, None, tip_um))
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


def _assemble(strands, seq, recs, cut_type, tip_um: float = 0.0) -> RawStrand:
    parts, strg, links = [], [], []
    members = 0
    n = len(seq)
    for k, (j, side) in enumerate(seq):
        st = strands[j]
        P = st.points if side == 0 else st.points[::-1]
        S = None if st.strength is None else (st.strength if side == 0 else st.strength[::-1])
        P, S = _trim_ends(P, S, front=k > 0, back=k < n - 1, tip_um=tip_um)
        parts.append(P)
        if S is not None:
            strg.append(S)
        links.extend(st.links)
        members += st.members
    links.extend({k: v for k, v in r.items() if k not in ("a", "b")} for r in recs)
    pts = np.vstack(parts)
    first_j, first_side = seq[0]
    last_j, last_side = seq[-1]
    e0 = strands[first_j].end_types[first_side]
    e1 = strands[last_j].end_types[1 - last_side]
    if cut_type is not None:
        e0 = e1 = cut_type
    base = strands[first_j]
    return RawStrand(points=pts, end_types=[e0, e1],
                     strength=np.concatenate(strg) if len(strg) == len(seq) else None,
                     links=links, members=members, nucleus_id=base.nucleus_id, source=base.source,
                     verified=base.verified, import_id=base.import_id)
