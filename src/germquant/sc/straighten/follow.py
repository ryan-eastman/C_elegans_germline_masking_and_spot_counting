"""Seeded ridge following (a predictor-corrector tracer) on the line fields of one nucleus crop.

From the strongest unclaimed crest seed (S >= ``follow.seed_min``, descending S, ties by voxel index), both
directions are followed: predict ``step_um`` along the current direction t; where the voxel is line-like
(lineness <= ``lineness_max``) correct with the Newton offset (normal-plane component, capped) and update
``t = normalise(m t + (1 - m) d_local)``; inside junction cores the momentum carries the trace through.
Stops and end types: weak ridge (S < ``lo``) for more than ``weak_run_um`` (trimmed to the last strong point;
``free``); leaving the territory (``territory_edge``); entering another strand's claim tube in parallel
(``contact``: with ``follow.contact_test: local``, the default and v1, the trace is compared with the local line
direction, which inside an unresolved junction is the trace's own; with ``claim_tangent`` (protocol changelog
entry 20, withdrawn by entry 21) with the tangent of the strand that claimed the voxel, kept with the claim);
re-entering its own path more than ``loop_back_um`` back along the strand (``loop``);
``max_length_um`` (``max_length``); within ``stack_edge_um`` of the first or last plane (``stack_edge``).
A claim tube (``claim_vox``) is stamped around every accepted strand; a claimed voxel never seeds.

An end that runs into another strand's interior (a T contact) is trimmed to where it enters that strand's tube
and becomes a ``contact`` end (``resolve_contacts``). The caller runs it after the axial shadows are removed
(germquant.sc.straighten.trace), so a shadow that is dropped later never leaves a ``contact`` end behind, and
``free_orphan_contacts`` turns a ``contact`` end back into a ``free`` one when no surviving strand lies within
reach of it (its partner was removed by the accounting). Positions are um in the pipeline frame
((crop voxel + org) x spacing), axis order (z, y, x).
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from .fields import Fields

LINKABLE_ENDS = ("free", "contact", "loop")


@dataclass
class RawStrand:
    points: np.ndarray                     # (n, 3) um, end0 -> end1
    end_types: list                        # [end0, end1]
    strength: np.ndarray | None = None     # S at each point (null units)
    links: list = field(default_factory=list)
    members: int = 1
    nucleus_id: int = 0
    source: str = "auto"
    verified: bool = False
    import_id: str = ""

    @property
    def length_um(self) -> float:
        if len(self.points) < 2:
            return 0.0
        return float(np.linalg.norm(np.diff(self.points, axis=0), axis=1).sum())

    def reversed(self) -> "RawStrand":
        return RawStrand(self.points[::-1].copy(), [self.end_types[1], self.end_types[0]],
                         None if self.strength is None else self.strength[::-1].copy(), list(self.links),
                         self.members, self.nucleus_id, self.source, self.verified, self.import_id)


def arclen(points: np.ndarray) -> np.ndarray:
    if len(points) < 2:
        return np.zeros(len(points))
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))])


class BudgetExceeded(Exception):
    pass


def follow_nucleus(F: Fields, terr: np.ndarray, spacing, org, p, *, z_planes: int | None = None,
                   budget_s: float | None = None, t_start: float | None = None) -> tuple[list[RawStrand], dict]:
    """Trace every strand of one crop. ``p`` is the resolved Params; ``z_planes`` the image's plane count
    (for the stack-edge stop). Returns (strands, info) with info: n_seeds, n_traced, budget_exceeded."""
    sp = [float(s) for s in spacing]
    o = [int(v) for v in org]
    S, off, D, LIN = F.S, F.off, F.D, F.LIN
    Z, Y, X = S.shape
    seed_min = float(p["follow.seed_min"])
    lo = float(p["follow.lo"])
    h = float(p["follow.step_um"])
    m = float(p["follow.momentum"])
    max_corr = float(p["follow.max_correction_um"])
    lin_max = float(p["follow.lineness_max"])
    weak_max = float(p["follow.weak_run_um"])
    cr = [int(v) for v in p["follow.claim_vox"]]
    contact_cos = float(p["follow.contact_cos"])
    loop_back = float(p["follow.loop_back_um"])
    min_len = float(p["follow.min_length_um"])
    max_len = float(p["follow.max_length_um"])
    edge_um = float(p["follow.stack_edge_um"])
    cos_turn = math.cos(math.radians(float(p["follow.max_turn_deg"])))
    # the contact stop against the claiming strand's own tangent (protocol changelog entry 20); "local" is v1
    claim_tangent = str(p.get("follow.contact_test", "local")) == "claim_tangent"
    z_lo_um = edge_um
    z_hi_um = ((z_planes - 1) * sp[0] - edge_um) if z_planes else math.inf
    t0_clock = time.perf_counter() if t_start is None else t_start
    budget = float(budget_s) if budget_s else math.inf
    claim = np.zeros(S.shape, np.int32)
    claim_t = np.zeros(S.shape + (3,), np.float32) if claim_tangent else None
    info = {"n_seeds": 0, "n_traced": 0, "budget_exceeded": False}
    steps_done = [0]

    def vox(q):
        return (int(round(q[0] / sp[0] - o[0])), int(round(q[1] / sp[1] - o[1])), int(round(q[2] / sp[2] - o[2])))

    def inside(v):
        return 0 <= v[0] < Z and 0 <= v[1] < Y and 0 <= v[2] < X and bool(terr[v])

    def phys(dv):
        a, b, c = float(dv[0]) * sp[0], float(dv[1]) * sp[1], float(dv[2]) * sp[2]
        n = math.sqrt(a * a + b * b + c * c)
        return (a / n, b / n, c / n) if n > 0 else (0.0, 0.0, 0.0)

    def trace(p0, t0, sid, own, sign, length_cap):
        pts = [p0]
        strg = [float(S[vox(p0)]) if inside(vox(p0)) else 0.0]
        q0, q1, q2 = p0
        t = list(t0)
        weak, last_strong = 0.0, 0
        stop = "free"
        s_here = 0.0
        n_max = int(length_cap / h) + 1
        for _ in range(n_max):
            steps_done[0] += 1
            if steps_done[0] % 200 == 0 and time.perf_counter() - t0_clock > budget:
                info["budget_exceeded"] = True
                stop = "budget"
                break
            a0, a1, a2 = q0 + h * t[0], q1 + h * t[1], q2 + h * t[2]
            v = vox((a0, a1, a2))
            if not inside(v):
                stop = "territory_edge"
                break
            lin = float(LIN[v])
            owner = int(claim[v])
            in_other = bool(owner) and owner != sid + 1
            if in_other and claim_t is not None:
                tc = claim_t[v]
                if abs(float(tc[0]) * t[0] + float(tc[1]) * t[1] + float(tc[2]) * t[2]) > contact_cos:
                    stop = "contact"
                    break
            elif in_other and lin <= lin_max:
                d2 = phys(D[v])
                if abs(d2[0] * t[0] + d2[1] * t[1] + d2[2] * t[2]) > contact_cos:
                    stop = "contact"
                    break
            # inside another strand's tube (a crossing) the field belongs to that strand: pure prediction; the same
            # when the local line turns away by more than max_turn_deg in one step (another ridge crossing ours:
            # an SC bends by at most a few degrees per step)
            dl = phys(D[v]) if lin <= lin_max and not in_other else None
            if dl is not None:
                cdl = dl[0] * t[0] + dl[1] * t[1] + dl[2] * t[2]
                if cdl < 0:
                    dl, cdl = (-dl[0], -dl[1], -dl[2]), -cdl
                if cdl < cos_turn:
                    dl = None
            if dl is not None:
                ov = off[v]
                c0 = (v[0] + max(-0.5, min(0.5, float(ov[0]))) + o[0]) * sp[0] - a0
                c1 = (v[1] + max(-0.5, min(0.5, float(ov[1]))) + o[1]) * sp[1] - a1
                c2 = (v[2] + max(-0.5, min(0.5, float(ov[2]))) + o[2]) * sp[2] - a2
                dot = c0 * t[0] + c1 * t[1] + c2 * t[2]
                c0, c1, c2 = c0 - dot * t[0], c1 - dot * t[1], c2 - dot * t[2]
                cn = math.sqrt(c0 * c0 + c1 * c1 + c2 * c2)
                if cn > max_corr:
                    f = max_corr / cn
                    c0, c1, c2 = c0 * f, c1 * f, c2 * f
                a0, a1, a2 = a0 + c0, a1 + c1, a2 + c2
                n0, n1, n2 = m * t[0] + (1 - m) * dl[0], m * t[1] + (1 - m) * dl[1], m * t[2] + (1 - m) * dl[2]
                nn = math.sqrt(n0 * n0 + n1 * n1 + n2 * n2)
                if nn > 0:
                    t = [n0 / nn, n1 / nn, n2 / nn]
                v = vox((a0, a1, a2))
                if not inside(v):
                    stop = "territory_edge"
                    break
            if a0 < z_lo_um or a0 > z_hi_um:
                stop = "stack_edge"
                break
            u_now = sign * (s_here + h)
            prev = own.get(v)
            if prev is not None and abs(u_now - prev) > loop_back:
                stop = "loop"
                break
            s_here += h
            if v not in own:
                own[v] = u_now
            pts.append((a0, a1, a2))
            sv = float(S[v])
            strg.append(sv)
            if sv >= lo:
                weak = 0.0
                last_strong = len(pts) - 1
            else:
                weak += h
                if weak > weak_max:
                    stop = "free"
                    break
            q0, q1, q2 = a0, a1, a2
        else:
            stop = "max_length"
        if stop in ("territory_edge", "stack_edge") and weak > 0:
            # the ridge had already faded when the trace reached the edge: the signal ended inside, a free end
            stop = "free"
        if stop == "free":
            pts = pts[: last_strong + 1]
            strg = strg[: last_strong + 1]
        return pts, strg, stop

    def stamp(pts, sid):
        if claim_t is not None and len(pts) > 1:
            P = np.asarray(pts, float)
            tg = np.gradient(P, axis=0)
            tg /= np.maximum(np.linalg.norm(tg, axis=1, keepdims=True), 1e-12)
        else:
            tg = None
        for i, q in enumerate(pts):
            v = vox(q)
            z0, z1 = max(v[0] - cr[0], 0), min(v[0] + cr[0] + 1, Z)
            y0, y1 = max(v[1] - cr[1], 0), min(v[1] + cr[1] + 1, Y)
            x0, x1 = max(v[2] - cr[2], 0), min(v[2] + cr[2] + 1, X)
            if z0 >= z1 or y0 >= y1 or x0 >= x1:
                continue
            sub = claim[z0:z1, y0:y1, x0:x1]
            free = sub == 0
            sub[free] = sid + 1
            if tg is not None:
                claim_t[z0:z1, y0:y1, x0:x1][free] = tg[i]

    cand = np.argwhere(F.crest & (S >= seed_min))
    if len(cand):
        sv = S[tuple(cand.T)]
        flat = np.ravel_multi_index(tuple(cand.T), S.shape)
        order = np.lexsort((flat, -sv))
        cand = cand[order]
    info["n_seeds"] = int(len(cand))
    strands: list[RawStrand] = []
    sid = 0
    for c in cand:
        c = (int(c[0]), int(c[1]), int(c[2]))
        if claim[c]:
            continue
        if info["budget_exceeded"]:
            break
        oc = off[c]
        p0 = tuple((c[k] + max(-0.5, min(0.5, float(oc[k]))) + o[k]) * sp[k] for k in range(3))
        t0 = phys(D[c])
        if t0 == (0.0, 0.0, 0.0):
            continue
        own: dict = {}
        fwd, sf, s1 = trace(p0, t0, sid, own, +1.0, max_len)
        len_f = (len(fwd) - 1) * h
        bwd, sb, s0 = trace(p0, (-t0[0], -t0[1], -t0[2]), sid, own, -1.0, max(max_len - len_f, h))
        pts = np.asarray(bwd[::-1] + fwd[1:], float)
        strg = np.asarray(sb[::-1] + sf[1:], float)
        stamp(pts, sid)
        info["n_traced"] += 1
        L = float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum()) if len(pts) > 1 else 0.0
        if L >= min_len:
            et = [("max_length" if s == "budget" else s) for s in (s0, s1)]
            strands.append(RawStrand(points=pts, end_types=et, strength=strg))
        sid += 1
    # T contacts are resolved by the caller after the axial shadows are removed (trace.trace_crop)
    return strands, info


def _scaled(points, zl):
    P = np.asarray(points, float)
    return np.stack([P[:, 0] / zl[0], P[:, 1] / zl[1], P[:, 2] / zl[1]], 1)


def resolve_contacts(strands: list[RawStrand], p, end_window_um: float = 0.5, other_end_guard_um: float = 0.5):
    """An end whose last ``end_window_um`` runs inside another strand's tube (``follow.contact_um``, (z, lateral)),
    away from that strand's own ends, is a T contact: trim it to the tube entry and call it ``contact``."""
    if len(strands) < 2:
        return
    zl = [float(v) for v in p["follow.contact_um"]]
    trees, arcs = [], []
    for s in strands:
        trees.append(cKDTree(_scaled(s.points, zl)))
        arcs.append(arclen(s.points))
    for i, st in enumerate(strands):
        for side in (0, 1):
            if st.end_types[side] == "contact":
                continue
            P = st.points if side == 1 else st.points[::-1]
            s = arclen(P)
            Ltot = s[-1]
            win = s >= Ltot - end_window_um
            idx = np.nonzero(win)[0]
            if len(idx) == 0:
                continue
            Q = _scaled(P[idx], zl)
            hit_first = None
            for j, tr in enumerate(trees):
                if j == i:
                    continue
                d, k = tr.query(Q, distance_upper_bound=1.0)
                ok = np.isfinite(d)
                if not ok.any():
                    continue
                aj = arcs[j]
                # ignore contacts near the other strand's own ends (an end-to-end meeting, for the linker)
                near_end = (aj[np.minimum(k[ok], len(aj) - 1)] < other_end_guard_um) | \
                           (aj[np.minimum(k[ok], len(aj) - 1)] > aj[-1] - other_end_guard_um)
                good = np.nonzero(ok)[0][~near_end]
                if len(good) == 0:
                    continue
                first = int(idx[good[0]])
                if hit_first is None or first < hit_first:
                    hit_first = first
            if hit_first is None or hit_first < 1:
                continue
            P2 = P[: hit_first + 1]
            strg = st.strength if st.strength is not None else None
            if side == 1:
                st.points = P2.copy()
                if strg is not None:
                    st.strength = strg[: hit_first + 1].copy()
            else:
                st.points = P2[::-1].copy()
                if strg is not None:
                    st.strength = strg[::-1][: hit_first + 1][::-1].copy()
            st.end_types[side] = "contact"
            trees[i] = cKDTree(_scaled(st.points, zl))
            arcs[i] = arclen(st.points)


def orphan_contact_ends(ends: list[tuple[np.ndarray, np.ndarray]], reach_um=(0.6, 0.4)) -> list[bool]:
    """For each (end point, other points) pair: True when no other point lies inside the ellipsoid of half axes
    ``reach_um`` (axial, lateral) around the end point (its contact partner is gone)."""
    ax, lat = (float(v) for v in reach_um)
    out = []
    for e, Q in ends:
        Q = np.asarray(Q, float).reshape(-1, 3)
        if not len(Q):
            out.append(True)
            continue
        d = Q - np.asarray(e, float)[None, :]
        out.append(not bool(np.any((d[:, 0] / ax) ** 2 + (d[:, 1] ** 2 + d[:, 2] ** 2) / lat ** 2 <= 1.0)))
    return out


def free_orphan_contacts(strands: list[RawStrand], p, others: np.ndarray | None = None) -> int:
    """A ``contact`` end with no other strand (of ``strands``, or a point of ``others``) within
    ``follow.contact_free_um`` (axial, lateral) becomes ``free``: the strand it stopped against was removed later
    (an axial shadow, a duplicate, a shadowed strand), so the end rule applies to it. Returns the number freed."""
    reach = tuple(float(v) for v in p["follow.contact_free_um"])
    n = 0
    for i, st in enumerate(strands):
        sides = [s for s in (0, 1) if st.end_types[s] == "contact"]
        if not sides or len(st.points) == 0:
            continue
        rest = [o.points for j, o in enumerate(strands) if j != i and len(o.points)]
        if others is not None and len(others):
            rest.append(np.asarray(others, float))
        Q = np.vstack(rest) if rest else np.zeros((0, 3))
        flags = orphan_contact_ends([(st.points[0 if s == 0 else -1], Q) for s in sides], reach)
        for s, orphan in zip(sides, flags):
            if orphan:
                st.end_types[s] = "free"
                n += 1
    return n
