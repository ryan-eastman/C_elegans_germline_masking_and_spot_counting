"""Decomposition changes of sc_straighten v2 (protocol changelog entries 20 on): the contact stop against the
claiming strand's tangent (entry 20, withdrawn by entry 21: a switch, off by default) and the end matching (entries
21 and 22; the chain bound of entry 21 is a switch, off by default)."""
from __future__ import annotations

import numpy as np
import pytest

from germquant.benchmark import sc_ribbons as B
from germquant.sc.straighten import params, synth
from germquant.sc.straighten import trace as TR
from germquant.sc.straighten.pipeline import straighten_image


def _run(sc, sets=()):
    return straighten_image(sc.channels["central_element"], sc.labels, sc.spacing, sc.nuclei,
                            params=params.resolve(None, list(sets)), image_id="t")


def _traced(monkeypatch, sc, sets):
    """The strands of the tracer itself (before T contacts, end matching and accounting)."""
    cap = []
    orig = TR.follow_nucleus

    def hook(*a, **k):
        out = orig(*a, **k)
        cap.append([(s.length_um, list(s.end_types)) for s in out[0]])
        return out

    monkeypatch.setattr(TR, "follow_nucleus", hook)
    res = _run(sc, sets)
    return cap[0], res


def _x_crossing(angle_deg: float, dz: float, seed: int):
    a = synth.arc(6.0, 8.0, 0.0, 0.3)
    b = synth.arc(6.0, 8.0, 0.0, 0.3 + np.radians(angle_deg))
    a[:, 0] -= dz / 2
    b[:, 0] += dz / 2
    return synth.custom_scene([{"points": a, "peak": 430.0}, {"points": b, "peak": 430.0}], seed=seed)


@pytest.mark.parametrize("angle", [64.0, 80.0])
def test_trace_passes_a_crossing_against_the_claiming_tangent(monkeypatch, angle):
    sc = _x_crossing(angle, 0.3, 2111)
    new, res = _traced(monkeypatch, sc, ["follow.contact_test=claim_tangent"])
    assert sum(1 for L, _e in new if L >= 6.0) >= 1               # the second ribbon is traced whole
    old, res_v1 = _traced(monkeypatch, sc, ["follow.contact_test=local"])
    long_old = [L for L, _e in old if L >= 6.0]
    assert sum(1 for L, e in old if "contact" in e and L < 3.0) >= 2 and len(long_old) <= 1   # v1 breaks it
    for r in (res, res_v1):                                         # the end matching joins v1's pieces again
        m = B.score_scene(sc, r).matches
        assert all(x["matched"] and x["cover"] >= 0.8 and x["purity"] >= 0.8 for x in m)


def test_parallel_pair_and_t_contact_are_unchanged():
    # a parallel pair 0.25 um apart in xy: no duplicated length
    a = synth.arc(4.0, 8.0, 0.0, 0.0)
    off = np.array([0.0, 0.25, 0.0])
    sc = synth.custom_scene([{"points": a - off / 2, "peak": 430.0}, {"points": a + off / 2, "peak": 430.0}],
                            seed=2121)
    res = _run(sc)
    assert res.strands["length_um"].sum() <= 1.15 * 2 * 4.0
    # the T contact still ends the stem at the bar with a contact end
    bar = synth.arc(6.0, 8.0, 0.0, 0.0)
    mid = bar[len(bar) // 2]
    s = np.linspace(0, 3.0, 601)
    stem = mid[None, :] + s[:, None] * np.array([0.0, 1.0, 0.0])[None, :]
    sc = synth.custom_scene([{"points": bar, "peak": 430.0}, {"points": stem, "peak": 430.0}], seed=104)
    res = _run(sc)
    assert "contact" in set(res.strands["end0_class"]) | set(res.strands["end1_class"])


# ------------------------------------------------------------------------------- the end matching (entries 21, 22)
def _poly(pts, n=30):
    pts = np.asarray(pts, float)
    out = [pts[0]]
    for a, b in zip(pts[:-1], pts[1:]):
        out.extend(a + np.linspace(0, 1, n)[1:, None] * (b - a))
    return np.asarray(out)


def test_ends_hooked_toward_a_crossing_ribbon_are_linked():
    from germquant.sc.straighten import link
    from germquant.sc.straighten.follow import RawStrand

    A = RawStrand(_poly([(2, 0, 0), (2, 0, 2.5), (2, 0.8, 2.65)]), ["free", "free"])
    Bs = RawStrand(_poly([(2, 0.8, 2.75), (2, 0, 2.9), (2, 0, 5.4)]), ["free", "free"])
    C = RawStrand(_poly([(2, -2.5, 2.7), (2, 2.5, 2.7)]), ["free", "free"])
    new = link.link_strands([A, C, Bs], params.resolve())
    old = link.link_strands([A, C, Bs], params.resolve(None, ["link.junction_dirs=false"]))
    assert sorted(len(s.links) for s in new) == [0, 1] and len(new) == 2
    assert len(old) == 3 and all(len(s.links) == 0 for s in old)


def test_pairs_beyond_a_junction_are_judged_from_the_tips():
    # entry 22: two ends measured beyond a junction whose tips lie 1.2 um apart: the anchor estimate of entry 21
    # gives them a gap of 0 and links them; the distance between the tips does not
    from germquant.sc.straighten import link
    from germquant.sc.straighten.follow import RawStrand

    A = RawStrand(_poly([(2, 0, 0), (2, 0, 2.5), (2, 0.8, 2.65)]), ["free", "free"])
    C = RawStrand(_poly([(2, -2.5, 2.7), (2, 2.5, 2.7)]), ["free", "free"])
    D = RawStrand(_poly([(2, 2.0, 2.75), (2, 1.6, 2.95), (2, 1.6, 5.4)]), ["free", "free"])
    new = link.link_strands([A, C, D], params.resolve())
    old = link.link_strands([A, C, D], params.resolve(None, ["link.junction_gap=anchors"]))
    assert len(new) == 3 and all(len(s.links) == 0 for s in new)
    assert len(old) == 2 and sorted(len(s.links) for s in old) == [0, 1]


@pytest.mark.parametrize("pieces,n_out", [([(0, 3.4), (3.7, 7.1), (7.4, 10.5)], 2), ([(0, 2.8), (3.1, 5.9), (6.2, 8.9)], 1)])
def test_chain_bound(pieces, n_out):
    from germquant.sc.straighten import link
    from germquant.sc.straighten.follow import RawStrand

    strands = [RawStrand(_poly([(2, 2, a), (2, 2, b)]), ["free", "free"]) for a, b in pieces]
    out = link.link_strands(strands, params.resolve(None, ["link.max_chain_um=9.0"]))
    assert len(out) == n_out and max(s.length_um for s in out) <= 9.0 + 1e-6
    if n_out == 2:
        assert all(e == "free" for s in out for e in s.end_types)
    # withdrawn by entry 22: no bound by default
    assert len(link.link_strands(strands, params.resolve())) == 1
