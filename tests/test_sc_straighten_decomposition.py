"""Decomposition changes of sc_straighten v2 (protocol changelog entries 20 on): the contact stop against the
claiming strand's tangent."""
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
    new, res = _traced(monkeypatch, sc, [])
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
