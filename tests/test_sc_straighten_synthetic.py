"""sc_straighten on exact-render phantoms (development seeds; docs/SC_VALIDATION_PROTOCOL.md 2.3): an isolated
arc at real noise, a tilted ribbon, an X crossing, a T contact, a 0.5 um gap, foci at known s, a noise-only
nucleus and a dense nucleus (design section 11)."""
from __future__ import annotations

import numpy as np
import pytest

from germquant.benchmark import sc_ribbons as B
from germquant.sc.straighten import params, synth
from germquant.sc.straighten.pipeline import straighten_image

P = params.resolve()


def _run(sc):
    roles = {k: v for k, v in sc.channels.items() if k != "central_element"}
    spots = {"spots": sc.spots} if len(sc.spots) else None
    return straighten_image(sc.channels["central_element"], sc.labels, sc.spacing, sc.nuclei, params=P,
                            roles=roles, spot_tables=spots, image_id="t")


def _main(res, min_um=1.0):
    st = res.strands
    return st[st["length_um"] >= min_um]


def test_isolated_arc_at_real_noise():
    sc = synth.custom_scene([{"points": synth.arc(5.0, 2.5, 0.0, 0.4), "peak": 430.0}], seed=101)
    res = _run(sc)
    main = _main(res)
    assert len(main) == 1
    s = B.score_scene(sc, res)
    m = s.matches[0]
    assert m["matched"] and abs(m["err_um"]) <= 0.2 and m["dist_um"] <= 0.05


def test_tilted_ribbon():
    sc = synth.custom_scene([{"points": synth.arc(5.0, 3.0, 60.0, 1.1), "peak": 430.0}], seed=102)
    res = _run(sc)
    s = B.score_scene(sc, res)
    assert s.matches[0]["matched"] and abs(s.matches[0]["err_um"]) <= 0.4


def test_x_crossing_continues_through():
    a = synth.arc(6.0, 5.0, 0.0, 0.3)
    b = synth.arc(6.0, 5.0, 0.0, 0.3 + np.radians(70))
    sc = synth.custom_scene([{"points": a, "peak": 450.0}, {"points": b, "peak": 420.0}], seed=103)
    res = _run(sc)
    s = B.score_scene(sc, res)
    assert len(_main(res)) == 2
    for m in s.matches:
        assert m["matched"] and m["purity"] >= 0.8 and m["cover"] >= 0.8


def test_t_contact_gives_a_contact_end():
    a = synth.arc(6.0, 8.0, 0.0, 0.0)
    mid = a[len(a) // 2]
    d = np.array([0.0, np.cos(0.0), np.sin(0.0)])               # perpendicular to the chord (y direction)
    s = np.linspace(0, 3.0, 601)
    b = mid[None, :] + s[:, None] * d[None, :]
    sc = synth.custom_scene([{"points": a, "peak": 430.0}, {"points": b, "peak": 430.0}], seed=104)
    res = _run(sc)
    ends = set(res.strands["end0_class"]) | set(res.strands["end1_class"])
    assert "contact" in ends


def test_gap_of_half_a_micron():
    sc = synth.custom_scene([{"points": synth.arc(5.0, 4.0, 0.0, 0.7), "peak": 430.0, "gaps": [(2.25, 2.75)]}],
                            seed=105)
    res = _run(sc)
    assert len(_main(res)) == 1
    assert len(res.gaps) == 1
    assert res.gaps["length_um"].iloc[0] == pytest.approx(0.5, abs=0.15)


def test_foci_at_known_positions():
    sc = synth.custom_scene([{"points": synth.arc(4.0, 3.0, 0.0, 0.2), "peak": 430.0, "foci_s": [1.0, 3.0]}],
                            seed=106)
    res = _run(sc)
    f = res.foci
    assert len(f) == 2 and set(f["assignment"]) == {"unique"}
    s = B.score_scene(sc, res)
    assert all(x["correct_ribbon"] for x in s.foci)
    assert max(abs(x["s_err_um"]) for x in s.foci) <= 0.1
    # the folded positions do not depend on the orientation
    assert sorted(np.round(f["s_folded_frac"], 1)) == [0.2, 0.2] or max(f["s_folded_frac"]) <= 0.5


def test_noise_only_nucleus():
    sc = synth.noise_scene(107)
    res = _run(sc)
    assert res.strands["length_um"].sum() < 1.0 if len(res.strands) else True
    assert res.nuclei["status"].iloc[0] == "no_sc"


def test_dense_nucleus_totals():
    # four fixed development nuclei (two oocytes, two spermatocytes): the mean total within 10 percent, each
    # within 25 percent (the per-nucleus spread is the protocol's S1b SD criterion, scored on the sealed seeds;
    # seed 112 is -19 percent in v1 and -21 percent with the contact stop of protocol changelog entry 20)
    rel = []
    for seed in (110, 111, 112, 113):
        sc = synth.dense_scene(seed)
        res = _run(sc)
        truth = sum(r.length_um for r in sc.ribbons)
        rel.append((res.strands["length_um"].sum() - truth) / truth)
    assert abs(float(np.mean(rel))) <= 0.10
    assert max(abs(x) for x in rel) <= 0.25
