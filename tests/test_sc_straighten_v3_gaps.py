"""Dim stretches that are not gaps (protocol changelog entry 33), off by default.

A gap candidate that passes every test of v1 is a dim stretch, not a gap, when no sample of its widened window lies
next to another strand or at a link between two traced members, its relative minimum lies above the minimum of an
empty gap of its missing length after the blur along s (with the image's measured axial sigma) by more than
gaps.dim_z standard errors, and its fill fraction exceeds gaps.dim_fill. The rule changes which candidates are gaps
and nothing else: the strands and every sample except in_gap are those of v1. An empty gap stays a gap (also between
two ribbons, and in a steep ribbon under a PSF wider than the fallback); a dim stretch of known length and
brightness is read as one, with its length. The known ambiguity (two ribbons with dim tips and a short empty core
that the tracer follows) is pinned at the top of the grid entry 33 chooses from.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
import sc_v3_phantoms as Ph

from germquant.sc.straighten import params, synth, synth_dev
from germquant.sc.straighten.measure import empty_gap_min, fill_fraction, find_gaps

WL = Ph.WL
P_V1 = params.resolve(None, params.v1_sets())
P_GAP = params.resolve(None, Ph.sets(params.V3_GAPS))
P_GAP_PSF = params.resolve(None, Ph.sets({**params.V3_GAPS, "gaps.dim_sigma_source": "psf"}))


# ------------------------------------------------------------------------------------ the model
def test_the_empty_gap_model_and_the_fill_fraction():
    sig = 0.2
    for g in (0.3, 0.5, 1.0):
        r = empty_gap_min(g, sig)
        assert r == pytest.approx(1.0 - math.erf(g / (2.0 * math.sqrt(2.0) * sig)))
        assert fill_fraction(g, r, sig) == pytest.approx(0.0, abs=1e-6)
    # a stretch at 40 percent over 1.2 um: missing length 0.72 um, its minimum after the blur
    f, ell = 0.4, 1.2
    r = f + (1 - f) * empty_gap_min(ell, sig)
    assert fill_fraction((1 - f) * ell, r, sig) == pytest.approx(f, abs=0.01)
    assert fill_fraction(0.5, 0.0, sig) == 0.0           # deeper than an empty gap: no fill


def _profile(s, events, sig, level=100.0):
    """A flat integrated profile with stretches [(s0, s1, fill)] blurred along s with ``sig``."""
    from scipy.special import erf

    y = np.full(len(s), level)
    for s0, s1, f in events:
        y -= level * (1 - f) * 0.5 * (erf((s - s0) / (math.sqrt(2) * sig)) - erf((s - s1) / (math.sqrt(2) * sig)))
    return y


def test_find_gaps_keeps_a_dim_stretch_and_an_empty_gap_apart():
    ds, L = 0.05, 8.0
    s = np.arange(-0.5, L + 0.5 + 1e-9, ds)
    ins = (s >= 0) & (s <= L)
    sig = 0.15
    rng = np.random.default_rng(0)
    y = _profile(s, [(2.0, 3.0, 0.4), (5.5, 6.0, 0.0)], sig) + rng.normal(0, 1.0, len(s))
    se = np.full(len(s), 1.5)
    clean = np.ones(len(s), bool)
    g_v1, _M, _sm = find_gaps(s, y, se, clean, ins, L, P_V1, ds, sigma_along=lambda _s: sig)
    assert len(g_v1) == 2                                 # v1 reads the dim stretch as a gap as well
    dims = []
    g3, _M, _sm = find_gaps(s, y, se, clean, ins, L, P_GAP, ds, sigma_along=lambda _s: sig, dim=dims)
    assert len(g3) == 1 and len(dims) == 1
    assert 5.4 < 0.5 * (g3[0]["s_start_um"] + g3[0]["s_end_um"]) < 6.1
    assert g3[0]["length_um"] == pytest.approx(0.5, abs=0.1)
    assert dims[0]["fill"] == pytest.approx(0.4, abs=0.1)
    assert 2.0 < 0.5 * (dims[0]["s_start_um"] + dims[0]["s_end_um"]) < 3.0
    # off: the list stays empty and the gaps are v1's
    d0 = []
    g_off, _M, _sm = find_gaps(s, y, se, clean, ins, L, P_V1, ds, sigma_along=lambda _s: sig, dim=d0)
    assert d0 == [] and [g["s_start_um"] for g in g_off] == [g["s_start_um"] for g in g_v1]


@pytest.mark.parametrize("g", [0.3, 0.5, 0.9])
@pytest.mark.parametrize("sig", [0.15, 0.25, 0.36])
def test_empty_gaps_are_never_dim_stretches(g, sig):
    ds, L = 0.05, 6.0
    s = np.arange(-0.5, L + 0.5 + 1e-9, ds)
    ins = (s >= 0) & (s <= L)
    rng = np.random.default_rng(int(100 * g + 1000 * sig))
    y = _profile(s, [(3.0 - g / 2, 3.0 + g / 2, 0.0)], sig) + rng.normal(0, 2.0, len(s))
    se = np.full(len(s), 2.5)
    kw = dict(sigma_along=lambda _s: sig)
    a, _, _ = find_gaps(s, y, se, np.ones(len(s), bool), ins, L, P_V1, ds, **kw)
    dims = []
    b, _, _ = find_gaps(s, y, se, np.ones(len(s), bool), ins, L, P_GAP, ds, dim=dims, **kw)
    assert len(b) == len(a) and dims == []


def test_no_dim_stretch_across_a_link_or_next_to_another_strand():
    """A candidate that the model would call a dim stretch stays a gap when one sample of its widened window is
    blocked (a link between traced members, or another strand within contamination_um)."""
    ds, L, sig = 0.05, 8.0, 0.15
    s = np.arange(-0.5, L + 0.5 + 1e-9, ds)
    ins = (s >= 0) & (s <= L)
    y = _profile(s, [(3.5, 4.5, 0.4)], sig) + np.random.default_rng(5).normal(0, 1.0, len(s))
    se = np.full(len(s), 1.5)
    clean = np.ones(len(s), bool)
    kw = dict(sigma_along=lambda _s: sig)
    free = []
    g_free, _, _ = find_gaps(s, y, se, clean, ins, L, P_GAP, ds, dim=free, **kw)
    assert len(free) == 1 and not g_free
    for where in (3.45, 4.0, 4.55):
        dims = []
        g, _, _ = find_gaps(s, y, se, clean, ins, L, P_GAP, ds, dim=dims, dim_block=np.abs(s - where) <= 0.1, **kw)
        assert dims == [] and len(g) == 1, where
    dims = []                                             # a block far from the stretch changes nothing
    g, _, _ = find_gaps(s, y, se, clean, ins, L, P_GAP, ds, dim=dims, dim_block=np.abs(s - 7.0) <= 0.1, **kw)
    assert len(dims) == 1 and not g


# ------------------------------------------------------------------------------------ phantoms
def test_the_rule_changes_gaps_and_nothing_else():
    a = Ph.line(-3.2, 3.2, 10.0)
    b = Ph.line(-2.6, 2.6, 35.0, origin=(0.0, 1.6, 0.4))
    sc = Ph.scene([{"points": a, "peak": 520.0, "weights": Ph.weights(a, [(1.6, 2.6, 0.35)])},
                   {"points": b, "peak": 480.0, "gaps": [(2.4, 2.9)]}], seed=31, R=3.4)
    v1 = Ph.run(sc, P_V1)
    v3 = Ph.run(sc, P_GAP)
    Ph.assert_same_strands(v1, v3)
    assert v3.record["dim_stretches"]["n"] >= 1 and len(v3.gaps) < len(v1.gaps)
    # every dim stretch is a candidate v1 read as a gap: same strand, overlapping s interval
    for _, d in v3.dim_stretches.iterrows():
        g = v1.gaps[v1.gaps["strand_id"] == d["strand_id"]]
        assert ((g["s_start_um"] < d["s_end_um"]) & (g["s_end_um"] > d["s_start_um"])).any()


def test_a_dim_stretch_of_known_length_and_brightness():
    """1.2 um at 40 percent of the ribbon's brightness in the image plane: one dim stretch whose missing length is
    the stretch's (0.6 x 1.2 = 0.72 um) and whose fill is about 0.4; the strand continues."""
    a = Ph.line(-3.4, 3.4, 5.0)
    sc = Ph.scene([{"points": a, "peak": 600.0, "weights": Ph.weights(a, [(2.8, 4.0, 0.4)])}], seed=81, R=3.4)
    v1 = Ph.run(sc, P_V1)
    v3 = Ph.run(sc, P_GAP)
    assert len(v1.gaps) == 1
    assert len(v3.gaps) == 0 and len(v3.dim_stretches) == 1
    d = v3.dim_stretches.iloc[0]
    assert d["length_um"] == pytest.approx(0.72, abs=0.15)
    assert d["fill"] == pytest.approx(0.4, abs=0.12)
    st = v3.strands.iloc[0]
    centre = 0.5 * (d["s_start_um"] + d["s_end_um"])
    assert min(abs(centre - 3.4), abs(centre - (float(st["length_um"]) - 3.4))) < 0.35
    assert v3.nuclei["piece_lw_mean_um"].iloc[0] == pytest.approx(float(st["length_um"]), rel=1e-9)


def test_a_dim_stretch_in_an_arc_continues_the_ribbon():
    arc = synth._densify(synth.arc(6.0, 6.0, 0.0, 0.3))
    sc = Ph.scene([{"points": arc, "peak": 500.0, "weights": Ph.weights(arc, [(2.5, 3.5, 0.35)])}], seed=11)
    v1 = Ph.run(sc, P_V1)
    assert len(v1.strands) == 1 and int(v1.strands["n_gaps"].iloc[0]) >= 1
    v3 = Ph.run(sc, P_GAP)
    st, nu = v3.strands, v3.nuclei.iloc[0]
    assert len(st) == 1 and int(st["n_gaps"].iloc[0]) == 0 and len(v3.gaps) == 0
    assert nu["piece_lw_mean_um"] == pytest.approx(float(st["length_um"].iloc[0]), rel=1e-9)
    assert v3.record["dim_stretches"]["n"] >= 1
    assert float(st["length_um"].iloc[0]) == pytest.approx(float(v1.strands["length_um"].iloc[0]), abs=1e-9)


@pytest.mark.parametrize("tilt", [0.0, 45.0])
def test_an_empty_gap_of_known_length_stays_a_gap(tilt):
    arc = synth._densify(synth.arc(6.0, 6.0, tilt, 0.3))
    sc = Ph.scene([{"points": arc, "peak": 500.0, "gaps": [(2.75, 3.35)]}], seed=12)
    v1 = Ph.run(sc, P_V1)
    v3 = Ph.run(sc, P_GAP)
    assert len(v1.gaps) >= 1
    assert len(v3.gaps) == len(v1.gaps)
    assert float(v3.gaps["length_um"].max()) == pytest.approx(0.6, abs=0.15)


@pytest.mark.parametrize("gap_um", [0.35, 0.6])
@pytest.mark.parametrize("tilt", [0.0, 40.0])
def test_an_empty_gap_between_two_ribbons_stays_a_gap(gap_um, tilt):
    a = Ph.line(-3.0 - gap_um / 2, -gap_um / 2, tilt)
    b = Ph.line(gap_um / 2, 3.0 + gap_um / 2, tilt)
    sc = Ph.scene([{"points": a, "peak": 520.0}, {"points": b, "peak": 520.0}], seed=int(40 + 100 * gap_um + tilt),
                  R=3.4)
    v1 = Ph.run(sc, P_V1)
    v3 = Ph.run(sc, P_GAP)
    Ph.assert_same_strands(v1, v3)
    assert v3.record["dim_stretches"]["n"] == 0
    assert len(v3.gaps) == len(v1.gaps)
    assert v3.nuclei["n_pieces"].iloc[0] == v1.nuclei["n_pieces"].iloc[0]


@pytest.mark.parametrize("flank", [0.4, 0.6])
@pytest.mark.parametrize("tilt", [0.0, 40.0])
def test_two_ribbons_with_dim_tips(tilt, flank):
    """The known ambiguity of the rule (entry 33): the tips of two different ribbons dimmed to ``flank`` of their
    brightness over their last 0.4 um with an empty 0.3 um between them, which the tracer follows as one strand.
    Behind the blur along s the empty core does not show, and the junction looks like a stretch filled to 0.1 to
    0.3 of its flanks: gaps.dim_fill 0.35, the top of the grid that entry 33 chooses from, keeps every such junction
    a gap, and their rate in the dense phantoms is a targeted metric of the entry."""
    g = 0.3
    a = Ph.line(-3.0 - g / 2, -g / 2, tilt)
    b = Ph.line(g / 2, 3.0 + g / 2, tilt)
    sc = Ph.scene([{"points": a, "peak": 520.0, "weights": Ph.weights(a, [(2.6, 3.0, flank)])},
                   {"points": b, "peak": 520.0, "weights": Ph.weights(b, [(0.0, 0.4, flank)])}], seed=int(60 + tilt),
                  R=3.4)
    v1 = Ph.run(sc, P_V1)
    v3 = Ph.run(sc, P_GAP)
    Ph.assert_same_strands(v1, v3)
    if len(v3.dim_stretches):
        assert float(v3.dim_stretches["fill"].max()) < 0.35
    strict = Ph.run(sc, params.resolve(None, Ph.sets({**params.V3_GAPS, "gaps.dim_fill": 0.35})))
    Ph.assert_same_strands(v1, strict)
    assert strict.record["dim_stretches"]["n"] == 0
    assert strict.nuclei["n_pieces"].iloc[0] == v1.nuclei["n_pieces"].iloc[0]


def test_recorded_dim_stretches_avoid_neighbours():
    """On a phantom with several ribbons, no recorded dim stretch overlaps a contaminated sample."""
    a = Ph.line(-3.2, 3.2, 15.0)
    b = Ph.line(-2.8, 2.8, 50.0, origin=(0.6, 0.25, 0.0))     # crosses a near its middle, 0.6 um higher
    c = Ph.line(-2.5, 2.5, 5.0, origin=(-0.8, -1.8, 0.3))
    sc = Ph.scene([{"points": a, "peak": 520.0, "weights": Ph.weights(a, [(1.2, 2.0, 0.35), (2.8, 3.6, 0.3)])},
                   {"points": b, "peak": 480.0, "weights": Ph.weights(b, [(2.2, 3.0, 0.3)])},
                   {"points": c, "peak": 450.0, "weights": Ph.weights(c, [(1.5, 2.4, 0.35)])}], seed=71, R=3.4)
    v3 = Ph.run(sc, P_GAP)
    smp = v3.samples
    for _, d in v3.dim_stretches.iterrows():
        m = smp[(smp["strand_id"] == d["strand_id"]) & (smp["s_um"] >= d["s_start_um"]) &
                (smp["s_um"] <= d["s_end_um"])]
        assert not m["contaminated"].astype(bool).any()


@pytest.mark.parametrize("tilt", [50.0, 70.0])
def test_an_empty_gap_in_a_steep_ribbon_under_a_wide_psf_stays_a_gap(tilt):
    """sigma_z 0.50 um (the fallback is 0.365): the empty-gap model with the image's measured axial sigma keeps a
    0.5 um empty gap of a ribbon tilted 50 or 70 degrees a gap."""
    a = Ph.line(-3.0, 3.0, tilt)
    sc = Ph.scene([{"points": a, "peak": 650.0, "gaps": [(2.75, 3.25)]}], sigma_z=0.50, foci=28, seed=int(90 + tilt),
                  R=3.4)
    v1 = Ph.run(sc, P_V1, WL)
    v3 = Ph.run(sc, P_GAP, WL)
    ref = v3.record["dim_reference"]
    assert ref["source"] in ("foci", "foci_depth") and ref["sigma_z_um"] == pytest.approx(0.50, rel=0.12)
    assert len(v1.gaps) >= 1
    assert len(v3.gaps) == len(v1.gaps) and v3.record["dim_stretches"]["n"] == 0
    # the image PSF itself is that of v1
    assert v3.record["psf"]["sigma_z_um"] == v1.record["psf"]["sigma_z_um"]


def test_the_gap_reference_without_foci_is_the_image_psf():
    a = Ph.line(-3.0, 3.0, 10.0)
    sc = Ph.scene([{"points": a, "peak": 500.0, "weights": Ph.weights(a, [(2.5, 3.5, 0.35)])}], seed=95, R=3.4)
    v3 = Ph.run(sc, P_GAP)
    assert v3.record["dim_reference"]["source"] == "psf"
    psf_src = Ph.run(sc, P_GAP_PSF)
    pd.testing.assert_frame_equal(v3.gaps, psf_src.gaps)
    pd.testing.assert_frame_equal(v3.dim_stretches, psf_src.dim_stretches)


# ------------------------------------------------------------------------------------ the development family
def test_the_gaps_psf_family():
    a = synth_dev.make_scene("gaps_psf", 2001)
    b = synth_dev.make_scene("gaps_psf", 2001)
    assert np.array_equal(a.channels["central_element"], b.channels["central_element"])
    assert np.array_equal(a.channels["foci"], b.channels["foci"])
    assert set(a.channels) == {"central_element", "foci"} and len(a.spots) > 20 and a.info["wavelengths"] == WL
    kinds = [r.kind for r in a.ribbons]
    assert kinds.count("gap_ribbon") == 6 and kinds.count("dim_ribbon") == 6
    for r in a.ribbons:
        if r.kind == "gap_ribbon":
            assert r.gaps and all(0.3 <= g1 - g0 <= 1.0 for g0, g1 in r.gaps)
            assert r.ribbon_id not in a.info["dim_stretches"]
        else:
            assert not r.gaps
            ev = a.info["dim_stretches"][r.ribbon_id]
            assert ev and all(0.4 <= e[1] - e[0] <= 1.5 and 0.15 <= e[2] <= 0.6 for e in ev)
            assert ev[0][0] >= 0.8 - 1e-9 and ev[-1][1] <= r.length_um - 0.8 + 1e-9
    assert a.info["sigma_z_level"] in synth_dev.PSF_LEVELS
    assert abs(a.info["sigma_z_gradient_per_um"]) <= synth_dev.PUNCTA_PSF_GRADIENT
    with pytest.raises(synth.SealedSeedError):
        synth_dev.make_scene("gaps_psf", 7_000_001)
    with pytest.raises(synth.SealedSeedError):
        synth_dev.make_scene("gaps_psf", 2001, final=True)
