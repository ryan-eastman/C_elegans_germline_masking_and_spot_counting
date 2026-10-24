"""The punctum class of sc_straighten v2 (protocol changelog entries 16 to 18): the axial blur of a SYP punctum is
classed and flagged punctum (by default a label, entry 18; under punctum.exclude, entry 16, it also leaves the
nucleus tier and the foci assignment), a tilted ribbon is not, and the rule before the end matching (entry 17)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from germquant.sc.straighten import params, synth, synth_dev
from germquant.sc.straighten.pipeline import straighten_image

P = params.resolve(None, ["punctum.exclude=true"])        # the tier exclusion of entry 16 (entry 18: off by default)


def _scene(ribbons=(), puncta=(), seed=0, foci=()):
    """One nucleus holding the given ribbons (dicts with points and peak) and point puncta (centre offsets from the
    nucleus centre, peaks); optional foci-channel spots at given positions (also offsets)."""
    rng = np.random.default_rng(seed)
    R = 2.4
    c0 = np.zeros(3)
    ribs = [synth.TruthRibbon(ribbon_id=k + 1, nucleus_id=1, points=synth._densify(np.asarray(r["points"], float)),
                              peak=float(r.get("peak", 430.0)), kind="custom") for k, r in enumerate(ribbons)]
    pts = [np.asarray(q["center"], float) for q in puncta]
    shape, shift = synth._place([r.points for r in ribs] + ([np.vstack(pts)] if pts else []), [c0], [R + 0.3])
    for r in ribs:
        r.points = r.points + shift
    center = c0 + shift
    labels, _ = synth._sphere_masks(shape, synth.SPACING, [center], [R + 0.3])
    nuc = synth._nucleoplasm(shape, synth.SPACING, [center], [R], [60.0])
    clean = np.full(shape, synth.OFFSET) + nuc
    for r in ribs:
        synth.render_ribbon(clean, r)
    for q, c in zip(puncta, pts):
        synth._render_points(clean, (c + shift)[None, :], synth.SPACING, synth.PSF_SIGMA, weight=float(q["peak"]))
    ce = synth.add_noise(clean, rng)
    channels = {"central_element": ce}
    spots = None
    if foci:
        fclean = np.full(shape, synth.OFFSET) + nuc
        rows = []
        for k, f in enumerate(foci):
            p = np.asarray(f, float) + shift
            synth._render_points(fclean, p[None, :], synth.SPACING, synth.PSF_SIGMA, weight=300.0)
            rows.append({"spot_id": k, "nucleus_id": 1, "z_um": p[0], "y_um": p[1], "x_um": p[2]})
        channels["foci"] = synth.add_noise(fclean, rng)
        spots = {"spots": pd.DataFrame(rows)}
    nuclei = pd.DataFrame([{"nucleus_id": 1, "germ_cell": "oocyte", "in_germline": True, "centroid_z_um": center[0],
                            "centroid_y_um": center[1], "centroid_x_um": center[2],
                            "volume_um3": float((labels == 1).sum() * np.prod(synth.SPACING))}])
    return ce, labels, nuclei, channels, spots, shift


def _run(scene, p=P):
    ce, labels, nuclei, channels, spots, _ = scene
    roles = {k: v for k, v in channels.items() if k != "central_element"}
    return straighten_image(ce, labels, synth.SPACING, nuclei, params=p, roles=roles, spot_tables=spots, image_id="t")


def test_a_blurred_point_is_a_punctum_left_out_of_the_tier():
    res = _run(_scene(puncta=[{"center": [0.0, 0.3, -0.2], "peak": 900.0}], seed=1))
    st, nu = res.strands, res.nuclei.iloc[0]
    assert len(st) == 1 and st["strand_class"].iloc[0] == "punctum"
    assert "punctum" in st["qc_flags"].iloc[0] and np.isfinite(st["punctum_sigma_z_um"].iloc[0])
    assert nu["n_puncta"] == 1 and nu["n_strands"] == 0 and nu["total_length_um"] == 0.0 and nu["status"] == "no_sc"
    assert nu["piece_lw_mean_um"] == 0.0 and 0.0 < nu["puncta_signal_frac"] <= 1.0
    # punctum.exclude off (the default of entry 18): the class and flag stay, the strand counts as in v1
    keep = _run(_scene(puncta=[{"center": [0.0, 0.3, -0.2], "peak": 900.0}], seed=1), params.resolve())
    assert keep.strands["strand_class"].iloc[0] == "punctum"
    assert keep.nuclei["n_strands"].iloc[0] == 1 and keep.nuclei["total_length_um"].iloc[0] > 0.5
    # punctum.enabled off: v1 (no class, no flag)
    v1 = _run(_scene(puncta=[{"center": [0.0, 0.3, -0.2], "peak": 900.0}], seed=1), params.resolve(None, params.v1_sets()))
    assert "punctum" not in set(v1.strands["strand_class"]) and v1.nuclei["n_puncta"].iloc[0] == 0


def test_a_ribbon_beside_a_punctum_keeps_its_own_totals():
    arc = synth.arc(4.0, 3.0, 0.0, 0.5)
    sc = _scene(ribbons=[{"points": arc, "peak": 450.0}], puncta=[{"center": [0.0, 1.6, 1.6], "peak": 900.0}], seed=2)
    res = _run(sc)
    st, nu = res.strands, res.nuclei.iloc[0]
    rib = st[st["strand_class"] != "punctum"]
    assert (st["strand_class"] == "punctum").sum() == 1 and len(rib) == 1
    assert nu["total_length_um"] == pytest.approx(float(rib["length_um"].sum()))
    assert nu["n_strands"] == 1 and nu["n_puncta"] == 1
    L = float(rib["length_um"].iloc[0])
    assert nu["piece_lw_mean_um"] == pytest.approx(L, rel=1e-6)          # one piece without gaps
    assert nu["on_strand_frac"] + nu["puncta_signal_frac"] <= 1.0 + 1e-9
    assert 0.0 <= nu["steep_length_frac"] <= 1.0


@pytest.mark.parametrize("L,tilt", [(2.0, 75.0), (1.2, 80.0)])
def test_steep_ribbons_are_not_puncta(L, tilt):
    res = _run(_scene(ribbons=[{"points": synth.arc(L, 20.0, tilt, 0.4), "peak": 600.0}], seed=3))
    assert "punctum" not in set(res.strands["strand_class"])


def test_under_exclusion_a_punctum_is_never_a_foci_candidate():
    c = [0.0, 0.3, -0.2]
    res = _run(_scene(puncta=[{"center": c, "peak": 900.0}], foci=[c], seed=4))
    assert (res.strands["strand_class"] == "punctum").all()
    assert res.foci["assignment"].iloc[0] == "no_strand" and pd.isna(res.foci["strand_id"].iloc[0])


def test_puncta_families_are_development_families():
    a = synth_dev.make_scene("dense_puncta", 2002)
    b = synth_dev.make_scene("dense_puncta_clean", 2002)
    assert len(a.info["puncta"]) == 6 and b.info["puncta"] == [] and len(a.ribbons) == len(b.ribbons)
    assert np.array_equal(a.labels, b.labels)
    assert synth_dev.make_scene("puncta_agg", 2001).info["puncta"][0]["diam"] > 0
    for fam in ("puncta_pt", "puncta_agg", "dense_puncta", "dense_puncta_clean", "psf_check"):
        with pytest.raises(synth.SealedSeedError):
            synth_dev.make_scene(fam, 7_000_001)


def test_the_rule_before_the_end_matching(monkeypatch):
    """Entry 17: a traced punctum takes no part in the end matching; a tilted ribbon's trace is not held out."""
    from germquant.sc.straighten import trace as TR
    from germquant.sc.straighten.noise import NoiseModel
    from germquant.sc.straighten.sampling import Sampler

    cap = []
    orig = TR.follow_nucleus

    def hook(*a, **k):
        out = orig(*a, **k)
        cap.append([s.points.copy() for s in out[0]])
        return out

    monkeypatch.setattr(TR, "follow_nucleus", hook)
    for scene, want in ((_scene(puncta=[{"center": [0.0, 0.3, -0.2], "peak": 900.0}], seed=1), True),
                        (_scene(ribbons=[{"points": synth.arc(2.0, 20.0, 75.0, 0.4), "peak": 600.0}], seed=3), False)):
        cap.clear()
        ce, labels, nuclei, channels, spots, _ = scene
        _run(scene)
        raw = ce.astype(np.float32)
        noise = NoiseModel(a=synth.NOISE_INTERCEPT, b=synth.NOISE_SLOPE)
        sampler = Sampler(raw, (0, 0, 0), synth.SPACING, cval=np.nan)
        bg = float(np.percentile(raw[labels == 1], 15))
        longest = max(cap[0], key=len)
        assert TR.traced_punctum(longest, sampler, bg, noise, synth.SPACING, P) is want
    assert params.resolve()["punctum.before_linking"] is False          # entry 18: the end matching of v1
