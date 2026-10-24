"""The per-image PSF of sc_straighten v2 (protocol changelog entries 12 and 13): the flat-stretch path and the foci
path each recover sigma_z, the wavelength ratio, the source order and flags, and the v1 rule under psf.foci off."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from germquant.sc.straighten import params, qc, synth, synth_dev
from germquant.sc.straighten import psf as PSF
from germquant.sc.straighten.noise import NoiseModel
from germquant.sc.straighten.pipeline import psf_foci_inputs, straighten_image
from germquant.sc.straighten.refine import psf_from_profiles

P = params.resolve(None, ['psf.foci="on"'])          # the per-image PSF of entries 12 to 15 (withdrawn: off by default)
NOISE = NoiseModel(a=synth.NOISE_INTERCEPT, b=synth.NOISE_SLOPE)
WL = {"central_element": 511.0, "foci": 595.0}


def _points_volume(n=36, sigma=(0.43 * 595 / 511, 0.15 * 595 / 511, 0.15 * 595 / 511), seed=0, peak=300.0):
    """Isolated point foci on a flat background in one volume: (raw counts, positions um)."""
    rng = np.random.default_rng(seed)
    sp = synth.SPACING
    k = int(np.ceil(np.sqrt(n)))
    shape = (24, 14 * k + 10, 14 * k + 10)
    clean = np.full(shape, synth.OFFSET + 60.0)
    pos = []
    for i in range(n):
        c = np.array([12 * sp[0], (7 + 14 * (i // k)) * sp[1], (7 + 14 * (i % k)) * sp[2]])
        c = c + rng.uniform(-0.5, 0.5, 3) * np.asarray(sp)
        synth._render_points(clean, c[None, :], sp, sigma, weight=peak)
        pos.append(c)
    return synth.add_noise(clean, rng), np.asarray(pos)


def test_foci_path_recovers_sigma_z_with_the_wavelength_ratio():
    vol, pos = _points_volume()
    fits = PSF.foci_fits_crop(vol, (0, 0, 0), synth.SPACING, NOISE,
                              {"pos": pos, "isolated": PSF.isolated_spots(pos, 1.0)}, P)
    assert sum(f["kept"] for f in fits) >= 30
    est = PSF.foci_estimate(fits, WL, P)
    assert est["ok"] and est["wavelength_ratio"] == pytest.approx(511 / 595)
    assert est["sigma_z_um"] == pytest.approx(0.43, rel=0.05)
    assert est["se_z_um"] > 0
    # without both wavelengths there is no foci estimate
    assert not PSF.foci_estimate(fits, {"central_element": 511.0}, P)["ok"]
    assert PSF.foci_estimate(fits, None, P)["why"] == "no_wavelengths"
    # too few kept fits
    assert PSF.foci_estimate(fits[:5], WL, P)["why"] == "too_few_kept"


def test_planar_background_absorbs_an_edge():
    from scipy.special import erf

    # a focus 0.4 um inside the edge of a nucleoplasm of 60 counts blurred like the phantoms' (sigma 0.43 um):
    # the planar background keeps it compact and on target
    sp = synth.SPACING
    rng = np.random.default_rng(3)
    shape = (24, 24, 24)
    zz = np.arange(shape[0])[:, None, None] * sp[0] * np.ones(shape)
    sig = (0.43 * 595 / 511, 0.15 * 595 / 511, 0.15 * 595 / 511)
    out = []
    for _ in range(20):
        clean = synth.OFFSET + 30.0 * (1.0 - erf((zz - 12 * sp[0] - 0.4) / (0.43 * np.sqrt(2.0))))
        c = np.array([12 * sp[0], 12 * sp[1], 12 * sp[2]]) + rng.uniform(-0.5, 0.5, 3) * np.asarray(sp)
        synth._render_points(clean, c[None, :], sp, sig, weight=250.0)
        vol = synth.add_noise(clean, rng)
        f = PSF.foci_fits_crop(vol, (0, 0, 0), sp, NOISE, {"pos": c[None, :], "isolated": np.ones(1, bool)}, P)[0]
        out.append(f)
    kept = [f for f in out if f["kept"]]
    assert len(kept) >= 18
    # within 10 percent so close to the edge (about 8 percent low at 0.4 um, unbiased without the edge)
    assert np.median([f["sz"] for f in kept]) == pytest.approx(sig[0], rel=0.10)


def test_window_outside_the_crop_and_crowded_spots_are_rejected():
    vol, pos = _points_volume(n=4)
    edge = pos.copy()
    edge[0, 0] = 0.1                                       # 0.1 um from the first plane: the window leaves the crop
    iso = PSF.isolated_spots(np.vstack([edge, edge[1] + [0.0, 0.0, 0.5]]), 1.0)[:4]
    fits = PSF.foci_fits_crop(vol, (0, 0, 0), synth.SPACING, NOISE, {"pos": edge, "isolated": iso}, P)
    assert fits[0]["reason"] == "window_outside_crop"
    assert fits[1]["reason"] == "not_isolated"
    assert PSF.isolated_spots(np.array([[0, 0, 0], [0, 0, 0.9], [5, 5, 5.0]]), 1.0).tolist() == [False, False, True]


def _stretch(sig_ax, sig_lat, n=40, amp=200.0):
    u = np.round(np.arange(-12, 13) * 0.05, 6)
    v = np.round(np.arange(-20, 21) * 0.05, 6)
    return {"lat": (amp * np.exp(-0.5 * (u / sig_lat) ** 2) + 5.0).tolist(),
            "ax": (amp * np.exp(-0.5 * (v / sig_ax) ** 2) + 5.0).tolist(), "n": n, "u": u.tolist(), "v": v.tolist()}


def test_stretch_path_and_its_bootstrap():
    rng = np.random.default_rng(5)
    st = []
    for _ in range(12):
        s = _stretch(0.43, 0.15)
        s["ax"] = (np.asarray(s["ax"]) + rng.normal(0, 8.0, len(s["ax"]))).tolist()
        s["lat"] = (np.asarray(s["lat"]) + rng.normal(0, 8.0, len(s["lat"]))).tolist()
        st.append(s)
    e = PSF.stretch_estimate(st, P, n_min=5)
    assert e["ok"] and e["sigma_z_um"] == pytest.approx(0.43, rel=0.05) and 0 < e["se_z_um"] < 0.05
    assert not PSF.stretch_estimate(st[:4], P, n_min=5)["ok"]
    # the vectorised grid fit used by the bootstrap agrees with refine.gauss_fit
    from germquant.sc.straighten.refine import gauss_fit

    v = np.asarray(st[0]["v"])
    assert PSF._gauss_grid_fit(v, np.asarray(st[0]["ax"])[None, :])[0] == pytest.approx(
        gauss_fit(v, np.asarray(st[0]["ax"]))[0], abs=1e-9)


def test_source_order_flags_and_record():
    good_st = [_stretch(0.43, 0.15) for _ in range(30)]
    vol, pos = _points_volume(n=36, sigma=(0.45 * 595 / 511, 0.15 * 595 / 511, 0.15 * 595 / 511))
    fits = PSF.foci_fits_crop(vol, (0, 0, 0), synth.SPACING, NOISE,
                              {"pos": pos, "isolated": np.ones(len(pos), bool)}, P)
    Pn = params.resolve(None, ['psf.foci="on"', "psf.min_stretches=20", "psf.foci_min=25"])
    # the foci estimate first (protocol changelog entry 14), the stretch estimate recorded beside it
    (sz, sxy), src, d = PSF.image_psf(good_st, fits, WL, Pn)
    assert src == "foci" and sz == pytest.approx(d["foci"]["sigma_z_um"]) and d["stretches"]["ok"]
    assert sxy == pytest.approx(0.15) and d["sigma_xy_from"] == "fallback"
    assert d["foci_minus_stretches_z_um"] == pytest.approx(d["foci"]["sigma_z_um"] - d["stretches"]["sigma_z_um"])
    assert PSF.image_psf(good_st, fits[:3], WL, Pn)[1] == "stretches"
    assert PSF.image_psf(good_st[:3], fits, WL, Pn)[1] == "foci"
    assert PSF.image_psf(good_st, fits, None, Pn)[1] == "stretches"
    assert PSF.image_psf(good_st[:3], fits, None, Pn)[1] == "fallback"
    (fz, fxy), src, _ = PSF.image_psf([], [], WL, Pn)
    assert src == "fallback" and (fz, fxy) == (0.365, 0.15)
    for source, want in (("foci", "psf_from_foci"), ("fallback", "psf_fallback")):
        assert want in qc.image_flags({"psf_source": source}, P)
    assert "psf_from_foci" not in qc.image_flags({"psf_source": "stretches"}, P)


def test_foci_off_gives_the_v1_rule():
    st = [_stretch(0.40, 0.15) for _ in range(25)]
    v1 = params.resolve(None, ["psf.foci=\"off\"", "psf.min_stretches=20"])
    assert PSF.image_psf(st, [], WL, v1) == psf_from_profiles(st, v1)
    assert PSF.image_psf(st[:19], [], WL, v1) == psf_from_profiles(st[:19], v1)
    assert PSF.image_psf(st, [], WL, v1)[1] == "estimated"
    # and phase A gets no foci inputs at all
    spots = pd.DataFrame({"spot_id": [0], "nucleus_id": [1], "z_um": [1.0], "y_um": [1.0], "x_um": [1.0]})
    assert psf_foci_inputs({"spots": spots}, {"foci": np.zeros((2, 2, 2))}, WL, [1], v1) is None
    assert psf_foci_inputs({"spots": spots}, {"foci": np.zeros((2, 2, 2))}, None, [1], P) is None
    assert psf_foci_inputs({"spots": spots}, {"foci": np.zeros((2, 2, 2))}, WL, [1], P)[1]["isolated"].tolist() == [True]
    # the default (entry 16: withdrawn) is the v1 rule
    assert psf_foci_inputs({"spots": spots}, {"foci": np.zeros((2, 2, 2))}, WL, [1], params.resolve()) is None


def test_psf_check_scene_and_the_pipeline_record():
    sc = synth_dev.make_scene("psf_check", 2002)
    assert sc.info["sigma_z_level"] == 0.43 and len(sc.ribbons) == 6 * 6 - 3 + 3 or len(sc.ribbons) >= 30
    assert set(sc.nuclei["kind"]) == {"dense", "arc"} and len(sc.spots) > 20
    res = straighten_image(sc.channels["central_element"], sc.labels, sc.spacing, sc.nuclei,
                           params=params.resolve(None, ['psf.foci="on"', "psf.foci_min=10"]),
                           roles={"foci": sc.channels["foci"]},
                           spot_tables={"spots": sc.spots}, image_id="t", wavelengths=sc.info["wavelengths"])
    rec = res.record["psf"]
    assert rec["source"] in ("foci", "both") and rec["sigma_z_um"] == pytest.approx(0.43, rel=0.08)
    assert rec["foci"]["n_kept"] >= 10 and rec["wavelengths_nm"] == WL
    assert "psf_from_foci" in res.flags and res.psf_items["foci"]
    with pytest.raises(synth.SealedSeedError):
        synth_dev.make_scene("psf_check", 7_000_000)
