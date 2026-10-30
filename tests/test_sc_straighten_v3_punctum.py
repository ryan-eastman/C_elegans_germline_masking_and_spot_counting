"""The punctum class B1 with the image's measured axial sigma (protocol changelog entry 32), off by default.

A traced steep strand (frac_steep60 above 0.5) at most punctum.max_xy_um (0.5 um) long in the image plane whose
core profile along s is one Gaussian (centre on the strand, residual RMS at most 0.2 of the amplitude) no wider in z
than punctum.sigma_ratio (1.3) times the reference axial sigma is of class punctum and leaves the ribbon tier. The
reference is the axial sigma of the image's compact isolated RAD-51 foci (the foci estimate of entry 12, scaled to
the SYP wavelength) with a depth term, or the image PSF's without foci. Real steep ribbons longer than the axial blur
stay ribbons at every PSF width; puncta leave the totals while ribbons stay.
"""
from __future__ import annotations

import numpy as np
import pytest
import sc_v3_phantoms as Ph

from germquant.sc.straighten import params, synth, synth_dev
from germquant.sc.straighten import psf as PSF
from germquant.sc.straighten.refine import along_profile_fit, along_profile_sigma, one_gaussian

WL = Ph.WL
P_B1 = params.resolve(None, Ph.sets(params.V3_PUNCTUM))
P_B1_V2 = params.resolve(None, ["punctum.enabled=true", "punctum.exclude=true"])   # entry 16: the image PSF
WIDE = 0.55                                                 # an axial PSF wider than the fallback 0.365 um


# ------------------------------------------------------------------------------------ the parts
def test_the_one_gaussian_profile_test():
    s = np.arange(-0.5, 1.5 + 1e-9, 0.05)
    rng = np.random.default_rng(3)
    one = 20.0 + 300.0 * np.exp(-0.5 * ((s - 0.5) / 0.4) ** 2) + rng.normal(0, 5.0, len(s))
    f = along_profile_fit(s, one)
    assert f["sigma"] == along_profile_sigma(s, one)          # the sigma of v2, bit for bit
    assert f["sigma"] == pytest.approx(0.4, abs=0.05) and one_gaussian(f, 1.0, 0.2)
    two = 20.0 + 300.0 * (np.exp(-0.5 * ((s + 0.1) / 0.25) ** 2) + np.exp(-0.5 * ((s - 1.1) / 0.25) ** 2))
    assert not one_gaussian(along_profile_fit(s, two), 1.0, 0.2)
    ramp = 20.0 + 300.0 * np.clip(s + 0.5, 0, None) / 2.0     # brightest beyond the end: no peak on the strand
    assert not one_gaussian(along_profile_fit(s, ramp), 1.0, 0.2)
    assert not one_gaussian(along_profile_fit(s[:5], one[:5]), 1.0, 0.2)   # too few samples: no fit


def _fits(cz, sz):
    return [{"kept": True, "sz": float(b), "cz": float(a), "sy": 0.15, "sx": 0.15} for a, b in zip(cz, sz)]


def test_the_measured_axial_sigma_with_its_depth_term():
    cz = np.linspace(1.0, 11.0, 60)
    sz = 0.40 + 0.01 * (cz - 6.0)                            # 0.35 to 0.45 um over 10 um of depth (foci channel)
    ratio = WL["central_element"] / WL["foci"]
    ref = PSF.punctum_reference(_fits(cz, sz), WL, P_B1, (0.365, 0.15))
    assert ref["source"] == "foci_depth"
    assert ref["slope_per_um"] == pytest.approx(0.01 * ratio, rel=1e-6)
    assert PSF.punctum_sigma_at(ref, 6.0) == pytest.approx(float(np.median(sz)) * ratio, rel=1e-6)
    assert PSF.punctum_sigma_at(ref, 9.0) == pytest.approx((0.40 + 0.03) * ratio, rel=0.02)
    hi = PSF.punctum_sigma_at(ref, 40.0)                    # beyond the foci: clipped to their range
    assert hi <= float(np.percentile(sz, 90)) * ratio + 1e-9
    p_const = params.resolve(None, Ph.sets({**params.V3_PUNCTUM, "punctum.sigma_source": "foci"}))
    ref_c = PSF.punctum_reference(_fits(cz, sz), WL, p_const, (0.365, 0.15))
    assert ref_c["source"] == "foci" and PSF.punctum_sigma_at(ref_c, 9.0) == pytest.approx(0.40 * ratio, rel=1e-6)
    # too few fits for a depth term: the constant estimate; too few for any estimate, or no wavelengths: the PSF
    assert PSF.punctum_reference(_fits(cz[:20], sz[:20]), WL, P_B1, (0.365, 0.15))["source"] == "foci"
    narrow = PSF.punctum_reference(_fits(np.full(60, 5.0) + np.linspace(0, 1, 60), sz), WL, P_B1, (0.365, 0.15))
    assert narrow["source"] == "foci"                        # the fits span less than psf.ref_depth_min_span_um
    assert PSF.punctum_reference(_fits(cz[:10], sz[:10]), WL, P_B1, (0.365, 0.15))["source"] == "psf"
    assert PSF.punctum_reference(_fits(cz, sz), None, P_B1, (0.365, 0.15))["source"] == "psf"
    v2 = PSF.punctum_reference(_fits(cz, sz), WL, P_B1_V2, (0.365, 0.15))
    assert v2["source"] == "psf" and PSF.punctum_sigma_at(v2, 9.0) == 0.365


# ------------------------------------------------------------------------------------ phantoms
def test_a_point_in_a_wide_psf_is_a_punctum_with_the_measured_sigma():
    sc = Ph.scene(puncta=[{"center": [0.0, 0.4, -0.3], "peak": 900.0}], sigma_z=WIDE, foci=24, seed=21, R=3.0)
    v3 = Ph.run(sc, P_B1, WL)
    ref = v3.record["punctum_reference"]
    assert ref["source"] in ("foci", "foci_depth")
    assert ref["sigma_z_um"] == pytest.approx(WIDE, rel=0.12)
    assert "punctum" in set(v3.strands["strand_class"])
    nu = v3.nuclei.iloc[0]
    assert nu["n_puncta"] >= 1 and nu["total_length_um"] == pytest.approx(
        float(v3.strands.loc[v3.strands["strand_class"] != "punctum", "length_um"].sum()))
    # the rule of v2 (the image PSF, here the fallback 0.365 um) does not see it: too wide for 1.3 x 0.365
    v2 = Ph.run(sc, P_B1_V2, WL)
    assert "punctum" not in set(v2.strands["strand_class"])
    # the image PSF itself stays that of v1 (psf.foci off): only the punctum reference is measured
    assert v3.record["psf"]["sigma_z_um"] == v2.record["psf"]["sigma_z_um"]


@pytest.mark.parametrize("L,tilt", [(2.6, 75.0), (3.5, 85.0)])
def test_steep_ribbons_longer_than_the_blur_stay_ribbons(L, tilt):
    rib = synth._densify(synth.arc(L, 20.0, tilt, 0.4))
    sc = Ph.scene(ribbons=[{"points": rib, "peak": 600.0}], sigma_z=WIDE, foci=24, seed=22, R=3.0)
    v3 = Ph.run(sc, P_B1, WL)
    assert "punctum" not in set(v3.strands["strand_class"])
    assert v3.nuclei["n_puncta"].iloc[0] == 0


@pytest.mark.parametrize("sigma_z", [0.33, 0.48])
@pytest.mark.parametrize("L,tilt", [(2.0, 80.0), (3.0, 70.0), (3.0, 88.0)])
def test_steep_ribbons_survive_at_every_psf_width(sigma_z, L, tilt):
    rib = Ph.line(-L / 2, L / 2, tilt)
    sc = Ph.scene([{"points": rib, "peak": 600.0}], sigma_z=sigma_z, foci=28, seed=int(100 * sigma_z + L + tilt),
                  R=3.4)
    v3 = Ph.run(sc, P_B1, WL)
    assert v3.record["punctum_reference"]["source"] in ("foci", "foci_depth")
    assert "punctum" not in set(v3.strands["strand_class"])
    assert v3.nuclei["n_puncta"].iloc[0] == 0


def test_puncta_leave_the_tier_and_ribbons_stay():
    """A nucleus with two in-plane ribbons and three point puncta under a wide PSF (0.48 um) with foci: the
    puncta's strands are of class punctum and out of the totals, the ribbons' are not."""
    a = Ph.line(-2.6, 2.6, 8.0, origin=(0.0, -1.0, 0.0))
    b = Ph.line(-2.2, 2.2, 25.0, origin=(0.3, 1.1, 0.2))
    puncta = [{"center": c, "peak": p} for c, p in (([1.6, 0.0, -1.6], 820.0), ([-1.6, 0.2, 1.4], 760.0),
                                                     ([1.4, 0.1, 1.8], 880.0))]
    sc = Ph.scene([{"points": a, "peak": 520.0}, {"points": b, "peak": 480.0}], puncta, sigma_z=0.48, foci=24,
                  seed=111, R=3.4)
    v3 = Ph.run(sc, P_B1, WL)
    st = v3.strands
    nu = v3.nuclei.iloc[0]
    assert nu["n_puncta"] >= 1
    assert nu["total_length_um"] == pytest.approx(float(st.loc[st["strand_class"] != "punctum", "length_um"].sum()))
    long_ = st[st["length_um"] >= 3.0]
    assert len(long_) >= 2 and "punctum" not in set(long_["strand_class"])


def test_without_foci_the_reference_is_the_image_psf():
    sc = Ph.scene(puncta=[{"center": [0.0, 0.3, -0.2], "peak": 900.0}], seed=1)
    res = Ph.run(sc, P_B1)
    assert res.record["punctum_reference"]["source"] == "psf"
    assert (res.strands["strand_class"] == "punctum").all() and len(res.strands) == 1


# ------------------------------------------------------------------------------------ the development family
def test_the_puncta_psf_family():
    a = synth_dev.make_scene("puncta_psf", 2001)
    b = synth_dev.make_scene("puncta_psf", 2001)
    assert np.array_equal(a.channels["central_element"], b.channels["central_element"])
    assert set(a.channels) == {"central_element", "foci"} and len(a.spots) > 20
    assert a.info["wavelengths"] == WL and len(a.info["puncta"]) > 0
    kinds = {r.kind for r in a.ribbons}
    assert kinds == {"dense", "steep_arc"}
    lv = a.info["sigma_z_level"]
    assert lv in synth_dev.PSF_LEVELS and synth_dev.puncta_psf_sigma_z(a.info, a.info["z_mid_um"]) == pytest.approx(lv)
    assert abs(a.info["sigma_z_gradient_per_um"]) <= synth_dev.PUNCTA_PSF_GRADIENT
    # the two layers sit at different depths, so the PSF differs between them by the gradient
    zs = sorted({round(float(p["center"][0]), 0) for p in a.info["puncta"]})
    assert zs[-1] - zs[0] > 3.0
    with pytest.raises(synth.SealedSeedError):
        synth_dev.make_scene("puncta_psf", 7_000_001)
    with pytest.raises(synth.SealedSeedError):
        synth_dev.make_scene("puncta_psf", 2001, final=True)
