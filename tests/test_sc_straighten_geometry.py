"""sc_straighten geometry: P-spline length, frames, the strip round trip, the end rule, the noise transform,
the crop frame and the end matching (docs/SC_STRAIGHTEN.md; design section 11)."""
from __future__ import annotations

import numpy as np
import pytest

from germquant.sc.straighten import frames as F
from germquant.sc.straighten import link, noise, params, refine, synth
from germquant.sc.straighten.follow import RawStrand
from germquant.sc.straighten.sampling import Sampler

P = params.resolve()


def _arc_points(r=1.2, L=5.0, n=400):
    th = np.linspace(0, L / r, n)
    return np.stack([np.zeros(n), r * np.sin(th), r * (1 - np.cos(th))], 1)


def _helix(r=0.8, pitch=1.5, turns=1.5, n=600):
    th = np.linspace(0, 2 * np.pi * turns, n)
    return np.stack([pitch * th / (2 * np.pi), r * np.cos(th), r * np.sin(th)], 1)


def _line_from_points(Pts, sig=(1e-3, 1e-3, 1e-3), knot=0.3, lam=1e-6):
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(Pts, axis=0), axis=1))])
    W = np.tile(1.0 / np.asarray(sig) ** 2, (len(Pts), 1))
    spl = refine.pspline(s, Pts, W, knot, lam * float(np.mean(W[0])))
    return refine.Centreline(spl, s[0], s[-1])


def test_pspline_arc_length_noise_free_and_noisy():
    Pts = _arc_points()
    line = _line_from_points(Pts)
    assert abs(line.L - 5.0) < 0.005
    rng = np.random.default_rng(3)
    errs = []
    for _ in range(20):
        noisy = Pts[::4] + rng.normal(0, 0.02, (len(Pts[::4]), 3))
        s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(noisy, axis=0), axis=1))])
        W = np.tile(1.0 / np.array([0.02, 0.02, 0.02]) ** 2, (len(noisy), 1))
        spl = refine.pspline(s, noisy, W, 0.3, 1.0 * float(np.mean(W[0])))
        errs.append(refine.Centreline(spl, s[0], s[-1]).L - 5.0)
    assert abs(float(np.mean(errs))) < 0.05


def test_helix_arc_length_under_anisotropic_sampling():
    Pts = _helix()
    true = float(np.linalg.norm(np.diff(_helix(n=20000), axis=0), axis=1).sum())
    # sampled unevenly (denser where it is steep in z), the spline arc length still converges
    idx = np.unique(np.r_[np.arange(0, 600, 3), np.arange(0, 600, 7)])
    line = _line_from_points(Pts[idx])
    assert abs(line.L - true) / true < 0.005


def test_rmf_frames_orthonormal_and_twist_free_on_a_helix():
    line = _line_from_points(_helix())
    fr = F.StrandFrame(line, ds=0.02, extend=0.0)
    T, nl, na = fr.T, fr.n_lat, fr.n_ax
    for a, b in ((T, nl), (T, na), (nl, na)):
        assert np.max(np.abs((a * b).sum(1))) < 1e-6
    for a in (T, nl, na):
        assert np.max(np.abs(np.linalg.norm(a, axis=1) - 1)) < 1e-9
    # rotation minimising: the normal turns only toward the tangent (no rotation about t)
    dn = np.diff(fr.n1, axis=0) / 0.02
    n2 = np.cross(fr.T, fr.n1)
    twist = np.abs((dn * 0.5 * (n2[1:] + n2[:-1])).sum(1))
    assert np.median(twist) < 0.02


def test_flat_curve_gives_axial_normal_along_z():
    line = _line_from_points(_arc_points())
    fr = F.StrandFrame(line)
    assert np.allclose(np.abs(fr.n_ax[:, 0]), 1.0, atol=1e-6)
    assert np.allclose(fr.n_lat[:, 0], 0.0, atol=1e-6)
    # and the roll from the optical axis is zero there
    assert np.nanmax(np.abs(fr.roll_from_axial_deg())) < 1e-3


def test_image_to_strip_to_image_round_trip():
    line = _line_from_points(_helix(turns=0.8))
    fr = F.StrandFrame(line)
    rng = np.random.default_rng(0)
    s = rng.uniform(0.2, line.L - 0.2, 200)
    u = rng.uniform(-0.5, 0.5, 200)
    v = rng.uniform(-0.5, 0.5, 200)
    pts = fr.from_local(s, u, v)
    s2, u2, v2, w2 = fr.to_local(pts)
    back = fr.from_local(s2, u2, v2)
    assert np.max(np.linalg.norm(back - pts, axis=1)) < 1e-6
    assert np.max(np.abs(w2)) < 1e-6
    assert np.max(np.abs(s2 - s)) < 1e-6
    # the multi-strand map picks the right strand
    other = _line_from_points(_arc_points() + np.array([5.0, 5.0, 5.0]))
    res = F.to_strip([F.StrandFrame(other), fr], pts[:20])
    assert (res["strand"] == 1).all()
    assert np.max(np.abs(res["s"] - s[:20])) < 1e-6


def test_reversed_centreline_is_the_same_curve():
    line = _line_from_points(_helix(turns=0.6))
    rev = line.reversed()
    assert abs(rev.L - line.L) < 1e-9
    s = np.linspace(0, line.L, 37)
    assert np.max(np.abs(rev.point(line.L - s) - line.point(s))) < 1e-9


def test_erf_end_fit_bias_at_snr_10():
    rng = np.random.default_rng(1)
    sg = np.arange(-0.6, 1.4001, 0.05)
    errs = []
    for _ in range(300):
        e = rng.uniform(0.2, 0.9)
        sig = 0.16
        a = 10.0
        y = 2.0 + a * (0.5 * (1 + __import__("scipy").special.erf((e - sg) / (sig * np.sqrt(2))))) \
            + rng.normal(0, 1.0, len(sg))
        r = refine.fit_erf(sg, y, 0.6 * 0.15, 1.6 * 0.15)
        errs.append(r[0] - e)
    assert abs(float(np.mean(errs))) < 0.02


def test_gat_gives_unit_noise():
    rng = np.random.default_rng(2)
    planes, masks = [], []
    base = np.linspace(110, 900, 300)[None, :].repeat(400, 0)
    for _ in range(6):
        var = np.maximum(4.04 * base - 394, 25)
        planes.append(base + np.sqrt(var) * rng.standard_normal(base.shape))
        masks.append(np.ones(base.shape, bool))
    m = noise.fit_noise_model(planes, masks)
    assert abs(m.b - 4.04) / 4.04 < 0.1
    y = noise.gat(np.stack(planes), m)
    d = y[:, :, 2:] - 2 * y[:, :, 1:-1] + y[:, :, :-2]
    assert abs(float(np.std(d) / np.sqrt(6)) - 1.0) < 0.05


def test_crop_origin_frame_mapping():
    rng = np.random.default_rng(4)
    vol = rng.normal(size=(12, 20, 22)).astype(np.float32)
    sp = (0.2, 0.108333, 0.108333)
    full = Sampler(vol, (0, 0, 0), sp)
    org = (3, 5, 7)
    crop = Sampler(vol[3:, 5:, 7:], org, sp)
    pts = (np.array([[5.3, 9.1, 11.7], [6.0, 12.2, 15.5]]) * np.asarray(sp))
    assert np.allclose(full(pts), crop(pts), atol=1e-6)


def test_end_matching_is_invariant_to_input_order():
    # three collinear pieces with gaps and one unrelated strand
    xs = [np.linspace(0, 1.5, 16), np.linspace(1.9, 3.4, 16), np.linspace(3.8, 5.0, 13)]
    pieces = [RawStrand(np.stack([np.zeros_like(x), np.zeros_like(x), x], 1), ["free", "free"]) for x in xs]
    other = RawStrand(np.stack([np.zeros(10), np.linspace(2, 3, 10), np.full(10, 9.0)], 1), ["free", "free"])
    a = link.link_strands(pieces + [other], P)
    b = link.link_strands([other, pieces[2], pieces[0], pieces[1]], P)

    def key(strands):
        return sorted((round(st.length_um, 6), len(st.links)) for st in strands)

    assert key(a) == key(b)
    assert max(len(st.links) for st in a) == 2           # the three pieces form one chain of two links


def test_synthetic_ribbons_obey_the_curvature_bound():
    sc = synth.dense_scene(11)
    for rib in sc.ribbons:
        k = synth.curvature(rib.points)
        assert np.max(k) <= 1.0 + 1e-3
    with pytest.raises(synth.SealedSeedError):
        synth.dense_scene(7_200_000)


def test_axial_ghosts_and_shadowed_strands_are_dropped():
    from germquant.sc.straighten.account import drop_axial_ghosts, drop_shadowed

    def seg(a, b, n=60):
        a, b = np.asarray(a, float), np.asarray(b, float)
        return RawStrand(points=a + np.linspace(0, 1, n)[:, None] * (b - a), end_types=["free", "free"])

    main = seg((2, 2, 2), (2, 2, 6))
    ghost = seg((2.05, 2, 6.0), (3.0, 2, 6.05))            # steep, standing on the main strand's end
    lone = seg((2, 8, 8), (3.0, 8, 8))                      # steep but alone
    def has(lst, x):
        return any(k is x for k in lst)

    kept, n = drop_axial_ghosts([main, ghost, lone], P)
    assert n == 1 and not has(kept, ghost) and has(kept, lone) and has(kept, main)
    near = seg((2, 2.2, 2.2), (2, 2.2, 5.8))                # 0.2 um beside the main strand
    far = seg((2, 2.6, 2.2), (2, 2.6, 5.8))                 # 0.6 um beside: resolvable
    kept, n = drop_shadowed([main, near, far], P)
    assert n == 1 and not has(kept, near) and has(kept, far)
