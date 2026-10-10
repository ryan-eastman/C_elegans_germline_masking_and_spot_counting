"""The EXPERIMENTAL net3d foci detector (germquant.spots.net3d), on the CPU with tiny synthetic data:

- the network's shapes, size divisor and receptive radius (a perturbation never reaches farther);
- tiled inference equals untiled inference (xy and z tiling, the occupancy skip);
- peak finding on a synthetic heatmap: ellipsoidal suppression, sub-voxel refinement, the threshold
  property the validation sweep relies on, nucleus assignment with the nearest-label fallback and the
  contested flag;
- the loss region and weights (exact dilations, Class A targets, Class B ignored, the masked loss);
- the backend contract with a randomly initialised model (tables, the DNA channel, the provenance, the
  pipeline's channel passing and sub-hash, the default path untouched);
- ``germquant train-spots`` for one epoch on a tiny synthetic bench cache, never opening the test
  gonad, and the bench harness scoring the result.
Nothing here needs microscopy data or a GPU.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

torch = pytest.importorskip("torch")

from germquant.spots.net3d import data as nd  # noqa: E402
from germquant.spots.net3d import peaks as pk  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SP = (0.2, 0.108, 0.108)


def _model(levels=3, base=4, in_channels=2, seed=0):
    from germquant.spots.net3d.model import build_model

    torch.manual_seed(seed)
    m = build_model({"in_channels": in_channels, "base_channels": base, "levels": levels})
    # non-trivial batch-norm statistics, so the eval-mode affine maps are not identities
    g = torch.Generator().manual_seed(seed + 1)
    for mod in m.modules():
        if isinstance(mod, torch.nn.BatchNorm3d):
            mod.running_mean.copy_(torch.randn(mod.num_features, generator=g) * 0.1)
            mod.running_var.copy_(torch.rand(mod.num_features, generator=g) + 0.5)
    return m.eval()


# ------------------------------------------------------------------------------------ model
def test_model_shapes_divisor_and_radius():
    from germquant.spots.net3d.model import (max_pool, pools_for, receptive_radius, size_divisor,
                                             tile_margin)

    assert pools_for(4) == [(1, 2, 2), (2, 2, 2), (2, 2, 2)]
    assert size_divisor(4) == (4, 8, 8) and size_divisor(3) == (2, 4, 4)
    assert receptive_radius(4) == (30, 58, 58) and receptive_radius(3) == (14, 26, 26)
    assert tile_margin(4) == (32, 64, 64) and tile_margin(3) == (14, 28, 28)
    with pytest.raises(ValueError):
        pools_for(1)
    m = _model(levels=4, base=4)
    out = m(torch.zeros(2, 2, 8, 32, 24))
    assert tuple(out.shape) == (2, 1, 8, 32, 24)
    assert float(m.head.bias.detach()) == -4.0
    with pytest.raises(ValueError, match="multiples"):
        max_pool(torch.zeros(1, 1, 3, 4, 4), (2, 2, 2))
    x = torch.arange(16.0).reshape(1, 1, 1, 4, 4)
    assert max_pool(x, (1, 2, 2)).flatten().tolist() == [5.0, 7.0, 13.0, 15.0]


@pytest.mark.parametrize("levels,shape", [(3, (32, 64, 64)), (4, (64, 128, 128))])
def test_a_perturbation_stays_inside_the_receptive_radius(levels, shape):
    from germquant.spots.net3d.model import receptive_radius

    m = _model(levels=levels, base=4, in_channels=1)
    x = torch.zeros((1, 1) + shape)
    c = tuple(s // 2 for s in shape)
    x2 = x.clone()
    x2[(0, 0) + c] = 10.0
    with torch.no_grad():
        d = (m(x2) - m(x)).abs()[0, 0].numpy()
    hit = np.argwhere(d > 0)
    assert len(hit) > 0
    reach = np.abs(hit - np.array(c)).max(axis=0)
    assert np.all(reach <= np.array(receptive_radius(levels))), (reach, receptive_radius(levels))


# -------------------------------------------------------------------------------- inference
def test_tiled_inference_equals_untiled_in_xy_and_z():
    from germquant.spots.net3d.infer import channel_ranges, predict_heatmap, tile_plan

    rng = np.random.default_rng(3)
    m = _model(levels=3, base=4)
    foci = rng.integers(90, 400, size=(20, 100, 90)).astype(np.uint16)
    dna = rng.integers(90, 900, size=(20, 100, 90)).astype(np.uint16)
    ranges = channel_ranges([foci, dna])
    whole = predict_heatmap(m, [foci, dna], ranges, tile_zyx=None)
    assert whole.shape == foci.shape and whole.dtype == np.float32
    assert len(tile_plan(foci.shape, 3, (96, 64, 64))) > 4
    tiled = predict_heatmap(m, [foci, dna], ranges, tile_zyx=(96, 64, 64))
    np.testing.assert_allclose(tiled, whole, rtol=0, atol=2e-6)
    # z tiling as well (40 planes, tiles of 32 with a 14-plane overlap)
    f2 = rng.integers(90, 400, size=(40, 24, 24)).astype(np.uint16)
    d2 = rng.integers(90, 900, size=(40, 24, 24)).astype(np.uint16)
    r2 = channel_ranges([f2, d2])
    np.testing.assert_allclose(predict_heatmap(m, [f2, d2], r2, tile_zyx=(32, None, None)),
                               predict_heatmap(m, [f2, d2], r2, tile_zyx=None), rtol=0, atol=2e-6)
    with pytest.raises(ValueError, match="no core"):
        predict_heatmap(m, [foci, dna], ranges, tile_zyx=(96, 56, 56))


def test_occupancy_skip_changes_nothing_near_the_nuclei():
    from germquant.spots.net3d.infer import channel_ranges, predict_heatmap

    rng = np.random.default_rng(4)
    m = _model(levels=3, base=4)
    foci = rng.integers(90, 400, size=(12, 96, 96)).astype(np.uint16)
    dna = rng.integers(90, 900, size=(12, 96, 96)).astype(np.uint16)
    labels = np.zeros(foci.shape, np.int32)
    labels[3:9, 10:20, 12:22] = 1
    ranges = channel_ranges([foci, dna])
    full = predict_heatmap(m, [foci, dna], ranges, tile_zyx=(96, 64, 64))
    skip = predict_heatmap(m, [foci, dna], ranges, tile_zyx=(96, 64, 64), labels=labels)
    near = np.zeros(foci.shape, bool)
    near[:, 6:24, 8:26] = True                         # the nucleus plus more than 0.3 um around it
    np.testing.assert_array_equal(skip[near], full[near])
    assert (skip[:, 60:, 60:] == 0).all() and (full[:, 60:, 60:] > 0).all()


def test_normalisation_is_robust_percentiles():
    from germquant.spots.net3d.infer import normalize_block, robust_range

    vol = np.arange(4 * 8 * 8, dtype=np.uint16).reshape(4, 8, 8)
    lo, hi = robust_range(vol, 1.0, 99.9, stride_yx=1)
    assert lo == pytest.approx(np.percentile(vol.astype(np.float32), 1.0))
    assert hi == pytest.approx(np.percentile(vol.astype(np.float32), 99.9))
    assert robust_range(np.full((2, 4, 4), 7, np.uint16)) == (7.0, 8.0)          # flat: finite scale
    out = normalize_block(np.array([lo, hi]), lo, hi)
    assert out.dtype == np.float32
    np.testing.assert_allclose(out, [0.0, 1.0], atol=1e-6)


# ------------------------------------------------------------------------------------ peaks
def _gauss_heat(shape, centres, sigma_vox=(1.2, 1.1, 1.1), amp=0.9):
    zz, yy, xx = np.indices(shape).astype(np.float64)
    h = np.zeros(shape)
    for c in centres:
        g = amp * np.exp(-0.5 * (((zz - c[0]) / sigma_vox[0]) ** 2 + ((yy - c[1]) / sigma_vox[1]) ** 2
                                 + ((xx - c[2]) / sigma_vox[2]) ** 2))
        h = np.maximum(h, g)
    return h.astype(np.float32)


def test_peaks_refine_suppress_and_threshold_property():
    shape = (16, 40, 40)
    centres = [(5.3, 10.2, 10.7), (5.0, 10.0, 14.6),       # 0.42 um apart in xy: both kept
               (8.0, 30.0, 30.0), (9.6, 30.0, 30.0)]       # 0.32 um apart in z, same xy: suppressed
    heat = _gauss_heat(shape, centres)
    heat[12, 30, 8] = heat[12, 30, 9] = 0.6                # a two-voxel plateau: one peak survives
    labels = np.ones(shape, np.int32)
    f = pk.find_spots(heat, labels, SP, 0.2)
    assert len(f) == 4, f
    pos = f[["zr", "yr", "xr"]].to_numpy()
    for c in centres[:2]:
        d = np.abs(pos - np.array(c)).max(axis=1).min()
        assert d < 0.05, (c, d)                            # log-parabola: exact for a Gaussian blob
    # the z neighbour (0.32 um < 0.6 um away) left no second peak at (y, x) = (30, 30)
    assert int(((np.abs(pos[:, 1] - 30) < 0.5) & (np.abs(pos[:, 2] - 30) < 0.5)).sum()) == 1
    plateau = f[(f["z"] == 12) & (f["y"] == 30)]
    assert len(plateau) == 1 and int(plateau["x"].iloc[0]) == 8          # ties: lower x first
    assert float(plateau["xr"].iloc[0]) == pytest.approx(8.5)            # refined to the plateau's middle
    # kept set at a higher threshold = kept set at a lower one, filtered
    lo = pk.find_spots(heat, labels, SP, 0.05)
    for t in (0.3, 0.55, 0.7, 0.85):
        hi = pk.find_spots(heat, labels, SP, t).sort_values(["z", "y", "x"]).reset_index(drop=True)
        cut = lo[lo["prob"] >= t].sort_values(["z", "y", "x"]).reset_index(drop=True)
        pd.testing.assert_frame_equal(hi, cut)
    assert pk.find_spots(np.zeros(shape, np.float32), labels, SP, 0.1).empty


def test_ellipsoid_offsets_and_suppression_radius():
    off = pk.ellipsoid_offsets((0.6, 0.3), SP)
    assert {tuple(o) for o in off} >= {(0, 0, 0), (3, 0, 0), (-3, 0, 0), (0, 2, 0), (0, 0, -2)}
    assert (4, 0, 0) not in {tuple(o) for o in off} and (0, 3, 0) not in {tuple(o) for o in off}
    cand = pd.DataFrame({"z": [5, 5, 5], "y": [10, 10, 10], "x": [10, 12, 14], "prob": [0.9, 0.8, 0.7]})
    kept = pk.suppress(cand, SP, (0.6, 0.3))
    assert kept["x"].tolist() == [10, 14]                 # 0.216 um is inside, 0.43 um outside


def test_assignment_nearest_label_and_contested():
    labels = np.zeros((10, 30, 30), np.int32)
    labels[3:7, 5:12, 5:12] = 1
    labels[3:7, 5:12, 15:22] = 2
    pos = np.array([[5.0, 8.0, 8.0],          # inside 1
                    [5.0, 8.0, 13.6],         # in the 3-voxel gap: 0.28 um from 1, 0.15 um from 2
                    [5.0, 8.0, 23.5],         # 0.27 um outside 2
                    [5.0, 25.0, 25.0],        # far from everything
                    [5.0, 8.0, 10.6]])        # inside 1 and 0.48 um from 2: not contested
    nid, rule, contested = pk.assign_nuclei(pos, labels, SP, assign_max_um=0.3)
    assert nid.tolist() == [1, 2, 2, 0, 1]
    assert rule.tolist() == ["label", "nearest", "nearest", "none", "label"]
    assert contested.tolist() == [False, True, False, False, False]


# ------------------------------------------------------------------------- data and loss
def _brute_dilate(mask, radius_um, sp):
    out = np.zeros_like(mask)
    pts = np.argwhere(mask)
    for idx in np.ndindex(mask.shape):
        d2 = (((pts - np.array(idx)) * np.array(sp)) ** 2).sum(axis=1)
        out[idx] = bool(len(d2)) and d2.min() <= radius_um ** 2 + 1e-9
    return out


def test_dilations_are_exact():
    rng = np.random.default_rng(5)
    mask = rng.random((7, 12, 12)) > 0.97
    np.testing.assert_array_equal(nd.dilate_ellipsoid(mask, 0.3, SP), _brute_dilate(mask, 0.3, SP))
    cyl = nd.dilate_cylinder(mask, 0.25, 1, SP)
    ref = np.zeros_like(mask)
    for k in (-1, 0, 1):
        shifted = np.zeros_like(mask)
        if k < 0:
            shifted[:k] = mask[-k:]
        elif k > 0:
            shifted[k:] = mask[:-k]
        else:
            shifted = mask.copy()
        ref |= _brute_dilate(shifted, 0.25, (1e6, SP[1], SP[2]))       # a huge dz keeps it per plane
    np.testing.assert_array_equal(cyl, ref)
    assert not nd.dilate_ellipsoid(np.zeros((3, 4, 4), bool), 0.3, SP).any()


def _toy_gonad():
    shape = (16, 48, 48)
    labels = np.zeros(shape, np.int32)
    labels[4:12, 8:24, 8:24] = 1                    # germline nucleus, traced
    labels[4:12, 28:44, 28:44] = 2                  # germline nucleus, not traced
    labels[4:12, 30:44, 2:14] = 3                   # off-gonad nucleus
    truth_labels = np.zeros(shape, np.int32)
    truth_labels[4:12, 8:24, 8:24] = 7
    region = nd.build_region(labels, [1, 2], truth_labels, SP)
    rng = np.random.default_rng(0)
    foci = rng.integers(90, 120, size=shape).astype(np.uint16)
    dna = np.where(labels > 0, 600, 100).astype(np.uint16)
    a = np.array([[8.0, 16.0, 16.0], [8.0, 36.0, 36.0], [8.0, 37.0, 8.0]])
    b = np.array([[8.0, 12.0, 20.0]])
    g = nd.TrainGonad(tag="toy", spacing=SP, channels=[foci, dna], ranges=[(90.0, 120.0), (100.0, 600.0)],
                      region=region, spots_a=a, spots_b=b, centers=np.array([[8.0, 16.0, 16.0]]))
    return g, labels


def test_region_targets_and_weights():
    g, labels = _toy_gonad()
    r = g.region
    assert r[8, 16, 16] == 2 and r[8, 36, 36] == 1 and r[8, 37, 8] == 0 and r[0, 0, 0] == 0
    assert r[8, 16, 25] == 2                      # germline dilated by 0.3 um (2 voxels in x) and traced
    assert r[8, 16, 27] == 0                      # beyond 0.3 um of the germline
    assert list(g.spot_pool) == [0, 1]            # the off-gonad Class A spot cannot centre a patch
    x, y, w = nd.make_patch(g, (0, 0, 0), g.shape, {}, rng=None, augment=False)
    assert x.shape == (2,) + g.shape and y.shape == w.shape == (1,) + g.shape
    assert y[0, 8, 16, 16] == pytest.approx(1.0) and y[0, 8, 36, 36] == pytest.approx(1.0)
    assert w[0, 8, 16, 16] == pytest.approx(10.0)                # traced: 1 x (1 + 9 x 1)
    assert w[0, 8, 36, 36] == pytest.approx(0.5 * 10.0)          # untraced germline: 0.5 x 10
    assert w[0, 8, 37, 8] == 0.0                                 # off-gonad spot: outside the region
    assert w[0, 8, 12, 20] == 0.0 and w[0, 8, 12, 22] == 0.0     # inside the Class B ball (0.3 um in xy)
    assert w[0, 8, 12, 24] > 0.0                                 # 0.43 um away: weighted again
    assert w[0, 8, 20, 12] == pytest.approx(1.0 * (1 + 9 * float(y[0, 8, 20, 12])))
    # the masked loss ignores what the weight zeroes
    from torch.nn import functional as F

    logits = torch.zeros((1,) + y.shape)
    yt, wt = torch.from_numpy(y[None]), torch.from_numpy(w[None])

    def loss(lg):
        return float((F.binary_cross_entropy_with_logits(lg, yt, reduction="none") * wt).sum() / wt.sum())

    base = loss(logits)
    junk = logits.clone()
    junk[0, 0][torch.from_numpy(w[0] == 0)] = 25.0
    assert loss(junk) == pytest.approx(base, abs=1e-6)
    junk[0, 0, 8, 16, 16] = -25.0
    assert loss(junk) > base
    # a patch hanging over the edge: zero input and zero weight outside the volume
    x2, _, w2 = nd.make_patch(g, (-4, -8, -8), (16, 32, 32), {}, rng=None, augment=False)
    assert (x2[:, :4] == 0).all() and (w2[:, :4] == 0).all()


def test_sampling_is_reproducible_and_stays_in_bounds():
    g, _ = _toy_gonad()
    a = nd.sample([g], (8, 16, 16), {}, seed=1, epoch=2, index=3)
    b = nd.sample([g], (8, 16, 16), {}, seed=1, epoch=2, index=3)
    c = nd.sample([g], (8, 16, 16), {}, seed=1, epoch=2, index=4)
    for u, v in zip(a, b):
        np.testing.assert_array_equal(u, v)
    assert not all(np.array_equal(u, v) for u, v in zip(a, c))
    rng = np.random.default_rng(0)
    for _ in range(50):
        o = np.array(nd.sample_origin(g, (8, 16, 16), {}, rng))
        assert np.all(o >= -np.array([2, 4, 4])) and np.all(o + [8, 16, 16] <= np.array(g.shape) + [2, 4, 4])
    batches = list(nd.Prefetch(nd.batches([g], (8, 16, 16), {}, 0, 0, 3, 2)))
    assert len(batches) == 3 and batches[0][0].shape == (2, 2, 8, 16, 16)


def test_truth_positions_use_the_imaris_voxel_centres():
    spots = pd.DataFrame({"imaris_x_um": [10.0 + 0.108 * 5.75], "imaris_y_um": [20.0 + 0.108 * 3.5],
                          "imaris_z_um": [1.0 + 0.2 * 2.5], "xi": [5], "yi": [3], "zi": [2]})
    params = {"ext_min_xyz_um": [10.0, 20.0, 1.0], "voxel_imaris_xyz_um": [0.108, 0.108, 0.2]}
    p = nd.truth_positions(spots, params)
    np.testing.assert_allclose(p, [[2.0, 3.0, 5.25]], atol=1e-9)
    assert np.array_equal(np.floor(p + 0.5).astype(int)[0], [2, 3, 5])     # the truth files' indices
    np.testing.assert_array_equal(nd.truth_positions(spots, None), [[2.0, 3.0, 5.0]])


# ------------------------------------------------------------------------------- backend
def _planted_stack(shape=(16, 64, 64), seed=2):
    rng = np.random.default_rng(seed)
    zz, yy, xx = np.indices(shape).astype(np.float32)
    labels = np.zeros(shape, np.int32)
    labels[3:13, 6:30, 6:30] = 1
    labels[3:13, 34:58, 34:58] = 2
    img = rng.normal(100, 5, size=shape).astype(np.float32)
    for c in [(8, 15, 15), (8, 22, 24), (7, 45, 40), (9, 50, 50)]:
        img += 800 * np.exp(-0.5 * (((zz - c[0]) / 1.5) ** 2 + ((yy - c[1]) / 1.4) ** 2 + ((xx - c[2]) / 1.4) ** 2))
    dna = np.where(labels > 0, 800, 100).astype(np.uint16)
    return np.clip(img, 0, 65535).astype(np.uint16), dna, labels


def _save_random_model(path, *, in_channels=2, threshold=0.5, levels=3, base=4):
    from germquant.spots.net3d.infer import save_model_file

    m = _model(levels=levels, base=base, in_channels=in_channels)
    meta = {"input_roles": ["foci", "dna"][:in_channels], "norm": {"low_pct": 1.0, "high_pct": 99.9, "stride_yx": 4},
            "peaks": {"threshold": threshold, "separation_um": [0.6, 0.3], "assign_max_um": 0.3}}
    return save_model_file(path, m, meta)


def test_backend_contract_with_a_random_model(tmp_path, caplog):
    from germquant import provenance
    from germquant.spots import detect_spots
    from germquant.spots.backends import (PER_NUC_COLS, PER_SPOT_COLS, backend_input_roles, backend_provenance,
                                          list_backends)
    from germquant.spots.net3d.infer import read_model_file

    assert "net3d" in list_backends() and backend_input_roles("net3d") == ("dna",)
    assert backend_input_roles("spotmax") == () and backend_input_roles("log") == ()
    img, dna, labels = _planted_stack()
    path = tmp_path / "model.pt"
    sha = _save_random_model(path)
    assert read_model_file(path)["input_roles"] == ["foci", "dna"]
    section = {"net3d": {"model": str(path), "threshold": 0.0, "device": "cpu"}}
    per_spot, per_nuc = detect_spots(img, labels, SP, backend="net3d", cfg_section=section, channels={"dna": dna})
    assert list(per_spot.columns[:9]) == PER_SPOT_COLS and list(per_spot.columns[9:]) == ["prob", "assignment", "contested"]
    assert list(per_nuc.columns) == PER_NUC_COLS
    assert len(per_spot) > 0 and (per_spot["detector"] == "net3d").all() and (per_nuc["detector"] == "net3d").all()
    assert per_nuc["nucleus_id"].tolist() == [1, 2] and int(per_nuc["n_spots"].sum()) == len(per_spot)
    assert per_spot["nucleus_id"].isin([1, 2]).all() and per_spot["spot_id"].tolist() == list(range(len(per_spot)))
    assert per_spot["assignment"].isin(["label", "nearest"]).all() and (per_spot["prob"] >= 0.0).all()
    # deterministic: the same call gives the same table
    again, _ = detect_spots(img, labels, SP, backend="net3d", cfg_section=section, channels={"dna": dna})
    pd.testing.assert_frame_equal(per_spot, again)
    # a high threshold keeps a subset; a threshold above every probability keeps nothing (a real zero per nucleus)
    hi, _ = detect_spots(img, labels, SP, backend="net3d", channels={"dna": dna},
                         cfg_section={"net3d": {"model": str(path), "threshold": float(per_spot["prob"].median()),
                                                "device": "cpu"}})
    assert 0 < len(hi) <= len(per_spot)
    none, none_nuc = detect_spots(img, labels, SP, backend="net3d", channels={"dna": dna},
                                  cfg_section={"net3d": {"model": str(path), "threshold": 1.01, "device": "cpu"}})
    assert none.empty and list(none.columns[:9]) == PER_SPOT_COLS and none_nuc["n_spots"].tolist() == [0, 0]
    # a model trained with DAPI refuses an image without it; a foci-only model needs nothing
    with pytest.raises(ValueError, match="dna"):
        detect_spots(img, labels, SP, backend="net3d", cfg_section=section)
    p1 = tmp_path / "foci_only.pt"
    _save_random_model(p1, in_channels=1)
    s1, _ = detect_spots(img, labels, SP, backend="net3d",
                         cfg_section={"net3d": {"model": str(p1), "threshold": 0.0, "device": "cpu"}})
    assert len(s1) > 0
    with pytest.raises(ValueError, match="spots.net3d.model"):
        detect_spots(img, labels, SP, backend="net3d", cfg_section={}, channels={"dna": dna})
    # provenance: the model file's sha256, also through the registry hook; never fatal
    rec = backend_provenance("net3d", section)
    assert rec["model_sha256"] == sha == provenance.file_sha256(path) and rec["backend"] == "net3d"
    assert backend_provenance("net3d", {})["model_sha256"] == "absent"
    assert backend_provenance("spotmax", {}) == {} and backend_provenance("nope", {}) == {}
    caplog.set_level("WARNING")
    detect_spots(img, labels, SP, backend="net3d", channels={"dna": dna},
                 cfg_section={"net3d": {"model": str(path), "device": "cpu", "threshold": 2.0, "bogus": 1}})
    assert any("unknown keys ['bogus']" in r.getMessage() for r in caplog.records)


def test_channels_reach_only_backends_that_ask(monkeypatch):
    from germquant.spots import backends, detect_spots
    from germquant.spots.backends import PER_NUC_COLS, PER_SPOT_COLS

    seen = {}

    def plain(image, labels, spacing, *, cfg_section, marker):
        seen["plain"] = True
        return (pd.DataFrame(columns=PER_SPOT_COLS),
                pd.DataFrame({"nucleus_id": [1], "marker": marker, "n_spots": [0], "mean_spot_intensity": [np.nan],
                              "detector": "plain_t"}, columns=PER_NUC_COLS))

    def asks(image, labels, spacing, *, cfg_section, marker, channels):
        seen["asks"] = sorted(channels)
        return plain(image, labels, spacing, cfg_section=cfg_section, marker=marker)

    asks.input_roles = ("dna",)
    monkeypatch.setitem(backends._REGISTRY, "plain_t", plain)
    monkeypatch.setitem(backends._REGISTRY, "asks_t", asks)
    img = np.zeros((2, 4, 4), np.uint16)
    lab = np.ones((2, 4, 4), np.int32)
    detect_spots(img, lab, SP, backend="plain_t", channels={"dna": img})          # silently not passed
    detect_spots(img, lab, SP, backend="asks_t", channels={"dna": img})
    assert seen == {"plain": True, "asks": ["dna"]}


def test_default_path_hashes_and_manifest_are_unchanged(tmp_path):
    from germquant import pipeline
    from germquant.config import load_config
    from germquant.stages import stage_hashes

    cfg = load_config(ROOT / "config" / "config.yaml")
    roles = {"dna": 0, "central_element": 1, "foci": 2, "granule": None}
    assert pipeline.spot_backend_provenance(cfg, roles) == {}
    tools = {"spotmax": "x"}
    assert stage_hashes(cfg, roles, tools, "m") == stage_hashes(cfg, roles, tools, "m", spot_backend={})
    path = tmp_path / "m.pt"
    sha = _save_random_model(path)
    cfg.set("spots.backend", "net3d")
    cfg.set("spots.net3d", {"model": str(path)})
    prov = pipeline.spot_backend_provenance(cfg, roles)
    assert prov == {"spots": {"backend": "net3d", "model": str(path), "model_sha256": sha, "input_roles": {"dna": 0}}}
    h1 = stage_hashes(cfg, roles, tools, "m", spot_backend=prov)
    h0 = stage_hashes(cfg, roles, tools, "m")
    assert h1["spots"] != h0["spots"] and {k: v for k, v in h1.items() if k != "spots"} == \
        {k: v for k, v in h0.items() if k != "spots"}


def test_pipeline_runs_net3d_with_the_dna_channel(tmp_path, monkeypatch):
    from test_end_to_end import CHANNELS, SPACING, _synthetic_stack

    from germquant import pipeline
    from germquant.config import load_config
    from germquant.spots import backends

    stack = _synthetic_stack()
    meta = {"sizes": {"C": 3, "Z": 20, "Y": 80, "X": 320}, "spacing": SPACING, "spacing_ok": True,
            "channel_names": CHANNELS, "dtype": "float32", "is_2d": False}
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    monkeypatch.setattr(pipeline, "read_stack", lambda *a, **k: stack)
    calls = []
    real = backends._REGISTRY["net3d"]

    def spy(image, labels, spacing, *, cfg_section, marker, channels=None):
        calls.append(channels)
        return real(image, labels, spacing, cfg_section=cfg_section, marker=marker, channels=channels)

    spy.input_roles = real.input_roles
    spy.provenance = real.provenance
    monkeypatch.setitem(backends._REGISTRY, "net3d", spy)
    path = tmp_path / "m.pt"
    sha = _save_random_model(path)
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    cfg.set("spots.backend", "net3d")
    cfg.set("spots.net3d", {"model": str(path), "threshold": 0.0, "device": "cpu"})
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / "r", prov=None)
    assert res["stages"]["spots"]["status"] == "ran", res["stages"]["spots"]
    assert len(calls) == 1 and np.array_equal(calls[0]["dna"], stack.data[0])
    spots = res["tables"]["spots"]
    assert (spots["detector"] == "net3d").all() and "spots:unvalidated_backend=net3d" in res["qc_flags"]
    man = json.loads((tmp_path / "r" / "run_manifest.json").read_text(encoding="utf-8"))
    assert man["spot_backend"]["spots"]["model_sha256"] == sha
    written = pd.read_csv(tmp_path / "r" / f"{res['image_id']}__spots.csv")
    assert {"prob", "assignment", "contested"} <= set(written.columns)


# ------------------------------------------------------------------- the CLI and the harness
@pytest.fixture
def tiny_bench(tmp_path):
    from test_bench_detector import _write_gonad

    root = tmp_path
    gonads = [_write_gonad(root, "synth_noHS_a", "train", "noHS", (2, 4, 6, 3), 11),
              _write_gonad(root, "synth_HS_b", "train", "HS", (5, 1, 3, 7), 12),
              _write_gonad(root, "synth_val_c", "val", "noHS", (3, 3, 5, 2), 13),
              _write_gonad(root, "synth_test_d", "test", "HS", (4, 2, 2, 6), 14)]
    unusable = dict(gonads[0]["entry"], tag="synth_rg_e", spot_truth_usable=False)
    manifest = root / "manifest.json"
    manifest.write_text(json.dumps({"gonads": [g["entry"] for g in gonads] + [unusable]}), encoding="utf-8")
    split = root / "split.json"
    split.write_text(json.dumps({"roles": {"train": ["synth_noHS_a", "synth_HS_b", "synth_rg_e"],
                                           "val": ["synth_val_c"], "test": ["synth_test_d"]}}), encoding="utf-8")
    return {"root": root, "manifest": manifest, "split": split}


def test_train_spots_cli_one_epoch_on_a_tiny_cache(tiny_bench, monkeypatch, capsys):
    from germquant import cli
    from germquant.benchmark import detector as bd
    from germquant.spots.net3d.infer import file_sha256, read_model_file

    b = tiny_bench
    opened = []
    real_load = bd.load_cache

    def spy_load(tag, *a, **k):
        opened.append(tag)
        return real_load(tag, *a, **k)

    monkeypatch.setattr(bd, "load_cache", spy_load)
    out = b["root"] / "model"
    rc = cli.main(["train-spots", "--cache", str(b["root"] / "cache"), "--split-file", str(b["split"]),
                   "--out", str(out), "--epochs", "1", "--patch", "16,32,32", "--seed", "3",
                   "--manifest", str(b["manifest"]), "--batch-size", "2", "--steps-per-epoch", "2",
                   "--device", "cpu", "--base-channels", "4", "--levels", "3", "--thresholds", "0.02,0.3,0.6"])
    text = capsys.readouterr().out
    assert rc == 0, text
    assert "synth_test_d" not in opened and "synth_rg_e" not in opened
    assert set(opened) == {"synth_noHS_a", "synth_HS_b", "synth_val_c"}
    for f in ("model.pt", "training_log.json", "model_card.json", "checkpoints/best.pt", "checkpoints/last.pt"):
        assert (out / f).is_file(), f
    logd = json.loads((out / "training_log.json").read_text(encoding="utf-8"))
    card = json.loads((out / "model_card.json").read_text(encoding="utf-8"))
    assert logd["status"] == "finished" and logd["best_epoch"] == 1 and len(logd["epochs"]) == 1
    assert logd["model_sha256"] == card["model_sha256"] == file_sha256(out / "model.pt")
    assert logd["checkpoints"][0]["sha256"] == file_sha256(out / "checkpoints" / "best.pt")
    assert logd["epochs"][0]["last_sha256"] == file_sha256(out / "checkpoints" / "last.pt")
    assert all(c["max_abs_diff"] <= 1e-9 for c in logd["harness_check"])
    assert [g["tag"] for g in card["training_gonads"]] == ["synth_noHS_a", "synth_HS_b"]
    assert [g["tag"] for g in card["val_gonads"]] == ["synth_val_c"]
    assert card["excluded_gonads"][0]["tag"] == "synth_rg_e" and card["test_gonads_never_opened"] == ["synth_test_d"]
    assert card["truth_classes_used"]["targets"] == ["Class A"] and card["epochs"]["run"] == 1
    assert card["threshold"] in (0.02, 0.3, 0.6) and card["git_sha"]
    sweep = logd["epochs"][0]["val"]["sweep"]
    assert [s["threshold"] for s in sweep] == [0.02, 0.3, 0.6] and "ccc_traced" in sweep[0]
    doc = read_model_file(out / "model.pt")
    assert doc["peaks"]["threshold"] == card["threshold"] and doc["input_roles"] == ["foci", "dna"]
    assert doc["training"]["train_gonads"] == ["synth_noHS_a", "synth_HS_b"]
    # the same seed gives the same weights
    out2 = b["root"] / "model2"
    assert cli.main(["train-spots", "--cache", str(b["root"] / "cache"), "--split-file", str(b["split"]),
                     "--out", str(out2), "--epochs", "1", "--patch", "16,32,32", "--seed", "3",
                     "--manifest", str(b["manifest"]), "--batch-size", "2", "--steps-per-epoch", "2",
                     "--device", "cpu", "--base-channels", "4", "--levels", "3", "--thresholds", "0.02,0.3,0.6",
                     "--prep-dir", str(out / "prep")]) == 0
    w1, w2 = read_model_file(out / "model.pt")["state_dict"], read_model_file(out2 / "model.pt")["state_dict"]
    assert all(torch.equal(w1[k], w2[k]) for k in w1)
    # the bench harness scores the model on the val gonad (the DNA channel comes from the cache)
    import importlib.util

    spec = importlib.util.spec_from_file_location("bench_detector", ROOT / "scripts" / "bench_detector.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    params = json.dumps({"net3d": {"model": str(out / "model.pt"), "device": "cpu"}})
    rc = mod.main(["--backend", "net3d", "--params", params, "--split", "val", "--out", str(b["root"] / "bench"),
                   "--manifest", str(b["manifest"]), "--cache-root", str(b["root"] / "cache"),
                   "--trials-log", str(b["root"] / "trials.jsonl"), "--n-boot", "20"])
    assert rc == 0, capsys.readouterr().out
    pg = pd.read_csv(b["root"] / "bench" / "per_gonad.csv")
    assert pg["tag"].tolist() == ["synth_val_c"] and pg["backend"].tolist() == ["net3d"]
    assert bd.read_jsonl(b["root"] / "trials.jsonl")[0]["backend"] == "net3d"
    # the refusal of a split file that moves a gonad
    bad = b["root"] / "bad_split.json"
    bad.write_text(json.dumps({"roles": {"train": ["synth_noHS_a"], "val": ["synth_HS_b"], "test": []}}),
                   encoding="utf-8")
    from germquant.spots.net3d.train import resolve_gonads

    with pytest.raises(ValueError, match="manifest says"):
        resolve_gonads(bad, b["manifest"])
    with pytest.raises(ValueError, match="pre-registered"):
        moved = b["root"] / "moved.json"
        moved.write_text(json.dumps({"roles": {"train": ["n2rad51_HS_male_001"]}}), encoding="utf-8")
        resolve_gonads(moved, b["manifest"])
    assert math.isfinite(float(card["val_metrics_at_best_epoch"]["n_expert"]))
