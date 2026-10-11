"""The pack knobs of package P2a (docs/DESIGN_system_packs.md sections 3.10, 3.12, 4.2, 4.5 and 4.7), each off by
default and each pinned to today's call when off:

* Cellpose, mocked: the default 3D call is exactly today's (do_3D=True, z_axis=0, channel_axis=None,
  anisotropy=dz/dy, no diameter); segmentation.nuclei.cellpose.object_diameter_um adds only diameter = d / dy;
  the 2D call of a projected run is model.eval(img2d, channel_axis=None), with the same model and CPU retry;
* SpotMAX, spied: without the optional keys (spot_radius_z_um, remove_hot_pixels, sharpen) the backend section
  and the preprocess, segmentation, detection and feature calls are exactly today's; with sharpen the detection
  gets the sharpened image and the feature call the unsharpened one plus sharp_spots_image, and the effect size
  is the one SpotMAX computes on the sharpened image; spot_radius_z_um makes the radii (rz, r, r) / spacing;
* segmentation.nuclei.require_channel fails the image before segmentation when the segmentation role is missing;
* qc.min_nuclei moves the qc:few_nuclei_<n> threshold only when set.
"""
from __future__ import annotations

import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_end_to_end import patched_reader  # noqa: F401  (fixture: the 3-channel synthetic gonad)
from test_segment_backends import registry  # noqa: F401  (fixture: a private segmentation backend registry)

from germquant import device as D
from germquant import pipeline, planar, qc
from germquant.config import load_config
from germquant.segment import backends as B
from germquant.segment import nuclei as N
from germquant.spots import detect as SD
from germquant.spots import detect_spots
from germquant.spots.detect import SPOTMAX_DEFAULTS, SPOTMAX_KEYS, SPOTMAX_OPTIONAL_KEYS, detect_spots_spotmax

ROOT = Path(__file__).resolve().parents[1]
ND2 = "20251105_n2_nohs_HERM_001.nd2"


# ---- Cellpose, mocked ------------------------------------------------------------------------------------

class _FakeCellposeModel:
    """Records how the model is built and every eval call; returns the voxels above the mean as one object."""
    inits: list = []   # noqa: RUF012 - test spy, reset by the fixture
    evals: list = []   # noqa: RUF012

    def __init__(self, **kw):
        _FakeCellposeModel.inits.append(kw)

    def eval(self, x, **kw):
        _FakeCellposeModel.evals.append({"x": x, "kw": kw})
        a = np.asarray(x, dtype=np.float64)
        return [(a > a.mean()).astype(np.int32), None, None]


@pytest.fixture
def fake_cellpose(monkeypatch):
    fake_models = types.SimpleNamespace(CellposeModel=_FakeCellposeModel)
    monkeypatch.setitem(sys.modules, "cellpose", types.SimpleNamespace(models=fake_models))
    monkeypatch.setitem(sys.modules, "cellpose.models", fake_models)
    monkeypatch.setattr(D, "available", lambda: {"cuda": False, "mps": False})
    monkeypatch.setattr(D, "torch_device", lambda name: f"device:{name}")
    monkeypatch.delenv("GERMQUANT_DEVICE", raising=False)
    _FakeCellposeModel.inits.clear()
    _FakeCellposeModel.evals.clear()
    return _FakeCellposeModel


def _volume(shape=(6, 16, 16), seed=1):
    return np.random.default_rng(seed).normal(100.0, 10.0, shape).astype(np.float32)


TODAY_3D = {"do_3D": True, "z_axis": 0, "channel_axis": None}


def test_the_default_3d_cellpose_call_is_exactly_todays(fake_cellpose, registry):  # noqa: F811
    section = load_config(ROOT / "config" / "config.yaml").get("segmentation.nuclei")
    dna = _volume()
    _labels, used = N.segment_nuclei(dna, (0.4, 0.2, 0.2), method="cellpose", cellpose_model="models/m",
                                     cfg_section=section, min_volume_um3=0.0)
    assert used == "cellpose"
    (call,) = fake_cellpose.evals
    assert call["x"] is dna
    assert call["kw"] == {**TODAY_3D, "anisotropy": 0.4 / 0.2}          # no diameter, nothing else
    assert fake_cellpose.inits == [{"pretrained_model": "models/m", "device": "device:cpu", "gpu": False,
                                    "use_bfloat16": False}]
    # the built-in model name builds the default model, the same call
    N.segment_nuclei(dna, (0.5, 0.25, 0.25), method="cellpose", cellpose_model="cpsam", min_volume_um3=0.0)
    assert fake_cellpose.inits[-1] == {"device": "device:cpu", "gpu": False, "use_bfloat16": False}
    assert fake_cellpose.evals[-1]["kw"] == {**TODAY_3D, "anisotropy": 2.0}


def test_the_default_pipeline_run_makes_todays_cellpose_call(tmp_path, fake_cellpose, patched_reader):  # noqa: F811
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg.set("segmentation.nuclei.method", "cellpose")
    cfg.set("segmentation.nuclei.cellpose_model", "cpsam")
    cfg.set("spots.enabled", False)
    res = pipeline.process_image(ND2, cfg, tmp_path / "r", prov=None)
    assert res["stages"]["segment"]["status"] == "ran"
    (call,) = fake_cellpose.evals
    assert call["x"].ndim == 3 and call["kw"] == {**TODAY_3D, "anisotropy": 0.4 / 0.2}
    assert not any(f.startswith("segment:projected") for f in res["qc_flags"])


def test_object_diameter_adds_only_the_diameter_in_pixels(fake_cellpose, registry):  # noqa: F811
    dna = _volume()
    section = {"method": "cellpose", "cellpose": {"object_diameter_um": 6.0}}
    N.segment_nuclei(dna, (0.4, 0.4, 0.4), method="cellpose", cellpose_model="cpsam", cfg_section=section,
                     min_volume_um3=0.0)
    assert fake_cellpose.evals[-1]["kw"] == {**TODAY_3D, "anisotropy": 1.0, "diameter": 15.0}
    # dy is the pixel size after xy_stride: a stride of 2 doubles it and halves the diameter in pixels
    N.segment_nuclei(dna, (0.4, 0.8, 0.8), method="cellpose", cellpose_model="cpsam", cfg_section=section,
                     min_volume_um3=0.0)
    assert fake_cellpose.evals[-1]["kw"] == {**TODAY_3D, "anisotropy": 0.5, "diameter": 7.5}
    # unset (absent, an empty cellpose mapping, null): exactly today's call
    for sec in ({}, {"cellpose": {}}, {"cellpose": {"object_diameter_um": None}}, {"cellpose": None}):
        N.segment_nuclei(dna, (0.4, 0.2, 0.2), method="cellpose", cfg_section=sec, min_volume_um3=0.0)
        assert fake_cellpose.evals[-1]["kw"] == {**TODAY_3D, "anisotropy": 2.0}
    assert B.object_diameter_px({"cellpose": {"object_diameter_um": 3}}, (1.0, 0.2, 0.2)) == pytest.approx(15.0)
    for bad in ({"cellpose": {"object_diameter_um": 0}}, {"cellpose": {"object_diameter_um": -2.0}},
                {"cellpose": 4.0}):
        with pytest.raises(ValueError, match="segmentation.nuclei.cellpose"):
            B.object_diameter_px(bad, (1.0, 0.2, 0.2))
        with pytest.raises(ValueError):
            N.segment_nuclei(dna, (0.4, 0.2, 0.2), method="cellpose", cfg_section=bad)


def test_the_2d_call_is_eval_of_the_plane_with_no_3d_arguments(fake_cellpose, registry):  # noqa: F811
    img = _volume(shape=(1, 20, 24))
    labels, used = planar.segment_planar(img, (0.4, 0.2, 0.2), method="cellpose", cellpose_model="models/m",
                                         cfg_section={"method": "cellpose"})
    assert used == "cellpose" and labels.shape == (1, 20, 24) and labels.dtype == np.int32 and labels.max() == 1
    (call,) = fake_cellpose.evals
    assert call["x"].shape == (20, 24) and np.array_equal(call["x"], img[0])
    assert call["kw"] == {"channel_axis": None}                       # no do_3D, z_axis or anisotropy
    planar.segment_planar(img, (0.4, 0.25, 0.25), method="cellpose", cellpose_model="cpsam",
                          cfg_section={"cellpose": {"object_diameter_um": 5.0}})
    assert fake_cellpose.evals[-1]["kw"] == {"channel_axis": None, "diameter": 20.0}
    with pytest.raises(ValueError, match="one \\(Y, X\\) plane"):
        N._cellpose_2d(img, "cpsam")


@pytest.mark.parametrize("have", [{"cuda": False, "mps": False}, {"cuda": True, "mps": False},
                                  {"cuda": False, "mps": True}])
def test_the_2d_model_is_built_exactly_as_the_3d_model(fake_cellpose, monkeypatch, have):
    monkeypatch.setattr(D, "available", lambda: have)
    for model in ("cpsam", "models/m"):
        N._cellpose(_volume(), (0.4, 0.2, 0.2), model, 3.0)
        N._cellpose_2d(_volume(shape=(16, 16)), model)
        assert fake_cellpose.inits[-1] == fake_cellpose.inits[-2]


def test_the_planar_and_diameter_branches_keep_the_cpu_retry(monkeypatch):
    seen = []

    def fake_2d(img2d, model, device_pref=None, diameter_px=None):
        seen.append(("2d", img2d.shape, device_pref, diameter_px))
        if device_pref != "cpu":
            raise RuntimeError("CUDA out of memory")
        return np.ones((1, *img2d.shape), np.int32)

    def fake_3d(dna, spacing, model, diam, device_pref=None, diameter_px=None):
        seen.append(("3d", dna.shape, device_pref, diameter_px))
        if device_pref != "cpu":
            raise RuntimeError("CUDA out of memory")
        return np.ones(dna.shape, np.int32)

    monkeypatch.setattr(N, "_cellpose_2d", fake_2d)
    monkeypatch.setattr(N, "_cellpose", fake_3d)
    monkeypatch.setattr(D, "available", lambda: {"cuda": True, "mps": False})
    monkeypatch.delenv("GERMQUANT_DEVICE", raising=False)
    out = B.cellpose_backend(np.zeros((1, 8, 9), np.float32), (0.4, 0.2, 0.2), model="m", diameter_um=3.0,
                             device="cuda", cfg_section={}, planar=True)
    assert out.shape == (1, 8, 9)
    assert seen == [("2d", (8, 9), "cuda", None), ("2d", (8, 9), "cpu", None)]
    seen.clear()
    B.cellpose_backend(np.zeros((3, 8, 9), np.float32), (0.4, 0.2, 0.2), model="m", diameter_um=3.0, device=None,
                       cfg_section={"cellpose": {"object_diameter_um": 2.0}})
    assert seen == [("3d", (3, 8, 9), None, 10.0), ("3d", (3, 8, 9), "cpu", 10.0)]
    seen.clear()
    monkeypatch.setattr(D, "available", lambda: {"cuda": False, "mps": False})   # already on the CPU: no retry
    with pytest.raises(RuntimeError, match="out of memory"):
        B.cellpose_backend(np.zeros((1, 8, 9), np.float32), (0.4, 0.2, 0.2), model="m", diameter_um=3.0,
                           device=None, cfg_section={}, planar=True)
    assert seen == [("2d", (8, 9), None, None)]
    with pytest.raises(ValueError, match="one plane as \\(1, Y, X\\)"):
        B.cellpose_backend(np.zeros((2, 8, 9), np.float32), (0.4, 0.2, 0.2), model="m", diameter_um=3.0,
                           device=None, cfg_section={}, planar=True)


# ---- SpotMAX: the optional keys --------------------------------------------------------------------------

def _two_nuclei():
    labels = np.zeros((4, 10, 10), np.int32)
    labels[1:3, 1:4, 1:4] = 1
    labels[1:3, 5:8, 5:8] = 2
    img = np.random.default_rng(5).normal(50.0, 5.0, labels.shape).astype(np.float32)
    return img, labels


class _FakeSpotmax:
    """spotmax.pipe and spotmax.filters as call recorders: each step returns a new array, so the test can tell
    which image went where; detection finds one spot in each nucleus."""

    def __init__(self):
        self.calls: list[tuple[str, tuple, dict]] = []
        self.out: dict[str, object] = {}
        self.pipe = types.SimpleNamespace(
            preprocess_image=self._rec("preprocess_image", lambda image, **kw: np.asarray(image) + 1.0),
            spots_semantic_segmentation=self._rec("spots_semantic_segmentation",
                                                  lambda image, **kw: np.ones(np.shape(image), bool)),
            spot_detection=self._rec("spot_detection", self._detect),
            spots_calc_features_and_filter=self._rec("spots_calc_features_and_filter", self._features),
        )
        self.filters = types.SimpleNamespace(
            DoG_spots=self._rec("DoG_spots", lambda image, radii, **kw: np.asarray(image) * 2.0))

    def _rec(self, name, fn):
        def wrapped(*args, **kw):
            self.calls.append((name, args, kw))
            out = fn(*args, **kw)
            self.out[name] = out
            return out
        return wrapped

    @staticmethod
    def _detect(image, **kw):
        return pd.DataFrame({"z": [1, 2], "y": [2, 6], "x": [2, 6]}), None

    @staticmethod
    def _features(image, radii, df, **kw):
        feat = pd.DataFrame({"z": df["z"].to_numpy(), "y": df["y"].to_numpy(), "x": df["x"].to_numpy(),
                             SD._EFFECT: [4.0, 5.0], "spot_raw_mean_in_spot": [70.0, 80.0]})
        return None, [feat]

    def call(self, name):
        hits = [c for c in self.calls if c[0] == name]
        assert len(hits) == 1, f"{name} called {len(hits)} times"
        return hits[0]


@pytest.fixture
def fake_spotmax(monkeypatch):
    fake = _FakeSpotmax()
    monkeypatch.setitem(sys.modules, "spotmax", types.SimpleNamespace(pipe=fake.pipe, filters=fake.filters))
    monkeypatch.setitem(sys.modules, "spotmax.pipe", fake.pipe)
    monkeypatch.setitem(sys.modules, "spotmax.filters", fake.filters)
    return fake


def _same(a, b) -> bool:
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        return np.array_equal(np.asarray(a), np.asarray(b)) and np.asarray(a).dtype == np.asarray(b).dtype
    return a == b


def _kw_equal(got: dict, want: dict) -> None:
    assert sorted(got) == sorted(want), f"keywords {sorted(got)} != today's {sorted(want)}"
    for k, v in want.items():
        assert _same(got[k], v), f"keyword {k}: {got[k]!r} != {v!r}"


PIPELINE_SPOT_ARGS = dict(spot_radius_um=0.3, gauss_sigma_um=0.08, thresholding_method="threshold_triangle",
                          effect_size_metric="spot_vs_backgr_effect_size_glass", effect_size_min=3.0,
                          merge_z_columns=True, z_merge_gap_um=0.8, z_merge_valley_frac=0.8,
                          max_spot_candidates=30000)


def test_without_the_optional_keys_the_backend_section_is_todays(monkeypatch):
    """SPOTMAX_DEFAULTS is not extended, and the spotmax backend hands detect_spots_spotmax exactly the
    resolved SPOTMAX_KEYS for the default config's spots: section (the pipeline's call)."""
    assert set(SPOTMAX_DEFAULTS) == set(SPOTMAX_KEYS)
    assert not set(SPOTMAX_OPTIONAL_KEYS) & set(SPOTMAX_DEFAULTS)
    assert SPOTMAX_OPTIONAL_KEYS == ("spot_radius_z_um", "remove_hot_pixels", "sharpen")
    seen = []

    def spy(image, labels, spacing, **kw):
        seen.append(kw)
        return SD._empty(labels, kw["marker"])

    monkeypatch.setattr(SD, "detect_spots_spotmax", spy)
    img, labels = _two_nuclei()
    section = load_config(ROOT / "config" / "config.yaml").get("spots")
    detect_spots(img, labels, (0.4, 0.2, 0.2), marker="RAD-51", cfg_section=section, **PIPELINE_SPOT_ARGS)
    assert seen == [{"marker": "RAD-51", **PIPELINE_SPOT_ARGS}]
    # with the optional keys in the raw section they are forwarded, and only they
    opt = {**section, "spot_radius_z_um": 0.6, "remove_hot_pixels": True, "sharpen": True}
    detect_spots(img, labels, (0.4, 0.2, 0.2), marker="RAD-51", cfg_section=opt, **PIPELINE_SPOT_ARGS)
    assert seen[-1] == {"marker": "RAD-51", **PIPELINE_SPOT_ARGS, "spot_radius_z_um": 0.6,
                        "remove_hot_pixels": True, "sharpen": True}
    detect_spots(img, labels, (0.4, 0.2, 0.2), marker="RAD-51", cfg_section={**section, "sharpen": None},
                 **PIPELINE_SPOT_ARGS)                                # a null value is no value
    assert seen[-1] == {"marker": "RAD-51", **PIPELINE_SPOT_ARGS}


def test_without_the_optional_keys_the_spotmax_calls_are_todays(fake_spotmax):
    img, labels = _two_nuclei()
    sp = (0.4, 0.2, 0.25)
    per_spot, per_nuc = detect_spots_spotmax(img, labels, sp, marker="RAD-51", **PIPELINE_SPOT_ARGS)
    spa = np.asarray(sp, float)
    radii = 0.3 / spa
    names = [c[0] for c in fake_spotmax.calls]
    assert names == ["preprocess_image", "spots_semantic_segmentation", "spot_detection",
                     "spots_calc_features_and_filter"]                 # no DoG_spots
    _, args, kw = fake_spotmax.call("preprocess_image")
    assert len(args) == 1 and _same(args[0], img.astype("float32"))
    _kw_equal(kw, {"lab": labels, "gauss_sigma": 0.08 / 0.2, "use_gpu": False, "spots_zyx_radii_pxl": radii})
    pre = fake_spotmax.out["preprocess_image"]
    _, args, kw = fake_spotmax.call("spots_semantic_segmentation")
    assert len(args) == 1 and args[0] is pre
    _kw_equal(kw, {"lab": labels, "spots_zyx_radii_pxl": radii, "do_try_all_thresholds": False,
                   "return_only_segm": True, "thresholding_method": "threshold_triangle", "use_gpu": False})
    _, args, kw = fake_spotmax.call("spot_detection")
    assert len(args) == 1 and args[0] is pre
    _kw_equal(kw, {"spots_segmantic_segm": fake_spotmax.out["spots_semantic_segmentation"],
                   "spots_zyx_radii_pxl": radii, "lab": labels, "return_df": True})
    _, args, kw = fake_spotmax.call("spots_calc_features_and_filter")
    assert len(args) == 3 and args[0] is pre and _same(args[1], radii)
    assert list(args[2].index) == [1, 2] and args[2].index.name == "Cell_ID"
    _kw_equal(kw, {"lab": labels, "raw_image": img.astype("float32"), "zyx_voxel_size": tuple(spa),
                   "gop_filtering_thresholds": None})
    assert per_spot["effect_size"].tolist() == [4.0, 5.0] and per_nuc["n_spots"].tolist() == [1, 1]


def test_sharpen_feeds_detection_and_gives_the_feature_call_both_images(fake_spotmax):
    img, labels = _two_nuclei()
    sp = (0.4, 0.2, 0.2)
    detect_spots_spotmax(img, labels, sp, marker="RAD-51", sharpen=True, **PIPELINE_SPOT_ARGS)
    radii = 0.3 / np.asarray(sp, float)
    pre = fake_spotmax.out["preprocess_image"]
    sharp = fake_spotmax.out["DoG_spots"]
    _, args, kw = fake_spotmax.call("preprocess_image")
    _kw_equal(kw, {"lab": labels, "gauss_sigma": 0.08 / 0.2, "use_gpu": False, "spots_zyx_radii_pxl": radii})
    _, args, kw = fake_spotmax.call("DoG_spots")                      # SpotMAX core's sharpen_spots on pre
    assert args[0] is pre and _same(args[1], radii)
    _kw_equal(kw, {"use_gpu": False, "lab": labels})
    assert fake_spotmax.call("spots_semantic_segmentation")[1][0] is sharp
    assert fake_spotmax.call("spot_detection")[1][0] is sharp
    _, args, kw = fake_spotmax.call("spots_calc_features_and_filter")
    assert args[0] is pre and kw["sharp_spots_image"] is sharp        # the unsharpened image plus the sharp one
    _kw_equal({k: v for k, v in kw.items() if k != "sharp_spots_image"},
              {"lab": labels, "raw_image": img.astype("float32"), "zyx_voxel_size": tuple(np.asarray(sp, float)),
               "gop_filtering_thresholds": None})


def test_spot_radius_z_and_hot_pixels_reach_every_spotmax_call(fake_spotmax):
    img, labels = _two_nuclei()
    sp = (0.4, 0.2, 0.25)
    detect_spots_spotmax(img, labels, sp, marker="RAD-51", spot_radius_z_um=0.8, remove_hot_pixels=True,
                         **PIPELINE_SPOT_ARGS)
    radii = np.array([0.8, 0.3, 0.3]) / np.asarray(sp, float)        # (rz, r, r) / spacing
    np.testing.assert_allclose(radii, [2.0, 1.5, 1.2])
    _, _, kw = fake_spotmax.call("preprocess_image")
    assert kw["do_remove_hot_pixels"] is True and _same(kw["spots_zyx_radii_pxl"], radii)
    for name in ("spots_semantic_segmentation", "spot_detection"):
        assert _same(fake_spotmax.call(name)[2]["spots_zyx_radii_pxl"], radii)
    assert _same(fake_spotmax.call("spots_calc_features_and_filter")[1][1], radii)
    assert "DoG_spots" not in [c[0] for c in fake_spotmax.calls]


def _planted_foci(seed=4):
    """Two spherical nuclei with Gaussian foci (a small stand-in for a real foci channel), uint16."""
    sp = np.array([0.25, 0.1, 0.1])
    shape = (14, 48, 80)
    rng = np.random.default_rng(seed)
    zz, yy, xx = (np.indices(shape).astype(np.float32) * sp[:, None, None, None])
    labels = np.zeros(shape, np.int32)
    img = np.full(shape, 40.0, np.float32)
    for k, (cy, cx) in enumerate(((2.4, 2.0), (2.4, 5.8)), start=1):
        inside = (zz - 1.75) ** 2 + (yy - cy) ** 2 + (xx - cx) ** 2 <= 1.6 ** 2
        labels[inside] = k
        img[inside] = 150.0
        for dz, dy, dx in ((0.0, -0.6, -0.5), (0.3, 0.5, 0.4), (-0.3, 0.2, -0.7)):
            d2 = (zz - 1.75 - dz) ** 2 / 0.3 ** 2 + (yy - cy - dy) ** 2 / 0.12 ** 2 + (xx - cx - dx) ** 2 / 0.12 ** 2
            img += 900.0 * np.exp(-d2 / 2)
    img += rng.normal(0, 8, shape).astype(np.float32)
    return np.clip(np.round(img), 0, 65535).astype(np.uint16), labels, tuple(sp)


def test_sharpen_puts_the_effect_size_on_the_sharpened_scale(monkeypatch):
    """With sharpen, germquant's effect_size is the value SpotMAX computes on the sharpened image: the feature
    call given (pre, sharp_spots_image=sharp) agrees with the same call given the sharpened image itself, and
    differs from the value on the unsharpened image."""
    P = pytest.importorskip("spotmax.pipe")
    img, labels, sp = _planted_foci()
    real = P.spots_calc_features_and_filter
    seen: dict = {}

    def wrapped(image, radii, df, **kw):
        seen.update(image=image, radii=radii, df=df.copy(), kw=dict(kw))
        return real(image, radii, df, **kw)

    monkeypatch.setattr(P, "spots_calc_features_and_filter", wrapped)
    params = dict(spot_radius_um=0.25, spot_radius_z_um=0.5, gauss_sigma_um=0.08,
                  thresholding_method="threshold_triangle", effect_size_min=0.0, merge_z_columns=False)
    per_spot, _ = detect_spots_spotmax(img, labels, sp, sharpen=True, **params)
    assert len(per_spot) >= 4
    sharp = seen["kw"].pop("sharp_spots_image")
    assert sharp is not None and not np.array_equal(sharp, seen["image"])

    def effect(image, **extra):
        out = real(image, seen["radii"], seen["df"].copy(), **seen["kw"], **extra)
        feat = pd.concat(out[1]).reset_index(drop=True)
        return feat.sort_values(["z", "y", "x"])[SD._EFFECT].to_numpy(float)

    on_both = effect(seen["image"], sharp_spots_image=sharp)
    on_sharp = effect(sharp)
    on_pre = effect(seen["image"])
    np.testing.assert_allclose(on_both, on_sharp, rtol=1e-12)
    assert not np.allclose(on_both, on_pre)
    got = per_spot.sort_values(["z_um", "y_um", "x_um"])["effect_size"].to_numpy(float)
    np.testing.assert_allclose(np.sort(got), np.sort(on_sharp), rtol=1e-12)


def test_neutral_options_give_the_default_tables():
    """The options branch mirrors the default body: with spot_radius_z_um equal to spot_radius_um and
    remove_hot_pixels and sharpen false, its tables equal the default call's."""
    pytest.importorskip("spotmax")
    img, labels, sp = _planted_foci(seed=9)
    params = dict(spot_radius_um=0.25, gauss_sigma_um=0.08, thresholding_method="threshold_triangle",
                  effect_size_min=1.0, merge_z_columns=True)
    ref = detect_spots_spotmax(img, labels, sp, **params)
    got = detect_spots_spotmax(img, labels, sp, spot_radius_z_um=0.25, remove_hot_pixels=False, sharpen=False,
                               **params)
    assert len(ref[0]) > 0
    pd.testing.assert_frame_equal(got[0], ref[0])
    pd.testing.assert_frame_equal(got[1], ref[1])


# ---- segmentation.nuclei.require_channel ----------------------------------------------------------------

def _spy_backend(calls):
    def spy(dna, spacing, *, model, diameter_um, device, cfg_section):
        calls.append(dna.shape)
        return N._classical(dna, spacing, diameter_um)
    return spy


def test_require_channel_fails_the_image_before_segmentation(tmp_path, patched_reader, registry):  # noqa: F811
    calls: list = []
    B.register_backend("spy", _spy_backend(calls))

    def cfg(require):
        c = load_config(ROOT / "config" / "config.yaml")
        c.set("segmentation.nuclei.method", "auto")
        c.set("segmentation.nuclei.backend", "spy")
        c.set("segmentation.nuclei.channel", "granule")            # no PGL-1 channel in a 3-channel image
        c.set("spots.enabled", False)
        if require is not None:
            c.set("segmentation.nuclei.require_channel", require)
        return c

    with pytest.raises(ValueError, match="require_channel is true .*segmentation role 'granule'"):
        pipeline.process_image(ND2, cfg(True), tmp_path / "required", prov=None)
    assert not calls                                                  # nothing was segmented
    assert not (tmp_path / "required" / "20251105_n2_nohs_HERM_001__done.json").exists()
    for require in (None, False):                                     # unset: today's fallback to channel 0
        res = pipeline.process_image(ND2, cfg(require), tmp_path / f"r{require}", prov=None)
        assert "segment:no_granule_channel_using_first_channel" in res["qc_flags"]
        assert res["stages"]["segment"]["status"] == "ran"
    assert len(calls) == 2
    # a role that resolves: the key changes nothing
    c = cfg(True)
    c.set("segmentation.nuclei.channel", "dna")
    res = pipeline.process_image(ND2, c, tmp_path / "dna", prov=None)
    assert not any(f.startswith("segment:no_") for f in res["qc_flags"]) and len(calls) == 3


# ---- qc.min_nuclei ---------------------------------------------------------------------------------------

def test_min_nuclei_moves_the_few_nuclei_threshold():
    base = dict(channel_flags=[], axis_flags=[], spots_found=True)
    assert "qc:few_nuclei_5" in qc.qc_flags(n_nuclei=5, **base)[1]           # the validated rule: below 10
    assert not any("few_nuclei" in f for f in qc.qc_flags(n_nuclei=10, **base)[1])
    assert qc.qc_flags(n_nuclei=5, **base) == qc.qc_flags(n_nuclei=5, min_nuclei=10, **base)
    assert not any("few_nuclei" in f for f in qc.qc_flags(n_nuclei=5, min_nuclei=3, **base)[1])
    assert not any("few_nuclei" in f for f in qc.qc_flags(n_nuclei=1, min_nuclei=1, **base)[1])
    assert "qc:few_nuclei_12" in qc.qc_flags(n_nuclei=12, min_nuclei=20, **base)[1]
    ok, flags = qc.qc_flags(n_nuclei=0, min_nuclei=1, **base)
    assert flags == ["qc:NO_NUCLEI"] and ok is False                   # no nuclei is still a hard fail


def test_process_image_passes_min_nuclei_only_when_set(tmp_path, patched_reader, monkeypatch):  # noqa: F811
    calls = []
    real = qc.qc_flags

    def spy(**kw):
        calls.append(kw)
        return real(**kw)

    monkeypatch.setattr(qc, "qc_flags", spy)
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg.set("segmentation.nuclei.method", "classical")
    cfg.set("spots.enabled", False)
    res = pipeline.process_image(ND2, cfg, tmp_path / "default", prov=None)
    assert "min_nuclei" not in calls[-1]                              # the default call is today's
    assert res["n_nuclei"] >= 10 and not any("few_nuclei" in f for f in res["qc_flags"])
    cfg.set("qc.min_nuclei", 40)
    res = pipeline.process_image(ND2, cfg, tmp_path / "forty", prov=None)
    assert calls[-1]["min_nuclei"] == 40
    assert f"qc:few_nuclei_{res['n_nuclei']}" in res["qc_flags"]
