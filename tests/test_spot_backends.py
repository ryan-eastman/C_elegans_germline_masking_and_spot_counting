"""Pluggable spot detection backends (germquant.spots.backends):
- the registry resolves the built-ins, rejects duplicates and names an unknown backend clearly;
- detect_spots dispatches on `backend`, hands the backend the spots section and checks its tables;
- entry-point discovery is never fatal and can not replace a registered backend;
- the spotmax wrapper is the reference implementation, byte for byte;
- the experimental `log` baseline recovers planted Gaussian foci inside synthetic nuclei and assigns
  them to the right nuclei;
- the pipeline reads spots.backend with the code default and records a non-default backend.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_end_to_end import _synthetic_stack, patched_reader  # noqa: F401  (fixture: synthetic 3-channel gonad)

from germquant import pipeline
from germquant.config import load_config
from germquant.spots import backends, detect_spots, get_backend, list_backends, register_backend
from germquant.spots.backends import PER_NUC_COLS, PER_SPOT_COLS
from germquant.spots.detect import detect_spots_spotmax
from germquant.spots.log_detect import detect_spots_log

ROOT = Path(__file__).resolve().parents[1]
SP = (0.2, 0.108, 0.108)          # dz, dy, dx: the real anisotropy
SHAPE = (36, 128, 128)            # 7.2 x 13.8 x 13.8 um


def _tables(marker, detector, n_nuc=2):
    per_spot = pd.DataFrame({
        "spot_id": [0], "nucleus_id": [1], "marker": [marker], "z_um": [0.2], "y_um": [0.3], "x_um": [0.4],
        "effect_size": [5.0], "spot_mean_intensity": [100.0], "detector": [detector],
    }, columns=PER_SPOT_COLS)
    per_nuc = pd.DataFrame({
        "nucleus_id": list(range(1, n_nuc + 1)), "marker": marker, "n_spots": [1] + [0] * (n_nuc - 1),
        "mean_spot_intensity": [100.0] + [np.nan] * (n_nuc - 1), "detector": detector,
    }, columns=PER_NUC_COLS)
    return per_spot, per_nuc


def _two_nuclei():
    labels = np.zeros((4, 8, 8), np.int32)
    labels[1:3, 1:4, 1:4] = 1
    labels[1:3, 5:7, 5:7] = 2
    return np.zeros(labels.shape, np.float32), labels


# ---- registry ----

def test_registry_resolves_builtins():
    names = list_backends()
    assert "spotmax" in names and "log" in names
    assert backends.DEFAULT_BACKEND == "spotmax"
    assert get_backend("spotmax") is backends._spotmax_backend
    assert get_backend("log") is backends._log_backend
    assert list(PER_SPOT_COLS) == ["spot_id", "nucleus_id", "marker", "z_um", "y_um", "x_um",
                                   "effect_size", "spot_mean_intensity", "detector"]


def test_unknown_backend_is_a_clear_error():
    with pytest.raises(ValueError, match="unknown spot backend 'nope'; known backends: .*log.*spotmax"):
        get_backend("nope")
    img, labels = _two_nuclei()
    with pytest.raises(ValueError, match="unknown spot backend 'nope'"):
        detect_spots(img, labels, SP, backend="nope")
    with pytest.raises(TypeError, match="must be a string"):
        get_backend(None)


def test_register_backend_rejects_duplicates_and_bad_names():
    def fn(image, labels, spacing, *, cfg_section, marker):
        return _tables(marker, "dup_test")

    register_backend("dup_test", fn)
    try:
        with pytest.raises(ValueError, match="already registered"):
            register_backend("dup_test", fn)
        register_backend("dup_test", fn, replace=True)              # explicit override is allowed
        assert get_backend("dup_test") is fn
    finally:
        backends._REGISTRY.pop("dup_test", None)
    for bad in ("", "with space", "a;b", None, 3):
        with pytest.raises(ValueError, match="spot backend name"):
            register_backend(bad, fn)
    with pytest.raises(TypeError, match="must be callable"):
        register_backend("not_callable_test", object())
    assert "not_callable_test" not in list_backends()

    @register_backend("deco_test")
    def deco(image, labels, spacing, *, cfg_section, marker):
        return _tables(marker, "deco_test")

    try:
        assert get_backend("deco_test") is deco                     # decorator form registers and returns fn
    finally:
        backends._REGISTRY.pop("deco_test", None)


# ---- dispatch ----

def test_fake_backend_is_dispatched_and_named():
    seen = {}

    def fake(image, labels, spacing, *, cfg_section, marker):
        seen.update(shape=image.shape, spacing=tuple(spacing), marker=marker, section=dict(cfg_section))
        return _tables(marker, "fake_test")

    register_backend("fake_test", fake)
    try:
        img, labels = _two_nuclei()
        ps, pn = detect_spots(img, labels, SP, marker="COSA-1", backend="fake_test",
                              cfg_section={"log": {"diameter_um": 0.4}, "custom": 1, "spot_radius_um": 9.9})
        assert (ps["detector"] == "fake_test").all() and (pn["detector"] == "fake_test").all()
        assert list(ps.columns) == PER_SPOT_COLS and list(pn.columns) == PER_NUC_COLS
        assert seen["marker"] == "COSA-1" and seen["shape"] == (4, 8, 8) and seen["spacing"] == SP
        # the backend sees the raw section plus the RESOLVED SpotMAX parameters (explicit arguments win)
        assert seen["section"]["spot_radius_um"] == 0.3 and seen["section"]["effect_size_min"] == 0.0
        assert seen["section"]["custom"] == 1 and seen["section"]["log"] == {"diameter_um": 0.4}
    finally:
        backends._REGISTRY.pop("fake_test", None)


def test_backend_result_contract_is_checked(caplog):
    img, labels = _two_nuclei()

    def missing(image, labels, spacing, *, cfg_section, marker):
        return pd.DataFrame(columns=["spot_id", "nucleus_id"]), _tables(marker, "missing_test")[1]

    def extra(image, labels, spacing, *, cfg_section, marker):
        ps, pn = _tables(marker, "extra_test")
        return ps.assign(size_um=0.3)[["size_um"] + PER_SPOT_COLS], pn

    def misnamed(image, labels, spacing, *, cfg_section, marker):
        return _tables(marker, "somebody_else")

    def not_a_pair(image, labels, spacing, *, cfg_section, marker):
        return _tables(marker, "pair_test")[0]

    for name, fn in (("missing_test", missing), ("extra_test", extra), ("misnamed_test", misnamed),
                     ("pair_test", not_a_pair)):
        register_backend(name, fn)
    try:
        with pytest.raises(ValueError, match="'missing_test' per-spot table lacks the columns"):
            detect_spots(img, labels, SP, backend="missing_test")
        ps, _ = detect_spots(img, labels, SP, backend="extra_test")
        assert list(ps.columns) == PER_SPOT_COLS + ["size_um"]     # declared columns first, extras kept after
        with caplog.at_level(logging.WARNING, logger="germquant.spots.backends"):
            ps, _ = detect_spots(img, labels, SP, backend="misnamed_test")
        assert "the contract is the backend name" in caplog.text and (ps["detector"] == "somebody_else").all()
        with pytest.raises(TypeError, match="must return \\(per_spot, per_nucleus\\)"):
            detect_spots(img, labels, SP, backend="pair_test")
    finally:
        for name in ("missing_test", "extra_test", "misnamed_test", "pair_test"):
            backends._REGISTRY.pop(name, None)


class _EP:
    """A stand-in for importlib.metadata.EntryPoint."""

    def __init__(self, name, value, fn=None, exc=None):
        self.name, self.value, self._fn, self._exc = name, value, fn, exc

    def load(self):
        if self._exc is not None:
            raise self._exc
        return self._fn


def test_entry_point_discovery_is_never_fatal(monkeypatch, caplog):
    def good(image, labels, spacing, *, cfg_section, marker):
        return _tables(marker, "ep_good")

    eps = [_EP("ep_good", "pkg.mod:good", fn=good),
           _EP("ep_broken", "pkg.mod:broken", exc=ImportError("no module named pkg")),
           _EP("spotmax", "evil.pkg:replace", fn=good),              # must not replace the built-in
           _EP("ep_notcallable", "pkg.mod:obj", fn=object()),
           _EP("bad name", "pkg.mod:x", fn=good)]
    monkeypatch.setattr(backends, "_iter_entry_points", lambda: eps)
    try:
        with caplog.at_level(logging.WARNING, logger="germquant.spots.backends"):
            backends.discover_backends(force=True)
        assert get_backend("ep_good") is good
        assert get_backend("spotmax") is backends._spotmax_backend
        assert "ep_broken" not in list_backends() and "ep_notcallable" not in list_backends()
        assert "failed to load" in caplog.text and "would replace" in caplog.text
        assert "not callable" in caplog.text and "unusable name" in caplog.text
    finally:
        backends._REGISTRY.pop("ep_good", None)

    def broken_metadata():
        raise RuntimeError("corrupt dist-info")

    monkeypatch.setattr(backends, "_iter_entry_points", broken_metadata)
    with caplog.at_level(logging.WARNING, logger="germquant.spots.backends"):
        backends.discover_backends(force=True)                     # discovery itself failing is logged
    assert "discovery failed" in caplog.text and "spotmax" in list_backends()


# ---- spotmax wrapper ----

def test_spotmax_wrapper_matches_the_reference_implementation():
    pytest.importorskip("spotmax")
    from skimage.measure import label as cc_label

    stack = _synthetic_stack()
    labels = cc_label(stack.data[0] > 1000).astype(np.int32)       # the 12 synthetic nuclei (+ junk blobs)
    rad = stack.data[2]
    params = {"thresholding_method": "threshold_triangle", "effect_size_min": 3.0}   # the pipeline's values
    ref = detect_spots_spotmax(rad, labels, stack.spacing, **params)
    via_default = detect_spots(rad, labels, stack.spacing, **params)               # backend defaults to spotmax
    via_name = detect_spots(rad, labels, stack.spacing, backend="spotmax", cfg_section={"log": {}}, **params)
    for got in (via_default, via_name):
        pd.testing.assert_frame_equal(got[0], ref[0])
        pd.testing.assert_frame_equal(got[1], ref[1])
    assert len(ref[0]) > 0 and (ref[0]["detector"] == "spotmax").all() and (ref[1]["detector"] == "spotmax").all()


# ---- the experimental LoG baseline ----

def _planted_stack(seed=7, n_per_nucleus=5, sigma_um=0.17):
    """Four spherical nuclei (radius 1.5 um, nucleoplasm 120 over a 30 background) with Gaussian foci
    (FWHM 0.4 um, amplitude 400 to 800) planted at known positions, Gaussian noise sd 10, uint16."""
    rng = np.random.default_rng(seed)
    sp = np.array(SP)
    zz, yy, xx = np.indices(SHAPE).astype(np.float32)
    zum, yum, xum = zz * sp[0], yy * sp[1], xx * sp[2]
    labels = np.zeros(SHAPE, np.int32)
    img = np.full(SHAPE, 30.0, np.float32)
    centers = [(3.6, 3.6, 3.6), (3.6, 3.6, 10.2), (3.6, 10.2, 3.6), (3.6, 10.2, 10.2)]
    planted = []
    for k, c in enumerate(centers, start=1):
        c = np.array(c)
        inside = ((zum - c[0]) ** 2 + (yum - c[1]) ** 2 + (xum - c[2]) ** 2) <= 1.5 ** 2
        labels[inside] = k
        img[inside] = 120.0
        pts: list[np.ndarray] = []
        while len(pts) < n_per_nucleus:                      # inside 0.9 um of the centre, 0.7 um apart
            v = rng.uniform(-0.9, 0.9, size=3)
            if np.linalg.norm(v) > 0.9 or any(np.linalg.norm(c + v - q) < 0.7 for q in pts):
                continue
            pts.append(c + v)
        for p in pts:
            amp = rng.uniform(400, 800)
            d2 = (zum - p[0]) ** 2 + (yum - p[1]) ** 2 + (xum - p[2]) ** 2
            img += amp * np.exp(-d2 / (2 * sigma_um ** 2))
            planted.append((k, *p))
    img += rng.normal(0, 10, size=SHAPE).astype(np.float32)
    img = np.clip(np.round(img), 0, 65535).astype(np.uint16)
    return img, labels, pd.DataFrame(planted, columns=["nucleus_id", "z_um", "y_um", "x_um"])


def _match(detected, planted, tol_um=0.3):
    """Greedy nearest matching: each planted focus takes the closest unused detection within tol_um."""
    det = detected[["z_um", "y_um", "x_um"]].to_numpy(dtype=float)
    used = np.zeros(len(det), bool)
    pairs = []
    for i, p in enumerate(planted[["z_um", "y_um", "x_um"]].to_numpy(dtype=float)):
        if not len(det):
            break
        d = np.linalg.norm(det - p, axis=1)
        d[used] = np.inf
        j = int(np.argmin(d))
        if d[j] <= tol_um:
            used[j] = True
            pairs.append((i, j, float(d[j])))
    return pairs


def test_log_backend_recovers_planted_foci():
    img, labels, planted = _planted_stack()
    per_spot, per_nuc = detect_spots(img, labels, SP, marker="RAD-51", backend="log")   # code defaults
    assert list(per_spot.columns) == PER_SPOT_COLS and list(per_nuc.columns) == PER_NUC_COLS
    assert (per_spot["detector"] == "log").all() and (per_nuc["detector"] == "log").all()
    assert (per_spot["marker"] == "RAD-51").all() and list(per_spot["spot_id"]) == list(range(len(per_spot)))

    pairs = _match(per_spot, planted)
    recall = len(pairs) / len(planted)
    extras = len(per_spot) - len(pairs)
    assert recall >= 0.9, f"recovered {len(pairs)} of {len(planted)} planted foci"
    assert extras <= 0.1 * len(planted), f"{extras} extra spots for {len(planted)} planted"
    # every recovered focus is assigned to the nucleus it was planted in, and lands close to it
    for i, j, _ in pairs:
        assert int(per_spot.at[j, "nucleus_id"]) == int(planted.at[i, "nucleus_id"])
    assert np.mean([d for _, _, d in pairs]) < 0.1        # sub-voxel refinement: well under one xy voxel
    # the per-nucleus table has one row per label and reconciles with the per-spot table
    assert list(per_nuc["nucleus_id"]) == [1, 2, 3, 4]
    counts = per_spot.groupby("nucleus_id").size()
    assert all(int(per_nuc.set_index("nucleus_id").at[k, "n_spots"]) == int(counts.get(k, 0)) for k in (1, 2, 3, 4))
    assert per_nuc["n_spots"].sum() == len(per_spot)
    # the spot-level readouts are the documented ones: bright core, strong effect over the annulus
    assert (per_spot["spot_mean_intensity"] > 120).all() and (per_spot["effect_size"] > 3).all()


def test_log_backend_reads_its_keys_from_the_spots_section(caplog):
    img, labels, planted = _planted_stack()
    section = {"log": {"diameter_um": 0.4, "min_snr": 2.5, "typo_key": 1}, "merge_z_columns": True}
    with caplog.at_level(logging.WARNING, logger="germquant.spots.log_detect"):
        per_spot, _ = detect_spots(img, labels, SP, backend="log", cfg_section=section)
    assert "unknown keys ['typo_key']" in caplog.text
    assert len(_match(per_spot, planted)) >= 0.9 * len(planted)
    # a gate nothing passes gives one zero row per nucleus, not an empty per-nucleus table
    ps, pn = detect_spots(img, labels, SP, backend="log", cfg_section={"log": {"min_snr": 1e9}})
    assert ps.empty and list(ps.columns) == PER_SPOT_COLS
    assert list(pn["nucleus_id"]) == [1, 2, 3, 4] and (pn["n_spots"] == 0).all() and (pn["detector"] == "log").all()


def test_log_backend_edge_cases_and_validation():
    img, labels, _ = _planted_stack(n_per_nucleus=0)            # nuclei with noise only: no foci
    ps, pn = detect_spots_log(img, labels, SP)
    assert ps.empty and list(pn["nucleus_id"]) == [1, 2, 3, 4] and (pn["n_spots"] == 0).all()
    ps, pn = detect_spots_log(img, np.zeros_like(labels), SP)    # no labels at all
    assert ps.empty and pn.empty and list(pn.columns) == PER_NUC_COLS
    with pytest.raises(ValueError, match="same shape"):
        detect_spots_log(img, labels[:, :, :-1], SP)
    with pytest.raises(ValueError, match="diameter_um must be positive"):
        detect_spots_log(img, labels, SP, diameter_um=0)
    with pytest.raises(ValueError, match="annulus_outer_um must exceed"):
        detect_spots_log(img, labels, SP, annulus_inner_um=1.0, annulus_outer_um=0.5)
    with pytest.raises(TypeError, match="spots.log must be a mapping"):
        detect_spots(img, labels, SP, backend="log", cfg_section={"log": 5})


# ---- pipeline wiring ----

def test_pipeline_dispatches_spots_backend_from_config(tmp_path, patched_reader):  # noqa: F811
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    assert cfg.get("spots.backend") is None                       # the config files do not name a backend
    cfg.set("spots.backend", "log")
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / "r", prov=None)
    assert res["stages"]["spots"]["status"] == "ran", res["stages"]["spots"]
    spots = res["tables"]["spots"]
    assert len(spots) > 0 and (spots["detector"] == "log").all()
    assert any(f.startswith("spots:log_n=") for f in res["qc_flags"])
    assert "spots:unvalidated_backend=log" in res["qc_flags"]
    assert not any(f.startswith("spots:spotmax") for f in res["qc_flags"])
    nuclei = res["tables"]["nuclei"]
    assert int(nuclei["n_spots"].sum()) == len(spots)
    written = pd.read_csv(tmp_path / "r" / f"{res['image_id']}__spots.csv")
    assert (written["detector"] == "log").all() and list(written.columns[:9]) == PER_SPOT_COLS


def test_pipeline_unknown_backend_fails_the_stage_not_the_run(tmp_path, patched_reader):  # noqa: F811
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    cfg.set("spots.backend", "nope")
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / "r", prov=None)
    assert res["stages"]["spots"]["status"] == "failed" and "unknown spot backend 'nope'" in res["stages"]["spots"]["error"]
    assert "spots:FAILED_ValueError" in res["qc_flags"]
    assert res["tables"]["spots"].empty and len(res["tables"]["nuclei"]) > 0
