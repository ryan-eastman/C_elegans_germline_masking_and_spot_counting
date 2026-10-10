"""The nucleus segmentation backend registry (germquant.segment.backends): the built-ins, a backend
registered in-process, entry-point discovery, the unknown-name error, and the wiring into
segment_nuclei and the pipeline's segmentation_method column. No test here calls the real Cellpose."""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import pytest
from test_end_to_end import patched_reader  # noqa: F401  (fixture: synthetic 3-channel gonad)

from germquant import device as D
from germquant import pipeline
from germquant.config import load_config
from germquant.segment import backends as B
from germquant.segment import nuclei as N

SPACING = (0.4, 0.2, 0.2)


def _dna():
    Z, Y, X = 8, 24, 24
    zz, yy, xx = np.indices((Z, Y, X)).astype(np.float32)
    dna = np.zeros((Z, Y, X), np.float32)
    for cx in (6, 17):
        r2 = (((zz - 4) * SPACING[0]) ** 2 + ((yy - 12) * SPACING[1]) ** 2
              + ((xx - cx) * SPACING[2]) ** 2)
        dna += np.exp(-r2 / (2 * 1.0 ** 2)) * 1000.0
    return dna


@pytest.fixture
def registry(monkeypatch):
    """A private copy of the registry for the test (restored afterwards), with entry-point discovery
    marked done so no plugin installed on the machine leaks into the assertions."""
    monkeypatch.setattr(B, "_REGISTRY", dict(B._REGISTRY))
    monkeypatch.setattr(B, "_SOURCES", dict(B._SOURCES))
    monkeypatch.setitem(B._STATE, "discovered", True)
    return B


def test_registry_resolves_the_builtins(registry):
    assert {"cellpose", "classical"} <= set(B.list_backends())
    assert B.get_backend("cellpose") is B.cellpose_backend
    assert B.get_backend("classical") is B.classical_backend
    assert B.get_backend(None) is B.cellpose_backend          # None / '' = the default
    assert B.get_backend("") is B.cellpose_backend
    assert B.backend_sources()["cellpose"] == "built-in"
    dna = _dna()
    labels = B.get_backend("classical")(dna, SPACING, model="cpsam", diameter_um=3.0, device=None,
                                        cfg_section={})
    assert labels.shape == dna.shape and labels.dtype == np.int32 and labels.max() >= 1


def test_classical_backend_is_the_existing_watershed(registry):
    """backend='classical' through segment_nuclei equals method='classical' voxel for voxel, and both
    equal the private functions they wrap."""
    dna = _dna()
    ref, m_ref = N.segment_nuclei(dna, SPACING, method="classical", min_volume_um3=0.1)
    got, m_got = N.segment_nuclei(dna, SPACING, method="auto", backend="classical", min_volume_um3=0.1)
    assert m_ref == "classical" and m_got == "classical"
    np.testing.assert_array_equal(got, ref)
    direct = N._drop_small(N._classical(dna, SPACING, 3.0), SPACING, 0.1)
    np.testing.assert_array_equal(ref, direct)


def test_fake_backend_gets_the_documented_arguments_and_names_the_method(registry):
    calls = []

    def fake(dna, spacing, *, model, diameter_um, device, cfg_section):
        calls.append({"dna": dna, "spacing": spacing, "model": model, "diameter_um": diameter_um,
                      "device": device, "cfg_section": cfg_section})
        out = np.zeros(dna.shape, np.int64)              # any integer dtype is accepted
        out[:, :12, :] = 1
        out[:, 12:, :] = 2
        return out

    B.register_backend("fake3d", fake)
    assert "fake3d" in B.list_backends() and B.backend_sources()["fake3d"] == "registered"
    dna = _dna()
    section = {"method": "auto", "backend": "fake3d", "fake3d": {"threshold": 0.5}}
    labels, used = N.segment_nuclei(dna, SPACING, method="auto", backend="fake3d",
                                    cellpose_model="models/x", diameter_um=2.5, device="cpu",
                                    cfg_section=section, min_volume_um3=0.0)
    assert used == "fake3d"
    assert labels.dtype == np.int32 and labels.shape == dna.shape and labels.max() == 2
    (c,) = calls
    assert c["dna"] is dna and c["spacing"] == SPACING
    assert c["model"] == "models/x" and c["diameter_um"] == 2.5 and c["device"] == "cpu"
    assert c["cfg_section"] == section and c["cfg_section"] is not section   # a copy: the config is safe
    # the min-volume filter is segment_nuclei's job, applied after the backend
    small, _ = N.segment_nuclei(dna, SPACING, method="auto", backend="fake3d", min_volume_um3=1e9)
    assert small.max() == 0


def test_pipeline_reports_the_backend_in_segmentation_method(tmp_path, registry, patched_reader):  # noqa: F811
    seen = {}

    def fake(dna, spacing, *, model, diameter_um, device, cfg_section):
        seen.update(model=model, diameter_um=diameter_um, device=device, cfg_section=cfg_section)
        return N._classical(dna, spacing, diameter_um)      # any instance labels will do

    B.register_backend("fake3d", fake)
    cfg = load_config("config/config.yaml")
    h0 = cfg.hash
    cfg.set("segmentation.nuclei.backend", "fake3d")       # an override, so it enters the config hash
    assert cfg.hash != h0
    cfg.set("spots.enabled", False)                         # keep the run cheap
    out = tmp_path / "r"
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, out, prov=None)

    summary = pd.read_csv(out / f"{res['image_id']}__image_summary.csv")
    assert summary["segmentation_method"].iloc[0] == "fake3d"
    assert res["stages"]["segment"]["status"] == "ran" and res["n_nuclei"] >= 1
    # the backend saw the config section and the resolved Cellpose arguments
    assert seen["cfg_section"]["backend"] == "fake3d" and seen["cfg_section"]["method"] == "auto"
    assert seen["diameter_um"] == 3.0 and seen["device"] is None
    assert "germline_nuclei_combined" in str(seen["model"])


def test_default_pipeline_path_still_says_classical(tmp_path, registry, patched_reader):  # noqa: F811
    """No backend key in the config (the default) + method classical: the column keeps saying
    'classical', as the goldens require."""
    cfg = load_config("config/config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    cfg.set("spots.enabled", False)
    out = tmp_path / "r"
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, out, prov=None)
    summary = pd.read_csv(out / f"{res['image_id']}__image_summary.csv")
    assert summary["segmentation_method"].iloc[0] == "classical"


def test_unknown_backend_raises_a_clear_error(registry):
    with pytest.raises(B.UnknownBackendError, match="unknown nucleus segmentation backend 'nope'") as ei:
        B.get_backend("nope")
    msg = str(ei.value)
    assert "cellpose" in msg and "classical" in msg and B.ENTRY_POINT_GROUP in msg
    assert isinstance(ei.value, LookupError)
    dna = _dna()
    for method in ("auto", "cellpose"):      # a misspelt name is never silently replaced by the fallback
        with pytest.raises(B.UnknownBackendError, match="nope"):
            N.segment_nuclei(dna, SPACING, method=method, backend="nope")
    _labels, used = N.segment_nuclei(dna, SPACING, method="classical", backend="nope", min_volume_um3=0.1)
    assert used == "classical"               # method=classical never consults the backend


def test_register_backend_validates_and_replaces(registry):
    def fn(dna, spacing, **kw):
        return np.zeros(dna.shape, np.int32)

    def other(dna, spacing, **kw):
        return np.ones(dna.shape, np.int32)

    with pytest.raises(TypeError, match="callable"):
        B.register_backend("notcallable", "x")
    with pytest.raises(ValueError, match="empty"):
        B.register_backend("  ", fn)
    with pytest.raises(TypeError, match="string"):
        B.register_backend(3, fn)
    B.register_backend("mine", fn)
    with pytest.raises(ValueError, match="already registered"):
        B.register_backend("mine", other)
    assert B.get_backend("mine") is fn
    B.register_backend("mine", other, replace=True)
    assert B.get_backend(" mine ") is other               # names are stripped
    B.unregister_backend("mine")
    assert "mine" not in B.list_backends()
    with pytest.raises(B.UnknownBackendError):
        B.unregister_backend("mine")
    # a test may stand in for a built-in; unregistering restores the built-in, never removes it
    B.register_backend("cellpose", other, replace=True)
    assert B.get_backend("cellpose") is other
    B.unregister_backend("cellpose")
    assert B.get_backend("cellpose") is B.cellpose_backend


class _EP:
    """Stands in for importlib.metadata.EntryPoint."""

    def __init__(self, name, value, obj=None, error=None):
        self.name, self.value, self._obj, self._error = name, value, obj, error

    def load(self):
        if self._error is not None:
            raise self._error
        return self._obj


def test_entry_points_are_discovered_and_failures_are_logged_not_raised(registry, monkeypatch, caplog):
    def good(dna, spacing, *, model, diameter_um, device, cfg_section):
        return np.zeros(dna.shape, np.int32)

    eps = [
        _EP("plugin_ok", "plugin_pkg.seg:segment", obj=good),
        _EP("plugin_broken", "plugin_broken:segment", error=ImportError("No module named 'plugin_broken'")),
        _EP("plugin_notcallable", "plugin_x:thing", obj=42),
        _EP("cellpose", "evil_pkg:replace_cellpose", obj=good),      # may not shadow a built-in
    ]
    monkeypatch.setattr(B, "_iter_entry_points", lambda: list(eps))
    monkeypatch.setitem(B._STATE, "discovered", False)
    with caplog.at_level(logging.WARNING, logger="germquant.segment.backends"):
        loaded = B.discover_entry_points()
    assert loaded == {"plugin_ok": "plugin_pkg.seg:segment"}
    assert B.get_backend("plugin_ok") is good
    assert B.backend_sources()["plugin_ok"] == "entry point plugin_pkg.seg:segment"
    assert B.get_backend("cellpose") is B.cellpose_backend
    assert "plugin_notcallable" not in B.list_backends() and "plugin_broken" not in B.list_backends()
    text = caplog.text
    assert "plugin_broken" in text and "ImportError" in text
    assert "plugin_notcallable" in text and "not callable" in text
    assert "evil_pkg" in text and "already registered" in text
    # the plugin runs through segment_nuclei like a built-in and names itself in the method
    _labels, used = N.segment_nuclei(_dna(), SPACING, method="cellpose", backend="plugin_ok")
    assert used == "plugin_ok"
    # discovery runs once per process; a forced rescan does not warn about its own earlier loads
    assert B.discover_entry_points() == {}
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="germquant.segment.backends"):
        assert B.discover_entry_points(force=True) == {}
    assert "plugin_ok" not in caplog.text

    # a damaged installation (reading the entry points itself fails) leaves the built-ins usable
    def boom():
        raise RuntimeError("bad metadata")

    monkeypatch.setattr(B, "_iter_entry_points", boom)
    with caplog.at_level(logging.WARNING, logger="germquant.segment.backends"):
        assert B.discover_entry_points(force=True) == {}
    assert "bad metadata" in caplog.text
    assert {"cellpose", "classical"} <= set(B.list_backends())


def test_auto_falls_back_to_classical_when_the_backend_fails(registry, caplog):
    def boom(dna, spacing, **kw):
        raise RuntimeError("no model for you")

    B.register_backend("boom", boom)
    dna = _dna()
    with caplog.at_level(logging.WARNING, logger="germquant.segment.nuclei"):
        labels, used = N.segment_nuclei(dna, SPACING, method="auto", backend="boom", min_volume_um3=0.1)
    assert used == "classical"
    assert "'boom'" in caplog.text and "classical watershed fallback" in caplog.text
    ref, _ = N.segment_nuclei(dna, SPACING, method="classical", min_volume_um3=0.1)
    np.testing.assert_array_equal(labels, ref)
    with pytest.raises(RuntimeError, match="no model for you"):     # method=cellpose: no fallback
        N.segment_nuclei(dna, SPACING, method="cellpose", backend="boom")


def test_backend_output_is_checked_against_the_contract(registry):
    def wrong_shape(dna, spacing, **kw):
        return np.zeros(dna.shape[1:], np.int32)

    def wrong_dtype(dna, spacing, **kw):
        return np.zeros(dna.shape, np.float32)

    B.register_backend("wrong_shape", wrong_shape)
    B.register_backend("wrong_dtype", wrong_dtype)
    dna = _dna()
    with pytest.raises(ValueError, match="'wrong_shape'.*shape"):
        N.segment_nuclei(dna, SPACING, method="cellpose", backend="wrong_shape")
    with pytest.raises(TypeError, match="'wrong_dtype'.*integer"):
        N.segment_nuclei(dna, SPACING, method="cellpose", backend="wrong_dtype")
    _, used = N.segment_nuclei(dna, SPACING, method="auto", backend="wrong_shape")   # auto: the fallback rules apply
    assert used == "classical"


def test_cellpose_backend_keeps_the_explicit_device_and_the_cpu_retry(registry, monkeypatch):
    seen = []

    def fake_cellpose(dna, spacing, model, diam, device_pref=None):
        seen.append((model, diam, device_pref))
        if device_pref != "cpu":
            raise RuntimeError("CUDA out of memory")
        return np.ones(dna.shape, np.int32)

    monkeypatch.setattr(N, "_cellpose", fake_cellpose)
    monkeypatch.setattr(D, "available", lambda: {"cuda": True, "mps": False})
    monkeypatch.delenv("GERMQUANT_DEVICE", raising=False)
    dna = _dna()
    out = B.get_backend("cellpose")(dna, SPACING, model="models/m", diameter_um=3.0, device="cuda",
                                    cfg_section={})
    assert out.dtype == np.int32 and out.max() == 1
    assert seen == [("models/m", 3.0, "cuda"), ("models/m", 3.0, "cpu")]
    # through segment_nuclei's default path the method is still reported as 'cellpose'
    seen.clear()
    labels, used = N.segment_nuclei(dna, SPACING, method="cellpose", cellpose_model="models/m",
                                    min_volume_um3=0.0)
    assert used == "cellpose" and labels.max() == 1 and [s[2] for s in seen] == [None, "cpu"]

    # already on the CPU: no retry, the error propagates (and method=auto then takes the fallback)
    seen.clear()
    monkeypatch.setattr(D, "available", lambda: {"cuda": False, "mps": False})

    def always_fails(dna, spacing, model, diam, device_pref=None):
        seen.append(device_pref)
        raise RuntimeError("boom")

    monkeypatch.setattr(N, "_cellpose", always_fails)
    with pytest.raises(RuntimeError, match="boom"):
        B.get_backend("cellpose")(dna, SPACING, model="cpsam", diameter_um=3.0, device=None, cfg_section={})
    assert seen == [None]
    labels, used = N.segment_nuclei(dna, SPACING, method="auto", min_volume_um3=0.1)
    assert used == "classical" and labels.max() >= 1
