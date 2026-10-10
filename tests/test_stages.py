"""The stage registry (germquant.stages) and its wiring into process_image / the CLI."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
from test_end_to_end import patched_reader  # noqa: F401  (fixture: synthetic 3-channel gonad)

from germquant import pipeline, schema, stages
from germquant.config import load_config
from germquant.stages import STAGES, Outcome, apply_cli_switches, run_stage, stage_enabled


def test_registry_is_ordered_and_unique():
    names = [s.name for s in STAGES]
    assert names == ["read", "acquisition", "segment", "measure", "signal_qc", "germline", "envelope", "audit",
                     "axis", "staging", "spots", "sc_trace", "granule", "coloc", "partition", "granule_tail",
                     "qc", "render", "write"]
    assert len(set(names)) == len(names)
    keys = [s.config_key for s in STAGES if s.config_key]
    assert len(set(keys)) == len(keys)
    for s in STAGES:
        if s.cli_switchable:
            assert s.config_key, f"{s.name}: a CLI switch needs a config key"


def test_stage_enabled_reads_code_defaults_without_touching_the_config(tmp_path):
    cfg = load_config("config/config.yaml")
    h = cfg.hash
    assert stage_enabled(cfg, "axis") is True          # key absent from config.yaml -> code default
    assert stage_enabled(cfg, "granule") is True       # follows coloc.enabled when absent
    assert stage_enabled(cfg, "spots") is True
    assert cfg.hash == h                                # reading switches never changes the hash
    cfg.set("coloc.enabled", False)
    assert stage_enabled(cfg, "granule") is False      # legacy --no-coloc semantics: no granules either
    cfg.set("granule.enabled", True)
    assert stage_enabled(cfg, "granule") is True       # but a profile can keep granules with coloc off


def test_cli_switches_are_overrides_that_enter_the_hash():
    class A:
        no_spots = True
        no_coloc = False
        no_germline = False
        no_axis = False
        no_granule = False

    cfg = load_config("config/config.yaml")
    h0 = cfg.hash
    msgs = apply_cli_switches(cfg, A())
    assert cfg.get("spots.enabled") is False and cfg.hash != h0
    assert msgs == ["segmentation only: skipping spot detection (--no-spots)"]
    assert cfg.get("coloc.enabled") is True            # untouched


def test_run_stage_records_and_flags():
    flags, outcomes = [], {}

    def ok():
        flags.append("x:done")
        return 42

    def boom():
        raise ValueError("nope")

    assert run_stage("x", ok, flags=flags, outcomes=outcomes) == 42
    assert outcomes["x"].status == "ran" and outcomes["x"].flags == ["x:done"]
    assert run_stage("y", boom, flags=flags, outcomes=outcomes, fail_value="fb") == "fb"
    assert outcomes["y"].status == "failed" and "ValueError" in outcomes["y"].error
    assert flags[-1] == "y:FAILED_ValueError"
    assert run_stage("z", ok, flags=flags, outcomes=outcomes, enabled=False, skip_reason="off") is None
    assert outcomes["z"].status == "skipped" and outcomes["z"].reason == "off"
    with pytest.raises(ValueError):
        run_stage("w", boom, flags=flags, outcomes=outcomes, fatal=True)
    assert outcomes["w"].status == "failed"
    assert isinstance(Outcome("ran").as_dict(), dict)


def test_process_image_writes_the_stage_record(tmp_path, patched_reader):  # noqa: F811
    cfg = load_config("config/config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / "r", prov=None)
    rec = res["stages"]
    assert rec["segment"]["status"] == "ran" and rec["germline"]["status"] == "ran"
    assert rec["granule"]["status"] == "skipped" and "granule" in rec["granule"]["reason"]   # 3-channel image
    assert rec["coloc"]["status"] == "skipped"
    f = next((tmp_path / "r").rglob("*__stages.json"))
    on_disk = json.loads(f.read_text(encoding="utf-8"))
    assert on_disk["stages"].keys() == rec.keys()


def test_failed_granule_stage_flags_and_continues(tmp_path, patched_reader, monkeypatch):  # noqa: F811
    """A granule failure must flag, skip coloc, and still write every table (batch keeps going)."""
    cfg = load_config("config/config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"

    def boom(*a, **k):
        raise RuntimeError("granules exploded")

    monkeypatch.setattr(pipeline, "_run_granule", boom)
    _give_granule_role(cfg, monkeypatch)
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / "r", prov=None)
    assert res["stages"]["granule"]["status"] == "failed"
    assert res["stages"]["coloc"]["status"] == "skipped"
    assert res["stages"]["coloc"]["reason"] == "granule stage failed"
    assert any(f.startswith("granule:FAILED_RuntimeError") for f in res["qc_flags"])
    assert set(res["tables"]) == set(schema.TABLES)
    assert res["tables"]["granules"].empty and res["tables"]["coloc"].empty
    assert set(res["stages"]) == {s.name for s in STAGES}       # every registered stage is recorded


def _give_granule_role(cfg, monkeypatch):
    """Make the 3-channel synthetic image resolve a `granule` role (aliased to the RAD-51 channel) so
    the granule / coloc stages are attempted."""
    orig_resolve = cfg.channel_map.resolve

    def resolve(names):
        r, fl = orig_resolve(names)
        r = dict(r)
        r["granule"] = r.get("foci")
        return r, fl

    monkeypatch.setattr(cfg.channel_map, "resolve", resolve)
    assert resolve(["405", "477", "545"])[0]["granule"] is not None


def test_failed_coloc_stage_keeps_granule_products(tmp_path, patched_reader, monkeypatch):  # noqa: F811
    """Contract (stages.py docstring): a coloc failure keeps the granule outputs, flags coloc:FAILED,
    and adds no coloc fields to image_summary."""
    cfg = load_config("config/config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    _give_granule_role(cfg, monkeypatch)

    def boom(*a, **k):
        raise RuntimeError("coloc exploded")

    monkeypatch.setattr(pipeline, "_run_coloc", boom)
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / "r", prov=None)
    assert res["stages"]["granule"]["status"] == "ran"
    assert res["stages"]["coloc"]["status"] == "failed"
    assert any(f.startswith("coloc:FAILED_RuntimeError") for f in res["qc_flags"])
    assert not any(f.startswith("coloc:granules_n") for f in res["qc_flags"])
    assert "n_granules" in res["tables"]["nuclei"].columns
    assert res["tables"]["coloc"].empty
    isum = res["tables"]["image_summary"]
    assert "n_granules" not in isum.columns or isum["n_granules"].isna().all()
    assert next((tmp_path / "r").rglob("*__stages.json")).exists()


def test_missing_roles():
    assert stages.missing_roles({"dna": 0, "central_element": 1, "foci": 2}, "granule") == ["granule"]
    assert stages.missing_roles({"dna": 0, "central_element": 1, "granule": 3}, "granule") == []
    assert np.isnan(float("nan"))  # keep numpy imported for future numeric checks
    assert isinstance(pd.DataFrame(), pd.DataFrame)
