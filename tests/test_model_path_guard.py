"""How pipeline.resolve_model_path resolves segmentation.nuclei.cellpose_model, pinned before the model path guard.

Cellpose's own names (cpsam, cpsam_v2) pass through to Cellpose, every value that resolves is returned exactly
as it is today (the same value and type: the joined Path under $GERMQUANT_MODELS or next to the config, the
config's own string for an absolute file or for a file found from the working directory), and runs that never
load the model (method classical, another backend) get the value back unresolved. Nothing here needs model
weights: the Cellpose call is replaced by a spy where a run is made.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
from test_end_to_end import patched_reader  # noqa: F401  (fixture: the 3-channel fixed gonad)

from germquant import pipeline, provenance
from germquant import models_registry as R
from germquant.config import load_config
from germquant.segment import nuclei as N

ROOT = Path(__file__).resolve().parents[1]
REL = "models/models/germline_nuclei_combined"            # the value config/config.yaml ships
ND2 = "20251105_n2_nohs_HERM_001.nd2"


@pytest.fixture
def no_models(tmp_path, monkeypatch):
    """No model anywhere: GERMQUANT_MODELS unset, an empty working directory, and an empty stand-in for the
    checkout's models/models (models_registry.models_dir), so no test depends on what this machine holds."""
    monkeypatch.delenv("GERMQUANT_MODELS", raising=False)
    real_models_dir = R.models_dir
    checkout = tmp_path / "checkout_root"
    monkeypatch.setattr(R, "models_dir", lambda base_dir=None: real_models_dir(checkout if base_dir is None else base_dir))
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    return cwd


def _study(folder: Path, model: str | None = REL, extra: str = "",
           base: Path = ROOT / "config" / "config.yaml") -> Path:
    """A study config outside the checkout that extends ``base`` and names ``model`` itself (None: inherits it),
    with ``extra`` lines under segmentation.nuclei."""
    folder.mkdir(parents=True, exist_ok=True)
    text = f"extends: {base.as_posix()}\n"
    if model is not None or extra:
        text += "segmentation:\n  nuclei:\n" + (f"    cellpose_model: {model}\n" if model is not None else "") + extra
    p = folder / "study.yaml"
    p.write_text(text, encoding="utf-8")
    return p


# ---- what is not refused, pinned before the guard ------------------------------------------------

@pytest.mark.parametrize("name", ["cpsam", "cpsam_v2"])
def test_cellpose_stock_weights_named_explicitly_pass_through(tmp_path, no_models, name):
    assert pipeline.resolve_model_path(load_config(_study(tmp_path / "study", model=name))) == name


def test_the_stock_name_reaches_cellpose_unchanged(tmp_path, no_models, patched_reader, monkeypatch):  # noqa: F811
    seen = []

    def spy(dna, spacing, model_name, diameter_um, device_pref=None, diameter_px=None):
        seen.append(model_name)
        return N._classical(dna, spacing, diameter_um)

    monkeypatch.setattr(N, "_cellpose", spy)
    cfg = load_config(_study(tmp_path / "study", model="cpsam_v2"))
    cfg._data["spots"]["enabled"] = False
    res = pipeline.process_image(ND2, cfg, tmp_path / "run", prov=None)
    assert seen == ["cpsam_v2"]
    assert res["tables"]["image_summary"]["segmentation_method"].iloc[0] == "cellpose"
    man = json.loads((tmp_path / "run" / "run_manifest.json").read_text(encoding="utf-8"))
    assert man["model_sha256"] == "absent"                         # a name, not a file: as it always recorded


@pytest.mark.parametrize("extra", ["    method: classical\n", "    backend: classical\n", "    backend: nuclei3d\n"])
def test_runs_that_never_load_the_cellpose_model_are_not_refused(tmp_path, no_models, extra):
    assert pipeline.resolve_model_path(load_config(_study(tmp_path / "study", extra=extra))) == REL


def test_values_that_resolve_are_returned_exactly_as_before(tmp_path, no_models, monkeypatch):
    # next to the config (base_dir): the joined Path
    model = tmp_path / "study" / REL
    model.parent.mkdir(parents=True)
    model.write_bytes(b"weights")
    cfg = load_config(_study(tmp_path / "study"))
    assert pipeline.resolve_model_path(cfg) == cfg.base_dir / REL
    # under $GERMQUANT_MODELS, with the leading models/ dropped: that Path, ahead of base_dir as before
    store = tmp_path / "store"
    (store / "models").mkdir(parents=True)
    (store / "models" / "germline_nuclei_combined").write_bytes(b"w")
    monkeypatch.setenv("GERMQUANT_MODELS", str(store))
    assert pipeline.resolve_model_path(cfg) == store / "models" / "germline_nuclei_combined"
    monkeypatch.delenv("GERMQUANT_MODELS")
    # an absolute path to a file: the config's own string, untouched
    assert pipeline.resolve_model_path(load_config(_study(tmp_path / "abs", model=model.as_posix()))) == model.as_posix()
    # a relative path that exists only from the working directory (Cellpose's own lookup): the string itself,
    # and the run records the file's digest as it always did
    (no_models / REL).parent.mkdir(parents=True)
    (no_models / REL).write_bytes(b"w")
    cfg_cwd = load_config(_study(tmp_path / "cwd_study"))
    assert pipeline.resolve_model_path(cfg_cwd) == REL
    assert provenance.file_sha256(pipeline.resolve_model_path(cfg_cwd)) == provenance.file_sha256(no_models / REL)
