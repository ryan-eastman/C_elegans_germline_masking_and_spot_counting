"""The model path guard of pipeline.resolve_model_path.

A configured Cellpose model that resolves to no file must stop the run: handed to Cellpose, a missing path is
taken for an unknown model name, logged as "pretrained model ... not found, using default model", and the
image is segmented with Cellpose's stock weights while the manifest records model_sha256 "absent". That
happened on 2026-10-23 for a study config outside the checkout run from a worktree without the models folder
(live_results/kate_movies/readiness.md, gap G1). The refusal names every place looked in and the fixes.
Cellpose's own names (cpsam, cpsam_v2) still pass through, every value that resolves is returned exactly as
before (the same value and type), and runs that never load the model (method classical, another backend) are
not refused. Nothing here needs model weights: the Cellpose call is replaced by a spy where a run is made.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest
from test_end_to_end import patched_reader  # noqa: F401  (fixture: the 3-channel fixed gonad)
from test_live_profile import ND2 as LIVE_ND2
from test_live_profile import patch_live_series
from test_segment_backends import registry  # noqa: F401  (fixture: a private backend registry)

from germquant import cli, pipeline, provenance
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


def _refusal(cfg) -> str:
    with pytest.raises(FileNotFoundError) as e:
        pipeline.resolve_model_path(cfg)
    return str(e.value)


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


@pytest.mark.parametrize("extra", ["    method: classical\n", "    backend: classical\n", "    backend: nuclei3d\n",
                                   "    backend: no_such_backend\n"])
def test_runs_that_never_load_the_cellpose_model_are_not_refused(tmp_path, no_models, extra):
    # (an unknown backend name is the segment stage's own error, UnknownBackendError, raised there)
    assert pipeline.resolve_model_path(load_config(_study(tmp_path / "study", extra=extra))) == REL


def test_a_stand_in_registered_under_the_cellpose_name_is_not_refused(tmp_path, no_models, registry):  # noqa: F811
    # tests/test_packs_run.py replaces the "cellpose" name with the watershed; only the built-in loads the model
    registry.register_backend("cellpose", registry.classical_backend, replace=True)
    assert pipeline.resolve_model_path(load_config(_study(tmp_path / "study"))) == REL


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


# ---- the refusal ---------------------------------------------------------------------------------------

def test_a_relative_model_path_found_nowhere_is_refused_with_the_places_looked_in(tmp_path, no_models):
    cfg = load_config(_study(tmp_path / "study"))
    msg = _refusal(cfg)
    assert f"segmentation.nuclei.cellpose_model {REL!r} is not a file" in msg
    assert str(cfg.base_dir / REL) in msg                          # next to the config
    assert str(Path.cwd() / REL) in msg                            # the working directory, where Cellpose looks
    assert "GERMQUANT_MODELS is not set" in msg
    assert "stock weights" in msg and "using default model" in msg  # what would have happened
    # how to fix it: the model store, the download, and how to ask for the stock weights on purpose
    assert "set GERMQUANT_MODELS to a folder that holds models/germline_nuclei_combined" in msg
    dest = R.models_dir() / "germline_nuclei_combined"            # the checkout's models/models
    assert f"germquant get-model germline_nuclei_combined` downloads it to {dest}" in msg
    assert "segmentation.nuclei.cellpose_model: cpsam" in msg
    # when this installation already holds the model, the message names the folder to point GERMQUANT_MODELS at
    dest.parent.mkdir(parents=True)
    dest.write_bytes(b"weights")
    msg = _refusal(cfg)
    assert f"holds models/germline_nuclei_combined ({dest.parent.parent} does)" in msg and "get-model" not in msg


def test_the_refusal_names_the_GERMQUANT_MODELS_candidates(tmp_path, no_models, monkeypatch):
    store = tmp_path / "empty_store"
    store.mkdir()
    monkeypatch.setenv("GERMQUANT_MODELS", str(store))
    msg = _refusal(load_config(_study(tmp_path / "study")))
    assert str(store / REL) in msg and str(store / "models" / "germline_nuclei_combined") in msg
    assert "GERMQUANT_MODELS is not set" not in msg
    assert str(store / "models" / "germline_nuclei_combined") in msg.split("get-model")[1]   # where it would download


def test_an_absolute_model_path_that_is_not_a_file_is_refused(tmp_path, no_models):
    missing = tmp_path / "nowhere" / "my_model"
    msg = _refusal(load_config(_study(tmp_path / "study", model=missing.as_posix())))
    assert "is not a file" in msg and str(missing) in msg
    assert "get-model" not in msg                                  # not a release or registered model name
    (tmp_path / "nowhere").mkdir()
    missing.mkdir()                                                # a folder is not a model file either
    assert "is not a file" in _refusal(load_config(_study(tmp_path / "study2", model=missing.as_posix())))
    canonical = R.models_dir() / "germline_nuclei_combined"       # the absolute path where get-model keeps it
    msg = _refusal(load_config(_study(tmp_path / "study3", model=canonical.as_posix())))
    assert f"at that path (`germquant get-model germline_nuclei_combined` downloads it to {canonical})" in msg


def test_an_inherited_value_anchored_to_a_folder_without_the_model_is_refused(tmp_path, no_models):
    # a base config with the shipped relative value, in a folder without models/ (a checkout without the models
    # junction), and a study profile elsewhere that inherits it: load_config makes the value absolute against the
    # base's folder, where GERMQUANT_MODELS is never consulted
    base = _study(tmp_path / "checkout")
    cfg = load_config(_study(tmp_path / "study", model=None, base=base))
    anchored = Path(cfg.get("segmentation.nuclei.cellpose_model"))
    assert anchored.is_absolute() and anchored == (tmp_path / "checkout" / REL).resolve()
    msg = _refusal(cfg)
    assert str(anchored) in msg
    assert f"inherited as {REL!r}" in msg and f"cellpose_model: {REL}" in msg


# ---- the run stops before any segmentation -------------------------------------------------------------

def test_process_image_stops_before_segmenting(tmp_path, no_models, patched_reader, monkeypatch):  # noqa: F811
    called = []
    monkeypatch.setattr(N, "_cellpose", lambda *a, **k: called.append(a))
    out = tmp_path / "run"
    with pytest.raises(FileNotFoundError, match="is not a file"):
        pipeline.process_image(ND2, load_config(_study(tmp_path / "study")), out, prov=None)
    assert not called
    assert not (out / "run_manifest.json").exists() and not list(out.glob("*__*"))


def test_batch_stops_before_its_first_image(tmp_path, no_models, monkeypatch):
    data = tmp_path / "data"
    data.mkdir()
    (data / ND2).write_bytes(b"")
    called = []
    monkeypatch.setattr(pipeline, "process_image", lambda *a, **k: called.append(a))
    with pytest.raises(FileNotFoundError, match="is not a file"):
        cli.main(["batch", str(data), "--config", str(_study(tmp_path / "study")), "--out", str(tmp_path / "out")])
    assert not called and not (tmp_path / "out").exists()


def test_a_time_series_stops_before_its_first_timepoint(tmp_path, no_models, monkeypatch):
    # the route of gap G1: `germquant run --t` with a study config that extends the live profile
    calls = patch_live_series(monkeypatch)
    segmented = []
    monkeypatch.setattr(N, "_cellpose", lambda *a, **k: segmented.append(a))     # never the real Cellpose
    study = _study(tmp_path / "study", base=ROOT / "config" / "profiles" / "live_syp2.yaml")
    with pytest.raises(FileNotFoundError, match="is not a file"):
        cli.main(["run", LIVE_ND2, "--config", str(study), "--out", str(tmp_path / "series"), "--t", "all"])
    assert calls == [] and segmented == []                          # no timepoint was read or segmented
