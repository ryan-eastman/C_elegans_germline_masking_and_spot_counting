"""germquant.api: the re-exports are the same objects as their sources, and run() is the path
`germquant run` takes (the CLI calls it), so both write the same files."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from test_end_to_end import patched_reader  # noqa: F401  (fixture: synthetic 3-channel gonad)

from germquant import api, batch, cli, models_registry, pipeline, provenance, stages
from germquant.config import Config, load_config
from germquant.io import read_image, read_image_metadata, read_stack
from germquant.segment import segment_nuclei
from germquant.spots import detect_spots

ROOT = Path(__file__).resolve().parents[1]
NAME = "20251105_n2_nohs_HERM_001.nd2"
IMAGE_ID = "20251105_n2_nohs_HERM_001"


def test_reexports_are_the_same_objects():
    assert api.process_image is pipeline.process_image
    assert api.load_config is load_config
    assert api.read_image is read_image and api.read_stack is read_stack
    assert api.read_image_metadata is read_image_metadata
    assert api.segment_nuclei is segment_nuclei and api.detect_spots is detect_spots
    assert api.STAGES is stages.STAGES and api.stage_enabled is stages.stage_enabled
    assert api.cli_switches is stages.cli_switches and api.apply_cli_switches is stages.apply_cli_switches
    assert api.resolve_model_path is pipeline.resolve_model_path and api.spot_instances is pipeline.spot_instances
    assert api.get_model is models_registry.get_model
    assert api.read_done_marker is provenance.read_done_marker
    assert api.discover_files is batch.discover_files and api.collect is batch.collect
    assert api.is_done is batch.is_done and api.load_exclusions is batch.load_exclusions
    missing = [n for n in api.__all__ if not hasattr(api, n)]
    assert not missing, missing
    assert callable(api.run) and "stability promise" in api.__doc__.lower()
    assert "run" in api.__all__ and "read_image" in api.__all__


def _classical_config(tmp_path) -> Path:
    """A profile outside the repo that forces the classical segmenter (no GPU, no Cellpose in tests)."""
    p = tmp_path / "study.yaml"
    p.write_text(f"extends: {(ROOT / 'config' / 'config.yaml').as_posix()}\n"
                 f"segmentation:\n  nuclei:\n    method: classical\n")
    return p


def _tables(out: Path) -> dict[str, pd.DataFrame]:
    return {f.name: pd.read_csv(f, low_memory=False, float_precision="round_trip").drop(columns=["run_timestamp"])
            for f in sorted(out.glob(f"{IMAGE_ID}__*.csv"))}


def test_run_writes_the_same_files_as_the_cli(tmp_path, patched_reader):  # noqa: F811
    cfg_path = _classical_config(tmp_path)
    out_cli, out_api = tmp_path / "cli", tmp_path / "api"
    assert cli.main(["run", NAME, "--config", str(cfg_path), "--out", str(out_cli), "--no-coloc"]) == 0
    res = api.run(NAME, config=cfg_path, out=out_api, overrides={"coloc.enabled": False})

    assert res["image_id"] == IMAGE_ID and res["n_nuclei"] >= 1 and res["out_dir"] == str(out_api)
    assert set(res["tables"]) >= {"nuclei", "spots", "image_summary"} and set(res["stages"])
    assert sorted(f.name for f in out_cli.iterdir()) == sorted(f.name for f in out_api.iterdir())
    a, b = _tables(out_cli), _tables(out_api)
    assert a.keys() == b.keys() and len(a) >= 5
    for name in a:
        pd.testing.assert_frame_equal(a[name], b[name], check_exact=True, obj=name)
    ma = provenance.read_done_marker(out_cli, IMAGE_ID)
    mb = provenance.read_done_marker(out_api, IMAGE_ID)
    for k in ("config_hash", "enabled_stages", "stage_hashes", "xy_stride", "z_range", "failed_stages"):
        assert ma[k] == mb[k], k
    assert "coloc" not in ma["enabled_stages"]            # the override took, on both paths alike
    assert a[f"{IMAGE_ID}__nuclei.csv"]["config_hash"].iloc[0] == mb["config_hash"]


def test_run_resolves_config_profile_geometry_and_default_out(tmp_path, monkeypatch):
    seen = []

    def fake_process_image(path, cfg, out_dir, *, xy_stride=1, z_range=None, prov=None):
        seen.append({"path": Path(path), "cfg": cfg, "out": Path(out_dir), "xy_stride": xy_stride,
                     "z_range": z_range, "prov": prov})
        return {"image_id": "x", "n_nuclei": 0, "qc_pass": True, "qc_flags": [], "out_dir": str(out_dir),
                "tables": {}, "stages": {}}

    monkeypatch.setattr(api, "process_image", fake_process_image)

    # a profile name, every geometry argument, and an override that enters the hash like --no-spots
    base = load_config(ROOT / "config" / "profiles" / "rad51_foci.yaml")
    api.run("g.nd2", profile="rad51_foci", out=tmp_path / "o", xy_stride=2, z_range=[3, 9],
            overrides={"spots.enabled": False})
    s = seen[-1]
    assert s["path"] == Path("g.nd2") and s["out"] == tmp_path / "o" and s["prov"] is None
    assert s["xy_stride"] == 2 and s["z_range"] == (3, 9)
    assert s["cfg"].get("spots.enabled") is False and s["cfg"].hash != base.hash
    assert s["cfg"].as_dict()["io"] == base.as_dict()["io"]

    # nothing but the image: the checkout's config/config.yaml and <stem>_results next to the image
    api.run(tmp_path / "gonad.ome.tif")
    s = seen[-1]
    assert s["out"] == tmp_path / "gonad_results" and s["xy_stride"] == 1 and s["z_range"] is None
    assert s["cfg"].hash == load_config(ROOT / "config" / "config.yaml").hash
    assert api.default_out_dir(Path("d") / "g.nd2") == Path("d") / "g_results"

    # a loaded Config passes through untouched (the CLI's own route into run)
    cfg = load_config(ROOT / "config" / "config.yaml")
    api.run("g.nd2", config=cfg, out=tmp_path / "c")
    assert seen[-1]["cfg"] is cfg and isinstance(seen[-1]["cfg"], Config)

    with pytest.raises(ValueError, match="not both"):
        api.run("g.nd2", config=cfg, profile="rad51_foci")
    with pytest.raises(FileNotFoundError, match="does_not_exist"):
        api.run("g.nd2", profile="does_not_exist")
    assert len(seen) == 3
