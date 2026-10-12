"""Running a pack (design 3.10, 3.11, 6.9 and 9; P1 acceptance 2 and 9, gate G8 for the test pack).

* api.resolve_config(pack="celegans_germline").hash is 84bdda5a5230, and for rad51_foci, segmentation_only,
  n2_sc, ccw77_partition and live_syp2 the pack route's hash equals the --profile hash.
* On the synthetic gonad, `run --pack celegans_germline` and `run --config config/config.yaml` write the same
  file list, cell-identical tables and equal done markers (every key but finished); the manifests differ only
  by the pack and runtime keys (and the time stamp), the stage records only by pack, runtime and post_run.
* `batch --pack` writes the pack and runtime records into the batch-root manifest; per-image __stages.json of
  pack runs carries post_run.segment_device_used; a default batch carries none of this.
* The test pack's recipe goldens pass with an empty model store, and read their OME-TIFF example with
  tifffile through process_image's read_fn seam (equal to the bioio reader where bioio is installed).
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
import pytest
from packs_helpers import (  # noqa: F401  (pack_env is a fixture)
    TEST_PACK,
    pack_env,
    patch_nd2_reader,
    small_stack,
    use_fixture_packs,
)
from test_end_to_end import patched_reader  # noqa: F401  (fixture: the synthetic 3-channel gonad)

from germquant import api, cli, packs, packs_golden
from germquant.config import load_config
from germquant.segment import backends as SB

ROOT = Path(__file__).resolve().parents[1]
NAME = "20251105_n2_nohs_HERM_001.nd2"
IID = "20251105_n2_nohs_HERM_001"


def test_celegans_pack_hashes_equal_the_frozen_files():
    assert api.resolve_config(pack="celegans_germline").hash == "84bdda5a5230"
    assert api.resolve_config(pack="celegans_germline@0.1.0", recipe="default").hash == "84bdda5a5230"
    for recipe in ("rad51_foci", "segmentation_only", "n2_sc", "ccw77_partition", "live_syp2"):
        assert api.resolve_config(pack="celegans_germline", recipe=recipe).hash == api.resolve_config(profile=recipe).hash
    cfg = api.resolve_config(pack="celegans_germline")
    assert cfg.pack_record["id"] == "celegans_germline" and cfg.pack_record["source"] == "builtin"
    assert cfg.as_dict() == load_config(ROOT / "config" / "config.yaml").as_dict()     # not data: never in as_dict
    assert load_config(ROOT / "config" / "config.yaml").pack_record is None


def test_resolve_config_takes_exactly_one_source(pack_env, monkeypatch):  # noqa: F811
    use_fixture_packs(monkeypatch)
    with pytest.raises(ValueError, match="not both"):
        api.resolve_config(config="config/config.yaml", pack="test_pack")
    with pytest.raises(ValueError, match="not both"):
        api.resolve_config(profile="generic", pack="test_pack")
    with pytest.raises(ValueError, match="recipe needs pack"):
        api.resolve_config(recipe="basic")
    with pytest.raises(FileNotFoundError, match="Packs available"):
        api.resolve_config(pack="nope")
    assert api.resolve_config(pack="test_pack").pack_record["recipe"] == "basic"
    for name in ("list_packs", "find_pack", "load_pack", "pack_status", "fetch_pack_models", "install_pack",
                 "list_segment_backends", "list_spot_backends", "known_models", "register_local_model"):
        assert name in api.__all__ and callable(getattr(api, name))
    assert {"cellpose", "classical"} <= set(api.list_segment_backends())
    assert {"spotmax", "log"} <= set(api.list_spot_backends())


def _tables(out: Path) -> dict[str, pd.DataFrame]:
    return {f.name: pd.read_csv(f, low_memory=False, float_precision="round_trip").drop(columns=["run_timestamp"])
            for f in sorted(out.glob(f"{IID}__*.csv"))}


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


@pytest.fixture
def classical_cellpose(monkeypatch):
    """The cellpose backend replaced by the classical watershed (no Cellpose, no weights in tests); both
    routes then segment identically, and the model file is the same (absent) for both."""
    monkeypatch.setattr(SB, "_REGISTRY", dict(SB._REGISTRY))
    monkeypatch.setattr(SB, "_SOURCES", dict(SB._SOURCES))     # restored afterwards, as both dicts must be
    SB.register_backend("cellpose", SB.classical_backend, replace=True)


def test_pack_and_config_routes_write_the_same_files(tmp_path, patched_reader, classical_cellpose,  # noqa: F811
                                                     pack_env):  # noqa: F811
    out_pack, out_cfg = tmp_path / "pack", tmp_path / "cfg"
    assert cli.main(["run", NAME, "--pack", "celegans_germline", "--no-pack-check", "--out", str(out_pack)]) == 0
    assert cli.main(["run", NAME, "--config", str(ROOT / "config" / "config.yaml"), "--out", str(out_cfg)]) == 0
    assert sorted(f.name for f in out_pack.iterdir()) == sorted(f.name for f in out_cfg.iterdir())
    a, b = _tables(out_pack), _tables(out_cfg)
    assert a.keys() == b.keys() and len(a) >= 12
    for name in a:
        pd.testing.assert_frame_equal(a[name], b[name], check_exact=True, obj=name)
    ma, mb = _json(out_pack / f"{IID}__done.json"), _json(out_cfg / f"{IID}__done.json")
    assert sorted(ma) == sorted(mb)
    assert {k: v for k, v in ma.items() if k != "finished"} == {k: v for k, v in mb.items() if k != "finished"}
    pa, pb = _json(out_pack / "run_manifest.json"), _json(out_cfg / "run_manifest.json")
    assert set(pa) - set(pb) == {"pack", "runtime"} and not set(pb) - set(pa)
    assert {k: v for k, v in pa.items() if k not in ("pack", "runtime", "run_timestamp")} == \
           {k: v for k, v in pb.items() if k != "run_timestamp"}
    assert pa["config_hash"] == "84bdda5a5230" and pa["pack"]["id"] == "celegans_germline"
    assert set(pa["runtime"]) == {"device_requested", "device_selected", "torch_version", "cuda_available",
                                  "python", "platform"}
    sa, sb = _json(out_pack / f"{IID}__stages.json"), _json(out_cfg / f"{IID}__stages.json")
    assert set(sa) - set(sb) == {"pack", "runtime", "post_run"} and not set(sb) - set(sa)
    assert sa["stage_hashes"] == sb["stage_hashes"]


def test_batch_pack_records_and_default_batch_without_them(tmp_path, pack_env, monkeypatch):  # noqa: F811
    use_fixture_packs(monkeypatch)
    patch_nd2_reader(monkeypatch, small_stack())
    folder = tmp_path / "in"
    (folder / "sub").mkdir(parents=True)
    for n in ("a_001.nd2", "sub/b_002.nd2"):
        (folder / n).write_bytes(b"x")
    out = tmp_path / "out_pack"
    assert cli.main(["batch", str(folder), "--pack", "test_pack", "--out", str(out)]) == 0
    man = _json(out / "run_manifest.json")
    assert man["pack"]["id"] == "test_pack" and man["pack"]["recipe"] == "basic"
    assert man["runtime"]["device_requested"] == "auto" and "python" in man["runtime"]
    for rel, iid in (("a_001", "a_001"), ("sub/b_002", "b_002")):
        rec = _json(out / rel / f"{iid}__stages.json")
        assert rec["post_run"] == {"segment_device_used": None}          # the classical recipe: Cellpose never ran
        assert rec["pack"]["id"] == "test_pack" and rec["runtime"]
        assert not (out / rel / "run_manifest.json").exists()             # batch writes the batch-root one only
    assert cli.main(["batch", str(folder), "--pack", "test_pack", "--out", str(out), "--resume"]) == 0
    # the same recipe by path: a default batch writes none of the pack records
    out2 = tmp_path / "out_cfg"
    assert cli.main(["batch", str(folder), "--config", str(TEST_PACK / "config" / "profiles" / "basic.yaml"),
                     "--out", str(out2)]) == 0
    man2 = _json(out2 / "run_manifest.json")
    assert "pack" not in man2 and "runtime" not in man2
    rec2 = _json(out2 / "a_001" / "a_001__stages.json")
    assert not {"pack", "runtime", "post_run"} & set(rec2)
    assert man2["config_hash"] == man["config_hash"]
    t1 = pd.read_csv(out / "a_001" / "a_001__nuclei.csv").drop(columns=["run_timestamp"])
    t2 = pd.read_csv(out2 / "a_001" / "a_001__nuclei.csv").drop(columns=["run_timestamp"])
    pd.testing.assert_frame_equal(t1, t2, check_exact=True)


def test_test_pack_recipe_goldens_pass_with_an_empty_model_store(pack_env, monkeypatch):  # noqa: F811
    use_fixture_packs(monkeypatch)
    seen = []
    real = packs_golden._empty_models_dir

    class _Watch:
        def __init__(self):
            self.cm = real()

        def __enter__(self):
            self.dir = self.cm.__enter__()
            assert os.environ["GERMQUANT_MODELS"] == str(self.dir) and not any(self.dir.iterdir())
            return self.dir

        def __exit__(self, *exc):
            seen.append(sorted(p.name for p in self.dir.rglob("*")))       # nothing resolved, nothing fetched
            return self.cm.__exit__(*exc)

    monkeypatch.setattr(packs_golden, "_empty_models_dir", _Watch)
    results = packs_golden.run_pack_goldens(packs.find_pack("test_pack"))
    assert [r["recipe"] for r in results] == ["basic", "with_model"]
    assert all(r["ok"] for r in results), [r["problems"] for r in results]
    assert seen == [[], []]
    assert os.environ["GERMQUANT_MODELS"] == str(pack_env.models)                # restored
    assert cli.main(["pack", "golden", "test_pack"]) == 0
    assert cli.main(["pack", "golden", "test_pack", "--tier", "model"]) == 2


def test_golden_comparison_catches_a_change():
    ref = pd.DataFrame({"n_spots": [2, 2], "marker": ["foci", "foci"]})
    assert packs_golden.compare_table(ref.copy(), ref, "nuclei") == []
    changed = ref.copy()
    changed.loc[0, "n_spots"] = 3
    assert packs_golden.compare_table(changed, ref, "nuclei") == ["nuclei.n_spots differs"]
    assert packs_golden.compare_table(ref.drop(columns=["marker"]), ref, "nuclei")


def test_tifffile_example_reader_matches_the_bioio_reader():
    example = TEST_PACK / "example" / "test_pack_synthetic.ome.tif"
    meta, st = packs_golden.read_example(example)
    assert st.data.shape == (2, 12, 64, 96) and st.spacing == (0.5, 0.2, 0.2) and st.spacing_ok
    assert st.channel_names == ["DAPI", "foci"] and meta["channel_names"] == ["DAPI", "foci"]
    m2, s2 = packs_golden.read_example(example, xy_stride=2, z_range=(2, 8))
    assert s2.data.shape == (2, 6, 32, 48) and s2.spacing == (0.5, 0.4, 0.4) and m2["spacing"] == (0.5, 0.2, 0.2)
    pytest.importorskip("bioio")
    pytest.importorskip("bioio_ome_tiff")
    from germquant.io.bioio_reader import read_bioio, read_bioio_metadata

    b = read_bioio(example)
    assert (b.data == st.data).all() and b.data.dtype == st.data.dtype
    assert tuple(b.spacing) == st.spacing and b.channel_names == st.channel_names and b.spacing_ok == st.spacing_ok
    bm = read_bioio_metadata(example)
    assert tuple(bm["spacing"]) == meta["spacing"] and bm["channel_names"] == meta["channel_names"]
