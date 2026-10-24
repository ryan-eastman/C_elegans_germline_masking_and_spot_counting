"""Single-plane runs through process_image (tranche 5 package P2b; design 4.3 to 4.9 and 12.3, P2b acceptance 2, 3,
5 and 6; docs/PLANAR.md).

The calibrated plane of tests/planar_helpers.py (dy = dx = 0.8 um, nine nuclei, two specks of debris, 30 puncta)
runs through the .nd2 reader seam under config/profiles/generic_2d.yaml with the classical segmenter forced, as
tests/test_golden.py forces it:

* test_synthetic_2d_golden: the tables equal tests/golden/synthetic_2d/ (sentinel: generic_2d's own hash), and
  the run meets acceptance 2: every planted nucleus is found and the debris is gated away; at least 90 percent of
  the planted puncta are found with at most 10 percent extras; volume_um3 NA; area_um2 within 2 percent of the
  planted area; every z_um 0.0; voxel_dz_um NaN; read:planar raised; the label TIF reads back as (Y, X) with
  PhysicalSizeX and Y and no Z, and the spots TIF is YX too.
* test_planar_tables_ignore_file_dz (G10): the same plane with a file dz of 0.2 and of 5.0 gives identical tables.
* test_the_z_merge_a_recipe_asks_for_is_flagged: spots.merge_z_columns true raises spots:z_merge_off_planar and
  changes no spot.
* test_the_yx_writers: write_labels_yx, write_spots_yx and write_mask_yx write YX OME and ImageJ TIFFs with the
  run's dy and dx (non-square pixels, so a swap shows).
* test_log_net3d_nuclei3d_refused_before_pixels: log, net3d (main spots and a spots instance) and nuclei3d stop a
  single-plane run before the segmenter or the detector sees a pixel (spies on both) and before any file.
* test_stage_hashes_carry_planar_only_in_planar_runs: "planar": true enters the read, segment, measure and spots
  payloads of a planar run only; the default hashes stay the fence's pins, and a 3D run under generic_2d's auto
  hashes as it would without the key.
* test_batch_resume_planar_folder and test_a_default_batch_resume_decides_as_before: `batch --resume` skips the
  finished planar images, reprocesses after a recipe edit, and makes today's is_done call on a default batch.

Regenerate the golden only for a deliberate change of this path: ``GERMQUANT_UPDATE_GOLDEN=1 pytest
tests/test_planar_run.py -k golden``, then check the run against the planted counts by hand.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
import planar_helpers as H
import pytest

from germquant import batch, cli, pipeline, planar, provenance, schema, stages
from germquant.config import load_config

ROOT = Path(__file__).resolve().parents[1]
GOLD = ROOT / "tests" / "golden" / "synthetic_2d"
PINS = json.loads((ROOT / "tests" / "golden" / "default_invariants.json").read_text(encoding="utf-8"))
TABLES = tuple(schema.TABLES)
ND2 = f"{H.STEM}.nd2"
PLANAR_STAGES = ("read", "segment", "measure", "spots")


def _run(monkeypatch, out: Path, *, dz: float = H.DZ, cfg=None, **kw) -> dict:
    with monkeypatch.context() as mp:
        H.patch_reader(mp, dz=dz)
        return pipeline.process_image(ND2, cfg if cfg is not None else H.generic_2d(), out, prov=None, **kw)


@pytest.fixture(scope="module")
def golden_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("synthetic_2d") / "results"
    cfg = H.generic_2d()
    with pytest.MonkeyPatch.context() as mp:
        H.patch_reader(mp)
        res = pipeline.process_image(ND2, cfg, out, prov=None)
    return {"res": res, "out": out, "cfg": cfg}


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def test_synthetic_2d_golden(golden_run):
    out, cfg, res = golden_run["out"], golden_run["cfg"], golden_run["res"]
    iid = res["image_id"]
    if os.environ.get("GERMQUANT_UPDATE_GOLDEN"):
        H.write_golden(out, iid, GOLD, TABLES, cfg.hash)
        pytest.skip("synthetic_2d golden regenerated: check it against the planted counts")
    sentinel = (GOLD / "config_hash.txt").read_text(encoding="utf-8").strip()
    assert cfg.hash == sentinel == PINS["golden_sentinels"]["synthetic_2d"] == \
        PINS["config_hashes"]["config/profiles/generic_2d.yaml"], "generic_2d.yaml or a base it extends changed"
    H.compare_tables(out, iid, GOLD, TABLES)


def test_synthetic_2d_meets_acceptance_2(golden_run):
    import tifffile

    out, res = golden_run["out"], golden_run["res"]
    iid = res["image_id"]
    info = H.planted()
    labels = tifffile.imread(out / f"{iid}__nuclei_labels.tif")
    assert labels.shape == (H.SIZE, H.SIZE) and labels.dtype == np.int32
    # every planted nucleus is found, pixel for pixel; the debris falls to the area gate (min_area_um2 20)
    yy, xx = np.indices(labels.shape)
    assert res["n_nuclei"] == len(info["nuclei"]) == int(labels.max())
    for n in info["nuclei"]:
        disk = (yy - n["y"]) ** 2 + (xx - n["x"]) ** 2 <= n["r_px"] ** 2
        (lab,) = set(np.unique(labels[disk]))
        assert lab > 0 and np.array_equal(labels == lab, disk)
    assert all(labels[d["y"], d["x"]] == 0 for d in info["debris"])
    nuclei = pd.read_csv(out / f"{iid}__nuclei.csv").set_index("nucleus_id")
    cols = list(pd.read_csv(out / f"{iid}__nuclei.csv", nrows=0).columns)
    assert cols[:len(schema.NUCLEI)] == schema.NUCLEI
    assert cols[len(schema.NUCLEI):len(schema.NUCLEI) + 2] == ["area_um2", "equivalent_diameter_um"]
    assert nuclei["volume_um3"].isna().all() and (nuclei["centroid_z_um"] == 0.0).all()
    for n in info["nuclei"]:
        row = nuclei.loc[int(labels[n["y"], n["x"]])]
        assert math.isclose(row["area_um2"], n["area_um2"], rel_tol=1e-12)          # the planted pixels
        assert abs(row["area_um2"] - n["disk_um2"]) <= 0.02 * n["disk_um2"]          # the planted disk, 2 percent
        assert math.isclose(row["equivalent_diameter_um"], 2 * math.sqrt(row["area_um2"] / math.pi))
        assert math.isclose(row["centroid_y_um"], n["y"] * H.DY) and math.isclose(row["centroid_x_um"], n["x"] * H.DX)
    assert nuclei["voxel_dz_um"].isna().all() and (nuclei["voxel_dy_um"] == H.DY).all()
    # at least 90 percent of the planted puncta, at most 10 percent extras, every spot at z 0
    spots = pd.read_csv(out / f"{iid}__spots.csv")
    found, extras = H.matched_puncta(spots, labels)
    n_planted = len(info["puncta"])
    assert found >= 0.9 * n_planted and extras <= 0.1 * n_planted, (found, extras, n_planted)
    assert (spots["z_um"] == 0.0).all() and spots["nucleus_id"].gt(0).all()
    assert "stage_flags" not in res["tables"]["spots"].attrs
    assert nuclei["n_spots"].sum() == len(spots)
    # flags and records of a single-plane run
    flags = res["qc_flags"]
    assert "read:planar" in flags and "spots:z_merge_off_planar" not in flags      # generic_2d asks for no merge
    assert not [f for f in flags if f.startswith("voxel:")]
    manifest, record = _json(out / "run_manifest.json"), _json(out / f"{iid}__stages.json")
    for rec in (manifest, record):
        assert rec["planar"] is True and rec["spacing_source"] == "file"
    assert "read:planar" in record["stages"]["read"]["flags"]
    done = _json(out / f"{iid}__done.json")
    assert done["planar"] is True and sorted(set(done) - set(PINS["runs"]["synthetic"]["done_marker_keys"])) == ["planar"]
    # the label TIF reads back as (Y, X) with PhysicalSizeX and Y and no Z; the spots TIF is YX too
    with tifffile.TiffFile(out / f"{iid}__nuclei_labels.tif") as tif:
        assert tif.series[0].axes == "YX"
        px = tif.ome_metadata
    assert f'PhysicalSizeX="{H.DX}"' in px and f'PhysicalSizeY="{H.DY}"' in px and "PhysicalSizeZ" not in px
    spots_img = tifffile.imread(out / f"{iid}__spots.tif")
    assert spots_img.shape == (H.SIZE, H.SIZE) and spots_img.dtype == np.uint16 and spots_img.max() == 65535
    status = {k: v["status"] for k, v in res["stages"].items()}
    assert status["segment"] == status["measure"] == status["spots"] == status["signal_qc"] == "ran"
    assert not [s for s, v in status.items() if v == "failed"]
    assert res["tables"]["image_summary"]["segmentation_method"].iloc[0] == "classical"
    assert (out / f"{iid}__montage.png").exists()
    for t in schema.TABLES:          # the 2D columns appear in the nuclei table only
        if t != "nuclei":
            assert "area_um2" not in pd.read_csv(out / f"{iid}__{t}.csv", nrows=0).columns


def test_planar_tables_ignore_file_dz(golden_run, tmp_path, monkeypatch):
    for dz in (0.2, 5.0):
        res = _run(monkeypatch, tmp_path / f"dz{dz}", dz=dz)
        for t in TABLES:
            a = H.read_table(Path(res["out_dir"]) / f"{res['image_id']}__{t}.csv")
            b = H.read_table(golden_run["out"] / f"{res['image_id']}__{t}.csv")
            pd.testing.assert_frame_equal(a, b, check_exact=True, obj=f"{t} at file dz {dz}")


def test_the_z_merge_a_recipe_asks_for_is_flagged(golden_run, tmp_path, monkeypatch):
    """spots.merge_z_columns true (config.yaml's value) cannot run on one plane: the spots stage raises
    spots:z_merge_off_planar, and the spots are those of the run without the merge."""
    res = _run(monkeypatch, tmp_path / "merge", cfg=H.generic_2d(**{"spots.merge_z_columns": True}))
    out = Path(res["out_dir"])
    assert "spots:z_merge_off_planar" in res["qc_flags"]
    assert "spots:z_merge_off_planar" in _json(out / f"{res['image_id']}__stages.json")["stages"]["spots"]["flags"]
    # the flags come back through attrs, which the stage removes: no table and no parquet file carries them
    assert "stage_flags" not in res["tables"]["spots"].attrs
    assert "stage_flags" not in pd.read_parquet(out / f"{res['image_id']}__spots.parquet").attrs
    pd.testing.assert_frame_equal(H.read_table(out / f"{res['image_id']}__spots.csv"),
                                  H.read_table(golden_run["out"] / f"{res['image_id']}__spots.csv"), check_exact=True)


def test_the_yx_writers(tmp_path):
    import tifffile

    labels = np.zeros((1, 20, 30), np.int32)
    labels[0, 2:6, 3:9] = 1
    sp = (0.5, 0.7)                                     # dy, dx: not square, so a swap would show
    planar.write_labels_yx(labels, tmp_path / "l.tif", sp)
    with tifffile.TiffFile(tmp_path / "l.tif") as tif:
        arr, axes, ome = tif.asarray(), tif.series[0].axes, tif.ome_metadata
    assert arr.shape == (20, 30) and arr.dtype == np.int32 and axes == "YX" and np.array_equal(arr, labels[0])
    assert 'PhysicalSizeX="0.7"' in ome and 'PhysicalSizeY="0.5"' in ome and "PhysicalSizeZ" not in ome
    spots = pd.DataFrame({"z_um": [0.0, 0.0], "y_um": [1.0, 7.0], "x_um": [7.0, 14.0]})
    planar.write_spots_yx(spots, (20, 30), sp, tmp_path / "s.tif", radius_um=1.0)
    with tifffile.TiffFile(tmp_path / "s.tif") as tif:
        img, ij, res = tif.asarray(), tif.imagej_metadata, tif.pages[0].resolution
    assert img.shape == (20, 30) and img.dtype == np.uint16 and ij.get("unit") == "um" and "spacing" not in ij
    assert np.allclose(res, (1 / 0.7, 1 / 0.5))
    assert img[2, 10] == 65535 and img[14, 20] == 65535 and img[2, 20] == 0       # (y, x) = (y_um / dy, x_um / dx)
    assert img[2, 10 + 2] == 0 and img[2 + 2, 10] == 65535                     # radius 1 um: 1 px in x, 2 px in y
    mask = np.zeros((1, 20, 30), bool)
    mask[0, 5, 5] = True
    planar.write_mask_yx(mask, sp, tmp_path / "m.tif")
    m = tifffile.imread(tmp_path / "m.tif")
    assert m.shape == (20, 30) and m.dtype == np.uint8 and m[5, 5] == 255 and m.sum() == 255


def test_log_net3d_nuclei3d_refused_before_pixels(tmp_path, monkeypatch):
    from germquant import spots as spots_pkg

    seen = []
    monkeypatch.setattr(planar, "segment_planar", lambda *a, **k: seen.append("segment_planar"))
    monkeypatch.setattr(pipeline, "segment_nuclei", lambda *a, **k: seen.append("segment_nuclei"))
    monkeypatch.setattr(spots_pkg, "detect_spots", lambda *a, **k: seen.append("detect_spots"))
    cases = {
        "log": ({"spots.backend": "log"}, "the spot backend 'log' does not declare planar_ok"),
        "net3d": ({"spots.backend": "net3d"}, "the spot backend 'net3d' does not declare planar_ok"),
        "instance": ({"spots.instances": [{"name": "second", "role": "foci", "backend": "log"}]},
                     "spots_second: the spot backend 'log' does not declare planar_ok"),
        "nuclei3d": ({"segmentation.nuclei.method": "cellpose", "segmentation.nuclei.backend": "nuclei3d"},
                     "'nuclei3d' does not declare planar_ok"),
    }
    for name, (overrides, message) in cases.items():
        out = tmp_path / name
        with pytest.raises(planar.PlanarConfigError, match="single-plane run refused") as err:
            _run(monkeypatch, out, cfg=H.generic_2d(**overrides))
        assert message in str(err.value), (name, str(err.value))
        assert out.is_dir() and not list(out.iterdir()), f"{name}: files written before the refusal"
    assert not seen


def test_stage_hashes_carry_planar_only_in_planar_runs(golden_run, tmp_path, monkeypatch):
    cfg = load_config(ROOT / "config" / "config.yaml")
    tools = {t: "PINNED" for ts in stages._STAGE_TOOLS.values() for t in ts}
    for case, pin in PINS["default_stage_hashes"].items():
        default = stages.stage_hashes(cfg, pin["roles"], tools, "absent")
        assert default == pin["hashes"], case
        flat = stages.stage_hashes(cfg, pin["roles"], tools, "absent", planar=True)
        assert sorted(s for s in default if default[s] != flat[s]) == sorted(PLANAR_STAGES), case
    # what the planar run wrote: the planar payload, in all three records
    out, res = golden_run["out"], golden_run["res"]
    run_cfg = golden_run["cfg"]
    roles, _ = run_cfg.channel_map.resolve(list(H.CHANNELS))
    want = stages.stage_hashes(run_cfg, roles, provenance.tool_versions(), "unused", planar=True)
    for name in ("run_manifest.json", f"{res['image_id']}__stages.json", f"{res['image_id']}__done.json"):
        assert _json(out / name)["stage_hashes"] == want, name
    # a 3D run under generic_2d (auto on a five-plane file): today's payload, no planar record
    with monkeypatch.context() as mp:
        H.patch_reader(mp, H.stack_3d(), two_d=False)
        res3 = pipeline.process_image(ND2, H.generic_2d(**{"spots.enabled": False}), tmp_path / "3d", prov=None)
    cfg3 = H.generic_2d(**{"spots.enabled": False})
    roles3, _ = cfg3.channel_map.resolve(list(H.CHANNELS))
    got = _json(tmp_path / "3d" / "run_manifest.json")
    assert got["stage_hashes"] == stages.stage_hashes(cfg3, roles3, provenance.tool_versions(), "unused")
    assert "planar" not in got and "planar" not in _json(tmp_path / "3d" / f"{res3['image_id']}__done.json")


def _study(tmp_path: Path, **nuclei) -> Path:
    keys = "".join(f"    {k}: {v}\n" for k, v in nuclei.items())
    p = tmp_path / "study_2d.yaml"
    p.write_text(f"extends: {(ROOT / 'config' / 'profiles' / 'generic_2d.yaml').as_posix()}\n"
                 f"segmentation:\n  nuclei:\n    method: classical\n{keys}spots:\n  enabled: false\n", encoding="utf-8")
    return p


def test_batch_resume_planar_folder(tmp_path, monkeypatch):
    src = tmp_path / "in"
    src.mkdir()
    for n in ("a_001.nd2", "b_002.nd2"):
        (src / n).write_bytes(b"")
    H.patch_reader(monkeypatch)
    calls, decisions = [], []
    real_process, real_done = pipeline.process_image, batch.is_done

    def process(path, *a, **k):
        calls.append(Path(path).name)
        return real_process(path, *a, **k)

    def is_done(*a, **k):
        decisions.append(k)
        return real_done(*a, **k)

    monkeypatch.setattr(pipeline, "process_image", process)
    monkeypatch.setattr(batch, "is_done", is_done)
    study, out = _study(tmp_path), tmp_path / "out"
    assert cli.main(["batch", str(src), "--config", str(study), "--out", str(out)]) == 0
    assert calls == ["a_001.nd2", "b_002.nd2"]
    for iid in ("a_001", "b_002"):
        assert _json(out / iid / f"{iid}__done.json")["planar"] is True
    calls.clear()
    assert cli.main(["batch", str(src), "--config", str(study), "--out", str(out), "--resume"]) == 0
    assert calls == []                                       # both planar images are done
    assert all(set(d) == {"xy_stride", "segment_hash", "segment_hash_planar"} for d in decisions)
    assert decisions[0]["segment_hash"] != decisions[0]["segment_hash_planar"]
    study.write_text(study.read_text(encoding="utf-8").replace("method: classical\n",
                                                               "method: classical\n    min_area_um2: 25.0\n"),
                     encoding="utf-8")
    assert cli.main(["batch", str(src), "--config", str(study), "--out", str(out), "--resume"]) == 0
    assert calls == ["a_001.nd2", "b_002.nd2"]               # a recipe edit reprocesses both


def test_a_default_batch_resume_decides_as_before(tmp_path, monkeypatch):
    """A recipe without io.planar: cli._batch makes exactly today's is_done call (no planar hash)."""
    from test_batch import _fake_process, _study as _batch_study

    src, cfg_path = _batch_study(tmp_path)
    calls, decisions = [], []
    real_done = batch.is_done
    monkeypatch.setattr(pipeline, "process_image", _fake_process(calls))
    monkeypatch.setattr(batch, "is_done", lambda *a, **k: decisions.append(k) or real_done(*a, **k))
    out = tmp_path / "out"
    assert cli.main(["batch", str(src), "--config", str(cfg_path), "--out", str(out)]) == 0
    calls.clear()
    assert cli.main(["batch", str(src), "--config", str(cfg_path), "--out", str(out), "--resume"]) == 0
    assert calls == [] and len(decisions) == 3
    cfg = load_config(cfg_path)
    want = stages.stage_hashes(cfg, {}, provenance.tool_versions(), "unused")["segment"]
    assert decisions == [{"xy_stride": 1, "segment_hash": want}] * 3
