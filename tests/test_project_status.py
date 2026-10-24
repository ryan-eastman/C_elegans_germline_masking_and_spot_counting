"""Project status (docs/PROJECTS.md, "Status"; design 7.7; PLAN.md 9.3 and 4.4 P5-RUN; P5 acceptance 3, run half).

* Every status of section 9.3 (done, stale, failed, missing, excluded, unreachable) and every reason code
  (config_hash, stages, geometry, planar_flag, segment_hash, failed_stages; sample_row comes with the sample
  sheet), with qc_pass, the last run id, the device used and the error.
* An edit marks every finished image stale with the reason (acceptance 3).
* The environment-mix warning, the lock report, the pack and model reports, the runs merged from every log, and
  `germquant project status --json` equal to api.project_status.
"""
from __future__ import annotations

import json
import os
import shutil
from datetime import datetime, timedelta, timezone

import pytest
from packs_helpers import pack_env  # noqa: F401
from project_helpers import BASIC_HASH, TEST_PACK_REF, project_env  # noqa: F401
from test_project_run import (
    _project_model,
    hand_edit_recipe,
    log_events,
    need_p2b,
    need_p2c,
    need_readers,
    project_with_images,
)

import germquant.project as P
from germquant import api, cli, cli_projects
from germquant.project import manifest as M
from germquant.project import runlog
from germquant.project import status as S

IMAGE_KEYS = ["id", "path", "status", "reason_code", "reason", "out_dir", "qc_pass", "last_run_id", "device_used",
              "error"]


def _fabricate(proj, src: str, dest: str, **marker_changes) -> None:
    """results/<dest>/ holding the marker and image summary of results/<src>/, re-named for dest, the marker's
    fields changed as given."""
    a, d = proj.root / "results" / src, proj.root / "results" / dest
    d.mkdir(parents=True)
    marker = json.loads((a / f"{src}__done.json").read_text(encoding="utf-8"))
    marker.update(image_id=dest, **marker_changes)
    (d / f"{dest}__done.json").write_text(json.dumps(marker, indent=1), encoding="utf-8")
    text = (a / f"{src}__image_summary.csv").read_text(encoding="utf-8")
    (d / f"{dest}__image_summary.csv").write_text(text.replace(f",{src},", f",{dest},"), encoding="utf-8")


def test_every_status_and_reason_code(project_env, tmp_path):  # noqa: F811
    need_readers()
    names = ("a", "cfg", "stg", "geo", "seg", "fst", "mis", "fld", "exc", "unr")
    proj = project_with_images(tmp_path, names=names)
    done = P.run_project(proj, images=["a"])
    assert done["exit_code"] == 0
    marker = json.loads((proj.root / "results" / "a" / "a__done.json").read_text(encoding="utf-8"))
    _fabricate(proj, "a", "cfg", config_hash="000000000000")
    _fabricate(proj, "a", "stg", enabled_stages=[s for s in marker["enabled_stages"] if s != "spots"])
    _fabricate(proj, "a", "geo", xy_stride=2)
    _fabricate(proj, "a", "seg", stage_hashes={**marker["stage_hashes"], "segment": "000000000000"})
    _fabricate(proj, "a", "fst", failed_stages=["spots"])
    (proj.root / "images" / "fld.ome.tif").write_bytes(b"no image")
    failed = P.run_project(proj, images=["fld"])
    assert failed["exit_code"] == 1
    m = M.read_manifest(proj.root)
    next(e for e in m["images"] if e["id"] == "exc").update(excluded=True, reason="out of focus")
    M.write_manifest(proj.root, m)
    os.remove(proj.root / "images" / "unr.ome.tif")

    st = P.project_status(proj.root)
    assert st["format"] == 1 and st["root"] == str(proj.root) and st["name"] == "proj"
    assert st["config_hash"] == st["expected_config_hash"] == BASIC_HASH and st["edited"] is False
    assert st["validated"]["value"] is False and "provisional" in st["validated"]["why"]
    assert st["pack"] == {"id": "test_pack", "pinned": TEST_PACK_REF, "installed": ["0.1.0"], "project_copy": True}
    rows = {r["id"]: r for r in st["images"]}
    assert [r["id"] for r in st["images"]] == list(names) and all(list(r) == IMAGE_KEYS for r in st["images"])
    want = {"a": ("done", None), "cfg": ("stale", "config_hash"), "stg": ("stale", "stages"),
            "geo": ("stale", "geometry"), "seg": ("stale", "segment_hash"), "fst": ("failed", "failed_stages"),
            "mis": ("missing", None), "fld": ("failed", None), "exc": ("excluded", None),
            "unr": ("unreachable", None)}
    assert {i: (r["status"], r["reason_code"]) for i, r in rows.items()} == want
    a = rows["a"]
    assert a["reason"] == "done" and a["qc_pass"] is True and a["last_run_id"] == done["run_id"]
    assert a["out_dir"] == str(proj.root / "results" / "a") and a["path"] == str(proj.root / "images" / "a.ome.tif")
    assert a["device_used"] is None and a["error"] is None                # the classical segmenter: no Cellpose
    assert "config_hash 000000000000" in rows["cfg"]["reason"] and rows["cfg"]["qc_pass"] is True
    assert "xy_stride" in rows["geo"]["reason"] and "spots" in rows["fst"]["reason"]
    assert rows["fld"]["last_run_id"] == failed["run_id"] and rows["fld"]["error"]
    assert "the last run failed" in rows["fld"]["reason"]
    assert rows["exc"]["reason"] == "excluded: out of focus"
    assert "not reachable" in rows["unr"]["reason"] and "NAS or VPN" in rows["unr"]["reason"]
    assert rows["mis"]["reason"] == "no results yet" and rows["mis"]["last_run_id"] is None
    assert st["lock"] == {"held": False, "host": None, "pid": None, "since": None, "stale": False}
    assert [r["run_id"] for r in st["runs"]] == [done["run_id"], failed["run_id"]]
    assert st["runs"][0] == {"run_id": done["run_id"], "via": "cli", "started": st["runs"][0]["started"],
                             "ended": st["runs"][0]["ended"], "counts": {"done": 1, "failed": 0, "skipped": 0},
                             "interrupted": False, "finished": True}
    # the reason codes come from the reasons batch.is_done gives
    assert S.reason_code("config_hash abc != def") == "config_hash"
    assert S.reason_code("enabled stages differ") == "stages"
    for why in ("xy_stride 1 != 2", "z_range differs", "t_index 3 != 4"):
        assert S.reason_code(why) == "geometry"
    assert S.reason_code("stage(s) failed last run: spots") == "failed_stages"
    assert S.reason_code("segment stage hash differs (model file or Cellpose version changed)") == "segment_hash"
    assert S.reason_code("planar flag: a single-plane marker, and no planar segment hash was given") == "planar_flag"
    assert S.reason_code("done") is None and S.reason_code("no completion marker") is None
    # a run that only skips an image does not become its last run: last_run_id names the run that processed it
    again = P.run_project(proj, images=["a", "fld"])
    assert again["counts"] == {"done": 0, "failed": 1, "skipped": 1}
    rows = {r["id"]: r for r in P.project_status(proj.root)["images"]}
    assert rows["a"]["last_run_id"] == done["run_id"] and rows["fld"]["last_run_id"] == again["run_id"]


def test_a_planar_marker_reads_stale_with_the_planar_flag(project_env, tmp_path):  # noqa: F811
    need_readers()
    need_p2b()
    from germquant import pipeline, provenance, stages

    proj = project_with_images(tmp_path, names=("a", "pla"))
    assert P.run_project(proj, images=["a"])["exit_code"] == 0
    # what a single-plane run of this recipe records: the marker's extra and the planar segment sub-hash
    cfg = P.project_config(proj)
    planar_seg = stages.stage_hashes(cfg, {}, provenance.tool_versions(), "unused",
                                     segment_backend=pipeline.segment_backend_provenance(cfg), planar=True)["segment"]
    marker = json.loads((proj.root / "results" / "a" / "a__done.json").read_text(encoding="utf-8"))
    _fabricate(proj, "a", "pla", planar=True, stage_hashes={**marker["stage_hashes"], "segment": planar_seg})
    rows = {r["id"]: r for r in P.project_status(proj)["images"]}
    assert (rows["pla"]["status"], rows["pla"]["reason_code"]) == ("stale", "planar_flag")


def test_status_reports_a_model_that_is_not_here(project_env, tmp_path, capsys):  # noqa: F811
    """A recipe whose segmentation model is not on this computer (the state of a new computer before `pack
    fetch-models`): status still answers, with the model missing and needed and a warning, while a run is refused
    (design 7.7: status reports whether models are present and verified)."""
    proj = project_with_images(tmp_path, names=("a",), recipe="with_model")       # the pack's tinynet, not fetched
    st = P.project_status(proj)
    (tiny,) = [m for m in st["models"] if m["ref"] == "tinynet@f5055f8a"]
    assert tiny["present"] is False and tiny["needed"] is True
    assert [(r["id"], r["status"]) for r in st["images"]] == [("a", "missing")]
    assert any("tinynet@f5055f8a" in w and "missing" in w for w in st["warnings"]), st["warnings"]
    assert not any("cannot be resolved here" in w for w in st["warnings"]), st["warnings"]   # said once
    assert cli.main(["project", "status", str(proj.root), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["models"] == json.loads(json.dumps(st["models"]))
    with pytest.raises(P.ProjectError, match="missing"):
        P.run_project(proj)
    # a model name the engine cannot resolve at all (no file, not registered, not a Cellpose built-in)
    hand_edit_recipe(proj, {"segmentation": {"nuclei": {"cellpose_model": "no_such_model"}}})
    st = P.project_status(proj)
    assert [r["status"] for r in st["images"]] == ["missing"]
    assert any("cannot be resolved here" in w and "no_such_model" in w for w in st["warnings"]), st["warnings"]
    with pytest.raises(P.ProjectError, match="cannot be resolved"):
        P.run_project(proj)
    # a project model removed after use-model: reported missing, never raised
    model, _sha = _project_model(project_with_images(tmp_path, names=("b",), recipe="with_model", folder="p2"))
    model.unlink()
    st = P.project_status(tmp_path / "p2")
    (rec,) = [m for m in st["models"] if m["key"] == "segmentation.nuclei.cellpose_model"]
    assert (rec["present"], rec["needed"]) == (False, True)
    assert [r["status"] for r in st["images"]] == ["missing"]
    assert any("training/r01/models/m1/m1" in w and "missing" in w for w in st["warnings"]), st["warnings"]


def test_status_names_the_refusal_of_a_series(project_env, tmp_path):  # noqa: F811
    """Before any run, a time series that the recipe cannot run as one (here labels_from, D15) reads missing with
    the refusal every run fails it with, in the run's own words."""
    need_p2c()
    from test_project_run import _series_image

    from germquant.io.labels import SERIES_REFUSAL

    proj = project_with_images(tmp_path, names=("a",))
    P.add_images(proj, [_series_image(proj.root / "images")])
    hand_edit_recipe(proj, {"segmentation": {"nuclei": {"labels_from": "{image_dir}/{image_stem}_labels.tif"}}})
    rows = {r["id"]: r for r in P.project_status(proj)["images"]}
    assert rows["a"]["reason"] == "no results yet"                              # a fixed image: no such note
    assert rows["movie"]["status"] == "missing" and SERIES_REFUSAL in rows["movie"]["reason"]
    res = P.run_project(proj, images=["movie"])
    (e,) = [e for e in log_events(res["log"]) if e["event"] == "image_failed"]
    assert e["error"].startswith("REFUSED: a time series of 2 timepoints") and e["error"] in rows["movie"]["reason"]


def test_edit_marks_every_image_stale_with_reason(project_env, tmp_path):  # noqa: F811
    """P5 acceptance 3, the run half: an edit changes the hash and every finished image reads stale."""
    need_readers()
    proj = project_with_images(tmp_path, names=("a", "b"))
    assert P.run_project(proj)["counts"]["done"] == 2
    assert {r["status"] for r in P.project_status(proj)["images"]} == {"done"}
    assert cli_projects.main(["edit", str(proj.root), "--set", "spots.log.min_snr=3", "--json"]) == 0
    st = P.project_status(proj.root)
    assert st["edited"] is True and st["config_hash"] != BASIC_HASH
    for r in st["images"]:
        assert (r["status"], r["reason_code"]) == ("stale", "config_hash")
        assert f"config_hash {BASIC_HASH} != {st['config_hash']}" in r["reason"]
    res = P.run_project(proj.root)
    assert res["counts"] == {"done": 2, "failed": 0, "skipped": 0}
    assert {r["status"] for r in P.project_status(proj)["images"]} == {"done"}


def _log(proj, run_id: str, *, device="cpu", interpreter="/env/a/python", sha="aaaa", end=True, t="2026-10-24T10"):
    runs = proj.root / "runs"
    runs.mkdir(exist_ok=True)
    start = {"event": "run_start", "run_id": run_id, "time": f"{t}:00:00.000+00:00", "via": "cli",
             "interpreter": interpreter, "engine": {"germquant": "0.2.0.dev0", "git_sha": sha},
             "device_requested": "auto", "device_selected": device, "images": []}
    lines = [start]
    if end:
        lines.append({"event": "run_end", "run_id": run_id, "time": f"{t}:00:01.000+00:00",
                      "counts": {"done": 0, "failed": 0, "skipped": 0}, "interrupted": False})
    (runs / f"{run_id}.jsonl").write_bytes("".join(json.dumps(x) + "\n" for x in lines).encode("ascii"))


def test_environment_mix_and_the_runs(project_env, tmp_path):  # noqa: F811
    proj = project_with_images(tmp_path, names=("a",))
    _log(proj, "20261024T100000Z_00000001")
    assert not [w for w in P.project_status(proj)["warnings"] if "environment" in w]
    _log(proj, "20261024T110000Z_00000002", device="cuda", interpreter="/env/b/python", sha="bbbb", end=False,
         t="2026-10-24T11")
    st = P.project_status(proj)
    (mix,) = [w for w in st["warnings"] if "environment" in w]
    for word in ("cpu", "cuda", "/env/a/python", "/env/b/python", "aaaa", "bbbb"):
        assert word in mix
    assert [(r["run_id"], r["finished"], r["interrupted"], r["ended"]) for r in st["runs"]] == [
        ("20261024T100000Z_00000001", True, False, "2026-10-24T10:00:01.000+00:00"),
        ("20261024T110000Z_00000002", False, True, None)]          # no run_end: stopped, or still running
    # a line still being written (no newline yet, even when it already reads as JSON) and a line that is not JSON
    # are left out; a file that is no run log is not read
    log2 = proj.root / "runs" / "20261024T110000Z_00000002.jsonl"
    with open(log2, "ab") as fh:
        fh.write(b"not json\n")
        fh.write(json.dumps({"event": "image_done", "run_id": "20261024T110000Z_00000002", "status": "done",
                             "time": "2026-10-24T11:00:02.000+00:00", "image_id": "a"}).encode("ascii"))
    assert [e["event"] for e in runlog.read_log(log2)] == ["run_start"]
    (proj.root / "runs" / "notes.jsonl").write_bytes(b"not a log\n")
    runs = P.project_status(proj)["runs"]
    assert len(runs) == 2 and runs[1]["counts"] == {"done": 0, "failed": 0, "skipped": 0}


def test_lock_report(project_env, tmp_path):  # noqa: F811
    proj = project_with_images(tmp_path, names=("a",))
    lock = proj.root / "germquant_project.lock"
    now = datetime.now(timezone.utc)
    lock.write_text(json.dumps({"host": "labpc", "pid": 7, "time": now.isoformat(timespec="seconds"),
                                "command": "edit"}), encoding="ascii")
    st = P.project_status(proj)
    assert st["lock"] == {"held": True, "host": "labpc", "pid": 7, "since": now.isoformat(timespec="seconds"),
                          "stale": False}
    assert any("locked by labpc" in w for w in st["warnings"])
    old = (now - timedelta(hours=13)).isoformat(timespec="seconds")
    lock.write_text(json.dumps({"host": "labpc", "pid": 7, "time": old, "command": "edit"}), encoding="ascii")
    st = P.project_status(proj)
    assert st["lock"]["stale"] is True
    assert any("older than 12 h" in w and "unlock" in w for w in st["warnings"]), st["warnings"]


def test_status_json_equals_the_api(project_env, tmp_path, capsys):  # noqa: F811
    need_readers()
    proj = project_with_images(tmp_path, names=("a", "b"))
    P.run_project(proj, images=["a"])
    assert cli.main(["project", "status", str(proj.root), "--json"]) == 0
    printed = json.loads(capsys.readouterr().out)
    assert printed == json.loads(json.dumps(api.project_status(proj.root)))
    assert cli.main(["project", "status", str(proj.root)]) == 0
    text = capsys.readouterr().out
    assert "a" in text and "done" in text and "missing" in text and BASIC_HASH in text
    assert cli.main(["project", "status", str(tmp_path / "nothing"), "--json"]) == 2
    assert json.loads(capsys.readouterr().out)["kind"] == "FileNotFoundError"
    # the status of a project moved to another folder reads its results where they now are
    moved = tmp_path / "moved"
    shutil.copytree(proj.root, moved)
    rows = {r["id"]: r["status"] for r in P.project_status(moved)["images"]}
    assert rows == {"a": "done", "b": "missing"}
    assert log_events(next((moved / "runs").iterdir()))[0]["project"]["root"] == str(proj.root)
