"""Project runs: run_project, its run log, models, devices, time series, collect and export-commands
(docs/PROJECTS.md; design 7.5 to 7.8; PLAN.md 4.4 P5-RUN and 10.6, decisions D12, D15, D17, D28, D33, D47 to D49).

* The run log runs/<run_id>.jsonl holds run_start, one event per image and run_end, with the plugin's envelope and
  keys (napari_germquant/_runner.py, RunLog and execute_run) plus project and shard in run_start and interrupted
  in run_end (D17); one ASCII line per event; done_sha256 on LF-normalised bytes (D48).
* A malformed run id is refused, --via is recorded, the exit codes are 0, 1 and 2, should_stop is honoured
  between images, and a run never writes the manifest.
* GERMQUANT_DEVICE is set from run.device (or --device, auto included) for the run and restored afterwards.
* A project model whose resolved file has another digest (a decoy under $GERMQUANT_MODELS/training/) is refused
  before any pixel is read; a project made on test_pack 0.1.0 still runs after its root moves test_pack to 0.2.0;
  a labels_from project skips the segmentation model's digest check with a warning (D12).
* An image with timepoints goes through process_series; a planar or labels_from recipe is refused for it before
  process_series, with the D9 and D15 messages, and nothing resolves a model or reads a timepoint (D33).
* In a fresh interpreter run_project imports germquant.pipeline (torch) before any reader or backend (D28).
* collect_project uses the recipe's qc.germline_outlier_factor (P5 acceptance 14); export-commands writes the
  commands that run the project.

The helpers at the top are shared with test_project_status.py, test_project_shards.py and test_project_parity.py.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pandas as pd
import pytest
import yaml
from packs_helpers import PYTHON, ROOT, TEST_PACK, copy_pack, pack_env, subprocess_env  # noqa: F401
from project_helpers import (  # noqa: F401
    BASIC_HASH,
    EXAMPLE,
    TEST_PACK_REF,
    new_test_project,
    project_env,
    read_manifest_file,
    write_ome_tiff,
)

import germquant.project as P
from germquant import api, cli, packs
from germquant.project import lock as L
from germquant.project import run as R
from germquant.project import runlog

# The plugin's run log (napari-germquant main 718663f, src/napari_germquant/_runner.py): the envelope of
# RunLog.event and the fields execute_run writes. The engine writes the same keys, plus project and shard in
# run_start and interrupted in run_end (PLAN.md D17); every image event carries the image_done keys.
PLUGIN_ENVELOPE = ["event", "run_id", "time"]
PLUGIN_RUN_START = ["via", "host", "os", "interpreter", "interpreter_source", "engine", "python", "constraints",
                    "device_requested", "device_selected", "config_hash", "config_hash_source", "recipe",
                    "xy_stride", "disabled_stages", "argv", "images", "plugin", "log_location"]
PLUGIN_RUN_START_IMAGE = ["path", "image_id", "out_dir"]
PLUGIN_ENGINE = ["germquant", "git_sha"]
PLUGIN_IMAGE_DONE = ["image_id", "out_dir", "status", "done_sha256", "config_hash", "stage_hashes", "failed_stages",
                     "qc_pass", "segment_device_used", "elapsed_s", "error"]
PLUGIN_RUN_END = ["counts"]
TIME_RE = re.compile(r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d{3}\+00:00$")


# ---- helpers shared with the other project test files ------------------------------------------------------

def example_copy(folder: Path, name: str = "test_pack_synthetic") -> Path:
    """A copy of the test pack's example OME-TIFF (two channels, five nuclei) as ``folder/<name>.ome.tif``."""
    folder.mkdir(parents=True, exist_ok=True)
    dest = folder / f"{name}.ome.tif"
    shutil.copyfile(EXAMPLE, dest)
    return dest


def project_with_images(tmp_path: Path, names=("a", "b"), recipe=None, *, folder: str = "proj"):
    """A test_pack project holding copies of the example under images/, one per name, added by relative path."""
    proj = new_test_project(tmp_path / folder, recipe)
    imgs = [example_copy(proj.root / "images", n) for n in names]
    res = P.add_images(proj, imgs)
    assert res["refused"] == [] and res["added"] == list(names), res
    return proj


def log_events(path) -> list[dict]:
    raw = Path(path).read_bytes()
    raw.decode("ascii")                                     # ASCII only
    lines = raw.decode("ascii").split("\n")
    assert lines[-1] == "", "the log ends with a complete line"
    return [json.loads(line) for line in lines[:-1]]


def marker_lf_sha(out_dir: Path, image_id: str) -> str:
    return hashlib.sha256((out_dir / f"{image_id}__done.json").read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def fake_run_loaded(monkeypatch, seen: list, *, fail=()):
    """api.run_loaded replaced by a stand-in that records each call (image name, out_dir, GERMQUANT_DEVICE) and
    writes nothing; images whose stem is in ``fail`` raise."""
    def fake(image, cfg, out=None, *, xy_stride=None, z_range=None, t=None):
        stem = Path(image).name.split(".")[0]
        seen.append({"image": Path(image).name, "out": Path(out), "device": os.environ.get("GERMQUANT_DEVICE"),
                     "xy_stride": xy_stride, "z_range": z_range})
        if stem in fail:
            raise RuntimeError(f"stand-in failure for {stem}")
        return {"image_id": stem, "n_nuclei": 0, "qc_pass": True, "qc_flags": [], "out_dir": str(out), "tables": {},
                "stages": {}}

    monkeypatch.setattr(api, "run_loaded", fake)


def need_readers():
    """Skip, with the reason, where the OME-TIFF readers are not installed (the readers extra): CI and the
    workstation have them (PLAN.md Q3), so these tests run there."""
    import importlib.util

    missing = [m for m in ("bioio", "bioio_ome_tiff") if importlib.util.find_spec(m) is None]
    if missing:
        pytest.skip(f"the OME-TIFF readers are not installed here ({', '.join(missing)}: the readers extra)")


def need_p2b():
    from germquant import planar

    if not hasattr(planar, "PLANAR_SERIES_REFUSED"):
        pytest.skip("needs P2b's single-plane names (planar.PLANAR_SERIES_REFUSED, one_plane_header): runs once "
                    "feat/t5-p5 has synced with merge #1")


def need_p2c():
    import importlib.util

    if importlib.util.find_spec("germquant.io.labels") is None:
        pytest.skip("needs P2c's labels_from (germquant.io.labels): runs once feat/t5-p5 has synced with merge #2")


def hand_edit_recipe(proj, data: dict) -> None:
    """Write recipe.yaml by hand (origin manual) and record it in the manifest, as a lab's hand edit is kept."""
    from germquant.project import manifest as M

    text = "extends: " + f"pack:test_pack@0.1.0/{proj.pack['recipe']}\n" + yaml.safe_dump(data, sort_keys=True)
    (proj.root / "recipe.yaml").write_bytes(text.encode("ascii"))
    m = M.read_manifest(proj.root)
    m["recipe"]["file"] = "recipe.yaml"
    M.write_manifest(proj.root, m)
    proj.manifest = m


# ---- the run log (P5 acceptance 2, second half; D17, D48) ---------------------------------------------------

def test_run_log_events(project_env, tmp_path):  # noqa: F811
    need_readers()
    proj = project_with_images(tmp_path)
    before = (proj.root / "germquant_project.yaml").read_bytes()
    seen = []
    res = P.run_project(proj, progress=seen.append)
    assert res["exit_code"] == 0 and res["counts"] == {"done": 2, "failed": 0, "skipped": 0}
    assert res["interrupted"] is False and runlog.RUN_ID_RE.fullmatch(res["run_id"])
    log = Path(res["log"])
    assert log == proj.root / "runs" / f"{res['run_id']}.jsonl"
    ev = log_events(log)
    assert [e["event"] for e in ev] == ["run_start", "image_done", "image_done", "run_end"]
    assert seen == ev                                         # progress gets every event as written
    for e in ev:
        assert list(e)[:3] == PLUGIN_ENVELOPE and e["run_id"] == res["run_id"] and TIME_RE.fullmatch(e["time"])
    start, a, b, end = ev
    assert list(start)[3:] == PLUGIN_RUN_START + ["project", "shard"]
    assert start["via"] == "cli" and start["shard"] is None and start["config_hash"] == BASIC_HASH
    assert list(start["engine"]) == PLUGIN_ENGINE and start["engine"]["germquant"] == api.__version__
    assert start["project"] == {"root": str(proj.root), "name": "proj", "id": proj.manifest["id"]}
    assert start["recipe"] == {"kind": "project", "pack": TEST_PACK_REF, "recipe": "basic", "edited": False,
                               "file": None}
    assert start["xy_stride"] == 1 and start["disabled_stages"] == [] and start["argv"] is None
    assert [sorted(i) for i in start["images"]] == [sorted(PLUGIN_RUN_START_IMAGE)] * 2
    assert [i["image_id"] for i in start["images"]] == ["a", "b"]
    assert start["images"][0]["out_dir"] == str(proj.root / "results" / "a")
    for e, iid in ((a, "a"), (b, "b")):
        out = proj.root / "results" / iid
        marker = json.loads((out / f"{iid}__done.json").read_text(encoding="utf-8"))
        assert list(e)[3:] == PLUGIN_IMAGE_DONE
        assert e["image_id"] == iid and e["out_dir"] == str(out) and e["status"] == "done" and e["error"] is None
        assert e["done_sha256"] == marker_lf_sha(out, iid)
        assert e["config_hash"] == marker["config_hash"] == BASIC_HASH
        assert e["stage_hashes"] == marker["stage_hashes"] and e["failed_stages"] == []
        assert e["qc_pass"] is True and e["elapsed_s"] >= 0
    assert list(end)[3:] == PLUGIN_RUN_END + ["interrupted"]
    assert end["counts"] == {"done": 2, "failed": 0, "skipped": 0} and end["interrupted"] is False
    # the run wrote no manifest and left no lock; the images' results are where design 7.1 says
    assert (proj.root / "germquant_project.yaml").read_bytes() == before
    assert not (proj.root / "germquant_project.lock").exists()

    # the second run resumes: every image skipped, with its marker's digest and record
    res2 = P.run_project(proj)
    ev2 = log_events(res2["log"])
    assert res2["exit_code"] == 0 and res2["counts"] == {"done": 0, "failed": 0, "skipped": 2}
    assert [e["event"] for e in ev2] == ["run_start", "image_skipped", "image_skipped", "run_end"]
    for e, old in zip(ev2[1:3], (a, b)):
        assert list(e)[3:] == PLUGIN_IMAGE_DONE and e["status"] == "skipped"
        assert {k: e[k] for k in ("done_sha256", "config_hash", "stage_hashes", "qc_pass")} == \
               {k: old[k] for k in ("done_sha256", "config_hash", "stage_hashes", "qc_pass")}
    # --no-resume reruns them
    res3 = P.run_project(proj, resume=False)
    assert res3["counts"] == {"done": 2, "failed": 0, "skipped": 0}
    assert sorted(p.name for p in (proj.root / "runs").iterdir()) == sorted(
        f"{r['run_id']}.jsonl" for r in (res, res2, res3))
    assert runlog.read_logs(proj.root) == sorted(ev + ev2 + log_events(res3["log"]), key=lambda e: e["time"])


def test_run_id_via_and_exit_codes(project_env, tmp_path, capsys):  # noqa: F811
    need_readers()
    proj = project_with_images(tmp_path, names=("a", "b"))
    runs = proj.root / "runs"
    # 2: refused before any image, nothing written
    for bad in ("20261024T120000Z", "20261024T120000Z_ABCDEF12", "2026-10-24_0123abcd", "x"):
        with pytest.raises(P.ProjectError, match="run id"):
            P.run_project(proj, run_id=bad)
        assert cli.main(["project", "run", str(proj.root), "--run-id", bad]) == 2
    err = capsys.readouterr().err
    assert "REFUSED:" in err and "run id" in err and not runs.exists()
    with pytest.raises(SystemExit) as exc:
        cli.main(["project", "run", str(proj.root), "--via", "gui"])
    assert exc.value.code == 2
    with pytest.raises(P.ProjectError, match="via"):
        P.run_project(proj, via="gui")
    assert cli.main(["project", "run", str(proj.root), "--images", "nope", "--json"]) == 2
    out = json.loads(capsys.readouterr().out)
    assert out["kind"] == "ProjectError" and "nope" in out["error"] and not runs.exists()
    # --run-id and --via are recorded; the same run id again is refused and the first log is kept
    rid = "20261024T120000Z_0123abcd"
    assert cli.main(["project", "run", str(proj.root), "--run-id", rid, "--via", "snakemake", "--images", "a",
                     "--json"]) == 0
    out = json.loads(capsys.readouterr().out)
    assert out["run_id"] == rid and out["exit_code"] == 0 and out["counts"] == {"done": 1, "failed": 0, "skipped": 0}
    first = (runs / f"{rid}.jsonl").read_bytes()
    start = log_events(runs / f"{rid}.jsonl")[0]
    assert start["via"] == "snakemake" and [i["image_id"] for i in start["images"]] == ["a"]
    assert start["argv"] == ["germquant", "project", "run", str(proj.root), "--images", "a", "--run-id", rid,
                             "--via", "snakemake", "--json"]
    assert cli.main(["project", "run", str(proj.root), "--run-id", rid]) == 2
    assert "already" in capsys.readouterr().err and (runs / f"{rid}.jsonl").read_bytes() == first
    # 1: some image failed (its file is no image any more); the others still run
    (proj.root / "images" / "b.ome.tif").write_bytes(b"not a TIFF any more")
    res = P.run_project(proj, resume=False)
    assert res["exit_code"] == 1 and res["counts"] == {"done": 1, "failed": 1, "skipped": 0}
    failed = [e for e in log_events(res["log"]) if e["event"] == "image_failed"]
    assert [e["image_id"] for e in failed] == ["b"] and failed[0]["status"] == "failed"
    assert list(failed[0])[3:] == PLUGIN_IMAGE_DONE and failed[0]["error"] and failed[0]["done_sha256"] is None
    assert cli.main(["project", "run", str(proj.root), "--images", "b"]) == 1
    # an image that cannot be reached says so (the NAS or VPN)
    os.remove(proj.root / "images" / "b.ome.tif")
    res = P.run_project(proj, images=["b"])
    (e,) = [e for e in log_events(res["log"]) if e["event"] == "image_failed"]
    assert "not reachable" in e["error"] and res["exit_code"] == 1
    # 2: a folder that is no project
    assert cli.main(["project", "run", str(tmp_path / "nothing")]) == 2


def test_lock_held_refuses_and_a_stale_lock_warns(project_env, tmp_path, monkeypatch):  # noqa: F811
    proj = project_with_images(tmp_path, names=("a",))
    with L.project_lock(proj.root, "edit"):
        with pytest.raises(P.ProjectLockedError, match="locked by"):
            P.run_project(proj)
    assert not (proj.root / "runs").exists()
    old = (datetime.now(timezone.utc) - timedelta(hours=13)).isoformat(timespec="seconds")
    (proj.root / "germquant_project.lock").write_text(json.dumps({"host": "gone", "pid": 4242, "time": old,
                                                                  "command": "edit"}), encoding="ascii")
    seen = []
    fake_run_loaded(monkeypatch, seen)
    res = P.run_project(proj)
    assert res["exit_code"] == 0 and [s["image"] for s in seen] == ["a.ome.tif"]
    assert any("older than 12 h" in w and "unlock" in w for w in res["warnings"]), res["warnings"]
    assert (proj.root / "germquant_project.lock").exists()            # reported, never removed


def test_an_unusable_pack_is_refused(project_env, tmp_path, monkeypatch, capsys):  # noqa: F811
    """I2: a pack that `pack check` finds unusable (here the project's copy no longer matches its SHA256SUMS) is a
    refusal before any image, with the fix."""
    proj = project_with_images(tmp_path, names=("a",))
    readme = proj.root / "packs" / "test_pack" / "0.1.0" / "README.md"
    readme.write_bytes(readme.read_bytes() + b"\nedited by hand\n")
    seen = []
    fake_run_loaded(monkeypatch, seen)
    with pytest.raises(P.ProjectError, match="SHA256SUMS") as exc:
        P.run_project(proj)
    assert "not usable" in str(exc.value) and "fix:" in str(exc.value)
    assert seen == [] and not (proj.root / "runs").exists()
    assert cli.main(["project", "run", str(proj.root)]) == 2
    assert "REFUSED: pack test_pack@0.1.0 recipe basic is not usable here" in capsys.readouterr().err


def test_should_stop_between_images(project_env, tmp_path, monkeypatch):  # noqa: F811
    proj = project_with_images(tmp_path, names=("a", "b", "c"))
    seen = []
    fake_run_loaded(monkeypatch, seen)
    asked = []

    def stop():
        asked.append(len(seen))
        return len(seen) >= 1

    res = P.run_project(proj, should_stop=stop)
    assert [s["image"] for s in seen] == ["a.ome.tif"] and asked[:2] == [0, 1]
    assert res["interrupted"] is True and res["exit_code"] == 0
    ev = log_events(res["log"])
    assert [e["event"] for e in ev] == ["run_start", "image_done", "run_end"]
    assert ev[-1]["interrupted"] is True and ev[-1]["counts"] == {"done": 1, "failed": 0, "skipped": 0}


def test_device_is_set_for_the_run_and_restored(project_env, tmp_path, monkeypatch):  # noqa: F811
    proj = project_with_images(tmp_path, names=("a",))
    seen = []
    fake_run_loaded(monkeypatch, seen)
    monkeypatch.delenv("GERMQUANT_DEVICE", raising=False)
    P.run_project(proj, resume=False)
    assert seen[-1]["device"] is None and "GERMQUANT_DEVICE" not in os.environ   # no default sets it
    P.save_project(proj, device="cpu")
    res = P.run_project(proj, resume=False)
    assert seen[-1]["device"] == "cpu" and "GERMQUANT_DEVICE" not in os.environ
    assert log_events(res["log"])[0]["device_requested"] == "cpu"
    assert log_events(res["log"])[0]["device_selected"] == "cpu"
    monkeypatch.setenv("GERMQUANT_DEVICE", "mps")
    for dev in ("auto", "cuda", "cpu"):                                   # the computer's choice wins (D47)
        P.run_project(proj, resume=False, device=dev)
        assert seen[-1]["device"] == dev and os.environ["GERMQUANT_DEVICE"] == "mps"
    assert cli.main(["project", "run", str(proj.root), "--no-resume", "--device", "auto"]) == 0
    assert seen[-1]["device"] == "auto" and os.environ["GERMQUANT_DEVICE"] == "mps"
    with pytest.raises(P.ProjectError, match="device"):
        P.run_project(proj, device="tpu")
    # a run that raises still restores the variable
    fake_run_loaded(monkeypatch, seen, fail=("a",))
    assert P.run_project(proj, resume=False, device="cpu")["exit_code"] == 1
    assert os.environ["GERMQUANT_DEVICE"] == "mps"
    # the stride and z range of the manifest reach every image
    P.save_project(proj, xy_stride=2, z_range=(1, 5))
    fake_run_loaded(monkeypatch, seen)
    P.run_project(proj, resume=False)
    assert (seen[-1]["xy_stride"], seen[-1]["z_range"]) == (2, (1, 5))


# ---- models (P5 acceptance 13, D12) -----------------------------------------------------------------------

def _project_model(proj, content=b"fine-tuned weights stand-in\n" * 8):
    """training/r01/models/m1/m1 with its card, put into the recipe with use_model."""
    P.new_round(proj)
    mdir = proj.root / "training" / "r01" / "models" / "m1"
    mdir.mkdir(parents=True, exist_ok=True)
    model = mdir / "m1"
    model.write_bytes(content)
    sha = hashlib.sha256(content).hexdigest()
    (mdir / "m1.card.json").write_text(json.dumps({"name": "m1", "model_file": str(model), "model_sha256": sha,
                                                   "model_size": len(content), "base_model": {"name": "cpsam"}}),
                                       encoding="ascii")
    P.use_model(proj, "training/r01/models/m1/m1")
    return model, sha


def test_decoy_model_refused(project_env, tmp_path, monkeypatch, capsys):  # noqa: F811
    proj = project_with_images(tmp_path, names=("a",), recipe="with_model")
    model, sha = _project_model(proj)
    decoy = project_env.models / "training" / "r01" / "models" / "m1" / "m1"
    decoy.parent.mkdir(parents=True)
    decoy.write_bytes(b"a decoy in the shared model store\n")
    from germquant import pipeline

    def no_pixels(*a, **k):
        raise AssertionError("a pixel was read")

    monkeypatch.setattr(pipeline, "_read_input", no_pixels)
    with pytest.raises(P.ProjectError) as exc:
        P.run_project(proj)
    msg = str(exc.value)
    assert str(decoy) in msg and sha[:12] in msg and "GERMQUANT_MODELS" in msg
    assert not (proj.root / "runs").exists()
    assert cli.main(["project", "run", str(proj.root), "--json"]) == 2
    out = json.loads(capsys.readouterr().out)
    assert out["kind"] == "ProjectError" and str(decoy) in out["error"]
    st = P.project_status(proj)
    (rec,) = [m for m in st["models"] if m["key"] == "segmentation.nuclei.cellpose_model"]
    assert rec["ref"] == "training/r01/models/m1/m1" and rec["sha256"] == sha and rec["needed"] is True
    assert rec["present"] is True and rec["verified"] is False and str(decoy) in rec["mismatch"]
    assert any("training/r01/models/m1/m1" in w for w in st["warnings"])
    # without the decoy the digest check passes (the pack's tinynet is no longer needed, so its missing blob is
    # no refusal either), and the run starts
    decoy.unlink()
    seen = []
    fake_run_loaded(monkeypatch, seen)
    res = P.run_project(proj)
    assert res["exit_code"] == 0 and [s["image"] for s in seen] == ["a.ome.tif"]
    st = P.project_status(proj)
    (rec,) = [m for m in st["models"] if m["key"] == "segmentation.nuclei.cellpose_model"]
    assert rec["verified"] is True and rec["mismatch"] is None
    (tiny,) = [m for m in st["models"] if m["ref"] == "tinynet@f5055f8a"]
    assert tiny["needed"] is False and tiny["present"] is False
    # the model file itself changed after use-model: refused the same way
    model.write_bytes(b"retrained in place\n")
    with pytest.raises(P.ProjectError, match="sha256"):
        P.run_project(proj, resume=False)


def test_a_model_that_cannot_be_resolved_is_refused(project_env, tmp_path, monkeypatch):  # noqa: F811
    """A recipe whose segmentation model the engine cannot resolve here (neither a file, nor registered, nor a
    Cellpose built-in) is refused before any image, instead of failing on each."""
    proj = project_with_images(tmp_path, names=("a",), recipe="with_model")
    hand_edit_recipe(proj, {"segmentation": {"nuclei": {"cellpose_model": "no_such_model"}}})
    seen = []
    fake_run_loaded(monkeypatch, seen)
    with pytest.raises(P.ProjectError, match="cannot be resolved") as exc:
        P.run_project(proj)
    assert "no_such_model" in str(exc.value) and seen == [] and not (proj.root / "runs").exists()
    # process_image resolves the name under the classical method too (it only skips hashing it), so the same holds
    hand_edit_recipe(proj, {"segmentation": {"nuclei": {"method": "classical", "cellpose_model": "no_such_model"}}})
    with pytest.raises(P.ProjectError, match="cannot be resolved"):
        P.run_project(proj)


def test_vendored_pack_survives_root_upgrade(project_env, tmp_path, monkeypatch):  # noqa: F811
    need_readers()
    from germquant import config as C

    root = tmp_path / "engine"
    (root / "src" / "germquant").mkdir(parents=True)
    copy_pack(TEST_PACK, root / "packs" / "test_pack", rebuild=False)
    monkeypatch.setattr(C, "builtin_roots", lambda: [root.resolve()])
    monkeypatch.delenv("GERMQUANT_PACKS", raising=False)
    packs.refresh()
    assert packs.find_pack("test_pack").source == "builtin"
    proj = new_test_project(tmp_path / "p")
    P.add_images(proj, [example_copy(proj.root / "images")])

    def to_020(d):
        d["version"] = "0.2.0"

    shutil.rmtree(root / "packs" / "test_pack")
    copy_pack(TEST_PACK, root / "packs" / "test_pack", edit=to_020)
    packs.refresh()
    assert [p.version for p in packs.discover().packs if p.id == "test_pack"] == ["0.2.0"]
    res = P.run_project(P.open_project(proj.root))
    assert res["exit_code"] == 0 and res["counts"]["done"] == 1
    man = json.loads((proj.root / "results" / "test_pack_synthetic" / "run_manifest.json").read_text("utf-8"))
    assert man["config_hash"] == BASIC_HASH
    assert man["pack"]["version"] == "0.1.0" and man["pack"]["source"] == "project" and man["pack"]["edited"] is False


def test_labels_project_skips_the_segmentation_model(project_env, tmp_path, monkeypatch):  # noqa: F811
    """D12: a project whose recipe sets labels_from skips the segmentation model's digest check (with a warning)
    and status reports that model as not needed; the run itself reads the labels files (P2c)."""
    from germquant.project.recipe import resolve_pack

    proj = project_with_images(tmp_path, names=("a",), recipe="with_model")
    _project_model(proj)
    decoy = project_env.models / "training" / "r01" / "models" / "m1" / "m1"
    decoy.parent.mkdir(parents=True)
    decoy.write_bytes(b"a decoy in the shared model store\n")
    with pytest.raises(P.ProjectError, match="sha256"):
        P.run_project(proj)
    recipe = yaml.safe_load((proj.root / "recipe.yaml").read_text(encoding="ascii"))
    recipe.pop("extends")
    recipe["segmentation"]["nuclei"]["labels_from"] = "{image_dir}/{image_stem}_labels.tif"
    hand_edit_recipe(proj, recipe)
    cfg = P.project_config(proj)
    pack, rec = resolve_pack(proj.root, proj.manifest)
    records, warnings = R.model_records(proj.root, proj.manifest, cfg, pack, rec)
    (seg,) = [m for m in records if m["key"] == "segmentation.nuclei.cellpose_model"]
    assert seg["needed"] is False and seg["mismatch"] is not None and str(decoy) in seg["mismatch"]
    assert R.model_refusal(records, pack) is None
    assert any("labels_from" in w and "not checked" in w for w in warnings), warnings
    need_p2c()
    need_readers()
    import numpy as np
    import tifffile

    from germquant.io.labels import labels_path

    st = P.project_status(proj)
    (rec,) = [m for m in st["models"] if m["key"] == "segmentation.nuclei.cellpose_model"]
    assert rec["needed"] is False and any("labels_from" in w and "not checked" in w for w in st["warnings"])
    lab = np.zeros((12, 64, 96), dtype=np.int32)
    lab[:, 10:20, 10:20] = 1
    tifffile.imwrite(labels_path(cfg, proj.root / "images" / "a.ome.tif"), lab)
    res = P.run_project(proj)
    assert res["exit_code"] == 0 and res["counts"]["done"] == 1
    assert any("labels_from" in w for w in res["warnings"])
    marker = json.loads((proj.root / "results" / "a" / "a__done.json").read_text(encoding="utf-8"))
    assert marker["config_hash"] == cfg.hash
    # the resume uses the per-image labels hash (D46): a second run skips the image, and status says done
    assert P.run_project(proj)["counts"] == {"done": 0, "failed": 0, "skipped": 1}
    assert [r["status"] for r in P.project_status(proj)["images"]] == ["done"]


# ---- time series (D9, D15, D33) ------------------------------------------------------------------------

def _series_image(folder: Path, name: str = "movie", t: int = 2) -> Path:
    """An OME-TIFF with ``t`` timepoints of the example field (axes TCZYX)."""
    import numpy as np
    import tifffile

    with tifffile.TiffFile(EXAMPLE) as tf:
        data = tf.series[0].asarray()                       # (C, Z, Y, X)
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"{name}.ome.tif"
    meta = {"axes": "TCZYX", "Channel": {"Name": ["DAPI", "foci"]}, "PhysicalSizeZ": 0.5, "PhysicalSizeY": 0.2,
            "PhysicalSizeX": 0.2}
    tifffile.imwrite(path, np.stack([data] * t), ome=True, metadata=meta)
    return path


def _spy_series(monkeypatch):
    from germquant import pipeline, timelapse

    calls = {"series": [], "model": [], "read": []}
    real_series = timelapse.process_series

    def series(*a, **k):
        calls["series"].append((Path(a[0]).name, k.get("timepoints")))
        return real_series(*a, **k)

    real_model, real_read = pipeline.resolve_model_path, pipeline._read_input
    monkeypatch.setattr(timelapse, "process_series", series)
    monkeypatch.setattr(pipeline, "resolve_model_path", lambda *a, **k: calls["model"].append(1) or real_model(*a, **k))
    monkeypatch.setattr(pipeline, "_read_input",
                        lambda *a, **k: calls["read"].append(k.get("t")) or real_read(*a, **k))
    return calls


def test_an_image_with_timepoints_goes_through_process_series(project_env, tmp_path, monkeypatch):  # noqa: F811
    need_readers()
    proj = new_test_project(tmp_path / "proj")
    P.add_images(proj, [_series_image(proj.root / "images")])
    assert proj.images[0]["timepoints"] == 2
    calls = _spy_series(monkeypatch)
    res = P.run_project(proj)
    assert res["exit_code"] == 0 and calls["series"] == [("movie.ome.tif", "all")] and calls["read"] == [0, 1]
    (e,) = [e for e in log_events(res["log"]) if e["event"].startswith("image_")]
    assert e["event"] == "image_done" and e["image_id"] == "movie" and e["config_hash"] == BASIC_HASH
    out = proj.root / "results" / "movie"
    assert (out / "movie__series.json").is_file() and (out / "t001" / "movie_t001__done.json").is_file()
    st = P.project_status(proj)
    assert [(i["id"], i["status"]) for i in st["images"]] == [("movie", "done")]
    res = P.run_project(proj)
    assert res["counts"] == {"done": 0, "failed": 0, "skipped": 1}


def test_series_refused_for_planar_recipes(project_env, tmp_path, monkeypatch):  # noqa: F811
    need_p2b()
    from germquant import planar

    proj = new_test_project(tmp_path / "proj")
    P.add_images(proj, [_series_image(proj.root / "images"), example_copy(proj.root / "images", "fixed")])
    for value in (True, "auto"):
        hand_edit_recipe(proj, {"io": {"planar": value}})
        calls = _spy_series(monkeypatch)
        if value == "auto":                    # a one-plane header for the series, so auto runs it as a plane
            monkeypatch.setattr(planar, "one_plane_header", lambda meta: True)
        res = P.run_project(proj, images=["movie"])
        (e,) = [e for e in log_events(res["log"]) if e["event"].startswith("image_")]
        assert e["event"] == "image_failed" and planar.PLANAR_SERIES_REFUSED in e["error"]
        assert res["exit_code"] == 1 and calls == {"series": [], "model": [], "read": []}


def test_a_series_whose_header_cannot_be_read_fails_that_image_only(project_env, tmp_path, monkeypatch):  # noqa: F811
    """io.planar auto reads a series' header to decide; a series file that cannot be reached, or whose header
    cannot be read, fails as that image, and the next image still runs (run.py's docstring, step 3)."""
    need_p2b()
    proj = new_test_project(tmp_path / "proj")
    series = _series_image(proj.root / "images")
    P.add_images(proj, [series, example_copy(proj.root / "images", "fixed")])
    hand_edit_recipe(proj, {"io": {"planar": "auto"}})
    seen = []
    fake_run_loaded(monkeypatch, seen)
    for broken in ("removed", "no TIFF"):
        if broken == "removed":
            os.remove(series)
        else:
            series.write_bytes(b"not a TIFF any more")
        res = P.run_project(proj, resume=False)
        assert res["exit_code"] == 1 and res["counts"] == {"done": 1, "failed": 1, "skipped": 0}, res
        assert [s["image"] for s in seen] == ["fixed.ome.tif"] and res["interrupted"] is False
        seen.clear()
        ev = log_events(res["log"])
        assert ev[-1]["event"] == "run_end" and ev[-1]["interrupted"] is False
        (e,) = [e for e in ev if e["event"] == "image_failed"]
        assert e["image_id"] == "movie" and e["error"]
        assert ("not reachable" in e["error"]) == (broken == "removed"), e["error"]


def test_series_refused_for_labels_recipes(project_env, tmp_path, monkeypatch):  # noqa: F811
    need_p2c()
    from germquant.io.labels import SERIES_REFUSAL

    proj = new_test_project(tmp_path / "proj")
    P.add_images(proj, [_series_image(proj.root / "images")])
    hand_edit_recipe(proj, {"segmentation": {"nuclei": {"labels_from": "{image_dir}/{image_stem}_labels.tif"}}})
    calls = _spy_series(monkeypatch)
    res = P.run_project(proj)
    (e,) = [e for e in log_events(res["log"]) if e["event"].startswith("image_")]
    assert e["event"] == "image_failed" and SERIES_REFUSAL in e["error"]
    assert res["exit_code"] == 1 and calls == {"series": [], "model": [], "read": []}


# ---- imports (D28, R5) -------------------------------------------------------------------------------

def test_run_project_imports_the_pipeline_first(project_env, tmp_path):  # noqa: F811
    need_readers()
    import importlib.util

    if importlib.util.find_spec("torch") is None:
        pytest.skip("torch is not installed here")
    proj = project_with_images(tmp_path, names=("a",))
    # the readers as germquant.io offers them (pipeline binds these names when it is imported) and the backend
    # lookups as the segment and spots stages call them (segment.nuclei and spots.detect bind get_backend)
    code = ("import json, sys\n"
            "import germquant.project as P\n"
            "from germquant import io as gio\n"
            "from germquant.segment import nuclei as SN\n"
            "from germquant.spots import detect as SD\n"
            "assert 'torch' not in sys.modules and 'germquant.pipeline' not in sys.modules\n"
            "calls = []\n"
            "def wrapped(owner, name):\n"
            "    orig = getattr(owner, name)\n"
            "    def w(*a, **k):\n"
            "        calls.append([name, 'germquant.pipeline' in sys.modules, 'torch' in sys.modules])\n"
            "        return orig(*a, **k)\n"
            "    setattr(owner, name, w)\n"
            "for name in ('read_image', 'read_image_metadata'):\n"
            "    wrapped(gio, name)\n"
            "for mod in (SN, SD):\n"
            "    wrapped(mod, 'get_backend')\n"
            "res = P.run_project(sys.argv[1])\n"
            "order = list(sys.modules)\n"
            "print(json.dumps({'exit_code': res['exit_code'], 'calls': calls,\n"
            "                  'pos': {m: order.index(m) for m in ('torch', 'cupy', 'spotmax', 'cellpose', 'dask')\n"
            "                          if m in sys.modules}}))\n")
    r = subprocess.run([PYTHON, "-c", code, str(proj.root)], capture_output=True, text=True, timeout=600,
                       env=subprocess_env(GERMQUANT_PACKS=TEST_PACK.parent), cwd=ROOT, encoding="utf-8",
                       errors="replace")
    assert r.returncode == 0, r.stderr[-3000:]
    out = json.loads(r.stdout.strip().splitlines()[-1])
    assert out["exit_code"] == 0
    names = {c[0] for c in out["calls"]}
    assert {"read_image", "read_image_metadata", "get_backend"} <= names, out["calls"]
    assert all(pipe and torch for _, pipe, torch in out["calls"]), out["calls"]
    pos = out["pos"]
    assert "torch" in pos and all(pos["torch"] < p for m, p in pos.items() if m != "torch"), pos


# ---- collect and export-commands (P5 acceptance 14) -----------------------------------------------------------

def _fake_results(results: Path, ids_counts: dict) -> None:
    """Finished image folders holding only what collect reads: an image summary and a completion marker."""
    from germquant import provenance

    for iid, n in ids_counts.items():
        d = results / iid
        d.mkdir(parents=True)
        pd.DataFrame([{"image_id": iid, "n_nuclei": n, "n_germline_nuclei": n, "qc_pass": True, "qc_flags": ""}]
                     ).to_csv(d / f"{iid}__image_summary.csv", index=False)
        provenance.write_done_marker(d, iid, config_hash=BASIC_HASH, enabled_stages=[], stage_hashes={})


def test_collect_project_uses_recipe_outlier_factor(project_env, tmp_path, capsys):  # noqa: F811
    counts = {"a": 10, "b": 11, "c": 9, "d": 10, "e": 60}           # e: six times the median, a dense field
    proj = project_with_images(tmp_path, names=tuple(counts))
    _fake_results(proj.root / "results", counts)
    _fake_results(proj.root / "results", {"stray": 10})              # a folder of no image of the project
    rec = P.collect_project(proj)
    assert rec["outlier_factor"] == 1.0e9 and rec["counts"]["batch_summary"] == 5
    summary = pd.read_csv(proj.root / "results" / "batch_summary.csv")
    assert sorted(summary["image_id"]) == sorted(counts)
    assert not summary["qc_flags"].fillna("").str.contains("germline_count_outlier").any()
    assert cli.main(["project", "collect", str(proj.root), "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["counts"] == rec["counts"]
    # the plain `germquant collect` of the same folder keeps 1.8
    assert cli.main(["collect", str(proj.root / "results")]) == 0
    plain = pd.read_csv(proj.root / "results" / "batch_summary.csv").set_index("image_id")
    assert "qc:germline_count_outlier_60" in plain.loc["e", "qc_flags"]
    # nothing to collect yet: refused
    empty = new_test_project(tmp_path / "empty")
    with pytest.raises(P.ProjectError, match="no results"):
        P.collect_project(empty)


def test_export_commands(project_env, tmp_path, capsys):  # noqa: F811
    proj = project_with_images(tmp_path, names=("a",), folder="my study")
    root = str(proj.root)
    sh = P.export_commands(proj)
    lines = [ln for ln in sh.splitlines() if ln and not ln.startswith("#")]
    assert sh.startswith("#!/bin/sh\n") and len(lines) == 2
    exe, q = f"{shlex.quote(sys.executable)} -m germquant", shlex.quote(root)
    assert q == f"'{root}'"                                           # the folder name holds a space
    assert lines[-2:] == [f"{exe} project run {q}", f"{exe} project collect {q}"]
    assert f"config_hash {BASIC_HASH}" in sh and "test_pack@0.1.0" in sh
    sh3 = P.export_commands(proj, shards=3)
    runs = [ln for ln in sh3.splitlines() if " project run " in ln and not ln.startswith("#")]
    assert runs == [f"{exe} project run {q} --shard {k}/3" for k in (1, 2, 3)]
    assert "SLURM_ARRAY_TASK_ID" in sh3
    ps = P.export_commands(proj, shell="ps1", shards=2)
    runs = [ln for ln in ps.splitlines() if " project run " in ln and not ln.startswith("#")]
    assert runs == [f"& '{sys.executable}' -m germquant project run '{root}' --shard {k}/2" for k in (1, 2)]
    # a recorded engine interpreter for the recipe is the one the commands use
    P.save_project(proj, interpreters={"basic": {"python": "/opt/lab env/bin/python", "germquant": "0.2.0",
                                                 "git_sha": "abc"}})
    assert "'/opt/lab env/bin/python' -m germquant project run" in P.export_commands(proj)
    for bad in ({"shell": "bat"}, {"shards": 0}):
        with pytest.raises(P.ProjectError):
            P.export_commands(proj, **bad)
    out = tmp_path / "run all.ps1"
    assert cli.main(["project", "export-commands", root, "--shell", "ps1", "--out", str(out)]) == 0
    assert capsys.readouterr().out.strip() == f"wrote {out}"
    data = out.read_bytes()
    assert data.startswith(b"\xef\xbb\xbf") and data[3:].decode("utf-8") == P.export_commands(proj, shell="ps1")
    assert cli.main(["project", "export-commands", root]) == 0
    assert capsys.readouterr().out == P.export_commands(proj)
    assert cli.main(["project", "export-commands", root, "--shards", "2", "--json"]) == 0
    assert json.loads(capsys.readouterr().out) == {"text": P.export_commands(proj, shards=2), "out": None,
                                                   "shell": "sh", "shards": 2}
    # nothing the project records can end a comment line and run as a command: a name with line breaks (also a
    # lone CR, which PowerShell reads as one) stays one comment, written with escapes
    odd = new_test_project(tmp_path / "odd", name="study\necho INJECTED\rWrite-Host INJECTED\u2028x")
    for shell in ("sh", "ps1"):
        text = P.export_commands(odd, shell=shell)
        assert all(ln.startswith("#") or " -m germquant project " in ln for ln in text.splitlines()), text
        assert "study\\necho INJECTED\\rWrite-Host INJECTED\\u2028x" in text
    # PowerShell also ends a single-quoted string at a typographic single quote (U+2018 to U+201B): a folder named
    # with one, such as Ryan's study typed with a curly apostrophe, has it doubled, as the ASCII quote is
    for q in (chr(0x2018), chr(0x2019), chr(0x201A), chr(0x201B), "'"):
        curly = new_test_project(tmp_path / f"Ryan{q}s study {ord(q)}")
        runs = [ln for ln in P.export_commands(curly, shell="ps1").splitlines() if not ln.startswith("#")]
        doubled = str(curly.root).replace(q, q * 2)
        assert runs == [f"& '{sys.executable}' -m germquant project run '{doubled}'",
                        f"& '{sys.executable}' -m germquant project collect '{doubled}'"], runs


# ---- the api re-exports (PLAN.md section 9.1) ----------------------------------------------------------------

PROJECT_NAMES = ["Project", "ProjectError", "ProjectLockedError", "new_project", "open_project", "save_project",
                 "add_images", "read_header", "project_config", "preview_roles", "write_recipe", "revert_recipe",
                 "use_model", "run_project", "project_status", "collect_project", "export_commands", "break_lock",
                 "list_rounds", "new_round"]
# the names of section 9.1 that other P5 units bring: (name, module, unit); each joins api once its unit is merged
LATER_NAMES = [("tune_spots", "germquant.tune", "P5-TUNE"), ("apply_tuning", "germquant.tune", "P5-TUNE"),
               ("import_labels", "germquant.import_labels", "P5-TRAIN"),
               ("SampleRow", "germquant.io.sample_sheet", "P5-SHEET"),
               ("SampleSheet", "germquant.io.sample_sheet", "P5-SHEET"),
               ("read_sample_sheet", "germquant.io.sample_sheet", "P5-SHEET"),
               ("write_sample_sheet", "germquant.io.sample_sheet", "P5-SHEET")]


def test_api_re_exports_the_project_names():
    from germquant.project import headers

    for name in PROJECT_NAMES:
        assert name in api.__all__, name
        assert getattr(api, name) is getattr(headers if name == "read_header" else P, name), name
    assert len(set(api.__all__)) == len(api.__all__)


@pytest.mark.parametrize("name,module,unit", LATER_NAMES)
def test_api_re_exports_the_other_p5_names(name, module, unit):
    import importlib
    import importlib.util

    if importlib.util.find_spec(module) is None:
        pytest.skip(f"{module} comes with {unit}, which is not merged into this branch yet")
    assert name in api.__all__ and getattr(api, name) is getattr(importlib.import_module(module), name)
