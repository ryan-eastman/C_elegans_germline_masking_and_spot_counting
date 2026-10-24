"""The project manifest, its lock, new, open and add, and the light header reads (docs/PROJECTS.md; design 7.1,
7.2, 7.8; PLAN.md 4.4 P5-PRJ, decisions D18, D27, D28; P5 acceptance 7 and 13).

* The manifest round-trips every key but the saved stamp; a non-ASCII project name is written as an escape; a
  failure mid-write leaves the old file (PermissionError is retried five times over 2 s).
* The lock is an exclusive create; a second writer is refused and the holder is named; a lock older than 12 h is
  reported; `project unlock` works only when given the holder's host and pid.
* new_project vendors the pinned pack version and writes no recipe.yaml; open after a pack version change warns
  and names the change; a project made from a fake built-in root keeps resolving its own copy after that root
  moves on.
* add_images refuses multi-scene files, nd2 files with P > 1, duplicate ids and unreadable headers, and still
  adds the rest; a folder is discovered with io.formats.SUPPORTED_SUFFIXES when the recipe sets no
  io.input_suffixes, skipping germquant outputs and the project's own folders.
* project.headers reads sizes, the OME image count, channel names, single plane and the yx spacing without
  bioio, and agrees with the run's reader where the readers are installed.
"""
from __future__ import annotations

import importlib.util
import json
import os
import socket
import uuid
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest
import yaml
from packs_helpers import FIXTURE_PACKS, ROOT, TEST_PACK, copy_pack, pack_env  # noqa: F401
from project_helpers import (
    BASIC_HASH,
    EXAMPLE,
    GAMMA,
    fake_nd2_headers,
    last_json,
    new_test_project,
    project_env,  # noqa: F401
    read_manifest_file,
    write_ome_tiff,
)

import germquant
from germquant import cli_projects, packs
from germquant import config as C
from germquant.io.ome_header import read_tiff_header
from germquant.project import (
    ProjectError,
    ProjectLockedError,
    add_images,
    break_lock,
    new_project,
    open_project,
    project_config,
    save_project,
    write_recipe,
)
from germquant.project import lock as L
from germquant.project import manifest as M
from germquant.project.headers import read_header

STAMP_KEYS = {"at", "germquant", "git_sha", "plugin"}


def _without_saved(m: dict) -> dict:
    return {k: v for k, v in m.items() if k != "saved"}


def test_the_manifest_round_trips_every_key_but_the_saved_stamp(project_env, tmp_path):  # noqa: F811
    proj = new_test_project(tmp_path / "p1", name="study one", plugin_version="0.2.0.dev1")
    m1 = read_manifest_file(proj.root)
    assert list(m1) == ["format", "name", "id", "created", "saved", "pack", "recipe", "run", "models", "images",
                        "samples", "validated"]
    assert m1["format"] == 1 and m1["name"] == "study one" and uuid.UUID(m1["id"]).version == 4
    for key in ("created", "saved"):
        assert set(m1[key]) == STAMP_KEYS
    assert m1["created"]["plugin"] == "0.2.0.dev1" and m1["created"]["germquant"] == germquant.__version__
    assert m1["pack"] == {"id": "test_pack", "version": "0.1.0", "recipe": "basic",
                          "sha256": packs.sums_sha256(TEST_PACK)}
    assert m1["recipe"] == {"file": None, "expected_config_hash": BASIC_HASH, "overrides": {}, "roles": {},
                            "origins": {}}
    assert m1["run"] == {"xy_stride": 1, "z_range": None, "device": None, "runner": None, "interpreters": {}}
    assert m1["models"] == [] and m1["images"] == [] and m1["samples"] is None
    assert m1["validated"]["value"] is False and m1["validated"]["record"] is None
    assert "provisional" in m1["validated"]["why"]
    opened = open_project(proj.root)
    assert opened.manifest == m1 and opened.warnings == [] and opened.config_hash == BASIC_HASH
    assert opened.pack is opened.manifest["pack"] and opened.images is opened.manifest["images"]
    save_project(opened)                                   # nothing to change: only the saved stamp moves
    m2 = read_manifest_file(proj.root)
    assert _without_saved(m2) == _without_saved(m1)
    assert m2["saved"]["plugin"] is None and m2["created"] == m1["created"]
    text = (proj.root / "germquant_project.yaml").read_bytes()
    assert b"\r" not in text and text == M.dump_yaml(m2)          # LF only, exactly what the writer makes


def test_a_non_ascii_project_name_is_written_as_an_escape(project_env, tmp_path):  # noqa: F811
    name = "Caf" + chr(0xE9) + " " + GAMMA + "H2AX"
    proj = new_test_project(tmp_path / "p", name=name)
    raw = (proj.root / "germquant_project.yaml").read_bytes()
    raw.decode("ascii")
    assert b"\\xE9" in raw and b"\\u03B3" in raw
    assert open_project(proj.root).manifest["name"] == name


def test_a_failure_mid_write_leaves_the_old_file(project_env, tmp_path, monkeypatch):  # noqa: F811
    proj = new_test_project(tmp_path / "p")
    path = proj.root / "germquant_project.yaml"
    before = path.read_bytes()

    def boom(src, dst):
        raise OSError("disk full")

    monkeypatch.setattr(M, "_replace", boom)
    with pytest.raises(OSError, match="disk full"):
        save_project(proj, device="cpu")
    assert path.read_bytes() == before
    assert [p.name for p in proj.root.iterdir() if p.name.startswith(".")] == []      # no temporary file left
    calls, sleeps = [], []

    def denied(src, dst):
        calls.append(dst)
        raise PermissionError("in use")

    monkeypatch.setattr(M, "_replace", denied)
    monkeypatch.setattr(M, "_sleep", sleeps.append)
    with pytest.raises(PermissionError):
        save_project(proj, device="cpu")
    assert len(calls) == 6 and sleeps == [0.4] * 5 and path.read_bytes() == before     # five retries over 2 s
    real, seen = os.replace, []

    def flaky(src, dst):
        seen.append(dst)
        if len(seen) <= 2:
            raise PermissionError("an indexer holds the file")
        return real(src, dst)

    monkeypatch.setattr(M, "_replace", flaky)
    save_project(proj, device="cpu")
    assert read_manifest_file(proj.root)["run"]["device"] == "cpu" and len(seen) == 3
    assert not (proj.root / "germquant_project.lock").exists()                          # released every time


def _write_lock(root: Path, *, host="otherhost", pid=4242, hours_ago=0.0, command="germquant project add"):
    at = (datetime.now(timezone.utc) - timedelta(hours=hours_ago)).isoformat()
    (root / "germquant_project.lock").write_text(json.dumps({"host": host, "pid": pid, "time": at,
                                                             "command": command}), encoding="ascii")
    return at


def test_the_lock_is_an_exclusive_create_and_a_second_writer_is_refused(project_env, tmp_path):  # noqa: F811
    proj = new_test_project(tmp_path / "p")
    lockf = proj.root / "germquant_project.lock"
    with L.project_lock(proj.root, "test writer"):
        rec = json.loads(lockf.read_text(encoding="ascii"))
        assert set(rec) == {"host", "pid", "time", "command"}
        assert rec["host"] == socket.gethostname() and rec["pid"] == os.getpid() and rec["command"] == "test writer"
        with pytest.raises(ProjectLockedError, match=str(os.getpid())):
            save_project(proj, device="cpu")
        rep = L.lock_report(proj.root)
        assert rep["held"] is True and rep["stale"] is False and rep["pid"] == os.getpid()
    assert not lockf.exists() and L.lock_report(proj.root) == {"held": False, "host": None, "pid": None,
                                                                "since": None, "stale": False}
    before = (proj.root / "germquant_project.yaml").read_bytes()
    _write_lock(proj.root)
    with pytest.raises(ProjectLockedError) as e:
        add_images(proj, [EXAMPLE])
    msg = str(e.value)
    assert "otherhost" in msg and "4242" in msg and "germquant project add" in msg and "older than" not in msg
    assert e.value.holder["host"] == "otherhost" and e.value.holder["pid"] == 4242
    assert (proj.root / "germquant_project.yaml").read_bytes() == before and lockf.exists()
    lockf.write_bytes(b"")                     # a writer between its exclusive create and its write: not an old lock
    with pytest.raises(ProjectLockedError, match="cannot be read") as e:
        save_project(proj, device="cpu")
    assert "older than" not in str(e.value) and "delete" in str(e.value) and lockf.read_bytes() == b""
    lockf.unlink()
    with L.project_lock(proj.root, "test writer"):     # the lock is removed only by its maker
        other = {"host": "otherhost", "pid": 4242, "time": "2026-10-24T00:00:00+00:00", "command": "x"}
        lockf.write_text(json.dumps(other), encoding="ascii")
    assert json.loads(lockf.read_text(encoding="ascii")) == other


def test_an_old_lock_is_reported_and_only_unlock_with_its_holder_removes_it(project_env, tmp_path,  # noqa: F811
                                                                           capsys):
    proj = new_test_project(tmp_path / "p")
    lockf = proj.root / "germquant_project.lock"
    since = _write_lock(proj.root, hours_ago=13)
    assert L.lock_report(proj.root) == {"held": True, "host": "otherhost", "pid": 4242, "since": since,
                                        "stale": True}
    with pytest.raises(ProjectLockedError, match="older than 12 h") as e:
        save_project(proj, device="cpu")
    assert "germquant project unlock" in str(e.value) and "--host otherhost --pid 4242" in str(e.value)
    assert open_project(proj.root).warnings and any("older than 12 h" in w for w in open_project(proj.root).warnings)
    for host, pid in (("otherhost", 1), ("thishost", 4242)):
        with pytest.raises(ProjectError, match="otherhost"):
            break_lock(proj, host=host, pid=pid)
        assert lockf.exists()
    assert break_lock(proj, host="otherhost", pid=4242) is True and not lockf.exists()
    assert break_lock(proj, host="otherhost", pid=4242) is False
    _write_lock(proj.root, hours_ago=13)
    assert cli_projects.main(["unlock", str(proj.root), "--host", "otherhost", "--pid", "1", "--json"]) == 2
    out = last_json(capsys.readouterr().out)
    assert out["removed"] is False and out["lock"]["pid"] == 4242 and lockf.exists()
    assert out["kind"] == "ProjectError" and "not by otherhost pid 1" in out["error"]   # like every refusal
    assert cli_projects.main(["unlock", str(proj.root), "--host", "otherhost", "--pid", "4242", "--json"]) == 0
    out = last_json(capsys.readouterr().out)
    assert out == {"removed": True, "lock": {"host": "otherhost", "pid": 4242, "since": out["lock"]["since"],
                                             "command": "germquant project add"}, "stale": True}
    assert not lockf.exists()


def test_new_vendors_the_pinned_pack_and_writes_no_recipe(project_env, tmp_path):  # noqa: F811
    proj = new_test_project(tmp_path / "p")
    vend = proj.root / "packs" / "test_pack" / "0.1.0"
    assert (vend / "pack.yaml").is_file() and packs.verify_sums(vend) == []
    assert not (proj.root / "recipe.yaml").exists() and not (proj.root / "channel_map.yaml").exists()
    cfg = project_config(proj)
    assert cfg.hash == BASIC_HASH and cfg.pack_record["edited"] is False
    assert cfg.pack_record["source"] == "project" and cfg.pack_record["id"] == "test_pack"
    assert sorted(p.name for p in proj.root.iterdir()) == ["germquant_project.yaml", "packs"]
    with pytest.raises(ProjectError, match="already a germquant project"):
        new_test_project(tmp_path / "p")
    with pytest.raises(FileNotFoundError, match="nope"):
        new_project(tmp_path / "q", "nope@0.1.0")
    assert not (tmp_path / "q" / "germquant_project.yaml").exists()
    with pytest.raises(ProjectError, match="config"):
        new_test_project(tmp_path / "config")                 # load_config would anchor its recipe to the parent


def test_open_after_a_pack_version_change_warns_and_names_it(project_env, tmp_path, monkeypatch):  # noqa: F811
    proj = new_test_project(tmp_path / "p")
    saved = write_recipe(proj, edits={"spots.log.min_snr": 3.0})
    assert open_project(proj.root).warnings == []

    def to_020(d):
        d["version"] = "0.2.0"

    v2 = copy_pack(TEST_PACK, tmp_path / "env" / "test_pack" / "0.2.0", edit=to_020, rebuild=False)
    prof = v2 / "config" / "profiles" / "basic.yaml"
    prof.write_bytes(prof.read_bytes().replace(b"min_volume_um3: 2.0", b"min_volume_um3: 2.5"))
    packs.write_sums(v2)
    monkeypatch.setenv("GERMQUANT_PACKS", f"{FIXTURE_PACKS}{os.pathsep}{tmp_path / 'env'}")
    packs.refresh()
    rp = proj.root / "recipe.yaml"
    rp.write_bytes(rp.read_bytes().replace(b"pack:test_pack@0.1.0/", b"pack:test_pack@0.2.0/"))
    warnings = open_project(proj.root).warnings
    hit = [w for w in warnings if saved in w]
    assert hit and "test_pack@0.2.0" in hit[0] and "test_pack@0.1.0" in hit[0], warnings
    # the same version with other content (the project's copy replaced) is named as such
    rp.write_bytes(rp.read_bytes().replace(b"pack:test_pack@0.2.0/", b"pack:test_pack@0.1.0/"))
    vend = proj.root / "packs" / "test_pack" / "0.1.0"
    vprof = vend / "config" / "profiles" / "basic.yaml"
    vprof.write_bytes(vprof.read_bytes().replace(b"min_volume_um3: 2.0", b"min_volume_um3: 2.5"))
    packs.write_sums(vend)
    warnings = open_project(proj.root).warnings
    assert any("differs from the copy the project was made with" in w for w in warnings), warnings
    assert any(saved in w for w in warnings), warnings


def test_a_project_keeps_its_pack_copy_when_the_built_in_root_moves_on(project_env, tmp_path,  # noqa: F811
                                                                      monkeypatch):
    """P5 acceptance 13, the open half: a project made on test_pack 0.1.0 from a fake built-in root (no version
    folder) still opens with the same hash after that root moves test_pack to 0.2.0."""
    root = tmp_path / "engine"
    (root / "src" / "germquant").mkdir(parents=True)
    copy_pack(TEST_PACK, root / "packs" / "test_pack", rebuild=False)
    monkeypatch.setattr(C, "builtin_roots", lambda: [root.resolve()])
    monkeypatch.delenv("GERMQUANT_PACKS", raising=False)
    packs.refresh()
    assert packs.find_pack("test_pack").source == "builtin"
    proj = new_test_project(tmp_path / "p")
    before = project_config(proj).hash

    def to_020(d):
        d["version"] = "0.2.0"

    import shutil

    shutil.rmtree(root / "packs" / "test_pack")
    copy_pack(TEST_PACK, root / "packs" / "test_pack", edit=to_020)
    packs.refresh()
    assert [p.version for p in packs.discover().packs if p.id == "test_pack"] == ["0.2.0"]
    opened = open_project(proj.root)
    assert opened.warnings == [] and opened.config_hash == before == BASIC_HASH
    assert project_config(opened).pack_record["version"] == "0.1.0"


def test_add_images_refusals(project_env, tmp_path, monkeypatch):  # noqa: F811
    """P5 acceptance 7 (the add half): a multi-scene file, an nd2 with P > 1, a duplicate id and an unreadable
    header are refused with their reasons, and the other files are still added."""
    proj = new_test_project(tmp_path / "p")
    data = tmp_path / "data"
    good = write_ome_tiff(data / "good.ome.tif")
    two = write_ome_tiff(data / "two_scenes.ome.tif", n_images=2)
    broken = data / "broken.ome.tif"
    broken.write_bytes(b"this is not a TIFF file at all")
    plate, single = data / "plate.nd2", data / "single.nd2"
    for f in (plate, single):
        f.write_bytes(b"\0" * 32)
    fake_nd2_headers(monkeypatch, {"plate.nd2": {"sizes": {"P": 3, "Z": 5, "C": 2, "Y": 32, "X": 32}},
                                   "single.nd2": {"sizes": {"T": 4, "Z": 5, "C": 2, "Y": 32, "X": 32}}})
    dup = write_ome_tiff(data / "other" / "good.ome.tif")
    case = write_ome_tiff(data / "case" / "GOOD.ome.tif")
    notes = data / "notes.txt"
    notes.write_text("not an image", encoding="ascii")
    output = write_ome_tiff(data / "good__nuclei_labels.tif")
    res = add_images(proj, [good, two, broken, plate, single, dup, case, data / "missing.ome.tif", notes, output])
    assert res["added"] == ["good", "single"]
    reasons = {r["path"]: r["reason"] for r in res["refused"]}
    assert set(reasons) == {str(p) for p in (two, broken, plate, dup, case, data / "missing.ome.tif", notes, output)}
    assert "2 images" in reasons[str(two)] and "one image per file" in reasons[str(two)]
    assert "3 positions" in reasons[str(plate)]
    assert reasons[str(broken)].startswith("header unreadable")
    assert "duplicate image id good" in reasons[str(dup)] and "letter case" in reasons[str(case)]
    assert "not found" in reasons[str(data / "missing.ome.tif")]
    assert "unsupported" in reasons[str(notes)] and "germquant output" in reasons[str(output)]
    entries = {e["id"]: e for e in open_project(proj.root).images}
    assert set(entries) == {"good", "single"}
    g = entries["good"]
    assert Path(g["path"]) == good and Path(g["path"]).is_absolute()
    assert g["size"] == good.stat().st_size and g["mtime"] == pytest.approx(good.stat().st_mtime)
    assert g["excluded"] is False and g["reason"] is None and "timepoints" not in g
    assert entries["single"]["timepoints"] == 4
    again = add_images(proj, [good])
    assert again["added"] == [] and "duplicate image id good" in again["refused"][0]["reason"]
    dot = write_ome_tiff(data / "dot..ome.tif")                       # the id would be "dot.", a folder Windows trims
    assert "ends with a dot" in add_images(proj, [dot])["refused"][0]["reason"]
    from germquant.project.images import id_reason

    for device in ("CON", "nul", "Com1", "lpt9"):                     # no such file can even exist on Windows
        assert "Windows device name" in id_reason(device)
    assert id_reason("console") is None and id_reason("COM10") is None and id_reason("gonad 01") is None


def test_add_folder_finds_ome_tif_without_input_suffixes(project_env, tmp_path):  # noqa: F811
    """P5 acceptance 7 and D27: the test pack's recipe sets no io.input_suffixes, so a folder is discovered with
    io.formats.SUPPORTED_SUFFIXES; germquant outputs and the project's own folders are skipped."""
    proj = new_test_project(tmp_path / "p")
    assert project_config(proj).get("io.input_suffixes") is None
    data = tmp_path / "data folder"
    write_ome_tiff(data / "a.ome.tif")
    write_ome_tiff(data / "sub" / "b.ome.tif")
    write_ome_tiff(data / "res" / "a__nuclei_labels.tif")
    write_ome_tiff(data / "c__spots.tif")
    (data / "res" / "a__nuclei.csv").write_text("x\n", encoding="ascii")
    assert add_images(proj, [data]) == {"added": ["a", "b"], "refused": []}
    for d in ("results", "runs", "training", "tuning", "annotations"):
        write_ome_tiff(proj.root / d / f"in_{d}.ome.tif")
    write_ome_tiff(proj.root / "images" / "c.ome.tif")
    assert add_images(proj, [proj.root]) == {"added": ["c"], "refused": []}    # packs/ holds the example: skipped
    entries = {e["id"]: e for e in open_project(proj.root).images}
    assert entries["c"]["path"] == "images/c.ome.tif"                          # inside the project: relative
    assert Path(entries["a"]["path"]).is_absolute() and " " in entries["a"]["path"]
    write_ome_tiff(data / "g1.ome.tif")
    write_ome_tiff(data / "g2.ome.tif")
    pattern, nothing = str(data / "g*.ome.tif"), str(data / "zz*.tif")
    assert add_images(proj, [pattern, nothing]) == {"added": ["g1", "g2"], "refused": [
        {"path": nothing, "reason": "no file matches the pattern"}]}
    for d in ("results", "runs"):
        res = add_images(proj, [proj.root / d / f"in_{d}.ome.tif"])
        assert res["added"] == [] and f"the project's {d}/ folder" in res["refused"][0]["reason"]
    one = add_images(proj, str(data / "g1.ome.tif"))                  # one path, not a string of letters
    assert one["added"] == [] and [r["reason"][:22] for r in one["refused"]] == ["duplicate image id g1:"]


def test_a_recipe_with_input_suffixes_discovers_only_those(project_env, tmp_path, monkeypatch):  # noqa: F811
    def to_030(d):
        d["version"] = "0.3.0"

    v3 = copy_pack(TEST_PACK, tmp_path / "env" / "test_pack" / "0.3.0", edit=to_030, rebuild=False)
    prof = v3 / "config" / "profiles" / "basic.yaml"
    prof.write_bytes(prof.read_bytes() + b"io:\n  input_suffixes: ['.nd2']\n")
    packs.write_sums(v3)
    monkeypatch.setenv("GERMQUANT_PACKS", f"{FIXTURE_PACKS}{os.pathsep}{tmp_path / 'env'}")
    packs.refresh()
    proj = new_project(tmp_path / "p", "test_pack@0.3.0")
    data = tmp_path / "data"
    write_ome_tiff(data / "a.ome.tif")
    (data / "x.nd2").write_bytes(b"\0" * 32)
    fake_nd2_headers(monkeypatch, {"x.nd2": {"sizes": {"Z": 3, "C": 2, "Y": 16, "X": 16}}})
    assert add_images(proj, [data]) == {"added": ["x"], "refused": []}
    assert add_images(proj, [data / "a.ome.tif"])["added"] == ["a"]          # a file named explicitly is added


def test_headers_are_read_without_bioio_and_agree_with_the_run_reader(tmp_path, monkeypatch):
    yx = write_ome_tiff(tmp_path / "yx.ome.tif", names=("DAPI", GAMMA + "H2AX"), zyx=(0, 32, 40),
                        spacing=(0, 0.12, 0.11))
    h = read_header(yx)
    assert h["sizes"] == {"T": 1, "C": 2, "Z": 1, "Y": 32, "X": 40} and h["n_images"] == 1
    assert h["channel_names"] == ["DAPI", GAMMA + "H2AX"] and h["single_plane"] is True
    assert tuple(h["spacing_yx_um"]) == (0.12, 0.11)
    assert tuple(h["spacing"]) == (1.0, 1.0, 1.0) and h["spacing_ok"] is False
    assert read_tiff_header(yx)["spacing"] == (1.0, 1.0, 1.0)            # where io.ome_header gives the guess
    ex = read_header(EXAMPLE)
    assert ex["sizes"] == {"T": 1, "C": 2, "Z": 12, "Y": 64, "X": 96} and ex["single_plane"] is False
    assert tuple(ex["spacing"]) == (0.5, 0.2, 0.2) and ex["spacing_ok"] is True
    assert tuple(ex["spacing_yx_um"]) == (0.2, 0.2) and ex["channel_names"] == ["DAPI", "foci"]
    two = write_ome_tiff(tmp_path / "two.ome.tif", n_images=2)
    assert read_header(two)["n_images"] == 2
    import numpy as np
    import tifffile

    plain = tmp_path / "plain.tif"
    tifffile.imwrite(plain, np.zeros((5, 16, 16), np.uint16))
    hp = read_header(plain)
    assert hp["channel_names"] == ["Channel:0:0"] and hp["sizes"]["Z"] == 5 and hp["spacing_ok"] is False
    ij = tmp_path / "ij.tif"
    tifffile.imwrite(ij, np.zeros((4, 2, 16, 16), np.uint16), imagej=True, resolution=(5.0, 5.0),
                     metadata={"axes": "ZCYX", "spacing": 0.5, "unit": "um"})
    hij = read_header(ij)
    assert hij["channel_names"] == ["Channel:0:0", "Channel:0:1"] and tuple(hij["spacing"]) == (0.5, 0.2, 0.2)
    fake_nd2_headers(monkeypatch, {"x.nd2": {"sizes": {"P": 2, "Z": 1, "C": 3, "Y": 8, "X": 8},
                                             "channel_names": ["a", "b", "c"]}})
    hn = read_header(tmp_path / "x.nd2")
    assert hn["n_images"] == 2 and hn["single_plane"] is True and hn["reader"] == "nd2"
    assert hn["channel_names"] == ["a", "b", "c"] and tuple(hn["spacing_yx_um"]) == (0.1, 0.1)
    if importlib.util.find_spec("bioio") is None:
        pytest.skip("the readers extra is not installed here: the comparison with the run's reader needs bioio")
    import germquant.pipeline  # noqa: F401  (torch before the readers' dask can load CuPy, rule R5)
    from germquant.io import read_image_metadata

    for f in (yx, EXAMPLE, two, plain, ij):
        light, full = read_header(f), read_image_metadata(f)
        assert light["channel_names"] == full["channel_names"], f
        assert light["sizes"] == full["sizes"], f
        assert light["n_images"] == full["n_scenes"], f
        assert tuple(light["spacing"]) == tuple(full["spacing"]) and light["spacing_ok"] == full["spacing_ok"], f


def test_cli_new_add_and_rounds_print_json_and_stamp_the_plugin_version(project_env, tmp_path,  # noqa: F811
                                                                       capsys):
    root = tmp_path / "cli project"
    img = write_ome_tiff(tmp_path / "in" / "a.ome.tif")
    rc = cli_projects.main(["new", str(root), "--pack", "test_pack@0.1.0", "--images", str(img), "--name", "cli",
                            "--plugin-version", "0.2.0.dev1", "--json"])
    out = last_json(capsys.readouterr().out)
    assert rc == 0, out
    assert out["root"] == str(root) and out["name"] == "cli" and out["config_hash"] == BASIC_HASH
    assert out["pack"]["id"] == "test_pack" and out["edited"] is False and out["added"] == ["a"]
    assert out["refused"] == [] and [e["id"] for e in out["images"]] == ["a"]
    m = read_manifest_file(root)
    assert m["created"]["plugin"] == m["saved"]["plugin"] == "0.2.0.dev1"
    b = write_ome_tiff(tmp_path / "in" / "b.ome.tif")
    rc = cli_projects.main(["add", str(root), str(b), str(tmp_path / "in" / "nope.ome.tif"), "--plugin-version",
                            "0.2.0.dev2", "--json"])
    out = last_json(capsys.readouterr().out)
    assert rc == 1 and out["added"] == ["b"] and out["refused"][0]["reason"].startswith("not found")
    assert out["headers"]["b"]["channel_names"] == ["DAPI", "foci"]
    assert read_manifest_file(root)["saved"]["plugin"] == "0.2.0.dev2"
    rc = cli_projects.main(["rounds", str(root), "--json"])
    assert rc == 0 and last_json(capsys.readouterr().out) == {"rounds": [], "created": None, "path": None}
    rc = cli_projects.main(["rounds", str(root), "--new", "--plugin-version", "0.2.0.dev3", "--json"])
    out = last_json(capsys.readouterr().out)
    assert rc == 0 and out["created"] == "r01" and out["rounds"] == ["r01"] and Path(out["path"]).is_dir()
    assert read_manifest_file(root)["saved"]["plugin"] == "0.2.0.dev3"
    rc = cli_projects.main(["rounds", str(root), "--new", "bad name!", "--json"])
    assert rc == 2 and "error" in last_json(capsys.readouterr().out)
    rc = cli_projects.main(["new", str(root), "--pack", "test_pack@0.1.0", "--json"])
    assert rc == 2 and "already a germquant project" in last_json(capsys.readouterr().out)["error"]
    plain = tmp_path / "not a project"
    plain.mkdir()
    for argv in (["add", str(plain), str(img)], ["rounds", str(plain), "--new"], ["edit", str(plain), "--device", "cpu"],
                 ["use-model", str(plain), "--model", "m"], ["unlock", str(plain), "--host", "h", "--pid", "1"]):
        assert cli_projects.main([*argv, "--json"]) == 2, argv
        assert "not a germquant project" in last_json(capsys.readouterr().out)["error"], argv
    assert list(plain.iterdir()) == []                            # no lock file was made there


def test_training_rounds_are_named_and_numbered(project_env, tmp_path):  # noqa: F811
    from germquant.project import list_rounds, new_round

    proj = new_test_project(tmp_path / "p")
    assert list_rounds(proj) == []
    assert new_round(proj) == proj.root / "training" / "r01"
    assert new_round(proj, "fine_tune-2") == proj.root / "training" / "fine_tune-2"
    (proj.root / "training" / "not a round").mkdir()
    assert new_round(proj) == proj.root / "training" / "r02"
    assert list_rounds(proj) == ["fine_tune-2", "r01", "r02"]
    for bad in ("", "a/b", "x" * 41, "..", "r 1", "CON", "aux"):
        with pytest.raises(ProjectError, match="round name"):
            new_round(proj, bad)
    with pytest.raises(ProjectError, match="already exists"):
        new_round(proj, "r01")


def test_the_full_command_line_has_the_project_group(project_env, tmp_path, capsys):  # noqa: F811
    from germquant import cli

    root = tmp_path / "p"
    assert cli.main(["project", "new", str(root), "--pack", "test_pack@0.1.0", "--json"]) == 0
    assert last_json(capsys.readouterr().out)["config_hash"] == BASIC_HASH
    assert yaml.safe_load((root / "germquant_project.yaml").read_text(encoding="ascii"))["format"] == 1
    assert ROOT.is_dir()
