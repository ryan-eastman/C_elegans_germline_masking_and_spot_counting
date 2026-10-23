"""sc-straighten on a finished synthetic run (tests/test_end_to_end.py's gonad): the fence (every pre-existing file
byte-identical, only declared new files), resume, the parameter hash, workers 1 == 2, the refusals, ids that end
in an underscore and hold a space, and the run_manifest fallback for the run geometry."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd
import pytest
import test_end_to_end as e2e

from germquant import pipeline
from germquant.config import load_config
from germquant.sc.straighten import inputs, params, runner
from germquant.sc.straighten import schema as S

IMAGE = "20251105_N2_nohs_HERM _.nd2"          # id "20251105_N2_nohs_HERM _": trailing underscore and a space


@pytest.fixture(scope="module")
def finished_run(tmp_path_factory):
    mp = pytest.MonkeyPatch()
    stack = e2e._synthetic_stack()
    meta = {"sizes": {"C": 3, "Z": 20, "Y": 80, "X": 320}, "spacing": e2e.SPACING, "spacing_ok": True,
            "channel_names": e2e.CHANNELS, "dtype": "float32", "is_2d": False}
    mp.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    mp.setattr(pipeline, "read_stack", lambda *a, **k: stack)
    cfg = load_config("config/config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    root = tmp_path_factory.mktemp("sc") / "results"
    res = pipeline.process_image(IMAGE, cfg, root / "img", prov=None)
    mp.undo()
    return root, res["image_id"], stack


def _reader(stack):
    def read(path, roles, xy_stride=1, z_range=None):
        return {r: stack.data[int(i)].copy() for r, i in roles.items()}
    return read


def _copy(finished_run, tmp_path) -> Path:
    root, _iid, _st = finished_run
    dst = tmp_path / "results"
    shutil.copytree(root, dst)
    return dst


def _digests(root: Path) -> dict:
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}


def _run(root, stack, **kw):
    kw.setdefault("params", params.resolve())
    return runner.run_results(root, reader=_reader(stack), echo=lambda *_: None, **kw)


def test_fence_new_files_only_and_declared_columns(finished_run, tmp_path):
    _root, iid, stack = finished_run
    assert iid.endswith(" _")
    root = _copy(finished_run, tmp_path)
    before = _digests(root)
    out = _run(root, stack)
    assert [o.status for o in out] == ["done"], [o.line for o in out]
    after = _digests(root)
    for k, v in before.items():
        assert after[k] == v, f"{k} changed"
    new = sorted(set(after) - set(before))
    want = sorted(f"img/{iid}{suf}" for suf in S.SUFFIXES.values())
    assert new == want
    for key, cols in S.TABLES.items():
        df = pd.read_csv(root / "img" / f"{iid}{S.SUFFIXES[key]}")
        assert list(df.columns[:len(cols)]) == cols, key
    strands = pd.read_csv(root / "img" / f"{iid}{S.SUFFIXES['strands']}")
    assert len(strands) >= 5                        # one short SC filament per germline nucleus
    assert (strands["image_id"] == iid).all() and "sc_straighten_hash" in strands
    rec = json.loads((root / "img" / f"{iid}__sc_straighten.json").read_text())
    assert rec["done"] and rec["status"] == "done" and "sc_straighten:uncalibrated" in rec["flags"]


def test_resume_param_change_and_force(finished_run, tmp_path):
    _r, iid, stack = finished_run
    root = _copy(finished_run, tmp_path)
    assert _run(root, stack)[0].status == "done"
    assert _run(root, stack)[0].status == "skipped"
    p2 = params.resolve(None, ["follow.seed_min=3.9"])
    assert p2.hash != params.resolve().hash
    o = _run(root, stack, params=p2)[0]
    assert o.status == "done"
    assert _run(root, stack, params=p2, force=True)[0].status == "done"


def test_workers_one_and_two_give_identical_tables(finished_run, tmp_path):
    _r, iid, stack = finished_run
    a = _copy(finished_run, tmp_path / "a")
    b = _copy(finished_run, tmp_path / "b")
    _run(a, stack, workers=1)
    _run(b, stack, workers=2)
    for key in ("strands", "samples", "gaps", "foci", "nuclei", "strips_flat", "strips_axial"):
        fa = a / "img" / f"{iid}{S.SUFFIXES[key]}"
        fb = b / "img" / f"{iid}{S.SUFFIXES[key]}"
        assert fa.read_bytes() == fb.read_bytes(), key


def test_refusals(finished_run, tmp_path, monkeypatch):
    _r, iid, stack = finished_run
    # missing source without a reader override: refused with the --source-map message
    root = _copy(finished_run, tmp_path / "src")
    o = runner.run_results(root, params=params.resolve(), echo=lambda *_: None)[0]
    assert o.status == "refused" and o.code == "source_missing" and "--source-map" in o.message
    # stride 2
    root = _copy(finished_run, tmp_path / "stride")
    m = root / "img" / f"{iid}__done.json"
    d = json.loads(m.read_text())
    d["xy_stride"] = 2
    m.write_text(json.dumps(d))
    o = _run(root, stack)[0]
    assert o.status == "refused" and o.code == "undersampled_xy"
    # a series folder
    root = _copy(finished_run, tmp_path / "series")
    (root / "img" / f"{iid}__series.json").write_text("{}")
    assert _run(root, stack)[0].code == "series_folder"
    # no central_element role in the run's channel map
    root = _copy(finished_run, tmp_path / "noce")
    f = root / "img" / f"{iid}__nuclei.csv"
    df = pd.read_csv(f)
    df["channel_map_json"] = json.dumps({"dna": 0, "central_element": None, "foci": 2})
    df.to_csv(f, index=False)
    assert _run(root, stack)[0].code == "no_central_element"
    # a golden tree: refused without --out, fine with it
    root = _copy(finished_run, tmp_path / "golden")
    (root / "golden_manifest.json").write_text("{}")
    assert _run(root, stack)[0].code == "golden_folder"
    out = _run(root, stack, out=tmp_path / "golden_out")[0]
    assert out.status == "done" and (tmp_path / "golden_out" / "img" / f"{iid}__sc_strands.csv").exists()
    # a sealed image: refused unless --allow-sealed, then recorded
    root = _copy(finished_run, tmp_path / "sealed")
    monkeypatch.setattr(inputs, "sealed_role", lambda *a, **k: "test")
    o = _run(root, stack)[0]
    assert o.status == "refused" and o.code == "sealed"
    o = _run(root, stack, allow_sealed=True)[0]
    assert o.status == "done" and o.record["sealed_override"] is True


def test_run_manifest_fallback_for_geometry(finished_run, tmp_path):
    _r, iid, stack = finished_run
    root = _copy(finished_run, tmp_path)
    (root / "img" / f"{iid}__done.json").unlink()
    o = _run(root, stack)[0]
    assert o.status == "done" and "geometry_from_manifest" in o.record["flags"]


def test_cli_exit_codes(finished_run, tmp_path, monkeypatch):
    from germquant import cli

    _r, iid, stack = finished_run
    root = _copy(finished_run, tmp_path)
    assert cli.main(["sc-straighten", str(root), "--set", "nonsense.key=1"]) == 2
    assert cli.main(["sc-straighten", str(root)]) == 1                       # source missing -> none finished
    monkeypatch.setattr(inputs, "read_roles", _reader(stack))
    assert cli.main(["sc-straighten", str(root), "--source-map", IMAGE + "=" + str(root / "img" / "run_manifest.json"),
                     "--workers", "1"]) == 0
