"""Post-stages in the run commands (docs/POST_STAGES.md, "Running post-stages"; design 5.6 and 9; PLAN.md 4.5 P5b
phase WIRE, decisions D12, D33 to D35 and D40; P5b acceptance 1 to 3, 5 and 7).

* Without post-stages nothing new is called: the default routes write no post file and never import
  germquant.post_stages (a fresh interpreter), and process_image gets exactly the keywords it got before (no
  on_stack); with post-stages on_stack joins them, and only then (D34, D40).
* `germquant batch --post-stages` calls the fixture post-stage once per image after process_image, on the stack the
  run read (the image is read once), and adds only the declared table and __post.json; `--resume` with the same
  file calls nothing and changes nothing; with changed parameters it reruns the post-stage only, while config_hash,
  the completion marker and is_done stay as they are (acceptance 2 and 3).
* collect writes batch_post_<name>_<table>.csv, the image id first, from the post tables in force only.
* `germquant run`, `germquant batch` and `germquant project run` each call the fixture once per image (a project's
  post_stages.yaml replaces its pack's), and a time series is refused at every entry point (acceptance 7).
* pack status greys a recipe whose needed post-stage is missing, from metadata only (its fresh-interpreter half is
  in tests/test_packs_isolation.py), and does not require the segmentation model of a recipe that reads its nuclei
  from labels files, warning instead (D12).
* A single plane: the post-stage sees (PLANAR_DZ, dy, dx) through a batch that keeps the stack and through a
  post-only rerun alike (germquant.planar.PLANAR_DZ, review S1 of P5b-CORE).
* A list that cannot run here is refused before any image (exit 2, or ProjectError), except a pack's own list under
  --no-pack-check, which each image then records as failed.
* The synthetic_post golden (tests/golden/synthetic_post/, its sentinel pinned in tests/golden/default_invariants.json):
  the fixture's table and the stable fields of __post.json on the synthetic gonad. Regenerate it only for a
  deliberate change of this path: ``GERMQUANT_UPDATE_GOLDEN=1 pytest tests/test_post_stages_wiring.py -k golden``.
"""
from __future__ import annotations

import json
import os
import shutil
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import planar_helpers as H
import pytest
from packs_helpers import TEST_PACK, copy_pack, pack_env  # noqa: F401  (pack_env is a fixture)
from post_helpers import (  # noqa: F401  (counting_fixture and post_pack_env are fixtures)
    POST_PACK_PARAMS,
    POST_PACK_REF,
    changed_files,
    counting_fixture,
    post_pack_env,
    post_record,
    register_counting,
    run_child,
    tree,
    write_post_file,
)
from project_helpers import fake_nd2_headers, new_test_project
from test_end_to_end import CHANNELS, SPACING, _synthetic_stack
from test_project_run import example_copy, need_readers

import germquant.post_stages as PS
from germquant import api, batch, cli, packs, pipeline, provenance, timelapse
from germquant import project as P
from germquant.config import load_config
from germquant.post_stages.testing import FIXTURE_COLUMNS, FIXTURE_NAME, FIXTURE_TABLE

ROOT = Path(__file__).resolve().parents[1]
GOLD = ROOT / "tests" / "golden" / "synthetic_post"
UPDATE = bool(os.environ.get("GERMQUANT_UPDATE_GOLDEN"))
_RTOL = os.environ.get("GERMQUANT_GOLDEN_RTOL")        # CI: another CPU's summation order (tests/README.md)
IN = "20251105_n2_nohs_HERM_001.nd2"                   # the synthetic golden's file name
OTHER = "20251021_n2_hs_MALE_2.nd2"
IN_ID, OTHER_ID = Path(IN).stem, Path(OTHER).stem
TABLE = f"{FIXTURE_NAME}_{FIXTURE_TABLE}"              # <id>__post_<TABLE>.csv; collect: batch_post_<TABLE>.csv
EXAMPLE_ID = "test_pack_synthetic"


def _study(folder: Path, *, spots: bool = True) -> Path:
    """config.yaml with the classical segmenter and no montage (and no spots when ``spots`` is false), outside the
    repository: the gonad runs in about a second. It extends config.yaml by name, found under the engine's
    built-in root, so its text and config_hash (the golden's sentinel) are the same in every checkout; an
    absolute extends path would enter the hash with the checkout's folder."""
    lines = ["extends: config.yaml", "segmentation:", "  nuclei:", "    method: classical", "render:",
             "  montage: false"]
    if not spots:
        lines += ["spots:", "  enabled: false"]
    path = Path(folder) / ("study.yaml" if spots else "study_fast.yaml")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _patch_gonad(mp, reads: list | None = None) -> None:
    """process_image reads the 3-channel synthetic gonad for any .nd2 path (tests/test_end_to_end.py); ``reads``
    collects one entry per pixel read."""
    stack = _synthetic_stack()
    meta = {"sizes": {"C": 3, "Z": 20, "Y": 80, "X": 320}, "spacing": SPACING, "spacing_ok": True,
            "channel_names": CHANNELS, "dtype": "float32", "is_2d": False}

    def read_stack(path, **kw):
        if reads is not None:
            reads.append(Path(path).name)
        return stack

    mp.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: dict(meta))
    mp.setattr(pipeline, "read_stack", read_stack)


def _spy_post(mp, posts: list) -> None:
    """Wrap germquant.post_stages.run_post_stages: each call records its image id, whether a stack was handed over,
    and the results folder before and after."""
    real = PS.run_post_stages

    def spy(out_dir, image_id, entries, **kw):
        before = tree(out_dir)
        res = real(out_dir, image_id, entries, **kw)
        posts.append({"image_id": image_id, "stack": kw.get("stack") is not None, "before": before,
                      "after": tree(out_dir), "res": res})
        return res

    mp.setattr(PS, "run_post_stages", spy)


def _image_folders(out: Path) -> dict[str, dict[str, bytes]]:
    return {iid: tree(out / iid) for iid in (IN_ID, OTHER_ID)}


def _stable(record: dict) -> dict:
    """The fields of __post.json that do not follow the run's time stamp (the golden's and the parity's): every
    stage's version, params, done, error and tables, and of its inputs the spacing, roles and image geometry."""
    keep = ("version", "params", "done", "error", "tables")
    stages = {}
    for name, st in record["stages"].items():
        inputs = st.get("inputs") or {}
        stages[name] = {**{k: st[k] for k in keep},
                        "inputs": {k: inputs[k] for k in ("spacing", "roles", "image") if k in inputs}}
    return {"format": record["format"], "image_id": record["image_id"],
            "post_stages_file_sha256": record["post_stages_file_sha256"], "stages": stages}


# ---------------------------------------------------------------------------------------------------
# one batch study through its post-stage steps (acceptance 2, 3, 5 and 7, and the golden)
# ---------------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def batch_runs(tmp_path_factory):
    """Two synthetic gonads through `germquant batch --post-stages P1` (main run and post-stage), then `--resume`
    with P1 (nothing to do), with P2 (other parameters: the post-stage only), without --post-stages, and collect."""
    work = tmp_path_factory.mktemp("post_batch")
    src = work / "in"
    src.mkdir()
    for name in (IN, OTHER):
        (src / name).write_bytes(b"")
    study = _study(work)
    p1 = write_post_file(work / "p1.yaml", {"factor": 1.5})
    p2 = write_post_file(work / "p2.yaml", {"factor": 3.0, "statistic": "max"})
    out = work / "out"
    calls, reads, posts = [], [], []
    steps = {}
    with pytest.MonkeyPatch.context() as mp:
        _patch_gonad(mp, reads)
        _spy_post(mp, posts)
        name = register_counting(calls)
        try:
            def step(key, *extra):
                calls.clear(), reads.clear(), posts.clear()
                code = cli.main(["batch", str(src), "--config", str(study), "--out", str(out), *extra])
                steps[key] = SimpleNamespace(code=code, calls=list(calls), reads=list(reads), posts=list(posts),
                                             folders=_image_folders(out))

            step("first", "--post-stages", str(p1))
            step("same", "--resume", "--post-stages", str(p1))
            step("changed", "--resume", "--post-stages", str(p2))
            step("without", "--resume")
            counts = batch.collect(out)
        finally:
            PS.unregister_post_stage(name)
    return SimpleNamespace(work=work, out=out, study=study, cfg=load_config(study), p1=p1, p2=p2, steps=steps,
                           counts=counts)


def test_fixture_writes_only_declared_table_and_record(batch_runs):
    """Acceptance 2: after process_image, the post-stage adds its declared table and __post.json to each folder and
    changes nothing that was there: the tables, label images, stage record and completion marker stay byte for
    byte as process_image left them. The stack the run read is handed over, so the image is read once."""
    first = batch_runs.steps["first"]
    assert first.code == 0
    assert [p["image_id"] for p in first.posts] == [OTHER_ID, IN_ID]          # sorted files: 20251021 first
    for post in first.posts:
        iid = post["image_id"]
        assert changed_files(post["before"], post["after"]) == {
            "added": sorted([f"{iid}__post.json", f"{iid}__post_{TABLE}.csv"]), "removed": [], "changed": []}
        assert post["stack"] is True and post["res"]["ran"] == [FIXTURE_NAME] and post["res"]["failed"] == []
        assert f"{iid}__done.json" in post["before"] and f"{iid}__nuclei_labels.tif" in post["before"]
        table = pd.read_csv(batch_runs.out / iid / f"{iid}__post_{TABLE}.csv")
        assert list(table.columns) == list(FIXTURE_COLUMNS) and len(table) > 0
    assert first.calls == [OTHER_ID, IN_ID] and first.reads == [OTHER, IN]


def test_param_change_reruns_only_post_stage(batch_runs):
    """Acceptance 3: `batch --resume` with the same post_stages.yaml calls nothing and leaves every file as it was;
    with other parameters it reruns the post-stage only (the main run is skipped), so only the post table and
    __post.json change, while config_hash, the completion marker and is_done stay as they are. A resumed batch
    without --post-stages calls nothing and keeps the post files."""
    s = batch_runs.steps
    assert s["same"].code == 0 and s["same"].calls == [] and s["same"].reads == []
    assert s["same"].folders == s["first"].folders
    changed = s["changed"]
    assert changed.code == 0 and changed.calls == [OTHER_ID, IN_ID]
    assert changed.reads == [OTHER, IN]                    # a post-only rerun reads the image again (D34)
    cfg = batch_runs.cfg
    for iid in (IN_ID, OTHER_ID):
        diff = changed_files(s["first"].folders[iid], changed.folders[iid])
        assert diff == {"added": [], "removed": [], "changed": sorted([f"{iid}__post.json", f"{iid}__post_{TABLE}.csv"])}
        assert s["first"].folders[iid][f"{iid}__done.json"] == changed.folders[iid][f"{iid}__done.json"]
        assert batch.is_done(batch_runs.out / iid, iid, cfg) == (True, "done")
        marker = provenance.read_done_marker(batch_runs.out / iid, iid)
        assert marker["config_hash"] == cfg.hash == load_config(batch_runs.study).hash
        old = json.loads(s["first"].folders[iid][f"{iid}__post.json"])["stages"][FIXTURE_NAME]
        new = post_record(batch_runs.out / iid, iid)["stages"][FIXTURE_NAME]
        assert old["params"]["factor"] == 1.5 and new["params"] == {"factor": 3.0, "min_voxels": 0, "statistic": "max",
                                                                    "fail": False}
        assert old["sub_hash"] != new["sub_hash"] and new["done"] is True
        assert old["inputs"]["main_run"] == new["inputs"]["main_run"]          # the main run was not redone
    without = s["without"]
    assert without.code == 0 and without.calls == [] and without.posts == [] and without.folders == changed.folders


def test_collect_gathers_post_tables(batch_runs, counting_fixture, tmp_path):  # noqa: F811
    """Acceptance 5, collect: batch_post_<name>_<table>.csv stacks the images' post tables, the image id first. Only
    the tables in force are gathered: a post table whose record lists it for a done stage run on the folder's
    current main run. Tables of a post-stage switched off since, or made before the main run was redone, stay on
    disk and are left out."""
    stacked = pd.read_csv(batch_runs.out / f"batch_post_{TABLE}.csv")
    assert list(stacked.columns) == ["image_id", *FIXTURE_COLUMNS]
    per_image = {iid: pd.read_csv(batch_runs.out / iid / f"{iid}__post_{TABLE}.csv") for iid in (OTHER_ID, IN_ID)}
    assert stacked["image_id"].tolist() == [iid for iid, t in per_image.items() for _ in range(len(t))]
    assert batch_runs.counts[f"post_{TABLE}"] == len(stacked) == sum(len(t) for t in per_image.values())
    assert "batch_summary" in batch_runs.counts and batch_runs.counts["nuclei"] > 0
    # the rule, on copies of the two folders: a main run redone without post-stages (its run_timestamp moves), and
    # a post-stage switched off (an empty list leaves a record with no stages)
    root = tmp_path / "copy"
    for iid in (IN_ID, OTHER_ID):
        shutil.copytree(batch_runs.out / iid, root / iid)
    summary = root / IN_ID / f"{IN_ID}__image_summary.csv"
    df = pd.read_csv(summary)
    df["run_timestamp"] = "2030-01-01T00:00:00+00:00"
    df.to_csv(summary, index=False)
    res = PS.run_post_stages(root / OTHER_ID, OTHER_ID, [], cfg=batch_runs.cfg, image_path=None)
    assert res["post"]["stages"] == {}
    counts = batch.collect(root)
    assert counts[f"post_{TABLE}"] == 0
    left = pd.read_csv(root / f"batch_post_{TABLE}.csv")
    assert list(left.columns) == ["image_id", *FIXTURE_COLUMNS] and len(left) == 0
    assert (root / IN_ID / f"{IN_ID}__post_{TABLE}.csv").is_file()      # nothing of a user's is deleted
    assert counting_fixture == []


def test_synthetic_post_golden(batch_runs):
    """The synthetic_post golden: the fixture's table on the synthetic gonad (study.yaml's recipe, parameters of
    P2) and the stable fields of its __post.json."""
    out = batch_runs.out / IN_ID
    table = pd.read_csv(out / f"{IN_ID}__post_{TABLE}.csv", low_memory=False, float_precision="round_trip")
    record = _stable(post_record(out, IN_ID))
    if UPDATE:
        GOLD.mkdir(parents=True, exist_ok=True)
        table.to_csv(GOLD / f"{FIXTURE_TABLE}.csv", index=False)
        (GOLD / "post_record.json").write_bytes((json.dumps(record, indent=1, sort_keys=True) + "\n").encode("ascii"))
        (GOLD / "config_hash.txt").write_text(batch_runs.cfg.hash + "\n", encoding="utf-8")
        pytest.skip("synthetic_post golden regenerated")
    sentinel = (GOLD / "config_hash.txt").read_text(encoding="utf-8").strip()
    assert batch_runs.cfg.hash == sentinel, "config/config.yaml or its channel map changed: the golden's hash moved"
    ref = pd.read_csv(GOLD / f"{FIXTURE_TABLE}.csv", low_memory=False, float_precision="round_trip")
    tol = {"check_exact": False, "rtol": float(_RTOL), "atol": 0.0} if _RTOL else {"check_exact": True}
    pd.testing.assert_frame_equal(table, ref, check_dtype=False, **tol)
    assert record == json.loads((GOLD / "post_record.json").read_text(encoding="ascii"))
    assert record["stages"][FIXTURE_NAME]["params"] == {"factor": 3.0, "min_voxels": 0, "statistic": "max",
                                                        "fail": False}
    assert len(ref) >= 3 and (ref["n_spots"] > 0).any() and (ref["foci_value_scaled"] == 3.0 * ref["foci_value"]).all()


# ---------------------------------------------------------------------------------------------------
# acceptance 1 and the call shape (D34, D40)
# ---------------------------------------------------------------------------------------------------

_DEFAULT_ROUTES = """
import sys
from pathlib import Path
sys.path.insert(0, {tests!r})
from test_end_to_end import CHANNELS, SPACING, _synthetic_stack
from germquant import api, cli, pipeline
from germquant.config import load_config
work, study = Path({work!r}), Path({study!r})
stack = _synthetic_stack()
meta = {{"sizes": {{"C": 3, "Z": 20, "Y": 80, "X": 320}}, "spacing": SPACING, "spacing_ok": True,
        "channel_names": CHANNELS, "dtype": "float32", "is_2d": False}}
pipeline.read_nd2_metadata = lambda *a, **k: dict(meta)
pipeline.read_stack = lambda *a, **k: stack
image = work / "in" / {image!r}
assert cli.main(["run", str(image), "--config", str(study), "--out", str(work / "run")]) == 0
assert cli.main(["batch", str(work / "in"), "--config", str(study), "--out", str(work / "batch")]) == 0
assert cli.main(["batch", str(work / "in"), "--config", str(study), "--out", str(work / "batch"), "--resume"]) == 0
assert cli.main(["collect", str(work / "batch")]) == 0
api.run(image, config=study, out=work / "api_run")
api.run_loaded(image, load_config(study), work / "api_loaded")
print("POST_FILES", sorted(p.name for p in work.rglob("*") if "__post" in p.name or "batch_post" in p.name))
print("LOADED", sorted(m for m in sys.modules if m.startswith("germquant.post_stages")))
"""


def test_no_post_stages_no_post_file_no_call(tmp_path):
    """Acceptance 1: `germquant run`, `batch`, `batch --resume`, `collect`, api.run and api.run_loaded without
    post-stages, in a fresh interpreter (the synthetic gonad through the .nd2 reader seam): no post file is written
    and germquant.post_stages is never imported, so nothing of it can be called."""
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / IN).write_bytes(b"")
    code = _DEFAULT_ROUTES.format(tests=str(ROOT / "tests"), work=str(tmp_path), study=str(_study(tmp_path, spots=False)),
                                  image=IN)
    r = run_child(["-c", code])
    assert r.returncode == 0, (r.stdout + r.stderr)[-3000:]
    assert "POST_FILES []" in r.stdout, r.stdout[-2000:]
    assert "LOADED []" in r.stdout, r.stdout[-2000:]
    assert (tmp_path / "batch" / Path(IN).stem / f"{Path(IN).stem}__done.json").is_file()


def _fake_process(seen: list):
    def fake(path, cfg, out_dir, **kw):
        seen.append(sorted(kw))
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        return {"image_id": Path(path).name.split(".")[0], "n_nuclei": 0, "qc_pass": True, "qc_flags": [],
                "out_dir": str(out_dir), "tables": {"image_summary": pd.DataFrame([{"n_nuclei": 0}])}, "stages": {}}
    return fake


def test_process_image_keywords_without_post_stages(post_pack_env, counting_fixture, monkeypatch, tmp_path):  # noqa: F811
    """D34 and D40: without post-stages every route calls process_image with exactly the keywords it passed before
    (run, run_loaded and project run: xy_stride, z_range and prov; batch: xy_stride and prov) and never calls the
    runner; with post-stages on_stack joins them, and only then, and the runner is called once per image."""
    seen_api, seen_batch, runner = [], [], []
    monkeypatch.setattr(api, "process_image", _fake_process(seen_api))
    monkeypatch.setattr(pipeline, "process_image", _fake_process(seen_batch))
    monkeypatch.setattr(PS, "run_post_stages", lambda out_dir, image_id, entries, **kw: runner.append(image_id) or
                        {"image_id": image_id, "ran": [], "skipped": [], "failed": [], "post": None})
    src = tmp_path / "in"
    src.mkdir()
    image = src / IN
    image.write_bytes(b"")
    study = _study(tmp_path, spots=False)
    cfg = load_config(study)
    proj = new_test_project(tmp_path / "proj")
    P.add_images(proj, [example_copy(proj.root / "images")])
    post = write_post_file(tmp_path / "post.yaml", {"factor": 2.0})

    def routes(with_post: bool):
        opt = ["--post-stages", str(post)] if with_post else []
        kw = {"post_stages": post} if with_post else {}
        assert cli.main(["run", str(image), "--config", str(study), "--out", str(tmp_path / "r"), *opt]) == 0
        assert cli.main(["batch", str(src), "--config", str(study), "--out", str(tmp_path / "b"), *opt]) == 0
        api.run(image, config=study, out=tmp_path / "a", **kw)
        api.run_loaded(image, cfg, tmp_path / "l", **kw)
        if with_post:
            write_post_file(proj.root / PS.POST_STAGES_FILE, {"factor": 2.0})
        assert P.run_project(proj, resume=False)["exit_code"] == 0

    routes(False)
    today_api, today_batch = ["prov", "xy_stride", "z_range"], ["prov", "xy_stride"]
    assert seen_api == [today_api] * 4 and seen_batch == [today_batch] and runner == []
    seen_api.clear(), seen_batch.clear()
    routes(True)
    assert seen_api == [sorted([*today_api, "on_stack"])] * 4 and seen_batch == [sorted([*today_batch, "on_stack"])]
    assert runner == [IN_ID, IN_ID, IN_ID, IN_ID, EXAMPLE_ID]


# ---------------------------------------------------------------------------------------------------
# acceptance 5: pack status (needs.post_stages, from metadata) and the labels_from model rule (D12)
# ---------------------------------------------------------------------------------------------------

def test_pack_status_greys_missing_post_stage(post_pack_env, monkeypatch, tmp_path):  # noqa: F811
    """A recipe that needs a post-stage this environment does not offer is unusable, with the post-stage named and
    the fix; offered (here registered in process) it is usable; installed but broken, the blocker says why. The
    check reads metadata only: the fixture's function is never called."""
    from germquant.post_stages.testing import declaration_text, write_distribution

    st = packs.pack_status(POST_PACK_REF, recipe="basic")
    assert st["usable"] is False
    (blocker,) = [b for b in st["blockers"] if FIXTURE_NAME in b["what"]]
    assert "post-stage" in blocker["what"] and "not installed" in blocker["what"]
    assert "germquant.post_stages" in blocker["fix"] and "POST_STAGES.md" in blocker["fix"]
    assert packs.pack_status("test_pack@0.1.0", recipe="basic")["usable"] is True      # needs none
    calls: list = []
    name = register_counting(calls)
    try:
        st = packs.pack_status(POST_PACK_REF, recipe="basic")
        assert st["usable"] is True, st["blockers"]
    finally:
        PS.unregister_post_stage(name)
    assert calls == []
    # installed but unusable: a distribution whose declaration file does not declare the entry point's name
    dist = tmp_path / "dist"
    write_distribution(dist, name="broken-post", version="1.0.0",
                       entry_points={FIXTURE_NAME: "broken_post:measure"},
                       files={"broken_post/__init__.py": "def measure(inputs, params):\n    return {}\n",
                              "broken_post/germquant_post_stage.yaml": declaration_text("other_name")})
    monkeypatch.syspath_prepend(str(dist))
    st = packs.pack_status(POST_PACK_REF, recipe="basic")
    assert st["usable"] is False
    (blocker,) = [b for b in st["blockers"] if FIXTURE_NAME in b["what"]]
    assert "cannot be used" in blocker["what"] and "list --all" in blocker["fix"]


def _labels_pack(folder: Path, *, labels: bool) -> str:
    """A copy of the test pack (id labels_pack) whose recipe with_model reads its nuclei from labels files when
    ``labels`` is true; its expected_config_hash follows the recipe text. Returns the pack ref."""
    dest = copy_pack(TEST_PACK, folder / "labels_pack", rebuild=False)
    profile = dest / "config" / "profiles" / "with_model.yaml"
    if labels:
        text = profile.read_text(encoding="ascii")
        profile.write_text(text.replace("    cellpose_model: tinynet@f5055f8a\n",
                                        "    cellpose_model: tinynet@f5055f8a\n"
                                        "    labels_from: '{image_dir}/{image_stem}_masks.tif'\n"), encoding="ascii")
    import yaml

    data = yaml.safe_load((dest / "pack.yaml").read_text(encoding="ascii"))
    data["id"] = "labels_pack"
    for r in data["recipes"]:
        if r["id"] == "with_model":
            r["expected_config_hash"] = load_config(profile).hash
    (dest / "pack.yaml").write_bytes(yaml.safe_dump(data, sort_keys=False, allow_unicode=False).encode("ascii"))
    packs.write_sums(dest)
    return "labels_pack@0.1.0"


def test_pack_status_labels_from_does_not_require_the_segmentation_model(pack_env, monkeypatch, tmp_path):  # noqa: F811
    """D12: a pack recipe that sets segmentation.nuclei.labels_from segments nothing, so pack status does not
    require its segmentation model: a missing one is a warning that says so, not a blocker. The same recipe without
    labels_from is blocked on the missing model, as before."""
    for labels in (False, True):
        root = tmp_path / ("labels" if labels else "model")
        ref = _labels_pack(root, labels=labels)
        monkeypatch.setenv("GERMQUANT_PACKS", str(root))
        packs.refresh()
        st = packs.pack_status(ref, recipe="with_model")
        model = [b for b in st["blockers"] if "tinynet@f5055f8a" in b["what"]]
        notes = [w for w in st["warnings"] if "tinynet@f5055f8a" in w["what"]]
        if labels:
            assert st["usable"] is True and model == [], st["blockers"]
            (note,) = notes
            assert "labels_from" in note["what"] and "not needed" in note["what"]
        else:
            assert st["usable"] is False and len(model) == 1 and notes == []
        assert [m["name"] for m in st["models"]] == ["tinynet@f5055f8a"]


# ---------------------------------------------------------------------------------------------------
# acceptance 7: series refused at every entry point; once per image on three routes
# ---------------------------------------------------------------------------------------------------

class _Reached(Exception):
    pass


@pytest.mark.parametrize("route", ["cli_t", "api_t", "project"])
def test_series_refused(route, post_pack_env, counting_fixture, monkeypatch, tmp_path, capsys):  # noqa: F811
    """D33: with post-stages switched on a time series is refused before any timepoint, in `germquant run --t`, in
    api.run and api.run_loaded given t, and in `germquant project run` (before process_series). An empty list
    switches them off, so the series is not refused for them."""
    entered = []

    def fake_series(*a, **k):
        entered.append(a[0])
        raise _Reached()

    monkeypatch.setattr(timelapse, "process_series", fake_series)
    post = write_post_file(tmp_path / "post.yaml", {"factor": 2.0})
    empty = write_post_file(tmp_path / "empty.yaml")
    image = tmp_path / "series.nd2"
    image.write_bytes(b"")
    study = _study(tmp_path, spots=False)
    if route == "cli_t":
        out = tmp_path / "out"
        base = ["run", str(image), "--config", str(study), "--out", str(out), "--t", "all"]
        assert cli.main([*base, "--post-stages", str(post)]) == 2
        err = capsys.readouterr().err
        assert f"REFUSED: {PS.SERIES_REFUSAL}" in err and not out.exists() and entered == []
        with pytest.raises(_Reached):
            cli.main([*base, "--post-stages", str(empty)])            # switched off: the series loop is reached
        # a pack recipe's own post_stages/<recipe>.yaml switches them on too, and the message says how to switch off
        assert cli.main(["run", str(image), "--pack", POST_PACK_REF, "--recipe", "basic", "--out", str(out), "--t",
                         "all"]) == 2
        err = capsys.readouterr().err
        assert f"REFUSED: {PS.SERIES_REFUSAL}" in err and "post_stages" in err and "[]" in err
        assert entered == [str(image)]
    elif route == "api_t":
        seen = []
        monkeypatch.setattr(api, "process_image", _fake_process(seen))
        with pytest.raises(ValueError, match="refused for time series"):
            api.run(image, config=study, out=tmp_path / "a", t=0, post_stages=post)
        with pytest.raises(ValueError, match="refused for time series"):
            api.run_loaded(image, load_config(study), tmp_path / "b", t=0, post_stages=[{"name": FIXTURE_NAME}])
        with pytest.raises(ValueError, match="refused for time series"):
            api.run(image, pack=POST_PACK_REF, recipe="basic", out=tmp_path / "c", t=0)
        assert seen == []
        api.run_loaded(image, load_config(study), tmp_path / "d", t=0, post_stages=[])
        assert seen == [["prov", "t", "xy_stride", "z_range"]]
    else:
        fake_nd2_headers(monkeypatch, {"series.nd2": {"sizes": {"T": 3, "C": 2, "Z": 4, "Y": 16, "X": 16}}})
        proj = new_test_project(tmp_path / "proj")
        P.add_images(proj, [image])
        write_post_file(proj.root / PS.POST_STAGES_FILE, {"factor": 2.0})
        res = P.run_project(proj)
        assert res["exit_code"] == 1 and res["counts"]["failed"] == 1 and entered == []
        assert "REFUSED: a time series of 3 timepoints" in res["message"] and PS.SERIES_REFUSAL in res["message"]
        write_post_file(proj.root / PS.POST_STAGES_FILE)                 # [] switches them off for every image
        res = P.run_project(proj)                                        # the series loop is reached (and fails)
        assert res["exit_code"] == 1 and "_Reached" in res["message"]
        assert [Path(p) for p in entered] == [image]
    assert counting_fixture == []


@pytest.mark.parametrize("route", ["run", "batch", "project"])
def test_called_once_per_image(route, batch_runs, post_pack_env, counting_fixture, monkeypatch, tmp_path):  # noqa: F811
    """Acceptance 7: `germquant run`, `germquant batch` and `germquant project run` each call the fixture
    post-stage exactly once per image. In the project, its own post_stages.yaml replaces the pack's
    post_stages/<recipe>.yaml."""
    if route == "batch":
        assert batch_runs.steps["first"].calls == [OTHER_ID, IN_ID]
        return
    if route == "run":
        reads: list = []
        _patch_gonad(monkeypatch, reads)
        image = tmp_path / IN
        image.write_bytes(b"")
        post = write_post_file(tmp_path / "post.yaml", {"factor": 1.5})
        out = tmp_path / "out"
        assert cli.main(["run", str(image), "--config", str(_study(tmp_path, spots=False)), "--out", str(out),
                         "--post-stages", str(post)]) == 0
        assert counting_fixture == [IN_ID] and reads == [IN]
        assert post_record(out, IN_ID)["stages"][FIXTURE_NAME]["done"] is True
        return
    need_readers()
    proj = P.new_project(tmp_path / "proj", POST_PACK_REF, "basic")
    P.add_images(proj, [example_copy(proj.root / "images")])
    write_post_file(proj.root / PS.POST_STAGES_FILE, {"factor": 4.0})
    assert P.run_project(proj)["exit_code"] == 0
    assert counting_fixture == [EXAMPLE_ID]
    rec = post_record(proj.root / "results" / EXAMPLE_ID, EXAMPLE_ID)["stages"][FIXTURE_NAME]
    assert rec["done"] is True and rec["params"]["factor"] == 4.0 and rec["params"]["statistic"] == "mean"
    # resumed: the main run and the post-stage are both up to date
    assert P.run_project(proj)["counts"] == {"done": 0, "failed": 0, "skipped": 1}
    assert counting_fixture == [EXAMPLE_ID]
    # without the project's file the pack's own list is in force: the resumed image runs its post-stage only
    (proj.root / PS.POST_STAGES_FILE).unlink()
    assert P.run_project(proj)["counts"] == {"done": 0, "failed": 0, "skipped": 1}
    assert counting_fixture == [EXAMPLE_ID, EXAMPLE_ID]
    rec = post_record(proj.root / "results" / EXAMPLE_ID, EXAMPLE_ID)["stages"][FIXTURE_NAME]
    assert {k: rec["params"][k] for k in POST_PACK_PARAMS} == POST_PACK_PARAMS


# ---------------------------------------------------------------------------------------------------
# a single plane (P5b-CORE review S1, with germquant.planar.PLANAR_DZ)
# ---------------------------------------------------------------------------------------------------

def test_planar_spacing_through_the_stack_and_a_post_only_rerun(counting_fixture, monkeypatch, tmp_path):  # noqa: F811
    """A single-plane run records voxel_dz_um as NaN and hands its array code (PLANAR_DZ, dy, dx). The post-stage
    sees that spacing through a batch that keeps the stack (one read) and through a post-only rerun (one more read)
    alike, so its sub-hash inputs and its volumes agree."""
    from germquant import planar

    reads: list = []
    H.patch_reader(monkeypatch, calls=reads)
    src = tmp_path / "in"
    src.mkdir()
    (src / f"{H.STEM}.nd2").write_bytes(b"")
    study = tmp_path / "study_2d.yaml"
    study.write_text(f"extends: {(ROOT / 'config' / 'profiles' / 'generic_2d.yaml').as_posix()}\n"
                     "segmentation:\n  nuclei:\n    method: classical\nspots:\n  enabled: false\n", encoding="utf-8")
    out = tmp_path / "out"
    run = ["batch", str(src), "--config", str(study), "--out", str(out)]
    seen = {}
    for key, factor, extra in (("kept", 1.0, []), ("rerun", 2.0, ["--resume"])):
        post = write_post_file(tmp_path / f"{key}.yaml", {"factor": factor})
        assert cli.main([*run, *extra, "--post-stages", str(post)]) == 0
        rec = post_record(out / H.STEM, H.STEM)["stages"][FIXTURE_NAME]
        table = pd.read_csv(out / H.STEM / f"{H.STEM}__post_{TABLE}.csv")
        seen[key] = (rec, table, len(reads))
    want = [planar.PLANAR_DZ, H.DY, H.DX]
    assert seen["kept"][0]["inputs"]["spacing"] == seen["rerun"][0]["inputs"]["spacing"] == want
    assert seen["kept"][2] == 1 and seen["rerun"][2] == 2
    assert counting_fixture == [H.STEM, H.STEM]
    for rec, table, _n in seen.values():
        assert rec["done"] is True
        assert table["volume_um3"].tolist() == pytest.approx((table["n_voxels"] * want[0] * want[1] * want[2]).tolist())
    assert seen["kept"][1]["volume_um3"].tolist() == seen["rerun"][1]["volume_um3"].tolist()
    nuclei = pd.read_csv(out / H.STEM / f"{H.STEM}__nuclei.csv")
    assert nuclei["voxel_dz_um"].isna().all()


# ---------------------------------------------------------------------------------------------------
# a list that cannot run here
# ---------------------------------------------------------------------------------------------------

def test_refused_before_any_image(post_pack_env, counting_fixture, monkeypatch, tmp_path, capsys):  # noqa: F811
    """A post-stage list that cannot be read, or holds an entry that cannot run here, is refused before any image:
    `run` and `batch` exit 2 with one REFUSED line and write nothing, and `project run` raises ProjectError before
    its run log. Under --no-pack-check a pack's own list that cannot run here is a warning instead, and the image
    records the failure in its __post.json."""
    real_process = api.process_image
    seen = []
    monkeypatch.setattr(api, "process_image", _fake_process(seen))
    monkeypatch.setattr(pipeline, "process_image", _fake_process(seen))
    src = tmp_path / "in"
    src.mkdir()
    image = src / IN
    image.write_bytes(b"")
    study = _study(tmp_path, spots=False)
    missing = tmp_path / "no_such.yaml"
    unknown = write_post_file(tmp_path / "unknown.yaml", {}, name=FIXTURE_NAME)
    unknown.write_text(unknown.read_text(encoding="ascii").replace(FIXTURE_NAME, "not_installed"), encoding="ascii")
    out = tmp_path / "out"
    assert cli.main(["run", str(image), "--config", str(study), "--out", str(out), "--post-stages", str(missing)]) == 2
    err = capsys.readouterr().err
    assert err.startswith(f"REFUSED: post-stages {missing.absolute()}: ") and not out.exists()
    assert cli.main(["batch", str(src), "--config", str(study), "--out", str(out), "--post-stages", str(unknown)]) == 2
    err = capsys.readouterr().err
    assert err.count("REFUSED:") == 1 and "not_installed" in err and "not installed" in err and not out.exists()
    proj = new_test_project(tmp_path / "proj")
    P.add_images(proj, [example_copy(proj.root / "images")])
    shutil.copyfile(unknown, proj.root / PS.POST_STAGES_FILE)
    with pytest.raises(P.ProjectError, match="not_installed"):
        P.run_project(proj)
    assert not (proj.root / "runs").exists() or list((proj.root / "runs").iterdir()) == []
    assert seen == []
    # --no-pack-check: the pack's own list names a post-stage this environment lacks (the fixture, no longer
    # registered); the run goes on, and the image records the failure
    need_readers()
    monkeypatch.setattr(api, "process_image", real_process)
    PS.unregister_post_stage(FIXTURE_NAME)
    try:
        example = example_copy(tmp_path / "ex")
        out = tmp_path / "pack_out"
        assert cli.main(["run", str(example), "--pack", POST_PACK_REF, "--recipe", "basic", "--out", str(out),
                         "--no-pack-check"]) == 0
        captured = capsys.readouterr()
        assert "warning:" in captured.out and FIXTURE_NAME in captured.out
        rec = post_record(out, EXAMPLE_ID)["stages"][FIXTURE_NAME]
        assert rec["done"] is False and "not installed" in rec["error"]
        assert provenance.read_done_marker(out, EXAMPLE_ID) is not None
    finally:
        register_counting(counting_fixture)
    assert counting_fixture == []


# ---------------------------------------------------------------------------------------------------
# reviews/final_partA.md N2: one series refusal for `run --t` and `project run`
# ---------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("case", ["plain", "labels", "planar_true", "planar_auto_plane", "planar_auto_stack",
                                  "planar_malformed", "post_stages"])
def test_one_series_refusal_for_run_t_and_project_run(case, monkeypatch, tmp_path):
    """`germquant run --t` (cli._series_refusal) and `germquant project run` (project.run.series_refusal) refuse a
    series for the same reasons with the same message: a recipe that reads labels from files, a single plane
    (io.planar true, or auto on a one-plane header), a malformed io.planar, and post-stages switched on; a plain
    recipe and io.planar auto on a stack run. Only the sample sheet differs, by design (D21): the command line
    refuses any --samples, a project only a series that has a row."""
    from germquant import planar
    from germquant.io import labels as L
    from germquant.io import sample_sheet as S
    from germquant.project import run as R

    keys = {"plain": "", "labels": "segmentation:\n  nuclei:\n    labels_from: '{image_dir}/{image_stem}_masks.tif'\n",
            "planar_true": "io:\n  planar: true\n", "planar_auto_plane": "io:\n  planar: auto\n",
            "planar_auto_stack": "io:\n  planar: auto\n", "planar_malformed": "io:\n  planar: 0\n",
            "post_stages": ""}[case]
    study = tmp_path / "study.yaml"
    study.write_text(f"extends: {(ROOT / 'config' / 'profiles' / 'generic.yaml').as_posix()}\n{keys}", encoding="utf-8")
    cfg = load_config(study)
    sizes = {"C": 2, "Z": 5, "Y": 16, "X": 16} if case == "planar_auto_stack" else {"C": 2, "Y": 16, "X": 16}
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: {"sizes": dict(sizes)})
    post = [{"name": FIXTURE_NAME}] if case == "post_stages" else None
    want = {"plain": None, "labels": L.SERIES_REFUSAL, "planar_true": planar.PLANAR_SERIES_REFUSED,
            "planar_auto_plane": planar.PLANAR_SERIES_REFUSED, "planar_auto_stack": None,
            "planar_malformed": "io.planar must be false, auto or true (or null for false); got 0",
            "post_stages": PS.SERIES_REFUSAL}[case]
    path = tmp_path / "series.nd2"
    assert cli._series_refusal(cfg, path, post_stages=post) == R.series_refusal(cfg, path, post_stages=post) == want
    # the sheet, after the recipe's own refusals and before the post-stages': any --samples on the command line, in
    # a project only a row the sheet has
    intrinsic = None if case in ("plain", "planar_auto_stack", "post_stages") else want
    assert cli._series_refusal(cfg, path, samples="samples.csv", post_stages=post) == (intrinsic or S.SERIES_REFUSAL)
    no_row = S.SampleRow(fields={}, source="samples.csv", row_sha256=None)
    row = S.SampleRow(fields={"sex": "male"}, source="samples.csv", row_sha256="0" * 64)
    assert R.series_refusal(cfg, path, row=no_row, post_stages=post) == want
    assert R.series_refusal(cfg, path, row=row, post_stages=post) == (intrinsic or S.SERIES_REFUSAL)
