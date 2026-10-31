"""Headless parity of the project route (docs/PROJECTS.md, "Portability"; design 9 routes 3 to 5, 7.7 and 3.4;
PLAN.md 4.4 P5-RUN and D25; P5 acceptance 2, 4 and 12).

* `germquant project run`, `germquant run --pack` and api.run(pack=) on the test pack's example write the same
  file list, cell-identical tables (run_timestamp aside), the same config_hash and stage hashes, and done markers
  equal in every key but finished; the project route also equals the test pack's committed recipe golden (G8).
* A project made under root A and copied to root B (one image by relative path, one by an absolute path holding
  a space) keeps its config_hash and stage hashes and reads done; an edited celegans_germline project keeps its
  stage hashes under both roots and with the engine config root moved.
* Where a snakemake executable is available (D25: a throwaway environment, never the .venv), a dry run of the
  Snakefile's project mode lists one job per image under its image_stem id; it skips everywhere else.
* With the fixture post-stage switched on (P5b acceptance 6): a project of the post-stage pack (its recipe basic
  switches the fixture on in post_stages/basic.yaml), `germquant run --pack` and api.run(pack=) write the same file
  list, post table and __post.json included, the same tables, markers, stage hashes and label pixels, and post
  records equal in every field but those that follow the run (the time stamp, and the digests of files that hold
  it); config_hash stays the test pack's basic hash.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tifffile
from packs_helpers import CELEGANS_PACK, ROOT, TEST_PACK, pack_env, subprocess_env  # noqa: F401
from project_helpers import BASIC_HASH, TEST_PACK_REF, new_test_project, project_env  # noqa: F401
from post_helpers import POST_PACK_PARAMS, POST_PACK_REF, counting_fixture, post_pack_env  # noqa: F401
from test_project_run import example_copy, need_readers

import germquant.project as P
from germquant import api, cli, packs, packs_golden, stages
from germquant import config as C

IID = "test_pack_synthetic"
FIXED = {tool: "PINNED" for tools in stages._STAGE_TOOLS.values() for tool in tools}
ROLES = {"dna": 0, "central_element": 1, "foci": 2, "granule": None, "lamin": None}


def _tables(out: Path) -> dict[str, pd.DataFrame]:
    # a post table carries no run_timestamp column (only the columns its post-stage declares)
    return {f.name: pd.read_csv(f, low_memory=False, float_precision="round_trip").drop(columns=["run_timestamp"],
                                                                                      errors="ignore")
            for f in sorted(out.glob(f"{IID}__*.csv"))}


def _marker(out: Path) -> dict:
    m = json.loads((out / f"{IID}__done.json").read_text(encoding="utf-8"))
    m.pop("finished")
    return m


def test_project_route_equals_pack_routes(project_env, tmp_path):  # noqa: F811
    need_readers()
    proj = new_test_project(tmp_path / "proj")
    image = example_copy(proj.root / "images")
    P.add_images(proj, [image])
    res = P.run_project(proj)
    assert res["exit_code"] == 0
    routes = {"project": proj.root / "results" / IID, "cli": tmp_path / "cli", "api": tmp_path / "api"}
    assert cli.main(["run", str(image), "--pack", TEST_PACK_REF, "--recipe", "basic", "--out",
                     str(routes["cli"])]) == 0
    api.run(image, pack=TEST_PACK_REF, recipe="basic", out=routes["api"])
    names = {k: sorted(f.name for f in v.iterdir()) for k, v in routes.items()}
    assert names["project"] == names["cli"] == names["api"] and len(names["project"]) > 20
    tables = {k: _tables(v) for k, v in routes.items()}
    assert len(tables["project"]) >= 12
    for other in ("cli", "api"):
        assert tables["project"].keys() == tables[other].keys()
        for name in tables["project"]:
            pd.testing.assert_frame_equal(tables["project"][name], tables[other][name], check_exact=True, obj=name)
    markers = {k: _marker(v) for k, v in routes.items()}
    assert markers["project"] == markers["cli"] == markers["api"]
    assert markers["project"]["config_hash"] == BASIC_HASH
    manifests = {k: json.loads((v / "run_manifest.json").read_text("utf-8")) for k, v in routes.items()}
    assert {m["config_hash"] for m in manifests.values()} == {BASIC_HASH}
    assert manifests["project"]["stage_hashes"] == manifests["cli"]["stage_hashes"] == manifests["api"]["stage_hashes"]
    assert manifests["project"]["pack"]["edited"] is False and "edited" not in manifests["cli"]["pack"]
    # G8: the project route equals the test pack's committed recipe golden
    gold = packs_golden.expected_dir(packs.find_pack(TEST_PACK_REF), "basic")
    for t in packs_golden.TABLES:
        got = packs_golden._read(routes["project"] / f"{IID}__{t}.csv")
        ref = pd.read_csv(gold / f"{t}.csv", low_memory=False, float_precision="round_trip")
        assert packs_golden.compare_table(got, ref, t) == []
    assert (gold / "config_hash.txt").read_text(encoding="ascii").strip() == BASIC_HASH


def _post_fields(out: Path) -> dict:
    """__post.json without what follows the run's time stamp: each stage's sub_hash and finished, and of its inputs
    main_run and the file digests of the tables (every table holds run_timestamp) and of the label image (its TIFF
    tags record the run; its pixels are compared on their own)."""
    rec = json.loads((out / f"{IID}__post.json").read_text(encoding="ascii"))
    for st in rec["stages"].values():
        st.pop("sub_hash"), st.pop("finished")
        st["inputs"].pop("main_run"), st["inputs"].pop("tables"), st["inputs"].pop("labels")
    return rec


def test_project_parity_with_fixture_post_stage(post_pack_env, counting_fixture, tmp_path):  # noqa: F811
    """P5b acceptance 6: the parity of the first test with the fixture post-stage switched on, here by the pack
    itself (post_stages/basic.yaml of the post-stage pack): the project route, `germquant run --pack` and
    api.run(pack=) each call it once and write the same files."""
    need_readers()
    proj = P.new_project(tmp_path / "proj", POST_PACK_REF, "basic")
    image = example_copy(proj.root / "images")
    P.add_images(proj, [image])
    assert P.run_project(proj)["exit_code"] == 0
    routes = {"project": proj.root / "results" / IID, "cli": tmp_path / "cli", "api": tmp_path / "api"}
    assert cli.main(["run", str(image), "--pack", POST_PACK_REF, "--recipe", "basic", "--out", str(routes["cli"])]) == 0
    api.run(image, pack=POST_PACK_REF, recipe="basic", out=routes["api"])
    assert counting_fixture == [IID, IID, IID]
    names = {k: sorted(f.name for f in v.iterdir()) for k, v in routes.items()}
    assert names["project"] == names["cli"] == names["api"]
    assert {f"{IID}__post.json", f"{IID}__post_fixture_nuclei_per_nucleus.csv"} <= set(names["project"])
    tables = {k: _tables(v) for k, v in routes.items()}
    assert f"{IID}__post_fixture_nuclei_per_nucleus.csv" in tables["project"]
    for other in ("cli", "api"):
        assert tables["project"].keys() == tables[other].keys()
        for name in tables["project"]:
            pd.testing.assert_frame_equal(tables["project"][name], tables[other][name], check_exact=True, obj=name)
    markers = {k: _marker(v) for k, v in routes.items()}
    assert markers["project"] == markers["cli"] == markers["api"]
    assert markers["project"]["config_hash"] == BASIC_HASH                 # the post-stage is never hashed
    manifests = {k: json.loads((v / "run_manifest.json").read_text("utf-8")) for k, v in routes.items()}
    assert manifests["project"]["stage_hashes"] == manifests["cli"]["stage_hashes"] == manifests["api"]["stage_hashes"]
    posts = {k: _post_fields(v) for k, v in routes.items()}
    assert posts["project"] == posts["cli"] == posts["api"]
    labels = {k: tifffile.imread(v / f"{IID}__nuclei_labels.tif") for k, v in routes.items()}
    assert np.array_equal(labels["project"], labels["cli"]) and np.array_equal(labels["project"], labels["api"])
    stage = posts["project"]["stages"]["fixture_nuclei"]
    assert stage["done"] is True and {k: stage["params"][k] for k in POST_PACK_PARAMS} == POST_PACK_PARAMS


def test_project_moves_between_roots(project_env, tmp_path):  # noqa: F811
    """P5 acceptance 4: root A to root B, one image by relative path, one by an absolute path with a space."""
    need_readers()
    a = tmp_path / "root A" / "study"
    proj = new_test_project(a)
    inside = example_copy(a / "images", "inside")
    outside = example_copy(tmp_path / "data dir", "outside image")
    P.add_images(proj, [inside, outside])
    paths = {e["id"]: e["path"] for e in proj.images}
    assert paths == {"inside": "images/inside.ome.tif", "outside image": str(outside)}
    hashes_a = stages.stage_hashes(P.project_config(proj), ROLES, FIXED, model_sha="fixed-model-digest")
    assert P.run_project(proj)["counts"]["done"] == 2
    markers_a = {iid: _marker_of(a / "results" / iid, iid) for iid in paths}
    b = tmp_path / "root B" / "deeper" / "study"
    shutil.copytree(a, b)
    shutil.rmtree(a)
    moved = P.open_project(b)
    assert moved.config_hash == BASIC_HASH and moved.warnings == []
    assert stages.stage_hashes(P.project_config(moved), ROLES, FIXED, model_sha="fixed-model-digest") == hashes_a
    st = P.project_status(b)
    assert {r["id"]: r["status"] for r in st["images"]} == {"inside": "done", "outside image": "done"}
    assert {r["id"]: r["path"] for r in st["images"]} == {"inside": str(b / "images" / "inside.ome.tif"),
                                                          "outside image": str(outside)}
    assert P.run_project(moved)["counts"] == {"done": 0, "failed": 0, "skipped": 2}
    # run again under root B: the same config_hash and stage hashes, so the same marker but finished. The image
    # by relative path is the one whose path moved (the skip above already showed both images' hashes hold)
    assert P.run_project(moved, images=["inside"], resume=False)["counts"]["done"] == 1
    assert _marker_of(b / "results" / "inside", "inside") == markers_a["inside"]


def _marker_of(out: Path, iid: str) -> dict:
    m = json.loads((out / f"{iid}__done.json").read_text(encoding="utf-8"))
    m.pop("finished")
    return m


def _engine_root_copy(dest: Path) -> Path:
    shutil.copytree(ROOT / "config", dest / "config")
    shutil.copytree(CELEGANS_PACK, dest / "packs" / "celegans_germline")
    (dest / "src" / "germquant").mkdir(parents=True)
    return dest


def test_edited_celegans_project_hashes_root_independent(pack_env, monkeypatch, tmp_path):  # noqa: F811
    """P5 acceptance 4, second half: one parameter edit through recipe.yaml on celegans_germline; the stage hashes
    agree under two project roots and with the engine config root moved (design 3.4's portable payload)."""
    seen = []
    for engine, project in (("engine_a", "root A/study"), ("deeper/engine b", "root B/x/study")):
        root = _engine_root_copy(tmp_path / engine)
        monkeypatch.setattr(C, "builtin_roots", lambda root=root: [root.resolve()])
        packs.refresh()
        proj = P.new_project(tmp_path / project, "celegans_germline")
        h = P.write_recipe(proj, edits={"spots.effect_size_min": 3.5})
        cfg = P.project_config(proj)
        assert cfg.hash == h and cfg.get("spots.effect_size_min") == 3.5
        assert str(root.resolve()) in cfg.get("io.channel_map")            # the absolute path differs per root ...
        seen.append((h, stages.stage_hashes(cfg, ROLES, FIXED, model_sha="fixed-model-digest")))
    (ha, sa), (hb, sb) = seen
    assert ha == hb                                                        # ... the hashes do not
    assert sa == sb, {k: (sa[k], sb[k]) for k in sa if sa[k] != sb[k]}


def _snakemake() -> str:
    exe = os.environ.get("GERMQUANT_SNAKEMAKE") or shutil.which("snakemake")
    if not exe:
        pytest.skip("snakemake is not installed here: it is checked in a throwaway environment, never in the .venv "
                    "(PLAN.md D25); set GERMQUANT_SNAKEMAKE to its snakemake executable to run this test")
    return exe


def test_snakemake_dry_run_project_mode(project_env, tmp_path):  # noqa: F811
    """P5 acceptance 12: `snakemake -n` lists the .ome.tif jobs of a project under their image_stem ids. It runs
    from the repository folder, as the docs say, and leaves nothing there."""
    exe = _snakemake()
    proj = new_test_project(tmp_path / "proj")
    P.add_images(proj, [example_copy(proj.root / "images", "g1"), example_copy(proj.root / "images", "g 2")])
    # Snakemake keeps its cache under the user's local application data, which Windows finds under USERPROFILE
    # (the fixture's home): without that folder the cache path comes out relative to the working folder
    (project_env.home / "AppData" / "Local").mkdir(parents=True, exist_ok=True)
    left = {name: (ROOT / name).exists() for name in (".snakemake", "snakemake")}
    try:
        r = subprocess.run([exe, "-n", "-p", "-s", str(ROOT / "workflow" / "Snakefile"), "--cores", "1", "--config",
                            f"project={proj.root}"], capture_output=True, text=True, cwd=ROOT,
                           env=subprocess_env(GERMQUANT_PACKS=TEST_PACK.parent), timeout=600)
        stray = [name for name, was in left.items() if name != ".snakemake" and not was and (ROOT / name).exists()]
    finally:
        for name, was in left.items():           # only what this run made: .snakemake is Snakemake's own record
            if not was:
                shutil.rmtree(ROOT / name, ignore_errors=True)
    assert r.returncode == 0, (r.stdout + r.stderr)[-3000:]
    assert stray == [], f"snakemake left {stray} in {ROOT}"
    out = r.stdout + r.stderr
    for iid in ("g1", "g 2"):
        assert f'project run "{proj.root}" --images "{iid}" --via snakemake' in out, out[-3000:]
    assert f'germquant project collect "{proj.root}"' in out
    assert ".ome" not in "".join(line for line in out.splitlines() if "--images" in line)
