"""Shards (docs/PROJECTS.md, "Shards and SLURM"; design 7.6 and 7.8; P5 acceptance 7, the shard half).

* --shard 1/3, 2/3 and 3/3 cover the image list exactly once, also after --images; malformed shards are refused.
* Three `germquant project run --shard K/3` processes started at once on one project write three run logs,
  lose no event, never write the manifest, and project_status merges their logs.
"""
from __future__ import annotations

import json
import subprocess

import pytest
from packs_helpers import PYTHON, ROOT, TEST_PACK, pack_env, subprocess_env  # noqa: F401
from project_helpers import project_env  # noqa: F401
from test_project_run import log_events, project_with_images

import germquant.project as P
from germquant.project import run as R


def _manifest(ids, excluded=()):
    return {"images": [{"id": i, "path": f"images/{i}.ome.tif", "excluded": i in excluded} for i in ids]}


def test_shards_cover_each_image_once():
    for n_images in range(0, 8):
        ids = [f"i{k}" for k in range(n_images)]
        for n in (1, 2, 3, 5):
            parts = [[e["id"] for e in R.select_images(_manifest(ids), shard=f"{k}/{n}")] for k in range(1, n + 1)]
            assert sorted(sum(parts, [])) == sorted(ids)                          # every image once
            assert all(len(p) in (n_images // n, -(-n_images // n)) for p in parts)  # balanced
            assert parts == [[e["id"] for e in R.select_images(_manifest(ids), shard=(k, n))]
                             for k in range(1, n + 1)]                            # deterministic
    ids = ["a", "b", "c", "d", "e", "x"]
    m = _manifest(ids, excluded=("x",))
    assert [e["id"] for e in R.select_images(m)] == ["a", "b", "c", "d", "e"]
    picked = ["e", "a", "c"]
    parts = [[e["id"] for e in R.select_images(m, images=picked, shard=f"{k}/2")] for k in (1, 2)]
    assert sorted(sum(parts, [])) == sorted(picked) and parts[0] == ["a", "e"]   # in manifest order
    assert R.parse_shard(" 2 / 3 ") == (2, 3) and R.parse_shard(None) is None
    for bad in ("0/3", "4/3", "1/0", "3", "a/b", "1/3/2", "-1/3", ""):
        with pytest.raises(P.ProjectError, match="shard"):
            R.parse_shard(bad)
    with pytest.raises(P.ProjectError, match="not in the project"):
        R.select_images(m, images=["a", "nope"])
    with pytest.raises(P.ProjectError, match="excluded"):
        R.select_images(m, images=["x"])


def test_three_parallel_shards(project_env, tmp_path):  # noqa: F811
    ids = ("a", "b", "c", "d", "e")
    proj = project_with_images(tmp_path, names=ids)
    manifest = (proj.root / "germquant_project.yaml").read_bytes()
    run_ids = [f"20261024T12000{k}Z_0000000{k}" for k in (1, 2, 3)]
    env = subprocess_env(GERMQUANT_PACKS=TEST_PACK.parent)
    procs = [subprocess.Popen([PYTHON, "-m", "germquant", "project", "run", str(proj.root), "--shard", f"{k}/3",
                               "--run-id", rid, "--json"], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              env=env, cwd=ROOT)
             for k, rid in zip((1, 2, 3), run_ids)]
    outs = [p.communicate(timeout=600) for p in procs]
    for p, (out, err) in zip(procs, outs):
        assert p.returncode == 0, err.decode("utf-8", "replace")[-3000:]
    results = [json.loads(out.decode("utf-8")) for out, _ in outs]
    assert [r["run_id"] for r in results] == run_ids
    logs = [log_events(proj.root / "runs" / f"{rid}.jsonl") for rid in run_ids]
    seen = []
    for k, (ev, rid) in enumerate(zip(logs, run_ids), 1):
        assert ev[0]["event"] == "run_start" and ev[0]["shard"] == f"{k}/3" and ev[-1]["event"] == "run_end"
        mine = [i["image_id"] for i in ev[0]["images"]]
        assert mine == [e["id"] for e in R.select_images(proj.manifest, shard=f"{k}/3")]
        assert [e["image_id"] for e in ev[1:-1]] == mine and all(e["event"] == "image_done" for e in ev[1:-1])
        assert ev[-1]["counts"] == {"done": len(mine), "failed": 0, "skipped": 0} and all(e["run_id"] == rid
                                                                                          for e in ev)
        seen += mine
    assert sorted(seen) == sorted(ids)                                     # no image twice, no event lost
    assert sorted(p.name for p in (proj.root / "runs").iterdir()) == sorted(f"{r}.jsonl" for r in run_ids)
    assert (proj.root / "germquant_project.yaml").read_bytes() == manifest
    assert not (proj.root / "germquant_project.lock").exists()
    st = P.project_status(proj)
    assert [r["run_id"] for r in st["runs"]] == run_ids and all(r["finished"] for r in st["runs"])
    rows = {r["id"]: r for r in st["images"]}
    assert all(r["status"] == "done" for r in rows.values())
    by_run = {i: rid for ev, rid in zip(logs, run_ids) for i in (x["image_id"] for x in ev[0]["images"])}
    assert {i: r["last_run_id"] for i, r in rows.items()} == by_run
