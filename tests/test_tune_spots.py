"""Spot preset tuning (germquant.tune; design 6.8 and 10 step 9; plan D24 and section 9; P5 acceptance 11).

On the planted examples of tests/tune_helpers.py, finished runs of the test pack's basic recipe:

* the planted optimum is recovered in 3D, by the LoG grid (spots.log.min_snr 3 finds exactly the true foci)
  and by SpotMAX (an effect_size_min between the brightest decoy and the faintest true focus); the 2D case
  needs package P2b and skips, naming it, until merge #1 is synced into feat/t5-p5;
* run once, then filter, equals direct SpotMAX runs at each effect size, table for table;
* a LoG grid makes one detection per grid point and image (D24), SpotMAX one per effect-size group;
* the run folder holds grid.yaml, trials.jsonl, ranked.csv and best.yaml as section 9.3 describes, plus run.json;
* --max-runs refuses an oversized grid before any detection;
* a detection at a run's own settings reproduces that run's spots table (the tuner makes the spots stage's
  call), also for a run with a stride and a z range;
* the synthetic_tune golden: ranked.csv and best.yaml of the LoG grid, values to 3 decimals, with its
  config-hash sentinel (tests/golden/synthetic_tune/; GERMQUANT_UPDATE_GOLDEN=1 rewrites it);
* tune_spots on a project, `germquant tune-spots` and apply_tuning, through a stand-in of germquant.project, and
  once more through the real germquant.project (package P5-PRJ): --apply refuses a best.yaml that is not YAML
  (exit 2), exits 1 for a write the project stops part way, and prints the hash the write leaves.
"""
from __future__ import annotations

import json
import os
import shutil
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import tune_helpers as H
import yaml

from germquant import cli, planar
from germquant.tune import apply as A
from germquant.tune import grid as G
from germquant.tune import spots as S
from germquant.tune.grid import TuneRefused, TuningWarning

ROOT = Path(__file__).resolve().parents[1]
GOLD = ROOT / "tests" / "golden" / "synthetic_tune"
LOG_GRID = {"spots.log.min_snr": [1.5, 2.0, 2.5, 3.0, 4.0, 5.0]}         # the default LoG grid (D24)
ES_KEY = "spots.effect_size_min"
SPOTMAX_GRID = {ES_KEY: [float(v) for v in range(0, 61)]}                 # fine enough to land in the gap
KEEP = (0, 20, 29)                                                        # sets whose SpotMAX tables are kept
SPOT_COLS = ["spot_id", "nucleus_id", "z_um", "y_um", "x_um", "effect_size", "spot_mean_intensity", "detector"]


def _require_p2b():
    if not hasattr(planar, "decide"):
        pytest.skip("the 2D case needs package P2b (single-plane runs: germquant.planar.decide); it runs once "
                    "merge #1 is synced into feat/t5-p5")


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    """Finished runs of the two planted fields (test pack basic: LoG at min_snr 2.0), their truth points, and
    the LoG grid tuned on them once, with every backend call counted and every table kept."""
    root = tmp_path_factory.mktemp("planted")
    with pytest.MonkeyPatch.context() as mp:
        stacks, truth, decoys = {}, {}, {}
        for name in H.FIELDS:
            stacks[f"{name}.nd2"], truth[name], decoys[name] = H.planted_stack(name)
        H.patch_planted_reader(mp, stacks)
        items = []
        for name in H.FIELDS:
            rd = root / "results" / name
            H.finished_run(root / f"{name}.nd2", H.basic_recipe(), rd)
            pts = H.write_truth(root / "annotations" / f"{name}__points_ryan.csv", truth[name], rd, name)
            items.append({"image_id": name, "image": root / f"{name}.nd2", "result_dir": rd, "points": pts})
        tables: dict = {}
        with pytest.MonkeyPatch.context() as spy:
            calls = H.counting_backend(spy, "log")
            res = S._tune(items, LOG_GRID, S.OBJECTIVE, "project", backend="log", out_dir=root / "tuning" / "log",
                          max_runs=50, shortcut=False, on_tables=lambda i, k, ps, pn: tables.__setitem__((i, k), ps))
        yield SimpleNamespace(root=root, items=items, truth=truth, decoys=decoys, stacks=stacks, log=res,
                              log_calls=calls, log_tables=tables)


@pytest.fixture(scope="module")
def planted_spotmax(planted):
    """The fine effect-size grid tuned with SpotMAX on the planted runs (run once per image, then filtered)."""
    pytest.importorskip("spotmax")
    tables: dict = {}

    def keep(i, k, ps, pn):
        if i in KEEP:
            tables[(i, k)] = (ps, pn)

    with pytest.MonkeyPatch.context() as spy:
        calls = H.counting_backend(spy, "spotmax")
        res = S._tune(planted.items, SPOTMAX_GRID, S.OBJECTIVE, "project", backend="spotmax",
                      out_dir=planted.root / "tuning" / "spotmax", max_runs=50, shortcut=True, on_tables=keep)
    return SimpleNamespace(res=res, calls=calls, tables=tables)


def _gap(tables, planted) -> tuple[float, float]:
    """(largest decoy effect size, smallest true one) of the unfiltered SpotMAX run, from the planted positions."""
    decoy, true = [], []
    for k, item in enumerate(planted.items):
        ps = tables[(0, k)][0]
        tr = np.asarray(planted.truth[item["image_id"]])
        de = np.asarray(planted.decoys[item["image_id"]])
        assert len(ps) >= len(tr), "the unfiltered run missed a true focus"
        for _, row in ps.iterrows():
            p = np.array([row["z_um"], row["y_um"], row["x_um"]])
            dt = np.min(np.linalg.norm(tr - p, axis=1))
            dd = np.min(np.linalg.norm(de - p, axis=1)) if len(de) else np.inf
            assert min(dt, dd) < 0.3, f"a detection at {p} is neither a planted focus nor a decoy"
            (true if dt < dd else decoy).append(float(row["effect_size"]))
    assert decoy and max(decoy) < min(true), (sorted(decoy), sorted(true))
    return max(decoy), min(true)


# ------------------------------------------------------------------------------ acceptance 11
@pytest.mark.parametrize("dims", ["3d", "2d"])
def test_recovers_planted_optimum(dims, planted, request, tmp_path):
    if dims == "2d":
        _require_p2b()
        _planar_case(tmp_path)
        return
    res = planted.log
    assert res["best"] == {"objective": S.OBJECTIVE_TEXT, "backend": "log", "keys": {"spots.log.min_snr": H.LOG_OPTIMUM}}
    top = res["ranked"][0]
    assert (top["objective_value"], top["tie_break_value"]) == (1.0, 1.0)
    assert all(row["objective_value"] < 1.0 for row in res["ranked"][1:]), res["ranked"]     # a unique optimum
    sm = request.getfixturevalue("planted_spotmax")
    lo, hi = _gap(sm.tables, planted)
    e = sm.res["best"]["keys"][ES_KEY]
    assert lo < e <= hi, (lo, e, hi)
    assert (sm.res["ranked"][0]["objective_value"], sm.res["ranked"][0]["tie_break_value"]) == (1.0, 1.0)
    assert sm.res["best"]["backend"] == "spotmax"


def _planar_case(tmp_path):
    """The 2D planted optimum: single-plane fields run with io.planar (package P2b) and SpotMAX (the LoG
    backend is refused in planar runs)."""
    pytest.importorskip("spotmax")
    with pytest.MonkeyPatch.context() as mp:
        stacks, truth, decoys, items = {}, {}, {}, []
        for name, field, seed in (("flat_a", H.FIELD_2D_A, 5), ("flat_b", H.FIELD_2D_B, 6)):
            pix, truth[name], decoys[name] = H.planted_pixels(field, seed, shape=H.SHAPE_2D)
            stacks[f"{name}.nd2"] = H.Stack(data=pix.astype(np.float32), spacing=H.SPACING, channel_names=H.CHANNELS,
                                            path=Path(f"{name}.nd2"), spacing_ok=True)
        H.patch_planted_reader(mp, stacks)
        for name in truth:
            rd = tmp_path / "results" / name
            H.finished_run(tmp_path / f"{name}.nd2", H.basic_recipe(**{"io.planar": True, "spots.backend": "spotmax"}), rd)
            pts = H.write_truth(tmp_path / f"{name}__points_ryan.csv", truth[name], rd, name)
            items.append({"image_id": name, "image": tmp_path / f"{name}.nd2", "result_dir": rd, "points": pts})
        tables: dict = {}
        res = S._tune(items, SPOTMAX_GRID, S.OBJECTIVE, "project", backend="spotmax", out_dir=tmp_path / "tuning",
                      max_runs=50, shortcut=True,
                      on_tables=lambda i, k, ps, pn: tables.__setitem__((i, k), (ps, pn)) if i == 0 else None)
    lo, hi = _gap(tables, SimpleNamespace(items=items, truth=truth, decoys=decoys))
    assert lo < res["best"]["keys"][ES_KEY] <= hi
    assert (res["ranked"][0]["objective_value"], res["ranked"][0]["tie_break_value"]) == (1.0, 1.0)


def test_a_batch_made_single_plane_run_is_tuned_as_one(tmp_path):
    """`germquant batch` writes one run_manifest.json for its whole folder, without per-image geometry: the tuner
    reads a single-plane run from the done marker ({"planar": true}), so it detects with planar=True, and it
    refuses a backend without planar_ok (log) before any detection (review P5-TUNE S1). The run keeps its spots
    stage off, so the test needs no SpotMAX."""
    _require_p2b()
    from germquant import pipeline, provenance

    name = "flat_a"
    cfg = H.basic_recipe(**{"io.planar": True, "spots.enabled": False})
    pix, truth, _ = H.planted_pixels(H.FIELD_2D_A, 5, shape=H.SHAPE_2D)
    stack = H.Stack(data=pix.astype(np.float32), spacing=H.SPACING, channel_names=H.CHANNELS,
                    path=Path(f"{name}.nd2"), spacing_ok=True)
    with pytest.MonkeyPatch.context() as mp:
        H.patch_planted_reader(mp, {f"{name}.nd2": stack})
        out_root = tmp_path / "batch"                                     # as cli._batch: one manifest at the root
        prov = provenance.write_manifest(out_root, config_hash=cfg.hash, config=cfg.as_dict(), extra={"n_files": 1})
        rd = out_root / name
        pipeline.process_image(tmp_path / f"{name}.nd2", cfg, rd, prov=prov)
        assert not (rd / "run_manifest.json").exists() and "planar" not in prov
        pts = H.write_truth(tmp_path / f"{name}__points_ryan.csv", truth, rd, name)
        item = {"image_id": name, "image": tmp_path / f"{name}.nd2", "result_dir": rd, "points": pts}
        ctx = S._prepare(item, 0, "spotmax", [], None)
        assert ctx["planar"] is True and ctx["manifest_path"] == out_root / "run_manifest.json"
        assert ctx["spacing"][0] == planar.PLANAR_DZ
        with pytest.MonkeyPatch.context() as spy:
            calls = H.counting_backend(spy, "log")
            with pytest.raises(TuneRefused, match="log spot backend does not declare planar_ok"):
                S.tune_grid([item], {"spots.log.min_snr": [2.0]}, backend="log", out_dir=tmp_path / "t")
        assert calls == [] and not (tmp_path / "t").exists()


def test_run_once_then_filter_equals_direct_runs(planted, planted_spotmax):
    sm = planted_spotmax
    # one SpotMAX run per image, with the effect-size filter at 0 (design 6.8)
    assert len(sm.calls) == len(planted.items) and sm.res["n_runs"] == len(planted.items)
    assert all(c["effect_size_min"] == 0.0 for c in sm.calls)
    grid = {ES_KEY: [SPOTMAX_GRID[ES_KEY][i] for i in KEEP]}
    direct: dict = {}
    with pytest.MonkeyPatch.context() as spy:
        calls = H.counting_backend(spy, "spotmax")
        res = S._tune(planted.items[:1], grid, S.OBJECTIVE, "project", backend="spotmax",
                      out_dir=planted.root / "tuning" / "spotmax_direct", max_runs=50, shortcut=False,
                      on_tables=lambda i, k, ps, pn: direct.__setitem__(i, (ps, pn)))
    assert [c["effect_size_min"] for c in calls] == grid[ES_KEY]          # one direct run per effect size
    for j, i in enumerate(KEEP):
        pd.testing.assert_frame_equal(sm.tables[(i, 0)][0], direct[j][0])
        pd.testing.assert_frame_equal(sm.tables[(i, 0)][1], direct[j][1])
    assert len({len(direct[j][0]) for j in direct}) >= 3                  # the minima filtered different spots
    shortcut_trials = _trials(Path(sm.res["run_dir"]))
    record = {name: json.loads((Path(r["run_dir"]) / "run.json").read_text(encoding="ascii"))["shortcut"]
              for name, r in (("shortcut", sm.res), ("direct", res))}
    assert record == {"shortcut": "run once, then filter", "direct": "one detection per set"}   # review N1
    for j, t in enumerate(_trials(Path(res["run_dir"]))):
        s = shortcut_trials[KEEP[j]]
        assert s["params"] == t["params"] and s["shared_detection"] == {ES_KEY: 0.0} and t["shared_detection"] is None
        assert _strip(s["images"][0]) == _strip(t["images"][0])


def _trials(run_dir: Path) -> list[dict]:
    return [json.loads(line) for line in (run_dir / "trials.jsonl").read_text(encoding="ascii").splitlines()]


def _strip(image_record: dict) -> dict:
    return {k: v for k, v in image_record.items() if k != "detect_s"}


def test_log_grid_runs_one_detection_per_point(planted):
    """D24: a backend other than spotmax runs once per grid point and image, with that point's value."""
    calls = planted.log_calls
    assert len(calls) == len(LOG_GRID["spots.log.min_snr"]) * len(planted.items) == planted.log["n_runs"]
    assert [c["log"]["min_snr"] for c in calls] == [v for _ in planted.items for v in LOG_GRID["spots.log.min_snr"]]
    assert all(c["backend"] == "log" for c in calls)


# ------------------------------------------------------------------------------ the run folder (section 9.3)
def test_run_folder_matches_section_9(planted):
    run_dir = Path(planted.log["run_dir"])
    names = sorted(p.name for p in run_dir.iterdir())
    assert names == ["best.yaml", "grid.yaml", "ranked.csv", "run.json", "trials.jsonl"]
    for name in names:
        raw = (run_dir / name).read_bytes()
        assert raw.isascii() and b"\r" not in raw, name
    assert G.load_grid(run_dir / "grid.yaml") == LOG_GRID                      # the search can be repeated
    trials = _trials(run_dir)
    assert len(trials) == 6 and [t["trial"] for t in trials] == list(range(6))
    for t, value in zip(trials, LOG_GRID["spots.log.min_snr"]):
        assert t["params"] == {"spots.log.min_snr": value} and t["backend"] == "log" and t["split"] == "project"
        assert set(t["pooled"]) >= {"object_ccc", "point_f1", "point_precision", "point_recall", "n_objects",
                                    "n_points", "n_matched", "n_ours"}
        assert [im["image_id"] for im in t["images"]] == ["planted_a", "planted_b"]
        for im in t["images"]:
            assert set(im) >= {"object_ccc", "point_f1", "n_points", "n_objects", "n_matched", "n_ours"}
        assert t["pooled"]["n_points"] == sum(im["n_points"] for im in t["images"]) == 26
        assert t["pooled"]["n_objects"] == sum(im["n_objects"] for im in t["images"]) == 12
        assert t["objective"]["name"] == "object_ccc" and t["objective"]["tie_break"] == "point_f1"
    ranked = pd.read_csv(run_dir / "ranked.csv")
    assert list(ranked.columns) == ["rank", "objective", "objective_value", "tie_break", "tie_break_value",
                                    "spots.log.min_snr"]
    assert ranked["rank"].tolist() == [1, 2, 3, 4, 5, 6]
    assert ranked["spots.log.min_snr"].tolist() == [3.0, 2.5, 1.5, 2.0, 4.0, 5.0]   # ties keep the grid's order
    by_rank = {t["rank"]: t for t in trials}
    for _, row in ranked.iterrows():
        t = by_rank[int(row["rank"])]
        assert row["objective_value"] == round(t["pooled"]["object_ccc"], 3)
        assert row["tie_break_value"] == round(t["pooled"]["point_f1"], 3)
    best = yaml.safe_load((run_dir / "best.yaml").read_text(encoding="ascii"))
    assert best == planted.log["best"] and list(best) == ["objective", "backend", "keys"]
    run = json.loads((run_dir / "run.json").read_text(encoding="ascii"))
    assert run["backend"] == "log" and run["n_runs"] == 12 and run["grid"] == LOG_GRID
    assert run["match"] == {"tol_xy_um": 0.3, "tol_z_um": 0.4, "scored": "annotated objects"}
    assert run["shortcut"] == "one detection per set"                        # review P5-TUNE N1
    import hashlib

    def lf(path):                                     # digests are taken on LF-normalised bytes (review P5-TUNE N2)
        return hashlib.sha256(Path(path).read_bytes().replace(b"\r\n", b"\n")).hexdigest()

    for it, rec in zip(planted.items, run["items"]):
        assert rec["points_sha256"] == lf(it["points"])
        assert rec["done_sha256"] == lf(it["result_dir"] / f"{it['image_id']}__done.json")
        assert rec["config_hash"] == "02bb83cf74f9" and rec["n_points_outside_run"] == 0


def test_digests_ignore_line_endings(planted, tmp_path):
    """A points file re-saved with CRLF (Excel) keeps its digest in run.json, as does a done marker written on
    Windows (review P5-TUNE N2)."""
    src = Path(planted.items[0]["points"])
    crlf = tmp_path / src.name
    crlf.write_bytes(src.read_bytes().replace(b"\r\n", b"\n").replace(b"\n", b"\r\n"))
    assert crlf.read_bytes() != src.read_bytes().replace(b"\r\n", b"\n")
    assert S._lf_sha256(crlf) == S._lf_sha256(src)


def test_an_image_of_another_name_gets_a_warning(planted, tmp_path):
    """Only the image's shape is checked against its run, so an image whose file name is not the one the run read
    gets a warning; a root move keeps the name and gets none (review P5-TUNE N8)."""
    item = planted.items[0]
    warned: list = []
    S._prepare({**item, "image": tmp_path / "moved" / "PLANTED_A.nd2"}, 0, "log", warned, None)
    assert warned == []
    S._prepare({**item, "image": tmp_path / "planted_b.nd2"}, 0, "log", warned, None)
    assert len(warned) == 1, warned
    assert "its run read a file named planted_a.nd2, and the image given is planted_b.nd2" in warned[0]


def test_a_nested_batch_folder_finds_its_manifest(planted, tmp_path):
    """`germquant batch` writes <out>/<rel>/<id>/ for images in subfolders and one run_manifest.json at <out>: the
    tuner takes the nearest manifest above a run, and refuses one of another config (review P5-TUNE N3)."""
    item = planted.items[0]
    out = tmp_path / "out"
    rd = out / "sub" / "deeper" / item["image_id"]
    shutil.copytree(item["result_dir"], rd)
    (rd / "run_manifest.json").replace(out / "run_manifest.json")
    assert S._prepare({**item, "result_dir": rd}, 0, "log", [], None)["manifest_path"] == out / "run_manifest.json"
    other = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    other["config_hash"] = "000000000000"
    (out / "sub" / "run_manifest.json").write_text(json.dumps(other), encoding="utf-8")
    with pytest.raises(TuneRefused, match="the manifest is not this run's"):
        S._prepare({**item, "result_dir": rd}, 0, "log", [], None)


def test_max_runs_refuses_an_oversized_grid(planted, tmp_path):
    grid = {**LOG_GRID, "spots.log.diameter_um": [0.4, 0.45, 0.5, 0.55, 0.6]}      # 30 sets on 2 images
    assert G.n_runs(grid, 2, "log") == 60
    with pytest.MonkeyPatch.context() as spy:
        calls = H.counting_backend(spy, "log")
        with pytest.raises(TuneRefused, match=r"60 log runs .* above --max-runs 50"):
            S.tune_grid(planted.items, grid, backend="log", out_dir=tmp_path / "big")
        with pytest.raises(TuneRefused, match="above --max-runs 11"):
            S.tune_grid(planted.items, LOG_GRID, backend="log", out_dir=tmp_path / "big", max_runs=11)
    assert calls == [] and not (tmp_path / "big").exists()
    # SpotMAX counts its runs once per effect-size group: 61 minima and two sharpen values are 2 runs per image
    assert G.n_runs({**SPOTMAX_GRID, "spots.sharpen": [False, True]}, 2, "spotmax") == 4


def test_refusals_before_any_detection(planted, tmp_path):
    items = planted.items
    with pytest.MonkeyPatch.context() as spy:
        calls = H.counting_backend(spy, "log")
        with pytest.raises(TuneRefused, match="never runs on a test split"):
            S.tune_grid(items, LOG_GRID, split="test", backend="log", out_dir=tmp_path / "a")
        with pytest.raises(TuneRefused, match="objective"):
            S.tune_grid(items, LOG_GRID, objective=("object_ccc", "bias"), backend="log", out_dir=tmp_path / "a")
        with pytest.raises(TuneRefused, match="not a key the log backend reads"):
            S.tune_grid(items, {ES_KEY: [1.0]}, backend="log", out_dir=tmp_path / "a")
        with pytest.raises(TuneRefused, match="unknown spot backend 'nosuch'"):
            S.tune_grid(items, {"spots.nosuch.knob": [1.0]}, backend="nosuch", out_dir=tmp_path / "a")
        with pytest.raises(TuneRefused, match="no finished run"):
            S.tune_grid([{**items[0], "result_dir": tmp_path}], LOG_GRID, backend="log", out_dir=tmp_path / "a")
        (tmp_path / "busy").mkdir()
        (tmp_path / "busy" / "x.txt").write_text("x", encoding="ascii")
        with pytest.raises(TuneRefused, match="never overwrites"):
            S.tune_grid(items, LOG_GRID, backend="log", out_dir=tmp_path / "busy")
        mixed = pd.read_csv(items[0]["points"]).assign(projected=lambda d: [i % 2 == 0 for i in range(len(d))])
        with pytest.raises(TuneRefused, match="mix points placed on a projection"):
            S.tune_grid([{**items[0], "points": mixed}], LOG_GRID, backend="log", out_dir=tmp_path / "a")
        nowhere = pd.read_csv(items[0]["points"]).assign(x_um=500.0)
        with pytest.raises(TuneRefused, match="nothing can be scored"):
            S.tune_grid([{**items[0], "points": nowhere}], LOG_GRID, backend="log", out_dir=tmp_path / "a")
        # a malformed points table or file is a refusal too, naming the image and the row (review P5-TUNE S2)
        maybe = pd.read_csv(items[0]["points"], dtype={"projected": str}).assign(projected="maybe")
        with pytest.raises(TuneRefused, match=r"planted_a: the points of planted_a: projected must be true or false"):
            S.tune_grid([{**items[0], "points": maybe}], LOG_GRID, backend="log", out_dir=tmp_path / "a")
        bad = tmp_path / "planted_a__points_bad.csv"
        bad.write_text(Path(items[0]["points"]).read_text(encoding="ascii") + "0,1,2,3,4,5,6,7,false\n",
                       encoding="ascii")
        with pytest.raises(TuneRefused, match="point_id 0 appears twice"):
            S.tune_grid([{**items[0], "points": bad}], LOG_GRID, backend="log", out_dir=tmp_path / "a")
    assert calls == []


def test_only_annotated_objects_are_scored(planted, tmp_path):
    """A nucleus without truth points is not scored: its detections are neither false positives nor misses."""
    item = planted.items[0]
    pts = pd.read_csv(item["points"])
    dropped = int(pts["object_id"].value_counts().idxmax())             # the nucleus with the most points
    kept = pts[pts["object_id"] != dropped]
    res = S._tune([{**item, "points": kept}], {"spots.log.min_snr": [3.0]}, S.OBJECTIVE, "project", backend="log",
                  out_dir=tmp_path / "t", max_runs=1, shortcut=False)
    image = _trials(Path(res["run_dir"]))[0]["images"][0]
    assert image["n_points"] == image["n_matched"] == image["n_ours"] == len(kept)
    assert image["n_objects"] == pts["object_id"].nunique() - 1
    assert image["n_spots_all"] == len(pts)                           # the dropped nucleus's spots were found
    assert (image["object_ccc"], image["point_f1"]) == (1.0, 1.0)


def test_equal_counts_rank_the_exact_set_first(planted, tmp_path):
    """When every annotated object holds the same number of points, Lin's CCC is 0 / 0: the set whose counts equal
    the truth scores 1.0 and every other set 0.0, so the exact set ranks first, not last (review P5-TUNE B1)."""
    assert S._ccc([1, 1, 1], [1, 1, 1]) == 1.0 and S._ccc([1, 1, 1], [1, 2, 1]) == 0.0
    assert S._ccc([2, 2], [3, 3]) == 0.0 and S._ccc([1, 3, 2], [1, 3, 2]) == 1.0
    assert np.isnan(S._ccc([2], [2])) and np.isnan(S._ccc([], []))       # one object defines no concordance
    item = planted.items[0]
    pts = pd.read_csv(item["points"])
    per_object = pts["object_id"].value_counts()
    single = per_object[per_object == 1].index                             # the nuclei with one true focus each
    sub = pts[pts["object_id"].isin(single)].reset_index(drop=True).assign(point_id=lambda d: np.arange(len(d)))
    assert len(sub) == 3
    # 1.5 also finds the faint decoys of those nuclei, 2.5 exactly their foci, 4.0 none of them
    res = S._tune([{**item, "points": sub}], {"spots.log.min_snr": [1.5, 2.5, 4.0]}, S.OBJECTIVE, "project",
                  backend="log", out_dir=tmp_path / "t", max_runs=3, shortcut=False)
    assert res["best"]["keys"] == {"spots.log.min_snr": 2.5}
    rows = [(r["spots.log.min_snr"], r["objective_value"], r["tie_break_value"]) for r in res["ranked"]]
    assert rows == [(2.5, 1.0, 1.0), (1.5, 0.0, 0.75), (4.0, 0.0, None)], rows
    pooled = _trials(Path(res["run_dir"]))[1]["pooled"]
    assert pooled["object_ccc"] == 1.0 and pooled["truth_count"] == pooled["ours_count"] == 3


def test_projected_points_are_matched_in_2d(planted, tmp_path):
    """Points placed on a projection carry no depth: they are matched in 2D (design 6.7), 3D points in 3D."""
    item = planted.items[0]
    shifted = pd.read_csv(item["points"]).assign(z_um=lambda d: d["z_um"] + 0.6)     # beyond the 0.4 um in z
    f1 = {}
    for flag in (True, False):
        res = S._tune([{**item, "points": shifted.assign(projected=flag)}], {"spots.log.min_snr": [3.0]},
                      S.OBJECTIVE, "project", backend="log", out_dir=tmp_path / str(flag), max_runs=1, shortcut=False)
        image = _trials(Path(res["run_dir"]))[0]["images"][0]
        f1[flag] = (image["point_f1"], image["n_matched"], image["n_points"])
    assert f1[True] == (1.0, len(shifted), len(shifted))
    assert f1[False][1] == 0


@pytest.mark.parametrize("backend", ["log", "spotmax"])
def test_detection_reproduces_the_run(backend, planted, tmp_path):
    """The tuner's detection at a run's own settings equals the spots table that run wrote (the tuner makes the
    call the spots stage makes, on the run's labels and geometry)."""
    if backend == "log":
        pairs = [(planted.log_tables[(1, k)], it) for k, it in enumerate(planted.items)]      # min_snr 2.0
    else:
        pytest.importorskip("spotmax")
        stack = H.planted_stack("tiny")[0]
        item = {"image_id": "tiny", "image": tmp_path / "tiny.nd2", "result_dir": tmp_path / "results" / "tiny"}
        got: dict = {}
        with pytest.MonkeyPatch.context() as mp:
            H.patch_planted_reader(mp, {"tiny.nd2": stack})
            H.finished_run(item["image"], H.basic_recipe(**{"spots.backend": "spotmax"}), item["result_dir"])
            item["points"] = H.write_truth(tmp_path / "tiny__points_ryan.csv", H.planted_stack("tiny")[1],
                                           item["result_dir"], "tiny")
            S._tune([item], {ES_KEY: [3.0]}, S.OBJECTIVE, "project", backend="spotmax", out_dir=tmp_path / "t",
                    max_runs=1, shortcut=False, on_tables=lambda i, k, ps, pn: got.__setitem__(k, ps))
        pairs = [(got[0], item)]
    for ps, item in pairs:
        run = pd.read_csv(item["result_dir"] / f"{item['image_id']}__spots.csv", float_precision="round_trip")
        assert len(run) > 0
        pd.testing.assert_frame_equal(ps[SPOT_COLS].reset_index(drop=True), run[SPOT_COLS], check_dtype=False,
                                      check_column_type=False)


def test_runs_with_a_stride_and_a_z_range(tmp_path):
    """A run made with xy_stride 2 on planes 2 to 10: the points saved on it lie on the file's grid, and the
    tuner reads the image, the labels and the points at that run's geometry."""
    name = "tiny"
    stack, truth, _ = H.planted_stack(name)
    with pytest.MonkeyPatch.context() as mp:
        H.patch_planted_reader(mp, {f"{name}.nd2": stack})
        rd = tmp_path / "results" / name
        H.finished_run(tmp_path / f"{name}.nd2", H.basic_recipe(), rd, xy_stride=2, z_range=(2, 11))
        pts = H.write_truth(tmp_path / f"{name}__points_ryan.csv", truth, rd, name, xy_stride=2, z_range=(2, 11))
        df = pd.read_csv(pts)
        want = np.rint(np.asarray(truth) / np.asarray(H.SPACING)).astype(int)
        assert (df[["zi", "yi", "xi"]].to_numpy() == want).all()
        got: dict = {}
        item = {"image_id": name, "image": tmp_path / f"{name}.nd2", "result_dir": rd, "points": pts}
        res = S._tune([item], {"spots.log.min_snr": [2.0]}, S.OBJECTIVE, "project", backend="log",
                      out_dir=tmp_path / "t", max_runs=1, shortcut=False,
                      on_tables=lambda i, k, ps, pn: got.__setitem__(k, ps))
    run = pd.read_csv(rd / f"{name}__spots.csv", float_precision="round_trip")
    assert len(run) > 0
    pd.testing.assert_frame_equal(got[0][SPOT_COLS].reset_index(drop=True), run[SPOT_COLS], check_dtype=False,
                                  check_column_type=False)
    rec = json.loads((tmp_path / "t" / "run.json").read_text(encoding="ascii"))["items"][0]
    assert rec["xy_stride"] == 2 and rec["z_range"] == [2, 11] and rec["n_points_outside_run"] == 0
    assert rec["spacing_um"] == [0.5, 0.4, 0.4]
    assert rec["n_points_scored"] == len(truth) and not res["warnings"]


# ------------------------------------------------------------------------------ the synthetic_tune golden
def test_synthetic_tune_golden(planted):
    run_dir = Path(planted.log["run_dir"])
    sentinel = H.basic_recipe().hash
    if os.environ.get("GERMQUANT_UPDATE_GOLDEN"):
        GOLD.mkdir(parents=True, exist_ok=True)
        for name in ("ranked.csv", "best.yaml"):
            (GOLD / name).write_bytes((run_dir / name).read_bytes())
        (GOLD / "config_hash.txt").write_text(sentinel + "\n", encoding="ascii", newline="\n")
        pytest.skip("synthetic_tune golden regenerated")
    assert GOLD.is_dir(), "no synthetic_tune golden checked in; run with GERMQUANT_UPDATE_GOLDEN=1 once"
    assert (GOLD / "config_hash.txt").read_text(encoding="ascii").strip() == sentinel == "02bb83cf74f9"
    got = pd.read_csv(run_dir / "ranked.csv")
    ref = pd.read_csv(GOLD / "ranked.csv")
    assert list(got.columns) == list(ref.columns) and len(got) == len(ref)
    for col in ref.columns:
        if pd.api.types.is_numeric_dtype(ref[col]):
            np.testing.assert_array_equal(got[col].astype(float).round(3), ref[col].astype(float).round(3), err_msg=col)
        else:
            assert got[col].astype(str).tolist() == ref[col].astype(str).tolist(), col
    best = yaml.safe_load((run_dir / "best.yaml").read_text(encoding="ascii"))
    ref_best = yaml.safe_load((GOLD / "best.yaml").read_text(encoding="ascii"))
    assert {k: best[k] for k in ("objective", "backend")} == {k: ref_best[k] for k in ("objective", "backend")}
    assert {k: round(float(v), 3) for k, v in best["keys"].items()} == \
        {k: round(float(v), 3) for k, v in ref_best["keys"].items()}


# ------------------------------------------------------------------------------ grids
def test_grid_rules(tmp_path):
    assert G.default_grid(H.basic_recipe()) == LOG_GRID
    assert G.default_grid(H.basic_recipe(**{"spots.backend": "spotmax"})) == \
        {ES_KEY: [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]}
    with pytest.raises(G.GridError, match="no default tuning grid"):
        G.default_grid(H.basic_recipe(**{"spots.backend": "net3d"}))
    assert G.check_grid({ES_KEY: [0, 2, 3.5]}, "spotmax") == {ES_KEY: [0.0, 2.0, 3.5]}   # the key's type
    assert G.check_grid({"spots.max_spot_candidates": [100, 200]}, "spotmax") == {"spots.max_spot_candidates": [100, 200]}
    assert G.check_grid({"spots.custom.knob": ["a", "b"]}) == {"spots.custom.knob": ["a", "b"]}   # an entry point's
    bad = [({}, "non-empty mapping"), ({"segmentation.nuclei.diameter_um": [1.0]}, "spot keys only"),
           ({"spots.backend": ["log"]}, "never tuned"), ({"spots.instances": [1]}, "never tuned"),
           ({"spots.log.min_snr": []}, "non-empty list"), ({"spots.log.min_snr": 2.0}, "non-empty list"),
           ({"spots.log.min_snr": [2, 2.0]}, "appears twice"), ({"spots.log.min_snr": [True]}, "finite numbers"),
           ({"spots.log.min_snr": [float("nan")]}, "finite numbers"), ({"spots.merge_z_columns": [1]}, "true or false"),
           ({"spots.max_spot_candidates": [1.5]}, "whole numbers"), ({"spots.log.min_snr": [[1, 2]]}, "scalar"),
           ({"spots.lg.min_snr": [2.0]}, "not a key the log backend reads")]
    for grid, message in bad:
        with pytest.raises(G.GridError, match=message):
            G.check_grid(grid, "log")
    f = tmp_path / "grid.yaml"
    f.write_bytes(b"\xef\xbb\xbfspots.log.min_snr: [2, 3]\nspots.log.diameter_um: [0.4, 0.5]\n")
    assert G.load_grid(f) == {"spots.log.min_snr": [2, 3], "spots.log.diameter_um": [0.4, 0.5]}
    assert [p["spots.log.min_snr"] for p in G.parameter_sets(G.load_grid(f))] == [2, 2, 3, 3]   # first key slowest
    f.write_text("- not a mapping\n", encoding="ascii")
    with pytest.raises(G.GridError, match="a grid file is a mapping"):
        G.load_grid(f)
    f.write_bytes(b"spots.log.min_snr: [2, 3]  # 0.5 \xb5m, saved as cp1252\n")          # not UTF-8 text
    with pytest.raises(G.GridError, match="not a YAML file"):
        G.load_grid(f)
    a = {ES_KEY: 1.0, "spots.sharpen": True}
    assert G.group_key(a) == G.group_key({ES_KEY: 4.0, "spots.sharpen": True}) != G.group_key({**a, "spots.sharpen": False})


# ------------------------------------------------------------------------------ projects, the command, apply
def _project(planted, root: Path, cfg=None, images=None) -> H.FakeProject:
    """A project folder holding copies of the planted runs and their truth points."""
    for it in planted.items:
        shutil.copytree(it["result_dir"], root / "results" / it["image_id"])
        (root / "annotations").mkdir(parents=True, exist_ok=True)
        shutil.copy2(it["points"], root / "annotations" / Path(it["points"]).name)
    images = images or [{"id": "planted_a", "path": "planted_a.nd2"},
                        {"id": "planted_b", "path": str(root / "elsewhere" / "planted_b.nd2")}]
    return H.FakeProject(root, images, cfg or H.basic_recipe())


def test_tune_spots_on_a_project(planted, tmp_path, monkeypatch):
    root = tmp_path / "proj"
    edited = H.basic_recipe(**{"spots.log.min_snr": 2.5})            # the recipe moved on since the runs
    project = _project(planted, root, cfg=edited, images=[
        {"id": "planted_a", "path": "planted_a.nd2"},
        {"id": "planted_b", "path": str(root / "elsewhere" / "planted_b.nd2")},
        {"id": "gone_c", "path": "gone_c.nd2", "excluded": True}])
    (root / "annotations" / "gone_c__points_ryan.csv").write_text("x", encoding="ascii")
    (root / "annotations" / "stranger__points_ryan.csv").write_text("x", encoding="ascii")
    (root / "annotations" / "notes.csv").write_text("x", encoding="ascii")
    (root / "grid.yaml").write_text("spots.log.min_snr: [2.0, 3.0]\n", encoding="ascii")
    monkeypatch.setattr(S, "_project_api", lambda: H.fake_project_api(project)[:2])
    with pytest.MonkeyPatch.context() as mp:
        H.patch_planted_reader(mp, planted.stacks)
        res = S.tune_spots(root, grid="grid.yaml")
    run_dir = Path(res["run_dir"])
    assert run_dir.parent == root / "tuning" and S.RUN_ID_RE.fullmatch(run_dir.name)
    assert res["best"]["keys"] == {"spots.log.min_snr": 3.0} and [r["spots.log.min_snr"] for r in res["ranked"]] == [3.0, 2.0]
    text = "\n".join(res["warnings"])
    assert "gone_c is excluded" in text and "stranger is not an image of this project" in text
    assert "notes.csv: not a truth-point file name" in text
    assert "was made with config 02bb83cf74f9, the recipe is now" in text and "with other spots: settings" in text
    run = json.loads((run_dir / "run.json").read_text(encoding="ascii"))
    assert [r["result_dir"] for r in run["items"]] == ["results/planted_a", "results/planted_b"]
    assert run["items"][0]["points"] == "annotations/planted_a__points_ryan.csv" and run["warnings"] == res["warnings"]
    # two annotators for one image, and an empty folder, are refused before any detection
    shutil.copy2(root / "annotations" / "planted_a__points_ryan.csv", root / "annotations" / "planted_a__points_kate.csv")
    with pytest.raises(TuneRefused, match="several annotators"):
        S.tune_spots(project)
    (root / "empty").mkdir()
    with pytest.raises(TuneRefused, match="belongs to an image of this project"):
        S.tune_spots(project, points="empty")


def test_tune_spots_command_and_apply(planted, tmp_path, monkeypatch, capsys):
    root = tmp_path / "proj"
    project = _project(planted, root, images=[{"id": "planted_a", "path": "planted_a.nd2"}])
    edits: list = []
    api = H.fake_project_api(project, edits)
    monkeypatch.setattr(S, "_project_api", lambda: api[:2])
    monkeypatch.setattr(A, "_project_api", lambda: api)
    (root / "g.yaml").write_text("spots.log.min_snr: [2.0, 3.0]\n", encoding="ascii")
    with pytest.MonkeyPatch.context() as mp:
        H.patch_planted_reader(mp, planted.stacks)
        orig_detect = S._detect

        def chatty(*a, **k):                     # a backend that prints, as SpotMAX does (review P5-TUNE N6)
            print("Quantifying spots...")
            return orig_detect(*a, **k)

        mp.setattr(S, "_detect", chatty)
        assert cli.main(["tune-spots", "--project", str(root), "--grid", "g.yaml", "--json"]) == 0
        captured = capsys.readouterr()
        out = json.loads(captured.out)                                # stdout is one JSON document
        assert "Quantifying spots..." in captured.err and "Quantifying" not in captured.out
        assert out["best"]["keys"] == {"spots.log.min_snr": 3.0} and len(out["ranked"]) == 2
        rel = Path(out["run_dir"]).relative_to(root).as_posix()
        assert cli.main(["tune-spots", "--project", str(root), "--grid", "g.yaml", "--max-runs", "1"]) == 2
        assert "above --max-runs 1" in capsys.readouterr().err
    assert cli.main(["tune-spots", "--project", str(root), "--apply", rel, "--json"]) == 0
    assert json.loads(capsys.readouterr().out) == {"config_hash": "c0ffee000000", "run_dir": rel, "warnings": []}
    assert edits == [{"edits": {"spots.log.min_snr": 3.0}, "origin": "tuning", "plugin_version": None}]
    # the plugin's Apply records its version in the saved stamp (D39, D42; review P5-TUNE S3)
    assert cli.main(["tune-spots", "--project", str(root), "--apply", rel, "--plugin-version", "0.2.0.dev1"]) == 0
    assert edits[-1] == {"edits": {"spots.log.min_snr": 3.0}, "origin": "tuning", "plugin_version": "0.2.0.dev1"}
    assert capsys.readouterr().out == f"applied {rel}: config_hash c0ffee000000\n"    # no hash at open: no change said
    assert cli.main(["tune-spots", "--project", str(root), "--plugin-version", "0.2.0.dev1"]) == 2
    assert "--plugin-version goes with --apply" in capsys.readouterr().err
    # --apply takes none of the search's options, however they are spelled (review P5-TUNE N5)
    for extra in (["--grid", "g.yaml"], ["--points", "annotations/"], ["--points", "annotations"],
                  ["--max-runs", "50"], ["--out", "tuning/x"]):
        assert cli.main(["tune-spots", "--project", str(root), "--apply", rel, *extra]) == 2, extra
        assert f"it takes no {extra[0]}" in capsys.readouterr().err
    assert len(edits) == 2
    # a recipe edited since the tuning run's images were run: applied, with a warning (review P5-TUNE N9)
    project.cfg = H.basic_recipe(**{"spots.log.diameter_um": 0.5})
    assert cli.main(["tune-spots", "--project", str(root), "--apply", rel, "--json"]) == 0
    captured = capsys.readouterr()
    notes = json.loads(captured.out)["warnings"]
    assert len(notes) == 1 and "the recipe is now" in notes[0] and f"WARNING: {notes[0]}" in captured.err
    project.cfg = H.basic_recipe()
    assert cli.main(["tune-spots", "--project", str(tmp_path / "nowhere"), "--json"]) == 2
    assert "cannot open the project" in capsys.readouterr().err

    def broken_config(p):
        raise ValueError("the pinned pack version is not installed")

    monkeypatch.setattr(S, "_project_api", lambda: (api[0], broken_config))
    assert cli.main(["tune-spots", "--project", str(root)]) == 2              # refused, before any detection
    assert "the project's recipe cannot be loaded" in capsys.readouterr().err
    # a malformed truth-point file is a refusal naming its row (exit 2); a failing detection exits 1 with one
    # ERROR line (review P5-TUNE S2)
    monkeypatch.setattr(S, "_project_api", lambda: api[:2])
    pfile = root / "annotations" / "planted_a__points_ryan.csv"
    good = pfile.read_text(encoding="ascii")
    pfile.write_text(good.replace(",false\n", ",maybe\n", 1), encoding="ascii")
    assert cli.main(["tune-spots", "--project", str(root), "--grid", "g.yaml"]) == 2
    assert "projected must be true or false (row 1: 'maybe')" in capsys.readouterr().err
    pfile.write_text(good, encoding="ascii")
    (root / "cp1252.yaml").write_bytes(b"spots.log.min_snr: [2.0, 3.0]  # 0.5 \xb5m\n")    # a grid that is not UTF-8
    assert cli.main(["tune-spots", "--project", str(root), "--grid", "cp1252.yaml"]) == 2
    assert "cp1252.yaml: not a YAML file" in capsys.readouterr().err

    def crash(*a, **k):
        raise RuntimeError("the backend crashed")

    monkeypatch.setattr(S, "_detect", crash)
    with pytest.MonkeyPatch.context() as mp:
        H.patch_planted_reader(mp, planted.stacks)
        assert cli.main(["tune-spots", "--project", str(root), "--grid", "g.yaml"]) == 1
    assert "ERROR: RuntimeError: the backend crashed" in capsys.readouterr().err


def test_apply_tuning_checks_the_run(planted, tmp_path, monkeypatch):
    root = tmp_path / "proj"
    project = _project(planted, root)
    edits: list = []
    monkeypatch.setattr(A, "_project_api", lambda: H.fake_project_api(project, edits))
    run = root / "tuning" / "r1"
    shutil.copytree(planted.log["run_dir"], run)
    with warnings.catch_warnings():
        warnings.simplefilter("error", TuningWarning)                    # the recipe is the runs' config: no warning
        assert A.apply_tuning(project, "tuning/r1", plugin_version="0.2.0.dev1") == "c0ffee000000"
    assert edits == [{"edits": {"spots.log.min_snr": 3.0}, "origin": "tuning", "plugin_version": "0.2.0.dev1"}]
    # a recipe edited after the tuning run's images were run: written, with a warning (review P5-TUNE N9)
    project.cfg = H.basic_recipe(**{"spots.log.diameter_um": 0.5})
    with pytest.warns(TuningWarning, match="r1: its images were run with config 02bb83cf74f9, and the recipe is now"):
        A.apply_tuning(project, run)
    assert len(edits) == 2
    project.cfg = H.basic_recipe()
    best = yaml.safe_load((run / "best.yaml").read_text(encoding="ascii"))
    (run / "best.yaml").write_text(yaml.safe_dump({**best, "keys": {**best["keys"], "spots.log.diameter_um": 0.4}}),
                                   encoding="ascii")
    with pytest.raises(TuneRefused, match="which the run's grid does not search"):
        A.apply_tuning(project, run)
    (run / "best.yaml").write_text(yaml.safe_dump(best), encoding="ascii")
    project.cfg = H.basic_recipe(**{"spots.backend": "spotmax"})
    with pytest.raises(TuneRefused, match="tuned the log backend; the recipe now uses spotmax"):
        A.apply_tuning(project, run)
    (run / "best.yaml").unlink()
    with pytest.raises(TuneRefused, match="no best.yaml"):
        A.apply_tuning(project, run)
    assert len(edits) == 2


def test_apply_refuses_a_best_yaml_that_is_not_yaml(planted, tmp_path, monkeypatch, capsys):
    """A best.yaml that is not YAML (a hand edit gone wrong) or not UTF-8 text is a refusal naming the file, in
    Python (TuneRefused) and from the command (exit 2), and nothing is written (review P5-TUNE-sync N4)."""
    root = tmp_path / "proj"
    project = _project(planted, root)
    edits: list = []
    monkeypatch.setattr(A, "_project_api", lambda: H.fake_project_api(project, edits))
    run = root / "tuning" / "r1"
    shutil.copytree(planted.log["run_dir"], run)
    for text, why in ((b"keys: {spots.log.min_snr: [3.0\n", "but got '<stream end>' at line 2, column 1"),
                      (b"keys:\n\tspots.log.min_snr: 3.0\n", "cannot start any token at line 2, column 1"),
                      (b"objective: object_ccc\nbackend: log\nkeys: {spots.log.min_snr: 3.0}  # 0.5 \xb5m\n",
                       "can't decode byte 0xb5")):
        (run / "best.yaml").write_bytes(text)
        with pytest.raises(TuneRefused, match="best.yaml: not a YAML file"):
            A.apply_tuning(project, run)
        assert cli.main(["tune-spots", "--project", str(root), "--apply", "tuning/r1"]) == 2, text
        err = capsys.readouterr().err
        assert "best.yaml: not a YAML file" in err and why in err and err.count("\n") == 1    # one ERROR line
    assert edits == []


# ------------------------------------------------------------------------------ with germquant.project (P5-PRJ)
def test_apply_tuning_writes_recipe_yaml_with_origin_tuning(planted, tmp_path, monkeypatch):
    """apply_tuning through the real project API: recipe.yaml extends the pack's recipe with the best keys,
    the config hash moves and the project reads the tuned value."""
    from packs_helpers import use_fixture_packs

    from germquant import project as P

    use_fixture_packs(monkeypatch)
    project = P.new_project(tmp_path / "proj", "test_pack")
    before = P.project_config(project).hash
    assert before == "02bb83cf74f9" and not (project.root / "recipe.yaml").exists()
    shutil.copytree(planted.log["run_dir"], Path(project.root) / "tuning" / "r1")
    h = A.apply_tuning(project, "tuning/r1", plugin_version="0.2.0.dev1")
    text = (Path(project.root) / "recipe.yaml").read_text(encoding="ascii")
    assert "extends: pack:test_pack@" in text and "min_snr" in text
    reopened = P.open_project(project.root)
    cfg = P.project_config(reopened)
    assert h == cfg.hash != before and cfg.get("spots.log.min_snr") == 3.0
    # the manifest records the key as a tuning edit (D19) and the plugin's version in the saved stamp (D42)
    assert reopened.manifest["recipe"]["origins"] == {"spots.log.min_snr": "tuning"}
    assert reopened.manifest["saved"]["plugin"] == "0.2.0.dev1"


def test_tune_spots_on_a_real_project(planted, tmp_path, monkeypatch):
    """tune_spots on a germquant.project project whose image is an OME-TIFF with a finished run under
    results/<id>/ (the layout run_project writes)."""
    from packs_helpers import use_fixture_packs

    from germquant import project as P

    pytest.importorskip("bioio")                      # the OME-TIFF run reads through bioio (the readers extra)
    use_fixture_packs(monkeypatch)
    pix = H.planted_pixels(H.FIELD_A, 1)[0]
    image = H.write_ome_tiff(tmp_path / "data" / "planted_a.ome.tif", pix)
    project = P.new_project(tmp_path / "proj", "test_pack", images=[image])
    iid = project.images[0]["id"]
    rd = Path(project.root) / "results" / iid
    H.finished_run(image, P.project_config(project), rd)
    H.write_truth(Path(project.root) / "annotations" / f"{iid}__points_ryan.csv", planted.truth["planted_a"], rd, iid)
    res = S.tune_spots(project, grid={"spots.log.min_snr": [2.0, 3.0]})
    assert res["best"]["keys"] == {"spots.log.min_snr": 3.0} and Path(res["run_dir"]).parent.name == "tuning"


def test_apply_that_stops_part_way_exits_1(planted, tmp_path, monkeypatch, capsys):
    """A write the real project stops part way (recipe.yaml replaced, then the manifest's os.replace refused
    after its retries, as a virus scanner can on Windows) is a failed write: exit 1, with the ERROR line naming
    what was written (9.4 #7; review P5-TUNE-sync S1). The same command run again finishes the write."""
    from packs_helpers import use_fixture_packs

    from germquant import project as P
    from germquant.project import manifest as M

    use_fixture_packs(monkeypatch)
    project = P.new_project(tmp_path / "proj", "test_pack")
    root = Path(project.root)
    shutil.copytree(planted.log["run_dir"], root / "tuning" / "r1")
    apply = ["tune-spots", "--project", str(root), "--apply", "tuning/r1"]
    real_replace = M._replace

    def replace(src, dst):
        if Path(str(dst)).name == M.MANIFEST_NAME:
            raise PermissionError(13, "held by a scanner", str(dst))
        return real_replace(src, dst)

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(M, "_replace", replace)
        mp.setattr(M, "_sleep", lambda s: None)
        assert cli.main(apply) == 1
    assert "the write stopped part way (wrote recipe.yaml)" in capsys.readouterr().err
    assert (root / "recipe.yaml").is_file() and P.open_project(root).manifest["recipe"]["origins"] == {}
    assert cli.main(apply) == 0                                        # run again: the write finishes
    assert P.open_project(root).manifest["recipe"]["origins"] == {"spots.log.min_snr": "tuning"}


def test_apply_reports_the_hash_it_leaves(planted, tmp_path, monkeypatch, capsys):
    """Without --json, --apply prints the project's hash after the write first on its line (the napari plugin's
    runner reads it there) and says what changed. A best value equal to the pack recipe's own is no edit
    (docs/PROJECTS.md): on a fresh project nothing changes, and after a tuning edit it drops that edit, so
    recipe.yaml goes and the project hashes as the pack recipe again (review P5-TUNE-sync N1)."""
    import re

    from packs_helpers import use_fixture_packs

    from germquant import project as P

    use_fixture_packs(monkeypatch)
    project = P.new_project(tmp_path / "proj", "test_pack")
    root = Path(project.root)
    pack_hash = project.config_hash                              # test_pack basic sets spots.log.min_snr 2.0
    for name, value in (("r3", 3.0), ("r2", 2.0)):
        run = root / "tuning" / name
        shutil.copytree(planted.log["run_dir"], run)
        best = yaml.safe_load((run / "best.yaml").read_text(encoding="ascii"))
        (run / "best.yaml").write_text(yaml.safe_dump({**best, "keys": {"spots.log.min_snr": value}}),
                                       encoding="ascii")

    def apply(name):
        assert cli.main(["tune-spots", "--project", str(root), "--apply", f"tuning/{name}"]) == 0
        captured = capsys.readouterr()
        now = P.open_project(root).config_hash
        # the plugin's reading (napari_germquant/_runner.py, tune_apply): the first such match in the merged output
        assert re.search(r"config_hash\W+([0-9a-f]{12})\b", captured.err + captured.out).group(1) == now
        return captured.out, now

    out, h = apply("r2")                                         # the pack's own value on a fresh project
    assert h == pack_hash and "unchanged: the recipe already holds this setting" in out
    assert not (root / "recipe.yaml").exists()
    out, tuned = apply("r3")
    assert tuned != pack_hash and f"(was {pack_hash})" in out and "stale until they are run again" in out
    out, h = apply("r2")                                         # back to the pack's value: the last edit goes
    assert h == pack_hash and f"(was {tuned})" in out and not (root / "recipe.yaml").exists()
    assert P.open_project(root).manifest["recipe"]["origins"] == {}
