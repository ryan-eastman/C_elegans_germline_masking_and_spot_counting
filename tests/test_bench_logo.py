"""Leave-one-gonad-out (germquant.benchmark.logo, scripts/bench_logo.py; docs/BENCHMARK_PROTOCOL.md
section 6.2): the fold rule (designated val gonads, a same-condition substitute by tag order for a held-out
val gonad, never a test gonad or a gonad without usable spot truth), the D3 stage generators, the staged
replay of the D3 search on a synthetic score index (including the parameter sets a fold still has to run),
and the driver end to end on the planted-foci gonads of test_bench_detector with the log backend
(held-out scoring, the trial records' notes, the pooled report and the paired comparison). The
pre-registered SpotMAX grid of protocol section 8 item 2: its 60 sets (values, order, names, the deployed
recipe among them, only spots: keys); the primary objective against bench_detector.py's own pooling (a
gonad without truth nuclei or with NaN counts inside R included) and against bench_tune.py's own run
folders; the truth-mask counts in the run index; the grid rule (the rounding, both tie-breaks, NaN and
missing values, a score without the truth-mask counts held as missing, no test gonad); each fold's own
selection gonads; and `bench_logo.py grid-configs` and `select --search prereg-grid` end to end (a file
bench_tune.py accepts, exit code 3 while pairs are missing, a selection that `heldout --selection` reads,
the D3 replay unchanged, the real split's folds with `to_run` merged over both val folds, the refusals).
Synthetic only.
"""
from __future__ import annotations

import importlib.util
import json
import math
import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import yaml
from test_bench_detector import planted_bench  # noqa: F401  (the fixture)

from germquant.benchmark import detector as bd
from germquant.benchmark import logo
from germquant.benchmark import stats
from germquant.benchmark import tuning as bt
from germquant.benchmark.manifest import Entry, load_manifest
from germquant.benchmark.report import write_outputs

ROOT = Path(__file__).resolve().parents[1]


def _load_script(name: str):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), ROOT / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------------- folds
def test_folds_follow_the_protocol_rule():
    def e(tag, role, cond, usable=True):
        return Entry(tag=tag, results_dir=Path("."), split=role,
                     extra={"role": role, "condition": cond, "spot_truth_usable": usable})

    entries = [e("g_noHS_t1", "train", "noHS"), e("g_noHS_t2", "train", "noHS"), e("g_HS_t1", "train", "HS"),
               e("g_HS_rg", "train", "HS", usable=False), e("g_noHS_v", "val", "noHS"), e("g_HS_v", "val", "HS"),
               e("g_noHS_test", "test", "noHS")]
    folds = logo.logo_folds(entries, require_cache=False)
    by = {f["held_out"]: f for f in folds}
    assert sorted(by) == ["g_HS_t1", "g_HS_v", "g_noHS_t1", "g_noHS_t2", "g_noHS_v"]     # no test, no unusable
    # a held-out train gonad: the designated val gonads select, the other train gonads tune
    assert by["g_noHS_t1"]["val"] == ["g_HS_v", "g_noHS_v"] and by["g_noHS_t1"]["train"] == ["g_HS_t1", "g_noHS_t2"]
    assert by["g_noHS_t1"]["substitute"] is None
    # a held-out val gonad: the other val gonad plus the first train gonad of its condition by tag order
    assert by["g_noHS_v"]["val"] == ["g_HS_v", "g_noHS_t1"] and by["g_noHS_v"]["substitute"] == "g_noHS_t1"
    assert by["g_noHS_v"]["train"] == ["g_HS_t1", "g_noHS_t2"]
    assert by["g_HS_v"]["val"] == ["g_noHS_v", "g_HS_t1"] and by["g_HS_v"]["train"] == ["g_noHS_t1", "g_noHS_t2"]
    for f in folds:
        assert f["held_out"] not in f["train"] + f["val"] and len(f["train"]) + len(f["val"]) == 4
    # the real split: the six usable train and val gonads, the two region-growing males and the test set left out
    real = [Entry(tag=t, results_dir=Path("."), split="", extra={"condition": "HS" if "_HS_" in t else "noHS",
                                                               "spot_truth_usable": t not in ("n2rad51_HS_male_000",
                                                                                              "n2rad51_HS_male_002")})
            for t in ("n2rad51_noHS_herm_000", "n2rad51_noHS_herm_001", "n2rad51_noHS_herm_002", "n2rad51_noHS_male_000",
                      "n2rad51_noHS_male_001", "n2rad51_noHS_male_002", "n2rad51_HS_herm_000", "n2rad51_HS_herm_001",
                      "n2rad51_HS_herm_002", "n2rad51_HS_male_000", "n2rad51_HS_male_001", "n2rad51_HS_male_002")]
    rf = {f["held_out"]: f for f in logo.logo_folds(real, require_cache=False)}
    assert sorted(rf) == ["n2rad51_HS_herm_000", "n2rad51_HS_herm_001", "n2rad51_noHS_herm_000",
                          "n2rad51_noHS_herm_001", "n2rad51_noHS_male_000", "n2rad51_noHS_male_001"]
    assert rf["n2rad51_HS_herm_001"]["substitute"] == "n2rad51_HS_herm_000"
    assert rf["n2rad51_noHS_herm_001"]["substitute"] == "n2rad51_noHS_herm_000"


# ----------------------------------------------------------------------- stage generators
def test_stage_generators():
    a = logo.log_stage_a()
    assert len(a) == 30 and a[0]["name"] == "logA_d0.35_s1.5" and a[0]["params"] == {"log": {"diameter_um": 0.35,
                                                                                             "min_snr": 1.5}}
    b = logo.log_stage_b(0.5, 2.5)
    assert len(b) == 17 and "logB_d0.5_s2.5_k0.5_g3_m1" not in [x["name"] for x in b]
    k75 = [x for x in b if x["name"] == "logB_d0.5_s2.5_k0.75_g4_m0"][0]
    assert k75["params"] == {"log": {"diameter_um": 0.5, "min_snr": 2.5, "min_separation_um": 0.375,
                                     "pre_gate_sigmas": 4.0}, "merge_z_columns": False}
    assert logo.log_point(k75["params"]) == {"d": 0.5, "s": 2.5, "k": 0.75, "g": 4.0, "m": False}
    c = logo.log_stage_c(logo.log_point(k75["params"]))
    assert len(c) == 8 and c[0]["name"] == "logC_d0.46_s2.25_k0.75_g4_m0"
    assert c[0]["params"]["log"]["min_separation_um"] == 0.345
    # at the edge of the stage-A ranges the neighbourhood is cut
    assert len(logo.log_stage_c({"d": 0.35, "s": 1.5, "k": 0.5, "g": 3.0, "m": True})) == 3
    s2 = logo.spotmax_stage_s2("threshold_otsu", 4.0)
    assert len(s2) == 12 and s2[0]["name"] == "smS2_otsu_r0.25_g0.08_e3.5"
    assert len(logo.spotmax_stage_s1()) == 12


# --------------------------------------------------------------- the staged replay, synthetic
_EXPERT = np.tile(np.arange(10, dtype=float), 5)


def _score(dist: float) -> logo.GonadScore:
    """A gonad score whose CCC and F1 fall with ``dist``."""
    return logo.GonadScore(expert=_EXPERT, ours=_EXPERT + 0.5 * dist, n_matched=max(1, int(100 - 10 * dist)),
                           n_ours=100, n_expert=100)


def _dist(params: dict, target: dict) -> float:
    p = logo.log_point(params)
    return (abs(p["d"] - target["d"]) / 0.04 + abs(p["s"] - target["s"]) / 0.25 + (p["k"] != target["k"])
            + (p["g"] != target["g"]) + (p["m"] != target["m"]))


def _add(idx: logo.RunIndex, configs, tags, target):
    for c in configs:
        sha = bd.params_sha256(c["params"])
        idx.params[sha] = c["params"]
        idx.names[sha] = c["name"]
        for t in tags:
            idx.scores.setdefault(sha, {})[t] = _score(_dist(c["params"], target))


def test_replay_log_search_runs_the_stages_and_names_what_is_missing():
    train, val = ["t1", "t2"], ["v1"]
    target = {"d": 0.47, "s": 2.6, "k": 0.75, "g": 2.0, "m": True}
    idx = logo.RunIndex()
    _add(idx, logo.log_stage_a(), train, target)
    r = logo.replay_log_search(idx, train, val)
    assert not r["complete"] and r["stage"] == "B" and r["stage_a_best"] == "logA_d0.5_s2.5"
    assert {m["name"] for m in r["missing"]} == {x["name"] for x in logo.log_stage_b(0.5, 2.5)}
    _add(idx, logo.log_stage_b(0.5, 2.5), train, target)
    r = logo.replay_log_search(idx, train, val)
    assert r["stage"] == "C" and r["stage_b_best"] == "logB_d0.5_s2.5_k0.75_g2_m1"
    centre = logo.log_point(logo.log_config("B", 0.5, 2.5, 0.75, 2.0, True)["params"])
    _add(idx, logo.log_stage_c(centre), train, target)
    r = logo.replay_log_search(idx, train, val)
    assert r["stage"] == "val" and len(r["finalists"]) == 5 and r["finalists"][0] == "logC_d0.46_s2.5_k0.75_g2_m1"
    assert r["n_sets"] == 55 and all(m["tags"] == ["v1"] for m in r["missing"])
    # val prefers another finalist: the selection follows val, among the finalists only
    fin = [x for x in logo.log_stage_a() + logo.log_stage_b(0.5, 2.5) + logo.log_stage_c(centre)
           if x["name"] in r["finalists"]]
    vt = logo.log_point([x for x in fin if x["name"] == r["finalists"][3]][0]["params"])
    _add(idx, fin, val, vt)
    r = logo.replay_log_search(idx, train, val)
    assert r["complete"] and r["selected"] == r["finalists"][3] and r["missing"] == []
    assert r["selected_params"] == [x for x in fin if x["name"] == r["selected"]][0]["params"]


def test_index_runs_refuses_disagreeing_duplicates(tmp_path):
    def write(d: Path, n_matched: int):
        d.mkdir(parents=True)
        (d / "pooled.json").write_text(json.dumps({"params": {"log": {"min_snr": 2.0}}, "backend": "log"}),
                                       encoding="utf-8")
        pd.DataFrame({"tag": ["g1"], "n_matched": [n_matched], "n_ours": [10], "n_expert": [12]}).to_csv(
            d / "per_gonad.csv", index=False)
        pd.DataFrame({"tag": ["g1", "g1"], "expert_n": [1, 2], "ours_n": [1, 3]}).to_csv(d / "per_nucleus.csv",
                                                                                      index=False)

    write(tmp_path / "a" / "x", 5)
    write(tmp_path / "b" / "x", 5)
    idx = logo.index_runs([tmp_path / "a", tmp_path / "b"], backend="log")
    sha = bd.params_sha256({"log": {"min_snr": 2.0}})
    assert idx.has(sha, "g1") and idx.scores[sha]["g1"].n_matched == 5
    assert logo.index_runs([tmp_path / "a"], backend="spotmax").scores == {}
    write(tmp_path / "c" / "x", 6)
    with pytest.raises(ValueError, match="disagree"):
        logo.index_runs([tmp_path / "a", tmp_path / "c"])


# ------------------------------------------------------------------ the driver end to end
def test_bench_logo_heldout_and_report(planted_bench, capsys):  # noqa: F811
    b = planted_bench
    mod = _load_script("bench_logo.py")
    common = ["--manifest", str(b["manifest"]), "--cache-root", str(b["root"] / "cache"), "--trials-log",
              str(b["root"] / "trials.jsonl"), "--n-boot", "50"]
    assert mod.main([*common, "folds"]) == 0
    folds = json.loads(capsys.readouterr().out)
    assert [f["held_out"] for f in folds] == ["synth_HS_b", "synth_noHS_a", "synth_val_c"]
    assert folds[2]["val"] == ["synth_noHS_a"] and folds[2]["substitute"] == "synth_noHS_a"
    for name, params in (("log_defaults", "{}"), ("log_snr3", json.dumps({"log": {"min_snr": 3.0}}))):
        rc = mod.main([*common, "heldout", "--detector", name, "--backend", "log", "--params", params,
                       "--out-root", str(b["root"] / "heldout" / name)])
        assert rc == 0, capsys.readouterr().out
    recs = bd.read_jsonl(b["root"] / "trials.jsonl")
    assert len(recs) == 6 and all(len(r["tags"]) == 1 for r in recs)
    assert recs[2]["note"]["logo"]["held_out"] == "synth_val_c" and recs[2]["note"]["logo"]["fold_val"] == ["synth_noHS_a"]
    assert recs[0]["split"] == "train" and recs[2]["split"] == "val" and "per_nucleus_truth_masks_detector_only" in recs[0]["pooled"]
    out = b["root"] / "report"
    rc = mod.main([*common, "report", "--heldout", f"log_defaults={b['root'] / 'heldout' / 'log_defaults'}",
                   f"log_snr3={b['root'] / 'heldout' / 'log_snr3'}", "--out", str(out)])
    assert rc == 0, capsys.readouterr().out
    doc = json.loads((out / "logo_report.json").read_text(encoding="utf-8"))
    assert doc["held_out_gonads"] == ["synth_HS_b", "synth_noHS_a", "synth_val_c"]
    tm = doc["pooled"]["log_defaults"]["truth_masks"]
    ptn = pd.concat([pd.read_csv(b["root"] / "heldout" / "log_defaults" / t / "per_truth_nucleus.csv")
                     for t in doc["held_out_gonads"]])
    assert tm["n_nuclei"] == len(ptn) == 12
    assert tm["stats"]["ccc"]["point"] == pytest.approx(stats.lin_ccc(ptn["expert_n"], ptn["ours_n"]))
    for k in ("germline", "traced", "truth_masks_class_a", "truth_masks_detector_only", "spots_all",
              "spots_without_class_b"):
        assert k in doc["pooled"]["log_snr3"], k
    pair = doc["paired"]["log_defaults - log_snr3"]["truth_masks"]["ccc"]
    assert pair["diff"]["point"] == pytest.approx(pair["a"] - pair["b"])
    pg = pd.read_csv(out / "logo_per_gonad.csv")
    assert len(pg) == 6 and set(pg["detector"]) == {"log_defaults", "log_snr3"} and "tm_ccc" in pg.columns
    assert {"pair", "block", "stat", "diff", "lo", "hi", "excludes_zero"} <= set(pd.read_csv(out / "logo_paired.csv").columns)
    # a missing held-out folder is refused
    assert mod.main([*common, "report", "--heldout", f"x={b['root'] / 'nowhere'}", "--out", str(out)]) == 2
    # --tags pools given gonads, from one bench_detector.py folder or from per-gonad folders: the numbers of
    # the two layouts agree for the same gonads
    det = _load_script("bench_detector.py")
    one = b["root"] / "one_folder"
    assert det.main(["--backend", "log", "--split", "train", "--detector-only", "--no-trial-log", "--out", str(one),
                     "--manifest", str(b["manifest"]), "--cache-root", str(b["root"] / "cache"), "--n-boot", "50"]) == 0
    out2 = b["root"] / "report_tags"
    rc = mod.main([*common, "report", "--tags", "synth_HS_b", "synth_noHS_a", "--heldout", f"single={one}",
                   f"per_gonad={b['root'] / 'heldout' / 'log_defaults'}", "--out", str(out2)])
    assert rc == 0, capsys.readouterr().out
    d2 = json.loads((out2 / "logo_report.json").read_text(encoding="utf-8"))
    assert d2["held_out_gonads"] == ["synth_HS_b", "synth_noHS_a"] and d2["folds"] == []
    s_one, s_per = d2["pooled"]["single"]["truth_masks"]["stats"], d2["pooled"]["per_gonad"]["truth_masks"]["stats"]
    assert s_one["ccc"]["point"] == pytest.approx(s_per["ccc"]["point"])
    assert d2["paired"]["single - per_gonad"]["spots_all"]["f1"]["diff"]["point"] == pytest.approx(0.0)
    assert mod.main([*common, "report", "--tags", "synth_val_c", "--heldout", f"single={one}", "--out", str(out2)]) == 2
    # the test gonad never takes part
    assert "synth_test_d" not in json.dumps(doc["folds"])
    assert load_manifest(b["manifest"])           # the manifest is untouched


# ---------------------------------------------------------- the pre-registered SpotMAX grid
_ABBR = {"threshold_otsu": "otsu", "threshold_li": "li", "threshold_triangle": "tri", "threshold_yen": "yen"}
_PG = {"thresholding_method": "threshold_li", "effect_size_min": 2.5, "z_merge_gap_um": 1.0}
_KEYS = (logo.PRIMARY_OBJECTIVE, logo.PRIMARY_TIE_BREAK)


def test_prereg_grid_is_the_protocol_grid():
    assert logo.PREREG_SPOTMAX_METHODS == ("threshold_otsu", "threshold_li", "threshold_triangle", "threshold_yen")
    assert logo.PREREG_SPOTMAX_EFFECT_SIZES == (2.0, 2.5, 3.0, 3.5, 4.0)
    assert logo.PREREG_SPOTMAX_Z_GAPS == (0.6, 0.8, 1.0)
    grid = logo.prereg_spotmax_grid()
    want = [(m, e, z) for m in ("threshold_otsu", "threshold_li", "threshold_triangle", "threshold_yen")
            for e in (2.0, 2.5, 3.0, 3.5, 4.0) for z in (0.6, 0.8, 1.0)]
    assert len(grid) == 60
    for c in grid:
        assert list(c["params"]) == ["thresholding_method", "effect_size_min", "z_merge_gap_um"]
        assert type(c["params"]["effect_size_min"]) is float and type(c["params"]["z_merge_gap_um"]) is float
    assert [tuple(c["params"].values()) for c in grid] == want
    assert [c["name"] for c in grid] == [f"pg_{_ABBR[m]}_e{e:.1f}_z{z:.1f}" for m, e, z in want]
    assert grid[0]["name"] == "pg_otsu_e2.0_z0.6" and grid[-1]["name"] == "pg_yen_e4.0_z1.0"
    assert len({c["name"] for c in grid}) == 60 and len({bd.params_sha256(c["params"]) for c in grid}) == 60
    assert bt.check_configs(grid) == grid                     # names bench_tune.py takes as folder names
    # keys of the spots: section only; every other key keeps the golden run's (deployed) value
    spots = yaml.safe_load((ROOT / "config" / "config.yaml").read_text(encoding="utf-8"))["spots"]
    for c in grid:
        assert set(c["params"]) <= set(spots)
        section = bd.effective_section(spots, c["params"])
        assert {k for k in section if section[k] != spots[k]} <= set(c["params"])
    deployed = [c for c in grid if c["name"] == "pg_tri_e3.0_z0.8"]
    assert len(deployed) == 1 and deployed[0] is grid[37]
    assert bd.effective_section(spots, deployed[0]["params"]) == spots          # the deployed recipe itself
    kw = bd.pipeline_spot_kwargs(spots)
    assert (kw["spot_radius_um"], kw["gauss_sigma_um"], kw["merge_z_columns"], kw["z_merge_valley_frac"],
            kw["effect_size_metric"], kw["max_spot_candidates"]) == (0.3, 0.08, True, 0.8,
                                                                     "spot_vs_backgr_effect_size_glass", 30000)
    # bench_tune.py shares one SpotMAX detection among the five effect sizes: 12 detections per gonad
    assert len({bt.effect_size_group_key(spots, c["params"]) for c in grid}) == 12


def _bench_tables(tags, seed=3):
    """Per-gonad rows, per-nucleus and per-truth-nucleus tables with every column bench_detector.py pools
    (random counts; two conditions, so the pooling also makes by_split)."""
    rng = np.random.default_rng(seed)
    rows, pn, ptn = [], [], []
    for i, t in enumerate(tags):
        cond = "HS" if i % 2 else "noHS"
        e = rng.integers(0, 9, 40 + 7 * i)
        pn.append(pd.DataFrame({"tag": t, "split": cond, "expert_n": e,
                                "ours_n": np.clip(e + rng.integers(-2, 4, len(e)), 0, None),
                                "traced": rng.random(len(e)) < 0.4}))
        te = rng.integers(0, 12, 11 + 4 * i)
        ta = np.minimum(te, rng.integers(0, 12, len(te)))
        ptn.append(pd.DataFrame({"tag": t, "split": cond, "role": "train", "truth_label": np.arange(1, len(te) + 1),
                                 "n_voxels": 4000, "expert_n": te, "expert_n_class_a": ta, "expert_n_class_b": te - ta,
                                 "ours_n": np.clip(te + rng.integers(-3, 5, len(te)), 0, None),
                                 "ours_n_detector_only": np.nan}))
        row = {"tag": t, "split": cond, "n_expert_spots": 300 + i, "n_our_spots": 280 + 3 * i,
               "matching_method": "hungarian", "truth_source": "synthetic"}
        for prefix in ("", "no_b_", "det_", "tm_spots_", "tm_spots_no_b_", "tm_det_spots_"):
            n_o, n_e = (0, 0) if prefix == "tm_det_spots_" else (int(rng.integers(50, 300)), int(rng.integers(50, 300)))
            m = int(rng.integers(0, min(n_o, n_e) + 1))
            row.update({f"{prefix}n_matched": m, f"{prefix}n_ours": n_o, f"{prefix}n_expert": n_e,
                        f"{prefix}f1": stats.precision_recall_f1(m, n_o, n_e)["f1"]})
        rows.append(row)
    return pd.DataFrame(rows), pd.concat(pn, ignore_index=True), pd.concat(ptn, ignore_index=True)


def _write_run(folder: Path, params: dict, rows, pn, ptn) -> dict:
    """A run folder as bench_tune.py and bench_detector.py write it (pooled.json from detector.pool_scores)."""
    pooled = bd.pool_scores(rows, pn, n_boot=20, seed=1, per_truth_nucleus=ptn)
    pooled.update({"backend": "spotmax", "params": params})
    write_outputs(folder, rows, pn, pooled, {"per_truth_nucleus": ptn})
    return pooled


def test_primary_objective_pools_as_bench_detector(tmp_path):
    tags = ["g_a", "g_b", "g_c", "g_d"]
    rows, pn, ptn = _bench_tables(tags)
    # g_d has no truth nuclei and NaN counts inside R (an older row; the scorer now writes zeros there), so
    # bench_detector.py leaves it out of both blocks
    ptn = ptn[ptn["tag"] != "g_d"].reset_index(drop=True)
    in_r = [*logo.TM_SPOTS_COLUMNS, "tm_spots_f1"]
    rows[in_r] = rows[in_r].astype(float)
    rows.loc[rows["tag"] == "g_d", in_r] = np.nan
    _write_run(tmp_path / "runs" / "x", _PG, rows, pn, ptn)
    sha = bd.params_sha256(_PG)
    idx = logo.index_runs([tmp_path / "runs"], backend="spotmax", truth_masks=True)
    doc = json.loads((tmp_path / "runs" / "x" / "pooled.json").read_text(encoding="utf-8"))
    got = logo.primary_objective_on(idx, sha, tags)
    assert got == bt.objective(doc, keys=_KEYS) and all(math.isfinite(v) for v in got)
    assert got == (round(doc["per_nucleus_truth_masks"]["ccc"]["point"], 3), doc["spots_truth_masks"]["f1"]["point"])
    assert doc["spots_truth_masks"]["n_gonads"] == 3
    d = idx.scores[sha]["g_d"]
    assert len(d.tm_expert) == len(d.tm_ours) == 0 and d.tm_spots == (0, 0, 0) and d.has_truth_masks()
    # any subset of the gonads, in any order: bench_detector.py's point estimates over exactly those gonads
    for sub in (["g_a"], ["g_c", "g_a"], ["g_b", "g_c"], ["g_d", "g_b"], ["g_c", "g_d", "g_a"], ["g_d"]):
        r, n, t = (df[df["tag"].isin(sub)].reset_index(drop=True) for df in (rows, pn, ptn))
        pooled = bd.pool_scores(r, n, n_boot=5, seed=1, per_truth_nucleus=t)
        assert logo.primary_objective_on(idx, sha, sub) == bt.objective(pooled, keys=_KEYS)
        if len(t):
            assert logo.primary_objective_on(idx, sha, sub)[0] == round(stats.lin_ccc(t["expert_n"], t["ours_n"]), 3)
    assert logo.primary_objective_on(idx, sha, ["g_d", "g_b"]) == logo.primary_objective_on(idx, sha, ["g_b"])
    assert logo.primary_objective_on(idx, sha, ["g_d"]) == (-math.inf, -math.inf)        # nothing to pool
    # the D3 objective reads the same folder as before and pools the germline block as bench_detector.py does
    plain = logo.index_runs([tmp_path / "runs"], backend="spotmax")
    assert logo.objective_on(plain, sha, tags) == logo.objective_on(idx, sha, tags) == bt.objective(doc)
    g = plain.scores[sha]["g_a"]
    assert g.tm_expert is None and g.tm_ours is None and g.tm_spots is None
    # missing values rank last: an index without the truth-mask counts, a gonad or a set without a score
    assert logo.primary_objective_on(plain, sha, tags) == (-math.inf, -math.inf)
    assert logo.primary_objective_on(idx, sha, ["g_a", "g_zz"]) == (-math.inf, -math.inf)
    assert logo.primary_objective_on(idx, "no-such-set", ["g_a"]) == (-math.inf, -math.inf)
    assert bt.objective({"per_nucleus_truth_masks": {"ccc": {"point": 0.12345}},
                         "spots_truth_masks": {"f1": {"point": 0.5}}}, keys=_KEYS) == (0.123, 0.5)


def test_index_runs_truth_masks_agree_or_refuse(tmp_path):
    rows, pn, ptn = _bench_tables(["g_a", "g_b"])
    sha = bd.params_sha256(_PG)
    _write_run(tmp_path / "a" / "x", _PG, rows, pn, ptn)
    _write_run(tmp_path / "b" / "x", _PG, rows, pn, ptn)                         # the same scores again
    _write_run(tmp_path / "c" / "x", _PG, rows.drop(columns=list(logo.TM_SPOTS_COLUMNS)), pn, ptn)
    (tmp_path / "c" / "x" / "per_truth_nucleus.csv").unlink()                    # a folder without them
    _write_run(tmp_path / "d" / "x", _PG, rows, pn, ptn.assign(ours_n=ptn["ours_n"] + 1))    # other counts
    alone = logo.index_runs([tmp_path / "c"], truth_masks=True).scores[sha]["g_a"]
    assert alone.tm_expert is None and alone.tm_ours is None and alone.tm_spots is None
    # a folder written before the section 5.2 columns: the truth-mask table alone is not enough for the
    # primary objective, so the pair counts as missing there (and only there)
    _write_run(tmp_path / "e" / "x", _PG, rows.drop(columns=list(logo.TM_SPOTS_COLUMNS)), pn, ptn)
    old = logo.index_runs([tmp_path / "e"], truth_masks=True)
    s = old.scores[sha]["g_a"]
    assert s.tm_expert is not None and s.tm_spots is None and not s.has_truth_masks()
    one = [{"name": "x", "params": _PG}]
    assert logo.missing_scores(old, one, ["g_a", "g_b"], truth_masks=True) == [{**one[0], "tags": ["g_a", "g_b"]}]
    assert logo.missing_scores(old, one, ["g_a", "g_b"]) == []
    g = logo.index_runs([tmp_path / "c", tmp_path / "a", tmp_path / "b"], truth_masks=True).scores[sha]["g_a"]
    assert np.array_equal(g.tm_ours, ptn.loc[ptn["tag"] == "g_a", "ours_n"].to_numpy(float))    # taken from a
    assert g.tm_spots == tuple(int(rows.loc[0, c]) for c in logo.TM_SPOTS_COLUMNS)
    with pytest.raises(ValueError, match="truth-mask"):
        logo.index_runs([tmp_path / "a", tmp_path / "d"], truth_masks=True)
    # the default index, the D3 replay's, reads none of it, so there the two folders still agree
    d3 = logo.index_runs([tmp_path / "a", tmp_path / "d"]).scores[sha]["g_a"]
    assert d3.tm_expert is None and d3.tm_spots is None and d3.n_matched == int(rows.loc[0, "n_matched"])


def test_primary_objective_reads_bench_tune_output(planted_bench):  # noqa: F811
    """The run folders bench_tune.py itself writes (log backend, so no SpotMAX run): every score carries the
    truth-mask counts, and the primary objective equals bench_tune.py's own pooling and summary."""
    b = planted_bench
    tune = _load_script("bench_tune.py")
    configs = [{"name": "defaults", "params": {}}, {"name": "snr3", "params": {"log": {"min_snr": 3.0}}},
               {"name": "gate2_nomerge", "params": {"log": {"pre_gate_sigmas": 2.0}, "merge_z_columns": False}}]
    cfile = b["root"] / "configs.json"
    cfile.write_text(json.dumps({"configs": configs}), encoding="utf-8")
    out = b["root"] / "tune"
    assert tune.main(["--backend", "log", "--configs", str(cfile), "--split", "train", "--out-root", str(out),
                      "--no-trial-log", "--manifest", str(b["manifest"]), "--cache-root", str(b["root"] / "cache"),
                      "--n-boot", "50"]) == 0
    tags = ["synth_HS_b", "synth_noHS_a"]
    idx = logo.index_runs([out], backend="log", truth_masks=True)
    summary = pd.read_csv(out / "summary.csv").set_index("name")
    for c in configs:
        sha = bd.params_sha256(c["params"])
        pooled = json.loads((out / c["name"] / "pooled.json").read_text(encoding="utf-8"))
        pg = pd.read_csv(out / c["name"] / "per_gonad.csv").set_index("tag")
        got = logo.primary_objective_on(idx, sha, tags)
        assert got == bt.objective(pooled, keys=_KEYS) and all(math.isfinite(v) for v in got)
        assert got[0] == round(summary.loc[c["name"], "tm_ccc"], 3)
        for t in tags:
            s = idx.scores[sha][t]
            assert s.has_truth_masks() and len(s.tm_expert) > 0
            assert s.tm_spots == tuple(int(pg.loc[t, k]) for k in logo.TM_SPOTS_COLUMNS)
    assert logo.missing_scores(idx, configs, tags, truth_masks=True) == []
    best = max(configs, key=lambda c: logo.primary_objective_on(idx, bd.params_sha256(c["params"]), tags))
    assert logo.select_grid(idx, configs, tags)["selected"]["name"] == best["name"]


def _tm(ccc, f1_hits=50, germline_shift=1.0) -> logo.GonadScore:
    """A gonad score with a chosen truth-mask CCC (the method's count is the truth count plus a constant c,
    so CCC = 1 / (1 + c^2 / (2 var)); None gives NaN, both counts constant) and an F1 inside R of
    f1_hits / 100 (None: no spot in R, F1 NaN). The germline counts are off by ``germline_shift``."""
    if ccc is None:
        te, to = np.full(6, 2.0), np.full(6, 2.0)
    else:
        te, to = _EXPERT, _EXPERT + math.sqrt(2 * _EXPERT.var() * (1 / ccc - 1))
    return logo.GonadScore(expert=_EXPERT, ours=_EXPERT + germline_shift, n_matched=50, n_ours=100, n_expert=100,
                           tm_expert=te, tm_ours=to, tm_spots=(0, 0, 0) if f1_hits is None else (f1_hits, 100, 100))


def _put(idx: logo.RunIndex, config: dict, tag: str, score: logo.GonadScore) -> None:
    sha = bd.params_sha256(config["params"])
    idx.params[sha], idx.names[sha] = config["params"], config["name"]
    idx.scores.setdefault(sha, {})[tag] = score


def test_select_grid_rounding_tie_breaks_and_missing():
    configs = [{"name": n, "params": {"effect_size_min": float(i)}} for i, n in enumerate("adbcgefh")]
    by = {c["name"]: c for c in configs}
    val = ["v1", "v2"]
    plan = {"a": _tm(0.8504, 40, germline_shift=0.0),      # the best germline (D3) objective, not the best here
            "d": _tm(0.8496, 50), "b": _tm(0.8496, 50),    # equal in every number: the grid's order decides
            "c": _tm(0.90, 10), "g": _tm(0.8501, None), "e": _tm(None, 95), "h": _tm(0.80, 99),
            "f": logo.GonadScore(expert=_EXPERT, ours=_EXPERT + 1.0, n_matched=50, n_ours=100, n_expert=100)}
    idx = logo.RunIndex()
    for name, score in plan.items():
        for t in (["v1"] if name == "h" else val):            # h has no score on v2 yet
            _put(idx, by[name], t, score)
    sel = logo.select_grid(idx, configs, val)
    assert [x["name"] for x in sel["ranking"]] == ["c", "d", "b", "a", "g", "e", "f", "h"]
    # f's scores lack the truth-mask counts (a run folder without them) and h has no score on v2: both are
    # pairs still to run, so nothing is selected yet
    assert sel["selected"] is None and sel["missing"] == [{"name": "f", "params": by["f"]["params"], "tags": val},
                                                          {"name": "h", "params": by["h"]["params"], "tags": ["v2"]}]
    assert all(set(x) == {"name", "params", "objective"} for x in sel["ranking"])
    obj = {x["name"]: x["objective"] for x in sel["ranking"]}
    # the rounding: a's CCC is above b's (0.8504 against 0.8496), both round to 0.850 and the F1 decides
    assert stats.lin_ccc(plan["a"].tm_expert, plan["a"].tm_ours) > stats.lin_ccc(plan["b"].tm_expert, plan["b"].tm_ours)
    assert obj["a"][0] == obj["b"][0] == obj["d"][0] == 0.85 and obj["c"][0] == 0.9
    assert obj["a"][1] == pytest.approx(0.4) and obj["b"] == obj["d"] and obj["b"][1] == pytest.approx(0.5)
    # NaN or missing: below every finite value, the other value still ordering them
    assert obj["g"][0] == 0.85 and obj["g"][1] == -math.inf and obj["e"][0] == -math.inf
    assert obj["e"][1] == pytest.approx(0.95) and obj["f"] == obj["h"] == (-math.inf, -math.inf)
    _put(idx, by["h"], "v2", plan["h"])
    for t in val:                       # counts that give a NaN CCC and no spot in R: f is complete but last
        _put(idx, by["f"], t, _tm(None, None))
    sel = logo.select_grid(idx, configs, val)
    assert sel["selected"] == {"name": "c", "params": by["c"]["params"]} and sel["missing"] == []
    assert [x["name"] for x in sel["ranking"]] == ["c", "d", "b", "a", "g", "h", "e", "f"]
    assert logo.rank_on(idx, configs, val)[0]["name"] == "a"            # the D3 objective would pick a
    # replay_grid: the rule on a fold's val gonads, with the returns bench_logo.py reads; train is not scored
    r = logo.replay_grid(idx, configs, ["t_unscored"], val)
    assert r["complete"] and r["missing"] == [] and r["selected"] == "c" and r["selected_params"] == by["c"]["params"]
    assert list(r["val_objective"]) == [x["name"] for x in sel["ranking"]]
    assert r["val_objective"]["c"] == pytest.approx([0.9, 0.1])
    assert r["finalists"] == [c["name"] for c in configs] and r["n_sets"] == 8 and "stage" not in r
    del idx.scores[bd.params_sha256(by["e"]["params"])]["v1"]
    r = logo.replay_grid(idx, configs, [], val)
    assert not r["complete"] and r["stage"] == "grid" and "selected" not in r and "selected_params" not in r
    assert r["missing"] == [{"name": "e", "params": by["e"]["params"], "tags": ["v1"]}]
    # the best set loses the truth-mask counts on one gonad: the selection waits for that pair instead of
    # passing the set over for d
    _put(idx, by["e"], "v1", plan["e"])
    _put(idx, by["c"], "v2", logo.GonadScore(expert=_EXPERT, ours=_EXPERT, n_matched=50, n_ours=100, n_expert=100))
    sel = logo.select_grid(idx, configs, val)
    assert sel["selected"] is None and sel["missing"] == [{"name": "c", "params": by["c"]["params"], "tags": ["v2"]}]
    assert [x["name"] for x in sel["ranking"]] == ["d", "b", "a", "g", "h", "e", "c", "f"]
    # a test gonad is refused
    with pytest.raises(bd.TestSetGuardError):
        logo.select_grid(idx, configs, ["v1", "n2rad51_HS_male_001"])
    with pytest.raises(bd.TestSetGuardError):
        logo.replay_grid(idx, configs, ["n2rad51_noHS_herm_002"], val)
    # an index without the truth-mask counts (built without truth_masks) has every pair still to run
    bare = logo.RunIndex()
    for c in configs:
        for t in val:
            _put(bare, c, t, plan["f"])
    sel = logo.select_grid(bare, configs, val)
    assert sel["selected"] is None and [m["tags"] for m in sel["missing"]] == [val] * len(configs)
    # nothing to select on: no set with a finite CCC, or no selection gonad
    flat = logo.RunIndex()
    for c in configs:
        for t in val:
            _put(flat, c, t, _tm(None, 50))
    with pytest.raises(ValueError, match="finite truth-mask CCC"):
        logo.select_grid(flat, configs, val)
    with pytest.raises(ValueError, match="selection gonad"):
        logo.replay_grid(idx, configs, ["v1"], [])


def _real_split_entries() -> list[Entry]:
    rg = ("n2rad51_HS_male_000", "n2rad51_HS_male_002")
    tags = [f"n2rad51_{c}_{s}_{i:03d}" for c in ("noHS", "HS") for s in ("herm", "male") for i in range(3)]
    return [Entry(tag=t, results_dir=Path("."), split="",
                  extra={"condition": "HS" if "_HS_" in t else "noHS", "spot_truth_usable": t not in rg}) for t in tags]


def test_replay_grid_selects_on_each_folds_own_val_gonads():
    folds = logo.logo_folds(_real_split_entries(), require_cache=False)
    grid = logo.prereg_spotmax_grid()
    # each gonad favours one set (truth-mask CCC 0.9 against 0.5); the train gonads' favourites come first in
    # the grid, so a rule that also pooled a fold's train gonads would pick one of them
    fav = {"n2rad51_noHS_male_000": 0, "n2rad51_noHS_male_001": 1, "n2rad51_HS_herm_000": 5,
           "n2rad51_noHS_herm_000": 6, "n2rad51_HS_herm_001": 10, "n2rad51_noHS_herm_001": 11,
           "n2rad51_noHS_herm_002": 2}                         # a test gonad: in the index, in no fold
    idx = logo.RunIndex()
    for i, c in enumerate(grid):
        for t, k in fav.items():
            _put(idx, c, t, _tm(0.9 if i == k else 0.5))
    d = logo.designated_split(folds)
    assert d == {"train": ["n2rad51_HS_herm_000", "n2rad51_noHS_herm_000", "n2rad51_noHS_male_000",
                           "n2rad51_noHS_male_001"], "val": ["n2rad51_HS_herm_001", "n2rad51_noHS_herm_001"]}
    # the designated val gonads select where they are among the five; for a held-out val gonad the first train
    # gonad of its condition takes its place (the two favourites then tie and the grid's order decides)
    want = {"n2rad51_HS_herm_000": "pg_otsu_e3.5_z0.8", "n2rad51_HS_herm_001": "pg_otsu_e2.5_z1.0",
            "n2rad51_noHS_herm_000": "pg_otsu_e3.5_z0.8", "n2rad51_noHS_herm_001": "pg_otsu_e3.0_z0.6",
            "n2rad51_noHS_male_000": "pg_otsu_e3.5_z0.8", "n2rad51_noHS_male_001": "pg_otsu_e3.5_z0.8"}
    got = {f["held_out"]: logo.replay_grid(idx, grid, f["train"], f["val"]) for f in folds}
    assert {t: r["selected"] for t, r in got.items()} == want
    assert all(r["complete"] and r["n_sets"] == 60 and len(r["val_objective"]) == 60 for r in got.values())
    assert logo.replay_grid(idx, grid, d["train"], d["val"])["selected"] == "pg_otsu_e3.5_z0.8"
    hs = next(f for f in folds if f["held_out"] == "n2rad51_HS_herm_000")
    assert logo.select_grid(idx, grid, hs["train"] + hs["val"])["selected"]["name"] == grid[0]["name"]
    for f in folds:
        assert "n2rad51_noHS_herm_002" not in f["train"] + f["val"] and f["held_out"] not in f["train"] + f["val"]
    # a train gonad's score is never needed, a substitute's is
    del idx.scores[bd.params_sha256(grid[40]["params"])]["n2rad51_noHS_male_001"]
    del idx.scores[bd.params_sha256(grid[41]["params"])]["n2rad51_HS_herm_000"]
    got = {f["held_out"]: logo.replay_grid(idx, grid, f["train"], f["val"]) for f in folds}
    assert [t for t, r in got.items() if not r["complete"]] == ["n2rad51_HS_herm_001"]
    assert got["n2rad51_HS_herm_001"]["missing"] == [{"name": grid[41]["name"], "params": grid[41]["params"],
                                                       "tags": ["n2rad51_HS_herm_000"]}]
    with pytest.raises(bd.TestSetGuardError):
        logo.replay_grid(idx, grid, [], ["n2rad51_noHS_herm_002"])


def _grid_runs(root: Path, tags, special: dict, grid) -> None:
    """Run folders of every grid set on ``tags`` holding what ``index_runs`` reads. Truth-mask CCC 0.5 and F1
    inside R 0.5 except ``special[tag][grid index] = (CCC, F1 hits of 100)``. The germline counts and the
    whole-field matching favour the first set of the grid, and set 20 over set 21, so a rule reading them
    would choose otherwise."""
    for i, c in enumerate(grid):
        pg, pn, ptn = [], [], []
        for t in tags:
            ccc, hits = special.get(t, {}).get(i, (0.5, 50))
            s = _tm(ccc, hits, germline_shift=0.0 if i == 0 else 1.0)
            pg.append({"tag": t, "n_matched": {0: 90, 20: 60}.get(i, 40), "n_ours": 100, "n_expert": 100,
                       **dict(zip(logo.TM_SPOTS_COLUMNS, s.tm_spots))})
            pn.append(pd.DataFrame({"tag": t, "expert_n": s.expert, "ours_n": s.ours}))
            ptn.append(pd.DataFrame({"tag": t, "truth_label": np.arange(1, len(s.tm_expert) + 1),
                                     "expert_n": s.tm_expert, "ours_n": s.tm_ours}))
        write_outputs(root / c["name"], pd.DataFrame(pg), pd.concat(pn, ignore_index=True),
                      {"backend": "spotmax", "params": c["params"]},
                      {"per_truth_nucleus": pd.concat(ptn, ignore_index=True)})


def test_bench_logo_grid_configs_and_select_prereg_grid(planted_bench, capsys, monkeypatch):  # noqa: F811
    b = planted_bench
    mod = _load_script("bench_logo.py")
    tune = _load_script("bench_tune.py")
    common = ["--manifest", str(b["manifest"]), "--cache-root", str(b["root"] / "cache"), "--trials-log",
              str(b["root"] / "trials.jsonl"), "--n-boot", "50"]
    grid = logo.prereg_spotmax_grid()
    g = b["root"] / "grid"
    # grid-configs: a configuration file that bench_tune.py takes as it is
    assert mod.main([*common, "grid-configs", "--out", str(g / "configs.json")]) == 0
    assert tune.load_configs(g / "configs.json") == bt.check_configs(grid) == grid
    capsys.readouterr()
    assert tune.main(["--backend", "spotmax", "--configs", str(g / "configs.json"), "--split", "val", "--chunks", "4",
                      "--list", "--out-root", str(g / "val"), "--manifest", str(b["manifest"]),
                      "--cache-root", str(b["root"] / "cache")]) == 0
    assert not (g / "val").exists()                               # --list runs nothing and writes nothing
    tasks = [ln for ln in capsys.readouterr().out.splitlines() if ln.startswith("  task synth_val_c")]
    assert len(tasks) == 4 and all(ln.count("'pg_") == 15 for ln in tasks)     # three shared detections each

    # the val gonad favours sets 20 and 21 equally on the CCC and 21 on the F1 inside R; the train gonads
    # favour sets 7 (synth_noHS_a, the substitute for the held-out val gonad) and 3 (synth_HS_b)
    special = {"synth_val_c": {20: (0.9, 50), 21: (0.9, 60)}, "synth_noHS_a": {7: (0.9, 50)},
               "synth_HS_b": {3: (0.9, 50)}}
    _grid_runs(g / "val", ["synth_val_c"], special, grid)
    sel = g / "selection.json"
    select = [*common, "select", "--backend", "spotmax", "--search", "prereg-grid", "--out", str(sel), "--runs"]
    assert mod.main([*select, str(g / "val")]) == 3
    doc = json.loads(sel.read_text(encoding="utf-8"))
    by = {f["held_out"]: f for f in doc["folds"]}
    assert not doc["complete"] and doc["search"] == "prereg-grid" and doc["configs"] == grid
    assert doc["objective"] == logo.PRIMARY_OBJECTIVE_TEXT
    assert by["synth_HS_b"]["selected"] == by["synth_noHS_a"]["selected"] == grid[21]["name"]
    assert doc["designated"]["train"] == ["synth_HS_b", "synth_noHS_a"] and doc["designated"]["val"] == ["synth_val_c"]
    assert doc["designated"]["selected"] == grid[21]["name"]
    lack = by["synth_val_c"]                  # its selection gonad is the substitute synth_noHS_a, not run yet
    assert not lack["complete"] and lack["stage"] == "grid" and lack["val"] == ["synth_noHS_a"]
    assert len(doc["missing"]) == 60
    assert {(m["fold"], m["stage"], tuple(m["tags"])) for m in doc["missing"]} == {
        ("synth_val_c", "grid", ("synth_noHS_a",))}
    assert doc["to_run"] == [{"name": c["name"], "params": c["params"], "tags": ["synth_noHS_a"]} for c in grid]
    heldout = [*common, "heldout", "--detector", "spotmax_grid", "--backend", "spotmax", "--selection", str(sel),
               "--out-root", str(g / "heldout")]
    assert mod.main(heldout) == 2                                 # an incomplete selection is refused

    _grid_runs(g / "train", ["synth_HS_b", "synth_noHS_a"], special, grid)
    runs = [str(g / "val"), str(g / "train")]
    assert mod.main([*select, *runs]) == 0
    doc = json.loads(sel.read_text(encoding="utf-8"))
    by = {f["held_out"]: f for f in doc["folds"]}
    assert set(doc) == {"backend", "search", "runs", "objective", "configs", "folds", "designated", "complete",
                        "missing", "to_run"}
    assert doc["complete"] and doc["missing"] == [] and doc["to_run"] == []
    assert {t: f["selected"] for t, f in by.items()} == {"synth_HS_b": grid[21]["name"],
                                                          "synth_noHS_a": grid[21]["name"],
                                                          "synth_val_c": grid[7]["name"]}
    assert by["synth_val_c"]["substitute"] == "synth_noHS_a"
    assert by["synth_val_c"]["selected_params"] == grid[7]["params"]
    assert doc["designated"]["selected"] == grid[21]["name"]
    assert doc["designated"]["selected_params"] == grid[21]["params"]
    v = by["synth_noHS_a"]["val_objective"]
    assert len(v) == 60 and list(v)[:2] == [grid[21]["name"], grid[20]["name"]]
    assert v[grid[21]["name"]] == pytest.approx([0.9, 0.6])
    assert by["synth_HS_b"]["finalists"] == [c["name"] for c in grid]
    assert "synth_test_d" not in json.dumps(doc["folds"]) + json.dumps(doc["designated"])
    assert mod.main([*common, "select", "--backend", "log", "--search", "prereg-grid", "--out", str(sel),
                     "--runs", *runs]) == 2                       # the grid is SpotMAX's

    # heldout --selection reads it as it reads a D3 selection (bench_detector.py stubbed: no detector runs)
    calls = []
    monkeypatch.setattr(mod, "_bench_detector", lambda: SimpleNamespace(main=lambda argv: calls.append(argv) or 0))
    assert mod.main(heldout) == 0
    assert sorted(a[a.index("--tags") + 1] for a in calls) == ["synth_HS_b", "synth_noHS_a", "synth_val_c"]
    for a in calls:
        f = by[a[a.index("--tags") + 1]]
        assert json.loads(a[a.index("--params") + 1]) == f["selected_params"] and a[a.index("--split") + 1] == f["role"]
        note = json.loads(a[a.index("--note") + 1])["logo"]
        assert note["selected"] == f["selected"] and note["finalists"] == [c["name"] for c in grid]
        assert note["fold_val"] == f["val"] and note["selection_file"] == str(sel)

    # the D3 replay is still the default and writes its own layout (here: its S1 sets were never run)
    d3 = g / "selection_d3.json"
    assert mod.main([*common, "select", "--backend", "spotmax", "--out", str(d3), "--runs", *runs]) == 3
    dd = json.loads(d3.read_text(encoding="utf-8"))
    assert set(dd) == {"backend", "runs", "objective", "folds", "complete", "missing"}
    assert dd["objective"] == bt.OBJECTIVE_TEXT and {f["stage"] for f in dd["folds"]} == {"S1"}


def _real_split_manifest(root: Path, *, with_val: bool = True) -> tuple[Path, Path]:
    """A manifest of the real split's train and val gonads (no test gonad; roles from the pre-registration)
    and a cache folder holding only the file the fold rule checks, so `bench_logo.py select` builds the real
    folds without image data. ``with_val=False`` leaves the val gonads out."""
    root.mkdir(parents=True)
    cache = root / "cache"
    entries = []
    for e in _real_split_entries():
        role = bd.preregistered_role(e.tag)
        if role == "test" or (role == "val" and not with_val):
            continue
        (cache / e.tag).mkdir(parents=True)
        (cache / e.tag / bd.CACHE_FILES["meta"]).write_text("{}", encoding="utf-8")
        entries.append({"tag": e.tag, "results_dir": str(root / "results" / e.tag), **e.extra})
    manifest = root / "manifest.json"
    manifest.write_text(json.dumps(entries), encoding="utf-8")
    return manifest, cache


def test_bench_logo_select_prereg_grid_on_the_real_folds(tmp_path, capsys):
    """`select --search prereg-grid` on the real split's six folds (synthetic run folders, no image data). The
    two folds that hold out a val gonad each need their own substitute, `to_run` merges both, and the two
    train males are never a selection gonad. A score without the truth-mask counts holds the selection up
    until a folder with them is added. Disagreeing folders and a split without val gonads are refused."""
    mod = _load_script("bench_logo.py")
    grid = logo.prereg_spotmax_grid()
    manifest, cache = _real_split_manifest(tmp_path / "split")
    common = ["--manifest", str(manifest), "--cache-root", str(cache), "--trials-log", str(tmp_path / "trials.jsonl")]
    hs_v, nohs_v = "n2rad51_HS_herm_001", "n2rad51_noHS_herm_001"
    hs_sub, nohs_sub = "n2rad51_HS_herm_000", "n2rad51_noHS_herm_000"
    # both val gonads favour set 30, each substitute its own set (two gonads pool to the harmonic mean of
    # their CCCs here: 0.655 for 0.5 and 0.95, against 0.643 for 0.9 and 0.5)
    special = {hs_v: {30: (0.9, 50)}, nohs_v: {30: (0.9, 50)}, hs_sub: {12: (0.95, 50)}, nohs_sub: {44: (0.95, 50)}}
    _grid_runs(tmp_path / "val", [hs_v, nohs_v], special, grid)
    sel = tmp_path / "selection.json"
    select = [*common, "select", "--backend", "spotmax", "--search", "prereg-grid", "--out", str(sel), "--runs"]
    assert mod.main([*select, str(tmp_path / "val")]) == 3
    doc = json.loads(sel.read_text(encoding="utf-8"))
    by = {f["held_out"]: f for f in doc["folds"]}
    assert list(by) == [hs_sub, hs_v, nohs_sub, nohs_v, "n2rad51_noHS_male_000", "n2rad51_noHS_male_001"]
    assert [t for t, f in by.items() if not f["complete"]] == [hs_v, nohs_v]
    assert by[hs_v]["substitute"] == hs_sub and by[nohs_v]["substitute"] == nohs_sub
    assert len(doc["missing"]) == 120
    assert {(m["fold"], tuple(m["tags"])) for m in doc["missing"]} == {(hs_v, (hs_sub,)), (nohs_v, (nohs_sub,))}
    assert doc["to_run"] == [{"name": c["name"], "params": c["params"], "tags": [hs_sub, nohs_sub]} for c in grid]
    assert doc["designated"]["complete"] and doc["designated"]["selected"] == grid[30]["name"]
    # run folders of the two substitutes alone complete every fold
    _grid_runs(tmp_path / "sub", [hs_sub, nohs_sub], special, grid)
    runs = [str(tmp_path / "val"), str(tmp_path / "sub")]
    assert mod.main([*select, *runs]) == 0
    doc = json.loads(sel.read_text(encoding="utf-8"))
    want = {t: grid[30]["name"] for t in by} | {hs_v: grid[12]["name"], nohs_v: grid[44]["name"]}
    assert {f["held_out"]: f["selected"] for f in doc["folds"]} == want
    assert doc["complete"] and doc["to_run"] == [] and doc["designated"]["selected"] == grid[30]["name"]

    # a run folder without per_truth_nucleus.csv: its pairs are to run again and nothing is selected, until the
    # same scores with the counts come from another folder
    name = grid[30]["name"]
    shutil.copytree(tmp_path / "val" / name, tmp_path / "rerun" / name)
    (tmp_path / "val" / name / "per_truth_nucleus.csv").unlink()
    assert mod.main([*select, *runs]) == 3
    doc = json.loads(sel.read_text(encoding="utf-8"))
    assert doc["to_run"] == [{"name": name, "params": grid[30]["params"], "tags": [hs_v, nohs_v]}]
    assert not doc["designated"]["complete"] and not any(f["complete"] for f in doc["folds"])
    runs.append(str(tmp_path / "rerun"))
    assert mod.main([*select, *runs]) == 0
    assert {f["held_out"]: f["selected"] for f in json.loads(sel.read_text(encoding="utf-8"))["folds"]} == want

    # refused (exit code 2, nothing written): a folder that disagrees on the truth-mask counts, and a split
    # without val gonads (no fold has a selection gonad)
    shutil.copytree(tmp_path / "rerun" / name, tmp_path / "bad" / name)
    ptn = pd.read_csv(tmp_path / "bad" / name / "per_truth_nucleus.csv")
    ptn.assign(ours_n=ptn["ours_n"] + 1).to_csv(tmp_path / "bad" / name / "per_truth_nucleus.csv", index=False)
    before = sel.read_text(encoding="utf-8")
    capsys.readouterr()
    assert mod.main([*select, *runs, str(tmp_path / "bad")]) == 2
    assert "REFUSED" in capsys.readouterr().out and sel.read_text(encoding="utf-8") == before
    manifest2, cache2 = _real_split_manifest(tmp_path / "no_val", with_val=False)
    assert mod.main(["--manifest", str(manifest2), "--cache-root", str(cache2), "select", "--backend", "spotmax",
                     "--search", "prereg-grid", "--out", str(tmp_path / "sel2.json"), "--runs", *runs]) == 2
    assert "selection gonad" in capsys.readouterr().out and not (tmp_path / "sel2.json").exists()
