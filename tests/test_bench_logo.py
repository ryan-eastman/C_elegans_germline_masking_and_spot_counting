"""Leave-one-gonad-out (germquant.benchmark.logo, scripts/bench_logo.py; docs/BENCHMARK_PROTOCOL.md
section 6.2): the fold rule (designated val gonads, a same-condition substitute by tag order for a held-out
val gonad, never a test gonad or a gonad without usable spot truth), the D3 stage generators, the staged
replay of the D3 search on a synthetic score index (including the parameter sets a fold still has to run),
and the driver end to end on the planted-foci gonads of test_bench_detector with the log backend
(held-out scoring, the trial records' notes, the pooled report and the paired comparison). Synthetic only.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_bench_detector import planted_bench  # noqa: F401  (the fixture)

from germquant.benchmark import detector as bd
from germquant.benchmark import logo
from germquant.benchmark import stats
from germquant.benchmark.manifest import Entry, load_manifest

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
