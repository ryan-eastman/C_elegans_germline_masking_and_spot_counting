"""The parameter search runner on the detector bench harness (germquant.benchmark.tuning,
scripts/bench_tune.py), on synthetic data:

- the SpotMAX effect-size filter applied after a run with effect_size_min 0 reproduces a run made with
  that minimum, table for table (directly and through run_configs);
- the objective and the ranking (rounding, the F1 tie-break, a missing value last) and the chunking;
- bench_tune.py on the planted-foci gonads of test_bench_detector: one trial record per parameter set
  with its tuning block, the same scores and detections as bench_detector.py for the same parameters,
  the process pool, the budget guard, the top-K selection for val and the test guard;
- config/profiles/bench_log_tuned.yaml and bench_spotmax_tuned.yaml extend config/config.yaml and
  change only spots: keys, and as bench parameters they give the section the profile gives.
Nothing here needs microscopy data or a GPU.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_bench_detector import _load_script, planted_bench  # noqa: F401  (fixture: planted-foci gonads)
from test_end_to_end import _synthetic_stack

from germquant.benchmark import detector as bd
from germquant.benchmark import tuning as bt

ROOT = Path(__file__).resolve().parents[1]
TIMING = ["detect_s", "score_s", "elapsed_s"]


def _args(b, *extra):
    r = b["root"]
    return ["--manifest", str(b["manifest"]), "--cache-root", str(r / "cache"), "--trials-log",
            str(r / "trials.jsonl"), "--n-boot", "50", *extra]


def _write_configs(path, configs):
    path.write_text(json.dumps({"configs": configs}), encoding="utf-8")
    return path


# ------------------------------------------------------------------ the SpotMAX effect-size filter
def test_filter_effect_size_reproduces_a_spotmax_run():
    pytest.importorskip("spotmax")
    from skimage.measure import label as cc_label

    from germquant.spots.detect import detect_spots_spotmax

    stack = _synthetic_stack()
    labels = cc_label(stack.data[0] > 1000).astype(np.int32)
    rad = stack.data[2]
    kw = {"thresholding_method": "threshold_triangle"}
    base = detect_spots_spotmax(rad, labels, stack.spacing, effect_size_min=0.0, **kw)
    es = base[0]["effect_size"].dropna().to_numpy()
    assert len(base[0]) > 3 and len(es) > 3
    for e in (0.0, 1.0, float(np.median(es)), 3.0, float(es.max()) + 1.0):
        direct = detect_spots_spotmax(rad, labels, stack.spacing, effect_size_min=e, **kw)
        got = bt.filter_effect_size(*base, e)
        pd.testing.assert_frame_equal(got[0], direct[0])
        pd.testing.assert_frame_equal(got[1], direct[1])
    # the filter keeps a missing effect size and filters nothing when no spot has one
    ps, pn = base
    ps_nan = ps.assign(effect_size=np.nan)
    got = bt.filter_effect_size(ps_nan, pn, 5.0)
    assert len(got[0]) == len(ps)
    one_nan = ps.copy()
    one_nan.loc[0, "effect_size"] = np.nan
    got = bt.filter_effect_size(one_nan, pn, float(es.max()) + 1.0)
    assert len(got[0]) == 1 and got[0]["spot_id"].tolist() == [0]
    assert int(got[1]["n_spots"].sum()) == 1 and len(got[1]) == len(pn)


def test_run_configs_shares_the_spotmax_detection_without_changing_a_result(planted_bench):  # noqa: F811
    pytest.importorskip("spotmax")
    b = planted_bench
    configs = [{"name": "es2", "params": {"effect_size_min": 2.0}}, {"name": "es3", "params": {}},
               {"name": "es4_r25", "params": {"effect_size_min": 4.0, "spot_radius_um": 0.25}}]
    kw = {"cache_root": b["root"] / "cache", "manifest": b["manifest"], "verbose": False}
    shared = bt.run_configs("synth_noHS_a", "spotmax", configs, **kw)
    alone = bt.run_configs("synth_noHS_a", "spotmax", configs, effect_size_shortcut=False, **kw)
    assert [r["name"] for r in shared] == ["es2", "es3", "es4_r25"]
    assert [r["shared_detection"] is not None for r in shared] == [True, True, False]
    assert shared[0]["shared_detection"]["effect_size_min"] == 0.0
    assert all(r["shared_detection"] is None for r in alone)
    for s, a in zip(shared, alone):
        pd.testing.assert_frame_equal(s["per_spot"], a["per_spot"])
        pd.testing.assert_frame_equal(s["per_nucleus"], a["per_nucleus"])
        rs = {k: v for k, v in s["row"].items() if k not in TIMING}
        ra = {k: v for k, v in a["row"].items() if k not in TIMING}
        assert rs.keys() == ra.keys()
        for k in rs:
            same = rs[k] == ra[k] or (isinstance(rs[k], float) and math.isnan(rs[k]) and math.isnan(ra[k]))
            assert same, (s["name"], k, rs[k], ra[k])


# ------------------------------------------------------------------------- objective and chunks
def test_objective_rank_and_chunks():
    def pooled(ccc, f1):
        return {"per_nucleus_all": {"ccc": {"point": ccc}}, "spots_all": {"f1": {"point": f1}}}

    items = [{"name": "a", "pooled": pooled(0.8504, 0.40)}, {"name": "b", "pooled": pooled(0.8496, 0.50)},
             {"name": "c", "pooled": pooled(0.9, 0.10)}, {"name": "d", "pooled": pooled(None, 0.90)},
             {"name": "e", "pooled": pooled(0.8501, 0.45)}, {"name": "f", "pooled": pooled(float("nan"), 0.95)}]
    # a missing CCC ranks below every scored one; among those the F1 still orders them
    assert [it["name"] for it in bt.rank(items)] == ["c", "b", "e", "a", "f", "d"]
    assert bt.objective(pooled(0.12345, 0.5)) == (0.123, 0.5)
    assert bt.objective({}) == (-math.inf, -math.inf)

    configs = [{"name": f"c{i}", "params": {"effect_size_min": e, "spot_radius_um": r}}
               for i, (r, e) in enumerate([(0.3, 2), (0.3, 3), (0.3, 4), (0.25, 2), (0.25, 3), (0.35, 3)])]
    golden = {"spot_radius_um": 0.3, "effect_size_min": 3.0}
    chunks = bt.chunk_configs(configs, 2, backend="spotmax", spots_config=golden)
    assert [[c["name"] for c in ch] for ch in chunks] == [["c0", "c1", "c2"], ["c3", "c4", "c5"]]
    flat = bt.chunk_configs(configs, 4, backend="log")
    assert sum(len(ch) for ch in flat) == 6 and len(flat) == 4
    assert [c["name"] for ch in flat for c in ch] == [c["name"] for c in configs]
    assert bt.chunk_configs(configs, 10, backend="spotmax", spots_config=golden)[0][0]["name"] == "c0"
    assert len(bt.chunk_configs(configs, 10, backend="spotmax", spots_config=golden)) == 3
    assert bt.effect_size_group_key(golden, {}) == bt.effect_size_group_key(golden, {"effect_size_min": 2.0})
    assert bt.effective_effect_size_min(golden, {}) == 3.0
    for bad in ([], [{"params": {}}], [{"name": "a b"}], [{"name": "x"}, {"name": "x"}],
                [{"name": "x", "params": {"backend": "log"}}]):
        with pytest.raises(ValueError):
            bt.check_configs(bad)


# ------------------------------------------------------------------------------- the runner
def test_bench_tune_log_on_planted_stacks(planted_bench, capsys):  # noqa: F811
    b = planted_bench
    tune = _load_script("bench_tune.py")
    det = _load_script("bench_detector.py")
    configs = [{"name": "defaults", "params": {}},
               {"name": "snr_huge", "params": {"log": {"min_snr": 1e9}}},   # nothing passes
               {"name": "gate2_nomerge", "params": {"log": {"pre_gate_sigmas": 2.0}, "merge_z_columns": False}}]
    cfile = _write_configs(b["root"] / "configs.json", configs)
    out = b["root"] / "tune"
    rc = tune.main(["--backend", "log", "--configs", str(cfile), "--split", "train", "--out-root", str(out),
                    "--search", "synth_log", "--stage", "A", *_args(b)])
    text = capsys.readouterr().out
    assert rc == 0, text
    lines = bd.read_jsonl(b["root"] / "trials.jsonl")
    assert [rec["tuning"]["name"] for rec in lines] == ["defaults", "snr_huge", "gate2_nomerge"]
    for rec, c in zip(lines, configs):
        assert rec["params"] == c["params"] and rec["params_sha256"] == bd.params_sha256(c["params"])
        assert rec["split"] == "train" and rec["tags"] == ["synth_noHS_a", "synth_HS_b"] and rec["final"] is False
        assert rec["tuning"]["search"] == "synth_log" and rec["tuning"]["stage"] == "A"
        assert rec["tuning"]["effect_size_shortcut"] is False and rec["tuning"]["complete"] is True
        assert "f1" in rec["pooled"]["spots_all"] and set(rec["per_gonad"]) == {"synth_noHS_a", "synth_HS_b"}
    # the same scores, pooled blocks and detections as bench_detector.py for the same parameter set
    assert det.main(["--backend", "log", "--params", json.dumps(configs[2]["params"]), "--split", "train",
                     "--out", str(b["root"] / "det"), "--no-trial-log", "--manifest", str(b["manifest"]),
                     "--cache-root", str(b["root"] / "cache"), "--n-boot", "50"]) == 0
    a = pd.read_csv(out / "gate2_nomerge" / "per_gonad.csv").drop(columns=TIMING, errors="ignore")
    d = pd.read_csv(b["root"] / "det" / "per_gonad.csv").drop(columns=TIMING, errors="ignore")
    pd.testing.assert_frame_equal(a, d)
    pa = json.loads((out / "gate2_nomerge" / "pooled.json").read_text(encoding="utf-8"))
    pdet = json.loads((b["root"] / "det" / "pooled.json").read_text(encoding="utf-8"))
    for k in bd.POOLED_KEYS + ("by_split",):
        assert pa[k] == pdet[k], k
    for t in ("synth_noHS_a", "synth_HS_b"):
        pd.testing.assert_frame_equal(pd.read_csv(out / "gate2_nomerge" / "detections" / f"{t}__spots.csv"),
                                      pd.read_csv(b["root"] / "det" / "detections" / f"{t}__spots.csv"))
    # the parameters reach the detector: three different spot sets
    n = [pd.read_csv(out / c["name"] / "per_gonad.csv")["n_our_spots"].sum() for c in configs]
    assert len(set(n)) > 1
    # the summary ranks every set by the objective, best first
    summ = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    keys = [bt.objective(c["pooled"]) for c in summ["configs"]]
    assert keys == sorted(keys, reverse=True) and len(keys) == 3
    table = pd.read_csv(out / "summary.csv")
    assert table["rank"].tolist() == [1, 2, 3] and table["name"].tolist() == [c["name"] for c in summ["configs"]]
    assert json.loads((out / "configs.json").read_text(encoding="utf-8"))["configs"] == configs

    # the process pool with two tasks per gonad gives the same scores
    out2 = b["root"] / "tune_pool"
    assert tune.main(["--backend", "log", "--configs", str(cfile), "--split", "train", "--out-root", str(out2),
                      "--jobs", "2", "--chunks", "2", "--no-trial-log", *_args(b)]) == 0
    for c in configs:
        x = pd.read_csv(out / c["name"] / "per_gonad.csv").drop(columns=TIMING, errors="ignore")
        y = pd.read_csv(out2 / c["name"] / "per_gonad.csv").drop(columns=TIMING, errors="ignore")
        pd.testing.assert_frame_equal(x, y)
    assert len(bd.read_jsonl(b["root"] / "trials.jsonl")) == 3          # --no-trial-log left the log alone

    # the budget counts distinct parameter sets per search
    more = _write_configs(b["root"] / "more.json", configs[:1] + [{"name": "snr4", "params": {"log": {"min_snr": 4.0}}},
                                                                  {"name": "snr5", "params": {"log": {"min_snr": 5.0}}}])
    capsys.readouterr()
    assert tune.main(["--backend", "log", "--configs", str(more), "--split", "train", "--out-root",
                      str(b["root"] / "o"), "--search", "synth_log", "--budget", "4", "--list", *_args(b)]) == 2
    assert "budget" in capsys.readouterr().out
    assert tune.main(["--backend", "log", "--configs", str(more), "--split", "train", "--out-root",
                      str(b["root"] / "o"), "--search", "synth_log", "--budget", "5", "--list", *_args(b)]) == 0
    assert tune.main(["--backend", "log", "--configs", str(more), "--split", "train", "--out-root",
                      str(b["root"] / "o"), "--budget", "5", "--list", *_args(b)]) == 2       # no --search

    # the best two of the train search go to val
    capsys.readouterr()
    val = b["root"] / "val"
    assert tune.main(["--backend", "log", "--top-from", str(out / "summary.json"), "--top", "2", "--split", "val",
                      "--out-root", str(val), "--search", "synth_log", "--stage", "val", *_args(b)]) == 0
    best_two = [c["name"] for c in summ["configs"][:2]]
    vsum = json.loads((val / "summary.json").read_text(encoding="utf-8"))
    assert sorted(c["name"] for c in vsum["configs"]) == sorted(best_two) and vsum["tags"] == ["synth_val_c"]
    vlines = bd.read_jsonl(b["root"] / "trials.jsonl")[3:]
    assert [rec["split"] for rec in vlines] == ["val", "val"] and all(rec["tuning"]["stage"] == "val" for rec in vlines)


def test_bench_tune_refuses_test_gonads(planted_bench, capsys):  # noqa: F811
    b = planted_bench
    tune = _load_script("bench_tune.py")
    cfile = _write_configs(b["root"] / "c.json", [{"name": "defaults", "params": {}}])
    out = b["root"] / "never"
    assert tune.main(["--backend", "log", "--configs", str(cfile), "--split", "test", "--out-root", str(out),
                      *_args(b)]) == 2
    assert "REFUSED" in capsys.readouterr().out
    assert tune.main(["--backend", "log", "--configs", str(cfile), "--split", "train", "--tags", "synth_test_d",
                      "--out-root", str(out), *_args(b)]) == 2
    with pytest.raises(SystemExit):                                       # there is no --final option
        tune.main(["--backend", "log", "--configs", str(cfile), "--split", "test", "--final",
                   "--out-root", str(out), *_args(b)])
    assert not out.exists() and not (b["root"] / "trials.jsonl").exists()
    with pytest.raises(bd.TestSetGuardError):
        bt.run_configs("synth_test_d", "log", [{"name": "x", "params": {}}], cache_root=b["root"] / "cache",
                       manifest=b["manifest"])
    with pytest.raises(bd.TestSetGuardError):
        bt.run_configs("n2rad51_noHS_male_002", "log", [{"name": "x", "params": {}}], cache_root=b["root"] / "nowhere",
                       manifest=b["manifest"])
    with pytest.raises(ValueError, match="unusable"):
        bt.run_configs("synth_rg_e", "log", [{"name": "x", "params": {}}], cache_root=b["root"] / "cache",
                       manifest=b["manifest"])


# ------------------------------------------------------------------------------ the tuned profiles
@pytest.mark.parametrize("name, backend", [("bench_log_tuned", "log"), ("bench_spotmax_tuned", "spotmax")])
def test_tuned_profiles_change_only_spots_keys(name, backend):
    import yaml

    from germquant.config import load_config
    from germquant.spots.detect import SPOTMAX_KEYS
    from germquant.spots.log_detect import LOG_DEFAULTS, log_kwargs

    path = ROOT / "config" / "profiles" / f"{name}.yaml"
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    assert raw["extends"] == "config.yaml" and set(raw) == {"extends", "spots"}
    base = load_config(ROOT / "config" / "config.yaml").as_dict()
    prof = load_config(path).as_dict()
    assert {k: v for k, v in prof.items() if k != "spots"} == {k: v for k, v in base.items() if k != "spots"}
    spots = prof["spots"]
    assert spots.get("backend", "spotmax") == backend
    overlay = {k: v for k, v in raw["spots"].items() if k != "backend"}
    # as bench parameters over the base section, the overlay gives exactly the profile's section
    assert bd.effective_section(base["spots"], overlay) == {k: v for k, v in spots.items() if k != "backend"}
    changed = {k for k in set(spots) | set(base["spots"]) if spots.get(k) != base["spots"].get(k)}
    if backend == "log":
        assert set(raw["spots"]["log"]) <= set(LOG_DEFAULTS) and "log" in changed
        kw = log_kwargs(spots)                                    # every key valid for the detector
        from germquant.spots.log_detect import _resolve

        _resolve(kw["diameter_um"], kw["min_snr"], kw["min_separation_um"], kw["annulus_inner_um"],
                 kw["annulus_outer_um"], kw["pre_gate_sigmas"])
    else:
        assert changed <= set(SPOTMAX_KEYS)
        assert spots["thresholding_method"] in ("threshold_triangle", "threshold_otsu", "threshold_li")
