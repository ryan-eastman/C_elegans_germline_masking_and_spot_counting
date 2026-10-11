"""The frozen recipe of net3d (germquant.spots.net3d.train.recipe_from, ``germquant train-spots
--recipe-from``): protocol sections 6.1 and 6.2 retrain a selected configuration on other gonads with no
validation, its epoch count and its threshold. On the tiny planted cache of test_net3d (CPU): retraining by
the recipe on the same gonads reproduces the selected checkpoint bit for bit, the train-plus-val model trains
on a val gonad and never opens the test gonad or the region-growing one, and the refusals. Skips without
torch."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

torch = pytest.importorskip("torch")

from test_net3d import tiny_bench  # noqa: E402,F401  (the fixture)

TRAIN_ARGS = ["--patch", "16,32,32", "--seed", "3", "--batch-size", "2", "--steps-per-epoch", "2",
              "--device", "cpu", "--base-channels", "4", "--levels", "3"]


def _base(b, out):
    return ["train-spots", "--cache", str(b["root"] / "cache"), "--split-file", str(b["split"]), "--out", str(out),
            "--manifest", str(b["manifest"]), "--prep-dir", str(b["root"] / "prep")]


def test_recipe_reproduces_the_selected_checkpoint(tiny_bench, monkeypatch, capsys):  # noqa: F811
    from germquant import cli
    from germquant.benchmark import detector as bd
    from germquant.spots.net3d.infer import file_sha256, read_model_file
    from germquant.spots.net3d.train import recipe_from

    b = tiny_bench
    sel_dir = b["root"] / "selected"
    rc = cli.main([*_base(b, sel_dir), "--epochs", "3", *TRAIN_ARGS, "--thresholds", "0.02,0.3,0.6"])
    assert rc == 0, capsys.readouterr().out
    logd = json.loads((sel_dir / "training_log.json").read_text(encoding="utf-8"))
    k, thr = logd["selection"]["epoch"], logd["selection"]["threshold"]

    kw, recipe = recipe_from(sel_dir)
    assert kw["run_epochs"] == k and kw["fixed_threshold"] == thr and kw["epochs"] == 3
    assert kw["patch"] == (16, 32, 32) and kw["seed"] == 3 and kw["base_channels"] == 4 and kw["use_dna"] is True
    assert recipe["from_model_sha256"] == logd["model_sha256"] and recipe["selected_epoch"] == k
    assert recipe["from_train_gonads"] == ["synth_noHS_a", "synth_HS_b"]

    # retrained by the recipe on the same gonads, without validation: the selected checkpoint, bit for bit
    rep = b["root"] / "reproduced"
    rc = cli.main([*_base(b, rep), "--recipe-from", str(sel_dir), "--train-tags", "synth_noHS_a", "synth_HS_b"])
    text = capsys.readouterr().out
    assert rc == 0, text
    assert "frozen recipe" in text
    want = read_model_file(sel_dir / "checkpoints" / f"epoch_{k:03d}.pt")["state_dict"]
    got_doc = read_model_file(rep / "model.pt")
    assert set(want) == set(got_doc["state_dict"])
    assert all(torch.equal(want[n], got_doc["state_dict"][n]) for n in want)
    assert got_doc["peaks"]["threshold"] == thr and got_doc["epoch"] == k
    assert got_doc["training"]["frozen_recipe"] is True and got_doc["training"]["val_gonads"] == []
    rlog = json.loads((rep / "training_log.json").read_text(encoding="utf-8"))
    card = json.loads((rep / "model_card.json").read_text(encoding="utf-8"))
    assert len(rlog["epochs"]) == k and all("val" not in e for e in rlog["epochs"])
    assert rlog["selection"]["n_configurations"] == 0 and rlog["recipe"]["from_model_sha256"] == logd["model_sha256"]
    assert rlog["settings"]["run_epochs"] == k and rlog["settings"]["epochs"] == 3
    assert [e["lr_end"] for e in rlog["epochs"]] == [e["lr_end"] for e in logd["epochs"][:k]]   # the same schedule
    assert card["model_sha256"] == file_sha256(rep / "model.pt") and card["threshold"] == thr
    assert card["training_gonad_tags"] == ["synth_noHS_a", "synth_HS_b"] and card["val_gonads"] == []
    assert card["recipe"]["selected_threshold"] == thr and card["epochs"]["run"] == k
    assert any("frozen recipe" in s for s in card["limitations"])

    # train plus val: the val gonad is trained on, the test and region-growing gonads are never opened
    opened = []
    real_load = bd.load_cache

    def spy_load(tag, *a, **kw2):
        opened.append(tag)
        return real_load(tag, *a, **kw2)

    monkeypatch.setattr(bd, "load_cache", spy_load)
    tv = b["root"] / "trainval"
    rc = cli.main([*_base(b, tv), "--recipe-from", str(sel_dir), "--train-tags", "synth_noHS_a", "synth_HS_b",
                   "synth_val_c"])
    assert rc == 0, capsys.readouterr().out
    assert set(opened) == {"synth_noHS_a", "synth_HS_b", "synth_val_c"}
    tcard = json.loads((tv / "model_card.json").read_text(encoding="utf-8"))
    assert tcard["training_gonad_tags"] == ["synth_noHS_a", "synth_HS_b", "synth_val_c"]
    assert tcard["test_gonads_never_opened"] == ["synth_test_d"] and tcard["threshold"] == thr


def test_bench_logo_net3d_folds(tiny_bench, capsys):  # noqa: F811
    """scripts/bench_logo.py net3d-folds on the tiny cache: a fold trains on the recipe's other gonads (in
    the recipe's order: its train gonads, then its val gonads) and its held-out gonad is scored with the
    detector-only count, one trial record per fold with the fold in its note."""
    import importlib.util

    from germquant import cli
    from germquant.benchmark import detector as bd
    from germquant.spots.net3d.infer import file_sha256, read_model_file

    b = tiny_bench
    sel_dir = b["root"] / "selected"
    assert cli.main([*_base(b, sel_dir), "--epochs", "2", *TRAIN_ARGS, "--thresholds", "0.02,0.3"]) == 0
    capsys.readouterr()
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("bench_logo", root / "scripts" / "bench_logo.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    common = ["--manifest", str(b["manifest"]), "--cache-root", str(b["root"] / "cache"), "--trials-log",
              str(b["root"] / "trials.jsonl"), "--n-boot", "20"]
    rc = mod.main([*common, "net3d-folds", "--tags", "synth_val_c", "synth_HS_b", "--recipe-from", str(sel_dir),
                   "--models-root", str(b["root"] / "folds"), "--out-root", str(b["root"] / "heldout"),
                   "--log-dir", str(b["root"] / "logs"), "--split-file", str(b["split"]), "--device", "cpu",
                   "--prep-dir", str(b["root"] / "prep")])
    logs = "".join(p.read_text(encoding="utf-8") for p in sorted((b["root"] / "logs").glob("*.log")))
    assert rc == 0, logs
    for held, trained in (("synth_val_c", ["synth_noHS_a", "synth_HS_b"]), ("synth_HS_b", ["synth_noHS_a", "synth_val_c"])):
        card = json.loads((b["root"] / "folds" / f"fold_{held}" / "model_card.json").read_text(encoding="utf-8"))
        assert card["training_gonad_tags"] == trained and card["val_gonads"] == []
        assert card["recipe"]["from_model_sha256"] == file_sha256(sel_dir / "model.pt")
        pg = pd.read_csv(b["root"] / "heldout" / held / "per_gonad.csv")
        assert pg["tag"].tolist() == [held] and pg["backend"].tolist() == ["net3d"] and "tm_det_ccc" in pg.columns
    doc = read_model_file(b["root"] / "folds" / "fold_synth_val_c" / "model.pt")
    assert doc["peaks"]["threshold"] == json.loads((sel_dir / "training_log.json").read_text(encoding="utf-8"))["best_threshold"]
    recs = bd.read_jsonl(b["root"] / "trials.jsonl")
    assert [r["note"]["logo"]["held_out"] for r in recs] == ["synth_val_c", "synth_HS_b"]
    assert recs[0]["note"]["logo"]["trained_on"] == ["synth_noHS_a", "synth_HS_b"]
    assert recs[0]["note"]["logo"]["model_sha256"] == file_sha256(b["root"] / "folds" / "fold_synth_val_c" / "model.pt")


def test_recipe_refusals(tiny_bench, capsys):  # noqa: F811
    from germquant import cli
    from germquant.benchmark import detector as bd
    from germquant.spots.net3d.train import choose_gonads, resolve_gonads

    b = tiny_bench
    gon = resolve_gonads(b["split"], b["manifest"])
    assert choose_gonads(gon) == (["synth_noHS_a", "synth_HS_b"], ["synth_val_c"])
    assert choose_gonads(gon, ["synth_val_c", "synth_noHS_a"]) == (["synth_val_c", "synth_noHS_a"], [])
    assert choose_gonads(gon, ["synth_noHS_a"], []) == (["synth_noHS_a"], [])
    with pytest.raises(bd.TestSetGuardError):
        choose_gonads(gon, ["synth_noHS_a", "synth_test_d"])
    with pytest.raises(ValueError, match="cannot be used"):
        choose_gonads(gon, ["synth_rg_e"])
    with pytest.raises(ValueError, match="both a training and a validation"):
        choose_gonads(gon, ["synth_noHS_a", "synth_val_c"], ["synth_val_c"])
    with pytest.raises(ValueError, match="named twice"):
        choose_gonads(gon, ["synth_noHS_a", "synth_noHS_a"])
    # no val gonad and no fixed epoch count and threshold: refused before training
    with pytest.raises(ValueError, match="no val gonad"):
        cli.main([*_base(b, b["root"] / "x"), "--epochs", "1", *TRAIN_ARGS, "--val-tags"])
    # a threshold is fixed only without val gonads
    with pytest.raises(ValueError, match="fixed_threshold"):
        cli.main([*_base(b, b["root"] / "y"), "--epochs", "1", *TRAIN_ARGS, "--threshold", "0.4"])
    # --recipe-from fixes the epoch count and the threshold itself, and selects nothing on val gonads
    assert cli.main([*_base(b, b["root"] / "z"), "--recipe-from", str(b["root"]), "--threshold", "0.4"]) == 2
    assert "REFUSED" in capsys.readouterr().out
    assert cli.main([*_base(b, b["root"] / "z"), "--recipe-from", str(b["root"]), "--val-tags", "synth_val_c"]) == 2
    assert "drop --val-tags" in capsys.readouterr().out
