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


ADDED = "n2rad51_20260105a_HS_male_001"


def _plant_addition(b) -> Path:
    """One planted gonad under a tag of protocol changelog entry 9, in the manifest as an addition (train, HS
    male), and an additions record naming it (the record validates: its role is fixed in code)."""
    from test_bench_detector import _write_gonad

    from germquant.validate import benchmark_split as bs

    g = _write_gonad(b["root"], ADDED, "train", "HS", (3, 5, 2, 4), 15)
    entry = dict(g["entry"], sex="male", addition=True, session="2026-01-05")
    doc = json.loads(b["manifest"].read_text(encoding="utf-8"))
    doc["gonads"].append(entry)
    b["manifest"].write_text(json.dumps(doc), encoding="utf-8")
    rec = {"name": "synthetic additions", "changelog_entry": 9,
           "registered_split": {"roles_sha256": bs.roles_sha256()},
           "roles": {"train": [ADDED], "val": [], "test": []},
           "gonads": [{"tag": ADDED, "role": "train", "session": "2026-01-05", "condition": "HS", "sex": "male",
                       "stem": entry["stem"],
                       "nd2": {"path": entry["nd2"], "size_bytes": 10, "sha256": "a" * 64},
                       "ims": {"path": entry["ims"], "size_bytes": 10, "sha256": "b" * 64, "region_growing": False},
                       "spot_truth": {"source": "ims", "status": "usable"}}]}
    bs.validate_additions(rec)
    path = b["root"] / "additions.json"
    path.write_text(json.dumps(rec), encoding="utf-8")
    return path


def test_recipe_with_the_additions_of_entry_9(tiny_bench, capsys):  # noqa: F811
    """With --additions an addition is a train gonad after the split's own (the recipe's order, then the
    additions in record order); the card names both sessions and the heat-shock male. Without the record the
    addition cannot be named, and a record that calls an addition test or names a registered gonad is refused.
    bench_logo.py net3d-fold trains a fold on the recipe's gonads plus the addition and refuses the folds without
    --additions."""
    import importlib.util

    from germquant import cli
    from germquant.spots.net3d.train import resolve_gonads

    b = tiny_bench
    adds = _plant_addition(b)
    sel_dir = b["root"] / "selected"
    assert cli.main([*_base(b, sel_dir), "--epochs", "2", *TRAIN_ARGS, "--thresholds", "0.02,0.3"]) == 0
    capsys.readouterr()
    gon = resolve_gonads(b["split"], b["manifest"], adds)
    assert gon["train"] == ["synth_noHS_a", "synth_HS_b", ADDED] and gon["val"] == ["synth_val_c"]
    assert gon["sessions"][ADDED] == "2026-01-05" and gon["sessions"]["synth_noHS_a"] is None
    order = ["synth_noHS_a", "synth_HS_b", "synth_val_c", ADDED]
    tv = b["root"] / "trainval_entry9"
    rc = cli.main([*_base(b, tv), "--recipe-from", str(sel_dir), "--train-tags", *order, "--additions", str(adds)])
    assert rc == 0, capsys.readouterr().out
    card = json.loads((tv / "model_card.json").read_text(encoding="utf-8"))
    assert card["training_gonad_tags"] == order
    lim = " ".join(card["limitations"])
    assert "2 imaging sessions (2025-11-05, 2026-01-05" in lim and f"Heat-shock males among the training gonads: ['{ADDED}']" in lim
    logd = json.loads((tv / "training_log.json").read_text(encoding="utf-8"))
    assert logd["settings"]["additions_file"]["path"].endswith("additions.json")
    # without the record the addition is not a gonad of the split
    with pytest.raises(ValueError, match="cannot be used"):
        cli.main([*_base(b, b["root"] / "x"), "--recipe-from", str(sel_dir), "--train-tags", *order])
    bad = json.loads(adds.read_text(encoding="utf-8"))
    bad["gonads"][0]["role"] = "test"
    adds.with_name("bad.json").write_text(json.dumps(bad), encoding="utf-8")
    with pytest.raises(ValueError, match="never val or test"):
        resolve_gonads(b["split"], b["manifest"], adds.with_name("bad.json"))
    bad["gonads"][0].update(role="train", tag="n2rad51_HS_male_001")
    adds.with_name("bad.json").write_text(json.dumps(bad), encoding="utf-8")
    with pytest.raises(ValueError, match="registered gonad"):
        resolve_gonads(b["split"], b["manifest"], adds.with_name("bad.json"))

    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("bench_logo", root / "scripts" / "bench_logo.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    common = ["--manifest", str(b["manifest"]), "--cache-root", str(b["root"] / "cache"), "--trials-log",
              str(b["root"] / "trials.jsonl"), "--n-boot", "20"]
    fold = ["net3d-fold", "--held-out", "synth_HS_b", "--recipe-from", str(sel_dir), "--models-root",
            str(b["root"] / "folds"), "--out-root", str(b["root"] / "heldout"), "--split-file", str(b["split"]),
            "--device", "cpu", "--prep-dir", str(b["root"] / "prep")]
    assert mod.main([*common, *fold]) == 2
    assert "pass --additions" in capsys.readouterr().out
    assert mod.main([*common, *fold, "--additions", str(adds)]) == 0, capsys.readouterr().out
    fcard = json.loads((b["root"] / "folds" / "fold_synth_HS_b" / "model_card.json").read_text(encoding="utf-8"))
    assert fcard["training_gonad_tags"] == ["synth_noHS_a", "synth_val_c", ADDED]
    rec = bd_read(b["root"] / "trials.jsonl")[-1]
    assert rec["note"]["logo"]["trained_on"] == ["synth_noHS_a", "synth_val_c", ADDED]
    assert rec["note"]["logo"]["additions_file"] == str(adds)


def bd_read(path):
    from germquant.benchmark import detector as bd

    return bd.read_jsonl(path)
