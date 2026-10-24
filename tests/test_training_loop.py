"""The nucleus training loop (germquant.segment.training, docs/TIMELAPSE.md "Segmenting without DAPI"):
`germquant export-training`, `germquant train-nuclei` and `germquant register-model`, plus the registry lookup the
segment stage makes for a registered model name.

The export runs on the synthetic live series of test_live_profile (two timepoints) and on a fixed 3-channel run.
Training is exercised with the Cellpose seams replaced (`_build_net`, `_train_seg`, `_predict`), so the default
suite never loads weights, never trains the 300 M parameter network and never writes its 1.2 GB file; the
bookkeeping around it (the split that never separates the slices of one image or series, the losses, the sha256
of every pair, the hold-out F1 per orientation, the card) is what is checked. A real random-initialisation
training run on the CPU (about 90 s and 1.2 GB) is opt-in: GERMQUANT_TRAIN_SMOKE=1.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tifffile
from test_end_to_end import patched_reader  # noqa: F401  (fixture: the 3-channel fixed gonad)
from test_live_profile import ND2, STEM, live_cfg, live_series_stacks, patch_live_series
from test_segment_backends import registry  # noqa: F401  (fixture: a private backend registry)

from germquant import cli, models_registry as R, pipeline, provenance
from germquant.config import load_config
from germquant.segment import backends as B
from germquant.segment import training as TR
from germquant.timelapse import process_series

ROOT = Path(__file__).resolve().parents[1]
T2 = 2


@pytest.fixture
def live_series_results(tmp_path, monkeypatch):
    """A finished two-timepoint live series under tmp_path/series (the live profile, classical segmenter)."""
    patch_live_series(monkeypatch, T2)
    out = tmp_path / "series"
    res = process_series(ND2, live_cfg(), out, timepoints="all")
    assert res["done"] == [0, 1]
    return out


def _sha(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---- export-training ----------------------------------------------------------------------------------

def test_export_training_from_a_series(tmp_path, live_series_results):
    out = live_series_results
    runs = TR.find_runs(out)
    assert [r.image_id for r in runs] == [f"{STEM}_t000", f"{STEM}_t001"]
    assert all(r.series_id == STEM and r.channel_role == "central_element" and r.channel_index == 1 for r in runs)
    assert [r.t_index for r in runs] == [0, 1] and all(r.xy_stride == 1 and r.z_range is None for r in runs)
    assert all(r.source_file.endswith(ND2) and r.nuclei_path is not None for r in runs)

    pairs = tmp_path / "pairs"
    man = TR.export_training_pairs(out, pairs, orientations=("xy", "xz", "yz"), slices_per_image=2, seed=1)
    assert list(man.columns) == TR.MANIFEST_COLS
    assert len(man) == 2 * 3 * 2 and man.groupby(["image_id", "orientation"]).size().eq(2).all()
    assert (pairs / "manifest.csv").exists() and (pairs / "export.json").exists()
    assert set(man["zone"].fillna("")) == {""}                       # no staging stage: no zone tags
    assert (man["pixel_um"] == 0.2).all() and (man["channel_role"] == "central_element").all()
    assert set(man.loc[man["orientation"] == "xy", "z_zoom"]) == {1.0} and set(man.loc[man["orientation"] != "xy", "z_zoom"]) == {2.0}

    stacks = live_series_stacks(T2)
    for r in man.to_dict("records"):
        im = tifffile.imread(pairs / r["image"])
        lab = tifffile.imread(pairs / r["labels"])
        assert im.dtype == np.uint16 and lab.dtype == np.int32 and im.shape == lab.shape == (r["height"], r["width"])
        ids = np.unique(lab)
        assert int(np.count_nonzero(ids)) == r["n_labels"] >= 1
        k = int(r["t_index"])
        labels3d = tifffile.imread(out / f"t{k:03d}" / f"{r['image_id']}__nuclei_labels.tif")
        gfp = stacks[k].data[1]
        assert set(ids[ids > 0]) <= set(np.unique(labels3d))
        if r["orientation"] == "xy":                                   # the raw plane, as uint16, and the run's labels
            assert im.shape == (80, 320)
            assert np.array_equal(im, np.clip(np.rint(gfp[r["slice_index"]]), 0, 65535).astype(np.uint16))
            assert np.array_equal(lab, labels3d[r["slice_index"]])
        else:                                                          # z stretched by dz / dy = 2 to the xy pixel size
            assert im.shape[0] == 40 and im.shape[1] == (320 if r["orientation"] == "xz" else 80)
            src_lab = labels3d[:, r["slice_index"], :] if r["orientation"] == "xz" else labels3d[:, :, r["slice_index"]]
            assert set(np.unique(lab)) == set(np.unique(src_lab))      # nearest neighbour keeps every id, invents none
            assert np.array_equal(lab[::2], src_lab) or np.array_equal(lab[1::2], src_lab)
        assert Path(r["image"]).name == f"{r['image_id']}_{r['orientation']}_{r['slice_index']:04d}.tif"
    exp = json.loads((pairs / "export.json").read_text(encoding="utf-8"))
    assert exp["seed"] == 1 and exp["stopped_at_cap"] is False and len(exp["runs"]) == 2

    # the draw is seeded per image: the same seed picks the same slices, another seed may not
    again = TR.export_training_pairs(out, tmp_path / "pairs_again", orientations=("xy",), slices_per_image=2, seed=1)
    assert list(again["slice_index"]) == list(man.loc[man["orientation"] == "xy", "slice_index"])
    # exporting into the same folder again replaces that image's rows instead of duplicating them
    man2 = TR.export_training_pairs(out, pairs, orientations=("xy",), slices_per_image=1, seed=5)
    assert len(man2) == 2 and man2["image_id"].nunique() == 2 and man2.groupby("image_id").size().eq(1).all()
    # the size cap stops the export early and says so
    capped = TR.export_training_pairs(out, tmp_path / "capped", orientations=("xy", "xz"), slices_per_image=3, max_bytes=1)
    assert len(capped) == 1 and json.loads((tmp_path / "capped" / "export.json").read_text(encoding="utf-8"))["stopped_at_cap"] is True
    with pytest.raises(ValueError):
        TR.export_training_pairs(out, tmp_path / "bad", orientations=("zz",))
    with pytest.raises(FileNotFoundError):
        TR.export_training_pairs(tmp_path / "empty_results", tmp_path / "bad2")
    # the command line
    cli_out = tmp_path / "cli_pairs"
    assert cli.main(["export-training", str(out), "--out", str(cli_out), "--orientations", "xy,yz",
                     "--slices-per-image", "1", "--seed", "3"]) == 0
    cm = pd.read_csv(cli_out / "manifest.csv")
    assert len(cm) == 4 and set(cm["orientation"]) == {"xy", "yz"}


def test_export_training_from_a_fixed_run_with_zone_tags(tmp_path, patched_reader):  # noqa: F811
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    out = tmp_path / "fixed"
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, out, prov=None)
    iid = res["image_id"]
    runs = TR.find_runs(out)
    assert len(runs) == 1 and runs[0].image_id == iid == runs[0].series_id and runs[0].t_index is None
    assert runs[0].channel_role == "dna" and runs[0].channel_index == 0        # the fixed recipe segments DAPI

    # a zone column (what the staging stage writes) becomes the slice tag: the most frequent zone of its nuclei
    nuc_path = out / f"{iid}__nuclei.csv"
    nuc = pd.read_csv(nuc_path)
    nuc["zone"] = np.where(nuc["centroid_x_um"] < 20.0, "early", "mid")
    nuc.to_csv(nuc_path, index=False)
    man = TR.export_training_pairs(out, tmp_path / "pairs", orientations=("xy", "xz"), slices_per_image=2, seed=0)
    assert len(man) == 4 and man["t_index"].isna().all() and (man["series_id"] == iid).all()
    assert (man["channel_role"] == "dna").all() and (man["channel_index"] == 0).all()
    zones = dict(zip(nuc["nucleus_id"].astype(int), nuc["zone"]))
    for r in man.to_dict("records"):
        lab = tifffile.imread(tmp_path / "pairs" / r["labels"])
        ids = [int(v) for v in np.unique(lab) if v > 0]
        tally = pd.Series([zones[i] for i in ids if i in zones]).value_counts()
        assert r["zone"] == sorted(tally[tally == tally.max()].index)[0]
    # the dna plane of the source stack is what was written
    stack = pipeline.read_stack("x")
    xy = man[man["orientation"] == "xy"].iloc[0]
    assert np.array_equal(tifffile.imread(tmp_path / "pairs" / xy["image"]),
                          np.clip(np.rint(stack.data[0][int(xy["slice_index"])]), 0, 65535).astype(np.uint16))
    # load_pairs reads the manifest back; without one it pairs the files by name
    pairs = TR.load_pairs(tmp_path / "pairs")
    assert len(pairs) == 4 and {p["series_id"] for p in pairs} == {iid} and all(p["zone"] in ("early", "mid") for p in pairs)
    (tmp_path / "pairs" / "manifest.csv").unlink()
    by_name = TR.load_pairs(tmp_path / "pairs")
    assert len(by_name) == 4 and {(p["orientation"], p["slice_index"]) for p in by_name} == {(p["orientation"], p["slice_index"]) for p in pairs}
    assert all(p["image_id"] == iid for p in by_name)


# ---- the split -------------------------------------------------------------------------------------------

def _fake_pairs(n_series=3, n_t=3, orients=("xy", "xz", "yz")):
    pairs = []
    for s in range(n_series):
        for t in range(n_t):
            for o in orients:
                iid = f"series{s}_t{t:03d}"
                pairs.append({"image": f"image/{iid}_{o}_0001.tif", "labels": f"labels/{iid}_{o}_0001.tif",
                              "image_id": iid, "series_id": f"series{s}", "orientation": o, "slice_index": 1, "zone": ""})
    return pairs


def test_split_holds_out_whole_series_or_images_never_slices():
    pairs = _fake_pairs()
    train, val, rec = TR.split_pairs(pairs, by="series", val_frac=0.34, seed=0)
    assert rec["n_groups"] == 3 and len(rec["val_groups"]) == 1 and len(train) + len(val) == len(pairs)
    vg = set(rec["val_groups"])
    assert all(p["series_id"] in vg for p in val) and not any(p["series_id"] in vg for p in train)
    assert len(val) == 9                                                   # every timepoint and orientation of that series
    train_i, val_i, rec_i = TR.split_pairs(pairs, by="image", val_frac=0.2, seed=1)
    assert rec_i["n_groups"] == 9 and len(rec_i["val_groups"]) == 2 and len(val_i) == 6
    assert {p["image_id"] for p in val_i} == set(rec_i["val_groups"])
    assert not {p["image_id"] for p in train_i} & set(rec_i["val_groups"])
    # at least one group is held out when there are two or more, none when there is one, never all of them
    _, v, r = TR.split_pairs(_fake_pairs(n_series=2), by="series", val_frac=0.01, seed=0)
    assert len(r["val_groups"]) == 1 and len(v) == 9
    _, v, r = TR.split_pairs(_fake_pairs(n_series=1), by="series", val_frac=0.5, seed=0)
    assert r["val_groups"] == [] and v == []
    _, _, r = TR.split_pairs(_fake_pairs(n_series=2), by="series", val_frac=1.0, seed=0)
    assert len(r["val_groups"]) == 1
    assert TR.split_pairs(pairs, by="series", seed=3)[2] == TR.split_pairs(pairs, by="series", seed=3)[2]   # seeded
    with pytest.raises(ValueError):
        TR.split_pairs(pairs, by="slice")
    assert TR.series_id_of("a_t007", 7) == "a" and TR.series_id_of("a_t007", None) == "a_t007"
    assert TR.series_id_of("a_t007", float("nan")) == "a_t007" and TR.series_id_of("a_t007", 8) == "a_t007"
    assert TR.choose_slices([0, 1, 5, 6, 7], 2, 3, np.random.default_rng(0)) and TR.choose_slices([0, 0], 2, 3, np.random.default_rng(0)) == []
    assert TR.choose_slices([1, 2, 0], 5, 3, np.random.default_rng(0)) == [0, 1]       # fewer candidates than asked: all non-empty


# ---- train-nuclei with the Cellpose seams replaced -------------------------------------------------------

@pytest.fixture
def fake_cellpose(monkeypatch):
    """Replace the three Cellpose seams: the net is a token, training writes a small 'model' file and returns
    linear losses, and prediction answers with the labels of the pair whose image it is given (perfect
    recall) unless told to answer with empty images."""
    calls: dict = {"build": [], "train": [], "predict": [], "answer": "perfect"}
    lookup: dict = {}

    def build(base, device):
        calls["build"].append((dict(base), device))
        return "net-token"

    def train_seg(net, **kw):
        calls["train"].append(kw)
        assert net == "net-token" and len(kw["train_data"]) == len(kw["train_labels"]) >= 1
        p = Path(kw["save_path"]) / "models" / kw["model_name"]
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"fake cellpose weights " * 64)
        n = int(kw["n_epochs"])
        return p, np.linspace(1.0, 0.5, n), np.linspace(0.9, 0.6, n)

    def predict(model_path, images, device):
        calls["predict"].append((str(model_path), len(images), device))
        out = []
        for img in images:
            key = hashlib.sha256(np.ascontiguousarray(img, dtype=np.float32).tobytes()).hexdigest()
            out.append(np.zeros(img.shape, np.int32) if calls["answer"] == "empty" else lookup[key])
        return out

    monkeypatch.setattr(TR, "_build_net", build)
    monkeypatch.setattr(TR, "_train_seg", train_seg)
    monkeypatch.setattr(TR, "_predict", predict)
    calls["lookup"] = lookup
    return calls


def _index_pairs(pairs_dir, lookup: dict) -> list[dict]:
    pairs = TR.load_pairs(pairs_dir)
    for p in pairs:
        img = tifffile.imread(p["image"]).astype(np.float32)
        lookup[hashlib.sha256(np.ascontiguousarray(img).tobytes()).hexdigest()] = tifffile.imread(p["labels"]).astype(np.int32)
    return pairs


def test_train_nuclei_writes_model_log_and_card(tmp_path, live_series_results, fake_cellpose):
    pairs_dir = tmp_path / "pairs"
    TR.export_training_pairs(live_series_results, pairs_dir, orientations=("xy", "xz"), slices_per_image=2, seed=1)
    pairs = _index_pairs(pairs_dir, fake_cellpose["lookup"])
    out = tmp_path / "models" / "live_nuclei_test"
    res = TR.train_nuclei(pairs_dir, out, base="cpsam", epochs=2, split_by="image", val_frac=0.5, seed=0, cpu=True,
                          batch_size=2, min_train_masks=1)
    # the model file, where the card and the registry expect it (cellpose's models/ subfolder is gone)
    model = Path(res["model_path"])
    assert model == out / "live_nuclei_test" and model.is_file() and not (out / "models").exists()
    assert res["model_sha256"] == _sha(model) and res["device"] == "cpu"
    assert fake_cellpose["build"] == [({"name": "cpsam", "kind": "cellpose_builtin", "path": None, "sha256": None}, "cpu")]
    kw = fake_cellpose["train"][0]
    assert kw["n_epochs"] == 2 and kw["learning_rate"] == 1e-5 and kw["weight_decay"] == 0.1 and kw["batch_size"] == 2
    assert kw["model_name"] == "live_nuclei_test" and kw["min_train_masks"] == 1
    # the split: one of the two timepoints is held out, with every one of its slices
    sp = res["split"]
    assert sp["by"] == "image" and len(sp["val_groups"]) == 1 and len(sp["train_groups"]) == 1
    assert sp["n_train"] == sp["n_val"] == 4 and len(kw["train_data"]) == 4 and len(kw["test_data"]) == 4
    val_ids = {p["image_id"] for p in pairs if p["image_id"] in set(sp["val_groups"])}
    assert val_ids == set(sp["val_groups"])
    # training_log.json
    logrec = json.loads(Path(res["log_path"]).read_text(encoding="utf-8"))
    assert logrec["losses"]["train"] == [1.0, 0.5] and logrec["losses"]["test"] == [0.9, 0.6]
    assert logrec["split"] == sp and logrec["counts"]["n_pairs"] == 8 and logrec["counts"]["n_val"] == 4
    assert sorted(logrec["counts"]["per_orientation"]) == ["xy", "xz"] and logrec["settings"]["epochs"] == 2
    assert set(logrec["train_files"]) | set(logrec["val_files"]) == {Path(p["image"]).name for p in pairs}
    assert not set(logrec["train_files"]) & set(logrec["val_files"])
    # the model card: base, sha256 of the weights and of every pair, orientations, epochs, date, version, metrics
    card = json.loads(Path(res["card_path"]).read_text(encoding="utf-8"))
    assert Path(res["card_path"]).name == "live_nuclei_test.card.json"
    assert card["name"] == "live_nuclei_test" and card["model_sha256"] == _sha(model) and card["base_model"]["name"] == "cpsam"
    assert card["training"]["epochs"] == 2 and card["training"]["split_by"] == "image" and card["orientations"] == ["xy", "xz"]
    assert card["date"] and card["germquant_version"] and "cellpose_version" in card and card["licence_note"]
    ts = card["training_set"]
    assert len(ts) == 8 and {t["split"] for t in ts} == {"train", "val"}
    for t in ts:
        assert t["image_sha256"] == _sha(pairs_dir / t["image"]) and t["labels_sha256"] == _sha(pairs_dir / t["labels"])
        assert t["orientation"] in ("xy", "xz") and t["image_id"].startswith(STEM)
    assert sum(t["split"] == "val" for t in ts) == 4
    hm = card["holdout_metrics"]
    assert hm["iou_threshold"] == 0.5 and hm["overall"]["f1"] == 1.0 and hm["overall"]["miss_rate"] == 0.0
    assert set(hm["per_orientation"]) == {"xy", "xz"} and all(m["f1"] == 1.0 and m["n_pairs"] == 2 for m in hm["per_orientation"].values())
    assert hm["overall"]["n_gt"] == sum(m["n_gt"] for m in hm["per_orientation"].values()) > 0
    assert len(card["holdout_per_pair"]) == 4 and all(r["f1"] == 1.0 for r in card["holdout_per_pair"])
    assert fake_cellpose["predict"] == [(str(model), 4, "cpu")]

    # an empty prediction scores zero: the metrics are honest, the model is still written
    fake_cellpose["answer"] = "empty"
    res2 = TR.train_nuclei(pairs_dir, tmp_path / "models" / "m2", base="none", epochs=1, split_by="image", cpu=True, min_train_masks=1)
    assert res2["metrics"]["overall"]["f1"] == 0.0 and res2["metrics"]["overall"]["miss_rate"] == 1.0
    assert res2["base_model"]["kind"] == "random_init" and fake_cellpose["build"][-1][0]["kind"] == "random_init"
    # a single series: nothing held out, the card says so with no metrics
    fake_cellpose["answer"] = "perfect"
    res3 = TR.train_nuclei(pairs_dir, tmp_path / "models" / "m3", base="none", epochs=1, split_by="series", cpu=True, min_train_masks=1)
    assert res3["split"]["val_groups"] == [] and res3["metrics"]["overall"] is None
    assert fake_cellpose["train"][-1]["test_data"] is None
    card3 = json.loads(Path(res3["card_path"]).read_text(encoding="utf-8"))
    assert card3["holdout_metrics"]["overall"] is None and card3["training"]["n_val"] == 0
    # every training slice too thin for Cellpose is an error before anything is trained
    with pytest.raises(ValueError, match="min_train_masks"):
        TR.train_nuclei(pairs_dir, tmp_path / "models" / "m4", base="none", epochs=1, cpu=True, min_train_masks=10_000)
    with pytest.raises(FileNotFoundError):
        TR.train_nuclei(tmp_path / "nothing_here", tmp_path / "models" / "m5", base="none", cpu=True)
    # the command line
    assert cli.main(["train-nuclei", "--pairs", str(pairs_dir), "--out", str(tmp_path / "models" / "cli_model"), "--base", "none",
                     "--epochs", "2", "--split-by", "image", "--val-frac", "0.5", "--seed", "0", "--cpu",
                     "--min-train-masks", "1"]) == 0
    assert (tmp_path / "models" / "cli_model" / "cli_model").is_file()
    assert (tmp_path / "models" / "cli_model" / "cli_model.card.json").is_file()
    assert (tmp_path / "models" / "cli_model" / "training_log.json").is_file()


# ---- register-model and the registry lookup -----------------------------------------------------------

@pytest.fixture
def models_root(tmp_path, monkeypatch):
    root = tmp_path / "models_root"
    monkeypatch.setenv("GERMQUANT_MODELS", str(root))
    return root


def _write_model_and_card(folder: Path, name: str, payload: bytes = b"weights " * 500) -> tuple[Path, Path]:
    folder.mkdir(parents=True, exist_ok=True)
    model = folder / name
    model.write_bytes(payload)
    card = folder / f"{name}.card.json"
    card.write_text(json.dumps({"name": name, "model_sha256": hashlib.sha256(payload).hexdigest(),
                                "what": "a test model", "base_model": {"name": "cpsam"}, "date": "2026-10-09"}),
                    encoding="utf-8")
    return model, card


def test_register_model_resolves_by_name(tmp_path, models_root):
    model, card = _write_model_and_card(tmp_path / "trained", "live_nuclei_test")
    assert R.local_models() == {} and R.resolve_registered("live_nuclei_test") is None
    entry = R.register_local_model(model, card)
    assert entry["name"] == "live_nuclei_test" and entry["sha256"] == _sha(model) and entry["size"] == model.stat().st_size
    assert Path(entry["path"]) == model.resolve() and Path(entry["card"]) == card.resolve() and entry["what"] == "a test model"
    reg = R.local_registry_path()
    assert reg == models_root / "models" / "local_models.json" and reg.is_file()
    assert set(R.local_models()) == {"live_nuclei_test"} and "live_nuclei_test" in R.known_models() and "germline_nuclei_combined" in R.known_models()
    assert R.resolve_registered("live_nuclei_test") == model.resolve()
    assert R.resolve_registered("cpsam") is None and R.resolve_registered("germline_nuclei_combined") is None   # not downloaded here
    assert R.get_model("live_nuclei_test", quiet=True) == model.resolve()                                        # verified, no download
    # refused: a release name, a Cellpose built-in, a path-like name, a card for other weights, a missing file
    for bad in ("germline_nuclei_combined", "cpsam", "CPSAM_v2", "a/b"):
        with pytest.raises(ValueError):
            R.register_local_model(model, card, name=bad)
    other, _ = _write_model_and_card(tmp_path / "other", "other_weights", payload=b"different")
    with pytest.raises(ValueError, match="model_sha256"):
        R.register_local_model(other, card, name="mismatch")
    with pytest.raises(FileNotFoundError):
        R.register_local_model(tmp_path / "nope.bin", None, name="x")
    # a name can be registered without a card, and registering again replaces the entry
    e2 = R.register_local_model(other, None, name="live_nuclei_test")
    assert e2["card"] is None and R.local_models()["live_nuclei_test"]["sha256"] == _sha(other)
    e3 = R.register_local_model(model, card)
    assert R.local_models()["live_nuclei_test"]["sha256"] == e3["sha256"] == _sha(model)
    # a corrupted or moved file is refused by get-model with advice, never downloaded
    model.write_bytes(b"garbage")
    with pytest.raises(RuntimeError, match="sha256"):
        R.get_model("live_nuclei_test", quiet=True)
    model.unlink()
    with pytest.raises(RuntimeError, match="missing"):
        R.get_model("live_nuclei_test", quiet=True)
    assert R.resolve_registered("live_nuclei_test") is None
    with pytest.raises(KeyError, match="live_nuclei_test"):
        R.get_model("never_registered", quiet=True)


def test_registered_name_in_the_config_reaches_the_segment_stage(tmp_path, models_root, patched_reader, registry, monkeypatch):  # noqa: F811
    model, card = _write_model_and_card(tmp_path / "trained", "live_nuclei_test")
    assert cli.main(["register-model", "--model", str(model), "--card", str(card)]) == 0
    assert cli.main(["get-model", "live_nuclei_test"]) == 0
    assert cli.main(["get-model", "no_such_model"]) == 1
    cfg = load_config(ROOT / "config" / "config.yaml")
    shipped = cfg.base_dir / str(cfg.get("segmentation.nuclei.cellpose_model"))
    if shipped.is_file():                     # a checkout that holds the release model
        assert Path(str(pipeline.resolve_model_path(cfg))) == shipped
    else:                                     # no models folder (CI): refused, never handed to Cellpose
        from germquant import model_paths

        monkeypatch.setattr(model_paths, "cellpose_installed", lambda: True)   # (without it nothing is loaded)
        with pytest.raises(FileNotFoundError, match="is not a file"):
            pipeline.resolve_model_path(cfg)
    cfg.set("segmentation.nuclei.cellpose_model", "live_nuclei_test")
    assert Path(str(pipeline.resolve_model_path(cfg))) == model.resolve()
    cfg.set("segmentation.nuclei.cellpose_model", "cpsam")
    assert pipeline.resolve_model_path(cfg) == "cpsam"                    # Cellpose's own name passes through
    # the segment stage hands the registered file to the backend, and the manifest records its digest
    seen: dict = {}

    def spy(dna, spacing, *, model, diameter_um, device, cfg_section):
        from germquant.segment import nuclei as N

        seen["model"] = model
        return N._classical(dna, spacing, diameter_um)

    B.register_backend("spy", spy)
    cfg.set("segmentation.nuclei.cellpose_model", "live_nuclei_test")
    cfg._data["segmentation"]["nuclei"]["method"] = "auto"
    cfg._data["segmentation"]["nuclei"]["backend"] = "spy"
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / "run", prov=None)
    assert Path(seen["model"]) == model.resolve() and res["stages"]["segment"]["status"] == "ran"
    man = json.loads((tmp_path / "run" / "run_manifest.json").read_text(encoding="utf-8"))
    assert man["model_sha256"] == _sha(model) == provenance.file_sha256(model)
    assert res["tables"]["image_summary"]["segmentation_method"].iloc[0] == "spy"


# ---- the real thing, opt-in ---------------------------------------------------------------------------

@pytest.mark.skipif(not os.environ.get("GERMQUANT_TRAIN_SMOKE"),
                    reason="set GERMQUANT_TRAIN_SMOKE=1 to train a randomly initialised Cellpose-SAM network for "
                           "two epochs on the CPU (about 90 s, writes a 1.2 GB model file)")
def test_real_random_init_training_on_the_cpu(tmp_path, live_series_results):
    pytest.importorskip("cellpose")
    pairs_dir = tmp_path / "pairs"
    TR.export_training_pairs(live_series_results, pairs_dir, orientations=("xy",), slices_per_image=2, seed=1)
    out = tmp_path / "models" / "smoke"
    res = TR.train_nuclei(pairs_dir, out, base="none", epochs=2, split_by="image", val_frac=0.5, seed=0, cpu=True,
                          batch_size=2, min_train_masks=1)
    model = Path(res["model_path"])
    assert model.is_file() and model.stat().st_size > 1_000_000_000          # the whole network was saved
    log = json.loads(Path(res["log_path"]).read_text(encoding="utf-8"))
    assert len(log["losses"]["train"]) == 2 and all(np.isfinite(log["losses"]["train"]))
    card = json.loads(Path(res["card_path"]).read_text(encoding="utf-8"))
    assert card["base_model"]["kind"] == "random_init" and card["holdout_metrics"]["overall"] is not None
    assert "error" not in card["holdout_metrics"]
