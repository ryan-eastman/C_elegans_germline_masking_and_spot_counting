"""The model card fields and the learning curve (docs/TRAINING.md, steps 5 and 6; design sections 3.6 and 10;
P5 acceptance 10, the card and curve half).

* train-nuclei's card adds kind, dims, pixel_um, object_kind, system, pack {id, version, recipe}, segmented_role,
  source_url and, with --dataset, dataset {id, licence, split_sha256}. The CI card test trains with --base
  tinynet@f5055f8a, the test pack's 64-byte stand-in, resolved through the blob store as
  tests/test_packs_models.py::test_train_base_resolves_a_pack_model does; the pairs come from two test pack runs.
  The same seamed command with --base cpsam_v2@0f1cc3f7 and the workstation's model store is opt-in
  (GERMQUANT_TRAIN_CARD_BASE), because CI cannot hold the 1.2 GB blob.
* train-curve holds one validation set out, trains nested, seeded subsets of the remaining groups (sizes by
  group, repeats with their own seeded order), scores each model on the held-out pairs right after training
  (the hold-out metrics, or a scorer= seam) and keeps cards, metrics and digests; --keep-models none keeps no
  weights, best keeps one model.

The Cellpose seams (_build_net, _train_seg, _predict) are replaced throughout, as tests/test_training_loop.py does.
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
from packs_helpers import (  # noqa: F401  (pack_env is a fixture)
    TINYNET,
    TINYNET_BYTES,
    TINYNET_SHA,
    pack_env,
    patch_nd2_reader,
    small_stack,
    use_fixture_packs,
)
from test_training_loop import _index_pairs, fake_cellpose  # noqa: F401  (fake_cellpose is a fixture)

from germquant import api, cli, packs, pipeline
from germquant.segment import training as TR

NAMES = ("20261011_field_001.nd2", "20261011_field_002.nd2")


def _sha(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _card(path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _pack_pairs(tmp_path, monkeypatch) -> Path:
    """Two test pack runs (recipe basic: classical nuclei, LoG puncta) through the .nd2 seam, exported."""
    patch_nd2_reader(monkeypatch, small_stack())
    cfg = api.resolve_config(pack="test_pack", recipe="basic")
    for name in NAMES:
        pipeline.process_image(name, cfg, tmp_path / "results" / Path(name).stem, prov=None)
    pairs = tmp_path / "pairs"
    TR.export_training_pairs(tmp_path / "results", pairs, slices_per_image=2, seed=0)
    return pairs


# ---- the card ---------------------------------------------------------------------------------------------

def test_card_fields(tmp_path, pack_env, monkeypatch, fake_cellpose):  # noqa: F811
    use_fixture_packs(monkeypatch)
    good = tmp_path / "t.bin"
    good.write_bytes(TINYNET_BYTES)
    packs.fetch_pack_models("test_pack", "with_model", {TINYNET: good})
    pairs_dir = _pack_pairs(tmp_path, monkeypatch)
    _index_pairs(pairs_dir, fake_cellpose["lookup"])
    out = tmp_path / "models" / "tp_nuclei"
    assert cli.main(["train-nuclei", "--pairs", str(pairs_dir), "--out", str(out), "--base", TINYNET,
                     "--epochs", "1", "--split-by", "image", "--val-frac", "0.5", "--cpu",
                     "--min-train-masks", "1"]) == 0
    card = _card(out / "tp_nuclei.card.json")
    blob = pack_env.blobs / TINYNET_SHA / "tinynet"
    assert card["base_model"] == {"name": TINYNET, "kind": "pack", "path": str(blob), "sha256": TINYNET_SHA}
    version = packs.find_pack("test_pack").version
    assert card["kind"] == "cellpose" and card["dims"] == "3d" and card["pixel_um"] == 0.2
    assert card["object_kind"] == "nucleus" and card["system"] is None and card["source_url"] is None
    assert card["pack"] == {"id": "test_pack", "version": version, "recipe": "basic"}
    assert card["segmented_role"] == "dna" and "dataset" not in card and "pack_note" not in card
    assert card["what"].startswith("Cellpose nucleus model trained with germquant train-nuclei from " + TINYNET)
    assert card["orientations"] == ["xy", "xz", "yz"] and card["holdout_metrics"]["overall"]["f1"] == 1.0
    keys = list(card)
    assert keys[:2] == ["name", "what"] and keys.index("source_url") < keys.index("model_file")
    # the options: --kind fills object_kind, --system the system, --dataset the dataset record of
    # benchmarks/ext/<id>/ (licence from dataset.json, the digest of split.json on LF-normalised bytes)
    ext = tmp_path / "root" / "benchmarks" / "ext" / "tiny_ds"
    ext.mkdir(parents=True)
    (ext / "dataset.json").write_text(json.dumps({"id": "tiny_ds", "licence": "CC0"}), encoding="ascii")
    split_lf = b'{"format": 1, "seed": 0, "tune": [], "val": ["a"], "test": ["b"]}\n'
    (ext / "split.json").write_bytes(split_lf.replace(b"\n", b"\r\n"))
    monkeypatch.setenv("GERMQUANT_CONFIG_ROOT", str(tmp_path / "root"))
    out2 = tmp_path / "models" / "yeastlike"
    assert cli.main(["train-nuclei", "--pairs", str(pairs_dir), "--out", str(out2), "--base", TINYNET,
                     "--epochs", "1", "--split-by", "image", "--cpu", "--min-train-masks", "1",
                     "--kind", "cell", "--system", "budding yeast", "--dataset", "tiny_ds"]) == 0
    card2 = _card(out2 / "yeastlike.card.json")
    assert card2["object_kind"] == "cell" and card2["system"] == "budding yeast"
    assert card2["dataset"] == {"id": "tiny_ds", "licence": "CC0",
                                "split_sha256": hashlib.sha256(split_lf).hexdigest()}
    assert card2["what"].startswith("Cellpose cell model trained")
    res = TR.train_nuclei(pairs_dir, tmp_path / "models" / "nowhere", base="none", epochs=1, cpu=True,
                          min_train_masks=1, dataset="nowhere_ds")
    assert _card(res["card_path"])["dataset"] == {"id": "nowhere_ds", "licence": None, "split_sha256": None}
    for bad in ("../x", "a/b", "", "x" * 80):
        with pytest.raises(ValueError, match="dataset"):
            TR.train_nuclei(pairs_dir, tmp_path / "models" / "bad", base="none", cpu=True, dataset=bad)
    with pytest.raises(ValueError, match="kind"):
        TR.train_nuclei(pairs_dir, tmp_path / "models" / "bad", base="none", cpu=True, kind="organelle")
    with pytest.raises(SystemExit):
        cli.main(["train-nuclei", "--pairs", str(pairs_dir), "--out", str(tmp_path / "x"), "--kind", "organelle"])
    assert not (tmp_path / "models" / "bad").exists()           # refused before anything was trained or written
    # the source runs are unknown when the pairs' export.json names folders that are gone
    moved = tmp_path / "moved_pairs"
    moved.mkdir()
    for sub in ("image", "labels"):
        (moved / sub).mkdir()
        for f in (pairs_dir / sub).iterdir():
            (moved / sub / f.name).write_bytes(f.read_bytes())
    (moved / "manifest.csv").write_bytes((pairs_dir / "manifest.csv").read_bytes())
    exp = json.loads((pairs_dir / "export.json").read_text(encoding="utf-8"))
    for r in exp["runs"]:
        r["folder"] = str(tmp_path / "gone" / r["image_id"])
    (moved / "export.json").write_text(json.dumps(exp), encoding="utf-8")
    res = TR.train_nuclei(moved, tmp_path / "models" / "unknown", base="none", epochs=1, cpu=True, min_train_masks=1)
    c = _card(res["card_path"])
    assert c["pack"] is None and "unknown" in c["pack_note"] and c["object_kind"] == "nucleus"
    assert c["segmented_role"] == "dna" and c["pixel_um"] == 0.2     # the manifest still says these


def _square_labels(shape, n, size=5) -> np.ndarray:
    lab = np.zeros(shape, np.int32)
    for k in range(n):
        y0, x0 = 2 + (k // 2) * (size + 3), 2 + (k % 2) * (size + 3)
        lab[y0:y0 + size, x0:x0 + size] = k + 1
    return lab


def _write_pairs(folder: Path, groups, *, per_group=2, orientations=("xy",), pixel_um=0.25, role="dna",
                 n_objects=4, lookup=None, export=False) -> Path:
    """A hand-made pairs folder: per group and orientation, per_group pairs of a 24 x 24 image whose labels are
    n_objects squares; manifest.csv with every column; no export.json unless asked."""
    for sub in ("image", "labels"):
        (folder / sub).mkdir(parents=True, exist_ok=True)
    rows = []
    rng = np.random.default_rng(1)
    for g in groups:
        for orient in orientations:
            for i in range(per_group):
                lab = _square_labels((24, 24), n_objects)
                img = (100 + 900 * (lab > 0) + rng.integers(0, 50, lab.shape)).astype(np.uint16)
                stem = f"{g}_{orient}_{i:04d}.tif"
                tifffile.imwrite(folder / "image" / stem, img)
                tifffile.imwrite(folder / "labels" / stem, lab)
                if lookup is not None:
                    lookup[hashlib.sha256(np.ascontiguousarray(img.astype(np.float32)).tobytes()).hexdigest()] = lab
                px = pixel_um(g) if callable(pixel_um) else pixel_um
                rows.append({"image_id": g, "series_id": g, "t_index": None, "orientation": orient, "slice_index": i,
                             "n_labels": n_objects, "zone": "", "image": f"image/{stem}", "labels": f"labels/{stem}",
                             "source_file": f"{g}.nd2", "xy_stride": 1, "z_range": "",
                             "channel_role": role(g) if callable(role) else role, "channel_index": 0,
                             "pixel_um": px, "z_zoom": 1.0, "height": 24, "width": 24})
    pd.DataFrame(rows, columns=TR.MANIFEST_COLS).to_csv(folder / "manifest.csv", index=False)
    return folder


def test_card_dims_pixel_size_and_role_from_the_pairs(tmp_path, fake_cellpose):  # noqa: F811
    lookup = fake_cellpose["lookup"]
    proj = _write_pairs(tmp_path / "proj", ["a", "b"], orientations=("proj",), lookup=lookup)
    res = TR.train_nuclei(proj, tmp_path / "m_proj", base="none", epochs=1, cpu=True, split_by="image",
                          min_train_masks=1)
    c = _card(res["card_path"])
    assert c["dims"] == "2d" and c["orientations"] == ["proj"] and c["pixel_um"] == 0.25
    assert c["pack"] is None and "unknown" in c["pack_note"] and c["segmented_role"] == "dna"
    assert set(c["holdout_metrics"]["per_orientation"]) == {"proj"}
    mixed = _write_pairs(tmp_path / "mixed", ["a", "b"], orientations=("proj", "xy"),
                         pixel_um=lambda g: 0.25 if g == "a" else 0.5, role=lambda g: "dna" if g == "a" else "nuclei",
                         lookup=lookup)
    c2 = _card(TR.train_nuclei(mixed, tmp_path / "m_mixed", base="none", epochs=1, cpu=True, min_train_masks=1)
               ["card_path"])
    assert c2["dims"] == "2d+3d" and c2["pixel_um"] is None and c2["pixel_um_values"] == [0.25, 0.5]
    assert c2["segmented_role"] is None and c2["segmented_roles"] == ["dna", "nuclei"]
    # without a manifest the pairs still train; the card says what it cannot know
    (mixed / "manifest.csv").unlink()
    c3 = _card(TR.train_nuclei(mixed, tmp_path / "m_named", base="none", epochs=1, cpu=True, min_train_masks=1)
               ["card_path"])
    assert c3["dims"] == "2d+3d" and c3["pixel_um"] is None and c3["segmented_role"] is None
    assert "pixel_um_values" not in c3 and c3["object_kind"] == "nucleus"


# ---- the learning curve ---------------------------------------------------------------------------------------

@pytest.fixture
def curve_seams(monkeypatch):
    """Cellpose replaced: a model file that names its training size, and a predictor that finds one object per
    two training pairs (all of them from 8 pairs on), so a larger subset scores higher."""
    state: dict = {"lookup": {}, "written": {}, "trained": [], "present_at_train": []}

    def build(base, device):
        return "net"

    def train_seg(net, **kw):
        out = Path(kw["save_path"])
        state["present_at_train"].append(sum(1 for p in out.parent.glob("*/*") if p.is_file() and p.suffix == ""
                                             and p.name == p.parent.name))
        p = out / "models" / kw["model_name"]
        p.parent.mkdir(parents=True, exist_ok=True)
        content = f"{kw['model_name']} n={len(kw['train_data'])}\n".encode("ascii") * 4
        p.write_bytes(content)
        state["written"][kw["model_name"]] = hashlib.sha256(content).hexdigest()
        state["trained"].append((kw["model_name"], len(kw["train_data"]), len(kw["test_data"] or [])))
        n = int(kw["n_epochs"])
        return p, np.linspace(1.0, 0.5, n), np.linspace(0.9, 0.6, n)

    def predict(model_path, images, device):
        n = int(Path(model_path).read_bytes().split(b"n=")[1].split(b"\n")[0])
        out = []
        for img in images:
            key = hashlib.sha256(np.ascontiguousarray(img, dtype=np.float32).tobytes()).hexdigest()
            lab = state["lookup"][key].copy()
            lab[lab > n // 2] = 0
            out.append(lab)
        return out

    monkeypatch.setattr(TR, "_build_net", build)
    monkeypatch.setattr(TR, "_train_seg", train_seg)
    monkeypatch.setattr(TR, "_predict", predict)
    return state


GROUPS = [f"img{k}" for k in range(6)]


def _curve(pairs, out, **kw):
    args = {"sizes": [1, 2, 3, "all", 16], "repeats": 2, "seed": 5, "base": "none", "epochs": 1, "cpu": True,
            "split_by": "image", "val_frac": 0.34, "min_train_masks": 1}
    args.update(kw)
    return TR.train_curve(pairs, out, **args)


def test_subsets_nested_and_deterministic(tmp_path, curve_seams):
    pairs = _write_pairs(tmp_path / "pairs", GROUPS, lookup=curve_seams["lookup"])
    res = _curve(pairs, tmp_path / "curve")
    curve = json.loads((tmp_path / "curve" / "curve.json").read_text(encoding="utf-8"))
    _, val, split = TR.split_pairs(TR.load_pairs(pairs), by="image", val_frac=0.34, seed=5)
    assert curve["val_groups"] == split["val_groups"] == res["val_groups"] and len(split["val_groups"]) == 2
    pool = split["train_groups"]
    assert curve["train_pool"] == pool and len(pool) == 4
    assert curve["sizes"] == [1, 2, 3, 4] and curve["skipped_sizes"] == [16] and curve["repeats"] == 2
    subsets = curve["subsets"]
    assert sorted(subsets) == ["r0", "r1"]
    for r in ("r0", "r1"):
        s = subsets[r]
        assert sorted(s) == ["1", "2", "3", "4"] and s["4"] == pool
        for a, b in (("1", "2"), ("2", "3"), ("3", "4")):
            assert set(s[a]) < set(s[b]) and len(s[a]) == int(a)          # nested: each subset holds the smaller
    # every model trained on exactly its subset's pairs and was scored on the one validation set
    rows = pd.read_csv(tmp_path / "curve" / "curve.csv", float_precision="round_trip")
    assert len(rows) == 8 and list(rows[["size", "repeat"]].itertuples(index=False, name=None)) == \
        [(k, r) for k in (1, 2, 3, 4) for r in (0, 1)]
    assert list(rows["size_label"].astype(str)) == ["1", "1", "2", "2", "3", "3", "all", "all"]
    val_files = sorted(Path(p["image"]).name for p in val)
    for row in rows.to_dict("records"):
        folder = tmp_path / "curve" / Path(row["card"]).parent
        log = json.loads((folder / "training_log.json").read_text(encoding="utf-8"))
        groups = subsets[f"r{row['repeat']}"][str(row["size"])]
        assert row["groups"] == ";".join(groups) and row["n_groups"] == len(groups)
        assert sorted({f.split("_xy_")[0] for f in log["train_files"]}) == groups
        assert sorted(log["val_files"]) == val_files and row["n_val"] == len(val_files)
        assert row["n_train"] == 2 * len(groups) and row["status"] == "scored"
    assert [t[1:] for t in curve_seams["trained"]] == [(2 * k, 4) for k in (1, 2, 3, 4) for _r in (0, 1)]
    assert rows["score"].is_monotonic_increasing and rows["score"].iloc[-1] == 1.0
    # deterministic: the same seed gives the same subsets, scores and digests; the defaults keep no weights
    _curve(pairs, tmp_path / "again")
    again = json.loads((tmp_path / "again" / "curve.json").read_text(encoding="utf-8"))
    assert again["subsets"] == subsets and again["val_groups"] == curve["val_groups"]
    rows2 = pd.read_csv(tmp_path / "again" / "curve.csv", float_precision="round_trip")
    for col in ("size", "repeat", "groups", "score", "f1"):
        assert list(rows2[col]) == list(rows[col]), col
    assert res["keep_models"] == "none" and res["scorer"] == "holdout_metrics"
    # a curve needs a validation set and a size that fits the pool
    one = _write_pairs(tmp_path / "one", ["solo"], lookup=curve_seams["lookup"])
    with pytest.raises(ValueError, match="held out"):
        _curve(one, tmp_path / "c1")
    with pytest.raises(ValueError, match="sizes"):
        _curve(pairs, tmp_path / "c2", sizes=[16, 32])
    for bad in ([0], ["most"], [1.5]):
        with pytest.raises(ValueError, match="sizes"):
            _curve(pairs, tmp_path / "c3", sizes=bad)
    with pytest.raises(ValueError, match="repeats"):
        _curve(pairs, tmp_path / "c4", repeats=0)
    with pytest.raises(ValueError, match="keep_models"):
        _curve(pairs, tmp_path / "c5", keep_models="all")
    with pytest.raises(TypeError, match="sizes_per_image"):
        _curve(pairs, tmp_path / "c6", sizes_per_image=3)


def test_keep_models_none_and_best(tmp_path, curve_seams):
    pairs = _write_pairs(tmp_path / "pairs", GROUPS, lookup=curve_seams["lookup"])
    res = _curve(pairs, tmp_path / "none", sizes=[1, "all"], keep_models="none")
    rows = pd.read_csv(tmp_path / "none" / "curve.csv", float_precision="round_trip")
    assert len(rows) == 4 and not rows["weights_kept"].any()
    for row in rows.to_dict("records"):
        card_path = tmp_path / "none" / row["card"]
        card = _card(card_path)
        name = card["name"]
        assert not (tmp_path / "none" / row["model"]).exists()          # no weights
        assert (card_path.parent / "training_log.json").is_file()        # the log and the card stay
        assert card["model_sha256"] == row["model_sha256"] == curve_seams["written"][name]
        assert card["model_size"] == row["model_size"] > 0
        assert card["holdout_metrics"]["overall"]["f1"] == row["f1"] == row["score"]
        assert card["curve"] == {"size": row["size"], "size_label": str(row["size_label"]), "repeat": row["repeat"],
                                 "groups": row["groups"].split(";"), "scorer": "holdout_metrics",
                                 "score": row["score"], "weights_kept": False, "curve": "curve.json"}
    assert res["best"]["size"] == 4 and res["best"]["repeat"] == 0 and res["best"]["weights_kept"] is False
    # best: one model stays, the highest score (ties keep the earlier model), never more than two on disk
    curve_seams["present_at_train"].clear()
    res = _curve(pairs, tmp_path / "best", sizes=[1, 2, "all"], keep_models="best")
    rows = pd.read_csv(tmp_path / "best" / "curve.csv", float_precision="round_trip")
    kept = rows[rows["weights_kept"]]
    assert len(kept) == 1 and (kept["size"].iloc[0], kept["repeat"].iloc[0]) == (4, 0)
    weights = [p for p in (tmp_path / "best" / "models").glob("*/*") if p.name == p.parent.name]
    assert [p.parent.name for p in weights] == [kept["model"].iloc[0].split("/")[1]]
    assert _sha(weights[0]) == kept["model_sha256"].iloc[0]
    for row in rows.to_dict("records"):
        assert _card(tmp_path / "best" / row["card"])["curve"]["weights_kept"] is bool(row["weights_kept"])
    assert max(curve_seams["present_at_train"]) <= 1
    assert res["best"]["model"] == kept["model"].iloc[0] and res["best"]["weights_kept"] is True
    # a scorer seam decides instead, and its name is recorded (P7's end-to-end scorer plugs in here)
    seen = []

    def smallest_first(model_path, val, *, device):
        seen.append((Path(model_path).name, len(val), device))
        n = int(Path(model_path).read_bytes().split(b"n=")[1].split(b"\n")[0])
        return {"score": -float(n), "metrics": {"n_train_pairs": n}}

    smallest_first.name = "smallest_first"
    res = _curve(pairs, tmp_path / "scored", sizes=[1, "all"], repeats=1, keep_models="best", scorer=smallest_first)
    assert res["scorer"] == "smallest_first" and [s[1:] for s in seen] == [(4, "cpu"), (4, "cpu")]
    rows = pd.read_csv(tmp_path / "scored" / "curve.csv", float_precision="round_trip")
    assert list(rows["score"]) == [-2.0, -8.0] and list(rows["weights_kept"]) == [True, False]
    card = _card(tmp_path / "scored" / rows["card"].iloc[0])
    assert card["curve"]["scorer"] == "smallest_first" and card["curve"]["metrics"] == {"n_train_pairs": 2}
    assert json.loads((tmp_path / "scored" / "curve.json").read_text(encoding="utf-8"))["scorer"] == "smallest_first"


def test_train_curve_command(tmp_path, curve_seams, capsys):
    # six objects per pair: the command trains with train-nuclei's defaults (min_train_masks 5)
    pairs = _write_pairs(tmp_path / "pairs", GROUPS, n_objects=6, lookup=curve_seams["lookup"])
    capsys.readouterr()
    assert cli.main(["train-curve", "--pairs", str(pairs), "--out", str(tmp_path / "c"), "--sizes", "1,all",
                     "--repeats", "1", "--seed", "3", "--base", "none", "--keep-models", "best", "--json"]) == 0
    rec = json.loads(capsys.readouterr().out)
    assert [m["size_label"] for m in rec["models"]] == ["1", "all"] and rec["keep_models"] == "best"
    assert rec["curve_csv"] == str(tmp_path / "c" / "curve.csv") and Path(rec["curve_json"]).is_file()
    assert sum(bool(m["weights_kept"]) for m in rec["models"]) == 1
    assert cli.main(["train-curve", "--pairs", str(pairs), "--out", str(tmp_path / "d"), "--sizes", "1",
                     "--repeats", "1", "--base", "none"]) == 0
    text = capsys.readouterr().out
    assert "curve.csv" in text and "size 1" in text
    one = _write_pairs(tmp_path / "one", ["solo"], n_objects=6, lookup=curve_seams["lookup"])
    assert cli.main(["train-curve", "--pairs", str(one), "--out", str(tmp_path / "e"), "--base", "none"]) == 2
    assert "held out" in capsys.readouterr().err
    assert cli.main(["train-curve", "--pairs", str(pairs), "--out", str(tmp_path / "f"), "--sizes", "1,x",
                     "--base", "none"]) == 2
    with pytest.raises(SystemExit):
        cli.main(["train-curve", "--pairs", str(pairs), "--out", str(tmp_path / "g"), "--keep-models", "all"])


# ---- the real base, once on the workstation ------------------------------------------------------------------

@pytest.mark.skipif(not os.environ.get("GERMQUANT_TRAIN_CARD_BASE"),
                    reason="set GERMQUANT_TRAIN_CARD_BASE=cpsam_v2@0f1cc3f7 (and GERMQUANT_MODELS to a store holding "
                           "that blob) to write the card with a real pack base; the Cellpose seams stay replaced")
def test_card_with_a_real_pack_base(tmp_path, monkeypatch, fake_cellpose):  # noqa: F811
    base = os.environ["GERMQUANT_TRAIN_CARD_BASE"]
    use_fixture_packs(monkeypatch)                     # the test pack, beside the built-in packs and the store
    pairs_dir = _pack_pairs(tmp_path, monkeypatch)
    _index_pairs(pairs_dir, fake_cellpose["lookup"])
    out = tmp_path / "models" / "real_base"
    assert cli.main(["train-nuclei", "--pairs", str(pairs_dir), "--out", str(out), "--base", base, "--epochs", "1",
                     "--split-by", "image", "--val-frac", "0.5", "--cpu", "--min-train-masks", "1"]) == 0
    card = _card(out / "real_base.card.json")
    bm = card["base_model"]
    assert bm["name"] == base and bm["kind"] == "pack" and bm["sha256"].startswith(base.split("@")[1])
    assert Path(bm["path"]).is_file() and TR._sha256(bm["path"]) == bm["sha256"]
    assert card["dims"] == "3d" and card["pixel_um"] == 0.2 and card["object_kind"] == "nucleus"
    assert card["pack"] == {"id": "test_pack", "version": packs.find_pack("test_pack").version, "recipe": "basic"}
    print(json.dumps({k: card[k] for k in ("base_model", "kind", "dims", "pixel_um", "object_kind", "system", "pack",
                                           "segmented_role", "source_url", "model_sha256")}, indent=1))


# ---- the card's pack and object kind from the run records; the command's card options -------------------------

def test_card_pack_and_object_kind_from_the_run_records(tmp_path, pack_env, monkeypatch, fake_cellpose):  # noqa: F811
    use_fixture_packs(monkeypatch)
    pairs_dir = _pack_pairs(tmp_path, monkeypatch)
    _index_pairs(pairs_dir, fake_cellpose["lookup"])
    runs = json.loads((pairs_dir / "export.json").read_text(encoding="utf-8"))["runs"]

    def records(edit):
        # every record a source run's pack is read from: run_manifest.json and <id>__stages.json
        for i, r in enumerate(runs):
            for name in ("run_manifest.json", f"{r['image_id']}__stages.json"):
                f = Path(r["folder"]) / name
                data = json.loads(f.read_text(encoding="utf-8"))
                data["pack"] = edit(i, data["pack"])
                f.write_text(json.dumps(data), encoding="utf-8")

    def card(name):
        res = TR.train_nuclei(pairs_dir, tmp_path / "m" / name, base="none", epochs=1, cpu=True, split_by="image",
                              min_train_masks=1)
        return _card(res["card_path"])

    version = packs.find_pack("test_pack").version
    # a project's edited recipe (project_config records edited, section 9.1) is not the pack recipe itself
    records(lambda i, p: {**p, "edited": True})
    assert card("edited")["pack"] == {"id": "test_pack", "version": version, "recipe": "basic", "edited": True}
    records(lambda i, p: {**p, "edited": i == 0})
    c = card("half_edited")
    assert c["pack"] is None and "different" in c["pack_note"] and "basic (edited)" in c["pack_note"]
    # an unedited project run is the pack recipe itself: the card's pack is the same as for --pack runs
    records(lambda i, p: {**p, "edited": False})
    assert card("unedited")["pack"] == {"id": "test_pack", "version": version, "recipe": "basic"}
    # a cell pack: object_kind comes from the pack's record, also when the runs used two of its recipes
    cell = {"id": "cellpack", "version": "0.1.0", "object_kind": "cell", "edited": False}
    records(lambda i, p: {**p, **cell, "recipe": "cells"})
    c = card("one_recipe")
    assert c["object_kind"] == "cell" and c["what"].startswith("Cellpose cell model")
    assert c["pack"] == {"id": "cellpack", "version": "0.1.0", "recipe": "cells"}
    records(lambda i, p: {**p, **cell, "recipe": "cells" if i == 0 else "cells_2d"})
    c = card("two_recipes")
    assert c["object_kind"] == "cell" and c["pack"] is None
    assert c["pack_note"].endswith("cellpack@0.1.0/cells; cellpack@0.1.0/cells_2d")
    # runs of a cell pack and a nucleus pack: nucleus unless --kind says otherwise
    records(lambda i, p: {**p, **cell, "recipe": "cells"} if i == 0 else {**p, "id": "test_pack", "recipe": "basic",
                                                                         "object_kind": "nucleus"})
    assert card("mixed")["object_kind"] == "nucleus"
