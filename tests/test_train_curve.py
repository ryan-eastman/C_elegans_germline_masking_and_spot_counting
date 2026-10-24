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
