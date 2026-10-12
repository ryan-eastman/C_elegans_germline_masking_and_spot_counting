"""Cards, lint and the base-weight comparison (design 3.6, 3.7, 3.11; P1 acceptance 11).

validate_card accepts every card train_nuclei writes today (model_size; size only as an alias in a
pack-written card) and marks model_file and base_model.path as non-portable; lint refuses a recipe key that
is neither set by config.yaml nor on the allow-list, a deny-listed key such as spots.channel, a tunable the
recipe does not set or that is not a scalar, and the other outside-recipe faults; the C. elegans card records
the base-weight comparison, its method and its result (cpsam_v2, from the run on the workstation CPU); the
comparison itself is checked on small fake state dicts.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import shutil
from pathlib import Path

import numpy as np
import pytest
import tifffile
import yaml
from packs_helpers import CELEGANS_PACK, ROOT, TEST_PACK, copy_pack, pack_env  # noqa: F401

from germquant import models_registry as MR
from germquant import packs
from germquant.config import load_config
from germquant.packs import basecheck as BC
from germquant.packs import lint as L


# ---- validate_card ------------------------------------------------------------------------------------------

def _train_a_card(tmp_path, monkeypatch) -> dict:
    """A card exactly as germquant train-nuclei writes it, with the three Cellpose seams replaced."""
    from germquant.segment import training as TR

    pairs = tmp_path / "pairs"
    for d in ("image", "labels"):
        (pairs / d).mkdir(parents=True)
    rng = np.random.default_rng(0)
    for name in ("imgA__xy_0001", "imgB__xy_0002"):
        lab = np.zeros((32, 32), np.uint16)
        lab[4:12, 4:12], lab[18:28, 16:26] = 1, 2
        tifffile.imwrite(pairs / "image" / f"{name}.tif", (rng.random((32, 32)) * 100).astype(np.float32))
        tifffile.imwrite(pairs / "labels" / f"{name}.tif", lab)

    def train_seg(net, **kw):
        p = Path(kw["save_path"]) / "models" / kw["model_name"]
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"fake weights " * 32)
        return p, np.linspace(1, 0.5, 2), np.linspace(0.9, 0.6, 2)

    monkeypatch.setattr(TR, "_build_net", lambda base, device: "net")
    monkeypatch.setattr(TR, "_train_seg", train_seg)
    monkeypatch.setattr(TR, "_predict", lambda path, images, device: [np.zeros(i.shape, np.int32) for i in images])
    res = TR.train_nuclei(pairs, tmp_path / "out" / "my_model", base="cpsam", epochs=2, split_by="image",
                          val_frac=0.5, cpu=True, min_train_masks=1)
    return json.loads(Path(res["card_path"]).read_text(encoding="utf-8"))


def test_validate_card_accepts_a_real_train_nuclei_card(tmp_path, monkeypatch):
    card = _train_a_card(tmp_path, monkeypatch)
    res = MR.validate_card(card)
    assert res["ok"], res["problems"]
    assert res["non_portable"] == ["model_file"]                       # base cpsam: no base_model.path
    spec = {"name": f"my_model@{card['model_sha256'][:8]}", "sha256": card["model_sha256"],
            "size": card["model_size"]}
    assert MR.validate_card(card, spec)["ok"]
    assert not MR.validate_card(card, {**spec, "sha256": "0" * 64})["ok"]
    assert not MR.validate_card(card, {**spec, "size": 1})["ok"]
    card2 = dict(card, base_model={"name": "x", "path": str(tmp_path / "base_weights")})
    assert MR.validate_card(card2)["non_portable"] == ["model_file", "base_model.path"]


def test_validate_card_reads_model_size_and_size_only_as_a_pack_alias():
    card = {"name": "m", "model_sha256": "a" * 64, "size": 10}
    assert not MR.validate_card(card)["ok"]                            # a train_nuclei card carries model_size
    assert MR.validate_card(card, pack_written=True)["ok"]
    assert MR.validate_card({"name": "m", "model_sha256": "a" * 64, "model_size": 10})["ok"]
    assert MR.validate_card({"model_sha256": "x"})["problems"]
    assert MR.validate_card("not a dict")["ok"] is False


def test_shipped_cards_validate_against_their_specs():
    for pack_dir in (TEST_PACK, CELEGANS_PACK):
        pack = packs.load_pack(pack_dir)
        for spec in pack.models:
            card = json.loads((pack_dir / spec.card).read_text(encoding="utf-8"))
            res = MR.validate_card(card, spec.as_dict(), pack_written=True)
            assert res["ok"] and res["non_portable"] == [], (spec.name, res)


# ---- lint -----------------------------------------------------------------------------------------------------

def test_shipped_packs_lint_clean():
    for pack_dir in (TEST_PACK, CELEGANS_PACK):
        rep = packs.lint_pack(pack_dir)
        assert rep["problems"] == [], rep


def _edited(tmp_path, rel, transform, *, fix_hashes=True) -> Path:
    pack = copy_pack(TEST_PACK, tmp_path / "test_pack", rebuild=False)
    p = pack / rel
    p.write_bytes(transform(p.read_bytes().decode("ascii")).encode("utf-8"))
    if fix_hashes:                       # the expected hashes and the recipe goldens' sentinels follow the edit
        data = yaml.safe_load((pack / "pack.yaml").read_text(encoding="ascii"))
        for r in data["recipes"]:
            r["expected_config_hash"] = load_config(pack / r["profile"]).hash
            sentinel = pack / "example" / "expected" / "recipe" / r["id"] / "config_hash.txt"
            if sentinel.is_file():
                sentinel.write_bytes((r["expected_config_hash"] + "\n").encode("ascii"))
        (pack / "pack.yaml").write_bytes(yaml.safe_dump(data, sort_keys=False).encode("ascii"))
    packs.write_sums(pack)
    return pack


def _problems(pack) -> str:
    return " | ".join(packs.lint_pack(pack)["problems"])


def test_lint_refuses_a_misspelt_key_and_a_deny_listed_key(tmp_path):
    typo = _edited(tmp_path / "a", "config/profiles/basic.yaml", lambda t: t + "  efect_size_min: 3.0\n")
    assert "key spots.efect_size_min: the engine does not read it" in _problems(typo)
    deny = _edited(tmp_path / "b", "config/profiles/pack_base.yaml", lambda t: t.replace("spots:\n", "spots:\n  channel: foci\n"))
    assert "key spots.channel: nothing reads it" in _problems(deny)
    t5 = _edited(tmp_path / "c", "config/profiles/basic.yaml", lambda t: t + "  sharpen: true\n")
    assert "sharpen" not in _problems(t5)                               # a key tranche 5 adds is allowed
    assert L.key_readable("spots.log.min_snr", set()) and L.key_readable("segmentation.nuclei.backend", set())
    assert not L.key_readable("spots.channel", {"spots.channel"})


@pytest.mark.parametrize("rel, transform, needle", [
    ("config/profiles/pack_base.yaml", lambda t: t.replace("germline:\n  enabled: false\n", "germline:\n"),
     "germline.enabled is not set explicitly"),
    ("config/profiles/pack_base.yaml", lambda t: t.replace("filename_regex: (?!x)x", "filename_regex: '(?P<x>.+)'"),
     "metadata.filename_regex must be"),
    ("config/profiles/basic.yaml", lambda t: t.replace("method: classical", "method: auto"),
     "never auto"),
    ("config/profiles/pack_base.yaml", lambda t: t.replace("germline_outlier_factor: 1000000000.0",
                                                           "germline_outlier_factor: 1.8"),
     "switch the batch framing check off"),
    ("config/profiles/pack_base.yaml", lambda t: t.replace("  min_nuclei: 1\n", ""), "qc.min_nuclei must be set"),
    ("config/profiles/basic.yaml", lambda t: t + "  instances:\n  - {name: foci2, role: foci}\n",
     "expected_from_germ_cell: false"),
    ("config/profiles/basic.yaml", lambda t: t.replace("extends: pack_base.yaml", "extends: generic.yaml"),
     "never config/config.yaml"),
    ("config/profiles/basic.yaml", lambda t: t + "# caf\u00e9\n", "not ASCII"),
    ("README.md", lambda t: t + "\u2014 dash\n", "an en or em dash"),
    ("config/profiles/basic.yaml", lambda t: t + "qc:\n  exclusions_file: C:/Users/x/excl.json\n", "machine path"),
    ("pack.yaml", lambda t: t.replace("{key: spots.log.min_snr", "{key: spots.no_such_key", 1),
     "recipe basic: tunable spots.no_such_key is not set by the recipe"),
    ("pack.yaml", lambda t: t.replace("{key: spots.log.min_snr", "{key: spots.log", 1),
     "recipe basic: tunable spots.log is not a scalar"),
])
def test_lint_outside_recipe_rules(tmp_path, rel, transform, needle):
    pack = _edited(tmp_path, rel, transform)
    assert needle in _problems(pack)


def test_lint_hashes_sums_versioned_pack_lines_and_config_yaml(tmp_path):
    stale = _edited(tmp_path / "a", "config/profiles/basic.yaml", lambda t: t + "  extra_line: 1\n", fix_hashes=False)
    assert "expected_config_hash" in _problems(stale)
    unversioned = _edited(tmp_path / "b", "config/profiles/with_model.yaml",
                          lambda t: t.replace("extends: pack_base.yaml", "extends: pack:test_pack/basic"),
                          fix_hashes=False)
    assert "must name a pack version" in _problems(unversioned)
    shadow = copy_pack(TEST_PACK, tmp_path / "c" / "test_pack")
    shutil.copyfile(ROOT / "config" / "config.yaml", shadow / "config" / "config.yaml")
    assert "would shadow" in _problems(shadow)
    (shadow / "config" / "config.yaml").unlink()
    (shadow / "README.md").write_bytes(b"edited\n")
    assert "SHA256SUMS" in _problems(shadow)


def test_lint_previous_enforces_the_version_bump(tmp_path):
    old = copy_pack(TEST_PACK, tmp_path / "old" / "test_pack")
    doc = copy_pack(TEST_PACK, tmp_path / "doc" / "test_pack",
                    edit=lambda d: d.update(version="0.1.1"))
    (doc / "CHANGELOG.md").write_bytes(b"# changed\n")
    packs.write_sums(doc)
    assert packs.lint_pack(doc, previous=old)["problems"] == []          # a doc change, a patch bump
    same = copy_pack(TEST_PACK, tmp_path / "same" / "test_pack")
    (same / "CHANGELOG.md").write_bytes(b"# changed\n")
    packs.write_sums(same)
    assert "must be above the previous" in _problems_prev(same, old)
    hashed = _edited(tmp_path / "hashed", "config/profiles/basic.yaml", lambda t: t.replace("min_snr: 2.0", "min_snr: 2.5"))
    data = yaml.safe_load((hashed / "pack.yaml").read_text(encoding="ascii"))
    data["version"] = "0.1.1"
    (hashed / "pack.yaml").write_bytes(yaml.safe_dump(data, sort_keys=False).encode("ascii"))
    packs.write_sums(hashed)
    assert "need a minor version bump" in _problems_prev(hashed, old)  # before 1.0.0: hashed text -> minor
    data["version"] = "0.2.0"
    (hashed / "pack.yaml").write_bytes(yaml.safe_dump(data, sort_keys=False).encode("ascii"))
    packs.write_sums(hashed)
    assert _problems_prev(hashed, old) == ""
    old1 = copy_pack(TEST_PACK, tmp_path / "old1" / "test_pack", edit=lambda d: d.update(version="1.0.0"))
    data["version"] = "1.1.0"
    (hashed / "pack.yaml").write_bytes(yaml.safe_dump(data, sort_keys=False).encode("ascii"))
    packs.write_sums(hashed)
    assert "need a major version bump" in _problems_prev(hashed, old1)  # from 1.0.0: hashed text -> major


def _problems_prev(pack, previous) -> str:
    return " | ".join(packs.lint_pack(pack, previous=previous)["problems"])


# ---- the base-weight comparison -------------------------------------------------------------------------------

def _fake_bases(seed=0):
    rng = np.random.default_rng(seed)
    v2 = {f"blocks.{i}.w": rng.normal(size=(8, 8)).astype(np.float32) for i in range(10)}
    v1 = {k: rng.normal(size=v.shape).astype(np.float32) for k, v in v2.items()}
    return v1, v2


def test_comparison_names_the_base_a_fine_tune_started_from():
    v1, v2 = _fake_bases()
    rng = np.random.default_rng(1)
    tuned = {k: v + 0.01 * rng.normal(size=v.shape).astype(np.float32) for k, v in v2.items()}
    tuned["blocks.0.w"] = v2["blocks.0.w"].copy()                         # a frozen layer stays identical
    rec = BC.compare_state_dicts(tuned, {"cpsam_v2": v2, "cpsam": v1})
    assert rec["result"] == "cpsam_v2" and rec["separation_ratio"] > 10
    assert rec["bases"]["cpsam_v2"]["n_identical"] == 1 and rec["bases"]["cpsam_v2"]["coverage"] == 1.0
    assert "tensor-by-tensor" in rec["method"]


def test_comparison_says_not_determined_when_the_bases_do_not_separate():
    v1, v2 = _fake_bases()
    mid = {k: (v1[k] + v2[k]) / 2 for k in v1}
    assert BC.compare_state_dicts(mid, {"cpsam_v2": v2, "cpsam": v1})["result"] == "not determined"
    partial = {"other.w": np.zeros((2, 2), np.float32)}
    assert BC.compare_state_dicts(partial, {"cpsam_v2": v2, "cpsam": v1})["result"] == "not determined"


def test_compare_files_and_the_card_writer_on_tiny_torch_files(tmp_path):
    torch = pytest.importorskip("torch")
    v1, v2 = _fake_bases(3)
    tuned = {k: torch.from_numpy(v + 0.001) for k, v in v2.items()}
    paths = {}
    for name, sd in (("model", tuned), ("cpsam_v2", {k: torch.from_numpy(v) for k, v in v2.items()}),
                     ("cpsam", {k: torch.from_numpy(v) for k, v in v1.items()})):
        paths[name] = tmp_path / name
        torch.save(sd, paths[name])
    rec = BC.compare_files(paths["model"], {"cpsam_v2": paths["cpsam_v2"], "cpsam": paths["cpsam"]})
    assert rec["result"] == "cpsam_v2"
    spec = importlib.util.spec_from_file_location("compare_base_weights", ROOT / "scripts" / "compare_base_weights.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    card = tmp_path / "card.json"
    shutil.copyfile(CELEGANS_PACK / "models" / "germline_nuclei_combined.card.json", card)
    out = mod.update_card(card, rec, {"cpsam_v2": paths["cpsam_v2"], "cpsam": paths["cpsam"]})
    assert out["base_model"]["name"] == "cpsam_v2"
    assert out["base_model"]["sha256"] == hashlib.sha256(paths["cpsam_v2"].read_bytes()).hexdigest()
    assert out["base_model"]["comparison"]["result"] == "cpsam_v2"
    card.read_bytes().decode("ascii")                                     # still ASCII


def test_celegans_card_records_the_comparison_method_and_result():
    card = json.loads((CELEGANS_PACK / "models" / "germline_nuclei_combined.card.json").read_text(encoding="utf-8"))
    base = card["base_model"]
    comp = base["comparison"]
    assert "tensor-by-tensor" in comp["method"] and comp["result"] in ("not determined", "cpsam_v2", "cpsam")
    assert base["name"] == comp["result"]                     # the card names exactly what the comparison found
    if comp["result"] != "not determined":
        want = base["candidates"][comp["result"]].get("sha256")
        assert re.fullmatch(r"[0-9a-f]{64}", base["sha256"]) and want in (None, base["sha256"])
        assert comp["separation_ratio"] >= comp["separation_required"]
        assert comp["bases"][comp["result"]]["coverage"] >= comp["min_coverage"]
    # the run on the workstation (scripts/compare_base_weights.py --write, CPU, 2026-10-11) named cpsam_v2,
    # about 39 times closer than cpsam (v1)
    assert comp["result"] == "cpsam_v2"
    assert base["sha256"] == "0f1cc3f7ecdd8a037a57c6c48d9d8921391be4cbce3fa9f13c3e3a2e1253c667"
    assert "compare_base_weights.py --write" in comp["how"]
    assert "CC BY-NC" in card["licence_note"] and card["model_sha256"] == MR.MODELS["germline_nuclei_combined"]["sha256"]
