"""`germquant pack derive`, `--project` on the pack entry, and the version features (Train B unit P5-DERIVE;
docs/DESIGN_system_packs.md 3.5 to 3.8 and section 10, step 8; P5 acceptance 10, the pack derive half; plan
decision D26 and the plugin's interface items I1 and I4).

* `pack derive --from test_pack@0.1.0 --recipe with_model --model M --card C --out <project>/packs` writes the
  derived pack test_pack_derived_<sha8>@0.1.0 beside a copy of its parent: one recipe that extends the parent by
  versioned name and names the model <name>@<sha8>, a model spec with sources [], the rewritten card, a README and
  SHA256SUMS. It verifies, `pack list --project` lists it, `pack check --project` finds it usable, and a run of
  its recipe resolves the model from the store through the pin check and writes the parent recipe's tables.
* `pack lint` reports no ASCII, dash, machine-path or unversioned pack: problem on it. Its outside-recipe problems
  are the expected ones listed here: no recipe golden, and the worm switches and the segmentation method that the
  parent's files set (packs/lint.py counts only the keys set inside the pack folder, and a derived recipe sets
  only its extends line and the model).
* The card holds relative names plus sha256 and ASCII text only; the model sits in the blob store with its
  record; on another computer `pack fetch-models --from-file NAME=PATH` installs it from the training folder.
* Deriving again is a no-op. Refused before anything is written: a recipe that never runs a model, a card of
  other weights, a bad --id, a different pack under an id and version already there, an --out inside the
  engine's own packs folder.
* `python -m germquant.packs --version --json` lists "features", each exactly when its code is present.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import yaml
from packs_helpers import (  # noqa: F401  (pack_env is a fixture)
    ROOT,
    TEST_PACK,
    TINYNET,
    TINYNET_BYTES,
    pack_env,
    patch_nd2_reader,
    small_stack,
    use_fixture_packs,
)

from germquant import cli, cli_packs, packs
from germquant import models_registry as MR
from germquant.config import load_config
from germquant.segment import backends as SB

MODEL_BYTES = b"germquant test: a stand-in for a model trained in a project, not weights\n"
MODEL_SHA = hashlib.sha256(MODEL_BYTES).hexdigest()
DERIVED = f"test_pack_derived_{MODEL_SHA[:8]}"
MODEL_NAME = f"r01_nuclei@{MODEL_SHA[:8]}"
CARD_REL = "models/r01_nuclei.card.json"
NAME = "20261023_field_001.nd2"
IID = "20261023_field_001"
DASH, MICRO = chr(0x2013), chr(0x00B5)        # non-ASCII card text, built with chr() (the G7 style check)
# the features of the Train B contract (plan section 9.3), in its order
CONTRACT = ("packs", "planar", "labels_from", "projects", "sample_sheet", "import_labels", "train_curve",
            "pack_derive", "tune_spots", "post_stages")


LICENCE = f"Cellpose-SAM weights are BSD-3-Clause; Cellpose's training data is CC BY-NC (https://www.cellpose.org/ {MICRO})"


def _card(model: Path, base: Path, sha: str = MODEL_SHA, size: int = len(MODEL_BYTES)) -> dict:
    """A card as train-nuclei writes it: absolute machine paths in model_file and base_model, the base path inside
    the "what" text (as train-nuclei quotes --base) and another path quoted there that no key holds, a URL, and
    non-ASCII text."""
    pairs = model.parents[2] / "pairs"
    return {"name": "r01_nuclei",
            "what": f"Cellpose nucleus model trained with germquant train-nuclei from {base} on 3 slice pair(s) "
                    f"(xy) {DASH} pixel size 0.2 {MICRO}m; pairs exported to {pairs}.",
            "model_file": str(model), "model_sha256": sha, "model_size": size,
            "base_model": {"name": str(base), "kind": "file", "path": str(base), "sha256": "a" * 64},
            "date": "2026-10-23T12:00:00+00:00", "germquant_version": "0.2.0.dev0",
            "training": {"epochs": 1, "n_train": 2, "n_val": 1},
            "orientations": ["xy"],
            "training_set": [{"image": "image/a_xy_0.tif", "labels": "labels/a_xy_0.tif", "split": "train"}],
            "holdout_metrics": {"iou_threshold": 0.5},
            "licence_note": LICENCE}


@pytest.fixture
def trained(pack_env, monkeypatch):  # noqa: F811
    """A project folder holding a model trained in round r01 and its card, written as UTF-8 with raw non-ASCII
    text; the test pack is found through $GERMQUANT_PACKS."""
    use_fixture_packs(monkeypatch)
    proj = pack_env.tmp / "proj"
    mdir = proj / "training" / "r01" / "models" / "r01_nuclei"
    mdir.mkdir(parents=True)
    model = mdir / "r01_nuclei"
    model.write_bytes(MODEL_BYTES)
    base = pack_env.tmp / "weights" / "base_model_file"
    card = _card(model, base)
    cpath = mdir / "r01_nuclei.card.json"
    cpath.write_text(json.dumps(card, indent=1, ensure_ascii=False), encoding="utf-8")
    return SimpleNamespace(env=pack_env, proj=proj, packs=proj / "packs", model=model, card=cpath, card_data=card,
                           folder=proj / "packs" / DERIVED / "0.1.0")


def _derive(t, *extra, out=None, model=None, card=None) -> int:
    return cli.main(["pack", "derive", "--from", "test_pack@0.1.0", "--recipe", "with_model",
                     "--model", str(model or t.model), "--card", str(card or t.card), "--out", str(out or t.packs),
                     *extra])


def _leaves(obj, path=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from _leaves(v, f"{path}.{k}" if path else k)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _leaves(v, f"{path}[{i}]")
    elif isinstance(obj, str):
        yield path, obj


def _tables(out: Path) -> dict[str, pd.DataFrame]:
    return {f.name: pd.read_csv(f, low_memory=False, float_precision="round_trip").drop(
        columns=["run_timestamp", "config_hash"], errors="ignore") for f in sorted(out.glob(f"{IID}__*.csv"))}


# ---- acceptance 10: the derived pack verifies, lists and runs ----------------------------------------------------

def test_derived_pack_verifies_lists_and_runs(trained, monkeypatch, tmp_path, capsys):
    t = trained
    assert _derive(t, "--json") == 0
    rec = json.loads(capsys.readouterr().out)
    folder = t.folder
    assert rec["ok"] and rec["pack"] == f"{DERIVED}@0.1.0" and Path(rec["folder"]) == folder
    assert rec["status"] == "derived" and rec["model"]["name"] == MODEL_NAME and rec["card"] == CARD_REL
    assert sorted(p.relative_to(folder).as_posix() for p in folder.rglob("*") if p.is_file()) == [
        "README.md", "SHA256SUMS", "config/profiles/with_model.yaml", CARD_REL, "pack.yaml"]
    # the recipe is the extends line and the model, nothing else (D26)
    assert yaml.safe_load((folder / "config" / "profiles" / "with_model.yaml").read_text(encoding="ascii")) == {
        "extends": "pack:test_pack@0.1.0/with_model", "segmentation": {"nuclei": {"cellpose_model": MODEL_NAME}}}
    pack = packs.load_pack(folder)
    parent = packs.find_pack("test_pack@0.1.0")
    assert pack.problems == [] and pack.id == DERIVED and pack.version == "0.1.0" and pack.root == "pack"
    assert pack.object_kind == parent.object_kind and pack.data["requires"] == parent.data["requires"]
    assert pack.data["derived_from"] == {"pack": "test_pack@0.1.0", "recipe": "with_model",
                                         "pack_sha256": packs.sums_sha256(TEST_PACK)}
    r = pack.recipe()
    assert r.id == "with_model" and r.default and r.status == "provisional" and r.validation is None
    assert r.need("models") == [MODEL_NAME] and r.need("segment_backends") == ["cellpose"]
    assert r.tunables == parent.recipe("with_model").tunables
    assert r.expected_config_hash == rec["config_hash"] == load_config(folder / r.profile).hash
    assert r.expected_config_hash != parent.recipe("with_model").expected_config_hash
    assert pack.model(MODEL_NAME).as_dict() == {
        "name": MODEL_NAME, "kind": "cellpose", "file_name": "r01_nuclei", "sha256": MODEL_SHA,
        "size": len(MODEL_BYTES), "sources": [], "card": CARD_REL,
        "licence": t.card_data["licence_note"]}
    # the parent's version is copied beside it, as new_project copies it (design 3.8)
    assert packs.verify_sums(t.packs / "test_pack" / "0.1.0") == []
    assert sorted(p.name for p in t.packs.iterdir()) == sorted([DERIVED, "test_pack"])  # no staging folder left
    # verifies
    assert packs.verify_pack(folder) == [] and cli.main(["pack", "verify", str(folder)]) == 0
    # lists: through the torch-free entry with --project, and only with it
    capsys.readouterr()
    assert cli_packs.main(["list", "--json", "--all", "--project", str(t.proj)]) == 0
    listing = json.loads(capsys.readouterr().out)
    rows = {(p["id"], p["version"]): p for p in listing["packs"]}
    assert rows[(DERIVED, "0.1.0")]["source"] == "project" and rows[(DERIVED, "0.1.0")]["state"] == "ok"
    assert Path(rows[(DERIVED, "0.1.0")]["folder"]) == folder
    assert rows[("test_pack", "0.1.0")]["source"] == "project"                 # the copy wins, identical: quiet
    assert listing["roots"][0] == {"kind": "project", "path": str(t.packs)}
    assert DERIVED not in {p["id"] for p in packs.list_packs()}
    # usable here: the model is in the store and verified
    assert cli_packs.main(["check", f"{DERIVED}@0.1.0", "--project", str(t.proj), "--json"]) == 0
    st = json.loads(capsys.readouterr().out)
    assert st["usable"] and [m["name"] for m in st["models"]] == [MODEL_NAME]
    assert st["models"][0]["state"] in ("trusted", "verified")
    # runs: the derived recipe and its parent through `run --pack`, with the cellpose backend replaced by the
    # classical watershed (no Cellpose weights in tests); the model still resolves and is pin-checked
    monkeypatch.setattr(SB, "_REGISTRY", dict(SB._REGISTRY))
    monkeypatch.setattr(SB, "_SOURCES", dict(SB._SOURCES))     # restored afterwards, as both dicts must be
    SB.register_backend("cellpose", SB.classical_backend, replace=True)
    patch_nd2_reader(monkeypatch, small_stack())
    tiny = tmp_path / "tinynet"
    tiny.write_bytes(TINYNET_BYTES)
    packs.fetch_pack_models("test_pack", "with_model", {TINYNET: tiny})       # the parent's own (stand-in) model
    out_d, out_p = tmp_path / "run_derived", tmp_path / "run_parent"
    with packs.use_project(t.proj):
        assert cli.main(["run", NAME, "--pack", f"{DERIVED}@0.1.0", "--out", str(out_d)]) == 0
        assert cli.main(["run", NAME, "--pack", "test_pack@0.1.0", "--recipe", "with_model",
                         "--out", str(out_p)]) == 0
    man = json.loads((out_d / "run_manifest.json").read_text(encoding="utf-8"))
    assert man["model_sha256"] == MODEL_SHA and man["config_hash"] == r.expected_config_hash
    assert man["pack"]["id"] == DERIVED and man["pack"]["recipe"] == "with_model" and man["pack"]["source"] == "project"
    assert json.loads((out_p / "run_manifest.json").read_text(encoding="utf-8"))["model_sha256"] != MODEL_SHA
    a, b = _tables(out_d), _tables(out_p)
    assert a.keys() == b.keys() and len(a) >= 3
    for name in a:                     # the same recipe with another model: the same classical masks here
        pd.testing.assert_frame_equal(a[name], b[name], check_exact=True, obj=name)


def test_cards_portable_ascii(trained):
    t = trained
    assert _derive(t) == 0
    raw = (t.folder / CARD_REL).read_bytes()
    raw.decode("ascii")                                     # ASCII only: non-ASCII text is written as escapes
    card = json.loads(raw)
    assert card["name"] == "r01_nuclei" and card["model_file"] == "r01_nuclei"
    assert card["model_sha256"] == MODEL_SHA and card["model_size"] == len(MODEL_BYTES)
    assert card["base_model"] == {"name": "base_model_file", "kind": "file", "path": "base_model_file",
                                  "sha256": "a" * 64}
    assert DASH in card["what"] and MICRO in card["what"]                      # escaped, never dropped
    assert "train-nuclei from base_model_file on 3 slice" in card["what"]       # the paths in the text: their names
    assert card["what"].endswith("; pairs exported to pairs.")
    assert card["licence_note"] == LICENCE                                     # the URL is not a path
    assert card["training_set"] == t.card_data["training_set"]                 # relative paths stay
    manifest = (t.folder / "pack.yaml").read_bytes()
    manifest.decode("ascii")
    assert packs.load_pack(t.folder).model(MODEL_NAME).licence == LICENCE      # YAML escapes, read back
    for key, value in _leaves(card):
        assert str(t.env.tmp) not in value and t.env.tmp.as_posix() not in value, key
    res = MR.validate_card(card, packs.load_pack(t.folder).model(MODEL_NAME).as_dict(), pack_written=True)
    assert res == {"ok": True, "problems": [], "non_portable": []}
    src_sha = hashlib.sha256(t.card.read_bytes()).hexdigest()
    assert card["derived"] == {"pack": f"{DERIVED}@0.1.0", "from": "test_pack@0.1.0", "recipe": "with_model",
                               "card_sha256": src_sha,
                               "rewritten": ["base_model.name", "base_model.path", "model_file", "what"]}
    for key in ("training", "orientations", "holdout_metrics", "date"):
        assert card[key] == t.card_data[key]
    for f in t.folder.rglob("*"):                           # every file of the derived pack is ASCII
        if f.is_file():
            f.read_bytes().decode("ascii")
    # the training folder's card is left as it was
    assert t.card.read_text(encoding="utf-8") == json.dumps(t.card_data, indent=1, ensure_ascii=False)


def test_portable_card_rewrites_the_paths_of_every_platform_and_keeps_the_rest():
    card = {"name": "m", "model_file": "/home/lab/models/m",
            "base_model": {"path": "\\\\nas\\share\\cp\\cpsam_v2", "sha256": "c" * 64},
            "notes": ["C:/Users/lab/run 1/x.tif", "~/models/m", "D:\\data\\set"],
            "what": "from C:\\Users\\lab\\m (see https://example.org/a/b), slices xy/xz/yz, 1/2 of /srv/data/x, "
                    "D:\\data\\settings\\y and C:/Users/lab/run 1/x.tif.",
            "keep": "image/a_xy_0.tif"}
    out, changed = packs.portable_card(card, name="m", file_name="m", sha256="d" * 64, size=5)
    assert out["model_file"] == "m" and out["base_model"] == {"path": "cpsam_v2", "sha256": "c" * 64}
    assert out["notes"] == ["x.tif", "m", "set"]
    # a path with spaces goes whole (it is a value elsewhere in the card), a value never cuts a longer path
    assert out["what"] == "from m (see https://example.org/a/b), slices xy/xz/yz, 1/2 of x, y and x.tif."
    assert out["keep"] == "image/a_xy_0.tif" and out["model_sha256"] == "d" * 64 and out["model_size"] == 5
    assert changed == ["base_model.path", "model_file", "model_sha256", "model_size", "notes[0]", "notes[1]",
                       "notes[2]", "what"]
    assert card["model_file"] == "/home/lab/models/m"                       # the input is not changed


def test_lint_reports_only_the_outside_recipe_problems(trained):
    t = trained
    assert _derive(t) == 0
    rep = packs.lint_pack(t.folder)
    for p in rep["problems"]:
        assert not re.search(r"ASCII|dash|machine path|must name a pack version|does not load|SHA256SUMS|"
                             r"manifest|expected_config_hash|card", p), p
    # the outside rules count only the keys set inside the pack folder (packs/lint.py), and a derived recipe holds
    # only its extends line and the model: these problems are expected, and nothing else is reported
    switches = ("germline", "envelope", "audit", "axis", "staging", "sc.trace", "granule", "coloc", "partition",
                "granule_tail")
    expected = sorted(
        ["recipe with_model: no recipe golden: example/expected/recipe/with_model/ lacks nuclei.csv, spots.csv, "
         "image_summary.csv, config_hash.txt (every recipe ships one, rule R7; write it with `germquant pack golden "
         f"{DERIVED} --recipe with_model --update`)",
         "recipe with_model: segmentation.nuclei.method must be set to cellpose (any backend) or classical, never "
         "auto (got 'cellpose')"]
        + [f"recipe with_model: worm stage switch {s}.enabled is not set explicitly" for s in switches])
    assert sorted(rep["problems"]) == expected
    assert rep["warnings"] == []


# ---- the model store, and another computer ----------------------------------------------------------------------

def test_model_in_the_blob_store_and_fetched_from_the_training_folder_elsewhere(trained, monkeypatch, tmp_path,
                                                                                 capsys):
    t = trained
    assert _derive(t) == 0
    blob = t.env.blobs / MODEL_SHA / "r01_nuclei"
    assert blob.read_bytes() == MODEL_BYTES
    rec = json.loads((t.env.blobs / MODEL_SHA / f"{MODEL_SHA}.json").read_text(encoding="utf-8"))
    assert rec["sha256"] == MODEL_SHA and rec["size"] == len(MODEL_BYTES) and rec["name"] == MODEL_NAME
    assert rec["packs"] == [f"{DERIVED}@0.1.0"] and rec["source"].startswith("pack derive ")
    with packs.use_project(t.proj):
        assert MR.resolve_model_name(MODEL_NAME) == blob                      # the lookup and the pin check
    with pytest.raises(FileNotFoundError, match="no installed or built-in pack declares it"):
        MR.resolve_model_name(MODEL_NAME)                                     # the pack is the project's
    # another computer: an empty model store; the spec names no source, so --from-file brings the model
    other = tmp_path / "other_store"
    monkeypatch.setenv("GERMQUANT_MODELS", str(other))
    capsys.readouterr()
    assert cli_packs.main(["check", f"{DERIVED}@0.1.0", "--project", str(t.proj)]) == 2
    assert f"fix: germquant pack fetch-models {DERIVED}" in capsys.readouterr().out
    assert cli_packs.main(["fetch-models", DERIVED, "--project", str(t.proj), "--no-download"]) == 1
    assert cli_packs.main(["fetch-models", DERIVED, "--project", str(t.proj),
                           "--from-file", f"{MODEL_NAME}={t.model}"]) == 0
    assert (other / "models" / "blobs" / MODEL_SHA / "r01_nuclei").read_bytes() == MODEL_BYTES
    assert cli_packs.main(["check", f"{DERIVED}@0.1.0", "--project", str(t.proj)]) == 0


# ---- --project on list, show, check and fetch-models (I1) --------------------------------------------------------

def test_project_option_on_list_show_check_and_fetch_models(trained, capsys):
    t = trained
    assert _derive(t) == 0
    capsys.readouterr()
    assert cli_packs.main(["show", DERIVED]) == 1                 # without --project the project is not searched
    assert f"no pack '{DERIVED}'" in capsys.readouterr().err
    assert cli_packs.main(["show", DERIVED, "--project", str(t.proj), "--json"]) == 0
    shown = json.loads(capsys.readouterr().out)
    assert shown["source"] == "project" and [r["id"] for r in shown["recipes"]] == ["with_model"]
    assert [m["name"] for m in shown["models"]] == [MODEL_NAME] and shown["pack_sha256"] == packs.sums_sha256(t.folder)
    assert cli_packs.main(["show", DERIVED, "--project", str(t.proj)]) == 0
    assert f"model {MODEL_NAME} (cellpose, {len(MODEL_BYTES)} bytes" in capsys.readouterr().out
    assert cli_packs.main(["list", "--project", str(t.proj)]) == 0
    out = capsys.readouterr().out
    assert f"project {t.packs}" in out.splitlines()[0] and f"{DERIVED} 0.1.0" in out
    assert cli_packs.main(["check", DERIVED, "--project", str(t.proj)]) == 0
    assert "usable" in capsys.readouterr().out
    assert cli_packs.main(["fetch-models", DERIVED, "--project", str(t.proj), "--no-download"]) == 0
    assert "already in the store" in capsys.readouterr().out
    # the full command line takes the option too; a folder that does not exist is an error
    assert cli.main(["pack", "list", "--json", "--project", str(t.proj)]) == 0
    assert DERIVED in {p["id"] for p in json.loads(capsys.readouterr().out)["packs"]}
    assert cli_packs.main(["list", "--project", str(t.proj / "nope")]) == 1
    assert "no such folder" in capsys.readouterr().err


# ---- refusals and a second derive -----------------------------------------------------------------------------

def test_derive_again_is_a_no_op_and_refusals_write_nothing(trained, tmp_path, capsys):
    t = trained
    base = ["pack", "derive", "--model", str(t.model), "--card", str(t.card), "--out", str(t.packs)]

    def refused(argv, needle):
        capsys.readouterr()
        assert cli.main(argv) == 1, argv
        err = capsys.readouterr().err
        assert "REFUSED" in err and needle in err, err

    refused([*base, "--from", "test_pack@0.1.0", "--recipe", "basic"], "never runs a Cellpose model")
    refused([*base, "--from", "test_pack@0.1.0", "--recipe", "nope"], "has no recipe 'nope'")
    refused([*base, "--from", "no_such_pack@0.1.0", "--recipe", "with_model"], "no pack 'no_such_pack@0.1.0'")
    refused([*base, "--from", "test_pack@0.1.0", "--recipe", "with_model", "--id", "Bad-Id"], "--id 'Bad-Id'")
    refused([*base, "--from", "test_pack@0.1.0", "--recipe", "with_model", "--id", "test_pack"], "the parent's id")
    wrong = t.card.with_name("wrong.card.json")
    wrong.write_text(json.dumps(_card(t.model, t.env.tmp / "w", sha="0" * 64)), encoding="utf-8")
    refused(["pack", "derive", "--from", "test_pack@0.1.0", "--recipe", "with_model", "--model", str(t.model),
             "--card", str(wrong), "--out", str(t.packs)], "model_sha256")
    refused([*base[:-1], str(ROOT / "packs"), "--from", "test_pack@0.1.0", "--recipe", "with_model"],
            "the engine's own packs folder")
    assert not t.packs.exists() and not t.env.blobs.exists()                  # nothing was written
    assert not (ROOT / "packs" / DERIVED).exists()
    capsys.readouterr()
    assert cli.main([*base, "--from", "test_pack@0.1.0", "--recipe", "with_model", "--json", "--id", "bad id"]) == 1
    answer = json.loads(capsys.readouterr().out)
    assert answer["ok"] is False and answer["error"].startswith("--id 'bad id' does not match")
    # the first derive, then the same inputs again: nothing changes
    assert _derive(t) == 0
    files = {p: p.read_bytes() for p in t.packs.rglob("*") if p.is_file()}
    capsys.readouterr()
    assert _derive(t, "--json") == 0
    assert json.loads(capsys.readouterr().out)["status"] == "already present"
    assert {p: p.read_bytes() for p in t.packs.rglob("*") if p.is_file()} == files
    # other weights under an id and version that already hold a derived pack: refused, the copy kept
    other = tmp_path / "other_model"
    other.write_bytes(b"other weights\n")
    oc = tmp_path / "other.card.json"
    oc.write_text(json.dumps(_card(other, t.env.tmp / "w", sha=hashlib.sha256(b"other weights\n").hexdigest(),
                                   size=14)), encoding="utf-8")
    capsys.readouterr()
    assert _derive(t, "--id", DERIVED, model=other, card=oc) == 1
    assert "already holds" in capsys.readouterr().err
    assert {p: p.read_bytes() for p in t.packs.rglob("*") if p.is_file()} == files


# ---- the version features (I4) --------------------------------------------------------------------------------

def test_version_json_carries_the_features_present(capsys, tmp_path):
    assert cli_packs.main(["--version", "--json"]) == 0
    info = json.loads(capsys.readouterr().out)
    feats = info["features"]
    assert feats == packs.features() and feats[0] == "packs" and "pack_derive" in feats
    assert feats == [f for f in CONTRACT if f in feats]                       # the contract's names, in its order
    names = [f for f, _, _ in packs.FEATURE_PROBES]
    assert names == [f for f in CONTRACT if f in names]
    assert packs.version_info()["features"] == feats
    # each feature is listed exactly when its code is there: a stand-in package folder, filled step by step
    pkg = tmp_path / "germquant"
    assert packs.features(pkg) == []
    for name, rel, attr in packs.FEATURE_PROBES:
        f = pkg / rel
        f.parent.mkdir(parents=True, exist_ok=True)
        if attr:
            f.write_text(f"# {attr} = 1\ndef other():\n    {attr} = 1\n", encoding="ascii")    # not top level
            assert name not in packs.features(pkg), name
            f.write_text(f"def {attr}(x):\n    return x\n" if attr.islower() else f"{attr} = 1.0\n", encoding="ascii")
        else:
            f.write_text("", encoding="ascii")
        assert packs.features(pkg)[-1] == name
    assert packs.features(pkg) == names
