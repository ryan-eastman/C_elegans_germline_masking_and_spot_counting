"""The two outside packs, cultured_cells and budding_yeast 0.1.0 (design 3.12 and 6.9; package P3, acceptance 3, 6, 7
and 9).

* Both packs load, verify against their SHA256SUMS and lint clean; every recipe is provisional with
  validated: false, extends the pack's frozen pack_base.yaml and hashes to its expected_config_hash.
* The frozen bases carry the QC settings of a field that is not a gonad and no deny-listed key; every spots
  instance sets expected_from_germ_cell: false; the yeast transmitted role has no index fallback.
* The recipe goldens pass with an empty model store; an edit of a pack_base.yaml fails them; edits of the engine's
  generic.yaml or config.yaml move no outside pack hash; a misspelt key in a pack_base.yaml fails lint.
* A stack with no transmitted channel fails the yeast recipe with the require_channel message instead of
  segmenting channel 0.
* The cards say what is provisional and which scale effect_size is on; the real examples record what was checked
  on the first download; every CITATION.cff has the fields the CFF 1.2.0 schema requires.
* The yeast recipe golden sees the recipe's own spot threshold and smallest-cell gate (a weak spot and a small
  object in the synthetic example).
* Lint asks every outside recipe for require_channel and a current recipe golden, every real example for a
  permissive data licence, and keeps large files and stray TIFFs out of a pack; `pack golden` fails a recipe
  without its golden and exits 1 when no golden ran; a frozen base is never rewritten silently; the cpsam
  comparison script refuses to run on another device than the one asked for.
"""
from __future__ import annotations

import importlib.util
import json
import shutil
from pathlib import Path

import pytest
import yaml
from packs_helpers import ROOT, copy_pack, pack_env  # noqa: F401  (pack_env is a fixture)

from germquant import cli, packs, packs_golden
from germquant.config import load_config
from germquant.io.channel_map import ChannelMap
from germquant.models_registry import validate_card
from germquant.packs import lint as L
from germquant.pipeline import process_image

CULTURED, YEAST = ROOT / "packs" / "cultured_cells", ROOT / "packs" / "budding_yeast"
OUTSIDE = {"cultured_cells": CULTURED, "budding_yeast": YEAST}
RECIPES = {"cultured_cells": {"nuclei_foci": "projected", "nuclei_foci_3d": "3d"},
           "budding_yeast": {"cells_spots": "projected"}}
MAPS = {"cultured_cells": "config/channel_maps/cultured_dna_foci.yaml",
        "budding_yeast": "config/channel_maps/yeast_transmitted_foci.yaml"}
MIN_NUCLEI = {"cultured_cells": 1, "budding_yeast": 10}
MODEL = "cpsam_v2@0f1cc3f7"
CPSAM_V2_SHA = "0f1cc3f7ecdd8a037a57c6c48d9d8921391be4cbce3fa9f13c3e3a2e1253c667"


def _get(d, dotted, default=None):
    node = d
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return default
        node = node[part]
    return node


@pytest.mark.parametrize("pid", list(OUTSIDE))
def test_outside_packs_load_verify_lint_and_are_provisional(pid):
    folder = OUTSIDE[pid]
    pack = packs.load_pack(folder)
    assert pack.problems == [] and pack.id == pid and pack.version == "0.1.0" and pack.root == "pack"
    assert pack.object_kind == {"cultured_cells": "nucleus", "budding_yeast": "cell"}[pid]
    assert packs.validate_pack(pack) == [] and packs.verify_pack(folder) == []
    rep = packs.lint_pack(folder)
    assert rep["problems"] == [], rep
    assert {r.id: r.dims for r in pack.recipes} == RECIPES[pid]
    assert pack.default_recipe().id == next(iter(RECIPES[pid]))
    for r in pack.recipes:
        assert r.status == "provisional" and r.raw.get("validated") is False and r.validation is None
        assert r.need("models") == [MODEL] and r.need("segment_backends") == ["cellpose"]
        assert r.need("spot_backends") == ["spotmax"]
        cfg = load_config(packs.recipe_profile_path(pack, r))
        assert cfg.hash == r.expected_config_hash
        recipe_text = (folder / r.profile).read_text(encoding="ascii")
        assert yaml.safe_load(recipe_text)["extends"] == "pack_base.yaml"
        assert cfg.get("segmentation.nuclei.method") == "cellpose"
        assert cfg.get("segmentation.nuclei.cellpose_model") == MODEL
        assert cfg.get("segmentation.nuclei.require_channel") is True
    (model,) = pack.models
    assert model.name == MODEL and model.sha256 == CPSAM_V2_SHA and model.size == 1233586851
    assert any(s.get("type") == "cellpose_cache" for s in model.sources)
    urls = [s["url"] for s in model.sources if s.get("type") == "url"]
    assert urls and all("/resolve/7c61431b5fbb078f3296754bd15d9f51b320f837/" in u for u in urls)  # a pinned commit
    card = json.loads((folder / model.card).read_text(encoding="ascii"))
    res = validate_card(card, model.as_dict(), pack_written=True)
    assert res["ok"] and res["non_portable"] == [], res
    assert "CC BY-NC" in card["licence_note"]


@pytest.mark.parametrize("pid", list(OUTSIDE))
def test_frozen_base_holds_the_qc_settings_and_no_deny_listed_key(pid):
    folder = OUTSIDE[pid]
    text = (folder / "config" / "profiles" / "pack_base.yaml").read_bytes().decode("ascii")
    base = yaml.safe_load(text)
    assert "extends" not in base and not any(line.lstrip().startswith("#") for line in text.splitlines())
    for key in L.DENY_KEYS:
        assert _get(base, key) is None, key                       # spots.channel among them
    assert "channel" not in base["spots"]
    assert base["io"]["channel_map"] == MAPS[pid]
    assert base["metadata"]["filename_regex"] == L.NEUTRAL_REGEX
    assert base["qc"]["germline_outlier_factor"] >= L.FRAMING_OFF_MIN
    assert base["qc"]["min_nuclei"] == MIN_NUCLEI[pid]
    for sw in L.WORM_SWITCHES:
        assert _get(base, sw) is False, sw
    pack = packs.load_pack(folder)
    for r in pack.recipes:
        merged = load_config(packs.recipe_profile_path(pack, r)).as_dict()
        for inst in _get(merged, "spots.instances", []) or []:
            assert inst["expected_from_germ_cell"] is False
        assert "channel" not in merged["spots"]


def test_recipe_values_are_the_design_starting_values():
    y = load_config(YEAST / "config" / "profiles" / "cells_spots.yaml")
    assert {k: y.get(f"segmentation.nuclei.{k}") for k in ("channel", "project", "project_plane", "min_area_um2",
                                                            "diameter_um")} == \
        {"channel": "transmitted", "project": "plane", "project_plane": "sharpest", "min_area_um2": 3.0,
         "diameter_um": 4.0}
    assert {k: y.get(f"spots.{k}") for k in ("spot_radius_um", "spot_radius_z_um", "gauss_sigma_um",
                                             "thresholding_method", "effect_size_min", "remove_hot_pixels",
                                             "sharpen", "merge_z_columns")} == \
        {"spot_radius_um": 0.244, "spot_radius_z_um": 1.0, "gauss_sigma_um": 0.054,
         "thresholding_method": "threshold_triangle", "effect_size_min": 2.5, "remove_hot_pixels": True,
         "sharpen": True, "merge_z_columns": False}
    assert y.get("spots.instances") in (None, [])
    for recipe, extra in (("nuclei_foci", {"project": "max", "min_area_um2": 20.0}),
                          ("nuclei_foci_3d", {"min_volume_um3": 60.0})):
        c = load_config(CULTURED / "config" / "profiles" / f"{recipe}.yaml")
        assert c.get("segmentation.nuclei.channel") == "dna" and c.get("segmentation.nuclei.diameter_um") == 12.0
        for k, v in extra.items():
            assert c.get(f"segmentation.nuclei.{k}") == v
        assert {k: c.get(f"spots.{k}") for k in ("spot_radius_um", "spot_radius_z_um", "gauss_sigma_um",
                                                 "thresholding_method", "effect_size_min", "merge_z_columns")} == \
            {"spot_radius_um": 0.3, "spot_radius_z_um": 0.6, "gauss_sigma_um": 0.08,
             "thresholding_method": "threshold_triangle", "effect_size_min": 3.0, "merge_z_columns": True}
        assert c.get("spots.sharpen") in (None, False)
        assert c.get("spots.instances") == [{"name": "foci2", "role": "foci2", "expected_from_germ_cell": False}]
    assert load_config(CULTURED / "config" / "profiles" / "nuclei_foci_3d.yaml").get("segmentation.nuclei.project") \
        in (None, "none")


def test_channel_maps():
    ym = ChannelMap.from_yaml(YEAST / MAPS["budding_yeast"])
    t = ym.roles["transmitted"]
    assert t.index is None and t.required and t.names == ["phase", "contr", "PhC", "Bright", "DIC", "trans"]
    assert "BF" not in t.names
    assert ym.roles["foci"].index == 1 and ym.roles["foci"].required
    assert ym.roles["dna"].index is None and not ym.roles["dna"].required
    assert ym.resolve(["Bright", "Cy3"]) == ({"transmitted": 0, "foci": 1, "dna": None}, [])
    roles, flags = ym.resolve(["phase_contr", "mKate", "mNeon"])          # the mtDNA positions of S-BIAD3750
    assert roles["transmitted"] == 0 and roles["foci"] == 2 and flags == []
    roles, flags = ym.resolve(["BFP", "Cy3"])                              # BF is left out: BFP is no bright field
    assert roles["transmitted"] is None and "channel:transmitted:MISSING_REQUIRED" in flags
    roles, flags = ym.resolve(["DAPI", "Cy3"])
    assert roles["transmitted"] is None and "channel:transmitted:MISSING_REQUIRED" in flags
    cm = ChannelMap.from_yaml(CULTURED / MAPS["cultured_cells"])
    assert (cm.roles["dna"].index, cm.roles["foci"].index, cm.roles["foci2"].index) == (0, 1, 2)
    assert cm.roles["dna"].required and not cm.roles["foci2"].required
    assert cm.resolve(["53BP1", "gH2AX", "DAPI"]) == ({"dna": 2, "foci": 1, "foci2": 0}, [])     # DeepFoci order
    roles, flags = cm.resolve(["405", "488", "561"])
    assert roles == {"dna": 0, "foci": 1, "foci2": 2} and len(flags) == 3


@pytest.mark.parametrize("pid, scale", [("budding_yeast", "sharpened"), ("cultured_cells", "unsharpened")])
def test_cards_say_what_is_provisional_and_which_scale_effect_size_is_on(pid, scale):
    card = (OUTSIDE[pid] / "README.md").read_text(encoding="ascii")
    for needle in ("provisional", "not validated", "validated: false", "effect_size", scale, "Object kind",
                   "object_kind"):
        assert needle in card, needle
    assert "What is provisional" in card and "## Limitations" in card and "CITATION.cff" in card
    assert "role editor" in card
    if pid == "budding_yeast":
        assert "MDN1" in card and "YeaZ" in card and "| cell | a whole cell |" in card
    else:
        assert "object_diameter_um" in card and "| nucleus | a nucleus |" in card
    for r in packs.load_pack(OUTSIDE[pid]).recipes:                    # the recipe text says it too
        text = (OUTSIDE[pid] / r.profile).read_text(encoding="ascii")
        assert "provisional" in text and scale.replace("un", "") in text.replace("unsharpened", "sharpened")


@pytest.mark.parametrize("pid", list(OUTSIDE))
def test_recipe_goldens_pass_with_an_empty_model_store(pack_env, pid):  # noqa: F811
    pack = packs.find_pack(pid)
    results = packs_golden.run_pack_goldens(pack)
    assert [r["recipe"] for r in results] == list(RECIPES[pid])
    assert all(r["ok"] for r in results), [r["problems"] for r in results]
    for r in results:
        sentinel = (pack.folder / "example" / "expected" / "recipe" / r["recipe"] / "config_hash.txt").read_text(
            encoding="ascii").strip()
        assert sentinel == r["config_hash"] == pack.recipe(r["recipe"]).expected_config_hash


def test_an_edit_of_pack_base_fails_the_recipe_goldens(tmp_path, pack_env):  # noqa: F811
    edited = copy_pack(YEAST, tmp_path / "budding_yeast", rebuild=False)
    base = edited / "config" / "profiles" / "pack_base.yaml"
    base.write_bytes(base.read_bytes() + b"# an edit\n")                 # even a comment moves the hash (rule R1)
    packs.write_sums(edited)
    res = packs_golden.run_recipe_golden(packs.load_pack(edited), "cells_spots")
    assert not res["ok"]
    assert any("!= the pack's expected_config_hash" in p for p in res["problems"])
    assert any(p.startswith("sentinel config_hash.txt") for p in res["problems"])
    # a value edit, with the hash brought up to date, still fails on the tables
    edited2 = copy_pack(YEAST, tmp_path / "v" / "budding_yeast", rebuild=False)
    base2 = edited2 / "config" / "profiles" / "pack_base.yaml"
    base2.write_bytes(base2.read_bytes().replace(b"  min_nuclei: 10\n", b"  min_nuclei: 4\n"))
    data = yaml.safe_load((edited2 / "pack.yaml").read_text(encoding="ascii"))
    new_hash = load_config(edited2 / "config" / "profiles" / "cells_spots.yaml").hash
    data["recipes"][0]["expected_config_hash"] = new_hash
    (edited2 / "pack.yaml").write_bytes(yaml.safe_dump(data, sort_keys=False).encode("ascii"))
    (edited2 / "example" / "expected" / "recipe" / "cells_spots" / "config_hash.txt").write_bytes(
        (new_hash + "\n").encode("ascii"))
    packs.write_sums(edited2)
    res2 = packs_golden.run_recipe_golden(packs.load_pack(edited2), "cells_spots")
    assert not res2["ok"] and "image_summary.qc_flags differs" in res2["problems"]


def test_engine_config_edits_move_no_outside_pack_hash(tmp_path, monkeypatch):
    """A scratch copy of the engine's config/ and of both packs: editing generic.yaml and config.yaml there moves the
    hash of a profile that extends generic.yaml (the control) and no outside recipe hash, even with the scratch copy
    as the built-in root, so no fallback could bind to the edited files."""
    from germquant import config as C

    root = tmp_path / "engine"
    shutil.copytree(ROOT / "config", root / "config")
    for pid, folder in OUTSIDE.items():
        shutil.copytree(folder, root / "packs" / pid)
    control = tmp_path / "study.yaml"
    control.write_bytes(f"extends: {(root / 'config' / 'profiles' / 'generic.yaml').as_posix()}\n".encode("ascii"))
    before = load_config(control).hash
    for rel in ("config/profiles/generic.yaml", "config/config.yaml"):
        p = root / rel
        p.write_bytes(p.read_bytes() + b"\n# scratch edit: the planned hash bump\n")
    monkeypatch.setattr(C, "builtin_roots", lambda: [root])
    assert load_config(control).hash != before                           # the edits are real
    for pid in OUTSIDE:
        pack = packs.load_pack(root / "packs" / pid)
        for r in pack.recipes:
            assert load_config(root / "packs" / pid / r.profile).hash == r.expected_config_hash


def test_a_misspelt_key_in_pack_base_fails_lint(tmp_path):
    def edited(where: Path, line: bytes) -> Path:
        pack = copy_pack(YEAST, where / "budding_yeast", rebuild=False)
        base = pack / "config" / "profiles" / "pack_base.yaml"
        base.write_bytes(base.read_bytes().replace(b"spots:\n", b"spots:\n" + line))
        data = yaml.safe_load((pack / "pack.yaml").read_text(encoding="ascii"))
        for r in data["recipes"]:
            r["expected_config_hash"] = load_config(pack / r["profile"]).hash
        (pack / "pack.yaml").write_bytes(yaml.safe_dump(data, sort_keys=False).encode("ascii"))
        packs.write_sums(pack)
        return pack

    typo = packs.lint_pack(edited(tmp_path / "a", b"  efect_size_min: 2.5\n"))["problems"]
    assert any("pack_base.yaml: key spots.efect_size_min: the engine does not read it" in p for p in typo), typo
    deny = packs.lint_pack(edited(tmp_path / "b", b"  channel: foci\n"))["problems"]
    assert any("pack_base.yaml: key spots.channel: nothing reads it" in p for p in deny), deny


def test_lint_keeps_a_packs_own_example_present_and_small(tmp_path):
    pack = copy_pack(YEAST, tmp_path / "budding_yeast", rebuild=False)
    ex = pack / "example" / "yeast_synthetic.ome.tif"
    ex.write_bytes(ex.read_bytes() + bytes(1_100_000))                # past 2 MB: real images are fetched
    packs.write_sums(pack)
    assert any("stays under 2,000,000" in p for p in packs.lint_pack(pack)["problems"])
    ex.unlink()
    packs.write_sums(pack)
    assert any("example/yeast_synthetic.ome.tif is missing" in p for p in packs.lint_pack(pack)["problems"])


def test_a_stack_without_a_transmitted_channel_fails_the_yeast_recipe(tmp_path, pack_env):  # noqa: F811
    pack = packs.find_pack("budding_yeast")
    example = pack.folder / "example" / "yeast_synthetic.ome.tif"

    def renamed(path, **kw):
        meta, stack = packs_golden.read_example(path, **kw)
        names = ["DAPI", "Cy3"]                                       # a fluorescence channel where Bright was
        meta = {**meta, "channel_names": names}
        stack.channel_names = names
        return meta, stack

    cfg = load_config(packs.recipe_profile_path(pack, "cells_spots"))
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"      # as the recipe golden runs it, without a model
    out = tmp_path / "no_transmitted"
    with pytest.raises(ValueError, match=r"require_channel is true .*segmentation role 'transmitted'"):
        process_image(example, cfg, out, prov=None, read_fn=renamed)
    assert not list(out.glob("*__done.json")) and not list(out.glob("*__nuclei.csv"))


def _edited_yeast(where: Path, old: bytes, new: bytes, *, rel: str = "config/profiles/cells_spots.yaml") -> Path:
    """A copy of the yeast pack with one edit, its hash, sentinel and sums brought up to date (so only the golden's
    tables, or lint, can catch the edit)."""
    pack = copy_pack(YEAST, where / "budding_yeast", rebuild=False)
    f = pack / rel
    text = f.read_bytes()
    assert old in text
    f.write_bytes(text.replace(old, new))
    data = yaml.safe_load((pack / "pack.yaml").read_text(encoding="ascii"))
    new_hash = load_config(pack / "config" / "profiles" / "cells_spots.yaml").hash
    data["recipes"][0]["expected_config_hash"] = new_hash
    (pack / "pack.yaml").write_bytes(yaml.safe_dump(data, sort_keys=False).encode("ascii"))
    sentinel = pack / "example" / "expected" / "recipe" / "cells_spots" / "config_hash.txt"
    if sentinel.parent.is_dir():
        sentinel.write_bytes((new_hash + "\n").encode("ascii"))
    packs.write_sums(pack)
    return pack


@pytest.mark.parametrize("old, new, needle", [
    (b"effect_size_min: 2.5", b"effect_size_min: 4.0", "spots: 19 rows, expected 20"),
    (b"min_area_um2: 3.0", b"min_area_um2: 1.0", "nuclei: 7 rows, expected 6"),
])
def test_the_yeast_golden_sees_its_spot_threshold_and_smallest_cell(tmp_path, pack_env, old, new, needle):  # noqa: F811
    """The synthetic example holds a weak spot (sharpened effect size about 3.2) and a small object (about 1.9 um2),
    so raising effect_size_min to 4 or lowering min_area_um2 to 1 moves the tables, not only the sentinel."""
    res = packs_golden.run_recipe_golden(packs.load_pack(_edited_yeast(tmp_path, old, new)), "cells_spots")
    assert not res["ok"] and needle in res["problems"], res["problems"]
    assert not any("config_hash" in p for p in res["problems"])          # the hash was brought up to date


def test_lint_asks_for_require_channel_and_a_current_recipe_golden(tmp_path, pack_env):  # noqa: F811
    no_req = _edited_yeast(tmp_path / "a", b"    require_channel: true\n", b"")
    assert load_config(no_req / "config" / "profiles" / "cells_spots.yaml").get(
        "segmentation.nuclei.require_channel") is None
    assert any("segmentation.nuclei.require_channel must be true" in p
               for p in packs.lint_pack(no_req)["problems"])
    stale = copy_pack(YEAST, tmp_path / "b" / "budding_yeast", rebuild=False)
    (stale / "example" / "expected" / "recipe" / "cells_spots" / "config_hash.txt").write_bytes(b"0123456789ab\n")
    packs.write_sums(stale)
    assert any("the recipe golden is stale" in p for p in packs.lint_pack(stale)["problems"])
    gone = copy_pack(YEAST, tmp_path / "c" / "budding_yeast", rebuild=False)
    shutil.rmtree(gone / "example" / "expected" / "recipe" / "cells_spots")
    packs.write_sums(gone)
    assert any("no recipe golden: example/expected/recipe/cells_spots/ lacks" in p
               for p in packs.lint_pack(gone)["problems"])
    (res,) = packs_golden.run_pack_goldens(packs.load_pack(gone))            # a failure, never a silent skip
    assert res["recipe"] == "cells_spots" and not res["ok"] and "no recipe golden" in res["problems"][0]


def test_pack_golden_exits_1_when_no_golden_ran(pack_env, capsys):  # noqa: F811
    assert cli.main(["pack", "golden", "celegans_germline"]) == 1
    assert "has no recipe golden" in capsys.readouterr().err


def test_lint_keeps_examples_permissive_and_large_files_and_stray_tiffs_out(tmp_path):
    pack = copy_pack(CULTURED, tmp_path / "cultured_cells", rebuild=False)
    text = (pack / "pack.yaml").read_bytes()
    (pack / "pack.yaml").write_bytes(text.replace(b"licence_spdx: CC-BY-4.0", b"licence_spdx: CC-BY-NC-4.0"))
    shutil.copyfile(pack / "example" / "cultured_synthetic.ome.tif", pack / "example" / "fetched.tif")
    (pack / "notes").mkdir()
    (pack / "notes" / "big.bin").write_bytes(bytes(L.EXAMPLE_MAX_BYTES))
    packs.write_sums(pack)
    probs = packs.lint_pack(pack)["problems"]
    assert any("licence_spdx 'CC-BY-NC-4.0' is not one of" in p for p in probs), probs
    assert any("example/fetched.tif: a TIFF that is not a declared local example" in p for p in probs), probs
    assert any("notes/big.bin is 2,000,000 bytes" in p for p in probs), probs


def test_every_citation_file_has_the_fields_cff_requires():
    """CFF 1.2.0: authors, cff-version, message and title at the root; authors, title and type on every reference."""
    files = sorted((ROOT / "packs").glob("*/CITATION.cff")) + [ROOT / "tests" / "fixtures" / "packs" / "test_pack"
                                                               / "CITATION.cff"]
    assert len(files) >= 4
    for f in files:
        cff = yaml.safe_load(f.read_text(encoding="ascii"))
        for key in ("authors", "cff-version", "message", "title"):
            assert cff.get(key), (f, key)
        for i, ref in enumerate(cff.get("references") or []):
            assert ref.get("authors") and ref.get("title") and ref.get("type"), (f, i, ref)
            assert all(isinstance(a, dict) and (a.get("family-names") or a.get("name")) for a in ref["authors"])
    yeast = yaml.safe_load((YEAST / "CITATION.cff").read_text(encoding="ascii"))
    dois = {r.get("doi") for r in yeast["references"]}
    assert {"10.6019/S-BIAD3750", "10.1101/2024.10.22.619610"} <= dois


def test_a_frozen_base_is_never_rewritten_silently(tmp_path):
    from germquant.packs import authoring as AU

    generic = ROOT / "config" / "profiles" / "generic.yaml"
    data = AU.frozen_base_data(generic, channel_map="config/channel_maps/x.yaml", min_nuclei=10)
    base = tmp_path / "pack_base.yaml"
    AU.write_frozen_base(base, data)
    assert AU.base_unchanged(base, data)
    AU.write_frozen_base(base, data)                                      # the same text: left as it is
    other = AU.frozen_base_data(generic, channel_map="config/channel_maps/x.yaml", min_nuclei=4)
    with pytest.raises(AU.FrozenBaseError, match=r"(?s)is frozen.*-  min_nuclei: 10.*\+  min_nuclei: 4"):
        AU.write_frozen_base(base, other)
    assert "min_nuclei: 10" in base.read_text(encoding="ascii")
    AU.write_frozen_base(base, other, rebase=True)
    assert "min_nuclei: 4" in base.read_text(encoding="ascii")


def test_the_cpsam_comparison_refuses_another_device_than_the_one_asked_for(monkeypatch):
    """A CUDA check must never pass on the CPU: with no accelerator (or $GERMQUANT_DEVICE pointing elsewhere), the
    script stops with exit code 2 before anything is loaded or run."""
    from germquant import device

    spec = importlib.util.spec_from_file_location("compare_cpsam_names", ROOT / "scripts" / "compare_cpsam_names.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    monkeypatch.setattr(device, "available", lambda: {"cuda": False, "mps": False, "cpu": True})
    monkeypatch.delenv("GERMQUANT_DEVICE", raising=False)
    assert mod.resolve_device("cpu") == ("cpu", "") and mod.resolve_device("auto") == ("cpu", "")
    dev, why = mod.resolve_device("cuda")
    assert dev == "cpu" and "cuda is not available here" in why
    assert mod.main(["--device", "cuda"]) == 2
    monkeypatch.setattr(device, "available", lambda: {"cuda": True, "mps": False, "cpu": True})
    monkeypatch.setenv("GERMQUANT_DEVICE", "cpu")
    dev, why = mod.resolve_device("cuda")
    assert dev == "cpu" and "$GERMQUANT_DEVICE=cpu overrides it" in why


def test_real_examples_record_what_was_checked():
    cult = {e["id"]: e for e in yaml.safe_load((CULTURED / "pack.yaml").read_text(encoding="ascii"))["examples"]}
    df = cult["deepfoci_valid_258"]
    assert df["base_url"] == "https://zenodo.org/api/records/4067741/files/" and df["licence"] == "CC BY 4.0"
    (src,) = df["sources"]
    assert src["path"] == "foci_detection_validation_subset.zip/content" and src["member"] == "valid/data_258.tif"
    assert (src["zip_size"], src["size"], src["crc32"]) == (10214754062, 145468367, "f3aa48ab")
    assert df["assemble"] == {"axes": "ZYXC", "channels": ["53BP1", "gH2AX", "DAPI"],
                              "spacing_zyx_um": [0.3, 0.0647, 0.0647], "crop": {"y": [224, 704], "x": [128, 1024]}}
    assert {"member", "page_order", "channel_order", "voxel_size", "crop"} <= set(df["checked"])
    assert "Vicar" in df["attribution"] and df["recipe"] == "nuclei_foci" and df["expect"]["n_objects"] == [1, 3]
    yex = {e["id"]: e for e in yaml.safe_load((YEAST / "pack.yaml").read_text(encoding="ascii"))["examples"]}
    sb = yex["sbiad3750_sut509_p1"]
    assert sb["base_url"] == "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/750/S-BIAD3750/Files/"
    assert sb["licence"] == "CC0" and "S-BIAD3750" in sb["attribution"]
    assert sb["assemble"] == {"axes": "ZYX", "channels": ["Bright", "Cy3"],
                              "spacing_zyx_um": [0.24, 0.0720635, 0.0720635]}
    assert [s["size"] for s in sb["sources"]] == [14701884, 14701884]
    for spec in (df, sb):
        assert len(spec["pixels_sha256"]) == 64 and all(len(s["sha256"]) == 64 for s in spec["sources"])
    for pid, folder in OUTSIDE.items():                                 # the synthetic examples stay small
        for ex in packs.list_examples(packs.load_pack(folder)):
            if ex["kind"] == "local":
                assert (folder / ex["image"]).stat().st_size < 2_000_000
