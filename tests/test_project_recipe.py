"""Recipe edits, channel roles and project models (docs/PROJECTS.md; design 7.3 and 7.5; PLAN.md 4.4 P5-PRJ,
decisions D19 and D42; P5 acceptance 3 and 10).

* An unmodified project has no recipe.yaml and hashes as `--pack` and the pack's expected_config_hash; an edit
  writes recipe.yaml with a versioned ``extends: pack:<id>@<version>/<recipe>`` and a new hash.
* preview_roles returns exactly what write_recipe writes; a channel named with a Greek gamma round-trips as a
  YAML escape and still binds; write_recipe validates by origin; revert brings back no recipe.yaml and the
  pack's hash; use_model writes the relative path and records the digest and the card.
* `project edit --dry-run` writes nothing and prints what preview_roles returns; edit's --z-range, --runner,
  --interpreter, --device, --xy-stride, --samples and --no-<stage> write the manifest (D42).
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest
import yaml
from packs_helpers import ROOT, pack_env  # noqa: F401
from project_helpers import (
    BASIC_HASH,
    GAMMA,
    TEST_PACK_REF,
    last_json,
    new_test_project,
    project_env,  # noqa: F401
    read_manifest_file,
    snapshot,
    write_ome_tiff,
)

import germquant
from germquant import cli_projects, stages
from germquant.config import load_config
from germquant.project import (
    ProjectError,
    add_images,
    new_project,
    new_round,
    open_project,
    preview_roles,
    project_config,
    revert_recipe,
    use_model,
    write_recipe,
)
from germquant.project.headers import read_header

FIXED = {tool: "PINNED" for tools in stages._STAGE_TOOLS.values() for tool in tools}
ROLES = {"dna": 0, "foci": 1, "central_element": None, "granule": None, "lamin": None}


def _recipe(root: Path) -> dict:
    return yaml.safe_load((root / "recipe.yaml").read_bytes().decode("ascii"))


def test_unmodified_project_hash_equals_pack(project_env, tmp_path):  # noqa: F811
    """P5 acceptance 3, first half: no recipe.yaml, and the hash (and every stage hash) of `--pack`."""
    from germquant import api

    for recipe, expected in (("basic", BASIC_HASH), ("with_model", "81e38c5cf371")):
        proj = new_test_project(tmp_path / recipe, recipe=recipe)
        assert not (proj.root / "recipe.yaml").exists()
        cfg = project_config(proj)
        ref = api.resolve_config(pack=TEST_PACK_REF, recipe=recipe)
        assert cfg.hash == ref.hash == expected == proj.manifest["recipe"]["expected_config_hash"]
        assert stages.stage_hashes(cfg, ROLES, FIXED) == stages.stage_hashes(ref, ROLES, FIXED)
        assert {k: v for k, v in cfg.pack_record.items() if k not in ("source", "edited")} == \
            {k: v for k, v in ref.pack_record.items() if k != "source"}
    worm = new_project(tmp_path / "worm", "celegans_germline@0.1.0")
    assert not (worm.root / "recipe.yaml").exists()
    assert project_config(worm).hash == "84bdda5a5230" == load_config(ROOT / "config" / "config.yaml").hash


def test_edit_writes_versioned_extends(project_env, tmp_path):  # noqa: F811
    """P5 acceptance 3, second half: an edit writes recipe.yaml with a versioned pack: line and a new hash."""
    proj = new_test_project(tmp_path / "p")
    h = write_recipe(proj, edits={"spots.log.min_snr": 3.0, "segmentation.nuclei.min_volume_um3": 3})
    raw = (proj.root / "recipe.yaml").read_bytes()
    assert raw.decode("ascii").splitlines()[0] == "extends: pack:test_pack@0.1.0/basic" and b"\r" not in raw
    assert b"    min_volume_um3: 3.0\n" in raw                    # the int 3 is written as the declared float
    assert _recipe(proj.root) == {"extends": "pack:test_pack@0.1.0/basic",
                                  "segmentation": {"nuclei": {"min_volume_um3": 3.0}},   # the declared float type
                                  "spots": {"log": {"min_snr": 3.0}}}
    cfg = project_config(proj)
    assert h != BASIC_HASH and h == cfg.hash == load_config(proj.root / "recipe.yaml").hash
    assert cfg.get("spots.log.min_snr") == 3.0 and cfg.pack_record["edited"] is True
    m = read_manifest_file(proj.root)
    assert m["recipe"]["file"] == "recipe.yaml" and m["recipe"]["expected_config_hash"] == h
    assert m["recipe"]["origins"] == {"segmentation.nuclei.min_volume_um3": "form", "spots.log.min_snr": "form"}
    assert m["validated"]["value"] is False and "edited" in m["validated"]["why"]
    assert write_recipe(proj, edits={"spots.log.min_snr": 3.0}) == h                  # the same edit again
    h2 = write_recipe(proj, edits={"segmentation.nuclei.min_volume_um3": 2.0})         # back to the pack's value
    assert "segmentation" not in _recipe(proj.root) and h2 not in (h, BASIC_HASH)
    other = new_test_project(tmp_path / "another folder")
    assert write_recipe(other, edits={"spots.log.min_snr": 3.0}) == h2               # folders do not enter the hash
    assert stages.stage_hashes(project_config(other), ROLES, FIXED) == \
        stages.stage_hashes(project_config(proj), ROLES, FIXED)


def test_preview_roles_returns_exactly_what_write_recipe_writes(project_env, tmp_path):  # noqa: F811
    proj = new_test_project(tmp_path / "p")
    a = write_ome_tiff(tmp_path / "in" / "a.ome.tif", names=("Hoechst 405", "Cy3 RAD-51"))
    b = write_ome_tiff(tmp_path / "in" / "b.ome.tif", names=("DAPI", "Cy3 RAD-51"))
    add_images(proj, [a, b])
    roles = {"dna": "Hoechst 405", "foci": "Cy3 RAD-51"}
    pv = preview_roles(proj, roles)
    assert set(pv) == {"yaml", "bindings", "problems"} and pv["problems"] == []
    assert pv["bindings"]["a"] == {"dna": {"channel": "Hoechst 405", "how": "name"},
                                   "foci": {"channel": "Cy3 RAD-51", "how": "name"}}
    assert pv["bindings"]["b"]["dna"] == {"channel": "DAPI", "how": "index"}          # the pack's index fallback
    cm = yaml.safe_load(pv["yaml"])
    assert cm["match_by"] == "name" and list(cm["roles"]) == ["dna", "foci"]
    assert cm["roles"]["dna"] == {"names": ["Hoechst 405"], "index": 0, "required": True, "marker": "DNA"}
    assert cm["roles"]["foci"] == {"names": ["Cy3 RAD-51"], "index": 1, "required": False, "marker": "foci"}
    h = write_recipe(proj, roles=roles, origin="roles")
    assert (proj.root / "channel_map.yaml").read_bytes() == pv["yaml"].encode("ascii")
    assert _recipe(proj.root) == {"extends": "pack:test_pack@0.1.0/basic", "io": {"channel_map": "channel_map.yaml"}}
    cfg = project_config(proj)
    assert cfg.hash == h != BASIC_HASH and cfg.channel_map.resolve(["Hoechst 405", "Cy3 RAD-51"]) == (
        {"dna": 0, "foci": 1}, [])
    assert read_manifest_file(proj.root)["recipe"]["roles"] == roles
    one = preview_roles(proj, roles, image_id="b")
    assert list(one["bindings"]) == ["b"] and one["yaml"] == pv["yaml"]
    clash = {"dna": "Cy3 RAD-51", "foci": "Cy3 RAD-51"}
    assert any("same channel" in p for p in preview_roles(proj, clash)["problems"])
    with pytest.raises(ProjectError, match="same channel"):
        write_recipe(proj, roles=clash, origin="roles")
    assert any("not a role" in p for p in preview_roles(proj, {"rad51": "Cy3 RAD-51"})["problems"])
    with pytest.raises(ProjectError, match="image"):
        preview_roles(proj, roles, image_id="nope")
    assert (proj.root / "channel_map.yaml").read_bytes() == pv["yaml"].encode("ascii")   # refusals wrote nothing
    assert write_recipe(proj, roles={}, origin="roles") == BASIC_HASH                   # back to the pack map
    assert not (proj.root / "channel_map.yaml").exists() and not (proj.root / "recipe.yaml").exists()


def test_a_greek_gamma_channel_round_trips_as_an_escape_and_still_binds(project_env, tmp_path, capsys):  # noqa: F811
    name = GAMMA + "H2AX"
    img = write_ome_tiff(tmp_path / "in" / "g.ome.tif", names=("DAPI", name))
    proj = new_test_project(tmp_path / "p", images=[img])
    write_recipe(proj, roles={"foci": name}, origin="roles")
    raw = (proj.root / "channel_map.yaml").read_bytes()
    raw.decode("ascii")
    assert b"\\u03B3H2AX" in raw
    cfg = project_config(proj)
    assert cfg.channel_map.roles["foci"].names == [name]
    assert cfg.channel_map.resolve(read_header(img)["channel_names"]) == ({"dna": 0, "foci": 1}, [])
    assert preview_roles(proj, {"foci": name})["bindings"]["g"]["foci"] == {"channel": name, "how": "name"}
    rc = cli_projects.main(["edit", str(proj.root), "--role", f"foci={name}", "--dry-run", "--json"])
    out = last_json(capsys.readouterr().out)
    assert rc == 0 and out["preview"]["bindings"]["g"]["foci"]["channel"] == name


def test_write_recipe_validates_by_origin(project_env, tmp_path):  # noqa: F811
    """D19: form edits only declared tunables, scalar, of the declared type, inside min, max and choices; tuning
    only spot keys; model only segmentation.nuclei.cellpose_model; roles the channel map."""
    proj = new_test_project(tmp_path / "p")
    for edits, needle in (({"spots.effect_size_min": 2.0}, "not a tunable"),
                          ({"spots.log.min_snr": "high"}, "float"),
                          ({"spots.log.min_snr": True}, "float"),
                          ({"spots.log.min_snr": 99.0}, "at most 10"),
                          ({"spots.log.min_snr": 0.5}, "at least 1"),
                          ({"spots.log": {"min_snr": 3.0}}, "not a tunable")):
        with pytest.raises(ProjectError, match=needle):
            write_recipe(proj, edits=edits)
    with pytest.raises(ProjectError, match="tuning"):
        write_recipe(proj, edits={"segmentation.nuclei.min_volume_um3": 3.0}, origin="tuning")
    with pytest.raises(ProjectError, match="cellpose_model"):
        write_recipe(proj, edits={"spots.log.min_snr": 3.0}, origin="model")
    with pytest.raises(ProjectError, match="origin"):
        write_recipe(proj, edits={"spots.log.min_snr": 3.0}, origin="gui")
    with pytest.raises(ProjectError, match="does not set"):
        write_recipe(proj, unset=["spots.log.min_snr"])
    with pytest.raises(ProjectError, match="dotted"):
        write_recipe(proj, edits={"spots..min_snr": 3.0}, origin="tuning")
    assert not (proj.root / "recipe.yaml").exists()
    h = write_recipe(proj, edits={"spots.log.min_snr": 4.0, "spots.log.diameter_um": 0.6}, origin="tuning")
    assert read_manifest_file(proj.root)["recipe"]["origins"] == {"spots.log.diameter_um": "tuning",
                                                                  "spots.log.min_snr": "tuning"}
    assert h == project_config(proj).hash
    cc = new_project(tmp_path / "cc", "cultured_cells@0.1.0")
    with pytest.raises(ProjectError, match="choices"):
        write_recipe(cc, edits={"segmentation.nuclei.project": "median"})
    assert write_recipe(cc, edits={"segmentation.nuclei.project": "mean"}) != cc.manifest["pack"]["sha256"]
    assert _recipe(cc.root)["segmentation"] == {"nuclei": {"project": "mean"}}


def test_revert_brings_back_no_recipe_and_the_pack_hash(project_env, tmp_path):  # noqa: F811
    img = write_ome_tiff(tmp_path / "in" / "a.ome.tif")
    proj = new_test_project(tmp_path / "p", images=[img])
    write_recipe(proj, edits={"spots.log.min_snr": 3.0})
    write_recipe(proj, roles={"foci": "foci"}, origin="roles")
    assert (proj.root / "recipe.yaml").exists() and (proj.root / "channel_map.yaml").exists()
    assert revert_recipe(proj, plugin_version="0.2.0.dev1") == BASIC_HASH
    assert not (proj.root / "recipe.yaml").exists() and not (proj.root / "channel_map.yaml").exists()
    m = read_manifest_file(proj.root)
    assert m["recipe"] == {"file": None, "expected_config_hash": BASIC_HASH, "overrides": {}, "roles": {},
                           "origins": {}}
    assert m["saved"]["plugin"] == "0.2.0.dev1" and project_config(proj).pack_record["edited"] is False
    write_recipe(proj, edits={"spots.log.min_snr": 3.0})
    assert write_recipe(proj, unset=["spots.log.min_snr"]) == BASIC_HASH                 # the last edit removed
    assert not (proj.root / "recipe.yaml").exists()


def test_a_run_setting_never_rewrites_the_recipe_files(project_env, tmp_path):  # noqa: F811
    """A change that touches no recipe edit (device, stride, stage switches) keeps recipe.yaml byte for byte, a hand
    edit included; the next recipe edit rewrites it from the manifest's pin, keeping the hand-set key."""
    from germquant.project import save_project

    proj = new_test_project(tmp_path / "p")
    write_recipe(proj, edits={"spots.log.min_snr": 3.0})
    rp = proj.root / "recipe.yaml"
    hand = rp.read_bytes() + b"# a note\nqc:\n  min_nuclei: 2\n"
    rp.write_bytes(hand)
    h = project_config(open_project(proj.root)).hash
    save_project(proj, device="cpu", xy_stride=2)
    assert rp.read_bytes() == hand and project_config(open_project(proj.root)).hash == h
    assert read_manifest_file(proj.root)["recipe"]["expected_config_hash"] == h
    crlf = hand.replace(b"\n", b"\r\n")                                 # saved by a Windows editor
    rp.write_bytes(crlf)
    import io
    from contextlib import redirect_stdout

    for argv, changed in ((["--json"], False), (["--device", "auto", "--json"], True)):
        buf = io.StringIO()
        with redirect_stdout(buf):
            assert cli_projects.main(["edit", str(proj.root), *argv]) == 0
        assert json.loads(buf.getvalue())["changed"] is changed, argv
        assert rp.read_bytes() == crlf and project_config(open_project(proj.root)).hash == h
    write_recipe(proj, edits={"segmentation.nuclei.min_volume_um3": 3.0})
    data = _recipe(proj.root)
    assert data["qc"] == {"min_nuclei": 2} and b"# a note" not in rp.read_bytes()
    assert read_manifest_file(proj.root)["recipe"]["origins"]["qc.min_nuclei"] == "manual"


def test_recipe_files_the_manifest_does_not_record_are_kept(project_env, tmp_path, capsys):  # noqa: F811
    """A recipe.yaml the manifest does not record (copied from another project) is not used by runs, and open says
    so; a run setting or a revert never removes it, and an edit that would replace it is refused (a dry run warns),
    unless it holds exactly the bytes the edit writes. An unrecorded channel_map.yaml is kept from a role edit."""
    from germquant.project import save_project

    donor = new_test_project(tmp_path / "donor")
    write_recipe(donor, edits={"spots.log.min_snr": 3.0})
    copied = (donor.root / "recipe.yaml").read_bytes()
    proj = new_test_project(tmp_path / "p")
    mine = proj.root / "recipe.yaml"
    hand = copied.replace(b"min_snr: 3.0", b"min_snr: 4.0")
    mine.write_bytes(hand)
    opened = open_project(proj.root)
    assert opened.config_hash == BASIC_HASH and project_config(opened).get("spots.log.min_snr") == 2.0
    assert [w for w in opened.warnings if "recipe.yaml" in w and "do not use it" in w], opened.warnings
    save_project(opened, device="cpu")                                          # a run setting keeps it
    assert mine.read_bytes() == hand
    assert revert_recipe(opened) == BASIC_HASH and mine.read_bytes() == hand    # so does a revert
    with pytest.raises(ProjectError, match="not recorded in the manifest"):
        write_recipe(opened, edits={"spots.log.min_snr": 3.0})
    assert mine.read_bytes() == hand and read_manifest_file(proj.root)["recipe"]["file"] is None
    rc = cli_projects.main(["edit", str(proj.root), "--set", "spots.log.min_snr=3", "--dry-run", "--json"])
    out = last_json(capsys.readouterr().out)
    assert rc == 0 and any("not recorded in the manifest" in w for w in out["warnings"]), out
    rc = cli_projects.main(["edit", str(proj.root), "--device", "auto", "--json"])
    out = last_json(capsys.readouterr().out)
    assert rc == 0 and out["changed"] is True and out["removed"] == [] and mine.read_bytes() == hand
    assert cli_projects.main(["edit", str(proj.root), "--json"]) == 0
    assert last_json(capsys.readouterr().out)["changed"] is False                # the file is no change either
    mine.write_bytes(copied)              # the bytes the edit writes (what a write stopped before its manifest leaves)
    assert write_recipe(opened, edits={"spots.log.min_snr": 3.0}) == project_config(donor).hash
    assert read_manifest_file(proj.root)["recipe"]["file"] == "recipe.yaml"
    write_recipe(opened, unset=["spots.log.min_snr"])
    cm, own = proj.root / "channel_map.yaml", b"match_by: name\nroles: {}\n"
    cm.write_bytes(own)
    with pytest.raises(ProjectError, match="role editor"):
        write_recipe(opened, roles={"foci": "Cy3 RAD-51"}, origin="roles")
    revert_recipe(opened)
    assert cm.read_bytes() == own and not mine.exists()


def _model(proj, name="m1", content=b"fine-tuned weights stand-in\n"):
    mdir = proj.root / "training" / "r01" / "models" / name
    mdir.mkdir(parents=True, exist_ok=True)
    model = mdir / name
    model.write_bytes(content * 8)
    sha = hashlib.sha256(model.read_bytes()).hexdigest()
    card = mdir / f"{name}.card.json"
    card.write_text(json.dumps({"name": name, "model_file": str(model), "model_sha256": sha,
                                "model_size": model.stat().st_size, "base_model": {"name": "cpsam"}}, indent=1),
                    encoding="ascii")
    return model, sha, card


def test_use_model_writes_relative_path(project_env, tmp_path, capsys):  # noqa: F811
    """P5 acceptance 10 (use-model): the relative path in recipe.yaml, the digest and the card in the manifest."""
    proj = new_test_project(tmp_path / "p", recipe="with_model")
    new_round(proj)
    model, sha, card = _model(proj)
    h = use_model(proj, "training/r01/models/m1/m1", plugin_version="0.2.0.dev1")
    assert _recipe(proj.root) == {"extends": "pack:test_pack@0.1.0/with_model",
                                  "segmentation": {"nuclei": {"cellpose_model": "training/r01/models/m1/m1"}}}
    m = read_manifest_file(proj.root)
    card_sha = hashlib.sha256(card.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    assert m["models"] == [{"key": "segmentation.nuclei.cellpose_model", "ref": "training/r01/models/m1/m1",
                            "sha256": sha, "card": "training/r01/models/m1/m1.card.json", "card_sha256": card_sha,
                            "source": "project"}]
    assert m["recipe"]["origins"] == {"segmentation.nuclei.cellpose_model": "model"}
    assert m["saved"]["plugin"] == "0.2.0.dev1" and m["validated"]["value"] is False
    assert h == project_config(proj).hash != "81e38c5cf371"
    from germquant.pipeline import resolve_model_path

    assert Path(resolve_model_path(project_config(proj))) == model           # the project's own file
    assert use_model(proj, model) == h                                       # an absolute path inside: relative
    outside = tmp_path / "elsewhere" / "m2"
    outside.parent.mkdir()
    outside.write_bytes(b"x")
    with pytest.raises(ProjectError, match="inside the project"):
        use_model(proj, outside)
    with pytest.raises(ProjectError, match="not found"):
        use_model(proj, "training/r01/models/m1/nothing")
    other, _, other_card = _model(proj, "m2", b"other weights\n")
    other_card.write_text(json.dumps({"name": "m2", "model_sha256": "0" * 64, "model_size": 3}), encoding="ascii")
    with pytest.raises(ProjectError, match="card"):
        use_model(proj, other)
    rc = cli_projects.main(["use-model", str(proj.root), "--model", "training/r01/models/m1/m1", "--plugin-version",
                            "0.2.0.dev9", "--json"])
    out = last_json(capsys.readouterr().out)
    assert rc == 0 and out["config_hash"] == h and out["model"]["sha256"] == sha
    assert read_manifest_file(proj.root)["saved"]["plugin"] == "0.2.0.dev9"


def test_edit_dry_run_writes_nothing_and_prints_what_preview_roles_returns(project_env, tmp_path,  # noqa: F811
                                                                          capsys):
    img = write_ome_tiff(tmp_path / "in" / "a.ome.tif", names=("DAPI", "RAD-51 Cy3"))
    proj = new_test_project(tmp_path / "p", images=[img])
    before = snapshot(proj.root)
    rc = cli_projects.main(["edit", str(proj.root), "--role", "foci=RAD-51 Cy3", "--set", "spots.log.min_snr=3",
                            "--dry-run", "--json"])
    out = last_json(capsys.readouterr().out)
    assert rc == 0 and snapshot(proj.root) == before
    assert out["dry_run"] is True and out["changed"] is True and out["written"] == []
    assert out["preview"] == preview_roles(proj, {"foci": "RAD-51 Cy3"})
    assert out["previous_config_hash"] == BASIC_HASH
    h = write_recipe(proj, edits={"spots.log.min_snr": 3.0}, roles={"foci": "RAD-51 Cy3"})
    assert out["config_hash"] == h
    assert out["recipe_yaml"] == (proj.root / "recipe.yaml").read_text(encoding="ascii")
    assert out["channel_map_yaml"] == (proj.root / "channel_map.yaml").read_text(encoding="ascii")
    rc = cli_projects.main(["edit", str(proj.root), "--role", "foci=", "--unset", "spots.log.min_snr", "--json"])
    out = last_json(capsys.readouterr().out)
    assert rc == 0 and out["config_hash"] == BASIC_HASH and sorted(out["removed"]) == ["channel_map.yaml",
                                                                                          "recipe.yaml"]
    rc = cli_projects.main(["edit", str(proj.root), "--revert", "--set", "spots.log.min_snr=3", "--json"])
    assert rc == 2 and "--revert" in last_json(capsys.readouterr().out)["error"]
    # a role choice with a problem: the dry run shows it (exit 0, the preview with its problems); the edit refuses it
    clash = {"dna": "RAD-51 Cy3", "foci": "RAD-51 Cy3"}
    before = snapshot(proj.root)
    argv = ["edit", str(proj.root), "--role", "dna=RAD-51 Cy3", "--role", "foci=RAD-51 Cy3", "--json"]
    rc = cli_projects.main([*argv, "--dry-run"])
    out = last_json(capsys.readouterr().out)
    assert rc == 0 and out["dry_run"] is True and out["written"] == [] and snapshot(proj.root) == before
    assert out["preview"] == preview_roles(proj, clash) and any("same channel" in p for p in out["preview"]["problems"])
    assert out["preview"]["bindings"]["a"]["foci"] == {"channel": "RAD-51 Cy3", "how": "name"}
    assert cli_projects.main(argv[:-1] + ["--dry-run"]) == 0
    assert "PROBLEM: image a: the roles dna and foci bind the same channel" in capsys.readouterr().out
    rc = cli_projects.main(argv)
    out = last_json(capsys.readouterr().out)
    assert rc == 2 and "same channel" in out["error"] and out["kind"] == "ProjectError"
    assert snapshot(proj.root) == before


def test_edit_writes_the_provenance_options(project_env, tmp_path, capsys):  # noqa: F811
    """D42: --z-range, --runner and --interpreter RECIPE=PYTHON write run.z_range, run.runner and
    run.interpreters; --device, --xy-stride, --samples and --no-<stage> write the rest of the manifest."""
    from germquant import api
    from germquant.packs import version_info

    proj = new_test_project(tmp_path / "p")
    (proj.root / "samples.csv").write_text("image_id,genotype\na,N2\n", encoding="ascii")
    rc = cli_projects.main(["edit", str(proj.root), "--z-range", "1", "3", "--runner", "subprocess", "--interpreter",
                            f"basic={sys.executable}", "--device", "cpu", "--xy-stride", "2", "--samples",
                            "samples.csv", "--no-spots", "--plugin-version", "0.2.0.dev1", "--json"])
    out = last_json(capsys.readouterr().out)
    assert rc == 0, out
    m = read_manifest_file(proj.root)
    info = version_info()
    assert m["run"] == {"xy_stride": 2, "z_range": [1, 3], "device": "cpu", "runner": "subprocess",
                        "interpreters": {"basic": {"python": sys.executable, "germquant": germquant.__version__,
                                                   "git_sha": info["git_sha"]}}}
    assert m["samples"] == "samples.csv" and m["recipe"]["overrides"] == {"spots.enabled": False}
    assert m["saved"]["plugin"] == "0.2.0.dev1" and m["validated"]["value"] is False
    cfg = project_config(open_project(proj.root))
    ref = api.resolve_config(pack=TEST_PACK_REF)
    ref.set("spots.enabled", False)                                     # what `run --pack ... --no-spots` hashes
    assert cfg.get("spots.enabled") is False and cfg.hash == ref.hash == m["recipe"]["expected_config_hash"]
    assert out["config_hash"] == cfg.hash and not (proj.root / "recipe.yaml").exists()
    rc = cli_projects.main(["edit", str(proj.root), "--z-range", "none", "--runner", "none", "--device", "none",
                            "--interpreter", "basic=none", "--samples", "none", "--unset", "spots.enabled", "--json"])
    assert rc == 0, capsys.readouterr()
    m = read_manifest_file(proj.root)
    assert m["run"] == {"xy_stride": 2, "z_range": None, "device": None, "runner": None, "interpreters": {}}
    assert m["samples"] is None and m["recipe"]["overrides"] == {} and m["recipe"]["expected_config_hash"] == BASIC_HASH
    capsys.readouterr()
    for argv, needle in ((["--z-range", "3", "1"], "z range"), (["--xy-stride", "0"], "stride"),
                         (["--interpreter", f"nope={sys.executable}"], "recipe"),
                         (["--interpreter", f"basic={tmp_path / 'no_python.exe'}"], "interpreter"),
                         (["--samples", str(tmp_path / "elsewhere.csv")], "samples"),
                         (["--samples", "missing.csv"], "samples"), (["--device", "tpu"], "invalid choice")):
        with pytest.raises(SystemExit) if needle == "invalid choice" else _noop():
            rc = cli_projects.main(["edit", str(proj.root), *argv, "--json"])
            assert rc == 2, argv
            assert needle in last_json(capsys.readouterr().out)["error"], argv


class _noop:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False
