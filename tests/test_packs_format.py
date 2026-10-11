"""Pack format 1, SHA256SUMS, discovery and recipe resolution (docs/DESIGN_system_packs.md 3.2 to 3.8, package P1).

Covers: validate_manifest's refusals (format major, id and name patterns, model digests and pinned names,
release names, validated without a record), LF-normalised sums that agree for CRLF and LF copies, junk never
hashed, the search order of pack_roots with built-in roots limited to the checkout and $GERMQUANT_CONFIG_ROOT,
shadowing (a built-in cannot be shadowed, an identical copy is not reported), a broken pack listed as
unusable, find_pack and resolve_recipe errors that name the packs available, and pack_record.
"""
from __future__ import annotations

import copy
import os

import pytest
import yaml
from packs_helpers import CELEGANS_PACK, TEST_PACK, copy_pack, pack_env, use_fixture_packs  # noqa: F401

from germquant import packs
from germquant.config import builtin_roots
from germquant.packs import format as F


def _manifest(path=TEST_PACK) -> dict:
    return yaml.safe_load((path / "pack.yaml").read_text(encoding="ascii"))


def test_both_shipped_manifests_are_valid():
    assert F.validate_manifest(_manifest(TEST_PACK), TEST_PACK) == []
    assert F.validate_manifest(_manifest(CELEGANS_PACK), CELEGANS_PACK) == []


@pytest.mark.parametrize("edit, needle", [
    (lambda d: d.update(format=2), "format"),
    (lambda d: d.update(id="Bad-Id"), "pack id"),
    (lambda d: d.update(version="1.0"), "semantic version"),
    (lambda d: d["recipes"][0].update(id="Basic Recipe"), "recipe id"),
    (lambda d: d["recipes"][0].update(expected_config_hash="xyz"), "expected_config_hash"),
    (lambda d: d["recipes"][0].update(dims="4d"), "dims"),
    (lambda d: d["recipes"][0].update(profile="../outside.yaml"), "relative path"),
    (lambda d: d["recipes"][0].update(profile="C:/abs/x.yaml"), "relative path"),
    (lambda d: d["models"][0].pop("sha256"), "sha256 missing"),
    (lambda d: d["models"][0].pop("size"), "size missing"),
    (lambda d: d["models"][0].update(name="tinynet@00000000"), "is not the start of its sha256"),
    (lambda d: d["models"][0].update(name="cpsam_v2"), "pack model pattern"),                  # a Cellpose name
    (lambda d: d["models"][0].update(name="germline_nuclei_combined"), "pack model pattern"),  # a release name
    (lambda d: d["models"][0].update(name="germline_nuclei_combined", kind="release"), "release record"),
    (lambda d: d["recipes"][1]["needs"].update(models=["other@12345678"]), "does not declare"),
    (lambda d: d["recipes"][0].update(status="validated"), "validation record"),
    (lambda d: d["recipes"][0].update(validated=True), "validation record"),
    (lambda d: d.update(requires={"germquant": "not a range"}), "version range"),
])
def test_validate_manifest_refuses(edit, needle):
    data = copy.deepcopy(_manifest())
    edit(data)
    probs = F.validate_manifest(data, TEST_PACK)
    assert any(needle in p for p in probs), probs


def test_validated_needs_a_passing_record(tmp_path):
    pack = copy_pack(TEST_PACK, tmp_path / "test_pack")
    data = _manifest(pack)
    data["recipes"][0].update(status="validated", validation="validation/basic.json")
    (pack / "validation").mkdir()
    (pack / "validation" / "basic.json").write_text('{"verdict": "fail"}', encoding="ascii")
    assert any("validation record" in p for p in F.validate_manifest(data, pack))
    (pack / "validation" / "basic.json").write_text('{"verdict": {"pass": true}}', encoding="ascii")
    assert F.validate_manifest(data, pack) == []


def test_sums_are_lf_normalised_and_skip_junk(tmp_path):
    lf = copy_pack(TEST_PACK, tmp_path / "lf" / "test_pack", rebuild=False)
    crlf = copy_pack(TEST_PACK, tmp_path / "crlf" / "test_pack", rebuild=False)
    for folder, eol in ((lf, b"\n"), (crlf, b"\r\n")):
        for rel in F.pack_files(folder):
            if F.is_text(rel):
                p = folder / rel
                p.write_bytes(p.read_bytes().replace(b"\r\n", b"\n").replace(b"\n", eol))
    assert F.verify_sums(lf) == [] and F.verify_sums(crlf) == []
    assert F.sums_sha256(lf) == F.sums_sha256(crlf) == F.sums_sha256(TEST_PACK)
    assert (crlf / "README.md").read_bytes() != (lf / "README.md").read_bytes()
    # junk is never hashed and never makes verify fail
    (crlf / ".DS_Store").write_bytes(b"junk")
    (crlf / "__pycache__").mkdir()
    (crlf / "__pycache__" / "x.cpython-311.pyc").write_bytes(b"junk")
    (crlf / "config" / "y.pyc").write_bytes(b"junk")
    assert F.verify_sums(crlf) == []
    text = F.build_sums_text(crlf)
    assert ".DS_Store" not in text and "pyc" not in text and "SHA256SUMS" not in text
    # a changed, missing or unlisted file is reported
    (crlf / "README.md").write_bytes(b"changed\n")
    (crlf / "extra.txt").write_bytes(b"new\n")
    (crlf / "CHANGELOG.md").unlink()
    probs = F.verify_sums(crlf)
    assert any("README.md: content differs" in p for p in probs)
    assert any("extra.txt: not listed" in p for p in probs)
    assert any("CHANGELOG.md: listed" in p for p in probs)


def test_sums_text_format_is_sorted_two_space_lf():
    text = (TEST_PACK / "SHA256SUMS").read_bytes()
    assert b"\r" not in text.replace(b"\r\n", b"\n")
    lines = text.replace(b"\r\n", b"\n").decode("ascii").splitlines()
    paths = [ln.split("  ", 1)[1] for ln in lines]
    assert paths == sorted(paths) and all(len(ln.split("  ", 1)[0]) == 64 for ln in lines)
    assert all("\\" not in p for p in paths)


def test_builtin_roots_are_the_checkout_and_the_env_root_only(pack_env, monkeypatch, tmp_path):  # noqa: F811
    from germquant import provenance

    roots = builtin_roots()
    assert roots and (roots[0] / "src" / "germquant").is_dir()
    assert all(r != tmp_path.resolve() for r in roots)
    monkeypatch.chdir(tmp_path)                       # the working directory is never a root
    assert builtin_roots() == roots
    env_root = tmp_path / "container_root"
    env_root.mkdir()
    monkeypatch.setenv("GERMQUANT_CONFIG_ROOT", str(env_root))
    assert builtin_roots() == [*roots, env_root.resolve()]
    # a repo root without src/germquant (another project's folder) is not a built-in root
    other = tmp_path / "otherproject"
    other.mkdir()
    (other / "pyproject.toml").write_text("[project]\nname='x'\n", encoding="ascii")
    monkeypatch.setattr(provenance, "_repo_root", lambda: other)
    assert builtin_roots() == [env_root.resolve()]


def test_pack_roots_order_and_dev_layout(pack_env, monkeypatch, tmp_path):  # noqa: F811
    a, b = tmp_path / "pa", tmp_path / "pb"
    a.mkdir()
    b.mkdir()
    monkeypatch.setenv("GERMQUANT_PACKS", f"{a}{os.pathsep}{b}")
    roots = packs.pack_roots()
    kinds = [k for k, _ in roots]
    assert kinds[0] == "builtin" and kinds[-1] == "user" and kinds.count("GERMQUANT_PACKS") == 2
    assert packs.install_root() == a
    with packs.use_project(tmp_path / "proj"):
        assert packs.pack_roots()[0] == ("project", tmp_path / "proj" / "packs")
    # a pack under development (<root>/<id>/pack.yaml) and an installed one (<root>/<id>/<version>/pack.yaml)
    copy_pack(TEST_PACK, a / "test_pack")
    copy_pack(TEST_PACK, b / "test_pack" / "0.2.0", edit=lambda d: d.update(version="0.2.0"))
    packs.refresh()
    found = {p.ref: p.source for p in packs.discover().packs}
    assert found["test_pack@0.1.0"] == "GERMQUANT_PACKS" and found["test_pack@0.2.0"] == "GERMQUANT_PACKS"
    assert found["celegans_germline@0.1.0"] == "builtin"
    # without a version: the newest version compatible with this engine
    assert packs.find_pack("test_pack").version == "0.2.0"
    assert packs.find_pack("test_pack@0.1.0").version == "0.1.0"


def test_builtin_cannot_be_shadowed_and_identical_copies_are_quiet(pack_env, monkeypatch, tmp_path):  # noqa: F811
    root = tmp_path / "env_packs"
    copy_pack(CELEGANS_PACK, root / "celegans_germline" / "0.1.0")                       # identical
    copy_pack(CELEGANS_PACK, root / "celegans_germline_other" / "0.1.0")
    monkeypatch.setenv("GERMQUANT_PACKS", str(root))
    packs.refresh()
    disc = packs.discover()
    assert [p.source for p in disc.get("celegans_germline")] == ["builtin"]
    assert disc.shadowed == []                                   # the identical copy is not reported
    # a different copy under the same id and version is listed as shadowed and ignored
    changed = root / "celegans_germline" / "0.1.0"
    (changed / "README.md").write_bytes((changed / "README.md").read_bytes() + b"\nedited\n")
    packs.write_sums(changed)
    packs.refresh()
    disc = packs.discover()
    assert disc.get("celegans_germline")[0].source == "builtin"
    assert [p.folder for p in disc.shadowed] == [changed]
    rows = packs.list_packs(all=True)
    assert any(r["state"] == "shadowed" and r["id"] == "celegans_germline" for r in rows)


def test_broken_pack_is_listed_unusable_never_raised(pack_env, monkeypatch, tmp_path):  # noqa: F811
    root = tmp_path / "env_packs"
    (root / "broken" / "0.1.0").mkdir(parents=True)
    (root / "broken" / "0.1.0" / "pack.yaml").write_text("format: 1\nid: [unclosed\n", encoding="ascii")
    copy_pack(TEST_PACK, root / "badid" / "0.1.0", edit=lambda d: d.update(id="Bad-Id"))
    monkeypatch.setenv("GERMQUANT_PACKS", str(root))
    packs.refresh()
    rows = packs.list_packs(all=True)
    bad = [r for r in rows if r["state"] == "unusable"]
    assert len(bad) == 2 and all(r["problems"] for r in bad)
    assert "celegans_germline" in {r["id"] for r in packs.list_packs()}        # the good ones still list
    with pytest.raises(FileNotFoundError, match="unusable"):
        packs.find_pack("broken")


def test_find_pack_and_resolve_recipe_errors_name_the_packs(pack_env, monkeypatch):  # noqa: F811
    use_fixture_packs(monkeypatch)
    assert packs.resolve_recipe("test_pack@0.1.0/basic") == TEST_PACK / "config" / "profiles" / "basic.yaml"
    assert packs.resolve_recipe("test_pack/with_model").name == "with_model.yaml"
    with pytest.raises(FileNotFoundError, match=r"Packs available: .*test_pack@0\.1\.0"):
        packs.find_pack("nope")
    with pytest.raises(FileNotFoundError, match="its recipes: basic, with_model"):
        packs.resolve_recipe("test_pack@0.1.0/nope")
    with pytest.raises(FileNotFoundError, match="no pack 'test_pack@9.9.9'"):
        packs.resolve_recipe("test_pack@9.9.9/basic")
    assert packs.parse_ref("a@1.2.3/r") == ("a", "1.2.3", "r") and packs.parse_ref("a") == ("a", None, None)


def test_pack_record_and_listing(pack_env, monkeypatch):  # noqa: F811
    use_fixture_packs(monkeypatch)
    rec = packs.pack_record("test_pack", "with_model")
    assert rec == {"id": "test_pack", "version": "0.1.0", "recipe": "with_model", "source": "GERMQUANT_PACKS",
                   "pack_sha256": F.sums_sha256(TEST_PACK), "object_kind": "nucleus", "dims": "3d",
                   "status": "provisional"}
    assert packs.pack_record("celegans_germline")["recipe"] == "default"
    row = next(r for r in packs.list_packs() if r["id"] == "test_pack")
    assert row["state"] == "ok" and [r["id"] for r in row["recipes"]] == ["basic", "with_model"]
    assert packs.verify_pack("test_pack") == [] and packs.verify_pack(CELEGANS_PACK) == []
