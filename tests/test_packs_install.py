"""Install, vendor, line endings and git (design 3.7 and 3.8, P1 acceptance 7 and part of 12).

install refuses zip-slip entries, drive letters, symlinks and Windows junctions (which build-sums and verify
never follow), bad sums, non-ASCII hashed files anywhere in the pack, a recipe outside config/profiles/,
faulty tunables and a collision with a built-in id and version; a zip's junk is never copied; an installed
version is immutable, and a reinstall over a damaged copy is refused; the SHA256SUMS of celegans_germline
verifies on CRLF and LF copies; every path in each committed SHA256SUMS is tracked by git; in a scratch
repository with the real .gitattributes and .gitignore, a PNG under packs/ is byte-identical after a fresh
checkout, cards and examples are tracked, and .DS_Store and caches are not; packs.vendor copies a pack version
that then verifies.
"""
from __future__ import annotations

import os
import shutil
import stat
import subprocess
import zipfile
from pathlib import Path

import pytest
import yaml
from packs_helpers import CELEGANS_PACK, ROOT, TEST_PACK, copy_pack, pack_env  # noqa: F401

from germquant import packs
from germquant.packs import format as F
from germquant.packs.install import PackInstallError


def _zip_of(folder: Path, dest: Path, *, extra: dict | None = None, prefix: str = "") -> Path:
    with zipfile.ZipFile(dest, "w") as zf:
        for rel in F.pack_files(folder) + ["SHA256SUMS"]:
            zf.write(folder / rel, prefix + rel)
        for name, data in (extra or {}).items():
            if isinstance(data, zipfile.ZipInfo):
                zf.writestr(data, b"target")
            else:
                info = zipfile.ZipInfo("placeholder")
                info.filename = name                  # the raw name: ZipInfo() would turn a backslash into '/'
                zf.writestr(info, data)
    return dest


def test_install_from_folder_and_zip_is_immutable(pack_env, monkeypatch, tmp_path):  # noqa: F811
    target = tmp_path / "installed"
    monkeypatch.setenv("GERMQUANT_PACKS", str(target))
    p = packs.install_pack(TEST_PACK)
    assert p.folder == target / "test_pack" / "0.1.0" and F.verify_sums(p.folder) == []
    assert packs.find_pack("test_pack@0.1.0").folder == p.folder
    again = packs.install_pack(_zip_of(TEST_PACK, tmp_path / "tp.zip", prefix="test_pack/"))   # same content: no-op
    assert again.folder == p.folder
    changed = copy_pack(TEST_PACK, tmp_path / "changed" / "test_pack")
    (changed / "README.md").write_bytes(b"other text\n")
    packs.write_sums(changed)
    with pytest.raises(PackInstallError, match="immutable"):
        packs.install_pack(changed)
    replaced = packs.install_pack(changed, replace=True)
    assert (replaced.folder / "README.md").read_bytes() == b"other text\n"


@pytest.mark.parametrize("entry, needle", [
    ("../evil.txt", "zip-slip"),
    ("/etc/evil.txt", "absolute"),
    ("C:/evil.txt", "drive letter"),
    ("C:evil.txt", "drive letter"),
    # Python's zipfile turns a backslash into '/' when it reads on Windows (then the '..' rule refuses it);
    # elsewhere the backslash rule does
    ("sub\\..\\..\\evil.txt", "backslash|zip-slip"),
])
def test_install_refuses_unsafe_archive_entries(pack_env, tmp_path, entry, needle):  # noqa: F811
    z = _zip_of(TEST_PACK, tmp_path / "bad.zip", extra={entry: b"x"})
    with pytest.raises(PackInstallError, match=needle):
        packs.install_pack(z, target=tmp_path / "t")
    assert not (tmp_path / "evil.txt").exists() and not (tmp_path / "t" / "test_pack").exists()


def test_install_refuses_symlink_entries(pack_env, tmp_path):  # noqa: F811
    info = zipfile.ZipInfo("link_to_somewhere")
    info.external_attr = (stat.S_IFLNK | 0o777) << 16
    z = _zip_of(TEST_PACK, tmp_path / "link.zip", extra={"link_to_somewhere": info})
    with pytest.raises(PackInstallError, match="symlink"):
        packs.install_pack(z, target=tmp_path / "t")


def _make_dir_link(link: Path, target: Path) -> None:
    """A folder link that needs no privilege: a junction on Windows (cmd mklink /J), a symlink elsewhere."""
    if os.name == "nt":
        r = subprocess.run(["cmd", "/c", "mklink", "/J", str(link), str(target)], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
    else:
        os.symlink(target, link, target_is_directory=True)


def _remove_dir_link(link: Path) -> None:
    if os.name == "nt":
        os.rmdir(link)                               # removes the junction itself, never its target's files
    else:
        link.unlink()


def test_links_in_a_pack_folder_are_never_followed_and_install_refuses_them(pack_env, tmp_path):  # noqa: F811
    folder = copy_pack(TEST_PACK, tmp_path / "src" / "test_pack")
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "secret.txt").write_bytes(b"not the pack's\n")
    link = folder / "linked"
    _make_dir_link(link, outside)                    # Python 3.11 sees a junction as neither symlink nor link
    try:
        assert F.is_link(link) and F.pack_links(folder) == ["linked"]
        assert "secret" not in F.build_sums_text(folder)                    # never walked into, never hashed
        assert any(p.startswith("linked: a symlink or junction") for p in F.verify_sums(folder))
        with pytest.raises(PackInstallError, match="a symlink or junction; refused"):
            packs.install_pack(folder, target=tmp_path / "t")
        assert not (tmp_path / "t" / "test_pack").exists()
    finally:
        _remove_dir_link(link)
    assert (outside / "secret.txt").is_file() and F.verify_sums(folder) == []
    try:                                             # a symlinked file, where this account may create one
        os.symlink(folder / "README.md", folder / "link.md")
    except (OSError, NotImplementedError):
        return                                       # Windows without developer mode: the junction case ran
    assert F.pack_links(folder) == ["link.md"]
    with pytest.raises(PackInstallError, match="symlink"):
        packs.install_pack(folder, target=tmp_path / "t2")


def test_zip_install_leaves_junk_behind(pack_env, tmp_path):  # noqa: F811
    z = _zip_of(TEST_PACK, tmp_path / "junk.zip", extra={
        "__pycache__/": b"", "__pycache__/evil.cpython-311.pyc": b"junk", ".DS_Store": b"junk",
        "config/x.pyc": b"junk", "__MACOSX/._pack.yaml": b"junk"})
    p = packs.install_pack(z, target=tmp_path / "t")
    left = sorted(q.relative_to(p.folder).as_posix() for q in p.folder.rglob("*"))
    assert not [n for n in left if n.endswith((".pyc", ".DS_Store")) or "__MACOSX" in n or "__pycache__" in n], left
    assert F.verify_sums(p.folder) == [] and "pack.yaml" in left


def test_install_refuses_non_ascii_hashed_text_anywhere_in_the_pack(pack_env, tmp_path):  # noqa: F811
    from germquant.config import load_config

    accent = ("# caf" + chr(0xE9) + "\n").encode("utf-8")                       # e-acute, two UTF-8 bytes
    # a YAML file outside config/
    a = copy_pack(TEST_PACK, tmp_path / "a" / "test_pack", rebuild=False)
    (a / "extra").mkdir()
    (a / "extra" / "notes.yaml").write_bytes(b"note: 1\n" + accent)
    packs.write_sums(a)
    with pytest.raises(PackInstallError, match=r"extra/notes\.yaml: byte \d+ is not ASCII"):
        packs.install_pack(a, target=tmp_path / "t")
    # a recipe's channel map under another suffix, outside config/
    b = copy_pack(TEST_PACK, tmp_path / "b" / "test_pack", rebuild=False)
    (b / "maps").mkdir()
    (b / "maps" / "alt_map.txt").write_bytes((b / "config" / "channel_maps" / "test_dna_foci.yaml").read_bytes() + accent)
    (b / "config" / "profiles" / "alt.yaml").write_bytes(b"extends: basic.yaml\nio:\n  channel_map: maps/alt_map.txt\n")
    data = yaml.safe_load((b / "pack.yaml").read_text(encoding="ascii"))
    data["recipes"].append(dict(data["recipes"][0], id="alt", default=False, profile="config/profiles/alt.yaml",
                                expected_config_hash=load_config(b / "config" / "profiles" / "alt.yaml").hash))
    (b / "pack.yaml").write_bytes(yaml.safe_dump(data, sort_keys=False).encode("ascii"))
    packs.write_sums(b)
    with pytest.raises(PackInstallError, match=r"maps/alt_map\.txt: byte \d+ is not ASCII"):
        packs.install_pack(b, target=tmp_path / "t")
    # a recipe outside config/profiles/ is refused by the manifest check (its relative paths would not anchor)
    c = copy_pack(TEST_PACK, tmp_path / "c" / "test_pack", edit=lambda d: d["recipes"].append(
        dict(d["recipes"][0], id="alt", default=False, profile="recipes/alt.yaml")))
    with pytest.raises(PackInstallError, match=r"must be config/profiles/<recipe>\.yaml"):
        packs.install_pack(c, target=tmp_path / "t")
    assert not (tmp_path / "t" / "test_pack").exists()


def test_install_refuses_faulty_tunables(pack_env, tmp_path):  # noqa: F811
    for key, needle in (("spots.no_such_key", "is not set by the recipe"), ("spots.log", "is not a scalar")):
        bad = copy_pack(TEST_PACK, tmp_path / key.replace(".", "_") / "test_pack",
                        edit=lambda d, key=key: d["recipes"][0].update(tunables=[{"key": key, "label": "x"}]))
        assert any(needle in p for p in packs.validate_pack(bad))
        with pytest.raises(PackInstallError, match=needle):
            packs.install_pack(bad, target=tmp_path / "t")
    assert packs.validate_pack(TEST_PACK) == [] and packs.validate_pack(CELEGANS_PACK) == []
    assert not (tmp_path / "t" / "test_pack").exists()


def test_reinstall_over_a_damaged_copy_is_refused(pack_env, tmp_path):  # noqa: F811
    target = tmp_path / "installed"
    p = packs.install_pack(TEST_PACK, target=target)
    (p.folder / "README.md").write_bytes(b"edited in place\n")              # the same SHA256SUMS, other files
    with pytest.raises(PackInstallError, match="no longer matches its SHA256SUMS"):
        packs.install_pack(TEST_PACK, target=target)                       # not a silent "already installed"
    fixed = packs.install_pack(TEST_PACK, target=target, replace=True)
    assert F.verify_sums(fixed.folder) == []


def test_install_refuses_bad_sums_non_ascii_and_builtin_collisions(pack_env, tmp_path):  # noqa: F811
    bad = copy_pack(TEST_PACK, tmp_path / "a" / "test_pack")
    (bad / "README.md").write_bytes(b"changed after the sums\n")
    with pytest.raises(PackInstallError, match="SHA256SUMS does not match"):
        packs.install_pack(bad, target=tmp_path / "t")
    nosums = copy_pack(TEST_PACK, tmp_path / "b" / "test_pack")
    (nosums / "SHA256SUMS").unlink()
    with pytest.raises(PackInstallError, match="no SHA256SUMS"):
        packs.install_pack(nosums, target=tmp_path / "t")
    accent = copy_pack(TEST_PACK, tmp_path / "c" / "test_pack", rebuild=False)
    rec = accent / "config" / "profiles" / "basic.yaml"
    rec.write_bytes(rec.read_bytes() + "# caf\u00e9\n".encode("utf-8"))
    packs.write_sums(accent)
    with pytest.raises(PackInstallError, match="not ASCII"):
        packs.install_pack(accent, target=tmp_path / "t")
    with pytest.raises(PackInstallError, match="built-in pack"):
        packs.install_pack(CELEGANS_PACK, target=tmp_path / "t")
    assert not (tmp_path / "t" / "celegans_germline").exists()


def _to_eol(folder: Path, eol: bytes) -> None:
    for rel in F.pack_files(folder):
        if F.is_text(rel):
            p = folder / rel
            p.write_bytes(p.read_bytes().replace(b"\r\n", b"\n").replace(b"\n", eol))


def test_celegans_sums_verify_on_crlf_and_lf_checkouts(tmp_path):
    crlf = copy_pack(CELEGANS_PACK, tmp_path / "win" / "celegans_germline", rebuild=False)
    lf = copy_pack(CELEGANS_PACK, tmp_path / "linux" / "celegans_germline", rebuild=False)
    _to_eol(crlf, b"\r\n")
    _to_eol(lf, b"\n")
    assert F.verify_sums(crlf) == [] and F.verify_sums(lf) == []
    assert F.sums_sha256(crlf) == F.sums_sha256(lf)
    # the engine files it records, LF-normalised, match this checkout whatever its line endings
    assert packs.engine_files_problems(packs.find_pack("celegans_germline")) == []


def _git(*args, cwd, check=True):
    return subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True, check=check)


def _in_git_checkout() -> bool:
    if shutil.which("git") is None:
        return False
    r = subprocess.run(["git", "rev-parse", "--is-inside-work-tree"], cwd=ROOT, capture_output=True, text=True)
    return r.returncode == 0 and r.stdout.strip() == "true"


@pytest.mark.skipif(not _in_git_checkout(), reason="not a git checkout")
def test_every_path_in_committed_sums_is_tracked():
    sums_files = sorted(ROOT.glob("packs/*/SHA256SUMS")) + sorted(ROOT.glob("tests/fixtures/packs/*/SHA256SUMS"))
    assert sums_files, "no SHA256SUMS found"
    for s in sums_files:
        rel_pack = s.parent.relative_to(ROOT).as_posix()
        paths = [f"{rel_pack}/{rel}" for rel in F.read_sums(s.parent)] + [f"{rel_pack}/SHA256SUMS"]
        r = subprocess.run(["git", "ls-files", "--error-unmatch", "--", *paths], cwd=ROOT, capture_output=True,
                           text=True)
        assert r.returncode == 0, f"{rel_pack}: paths in SHA256SUMS are not tracked by git: {r.stderr.strip()}"


@pytest.mark.skipif(shutil.which("git") is None, reason="git not installed")
def test_git_rules_keep_pack_files_byte_identical(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    _git("init", "-q", cwd=repo)
    _git("config", "core.autocrlf", "true", cwd=repo)
    _git("config", "user.email", "t@example.invalid", cwd=repo)
    _git("config", "user.name", "t", cwd=repo)
    shutil.copyfile(ROOT / ".gitattributes", repo / ".gitattributes")
    shutil.copyfile(ROOT / ".gitignore", repo / ".gitignore")
    png = b"\x89PNG\r\n\x1a\n\r\n\r\nzz"                                  # CRLF pairs a text rule would mangle
    files = {"packs/y/img.png": png, "packs/y/models/a.card.json": b'{"name": "a"}\n',
             "packs/y/example/e.ome.tif": b"II*\x00\r\n\r\n", "packs/y/pack.yaml": b"id: y\r\n",
             "packs/y/.DS_Store": b"junk", "packs/y/__pycache__/j.pyc": b"junk", "packs/y/k.pyc": b"junk",
             "tests/fixtures/packs/z/example/e.ome.tif": b"II*\x00\r\n",
             "tests/fixtures/packs/z/.DS_Store": b"junk"}
    for rel, data in files.items():
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        (repo / rel).write_bytes(data)
    _git("add", "-A", cwd=repo)
    tracked = set(_git("ls-files", cwd=repo).stdout.split())
    assert {"packs/y/img.png", "packs/y/models/a.card.json", "packs/y/example/e.ome.tif",
            "tests/fixtures/packs/z/example/e.ome.tif"} <= tracked
    assert not any(t.endswith((".DS_Store", ".pyc")) for t in tracked)
    _git("commit", "-q", "-m", "x", cwd=repo)
    clone = tmp_path / "clone"
    _git("-c", "core.autocrlf=true", "clone", "-q", str(repo), str(clone), cwd=tmp_path)
    assert (clone / "packs/y/img.png").read_bytes() == png
    assert (clone / "packs/y/example/e.ome.tif").read_bytes() == files["packs/y/example/e.ome.tif"]
    assert (clone / "packs/y/pack.yaml").read_bytes() == b"id: y\n"        # text: LF in every checkout


def test_ds_store_is_not_listed_in_sums(tmp_path):
    p = copy_pack(TEST_PACK, tmp_path / "test_pack", rebuild=False)
    (p / ".DS_Store").write_bytes(b"junk")
    (p / "models" / ".DS_Store").write_bytes(b"junk")
    text = F.build_sums_text(p)
    assert ".DS_Store" not in text and F.verify_sums(p) == []


def test_vendor_copies_a_version_that_verifies(pack_env, monkeypatch, tmp_path):  # noqa: F811
    monkeypatch.setenv("GERMQUANT_PACKS", str(TEST_PACK.parent))
    packs.refresh()
    proj = tmp_path / "project"
    out = packs.vendor("test_pack", "0.1.0", proj / "packs")
    assert out == proj / "packs" / "test_pack" / "0.1.0" and F.verify_sums(out) == []
    # a root: checkout pack copies its own files only; its recipes still resolve to the engine files
    cel = packs.vendor("celegans_germline", "0.1.0", proj / "packs")
    assert F.verify_sums(cel) == [] and not (cel / "config").exists()
    with packs.use_project(proj):
        disc = packs.discover(refresh_cache=True)
        assert disc.get("celegans_germline")[0].source == "builtin" and disc.shadowed == []   # identical: quiet
        assert disc.get("test_pack")[0].source == "project"
        vend = packs.load_pack(cel)
        assert packs.recipe_profile_path(vend, "default") == ROOT.resolve() / "config" / "config.yaml"
