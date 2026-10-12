"""`germquant pack install` and packs.vendor (docs/DESIGN_system_packs.md, sections 3.7 and 3.8).

Install takes a folder or a zip and writes ``<packs root>/<id>/<version>/``. It refuses absolute paths, drive
letters, ``..`` components and backslashes in archive names (zip-slip), symlinks and Windows junctions, a
pack that packs.validate_pack faults (the manifest, and tunables the recipe does not set), a missing or wrong
SHA256SUMS (a listed file missing or changed, an unlisted file present), non-ASCII hashed files (every YAML
file of the pack, and each recipe's extends chain and channel map inside the pack, which the config hash
decodes as cp1252, rule R9) and a collision with a built-in id and version (a built-in cannot be shadowed).
Junk (caches, compiled Python, .DS_Store, a zip's __MACOSX folder) is never copied. An installed version is
immutable: the same version with a different SHA256SUMS is refused unless --replace (pack developers only);
the same content installed again is a no-op once the installed copy still verifies. Nothing from a pack is
ever executed.

vendor copies one pack version through the same checks into a folder of the caller's choice (a project's
packs/), so a project keeps the exact version it was made with when the engine moves on (section 3.8).
"""
from __future__ import annotations

import logging
import os
import re
import shutil
import stat
import tempfile
import zipfile
from pathlib import Path, PurePosixPath

from .format import (
    MANIFEST,
    SKIP_DIRS,
    Pack,
    PackFormatError,
    is_link,
    load_pack_folder,
    pack_files,
    read_sums,
    skipped,
    sums_sha256,
    verify_sums,
)

log = logging.getLogger(__name__)

_DRIVE_RE = re.compile(r"^[A-Za-z]:")
_MACOS_FORK_DIR = "__MACOSX"           # resource forks a macOS zip adds next to the files


class PackInstallError(ValueError):
    """An archive or folder that may not be installed (the message says why)."""


def _check_member_name(name: str) -> PurePosixPath:
    if not name or "\x00" in name:
        raise PackInstallError(f"archive entry {name!r}: empty name")
    if "\\" in name:
        raise PackInstallError(f"archive entry {name!r}: backslash in the name (a path trick); refused")
    if _DRIVE_RE.match(name):
        raise PackInstallError(f"archive entry {name!r}: a drive letter; refused")
    pp = PurePosixPath(name)
    if pp.is_absolute() or name.startswith("/"):
        raise PackInstallError(f"archive entry {name!r}: an absolute path; refused")
    if ".." in pp.parts:
        raise PackInstallError(f"archive entry {name!r}: a '..' component (zip-slip); refused")
    return pp


def _is_symlink_entry(info: zipfile.ZipInfo) -> bool:
    mode = (info.external_attr >> 16) & 0xFFFF
    return stat.S_ISLNK(mode)


def _junk_entry(pp: PurePosixPath, is_dir: bool) -> bool:
    """An archive entry that is never extracted: junk (format.skipped) or anything in a __MACOSX folder."""
    if _MACOS_FORK_DIR in pp.parts:
        return True
    if is_dir:
        return any(part in SKIP_DIRS for part in pp.parts)
    return skipped(pp.parts)


def extract_zip(src: Path, dest: Path) -> None:
    """Extract a pack archive into ``dest`` after checking every entry name; refuses symlink entries and leaves
    junk (caches, compiled Python, .DS_Store, __MACOSX) behind, as a folder install does."""
    with zipfile.ZipFile(src) as zf:
        infos = zf.infolist()
        for info in infos:                       # check everything before writing anything
            _check_member_name(info.filename)
            if _is_symlink_entry(info):
                raise PackInstallError(f"archive entry {info.filename!r} is a symlink; refused")
        root = dest.resolve()
        for info in infos:
            pp = PurePosixPath(info.filename)
            if _junk_entry(pp, info.is_dir()):
                continue
            target = (dest / Path(*pp.parts)).resolve()
            if root != target and root not in target.parents:
                raise PackInstallError(f"archive entry {info.filename!r} would land outside the pack folder")
            if info.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with zf.open(info) as fin, open(target, "wb") as fout:
                shutil.copyfileobj(fin, fout)


def _copy_folder(src: Path, dest: Path) -> None:
    """Copy a pack folder, refusing symlinks and Windows junctions anywhere in it (format.is_link; os.walk
    would walk into a junction) and leaving junk (caches, .DS_Store) behind."""
    for dirpath, dirnames, filenames in os.walk(src, followlinks=False):
        d = Path(dirpath)
        for name in list(dirnames):
            if is_link(d / name):
                raise PackInstallError(f"{d / name}: a symlink or junction; refused")
        for name in filenames:
            f = d / name
            if is_link(f):
                raise PackInstallError(f"{f}: a symlink or junction; refused")
            rel = f.relative_to(src)
            if skipped(rel.parts):
                continue
            out = dest / rel
            out.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(f, out)


def _pack_folder_in(tree: Path) -> Path:
    """The pack root inside an extracted tree: the tree itself, or its single top-level folder."""
    if (tree / MANIFEST).is_file():
        return tree
    subs = [p for p in tree.iterdir() if p.is_dir() and not p.name.startswith(".") and p.name != "__MACOSX"]
    if len(subs) == 1 and (subs[0] / MANIFEST).is_file():
        return subs[0]
    raise PackInstallError(f"no {MANIFEST} at the top of the archive or folder")


def _hashed_text_files(folder: Path, pack: Pack) -> list[Path]:
    """The files of the pack that a config hash reads (rule R9): every .yaml or .yml file it ships, and each
    loadable recipe's extends-chain files and channel map that lie inside the pack, whatever their suffix."""
    from .discovery import recipe_profile_path
    from .lint import _chain_files, _inside
    from .validate import load_recipe_config

    root = folder.resolve()
    files = {(root / rel) for rel in pack_files(root) if PurePosixPath(rel).suffix.lower() in (".yaml", ".yml")}
    for r in pack.recipes:
        try:
            chain = _chain_files(recipe_profile_path(pack, r))
            cfg = load_recipe_config(pack, r)
        except Exception:  # noqa: BLE001 - a recipe that does not load is pack status's and lint's to report
            continue
        cm = Path(str(cfg.get("io.channel_map")))
        for f in [*chain, cm if cm.is_absolute() else cfg.base_dir / cm]:
            if f.is_file() and _inside(f, root):
                files.add(f.resolve())
    return sorted(files)


def check_installable(folder: Path) -> Pack:
    """The Pack in ``folder`` if it may be installed; PackInstallError otherwise."""
    from .validate import validate_pack

    try:
        pack = load_pack_folder(folder, source="path")
    except PackFormatError as e:
        raise PackInstallError(str(e)) from e
    if pack.problems:
        raise PackInstallError("the manifest has problems: " + "; ".join(pack.problems))
    try:
        read_sums(folder)
    except PackFormatError as e:
        raise PackInstallError(f"{e} (a pack ships SHA256SUMS; `germquant pack build-sums` writes it)") from e
    bad = verify_sums(folder)
    if bad:
        raise PackInstallError("SHA256SUMS does not match: " + "; ".join(bad[:5]) + (" ..." if len(bad) > 5 else ""))
    for p in _hashed_text_files(folder, pack):
        try:
            p.read_bytes().decode("ascii")
        except UnicodeDecodeError as e:
            raise PackInstallError(f"{p.relative_to(folder.resolve()).as_posix()}: byte {e.start} is not ASCII; "
                                   f"hashed recipe files must be ASCII (rule R9)") from e
    faults = validate_pack(pack)
    if faults:
        raise PackInstallError("the pack has faults: " + "; ".join(faults))
    return pack


def install_pack(src: str | Path, *, replace: bool = False, target: str | Path | None = None,
                 _vendor: bool = False) -> Pack:
    """Install a pack folder or zip into ``target`` (default: the first $GERMQUANT_PACKS folder, else the
    per-user folder) as ``<target>/<id>/<version>/``. Returns the installed Pack. (``_vendor``: vendor's
    copy of a built-in into a project, which discovery treats as identical, not as shadowing.)"""
    from .discovery import ROOT_BUILTIN, discover, install_root, refresh

    src = Path(src)
    if not src.exists():
        raise PackInstallError(f"{src}: no such folder or archive")
    root = Path(target) if target is not None else install_root()
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".install-", dir=root) as tmp:
        tree = Path(tmp) / "src"
        tree.mkdir()
        if src.is_dir():
            if is_link(src):
                raise PackInstallError(f"{src}: a symlink or junction; refused")
            _copy_folder(src, tree)
        elif zipfile.is_zipfile(src):
            extract_zip(src, tree)
        else:
            raise PackInstallError(f"{src}: neither a folder nor a zip archive")
        folder = _pack_folder_in(tree)
        pack = check_installable(folder)
        builtin = [p for p in discover(refresh_cache=True).packs
                   if p.source == ROOT_BUILTIN and p.id == pack.id and p.version == pack.version]
        if builtin and not _vendor:
            raise PackInstallError(f"{pack.ref} is a built-in pack ({builtin[0].folder}); a built-in id and version "
                                   f"cannot be shadowed: give the pack a new version")
        dest = root / pack.id / pack.version
        if dest.exists():
            same = sums_sha256(dest) is not None and sums_sha256(dest) == sums_sha256(folder)
            damaged = verify_sums(dest) if same else []
            if same and not damaged:
                log.info("%s is already installed at %s", pack.ref, dest)
                refresh()
                return load_pack_folder(dest, source="installed")
            if not replace and damaged:
                raise PackInstallError(f"{pack.ref} is installed at {dest}, but that copy no longer matches its "
                                       f"SHA256SUMS ({damaged[0]}); an installed version is immutable: pass --replace "
                                       f"to install it again")
            if not replace:
                raise PackInstallError(f"{pack.ref} is already installed at {dest} with other content; an installed "
                                       f"version is immutable (give the pack a new version, or pass --replace)")
            log.warning("replacing the installed %s at %s (--replace)", pack.ref, dest)
            trash = Path(tmp) / "old"
            os.replace(dest, trash)
        dest.parent.mkdir(parents=True, exist_ok=True)
        staged = Path(tmp) / "staged"
        shutil.copytree(folder, staged)
        os.replace(staged, dest)
    refresh()
    return load_pack_folder(dest, source="installed")


def vendor(pack_id: str, version: str, dest: str | Path) -> Path:
    """Copy pack ``pack_id`` at ``version`` through the install checks into ``<dest>/<id>/<version>/`` and
    return that folder. A root: checkout pack copies only its own files; its recipes keep pointing at the
    frozen engine files, whose recorded digests are checked here. An identical copy already there is kept."""
    from .discovery import find_pack
    from .status import engine_files_problems

    pack = find_pack(f"{pack_id}@{version}")
    bad = engine_files_problems(pack)
    if bad:
        raise PackInstallError(f"{pack.ref}: " + "; ".join(bad))
    installed = install_pack(pack.folder, target=dest, _vendor=True)
    return installed.folder
