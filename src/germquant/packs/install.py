"""`germquant pack install` and packs.vendor (docs/DESIGN_system_packs.md, sections 3.7 and 3.8).

Install takes a folder or a zip and writes ``<packs root>/<id>/<version>/``. It refuses absolute paths, drive
letters, ``..`` components and backslashes in archive names (zip-slip), symlinks, a missing or wrong
SHA256SUMS (a listed file missing or changed, an unlisted file present), non-ASCII hashed files (the recipe
and channel-map YAML under config/, decoded as cp1252 by the config hash, rule R9) and a collision with a
built-in id and version (a built-in cannot be shadowed). An installed version is immutable: the same version
with a different SHA256SUMS is refused unless --replace (pack developers only); the same content installed
again is a no-op. Nothing from a pack is ever executed.

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

from .format import MANIFEST, Pack, PackFormatError, load_pack_folder, read_sums, skipped, sums_sha256, verify_sums

log = logging.getLogger(__name__)

_DRIVE_RE = re.compile(r"^[A-Za-z]:")


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


def extract_zip(src: Path, dest: Path) -> None:
    """Extract a pack archive into ``dest`` after checking every entry name; refuses symlink entries."""
    with zipfile.ZipFile(src) as zf:
        infos = zf.infolist()
        for info in infos:                       # check everything before writing anything
            _check_member_name(info.filename)
            if _is_symlink_entry(info):
                raise PackInstallError(f"archive entry {info.filename!r} is a symlink; refused")
        root = dest.resolve()
        for info in infos:
            pp = PurePosixPath(info.filename)
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
    """Copy a pack folder, refusing symlinks anywhere in it and leaving junk (caches, .DS_Store) behind."""
    for dirpath, dirnames, filenames in os.walk(src, followlinks=False):
        d = Path(dirpath)
        for name in list(dirnames):
            if (d / name).is_symlink():
                raise PackInstallError(f"{d / name}: a symlink; refused")
        for name in filenames:
            f = d / name
            if f.is_symlink():
                raise PackInstallError(f"{f}: a symlink; refused")
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


def check_installable(folder: Path) -> Pack:
    """The Pack in ``folder`` if it may be installed; PackInstallError otherwise."""
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
    for p in sorted((folder / "config").rglob("*")) if (folder / "config").is_dir() else []:
        if p.is_file() and p.suffix.lower() in (".yaml", ".yml"):
            try:
                p.read_bytes().decode("ascii")
            except UnicodeDecodeError as e:
                raise PackInstallError(f"{p.relative_to(folder).as_posix()}: byte {e.start} is not ASCII; hashed "
                                       f"recipe files must be ASCII (rule R9)") from e
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
            if src.is_symlink():
                raise PackInstallError(f"{src}: a symlink; refused")
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
            if same and not verify_sums(dest):
                log.info("%s is already installed at %s", pack.ref, dest)
                refresh()
                return load_pack_folder(dest, source="installed")
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
