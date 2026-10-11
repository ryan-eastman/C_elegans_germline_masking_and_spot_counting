"""Where packs live and how a name finds one (docs/DESIGN_system_packs.md, section 3.8).

packs.pack_roots() searches, in order:

1. ``<project>/packs/`` while a project is open or running (vendored and derived packs; ``use_project``);
2. the built-in packs, ``<root>/packs/<id>/pack.yaml``, for each built-in root (config.builtin_roots: the
   engine checkout that holds src/germquant, then ``$GERMQUANT_CONFIG_ROOT``; never the working directory);
3. each folder of ``$GERMQUANT_PACKS`` (an os.pathsep list; the first folder is the install target);
4. the per-user folder (``%LOCALAPPDATA%/germquant/packs``, ``~/Library/Application Support/germquant/packs``
   or ``~/.local/share/germquant/packs``).

Under 1, 3 and 4 installed packs sit at ``<root>/<id>/<version>/pack.yaml``; a folder of 3 or 4 may also hold
a pack under development at ``<root>/<id>/pack.yaml``. The first (id, version) found wins, except that a
built-in id and version cannot be shadowed: a different copy elsewhere is listed as shadowed and ignored, and
an identical copy (the same pack_sha256, such as a project's vendored copy) is not reported at all. A pack
whose manifest is broken is logged and listed as unusable; it never breaks a run that does not use it.
Discovery is cached per process and per set of roots; refresh() rescans.
"""
from __future__ import annotations

import contextlib
import contextvars
import logging
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

from .format import MANIFEST, Pack, PackFormatError, load_pack_folder, sums_sha256

log = logging.getLogger(__name__)

ROOT_PROJECT, ROOT_BUILTIN, ROOT_ENV, ROOT_USER, ROOT_PATH = "project", "builtin", "GERMQUANT_PACKS", "user", "path"

_ACTIVE_PROJECT: contextvars.ContextVar[Path | None] = contextvars.ContextVar("germquant_active_project",
                                                                              default=None)


@contextlib.contextmanager
def use_project(project: str | Path | None):
    """While open, ``<project>/packs/`` is searched first (package P5 opens projects with it)."""
    token = _ACTIVE_PROJECT.set(Path(project) if project is not None else None)
    try:
        yield
    finally:
        _ACTIVE_PROJECT.reset(token)


def user_packs_dir() -> Path:
    """The per-user install folder of this platform."""
    if os.name == "nt":
        base = os.environ.get("LOCALAPPDATA") or str(Path.home() / "AppData" / "Local")
        return Path(base) / "germquant" / "packs"
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "germquant" / "packs"
    base = os.environ.get("XDG_DATA_HOME") or str(Path.home() / ".local" / "share")
    return Path(base) / "germquant" / "packs"


def env_pack_dirs() -> list[Path]:
    raw = os.environ.get("GERMQUANT_PACKS", "")
    return [Path(p) for p in raw.split(os.pathsep) if p.strip()]


def install_root() -> Path:
    """Where `germquant pack install` writes: the first folder of $GERMQUANT_PACKS, else the per-user folder."""
    env = env_pack_dirs()
    return env[0] if env else user_packs_dir()


def pack_roots(project: str | Path | None = None) -> list[tuple[str, Path]]:
    """``[(kind, folder)]`` in search order (module docstring)."""
    from ..config import builtin_roots

    roots: list[tuple[str, Path]] = []
    proj = Path(project) if project is not None else _ACTIVE_PROJECT.get()
    if proj is not None:
        roots.append((ROOT_PROJECT, proj / "packs"))
    for r in builtin_roots():
        roots.append((ROOT_BUILTIN, r / "packs"))
    for r in env_pack_dirs():
        roots.append((ROOT_ENV, r))
    roots.append((ROOT_USER, user_packs_dir()))
    seen, out = set(), []
    for kind, folder in roots:
        key = str(Path(folder).resolve()) if Path(folder).exists() else str(folder)
        if key in seen:
            continue
        seen.add(key)
        out.append((kind, Path(folder)))
    return out


def _candidate_folders(kind: str, root: Path) -> list[Path]:
    """Pack folders under one root, by the layout of its kind."""
    if not root.is_dir():
        return []
    out = []
    try:
        children = sorted(p for p in root.iterdir() if p.is_dir() and not p.name.startswith("."))
    except OSError as e:
        log.warning("cannot list the pack folder %s: %s", root, e)
        return []
    for child in children:
        if kind == ROOT_BUILTIN:
            if (child / MANIFEST).is_file():
                out.append(child)
            continue
        if kind != ROOT_PROJECT and (child / MANIFEST).is_file():
            out.append(child)                      # a pack under development: <root>/<id>/pack.yaml
        try:
            versions = sorted(p for p in child.iterdir() if p.is_dir() and (p / MANIFEST).is_file())
        except OSError:
            versions = []
        out.extend(versions)                       # installed: <root>/<id>/<version>/pack.yaml
    return out


@dataclass
class Discovery:
    packs: list[Pack] = field(default_factory=list)        # usable manifests, winners only, in search order
    unusable: list[dict] = field(default_factory=list)     # {"folder", "source", "id", "version", "problems"}
    shadowed: list[Pack] = field(default_factory=list)
    roots: list[tuple[str, Path]] = field(default_factory=list)

    @property
    def builtin_found(self) -> bool:
        return any(kind == ROOT_BUILTIN for kind, _ in self.roots)

    def get(self, pack_id: str, version: str | None = None) -> list[Pack]:
        return [p for p in self.packs if p.id == pack_id and (version is None or p.version == version)]


_CACHE: dict[tuple, Discovery] = {}


def refresh() -> None:
    """Forget the cached scan; the next lookup rescans every root."""
    _CACHE.clear()


def discover(project: str | Path | None = None, *, refresh_cache: bool = False) -> Discovery:
    roots = pack_roots(project)
    key = tuple((k, str(r)) for k, r in roots)
    if not refresh_cache and key in _CACHE:
        return _CACHE[key]
    found: list[Pack] = []
    unusable: list[dict] = []
    for kind, root in roots:
        for folder in _candidate_folders(kind, root):
            try:
                pack = load_pack_folder(folder, source=kind)
            except PackFormatError as e:
                log.warning("pack folder %s is unusable: %s", folder, e)
                installed = folder.parent != root                 # <root>/<id>/<version>/: the names say which
                unusable.append({"folder": str(folder), "source": kind,
                                 "id": folder.parent.name if installed else folder.name,
                                 "version": folder.name if installed else None, "problems": [str(e)]})
                continue
            if pack.problems:
                log.warning("pack %s at %s is unusable: %s", pack.ref, folder, "; ".join(pack.problems))
                unusable.append({"folder": str(folder), "source": kind, "id": pack.id, "version": pack.version,
                                 "problems": list(pack.problems)})
                continue
            found.append(pack)
    winners: dict[tuple[str, str], Pack] = {}
    order: list[tuple[str, str]] = []
    shadowed: list[Pack] = []
    for pack in found:                              # built-ins first: they cannot be shadowed
        key2 = (pack.id, pack.version)
        if pack.source == ROOT_BUILTIN and key2 not in winners:
            winners[key2] = pack
    for pack in found:
        key2 = (pack.id, pack.version)
        if key2 not in winners:
            winners[key2] = pack
        if key2 not in order:
            order.append(key2)
        win = winners[key2]
        if win is pack:
            continue
        if sums_sha256(win.folder) and sums_sha256(win.folder) == sums_sha256(pack.folder):
            continue                                # an identical copy: not shadowed
        shadowed.append(pack)
    disc = Discovery(packs=[winners[k] for k in order], unusable=unusable, shadowed=shadowed, roots=roots)
    _CACHE[key] = disc
    return disc


def parse_ref(ref: str) -> tuple[str, str | None, str | None]:
    """``id``, ``id@version``, ``id/recipe`` or ``id@version/recipe`` -> (id, version or None, recipe or None)."""
    text = str(ref).strip()
    recipe = None
    if "/" in text:
        text, recipe = text.split("/", 1)
        recipe = recipe.strip() or None
    version = None
    if "@" in text:
        text, version = text.split("@", 1)
        version = version.strip() or None
    return text.strip(), version, recipe


def _version_key(v: str):
    try:
        from packaging.version import Version

        return (0, Version(v))
    except Exception:  # noqa: BLE001 - a non-PEP 440 semver tail sorts by its numbers
        nums = [int(x) if x.isdigit() else 0 for x in v.replace("-", ".").split(".")[:3]]
        return (1, tuple(nums))


def engine_compatible(pack: Pack) -> bool:
    """Whether this engine's version lies inside the pack's requires.germquant range (True without one)."""
    spec = (pack.data.get("requires") or {}).get("germquant")
    if not spec:
        return True
    try:
        from packaging.specifiers import SpecifierSet

        from .. import __version__

        return SpecifierSet(str(spec)).contains(__version__, prereleases=True)
    except Exception:  # noqa: BLE001 - an unparsable range is a manifest problem reported elsewhere
        return False


def available_refs(disc: Discovery | None = None) -> list[str]:
    disc = disc or discover()
    return sorted(p.ref for p in disc.packs)


def find_pack(ref: str, project: str | Path | None = None) -> Pack:
    """The pack ``id`` or ``id@version`` names. Without a version: the newest version compatible with this
    engine (interactive use; projects and recipes always pin a version). FileNotFoundError naming the packs
    available; an unusable pack of that name is named with its problems."""
    if isinstance(ref, Pack):
        return ref
    pid, version, _ = parse_ref(ref)
    disc = discover(project)
    hits = disc.get(pid, version)
    if version is None and hits:
        compatible = [p for p in hits if engine_compatible(p)]
        pool = compatible or hits
        return sorted(pool, key=lambda p: _version_key(p.version))[-1]
    if hits:
        return hits[0]
    bad = [u for u in disc.unusable if u["id"] == pid and (version is None or u["version"] == version)]
    why = ("; unusable: " + "; ".join(f"{u['folder']}: {'; '.join(u['problems'])}" for u in bad)) if bad else ""
    roots = ", ".join(f"{k} {r}" for k, r in disc.roots)
    raise FileNotFoundError(f"no pack {ref!r}{why}. Packs available: {', '.join(available_refs(disc)) or 'none'} "
                            f"(searched: {roots})")


def engine_root_for(pack: Pack) -> Path:
    """The engine config root a root: checkout pack's recipe paths are relative to: the built-in root the pack
    lives in, else (a vendored or installed copy) the first built-in root."""
    from ..config import builtin_roots

    if pack.source == ROOT_BUILTIN:
        return pack.folder.parent.parent
    roots = builtin_roots()
    if not roots:
        raise FileNotFoundError(f"pack {pack.ref} points at engine files (root: checkout), but no engine checkout was "
                                f"found: install the engine as an editable checkout or set GERMQUANT_CONFIG_ROOT")
    return roots[0]


def recipe_profile_path(pack: Pack, recipe) -> Path:
    """The profile file of a recipe (a Recipe or a recipe id): inside the pack (root: pack) or under the
    engine config root (root: checkout)."""
    r = pack.recipe(recipe) if not hasattr(recipe, "profile") else recipe
    base = engine_root_for(pack) if pack.root == "checkout" else pack.folder
    return base / r.profile


def resolve_recipe(ref: str, project: str | Path | None = None) -> Path:
    """``<id>@<version>/<recipe>`` (or the unversioned ``<id>/<recipe>``, newest compatible version) -> the
    recipe's profile file, which load_config reads like any other profile. FileNotFoundError naming the
    candidates and the packs available."""
    pid, version, recipe = parse_ref(ref)
    if not pid or recipe is None:
        raise FileNotFoundError(f"{ref!r} is not <id>@<version>/<recipe>")
    pack = find_pack(f"{pid}@{version}" if version else pid, project)   # its error names the packs available
    available = f"Packs available: {', '.join(available_refs(discover(project))) or 'none'}"
    try:
        r = pack.recipe(recipe)
    except FileNotFoundError as e:
        raise FileNotFoundError(f"{e}. {available}") from None
    path = recipe_profile_path(pack, r)
    if not path.is_file():
        raise FileNotFoundError(f"pack {pack.ref} recipe {recipe!r}: profile {path} is missing. {available}")
    return path
