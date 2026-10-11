"""Config loading. Thin wrapper over YAML with attribute access + the channel map."""
from __future__ import annotations

import hashlib
import json
import logging
import os
from pathlib import Path
from typing import Any

import yaml

from .fsutil import read_config_text
from .io.channel_map import ChannelMap

log = logging.getLogger(__name__)


class Config:
    """Dict-backed config with attribute access and a resolved ChannelMap.

    Access nested keys via attributes (``cfg.segmentation.nuclei.method``) or
    ``cfg.get("segmentation.nuclei.method", default)``.

    Two instance attributes are not config data and never enter ``as_dict()``, the manifest's config block
    or ``hash`` (docs/DESIGN_system_packs.md, section 3.4): ``pack_record`` (None, or the pack a run was
    started from: set by ``api.resolve_config(pack=...)``) and ``anchors`` (empty, or the inherited relative
    paths load_config made absolute, each with its original value, and whether the chain went through one
    of the two extends fallbacks; stages.stage_hashes reads it). Both are set here, so ``__getattr__`` is
    never consulted for them.
    """

    def __init__(self, data: dict, *, base_dir: Path, channel_map: ChannelMap, raw_text: str,
                 channel_map_text: str = "", anchors: dict | None = None):
        self._data = data
        self.base_dir = base_dir
        self.channel_map = channel_map
        self._raw_text = raw_text
        self._channel_map_text = channel_map_text
        self._overrides: dict[str, Any] = {}
        self.pack_record: dict | None = None
        self.anchors: dict = dict(anchors) if anchors else {}

    # ---- access helpers ----
    def __getattr__(self, name: str) -> Any:
        try:
            val = self._data[name]
        except KeyError as e:  # pragma: no cover
            raise AttributeError(name) from e
        return _Node(val) if isinstance(val, dict) else val

    def get(self, dotted: str, default: Any = None) -> Any:
        node: Any = self._data
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node

    def set(self, dotted: str, value: Any) -> None:
        """Override a (possibly nested) value at runtime — for CLI flags like ``--no-spots`` that
        flip a config switch. The override shows up in ``as_dict()`` AND changes ``hash``, so a
        ``--no-spots`` run is provenance-distinguishable from a normal one."""
        node = self._data
        parts = dotted.split(".")
        for part in parts[:-1]:
            nxt = node.setdefault(part, {})
            if not isinstance(nxt, dict):
                raise ValueError(f"cannot set {dotted!r}: {part!r} is not a section")
            node = nxt
        node[parts[-1]] = value
        self._overrides[dotted] = value

    @property
    def hash(self) -> str:
        """Stable hash of the resolved config — goes into the provenance manifest.

        Includes the referenced channel-map file (it sets which channel is DAPI/SYP/RAD-51, as
        much an effective parameter as config.yaml) and any runtime overrides (e.g. ``--no-spots``).
        """
        payload = self._raw_text + "\x00" + self._channel_map_text
        if self._overrides:
            payload += "\x00" + json.dumps(self._overrides, sort_keys=True)
        return hashlib.sha256(payload.encode()).hexdigest()[:12]

    def as_dict(self) -> dict:
        return json.loads(json.dumps(self._data))


class _Node:
    def __init__(self, d: dict):
        self._d = d

    def __getattr__(self, name: str) -> Any:
        try:
            val = self._d[name]
        except KeyError as e:
            raise AttributeError(name) from e
        return _Node(val) if isinstance(val, dict) else val

    def get(self, key: str, default: Any = None) -> Any:
        return self._d.get(key, default)

    def __repr__(self) -> str:  # pragma: no cover
        return f"_Node({self._d!r})"


def _base_dir_for(path: Path) -> Path:
    """The directory that relative paths in a config (channel map, model, exclusions file) resolve
    against. Two layouts hop to the repo root: ``config/x.yaml`` (the rule the pipeline always had) and
    ``config/profiles/x.yaml``; any other file resolves against its own directory."""
    p = path.resolve()
    if p.parent.name == "config":
        return p.parent.parent
    if p.parent.name == "profiles" and p.parent.parent.name == "config":
        return p.parent.parent.parent
    return p.parent


def _deep_merge(base: dict, overlay: dict, _path: str = "") -> dict:
    """Overlay wins; mappings merge recursively; lists and scalars replace the base value wholesale
    (so ``io.exclude_patterns`` in a profile replaces the base list). A section header with no keys
    (``coloc:`` followed only by comments) is an authoring mistake and raises rather than erasing the
    base section."""
    out = dict(base)
    for k, v in overlay.items():
        here = f"{_path}.{k}" if _path else k
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v, here)
        elif v is None and isinstance(out.get(k), dict):
            raise ValueError(f"profile section {here!r} is empty (only comments?) and would erase the base "
                             f"section; delete the header or give it keys")
        else:
            out[k] = v
    return out


# relative paths that a profile inherits from its base must keep resolving against the BASE's directory
_PATH_KEYS = (("io", "channel_map"), ("segmentation", "nuclei", "cellpose_model"), ("qc", "exclusions_file"))


def _is_path_like(value: str, base_dir: Path) -> bool:
    """Whether an inherited relative value is a path to anchor. A bare name with no directory part that
    exists nowhere next to the base (a built-in Cellpose model such as ``cpsam``, the generic profile's
    model) is a name, not a path: anchoring it would turn it into ``<repo>/cpsam``, which Cellpose
    cannot load (and ``method: auto`` would silently fall back to the classical watershed)."""
    return len(Path(value).parts) > 1 or (base_dir / value).exists()


def _get_in(d: dict, keys):
    for k in keys:
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def _set_in(d: dict, keys, value):
    for k in keys[:-1]:
        d = d.setdefault(k, {})
    d[keys[-1]] = value


PACK_EXTENDS_PREFIX = "pack:"


def _legacy_extends_candidates(path: Path, base_dir: Path, ext: str) -> list[Path]:
    """The files today's rule tries for ``extends: <ext>``, in order (unchanged since profiles existed)."""
    ext_p = Path(ext)
    if ext_p.is_absolute():
        return [ext_p]
    if len(ext_p.parts) > 1:                         # a relative path: the file's own directory first
        return [path.parent / ext_p, base_dir / "config" / ext_p, base_dir / ext_p]
    # a bare name: the repo config dir first, so a Snakemake profile cannot shadow it
    return [base_dir / "config" / ext_p, base_dir / ext_p, path.parent / ext_p]


def _is_file(p: Path) -> bool:
    """Path.is_file, False for a name the platform cannot even stat (a ``pack:`` value on Windows)."""
    try:
        return p.is_file()
    except (OSError, ValueError):
        return False


def builtin_roots() -> list[Path]:
    """The engine's own roots, searched for the built-in packs (``<root>/packs/<id>/pack.yaml``) and by the
    bare-name extends fallback: the checkout that holds this package (provenance._repo_root, taken only when
    it holds ``src/germquant``), then ``$GERMQUANT_CONFIG_ROOT``. Never the working directory. In a
    non-editable install inside another project's environment, _repo_root finds that project's folder, which
    holds no src/germquant, so it is not taken (docs/DESIGN_system_packs.md, section 3.8)."""
    from .provenance import _repo_root

    roots: list[Path] = []
    repo = _repo_root()
    if repo is not None and (repo / "src" / "germquant").is_dir():
        roots.append(repo.resolve())
    env = os.environ.get("GERMQUANT_CONFIG_ROOT")
    if env:
        p = Path(env).resolve()
        if p not in roots:
            roots.append(p)
    return roots


def _packs_available() -> str:
    """The packs discovery finds, for an error message; never raises (a broken pack must not hide the
    real error)."""
    try:
        from .packs import list_packs

        names = sorted(f"{p['id']}@{p['version']}" for p in list_packs())
        return ", ".join(names) if names else "none"
    except Exception as e:  # noqa: BLE001 - only a message
        return f"unknown ({type(e).__name__}: {e})"


def _resolve_extends_ex(path: Path, base_dir: Path, ext: str) -> tuple[Path, bool]:
    """``(file, new_form)``: the file ``extends: <ext>`` names. Today's candidates are tried first and in
    their old order, so an existing file resolves exactly as before (design rule R4). Only when all of them
    fail are the two new forms tried (``new_form`` True):

    * ``pack:<id>@<version>/<recipe>`` (or the unversioned ``pack:<id>/<recipe>``, for interactive use):
      the recipe file of an installed or built-in pack, through packs.resolve_recipe;
    * a bare name with no directory part: ``<root>/config/<name>``, then ``<root>/config/profiles/<name>``,
      for each built-in root (builtin_roots: the checkout that holds src/germquant, then
      ``$GERMQUANT_CONFIG_ROOT``; never the working directory).

    The error lists every candidate tried and the packs available."""
    cand = _legacy_extends_candidates(path, base_dir, ext)
    for c in cand:
        if _is_file(c) and c.resolve() != path.resolve():
            return c, False
    tried = [str(c) for c in cand]
    if ext.startswith(PACK_EXTENDS_PREFIX):
        from .packs import resolve_recipe    # lazy: packs imports this module

        try:
            found = resolve_recipe(ext[len(PACK_EXTENDS_PREFIX):])
        except FileNotFoundError as e:
            raise FileNotFoundError(f"{path}: extends {ext!r} not found (tried {tried}); {e}") from None
        if found.resolve() == path.resolve():
            raise ValueError(f"{path}: extends {ext!r} names this file itself")
        log.info("%s: extends %r resolved to %s", path, ext, found)
        return found, True
    ext_p = Path(ext)
    if not ext_p.is_absolute() and len(ext_p.parts) == 1:
        for root in builtin_roots():
            for c in (root / "config" / ext_p, root / "config" / "profiles" / ext_p):
                tried.append(str(c))
                if _is_file(c) and c.resolve() != path.resolve():
                    log.info("%s: extends %r resolved to %s (built-in root %s)", path, ext, c, root)
                    return c, True
    raise FileNotFoundError(f"{path}: extends {ext!r} not found (tried {tried}); packs available: "
                            f"{_packs_available()}")


def _resolve_extends(path: Path, base_dir: Path, ext: str) -> Path:
    return _resolve_extends_ex(path, base_dir, ext)[0]


def load_config(path: str | Path, _seen: frozenset[Path] = frozenset()) -> Config:
    """Load a config file. A file may name a base with a top-level ``extends: <path>`` (an absolute
    path, a path relative to the file, or a bare name looked up in the repo's config dir): the base is
    loaded first and the file's own sections are deep-merged over it (a "profile"). Relative paths the
    profile inherits from the base (channel map, Cellpose model, exclusions file) keep resolving against
    the base's directory, so a profile can live anywhere (a bare built-in model name such as ``cpsam``
    is not a path and passes through unchanged). Hash rule: a file WITHOUT ``extends`` hashes
    exactly as it always has (its raw text plus the channel-map text plus runtime overrides), so existing
    runs keep reproducing; a profile hashes the base text plus its own text plus the channel-map text.
    The text is decoded with ``fsutil.CONFIG_TEXT_ENCODING`` on every platform (see the note there), so
    the hash does not depend on the machine's locale.

    Where today's rule finds no base, two more forms resolve (``_resolve_extends_ex``): ``pack:<id>@<v>/<r>``
    and a bare name under the built-in roots. Each inherited path made absolute is recorded in
    ``Config.anchors`` with its original relative value; in a chain that went through one of the two new
    forms the stage hashes use the relative value, so the same edited recipe hashes alike on every machine
    (section 3.4). Every other chain records the same values but keeps today's stage-hash payload."""
    path = Path(path)
    if path.resolve() in _seen:
        raise ValueError(f"{path}: extends cycle")
    raw_text = read_config_text(path)
    data = yaml.safe_load(raw_text) or {}
    base_dir = _base_dir_for(path)
    hash_text = raw_text
    anchored: dict[str, dict] = {}
    portable = False

    ext = data.pop("extends", None)
    if ext:
        ext_path, new_form = _resolve_extends_ex(path, base_dir, str(ext))
        base = load_config(ext_path, _seen | {path.resolve()})
        base_data = base.as_dict()
        if base.base_dir != base_dir:
            for keys in _PATH_KEYS:            # inherited relative paths stay anchored to the base
                v = _get_in(base_data, keys)
                if (isinstance(v, str) and v and not Path(v).is_absolute() and _get_in(data, keys) is None
                        and _is_path_like(v, base.base_dir)):
                    _set_in(base_data, keys, str((base.base_dir / v).resolve()))
                    anchored[".".join(keys)] = {"relative": v, "absolute": _get_in(base_data, keys)}
        # anchors the base inherited itself, for the keys this file leaves alone
        for key, rec in (base.anchors.get("paths") or {}).items():
            if key not in anchored and _get_in(data, key.split(".")) is None:
                anchored[key] = dict(rec)
        portable = new_form or bool(base.anchors.get("portable"))
        data = _deep_merge(base_data, data)
        hash_text = base._raw_text + "\x00extends\x00" + raw_text

    cm_rel = data.get("io", {}).get("channel_map")
    if cm_rel is None:
        raise ValueError("config.io.channel_map is required")
    cm_path = (base_dir / cm_rel) if not Path(cm_rel).is_absolute() else Path(cm_rel)
    channel_map = ChannelMap.from_yaml(cm_path)
    cm_text = read_config_text(cm_path) if cm_path.exists() else ""

    anchors = {"paths": anchored, "portable": portable} if (anchored or portable) else None
    return Config(data, base_dir=base_dir, channel_map=channel_map, raw_text=hash_text,
                  channel_map_text=cm_text, anchors=anchors)


def config_roots() -> list[Path]:
    """Where a ``config/`` folder is looked for, in order: the germquant checkout, the working directory,
    then ``$GERMQUANT_CONFIG_ROOT`` (containers and non-editable installs)."""
    from .provenance import _repo_root

    roots: list[Path | None] = [_repo_root(), Path.cwd()]
    env = os.environ.get("GERMQUANT_CONFIG_ROOT")
    if env:
        roots.append(Path(env))
    return [r for r in roots if r is not None]


def find_profile(name: str) -> Path:
    """``config/profiles/<name>.yaml`` under the first config root that has it. FileNotFoundError naming
    the places tried and the profiles available (the message ``germquant --profile`` prints)."""
    tried = []
    for root in config_roots():
        p = root / "config" / "profiles" / f"{name}.yaml"
        tried.append(p)
        if p.is_file():
            return p
    have = sorted({x.stem for r in config_roots() for x in (r / "config" / "profiles").glob("*.yaml")})
    raise FileNotFoundError(f"no profile {name!r}; looked in {[str(t) for t in tried]}; available: {have} "
                            f"(set GERMQUANT_CONFIG_ROOT to the checkout that holds config/profiles)")


def default_config_path() -> Path:
    """``config/config.yaml`` under the first config root that has it: the config the drag-and-drop
    scripts and the README use when nothing else is named."""
    tried = []
    for root in config_roots():
        p = root / "config" / "config.yaml"
        tried.append(p)
        if p.is_file():
            return p
    raise FileNotFoundError(f"no config/config.yaml found; looked in {[str(t) for t in tried]} "
                            f"(pass a config path or set GERMQUANT_CONFIG_ROOT)")
