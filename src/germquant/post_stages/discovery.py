"""Finding post-stages without importing them (docs/POST_STAGES.md; decision D31).

An installed post-stage is an entry point in the group ``germquant.post_stages`` (``name = "package.module:fn"``)
plus a declaration file, germquant_post_stage.yaml, inside the top-level package that the entry point names.
Discovery reads the entry points from the installed distributions' metadata (importlib.metadata), finds that
package with importlib.util.find_spec on its top-level name (which locates the package without executing it, and
works for editable installs, whose metadata lists no package files) and reads the file with yaml.safe_load. No
post-stage module is imported until a run calls it (``load_post_stage``), so listing post-stages and checking
a pack's needs.post_stages stay free of torch, CuPy, SpotMAX and Cellpose (gate G6).

A post-stage can also be registered in process (``register_post_stage``), as the backend tests register
backends: for tests and notebooks. A registration shadows an installed post-stage of the same name.

An installed post-stage that cannot be used (no declaration file, a malformed one, a name declared by two
distributions, an entry point that names a module rather than a package) is never offered; it is listed, with
the reason, by ``list_post_stages(include_unusable=True)``.
"""
from __future__ import annotations

import importlib.util
import logging
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

import yaml

from ..fsutil import long_path
from .declaration import (
    DECLARATION_FILE,
    declaration_from_document,
    declaration_problems,
    name_problem,
    normalise_declaration,
)

log = logging.getLogger(__name__)

ENTRY_POINT_GROUP = "germquant.post_stages"

# name -> {"fn": callable, "declaration": the listing form}; in-process registrations (tests, notebooks)
_REGISTRY: dict[str, dict] = {}


class PostStageError(LookupError):
    """A post-stage that is not installed here, cannot be used, or fails to load."""


# ---- in-process registration -------------------------------------------------------------------------------

def register_post_stage(name: str, fn: Callable, declaration: dict, *, replace: bool = False) -> None:
    """Make ``fn`` available as post-stage ``name`` in this process, with ``declaration`` (the mapping a
    germquant_post_stage.yaml holds under post_stages.<name>). The declaration is checked by the rules of an
    installed one (ValueError naming every problem). Registering a name twice is an error unless
    ``replace=True``. A registration shadows an installed post-stage of the same name."""
    if not callable(fn):
        raise TypeError(f"post-stage {name!r} must be a callable fn(inputs, params), got {type(fn).__name__}")
    probs = declaration_problems(name, declaration)
    if probs:
        raise ValueError("; ".join(probs))
    if name in _REGISTRY and not replace:
        raise ValueError(f"post-stage {name!r} is already registered in this process; pass replace=True")
    _REGISTRY[name] = {"fn": fn, "declaration": normalise_declaration(name, declaration)}


def unregister_post_stage(name: str) -> None:
    """Undo register_post_stage (PostStageError for a name that is not registered in this process)."""
    if name not in _REGISTRY:
        raise PostStageError(f"post-stage {name!r} is not registered in this process")
    del _REGISTRY[name]


# ---- installed post-stages ---------------------------------------------------------------------------------

def _iter_entry_points() -> Iterable[Any]:
    """The raw entry points of the group (a seam the tests may replace)."""
    from importlib.metadata import entry_points

    return entry_points(group=ENTRY_POINT_GROUP)


def _dist_name(ep) -> str | None:
    dist = getattr(ep, "dist", None)
    name = getattr(dist, "name", None) if dist is not None else None
    return str(name) if name else None


def read_declaration(top: str) -> tuple[dict | None, Path | None, list[str]]:
    """(the parsed germquant_post_stage.yaml of top-level package ``top``, its path, problems), found through
    importlib.util.find_spec, which locates the package without executing it."""
    try:
        spec = importlib.util.find_spec(top)
    except Exception as e:  # noqa: BLE001 - a broken package must not stop a listing
        return None, None, [f"package {top!r} cannot be located ({type(e).__name__}: {e})"]
    if spec is None:
        return None, None, [f"package {top!r} is not importable in this environment"]
    locations = list(spec.submodule_search_locations or [])
    if not locations:
        return None, None, [f"{top!r} is a module, not a package: {DECLARATION_FILE} must sit in the top-level "
                            f"package the entry point names"]
    for loc in locations:
        path = Path(loc) / DECLARATION_FILE
        try:
            raw = Path(long_path(path)).read_bytes()
        except FileNotFoundError:
            continue
        except OSError as e:
            return None, path, [f"{path} cannot be read ({type(e).__name__}: {e})"]
        try:
            text = raw.decode("utf-8-sig")
        except UnicodeDecodeError as e:
            return None, path, [f"{path} is not UTF-8 text (byte {e.start})"]
        try:
            return yaml.safe_load(text), path, []
        except yaml.YAMLError as e:
            mark = getattr(e, "problem_mark", None)
            where = f" (line {mark.line + 1})" if mark is not None else ""
            return None, path, [f"{path} cannot be read as YAML{where}"]
    return None, None, [f"package {top!r} holds no {DECLARATION_FILE} (searched {locations})"]


def _installed() -> tuple[dict[str, dict], list[dict]]:
    """({name: row} of the usable installed post-stages, [rows of the unusable ones]); rows carry the entry
    point under "_ep" (left out of every listing)."""
    try:
        eps = list(_iter_entry_points())
    except Exception as e:  # noqa: BLE001 - a damaged installation must not stop a listing
        log.warning("could not read the %r entry points (%s: %s)", ENTRY_POINT_GROUP, type(e).__name__, e)
        return {}, []
    by_name: dict[str, list] = {}
    for ep in eps:
        by_name.setdefault(str(getattr(ep, "name", "")).strip(), []).append(ep)
    usable: dict[str, dict] = {}
    unusable: list[dict] = []
    for name, group in sorted(by_name.items()):
        ep = group[0]
        value = str(getattr(ep, "value", ""))
        base = {"name": name, "version": None, "dist": _dist_name(ep), "title": "", "reads": {}, "writes": {},
                "params": [], "source": f"entry point {value}", "_ep": ep}
        if len(group) > 1:
            dists = sorted({_dist_name(e) or "?" for e in group})
            unusable.append({**base, "problems": [f"post-stage {name!r} is declared by more than one distribution "
                                                  f"({', '.join(dists)})"]})
            continue
        bad_name = name_problem(name)
        if bad_name:
            unusable.append({**base, "problems": [f"entry point {name!r}: {bad_name}"]})
            continue
        top = value.split(":", 1)[0].strip().split(".", 1)[0]
        if not top:
            unusable.append({**base, "problems": [f"entry point {name!r} names no module ({value!r})"]})
            continue
        doc, _path, probs = read_declaration(top)
        decl = None
        if not probs:
            decl, probs = declaration_from_document(doc, name)
        if probs:
            unusable.append({**base, "problems": probs})
            continue
        usable[name] = {**normalise_declaration(name, decl), "dist": base["dist"], "source": base["source"],
                        "_ep": ep}
    return usable, unusable


def _public(row: dict, *, usable: bool | None = None) -> dict:
    out = {k: v for k, v in row.items() if k != "_ep" and k != "problems"}
    if usable is not None:
        out["usable"] = usable
        out["problems"] = list(row.get("problems") or [])
    return out


def list_post_stages(*, include_unusable: bool = False) -> list[dict]:
    """The post-stages this environment offers, sorted by name, read from metadata and the declaration files
    without importing any: [{"name", "version", "dist", "title", "reads", "writes", "params", "source"}].
    ``dist`` is the distribution's name (None for an in-process registration), ``source`` is
    "entry point <module:attr>" or "registered". With ``include_unusable`` every row also carries "usable" and
    "problems", and installed post-stages that cannot be used are listed too."""
    installed, broken = _installed()
    rows = {name: row for name, row in installed.items()}
    for name, reg in _REGISTRY.items():
        rows[name] = {**reg["declaration"], "dist": None, "source": "registered"}
    out = [_public(rows[n], usable=True if include_unusable else None) for n in sorted(rows)]
    if include_unusable:
        out += [_public(r, usable=False) for r in broken if r["name"] not in rows]
        out.sort(key=lambda r: r["name"])
    return out


def declarations() -> dict[str, dict]:
    """{name: declaration in the listing form} of every usable post-stage here (installed and registered)."""
    return {r["name"]: r for r in list_post_stages()}


def unusable_reasons() -> dict[str, list[str]]:
    """{name: problems} of the installed post-stages that cannot be used (and are not shadowed)."""
    installed, broken = _installed()
    return {r["name"]: list(r["problems"]) for r in broken if r["name"] not in _REGISTRY}


def load_post_stage(name: str) -> tuple[Callable, dict]:
    """(the callable, the declaration in listing form) of post-stage ``name``: a registration, else the
    installed entry point, whose module is imported here and only here. PostStageError when it is not
    installed, cannot be used, or fails to load."""
    if name in _REGISTRY:
        reg = _REGISTRY[name]
        return reg["fn"], dict(reg["declaration"])
    installed, broken = _installed()
    if name not in installed:
        for row in broken:
            if row["name"] == name:
                raise PostStageError(f"post-stage {name!r} is installed but not usable: {'; '.join(row['problems'])}")
        known = sorted(installed) + sorted(_REGISTRY)
        raise PostStageError(f"post-stage {name!r} is not installed in this environment (installed: "
                             f"{', '.join(known) or 'none'}; a post-stage is a distribution with an entry point in "
                             f"the {ENTRY_POINT_GROUP!r} group, docs/POST_STAGES.md)")
    row = installed[name]
    try:
        fn = row["_ep"].load()
    except (Exception, SystemExit) as e:  # noqa: BLE001 - a broken post-stage fails its own stage only
        raise PostStageError(f"post-stage {name!r} failed to load ({row['source']}: {type(e).__name__}: {e})") from e
    if not callable(fn):
        raise PostStageError(f"post-stage {name!r} ({row['source']}) is not callable")
    return fn, _public(row)
