"""The project manifest, germquant_project.yaml, format 1 (docs/PROJECTS.md, "Layout and manifest"; design 7.2).

The manifest is not hashed. It holds, in this order:

* ``format`` (1), ``name``, ``id`` (a uuid4);
* ``created`` and ``saved`` stamps, each ``{at, germquant, git_sha, plugin}``: when, the engine version and git sha
  that wrote it, and the plugin version when the plugin asked for the write (``plugin`` is null otherwise);
* ``pack``: ``{id, version, recipe, sha256}``, the pinned pack version (copied into packs/<id>/<version>/ at
  creation) and the pack_sha256 of that copy;
* ``recipe``: ``{file, expected_config_hash, overrides, roles, origins}``: recipe.yaml or null, the config hash
  the project had when it was last written (the pack recipe's own hash while nothing is edited), the stage
  switches turned off (``{"spots.enabled": false}``, what ``--no-spots`` sets; they enter the hash), the channel
  role choices behind channel_map.yaml (``{role: channel name}``), and where each recipe.yaml key came from
  (``{dotted key: form | tuning | model}``);
* ``run``: ``{xy_stride, z_range, device, runner, interpreters}``, the run settings (device and runner null: the
  engine's choice; interpreters ``{recipe: {python, germquant, git_sha}}``);
* ``models``: the models the project chose, ``[{key, ref, sha256, card, card_sha256, source}]``;
* ``images``: ``[{id, path, size, mtime, timepoints (more than one only), excluded, reason}]``, a path inside the
  project relative with forward slashes, any other absolute (UNC paths verbatim);
* ``samples``: the sample sheet's file name in the project, or null;
* ``validated``: ``{value, why, record}``, computed at every write.

It is written atomically (a temporary file next to it, fsync, then os.replace, all through fsutil.long_path,
with five retries over 2 s on PermissionError, which an indexer or a virus scanner holding the file raises on
Windows), as ASCII with escapes (yaml.safe_dump, allow_unicode False) and LF line ends. Only the interactive
writers write it, under the lock (project/lock.py); runs never do.
"""
from __future__ import annotations

import os
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

import yaml

from ..fsutil import long_path

FORMAT = 1
MANIFEST_NAME = "germquant_project.yaml"
LOCK_NAME = "germquant_project.lock"
RECIPE_NAME = "recipe.yaml"
CHANNEL_MAP_NAME = "channel_map.yaml"
# the project's own folders (design 7.1): never searched for images, and never a home for one
PROJECT_DIRS = ("results", "runs", "training", "tuning", "annotations", "packs")
KEY_ORDER = ("format", "name", "id", "created", "saved", "pack", "recipe", "run", "models", "images", "samples",
             "validated")
RECIPE_KEYS = ("file", "expected_config_hash", "overrides", "roles", "origins")
RUN_KEYS = ("xy_stride", "z_range", "device", "runner", "interpreters")
RETRIES = 5                      # after the first attempt: five more over RETRIES * RETRY_WAIT_S = 2 s
RETRY_WAIT_S = 0.4

# seams the tests replace (a failing or a briefly locked file system)
_replace = os.replace
_sleep = time.sleep


class ProjectError(ValueError):
    """A project request that is refused; the message says why and what to do."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def stamp(plugin_version: str | None = None) -> dict:
    """A created or saved stamp: ``{at, germquant, git_sha, plugin}``."""
    from .. import __version__
    from ..provenance import git_sha

    return {"at": utc_now(), "germquant": __version__, "git_sha": git_sha(),
            "plugin": None if plugin_version is None else str(plugin_version)}


def dump_yaml(data, *, sort_keys: bool = False) -> bytes:
    """YAML as ASCII bytes with LF line ends: every non-ASCII character becomes an escape (allow_unicode False)."""
    text = yaml.safe_dump(data, sort_keys=sort_keys, allow_unicode=False, default_flow_style=False, width=4096)
    return text.encode("ascii")


def ordered(manifest: dict) -> dict:
    """The manifest with its keys in KEY_ORDER (and the recipe and run blocks in theirs), unknown keys after."""
    out = {k: manifest[k] for k in KEY_ORDER if k in manifest}
    out.update({k: v for k, v in manifest.items() if k not in out})
    for block, keys in (("recipe", RECIPE_KEYS), ("run", RUN_KEYS)):
        if isinstance(out.get(block), dict):
            b = out[block]
            out[block] = {**{k: b[k] for k in keys if k in b}, **{k: v for k, v in b.items() if k not in keys}}
    return out


def manifest_path(root) -> Path:
    return Path(root) / MANIFEST_NAME


def new_manifest(*, name: str, pack: dict, expected_config_hash: str, validated: dict,
                 plugin_version: str | None = None) -> dict:
    st = stamp(plugin_version)
    return {"format": FORMAT, "name": str(name), "id": str(uuid.uuid4()), "created": st, "saved": dict(st),
            "pack": dict(pack),
            "recipe": {"file": None, "expected_config_hash": expected_config_hash, "overrides": {}, "roles": {},
                       "origins": {}},
            "run": {"xy_stride": 1, "z_range": None, "device": None, "runner": None, "interpreters": {}},
            "models": [], "images": [], "samples": None, "validated": dict(validated)}


def _normalise(data: dict) -> dict:
    """Fill the optional blocks a hand-edited manifest may lack, so every reader sees the same shape."""
    recipe = data.setdefault("recipe", {})
    if not isinstance(recipe, dict):
        raise ProjectError("the manifest's recipe entry is not a mapping")
    recipe.setdefault("file", None)
    recipe.setdefault("expected_config_hash", None)
    for key in ("overrides", "roles", "origins"):
        if recipe.get(key) is None:
            recipe[key] = {}
    run = data.setdefault("run", {})
    if not isinstance(run, dict):
        raise ProjectError("the manifest's run entry is not a mapping")
    run.setdefault("xy_stride", 1)
    for key in ("z_range", "device", "runner"):
        run.setdefault(key, None)
    if run.get("interpreters") is None:
        run["interpreters"] = {}
    for key in ("models", "images"):
        if data.get(key) is None:
            data[key] = []
    data.setdefault("samples", None)
    pack = data.get("pack")
    if not isinstance(pack, dict) or not pack.get("id") or not pack.get("version") or not pack.get("recipe"):
        raise ProjectError("the manifest names no pack {id, version, recipe}")
    return data


def read_manifest(root) -> dict:
    """The parsed manifest of the project in ``root``. FileNotFoundError when there is none; ProjectError when it
    cannot be read as format 1."""
    p = manifest_path(root)
    try:
        with open(long_path(p), "rb") as fh:
            raw = fh.read()
    except FileNotFoundError:
        raise FileNotFoundError(f"{root}: no {MANIFEST_NAME}, so this folder is not a germquant project (make one "
                                f"with `germquant project new DIR --pack ID`)") from None
    try:
        data = yaml.safe_load(raw.decode("utf-8"))
    except (UnicodeDecodeError, yaml.YAMLError) as e:
        raise ProjectError(f"{p}: not readable as YAML ({str(e).splitlines()[0]})") from e
    if not isinstance(data, dict):
        raise ProjectError(f"{p}: not a mapping")
    if data.get("format") != FORMAT:
        raise ProjectError(f"{p}: format {data.get('format')!r} is not format {FORMAT} (this engine reads format "
                           f"{FORMAT} projects)")
    return _normalise(data)


def atomic_write(path, data: bytes) -> Path:
    """Write ``data`` to ``path`` through a temporary file in the same folder, fsync, then os.replace; a failure at
    any point leaves the old file as it was and no temporary file behind."""
    path = Path(path)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex[:8]}.tmp")
    try:
        with open(long_path(tmp), "wb") as fh:
            fh.write(data)
            fh.flush()
            os.fsync(fh.fileno())
        for attempt in range(RETRIES + 1):
            try:
                _replace(long_path(tmp), long_path(path))
                return path
            except PermissionError:
                if attempt == RETRIES:
                    raise
                _sleep(RETRY_WAIT_S)
    except BaseException:
        try:
            os.remove(long_path(tmp))
        except OSError:
            pass
        raise
    return path                                   # not reached: the loop returns or raises


def remove_file(path) -> bool:
    """Remove ``path`` (True) or report that it was not there (False), with the retries of atomic_write."""
    for attempt in range(RETRIES + 1):
        try:
            os.remove(long_path(path))
            return True
        except FileNotFoundError:
            return False
        except PermissionError:
            if attempt == RETRIES:
                raise
            _sleep(RETRY_WAIT_S)
    return False


def write_manifest(root, manifest: dict) -> Path:
    return atomic_write(manifest_path(root), dump_yaml(ordered(manifest)))
