"""Recipe edits, channel roles and project models (docs/PROJECTS.md; design 7.3 and 7.5; PLAN.md decision D19).

An unmodified project has no recipe.yaml: project_config loads the pack recipe file itself, so the project's
config_hash equals `germquant run --pack` and the pack's expected_config_hash. Any edit writes recipe.yaml,

    extends: pack:<id>@<version>/<recipe>
    <the edits, nested, keys sorted>

with yaml.safe_dump (ASCII with escapes), and nothing else: no comment, no project name and no path, so the same
edit of the same pack version hashes alike in every project and on every machine. A role edit writes
channel_map.yaml, the pack recipe's channel map with each edited role's names replaced by the exact channel name
chosen (its index fallback, required flag and marker kept), and points recipe.yaml at it.

write_recipe validates by origin: ``form`` edits only the pack recipe's declared tunables, with a scalar of the
declared type inside min, max and choices; ``tuning`` edits spot keys only (spots.*; the tuner checks the keys
against its grid); ``model`` edits segmentation.nuclei.cellpose_model only (use_model); ``roles`` edits the
channel map (``roles=``). An edit that sets a key to the pack recipe's own value is no edit, and once no edit is
left recipe.yaml (and channel_map.yaml) are removed, so the project hashes as the pack recipe again.

Every change is computed first, then written under the project lock: channel_map.yaml, then recipe.yaml, then the
manifest, each atomically. plan_change loads new recipe files through load_config in a temporary folder (with a
channel map that recipe.yaml names by hand copied next to them), so the hash is the one the project will have; a
change that keeps the recipe files loads them where they are, as a run does. A recipe file the manifest does not
record (put there by hand) is never removed, and a change that would replace it is refused (project/__init__.py,
_apply). A dry run computes and writes nothing.
"""
from __future__ import annotations

import copy
import hashlib
import json
import logging
import os
import re
import sys
import tempfile
from pathlib import Path

import yaml

from .manifest import CHANNEL_MAP_NAME, RECIPE_NAME, ProjectError, dump_yaml

log = logging.getLogger(__name__)

MODEL_KEY = "segmentation.nuclei.cellpose_model"
CHANNEL_MAP_KEY = "io.channel_map"
ORIGINS = ("form", "tuning", "model", "roles")
DEVICES = ("auto", "cuda", "mps", "cpu")
RUNNERS = ("subprocess", "inprocess")
_KEY_RE = re.compile(r"^[A-Za-z0-9_]+(\.[A-Za-z0-9_]+)*$")
_MISSING = object()


# ---- the pack, its recipe and its channel map -------------------------------------------------------------------

def extends_ref(manifest: dict) -> str:
    p = manifest["pack"]
    return f"pack:{p['id']}@{p['version']}/{p['recipe']}"


def resolve_pack(root, manifest: dict, *, refresh: bool = False):
    """``(Pack, Recipe)`` of the pinned version, found with the project's packs/ searched first."""
    from ..packs import discover, find_pack

    p = manifest["pack"]
    if refresh:
        discover(root, refresh_cache=True)
    pack = find_pack(f"{p['id']}@{p['version']}", project=root)
    return pack, pack.recipe(p["recipe"])


def pack_recipe_config(root, pack, recipe):
    """The Config of the pack recipe file itself (no project edit, no override)."""
    from ..config import load_config
    from ..packs import recipe_profile_path, use_project

    with use_project(root):
        return load_config(recipe_profile_path(pack, recipe))


def pack_channel_map(cfg) -> dict:
    """The parsed channel map a recipe Config names (match_by and roles, in the file's order)."""
    from ..fsutil import read_config_text

    cm = str(cfg.get("io.channel_map"))
    path = Path(cm) if Path(cm).is_absolute() else Path(cfg.base_dir) / cm
    data = yaml.safe_load(read_config_text(path)) or {}
    if not isinstance(data.get("roles"), dict):
        raise ProjectError(f"{path}: the recipe's channel map has no roles")
    return data


# ---- reading and rendering ------------------------------------------------------------------------------------

def flatten(d: dict, prefix: str = "") -> dict:
    out = {}
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else str(k)
        if isinstance(v, dict) and v:
            out.update(flatten(v, key))
        else:
            out[key] = v
    return out


def _set_nested(d: dict, dotted: str, value) -> None:
    parts = dotted.split(".")
    node = d
    for i, part in enumerate(parts[:-1]):
        nxt = node.setdefault(part, {})
        if not isinstance(nxt, dict):
            raise ProjectError(f"the edits {'.'.join(parts[:i + 1])} and {dotted} overlap")
        node = nxt
    if isinstance(node.get(parts[-1]), dict):
        raise ProjectError(f"the edit {dotted} overlaps a deeper edit under it")
    node[parts[-1]] = value


def yaml_reason(e: yaml.YAMLError) -> str:
    """A YAML error in one line: what is wrong and where (line and column, counted from 1)."""
    problem, mark = getattr(e, "problem", None), getattr(e, "problem_mark", None)
    if problem and mark is not None:
        return f"{problem}, line {mark.line + 1} column {mark.column + 1}"
    lines = str(e).strip().splitlines()
    return lines[0] if lines else type(e).__name__


def read_edits(root, manifest: dict) -> dict:
    """``{dotted key: value}`` that recipe.yaml sets: the extends line left out, and io.channel_map left out when it
    is the role editor's channel_map.yaml. ProjectError when the file is missing or is not YAML (a hand edit)."""
    from ..fsutil import read_config_text

    if not manifest["recipe"].get("file"):
        return {}
    p = Path(root) / str(manifest["recipe"]["file"])
    if not p.is_file():
        raise ProjectError(f"{p}: the manifest names it but it is missing; restore it, or go back to the pack recipe "
                           f"with `germquant project edit DIR --revert`")
    try:
        data = yaml.safe_load(read_config_text(p)) or {}
    except yaml.YAMLError as e:
        raise ProjectError(f"{p}: not readable as YAML ({yaml_reason(e)}); fix the file, or go back to the pack "
                           f"recipe with `germquant project edit DIR --revert`") from None
    if not isinstance(data, dict):
        raise ProjectError(f"{p}: not a mapping")
    data.pop("extends", None)
    flat = flatten(data)
    if manifest["recipe"].get("roles") and flat.get(CHANNEL_MAP_KEY) == CHANNEL_MAP_NAME:
        flat.pop(CHANNEL_MAP_KEY)
    return flat


def render_recipe(manifest: dict, edits: dict, has_roles: bool) -> str | None:
    """recipe.yaml's text for these edits, or None when there is nothing to write."""
    if not edits and not has_roles:
        return None
    nested: dict = {}
    for key in sorted(edits):
        _set_nested(nested, key, edits[key])
    if has_roles:
        _set_nested(nested, CHANNEL_MAP_KEY, CHANNEL_MAP_NAME)
    return (dump_yaml({"extends": extends_ref(manifest)}) + dump_yaml(nested, sort_keys=True)).decode("ascii")


def render_channel_map(pack_map: dict, roles: dict) -> str:
    """The pack's channel map with each role of ``roles`` bound by its exact channel name."""
    data = copy.deepcopy(pack_map)
    for role, channel in roles.items():
        data["roles"][role]["names"] = [str(channel)]
    return dump_yaml(data).decode("ascii")


# ---- validation ---------------------------------------------------------------------------------------------

def check_key(key: str) -> str:
    if not isinstance(key, str) or not _KEY_RE.fullmatch(key):
        raise ProjectError(f"{key!r} is not a dotted config key such as spots.effect_size_min")
    return key


def _scalar(key: str, value):
    if isinstance(value, (dict, list, tuple, set)) or value is None:
        raise ProjectError(f"{key}: an edit sets one scalar value (a number, true or false, or a word), got {value!r}")
    if not isinstance(value, (bool, int, float, str)):
        raise ProjectError(f"{key}: {type(value).__name__} values cannot be written to a recipe")
    return value


def _number(key: str, value, kind: str, t: dict):
    is_num = isinstance(value, (int, float)) and not isinstance(value, bool)
    if kind == "int" and not (is_num and float(value).is_integer()):
        raise ProjectError(f"{key} takes a whole number (int), got {value!r}")
    if kind == "float" and not is_num:
        raise ProjectError(f"{key} takes a float (a number), got {value!r}")
    v = int(value) if kind == "int" else float(value)
    lo, hi = t.get("min"), t.get("max")
    for bound in (lo, hi):
        if bound is not None and (isinstance(bound, bool) or not isinstance(bound, (int, float))):
            raise ProjectError(f"{key}: the pack declares the bound {bound!r}, which is not a number; the form cannot "
                               f"check the value (a pack problem: `germquant pack lint` names it)")
    if lo is not None and v < lo:
        raise ProjectError(f"{key} must be at least {lo:g}, got {v:g}")
    if hi is not None and v > hi:
        raise ProjectError(f"{key} must be at most {hi:g}, got {v:g}")
    return v


def _same(a, b) -> bool:
    if b is _MISSING:
        return False
    if isinstance(a, bool) or isinstance(b, bool):
        return isinstance(a, bool) and isinstance(b, bool) and a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return float(a) == float(b)
    return type(a) is type(b) and a == b


def check_edit(origin: str, key: str, value, pack, recipe):
    """The value to write for ``key`` from an edit of ``origin`` (D19), or ProjectError naming the key and the rule."""
    check_key(key)
    if origin == "form":
        tun = {str(t.get("key")): t for t in recipe.tunables or [] if isinstance(t, dict)}
        if key not in tun:
            raise ProjectError(f"{key} is not a tunable of the pack recipe {pack.ref}/{recipe.id}: the parameter form "
                               f"edits only its tunables ({', '.join(sorted(tun)) or 'none declared'})")
        value = _scalar(key, value)
        t = tun[key]
        kind = str(t.get("type", "")).lower()
        if kind in ("float", "int"):
            return _number(key, value, kind, t)
        if kind == "bool":
            if not isinstance(value, bool):
                raise ProjectError(f"{key} takes true or false (bool), got {value!r}")
            return value
        if kind == "choice":
            choices = list(t.get("choices") or [])
            if not any(_same(value, c) for c in choices):
                raise ProjectError(f"{key} must be one of the choices {choices}, got {value!r}")
            return next(c for c in choices if _same(value, c))
        if kind in ("str", "string"):
            if not isinstance(value, str):
                raise ProjectError(f"{key} takes a word (str), got {value!r}")
            return value
        raise ProjectError(f"{key}: the form does not write tunables of type {t.get('type')!r}")
    if origin == "tuning":
        if not key.startswith("spots."):
            raise ProjectError(f"{key}: a tuning edit sets spot keys only (spots.*, the keys of the tuning grid)")
        return _scalar(key, value)
    if origin == "model":
        if key != MODEL_KEY or not isinstance(value, str):
            raise ProjectError(f"{key}: a model edit sets {MODEL_KEY} (a model path) only; use use_model")
        return value
    if origin == "roles":
        raise ProjectError("a roles edit changes the channel map: pass roles={role: channel name}, not edits")
    raise ProjectError(f"unknown origin {origin!r}: one of {', '.join(ORIGINS)}")


def validated_record(pack, recipe, *, edited: bool, overrides: dict) -> dict:
    """The manifest's validated entry: ``{value, why, record}``."""
    from ..packs.format import validation_record_passes

    if edited:
        return {"value": False, "why": "edited: recipe.yaml changes the pack recipe, so the pack's validation does "
                                       "not cover this project", "record": None}
    if overrides:
        return {"value": False, "why": f"edited: the stage switches {', '.join(sorted(overrides))} change the config "
                                       f"hash, so the pack's validation does not cover this project", "record": None}
    claims = recipe.status == "validated" or recipe.raw.get("validated") is True
    if not claims:
        return {"value": False, "why": f"not validated: the pack recipe {pack.ref}/{recipe.id} is provisional",
                "record": None}
    ok, why = validation_record_passes(pack.folder, recipe.validation)
    if not ok:
        return {"value": False, "why": f"not validated: {why}", "record": None}
    return {"value": True, "why": f"validated: {why}", "record": f"packs/{pack.id}/{pack.version}/{recipe.validation}"}


def _check_save(root, manifest: dict, pack, save: dict) -> None:
    """Apply the run, samples and override fields of ``save`` to ``manifest`` (save_project, `project edit`)."""
    run = manifest["run"]
    if "device" in save:
        d = str(save["device"] or "").strip().lower()
        if d and d not in DEVICES:
            raise ProjectError(f"device {save['device']!r} is not one of {', '.join(DEVICES)}")
        run["device"] = d or None
    if "runner" in save:
        r = str(save["runner"] or "").strip().lower()
        if r and r not in RUNNERS:
            raise ProjectError(f"runner {save['runner']!r} is not one of {', '.join(RUNNERS)}")
        run["runner"] = r or None
    if "xy_stride" in save:
        s = save["xy_stride"]
        if isinstance(s, bool) or not isinstance(s, int) or s < 1:
            raise ProjectError(f"xy stride {s!r} must be a whole number of 1 or more (1: full resolution)")
        run["xy_stride"] = int(s)
    if "z_range" in save:
        z = save["z_range"]
        if z is None or (isinstance(z, (list, tuple)) and len(z) == 0):
            run["z_range"] = None
        else:
            try:
                z0, z1 = (int(v) for v in z)
                ok = all(not isinstance(v, bool) for v in z) and 0 <= z0 < z1
            except (TypeError, ValueError):
                ok = False
            if not ok:
                raise ProjectError(f"z range {z!r} must be two planes Z0 Z1 with 0 <= Z0 < Z1 (or none)")
            run["z_range"] = [z0, z1]
    if "interpreters" in save:
        cur = dict(run.get("interpreters") or {})
        ids = [r.id for r in pack.recipes]
        for rid, py in dict(save["interpreters"]).items():
            if rid not in ids:
                raise ProjectError(f"--interpreter {rid}=...: {rid!r} is not a recipe of the pack {pack.ref} "
                                   f"(its recipes: {', '.join(ids)})")
            if py is None or (isinstance(py, str) and py.strip() in ("", "none")):
                cur.pop(rid, None)
            elif isinstance(py, dict):
                if not py.get("python"):
                    raise ProjectError(f"the interpreter record of {rid} names no python")
                cur[rid] = {"python": str(py["python"]), "germquant": str(py.get("germquant", "unknown")),
                            "git_sha": str(py.get("git_sha", "unknown"))}
            else:
                cur[rid] = interpreter_record(str(py))
        run["interpreters"] = cur
    if "samples" in save:
        s = str(save["samples"] or "").strip()
        if not s or s.lower() == "none":
            manifest["samples"] = None
        else:
            p = Path(s)
            full = Path(os.path.abspath(p if p.is_absolute() else Path(root) / p))
            try:
                rel = full.relative_to(Path(os.path.abspath(root))).as_posix()
            except ValueError:
                raise ProjectError(f"samples {s}: the sample sheet must be a file inside the project "
                                   f"({root}), such as samples.csv") from None
            if not full.is_file():
                raise ProjectError(f"samples {s}: no such file in the project ({full})")
            manifest["samples"] = rel
    if "overrides" in save:
        from ..stages import cli_switches

        keys = {st.config_key: st.name for st in cli_switches()}
        out = {}
        for key, value in dict(save["overrides"]).items():
            if key not in keys:
                raise ProjectError(f"override {key!r}: a project overrides only the stage switches of --no-<stage> "
                                   f"({', '.join(sorted(keys))})")
            if value is not False:
                raise ProjectError(f"override {key}: a stage switch is overridden to false only (--no-{keys[key]})")
            out[key] = False
        manifest["recipe"]["overrides"] = out


def interpreter_record(python: str) -> dict:
    """``{python, germquant, git_sha}`` of an engine interpreter: this one directly, another by running its
    torch-free ``-m germquant.packs --version --json``."""
    import subprocess

    p = Path(python)
    if not p.is_file():
        raise ProjectError(f"interpreter {python}: no such file")
    try:
        same = os.path.samefile(p, sys.executable)
    except OSError:
        same = False
    if same:
        from ..packs import version_info

        v = version_info()
        return {"python": str(python), "germquant": str(v["germquant"]), "git_sha": str(v["git_sha"])}
    try:
        r = subprocess.run([str(p), "-m", "germquant.packs", "--version", "--json"], capture_output=True, text=True,
                           timeout=120)
        info = json.loads(r.stdout.strip().splitlines()[-1])
        return {"python": str(python), "germquant": str(info["germquant"]), "git_sha": str(info["git_sha"])}
    except Exception as e:  # noqa: BLE001 - any failure means this is not an engine interpreter
        raise ProjectError(f"interpreter {python}: `-m germquant.packs --version --json` did not answer ({type(e).__name__}"
                           f": {e}); is germquant installed in that environment?") from None


# ---- planning a change ---------------------------------------------------------------------------------------

def staged_config(root, pack, recipe, recipe_bytes: bytes | None, cm_bytes: bytes | None, overrides: dict, *,
                  hand_map=None):
    """The Config the project will have with these files: loaded through load_config from a temporary folder
    (the hash reads the text only, never the folder), with the project's packs/ searched first. ``hand_map``: the
    channel map that recipe.yaml names by a hand edit of io.channel_map, copied from the project into that folder
    when ``cm_bytes`` (the role editor's map) is None."""
    from ..config import load_config
    from ..packs import recipe_profile_path, use_project

    with use_project(root):
        if recipe_bytes is None:
            cfg = load_config(recipe_profile_path(pack, recipe))
        else:
            with tempfile.TemporaryDirectory(prefix="germquant_recipe_") as tmp:
                t = Path(tmp)
                (t / RECIPE_NAME).write_bytes(recipe_bytes)
                if cm_bytes is not None:
                    (t / CHANNEL_MAP_NAME).write_bytes(cm_bytes)
                elif hand_map is not None:
                    _stage_hand_map(root, hand_map, t)
                cfg = load_config(t / RECIPE_NAME)
    for key, value in sorted(overrides.items()):
        cfg.set(key, value)
    return cfg


def _stage_hand_map(root, hand_map, folder: Path) -> None:
    """Copy the channel map recipe.yaml names by hand (io.channel_map: a lab's own map, several tokens per role or
    match_by index, which the role editor does not write) from the project into the staging folder, at the same
    relative path, so the staged recipe loads as the project's own does. An absolute path loads where it is; a
    relative path that leaves the project folder is refused."""
    if not isinstance(hand_map, str) or not hand_map.strip() or Path(hand_map).is_absolute() or Path(hand_map).anchor:
        return
    parts = Path(os.path.normpath(hand_map)).parts
    if not parts or parts[0] == os.pardir:
        raise ProjectError(f"{RECIPE_NAME} names its channel map {hand_map} (io.channel_map) by a path that leaves the "
                           f"project folder: name it by its absolute path, or copy it into the project")
    data = _file_bytes(root, hand_map)
    if data is None:
        raise ProjectError(f"{RECIPE_NAME} names its channel map {hand_map} (io.channel_map), which is not in the "
                           f"project folder ({Path(root) / hand_map})")
    dst = folder / hand_map
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_bytes(data)


def files_config(root, pack, recipe, overrides: dict, *, recipe_file: bool):
    """The Config of the project's recipe files as they are in its folder: recipe.yaml (with the channel map it
    names, the role editor's or a hand-written one) when ``recipe_file``, else the pack recipe file; then the
    stage switches. This is what project_config loads for a run."""
    from ..config import load_config
    from ..packs import recipe_profile_path, use_project

    with use_project(root):
        cfg = load_config(Path(root) / RECIPE_NAME) if recipe_file else load_config(recipe_profile_path(pack, recipe))
    for key, value in sorted(overrides.items()):
        cfg.set(key, value)
    return cfg


def _file_bytes(root, name: str) -> bytes | None:
    from ..fsutil import long_path

    try:
        with open(long_path(Path(root) / name), "rb") as fh:
            return fh.read()
    except FileNotFoundError:
        return None


def bindings(root, manifest: dict, channel_map, *, image_id=None) -> tuple[dict, list[str], list[str]]:
    """``(bindings, problems, warnings)`` of a channel map on the project's images (excluded images left out):
    for each image ``{role: {"channel", "how": "name" | "index" | "unbound"}}``; a required role that binds no
    channel and two roles on one channel are problems; an unreadable header is a warning (that image is left out
    of the bindings)."""
    from .headers import read_header
    from .images import resolve_path

    images = [e for e in manifest.get("images") or [] if not e.get("excluded")]
    if image_id is not None:
        images = [e for e in images if e.get("id") == image_id]
        if not images:
            raise ProjectError(f"no image {image_id!r} in the project")
    out: dict = {}
    problems: list[str] = []
    warnings: list[str] = []
    for e in images:
        try:
            names = list(read_header(resolve_path(root, e["path"]))["channel_names"])
        except Exception as ex:  # noqa: BLE001 - one unreadable image must not hide the others
            warnings.append(f"image {e['id']}: header unreadable ({type(ex).__name__}: {ex}); its bindings are not shown")
            continue
        mapping, flags = channel_map.resolve(names)
        b = {}
        for role, idx in mapping.items():
            if idx is None:
                b[role] = {"channel": None, "how": "unbound"}
            else:
                by_index = channel_map.match_by != "name" or f"channel:{role}:name_unmatched_used_index_{idx}" in flags
                b[role] = {"channel": names[idx] if idx < len(names) else None, "how": "index" if by_index else "name"}
        out[e["id"]] = b
        for f in flags:
            if f.endswith(":MISSING_REQUIRED"):
                problems.append(f"image {e['id']}: the required role {f.split(':')[1]} binds no channel")
            elif f.startswith("channel:role_collision:"):
                spec = f[len("channel:role_collision:"):]
                a, b2, idx = spec.split("=")
                i = int(idx.replace("idx", ""))
                problems.append(f"image {e['id']}: the roles {a} and {b2} bind the same channel {i} "
                                f"({names[i] if i < len(names) else '?'})")
    return out, problems, warnings


def plan_change(root, manifest: dict, *, edits=None, unset=(), roles=None, origin="form", revert=False,
                model_entry=None, save=None, image_id=None) -> dict:
    """The result of a change, computed on a copy of ``manifest`` and never written: ``{"manifest", "recipe_bytes",
    "channel_map_bytes", "recipe_text", "channel_map_text", "config", "preview", "problems", "warnings"}``
    (None for a file the project then lacks). A change that touches no recipe edit (run settings, stage switches)
    keeps recipe.yaml and channel_map.yaml byte for byte and loads them where they are. ``preview`` (with
    ``roles``) is what preview_roles returns."""
    m = copy.deepcopy(manifest)
    pack, recipe = resolve_pack(root, m)
    base = pack_recipe_config(root, pack, recipe)
    pack_map = pack_channel_map(base)
    if revert:
        cur, origins, cur_roles, models = {}, {}, {}, []
    else:
        cur = read_edits(root, m)
        origins = {k: v for k, v in (m["recipe"].get("origins") or {}).items() if k in cur}
        cur_roles = dict(m["recipe"].get("roles") or {})
        models = [dict(e) for e in m.get("models") or []]
    for key in unset or ():
        key = str(key)
        if key == CHANNEL_MAP_KEY and cur_roles:
            cur_roles = {}
            continue
        if key not in cur:
            sets = sorted(cur) + ([CHANNEL_MAP_KEY] if cur_roles else [])
            raise ProjectError(f"{key}: this project does not set it (it sets: {', '.join(sets) or 'nothing'})")
        del cur[key]
        origins.pop(key, None)
    if origin not in ORIGINS:
        raise ProjectError(f"unknown origin {origin!r}: one of {', '.join(ORIGINS)}")
    for key, value in (edits or {}).items():
        v = check_edit(origin, str(key), value, pack, recipe)
        if _same(v, base.get(str(key), _MISSING)):
            cur.pop(str(key), None)                      # the pack's own value: no edit
            origins.pop(str(key), None)
        else:
            cur[str(key)] = v
            origins[str(key)] = origin
    problems: list[str] = []
    if roles is not None:
        clean = {}
        for role, channel in dict(roles).items():
            role = str(role)
            if role not in pack_map["roles"]:
                problems.append(f"{role} is not a role of the pack's channel map (its roles: "
                                f"{', '.join(pack_map['roles'])})")
                continue
            if channel is None or str(channel) == "":
                continue                                 # no edit for this role
            if [str(channel)] == [str(n) for n in pack_map["roles"][role].get("names") or []]:
                continue                                 # the pack map already says exactly this
            clean[role] = str(channel)
        cur_roles = clean
    if model_entry is not None:
        models = [e for e in models if e.get("key") != model_entry["key"]] + [dict(model_entry)]
    models = [e for e in models if e.get("key") in cur]          # a model edit that went away takes its record
    if save:
        _check_save(root, m, pack, save)
    overrides = dict(m["recipe"].get("overrides") or {})
    touched = bool(edits) or bool(unset) or roles is not None or revert or model_entry is not None
    if touched:
        cm_text = render_channel_map(pack_map, cur_roles) if cur_roles else None
        recipe_text = render_recipe(m, cur, bool(cur_roles))
        recipe_bytes = None if recipe_text is None else recipe_text.encode("ascii")
        cm_bytes = None if cm_text is None else cm_text.encode("ascii")
        hand_map = cur.get(CHANNEL_MAP_KEY) if recipe_bytes is not None and cm_bytes is None else None
        cfg = staged_config(root, pack, recipe, recipe_bytes, cm_bytes, overrides, hand_map=hand_map)
    else:                       # the recipe files stay exactly as they are (a run setting never rewrites them)
        recipe_bytes = _file_bytes(root, RECIPE_NAME) if m["recipe"].get("file") else None
        cm_bytes = _file_bytes(root, CHANNEL_MAP_NAME) if recipe_bytes is not None and cur_roles else None
        cfg = files_config(root, pack, recipe, overrides, recipe_file=recipe_bytes is not None)
    preview, warnings = None, []
    if roles is not None:
        b, role_problems, warnings = bindings(root, m, cfg.channel_map, image_id=image_id)
        problems += role_problems
        preview = {"yaml": render_channel_map(pack_map, cur_roles), "bindings": b, "problems": list(problems)}
    m["recipe"].update({"file": RECIPE_NAME if recipe_bytes is not None else None, "expected_config_hash": cfg.hash,
                        "roles": cur_roles, "origins": {k: origins.get(k, "manual") for k in sorted(cur)}})
    m["models"] = models
    m["validated"] = validated_record(pack, recipe, edited=recipe_bytes is not None, overrides=overrides)

    def text(b):
        return None if b is None else b.decode("ascii", "replace")

    return {"manifest": m, "recipe_bytes": recipe_bytes, "channel_map_bytes": cm_bytes, "recipe_text": text(recipe_bytes),
            "channel_map_text": text(cm_bytes), "config": cfg, "preview": preview, "problems": problems,
            "warnings": warnings}


# ---- models -------------------------------------------------------------------------------------------------

def model_entry(root, model, card=None) -> dict:
    """The manifest record of a model file inside the project: ``{key, ref, sha256, card, card_sha256, source}``.
    Relative paths are relative to the project root. The card defaults to <file>.card.json next to the file (what
    germquant train-nuclei writes); a card that names another digest is refused."""
    from ..fsutil import long_path
    from ..models_registry import validate_card
    from ..provenance import file_sha256
    from .images import path_too_long

    root = Path(os.path.abspath(root))

    def inside(path, what):
        p = Path(path)
        full = Path(os.path.abspath(p if p.is_absolute() else root / p))
        try:
            return full, full.relative_to(root).as_posix()
        except ValueError:
            raise ProjectError(f"{path}: the {what} must be a file inside the project ({root}), such as "
                               f"training/<round>/models/<name>/<file>; copy it there first, or promote it into a "
                               f"pack with `germquant pack derive`") from None

    full, rel = inside(model, "model")
    too_long = path_too_long(full, "the model loader")
    if too_long is not None:
        raise ProjectError(f"{model}: {too_long}")
    if not full.is_file():
        raise ProjectError(f"{model}: model file not found ({full})")
    sha = file_sha256(full)
    card_full = None
    if card is not None:
        card_full, card_rel = inside(card, "card")
        if not os.path.isfile(long_path(card_full)):
            raise ProjectError(f"{card}: card not found ({card_full})")
    else:
        guess = full.with_name(full.name + ".card.json")
        if os.path.isfile(long_path(guess)):         # read here only, so through the long form (a deep project)
            card_full, card_rel = guess, guess.relative_to(root).as_posix()
    card_sha = None
    if card_full is not None:
        with open(long_path(card_full), "rb") as fh:
            raw = fh.read()
        try:
            data = json.loads(raw.decode("utf-8"))
        except ValueError as e:
            raise ProjectError(f"{card_full}: the card is not JSON ({e})") from None
        check = validate_card(data)
        if not check["ok"]:
            raise ProjectError(f"{card_full}: the card has problems: {'; '.join(check['problems'])}")
        if str(data.get("model_sha256", "")).lower() != sha:
            raise ProjectError(f"{card_full}: the card describes a file with sha256 {str(data.get('model_sha256'))[:12]}"
                               f"..., but {rel} has {sha[:12]}...: give the model's own card")
        card_sha = hashlib.sha256(raw.replace(b"\r\n", b"\n")).hexdigest()
    else:
        card_rel = None
    return {"key": MODEL_KEY, "ref": rel, "sha256": sha, "card": card_rel, "card_sha256": card_sha,
            "source": "project"}
