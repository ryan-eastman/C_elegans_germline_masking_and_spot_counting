"""Project folders: the manifest, recipe edits, images, models and runs of one study (docs/PROJECTS.md,
docs/DESIGN_system_packs.md, section 7).

A project folder holds germquant_project.yaml (the manifest; project/manifest.py), the pinned pack version copied
into packs/<id>/<version>/, and, only after an edit, recipe.yaml (and channel_map.yaml after a role edit;
project/recipe.py). Results go to results/<image_id>/, run logs to runs/, training rounds to training/<round>/.

Every function takes the project explicitly, as a Project or as its folder, and searches the project's packs/
first by itself, so a GUI calling from a worker thread never depends on a context it did not set (plugin need
I11). Every write takes the project lock (project/lock.py), reads the manifest fresh, writes atomically, and
updates the Project it was given; runs never write the manifest.

Importing this package loads nothing heavy (no torch, CuPy, SpotMAX, Cellpose, bioio or dask.array; gate G6,
rule R5), and neither does any function here on an .nd2 or TIFF image: the torch-free entry ``python -m
germquant.project`` runs them (new, add, edit, use-model, rounds, unlock). run_project, project_status,
collect_project and export_commands belong to unit P5-RUN.
"""
from __future__ import annotations

import logging
import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from ..fsutil import long_path
from .lock import ProjectLockedError, describe, lock_report, project_lock, read_lock, unlock_hint
from .manifest import (
    MANIFEST_NAME,
    RECIPE_NAME,
    ProjectError,
    atomic_write,
    new_manifest,
    read_manifest,
    remove_file,
    stamp,
    write_manifest,
)

log = logging.getLogger(__name__)

__all__ = [
    "Project", "ProjectError", "ProjectLockedError", "add_images", "break_lock", "collect_project",
    "export_commands", "list_rounds", "new_project", "new_round", "open_project", "preview_roles",
    "project_config", "project_status", "revert_recipe", "run_project", "save_project", "use_model",
    "write_recipe",
]


def _todo(name: str, unit: str) -> str:
    return f"germquant.project.{name} is not built yet: unit {unit} of tranche 5 Train B builds it"


@dataclass
class Project:
    """An open project: ``root`` (the folder, absolute), ``manifest`` (germquant_project.yaml as a dict),
    ``warnings`` (what open_project found), ``config_hash`` (the hash computed at open or after the last write),
    and the manifest's ``images`` and ``pack`` entries as attributes."""

    root: Path
    manifest: dict
    warnings: list[str] = field(default_factory=list)
    config_hash: str | None = None

    @property
    def images(self) -> list[dict]:
        return self.manifest.setdefault("images", [])

    @property
    def pack(self) -> dict:
        return self.manifest.setdefault("pack", {})


def _root(project) -> Path:
    if isinstance(project, Project):
        return project.root
    return Path(os.path.abspath(os.fspath(project)))


def _manifest(project) -> dict:
    return project.manifest if isinstance(project, Project) else read_manifest(_root(project))


def _as_project(project) -> Project:
    if isinstance(project, Project):
        return project
    root = _root(project)
    return Project(root=root, manifest=read_manifest(root))


# ---- create and open ------------------------------------------------------------------------------------------

def new_project(root, pack, recipe=None, images=(), *, name=None, plugin_version=None) -> Project:
    """Create a project in ``root`` (made if missing; it may hold other files, but not a project): copy the pinned
    pack version (``pack`` is ``id`` or ``id@version``; without a version the newest compatible one, which the
    project then pins) into <root>/packs/<id>/<version>/ through the install checks, write the manifest and no
    recipe.yaml, and add ``images`` (files, folders, globs; refusals become warnings). FileNotFoundError for an
    unknown pack or recipe, ProjectError when ``root`` is already a project."""
    from ..packs import find_pack, sums_sha256, vendor
    from .recipe import resolve_pack, staged_config, validated_record

    root = Path(os.path.abspath(os.fspath(root)))
    if (root / MANIFEST_NAME).exists():
        raise ProjectError(f"{root} is already a germquant project ({MANIFEST_NAME}); open it, or choose another folder")
    if root.name.lower() == "config" or (root.name.lower() == "profiles" and root.parent.name.lower() == "config"):
        raise ProjectError(f"{root}: a project folder may not be named config (or config/profiles): load_config would "
                           f"resolve the project's recipe.yaml against the folder above it")
    p = pack if hasattr(pack, "recipes") else find_pack(str(pack))
    r = p.recipe(recipe)
    os.makedirs(root, exist_ok=True)
    with project_lock(root, "new"):
        if (root / MANIFEST_NAME).exists():      # again under the lock: another `new` may have finished meanwhile
            raise ProjectError(f"{root} is already a germquant project ({MANIFEST_NAME}); open it, or choose another "
                               f"folder")
        made_packs = not (root / "packs").exists()
        try:
            folder = vendor(p.id, p.version, root / "packs")
            entry = {"id": p.id, "version": p.version, "recipe": r.id, "sha256": sums_sha256(folder)}
            vp, vr = resolve_pack(root, {"pack": entry}, refresh=True)
            cfg = staged_config(root, vp, vr, None, None, {})
            m = new_manifest(name=name if name is not None else root.name, pack=entry, expected_config_hash=cfg.hash,
                             validated=validated_record(vp, vr, edited=False, overrides={}),
                             plugin_version=plugin_version)
            write_manifest(root, m)
        except BaseException:
            if made_packs:
                shutil.rmtree(root / "packs", ignore_errors=True)
            raise
    proj = Project(root=root, manifest=m, config_hash=cfg.hash)
    if cfg.hash != r.expected_config_hash:
        proj.warnings.append(f"the pack recipe {p.ref}/{r.id} hashes to {cfg.hash}, not its expected_config_hash "
                             f"{r.expected_config_hash}: `germquant pack check {p.ref} --recipe {r.id}` says why")
    if images:
        res = add_images(proj, _paths(images), plugin_version=plugin_version)
        proj.warnings += [f"not added: {x['path']}: {x['reason']}" for x in res["refused"]]
    return proj


def _hash_change(root, m: dict, now: str, saved: str, pack, recipe, sha_changed: bool) -> str:
    from ..packs import parse_ref

    msg = f"the config_hash {now} differs from the {saved} this project was saved with"
    pinned = f"{m['pack']['id']}@{m['pack']['version']}"
    why = []
    if m["recipe"].get("file"):
        import yaml

        from ..fsutil import read_config_text

        try:
            ext = str((yaml.safe_load(read_config_text(Path(root) / m["recipe"]["file"])) or {}).get("extends", ""))
        except Exception:  # noqa: BLE001 - the warning is built from what can be read
            ext = ""
        if ext.startswith("pack:"):
            pid, ver, rid = parse_ref(ext[len("pack:"):])
            if (pid, ver, rid) != (m["pack"]["id"], m["pack"]["version"], m["pack"]["recipe"]):
                why.append(f"recipe.yaml extends {ext}, while the manifest pins {pinned} recipe {m['pack']['recipe']}")
    if sha_changed:
        why.append(f"the copy of {pinned} in use differs from the one the project was made with")
    if not m["recipe"].get("file") and not m["recipe"].get("overrides") and now != recipe.expected_config_hash:
        why.append(f"the pack recipe {pinned}/{recipe.id} now hashes to {now}, not its expected_config_hash "
                   f"{recipe.expected_config_hash}")
    return msg + (": " + "; ".join(why) if why else ": a recipe file was changed outside germquant")


def open_project(root) -> Project:
    """The project in ``root``, with its config hash recomputed. Warnings (never errors) name what differs from
    the manifest: the hash it was saved with (naming the pack version change), the project's copy of the pack
    missing or replaced, a held lock (stale after 12 h)."""
    from ..packs import sums_sha256
    from .recipe import resolve_pack

    root = _root(root)
    m = read_manifest(root)
    proj = Project(root=root, manifest=m)
    p = m["pack"]
    pinned = f"{p['id']}@{p['version']}"
    try:
        pack, recipe = resolve_pack(root, m, refresh=True)
    except FileNotFoundError as e:
        proj.warnings.append(f"the pinned pack {pinned} does not resolve: {e}")
        pack = recipe = None
    if pack is not None:
        if not (root / "packs" / p["id"] / p["version"] / "pack.yaml").is_file():
            proj.warnings.append(f"the project's copy of the pack {pinned} is missing (packs/{p['id']}/{p['version']}/);"
                                 f" it resolves to {pack.folder} ({pack.source}) instead")
        sha = sums_sha256(pack.folder)
        sha_changed = bool(p.get("sha256")) and sha != p.get("sha256")
        if sha_changed:
            proj.warnings.append(f"the pack {pinned} in use ({pack.folder}) differs from the copy the project was made "
                                 f"with (pack_sha256 {str(sha)[:12]}... against {str(p['sha256'])[:12]}...)")
        try:
            proj.config_hash = project_config(proj).hash
        except Exception as e:  # noqa: BLE001 - reported, so a GUI can show the project and the fix
            proj.warnings.append(f"the project's recipe does not load: {type(e).__name__}: {e}")
        saved = m["recipe"].get("expected_config_hash")
        if proj.config_hash is not None and saved and proj.config_hash != saved:
            proj.warnings.append(_hash_change(root, m, proj.config_hash, saved, pack, recipe, sha_changed))
    if not m["recipe"].get("file") and os.path.isfile(long_path(root / RECIPE_NAME)):
        proj.warnings.append(f"{RECIPE_NAME} is in the project folder, but the manifest records no recipe edit, so runs "
                             f"do not use it, and a recipe edit is refused while it is there: move it out of the "
                             f"folder, then make its settings with `germquant project edit DIR --set KEY=VALUE`")
    holder = read_lock(root)
    if holder is not None:
        rep = lock_report(root)
        proj.warnings.append(f"the project is locked by {describe(holder)}"
                             + (f"; the lock is older than 12 h: {unlock_hint(root, holder)}" if rep["stale"] else ""))
    return proj


def project_config(project):
    """The Config the project runs with: recipe.yaml when the project was edited, else the pack recipe file itself
    (so the hash equals `germquant run --pack`), then the stage switches the project turns off. cfg.pack_record
    is the pack record of `--pack` runs plus ``edited`` (whether recipe.yaml exists)."""
    from ..config import load_config
    from ..packs import pack_record, recipe_profile_path, use_project
    from .recipe import resolve_pack

    root, m = _root(project), _manifest(project)
    pack, recipe = resolve_pack(root, m)
    with use_project(root):
        if m["recipe"].get("file"):
            p = root / str(m["recipe"]["file"])
            if not p.is_file():
                raise ProjectError(f"{p}: the manifest names it but it is missing; restore it, or go back to the pack "
                                   f"recipe with `germquant project edit DIR --revert`")
            cfg = load_config(p)
        else:
            cfg = load_config(recipe_profile_path(pack, recipe))
    for key, value in sorted((m["recipe"].get("overrides") or {}).items()):
        cfg.set(key, value)
    cfg.pack_record = {**pack_record(pack, recipe.id), "edited": bool(m["recipe"].get("file"))}
    return cfg


# ---- writes -----------------------------------------------------------------------------------------------------

def _unrecorded(root, name: str) -> str:
    """Why a change may not replace ``name``, a recipe file in the project folder that the manifest does not
    record (put there by hand, or copied from another project)."""
    from .manifest import CHANNEL_MAP_NAME

    if name == CHANNEL_MAP_NAME:
        return (f"{Path(root) / name} was not written by the role editor (the manifest records no role choice), and "
                f"this role edit would replace it: move it out of the project folder first (and if recipe.yaml names "
                f"it in io.channel_map, drop that with --unset io.channel_map)")
    return (f"{Path(root) / name} is not recorded in the manifest, so runs do not use it, and this edit would replace "
            f"it: move it out of the project folder first; its settings can then be made with `germquant project "
            f"edit DIR --set KEY=VALUE`")


def _apply(project, *, command: str, plugin_version=None, dry_run: bool = False, **change) -> dict:
    """Plan a change (project/recipe.py plan_change) and, unless ``dry_run``, write it under the lock: the recipe
    files, then the manifest with a new saved stamp. Returns the plan plus ``to_write`` and ``to_remove`` (the
    recipe files the change writes and removes), ``written`` and ``removed``.

    Only the recipe files the manifest records are ever removed: recipe.yaml when it records a recipe edit,
    channel_map.yaml when it records role choices. A write is refused (ProjectError) while a role edit has
    problems, and when it would replace a recipe file the manifest does not record. A dry run refuses neither: it
    returns the plan with the role problems in its preview and the replacement in its warnings, so the preview
    shows what a write refuses."""
    from .manifest import CHANNEL_MAP_NAME
    from .recipe import _file_bytes, plan_change

    root = _root(project)
    read_manifest(root)                          # not a project: refused before any lock file is made

    def plan_for(m, *, refuse: bool) -> dict:
        plan = plan_change(root, m, **change)
        if refuse and plan["problems"] and change.get("roles") is not None:
            raise ProjectError("the role edit is refused: " + "; ".join(plan["problems"]))
        recorded = {RECIPE_NAME: bool(m["recipe"].get("file")), CHANNEL_MAP_NAME: bool(m["recipe"].get("roles"))}
        files = ((CHANNEL_MAP_NAME, plan["channel_map_bytes"]), (RECIPE_NAME, plan["recipe_bytes"]))
        on_disk = {name: _file_bytes(root, name) for name, _ in files}
        # the map is written first, so recipe.yaml never names a missing map; recipe.yaml is removed first, so it
        # never names a removed map
        plan["to_write"] = [n for n, data in files if data is not None and on_disk[n] != data]
        plan["to_remove"] = [n for n, data in reversed(files) if data is None and recorded[n] and on_disk[n] is not None]
        for name in plan["to_write"]:
            if on_disk[name] is not None and not recorded[name]:
                if refuse:
                    raise ProjectError(_unrecorded(root, name))
                plan["warnings"].append(_unrecorded(root, name))
        return plan

    if dry_run:
        plan = plan_for(read_manifest(root), refuse=False)
        plan.update(written=[], removed=[])
        return plan
    written, removed = [], []
    with project_lock(root, command):
        plan = plan_for(read_manifest(root), refuse=True)
        new = plan["manifest"]
        new["saved"] = stamp(plugin_version)
        data = {CHANNEL_MAP_NAME: plan["channel_map_bytes"], RECIPE_NAME: plan["recipe_bytes"]}
        try:
            for name in plan["to_write"]:
                atomic_write(root / name, data[name])
                written.append(name)
            for name in plan["to_remove"]:
                if remove_file(root / name):
                    removed.append(name)
            write_manifest(root, new)
        except OSError as e:
            if not (written or removed):
                raise                                # nothing was changed: the project is as it was
            done = ", ".join([f"wrote {n}" for n in written] + [f"removed {n}" for n in removed])
            raise ProjectError(f"{e}; the write stopped part way ({done}), and the manifest still records the project "
                               f"as it was: run the same command again to finish it") from e
        written.append(MANIFEST_NAME)
    if isinstance(project, Project):
        project.manifest = new
        project.config_hash = plan["config"].hash
    plan.update(written=written, removed=removed)
    return plan


def save_project(project, *, plugin_version=None, device=None, interpreters=None, runner=None, xy_stride=None,
                 z_range=None, overrides=None, samples=None) -> Project:
    """Write the run settings and the other fields that are not recipe edits; None leaves a field as it is.
    ``device`` auto|cuda|mps|cpu ("" clears it: the engine's choice); ``runner`` subprocess|inprocess ("" clears);
    ``xy_stride`` 1 or more; ``z_range`` (z0, z1) (() clears); ``interpreters`` {recipe: python path, a record
    {python, germquant, git_sha}, or None to remove}; ``overrides`` the stage switches turned off,
    {"spots.enabled": False}, replacing the current set (they enter the config hash, as ``--no-spots`` does);
    ``samples`` the sample sheet's file in the project ("" clears)."""
    save = {k: v for k, v in (("device", device), ("interpreters", interpreters), ("runner", runner),
                              ("xy_stride", xy_stride), ("z_range", z_range), ("overrides", overrides),
                              ("samples", samples)) if v is not None}
    proj = _as_project(project)
    _apply(proj, command="save", plugin_version=plugin_version, save=save)
    return proj


def _add(project, paths, *, plugin_version=None) -> tuple[dict, list[dict], dict]:
    """add_images, plus the manifest records and the headers of the images it added (for `project add --json`)."""
    from .images import candidates, check_candidate, duplicate_reason, image_entry, input_suffixes

    root = _root(project)
    _manifest(project)                           # not a project: refused before any lock file is made
    try:
        cfg = project_config(project)
    except Exception as e:  # noqa: BLE001 - a recipe that does not load must not stop adding images
        log.warning("the project's recipe does not load (%s: %s): folders are searched for every supported suffix",
                    type(e).__name__, e)
        cfg = {}
    suffixes = input_suffixes(cfg)
    refused, checked = [], []
    for i, (given, path, refusal) in enumerate(candidates(root, _paths(paths), suffixes)):
        header, reason = (None, refusal) if path is None else check_candidate(root, path)
        if reason is not None:
            refused.append((i, {"path": given, "reason": reason}))
        else:
            checked.append((i, given, path, header))
    added, entries, headers = [], [], {}
    new = None
    if checked:
        with project_lock(root, "add"):
            m = read_manifest(root)
            for i, given, path, header in checked:
                entry = image_entry(root, path, header)
                why = duplicate_reason(entry["id"], m["images"])
                if why is not None:
                    refused.append((i, {"path": given, "reason": why}))
                    continue
                m["images"].append(entry)
                added.append(entry["id"])
                entries.append(entry)
                headers[entry["id"]] = header
            if added:
                m["saved"] = stamp(plugin_version)
                write_manifest(root, m)
                new = m
    if new is not None and isinstance(project, Project):
        project.manifest = new
    return ({"added": added, "refused": [r for _, r in sorted(refused, key=lambda x: x[0])]}, entries, headers)


def add_images(project, paths, *, plugin_version=None) -> dict:
    """Add files, folders and globs: ``{"added": [image_id], "refused": [{"path": str, "reason": str}]}``. A folder
    is searched for the recipe's io.input_suffixes, or every supported suffix when the recipe sets none; germquant
    outputs and the project's own folders are skipped. Refused: multi-scene files, an .nd2 with P > 1, duplicate
    ids, unreadable headers, and the rest of project/images.py; the other files are still added."""
    return _add(project, paths, plugin_version=plugin_version)[0]


def _paths(paths) -> list:
    """``paths`` as a list: one path (a str or os.PathLike) is a list of one, never a string split into letters."""
    if isinstance(paths, (str, os.PathLike)):
        return [paths]
    return list(paths)


def preview_roles(project, roles, *, image_id=None) -> dict:
    """``{"yaml", "bindings", "problems"}``. ``yaml``: the channel map these role choices give, the pack's map with
    each chosen role bound by its exact channel name; ``write_recipe(project, roles=roles)`` writes exactly this
    text as channel_map.yaml when a role differs from the pack map, and no map otherwise. ``bindings``: how each
    role binds on each image (``{image_id: {role: {"channel", "how": "name" | "index" | "unbound"}}}``, all images
    or ``image_id`` only; an image whose header cannot be read is left out). ``problems``: what makes write_recipe
    refuse (an unknown role, a required role that binds no channel, two roles on one channel). Nothing is
    written."""
    from .recipe import plan_change

    root = _root(project)
    plan = plan_change(root, read_manifest(root), roles=dict(roles), origin="roles", image_id=image_id)
    for w in plan["warnings"]:
        log.warning("%s", w)
    return plan["preview"]


def write_recipe(project, *, edits=None, unset=(), roles=None, origin="form", plugin_version=None) -> str:
    """Write recipe.yaml (and channel_map.yaml for ``roles``), validated by origin (form, tuning, model, roles;
    project/recipe.py), and return the new config_hash. ``edits`` {dotted key: value}; ``unset`` keys to drop
    (io.channel_map drops the role edits); ``roles`` {role: exact channel name}, replacing the role edits ({}
    goes back to the pack map). Once no edit is left recipe.yaml is removed and the hash is the pack recipe's."""
    plan = _apply(project, command="write_recipe", plugin_version=plugin_version, edits=dict(edits or {}),
                  unset=tuple(unset or ()), roles=None if roles is None else dict(roles), origin=origin)
    return plan["config"].hash


def revert_recipe(project, *, plugin_version=None) -> str:
    """Remove every recipe edit, role edit and project model: no recipe.yaml, no channel_map.yaml, and the hash of
    the pack recipe (with the project's stage switches, if any). Returns it."""
    plan = _apply(project, command="revert", plugin_version=plugin_version, revert=True)
    return plan["config"].hash


def use_model(project, model, *, card=None, plugin_version=None) -> str:
    """Name a model file of the project (training/<round>/models/<name>/<file>; relative to the project root, or
    absolute inside it) by relative path in recipe.yaml (segmentation.nuclei.cellpose_model), and record its sha256
    and its card (default <file>.card.json next to it) in the manifest's models. Returns the new config_hash."""
    from .recipe import MODEL_KEY, model_entry

    root = _root(project)
    _manifest(project)                           # not a project: refused first
    entry = model_entry(root, model, card)
    plan = _apply(project, command="use-model", plugin_version=plugin_version, edits={MODEL_KEY: entry["ref"]},
                  origin="model", model_entry=entry)
    if str(plan["config"].get("segmentation.nuclei.method", "auto")) == "classical":
        log.warning("the recipe segments with method classical, which loads no model: %s is recorded but not used",
                    entry["ref"])
    return plan["config"].hash


# ---- locks and rounds -------------------------------------------------------------------------------------------

def break_lock(project, *, host, pid) -> bool:
    """Remove the project's lock when its holder is ``host`` and ``pid`` (True); False when it is not locked;
    ProjectError naming the holder otherwise. A stale lock is never removed any other way."""
    from .lock import break_lock as _break

    return _break(_root(project), host=host, pid=pid)


def list_rounds(project) -> list[str]:
    """The training rounds under training/ (names matching ^[A-Za-z0-9_-]{1,40}$), sorted."""
    from .rounds import rounds_in

    return rounds_in(_root(project))


def new_round(project, name=None, *, plugin_version=None) -> Path:
    """Create training/<name>/ (default: the next free of r01, r02 and so on) and return it; the manifest's saved
    stamp records the write."""
    from .rounds import check_round_name, make_round, next_round_name

    root = _root(project)
    if name is not None:
        check_round_name(name)
    read_manifest(root)                          # not a project: refused before any lock file is made
    with project_lock(root, "rounds --new"):
        m = read_manifest(root)
        path = make_round(root, name if name is not None else next_round_name(root))
        m["saved"] = stamp(plugin_version)
        write_manifest(root, m)
    if isinstance(project, Project):
        project.manifest = m
    return path


# ---- unit P5-RUN ------------------------------------------------------------------------------------------------

def run_project(project, images=None, resume=True, device=None, shard=None, progress=None, *, via="cli",
                run_id=None, should_stop=None) -> dict:
    """Run the project's images: ``{"run_id", "log", "counts": {"done", "failed", "skipped"}, "exit_code"}``."""
    raise NotImplementedError(_todo("run_project", "P5-RUN"))


def project_status(project) -> dict:
    """The status record of PLAN.md section 9.3 (what `germquant project status --json` prints)."""
    raise NotImplementedError(_todo("project_status", "P5-RUN"))


def collect_project(project) -> dict:
    """Collect the per-image tables of results/ with the recipe's qc.germline_outlier_factor."""
    raise NotImplementedError(_todo("collect_project", "P5-RUN"))


def export_commands(project, *, shell="sh", shards=None) -> str:
    """The commands that run the project headless (sh or ps1), optionally split into shards."""
    raise NotImplementedError(_todo("export_commands", "P5-RUN"))
