"""`germquant project ...` and the torch-free `python -m germquant.project ...` (docs/PROJECTS.md).

  germquant project new DIR --pack ID[@VER] [--recipe R] [--images GLOB ...] [--name N] [--plugin-version V] [--json]
  germquant project add DIR PATH|FOLDER|GLOB ... [--plugin-version V] [--json]
  germquant project edit DIR [--set KEY=VALUE ...] [--unset KEY ...] [--role ROLE=CHANNEL ...] [--revert]
                             [--device D] [--xy-stride N] [--z-range Z0 Z1] [--runner subprocess|inprocess]
                             [--interpreter RECIPE=PYTHON ...] [--samples FILE] [--no-<stage> ...]
                             [--plugin-version V] [--dry-run] [--json]
  germquant project use-model DIR --model training/<round>/models/<name>/<file> [--card CARD] [--plugin-version V]
                             [--json]
  germquant project rounds DIR [--new [NAME]] [--plugin-version V] [--json]
  germquant project unlock DIR --host H --pid P [--json]

Exit codes: 0 done; 1 done in part (`add` and `new --images` refused some files, each named with its reason);
2 refused before anything was written (no project, a locked project, an invalid edit, an unknown pack), with the
reason on stderr and, with --json, as {"error": ...} on stdout.

`edit` takes everything the project records: recipe edits (--set KEY=VALUE, a declared tunable of the pack
recipe; VALUE is read as YAML, so 3 is a number and true a boolean), --unset KEY (a recipe edit, a stage switch
or io.channel_map, which drops the role edits), --role ROLE=CHANNEL (bind a role by the image's exact channel
name; ROLE= drops that role's edit), --revert (back to the pack recipe), and the run settings (--device, the
word none clearing it; --xy-stride; --z-range Z0 Z1 or none; --runner; --interpreter RECIPE=PYTHON, PYTHON none
removing it; --samples FILE or none; --no-<stage>). --dry-run prints what would be written, the role preview
with the bindings per image and the new config_hash, and writes nothing. --plugin-version goes into the
manifest's stamps.

This module imports nothing heavy at module level, so the same parsers serve cli.main (one line there) and
germquant.project.__main__, whose commands never load torch (gate G6). The run, status, collect and
export-commands commands come with unit P5-RUN and are offered by the full command line only.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

LIGHT_COMMANDS = ("new", "add", "edit", "use-model", "rounds", "unlock")


def _common(p, *, plugin: bool = True) -> None:
    if plugin:
        p.add_argument("--plugin-version", default=None, metavar="V",
                       help="the napari-germquant version asking for this write (recorded in the manifest's stamps)")
    p.add_argument("--json", action="store_true", help="machine-readable output")


def add_project_commands(psub, *, light: bool = False) -> None:
    """Add the project subcommands to a subparsers object (``light``: only the torch-free writers)."""
    from .stages import cli_switches

    pn = psub.add_parser("new", help="create a project: copy the pinned pack version into DIR/packs/ and write the "
                                     "manifest (no recipe.yaml)")
    pn.add_argument("dir", metavar="DIR")
    pn.add_argument("--pack", required=True, metavar="ID[@VER]", help="the system pack (germquant pack list)")
    pn.add_argument("--recipe", default=None, help="one of the pack's recipes (default: its default recipe)")
    pn.add_argument("--images", nargs="+", default=[], metavar="PATH|FOLDER|GLOB", help="images to add at once")
    pn.add_argument("--name", default=None, help="the project's name (default: the folder's name)")
    _common(pn)

    pa = psub.add_parser("add", help="add image files, folders and glob patterns (refusals are named; the rest "
                                     "is added)")
    pa.add_argument("dir", metavar="DIR")
    pa.add_argument("paths", nargs="+", metavar="PATH|FOLDER|GLOB")
    _common(pa)

    pe = psub.add_parser("edit", help="edit the recipe (tunables, channel roles), the run settings and the stage "
                                      "switches; --dry-run previews")
    pe.add_argument("dir", metavar="DIR")
    pe.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="a tunable of the pack recipe")
    pe.add_argument("--unset", action="append", default=[], metavar="KEY",
                    help="drop a recipe edit, a stage switch, or io.channel_map (every role edit)")
    pe.add_argument("--role", action="append", default=[], metavar="ROLE=CHANNEL",
                    help="bind a role by an image's exact channel name (ROLE= drops that role's edit)")
    pe.add_argument("--revert", action="store_true", help="drop every recipe edit: back to the pack recipe")
    pe.add_argument("--device", default=None, choices=["auto", "cuda", "mps", "cpu", "none"])
    pe.add_argument("--xy-stride", type=int, default=None, metavar="N")
    pe.add_argument("--z-range", nargs="+", default=None, metavar="Z0 Z1", help="two planes, or the word none")
    pe.add_argument("--runner", default=None, choices=["subprocess", "inprocess", "none"])
    pe.add_argument("--interpreter", action="append", default=[], metavar="RECIPE=PYTHON",
                    help="the engine interpreter a recipe runs with (PYTHON none removes it)")
    pe.add_argument("--samples", default=None, metavar="FILE", help="the sample sheet in the project, or none")
    for st in cli_switches():
        pe.add_argument(f"--no-{st.name}", action="store_true", help=f"turn the {st.name} stage off in this project")
    pe.add_argument("--dry-run", action="store_true", help="print what would be written, and write nothing")
    _common(pe)

    pu = psub.add_parser("use-model", help="name a model of the project's training/ folder in recipe.yaml and "
                                           "record its digest and card")
    pu.add_argument("dir", metavar="DIR")
    pu.add_argument("--model", required=True, metavar="training/<round>/models/<name>/<file>")
    pu.add_argument("--card", default=None, help="its model card (default: <file>.card.json next to it)")
    _common(pu)

    pr = psub.add_parser("rounds", help="list the training rounds, or make a new one with --new")
    pr.add_argument("dir", metavar="DIR")
    pr.add_argument("--new", nargs="?", const="", default=None, metavar="NAME",
                    help="make training/<NAME>/ (default: the next of r01, r02, ...)")
    _common(pr)

    pk = psub.add_parser("unlock", help="remove a project lock, given its holder's host and pid")
    pk.add_argument("dir", metavar="DIR")
    pk.add_argument("--host", required=True)
    pk.add_argument("--pid", required=True, type=int)
    _common(pk, plugin=False)

    if not light:
        pass                       # unit P5-RUN adds run, status, collect and export-commands here


def add_parsers(sub) -> None:
    """The `germquant project` group (cli.main calls this once)."""
    pk = sub.add_parser("project", help="projects: new, add, edit, use-model, rounds, unlock (docs/PROJECTS.md)")
    psub = pk.add_subparsers(dest="project_cmd", required=True)
    add_project_commands(psub)


def _print(*a) -> None:
    print(*a, flush=True)


def _emit(args, payload: dict, lines: list[str]) -> None:
    if getattr(args, "json", False):
        _print(json.dumps(payload, indent=1, default=str))
    else:
        for line in lines:
            _print(line)


def _summary(proj) -> dict:
    m = proj.manifest
    return {"root": str(proj.root), "name": m.get("name"), "id": m.get("id"), "pack": m.get("pack"),
            "recipe": m.get("recipe"), "config_hash": proj.config_hash,
            "expected_config_hash": m["recipe"].get("expected_config_hash"), "edited": bool(m["recipe"].get("file")),
            "validated": m.get("validated"), "run": m.get("run"), "samples": m.get("samples"),
            "images": m.get("images"), "warnings": list(proj.warnings)}


def _new(args) -> int:
    from . import project as P

    proj = P.new_project(args.dir, args.pack, args.recipe, name=args.name, plugin_version=args.plugin_version)
    res = {"added": [], "refused": []}
    if args.images:
        res = P._add(proj, args.images, plugin_version=args.plugin_version)[0]
    out = {**_summary(proj), **res}
    lines = [f"project {proj.root}: {proj.pack['id']}@{proj.pack['version']} recipe {proj.pack['recipe']}, "
             f"config_hash {proj.config_hash}"] + [f"  warning: {w}" for w in proj.warnings]
    lines += [f"  added {i}" for i in res["added"]] + [f"  REFUSED {r['path']}: {r['reason']}" for r in res["refused"]]
    _emit(args, out, lines)
    return 1 if res["refused"] else 0


def _add(args) -> int:
    from . import project as P

    res, entries, headers = P._add(Path(args.dir), args.paths, plugin_version=args.plugin_version)
    out = {**res, "images": entries, "headers": headers}
    lines = [f"added {i}" for i in res["added"]] + [f"REFUSED {r['path']}: {r['reason']}" for r in res["refused"]]
    if not res["added"] and not res["refused"]:
        lines.append("nothing to add")
    _emit(args, out, lines)
    return 1 if res["refused"] else 0


def _pairs(items: list[str], what: str) -> list[tuple[str, str]]:
    out = []
    for item in items:
        if "=" not in item:
            raise ValueError(f"{what} wants {'KEY=VALUE' if what == '--set' else what[2:].upper() + '=...'}, got {item!r}")
        k, v = item.split("=", 1)
        out.append((k.strip(), v))
    return out


def _edit(args) -> int:
    import yaml

    from . import project as P
    from .project.manifest import CHANNEL_MAP_NAME, RECIPE_NAME, ProjectError, read_manifest
    from .project.recipe import _file_bytes
    from .stages import cli_switches

    root = Path(args.dir)
    if args.revert and (args.set or args.unset or args.role):
        raise ProjectError("--revert cannot be combined with --set, --unset or --role: it drops every recipe edit")
    m = read_manifest(root)
    edits = {}
    for k, v in _pairs(args.set, "--set"):
        try:
            edits[k] = yaml.safe_load(v) if v.strip() else None
        except yaml.YAMLError as e:
            raise ProjectError(f"--set {k}={v}: the value is not a YAML scalar ({str(e).splitlines()[0]})") from None
    roles = None
    if args.role:
        roles = dict(m["recipe"].get("roles") or {})
        for k, v in _pairs(args.role, "--role"):
            if v == "":
                roles.pop(k, None)
            else:
                roles[k] = v
    overrides = dict(m["recipe"].get("overrides") or {})
    over_changed = False
    for st in cli_switches():
        if getattr(args, f"no_{st.name.replace('-', '_')}", False):
            overrides[st.config_key] = False
            over_changed = True
    unset = []
    for key in args.unset:
        if key in overrides:
            overrides.pop(key)
            over_changed = True
        else:
            unset.append(key)
    save = {}
    if args.device is not None:
        save["device"] = "" if args.device == "none" else args.device
    if args.xy_stride is not None:
        save["xy_stride"] = args.xy_stride
    if args.z_range is not None:
        z = args.z_range
        if [s.lower() for s in z] == ["none"]:
            save["z_range"] = ()
        else:
            try:
                save["z_range"] = [int(s) for s in z]
            except ValueError:
                raise ProjectError(f"--z-range wants two planes Z0 Z1 or the word none, got {' '.join(z)}") from None
            if len(save["z_range"]) != 2:
                raise ProjectError(f"--z-range wants two planes Z0 Z1 or the word none, got {' '.join(z)}")
    if args.runner is not None:
        save["runner"] = "" if args.runner == "none" else args.runner
    if args.interpreter:
        save["interpreters"] = {k: (None if v.strip().lower() == "none" else v)
                                for k, v in _pairs(args.interpreter, "--interpreter")}
    if args.samples is not None:
        save["samples"] = "" if args.samples.strip().lower() == "none" else args.samples
    if over_changed:
        save["overrides"] = overrides
    before = {k: v for k, v in m.items() if k != "saved"}
    files_before = (_file_bytes(root, RECIPE_NAME), _file_bytes(root, CHANNEL_MAP_NAME))   # bytes, as written
    plan = P._apply(root, command="edit", plugin_version=args.plugin_version, dry_run=args.dry_run,
                    edits=edits or None, unset=tuple(unset), roles=roles, origin="form", revert=args.revert,
                    save=save or None)
    new = plan["manifest"]
    changed = ({k: v for k, v in new.items() if k != "saved"} != before
               or files_before != (plan["recipe_bytes"], plan["channel_map_bytes"]))
    cfg = plan["config"]
    out = {"dry_run": bool(args.dry_run), "changed": bool(changed), "config_hash": cfg.hash,
           "previous_config_hash": m["recipe"].get("expected_config_hash"), "edited": plan["recipe_text"] is not None,
           "validated": new["validated"], "recipe_yaml": plan["recipe_text"],
           "channel_map_yaml": plan["channel_map_text"], "preview": plan["preview"], "written": plan["written"],
           "removed": plan["removed"], "run": new["run"], "samples": new["samples"],
           "overrides": new["recipe"]["overrides"], "warnings": plan["warnings"]}
    head = "dry run, nothing written: " if args.dry_run else ""
    lines = [f"{head}config_hash {m['recipe'].get('expected_config_hash')} -> {cfg.hash}"
             + (" (edited)" if plan["recipe_text"] is not None else " (the pack recipe)")]
    if plan["preview"] is not None:
        for iid, b in plan["preview"]["bindings"].items():
            lines.append(f"  {iid}: " + ", ".join(f"{r}={v['channel']} ({v['how']})" for r, v in b.items()))
        lines += [f"  PROBLEM: {p}" for p in plan["preview"]["problems"]]
    lines += [f"  warning: {w}" for w in plan["warnings"]]
    lines += [f"  wrote {n}" for n in plan["written"]] + [f"  removed {n}" for n in plan["removed"]]
    _emit(args, out, lines)
    return 0


def _use_model(args) -> int:
    from . import project as P

    proj = P._as_project(Path(args.dir))
    h = P.use_model(proj, args.model, card=args.card, plugin_version=args.plugin_version)
    m = proj.manifest
    entry = next(e for e in m["models"] if e["key"] == "segmentation.nuclei.cellpose_model")
    _emit(args, {"config_hash": h, "model": entry, "edited": bool(m["recipe"].get("file")),
                 "validated": m["validated"]},
          [f"recipe.yaml now names {entry['ref']} (sha256 {entry['sha256'][:12]}...); config_hash {h}"])
    return 0


def _rounds(args) -> int:
    from . import project as P

    root = Path(args.dir)
    created = path = None
    if args.new is not None:
        p = P.new_round(root, args.new or None, plugin_version=args.plugin_version)
        created, path = p.name, str(p)
    else:
        P._as_project(root)                         # a folder that is not a project is refused
    rounds = P.list_rounds(root)
    _emit(args, {"rounds": rounds, "created": created, "path": path},
          ([f"made {path}"] if created else []) + [f"  {r}" for r in rounds] + ([] if rounds else ["no rounds"]))
    return 0


def _unlock(args) -> int:
    from .project import lock as L
    from .project.manifest import ProjectError, read_manifest

    root = Path(args.dir)
    read_manifest(root)                             # only a project's lock
    holder = L.read_lock(root)
    lock = None if holder is None else {"host": holder.get("host"), "pid": holder.get("pid"),
                                        "since": holder.get("time"), "command": holder.get("command")}
    stale = holder is not None and L.is_stale(root, holder)
    try:
        removed = L.break_lock(root, host=args.host, pid=args.pid)
    except ProjectError as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        _emit(args, {"removed": False, "lock": lock, "stale": stale, "error": str(e)}, [])
        return 2
    _emit(args, {"removed": removed, "lock": lock, "stale": stale},
          ["removed the lock of " + L.describe(holder) if removed else "the project is not locked"])
    return 0


HANDLERS = {"new": _new, "add": _add, "edit": _edit, "use-model": _use_model, "rounds": _rounds,
            "unlock": _unlock}


def dispatch(args) -> int:
    """Run the project subcommand parsed into ``args`` (``args.project_cmd``). A refusal (ProjectError, an unknown
    pack or file, a bad value, a folder that cannot be written) is exit code 2 with the reason on stderr, and as
    {"error": ...} on stdout with --json."""
    from .packs import PackInstallError
    from .project.manifest import ProjectError

    try:
        return HANDLERS[args.project_cmd](args)
    except (ProjectError, PackInstallError, ValueError, OSError) as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        if getattr(args, "json", False):
            _print(json.dumps({"error": str(e), "kind": type(e).__name__}))
        return 2


def main(argv: list[str] | None = None) -> int:
    """``python -m germquant.project``: the torch-free entry (new, add, edit, use-model, rounds, unlock)."""
    ap = argparse.ArgumentParser(prog="python -m germquant.project", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    psub = ap.add_subparsers(dest="project_cmd", required=True)
    add_project_commands(psub, light=True)
    args = ap.parse_args(sys.argv[1:] if argv is None else argv)
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass
    import logging

    logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(name)s: %(message)s")
    return dispatch(args)
