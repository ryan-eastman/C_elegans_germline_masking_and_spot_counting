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
  germquant project run DIR [--images ID ...] [--no-resume] [--device auto|cuda|mps|cpu] [--shard K/N]
                            [--run-id ID] [--via gui-subprocess|gui-inprocess|cli|snakemake] [--json]
  germquant project status DIR [--json]
  germquant project collect DIR [--json]
  germquant project export-commands DIR [--shell sh|ps1] [--shards N] [--out FILE] [--json]

Exit codes: 0 done; 1 done in part (`add` and `new --images` refused some files, each named with its reason;
`run`: an image failed, each named in the run's log, or the run stopped part way, with one ERROR line and, with
--json, {"error", "kind": "RunError", "run_id", "log"});
2 refused (no project, a locked project, an invalid edit, an unknown pack, a project file that is not YAML), with
nothing written, or a write that stopped part way (a file still denied after 2 s), whose message names what it
wrote: the same command again finishes it. The reason goes to stderr and, with --json, as {"error", "kind"} to
stdout.

`edit` takes everything the project records: recipe edits (--set KEY=VALUE, a declared tunable of the pack
recipe; VALUE is read as YAML, so 3 is a number and true a boolean), --unset KEY (a recipe edit, a stage switch
or io.channel_map, which drops the role edits), --role ROLE=CHANNEL (bind a role by the image's exact channel
name; ROLE= drops that role's edit), --revert (back to the pack recipe), and the run settings (--device, the
word none clearing it; --xy-stride; --z-range Z0 Z1 or none; --runner; --interpreter RECIPE=PYTHON, PYTHON none
removing it; --samples FILE or none; --no-<stage>). --dry-run prints what would be written, the role preview
with the bindings per image and the problems that make the edit itself refused (exit 0 all the same), and the new
config_hash, and writes nothing. --plugin-version goes into the manifest's stamps.

`run` runs the images into results/<image_id>/ and writes the run's log, runs/<run_id>.jsonl (germquant.project.run;
docs/PROJECTS.md, "Running a project"); `status` prints what each image needs; `collect` stacks the tables with
the recipe's QC settings; `export-commands` writes the commands that run the project headless. A refusal before any
image (a bad run id, shard or device, a model that is not the one the manifest records, an unusable pack, a lock
held by an edit) is exit 2, as for the commands above.

This module imports nothing heavy at module level, so the same parsers serve cli.main (one line there) and
germquant.project.__main__, whose commands never load torch (gate G6). run, status, collect and export-commands
are offered by the full command line only: run and status load the pipeline (torch first).
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
    pe.add_argument("--device", default=None, choices=["auto", "cuda", "mps", "cpu", "none"],
                    help="the device the project's runs ask for (auto = CUDA when available; none: the engine's "
                         "choice)")
    pe.add_argument("--xy-stride", type=int, default=None, metavar="N",
                    help="downsample xy by N for a quick run (1: full resolution)")
    pe.add_argument("--z-range", nargs="+", default=None, metavar="Z0 Z1",
                    help="two planes Z0 Z1 (0 <= Z0 < Z1), or the word none")
    pe.add_argument("--runner", default=None, choices=["subprocess", "inprocess", "none"],
                    help="how the plugin runs this project (none: its default)")
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
        _add_run_commands(psub)


def _add_run_commands(psub) -> None:
    """run, status, collect and export-commands: the full command line only (they load the pipeline)."""
    from .project.runlog import VIA

    prun = psub.add_parser("run", help="run the project's images into results/<image_id>/, one log per run under "
                                       "runs/ (exit 0 done, 1 an image failed, 2 refused)")
    prun.add_argument("dir", metavar="DIR")
    prun.add_argument("--images", nargs="+", default=None, metavar="ID",
                      help="run only these images (default: every image not excluded)")
    prun.add_argument("--no-resume", action="store_true", help="run every image again, also those already done")
    prun.add_argument("--device", default=None, choices=["auto", "cuda", "mps", "cpu"],
                      help="the device to ask for on this computer (default: the project's run.device, else the "
                           "engine's choice)")
    prun.add_argument("--shard", default=None, metavar="K/N",
                      help="run the K-th of N slices of the images (1 <= K <= N), for N jobs side by side")
    prun.add_argument("--run-id", default=None, metavar="ID",
                      help="the run's id, YYYYMMDDTHHMMSSZ_<8 hex> (default: a new one); names runs/<ID>.jsonl")
    prun.add_argument("--via", default="cli", choices=list(VIA), help="who started the run (recorded in its log)")
    prun.add_argument("--json", action="store_true", help="print the run's record as JSON (progress goes to stderr)")

    pst = psub.add_parser("status", help="each image's status (done, stale with the reason, failed, missing), the "
                                         "pack, models, lock and runs")
    pst.add_argument("dir", metavar="DIR")
    pst.add_argument("--json", action="store_true", help="machine-readable output (api.project_status)")

    pco = psub.add_parser("collect", help="stack the images' tables into results/batch_<table>.csv with the "
                                          "recipe's QC settings")
    pco.add_argument("dir", metavar="DIR")
    pco.add_argument("--json", action="store_true", help="machine-readable output")

    pex = psub.add_parser("export-commands", help="print (or write) the commands that run the project headless")
    pex.add_argument("dir", metavar="DIR")
    pex.add_argument("--shell", default="sh", choices=["sh", "ps1"],
                     help="sh (Linux, macOS, SLURM) or ps1 (Windows PowerShell); default sh")
    pex.add_argument("--shards", type=int, default=None, metavar="N", help="one run line per shard, --shard K/N")
    pex.add_argument("--out", default=None, metavar="FILE", help="write the commands to FILE instead of printing them")
    pex.add_argument("--json", action="store_true", help="print {text, out, shell, shards} as JSON")


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
    from .project.manifest import ProjectError, read_manifest
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
    plan = P._apply(root, command="edit", plugin_version=args.plugin_version, dry_run=args.dry_run,
                    edits=edits or None, unset=tuple(unset), roles=roles, origin="form", revert=args.revert,
                    save=save or None)
    new = plan["manifest"]
    # what the change writes: the manifest, and the recipe files whose bytes it replaces or removes
    changed = ({k: v for k, v in new.items() if k != "saved"} != before or bool(plan["to_write"])
               or bool(plan["to_remove"]))
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
        if plan["preview"]["problems"]:
            lines.append("  without --dry-run this edit is refused while a problem is left")
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
        _emit(args, {"removed": False, "lock": lock, "stale": stale, "error": str(e), "kind": type(e).__name__}, [])
        return 2
    _emit(args, {"removed": removed, "lock": lock, "stale": stale},
          ["removed the lock of " + L.describe(holder) if removed else "the project is not locked"])
    return 0


def _run_argv(args) -> list[str]:
    """The command line run_start records: the options given, in a fixed order."""
    argv = ["germquant", "project", "run", str(args.dir)]
    if args.images:
        argv += ["--images", *[str(i) for i in args.images]]
    if args.no_resume:
        argv.append("--no-resume")
    for opt, value in (("--device", args.device), ("--shard", args.shard), ("--run-id", args.run_id)):
        if value is not None:
            argv += [opt, str(value)]
    argv += ["--via", str(args.via)]
    if args.json:
        argv.append("--json")
    return argv


def _progress_printer(stream):
    state = {"n": 0, "i": 0}

    def show(ev: dict) -> None:
        kind = ev.get("event")
        if kind == "run_start":
            state["n"] = len(ev.get("images") or [])
            shard = f", shard {ev['shard']}" if ev.get("shard") else ""
            print(f"run {ev['run_id']}: {state['n']} image(s){shard}, config_hash {ev.get('config_hash')}, device "
                  f"{ev.get('device_selected')}", file=stream, flush=True)
        elif kind in ("image_done", "image_skipped", "image_failed"):
            state["i"] += 1
            what = {"image_done": f"done ({ev.get('elapsed_s')} s)", "image_skipped": "already done; skipped",
                    "image_failed": f"FAILED: {ev.get('error')}"}[kind]
            print(f"[{state['i']}/{state['n']}] {ev.get('image_id')}: {what}", file=stream, flush=True)
        elif kind == "run_end":
            c = ev.get("counts") or {}
            stop = "; interrupted" if ev.get("interrupted") else ""
            print(f"done {c.get('done', 0)}, skipped {c.get('skipped', 0)}, failed {c.get('failed', 0)}{stop}",
                  file=stream, flush=True)

    return show


def _run(args) -> int:
    from .project import run as R

    token = R.RUN_ARGV.set(_run_argv(args))
    try:
        res = R.run_project(Path(args.dir), images=args.images, resume=not args.no_resume, device=args.device,
                            shard=args.shard, progress=_progress_printer(sys.stderr if args.json else sys.stdout),
                            via=args.via, run_id=args.run_id)
    except R.RunError as e:                         # stopped after its start: an error, not a refusal
        print(f"ERROR: {e}", file=sys.stderr)
        if args.json:
            _print(json.dumps({"error": str(e), "kind": "RunError", "run_id": e.run_id, "log": e.log}))
        return 1
    finally:
        R.RUN_ARGV.reset(token)
    if args.json:                                   # the warnings were logged at the start, and are in the record
        _print(json.dumps(res, indent=1))
    else:
        _print(f"log: {res['log']}")
        if res["message"]:
            _print(res["message"])
    return int(res["exit_code"])


def _status(args) -> int:
    from .project.status import project_status

    st = project_status(Path(args.dir))
    lines = [f"project {st['root']}: {st['name']}, pack {st['pack']['pinned']} (installed: "
             f"{', '.join(st['pack']['installed']) or 'none'}; project copy: {'yes' if st['pack']['project_copy'] else 'no'}),"
             f" config_hash {st['config_hash']}" + (" (edited)" if st["edited"] else ""),
             f"validated: {st['validated'].get('value') if st['validated'] else None}"
             + (f" ({st['validated'].get('why')})" if st["validated"] else "")]
    for r in st["images"]:
        why = "" if r["status"] == "done" else f": {r['reason']}"
        qc = "" if r["qc_pass"] is None else f" (qc_pass {r['qc_pass']})"
        lines.append(f"  {r['id']}  {r['status']}{qc}{why}")
    lines += [f"  model {m['ref']}: " + ("not needed" if not m["needed"] else "verified" if m["verified"] else
                                        "MISSING" if not m["present"] else f"MISMATCH {m['mismatch']}")
              for m in st["models"]]
    lines += [f"warning: {w}" for w in st["warnings"]]
    lines.append(f"{len(st['runs'])} run(s) logged under runs/")
    _emit(args, st, lines)
    return 0


def _collect(args) -> int:
    from .project.collect import collect_project

    rec = collect_project(Path(args.dir))
    lines = [f"collected {rec['counts']} -> {rec['summary']} (qc.germline_outlier_factor {rec['outlier_factor']:g})"]
    if rec["sample_sheet"] is not None:
        joined = ", ".join(rec["sample_sheet"]["columns"]) or "none: the sheet has no column of its own"
        lines.append(f"sample sheet {rec['sample_sheet']['path']}: columns joined: {joined}")
    _emit(args, rec, lines)
    return 0


def _export_commands(args) -> int:
    from .project.commands import export_commands
    from .project.manifest import atomic_write

    text = export_commands(Path(args.dir), shell=args.shell, shards=args.shards)
    if args.out is not None:
        # PowerShell 5.1 reads a script without a byte-order mark in the ANSI code page, so ps1 gets one
        data = (b"\xef\xbb\xbf" if args.shell == "ps1" else b"") + text.encode("utf-8")
        atomic_write(Path(args.out), data)
    if args.json:
        _print(json.dumps({"text": text, "out": None if args.out is None else str(args.out), "shell": args.shell,
                           "shards": args.shards}))
    elif args.out is not None:
        _print(f"wrote {args.out}")
    else:
        sys.stdout.write(text)
        sys.stdout.flush()
    return 0


HANDLERS = {"new": _new, "add": _add, "edit": _edit, "use-model": _use_model, "rounds": _rounds,
            "unlock": _unlock, "run": _run, "status": _status, "collect": _collect,
            "export-commands": _export_commands}


def dispatch(args) -> int:
    """Run the project subcommand parsed into ``args`` (``args.project_cmd``). A refusal (ProjectError, an unknown
    pack or file, a bad value, a folder that cannot be written, a project file that is not YAML) is exit code 2
    with the reason on stderr, and as {"error": ...} on stdout with --json."""
    import yaml

    from .packs import PackInstallError
    from .project.manifest import ProjectError
    from .project.recipe import yaml_reason

    try:
        return HANDLERS[args.project_cmd](args)
    except (ProjectError, PackInstallError, ValueError, OSError, yaml.YAMLError) as e:
        msg = str(e)
        if isinstance(e, yaml.YAMLError):       # a hand-edited channel map, say: one line, and the way back
            msg = (f"a YAML file of the project does not load ({yaml_reason(e)}): if recipe.yaml or channel_map.yaml "
                   f"was edited by hand, fix it, or go back to the pack recipe with `germquant project edit DIR "
                   f"--revert`")
        print(f"REFUSED: {msg}", file=sys.stderr)
        if getattr(args, "json", False):
            _print(json.dumps({"error": msg, "kind": type(e).__name__}))
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
