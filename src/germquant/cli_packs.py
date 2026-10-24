"""`germquant pack ...` and the torch-free `python -m germquant.packs ...` (docs/PACKS.md).

  germquant pack list [--all] [--json] [--project DIR]
  germquant pack show ID[@VER] [--recipe R] [--json] [--project DIR]
  germquant pack check ID[@VER] [--recipe R] [--image FILE] [--deep] [--json] [--project DIR]
                                                                                  (exit 2 when unusable or unknown)
  germquant pack fetch-models ID[@VER] [--recipe R] [--from-file NAME=PATH ...] [--project DIR]
  germquant pack fetch-example ID[@VER] [--example NAME] [--out DIR] [--no-download] [--json]
  germquant pack install SRC [--replace]
  germquant pack verify ID|DIR
  germquant pack build-sums DIR
  germquant pack lint DIR [--previous ZIP]
  germquant pack derive --from ID@VER --recipe R --model PATH --card CARD --out DIR [--id NEW_ID] [--json]
                                                                                  (exit 1 when refused)
  germquant pack golden ID[@VER] [--recipe R] [--tier recipe|model] [--update]   (germquant only; exit 1 when a
                                                                                  golden fails or none ran)

`--project DIR` searches DIR/packs (a project's copy of its pinned pack and its derived packs) before the other
pack folders. This module imports nothing heavy at module level, so the same parsers serve cli.main (one line
there) and germquant.packs.__main__, whose commands never load torch (gate G6). `golden` runs the pipeline and is
offered by the full command line only; `check --deep` imports the pipeline first, as its docstring says.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

PROJECT_HELP = "search the project's packs folder (DIR/packs: its pinned pack and derived packs) first"


def add_pack_commands(psub, *, light: bool = False) -> None:
    """Add the pack subcommands to a subparsers object (``light``: leave out the ones that run the pipeline)."""
    pl = psub.add_parser("list", help="list the packs this engine finds (built-in, $GERMQUANT_PACKS, user folder)")
    pl.add_argument("--all", action="store_true", help="also list unusable and shadowed packs, with the reason")
    pl.add_argument("--json", action="store_true", help="machine-readable output")
    pl.add_argument("--project", default=None, metavar="DIR", help=PROJECT_HELP)

    ps = psub.add_parser("show", help="print a pack's card: recipes, models, citation")
    ps.add_argument("pack", metavar="ID[@VER]")
    ps.add_argument("--recipe", default=None)
    ps.add_argument("--json", action="store_true")
    ps.add_argument("--project", default=None, metavar="DIR", help=PROJECT_HELP)

    pc = psub.add_parser("check", help="can this recipe run here? exit code 2 when not, with every fix")
    pc.add_argument("pack", metavar="ID[@VER]")
    pc.add_argument("--recipe", default=None, help="a recipe of the pack (default: its default recipe)")
    pc.add_argument("--image", default=None, help="an image to check the readers extra and expected voxel size on")
    pc.add_argument("--deep", action="store_true",
                    help="also load the pipeline (torch first) and the backends, and report the device")
    pc.add_argument("--json", action="store_true")
    pc.add_argument("--project", default=None, metavar="DIR", help=PROJECT_HELP)

    pf = psub.add_parser("fetch-models", help="put the models a recipe needs into the model store (store, Cellpose "
                                              "cache, --from-file, then the pinned URL; every file checked)")
    pf.add_argument("pack", metavar="ID[@VER]")
    pf.add_argument("--recipe", default=None, help="only this recipe's models (default: every model of the pack)")
    pf.add_argument("--from-file", action="append", default=[], metavar="NAME=PATH",
                    help="a local copy of model NAME, checked by size and sha256 (repeatable)")
    pf.add_argument("--no-download", action="store_true", help="never download (the store, cache and files only)")
    pf.add_argument("--project", default=None, metavar="DIR", help=PROJECT_HELP)

    pe = psub.add_parser("fetch-example", help="fetch a pack's real example image (public data, checked by size, "
                                               "sha256 and pixels_sha256) and assemble it into one OME-TIFF")
    pe.add_argument("pack", metavar="ID[@VER]")
    pe.add_argument("--example", default=None, help="the example's id (default: the pack's first real example)")
    pe.add_argument("--out", default=None, help="the examples folder (default: $GERMQUANT_EXAMPLES, else the "
                                                "per-user folder); the image goes to <out>/<pack>/<example>/")
    pe.add_argument("--no-download", action="store_true", help="never download (use files already in --out)")
    pe.add_argument("--json", action="store_true", help="print the example's record as JSON")

    pi = psub.add_parser("install", help="install a pack folder or zip into the first $GERMQUANT_PACKS folder "
                                         "(else the per-user folder) as <id>/<version>/")
    pi.add_argument("src")
    pi.add_argument("--replace", action="store_true",
                    help="replace an installed copy of the same version with other content (pack developers only)")
    pi.add_argument("--target", default=None, help=argparse.SUPPRESS)

    pv = psub.add_parser("verify", help="check a pack's manifest and SHA256SUMS (and a root: checkout pack's "
                                        "engine files)")
    pv.add_argument("pack", metavar="ID|DIR")

    pb = psub.add_parser("build-sums", help="(re)write DIR/SHA256SUMS")
    pb.add_argument("dir")

    pn = psub.add_parser("lint", help="the pack rules: ASCII, hashes, sums, keys the engine reads, outside-recipe "
                                      "rules; --previous checks the version bump")
    pn.add_argument("dir")
    pn.add_argument("--previous", default=None, metavar="ZIP", help="the previous release (zip or folder)")
    pn.add_argument("--json", action="store_true")

    pd = psub.add_parser("derive", help="promote a model trained in a project into a derived pack: a recipe that "
                                        "extends the parent by versioned name and names the model, its portable "
                                        "card, and a copy of the parent (docs/PACKS.md, \"Derived packs\")")
    pd.add_argument("--from", dest="parent", required=True, metavar="ID@VER", help="the parent pack and its version")
    pd.add_argument("--recipe", required=True, help="the parent's recipe that the derived pack keeps")
    pd.add_argument("--model", required=True, metavar="PATH",
                    help="the model file (in a project: training/<round>/models/<name>/<file>)")
    pd.add_argument("--card", required=True, metavar="CARD",
                    help="its model card (<name>.card.json), written into the pack without machine paths")
    pd.add_argument("--out", required=True, metavar="DIR",
                    help="the packs folder to write into (<project>/packs, a $GERMQUANT_PACKS folder or the per-user "
                         "folder); the pack goes to DIR/<id>/0.1.0/ and the parent's copy beside it")
    pd.add_argument("--id", dest="new_id", default=None, metavar="NEW_ID",
                    help="the derived pack's id (default: <parent id>_derived_<first 8 hex of the model's sha256>)")
    pd.add_argument("--json", action="store_true", help="print the derived pack's record as JSON")

    if not light:
        pg = psub.add_parser("golden", help="run a recipe golden: the recipe tier on the pack's synthetic example "
                                            "with the classical segmenter and an empty model store")
        pg.add_argument("pack", metavar="ID[@VER]")
        pg.add_argument("--recipe", default=None, help="one recipe (default: every recipe; an outside recipe without "
                                                       "an expected folder fails)")
        pg.add_argument("--tier", choices=["recipe", "model"], default="recipe")
        pg.add_argument("--update", action="store_true",
                        help="write the expected tables instead of comparing (a pack maintainer's step)")


def add_parsers(sub) -> None:
    """The `germquant pack` group (cli.main calls this once)."""
    pk = sub.add_parser("pack", help="system packs: list, show, check, fetch-models, fetch-example, install, verify, "
                                     "build-sums, lint, derive, golden (docs/PACKS.md)")
    psub = pk.add_subparsers(dest="pack_cmd", required=True)
    add_pack_commands(psub)


def _print(*a) -> None:
    print(*a, flush=True)


def _list(args) -> int:
    from . import packs as P

    disc = P.discover()
    msg = P.no_builtin_root_message()
    if getattr(args, "json", False):
        _print(json.dumps({"version": P.version_info(), "builtin_root_found": msg is None, "message": msg,
                           "roots": [{"kind": k, "path": str(r)} for k, r in disc.roots],
                           "packs": P.list_packs(all=args.all)}, indent=1, default=str))
        return 0
    if msg:
        _print(msg)
    _print("search order: " + "; ".join(f"{k} {r}" for k, r in disc.roots))
    rows = P.list_packs(all=args.all)
    if not rows:
        _print("no packs found")
    for row in rows:
        state = "" if row["state"] == "ok" else f"  [{row['state'].upper()}]"
        recipes = ", ".join(r["id"] + ("*" if r.get("default") else "") for r in row.get("recipes", []))
        _print(f"  {row['id']} {row.get('version')}  {row.get('title', '')}  ({row['source']}){state}")
        if recipes:
            _print(f"      recipes: {recipes}")
        for prob in row.get("problems") or []:
            _print(f"      problem: {prob}")
    return 0


def _show(args) -> int:
    from . import packs as P

    p = P.find_pack(args.pack)
    recipes = [p.recipe(args.recipe)] if args.recipe else p.recipes
    if args.json:
        s = p.summary()
        s["recipes"] = [r.as_dict() for r in recipes]
        s["pack_sha256"] = P.sums_sha256(p.folder)
        _print(json.dumps(s, indent=1, default=str))
        return 0
    d = p.data
    _print(f"{p.ref}: {p.title}")
    if d.get("summary"):
        _print(f"  {d['summary']}")
    sysd = d.get("system") or {}
    if sysd:
        _print("  system: " + "; ".join(f"{k} {v}" for k, v in sysd.items()))
    _print(f"  object kind: {p.object_kind} (one row of the nuclei table is one {p.object_kind})")
    _print(f"  engine range: {(d.get('requires') or {}).get('germquant', 'any')}; root: {p.root}; source: {p.source} "
           f"({p.folder})")
    for r in recipes:
        flag = " (default)" if r is p.default_recipe() else ""
        _print(f"  recipe {r.id}{flag}: {r.title} [{r.dims}, {r.status}, config_hash {r.expected_config_hash}]")
        if r.need("models"):
            _print(f"      models: {', '.join(r.need('models'))}")
    for m in p.models:
        _print(f"  model {m.name} ({m.kind}, {m.size} bytes, sha256 {m.sha256[:12]}...)")
    for ex in P.list_examples(p):
        where = (f"{ex['image']}, in the pack" if ex["kind"] == "local"
                 else f"fetched with `germquant pack fetch-example {p.id} --example {ex['id']}`; licence "
                      f"{ex.get('licence')}; {ex.get('attribution')}")
        _print(f"  example {ex['id']} (for recipe {ex['recipe']}): {where}")
    cit = d.get("citation") or {}
    if cit.get("text"):
        _print(f"  cite: {cit['text']}")
    _print(f"  licence: {d.get('licence', 'not stated')}")
    return 0


def _check(args) -> int:
    from . import packs as P

    try:
        st = P.pack_status(args.pack, recipe=args.recipe, image=args.image, deep=args.deep)
    except FileNotFoundError as e:              # an unknown pack is unusable too: exit 2, as for any blocker
        if args.json:
            _print(json.dumps({"pack": args.pack, "recipe": args.recipe, "usable": False,
                               "blockers": [{"what": str(e), "fix": "germquant pack list"}], "warnings": []},
                              indent=1))
        else:
            print(f"ERROR: {e}\n      fix: germquant pack list", file=sys.stderr)
        return 2
    if args.json:
        _print(json.dumps(st, indent=1, default=str))
        return 0 if st["usable"] else 2
    head = f"{st['pack']}@{st['version']} recipe {st['recipe']}: " + ("usable" if st["usable"] else "NOT usable")
    _print(head)
    for b in st["blockers"]:
        _print(f"  BLOCKER: {b['what']}\n      fix: {b['fix']}")
    for w in st["warnings"]:
        _print(f"  warning: {w['what']}\n      {w['fix']}")
    for m in st["models"]:
        _print(f"  model {m['name']}: {m['state']} ({m['detail']})")
    if st.get("device"):
        _print(f"  device: {st['device']}")
    return 0 if st["usable"] else 2


def _parse_from_files(items: list[str]) -> dict[str, str]:
    out = {}
    for item in items or []:
        if "=" not in item:
            raise SystemExit(f"--from-file wants NAME=PATH, got {item!r}")
        name, path = item.split("=", 1)
        out[name.strip()] = path.strip()
    return out


def _fetch_models(args) -> int:
    from . import packs as P

    try:
        res = P.fetch_pack_models(args.pack, args.recipe, _parse_from_files(args.from_file),
                                  download=not args.no_download, say=_print)
    except P.PackFormatError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    missing = [r for r in res if r["status"] == "missing"]
    for r in res:
        _print(f"  {r['name']}: {r['status']} -> {r['path']}")
        for t in r["tried"] if r["status"] == "missing" else []:
            _print(f"      tried: {t}")
    if not res:
        _print("  nothing to fetch: the recipe needs no model")
    return 1 if missing else 0


def _fetch_example(args) -> int:
    from . import packs as P
    from .assemble import AssembleError
    from .datafetch import DownloadError, RangeNotSupported, ZipMemberError

    pack = P.find_pack(args.pack)
    say = (lambda *a: None) if args.json else _print
    try:
        rec = P.fetch_example(pack, args.example, args.out, download=not args.no_download, say=say)
    except (P.ExampleError, DownloadError, RangeNotSupported, ZipMemberError, AssembleError, OSError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    if args.json:
        _print(json.dumps(rec, indent=1, default=str))
        return 0
    _print(f"{pack.ref} example {rec['id']} ({rec['status']}): {rec['image']}")
    if rec.get("licence"):
        _print(f"  licence {rec['licence']}" + (f" ({rec['licence_url']})" if rec.get("licence_url") else "")
               + f"; {rec.get('attribution')}")
    recipe = f" --recipe {rec['recipe']}" if rec["recipe"] != pack.default_recipe().id else ""
    _print(f"  run it: germquant run \"{rec['image']}\" --pack {pack.ref}{recipe} --out RESULTS_FOLDER")
    if rec.get("expect"):
        _print("  expected: " + "; ".join(f"{k} in [{v[0]:g}, {v[1]:g}]" for k, v in rec["expect"].items()))
    return 0


def _install(args) -> int:
    from . import packs as P

    try:
        pack = P.install_pack(args.src, replace=args.replace, target=args.target)
    except (P.PackInstallError, FileNotFoundError) as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 1
    _print(f"installed {pack.ref} -> {pack.folder}")
    return 0


def _verify(args) -> int:
    from . import packs as P

    problems = P.verify_pack(args.pack)
    for prob in problems:
        _print(f"  {prob}")
    _print("OK: the pack matches its SHA256SUMS" if not problems else f"FAILED: {len(problems)} problem(s)")
    return 1 if problems else 0


def _build_sums(args) -> int:
    from . import packs as P

    p = P.build_sums(args.dir)
    n = len([ln for ln in p.read_text(encoding="ascii").splitlines() if ln.strip()])
    _print(f"wrote {p} ({n} files)")
    return 0


def _lint(args) -> int:
    from . import packs as P

    rep = P.lint_pack(args.dir, previous=args.previous)
    if args.json:
        _print(json.dumps(rep, indent=1))
        return 1 if rep["problems"] else 0
    for prob in rep["problems"]:
        _print(f"  problem: {prob}")
    for w in rep["warnings"]:
        _print(f"  warning: {w}")
    _print(f"{rep['pack']}: " + ("lint clean" if not rep["problems"] else f"{len(rep['problems'])} problem(s)"))
    return 1 if rep["problems"] else 0


def _derive(args) -> int:
    from . import packs as P

    try:
        rec = P.derive_pack(args.parent, args.recipe, args.model, args.card, args.out, pack_id=args.new_id)
    except (P.PackDeriveError, P.PackInstallError, P.PackFormatError, FileNotFoundError) as e:
        if args.json:
            _print(json.dumps({"ok": False, "error": str(e)}, indent=1))
        print(f"REFUSED: {e}", file=sys.stderr)
        return 1
    if args.json:
        _print(json.dumps(rec, indent=1))
        return 0
    m, par = rec["model"], rec["parent"]
    _print(("derived " if rec["status"] == "derived" else "already present, unchanged: ") + f"{rec['pack']} -> "
           f"{rec['folder']}")
    _print(f"  recipe {rec['recipe']}: extends {rec['extends']}, with the model {m['name']}; config_hash "
           f"{rec['config_hash']}")
    _print(f"  model in the store ({m['placed']}): {m['blob']}")
    _print(f"  parent {par['pack']}: " + ("copied to " if par["copied"] else "") + par["folder"])
    _print(f"  card {rec['card']}" + (f" (rewritten: {', '.join(rec['rewritten'])})" if rec["rewritten"] else ""))
    if rec["where"] == "project":
        _print(f"  a project's packs folder: list it with `germquant pack list --project "
               f"{os.path.dirname(os.path.abspath(args.out))}`")
    elif rec["where"] == "not searched":
        _print(f"  warning: {args.out} is not a folder this engine searches for packs: install the pack with "
               f"`germquant pack install {rec['folder']}`, or derive into a project's packs folder, a $GERMQUANT_PACKS "
               f"folder or the per-user folder")
    _print(f"  on another computer: germquant pack fetch-models {rec['id']} --from-file {m['name']}=PATH")
    return 0


def _golden(args) -> int:
    from . import packs as P
    from .packs_golden import run_pack_goldens

    pack = P.find_pack(args.pack)
    if args.tier == "model":
        print("REFUSED: the model tier runs a fetched public example on the RTX 5090 (CUDA, bfloat16) and compares "
              "it with example/expected/model_cuda/; it arrives with the model goldens of the 1.0.0 packs "
              "(docs/DESIGN_system_packs.md, 6.9). Use --tier recipe here.", file=sys.stderr)
        return 2
    results = run_pack_goldens(pack, recipe=args.recipe, update=args.update)
    bad = [r for r in results if not r["ok"]]
    for r in results:
        _print(f"  {pack.ref} {r['recipe']}: {'updated' if args.update else ('PASS' if r['ok'] else 'FAIL')}"
               + ("" if r["ok"] else f": {'; '.join(r['problems'][:5])}"))
    if not results:                  # nothing was checked: never a pass
        where = (" (its recipes are the engine's frozen config files; their goldens are tests/golden, run by "
                 "`pytest tests/test_golden.py`)" if pack.root == "checkout" else "")
        print(f"FAILED: {pack.ref} has no recipe golden (example/expected/recipe/<recipe>/) to run{where}",
              file=sys.stderr)
        return 1
    return 1 if bad else 0


HANDLERS = {"list": _list, "show": _show, "check": _check, "fetch-models": _fetch_models,
            "fetch-example": _fetch_example, "install": _install, "verify": _verify, "build-sums": _build_sums,
            "lint": _lint, "derive": _derive, "golden": _golden}


def dispatch(args) -> int:
    """Run the pack subcommand parsed into ``args`` (``args.pack_cmd``), with DIR/packs searched first when
    ``--project DIR`` is given (packs.use_project); FileNotFoundError is a clean exit 1."""
    try:
        project = getattr(args, "project", None)
        if project is None:
            return HANDLERS[args.pack_cmd](args)
        from .packs import use_project

        if not os.path.isdir(project):
            raise FileNotFoundError(f"--project {project}: no such folder")
        with use_project(project):
            return HANDLERS[args.pack_cmd](args)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


def main(argv: list[str] | None = None) -> int:
    """``python -m germquant.packs``: the torch-free entry (list, show, check, fetch-models, fetch-example,
    install, verify, build-sums, lint, derive, --version)."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--version" in argv:
        from . import packs as P

        if "--json" in argv:
            _print(json.dumps(P.version_info()))
        else:
            _print(P.version_line())
        return 0
    ap = argparse.ArgumentParser(prog="python -m germquant.packs", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="store_true", help="print the engine version, git sha and python (add "
                                                           "--json for JSON)")
    psub = ap.add_subparsers(dest="pack_cmd", required=True)
    add_pack_commands(psub, light=True)
    args = ap.parse_args(argv)
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass
    import logging

    logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(name)s: %(message)s")
    return dispatch(args)
