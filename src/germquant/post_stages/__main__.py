"""``python -m germquant.post_stages list [--json] [--all]`` and ``write DIR --entries JSON [--json]``.

The torch-free entry (docs/POST_STAGES.md; plugin needs I12 and I13; decision D39). It imports neither
germquant.cli nor germquant.pipeline, so it never loads torch, CuPy, SpotMAX, Cellpose or a post-stage module:

* ``list`` prints the post-stages this environment offers, read from metadata and their declaration files (with
  ``--all`` also the installed ones that cannot be used, with the reason); ``--json`` prints what
  germquant.post_stages.list_post_stages returns.
* ``write DIR --entries JSON`` checks the entries (a JSON list of {name, params}) against the post-stages
  installed in this environment and writes DIR/post_stages.yaml (ASCII, LF, atomic). The plugin writes a
  project's file this way, in the engine environment, because only that environment has the post-stages.

Exit codes: 0 done; 2 refused before anything was written (malformed JSON, no such folder, an entry that cannot
run here), with every problem on stderr; 1 the file could not be written.
"""
from __future__ import annotations

import argparse
import json
import sys


def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="python -m germquant.post_stages",
                                 description="List installed germquant post-stages and write post_stages.yaml "
                                             "(docs/POST_STAGES.md). Loads no torch.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    pl = sub.add_parser("list", help="the post-stages this environment offers, from metadata only")
    pl.add_argument("--json", action="store_true", help="machine-readable output")
    pl.add_argument("--all", action="store_true", help="also list installed post-stages that cannot be used, "
                                                       "with the reason")
    pw = sub.add_parser("write", help="check entries against this environment's post-stages and write "
                                      "DIR/post_stages.yaml")
    pw.add_argument("dir", metavar="DIR", help="the folder to write post_stages.yaml into (a project root)")
    pw.add_argument("--entries", required=True, metavar="JSON",
                    help='a JSON list of {"name": ..., "params": {...}}; [] switches every post-stage off')
    pw.add_argument("--json", action="store_true", help="print the outcome as JSON")
    return ap


def _list(args) -> int:
    from . import list_post_stages

    rows = list_post_stages(include_unusable=args.all)
    if args.json:
        print(json.dumps(rows, indent=1), flush=True)
        return 0
    if not rows:
        print("no post-stages installed in this environment", flush=True)
    for r in rows:
        state = "" if r.get("usable", True) else "  [UNUSABLE]"
        print(f"  {r['name']} {r.get('version') or '?'}  ({r.get('dist') or r.get('source')}){state}", flush=True)
        if r.get("title"):
            print(f"      {r['title']}", flush=True)
        for prob in r.get("problems") or []:
            print(f"      problem: {prob}", flush=True)
    return 0


def _refuse(args, problems: list[str]) -> int:
    for p in problems:
        print(f"refused: {p}", file=sys.stderr, flush=True)
    if args.json:
        print(json.dumps({"written": False, "problems": problems}, indent=1), flush=True)
    return 2


def _write(args) -> int:
    import os
    from pathlib import Path

    from ..fsutil import long_path
    from . import POST_STAGES_FILE, load_post_stages_file, post_stages_sha256, validate_post_stages
    from .config import write_post_stages_file

    try:
        entries = json.loads(args.entries)
    except ValueError as e:
        return _refuse(args, [f"--entries is not JSON ({e})"])
    folder = Path(args.dir)
    if not os.path.isdir(long_path(folder)):         # as the writer writes: past MAX_PATH on Windows too
        return _refuse(args, [f"{folder} is not a folder"])
    problems = validate_post_stages(entries)
    if problems:
        return _refuse(args, problems)
    target = folder / POST_STAGES_FILE
    try:
        write_post_stages_file(target, entries)
    except (OSError, ValueError) as e:
        print(f"error: {target} could not be written ({type(e).__name__}: {e})", file=sys.stderr, flush=True)
        return 1
    if args.json:
        print(json.dumps({"written": True, "path": str(target), "sha256": post_stages_sha256(target),
                          "entries": load_post_stages_file(target)}, indent=1), flush=True)
    else:
        print(f"wrote {target}", flush=True)
    return 0


def main(argv=None) -> int:
    # a redirected stdout on Windows is in the ANSI code page (cp1252), which cannot encode every title or path a
    # listing prints: escape such a character instead of failing (--json output is ASCII already)
    reconfigure = getattr(sys.stdout, "reconfigure", None)
    if reconfigure is not None:
        reconfigure(errors="backslashreplace")
    args = _parser().parse_args(argv)
    return _list(args) if args.cmd == "list" else _write(args)


if __name__ == "__main__":
    sys.exit(main())
