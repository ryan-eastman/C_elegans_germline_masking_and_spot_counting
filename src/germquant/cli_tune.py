"""`germquant tune-spots`: tune a project's spot keys on its truth points, and apply the result
(germquant.tune; docs/TUNING.md; design 6.8 and 10 step 9).

  germquant tune-spots --project DIR [--points annotations/] [--grid FILE] [--max-runs N] [--out DIR] [--json]
  germquant tune-spots --project DIR --apply tuning/<run> [--plugin-version V] [--json]

The first form runs the search and writes tuning/<run_id>/ (grid.yaml, trials.jsonl, ranked.csv, best.yaml
and run.json); the second writes the run's best.yaml into the project's recipe.yaml. Exit codes: 0 done,
1 the search ran but nothing could be scored, or a detection (or the write of --apply) failed, 2 refused
before any detection (the project, the grid, --max-runs, an image without a finished run, a truth-point file
that cannot be read, no truth point in an object of its run, two annotators for one image).
"""
from __future__ import annotations

import argparse
import json
import math
import sys


def add_parsers(sub) -> None:
    """The `germquant tune-spots` command (cli.main calls this once)."""
    p = sub.add_parser(
        "tune-spots", formatter_class=argparse.RawDescriptionHelpFormatter,
        help="tune the recipe's spot keys of a project on its truth points, or apply a tuning run "
             "(docs/TUNING.md)",
        description=__doc__)
    p.add_argument("--project", required=True, metavar="DIR", help="the project folder")
    # --points and --max-runs default to None, so that --apply refuses them however they are spelled
    p.add_argument("--points", default=None, metavar="DIR",
                   help="the folder of truth-point files <image_id>__points_<annotator>.csv, relative to the "
                        "project (default annotations/); one annotator per image")
    p.add_argument("--grid", default=None, metavar="FILE",
                   help="a YAML grid file of dotted spots keys and the values to try, relative to the project "
                        "(default: spots.effect_size_min for spotmax, spots.log.min_snr for log; "
                        "docs/TUNING.md)")
    p.add_argument("--max-runs", type=int, default=None, metavar="N",
                   help="refuse a grid that means more than N detections over all images (default 50)")
    p.add_argument("--out", default=None, metavar="DIR",
                   help="the run folder to write, relative to the project (default tuning/<run_id>/; it must "
                        "be absent or empty)")
    p.add_argument("--apply", default=None, metavar="RUN_DIR",
                   help="instead of tuning, write the best.yaml of this tuning run (for example "
                        "tuning/<run_id>, relative to the project) into the project's recipe.yaml")
    p.add_argument("--plugin-version", default=None, metavar="V",
                   help="with --apply: the napari-germquant version asking for this write (recorded in the "
                        "manifest's saved stamp)")
    p.add_argument("--json", action="store_true", help="print the result as JSON")


def _clean(obj):
    """NaN and infinities to None, for JSON."""
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    return obj


def _failed(what: str, e: BaseException) -> int:
    """Exit 1 for a failure after the checks: the traceback through logging, then one ERROR line."""
    import logging

    logging.getLogger("germquant.tune").error("%s failed", what, exc_info=e)
    print(f"ERROR: {type(e).__name__}: {e}", file=sys.stderr)
    return 1


def dispatch(args) -> int:
    from .tune.apply import _project_api as apply_api
    from .tune.grid import TuneRefused
    from .tune.spots import MAX_RUNS
    from .tune.spots import _project_api as tune_api

    if args.apply is not None:
        given = [opt for opt, used in (("--points", args.points is not None), ("--grid", args.grid is not None),
                                       ("--max-runs", args.max_runs is not None), ("--out", args.out is not None))
                 if used]
        if given:
            print(f"ERROR: --apply writes a finished run's best.yaml; it takes no {', '.join(given)}", file=sys.stderr)
            return 2
    elif args.plugin_version is not None:
        print("ERROR: --plugin-version goes with --apply (a search writes no manifest stamp)", file=sys.stderr)
        return 2
    try:
        open_project = (apply_api() if args.apply is not None else tune_api())[0]
        project = open_project(args.project)
    except Exception as e:  # noqa: BLE001 - any failure to open the project refuses the command
        print(f"ERROR: cannot open the project {args.project}: {type(e).__name__}: {e}", file=sys.stderr)
        return 2

    if args.apply is not None:
        from .tune import apply_tuning

        try:
            config_hash = apply_tuning(project, args.apply, plugin_version=args.plugin_version)
        except (TuneRefused, FileNotFoundError, ValueError) as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 2
        except Exception as e:  # noqa: BLE001 - a failed write: exit 1, the traceback through logging
            return _failed("tune-spots --apply", e)
        if args.json:
            print(json.dumps({"config_hash": config_hash, "run_dir": args.apply}, indent=1))
        else:
            print(f"applied {args.apply}: recipe.yaml written, config_hash {config_hash} (validated: false)")
        return 0

    import contextlib

    from .tune import tune_spots

    # with --json, stdout holds the one JSON document: what a backend prints during the search (SpotMAX's
    # progress report) goes to stderr, with the progress lines and the warnings
    quiet = contextlib.redirect_stdout(sys.stderr) if args.json else contextlib.nullcontext()
    try:
        with quiet:
            res = tune_spots(project, points="annotations" if args.points is None else args.points, grid=args.grid,
                             out=args.out, max_runs=MAX_RUNS if args.max_runs is None else args.max_runs)
    except (TuneRefused, FileNotFoundError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    except Exception as e:  # noqa: BLE001 - a detection failed after the search started: exit 1
        return _failed("tune-spots", e)
    for w in res.get("warnings") or []:
        print(f"WARNING: {w}", file=sys.stderr)
    if args.json:
        print(json.dumps(_clean(res), indent=1, ensure_ascii=True))
    else:
        print(f"tuning run {res['run_dir']}")
        for row in res["ranked"][:5]:
            keys = ", ".join(f"{k}={v}" for k, v in row.items() if k.startswith("spots."))
            print(f"  {row['rank']:>3}  {row['objective']} {row['objective_value']}  {row['tie_break']} "
                  f"{row['tie_break_value']}  {keys}")
        if res["best"] is not None:
            print(f"best: {res['best']['keys']} ({res['best']['backend']}); apply it with "
                  f"germquant tune-spots --project {args.project} --apply {res['run_dir']}")
    if res["best"] is None:
        print("ERROR: no parameter set could be scored; see the warnings and docs/TUNING.md", file=sys.stderr)
        return 1
    return 0
