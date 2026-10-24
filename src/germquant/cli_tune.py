"""`germquant tune-spots`: tune a project's spot keys on its truth points, and apply the result
(germquant.tune; docs/TUNING.md; design 6.8 and 10 step 9).

  germquant tune-spots --project DIR [--points annotations/] [--grid FILE] [--max-runs N] [--out DIR] [--json]
  germquant tune-spots --project DIR --apply tuning/<run> [--json]

The first form runs the search and writes tuning/<run_id>/ (grid.yaml, trials.jsonl, ranked.csv, best.yaml
and run.json); the second writes the run's best.yaml into the project's recipe.yaml. Exit codes: 0 done,
1 the search ran but nothing could be scored or a detection failed, 2 refused before any detection (the
project, the grid, --max-runs, an image without a finished run or without truth points).
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
    p.add_argument("--points", default="annotations", metavar="DIR",
                   help="the folder of truth-point files <image_id>__points_<annotator>.csv, relative to the "
                        "project (default annotations/); one annotator per image")
    p.add_argument("--grid", default=None, metavar="FILE",
                   help="a YAML grid of dotted spots keys and the values to try (default: "
                        "spots.effect_size_min for spotmax, spots.log.min_snr for log; docs/TUNING.md)")
    p.add_argument("--max-runs", type=int, default=50, metavar="N",
                   help="refuse a grid that means more than N detections over all images (default 50)")
    p.add_argument("--out", default=None, metavar="DIR",
                   help="the run folder to write (default tuning/<run_id>/ in the project; it must be absent "
                        "or empty)")
    p.add_argument("--apply", default=None, metavar="RUN_DIR",
                   help="instead of tuning, write the best.yaml of this tuning run (for example "
                        "tuning/<run_id>, relative to the project) into the project's recipe.yaml")
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


def dispatch(args) -> int:
    from .tune.apply import _project_api as apply_api
    from .tune.grid import TuneRefused
    from .tune.spots import _project_api as tune_api

    if args.apply is not None:
        given = [opt for opt, used in (("--grid", args.grid is not None), ("--out", args.out is not None),
                                       ("--points", args.points != "annotations"), ("--max-runs", args.max_runs != 50))
                 if used]
        if given:
            print(f"ERROR: --apply writes a finished run's best.yaml; it takes no {', '.join(given)}", file=sys.stderr)
            return 2
    open_project = (apply_api() if args.apply is not None else tune_api())[0]
    try:
        project = open_project(args.project)
    except Exception as e:  # noqa: BLE001 - any failure to open the project refuses the command
        print(f"ERROR: cannot open the project {args.project}: {type(e).__name__}: {e}", file=sys.stderr)
        return 2

    if args.apply is not None:
        from .tune import apply_tuning

        try:
            config_hash = apply_tuning(project, args.apply)
        except (TuneRefused, FileNotFoundError, ValueError) as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 2
        if args.json:
            print(json.dumps({"config_hash": config_hash, "run_dir": args.apply}, indent=1))
        else:
            print(f"applied {args.apply}: recipe.yaml written, config_hash {config_hash} (validated: false)")
        return 0

    from .tune import tune_spots

    try:
        res = tune_spots(project, points=args.points, grid=args.grid, out=args.out, max_runs=args.max_runs)
    except TuneRefused as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
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
