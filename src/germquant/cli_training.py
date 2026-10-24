"""The training loop's commands that live outside cli.py (docs/TRAINING.md).

  germquant train-curve --pairs DIR --out DIR [--sizes 1,2,4,8,16,all] [--repeats 3] [--seed S] [--base B]
                        [--keep-models none|best] [--json]          (exit 2 when refused before any training)

train-curve trains nested, seeded subsets of the training groups of a pairs folder against one held-out set and
scores each model right after training (germquant.segment.training.train_curve). Its models train with
train-nuclei's defaults (150 epochs, split by series, 20 percent held out, at least 5 nuclei per slice, the GPU
when one is visible); the API takes those settings as keywords.

This module imports nothing heavy at module level; cli.main adds its parsers with add_parsers and dispatches
every name of COMMANDS to dispatch.
"""
from __future__ import annotations

import json
import sys

COMMANDS = ("train-curve",)


def add_parsers(sub) -> None:
    """The commands of this module (cli.main calls this once)."""
    pc = sub.add_parser("train-curve",
                        help="a learning curve: train-nuclei on nested, seeded subsets of the training groups of a "
                             "pairs folder, every model scored on the same held-out pairs right after training; "
                             "keeps every card with its metrics and sha256, plus curve.csv (docs/TRAINING.md)")
    pc.add_argument("--pairs", required=True, help="the pairs folder (image/, labels/, manifest.csv)")
    pc.add_argument("--out", required=True,
                    help="the curve folder: curve.csv, curve.json and models/n<size>_r<repeat>/ (cards, logs and "
                         "the weights --keep-models keeps)")
    pc.add_argument("--sizes", default="1,2,4,8,16,all",
                    help="subset sizes in training groups (source images, or series), comma-separated; all is "
                         "every training group, and a size above that is skipped (default 1,2,4,8,16,all)")
    pc.add_argument("--repeats", type=int, default=3,
                    help="seeded subset orders per size; each order gives nested subsets (default 3)")
    pc.add_argument("--seed", type=int, default=0,
                    help="seed of the held-out split (the one train-nuclei --seed makes) and of the subsets "
                         "(default 0)")
    pc.add_argument("--base", default="cpsam",
                    help="starting weights of every model, as train-nuclei --base: a registered model, a model "
                         "file, a pinned pack model such as cpsam_v2@0f1cc3f7, cpsam (default) or none")
    pc.add_argument("--keep-models", choices=["none", "best"], default="none",
                    help="which weights stay: none (default; a Cellpose-SAM fine-tune is about 1.2 GB) or the best "
                         "one by its score")
    pc.add_argument("--json", action="store_true", help="print the curve record as JSON")


def dispatch(args) -> int:
    return {"train-curve": _train_curve}[args.cmd](args)


def _train_curve(args) -> int:
    from .segment.training import parse_sizes, train_curve

    try:
        res = train_curve(args.pairs, args.out, sizes=parse_sizes(args.sizes), repeats=args.repeats,
                          seed=args.seed, base=args.base, keep_models=args.keep_models)
    except (ValueError, FileNotFoundError) as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(res, indent=1, default=str))
        return 0
    print(f"\nlearning curve: {len(res['models'])} model(s), sizes {res['sizes']} of {len(res['train_pool'])} "
          f"training group(s), {len(res['val_groups'])} group(s) held out, scored by {res['scorer']}")
    if res["skipped_sizes"]:
        print(f"  skipped sizes {res['skipped_sizes']}: more than the training groups")
    for m in res["models"]:
        score = "failed: " + str(m["error"]) if m["status"] == "failed" else (
            "n/a" if m["score"] is None else f"{m['score']:.3f}")
        print(f"  size {m['size_label']} (repeat {m['repeat']}, {m['n_train']} pair(s) from {m['n_groups']} "
              f"group(s)): score {score}{'  [weights kept]' if m['weights_kept'] else ''}")
    best = res["best"]
    if best:
        print(f"  best: size {best['size_label']}, repeat {best['repeat']}, score {best['score']:.3f}")
    print(f"  curve -> {res['curve_csv']}")
    return 0
