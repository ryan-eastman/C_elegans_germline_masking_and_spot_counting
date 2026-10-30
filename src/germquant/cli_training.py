"""The training loop's commands that live outside cli.py (docs/TRAINING.md).

  germquant train-curve --pairs DIR --out DIR [--sizes 1,2,4,8,16,all] [--repeats 3] [--seed S] [--base B]
                        [--keep-models none|best] [--json]          (exit 2 when refused before any training)
  germquant import-labels SRC --format tiff|png|cellacdc-npz|cellpose-npy --image IMG --pairs DIR
                          [--allow-pickle] [--frame N] [--channel C] [--json]
                                    (exit 2 when refused before anything is written, 1 when a write fails)

train-curve trains nested, seeded subsets of the training groups of a pairs folder against one held-out set and
scores each model right after training (germquant.segment.training.train_curve). Its models train with
train-nuclei's defaults (150 epochs, split by series, 20 percent held out, at least 5 nuclei per slice, the GPU
when one is visible); the API takes those settings as keywords.

import-labels writes a label image made in another program, with the image it was drawn on, into a pairs folder
as export-training writes one (germquant.import_labels.import_labels): manifest rows with the source
import:<tool>. A refusal prints one REFUSED line on stderr and, with --json, {"error", "kind"} on stdout.

This module imports nothing heavy at module level; cli.main adds its parsers with add_parsers and dispatches
every name of COMMANDS to dispatch.
"""
from __future__ import annotations

import json
import sys

COMMANDS = ("train-curve", "import-labels")
IMPORT_FORMATS = ("tiff", "png", "cellacdc-npz", "cellpose-npy")     # germquant.import_labels.FORMATS


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
                    help="seeded orders of the training groups; each gives one model per size, its subsets "
                         "nested (default 3)")
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

    what = ("write a label image made in another program (the Cellpose GUI, Cell-ACDC, ilastik, Labkit, micro-sam) "
            "and the image it was drawn on into a pairs folder, as export-training writes one: image/ and labels/ "
            "TIFFs and manifest.csv rows with the source import:<tool>. Refuses a Cellpose _seg.npy without "
            "--allow-pickle and an image that a committed benchmark split holds as test data (docs/TRAINING.md)")
    pi = sub.add_parser("import-labels", help=what, description=what)
    pi.add_argument("src", help="the labels file: a TIFF or OME-TIFF, a PNG, a Cell-ACDC segm.npz or a Cellpose "
                                "_seg.npy")
    pi.add_argument("--format", required=True, choices=list(IMPORT_FORMATS),
                    help="the labels file's format: tiff (also OME-TIFF), png, cellacdc-npz or cellpose-npy (a "
                         "_seg.npy, read only with --allow-pickle)")
    pi.add_argument("--image", required=True,
                    help="the image the labels were drawn on, any format germquant reads; the labels hold its own "
                         "grid at full resolution, one plane or every plane")
    pi.add_argument("--pairs", required=True,
                    help="the pairs folder to write into (image/, labels/, manifest.csv); earlier rows of this image "
                         "are replaced, rows of other images stay")
    pi.add_argument("--allow-pickle", action="store_true",
                    help="read a Cellpose _seg.npy: it is a Python pickle, and loading it can run code, so allow "
                         "only files you made")
    pi.add_argument("--frame", type=int, default=0,
                    help="the frame of a T axis of the labels file and, when the image is a time series, its "
                         "timepoint (default 0); for a time series, a labels file without a T axis holds "
                         "the labels of that timepoint")
    pi.add_argument("--channel", default=None,
                    help="the image channel the labels were drawn on: its index or its name (default: the only "
                         "channel, else channel 0 with a warning)")
    pi.add_argument("--json", action="store_true",
                    help="print the import record as JSON; a refusal prints {\"error\", \"kind\"}")


def dispatch(args) -> int:
    return {"train-curve": _train_curve, "import-labels": _import_labels}[args.cmd](args)


def _import_labels(args) -> int:
    from .import_labels import import_labels

    channel = args.channel
    if isinstance(channel, str) and channel.strip().isdigit():
        channel = int(channel.strip())             # an index; any other text is a channel name
    try:
        res = import_labels(args.src, args.format, args.image, args.pairs, allow_pickle=args.allow_pickle,
                            frame=args.frame, channel=channel)
    except (ValueError, FileNotFoundError, ImportError) as e:
        # refused before anything was written: a bad file, image or option, a test-split image, a missing reader
        if args.json:
            print(json.dumps({"error": str(e), "kind": type(e).__name__}))
        print(f"REFUSED: {e}", file=sys.stderr)
        return 2
    except OSError as e:
        if args.json:
            print(json.dumps({"error": str(e), "kind": type(e).__name__}))
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    if args.json:
        print(json.dumps(res, indent=1, default=str))
        return 0
    print(f"\n{res['image_id']}: {res['n_pairs']} pair(s) imported ({res['dims']}, {res['source']}, channel "
          f"{res['channel_index']} {res['channel_name'] or 'unnamed'}, {res['n_objects']} object(s))")
    if res["replaced"]:
        print(f"  replaced {res['replaced']} earlier pair(s) of that image")
    print(f"  manifest -> {res['manifest']}")
    return 0


def _train_curve(args) -> int:
    from .models_registry import ModelPinError
    from .segment.training import parse_sizes, train_curve

    try:
        res = train_curve(args.pairs, args.out, sizes=parse_sizes(args.sizes), repeats=args.repeats,
                          seed=args.seed, base=args.base, keep_models=args.keep_models)
    except (ValueError, FileNotFoundError, ModelPinError) as e:
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
