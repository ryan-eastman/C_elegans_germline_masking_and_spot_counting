"""germquant command line.

  germquant info  STACK.nd2
  germquant run   STACK.nd2 --config config/config.yaml --out results/ [--xy-stride 4 --z-range 20 40]
  germquant batch /nas/folder --config config/config.yaml --out /nas/folder_results
"""
from __future__ import annotations

import argparse
import fnmatch
import logging
import sys
from pathlib import Path

from . import provenance
from .config import load_config
from .io import read_nd2_metadata


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="germquant", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    pi = sub.add_parser("info", help="print .nd2 channels + voxel size (no pixels loaded)")
    pi.add_argument("nd2")

    pr = sub.add_parser("run", help="process a single .nd2")
    pr.add_argument("nd2")
    pr.add_argument("--config", required=True)
    pr.add_argument("--out", required=True)
    pr.add_argument("--xy-stride", type=int, default=1, help="downsample xy for a quick test")
    pr.add_argument("--z-range", type=int, nargs=2, default=None, metavar=("Z0", "Z1"))
    pr.add_argument("--no-spots", action="store_true",
                    help="segmentation only: skip RAD-51/SpotMAX spot detection (fast, never wedges)")
    pr.add_argument("--no-coloc", action="store_true",
                    help="skip PGL-1 granule surfacing + SYP<->PGL-1 colocalization stage")

    pb = sub.add_parser("batch", help="process every .nd2 under a folder, mirroring the tree")
    pb.add_argument("folder")
    pb.add_argument("--config", required=True)
    pb.add_argument("--out", required=True)
    pb.add_argument("--xy-stride", type=int, default=1)
    pb.add_argument("--no-spots", action="store_true",
                    help="segmentation only: skip RAD-51/SpotMAX spot detection (fast, never wedges)")
    pb.add_argument("--no-coloc", action="store_true",
                    help="skip PGL-1 granule surfacing + SYP<->PGL-1 colocalization stage")

    pv = sub.add_parser("validate", help="compare pipeline output to hand-scored ground truth")
    pv.add_argument("--pred", help="pipeline CSV (counts/lengths mode)")
    pv.add_argument("--truth", help="ground-truth CSV")
    pv.add_argument("--key", nargs="+", default=["image_id", "nucleus_id"], help="join key columns")
    pv.add_argument("--pred-col")
    pv.add_argument("--truth-col")
    pv.add_argument("--seg-pred", help="predicted label image .tif (segmentation mode)")
    pv.add_argument("--seg-truth", help="ground-truth label image .tif")
    pv.add_argument("--iou", type=float, default=0.5)
    pv.add_argument("--out", default="validation")

    pp = sub.add_parser("prep-training", help="export DAPI z-slices from .nd2 for annotation")
    pp.add_argument("folder")
    pp.add_argument("--config", required=True)
    pp.add_argument("--out", required=True)
    pp.add_argument("--n-slices", type=int, default=3)
    pp.add_argument("--xy-stride", type=int, default=1)

    pf = sub.add_parser("finetune", help="fine-tune a germline Cellpose model (GPU)")
    pf.add_argument("labeled_dir")
    pf.add_argument("--out-model", required=True)
    pf.add_argument("--pretrained", default="cpsam")
    pf.add_argument("--epochs", type=int, default=100)
    pf.add_argument("--print-only", action="store_true", help="print the command, don't run")

    sub.add_parser("check-gpu", help="assert the GPU is a Blackwell sm_120 (RTX 5090) with cu128 torch")

    args = p.parse_args(argv)
    _force_utf8_stdio()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    return {
        "info": lambda: _info(args.nd2),
        "run": lambda: _run(args),
        "batch": lambda: _batch(args),
        "validate": lambda: _validate(args),
        "prep-training": lambda: _prep_training(args),
        "finetune": lambda: _finetune(args),
        "check-gpu": lambda: _check_gpu(),
    }[args.cmd]()


def _force_utf8_stdio() -> None:
    """Windows consoles default to cp1252; our status glyphs (✓ ⚠ ≥ µ) would raise
    UnicodeEncodeError *after* the work is done, reporting a success as a crash. Reconfigure
    stdout/stderr to UTF-8 (replacing anything truly unencodable) so output never aborts a run.
    """
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")  # py3.7+ TextIOWrapper
        except (AttributeError, ValueError):  # already wrapped / not reconfigurable
            pass


def _info(nd2: str) -> int:
    m = read_nd2_metadata(nd2)
    print(f"file       : {nd2}")
    print(f"sizes      : {m['sizes']}")
    print(f"dtype      : {m['dtype']}   2D={m['is_2d']}")
    dz, dy, dx = m["spacing"]
    print(f"voxel µm   : dz={dz:.4f} dy={dy:.4f} dx={dx:.4f}  (anisotropy z/xy={dz/dy:.2f})")
    print(f"channels   : {m['channel_names']}")
    return 0


def _run(args) -> int:
    from .pipeline import process_image

    cfg = load_config(args.config)
    if getattr(args, "no_spots", False):
        cfg.set("spots.enabled", False)
        print("segmentation only: skipping spot detection (--no-spots)")
    if getattr(args, "no_coloc", False):
        cfg.set("coloc.enabled", False)
        print("skipping PGL-1 granule surfacing + colocalization (--no-coloc)")
    out = Path(args.out)
    prov = provenance.write_manifest(out, config_hash=cfg.hash, config=cfg.as_dict())
    z_range = tuple(args.z_range) if args.z_range else None
    res = process_image(args.nd2, cfg, out, xy_stride=args.xy_stride, z_range=z_range, prov=prov)
    mark = "✓" if res["qc_pass"] else "⚠"
    print(f"\n{mark} {res['image_id']}: {res['n_nuclei']} nuclei, qc_pass={res['qc_pass']}")
    if res["qc_flags"]:
        print("  flags:", "; ".join(res["qc_flags"]))
    print(f"  -> {res['out_dir']}")
    # QC status is recorded in the output tables; a flagged image is NOT a process
    # failure (so unattended Snakemake batches complete). Only exceptions are errors.
    return 0


def _batch(args) -> int:
    import pandas as pd

    from .pipeline import process_image

    cfg = load_config(args.config)
    if getattr(args, "no_spots", False):
        cfg.set("spots.enabled", False)
        print("segmentation only: skipping spot detection (--no-spots)")
    if getattr(args, "no_coloc", False):
        cfg.set("coloc.enabled", False)
        print("skipping PGL-1 granule surfacing + colocalization (--no-coloc)")
    root = Path(args.folder)
    out_root = Path(args.out)
    glob = cfg.get("io.input_glob", "**/*.nd2")
    excludes = cfg.get("io.exclude_patterns", [])

    files = [
        f for f in sorted(root.glob(glob))
        if not any(fnmatch.fnmatch(f.name.lower(), pat.lower()) for pat in excludes)
    ]
    if not files:
        print(f"No .nd2 files matched {glob} under {root}", file=sys.stderr)
        return 1

    prov = provenance.write_manifest(out_root, config_hash=cfg.hash, config=cfg.as_dict(),
                                     extra={"n_files": len(files), "input_root": str(root)})
    print(f"Processing {len(files)} files -> {out_root}")
    summaries = []
    for i, f in enumerate(files, 1):
        rel = f.relative_to(root).parent
        out_dir = out_root / rel / f.stem
        print(f"[{i}/{len(files)}] {f.name}")
        try:
            res = process_image(f, cfg, out_dir, xy_stride=args.xy_stride, prov=prov)
            isum = res["tables"]["image_summary"]
            n_germ = int(isum["n_germline_nuclei"].iloc[0]) if "n_germline_nuclei" in isum else 0
            summaries.append({"image_id": res["image_id"], "n_nuclei": res["n_nuclei"],
                              "n_germline": n_germ, "qc_pass": res["qc_pass"],
                              "qc_flags": ";".join(res["qc_flags"]), "out_dir": res["out_dir"]})
        except Exception as e:  # noqa: BLE001
            logging.exception("FAILED %s", f.name)
            summaries.append({"image_id": f.stem, "n_nuclei": 0, "n_germline": 0, "qc_pass": False,
                              "qc_flags": f"EXCEPTION:{e}", "out_dir": str(out_dir)})

    # Framing QC: flag images whose germline count is a strong outlier vs the batch median (a robust,
    # threshold-free proxy for "two gonad arms / extra tissue / fuller distal capture in frame" — worth
    # eyeballing the montage; the axis/position readout for such gonads is unreliable). Spot COUNTS are
    # unaffected, so this is advisory, not a failure.
    germ = [s["n_germline"] for s in summaries if s["qc_pass"] and s["n_germline"] > 0]
    if len(germ) >= 4:
        import statistics
        med = statistics.median(germ)
        factor = float(cfg.get("qc.germline_outlier_factor", 1.8))
        for s in summaries:
            if med > 0 and s["n_germline"] > factor * med:
                flag = f"qc:germline_count_outlier_{s['n_germline']}_vs_median{med:.0f}_review_framing"
                s["qc_flags"] = f"{s['qc_flags']};{flag}" if s["qc_flags"] else flag

    from .fsutil import long_path
    pd.DataFrame(summaries).to_csv(long_path(out_root / "batch_summary.csv"), index=False)
    n_pass = sum(s["qc_pass"] for s in summaries)
    print(f"\nDone. {n_pass}/{len(files)} passed QC. Summary -> {out_root / 'batch_summary.csv'}")
    return 0


def _validate(args) -> int:
    import json

    import pandas as pd

    from .validate import bland_altman_plot, compare_table, segmentation_metrics

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if args.seg_pred and args.seg_truth:
        import tifffile

        pred = tifffile.imread(args.seg_pred)
        gt = tifffile.imread(args.seg_truth)
        m = segmentation_metrics(pred, gt, iou_threshold=args.iou)
        (out / "segmentation_metrics.json").write_text(json.dumps(m, indent=2))
        print(f"Segmentation @ IoU≥{args.iou}: F1={m['f1']:.3f}  precision={m['precision']:.3f}  "
              f"recall={m['recall']:.3f}  mean_IoU={m['mean_iou']:.3f}  (TP={m['tp']} FP={m['fp']} FN={m['fn']})")
        return 0

    if not (args.pred and args.truth and args.pred_col and args.truth_col):
        print("counts mode needs --pred --truth --pred-col --truth-col (and --key)", file=sys.stderr)
        return 1
    paired, stats = compare_table(
        pd.read_csv(args.pred), pd.read_csv(args.truth), args.key, args.pred_col, args.truth_col
    )
    paired.to_csv(out / "paired.csv", index=False)
    (out / "agreement.json").write_text(json.dumps(stats, indent=2))
    if stats.get("n", 0) >= 2:
        bland_altman_plot(paired["pred"], paired["truth"], out / "bland_altman.png",
                          title=f"{args.pred_col} vs {args.truth_col}")
    print(f"n={stats.get('n')}  CCC={stats.get('ccc', float('nan')):.3f}  "
          f"Pearson={stats.get('pearson_r', float('nan')):.3f}  bias={stats.get('bias_mean_diff', float('nan')):.2f}  "
          f"MAE={stats.get('mae', float('nan')):.2f}")
    print(f"  -> {out}")
    return 0


def _prep_training(args) -> int:
    from .segment.finetune import prep_training_data

    cfg = load_config(args.config)
    root = Path(args.folder)
    excludes = cfg.get("io.exclude_patterns", [])
    files = [f for f in sorted(root.rglob("*.nd2"))
             if not any(fnmatch.fnmatch(f.name.lower(), pat.lower()) for pat in excludes)]
    if not files:
        print(f"No .nd2 under {root}", file=sys.stderr)
        return 1
    written = prep_training_data(files, cfg.channel_map, args.out,
                                 n_slices=args.n_slices, xy_stride=args.xy_stride)
    print(f"Wrote {len(written)} DAPI slices to {args.out}. Annotate them in the Cellpose GUI "
          f"(see docs/ANNOTATION.md), then `germquant finetune`.")
    return 0


def _finetune(args) -> int:
    from .segment.finetune import finetune_cellpose

    finetune_cellpose(args.labeled_dir, args.out_model, pretrained=args.pretrained,
                      n_epochs=args.epochs, run=not args.print_only)
    return 0


def _check_gpu() -> int:
    """Assert the GPU is a Blackwell sm_120 (RTX 5090) on a cu128 torch wheel. Returns nonzero
    if torch/CUDA is missing or the capability isn't (12, 0) — wire into CI/containers before
    trusting a run (ARCHITECTURE.md §4)."""
    try:
        import torch
    except Exception as e:  # noqa: BLE001
        print(f"torch not importable ({e}); install germquant[gpu] on the 5090/HPC.", file=sys.stderr)
        return 1
    if not torch.cuda.is_available():
        print("CUDA not available to torch (CPU-only build or no GPU visible).", file=sys.stderr)
        return 1
    cap = torch.cuda.get_device_capability()
    name = torch.cuda.get_device_name(0)
    print(f"torch {torch.__version__}  device={name}  CUDA cap {tuple(cap)}")
    if tuple(cap) != (12, 0):
        print(f"WARNING: expected sm_120 (12, 0) for the RTX 5090; got {tuple(cap)}. "
              "If this isn't a 5090 that's fine; if it is, the torch wheel didn't match "
              "Blackwell — reinstall from the cu128 index (ARCHITECTURE.md §4).", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
