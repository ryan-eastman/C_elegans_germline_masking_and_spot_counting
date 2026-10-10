"""Benchmark our 3D nucleus segmentation against the expert's partial label image over a manifest.

    python scripts/benchmark_nuclei3d.py --manifest benchmark_manifest.json --out benchmark_nuclei3d_out
        [--split HS] [--coverage auto|roi|bbox|none] [--dilate-um 2.0] [--majority-frac 0.5]
        [--max-dense-cells N] [--n-boot 1000] [--seed 0]

The manifest is the one scripts/benchmark_manifest.py builds (tag, split, results_dir, truth_dir, nd2,
ims per entry). Every gonad needs <image>__nuclei_labels.tif in results_dir and the truth label image
<stem>__imaris_labels.tif in truth_dir (the truth file contract of the Imaris reader package); a gonad
without either is skipped and listed under "skipped" in pooled.json. The ROI box is read from
<stem>__imaris_params.json when present and used as the coverage region (--coverage auto); "bbox"
forces the dilated truth bounding box, "none" counts every unmatched prediction as a false positive.

Outputs in --out: per_gonad.csv (one row per gonad: F1 at IoU 0.5, AP over IoU 0.5 to 0.95, splits,
merges, ignored predictions, z-extent error), per_nucleus.csv (one row per truth nucleus with its
matched prediction, IoU and z-extent error) and pooled.json (pooled scores with cluster-bootstrap
intervals by gonad, overall and by split). docs/BENCHMARK.md explains every column and key. Exit code
0 when at least one gonad was scored, 1 otherwise.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from germquant.benchmark import nuclei3d as n3  # noqa: E402
from germquant.benchmark.manifest import (  # noqa: E402
    find_pipeline_outputs,
    load_labels,
    load_manifest,
    load_truth,
    our_voxel_zyx,
    roi_box_from_params,
)
from germquant.benchmark.matching import DEFAULT_MAX_DENSE_CELLS  # noqa: E402
from germquant.benchmark.report import write_outputs  # noqa: E402


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="benchmark manifest JSON")
    ap.add_argument("--out", required=True, help="output folder")
    ap.add_argument("--split", default=None, help="only entries with this split")
    ap.add_argument("--coverage", choices=["auto", "roi", "bbox", "none"], default="auto",
                    help="coverage region for the partial-annotation rule")
    ap.add_argument("--dilate-um", type=float, default=n3.DEFAULT_DILATE_UM,
                    help="dilation of the truth bounding box (um) when it is the coverage region")
    ap.add_argument("--majority-frac", type=float, default=n3.DEFAULT_MAJORITY_FRAC,
                    help="volume fraction that makes an object belong to another (splits and merges)")
    ap.add_argument("--max-dense-cells", type=int, default=DEFAULT_MAX_DENSE_CELLS,
                    help="largest dense matrix (cells) solved exactly per connected component")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    entries = load_manifest(args.manifest)
    if args.split is not None:
        entries = [e for e in entries if e.split == args.split]
    rows, per_truth, skipped = [], [], []
    for e in entries:
        t0 = time.time()
        try:
            outputs = find_pipeline_outputs(e.results_dir)
        except (FileNotFoundError, ValueError) as err:
            skipped.append({"tag": e.tag, "reason": f"pipeline outputs: {err}"})
            print(f"  skip {e.tag}: {err}")
            continue
        if outputs.labels_path is None:
            skipped.append({"tag": e.tag, "reason": "no __nuclei_labels.tif in results_dir"})
            print(f"  skip {e.tag}: no nucleus label image")
            continue
        truth = load_truth(e.truth_dir, e.truth_stem)
        if truth is None or not truth.has_labels:
            skipped.append({"tag": e.tag, "reason": f"no truth label image for stem {e.truth_stem!r} in {e.truth_dir}"})
            print(f"  skip {e.tag}: no truth label image")
            continue
        roi_box = roi_box_from_params(truth.params)
        if args.coverage == "roi" and roi_box is None:
            skipped.append({"tag": e.tag, "reason": "--coverage roi but the params JSON has no ROI box"})
            print(f"  skip {e.tag}: no ROI box")
            continue
        voxel, vsrc = our_voxel_zyx(outputs, e, truth)
        pred = load_labels(outputs.labels_path)
        tl = load_labels(truth.labels_path)
        if pred.shape != tl.shape:
            skipped.append({"tag": e.tag, "reason": f"label images differ in shape: ours {pred.shape}, truth {tl.shape}"})
            print(f"  skip {e.tag}: shape mismatch {pred.shape} vs {tl.shape}")
            continue
        row, pt = n3.benchmark_gonad(pred, tl, voxel, tag=e.tag, split=e.split, roi_box=roi_box,
                                     coverage=args.coverage, dilate_um=args.dilate_um,
                                     majority_frac=args.majority_frac, max_dense_cells=args.max_dense_cells)
        row["voxel_source"] = vsrc
        row["roi_box"] = "" if roi_box is None else (
            " ".join(str(v) for v in roi_box["min_xyz"]) + " to " + " ".join(str(v) for v in roi_box["max_xyz"]))
        row["elapsed_s"] = round(time.time() - t0, 1)
        rows.append(row)
        per_truth.append(pt)
        print(f"  {e.tag}: pred {row['n_pred']} (ignored {row['n_pred_ignored_outside_coverage']}, "
              f"{row['coverage_rule']}) truth {row['n_truth']} | IoU0.5 TP/FP/FN {row['tp_iou50']}/"
              f"{row['fp_iou50']}/{row['fn_iou50']} F1={row['f1_iou50']:.3f} AP={row['ap_iou50_95']:.3f} | "
              f"splits {row['n_splits']} merges {row['n_merges']} | z-extent err "
              f"{row['z_extent_error_mean_um']:+.2f} um | {row['elapsed_s']} s")
        del pred, tl
    if not rows:
        print("no gonad could be scored")
        write_outputs(args.out, pd.DataFrame(), pd.DataFrame(), {"n_gonads": 0, "skipped": skipped})
        return 1
    per_gonad = pd.DataFrame(rows)
    per_nucleus = pd.concat(per_truth, ignore_index=True)
    pooled = n3.pool(per_gonad, per_nucleus, n_boot=args.n_boot, seed=args.seed)
    pooled["skipped"] = skipped
    pooled["manifest"] = str(args.manifest)
    pooled["settings"] = {"coverage": args.coverage, "dilate_um": args.dilate_um,
                          "majority_frac": args.majority_frac, "max_dense_cells": args.max_dense_cells,
                          "iou_thresholds": list(n3.IOU_THRESHOLDS)}
    paths = write_outputs(args.out, per_gonad, per_nucleus, pooled)
    f1, ap = pooled["f1_iou50"], pooled["ap_iou50_95"]
    print(f"POOLED {len(rows)} gonads: F1@0.5={f1['point']:.3f} [{f1['lo']:.3f}, {f1['hi']:.3f}] "
          f"AP={ap['point']:.3f} [{ap['lo']:.3f}, {ap['hi']:.3f}] splits {pooled['n_splits']} "
          f"merges {pooled['n_merges']} ignored {pooled['n_pred_ignored_outside_coverage']}")
    print("outputs ->", paths["pooled"].parent)
    return 0


if __name__ == "__main__":
    sys.exit(main())
