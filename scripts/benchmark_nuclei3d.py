"""Benchmark a 3D nucleus segmentation against the expert's partial label image over a manifest.

    python scripts/benchmark_nuclei3d.py --manifest benchmark_manifest.json --out benchmark_nuclei3d_out
        [--split HS] [--tags T1 T2] [--coverage auto|roi|bbox|none|protocol] [--dilate-um 2.0]
        [--majority-frac 0.5] [--max-dense-cells N] [--n-boot 1000] [--seed 0]
        [--pred-dir DIR | --model MODEL.pt --cache-root DIR [--device cuda]]

The manifest is the one scripts/benchmark_manifest.py builds (tag, split, results_dir, truth_dir, nd2,
ims per entry). Every gonad needs the truth label image <stem>__imaris_labels.tif in truth_dir (the truth
file contract of the Imaris reader package) and a predicted label image; a gonad without either is
skipped and listed under "skipped" in pooled.json.

Which prediction is scored:

* by default the pipeline's own <image>__nuclei_labels.tif in results_dir (for the frozen goldens, the
  deployed Cellpose-SAM labels);
* ``--pred-dir DIR``: DIR/<tag>__labels.tif, a label image written by an earlier ``--model`` run or by
  any other segmenter;
* ``--model MODEL.pt``: the EXPERIMENTAL nuclei3d model (germquant train-nuclei3d) run through the
  pipeline's own path (``segment_nuclei`` with backend nuclei3d and the size filter
  ``--min-volume-um3``, 4 um3 as in config/config.yaml) on the DAPI channel of the bench cache
  (``--cache-root``); the labels are written to <out>/labels/<tag>__labels.tif so they can be scored
  again with ``--pred-dir``. ``--fg-threshold``, ``--h-um`` and ``--neighbour-weight`` override the
  model file's instance parameters. It needs torch, and a GPU to be quick.

Coverage. ``--coverage protocol`` is the pre-registered rule of docs/BENCHMARK_PROTOCOL.md section 5.3
(truth fragments under ``--min-truth-voxels`` removed; predictions in scope when at least
``--scope-frac`` of their volume lies within ``--coverage-xy-um`` in xy and ``--coverage-z-planes``
planes in z of a truth mask; splits and merges by ``--split-frac`` of the truth volume). The default,
``auto``, is the rule this script had before the protocol: the ROI box from <stem>__imaris_params.json
when present, else the truth bounding box dilated by ``--dilate-um``; "bbox" forces the latter, "none"
counts every unmatched prediction as a false positive. germquant.benchmark.nuclei3d documents both.

Test gonads. Under the protocol the test gonads are scored once, by a frozen method. A test gonad (the
pre-registered split's, or one the manifest calls test) is skipped and listed under "skipped" unless
``--final`` is given with a ``--reason``.

Outputs in --out: per_gonad.csv (one row per gonad: F1 at IoU 0.5, AP over IoU 0.5 to 0.95, splits,
merges, ignored predictions, z-extent error, the coverage rule used), per_nucleus.csv (one row per truth
nucleus with its matched prediction, IoU and z-extent error) and pooled.json (pooled scores with
cluster-bootstrap intervals by gonad, overall and by split). docs/BENCHMARK.md explains every column
and key. Exit code 0 when at least one gonad was scored, 1 otherwise.
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
    ap.add_argument("--tags", nargs="*", default=None, help="only these tags")
    ap.add_argument("--group-by", choices=["split", "condition"], default="split",
                    help="what the per-gonad 'split' column (and the pooling by split) holds: the manifest's "
                         "split (default) or the gonad's condition (noHS, HS)")
    ap.add_argument("--coverage", choices=list(n3.COVERAGE_CHOICES), default="auto",
                    help="coverage region for the partial-annotation rule (protocol = docs/BENCHMARK_PROTOCOL.md 5.3)")
    ap.add_argument("--dilate-um", type=float, default=n3.DEFAULT_DILATE_UM,
                    help="dilation of the truth bounding box (um) when it is the coverage region")
    ap.add_argument("--majority-frac", type=float, default=n3.DEFAULT_MAJORITY_FRAC,
                    help="volume fraction that makes an object belong to another (splits and merges of the box rules)")
    ap.add_argument("--min-truth-voxels", type=int, default=n3.PROTOCOL_RULES["min_truth_voxels"],
                    help="protocol: truth fragments below this many voxels are removed")
    ap.add_argument("--coverage-xy-um", type=float, default=n3.PROTOCOL_RULES["coverage_xy_um"],
                    help="protocol: region C is the truth dilated by this in xy (um)")
    ap.add_argument("--coverage-z-planes", type=int, default=n3.PROTOCOL_RULES["coverage_z_planes"],
                    help="protocol: ... and by this many planes in z")
    ap.add_argument("--scope-frac", type=float, default=n3.PROTOCOL_RULES["scope_frac"],
                    help="protocol: a prediction is in scope with at least this fraction of its volume inside C")
    ap.add_argument("--split-frac", type=float, default=n3.PROTOCOL_RULES["split_frac"],
                    help="protocol: splits and merges count predictions covering this fraction of a truth nucleus")
    ap.add_argument("--max-dense-cells", type=int, default=DEFAULT_MAX_DENSE_CELLS,
                    help="largest dense matrix (cells) solved exactly per connected component")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--pred-dir", default=None, help="score DIR/<tag>__labels.tif instead of the pipeline's labels")
    src.add_argument("--model", default=None, help="segment with this nuclei3d model file and score the result")
    ap.add_argument("--cache-root", default=None, help="bench cache root for --model (default: the harness's)")
    ap.add_argument("--device", default="auto", help="--model: auto, cuda or cpu")
    ap.add_argument("--tile", default="96,384,384", help="--model: inference tile Z,Y,X")
    ap.add_argument("--min-volume-um3", type=float, default=4.0,
                    help="--model: the size filter (segmentation.nuclei.min_volume_um3, 4.0 in config/config.yaml)")
    ap.add_argument("--fg-threshold", type=float, default=None, help="--model: override the model's foreground threshold")
    ap.add_argument("--h-um", type=float, default=None, help="--model: override the model's h (um)")
    ap.add_argument("--neighbour-weight", type=float, default=None, help="--model: override the model's neighbour weight")
    ap.add_argument("--final", action="store_true", help="allow test gonads (needs --reason)")
    ap.add_argument("--reason", default=None, help="why test gonads are scored (with --final)")
    return ap.parse_args(argv)


def _model_labels(args, entry, out_dir: Path) -> tuple[object, dict]:
    """Segment the cached DAPI of ``entry`` with the nuclei3d model through segment_nuclei; write and
    return the labels and a record."""
    from germquant.benchmark import detector as bd
    from germquant.pipeline import _save_labels
    from germquant.segment import segment_nuclei
    from germquant.segment.nuclei3d.infer import parse_zyx

    cache = bd.load_cache(entry.tag, args.cache_root, arrays=("dna",), mmap=True)
    if cache.dna is None:
        raise FileNotFoundError(f"{entry.tag}: the bench cache holds no DNA channel")
    section = {"nuclei3d": {"model": str(args.model), "device": args.device, "tile_zyx": list(parse_zyx(args.tile)),
                            "fg_threshold": args.fg_threshold, "h_um": args.h_um,
                            "neighbour_weight": args.neighbour_weight}}
    t0 = time.time()
    labels, method = segment_nuclei(cache.dna, cache.spacing, method="cellpose", backend="nuclei3d",
                                    min_volume_um3=args.min_volume_um3, cfg_section=section)
    seg_s = time.time() - t0
    out = out_dir / "labels" / f"{entry.tag}__labels.tif"
    out.parent.mkdir(parents=True, exist_ok=True)
    _save_labels(labels, out, cache.spacing)
    return labels, {"segment_s": round(seg_s, 1), "labels_file": str(out), "segmentation_method": method}


def main(argv=None) -> int:
    from germquant.benchmark import detector as bd

    args = parse_args(argv)
    if args.final and not args.reason:
        print("--final needs --reason", file=sys.stderr)
        return 2
    entries = load_manifest(args.manifest)
    if args.split is not None:
        entries = [e for e in entries if e.split == args.split]
    if args.tags:
        want = set(args.tags)
        entries = [e for e in entries if e.tag in want]
    model_rec = None
    if args.model:
        from germquant import provenance

        model_rec = {"path": str(args.model), "sha256": provenance.file_sha256(args.model)}
    out_dir = Path(args.out)
    rows, per_truth, skipped = [], [], []
    for e in entries:
        role = e.extra.get("role") or e.split
        if bd.is_test_tag(e.tag, role) and not args.final:
            skipped.append({"tag": e.tag, "reason": "test gonad (docs/BENCHMARK_PROTOCOL.md): refused without --final"})
            print(f"  skip {e.tag}: test gonad, refused without --final")
            continue
        t0 = time.time()
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
        outputs = None
        try:
            outputs = find_pipeline_outputs(e.results_dir)
        except (FileNotFoundError, ValueError) as err:
            if not (args.pred_dir or args.model):
                skipped.append({"tag": e.tag, "reason": f"pipeline outputs: {err}"})
                print(f"  skip {e.tag}: {err}")
                continue
        extra = {}
        if args.model:
            try:
                pred, extra = _model_labels(args, e, out_dir)
            except FileNotFoundError as err:
                skipped.append({"tag": e.tag, "reason": f"model run: {err}"})
                print(f"  skip {e.tag}: {err}")
                continue
            source = "nuclei3d model"
        elif args.pred_dir:
            path = Path(args.pred_dir) / f"{e.tag}__labels.tif"
            if not path.is_file():
                skipped.append({"tag": e.tag, "reason": f"no {path.name} in {args.pred_dir}"})
                print(f"  skip {e.tag}: no predicted label image")
                continue
            pred = load_labels(path)
            source = f"pred_dir {path}"
        else:
            if outputs.labels_path is None:
                skipped.append({"tag": e.tag, "reason": "no __nuclei_labels.tif in results_dir"})
                print(f"  skip {e.tag}: no nucleus label image")
                continue
            pred = load_labels(outputs.labels_path)
            source = "pipeline labels"
        voxel, vsrc = our_voxel_zyx(outputs, e, truth)
        tl = load_labels(truth.labels_path)
        if pred.shape != tl.shape:
            skipped.append({"tag": e.tag, "reason": f"label images differ in shape: ours {pred.shape}, truth {tl.shape}"})
            print(f"  skip {e.tag}: shape mismatch {pred.shape} vs {tl.shape}")
            continue
        group = str(e.extra.get("condition") or e.split) if args.group_by == "condition" else e.split
        row, pt = n3.benchmark_gonad(pred, tl, voxel, tag=e.tag, split=group, roi_box=roi_box,
                                     coverage=args.coverage, dilate_um=args.dilate_um,
                                     majority_frac=args.majority_frac, max_dense_cells=args.max_dense_cells,
                                     min_truth_voxels=args.min_truth_voxels, coverage_xy_um=args.coverage_xy_um,
                                     coverage_z_planes=args.coverage_z_planes, scope_frac=args.scope_frac,
                                     split_frac=args.split_frac)
        row["voxel_source"] = vsrc
        row["roi_box"] = "" if roi_box is None else (
            " ".join(str(v) for v in roi_box["min_xyz"]) + " to " + " ".join(str(v) for v in roi_box["max_xyz"]))
        row["elapsed_s"] = round(time.time() - t0, 1)
        row["pred_source"] = source
        row["role"] = bd.gonad_role(e.tag, role) or ""
        for k, v in extra.items():
            row[k] = v
        rows.append(row)
        per_truth.append(pt)
        print(f"  {e.tag}: pred {row['n_pred']} (ignored {row['n_pred_ignored_outside_coverage']}, "
              f"{row['coverage_rule']}) truth {row['n_truth']} | IoU0.5 TP/FP/FN {row['tp_iou50']}/"
              f"{row['fp_iou50']}/{row['fn_iou50']} F1={row['f1_iou50']:.3f} AP={row['ap_iou50_95']:.3f} | "
              f"splits {row['n_splits']} merges {row['n_merges']} | z-extent err "
              f"{row['z_extent_error_mean_um']:+.2f} um | {row['elapsed_s']} s")
        del pred, tl
    settings = {"coverage": args.coverage, "dilate_um": args.dilate_um, "majority_frac": args.majority_frac,
                "max_dense_cells": args.max_dense_cells, "iou_thresholds": list(n3.IOU_THRESHOLDS),
                "group_by": args.group_by, "final": bool(args.final), "reason": args.reason,
                "pred_dir": args.pred_dir, "model": model_rec}
    if args.coverage == "protocol":
        settings["protocol"] = {"min_truth_voxels": args.min_truth_voxels, "coverage_xy_um": args.coverage_xy_um,
                                "coverage_z_planes": args.coverage_z_planes, "scope_frac": args.scope_frac,
                                "split_frac": args.split_frac}
    if args.model:
        settings["model_run"] = {"device": args.device, "tile": args.tile, "min_volume_um3": args.min_volume_um3,
                                 "fg_threshold": args.fg_threshold, "h_um": args.h_um,
                                 "neighbour_weight": args.neighbour_weight, "cache_root": args.cache_root}
    if not rows:
        print("no gonad could be scored")
        write_outputs(args.out, pd.DataFrame(), pd.DataFrame(),
                      {"n_gonads": 0, "skipped": skipped, "settings": settings})
        return 1
    per_gonad = pd.DataFrame(rows)
    per_nucleus = pd.concat(per_truth, ignore_index=True)
    pooled = n3.pool(per_gonad, per_nucleus, n_boot=args.n_boot, seed=args.seed)
    pooled["skipped"] = skipped
    pooled["manifest"] = str(args.manifest)
    pooled["settings"] = settings
    paths = write_outputs(args.out, per_gonad, per_nucleus, pooled)
    f1, ap = pooled["f1_iou50"], pooled["ap_iou50_95"]
    print(f"POOLED {len(rows)} gonads: F1@0.5={f1['point']:.3f} [{f1['lo']:.3f}, {f1['hi']:.3f}] "
          f"AP={ap['point']:.3f} [{ap['lo']:.3f}, {ap['hi']:.3f}] splits {pooled['n_splits']} "
          f"merges {pooled['n_merges']} ignored {pooled['n_pred_ignored_outside_coverage']}")
    print("outputs ->", paths["pooled"].parent)
    return 0


if __name__ == "__main__":
    sys.exit(main())
