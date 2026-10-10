"""Export 2D training slices of the TRAIN gonads' DAPI channel with their recovered Imaris nucleus masks, in
the pairs format of ``germquant train-nuclei`` (image/, labels/, manifest.csv), for the fine-tuned Cellpose
competitor of the nuclei3d benchmark (docs/BENCHMARK.md, "nuclei3d").

    python scripts/export_imaris_training_pairs.py --out PAIRS [--split-file tests/golden/rad51_split.json]
        [--manifest M] [--cache-root C] [--prep-dir P] [--slices-per-image 8] [--orientations xy,xz,yz]
        [--seed 0] [--min-labels 3] [--crop-margin 48] [--no-mask-untraced]
    germquant train-nuclei --pairs PAIRS --out MODELDIR --base <a Cellpose model file> --split-by series
    python scripts/benchmark_nuclei3d.py --manifest M --split val --coverage protocol \\
        --cellpose-model MODELDIR/<name> --out DIR

Which gonads. The train gonads of the split file, all six (the nucleus surfaces of the two heat-shock males
are intact); the val gonads are for model selection and the test gonads are never opened (every gonad
passes the harness's test guard before anything is read). The inputs are those nuclei3d trains on: the
bench cache's DAPI channel and the truth files' ``<stem>__imaris_labels.tif``.

Partial annotation. Only the traced nuclei carry labels, and Cellpose's training loss has no ignore label,
so an untraced nucleus left in a slice would be taught as background. The voxels nuclei3d leaves out of its
loss (inside an untraced deployed Cellpose label or within 0.5 um of one, and outside every traced mask:
the ignore class of germquant.segment.nuclei3d.data, read from or written to ``--prep-dir``) are filled with
the gonad's 1st DAPI percentile, the background level, before slicing; ``--no-mask-untraced`` keeps the raw
image. The fill changes what the network sees (an untraced nucleus becomes a dark hole with a sharp edge),
a known limit of this competitor.

Slices. Per gonad and orientation, ``--slices-per-image`` slices are drawn at random (seeded per gonad)
among those that cut at least ``--min-labels`` traced nuclei (germquant.segment.training.choose_slices);
xz and yz slices are resampled along z to the xy pixel size (training.isotropic_slice, linear for the
image and nearest neighbour for the labels), and every slice is cropped to the bounding box of its traced
labels plus ``--crop-margin`` pixels. All pairs carry the series id ``train``, so ``train-nuclei
--split-by series`` trains on every pair (a single group leaves no hold-out); the evaluation is the 3D
scoring of the val gonads.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="the pairs folder (image/, labels/, manifest.csv, export.json)")
    ap.add_argument("--split-file", default=str(REPO / "tests" / "golden" / "rad51_split.json"),
                    help="the pre-registered split record (train / val / test roles)")
    ap.add_argument("--manifest", default=None, help="benchmark manifest (default: the harness's)")
    ap.add_argument("--cache-root", default=None, help="bench cache root (default: the harness's)")
    ap.add_argument("--prep-dir", default=None,
                    help="folder of the nuclei3d target maps (read, or written when absent; default <out>/prep)")
    ap.add_argument("--slices-per-image", type=int, default=8, help="slices per gonad and orientation (default 8)")
    ap.add_argument("--orientations", default="xy,xz,yz", help="comma-separated subset of xy,xz,yz")
    ap.add_argument("--seed", type=int, default=0, help="seed of the slice draw (default 0)")
    ap.add_argument("--min-labels", type=int, default=3, help="a slice must cut this many traced nuclei (default 3)")
    ap.add_argument("--crop-margin", type=int, default=48, help="pixels kept around the traced labels (default 48)")
    ap.add_argument("--no-mask-untraced", action="store_true", help="keep the raw image (no fill of untraced nuclei)")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    import tifffile

    from germquant.benchmark import detector as bd
    from germquant.benchmark.manifest import load_labels, load_manifest, load_truth
    from germquant.segment.nuclei3d import data as d3
    from germquant.segment.nuclei3d.infer import robust_range
    from germquant.segment.training import (MANIFEST_COLS, ORIENTATIONS, _rng_for, _to_uint16, choose_slices,
                                            isotropic_slice, slice_counts)

    args = parse_args(argv)
    orientations = tuple(o.strip().lower() for o in str(args.orientations).split(",") if o.strip())
    if not orientations or any(o not in ORIENTATIONS for o in orientations):
        print(f"--orientations must be among {ORIENTATIONS}", file=sys.stderr)
        return 2
    split = json.loads(Path(args.split_file).read_text(encoding="utf-8"))
    roles = {str(t): bd.gonad_role(str(t), role) for role, tags in (split.get("roles") or {}).items() for t in tags}
    entries = {e.tag: e for e in load_manifest(args.manifest or bd.DEFAULT_MANIFEST)}
    train = [t for t, r in roles.items() if r == "train"]
    out = Path(args.out)
    prep = Path(args.prep_dir) if args.prep_dir else out / "prep"
    for sub in ("image", "labels"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    rows, done = [], []
    for tag in train:
        e = entries.get(tag)
        if e is None:
            print(f"  skip {tag}: not in the manifest")
            continue
        man_role = str(e.extra.get("role") or e.split or "").strip()
        if man_role in ("val", "test"):
            raise ValueError(f"{tag}: the split file says train, the manifest says {man_role!r}")
        bd.guard_test([tag], allow_test=False, roles={tag: bd.entry_role(e)})
        cache = bd.load_cache(tag, args.cache_root, arrays=("dna", "labels"), mmap=True)
        truth = load_truth(e.truth_dir, e.truth_stem)
        if truth is None or not truth.has_labels or cache.dna is None:
            print(f"  skip {tag}: no truth label image or no DNA channel")
            continue
        tl = load_labels(truth.labels_path)
        img = np.array(cache.dna)
        n_filled = 0
        if not args.no_mask_untraced:
            tg = d3.load_targets(prep, tag, cache.meta, cache.labels, truth.labels_path, cache.spacing, {})
            fill = np.asarray(np.rint(robust_range(img)[0])).astype(img.dtype)
            for z in range(img.shape[0]):
                m = np.asarray(tg["cls"][z]) == d3.CLASS_IGNORE
                img[z][m] = fill
                n_filled += int(m.sum())
            del tg
        rng = _rng_for(args.seed, tag)
        dna_idx = (cache.meta.get("role_indices") or {}).get("dna")
        n_before = len(rows)
        for orient in orientations:
            counts = slice_counts(tl, orient)
            for i in choose_slices(counts, args.slices_per_image, args.min_labels, rng):
                im2, lab2, px, z_zoom = isotropic_slice(img, tl, orient, i, cache.spacing)
                ys, xs = np.nonzero(lab2 > 0)
                if ys.size == 0:
                    continue
                m = int(args.crop_margin)
                y0, y1 = max(0, ys.min() - m), min(lab2.shape[0], ys.max() + m + 1)
                x0, x1 = max(0, xs.min() - m), min(lab2.shape[1], xs.max() + m + 1)
                im2, lab2 = im2[y0:y1, x0:x1], lab2[y0:y1, x0:x1]
                stem = f"{tag}_{orient}_{i:04d}.tif"
                tifffile.imwrite(out / "image" / stem, _to_uint16(im2), compression="zlib")
                tifffile.imwrite(out / "labels" / stem, lab2.astype(np.int32), compression="zlib")
                ids = np.unique(lab2)
                rows.append({"image_id": tag, "series_id": "train", "t_index": "", "orientation": orient,
                             "slice_index": int(i), "n_labels": int(np.count_nonzero(ids)), "zone": "",
                             "image": f"image/{stem}", "labels": f"labels/{stem}",
                             "source_file": str(e.nd2 or ""), "xy_stride": 1, "z_range": "", "channel_role": "dna",
                             "channel_index": dna_idx, "pixel_um": px, "z_zoom": z_zoom,
                             "height": int(lab2.shape[0]), "width": int(lab2.shape[1])})
        done.append({"tag": tag, "n_pairs": len(rows) - n_before, "n_voxels_filled": n_filled,
                     "truth_labels": str(truth.labels_path)})
        print(f"  {tag}: {len(rows) - n_before} slice pairs ({n_filled} voxels filled as background)")
        del img, tl
    pd.DataFrame(rows, columns=MANIFEST_COLS).to_csv(out / "manifest.csv", index=False)
    (out / "export.json").write_text(json.dumps({
        "what": "train gonads' DAPI slices with the recovered Imaris masks (scripts/export_imaris_training_pairs.py)",
        "split_file": str(args.split_file), "gonads": done, "orientations": list(orientations),
        "slices_per_image": args.slices_per_image, "seed": args.seed, "min_labels": args.min_labels,
        "crop_margin": args.crop_margin, "mask_untraced": not args.no_mask_untraced,
        "test_never_opened": [t for t, r in roles.items() if r == "test"],
        "date": datetime.now(timezone.utc).isoformat(timespec="seconds")}, indent=1), encoding="utf-8")
    print(f"{len(rows)} slice pairs from {len(done)} gonads -> {out}")
    return 0 if rows else 1


if __name__ == "__main__":
    sys.exit(main())
