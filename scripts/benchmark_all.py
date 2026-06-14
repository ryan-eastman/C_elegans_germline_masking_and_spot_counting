#!/usr/bin/env python
"""Unified benchmark: evaluate every nucleus model on BOTH held-out real (Koehler) and synthetic
test sets. One table -> which training recipe actually generalizes to real germline data.

    python scripts/benchmark_all.py
"""
import glob
import os

import numpy as np
import tifffile
from skimage.io import imread

from cellpose import models
from germquant.validate import segmentation_metrics

KOHLER = "data/train_kohler/test"          # held-out REAL
SYNTH = "data/train_nuclei/test"           # held-out SYNTHETIC

MODELS = [
    ("STOCK cpsam", None),
    ("synthetic-only", "models/models/germline_nuclei"),
    ("real-Koehler", "models/models/germline_nuclei_real"),
]


def gray(a):
    a = np.squeeze(np.asarray(a))
    while a.ndim > 2:
        a = a[0] if a.shape[0] <= 4 else a[..., 0]
    return a


def eval_set(model, folder):
    imgs = sorted(f for f in glob.glob(folder + "/*.tif") if not f.endswith("_masks.tif"))
    f1s, ious, dets, gts = [], [], [], []
    for f in imgs:
        gtf = f[:-4] + "_masks.tif"
        if not os.path.exists(gtf):
            continue
        img = gray(tifffile.imread(f)).astype(np.float32)
        gt = gray(tifffile.imread(gtf)).astype(np.int32)
        pred = np.asarray(model.eval(img, channel_axis=None)[0])
        m = segmentation_metrics(pred, gt, iou_threshold=0.5)
        f1s.append(m["f1"]); ious.append(m["mean_iou"]); dets.append(pred.max()); gts.append(gt.max())
    return dict(f1=np.mean(f1s), iou=np.nanmean(ious), pred=np.mean(dets), gt=np.mean(gts), n=len(f1s))


def main():
    rows = []
    for name, path in MODELS:
        if path and not os.path.exists(path):
            print(f"skip {name}: {path} not found")
            continue
        model = models.CellposeModel(gpu=True, pretrained_model=path) if path \
            else models.CellposeModel(gpu=True)
        r = eval_set(model, KOHLER)
        s = eval_set(model, SYNTH)
        rows.append((name, r, s))
        print(f"[{name}] done")
    print("\n" + "=" * 78)
    print(f"{'model':<16} | {'REAL (Koehler held-out)':^28} | {'SYNTHETIC held-out':^26}")
    print(f"{'':<16} | {'F1@.5':>6} {'IoU':>6} {'pred':>5} {'gt':>5} | {'F1@.5':>6} {'IoU':>6} {'pred':>5} {'gt':>5}")
    print("-" * 78)
    for name, r, s in rows:
        print(f"{name:<16} | {r['f1']:6.3f} {r['iou']:6.3f} {r['pred']:5.0f} {r['gt']:5.0f} | "
              f"{s['f1']:6.3f} {s['iou']:6.3f} {s['pred']:5.0f} {s['gt']:5.0f}")
    print("=" * 78)
    print(f"(real n={rows[0][1]['n']} imgs, synth n={rows[0][2]['n']} imgs)")


if __name__ == "__main__":
    main()
