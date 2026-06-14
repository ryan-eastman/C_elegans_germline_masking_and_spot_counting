#!/usr/bin/env python
"""Hardened benchmark (post adversarial-review). Evaluates every model on:
  (1) LEAK-FREE 0.108 um/px held-out Koehler xy test (the trustworthy, Cahoon-scale number)
  (2) Koehler xz/yz orientation slices  -- probes the 'trained-on-xy, deployed-in-3D' gap
  (3) synthetic held-out
Reports mean F1 AND per-image spread (min-max) + predicted/GT nuclei. 2 sig figs.

    python scripts/benchmark2.py
"""
import glob
import os

import numpy as np
import tifffile
from skimage.io import imread

from cellpose import models
from germquant.validate import segmentation_metrics

SRC = r"C:/Users/ryane/Cellpose_germlineNuclei/training/trainingData"
GRP = "data/train_kohler_grp/test"      # leak-free 0.108 xy
SYNTH = "data/train_nuclei/test"

MODELS = [
    ("STOCK cpsam", None),
    ("synthetic-only", "models/models/germline_nuclei"),
    ("real (leaky split)", "models/models/germline_nuclei_real"),
    ("real (leak-free)", "models/models/germline_nuclei_real_grp"),
    ("real+synth combined", "models/models/germline_nuclei_combined"),
]


def gray(a):
    a = np.squeeze(np.asarray(a))
    while a.ndim > 2:
        a = a[0] if a.shape[0] <= 4 else a[..., 0]
    return a


def eval_pairs(model, pairs):
    f1s, dets, gts = [], [], []
    for imgf, gtf, is_png in pairs:
        img = gray(tifffile.imread(imgf)).astype(np.float32)
        gt = gray(imread(gtf) if is_png else tifffile.imread(gtf)).astype(np.int32)
        pred = np.asarray(model.eval(img, channel_axis=None)[0])
        m = segmentation_metrics(pred, gt, iou_threshold=0.5)
        f1s.append(m["f1"]); dets.append(pred.max()); gts.append(gt.max())
    return dict(f1=np.mean(f1s), lo=min(f1s), hi=max(f1s), pred=np.mean(dets), gt=np.mean(gts), n=len(f1s))


def tif_pairs(folder):
    out = []
    for f in sorted(glob.glob(folder + "/*.tif")):
        if f.endswith("_masks.tif"):
            continue
        gtf = f[:-4] + "_masks.tif"
        if os.path.exists(gtf):
            out.append((f, gtf, False))
    return out


def ortho_pairs():
    """Clean (non-saturated) xz + yz Koehler slices."""
    out = []
    for f in sorted(glob.glob(SRC + "/tifs/*_xz.tif") + glob.glob(SRC + "/tifs/*_yz.tif")):
        base = os.path.basename(f)[:-4]
        mp = f"{SRC}/masks/{base}_masks.png"
        if not os.path.exists(mp):
            continue
        raw = imread(mp)
        if raw.dtype == np.uint8 and gray(raw).max() >= 255:
            continue
        out.append((f, mp, True))
    return out


def main():
    grp = tif_pairs(GRP)
    ortho = ortho_pairs()
    synth = tif_pairs(SYNTH)
    print(f"sets: leak-free-xy n={len(grp)} | xz/yz n={len(ortho)} | synth n={len(synth)}\n")
    print(f"{'model':<20} | {'leak-free 0.108 xy':^26} | {'xz/yz ortho':^14} | {'synthetic':^14}")
    print(f"{'':<20} | {'F1 (min-max)':>16} {'pr/gt':>8} | {'F1':>6} {'pr/gt':>6} | {'F1':>6} {'pr/gt':>6}")
    print("-" * 92)
    for name, path in MODELS:
        if path and not os.path.exists(path):
            print(f"{name:<20} | (model not on disk -- skipped)")
            continue
        model = models.CellposeModel(gpu=True, pretrained_model=path) if path else models.CellposeModel(gpu=True)
        g = eval_pairs(model, grp); o = eval_pairs(model, ortho); s = eval_pairs(model, synth)
        print(f"{name:<20} | {g['f1']:.2f} ({g['lo']:.2f}-{g['hi']:.2f})  {g['pred']:.0f}/{g['gt']:.0f}"
              f" | {o['f1']:.2f}  {o['pred']:.0f}/{o['gt']:.0f} | {s['f1']:.2f}  {s['pred']:.0f}/{s['gt']:.0f}")


if __name__ == "__main__":
    main()
