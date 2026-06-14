#!/usr/bin/env python
"""Validate a nucleus-segmentation model against REAL hand-annotated germline masks from the
KoehlerLab Cellpose_germlineNuclei repo -- independent real ground truth. Instance F1 / IoU
(via germquant's Metrics-Reloaded-style segmentation_metrics).

    python scripts/validate_vs_real_masks.py            # stock cpsam
    python scripts/validate_vs_real_masks.py models/models/germline_nuclei   # our model
"""
import glob
import os
import sys

import numpy as np
import tifffile
from skimage.io import imread

from cellpose import models
from germquant.validate import segmentation_metrics

DATA = r"C:/Users/ryane/Cellpose_germlineNuclei/training/trainingData"


def gray(a):
    a = np.squeeze(np.asarray(a))
    while a.ndim > 2:                      # take the largest 2D plane / drop channels
        a = a[0] if a.shape[0] <= 4 else a[..., 0]
    return a


def main():
    model_path = sys.argv[1] if len(sys.argv) > 1 else None
    model = models.CellposeModel(gpu=True, pretrained_model=model_path) if model_path \
        else models.CellposeModel(gpu=True)
    tifs = sorted(f for f in glob.glob(DATA + "/tifs/*_xy.tif"))   # xy slices (match our training)
    print(f"model: {model_path or 'STOCK cpsam'}  |  {len(tifs)} real xy images")
    f1s, ious, dets, gts = [], [], [], []
    for f in tifs:
        base = os.path.basename(f)[:-4]
        mp = f"{DATA}/masks/{base}_masks.png"
        if not os.path.exists(mp):
            continue
        img = gray(tifffile.imread(f)).astype(np.float32)
        gt = gray(imread(mp)).astype(np.int32)
        pred = np.asarray(model.eval(img, channel_axis=None)[0])
        m = segmentation_metrics(pred, gt, iou_threshold=0.5)
        f1s.append(m["f1"]); ious.append(m["mean_iou"]); dets.append(pred.max()); gts.append(gt.max())
    print(f"  mean F1@0.5 = {np.mean(f1s):.3f}   mean matched-IoU = {np.nanmean(ious):.3f}")
    print(f"  mean nuclei: predicted {np.mean(dets):.0f}  vs  ground-truth {np.mean(gts):.0f}")
    print(f"  (n={len(f1s)} images; F1 range {min(f1s):.2f}-{max(f1s):.2f})")


if __name__ == "__main__":
    main()
