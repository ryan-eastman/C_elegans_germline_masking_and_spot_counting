#!/usr/bin/env python
"""Sanity check: evaluate STOCK cpsam and our fine-tuned model on the HELD-OUT SYNTHETIC test set
(same distribution the model trained on). Isolates 'did training work at all' from 'does it
generalize to a different microscope (Koehler)'.

    python scripts/validate_vs_synth_test.py            # stock
    python scripts/validate_vs_synth_test.py models/models/germline_nuclei   # our model
"""
import glob
import os
import sys

import numpy as np
import tifffile

from cellpose import models
from germquant.validate import segmentation_metrics

DATA = r"data/train_nuclei/test"


def main():
    model_path = sys.argv[1] if len(sys.argv) > 1 else None
    model = models.CellposeModel(gpu=True, pretrained_model=model_path) if model_path \
        else models.CellposeModel(gpu=True)
    imgs = sorted(f for f in glob.glob(DATA + "/*.tif") if not f.endswith("_masks.tif"))
    print(f"model: {model_path or 'STOCK cpsam'}  |  {len(imgs)} synthetic test images")
    f1s, ious, dets, gts = [], [], [], []
    for f in imgs:
        gtf = f[:-4] + "_masks.tif"
        if not os.path.exists(gtf):
            continue
        img = tifffile.imread(f).astype(np.float32)
        gt = tifffile.imread(gtf).astype(np.int32)
        pred = np.asarray(model.eval(img, channel_axis=None)[0])
        m = segmentation_metrics(pred, gt, iou_threshold=0.5)
        f1s.append(m["f1"]); ious.append(m["mean_iou"]); dets.append(pred.max()); gts.append(gt.max())
    print(f"  mean F1@0.5 = {np.mean(f1s):.3f}   mean matched-IoU = {np.nanmean(ious):.3f}")
    print(f"  mean nuclei: predicted {np.mean(dets):.0f}  vs  ground-truth {np.mean(gts):.0f}")
    print(f"  (n={len(f1s)} images; F1 range {min(f1s):.2f}-{max(f1s):.2f})")


if __name__ == "__main__":
    main()
