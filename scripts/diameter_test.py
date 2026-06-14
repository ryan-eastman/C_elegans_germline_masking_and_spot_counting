#!/usr/bin/env python
"""Empirical: which inference regime is MORE CORRECT on real ground truth -- native scale (no diameter,
how we benchmarked) or the pipeline's diameter rescale (~28px)? Score both on the leak-free 0.108 um/px
Koehler held-out GT and let the data decide.

    python scripts/diameter_test.py
"""
import glob
import os

import numpy as np
import tifffile

from cellpose import models
from germquant.validate import segmentation_metrics

GRP = "data/train_kohler_grp/test"
DIAM_PX = 3.0 / 0.108333          # what the pipeline passes (diameter_um / dy)
MODELS = [
    ("real (leak-free)", "models/models/germline_nuclei_real_grp"),
    ("real+synth combined", "models/models/germline_nuclei_combined"),
]


def pairs():
    out = []
    for f in sorted(glob.glob(GRP + "/*.tif")):
        if f.endswith("_masks.tif"):
            continue
        g = f[:-4] + "_masks.tif"
        if os.path.exists(g):
            out.append((f, g))
    return out


def score(model, ps, diameter):
    f1s = []
    for imgf, gtf in ps:
        img = np.squeeze(tifffile.imread(imgf)).astype(np.float32)
        gt = np.squeeze(tifffile.imread(gtf)).astype(np.int32)
        kw = {} if diameter is None else {"diameter": diameter}
        pred = np.asarray(model.eval(img, channel_axis=None, **kw)[0])
        f1s.append(segmentation_metrics(pred, gt, iou_threshold=0.5)["f1"])
    return np.mean(f1s), min(f1s), max(f1s)


def main():
    ps = pairs()
    print(f"{len(ps)} leak-free 0.108 um/px GT images;  pipeline diameter = {DIAM_PX:.1f} px\n")
    print(f"{'model':<20} {'NATIVE (no diam)':>22} {'DIAMETER rescale':>22}")
    for name, path in MODELS:
        m = models.CellposeModel(gpu=True, pretrained_model=path)
        nf, nlo, nhi = score(m, ps, None)
        df, dlo, dhi = score(m, ps, DIAM_PX)
        print(f"{name:<20} {nf:.3f} ({nlo:.2f}-{nhi:.2f}){'':>4} {df:.3f} ({dlo:.2f}-{dhi:.2f})")
    print("\n-> use whichever scores higher on this real GT as the pipeline default.")


if __name__ == "__main__":
    main()
