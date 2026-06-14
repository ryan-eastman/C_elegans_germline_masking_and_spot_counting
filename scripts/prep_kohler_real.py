#!/usr/bin/env python
"""Build a CLEAN Cellpose training set from the KoehlerLab real hand-annotated germline nuclei
(xy slices only). Excludes uint8-saturated masks (>=255 labels -> corrupted instance ids).
Held-out split BY IMAGE (no slice leakage). Writes Cellpose-style image.tif / image_masks.tif.

    python scripts/prep_kohler_real.py
"""
import glob
import os

import numpy as np
import tifffile
from skimage.io import imread
from skimage.segmentation import relabel_sequential

SRC = r"C:/Users/ryane/Cellpose_germlineNuclei/training/trainingData"
OUT = "data/train_kohler"
TEST_EVERY = 5   # every 5th clean image -> test (~20%)


def gray(a):
    a = np.squeeze(np.asarray(a))
    while a.ndim > 2:
        a = a[0] if a.shape[0] <= 4 else a[..., 0]
    return a


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(OUT + "/test", exist_ok=True)
    xy = sorted(glob.glob(SRC + "/tifs/*_xy.tif"))
    ntr = nte = skipped = 0
    for i, f in enumerate(xy):
        base = os.path.basename(f)[:-4]
        mp = f"{SRC}/masks/{base}_masks.png"
        if not os.path.exists(mp):
            continue
        raw = imread(mp)
        gt = gray(raw).astype(np.int64)
        # drop uint8-saturated masks: >=255 labels in a uint8 png == capped/corrupted ids
        if raw.dtype == np.uint8 and gt.max() >= 255:
            skipped += 1
            continue
        img = gray(tifffile.imread(f)).astype(np.float32)
        lbl = relabel_sequential(gt)[0].astype(np.uint16)
        if lbl.max() < 3:                      # near-empty field, skip
            skipped += 1
            continue
        dst = OUT + "/test" if (i % TEST_EVERY == 0) else OUT
        tifffile.imwrite(f"{dst}/{base}.tif", img)
        tifffile.imwrite(f"{dst}/{base}_masks.tif", lbl)
        if dst.endswith("test"):
            nte += 1
        else:
            ntr += 1
    print(f"wrote {ntr} train + {nte} test clean Koehler xy images  (skipped {skipped})")


if __name__ == "__main__":
    main()
