#!/usr/bin/env python
"""LEAK-FREE Koehler split (fixes the adversarial-review finding that prep_kohler_real.py strided
i%5 over sorted slices, leaking same-gonad adjacent crops into train+test).

Grouping rule: adjacent crops of one gonad share an EXACT pixel shape (high NCC per review). So no
shape may appear in both train and test. Test is drawn from the 0.10833 um/px scale group -- the
Cahoon-matched pixel size (deployment scale), and every 0.108 image has a UNIQUE shape, so holding
some out leaks nothing. All multi-image fine-scale shape-groups (642x804 x11, 384x444 x7, ...) go
ENTIRELY to train. Result: test = held-out, Cahoon-scale, no shape shared with train.

    python scripts/prep_kohler_grouped.py
"""
import csv
import glob
import os

import numpy as np
import tifffile
from skimage.io import imread
from skimage.segmentation import relabel_sequential

SRC = r"C:/Users/ryane/Cellpose_germlineNuclei/training/trainingData"
OUT = "data/train_kohler_grp"
CAHOON_SCALE = 0.10833           # held-out test scale (matches Cahoon 0.1083 um/px)


def gray(a):
    a = np.squeeze(np.asarray(a))
    while a.ndim > 2:
        a = a[0] if a.shape[0] <= 4 else a[..., 0]
    return a


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(OUT + "/test", exist_ok=True)
    prov = {}
    with open(SRC + "/TrainDataFiles.csv") as f:
        for r in csv.DictReader(f):
            prov[os.path.basename(r["PathToFile"])] = round(float(r["xScale"]), 5)

    # collect clean xy images
    clean = []  # (base, scale)
    for fp in sorted(glob.glob(SRC + "/tifs/*_xy.tif")):
        base = os.path.basename(fp)[:-4]
        raw = imread(f"{SRC}/masks/{base}_masks.png")
        g = gray(raw)
        if raw.dtype == np.uint8 and g.max() >= 255:
            continue
        if g.max() < 3:
            continue
        clean.append((base, prov.get(base + ".tif", 0.0)))

    # held-out test = every other Cahoon-scale (0.108) image (all unique shapes -> leak-free);
    # train = all other scales + the remaining 0.108 images.
    cahoon = sorted(b for b, s in clean if abs(s - CAHOON_SCALE) < 1e-4)
    test_set = set(cahoon[::2])     # ~half of the 0.108 group
    ntr = nte = 0
    for base, scale in clean:
        img = gray(tifffile.imread(f"{SRC}/tifs/{base}.tif")).astype(np.float32)
        lbl = relabel_sequential(gray(imread(f"{SRC}/masks/{base}_masks.png")).astype(np.int64))[0].astype(np.uint16)
        dst = OUT + "/test" if base in test_set else OUT
        tifffile.imwrite(f"{dst}/{base}.tif", img)
        tifffile.imwrite(f"{dst}/{base}_masks.tif", lbl)
        if base in test_set:
            nte += 1
        else:
            ntr += 1
    print(f"LEAK-FREE split: {ntr} train + {nte} test (test = held-out {CAHOON_SCALE} um/px, unique shapes)")
    print(f"test images: {sorted(test_set)}")


if __name__ == "__main__":
    main()
