#!/usr/bin/env python
"""Generate a varied synthetic training set for the germline-nucleus Cellpose model.
Produces 2D (image, instance-mask) slice pairs from many synthetic 3D volumes with different
density / size / chromosome-count, in the Cellpose convention (<name>.tif + <name>_masks.tif).

    python scripts/make_training_set.py
"""
import os
import sys

import numpy as np
import tifffile
from skimage.segmentation import relabel_sequential

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from synth_nuclei import P, degrade, synthesize  # noqa: E402

OUT = "data/train_nuclei"


def configs():
    cfgs = []
    for sp in (4.6, 5.4, 6.2):          # packing density: crowded -> sparse
        for nr in (2.0, 2.4, 2.8):      # nucleus size
            for ns in (6, 5):           # oocyte (6) / spermatocyte (5) chromosomes
                cfgs.append(dict(pack_spacing_um=sp, nuc_radius_um=nr, n_strings=ns))
    return cfgs


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(OUT + "/test", exist_ok=True)
    cfgs = configs()
    n_pairs = 0
    for ci, ov in enumerate(cfgs):
        p = {**P, **ov}
        clean, label = synthesize(p=p, seed=ci)
        img = degrade(clean, p=p, seed=1000 + ci)
        im16 = (np.clip(img / (img.max() + 1e-9), 0, 1) * 65535).astype(np.uint16)
        dest = OUT + "/test" if ci % 6 == 5 else OUT          # ~1 in 6 volumes held out
        nz = img.shape[0]
        for z in range(int(0.20 * nz), int(0.80 * nz), 4):    # skip edge slices (nuclei cut)
            lab = label[z]
            if np.unique(lab).size < 7:                       # need several nuclei in the slice
                continue
            lab2, _, _ = relabel_sequential(lab)
            base = f"{dest}/v{ci:02d}_z{z:03d}"
            tifffile.imwrite(base + ".tif", im16[z])
            tifffile.imwrite(base + "_masks.tif", lab2.astype(np.uint16))
            n_pairs += 1
        print(f"  vol {ci:02d} {ov} -> {'TEST' if dest.endswith('test') else 'train'}")
    ntest = len([f for f in os.listdir(OUT + "/test") if f.endswith("_masks.tif")])
    print(f"wrote {n_pairs} slice pairs total ({ntest} held out for test) under {OUT}/")


if __name__ == "__main__":
    main()
