#!/usr/bin/env python
"""Re-train the real-Koehler nucleus model on the LEAK-FREE group-by-shape split (data/train_kohler_grp)
so the held-out F1 is trustworthy (no same-gonad adjacent-crop leakage). Same recipe as the original
real model; only the split changed.  -> models/models/germline_nuclei_real_grp

    python scripts/train_real_grouped.py
"""
import glob

import numpy as np
import tifffile
from cellpose import models, train


def load(folder):
    imgs, masks = [], []
    for f in sorted(glob.glob(folder + "/*.tif")):
        if f.endswith("_masks.tif"):
            continue
        imgs.append(tifffile.imread(f).astype(np.float32))
        masks.append(tifffile.imread(f.replace(".tif", "_masks.tif")).astype(np.uint16))
    return imgs, masks


def main():
    trd, trl = load("data/train_kohler_grp")
    ted, tel = load("data/train_kohler_grp/test")
    print(f"LEAK-FREE train {len(trd)} imgs, test {len(ted)} imgs (held-out 0.108 um/px gonads)")
    model = models.CellposeModel(gpu=True)
    out = train.train_seg(
        model.net,
        train_data=trd, train_labels=trl,
        test_data=ted, test_labels=tel,
        n_epochs=150, learning_rate=1e-5, weight_decay=0.1, batch_size=8,
        save_path="models", model_name="germline_nuclei_real_grp",
    )
    print("MODEL SAVED:", out[0] if isinstance(out, (list, tuple)) else out)


if __name__ == "__main__":
    main()
