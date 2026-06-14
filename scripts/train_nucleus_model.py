#!/usr/bin/env python
"""Fine-tune Cellpose-SAM into a germline-nucleus model on the synthetic training set.
The model learns that '6 sharp chromosome threads in a ball = ONE nucleus' (the thing stock
cpsam fails at) from synthetic data with perfect whole-nucleus labels.

    python scripts/train_nucleus_model.py
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
    trd, trl = load("data/train_nuclei")
    ted, tel = load("data/train_nuclei/test")
    print(f"train {len(trd)} imgs, test {len(ted)} imgs")
    model = models.CellposeModel(gpu=True)  # cpsam
    out = train.train_seg(
        model.net,
        train_data=trd, train_labels=trl,
        test_data=ted, test_labels=tel,
        n_epochs=150, learning_rate=1e-5, weight_decay=0.1, batch_size=8,
        save_path="models", model_name="germline_nuclei",
    )
    print("MODEL SAVED:", out[0] if isinstance(out, (list, tuple)) else out)


if __name__ == "__main__":
    main()
