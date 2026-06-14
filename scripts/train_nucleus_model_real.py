#!/usr/bin/env python
"""Fine-tune Cellpose-SAM on REAL KoehlerLab germline-nucleus annotations (the sim-to-real test).
The synthetic-only model overfit (synth-test F1 0.98 but Koehler-real F1 0.27 < stock 0.37).
This trains on real hand labels to see whether real-data fine-tuning beats stock on held-out real.

    python scripts/train_nucleus_model_real.py
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
    trd, trl = load("data/train_kohler")
    ted, tel = load("data/train_kohler/test")
    print(f"train {len(trd)} real imgs, test {len(ted)} real imgs")
    model = models.CellposeModel(gpu=True)  # cpsam base
    out = train.train_seg(
        model.net,
        train_data=trd, train_labels=trl,
        test_data=ted, test_labels=tel,
        n_epochs=150, learning_rate=1e-5, weight_decay=0.1, batch_size=8,
        save_path="models", model_name="germline_nuclei_real",
    )
    print("MODEL SAVED:", out[0] if isinstance(out, (list, tuple)) else out)


if __name__ == "__main__":
    main()
