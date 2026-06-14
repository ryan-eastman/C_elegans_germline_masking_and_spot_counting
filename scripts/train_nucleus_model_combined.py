#!/usr/bin/env python
"""Ablation: train on REAL Koehler + a BALANCED subset of synthetic (~1:1). Tests whether synthetic
data ADDS Cahoon-matched pachytene morphology on top of real texture, or DRAGS the model toward the
synthetic-only failure mode (real held-out F1 0.36). Test set = real held-out (the Cahoon-relevant
generalization measure).

    python scripts/train_nucleus_model_combined.py
"""
import glob

import numpy as np
import tifffile
from cellpose import models, train


def load(folder, every=1):
    imgs, masks = [], []
    fs = [f for f in sorted(glob.glob(folder + "/*.tif")) if not f.endswith("_masks.tif")]
    for f in fs[::every]:
        imgs.append(tifffile.imread(f).astype(np.float32))
        masks.append(tifffile.imread(f.replace(".tif", "_masks.tif")).astype(np.uint16))
    return imgs, masks


def main():
    rd, rl = load("data/train_kohler")              # 48 real
    sd, sl = load("data/train_nuclei", every=4)     # ~49 synthetic (every 4th of 198) -> ~1:1
    trd, trl = rd + sd, rl + sl
    ted, tel = load("data/train_kohler/test")       # real held-out
    print(f"train {len(rd)} real + {len(sd)} synth = {len(trd)} imgs ; test {len(ted)} real imgs")
    model = models.CellposeModel(gpu=True)
    out = train.train_seg(
        model.net,
        train_data=trd, train_labels=trl,
        test_data=ted, test_labels=tel,
        n_epochs=150, learning_rate=1e-5, weight_decay=0.1, batch_size=8,
        save_path="models", model_name="germline_nuclei_combined",
    )
    print("MODEL SAVED:", out[0] if isinstance(out, (list, tuple)) else out)


if __name__ == "__main__":
    main()
