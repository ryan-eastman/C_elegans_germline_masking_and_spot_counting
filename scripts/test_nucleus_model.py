#!/usr/bin/env python
"""Moment of truth: run STOCK cpsam vs the fine-tuned germline model on a REAL pachytene crop
and overlay nucleus boundaries. The question: does the fine-tuned model enclose WHOLE nuclei
(all chromosomes) where stock cpsam grabs part and misses chromosomes.

    python scripts/test_nucleus_model.py
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from cellpose import models
from skimage.segmentation import find_boundaries

import germquant.io.nd2_reader as R

CY, CX, CZ, H = 1761, 1423, 47, 185
REAL = r"data/raw_examples/madeleline images/20251105_N2_noHS/20251105_N2_nohs_HERM _001.nd2"


def norm(a):
    lo, hi = np.percentile(a, 1), np.percentile(a, 99.5)
    return np.clip((a - lo) / (hi - lo + 1e-9), 0, 1)


def seg(model, img):
    out = model.eval(img, channel_axis=None)
    return np.asarray(out[0])


def main():
    st = R.read_stack(REAL)
    crop = st.data[0][CZ, CY - H:CY + H, CX - H:CX + H].astype(np.float32)
    stock = models.CellposeModel(gpu=True)
    custom = models.CellposeModel(gpu=True, pretrained_model="models/germline_nuclei")

    ms = seg(stock, crop)
    mc = seg(custom, crop)
    dn = norm(crop)
    print(f"stock: {ms.max()} nuclei | custom: {mc.max()} nuclei")

    fig, ax = plt.subplots(1, 3, figsize=(21, 7))
    ax[0].imshow(dn, cmap="gray"); ax[0].set_title("REAL DAPI (pachytene)"); ax[0].axis("off")
    for a, m, t in ((ax[1], ms, f"STOCK cpsam  (n={ms.max()})"), (ax[2], mc, f"FINE-TUNED germline  (n={mc.max()})")):
        a.imshow(dn, cmap="gray")
        b = find_boundaries(m, mode="outer"); ov = np.zeros((*b.shape, 4)); ov[b] = (1, 1, 0, 1)
        a.imshow(ov); a.set_title(t); a.axis("off")
    fig.tight_layout(); fig.savefig("results_fullres/model_compare.png", dpi=130, bbox_inches="tight")
    print("wrote results_fullres/model_compare.png")


if __name__ == "__main__":
    main()
