#!/usr/bin/env python
"""3-way qualitative comparison on the REAL Cahoon pachytene crop (our actual target data, no GT
masks yet -- Imaris arrives Monday). Stock vs synthetic-only vs real-Koehler nucleus model.
The two failure modes to watch: (a) hallucinating nuclei in the dark rachis, (b) under-segmenting
(merging real nuclei into few blobs).

    python scripts/compare_models_cahoon.py
"""
import os

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from cellpose import models
from skimage.segmentation import find_boundaries

import germquant.io.nd2_reader as R

CY, CX, CZ, H = 1761, 1423, 47, 185
REAL = r"data/raw_examples/madeleline images/20251105_N2_noHS/20251105_N2_nohs_HERM _001.nd2"

MODELS = [
    ("STOCK cpsam", None),
    ("synthetic-only", "models/models/germline_nuclei"),
    ("real-Koehler", "models/models/germline_nuclei_real"),
]


def norm(a):
    lo, hi = np.percentile(a, 1), np.percentile(a, 99.5)
    return np.clip((a - lo) / (hi - lo + 1e-9), 0, 1)


def main():
    st = R.read_stack(REAL)
    crop = st.data[0][CZ, CY - H:CY + H, CX - H:CX + H].astype(np.float32)
    dn = norm(crop)
    panels = [("REAL DAPI (pachytene)", None)]
    for name, path in MODELS:
        if path and not os.path.exists(path):
            print(f"skip {name}: not found")
            continue
        model = models.CellposeModel(gpu=True, pretrained_model=path) if path \
            else models.CellposeModel(gpu=True)
        m = np.asarray(model.eval(crop, channel_axis=None)[0])
        print(f"{name}: {m.max()} nuclei")
        panels.append((f"{name}  (n={m.max()})", m))

    fig, ax = plt.subplots(1, len(panels), figsize=(7 * len(panels), 7))
    for a, (t, m) in zip(ax, panels):
        a.imshow(dn, cmap="gray"); a.set_title(t); a.axis("off")
        if m is not None:
            b = find_boundaries(m, mode="outer"); ov = np.zeros((*b.shape, 4)); ov[b] = (1, 1, 0, 1)
            a.imshow(ov)
    fig.tight_layout(); fig.savefig("results_fullres/model_compare3.png", dpi=130, bbox_inches="tight")
    print("wrote results_fullres/model_compare3.png")


if __name__ == "__main__":
    main()
