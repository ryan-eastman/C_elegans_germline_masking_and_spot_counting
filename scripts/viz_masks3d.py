#!/usr/bin/env python
"""Visual check of the 3D nucleus masks on the Cahoon pachytene crop: thin-MIP DAPI with each
model's label boundaries overlaid. Rules out the 'few giant blobs' artifact -- confirms whether
real-Koehler masks follow individual whole nuclei (vs stock/synthetic fragmenting).

    python scripts/viz_masks3d.py
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from skimage.segmentation import find_boundaries

import germquant.io.nd2_reader as R

REAL = r"data/raw_examples/madeleline images/20251105_N2_noHS/20251105_N2_nohs_HERM _001.nd2"
CY, CX, HY, HX = 1761, 1423, 350, 350
SP = (0.2, 0.108333333333333, 0.108333333333333)
MODELS = [
    ("STOCK cpsam", None),
    ("synthetic-only", "models/models/germline_nuclei"),
    ("real (leak-free)", "models/models/germline_nuclei_real_grp"),
    ("real+synth combined", "models/models/germline_nuclei_combined"),
]


def norm(a):
    lo, hi = np.percentile(a, 1), np.percentile(a, 99.5)
    return np.clip((a - lo) / (hi - lo + 1e-9), 0, 1)


def seg3d(path, dna):
    from cellpose import models
    aniso = SP[0] / SP[1]
    m = models.CellposeModel(gpu=True, pretrained_model=path) if path else models.CellposeModel(gpu=True)
    return np.asarray(m.eval(dna, do_3D=True, z_axis=0, channel_axis=None, anisotropy=aniso)[0]).astype(np.int32)


def main():
    import os
    st = R.read_stack(REAL)
    dna = st.data[0][:, CY - HY:CY + HY, CX - HX:CX + HX].astype(np.float32)
    z = dna.shape[0] // 2
    mip = dna[max(0, z - 4):z + 5].max(0)
    dn = norm(mip)
    avail = [(n, p) for n, p in MODELS if p is None or os.path.exists(p)]
    panels = [("DAPI thin-MIP", None)]
    for name, path in avail:
        lab = seg3d(path, dna)
        labmip = lab[max(0, z - 4):z + 5].max(0)   # any label present in the slab
        print(f"{name}: {lab.max()} 3D nuclei")
        panels.append((f"{name} (n={lab.max()})", labmip))

    fig, ax = plt.subplots(1, len(panels), figsize=(6.5 * len(panels), 6.5))
    for a, (t, lm) in zip(ax, panels):
        a.imshow(dn, cmap="gray"); a.set_title(t, fontsize=11); a.axis("off")
        if lm is not None:
            b = find_boundaries(lm, mode="outer"); ov = np.zeros((*b.shape, 4)); ov[b] = (1, 0.9, 0, 1)
            a.imshow(ov)
    fig.tight_layout(); fig.savefig("results_fullres/masks3d_compare.png", dpi=120, bbox_inches="tight")
    print("wrote results_fullres/masks3d_compare.png")


if __name__ == "__main__":
    main()
