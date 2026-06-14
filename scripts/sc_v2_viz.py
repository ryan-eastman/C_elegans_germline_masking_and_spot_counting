"""Visual before/after for the v2 SC tracer: SYP vs OLD p90-threshold mask vs NEW hysteresis mask
vs NEW skeleton. Confirms the v2 mask is a continuous strand, not the old disconnected blobs."""
import os
import sys

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage as ndi
from skimage.filters import apply_hysteresis_threshold, sato
from skimage.morphology import skeletonize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from synth_nuclei import DX, DY, DZ, degrade, synthesize_scene  # noqa: E402

SP = np.array([DZ, DY, DX]); TARGET = float(SP.min()); ZOOM = SP / TARGET


def norm(a):
    m = a.max()
    return a / m if m > 0 else a


def masks(sub_sc, sub_lab):
    sc_iso = ndi.zoom(sub_sc, ZOOM, order=1)
    lab = ndi.zoom(sub_lab.astype(np.float32), ZOOM, order=0) > 0.5
    ridge = sato(sc_iso, sigmas=[s / TARGET for s in (0.15, 0.25, 0.40)], black_ridges=False)
    interior = ndi.binary_erosion(lab)
    if interior.sum() < 2:
        interior = lab
    pos = ridge[interior]
    pos = pos[pos > 0]
    old = (ridge >= np.percentile(pos, 90)) & interior
    new = apply_hysteresis_threshold(ridge, np.percentile(pos, 45), np.percentile(pos, 80)) & interior
    new = _smooth_mask(new, interior)        # production mask-smoothing (matches sc/skeleton._trace_one)
    return sc_iso, old, new, skeletonize(new)


def _smooth_mask(mask, interior):
    from skimage.morphology import disk
    mask = ndi.binary_closing(mask, disk(1)[None, :, :])
    holes = ndi.binary_fill_holes(mask) & ~mask
    hl, hn = ndi.label(holes)
    if hn:
        hsz = np.bincount(hl.ravel())
        small = np.where(hsz < 27)[0]
        mask = mask | np.isin(hl, small[small != 0])
    ol, on = ndi.label(mask)
    if on:
        osz = np.bincount(ol.ravel())
        big = np.where(osz >= 30)[0]
        mask = np.isin(ol, big[big != 0]) & interior
    return mask


def main():
    dapi, syp, label, gt = synthesize_scene(seed=3, frag_rate=0.0)
    sypimg = degrade(syp, seed=30)
    objs = ndi.find_objects(label)
    nids = [1, 2, 3]
    fig, ax = plt.subplots(3, 4, figsize=(16, 12))
    for r, nid in enumerate(nids):
        sl = objs[nid - 1]
        sc_iso, old, new, skel = masks(sypimg[sl], label[sl] == nid)
        z = sc_iso.shape[0] // 2

        def mip(v):
            return v[max(0, z - 4):z + 5].max(0)
        for c, (im, t) in enumerate([(norm(mip(sc_iso)), "SYP MIP"), (mip(old), "OLD p90 mask"),
                                     (mip(new), "NEW hysteresis mask"), (mip(skel), "NEW skeleton")]):
            ax[r, c].imshow(im, cmap="gray"); ax[r, c].axis("off"); ax[r, c].set_title(f"nuc{nid} | {t}", fontsize=9)
    fig.suptitle("SC tracer v2: hysteresis mask recovers continuous strands (old p90 = blobs)", fontsize=12)
    fig.tight_layout(); fig.savefig("results_fullres/sc_v2_compare.png", dpi=120, bbox_inches="tight")
    print("wrote results_fullres/sc_v2_compare.png")


if __name__ == "__main__":
    main()
