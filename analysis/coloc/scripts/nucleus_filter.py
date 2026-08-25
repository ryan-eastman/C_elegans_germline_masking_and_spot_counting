"""Per-nucleus filters for the lamin-masked analysis (Ryan, 2026-08-24; from the mask audit).

The pipeline's `in_germline` flag is DAPI-only and lets two kinds of non-germline objects through:
  1. objects with no nuclear envelope: sperm / spermatids, condensed debris, the sperm mass at the
     proximal end of males and in the herm spermatheca. Test: LMN-1 just outside the label is flat
     against the interior (ring_ratio > 0.97) AND well below the ring threshold (shell < 0.75 x Otsu
     ring level of the crop). Real germline nuclei score ratio ~0.85 and shell ~1.0 x threshold.
  2. ring-bearing nuclei that are not part of the traced gonad: a second gonad arm, somatic tissue,
     another worm. Test: the germline labels are split into 2D territories (max projection dilated
     by 4 um, connected components); only the territory that Ryan's hand-traced pachytene polyline
     runs through is kept.
Both are number-changing for the whole-gonad tables; the zone tables already use the trace and its
off-axis cutoff, so they are affected only by filter 1."""
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu

import chload

SP = chload.SP
RING_RATIO_MAX = 0.97      # shell / inside above this = no ring
SHELL_OVER_THR_MAX = 0.75  # shell below this fraction of the Otsu ring level = no ring
TERRITORY_DILATE_UM = 4.0


def ring_scores(lab_c, germ_ids, lamin):
    """DataFrame indexed by nucleus id: ring_ratio, shell_over_thr, plus the crop ring threshold."""
    ids = np.asarray(germ_ids, dtype=np.int64)
    gn = np.isin(lab_c, ids)
    sm = ndi.gaussian_filter(lamin.astype(np.float32), sigma=0.25 / SP)
    thr = float(threshold_otsu(sm[sm > np.percentile(sm, 50)]))
    dt_out, idx = ndi.distance_transform_edt(~gn, sampling=tuple(SP), return_indices=True)
    shell = (dt_out > 0) & (dt_out <= 0.4)
    near = lab_c[idx[0], idx[1], idx[2]]
    shell_lab = np.where(shell, near, 0)
    inside = np.asarray(ndi.mean(sm, np.where(gn, lab_c, 0), ids))
    outside = np.asarray(ndi.mean(sm, shell_lab, ids))
    df = pd.DataFrame({"ring_ratio": outside / np.maximum(inside, 1.0), "shell_over_thr": outside / thr}, index=ids)
    df.attrs["ring_thr"] = thr
    return df


def no_ring_ids(scores):
    m = (scores.ring_ratio > RING_RATIO_MAX) & (scores.shell_over_thr < SHELL_OVER_THR_MAX)
    return [int(i) for i in scores.index[m]]


def territory_ids(lab_c, germ_ids, poly_um, nc=None):
    """ids whose centroid lies in the 2D territory crossed by the trace polyline (crop-um x, y).
    Returns (kept_ids, n_territories, note). Falls back to the largest territory if the polyline
    misses every territory (note says so)."""
    ids = np.asarray(germ_ids, dtype=np.int64)
    gn = np.isin(lab_c, ids)
    mp = gn.max(0)
    r = int(round(TERRITORY_DILATE_UM / SP[1]))
    yy, xx = np.ogrid[-r:r + 1, -r:r + 1]
    comp, n = ndi.label(ndi.binary_dilation(mp, structure=(yy ** 2 + xx ** 2) <= r * r))
    # sample the polyline densely and collect the territories it crosses
    pts = np.asarray(poly_um, dtype=float)
    hit = set()
    for a, b in zip(pts[:-1], pts[1:]):
        for t in np.linspace(0, 1, 200):
            x, y = a + t * (b - a)
            yi, xi = int(round(y / SP[1])), int(round(x / SP[2]))
            if 0 <= yi < comp.shape[0] and 0 <= xi < comp.shape[1] and comp[yi, xi] > 0:
                hit.add(int(comp[yi, xi]))
    note = ""
    if not hit:
        sizes = np.bincount(comp.ravel())
        sizes[0] = 0
        hit = {int(np.argmax(sizes))}
        note = "polyline missed all territories; kept the largest"
    coms = ndi.center_of_mass(gn, lab_c, ids)
    keep = []
    for i, c in zip(ids, coms):
        yi, xi = int(round(c[1])), int(round(c[2]))
        yi = min(max(yi, 0), comp.shape[0] - 1)
        xi = min(max(xi, 0), comp.shape[1] - 1)
        if comp[yi, xi] in hit:
            keep.append(int(i))
    return keep, int(n), note
