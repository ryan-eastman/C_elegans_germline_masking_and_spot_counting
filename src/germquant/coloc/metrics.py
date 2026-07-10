"""Direct voxel/object colocalization between an SC (SYP) mask and the PGL-1 granule mask.

Everything is computed WITHIN a supplied region mask R (the germline dilated by a perinuclear shell —
see the pipeline), because coloc over the whole image is dominated by empty background and off-gonad
noise. The HEADLINE metrics are object/mask overlap plus a translation-null test; Manders M1/M2 are
mask-restricted intensity fractions; Pearson r is reported as a region-restricted DIAGNOSTIC only
(it is inflated by the large shared near-zero background of two sparse channels and must not be used
as the coincidence claim). Costes auto-thresholding is optional and off by default.

Pure function: no config parsing, no I/O. Physical volumes come from voxel counts × voxel volume.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

PER_GRANULE_COLS = ["granule_id", "overlap_frac", "overlaps", "nearest_um"]


def colocalize(
    sc_mask: np.ndarray,
    granule_mask: np.ndarray,
    granule_labels: np.ndarray,
    sc_img: np.ndarray,
    granule_img: np.ndarray,
    region_mask: np.ndarray,
    spacing: tuple[float, float, float],
    *,
    sc_operand: str = "syp_aggregate",
    region_name: str = "perinuclear_shell",
    region_dilation_um: float = 1.5,
    n_random: int = 100,
    object_overlap_min_frac: float = 0.0,
    costes: bool = False,
    rng_seed: int = 0,
) -> tuple[dict, pd.DataFrame]:
    """Return (coloc_row_dict, per_granule_df).

    `object_overlap_min_frac` is the fraction of a granule's voxels that must fall in the SC mask for
    it to count as "overlapping" (0 = any single shared voxel counts). `n_random` sets the size of the
    translation null used for the overlap p-value/z-score (0 disables it → NaN)."""
    sp = np.asarray(spacing, dtype=float)
    vox_vol = float(sp[0] * sp[1] * sp[2])
    region = np.asarray(region_mask, dtype=bool)
    sc = np.asarray(sc_mask, dtype=bool) & region
    gr = np.asarray(granule_mask, dtype=bool) & region
    glab = np.asarray(granule_labels)

    inter = sc & gr
    n_sc, n_gr, n_int = int(sc.sum()), int(gr.sum()), int(inter.sum())
    union = n_sc + n_gr - n_int

    row: dict = {
        "sc_operand": sc_operand,
        "region": region_name,
        "region_dilation_um": float(region_dilation_um),
        "region_voxels": int(region.sum()),
        "region_volume_um3": float(region.sum() * vox_vol),
        "overlap_volume_um3": float(n_int * vox_vol),
        "dice": float(2 * n_int / (n_sc + n_gr)) if (n_sc + n_gr) else float("nan"),
        "jaccard": float(n_int / union) if union else float("nan"),
        "frac_sc_in_granules": float(n_int / n_sc) if n_sc else float("nan"),
        "frac_granule_in_sc": float(n_int / n_gr) if n_gr else float("nan"),
        "manders_m1": _manders(sc_img, gr, region),      # frac of SYP intensity inside the PGL mask
        "manders_m2": _manders(granule_img, sc, region),  # frac of PGL intensity inside the SYP mask
        "pearson_r": _pearson(sc_img, granule_img, region),
        "sc_mask_voxels": n_sc,
        "granule_mask_voxels": n_gr,
        "n_random": int(n_random),
        "costes_threshold_syp": float("nan"),
        "costes_threshold_pgl": float("nan"),
    }

    # per-granule overlap fraction + nearest-SC distance
    per_g, n_over, dists = _per_granule(glab, sc, sp, object_overlap_min_frac)
    n_granules = int(len(per_g))
    row["n_granules"] = n_granules
    row["n_granules_overlapping_sc"] = int(n_over)
    row["frac_granules_overlapping_sc"] = float(n_over / n_granules) if n_granules else float("nan")
    finite = dists[np.isfinite(dists)] if dists.size else dists
    row["mean_granule_to_sc_um"] = float(finite.mean()) if finite.size else float("nan")
    row["median_granule_to_sc_um"] = float(np.median(finite)) if finite.size else float("nan")

    # translation null: is the observed overlap more than chance given the two masks' sizes?
    pval, z = _overlap_null(sc, gr, n_int, n_random, rng_seed)
    row["overlap_pvalue"] = pval
    row["overlap_zscore"] = z

    if costes:
        try:
            tsyp, tpgl = _costes_threshold(sc_img, granule_img, region)
            row["costes_threshold_syp"] = tsyp
            row["costes_threshold_pgl"] = tpgl
        except Exception:  # noqa: BLE001 - optional enrichment
            pass

    return row, per_g


def _manders(intensity: np.ndarray, other_mask: np.ndarray, region: np.ndarray) -> float:
    """Fraction of `intensity` (within `region`) that falls inside `other_mask`. Non-negative
    intensities assumed (fluorescence); a zero floor guards a rare negative background."""
    a = np.asarray(intensity, dtype=np.float64)
    a = np.clip(a, 0, None)
    tot = float(a[region].sum())
    if tot <= 0:
        return float("nan")
    return float(a[other_mask & region].sum() / tot)


def _pearson(a_img: np.ndarray, b_img: np.ndarray, region: np.ndarray) -> float:
    a = np.asarray(a_img, dtype=np.float64)[region]
    b = np.asarray(b_img, dtype=np.float64)[region]
    if a.size < 2:
        return float("nan")
    sa, sb = a.std(), b.std()
    if sa == 0 or sb == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _per_granule(granule_labels, sc, spacing, min_frac):
    """Per-granule overlap fraction (voxels in SC / granule voxels), overlap bool, and nearest-SC
    surface distance (µm). Distances via a physical EDT of the SC-mask complement."""
    from scipy import ndimage as ndi

    ids = [int(v) for v in np.unique(granule_labels) if v != 0]
    if not ids:
        return pd.DataFrame(columns=PER_GRANULE_COLS), 0, np.array([])

    if sc.any():
        dt = ndi.distance_transform_edt(~sc, sampling=tuple(float(s) for s in spacing))
    else:
        dt = None
    sizes = ndi.sum(np.ones_like(granule_labels, dtype=np.float32), granule_labels, index=ids)
    in_sc = ndi.sum(sc.astype(np.float32), granule_labels, index=ids)
    rows, dists = [], []
    n_over = 0
    for gid, sz, insc in zip(ids, np.atleast_1d(sizes), np.atleast_1d(in_sc)):
        frac = float(insc / sz) if sz else 0.0
        overlaps = (insc > 0) if min_frac <= 0 else (frac >= min_frac)
        n_over += int(bool(overlaps))
        if dt is not None:
            nearest = float(dt[granule_labels == gid].min())
        else:
            nearest = float("nan")
        rows.append({"granule_id": gid, "overlap_frac": frac,
                     "overlaps": bool(overlaps), "nearest_um": nearest})
        dists.append(nearest)
    return pd.DataFrame(rows, columns=PER_GRANULE_COLS), n_over, np.asarray(dists, dtype=float)


def _overlap_null(sc, gr, observed_vox, n_random, rng_seed):
    """Translation null: keep the SC mask fixed, circularly shift the granule mask to random
    positions N times, and compare the observed SC∩granule voxel count to that null distribution.
    Tests whether the coincidence exceeds chance given each mask's size within the volume. Returns
    (p_value, z_score); NaN if disabled or either mask is empty."""
    if n_random <= 0 or not sc.any() or not gr.any():
        return float("nan"), float("nan")
    rng = np.random.default_rng(rng_seed)
    shape = gr.shape
    null = np.empty(n_random, dtype=np.int64)
    for i in range(n_random):
        shifts = [int(rng.integers(0, s)) for s in shape]
        rolled = np.roll(gr, shift=shifts, axis=(0, 1, 2))
        null[i] = int((sc & rolled).sum())
    p = float((1 + int((null >= observed_vox).sum())) / (1 + n_random))
    mu, sd = float(null.mean()), float(null.std())
    z = float((observed_vox - mu) / sd) if sd > 0 else float("nan")
    return p, z


def _costes_threshold(a_img, b_img, region):
    """Costes automatic threshold (optional): the largest intensity cutoff along the a→b regression
    line below which the two channels are no longer positively correlated. Returns (T_a, T_b)."""
    a = np.asarray(a_img, dtype=np.float64)[region]
    b = np.asarray(b_img, dtype=np.float64)[region]
    if a.size < 3 or a.std() == 0 or b.std() == 0:
        return float("nan"), float("nan")
    slope = np.cov(a, b)[0, 1] / np.var(a)
    intercept = b.mean() - slope * a.mean()
    ta = float(a.max())
    step = (a.max() - a.min()) / 100.0 or 1.0
    while ta > a.min():
        tb = slope * ta + intercept
        below = (a < ta) & (b < tb)
        if below.sum() < 2 or a[below].std() == 0 or b[below].std() == 0:
            ta -= step
            continue
        if np.corrcoef(a[below], b[below])[0, 1] <= 0:
            break
        ta -= step
    return float(ta), float(slope * ta + intercept)
