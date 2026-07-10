"""Surface the synaptonemal-complex (SYP) signal into a binary ribbon mask.

This ports the *validated mask core* of the upstream per-nucleus SC tracer
(`C_elegans_ml/src/germquant/sc/skeleton.py::_trace_one`) — the headless-Python equivalent of the
Fiji "Ridge Detection" (Steger) plugin — and STOPS before skeletonization. Per nucleus we resample
the SYP crop to isotropic voxels (so the ridge filter isn't biased by z-anisotropy), enhance the
filament with a Sato (tubeness) ridge filter at physical scales, and build a CONTINUOUS mask with
hysteresis thresholding inside the nucleus (a single high percentile shreds the thin SC into
disconnected blobs — verified upstream on synthetic SC ground truth). The mask is resampled back to
the original anisotropic grid so it aligns voxel-for-voxel with the nucleus labels and the p-granule
mask for colocalization.

  *** WHAT IS / ISN'T RECOVERABLE (validated upstream, scripts/validate_sc_tracer.py) ***
  * SC ribbon MASK + total length — usable (continuous-mask skeleton length ~35 vs 35 µm truth).
  * per-nucleus SC FRAGMENT COUNT — NOT recoverable from light microscopy at pachytene density (the
    ~6 SCs overlap in 3D); use Imaris/SNT for absolute counts. This module never exposes it.

This ribbon is the INTRANUCLEAR SC. The cytoplasmic SYP-aggregate pool that would coincide with
perinuclear P-granules is segmented separately (granule.segment on SYP in the perinuclear shell).
`skan` is only needed for the optional length readout; the mask itself needs only scikit-image.
"""
from __future__ import annotations

import logging

import numpy as np

log = logging.getLogger(__name__)


def surface_sc_ribbon(
    sc_img: np.ndarray,
    labels: np.ndarray,
    spacing: tuple[float, float, float],
    *,
    ridge_sigmas_um: tuple[float, ...] = (0.15, 0.25, 0.40),
    ridge_hyst_low_pct: float = 45.0,
    ridge_hyst_high_pct: float = 80.0,
    keep_nucleus_ids: set[int] | None = None,
) -> np.ndarray:
    """Return a boolean (Z,Y,X) mask of the SC ribbon, aligned to `labels`.

    Runs per-nucleus-crop (percentiles computed inside each nucleus — a global filter is dominated by
    empty space) and confines each crop's mask to the nucleus interior. `keep_nucleus_ids` limits work
    to the germline nuclei; None does all.
    """
    from scipy import ndimage as ndi
    from skimage.filters import sato

    sc_img = np.asarray(sc_img).astype(np.float32)
    labels = np.asarray(labels)
    sp = np.asarray(spacing, dtype=float)
    target = float(sp.min())
    zoom = sp / target
    sigmas_iso = [float(s / target) for s in ridge_sigmas_um]

    out = np.zeros(labels.shape, dtype=bool)
    for nid_minus1, sl in enumerate(ndi.find_objects(labels)):
        nid = nid_minus1 + 1
        if sl is None:
            continue
        if keep_nucleus_ids is not None and nid not in keep_nucleus_ids:
            continue
        sub_lab = labels[sl] == nid
        if sub_lab.sum() < 2:
            continue
        crop_mask = _ribbon_one(
            sc_img[sl], sub_lab, zoom, sigmas_iso, ridge_hyst_low_pct, ridge_hyst_high_pct, ndi, sato)
        if crop_mask is not None:
            out[sl] |= crop_mask & sub_lab
    return out


def _ribbon_one(sub_sc, sub_lab, zoom, sigmas_iso, low_pct, high_pct, ndi, sato):
    """Ridge + hysteresis mask for one nucleus crop, returned on the crop's ORIGINAL grid (or None)."""
    try:
        from skimage.filters import apply_hysteresis_threshold
        from skimage.morphology import disk

        sc_iso = ndi.zoom(sub_sc.astype(np.float32), zoom, order=1)
        lab_iso = ndi.zoom(sub_lab.astype(np.float32), zoom, order=0) > 0.5
        if lab_iso.sum() < 2:
            return None
        ridge = sato(sc_iso, sigmas=sigmas_iso, black_ridges=False)
        interior = ndi.binary_erosion(lab_iso)
        if interior.sum() < 2:
            interior = lab_iso
        pos = ridge[interior]
        pos = pos[pos > 0]
        if pos.size == 0:
            return None
        lo, hi = np.percentile(pos, low_pct), np.percentile(pos, high_pct)
        if hi <= lo:
            mask = (ridge >= hi) & interior
        else:
            mask = apply_hysteresis_threshold(ridge, lo, hi) & interior
        if mask.sum() < 2:
            return None
        # light cleanup: fuse the Sato double-rail into one strand + drop orphan speckle (bincount,
        # version-robust vs skimage.remove_small_objects — same convention as segment/nuclei).
        mask = ndi.binary_closing(mask, disk(1)[None, :, :])
        ol, on = ndi.label(mask)
        if on:
            osz = np.bincount(ol.ravel())
            big = np.where(osz >= 10)[0]
            mask = np.isin(ol, big[big != 0]) & interior
        if mask.sum() < 2:
            return None
        # back to the crop's original (anisotropic) grid so it aligns with labels / granule mask
        return _resize_nn(mask, sub_lab.shape)
    except Exception as e:  # noqa: BLE001 - one bad nucleus must not kill the batch
        log.debug("SC ribbon failed for a nucleus (%s)", e)
        return None


def _resize_nn(mask: np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    """Nearest-neighbour resize of a boolean mask to an exact target shape (crop/pad after zoom to
    guarantee the shape matches even when rounding is off by one voxel)."""
    from scipy import ndimage as ndi

    if mask.shape == tuple(shape):
        return mask
    factor = np.asarray(shape, dtype=float) / np.asarray(mask.shape, dtype=float)
    z = ndi.zoom(mask.astype(np.float32), factor, order=0) > 0.5
    out = np.zeros(shape, dtype=bool)
    sl = tuple(slice(0, min(a, b)) for a, b in zip(shape, z.shape))
    out[sl] = z[sl]
    return out


def sc_ribbon_length(sc_mask: np.ndarray, spacing: tuple[float, float, float]) -> float:
    """Optional total SC ribbon length (µm) via skan skeletonization of the mask. Returns NaN if
    `skan` isn't installed (it's an optional dependency; the mask/coloc never depend on it)."""
    try:
        import skan
        from skimage.morphology import skeletonize
    except Exception:
        return float("nan")
    if sc_mask is None or not sc_mask.any():
        return 0.0
    try:
        skel = skeletonize(np.asarray(sc_mask, dtype=bool))
        if skel.sum() < 2:
            return 0.0
        sk = skan.Skeleton(skel, spacing=tuple(float(s) for s in spacing))
        summary = skan.summarize(sk, separator="_")
        col = "branch_distance" if "branch_distance" in summary else "branch-distance"
        return float(summary[col].sum())
    except Exception as e:  # noqa: BLE001
        log.debug("SC length (skan) failed (%s)", e)
        return float("nan")
