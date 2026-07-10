"""Segment p-granule (PGL-1) puncta as 3D objects inside a region mask.

P-granules are perinuclear condensates (bright, roughly-round puncta of ~0.2-1.5 µm) that dock on
the cytoplasmic face of the germ-cell nuclear envelope. We segment them the same spacing-aware way
the classical nucleus fallback works (`segment/nuclei.py::_classical`): light gaussian -> a global
threshold computed WITHIN the region -> connected components -> a physical size gate. The threshold
default is `threshold_triangle`, matching the cross-validated SpotMAX spots default (robust to
per-image brightness). The same function is reused on the SYP channel (restricted to the perinuclear
shell) to segment cytoplasmic SYP aggregates for colocalization.

Returns a label image on the input voxel grid (-> Imaris Surfaces) and a tidy per-granule table
(centroid µm, volume µm³, intensities) via the shared `measure.measure_objects`.
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

GRANULE_COLS = [
    "granule_id", "marker", "z_um", "y_um", "x_um",
    "volume_um3", "n_voxels", "granule_mean_intensity", "granule_max_intensity", "detector",
]

_THRESHOLDS = {
    "threshold_triangle": "threshold_triangle",
    "threshold_otsu": "threshold_otsu",
    "threshold_li": "threshold_li",
    "triangle": "threshold_triangle",
    "otsu": "threshold_otsu",
    "li": "threshold_li",
}


def segment_granules(
    img: np.ndarray,
    region_mask: np.ndarray,
    spacing: tuple[float, float, float],
    *,
    marker: str = "PGL-1",
    thresholding_method: str = "threshold_triangle",
    gauss_sigma_um: float = 0.1,
    min_volume_um3: float = 0.03,
    max_volume_um3: float = 8.0,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Return (granule_labels (Z,Y,X) int32, per_granule_df).

    Objects are found only where `region_mask` is True (the dilated germline / perinuclear shell), so
    off-gonad autofluorescence is excluded. The threshold is computed over the smoothed intensities
    inside the region only. Components outside [min_volume_um3, max_volume_um3] are dropped (the max
    gate rejects large blobs that are not discrete granules, e.g. bright autofluor or a merged sheet).
    """
    from scipy import ndimage as ndi

    img = np.asarray(img).astype(np.float32)
    region = np.asarray(region_mask, dtype=bool)
    sp = np.asarray(spacing, dtype=float)
    if not region.any():
        return np.zeros(img.shape, np.int32), pd.DataFrame(columns=GRANULE_COLS)

    sigma_vox = tuple(float(gauss_sigma_um / s) for s in sp)
    sm = ndi.gaussian_filter(img, sigma=sigma_vox)

    vals = sm[region]
    vals = vals[np.isfinite(vals)]
    thr = _threshold(vals, thresholding_method)
    if thr is None:
        return np.zeros(img.shape, np.int32), pd.DataFrame(columns=GRANULE_COLS)

    fg = (sm > thr) & region
    labels, _ = ndi.label(fg)
    labels = _size_gate(labels, sp, min_volume_um3, max_volume_um3)

    df = _measure(labels, img, sp, marker)
    return labels.astype(np.int32), df


def _threshold(vals: np.ndarray, method: str) -> float | None:
    """Global intensity threshold over the in-region values. Returns None if it can't be computed
    (empty / flat region), so the caller yields zero granules rather than crashing."""
    import skimage.filters as F

    if vals.size == 0 or float(vals.max()) <= float(vals.min()):
        return None
    fn_name = _THRESHOLDS.get(str(method).lower(), "threshold_triangle")
    try:
        return float(getattr(F, fn_name)(vals))
    except Exception as e:  # noqa: BLE001 - degenerate histogram; skip granules this image
        log.warning("granule threshold (%s) failed (%s); no granules this image.", method, e)
        return None


def _size_gate(labels: np.ndarray, spacing: np.ndarray, min_um3: float, max_um3: float) -> np.ndarray:
    """Keep connected components whose volume is within [min_um3, max_um3], relabel consecutively.

    Uses np.bincount directly (version-robust vs skimage.remove_small_objects, whose min_size
    semantics changed in 0.26 — same convention as segment/nuclei._drop_small)."""
    from skimage.segmentation import relabel_sequential

    lab = labels.astype(np.int32)
    vox_vol = float(spacing[0] * spacing[1] * spacing[2])
    min_vox = max(1, int(min_um3 / vox_vol))
    max_vox = int(max_um3 / vox_vol) if max_um3 and max_um3 > 0 else np.iinfo(np.int64).max
    counts = np.bincount(lab.ravel())
    drop = np.where((counts < min_vox) | (counts > max_vox))[0]
    drop = drop[drop != 0]
    if drop.size:
        lab[np.isin(lab, drop)] = 0
    lab, _, _ = relabel_sequential(lab)
    return lab.astype(np.int32)


def _measure(labels: np.ndarray, img: np.ndarray, spacing: np.ndarray, marker: str) -> pd.DataFrame:
    """Per-granule table via the shared measure_objects (µm³ volume, µm centroid, intensity)."""
    from ..measure import measure_objects

    m = measure_objects(labels, {"granule": img}, tuple(float(s) for s in spacing))
    if m.empty:
        return pd.DataFrame(columns=GRANULE_COLS)
    out = pd.DataFrame({
        "granule_id": m["label"].astype(int),
        "marker": marker,
        "z_um": m["centroid_z_um"].to_numpy(),
        "y_um": m["centroid_y_um"].to_numpy(),
        "x_um": m["centroid_x_um"].to_numpy(),
        "volume_um3": m["volume_um3"].to_numpy(),
        "n_voxels": m["n_voxels"].astype(int),
        "granule_mean_intensity": m["granule_mean_intensity"].to_numpy(),
        "granule_max_intensity": m["granule_max_intensity"].to_numpy(),
        "detector": "threshold_cc",
    }, columns=GRANULE_COLS)
    return out
