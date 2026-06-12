"""Generic 3D-object (granule) quantification — off the critical path for the N2 SC/foci
experiment, but a promised module (ARCHITECTURE.md stage 7). Threshold (Li/Otsu) -> label ->
spacing-aware volume / surface area / sphericity / intensity, one row per object.

This is the stage that exercises ``measure_objects(..., compute_surface=True)``: granule shape
(sphericity via marching-cubes surface area) is the field-standard readout, where it isn't for
nuclei. Enable via ``config.granules.enabled`` with a channel map that defines a ``granule`` role.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .measure import measure_objects
from .schema import GRANULES


def detect_granules(
    img: np.ndarray,
    spacing: tuple[float, float, float],
    *,
    labels: np.ndarray | None = None,
    marker: str = "granule",
    threshold_method: str = "li",
    min_volume_um3: float = 0.05,
) -> pd.DataFrame:
    """Return a tidy granules table (schema.GRANULES). ``labels`` (nucleus mask) is optional;
    when given, each granule is assigned the nucleus_id at its centroid."""
    from scipy import ndimage as ndi
    from skimage.filters import threshold_li, threshold_otsu

    imgf = img.astype(np.float32)
    thr_fn = {"li": threshold_li, "otsu": threshold_otsu}.get(threshold_method, threshold_li)
    if not np.isfinite(imgf).any() or float(imgf.max() - imgf.min()) <= 0:
        return pd.DataFrame(columns=GRANULES)
    fg = imgf > thr_fn(imgf)
    lab, n = ndi.label(fg)
    if n == 0:
        return pd.DataFrame(columns=GRANULES)

    df = measure_objects(lab, {"granule": img}, spacing, compute_surface=True)
    if df.empty:
        return pd.DataFrame(columns=GRANULES)
    df = df[df["volume_um3"] >= float(min_volume_um3)].copy()
    if df.empty:
        return pd.DataFrame(columns=GRANULES)

    sp = np.asarray(spacing, float)
    def _nucleus_of(row):
        if labels is None:
            return 0
        z = int(np.clip(round(row["centroid_z_um"] / sp[0]), 0, labels.shape[0] - 1))
        y = int(np.clip(round(row["centroid_y_um"] / sp[1]), 0, labels.shape[1] - 1))
        x = int(np.clip(round(row["centroid_x_um"] / sp[2]), 0, labels.shape[2] - 1))
        return int(labels[z, y, x])

    out = pd.DataFrame({
        "granule_id": np.arange(1, len(df) + 1),
        "nucleus_id": df.apply(_nucleus_of, axis=1).to_numpy() if len(df) else [],
        "marker": marker,
        "volume_um3": df["volume_um3"].to_numpy(),
        "surface_area_um2": df["surface_area_um2"].to_numpy(),
        "n_voxels": df["n_voxels"].to_numpy(),
        "sphericity": df["sphericity"].to_numpy(),
        "intensity_mean": df["granule_mean_intensity"].to_numpy(),
        "intensity_integrated": df["granule_integrated_intensity"].to_numpy(),
        "intensity_max": df["granule_max_intensity"].to_numpy(),
        "centroid_z_um": df["centroid_z_um"].to_numpy(),
        "centroid_y_um": df["centroid_y_um"].to_numpy(),
        "centroid_x_um": df["centroid_x_um"].to_numpy(),
        "threshold_method": threshold_method,
    })
    return out[GRANULES]
