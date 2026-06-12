"""3D foci detection (RAD-51, or crossover markers COSA-1/MSH-5/ZHP-3).

Laplacian-of-Gaussian blob detection with per-axis sigma derived from PHYSICAL size
(µm) via voxel spacing — so the detector looks for the same real spot size regardless
of z-anisotropy. Each focus is assigned to the nucleus label it falls in.

v1 baseline. For sub-pixel rigor on dense fields, swap in big-fish / RS-FISH (see
ARCHITECTURE.md stage 6) — the interface (returns a tidy DataFrame) stays the same.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def detect_foci(
    img: np.ndarray,
    spacing: tuple[float, float, float],
    *,
    marker: str = "RAD-51",
    min_sigma_um: float = 0.10,
    max_sigma_um: float = 0.35,
    threshold_rel: float = 0.10,
    labels: np.ndarray | None = None,
) -> pd.DataFrame:
    from skimage.feature import blob_log

    sp = np.asarray(spacing, dtype=float)
    min_sigma = (min_sigma_um / sp).tolist()      # per-axis (z, y, x) in voxels
    max_sigma = (max_sigma_um / sp).tolist()

    imgf = img.astype(np.float32)
    rng = float(imgf.max() - imgf.min())
    if rng <= 0:
        return pd.DataFrame(columns=_COLS)
    norm = (imgf - imgf.min()) / rng

    blobs = blob_log(norm, min_sigma=min_sigma, max_sigma=max_sigma, threshold_rel=threshold_rel)
    if blobs.size == 0:
        return pd.DataFrame(columns=_COLS)

    rows = []
    for i, b in enumerate(blobs):
        z, y, x = int(b[0]), int(b[1]), int(b[2])
        sig_vox = b[3:6] if blobs.shape[1] >= 6 else np.array([b[3]] * 3)
        radius_um = float(np.sqrt(3) * np.mean(sig_vox * sp))   # LoG radius ≈ sqrt(ndim)*sigma
        nid = int(labels[z, y, x]) if labels is not None else 0
        rows.append(
            {
                "focus_id": i,
                "nucleus_id": nid,
                "marker": marker,
                "z_um": z * sp[0],
                "y_um": y * sp[1],
                "x_um": x * sp[2],
                "sigma_um": float(np.mean(sig_vox * sp)),
                "intensity_mean": float(img[z, y, x]),
                "intensity_max": float(img[z, y, x]),
                "detector": "blob_log",
                "detection_radius_um": radius_um,
            }
        )
    return pd.DataFrame(rows, columns=_COLS)


_COLS = [
    "focus_id", "nucleus_id", "marker", "z_um", "y_um", "x_um", "sigma_um",
    "intensity_mean", "intensity_max", "detector", "detection_radius_um",
]
