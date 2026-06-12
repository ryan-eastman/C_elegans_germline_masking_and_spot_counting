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

    shape = np.asarray(img.shape)
    rows = []
    for i, b in enumerate(blobs):
        # round (not floor) the sub-voxel center so foci near a nucleus edge are assigned
        # to the nearest voxel rather than always the lower one.
        z, y, x = (int(v) for v in np.clip(np.round(b[:3]), 0, shape - 1))
        sig_vox = b[3:6] if blobs.shape[1] >= 6 else np.array([b[3]] * 3)
        radius_um = float(np.sqrt(3) * np.mean(sig_vox * sp))   # LoG radius ≈ sqrt(ndim)*sigma
        nid = int(labels[z, y, x]) if labels is not None else 0
        # sample intensity over the spot, not a single voxel: a box ~1 sigma in each axis.
        half = np.maximum(1, np.round(sig_vox)).astype(int)
        zsl = slice(max(0, z - half[0]), min(shape[0], z + half[0] + 1))
        ysl = slice(max(0, y - half[1]), min(shape[1], y + half[1] + 1))
        xsl = slice(max(0, x - half[2]), min(shape[2], x + half[2] + 1))
        patch = img[zsl, ysl, xsl]
        rows.append(
            {
                "focus_id": i,
                "nucleus_id": nid,
                "marker": marker,
                "z_um": z * sp[0],
                "y_um": y * sp[1],
                "x_um": x * sp[2],
                "sigma_um": float(np.mean(sig_vox * sp)),
                "intensity_mean": float(patch.mean()),
                "intensity_max": float(patch.max()),
                "detector": "blob_log",
                "detection_radius_um": radius_um,
            }
        )
    return pd.DataFrame(rows, columns=_COLS)


_COLS = [
    "focus_id", "nucleus_id", "marker", "z_um", "y_um", "x_um", "sigma_um",
    "intensity_mean", "intensity_max", "detector", "detection_radius_um",
]
