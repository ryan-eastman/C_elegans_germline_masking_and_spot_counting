"""Spacing-aware 3D object measurement.

Physical quantities are computed EXPLICITLY from voxel counts × voxel volume and from
spacing-scaled centroids, so results are correct regardless of scikit-image version
behaviour around the `spacing=` kwarg. Surface area uses marching_cubes(spacing=...).
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _voxel_volume(spacing: tuple[float, float, float]) -> float:
    return float(spacing[0] * spacing[1] * spacing[2])


def surface_area_um2(mask: np.ndarray, spacing: tuple[float, float, float]) -> float:
    """Marching-cubes surface area in µm² (physically correct only with spacing)."""
    from skimage.measure import marching_cubes, mesh_surface_area

    if mask.sum() == 0:
        return float("nan")
    try:
        verts, faces, _, _ = marching_cubes(mask.astype(np.uint8), level=0.5, spacing=spacing)
        return float(mesh_surface_area(verts, faces))
    except (ValueError, RuntimeError):
        return float("nan")


def measure_objects(
    label_img: np.ndarray,
    intensity: dict[str, np.ndarray],
    spacing: tuple[float, float, float],
    *,
    compute_surface: bool = False,
) -> pd.DataFrame:
    """One row per labelled object with µm³ volume, µm centroid, and per-channel intensity.

    Parameters
    ----------
    label_img : (Z, Y, X) int label image
    intensity : {channel_role: (Z, Y, X) image} for intensity stats
    spacing   : (dz, dy, dx) microns — REQUIRED
    """
    from skimage.measure import regionprops

    vox_vol = _voxel_volume(spacing)
    sp = np.asarray(spacing, dtype=float)
    rows = []
    for rp in regionprops(label_img):
        n_vox = int(rp.area)  # voxel count (unscaled) regardless of skimage version
        cz, cy, cx = (np.asarray(rp.centroid) * sp).tolist()
        row = {
            "label": int(rp.label),
            "n_voxels": n_vox,
            "volume_um3": n_vox * vox_vol,
            "centroid_z_um": cz,
            "centroid_y_um": cy,
            "centroid_x_um": cx,
        }
        coords = rp.coords  # (n_vox, 3) voxel indices
        for role, img in intensity.items():
            vals = img[coords[:, 0], coords[:, 1], coords[:, 2]]
            row[f"{role}_mean_intensity"] = float(vals.mean())
            row[f"{role}_max_intensity"] = float(vals.max())
            row[f"{role}_integrated_intensity"] = float(vals.sum())
        if compute_surface:
            sa = surface_area_um2(label_img[rp.slice] == rp.label, spacing)
            row["surface_area_um2"] = sa
            # sphericity = (pi^(1/3) (6V)^(2/3)) / A
            v = row["volume_um3"]
            row["sphericity"] = (
                float((np.pi ** (1 / 3)) * ((6 * v) ** (2 / 3)) / sa)
                if sa and not np.isnan(sa) and sa > 0
                else float("nan")
            )
        rows.append(row)
    return pd.DataFrame(rows)
