"""Read objects directly from an Imaris .ims (HDF5) project — Spots (RAD-51 foci) and Surfaces
(nuclei) — in IMAGE-RELATIVE microns, so they can be matched against our pipeline output for the
exact same image (our centroids are voxel*spacing, origin at the image corner).

Imaris stores stage coordinates offset by the image extent minimum (DataSetInfo/Image ExtMin{0,1,2}
for X/Y/Z). Subtracting ExtMin puts objects in the same origin-at-corner frame as our pipeline.
Surface VOXEL masks live in Imaris's proprietary MegaSurfaces block tree (not extracted here); we
read per-surface centroid + volume from SurfaceModelInfo, which is enough for count/centroid checks.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

_SPOTS = "Scene8/Content/Points0/Spot"
_SURF = "Scene8/Content/MegaSurfaces0/SurfaceModelInfo"


def _fattr(group, key: str) -> float:
    v = group.attrs[key]
    if isinstance(v, np.ndarray) and v.dtype.kind == "S":
        return float(v.tobytes().decode("latin1").strip("\x00"))
    if isinstance(v, np.ndarray):
        return float(v.ravel()[0])
    return float(v)


def read_imaris_ims(path: str | Path):
    """Return (spots_df, surfaces_df, geom). spots_df: x_um/y_um/z_um/radius_um (image-relative).
    surfaces_df: x_um/y_um/z_um/volume_um3 (image-relative). geom: extent + shape + voxel size."""
    import h5py

    with h5py.File(str(path), "r") as f:
        info = f["DataSetInfo"]["Image"]
        ext_min = np.array([_fattr(info, f"ExtMin{i}") for i in range(3)])  # x, y, z (µm)
        ext_max = np.array([_fattr(info, f"ExtMax{i}") for i in range(3)])
        shape_xyz = np.array([int(_fattr(info, k)) for k in ("X", "Y", "Z")])
        spot = f[_SPOTS][:] if _SPOTS in f else None
        smi = f[_SURF][:] if _SURF in f else None

    spots = pd.DataFrame(columns=["x_um", "y_um", "z_um", "radius_um"])
    if spot is not None and len(spot):
        spots = pd.DataFrame({
            "x_um": spot["PositionX"] - ext_min[0],
            "y_um": spot["PositionY"] - ext_min[1],
            "z_um": spot["PositionZ"] - ext_min[2],
            "radius_um": spot["Radius"],
        })
    surfaces = pd.DataFrame(columns=["x_um", "y_um", "z_um", "volume_um3"])
    if smi is not None and len(smi):
        surfaces = pd.DataFrame({
            "x_um": smi["CenterOfMassX"] - ext_min[0],
            "y_um": smi["CenterOfMassY"] - ext_min[1],
            "z_um": smi["CenterOfMassZ"] - ext_min[2],
            "volume_um3": smi["Volume"],
        })
    voxel_xyz = (ext_max - ext_min) / np.maximum(shape_xyz, 1)
    geom = {"ext_min_xyz": ext_min, "ext_max_xyz": ext_max,
            "shape_xyz": shape_xyz, "voxel_xyz_um": voxel_xyz}
    return spots, surfaces, geom


def spots_per_surface(spots: pd.DataFrame, surfaces: pd.DataFrame, max_dist_um: float = 3.0,
                      use_surface_radius: bool = True) -> np.ndarray:
    """Assign each spot to the nearest surface centroid and count per surface. By default the cutoff
    is the surface's own sphere-equivalent radius from `volume_um3` (capped at max_dist_um) rather
    than a flat max_dist_um — a flat 3 um exceeds the ~2 um nucleus radius and over-counts. This is
    still an approximation (centroid-based, no voxel masks) — the lab's coloc method (imaris_xlsx) is
    canonical for per-nucleus counts."""
    if surfaces.empty:
        return np.array([], dtype=int)
    sc = surfaces[["x_um", "y_um", "z_um"]].to_numpy()
    if use_surface_radius and "volume_um3" in surfaces.columns:
        cutoff = np.minimum((3 * surfaces["volume_um3"].to_numpy() / (4 * np.pi)) ** (1 / 3), max_dist_um)
    else:
        cutoff = np.full(len(surfaces), max_dist_um)
    counts = np.zeros(len(surfaces), dtype=int)
    if spots.empty:
        return counts
    sp = spots[["x_um", "y_um", "z_um"]].to_numpy()
    for p in sp:
        d = np.linalg.norm(sc - p, axis=1)
        j = int(d.argmin())
        if d[j] <= cutoff[j]:
            counts[j] += 1
    return counts
