"""Render detected spots as a 3D image for viewing / importing in Imaris — a channel of small blobs
on the SAME voxel grid as the nucleus-label TIF and the original image, so it overlays directly.
Load it in Imaris as an extra Channel, or run Imaris' Spots detection on it (each blob -> one Spot).
"""
from __future__ import annotations

import numpy as np


def spots_to_image(spots_df, shape, spacing, radius_um: float = 0.3, label: bool = False) -> np.ndarray:
    """uint16 (Z,Y,X) image with a small ellipsoid blob at each spot (anisotropy-aware, ~equal physical
    radius in z and xy). `spots_df` needs z_um/y_um/x_um columns. label=True writes a unique integer id
    per spot (capped at uint16); otherwise every blob is 65535 (for viewing / Imaris Spots detection)."""
    out = np.zeros(tuple(int(s) for s in shape), np.uint16)
    if spots_df is None or len(spots_df) == 0:
        return out
    sp = tuple(float(s) for s in spacing)
    rz, ry, rx = (max(1, int(round(radius_um / s))) for s in sp)
    dz, dy, dx = np.ogrid[-rz:rz + 1, -ry:ry + 1, -rx:rx + 1]
    offs = np.argwhere(((dz / rz) ** 2 + (dy / ry) ** 2 + (dx / rx) ** 2) <= 1.0) - np.array([rz, ry, rx])
    zc = np.clip(np.round(spots_df["z_um"].to_numpy() / sp[0]).astype(int), 0, out.shape[0] - 1)
    yc = np.clip(np.round(spots_df["y_um"].to_numpy() / sp[1]).astype(int), 0, out.shape[1] - 1)
    xc = np.clip(np.round(spots_df["x_um"].to_numpy() / sp[2]).astype(int), 0, out.shape[2] - 1)
    for i, (z, y, x) in enumerate(zip(zc, yc, xc), start=1):
        zz = np.clip(z + offs[:, 0], 0, out.shape[0] - 1)
        yy = np.clip(y + offs[:, 1], 0, out.shape[1] - 1)
        xx = np.clip(x + offs[:, 2], 0, out.shape[2] - 1)
        out[zz, yy, xx] = (i & 0xFFFF) if label else 65535
    return out
