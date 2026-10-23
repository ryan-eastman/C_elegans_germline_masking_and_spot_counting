"""Sampling a crop at whole-image micron positions (pipeline frame) and the normal-plane grids."""
from __future__ import annotations

import numpy as np
from scipy.ndimage import map_coordinates


class Sampler:
    """``vol`` is a crop whose voxel [0, 0, 0] is whole-image voxel ``org``; ``self(points_um)`` interpolates it
    (order 1 by default) at whole-image micron positions; outside the crop -> ``cval``."""

    def __init__(self, vol: np.ndarray, org, spacing, cval: float = np.nan):
        self.vol = np.asarray(vol, np.float32)
        self.org = np.asarray(org, float)
        self.sp = np.asarray(spacing, float)
        self.cval = cval

    def index(self, points_um: np.ndarray) -> np.ndarray:
        return np.asarray(points_um, float) / self.sp - self.org

    def __call__(self, points_um: np.ndarray, order: int = 1, cval: float | None = None) -> np.ndarray:
        P = np.asarray(points_um, float)
        shp = P.shape[:-1]
        idx = self.index(P.reshape(-1, 3)).T
        v = map_coordinates(self.vol, idx, order=order, mode="constant",
                            cval=self.cval if cval is None else cval, prefilter=order > 1)
        return v.reshape(shp)

    def nearest(self, points_um: np.ndarray) -> np.ndarray:
        """Nearest-voxel lookup (for labels and masks); -1 outside the crop."""
        P = np.asarray(points_um, float)
        shp = P.shape[:-1]
        idx = np.rint(self.index(P.reshape(-1, 3))).astype(np.int64)
        ok = np.all((idx >= 0) & (idx < np.asarray(self.vol.shape)), axis=1)
        out = np.full(len(idx), -1, np.int64)
        out[ok] = self.vol[tuple(idx[ok].T)].astype(np.int64)
        return out.reshape(shp)


class LabelLookup:
    """Nearest-voxel lookup on an integer crop (labels), -1 outside."""

    def __init__(self, labels: np.ndarray, org, spacing):
        self.lab = labels
        self.org = np.asarray(org, float)
        self.sp = np.asarray(spacing, float)

    def __call__(self, points_um: np.ndarray) -> np.ndarray:
        P = np.asarray(points_um, float)
        shp = P.shape[:-1]
        idx = np.rint(P.reshape(-1, 3) / self.sp - self.org).astype(np.int64)
        ok = np.all((idx >= 0) & (idx < np.asarray(self.lab.shape)), axis=1)
        out = np.full(len(idx), -1, np.int64)
        out[ok] = self.lab[tuple(idx[ok].T)].astype(np.int64)
        return out.reshape(shp)


def plane_points(C: np.ndarray, n_lat: np.ndarray, n_ax: np.ndarray, U: np.ndarray, V: np.ndarray) -> np.ndarray:
    """(m, nu, nv, 3) points C + u n_lat + v n_ax."""
    return (C[:, None, None, :] + U[None, :, None, None] * n_lat[:, None, None, :]
            + V[None, None, :, None] * n_ax[:, None, None, :])


def grid(half: float, step: float) -> np.ndarray:
    n = int(round(half / step))
    return np.arange(-n, n + 1) * step
