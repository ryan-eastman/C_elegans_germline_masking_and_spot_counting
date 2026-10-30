"""Hessian line fields on the GAT crop, with a strength calibrated against a unit-noise null.

Gaussian derivatives at sigma = ``fields.sigma_um`` / spacing (about one voxel per axis at stride 1, where the
PSF is nearly isotropic in voxel units). On territory voxels only: eigenvalues w0 <= w1 <= w2 and eigenvectors;
strength ``S = max(-w1, 0) * mean(sigma^2) / null_q`` (null units: about a fraction 1 - q of unit white-noise
crest voxels pass S >= 1); the Newton offset ``p = -sum_k (g . v_k / w_k) v_k`` over the two cross-sectional
directions; lineness ``|w2| / |w1|`` (0 for an ideal line, about 1 for a blob or a junction core); a crest voxel
has w1 < 0 and ``|p|_inf <= crest_offset_max`` voxel. Directions and offsets are in voxel units.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.ndimage import gaussian_filter


@dataclass
class Fields:
    S: np.ndarray            # (Z, Y, X) float32, null units, 0 outside the territory
    off: np.ndarray | None   # (Z, Y, X, 3) float32 Newton offsets (voxels)
    D: np.ndarray | None     # (Z, Y, X, 3) float32 line direction (voxel units, unit length)
    LIN: np.ndarray          # (Z, Y, X) float32 lineness (1 outside the territory)
    crest: np.ndarray        # (Z, Y, X) bool


def sigma_vox(sigma_um, spacing) -> tuple[float, float, float]:
    return tuple(float(s) / float(d) for s, d in zip(sigma_um, spacing))


def _hessian(y: np.ndarray, sig) -> tuple[np.ndarray, list[np.ndarray]]:
    g = [gaussian_filter(y, sig, order=o, mode="nearest") for o in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    h = {}
    for name, o in (("zz", (2, 0, 0)), ("yy", (0, 2, 0)), ("xx", (0, 0, 2)), ("zy", (1, 1, 0)),
                    ("zx", (1, 0, 1)), ("yx", (0, 1, 1))):
        h[name] = gaussian_filter(y, sig, order=o, mode="nearest")
    return h, g


def line_fields(y: np.ndarray, terr: np.ndarray, sig, null_q: float | None, *, max_off: float = 0.6,
                with_vectors: bool = True) -> Fields:
    y = np.asarray(y, np.float32)
    sig = tuple(float(s) for s in sig)
    h, g = _hessian(y, sig)
    idx = np.nonzero(terr)
    n = idx[0].size
    shp = y.shape
    S = np.zeros(shp, np.float32)
    LIN = np.ones(shp, np.float32)
    crest = np.zeros(shp, bool)
    off = np.zeros(shp + (3,), np.float32) if with_vectors else None
    D = np.zeros(shp + (3,), np.float32) if with_vectors else None
    if n == 0:
        return Fields(S, off, D, LIN, crest)
    H = np.empty((n, 3, 3), np.float64)
    H[:, 0, 0] = h["zz"][idx]
    H[:, 1, 1] = h["yy"][idx]
    H[:, 2, 2] = h["xx"][idx]
    H[:, 0, 1] = H[:, 1, 0] = h["zy"][idx]
    H[:, 0, 2] = H[:, 2, 0] = h["zx"][idx]
    H[:, 1, 2] = H[:, 2, 1] = h["yx"][idx]
    w, v = np.linalg.eigh(H)
    norm = float(np.mean(np.square(sig)))
    s_raw = np.clip(-w[:, 1], 0, None) * norm
    Sv = s_raw / null_q if null_q else s_raw
    gv = np.stack([g[0][idx], g[1][idx], g[2][idx]], -1).astype(np.float64)
    p = np.zeros((n, 3), np.float64)
    for k in (0, 1):
        den = np.where(w[:, k] < -1e-9, w[:, k], -1e-9)
        p -= ((gv * v[:, :, k]).sum(1) / den)[:, None] * v[:, :, k]
    lin = np.abs(w[:, 2]) / np.maximum(np.abs(w[:, 1]), 1e-12)
    S[idx] = Sv
    LIN[idx] = lin
    ok = (w[:, 1] < 0) & (np.abs(p).max(1) <= max_off)
    crest[tuple(a[ok] for a in idx)] = True
    if with_vectors:
        off[idx] = p
        D[idx] = v[:, :, 2]
    return Fields(S, off, D, LIN, crest)


@lru_cache(maxsize=16)
def _null_cached(sig: tuple, q: float, shape: tuple, seed: int, max_off: float) -> float:
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal(shape).astype(np.float32)
    terr = np.zeros(shape, bool)
    terr[4:-4, 8:-8, 8:-8] = True
    f = line_fields(noise, terr, sig, None, max_off=max_off, with_vectors=False)
    sel = f.crest & terr & (f.S > 0)
    return float(np.quantile(f.S[sel], q))


def calibrate_null(sig, q: float = 0.999, shape=(40, 96, 96), seed: int = 0, max_off: float = 0.6) -> float:
    """The q quantile of the raw strength on crest voxels of unit white noise (cached per arguments)."""
    key = tuple(round(float(s), 6) for s in sig)
    return _null_cached(key, float(q), tuple(int(v) for v in shape), int(seed), float(max_off))
