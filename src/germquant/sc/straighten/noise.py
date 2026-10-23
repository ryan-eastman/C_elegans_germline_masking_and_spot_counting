"""Per-image noise model and the generalised Anscombe transform (GAT).

The model ``var(I) = a + b * I`` (counts^2, I the raw count) is fitted on second differences along x
(variance 6 sigma^2 for white noise) inside the germline territory: a few planes, intensity bins up to a
quantile (ridges inflate the top), a robust MAD per bin, weighted least squares. The GAT
``y = (2 / b) * sqrt(max(b * I + 3/8 * b^2 + a, 0))`` then has unit noise, so ridge strengths can be
calibrated once against unit white noise. With no signal-dependent term (b about 0) the transform is the
plain scaling ``I / sqrt(a)``. The GAT is applied per nucleus crop, never to a whole channel.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np


@dataclass
class NoiseModel:
    a: float                       # intercept (counts^2)
    b: float                       # slope (counts^2 per count)
    rel_rms: float = 0.0           # relative RMS of the fit over the bins
    n_bins: int = 0
    table: list = field(default_factory=list)   # [(mean, var, n)] per bin
    fallback: bool = False

    def var(self, intensity) -> np.ndarray:
        return np.maximum(self.a + self.b * np.asarray(intensity, float), self.floor)

    @property
    def floor(self) -> float:
        return 1.0

    def sd(self, intensity) -> np.ndarray:
        return np.sqrt(self.var(intensity))

    def as_dict(self) -> dict:
        d = asdict(self)
        d["table"] = [[round(float(x), 4) for x in row] for row in self.table]
        return d


def _plane_sample(planes, masks, rng, max_n=3_000_000):
    ds, ms = [], []
    for p, mm in zip(planes, masks):
        p = np.asarray(p, dtype=np.float32)
        d = p[:, 2:] - 2.0 * p[:, 1:-1] + p[:, :-2]
        m = (p[:, 2:] + p[:, 1:-1] + p[:, :-2]) / 3.0
        if mm is not None:
            sel = mm[:, 1:-1] & mm[:, 2:] & mm[:, :-2]
            d, m = d[sel], m[sel]
        ds.append(d.ravel())
        ms.append(m.ravel())
    d = np.concatenate(ds) if ds else np.zeros(0, np.float32)
    m = np.concatenate(ms) if ms else np.zeros(0, np.float32)
    if d.size > max_n:
        idx = rng.choice(d.size, max_n, replace=False)
        idx.sort()
        d, m = d[idx], m[idx]
    return d, m


def fit_noise_model(planes, masks=None, *, bins: int = 30, qmax: float = 0.85, var_floor: float = 1.0,
                    seed: int = 0, min_per_bin: int = 400) -> NoiseModel:
    """Fit ``var = a + b I`` on a list of 2D planes (raw counts) restricted to the boolean ``masks``."""
    rng = np.random.default_rng(seed)
    masks = masks if masks is not None else [None] * len(planes)
    d, m = _plane_sample(planes, masks, rng)
    if d.size < 4 * min_per_bin:
        return NoiseModel(a=float(var_floor), b=0.0, rel_rms=float("nan"), fallback=True)
    qs = np.unique(np.quantile(m, np.linspace(0.0, qmax, bins + 1)))
    rows = []
    for lo, hi in zip(qs[:-1], qs[1:]):
        sel = (m >= lo) & (m < hi)
        n = int(sel.sum())
        if n < min_per_bin:
            continue
        dd = d[sel]
        mad = float(np.median(np.abs(dd - np.median(dd))))
        rows.append((float(m[sel].mean()), (1.4826 * mad) ** 2 / 6.0, n))
    if len(rows) < 2:
        v = float(np.var(d) / 6.0) if d.size else var_floor
        return NoiseModel(a=max(v, var_floor), b=0.0, rel_rms=float("nan"), n_bins=len(rows),
                          table=rows, fallback=True)
    t = np.asarray(rows, float)
    w = np.sqrt(t[:, 2])
    A = np.stack([np.ones(len(t)), t[:, 0]], 1)
    a, b = np.linalg.lstsq(A * w[:, None], t[:, 1] * w, rcond=None)[0]
    pred = a + b * t[:, 0]
    rel = float(np.sqrt(np.mean(((t[:, 1] - pred) / np.maximum(pred, 1e-9)) ** 2)))
    fallback = False
    if not np.isfinite(a) or not np.isfinite(b):
        a, b, fallback = float(np.median(t[:, 1])), 0.0, True
    if b < 0:
        # a variance that falls with intensity is not a camera; keep a constant model
        a, b, fallback = float(np.average(t[:, 1], weights=t[:, 2])), 0.0, True
    # the model must stay positive over the observed range
    lo_i = float(t[:, 0].min())
    if a + b * lo_i < var_floor:
        a = var_floor - b * lo_i
    return NoiseModel(a=float(a), b=float(b), rel_rms=rel, n_bins=len(rows), table=rows, fallback=fallback)


def gat(x, model: NoiseModel) -> np.ndarray:
    """Generalised Anscombe transform to unit noise (float32)."""
    x = np.asarray(x, dtype=np.float32)
    b, a = float(model.b), float(model.a)
    if b <= 1e-6:
        return (x / np.float32(np.sqrt(max(a, model.floor)))).astype(np.float32)
    return ((2.0 / b) * np.sqrt(np.maximum(b * x + 0.375 * b * b + a, 0.0))).astype(np.float32)


def nucleus_noise_sd(y_crop: np.ndarray, terr: np.ndarray) -> float:
    """Robust SD of the GAT second differences along x inside the territory (about 1 when the model fits)."""
    d = y_crop[:, :, 2:] - 2.0 * y_crop[:, :, 1:-1] + y_crop[:, :, :-2]
    sel = terr[:, :, 1:-1] & terr[:, :, 2:] & terr[:, :, :-2]
    v = d[sel]
    if v.size < 50:
        return float("nan")
    mad = float(np.median(np.abs(v - np.median(v))))
    return 1.4826 * mad / np.sqrt(6.0)
