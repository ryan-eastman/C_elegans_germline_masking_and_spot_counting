"""Rotation-minimising frames with one constant roll toward the optical axis, and the strand <-> strip maps.

Frames are computed by double reflection (Wang et al. 2008, ACM TOG 27:2) on samples every ``ds`` um over
s from -``extend`` to L + ``extend``; never Frenet frames (they flip at inflections). One roll per strand,
``phi = atan2(sum w (n2 . z_hat_perp), sum w (n1 . z_hat_perp))`` with ``w = |z_perp|^2`` (z_perp: the optical
axis projected into the normal plane), turns the frame so that ``n_ax = cos(phi) n1 + sin(phi) n2`` points along
the optical axis as far as one constant roll can; ``n_lat = t x n_ax``. A strand along z (sum w about 0) keeps
n_ax = n1. ``to_local`` / ``from_local`` map whole-image microns to (s, u, v, w) and back, where u runs along
n_lat, v along n_ax and w along the tangent (zero at the closest point).
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .refine import Centreline, unit

ZHAT = np.array([1.0, 0.0, 0.0])


def initial_normal(t: np.ndarray) -> np.ndarray:
    v = ZHAT - t * float(t @ ZHAT)
    if np.linalg.norm(v) < 0.15:
        yhat = np.array([0.0, 1.0, 0.0])
        v = yhat - t * float(t @ yhat)
    return v / np.linalg.norm(v)


def rmf(C: np.ndarray, T: np.ndarray, n0: np.ndarray | None = None) -> np.ndarray:
    """Double-reflection rotation-minimising normals along the samples (unit, perpendicular to T)."""
    m = len(C)
    R = np.empty((m, 3))
    R[0] = initial_normal(T[0]) if n0 is None else unit(n0 - T[0] * float(n0 @ T[0]))
    for i in range(m - 1):
        v1 = C[i + 1] - C[i]
        c1 = float(v1 @ v1)
        if c1 < 1e-24:
            R[i + 1] = R[i]
            continue
        rL = R[i] - (2.0 / c1) * float(v1 @ R[i]) * v1
        tL = T[i] - (2.0 / c1) * float(v1 @ T[i]) * v1
        v2 = T[i + 1] - tL
        c2 = float(v2 @ v2)
        r = rL - (2.0 / c2) * float(v2 @ rL) * v2 if c2 > 1e-24 else rL
        r = r - T[i + 1] * float(r @ T[i + 1])
        R[i + 1] = r / np.linalg.norm(r)
    return R


def roll_angle(n1: np.ndarray, T: np.ndarray) -> float:
    n2 = np.cross(T, n1)
    zp = ZHAT[None] - T * (T @ ZHAT)[:, None]
    w = (zp * zp).sum(1)
    zu = zp / np.maximum(np.sqrt(w), 1e-12)[:, None]
    a = float((w * (n2 * zu).sum(1)).sum())
    b = float((w * (n1 * zu).sum(1)).sum())
    if w.sum() < 1e-9 or (abs(a) < 1e-12 and abs(b) < 1e-12):
        return 0.0
    return float(np.arctan2(a, b))


class StrandFrame:
    """Samples, tangents and rolled rotation-minimising frames of one measured centreline."""

    def __init__(self, line: Centreline, ds: float = 0.05, extend: float = 0.5):
        self.line = line
        self.ds = float(ds)
        self.extend = float(extend)
        n = int(np.floor((line.L + 2 * extend) / ds + 1e-9)) + 1
        self.s = -extend + np.arange(n) * ds
        self.C = line.point(self.s)
        self.T = line.tangent(self.s)
        self.n1 = rmf(self.C, self.T)
        self.phi = roll_angle(self.n1, self.T)
        self.n_ax, self.n_lat = self._roll(self.n1, self.T)
        self._tree = None

    def _roll(self, n1, T):
        n2 = np.cross(T, n1)
        n_ax = np.cos(self.phi) * n1 + np.sin(self.phi) * n2
        n_ax = unit(n_ax - T * (n_ax * T).sum(1)[:, None])
        n_lat = unit(np.cross(T, n_ax))
        return n_ax, n_lat

    @property
    def in_strand(self) -> np.ndarray:
        return (self.s >= -1e-9) & (self.s <= self.line.L + 1e-9)

    def roll_from_axial_deg(self) -> np.ndarray:
        zp = ZHAT[None] - self.T * (self.T @ ZHAT)[:, None]
        nz = np.linalg.norm(zp, axis=1)
        zu = zp / np.maximum(nz, 1e-12)[:, None]
        ang = np.degrees(np.arctan2((np.cross(zu, self.n_ax) * self.T).sum(1), (zu * self.n_ax).sum(1)))
        ang[nz < 0.15] = np.nan
        return ang

    def frame_at(self, s) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """(C, T, n_lat, n_ax) at arbitrary s (normals interpolated between samples, re-orthogonalised)."""
        s = np.atleast_1d(np.asarray(s, float))
        C = self.line.point(s)
        T = self.line.tangent(s)
        f = np.clip((s - self.s[0]) / self.ds, 0, len(self.s) - 1)
        i0 = np.clip(np.floor(f).astype(int), 0, len(self.s) - 2) if len(self.s) > 1 else np.zeros(len(s), int)
        a = (f - i0)[:, None]
        if len(self.s) > 1:
            n1 = (1 - a) * self.n1[i0] + a * self.n1[i0 + 1]
        else:
            n1 = np.repeat(self.n1[:1], len(s), 0)
        n1 = unit(n1 - T * (n1 * T).sum(1)[:, None])
        n_ax, n_lat = self._roll(n1, T)
        return C, T, n_lat, n_ax

    def from_local(self, s, u, v, w=None) -> np.ndarray:
        C, T, n_lat, n_ax = self.frame_at(s)
        u = np.atleast_1d(np.asarray(u, float))[:, None]
        v = np.atleast_1d(np.asarray(v, float))[:, None]
        out = C + u * n_lat + v * n_ax
        if w is not None:
            out = out + np.atleast_1d(np.asarray(w, float))[:, None] * T
        return out

    def project(self, points: np.ndarray, s0: np.ndarray | None = None, iters: int = 12) -> np.ndarray:
        """s of the closest point on the (extended) centreline, by Newton from the nearest sample."""
        P = np.atleast_2d(np.asarray(points, float))
        if s0 is None:
            if self._tree is None:
                self._tree = cKDTree(self.C)
            _, k = self._tree.query(P)
            s = self.s[k].astype(float)
        else:
            s = np.asarray(s0, float).copy()
        lo, hi = self.s[0], self.s[-1]
        for _ in range(iters):
            C = self.line.point(s)
            T = self.line.tangent(s)
            K = self.line.curvature_vector(s)
            r = C - P
            g = (r * T).sum(1)
            gp = 1.0 + (r * K).sum(1)
            gp = np.where(gp > 0.2, gp, 1.0)
            step = g / gp
            s = np.clip(s - step, lo, hi)
            if np.all(np.abs(step) < 1e-12):
                break
        return s

    def to_local(self, points: np.ndarray, s0=None) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        P = np.atleast_2d(np.asarray(points, float))
        s = self.project(P, s0)
        C, T, n_lat, n_ax = self.frame_at(s)
        r = P - C
        return s, (r * n_lat).sum(1), (r * n_ax).sum(1), (r * T).sum(1)


def to_strip(frames: list[StrandFrame], points: np.ndarray, search_um: float = 1.5):
    """For each point the closest strand: (strand index, s, u, v, w, dist, second index, second dist).
    Candidates are the strands with a sample within ``search_um`` (else the strand of the nearest sample).
    -1 / NaN when ``frames`` is empty."""
    P = np.atleast_2d(np.asarray(points, float))
    n = len(P)
    out = {"strand": np.full(n, -1), "s": np.full(n, np.nan), "u": np.full(n, np.nan), "v": np.full(n, np.nan),
           "w": np.full(n, np.nan), "dist": np.full(n, np.inf), "second": np.full(n, -1),
           "second_dist": np.full(n, np.inf), "second_u": np.full(n, np.nan), "second_v": np.full(n, np.nan)}
    if not frames or n == 0:
        return out
    allC = np.vstack([f.C for f in frames])
    owner = np.concatenate([np.full(len(f.C), j) for j, f in enumerate(frames)])
    tree = cKDTree(allC)
    for i in range(n):
        near = tree.query_ball_point(P[i], r=search_um)
        if not near:
            _, j0 = tree.query(P[i])
            near = [int(j0)]
        cands = sorted({int(owner[j]) for j in near})
        res = []
        for j in cands:
            s, u, v, w = frames[j].to_local(P[i:i + 1])
            d = float(np.sqrt(u[0] ** 2 + v[0] ** 2 + w[0] ** 2))
            res.append((d, j, float(s[0]), float(u[0]), float(v[0]), float(w[0])))
        res.sort(key=lambda t: (t[0], t[1]))
        d, j, s, u, v, w = res[0]
        out["strand"][i], out["s"][i], out["u"][i], out["v"][i], out["w"][i], out["dist"][i] = j, s, u, v, w, d
        if len(res) > 1:
            out["second"][i] = res[1][1]
            out["second_dist"][i] = res[1][0]
            out["second_u"][i] = res[1][3]
            out["second_v"][i] = res[1][4]
    return out
