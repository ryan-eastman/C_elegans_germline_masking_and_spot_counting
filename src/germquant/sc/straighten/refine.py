"""Centreline refinement, arc length, the PSF estimate and the anchored half-maximum end rule.

* Re-centring: the traced polyline is resampled every ``refine.ds_um``, smoothed by a first penalised spline,
  then each sample moves in its normal plane by mean-shift on raw counts minus the territory background
  (Gaussian window ``window_um`` (axial, lateral) over a grid ``grid_um``, total shift capped at
  ``max_shift_um``); samples near another strand (contaminated) and in junction cores do not move.
* Final fit: a cubic B-spline with knots every ``knot_um``, second-difference penalty ``lambda * mean(W)``,
  weights ``W = 1 / sigma_loc^2`` per axis (z, y, x); contaminated and junction samples weigh 1e-3 of that.
  Arc length by 5-node Gauss-Legendre quadrature on sub-intervals of at most 0.02 um.
* End rule: at a ``free`` end, an anchor min(``anchor_um``, L / 3) inside the end; the profile p(s) of the core
  mean runs from -``inside_um`` (along the spline) to +``outside_um`` (straight along the anchor tangent, each
  sample re-centred by at most ``max_recentre_um``); ``p = b + a * Phi((e - s) / sigma)`` is fitted with sigma
  bounded to ``sigma_bounds`` x sigma_t, ``sigma_t = sqrt(sigma_xy^2 (1 - t_z^2) + sigma_z^2 t_z^2)``. The end
  is ``free_fit`` when a >= ``fit_min_amp_se`` x SE(a) and the residual RMS <= ``fit_max_rel_rms`` x a, else
  ``free_unfit`` (the traced end is kept). Non-free ends stay where the trace stopped.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import BSpline
from scipy.special import ndtr

from .account import resample
from .sampling import Sampler, grid, plane_points

_GL_X, _GL_W = np.polynomial.legendre.leggauss(5)


# ------------------------------------------------------------------------------------ P-spline
def pspline(s: np.ndarray, P: np.ndarray, W: np.ndarray, knot_um: float, lam_eff: float, k: int = 3) -> BSpline:
    """Penalised cubic B-spline per axis (second-difference penalty ``lam_eff``); ``W`` (m, 3) weights."""
    s = np.asarray(s, float)
    L = float(s[-1] - s[0])
    nint = max(1, int(np.ceil(L / knot_um - 1e-9)))
    inner = np.linspace(s[0], s[-1], nint + 1)
    t = np.r_[[s[0]] * k, inner, [s[-1]] * k]
    B = BSpline.design_matrix(s, t, k).toarray()
    nb = B.shape[1]
    D = np.diff(np.eye(nb), n=2, axis=0) if nb > 2 else np.zeros((0, nb))
    DtD = D.T @ D
    C = np.empty((nb, 3))
    for a in range(3):
        w = W[:, a]
        A = B.T @ (w[:, None] * B) + lam_eff * DtD + 1e-10 * np.eye(nb)
        C[:, a] = np.linalg.solve(A, B.T @ (w * P[:, a]))
    return BSpline(t, C, k)


def unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.maximum(n, 1e-12)


def axial_frames(T: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """n_ax = the optical axis projected into the normal plane (y axis when the tangent is near z),
    n_lat = T x n_ax. Used for re-centring and profiles; strips use rotation-minimising frames."""
    zhat = np.array([1.0, 0.0, 0.0])
    yhat = np.array([0.0, 1.0, 0.0])
    v = zhat[None] - T * (T @ zhat)[:, None]
    bad = np.linalg.norm(v, axis=1) < 0.15
    if bad.any():
        v[bad] = yhat[None] - T[bad] * (T[bad] @ yhat)[:, None]
    n_ax = unit(v)
    n_lat = unit(np.cross(T, n_ax))
    return n_lat, n_ax


# ------------------------------------------------------------------------------------ centreline
class Centreline:
    """Arc-length parametrised centreline: a straight start extension of ``ext0`` um ending at P(qa), the spline
    P(q) for q in [qa, qb], a straight end extension of ``ext1`` um. ``s`` runs over [0, L]; beyond either end
    the end tangent continues straight (used for the 0.5 um strip margins)."""

    def __init__(self, spl: BSpline, qa: float, qb: float, ext0: float = 0.0, ext1: float = 0.0):
        self.spl = spl
        self.d1 = spl.derivative(1)
        self.d2 = spl.derivative(2)
        self.qa, self.qb = float(qa), float(qb)
        self.ext0, self.ext1 = max(float(ext0), 0.0), max(float(ext1), 0.0)
        n = max(4, int(np.ceil((self.qb - self.qa) / 0.02)))
        self._qtab = np.linspace(self.qa, self.qb, n + 1)
        seg = self._gl(self._qtab[:-1], self._qtab[1:])
        self._stab = np.concatenate([[0.0], np.cumsum(seg)])
        self.Ls = float(self._stab[-1])
        self.L = self.ext0 + self.Ls + self.ext1
        self.Pa = spl(self.qa)
        self.Pb = spl(self.qb)
        self.Ta = unit(self.d1(self.qa))
        self.Tb = unit(self.d1(self.qb))

    def _speed(self, q):
        return np.linalg.norm(self.d1(q), axis=-1)

    def _gl(self, a, b):
        a = np.asarray(a, float)
        b = np.asarray(b, float)
        mid, half = 0.5 * (a + b), 0.5 * (b - a)
        x = mid[..., None] + half[..., None] * _GL_X
        return (self._speed(x) * _GL_W).sum(-1) * half

    def spline_arc(self, q) -> np.ndarray:
        """Arc length from qa to q along the spline (vectorised)."""
        q = np.clip(np.asarray(q, float), self.qa, self.qb)
        k = np.clip(np.searchsorted(self._qtab, q, side="right") - 1, 0, len(self._qtab) - 2)
        return self._stab[k] + self._gl(self._qtab[k], q)

    def q_of(self, sigma) -> np.ndarray:
        """Spline parameter at arc length ``sigma`` from qa (clipped to [0, Ls])."""
        sig = np.clip(np.atleast_1d(np.asarray(sigma, float)), 0.0, self.Ls)
        k = np.clip(np.searchsorted(self._stab, sig, side="right") - 1, 0, len(self._stab) - 2)
        s0, s1 = self._stab[k], self._stab[k + 1]
        q0, q1 = self._qtab[k], self._qtab[k + 1]
        q = q0 + (sig - s0) / np.maximum(s1 - s0, 1e-15) * (q1 - q0)
        for _ in range(4):
            f = self._stab[k] + self._gl(q0, q) - sig
            q = np.clip(q - f / np.maximum(self._speed(q), 1e-12), q0, q1)
        return q

    def point(self, s) -> np.ndarray:
        s = np.atleast_1d(np.asarray(s, float))
        out = np.empty((len(s), 3))
        a = s < self.ext0
        b = s > self.ext0 + self.Ls
        m = ~(a | b)
        if a.any():
            out[a] = self.Pa[None] + (s[a] - self.ext0)[:, None] * self.Ta[None]
        if b.any():
            out[b] = self.Pb[None] + (s[b] - self.ext0 - self.Ls)[:, None] * self.Tb[None]
        if m.any():
            out[m] = self.spl(self.q_of(s[m] - self.ext0))
        return out

    def tangent(self, s) -> np.ndarray:
        s = np.atleast_1d(np.asarray(s, float))
        out = np.empty((len(s), 3))
        a = s < self.ext0
        b = s > self.ext0 + self.Ls
        m = ~(a | b)
        out[a] = self.Ta
        out[b] = self.Tb
        if m.any():
            out[m] = unit(self.d1(self.q_of(s[m] - self.ext0)))
        return out

    def curvature_vector(self, s) -> np.ndarray:
        """dT/ds (1/um); zero on the straight extensions."""
        s = np.atleast_1d(np.asarray(s, float))
        out = np.zeros((len(s), 3))
        m = (s >= self.ext0) & (s <= self.ext0 + self.Ls)
        if m.any():
            q = self.q_of(s[m] - self.ext0)
            d1 = self.d1(q)
            d2 = self.d2(q)
            sp2 = np.maximum((d1 * d1).sum(1), 1e-24)
            T = d1 / np.sqrt(sp2)[:, None]
            out[m] = (d2 - (d2 * T).sum(1)[:, None] * T) / sp2[:, None]
        return out

    def reversed(self) -> "Centreline":
        t, c, k = self.spl.t, self.spl.c, self.spl.k
        T0, T1 = float(t[k]), float(t[len(t) - k - 1])
        cc = T0 + T1
        spl_r = BSpline(cc - t[::-1], c[::-1], k)
        return Centreline(spl_r, cc - self.qb, cc - self.qa, self.ext1, self.ext0)

    def to_state(self) -> dict:
        return {"t": self.spl.t, "c": self.spl.c, "k": int(self.spl.k), "qa": self.qa, "qb": self.qb,
                "ext0": self.ext0, "ext1": self.ext1}

    @classmethod
    def from_state(cls, st: dict) -> "Centreline":
        return cls(BSpline(st["t"], st["c"], st["k"]), st["qa"], st["qb"], st["ext0"], st["ext1"])


# ------------------------------------------------------------------------------------ re-centring
def mean_shift(P: np.ndarray, T: np.ndarray, sampler_bg: Sampler, *, window=(0.30, 0.15), grid_um=(0.6, 0.4),
               step=0.05, iterations=4, max_shift=0.3, frozen: np.ndarray | None = None) -> np.ndarray:
    """Normal-plane mean-shift of the samples ``P`` (tangents ``T``) on clipped background-subtracted counts."""
    n_lat, n_ax = axial_frames(T)
    U = grid(grid_um[1], step)
    V = grid(grid_um[0], step)
    Wn = np.exp(-0.5 * ((U[:, None] / window[1]) ** 2 + (V[None, :] / window[0]) ** 2))
    P0 = P.copy()
    Pr = P.copy()
    move = np.ones(len(P), bool) if frozen is None else ~frozen
    for _ in range(int(iterations)):
        vals = sampler_bg(plane_points(Pr, n_lat, n_ax, U, V), cval=0.0)
        vals = np.clip(np.nan_to_num(vals, nan=0.0), 0, None) * Wn[None]
        tot = vals.sum(axis=(1, 2))
        ok = (tot > 1e-9) & move
        du = np.where(ok, (vals * U[None, :, None]).sum(axis=(1, 2)) / np.maximum(tot, 1e-9), 0.0)
        dv = np.where(ok, (vals * V[None, None, :]).sum(axis=(1, 2)) / np.maximum(tot, 1e-9), 0.0)
        Pr = Pr + du[:, None] * n_lat + dv[:, None] * n_ax
        sh = Pr - P0
        nrm = np.linalg.norm(sh, axis=1)
        big = nrm > max_shift
        if big.any():
            Pr[big] = P0[big] + sh[big] * (max_shift / nrm[big])[:, None]
    return Pr


@dataclass
class Refined:
    line: Centreline
    fit_rms_um: float
    shifts_um: np.ndarray = field(default_factory=lambda: np.zeros(0))
    n_frozen: int = 0


def refine_points(points: np.ndarray, sampler_bg: Sampler, p, *, frozen_fn=None) -> Refined | None:
    """Re-centre and fit the final penalised spline. ``frozen_fn(P)`` -> bool per sample (contaminated or junction
    core: no shift, low weight)."""
    ds = float(p["refine.ds_um"])
    P0, _, q = resample(np.asarray(points, float), ds)
    if len(P0) < 4:
        P0, _, q = resample(np.asarray(points, float), max(q[-1] / 4.0, 1e-3)) if len(P0) >= 2 else (P0, None, q)
        if len(P0) < 4:
            return None
    knot = float(p["refine.knot_um"])
    lam = float(p["refine.lambda"])
    sig_loc = np.asarray(p["refine.sigma_loc_um"], float)
    W0 = np.ones_like(P0)
    spl0 = pspline(q, P0, W0, knot, lam)
    P1 = spl0(q)
    T1 = unit(spl0.derivative()(q))
    frozen = frozen_fn(P1) if frozen_fn is not None else np.zeros(len(P1), bool)
    Pr = mean_shift(P1, T1, sampler_bg, window=tuple(p["refine.window_um"]), grid_um=tuple(p["refine.grid_um"]),
                    iterations=int(p["refine.iterations"]), max_shift=float(p["refine.max_shift_um"]),
                    frozen=frozen)
    q2 = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(Pr, axis=0), axis=1))])
    keep = np.concatenate([[True], np.diff(q2) > 1e-6])
    Pr, q2, frozen_k = Pr[keep], q2[keep], frozen[keep]
    if len(Pr) < 4:
        return None
    w_nom = 1.0 / sig_loc ** 2
    W = np.tile(w_nom, (len(Pr), 1))
    W[frozen_k] *= 1e-3
    spl = pspline(q2, Pr, W, knot, lam * float(np.mean(w_nom)))
    res = spl(q2) - Pr
    clean = ~frozen_k
    rms = float(np.sqrt(np.mean((res[clean] ** 2).sum(1)))) if clean.any() else float(np.sqrt(np.mean((res ** 2).sum(1))))
    line = Centreline(spl, float(q2[0]), float(q2[-1]))
    return Refined(line=line, fit_rms_um=rms, shifts_um=np.linalg.norm(Pr - P1[keep], axis=1),
                   n_frozen=int(frozen_k.sum()))


def refine_imported(points: np.ndarray, p) -> Refined | None:
    """A drawn path measured as drawn: the final spline on the given points, no re-centring."""
    ds = float(p["refine.ds_um"])
    P0, _, q = resample(np.asarray(points, float), ds)
    if len(P0) < 4:
        return None
    sig_loc = np.asarray(p["refine.sigma_loc_um"], float)
    w_nom = 1.0 / sig_loc ** 2
    W = np.tile(w_nom, (len(P0), 1))
    spl = pspline(q, P0, W, float(p["refine.knot_um"]), float(p["refine.lambda"]) * float(np.mean(w_nom)))
    res = spl(q) - P0
    return Refined(line=Centreline(spl, float(q[0]), float(q[-1])),
                   fit_rms_um=float(np.sqrt(np.mean((res ** 2).sum(1)))))


# ------------------------------------------------------------------------------------ core values
def ellipse_mask(U: np.ndarray, V: np.ndarray, r_lat: float, r_ax: float) -> np.ndarray:
    return (U[:, None] / r_lat) ** 2 + (V[None, :] / r_ax) ** 2 <= 1.0 + 1e-9


def n_eff(r_lat: float, r_ax: float, spacing, n_ax_z: np.ndarray) -> np.ndarray:
    """Independent voxels in a cross-section ellipse (an approximation for the core standard error)."""
    dz, dxy = float(spacing[0]), float(spacing[1])
    e_ax = np.abs(n_ax_z) * dz + (1.0 - np.abs(n_ax_z)) * dxy
    return np.maximum(np.pi * r_lat * r_ax / (dxy * e_ax), 1.0)


def core_profile(sampler_raw: Sampler, C, n_lat, n_ax, bg: float, noise, spacing, r_lat=0.12, r_ax=0.25,
                 step=0.04):
    """Core mean of (I - bg) in the ellipse r_lat x r_ax at each sample, its SE from the noise model and the SNR."""
    U = grid(r_lat, step)
    V = grid(r_ax, 0.05)
    ell = ellipse_mask(U, V, r_lat, r_ax)
    vals = sampler_raw(plane_points(C, n_lat, n_ax, U, V))
    w = ell[None].astype(float) * np.isfinite(vals)
    vv = np.nan_to_num(vals, nan=0.0)
    cnt = np.maximum(w.sum(axis=(1, 2)), 1e-9)
    mean_abs = (vv * w).sum(axis=(1, 2)) / cnt
    core = mean_abs - bg
    se = noise.sd(mean_abs) / np.sqrt(n_eff(r_lat, r_ax, spacing, n_ax[:, 0]))
    snr = core / np.maximum(se, 1e-9)
    core[cnt < 0.5] = np.nan
    return core, se, snr


# ------------------------------------------------------------------------------------ end rule
@dataclass
class EndFit:
    end_class: str            # free_fit | free_unfit | contact | territory_edge | stack_edge | cycle_cut | max_length | ...
    e_um: float = 0.0         # edge position along the outward direction, from the anchor (um)
    e_se_um: float = 0.0
    anchor_arc_um: float = 0.0   # arc length of the anchor from this end (um, along the refined spline)
    sigma_um: float = float("nan")
    amp: float = float("nan")
    amp_se: float = float("nan")
    rel_rms: float = float("nan")
    profile_s: list = field(default_factory=list)
    profile: list = field(default_factory=list)


def fit_erf(sgrid: np.ndarray, prof: np.ndarray, sig_lo: float, sig_hi: float):
    """Least squares fit of p = b + a * Phi((e - s) / sigma), sigma in [sig_lo, sig_hi], a > 0.
    Returns (e, sigma, b, a, se_e, se_a, rms) or None."""
    from scipy.optimize import least_squares

    ok = np.isfinite(prof)
    s = sgrid[ok]
    y = prof[ok]
    if len(s) < 8 or not np.isfinite(y).all():
        return None
    e_grid = np.arange(s[0] + 0.1, s[-1] - 0.1 + 1e-9, 0.01)
    sig_grid = np.linspace(sig_lo, sig_hi, 9)
    E, SG = np.meshgrid(e_grid, sig_grid, indexing="ij")
    phi = ndtr((E.ravel()[:, None] - s[None, :]) / SG.ravel()[:, None])        # (G, n)
    n = len(s)
    sx = phi.sum(1)
    sxx = (phi * phi).sum(1)
    sy = y.sum()
    sxy = phi @ y
    det = n * sxx - sx * sx
    a = np.where(np.abs(det) > 1e-12, (n * sxy - sx * sy) / np.where(np.abs(det) > 1e-12, det, 1.0), 0.0)
    b = (sy - a * sx) / n
    sse = (y * y).sum() - 2 * b * sy - 2 * a * sxy + n * b * b + 2 * a * b * sx + a * a * sxx
    sse = np.where(a > 0, sse, np.inf)
    g = int(np.argmin(sse))
    if not np.isfinite(sse[g]):
        return None
    x0 = np.array([E.ravel()[g], SG.ravel()[g], b[g], max(a[g], 1e-6)])
    lo = np.array([s[0], sig_lo, -np.inf, 0.0])
    hi = np.array([s[-1], sig_hi, np.inf, np.inf])
    x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)

    def resid(x):
        return x[2] + x[3] * ndtr((x[0] - s) / x[1]) - y

    try:
        r = least_squares(resid, x0, bounds=(lo, hi), method="trf", x_scale="jac", max_nfev=200)
        x = r.x
        J = r.jac
        res = r.fun
    except Exception:  # noqa: BLE001 - fall back to the grid optimum
        x = x0
        res = resid(x)
        J = None
    dof = max(len(s) - 4, 1)
    s2 = float(res @ res) / dof
    se = np.full(4, np.nan)
    if J is not None:
        try:
            cov = np.linalg.pinv(J.T @ J) * s2
            se = np.sqrt(np.clip(np.diag(cov), 0, None))
        except Exception:  # noqa: BLE001
            pass
    rms = float(np.sqrt(np.mean(res ** 2)))
    return float(x[0]), float(x[1]), float(x[2]), float(x[3]), float(se[0]), float(se[3]), rms


def end_profile(line: Centreline, at_start: bool, sampler_bg: Sampler, p, *, anchor_um: float):
    """The core profile around one end: s <= 0 inside along the spline from the anchor, s > 0 straight outward,
    re-centred by mean-shift (capped). Returns (sgrid, profile, anchor_s, outward tangent at the anchor)."""
    inside = float(p["ends.inside_um"])
    outside = float(p["ends.outside_um"])
    ds = float(p["ends.step_um"])
    rmax = float(p["ends.max_recentre_um"])
    core_ax, core_lat = (float(v) for v in p["profile.core_um"])
    sg = np.arange(-round(inside / ds), round(outside / ds) + 1) * ds
    L = line.L
    sa = anchor_um if at_start else L - anchor_um
    ta = line.tangent(sa)[0]
    out_dir = -ta if at_start else ta
    pa = line.point(sa)[0]
    # s <= 0 runs inward along the spline from the anchor (larger s at the start end, smaller at the far end)
    s_in = np.clip(sa - sg if at_start else sa + sg, 0.0, L)
    P = np.where((sg <= 0)[:, None], line.point(s_in), pa[None] + sg[:, None] * out_dir[None])
    T = np.where((sg <= 0)[:, None], line.tangent(s_in), ta[None])
    if rmax > 0:
        out_side = sg > 0
        Pm = mean_shift(P, T, sampler_bg, window=(0.30, 0.15), grid_um=(0.5, 0.3), iterations=3, max_shift=rmax)
        P = np.where(out_side[:, None], Pm, P)
    n_lat, n_ax = axial_frames(T)
    U = grid(core_lat, 0.04)
    V = grid(core_ax, 0.05)
    ell = ellipse_mask(U, V, core_lat, core_ax)
    vals = sampler_bg(plane_points(P, n_lat, n_ax, U, V))
    w = ell[None].astype(float) * np.isfinite(vals)
    prof = (np.nan_to_num(vals) * w).sum(axis=(1, 2)) / np.maximum(w.sum(axis=(1, 2)), 1e-9)
    prof[w.sum(axis=(1, 2)) < 0.5] = np.nan
    return sg, prof, sa, out_dir, ta


def sigma_t(psf_sigma, t_z: float) -> float:
    sz, sxy = float(psf_sigma[0]), float(psf_sigma[1])
    tz2 = min(float(t_z) ** 2, 1.0)
    return float(np.sqrt(sxy ** 2 * (1 - tz2) + sz ** 2 * tz2))


def fit_end(line: Centreline, at_start: bool, end_type: str, sampler_bg: Sampler, psf_sigma, p) -> EndFit:
    if end_type != "free":
        return EndFit(end_class=end_type)
    anchor = min(float(p["ends.anchor_um"]), line.L / 3.0)
    sg, prof, sa, out_dir, ta = end_profile(line, at_start, sampler_bg, p, anchor_um=anchor)
    st = sigma_t(psf_sigma, ta[0])
    lo_f, hi_f = (float(v) for v in p["ends.sigma_bounds"])
    r = fit_erf(sg, prof, lo_f * st, hi_f * st)
    ef = EndFit(end_class="free_unfit", anchor_arc_um=anchor, profile_s=[float(v) for v in sg],
                profile=[float(v) if np.isfinite(v) else None for v in prof])
    if r is None:
        return ef
    e, sig, b0, a, se_e, se_a, rms = r
    ef.sigma_um, ef.amp, ef.amp_se = sig, a, se_a
    ef.rel_rms = rms / a if a > 0 else float("inf")
    good = (a > 0 and np.isfinite(se_a) and a >= float(p["ends.fit_min_amp_se"]) * se_a
            and rms <= float(p["ends.fit_max_rel_rms"]) * a)
    if good:
        ef.end_class = "free_fit"
        ef.e_um = e
        ef.e_se_um = se_e if np.isfinite(se_e) else 0.0
    return ef


def finalize_line(line: Centreline, fits: tuple[EndFit, EndFit]) -> tuple[Centreline, float, float]:
    """The measured centreline from the refined one and the two end fits. Returns (line, e0, e1) where e_k is
    the net change of length at end k against the traced end (positive: the end lies beyond it)."""
    L = line.L
    f0, f1 = fits
    # arc positions (on the refined line) of the new start and end, and the straight extensions
    if f0.end_class == "free_fit":
        a0 = f0.anchor_arc_um
        if f0.e_um >= 0:
            s_start, ext0 = a0, f0.e_um
        else:
            s_start, ext0 = a0 - f0.e_um, 0.0
    else:
        s_start, ext0 = 0.0, 0.0
    if f1.end_class == "free_fit":
        a1 = L - f1.anchor_arc_um
        if f1.e_um >= 0:
            s_end, ext1 = a1, f1.e_um
        else:
            s_end, ext1 = a1 + f1.e_um, 0.0
    else:
        s_end, ext1 = L, 0.0
    if s_end - s_start < 0.05:
        return Centreline(line.spl, line.qa, line.qb, line.ext0, line.ext1), 0.0, 0.0
    qa = float(line.q_of(s_start - line.ext0)[0]) if s_start > line.ext0 else line.qa
    qb = float(line.q_of(s_end - line.ext0)[0]) if s_end < line.ext0 + line.Ls else line.qb
    new = Centreline(line.spl, qa, qb, ext0, ext1)
    e0 = ext0 - s_start
    e1 = ext1 - (L - s_end)
    return new, float(e0), float(e1)


# ------------------------------------------------------------------------------------ PSF estimate
def gauss_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """y = b + a exp(-x^2 / 2 s^2) (centred): returns (s, a, b) by a grid over s then least squares in (a, b)."""
    ok = np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 5:
        return float("nan"), float("nan"), float("nan")
    best = None
    for s in np.arange(0.04, 1.2, 0.0025):
        g = np.exp(-0.5 * (x / s) ** 2)
        X = np.stack([np.ones_like(x), g], 1)
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        if coef[1] <= 0:
            continue
        r = y - X @ coef
        sse = float(r @ r)
        if best is None or sse < best[0]:
            best = (sse, s, coef[1], coef[0])
    if best is None:
        return float("nan"), float("nan"), float("nan")
    return float(best[1]), float(best[2]), float(best[3])


def psf_from_profiles(stretches: list[dict], p) -> tuple[tuple[float, float], str, dict]:
    """Pool the per-stretch mean cross-sections (lateral along u, axial along v) and fit Gaussians.
    Returns ((sigma_z, sigma_xy), source, detail)."""
    fb = tuple(float(v) for v in p["psf.fallback_sigma_um"])
    n = len(stretches)
    detail = {"n_stretches": n}
    if p["psf.mode"] == "fixed" or n < int(p["psf.min_stretches"]):
        return fb, "fallback", detail
    w = np.array([s["n"] for s in stretches], float)
    lat = np.average(np.vstack([s["lat"] for s in stretches]), axis=0, weights=w)
    ax = np.average(np.vstack([s["ax"] for s in stretches]), axis=0, weights=w)
    u = np.asarray(stretches[0]["u"], float)
    v = np.asarray(stretches[0]["v"], float)
    s_lat, a_lat, _ = gauss_fit(u, lat)
    s_ax, a_ax, _ = gauss_fit(v, ax)
    detail.update({"sigma_lat_fit": s_lat, "sigma_ax_fit": s_ax})
    if not (np.isfinite(s_lat) and np.isfinite(s_ax)) or s_lat <= 0.05 or s_ax <= 0.1:
        return fb, "fallback", detail
    return (float(s_ax), float(s_lat)), "estimated", detail
