"""The per-image PSF (protocol changelog entry 12): the flat-stretch estimate of v1 and an estimate from compact,
isolated RAD-51 foci of the run's ``spots`` table, scaled from the foci to the central_element emission
wavelength, combined by inverse variance.

* Stretches (DESIGN 4.9, v1): the mean lateral and axial cross-sections of flat (|t_z| <= ``psf.flat_tz_max``),
  isolated, bright stretches of at least ``psf.min_stretch_um``, pooled over the image and fitted with centred
  Gaussians; at least ``psf.min_stretches`` stretches are needed.
* Foci: each focus in a processed nucleus is fitted in the foci channel with a 3D Gaussian of diagonal covariance
  on a planar background, weighted by the noise model, in a window of +/- ``psf.foci_window_um`` (axial, lateral). A fit
  is kept when the window lies inside the nucleus crop, no other spot lies within ``psf.foci_isolation_um``, both
  lateral sigmas are at most ``psf.foci_max_sigma_xy_um`` with a ratio of at most ``psf.foci_max_ratio``, sigma_z
  lies inside its bounds, the relative residual RMS is at most ``psf.foci_max_rel_rms`` and the amplitude is at
  least ``psf.foci_min_snr`` standard errors. With ``psf.foci_min`` kept fits and both emission wavelengths, the
  estimate is the median sigma_z and the median sqrt(sigma_y sigma_x) times lambda_central_element / lambda_foci.
* The variance of each estimate is a bootstrap over its items (200 resamples, seed 0); with both, sigma_z is
  their inverse-variance weighted mean. The lateral sigma stays the fallback's (protocol changelog entry 13: both
  sources over-estimate it on the design phantoms); the lateral estimates are recorded. With ``psf.foci: off`` the
  v1 rule (``refine.psf_from_profiles``) is used unchanged.
* The foci fit has a planar background (entry 13), which absorbs the blurred edge of the nucleoplasm.
"""
from __future__ import annotations

import math

import numpy as np

from .refine import gauss_fit, psf_from_profiles

N_BOOT = 200
BOOT_SEED = 0
VAR_FLOOR = 25.0                       # counts^2: the floor of the fit weights (as in the entry-9 diagnostics)
SZ_BOUNDS = (0.08, 2.0)
SXY_BOUNDS = (0.05, 1.0)
CENTRE_SHIFT = (0.4, 0.3)              # um, axial and lateral
REJECT_ORDER = ("window_outside_crop", "not_isolated", "fit_failed", "lateral_sigma", "lateral_ratio",
                "sigma_z_bound", "rel_rms", "snr")
_GRID = np.arange(0.04, 1.2, 0.0025)  # the sigma grid of refine.gauss_fit


# ------------------------------------------------------------------------------------ foci fits
def fit_gauss3d_diag(win: np.ndarray, var: np.ndarray, spacing, c0_um) -> dict | None:
    """3D Gaussian with a diagonal covariance on a planar background (B + g . (x - window centre)) fitted to one
    window of raw counts, weighted by 1 / sqrt(var); ``c0_um`` the initial centre in window microns (voxel index x
    spacing). Analytic Jacobian. The plane absorbs the blurred edge of the nucleoplasm, which a constant cannot
    (protocol changelog entry 13: on design seeds most foci lie within 1 um of it)."""
    from scipy.optimize import least_squares

    sp = [float(s) for s in spacing]
    Z, Y, X = win.shape
    zz, yy, xx = np.meshgrid(np.arange(Z) * sp[0], np.arange(Y) * sp[1], np.arange(X) * sp[2], indexing="ij")
    zz, yy, xx = zz.ravel(), yy.ravel(), xx.ravel()
    zc, yc, xc = float(zz.mean()), float(yy.mean()), float(xx.mean())
    pz, py, px = zz - zc, yy - yc, xx - xc
    d = win.ravel().astype(np.float64)
    w = 1.0 / np.sqrt(np.maximum(var.ravel().astype(np.float64), VAR_FLOOR))
    bg0 = float(np.percentile(d, 10))
    ci = tuple(int(np.clip(round(c0_um[k] / sp[k]), 0, n - 1)) for k, n in enumerate(win.shape))
    A0 = max(float(win[ci]) - bg0, 1.0)
    z0, y0, x0 = (float(v) for v in c0_um)
    big = 1e6
    x_init = np.array([A0, bg0, 0.0, 0.0, 0.0, z0, y0, x0, 0.40, 0.17, 0.17])
    lo = np.array([0.0, -big, -big, -big, -big, z0 - CENTRE_SHIFT[0], y0 - CENTRE_SHIFT[1], x0 - CENTRE_SHIFT[1],
                   SZ_BOUNDS[0], SXY_BOUNDS[0], SXY_BOUNDS[0]])
    hi = np.array([10.0 * A0 + 10.0, big, big, big, big, z0 + CENTRE_SHIFT[0], y0 + CENTRE_SHIFT[1],
                   x0 + CENTRE_SHIFT[1], SZ_BOUNDS[1], SXY_BOUNDS[1], SXY_BOUNDS[1]])
    x_init = np.clip(x_init, lo + 1e-9, hi - 1e-9)

    def parts(q):
        A, B, gz, gy, gx, cz, cy, cx, sz, sy, sx = q
        dz, dy, dx = zz - cz, yy - cy, xx - cx
        E = np.exp(-0.5 * (dz * dz / (sz * sz) + dy * dy / (sy * sy) + dx * dx / (sx * sx)))
        return A, B + gz * pz + gy * py + gx * px, dz, dy, dx, sz, sy, sx, E

    def res(q):
        A, base, dz, dy, dx, sz, sy, sx, E = parts(q)
        return (base + A * E - d) * w

    def jac(q):
        A, base, dz, dy, dx, sz, sy, sx, E = parts(q)
        AE = A * E
        J = np.empty((len(d), 11))
        J[:, 0] = E
        J[:, 1] = 1.0
        J[:, 2] = pz
        J[:, 3] = py
        J[:, 4] = px
        J[:, 5] = AE * dz / (sz * sz)
        J[:, 6] = AE * dy / (sy * sy)
        J[:, 7] = AE * dx / (sx * sx)
        J[:, 8] = AE * dz * dz / sz ** 3
        J[:, 9] = AE * dy * dy / sy ** 3
        J[:, 10] = AE * dx * dx / sx ** 3
        return J * w[:, None]

    try:
        r = least_squares(res, x_init, jac=jac, bounds=(lo, hi), method="trf", x_scale="jac", max_nfev=200)
    except Exception:  # noqa: BLE001 - a failed fit is a rejected focus
        return None
    q = r.x
    A = float(q[0])
    B = float(q[1])
    cz, cy, cx, sz, sy, sx = (float(v) for v in q[5:])
    m = parts(q)[1] + A * parts(q)[-1]
    dof = max(len(d) - 11, 1)
    chi2 = float(np.sum(r.fun ** 2)) / dof
    try:
        J = r.jac
        cov = np.linalg.pinv(J.T @ J) * max(chi2, 1e-12)
        se_a = float(np.sqrt(max(cov[0, 0], 0.0)))
    except Exception:  # noqa: BLE001
        se_a = float("nan")
    return {"A": A, "B": B, "cz": cz, "cy": cy, "cx": cx, "sz": sz, "sy": sy, "sx": sx,
            "rel_rms": float(np.sqrt(np.mean((d - m) ** 2)) / max(A, 1e-9)), "se_A": se_a, "chi2_red": chi2}


def half_window(spacing, window_um) -> tuple[int, int, int]:
    ax, lat = (float(v) for v in window_um)
    return (int(round(ax / float(spacing[0]))), int(round(lat / float(spacing[1]))),
            int(round(lat / float(spacing[2]))))


def foci_fits_crop(vol: np.ndarray, org, spacing, noise, foci: dict, p) -> list[dict]:
    """Fits of the foci of one nucleus (``foci``: ``pos`` (n, 3) whole-image um and ``isolated`` (n,) bool) on the
    foci-channel crop ``vol`` whose voxel [0, 0, 0] is whole-image voxel ``org``. Each item says whether it is
    kept and, if not, the first failing rule (REJECT_ORDER)."""
    sp = [float(s) for s in spacing]
    hz, hy, hx = half_window(sp, p["psf.foci_window_um"])
    max_sxy = float(p["psf.foci_max_sigma_xy_um"])
    max_ratio = float(p["psf.foci_max_ratio"])
    max_rel = float(p["psf.foci_max_rel_rms"])
    min_snr = float(p["psf.foci_min_snr"])
    pos = np.asarray(foci.get("pos", np.zeros((0, 3))), float).reshape(-1, 3)
    iso = np.asarray(foci.get("isolated", np.ones(len(pos), bool)), bool)
    out = []
    Z, Y, X = vol.shape
    for k in range(len(pos)):
        rec = {"kept": False, "reason": ""}
        c = [int(round(pos[k, a] / sp[a] - float(org[a]))) for a in range(3)]
        if not (hz <= c[0] < Z - hz and hy <= c[1] < Y - hy and hx <= c[2] < X - hx):
            rec["reason"] = "window_outside_crop"
            out.append(rec)
            continue
        if not iso[k]:
            rec["reason"] = "not_isolated"
            out.append(rec)
            continue
        win = np.asarray(vol[c[0] - hz:c[0] + hz + 1, c[1] - hy:c[1] + hy + 1, c[2] - hx:c[2] + hx + 1], np.float64)
        w0 = [(c[a] - (hz, hy, hx)[a] + float(org[a])) * sp[a] for a in range(3)]
        c0 = [pos[k, a] - w0[a] for a in range(3)]
        f = fit_gauss3d_diag(win, noise.var(win), sp, c0)
        if f is None or not all(math.isfinite(f[x]) for x in ("A", "sz", "sy", "sx", "rel_rms")):
            rec["reason"] = "fit_failed"
            out.append(rec)
            continue
        rec.update({k_: f[k_] for k_ in ("A", "B", "sz", "sy", "sx", "rel_rms", "se_A")})
        rec["cz"], rec["cy"], rec["cx"] = f["cz"] + w0[0], f["cy"] + w0[1], f["cx"] + w0[2]
        smax, smin = max(f["sy"], f["sx"]), min(f["sy"], f["sx"])
        if smax > max_sxy:
            rec["reason"] = "lateral_sigma"
        elif smax / max(smin, 1e-9) > max_ratio:
            rec["reason"] = "lateral_ratio"
        elif not (SZ_BOUNDS[0] * 1.01 < f["sz"] < SZ_BOUNDS[1] * 0.99):
            rec["reason"] = "sigma_z_bound"
        elif f["rel_rms"] > max_rel:
            rec["reason"] = "rel_rms"
        elif not (math.isfinite(f["se_A"]) and f["A"] >= min_snr * f["se_A"]):
            rec["reason"] = "snr"
        else:
            rec["kept"] = True
        out.append(rec)
    return out


def isolated_spots(pos: np.ndarray, radius_um: float) -> np.ndarray:
    """True for the spots with no other spot within ``radius_um`` (3D)."""
    from scipy.spatial import cKDTree

    P = np.asarray(pos, float).reshape(-1, 3)
    if len(P) < 2:
        return np.ones(len(P), bool)
    cnt = np.array([len(x) for x in cKDTree(P).query_ball_point(P, float(radius_um))])
    return cnt <= 1


def _median_boot(x: np.ndarray, rng) -> float:
    idx = rng.integers(0, len(x), size=(N_BOOT, len(x)))
    return float(np.std(np.median(x[idx], axis=1), ddof=1))


def foci_estimate(fits: list[dict], wavelengths: dict | None, p) -> dict:
    """The foci estimate (``ok`` False when it cannot be made) and its bookkeeping."""
    kept = [f for f in fits if f.get("kept")]
    reasons = {r: sum(1 for f in fits if f.get("reason") == r) for r in REJECT_ORDER}
    out = {"ok": False, "n_fits": len(fits), "n_kept": len(kept), "rejected": reasons}
    lam_ce = (wavelengths or {}).get("central_element")
    lam_f = (wavelengths or {}).get("foci")
    out["wavelengths_nm"] = {"central_element": lam_ce, "foci": lam_f}
    if not (lam_ce and lam_f and math.isfinite(float(lam_ce)) and math.isfinite(float(lam_f))):
        out["why"] = "no_wavelengths"
        return out
    ratio = float(lam_ce) / float(lam_f)
    out["wavelength_ratio"] = ratio
    if len(kept) < int(p["psf.foci_min"]):
        out["why"] = "too_few_kept"
        return out
    sz = np.array([f["sz"] for f in kept], float)
    sxy = np.sqrt(np.array([f["sy"] for f in kept], float) * np.array([f["sx"] for f in kept], float))
    rng = np.random.default_rng(BOOT_SEED)
    out.update({"ok": True, "sigma_z_um": float(np.median(sz)) * ratio, "sigma_xy_um": float(np.median(sxy)) * ratio,
                "se_z_um": _median_boot(sz, rng) * ratio, "se_xy_um": _median_boot(sxy, rng) * ratio})
    return out


# ------------------------------------------------------------------------------------ stretches
def _gauss_grid_fit(x: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """refine.gauss_fit's sigma for each row of ``Y`` (profiles on ``x``), vectorised (closed-form least squares
    of y = b + a g over the same sigma grid, a > 0)."""
    ok = np.isfinite(Y).all(0)
    x, Y = x[ok], Y[:, ok]
    n = len(x)
    if n < 5:
        return np.full(len(Y), np.nan)
    G = np.exp(-0.5 * (x[None, :] / _GRID[:, None]) ** 2)            # (S, n)
    Sg, Sgg = G.sum(1), (G * G).sum(1)
    Sy = Y.sum(1)                                                      # (B,)
    Sgy = Y @ G.T                                                      # (B, S)
    det = n * Sgg - Sg * Sg
    a = (n * Sgy - Sg[None, :] * Sy[:, None]) / np.where(np.abs(det) > 1e-12, det, np.nan)[None, :]
    b = (Sy[:, None] - a * Sg[None, :]) / n
    Syy = (Y * Y).sum(1)
    sse = Syy[:, None] - 2 * b * Sy[:, None] - 2 * a * Sgy + n * b * b + 2 * a * b * Sg[None, :] + a * a * Sgg[None, :]
    sse = np.where(a > 0, sse, np.inf)
    k = np.argmin(sse, axis=1)
    s = _GRID[k]
    s[~np.isfinite(sse[np.arange(len(Y)), k])] = np.nan
    return s


def pooled_profiles(stretches: list[dict]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(u, v, lat (k, nu), ax (k, nv), n (k,)) of a list of stretches."""
    u = np.asarray(stretches[0]["u"], float)
    v = np.asarray(stretches[0]["v"], float)
    lat = np.vstack([s["lat"] for s in stretches]).astype(float)
    ax = np.vstack([s["ax"] for s in stretches]).astype(float)
    n = np.array([s["n"] for s in stretches], float)
    return u, v, lat, ax, n


def stretch_estimate(stretches: list[dict], p, n_min: int | None = None, with_se: bool = True) -> dict:
    """The v1 flat-stretch estimate with at least ``psf.min_stretches`` (or ``n_min``) stretches, and its bootstrap
    standard errors."""
    n_min = int(p["psf.min_stretches"]) if n_min is None else int(n_min)
    out = {"ok": False, "n_stretches": len(stretches)}
    if len(stretches) < n_min or not stretches:
        out["why"] = "too_few_stretches"
        return out
    u, v, lat, ax, w = pooled_profiles(stretches)
    s_lat, _, _ = gauss_fit(u, np.average(lat, axis=0, weights=w))
    s_ax, _, _ = gauss_fit(v, np.average(ax, axis=0, weights=w))
    out.update({"sigma_lat_fit": s_lat, "sigma_ax_fit": s_ax})
    if not (np.isfinite(s_lat) and np.isfinite(s_ax)) or s_lat <= 0.05 or s_ax <= 0.1:
        out["why"] = "fit_failed"
        return out
    out.update({"ok": True, "sigma_z_um": float(s_ax), "sigma_xy_um": float(s_lat)})
    if with_se:
        rng = np.random.default_rng(BOOT_SEED)
        idx = rng.integers(0, len(w), size=(N_BOOT, len(w)))
        W = w[idx]                                                     # (B, k)
        lat_b = np.einsum("bk,bkn->bn", W, lat[idx]) / W.sum(1, keepdims=True)
        ax_b = np.einsum("bk,bkn->bn", W, ax[idx]) / W.sum(1, keepdims=True)
        sl, sa = _gauss_grid_fit(u, lat_b), _gauss_grid_fit(v, ax_b)
        out["se_xy_um"] = float(np.nanstd(sl, ddof=1))
        out["se_z_um"] = float(np.nanstd(sa, ddof=1))
    return out


# ------------------------------------------------------------------------------------ the image PSF
def _ivw(a: float, sa: float, b: float, sb: float) -> float:
    """Inverse-variance weighted mean; a standard error that is not finite carries no weight (both: the mean)."""
    fa, fb_ = math.isfinite(sa), math.isfinite(sb)
    if not fa and not fb_:
        return float(0.5 * (a + b))
    wa = 1.0 / max(sa, 1e-4) ** 2 if fa else 0.0
    wb = 1.0 / max(sb, 1e-4) ** 2 if fb_ else 0.0
    return float((wa * a + wb * b) / (wa + wb))


def image_psf(stretches: list[dict], foci_fits: list[dict] | None, wavelengths: dict | None, p
              ) -> tuple[tuple[float, float], str, dict]:
    """((sigma_z, sigma_xy), source, detail) of one image. With ``psf.foci: off`` this is v1's
    ``psf_from_profiles`` (sources estimated and fallback). Otherwise sigma_z comes from the estimates (protocol
    changelog entries 12 and 13) and sigma_xy stays the fallback's lateral sigma (both sources over-estimate it on
    the design phantoms; the lateral estimates are recorded and choose nothing)."""
    from .params import switch_on

    if not switch_on(p["psf.foci"]) or p["psf.mode"] == "fixed":
        return psf_from_profiles(stretches, p)
    fb = tuple(float(v) for v in p["psf.fallback_sigma_um"])
    st = stretch_estimate(stretches, p)
    fo = foci_estimate(foci_fits or [], wavelengths, p)
    detail = {"n_stretches": len(stretches), "stretches": st, "foci": fo, "sigma_xy_from": "fallback"}
    if "sigma_lat_fit" in st:
        detail.update({"sigma_lat_fit": st["sigma_lat_fit"], "sigma_ax_fit": st["sigma_ax_fit"]})
    if st["ok"] and fo["ok"]:
        sz = _ivw(st["sigma_z_um"], st.get("se_z_um", np.nan), fo["sigma_z_um"], fo["se_z_um"])
        detail["foci_minus_stretches_z_um"] = fo["sigma_z_um"] - st["sigma_z_um"]
        return (sz, fb[1]), "both", detail
    if st["ok"]:
        return (st["sigma_z_um"], fb[1]), "stretches", detail
    if fo["ok"]:
        return (fo["sigma_z_um"], fb[1]), "foci", detail
    return fb, "fallback", detail
