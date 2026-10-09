"""EXPERIMENTAL classical spot baseline: anisotropy-aware Laplacian-of-Gaussian (LoG) blob detection.

This is the first in-house seed of a spot detector, registered as the ``log`` backend in
spots/backends.py. It is NEVER selected by default (``spots.backend`` must be set to ``log``), it has
not been validated against the lab's Imaris counts, and its numbers are not comparable with the
SpotMAX defaults. It exists so that a candidate detector can be run through the same tables, flags
and goldens as the validated path and compared head to head.

Chain, per nucleus label (each label is processed in its own padded bounding-box crop, so memory
stays small on a 2600 x 2600 x 90 stack):

1. LoG response ``-gaussian_laplace(crop, sigma)`` with a per-axis sigma in voxels from one physical
   sigma: the standard deviation of a Gaussian whose full width at half maximum equals the expected
   focus diameter (``spots.log.diameter_um``). Bright blobs give a positive response.
2. Local maxima: a voxel is a candidate when it is the maximum of the response inside an ellipsoid of
   radius ``spots.log.min_separation_um`` (physical units, so z and xy are treated alike), its response
   is positive and it lies inside the label. A greedy pass in physical coordinates settles ties.
3. Flood guard: ``_cap_candidates`` (the SpotMAX path's guard, same key ``spots.max_spot_candidates``).
4. Signal-to-background gate: the raw peak (the brightest voxel in a core ellipsoid of radius half the
   diameter) divided by the median of the raw image in a local annulus (``annulus_inner_um`` to
   ``annulus_outer_um``, restricted to the same nucleus label and excluding the cores of the other
   candidates when enough voxels are there) must reach ``spots.log.min_snr``.
5. ``_merge_z_columns`` (the SpotMAX path's z-split merge) on the integer peak positions, with the
   raw image for the valley test; the shared keys ``spots.merge_z_columns``, ``spots.z_merge_gap_um``
   and ``spots.z_merge_valley_frac`` apply.
6. Sub-voxel refinement: the intensity-weighted centroid of the background-subtracted core voxels.
7. Nucleus assignment: the label at the refined position; when that is 0 (the centroid drifted just
   outside the mask) the label of the peak voxel, which is never 0 because candidates are only taken
   inside labels. The SpotMAX path has no nearest-label fallback (spots outside every label are
   dropped), so none is added here either.

Columns match the SpotMAX tables exactly. ``spot_mean_intensity`` is the mean raw intensity of the
core voxels; ``effect_size`` is a Glass's delta of core mean over annulus mean in annulus standard
deviations (the same idea as SpotMAX's ``spot_vs_backgr_effect_size_glass``, not the same number).
"""
from __future__ import annotations

import logging
import math

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

from .detect import _cap_candidates, _empty, _finish_tables, _merge_z_columns

log = logging.getLogger(__name__)

DETECTOR = "log"
PRE_GATE_SIGMAS = 3.0   # candidate pre-gate: LoG response above this many robust SDs of the in-nucleus response

# spots.log.* keys and their code defaults (None = derived from the diameter, see _resolve)
LOG_DEFAULTS: dict = {
    "diameter_um": 0.5,          # expected focus diameter (sets the LoG sigma and the core radius)
    "min_snr": 2.0,              # raw peak over local annulus median
    "min_separation_um": None,   # minimum distance between two maxima; default half the diameter
    "annulus_inner_um": None,    # local background annulus; defaults: diameter to twice the diameter
    "annulus_outer_um": None,
}
# keys shared with the SpotMAX path (read from the top of the spots section, same code defaults)
SHARED_DEFAULTS: dict = {
    "merge_z_columns": True, "z_merge_gap_um": 0.8, "z_merge_valley_frac": 0.8, "max_spot_candidates": 30000,
}
_FWHM_PER_SIGMA = 2.0 * math.sqrt(2.0 * math.log(2.0))
_MIN_ANNULUS_VOXELS = 8      # below this the same-label restriction of the annulus is dropped
_warned = False


def log_kwargs(cfg_section: dict) -> dict:
    """The keyword arguments of ``detect_spots_log`` from the raw ``spots:`` section: its own keys under
    ``spots.log`` and the shared merge / flood-guard keys, every one with a code default."""
    sub = cfg_section.get("log", None)
    if sub is None:
        sub = {}
    if not isinstance(sub, dict):
        raise TypeError(f"spots.log must be a mapping of the log backend's parameters, got {type(sub).__name__}")
    unknown = sorted(set(sub) - set(LOG_DEFAULTS))
    if unknown:
        log.warning("spots.log: unknown keys %s ignored (known: %s)", unknown, sorted(LOG_DEFAULTS))
    kw = {k: sub.get(k, d) for k, d in LOG_DEFAULTS.items()}
    kw.update({k: cfg_section.get(k, d) for k, d in SHARED_DEFAULTS.items()})
    return kw


def detect_spots_log(
    spot_img: np.ndarray,
    labels: np.ndarray,
    spacing: tuple[float, float, float],
    *,
    marker: str = "RAD-51",
    diameter_um: float = 0.5,
    min_snr: float = 2.0,
    min_separation_um: float | None = None,
    annulus_inner_um: float | None = None,
    annulus_outer_um: float | None = None,
    merge_z_columns: bool = True,
    z_merge_gap_um: float = 0.8,
    z_merge_valley_frac: float = 0.8,
    max_spot_candidates: int = 30000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (per_spot_df, per_nucleus_df) with the SpotMAX columns and detector ``log``. See the
    module docstring for the chain; every length is in microns and converted with ``spacing``."""
    _warn_experimental()
    img = np.asarray(spot_img)
    lab = np.asarray(labels)
    if img.ndim != 3 or lab.shape != img.shape:
        raise ValueError("the log backend needs a (Z, Y, X) foci image and a label image of the same shape, "
                         f"got {img.shape} and {lab.shape}")
    sp = np.asarray(spacing, dtype=float)
    if sp.shape != (3,) or not np.all(sp > 0):
        raise ValueError(f"spacing must be three positive voxel sizes (dz, dy, dx) in microns, got {spacing!r}")
    if not np.issubdtype(lab.dtype, np.integer):
        lab = lab.astype(np.int64)
    p = _resolve(diameter_um, min_snr, min_separation_um, annulus_inner_um, annulus_outer_um)

    sigma_vox = p["sigma_um"] / sp
    footprint = _footprint(p["separation_um"], sp)
    core_off = _ellipsoid_offsets(p["core_um"], sp)
    shell_off = _shell_offsets(p["annulus_inner_um"], p["annulus_outer_um"], sp)
    pad = np.ceil(max(4.0 * p["sigma_um"], p["separation_um"]) / sp).astype(int) + 1

    cand = _candidates(img, lab, sigma_vox, footprint, pad)
    if cand.empty:
        return _empty(lab, marker, DETECTOR)
    # flood guard, exactly as the SpotMAX path: the brightest raw voxels survive, with a loud warning
    cand = _cap_candidates(cand, img, int(max_spot_candidates))
    cand = _suppress(cand, sp, p["separation_um"])
    cand = _measure(cand, img, lab, sp, core_off, shell_off, p["core_um"], p["annulus_outer_um"])
    cand = cand[cand["snr"] >= p["min_snr"]].reset_index(drop=True)
    if cand.empty:
        return _empty(lab, marker, DETECTOR)

    # _merge_z_columns IS applied here, on the integer peak positions (before refinement, so stacked
    # peaks share an exact (x, y) column) with the raw foci channel for the valley test. The anisotropic
    # non-maximum suppression already removes stacked duplicates closer than the separation radius;
    # this is the second guard against a PSF-elongated focus splitting along z, as in the SpotMAX path.
    if bool(merge_z_columns) and len(cand) > 1:
        cand = _merge_z_columns(cand, sp, float(z_merge_gap_um), img, float(z_merge_valley_frac))

    zr = np.clip(np.rint(cand["zr"].to_numpy()).astype(int), 0, img.shape[0] - 1)
    yr = np.clip(np.rint(cand["yr"].to_numpy()).astype(int), 0, img.shape[1] - 1)
    xr = np.clip(np.rint(cand["xr"].to_numpy()).astype(int), 0, img.shape[2] - 1)
    at_refined = lab[zr, yr, xr]
    nucleus_id = np.where(at_refined > 0, at_refined, cand["nucleus_id"].to_numpy())
    out = pd.DataFrame({
        "nucleus_id": nucleus_id.astype(int),
        "z": cand["zr"].to_numpy(), "y": cand["yr"].to_numpy(), "x": cand["xr"].to_numpy(),
        "effect_size": cand["effect_size"].to_numpy(),
        "spot_mean_intensity": cand["spot_mean_intensity"].to_numpy(),
    })
    out = out.sort_values(["nucleus_id", "z", "y", "x"], kind="stable").reset_index(drop=True)
    return _finish_tables(out, lab, sp, marker, DETECTOR)


# ---- parameters and geometry ----

def _resolve(diameter_um, min_snr, min_separation_um, annulus_inner_um, annulus_outer_um) -> dict:
    d = float(diameter_um)
    if not d > 0:
        raise ValueError(f"spots.log.diameter_um must be positive, got {diameter_um!r}")
    sep = d / 2 if min_separation_um is None else float(min_separation_um)
    r_in = d if annulus_inner_um is None else float(annulus_inner_um)
    r_out = 2 * d if annulus_outer_um is None else float(annulus_outer_um)
    if not sep > 0:
        raise ValueError(f"spots.log.min_separation_um must be positive, got {min_separation_um!r}")
    if not (r_out > r_in > 0):
        raise ValueError("spots.log.annulus_outer_um must exceed spots.log.annulus_inner_um and both must be "
                         f"positive, got {annulus_inner_um!r} and {annulus_outer_um!r}")
    return {"sigma_um": d / _FWHM_PER_SIGMA, "core_um": d / 2, "separation_um": sep,
            "annulus_inner_um": r_in, "annulus_outer_um": r_out, "min_snr": float(min_snr)}


def _ellipsoid_offsets(radius_um: float, sp: np.ndarray) -> np.ndarray:
    """Integer voxel offsets (n, 3) whose physical distance from the origin is at most ``radius_um``."""
    r = np.ceil(radius_um / sp).astype(int)
    grid = np.mgrid[-r[0]:r[0] + 1, -r[1]:r[1] + 1, -r[2]:r[2] + 1].reshape(3, -1).T
    d2 = ((grid * sp) ** 2).sum(axis=1)
    return grid[d2 <= radius_um ** 2 + 1e-12]


def _shell_offsets(inner_um: float, outer_um: float, sp: np.ndarray) -> np.ndarray:
    """Offsets of the annulus: farther than ``inner_um`` and at most ``outer_um`` away."""
    off = _ellipsoid_offsets(outer_um, sp)
    d2 = ((off * sp) ** 2).sum(axis=1)
    return off[d2 > inner_um ** 2 + 1e-12]


def _footprint(radius_um: float, sp: np.ndarray) -> np.ndarray:
    """Boolean maximum-filter footprint: the ellipsoid of ``radius_um`` (always includes the centre)."""
    off = _ellipsoid_offsets(radius_um, sp)
    r = np.abs(off).max(axis=0)
    fp = np.zeros(tuple(2 * r + 1), dtype=bool)
    fp[tuple((off + r).T)] = True
    return fp


# ---- the chain ----

def _candidates(img, lab, sigma_vox, footprint, pad) -> pd.DataFrame:
    """Local maxima of the LoG response inside every label, each label in its own padded crop."""
    objs = ndi.find_objects(lab) if lab.size and lab.max() > 0 else []
    rows = []
    for lab_id, sl in enumerate(objs, start=1):
        if sl is None:
            continue
        lo = [max(0, s.start - int(q)) for s, q in zip(sl, pad)]
        hi = [min(n, s.stop + int(q)) for s, q, n in zip(sl, pad, img.shape)]
        win = tuple(slice(a, b) for a, b in zip(lo, hi))
        inside = lab[win] == lab_id
        if not inside.any():
            continue
        crop = np.asarray(img[win], dtype=np.float32)
        resp = -ndi.gaussian_laplace(crop, sigma=sigma_vox, mode="nearest")
        # pre-gate: a local maximum must clear a noise-scaled floor (PRE_GATE_SIGMAS robust standard
        # deviations of the response inside the nucleus) before it can count as a candidate, so that on a
        # noisy or empty nucleus the flood cap is never filled with noise maxima that crowd out real foci.
        vals = resp[inside]
        noise = 1.4826 * float(np.median(np.abs(vals - np.median(vals)))) if vals.size else 0.0
        floor = PRE_GATE_SIGMAS * noise
        peaks = (resp >= ndi.maximum_filter(resp, footprint=footprint, mode="nearest")) & (resp > floor) & inside
        idx = np.argwhere(peaks)
        if len(idx) == 0:
            continue
        gpos = idx + np.array(lo)
        rows.append(pd.DataFrame({"nucleus_id": lab_id, "z": gpos[:, 0], "y": gpos[:, 1], "x": gpos[:, 2],
                                  "log_response": resp[peaks].astype(np.float64)}))
    if not rows:
        return pd.DataFrame(columns=["nucleus_id", "z", "y", "x", "log_response"])
    return pd.concat(rows, ignore_index=True)


def _suppress(cand: pd.DataFrame, sp: np.ndarray, sep_um: float) -> pd.DataFrame:
    """Greedy non-maximum suppression in physical coordinates: of two candidates closer than ``sep_um``
    the stronger response wins (ties: lower z, y, x). The footprint filter already enforces this for
    strict maxima; this settles plateaus deterministically."""
    if len(cand) < 2:
        return cand
    pts = cand[["z", "y", "x"]].to_numpy(dtype=float) * sp
    order = np.lexsort((cand["x"].to_numpy(), cand["y"].to_numpy(), cand["z"].to_numpy(),
                        -cand["log_response"].to_numpy()))
    tree = cKDTree(pts)
    keep = np.ones(len(cand), dtype=bool)
    for i in order:
        if not keep[i]:
            continue
        for j in tree.query_ball_point(pts[i], sep_um):
            if j != i:
                keep[j] = False
    return cand[keep].reset_index(drop=True)


def _measure(cand: pd.DataFrame, img, lab, sp, core_off, shell_off, core_um, outer_um) -> pd.DataFrame:
    """Per candidate: raw peak, local background (annulus median / mean / sd), core mean, Glass's delta,
    signal-to-background ratio and the background-subtracted intensity-weighted centroid (zr, yr, xr).
    The annulus is restricted to the candidate's own nucleus label and excludes the core of every other
    candidate, so a neighbouring focus in a dense nucleus does not inflate the background (the
    restriction is dropped when fewer than _MIN_ANNULUS_VOXELS voxels would remain)."""
    shape = np.array(img.shape)
    pos = cand[["z", "y", "x"]].to_numpy(dtype=int)
    pos_um = pos * sp
    nuc = cand["nucleus_id"].to_numpy()
    tree = cKDTree(pos_um)
    n = len(cand)
    peak = np.zeros(n)
    bg_med = np.zeros(n)
    bg_mean = np.zeros(n)
    bg_std = np.zeros(n)
    core_mean = np.zeros(n)
    refined = pos.astype(float)
    for i in range(n):
        p = pos[i]
        c = p + core_off
        ok = np.all((c >= 0) & (c < shape), axis=1)
        c, co = c[ok], core_off[ok]
        cv = img[c[:, 0], c[:, 1], c[:, 2]].astype(np.float64)
        s = p + shell_off
        ok = np.all((s >= 0) & (s < shape), axis=1)
        s = s[ok]
        sv = img[s[:, 0], s[:, 1], s[:, 2]].astype(np.float64)
        clean = lab[s[:, 0], s[:, 1], s[:, 2]] == nuc[i]
        others = [j for j in tree.query_ball_point(pos_um[i], outer_um + core_um) if j != i]
        if others:
            d2 = ((s[:, None, :] * sp - pos_um[others][None, :, :]) ** 2).sum(axis=2)
            clean &= (d2 > core_um ** 2).all(axis=1)
        if clean.sum() >= _MIN_ANNULUS_VOXELS:
            sv = sv[clean]
        if len(sv):
            bg_med[i], bg_mean[i] = np.median(sv), sv.mean()
            bg_std[i] = sv.std(ddof=1) if len(sv) > 1 else 0.0
        peak[i], core_mean[i] = cv.max(), cv.mean()
        w = np.maximum(cv - bg_med[i], 0.0)
        if w.sum() > 0:
            refined[i] = p + (w @ co) / w.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = np.where(bg_med > 0, peak / bg_med, np.where(peak > 0, np.inf, 0.0))
        effect = np.where(bg_std > 0, (core_mean - bg_mean) / bg_std,
                          np.where(core_mean > bg_mean, np.inf, np.nan))
    return cand.assign(peak=peak, bg_median=bg_med, snr=snr, spot_mean_intensity=core_mean,
                       effect_size=effect, zr=refined[:, 0], yr=refined[:, 1], xr=refined[:, 2])


def _warn_experimental() -> None:
    global _warned
    if not _warned:
        _warned = True
        log.warning("the 'log' spot backend is EXPERIMENTAL: a classical Laplacian-of-Gaussian baseline that "
                    "is not validated against Imaris; its counts are not comparable with the spotmax defaults.")
