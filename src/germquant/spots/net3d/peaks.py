"""Peak finding on a net3d heatmap and nucleus assignment, shared by the backend, the training-time
validation and the tests (numpy and scipy only, no torch).

1. Candidates: voxels whose probability reaches the threshold and equals the maximum of the heatmap
   over an ellipsoid with the separation radii (z, xy) in microns.
2. Ellipsoidal non-maximum suppression: candidates are visited from the most to the least probable
   (ties: lower z, y, x first) and every candidate closer than the separation ellipsoid to a kept one
   is dropped. The kept set at a higher threshold is exactly the kept set at a lower one restricted
   to the candidates above the higher threshold (a candidate's fate depends only on more probable
   candidates), which is what lets the validation sweep thresholds from one pass.
3. Sub-voxel refinement: per axis, the vertex of the parabola through the log probabilities at the
   peak and its two neighbours (exact for a Gaussian blob), clipped to half a voxel; 0 at the border or
   when the three values do not curve downwards.
4. Nucleus assignment: the label at the refined position (rounded to the nearest voxel); when that is
   background, the nearest labelled voxel within ``assign_max_um`` (physical distance from the refined
   position to voxel centres; ties go to the smaller label); no label within that distance leaves the
   spot unassigned (nucleus 0, dropped by the backend, as SpotMAX drops spots outside every nucleus).
   ``contested`` is true when two or more labels lie within ``assign_max_um`` of the refined position,
   whichever rule assigned the spot.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

DEFAULT_SEPARATION_UM = (0.6, 0.3)     # (z, xy) radii of the suppression ellipsoid
DEFAULT_ASSIGN_MAX_UM = 0.3
_EPS = 1e-12


def ellipsoid_offsets(radii_um, spacing, *, include_center: bool = True) -> np.ndarray:
    """Integer (dz, dy, dx) offsets whose physical position lies inside the ellipsoid with semi-axes
    ``radii_um`` = (rz, ry, rx) (inclusive of the boundary)."""
    r = np.asarray(radii_um, float)
    sp = np.asarray(spacing, float)
    if r.shape == (2,):
        r = np.array([r[0], r[1], r[1]])
    n = np.floor(r / sp + 1e-9).astype(int)
    grid = np.mgrid[-n[0]:n[0] + 1, -n[1]:n[1] + 1, -n[2]:n[2] + 1].reshape(3, -1).T
    q = (((grid * sp) / np.maximum(r, _EPS)) ** 2).sum(axis=1)
    off = grid[q <= 1.0 + 1e-9]
    if not include_center:
        off = off[np.any(off != 0, axis=1)]
    return off


def _gather(heat: np.ndarray, z, y, x, fill: float) -> np.ndarray:
    Z, Y, X = heat.shape
    ok = (z >= 0) & (z < Z) & (y >= 0) & (y < Y) & (x >= 0) & (x < X)
    out = np.full(len(z), fill, dtype=np.float64)
    if ok.any():
        out[ok] = heat[z[ok], y[ok], x[ok]]
    return out


def candidate_peaks(heat: np.ndarray, threshold: float, spacing, *,
                    separation_um=DEFAULT_SEPARATION_UM, region=None) -> pd.DataFrame:
    """Voxels with ``heat >= threshold`` that are the maximum over the separation ellipsoid (plateaus
    keep every tied voxel; ``suppress`` settles them). ``region`` is an optional tuple of slices: only
    voxels inside it can be candidates (the neighbourhood test still reads the whole array). Returns
    z, y, x (int) and prob."""
    heat = np.asarray(heat)
    sub = heat if region is None else heat[region]
    zz, yy, xx = np.nonzero(sub >= threshold)
    if region is not None and len(zz):
        zz = zz + (region[0].start or 0)
        yy = yy + (region[1].start or 0)
        xx = xx + (region[2].start or 0)
    if len(zz) == 0:
        return pd.DataFrame({"z": np.zeros(0, int), "y": np.zeros(0, int), "x": np.zeros(0, int),
                             "prob": np.zeros(0, float)})
    vals = heat[zz, yy, xx].astype(np.float64)
    keep = np.ones(len(zz), bool)
    for dz, dy, dx in ellipsoid_offsets(separation_um, spacing, include_center=False):
        idx = np.flatnonzero(keep)
        if idx.size == 0:
            break
        nb = _gather(heat, zz[idx] + dz, yy[idx] + dy, xx[idx] + dx, -np.inf)
        keep[idx[nb > vals[idx]]] = False
    return pd.DataFrame({"z": zz[keep].astype(int), "y": yy[keep].astype(int), "x": xx[keep].astype(int),
                         "prob": vals[keep]})


def suppress(cand: pd.DataFrame, spacing, separation_um=DEFAULT_SEPARATION_UM) -> pd.DataFrame:
    """Greedy ellipsoidal non-maximum suppression in physical coordinates (module docstring, step 2).
    Returns the kept rows in visiting order (most probable first)."""
    if len(cand) == 0:
        return cand.reset_index(drop=True)
    r = np.asarray(separation_um, float)
    if r.shape == (2,):
        r = np.array([r[0], r[1], r[1]])
    pts = cand[["z", "y", "x"]].to_numpy(float) * np.asarray(spacing, float) / r
    order = np.lexsort((cand["x"].to_numpy(), cand["y"].to_numpy(), cand["z"].to_numpy(),
                        -cand["prob"].to_numpy()))
    tree = cKDTree(pts)
    alive = np.ones(len(cand), bool)
    kept = []
    for i in order:
        if not alive[i]:
            continue
        kept.append(i)
        for j in tree.query_ball_point(pts[i], 1.0 - 1e-9):   # strictly inside the ellipsoid
            alive[j] = False
    return cand.iloc[kept].reset_index(drop=True)


def refine_subvoxel(heat: np.ndarray, cand: pd.DataFrame) -> np.ndarray:
    """(n, 3) refined positions in voxel units (module docstring, step 3)."""
    n = len(cand)
    pos = np.array(cand[["z", "y", "x"]].to_numpy(float), dtype=np.float64, copy=True).reshape(n, 3)
    if n == 0:
        return pos
    base = np.array(cand[["z", "y", "x"]].to_numpy(), dtype=np.int64, copy=True)
    c = np.log(np.maximum(heat[base[:, 0], base[:, 1], base[:, 2]].astype(np.float64), _EPS))
    for a in range(3):
        lo, hi = base.copy(), base.copy()
        lo[:, a] -= 1
        hi[:, a] += 1
        inside = (lo[:, a] >= 0) & (hi[:, a] < heat.shape[a])
        if not inside.any():
            continue
        vl = np.full(n, np.nan)
        vh = np.full(n, np.nan)
        vl[inside] = np.log(np.maximum(heat[lo[inside, 0], lo[inside, 1], lo[inside, 2]].astype(np.float64), _EPS))
        vh[inside] = np.log(np.maximum(heat[hi[inside, 0], hi[inside, 1], hi[inside, 2]].astype(np.float64), _EPS))
        denom = vl - 2.0 * c + vh
        with np.errstate(invalid="ignore", divide="ignore"):
            off = np.where(inside & (denom < -_EPS), 0.5 * (vl - vh) / denom, 0.0)
        pos[:, a] += np.clip(np.nan_to_num(off), -0.5, 0.5)
    return pos


def assign_nuclei(pos_vox: np.ndarray, labels: np.ndarray, spacing, *,
                  assign_max_um: float = DEFAULT_ASSIGN_MAX_UM) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(nucleus_id, rule, contested) for refined positions in voxel units (module docstring, step 4).
    ``rule`` is "label" (inside a label), "nearest" (the fallback) or "none" (unassigned, id 0)."""
    pos = np.asarray(pos_vox, float).reshape(-1, 3)
    n = len(pos)
    sp = np.asarray(spacing, float)
    shape = np.asarray(labels.shape)
    if n == 0:
        return np.zeros(0, np.int64), np.zeros(0, dtype=object), np.zeros(0, bool)
    r = np.clip(np.rint(pos).astype(np.int64), 0, shape - 1)
    at = np.asarray(labels[r[:, 0], r[:, 1], r[:, 2]]).astype(np.int64)
    # every voxel centre that can lie within assign_max_um of a position rounded to its voxel
    reach = np.ceil(assign_max_um / sp + 0.5).astype(int)
    grid = np.mgrid[-reach[0]:reach[0] + 1, -reach[1]:reach[1] + 1, -reach[2]:reach[2] + 1].reshape(3, -1).T
    best_lab = np.zeros(n, np.int64)
    best_d2 = np.full(n, np.inf)
    hit_k, hit_lab = [], []
    lim2 = assign_max_um ** 2 + 1e-12
    for off in grid:
        v = r + off
        ok = np.all((v >= 0) & (v < shape), axis=1)
        if not ok.any():
            continue
        lab = np.zeros(n, np.int64)
        lab[ok] = np.asarray(labels[v[ok, 0], v[ok, 1], v[ok, 2]])
        d2 = (((v - pos) * sp) ** 2).sum(axis=1)
        hit = (lab > 0) & (d2 <= lim2)
        if not hit.any():
            continue
        hit_k.append(np.flatnonzero(hit))
        hit_lab.append(lab[hit])
        better = hit & ((d2 < best_d2 - 1e-12) | ((np.abs(d2 - best_d2) <= 1e-12) & (lab < best_lab)))
        best_lab[better] = lab[better]
        best_d2[better] = d2[better]
    nucleus = np.where(at > 0, at, best_lab)
    rule = np.where(at > 0, "label", np.where(best_lab > 0, "nearest", "none")).astype(object)
    contested = np.zeros(n, bool)
    if hit_k:
        pairs = np.unique(np.stack([np.concatenate(hit_k), np.concatenate(hit_lab)], axis=1), axis=0)
        n_labels = np.bincount(pairs[:, 0], minlength=n)
        contested = n_labels > 1
    return nucleus, rule, contested


def find_spots(heat: np.ndarray, labels: np.ndarray, spacing, threshold: float, *,
               separation_um=DEFAULT_SEPARATION_UM, assign_max_um: float = DEFAULT_ASSIGN_MAX_UM,
               max_candidates: int | None = None) -> pd.DataFrame:
    """The whole chain at one threshold: candidates, suppression, refinement, assignment. Returns one
    row per kept peak (assigned or not): z, y, x (int peak voxel), zr, yr, xr (refined, voxel units),
    prob, nucleus_id (0 = unassigned), rule, contested, in suppression order."""
    cand = candidate_peaks(heat, threshold, spacing, separation_um=separation_um)
    if max_candidates is not None and len(cand) > int(max_candidates):
        order = np.lexsort((cand["x"].to_numpy(), cand["y"].to_numpy(), cand["z"].to_numpy(),
                            -cand["prob"].to_numpy()))
        cand = cand.iloc[np.sort(order[:int(max_candidates)])].reset_index(drop=True)
    kept = suppress(cand, spacing, separation_um)
    pos = refine_subvoxel(heat, kept)
    nucleus, rule, contested = assign_nuclei(pos, labels, spacing, assign_max_um=assign_max_um)
    return kept.assign(zr=pos[:, 0], yr=pos[:, 1], xr=pos[:, 2], nucleus_id=nucleus, rule=rule,
                       contested=contested)
