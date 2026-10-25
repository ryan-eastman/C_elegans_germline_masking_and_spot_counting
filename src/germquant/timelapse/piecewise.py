"""Piecewise registration of a time-lapse series along the gonad axis (``germquant register --mode piecewise``;
docs/TIMELAPSE.md "Piecewise registration").

Why. On the first real series (the October 2026 GFP::SYP-2 hermaphrodite) the worm does not drift, but the
germline slides inside it non-rigidly: the proximal half moves up to 3.6 um per 5 minutes while the distal end
barely moves. One translation per frame (``germquant register``) cannot describe that: registering on GFP
breaks the distal links, registering on brightfield the proximal ones.

What. At the origin (the first timepoint with a reference volume) the germline nuclei are split into K segments
of equal nucleus count along the fixed axis (track.fixed_axis_model: the origin's axis_position_um, smoothed by
a local linear fit over the 16 nearest nuclei; a fitted centerline when the axis stage did not run). For every
later timepoint and every segment, the translation of the segment between the previous registrable frame and
this one is estimated exactly as ``germquant register`` estimates a whole frame (register.displacement: the
same phase correlation, upsampled peak refinement, normalization and step cap), on a crop of the two reference
volumes: the y and x extent of the segment's origin nuclei, moved to where the segment sat in the previous frame
and widened by ``margin_um``, every plane in z, mean-subtracted and tapered (``taper``, a Tukey window, because
the crop always cuts the gonad where the segment ends). The steps are accumulated per segment. A segment whose
own estimate fails takes the step interpolated from the other segments (``filled``); a timepoint where every
segment fails is a row with an error and the chain continues from the last registrable frame.

The displacement of a nucleus. The cumulative displacements D_s(t) of the segments are interpolated linearly
along the axis between the segment centres (the median axis position of each segment's origin nuclei) and held
constant beyond the outermost centres. A nucleus at axis position a has the displacement D(t, a), and its
registered coordinate is its raw centroid minus D(t, a). Since a is the nucleus's position on the origin's axis
AFTER registration, it is found by a fixed-point iteration (a few steps: the displacement varies slowly along
the axis). Nothing is resampled; the per-timepoint tables and label images stay in their raw frame.

What it means for the tracks. The registered frame of the piecewise registration follows the GFP tissue itself,
so linking in it sees nearly stationary nuclei, which is its purpose; but a velocity measured in that frame is
the motion relative to the local flow, not the flow. The flow is the registration itself (the cumulative
displacement per segment over time), or the tracks' ``axis_velocity_raw_um_per_min`` (the raw, microscope
frame), which does not depend on the registration.

Writes ``<stem>__register_piecewise.csv`` and ``.parquet`` (schema.SERIES_TABLES["register_piecewise"], one row
per timepoint per segment) and ``<stem>__register_piecewise_report.json``. The rigid ``<stem>__register_shifts``
table is not touched. ``germquant track-link --registration piecewise`` reads the table (load_piecewise).

Opt-in (docs/TIMELAPSE.md "Seeding the piecewise registration", "Per-gap caps and the acquisition record"; every
default reproduces the registration above number for number):

* ``seed="rigid"`` (key ``timelapse.register.piecewise_seed``) starts every gap from the rigid step of that gap, the
  difference of the cumulative shifts of ``<stem>__register_shifts.csv`` (or ``seed_table``). The moving crop is the
  reference crop moved by the seed rounded to whole reference voxels in y and x (z is never seeded: every plane and
  the circular z estimate stay), so a segment is followed across a stage re-centring that one crop box for both
  volumes cannot follow. The estimate is the correlation of the two crops plus that whole-voxel offset; the cap
  applies to the departure from the seed, and when nothing is capped the step is the estimate itself. A seed below
  half a reference voxel moves no crop, so there the estimate is the unseeded one exactly. A seed table whose
  ``series_id`` names another series is refused; a gap whose rigid step the seed table capped is seeded with that
  shortened step and says so (``seed_source`` "rigid, capped in the seed table", a warning with the command to run).
* ``max_step_at`` and ``prior`` set the cap of a gap and flag it from the acquisition record as for the rigid
  registration (timelapse.priors); a seeded step counts only the recorded z move in its cap, and a seed that does not
  fit the recorded stage step is flagged ``seed_step_mismatch``.
* ``axis_guard`` (key ``timelapse.register.piecewise_axis_guard``) replaces a folded origin axis by the straight PC1
  axis of the origin's germline (guard_origin_axis), and warns when PC1 explains under PC1_MIN_VARIANCE_FRACTION of
  their yx variance.
"""
from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .. import __version__, provenance, schema
from ..fsutil import long_path
from .layout import _truthy
from .register import (
    CENTROID_COLS,
    METHODS,
    NORMALIZATIONS,
    _clean,
    _f,
    _l,
    _prior_option,
    _series_meta,
    cap_step,
    collect_reference_volumes,
    discover_series,
    displacement,
    nn_median,
)
from .register import DEFAULTS as RIGID_DEFAULTS
from .series import BRIGHTFIELD_ROLE, write_series_table
from .track import cfg_get, fixed_axis_model

log = logging.getLogger(__name__)

TABLE = "register_piecewise"
REPORT_SUFFIX = "__register_piecewise_report.json"
# code defaults of the timelapse.register.piecewise_* keys (read with cfg_get; no config file names them)
DEFAULTS = {"segments": 4, "margin_um": 6.0, "taper": 0.25}
KEYS = {"segments": "timelapse.register.piecewise_segments", "margin_um": "timelapse.register.piecewise_margin_um",
        "taper": "timelapse.register.piecewise_taper"}
MIN_NUCLEI_PER_SEGMENT = 5          # a segment needs at least this many origin nuclei
MIN_CROP_VOXELS = 4                 # a crop narrower than this many reference voxels in y or x is not estimated
FIXED_POINT_ITER = 25               # iterations of the registered axis position (it converges in a few)
FIXED_POINT_TOL_UM = 1e-6
# the opt-in keys (module docstring) and their code defaults
SEED_MODES = ("none", "rigid")
OPTION_DEFAULTS = {"seed": "none", "axis_guard": False, "axis_fold_ratio": 4.0}
OPTION_KEYS = {"seed": "timelapse.register.piecewise_seed", "axis_guard": "timelapse.register.piecewise_axis_guard",
               "axis_fold_ratio": "timelapse.register.piecewise_axis_fold_ratio"}
RIGID_TABLE_SUFFIX = "__register_shifts.csv"
RIGID_REPORT_SUFFIX = "__register_report.json"
SEED_COLUMNS = ("seed_dz_um", "seed_dy_um", "seed_dx_um", "seed_shift_y_vox", "seed_shift_x_vox", "mov_crop_y0",
                "mov_crop_y1", "mov_crop_x0", "mov_crop_x1", "seed_source", "cap_basis")


class PiecewiseRegistrationMissing(FileNotFoundError):
    """``<stem>__register_piecewise.csv`` is not in the series folder: run ``germquant register --mode piecewise``."""


def piecewise_settings(cfg=None, *, segments=None, margin_um=None, taper=None) -> dict:
    """``{"segments", "margin_um", "taper"}``: an explicit value wins, else the timelapse.register.piecewise_* key
    of ``cfg`` (a Config, a dict or None), else the code default. ValueError for a value out of range."""
    out = {}
    for name, given in (("segments", segments), ("margin_um", margin_um), ("taper", taper)):
        out[name] = given if given is not None else cfg_get(cfg, KEYS[name], DEFAULTS[name])
    try:
        k = int(out["segments"])
    except (TypeError, ValueError) as e:
        raise ValueError(f"{KEYS['segments']} must be a positive integer, got {out['segments']!r}") from e
    if isinstance(out["segments"], bool) or k < 1:
        raise ValueError(f"{KEYS['segments']} must be a positive integer, got {out['segments']!r}")
    m, tap = float(out["margin_um"]), float(out["taper"])
    if not m >= 0:
        raise ValueError(f"{KEYS['margin_um']} must be >= 0, got {out['margin_um']!r}")
    if not 0.0 <= tap <= 1.0:
        raise ValueError(f"{KEYS['taper']} must be a fraction between 0 and 1, got {out['taper']!r}")
    return {"segments": k, "margin_um": m, "taper": tap}


def piecewise_options(cfg=None, *, seed=None, axis_guard=None, axis_fold_ratio=None) -> dict:
    """``{"seed", "axis_guard", "axis_fold_ratio"}`` of the opt-in keys (module docstring): an explicit value wins,
    else the timelapse.register.piecewise_* key of ``cfg``, else the code default (no seed, no guard, ratio 4).
    ValueError for a value out of range."""
    out = {}
    for name, given in (("seed", seed), ("axis_guard", axis_guard), ("axis_fold_ratio", axis_fold_ratio)):
        out[name] = given if given is not None else cfg_get(cfg, OPTION_KEYS[name], OPTION_DEFAULTS[name])
    mode = str(out["seed"]).strip().lower()
    if mode not in SEED_MODES:
        raise ValueError(f"{OPTION_KEYS['seed']} (--seed) must be one of {SEED_MODES}, got {out['seed']!r}")
    try:
        ratio = float(out["axis_fold_ratio"])
    except (TypeError, ValueError) as e:
        raise ValueError(f"{OPTION_KEYS['axis_fold_ratio']} must be a number, got {out['axis_fold_ratio']!r}") from e
    if not ratio > 0:
        raise ValueError(f"{OPTION_KEYS['axis_fold_ratio']} must be > 0, got {out['axis_fold_ratio']!r}")
    return {"seed": mode, "axis_guard": _truthy(out["axis_guard"]), "axis_fold_ratio": ratio}


def load_seed_table(series_dir: Path, stem: str, seed_table=None) -> tuple[dict[int, np.ndarray], dict]:
    """The rigid cumulative shift per timepoint (``{t_index: (dz, dy, dx) um}``, finite rows only) of the seed table:
    ``seed_table``, else the folder's ``<stem>__register_shifts.csv``; and its record for the report (path, sha256,
    the parameters of its registration report when one sits next to it, ``capped_t_indices``: the timepoints whose
    rigid step that registration capped, so their seeds are shortened steps). ValueError naming the command to run
    when there is none, and for a table whose ``series_id`` column names another series."""
    p = Path(seed_table) if seed_table is not None else Path(series_dir) / f"{stem}{RIGID_TABLE_SUFFIX}"
    if not p.is_file():
        if seed_table is None:
            raise ValueError(f"--seed rigid starts every gap from the rigid step of {p.name}, which is not in "
                             f"{series_dir}: run `germquant register {series_dir}` first (with the cap the movie "
                             f"needs, for example --max-step-um 15 --prior stage), or give a rigid table with --seed-table")
        raise ValueError(f"--seed-table {p}: no such file")
    df = pd.read_csv(long_path(p), low_memory=False)
    cols = next((c for c in (("dz_um", "dy_um", "dx_um"), ("cum_shift_z_um", "cum_shift_y_um", "cum_shift_x_um"))
                 if set(c) <= set(df.columns)), None)
    if cols is None or "t_index" not in df.columns:
        raise ValueError(f"{p.name}: not a register_shifts table (no t_index and dz_um, dy_um, dx_um columns)")
    if "series_id" in df.columns:
        others = sorted({str(v) for v in df["series_id"].dropna()} - {str(stem)})
        if others:
            raise ValueError(f"--seed-table {p.name} is the rigid table of {', '.join(others)}, not of {stem}: seed "
                             f"from this series' own rigid table (`germquant register {series_dir}`)")
    vals = df[list(cols)].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    seeds = {}
    for t, v in zip(pd.to_numeric(df["t_index"], errors="coerce"), vals):
        if pd.notna(t) and np.all(np.isfinite(v)):
            seeds[int(t)] = v
    capped = []
    if "capped" in df.columns:
        capped = sorted({int(t) for t, c in zip(pd.to_numeric(df["t_index"], errors="coerce"), df["capped"])
                         if pd.notna(t) and _truthy(c)})
    params = None
    if p.name.endswith(RIGID_TABLE_SUFFIX):
        rp = p.with_name(p.name[: -len(RIGID_TABLE_SUFFIX)] + RIGID_REPORT_SUFFIX)
        try:
            with open(long_path(rp), encoding="utf-8") as fh:
                par = json.load(fh).get("parameters") or {}
            params = {k: par.get(k) for k in ("reference", "max_step_um", "max_step_at", "prior", "taper",
                                              "channel_role", "normalization", "stride")}
        except (OSError, ValueError, AttributeError):
            params = None
    with open(long_path(p), "rb") as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    return seeds, {"table": str(p), "sha256": sha, "rigid_parameters": params, "capped_t_indices": capped}


# ---- the axis, the segments and the interpolation ------------------------------------------------------------

def germline_points(df: pd.DataFrame | None, t_index: int) -> pd.DataFrame | None:
    """The germline nuclei of one nuclei table as the frame track.fixed_axis_model reads: ``t_index``,
    ``nucleus_id``, ``zr_um``, ``yr_um``, ``xr_um`` (the raw centroid; at the origin the displacement is zero),
    ``in_germline`` and ``axis_position_um`` (NaN without one). The in_germline rows when the column marks any,
    else every row; rows without a finite centroid are dropped. None for a missing or empty table."""
    if df is None or df.empty or any(c not in df.columns for c in CENTROID_COLS):
        return None
    if "in_germline" in df.columns:
        m = df["in_germline"].map(_truthy).to_numpy(bool)
        if m.any():
            df = df[m]
    c = df[list(CENTROID_COLS)].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    ok = np.isfinite(c).all(axis=1)
    if not ok.any():
        return None
    ap = (pd.to_numeric(df["axis_position_um"], errors="coerce").to_numpy(float) if "axis_position_um" in df.columns
          else np.full(len(df), np.nan))
    nid = (pd.to_numeric(df["nucleus_id"], errors="coerce").to_numpy(float) if "nucleus_id" in df.columns
           else np.arange(1, len(df) + 1, dtype=float))
    return pd.DataFrame({"t_index": int(t_index), "nucleus_id": nid[ok], "zr_um": c[ok, 0], "yr_um": c[ok, 1],
                         "xr_um": c[ok, 2], "in_germline": True, "axis_position_um": ap[ok]})


def read_nuclei_table(folder: str | Path, image_id: str) -> pd.DataFrame | None:
    f = Path(folder) / f"{image_id}__nuclei.csv"
    if not f.is_file():
        return None
    try:
        return pd.read_csv(long_path(f), low_memory=False)
    except Exception as e:  # noqa: BLE001 - a damaged table is a missing quality number
        log.warning("%s: unreadable nuclei table (%s)", image_id, e)
        return None


def origin_axis(points: pd.DataFrame) -> tuple:
    """``(model, source)``: the fixed-axis model of the origin's germline points (track.fixed_axis_model) and
    which rule built it (``axis_position_um`` or ``centerline``). ValueError when it cannot be built."""
    model = fixed_axis_model(points)
    if model is None:
        raise ValueError("the origin timepoint has fewer than 5 germline nuclei with a centroid; the piecewise "
                         "registration needs them to place its segments along the gonad")
    n_ax = int(np.isfinite(points["axis_position_um"].to_numpy(float)).sum())
    return model, ("axis_position_um" if n_ax >= 5 else "centerline")


# ---- the folded-axis guard (docs/TIMELAPSE.md "The folded-axis guard"; opt-in) ---------------------------------------

PC1_COLUMNS = {"vy": "axis_pc1_vy", "vx": "axis_pc1_vx", "cy_um": "axis_pc1_cy_um", "cx_um": "axis_pc1_cx_um",
               "offset_um": "axis_pc1_offset_um"}
PC1_MIN_VARIANCE_FRACTION = 0.85    # a tripped guard warns when PC1 explains less of the origin's yx variance


def _pc1_direction(yx: np.ndarray) -> np.ndarray:
    """The first principal direction (unit, (y, x)) of the centred points."""
    return np.linalg.svd(yx - yx.mean(axis=0), full_matrices=False)[2][0]


def pc1_variance_fraction(yx: np.ndarray) -> float:
    """The share of the points' yx variance along their first principal direction (1 for points on a line, 0.5 for
    an isotropic cloud): how well one straight axis can describe the origin's germline. The spread ratio cannot tell,
    since PC1 is straight by construction and always reads well on it."""
    s = np.linalg.svd(np.asarray(yx, dtype=float) - np.asarray(yx, dtype=float).mean(axis=0), compute_uv=False)
    total = float(np.sum(s ** 2))
    return float(s[0] ** 2 / total) if total > 0 else 1.0


def spread_ratios(points: pd.DataFrame, a0, k: int) -> dict:
    """The folded-axis measure: for each of the K equal-count segments of the axis positions ``a0`` of the origin's
    germline ``points``, the 5th to 95th percentile range of its nuclei's projection on the first principal direction
    of their yx centroids, over the segment's axis span. A sound axis, straight or curved, reads about 1 (a segment on
    a bend spreads less along the global direction than along the axis); a fold, two limbs side by side on the axis,
    or a collapsed stretch reads high. ``{"ratios", "spread_ratio" (the largest), "spans_um", "extents_um"}``."""
    yx = points[["yr_um", "xr_um"]].to_numpy(float)
    pc = (yx - yx.mean(axis=0)) @ _pc1_direction(yx)
    edges, _centres, seg = segment_bounds(a0, k)
    spans = np.diff(edges)
    ext = np.array([float(np.percentile(pc[seg == s], 95) - np.percentile(pc[seg == s], 5)) if np.any(seg == s)
                    else np.nan for s in range(int(k))])
    ratios = ext / np.maximum(spans, 1e-6)
    return {"ratios": ratios, "spread_ratio": float(np.nanmax(ratios)), "spans_um": spans, "extents_um": ext}


def pc1_model(params: dict):
    """The straight PC1 axis as a callable on (N, 3) points (z, y, x um): the projection of (y, x) minus the centre on
    the direction, minus the offset (so the origin's smallest value is 0). Built from the stored numbers only, so
    the registration and load_piecewise evaluate the same function."""
    v = np.array([float(params["vy"]), float(params["vx"])])
    mu = np.array([float(params["cy_um"]), float(params["cx_um"])])
    off = float(params["offset_um"])

    def model(x):
        x = np.asarray(x, dtype=float).reshape(-1, 3)
        return (x[:, 1:] - mu) @ v - off

    return model


def pc1_axis(points: pd.DataFrame) -> tuple:
    """``(model, params)``: the straight principal axis of the origin's germline yx centroids, the guard's fallback.
    Its sign makes the Spearman correlation with ``axis_position_um`` non-negative where the origin has five or
    more positions (segment 0 then sits as near the axis stage's distal end as the fold allows), else the SVD's sign
    is kept; it is zeroed at the smallest origin value."""
    from scipy.stats import spearmanr

    yx = points[["yr_um", "xr_um"]].to_numpy(float)
    mu = yx.mean(axis=0)
    v = _pc1_direction(yx)
    proj = (yx - mu) @ v
    ap = points["axis_position_um"].to_numpy(float) if "axis_position_um" in points.columns else np.full(len(yx), np.nan)
    ok = np.isfinite(ap)
    rule = "svd (fewer than 5 axis positions)"
    if int(ok.sum()) >= 5:
        rule = "spearman with axis_position_um, non-negative"
        rho = float(spearmanr(proj[ok], ap[ok])[0])
        if np.isfinite(rho) and rho < 0:
            v = -v
    params = {"vy": float(v[0]), "vx": float(v[1]), "cy_um": float(mu[0]), "cx_um": float(mu[1])}
    params["offset_um"] = float(np.min((yx - mu) @ np.array([params["vy"], params["vx"]])))
    return pc1_model(params), {**params, "sign_rule": rule}


def guard_origin_axis(points: pd.DataFrame, c0, a0, axis, source: str, k: int, threshold: float) -> tuple:
    """The guard on the origin axis: ``(axis, source, a0, block)``, the axis to use (the default one, or PC1 when the
    spread ratio is above ``threshold``), its positions of the origin nuclei ``c0`` and the report block. A default
    axis whose positions cannot be split into K segments at all counts as folded. The block also records
    ``pc1_variance_fraction`` (pc1_variance_fraction of the origin's yx centroids) and, when the guard tripped and
    that is under PC1_MIN_VARIANCE_FRACTION, a ``warning``: the straight axis it fell back on is a poor fit."""
    block = {"threshold": float(threshold), "axis_before": source}
    try:
        sr = spread_ratios(points, a0, k)
        block.update(spread_ratio_by_segment=[float(v) for v in sr["ratios"]], spread_ratio=sr["spread_ratio"],
                     segment_spans_before_um=[float(v) for v in sr["spans_um"]])
        tripped = bool(sr["spread_ratio"] > threshold)
    except ValueError as e:                 # the positions collapse into ties: as folded as an axis gets
        block.update(spread_ratio_by_segment=None, spread_ratio=None, segment_spans_before_um=None,
                     why=f"the default axis does not separate into {k} segments: {e}")
        tripped = True
    frac = pc1_variance_fraction(points[["yr_um", "xr_um"]].to_numpy(float))
    block.update(tripped=tripped, pc1_variance_fraction=frac, pc1_min_variance_fraction=PC1_MIN_VARIANCE_FRACTION,
                 warning=None)
    if not tripped:
        block.update(axis_used=source, pc1=None)
        return axis, source, a0, block
    model, params = pc1_axis(points)
    block.update(axis_used="pc1", pc1=params)
    if frac < PC1_MIN_VARIANCE_FRACTION:
        block["warning"] = (f"the axis guard fell back on PC1, the straight principal axis, but PC1 explains only "
                            f"{frac:.2f} of the origin's yx variance (under {PC1_MIN_VARIANCE_FRACTION:g}): the germline "
                            f"is not straight enough for one straight axis; check the origin's axis, or use fewer "
                            f"segments")
    return model, "pc1", np.asarray(model(c0), dtype=float), block


def _pc1_from_table(path: Path) -> dict:
    """The stored PC1 axis of a register_piecewise table, its numbers read back exactly (round-trip parser)."""
    cols = list(PC1_COLUMNS.values())
    df = pd.read_csv(long_path(path), usecols=cols, nrows=1, float_precision="round_trip")
    if df.empty or not np.isfinite(df[cols].to_numpy(float)).all():
        raise ValueError(f"{path.name}: axis_source is pc1 but the axis columns {', '.join(cols)} are empty")
    return {k: float(df[c].iloc[0]) for k, c in PC1_COLUMNS.items()}


def assign_segments(a, edges) -> np.ndarray:
    """The segment (0 .. K-1) of each axis position for the segment edges (K + 1 values); positions before the
    first edge belong to segment 0 and after the last to segment K-1; NaN gives -1."""
    a = np.asarray(a, dtype=float)
    seg = np.searchsorted(np.asarray(edges, dtype=float)[1:-1], a, side="right").astype(int)
    seg[~np.isfinite(a)] = -1
    return seg


def segment_bounds(a0, k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """K segments of equal nucleus count along the axis: ``(edges, centres, seg)`` with the K + 1 quantile
    edges of the finite positions ``a0``, the median position of each segment (strictly increasing) and the
    segment of every entry of ``a0``. ValueError with fewer than MIN_NUCLEI_PER_SEGMENT positions per segment."""
    a0 = np.asarray(a0, dtype=float)
    fin = a0[np.isfinite(a0)]
    if len(fin) < MIN_NUCLEI_PER_SEGMENT * int(k):
        raise ValueError(f"{len(fin)} germline nuclei with an axis position at the origin: too few for {k} segments "
                         f"of at least {MIN_NUCLEI_PER_SEGMENT} (use fewer segments)")
    edges = np.quantile(fin, np.linspace(0.0, 1.0, int(k) + 1))
    seg = assign_segments(a0, edges)
    centres = np.array([float(np.median(a0[seg == s])) if np.any(seg == s) else np.nan for s in range(int(k))])
    if not np.all(np.isfinite(centres)) or (len(centres) > 1 and np.any(np.diff(centres) <= 0)):
        raise ValueError(f"the origin's axis positions do not separate into {k} segments (ties?); use fewer segments")
    return edges, centres, seg


def interp_displacement(cum, centres, a) -> np.ndarray:
    """D(t, a): the per-segment cumulative displacements ``cum`` (K x 3, um) interpolated linearly along the axis
    between the segment ``centres`` and held constant beyond the outermost ones, at the positions ``a`` (N x 3
    result; NaN for a NaN position)."""
    cum = np.asarray(cum, dtype=float).reshape(-1, 3)
    c = np.asarray(centres, dtype=float).ravel()
    a = np.asarray(a, dtype=float).ravel()
    out = np.full((len(a), 3), np.nan)
    ok = np.isfinite(a)
    if not ok.any():
        return out
    if len(c) == 1:
        out[ok] = cum[0]
        return out
    for ax in range(3):
        out[ok, ax] = np.interp(a[ok], c, cum[:, ax])
    return out


def register_points(raw, cum, centres, axis, *, max_iter: int = FIXED_POINT_ITER,
                    tol_um: float = FIXED_POINT_TOL_UM) -> tuple[np.ndarray, np.ndarray]:
    """``(registered, a)`` for raw-frame points (N x 3, um) at a timepoint with the per-segment cumulative
    displacements ``cum``: registered = raw - D(t, a) with a = axis(registered), the point's position on the
    origin's axis after registration, solved by fixed-point iteration from the median segment displacement.
    Rows with a non-finite coordinate come back NaN."""
    raw = np.asarray(raw, dtype=float).reshape(-1, 3)
    cum = np.asarray(cum, dtype=float).reshape(-1, 3)
    reg = np.full(raw.shape, np.nan)
    a = np.full(len(raw), np.nan)
    ok = np.isfinite(raw).all(axis=1)
    if not ok.any():
        return reg, a
    x = raw[ok]
    if not np.any(cum):
        reg[ok] = x
        a[ok] = axis(x)
        return reg, a
    aa = np.asarray(axis(x - np.median(cum, axis=0)), dtype=float)
    for _ in range(int(max_iter)):
        a_new = np.asarray(axis(x - interp_displacement(cum, centres, aa)), dtype=float)
        done = bool(np.all(np.abs(a_new - aa) <= tol_um))
        aa = a_new
        if done:
            break
    reg[ok] = x - interp_displacement(cum, centres, aa)
    a[ok] = aa
    return reg, a


def crop_box(points_um, spacing, shape, margin_um: float) -> tuple[int, int, int, int]:
    """``(y0, y1, x0, x1)``: the reference-volume voxel range (stop exclusive, clipped to ``shape``) of the y and
    x extent of the points (um) widened by ``margin_um``. A reference voxel j sits at j * spacing."""
    p = np.asarray(points_um, dtype=float).reshape(-1, 3)
    p = p[np.isfinite(p).all(axis=1)]
    if not len(p):
        return 0, 0, 0, 0
    sy, sx = float(spacing[1]), float(spacing[2])
    y0 = int(np.floor((p[:, 1].min() - margin_um) / sy))
    y1 = int(np.ceil((p[:, 1].max() + margin_um) / sy)) + 1
    x0 = int(np.floor((p[:, 2].min() - margin_um) / sx))
    x1 = int(np.ceil((p[:, 2].max() + margin_um) / sx)) + 1
    Y, X = int(shape[1]), int(shape[2])
    return max(0, y0), min(Y, y1), max(0, x0), min(X, x1)


# ---- the command ----------------------------------------------------------------------------------------------

def register_series_piecewise(
    series_dir: str | Path,
    *,
    segments: int | None = None,
    margin_um: float | None = None,
    taper: float | None = None,
    upsample: int = RIGID_DEFAULTS["upsample"],
    max_step_um: float | None = RIGID_DEFAULTS["max_step_um"],
    channel: str | None = None,
    disagreement_um: float = RIGID_DEFAULTS["disagreement_um"],
    normalization: str = RIGID_DEFAULTS["normalization"],
    cfg=None,
    max_step_at: dict | None = None,
    prior: str | None = None,
    nd2: str | Path | None = None,
    seed: str | None = None,
    seed_table: str | Path | None = None,
    axis_guard: bool | None = None,
    axis_fold_ratio: float | None = None,
) -> dict:
    """Register a finished series folder piecewise along the gonad axis (module docstring) and write the
    register_piecewise table and its report next to the series record.

    ``segments`` (K), ``margin_um`` and ``taper`` default to the timelapse.register.piecewise_* keys of ``cfg``,
    else 4, 6.0 um and 0.25. ``upsample``, ``max_step_um`` (the cap on one segment's step, its direction kept),
    ``channel`` (the role whose volume is registered; the other kept role is the per-segment cross-check),
    ``disagreement_um`` and ``normalization`` mean what they mean for ``germquant register``. The reference frame
    is always the previous registrable timepoint.

    Opt-in (module docstring; piecewise_options for the keys): ``seed`` "rigid" seeds every gap with the rigid
    step of ``seed_table`` (default the folder's ``<stem>__register_shifts.csv``) and caps the departure from it;
    ``max_step_at``, ``prior`` and ``nd2`` as for register.register_series; ``axis_guard`` and
    ``axis_fold_ratio`` the folded-axis guard. With any of them the table gains its columns after the declared
    ones and the report its blocks; without them both are main's.

    Returns ``series_id``, ``series_dir``, ``table`` (declared columns), ``table_path`` (the CSV), ``report`` and
    ``report_path``. ReferenceVolumesMissing when the series kept no reference volumes; ValueError for a bad option
    or an origin without enough germline nuclei."""
    st = piecewise_settings(cfg, segments=segments, margin_um=margin_um, taper=taper)
    K, margin, tap = st["segments"], st["margin_um"], st["taper"]
    norm = "none" if normalization is None else str(normalization).strip().lower()
    if norm not in NORMALIZATIONS:
        raise ValueError(f"normalization must be one of {NORMALIZATIONS}, got {normalization!r}")
    method = METHODS[norm]
    up = int(upsample)
    if up < 1:
        raise ValueError(f"upsample must be a positive integer, got {upsample!r}")
    cap = float(max_step_um) if max_step_um is not None else 0.0
    thr = float(disagreement_um)
    opt = piecewise_options(cfg, seed=seed, axis_guard=axis_guard, axis_fold_ratio=axis_fold_ratio)
    seeding = opt["seed"] == "rigid"
    if seed_table is not None and not seeding:
        raise ValueError("seed_table is the rigid table of the seed: give it with seed='rigid' (--seed rigid "
                         "--seed-table PATH)")
    pri = _prior_option(prior, nd2)
    per_gap = bool(max_step_at) or pri is not None or seeding

    found = discover_series(series_dir)
    sdir, stem, record, tps = found["series_dir"], found["stem"], found["record"], found["timepoints"]
    vols = collect_reference_volumes(sdir, tps)
    stride, main_role, files, roles = vols["stride"], vols["role"], vols["files"], vols["roles"]
    spacing = np.asarray(vols["spacing_um"], dtype=float)
    primary = main_role if channel is None else str(channel).strip()
    if primary not in roles:
        raise ValueError(f"no reference volume for channel role {primary!r}; this series kept {sorted(roles)}")
    if primary != BRIGHTFIELD_ROLE and BRIGHTFIELD_ROLE in roles:
        cross: str | None = BRIGHTFIELD_ROLE
    elif primary == BRIGHTFIELD_ROLE and main_role in roles:
        cross = main_role
    else:
        cross = None
    if per_gap:
        from . import priors as PR

        msa = PR.check_max_step_at(max_step_at, [int(tp["t_index"]) for j, tp in enumerate(tps)
                                                 if primary in files[j]][1:])
        stage = (PR.load_stage_record(sdir, stem, record, nd2=nd2, t_indices=[int(tp["t_index"]) for tp in tps])
                 if pri == "stage" else None)
        z_range = (record or {}).get("z_range")
        gap_rows: dict[int, dict] = {}          # t_index -> the per-gap columns of its segment rows
        gap_list: list[dict] = []               # the registered gaps with the prior, for the report's prior block
        seeds, seed_info = load_seed_table(sdir, stem, seed_table) if seeding else ({}, None)
        seed_capped = set(seed_info["capped_t_indices"]) if seeding else set()     # shortened rigid steps
        seed_rows: dict[tuple[int, int], dict] = {}     # (t_index, segment) -> the seed columns of that row
        unseeded_gaps: list[dict] = []
        seed_warnings: list[str] = []
        prior_warnings: list[str] = []

    n = len(tps)
    with_vol = [i for i in range(n) if primary in files[i]]
    if not with_vol:
        raise ValueError(f"{sdir}: no timepoint has a {primary!r} reference volume")
    origin = with_vol[0]
    k_origin = int(tps[origin]["t_index"])
    pts0 = germline_points(read_nuclei_table(tps[origin]["folder"], tps[origin]["image_id"]), k_origin)
    if pts0 is None:
        raise ValueError(f"{tps[origin]['image_id']}: no nuclei table with germline centroids at the origin; the "
                         f"piecewise registration places its segments on the origin's germline nuclei")
    axis, axis_source = origin_axis(pts0)
    c0 = pts0[["zr_um", "yr_um", "xr_um"]].to_numpy(float)
    a0 = np.asarray(axis(c0), dtype=float)
    guard = None
    if opt["axis_guard"]:           # a folded origin axis is replaced by the straight PC1 axis (guard_origin_axis)
        axis, axis_source, a0, guard = guard_origin_axis(pts0, c0, a0, axis, axis_source, K, opt["axis_fold_ratio"])
    edges, centres, seg0 = segment_bounds(a0, K)
    n0 = np.array([int(np.sum(seg0 == s)) for s in range(K)])
    if guard is not None:
        guard["segment_spans_after_um"] = [float(v) for v in np.diff(edges)]

    cache: dict[tuple[int, str], np.ndarray] = {}

    def load(i: int, role: str) -> np.ndarray:
        key = (i, role)
        if key not in cache:
            cache[key] = np.load(long_path(files[i][role]))
        return cache[key]

    cents: dict[int, np.ndarray | None] = {}

    def centroids(i: int) -> np.ndarray | None:
        if i not in cents:
            p = germline_points(read_nuclei_table(tps[i]["folder"], tps[i]["image_id"]), int(tps[i]["t_index"]))
            cents[i] = None if p is None else p[["zr_um", "yr_um", "xr_um"]].to_numpy(float)
        return cents[i]

    cum: dict[int, np.ndarray] = {origin: np.zeros((K, 3))}
    cum_bf: dict[int, np.ndarray | None] = {}
    rows: list[dict] = []
    steps_rep: list[dict] = []
    warnings: list[str] = []
    if guard is not None and guard.get("warning"):
        warnings.append(guard["warning"])
        log.warning("%s: %s", stem, guard["warning"])
    errors: dict[str, str] = {}
    capped_list: list[dict] = []
    disagree_list: list[dict] = []
    filled_list: list[dict] = []
    last = origin

    for i, tp in enumerate(tps):
        k = int(tp["t_index"])
        if i < origin or primary not in files[i]:
            why = "no reference volume"
            if tp.get("status") not in (None, "done"):
                why += f" (timepoint {tp['status']})"
            rows += [_row(tp, s, K, edges, centres, n0, k_origin, primary, stride, up, method, tap, margin, error=why)
                     for s in range(K)]
            errors[str(k)] = why
            continue
        if i == origin:
            has_bf = bool(cross) and cross in files[i]
            cum_bf[i] = np.zeros((K, 3)) if has_bf else None
            for s in range(K):
                r = _row(tp, s, K, edges, centres, n0, k_origin, primary, stride, up, method, tap, margin)
                _fill(r, step=np.zeros(3), cum=np.zeros(3), bf_step=np.zeros(3) if has_bf else None,
                      bf_cum=np.zeros(3) if has_bf else None)
                y0, y1, x0, x1 = crop_box(c0[seg0 == s], spacing, load(i, primary).shape, margin)
                r.update(crop_y0=y0, crop_y1=y1, crop_x0=x0, crop_x1=x1)
                rows.append(r)
            continue
        ref_i = last
        ref_k = int(tps[ref_i]["t_index"])
        vr, vm = load(ref_i, primary), load(i, primary)
        if vr.shape != vm.shape:
            why = (f"ValueError: reference volume shapes differ: {vr.shape} vs {vm.shape} (was timelapse.reference_stride "
                   f"changed between runs of this series?)")
            rows += [_row(tp, s, K, edges, centres, n0, k_origin, primary, stride, up, method, tap, margin, error=why)
                     for s in range(K)]
            errors[str(k)] = why
            log.warning("t%03d: %s", k, why)
            continue
        # where the segments sat in the reference frame: the origin nuclei moved by the interpolated displacement
        pos_ref = c0 + interp_displacement(cum[ref_i], centres, a0)
        sd, sy, sx, seed_src = None, 0, 0, None     # the seed of this gap and its whole-voxel crop offset in y, x
        cap_k = cap
        if per_gap:
            if seeding:
                if ref_k in seeds and k in seeds:
                    sd, seed_src = seeds[k] - seeds[ref_k], "rigid"
                    sy, sx = int(np.rint(sd[1] / spacing[1])), int(np.rint(sd[2] / spacing[2]))
                    hit = [t for t in sorted(seed_capped) if ref_k < t <= k]
                    if hit:                 # the seed table clipped a step of this gap: the seed is not the move
                        seed_src = "rigid, capped in the seed table"
                        w = (f"t{ref_k}-t{k}: the seed table capped its rigid step into "
                             f"{', '.join(f't{t:03d}' for t in hit)}, so this gap's seed is a shortened step, not the "
                             f"move; register the rigid table with the cap the movie needs, for example `germquant "
                             f"register {sdir} --max-step-um 15 --prior stage`, then seed again")
                        seed_warnings.append(w)
                        warnings.append(w)
                        log.warning("%s", w)
                else:
                    seed_src = f"none: no rigid step for t{(k if k not in seeds else ref_k):03d}"
                    unseeded_gaps.append({"gap": f"t{ref_k}-t{k}", "why": seed_src})
            gp = PR.gap_prior(stage, ref_k, k, z_range) if stage is not None else None
        raw_steps = np.full((K, 3), np.nan)
        seg_err = [""] * K
        boxes = []
        for s in range(K):
            box = crop_box(pos_ref[seg0 == s], spacing, vr.shape, margin)
            if sy or sx:            # the reference box, clipped so that the box moved by the seed stays in the volume
                box = (max(box[0], -sy), min(box[1], vr.shape[1] - sy), max(box[2], -sx), min(box[3], vr.shape[2] - sx))
            boxes.append(box)
            y0, y1, x0, x1 = box
            if y1 - y0 < MIN_CROP_VOXELS or x1 - x0 < MIN_CROP_VOXELS:
                seg_err[s] = f"crop of {y1 - y0} x {x1 - x0} reference voxels (the segment left the field?)"
                continue
            try:
                if sy or sx:        # the moving crop follows the seed; the whole-voxel offset is added back
                    raw_steps[s] = (displacement(vr[:, y0:y1, x0:x1], vm[:, y0 + sy:y1 + sy, x0 + sx:x1 + sx], up, norm,
                                                 tap) * spacing + np.array([0.0, sy * spacing[1], sx * spacing[2]]))
                else:
                    raw_steps[s] = displacement(vr[:, y0:y1, x0:x1], vm[:, y0:y1, x0:x1], up, norm, tap) * spacing
            except Exception as e:  # noqa: BLE001 - one segment that cannot be estimated is filled from the others
                seg_err[s] = f"{type(e).__name__}: {e}"
        good = [s for s in range(K) if not seg_err[s]]
        if not good:
            why = "no segment could be estimated: " + "; ".join(f"segment {s}: {seg_err[s]}" for s in range(K))
            rows += [_row(tp, s, K, edges, centres, n0, k_origin, primary, stride, up, method, tap, margin, error=why)
                     for s in range(K)]
            errors[str(k)] = why
            log.warning("t%03d: %s", k, why)
            continue
        if per_gap:     # the gap's cap (priors.gap_cap); its focus term is what the segments' own z uses of it
            cap_k, cap_src = PR.gap_cap(k, cap, max_step_at=msa, gp=gp, seeded=sd is not None,
                                        est_dz_um=max(abs(float(raw_steps[s][0])) for s in good))
        steps = np.full((K, 3), np.nan)
        capped = [False] * K
        deps = np.full((K, 3), np.nan)              # seeded: the departure of each estimate from the seed
        for s in good:
            if sd is not None:      # the cap is on the departure; an uncapped step is the estimate itself
                sd_yx = np.array([0.0, sd[1], sd[2]])
                deps[s] = raw_steps[s] - sd_yx
                dep, capped[s] = cap_step(deps[s], cap_k)
                steps[s] = sd_yx + dep if capped[s] else raw_steps[s]
            else:
                steps[s], capped[s] = cap_step(raw_steps[s], cap_k)
        filled = [s for s in range(K) if s not in good]
        for s in filled:
            steps[s] = interp_displacement(steps[good], centres[good], [centres[s]])[0]
            filled_list.append({"t_index": k, "segment": s, "why": seg_err[s]})
            warnings.append(f"t{k:03d} segment {s}: {seg_err[s]}; step filled from the other segments")
        cum[i] = cum[ref_i] + steps
        # the cross-check: the same crops on the other kept channel, never capped
        bf_steps = np.full((K, 3), np.nan)
        if cross and cross in files[i] and cross in files[ref_i]:
            br, bm = load(ref_i, cross), load(i, cross)
            for s in range(K):
                y0, y1, x0, x1 = boxes[s]
                if br.shape != bm.shape or y1 - y0 < MIN_CROP_VOXELS or x1 - x0 < MIN_CROP_VOXELS:
                    continue
                try:
                    if sy or sx:    # the same moving box as the primary estimate
                        bf_steps[s] = (displacement(br[:, y0:y1, x0:x1], bm[:, y0 + sy:y1 + sy, x0 + sx:x1 + sx], up,
                                                    norm, tap) * spacing + np.array([0.0, sy * spacing[1], sx * spacing[2]]))
                    else:
                        bf_steps[s] = displacement(br[:, y0:y1, x0:x1], bm[:, y0:y1, x0:x1], up, norm, tap) * spacing
                except Exception as e:  # noqa: BLE001 - the cross-check is advisory
                    warnings.append(f"t{k:03d} segment {s}: {cross} cross-check failed: {type(e).__name__}: {e}")
            base = cum_bf.get(ref_i)
            base = cum[ref_i] if base is None else np.where(np.isfinite(base), base, cum[ref_i])
            cum_bf[i] = base + bf_steps
        else:
            cum_bf[i] = None
        # quality per segment: this timepoint's germline nuclei against the reference frame's, raw and registered
        c_mov, c_ref = centroids(i), centroids(ref_i)
        q = {s: (0, float("nan"), float("nan")) for s in range(K)}
        if c_mov is not None and c_ref is not None:
            reg_mov, a_mov = register_points(c_mov, cum[i], centres, axis)
            reg_ref, _a_ref = register_points(c_ref, cum[ref_i], centres, axis)
            seg_mov = assign_segments(a_mov, edges)
            for s in range(K):
                m = seg_mov == s
                q[s] = (int(m.sum()), nn_median(c_mov[m], c_ref), nn_median(reg_mov[m], reg_ref))
        seg_rep = []
        for s in range(K):
            r = _row(tp, s, K, edges, centres, n0, k_origin, primary, stride, up, method, tap, margin, error=seg_err[s])
            dis = float(np.linalg.norm(bf_steps[s] - raw_steps[s])) if np.all(np.isfinite(bf_steps[s])) else float("nan")
            if s in filled and np.all(np.isfinite(bf_steps[s])):
                dis = float(np.linalg.norm(bf_steps[s] - steps[s]))
            flag = bool(np.isfinite(dis) and dis > thr)
            bfc = None if cum_bf[i] is None else cum_bf[i][s]
            _fill(r, step=steps[s], cum=cum[i][s], capped=capped[s], filled=s in filled,
                  bf_step=bf_steps[s] if np.all(np.isfinite(bf_steps[s])) else None,
                  bf_cum=bfc if bfc is not None and np.all(np.isfinite(bfc)) else None,
                  disagreement=dis, flag=flag, n_nuc=q[s][0], q_before=q[s][1], q_after=q[s][2])
            r.update(reference_t_index=ref_k, crop_y0=boxes[s][0], crop_y1=boxes[s][1], crop_x0=boxes[s][2],
                     crop_x1=boxes[s][3])
            rows.append(r)
            if capped[s]:
                capped_list.append({"t_index": k, "segment": s})
            if flag:
                disagree_list.append({"t_index": k, "segment": s})
            seg_rep.append({"segment": s, "step_um": _l(steps[s]), "raw_step_um": _l(raw_steps[s]) if capped[s] else None,
                            "cumulative_um": _l(cum[i][s]), "capped": capped[s], "filled": s in filled,
                            "brightfield_step_um": _l(bf_steps[s]) if np.all(np.isfinite(bf_steps[s])) else None,
                            "disagreement_um": _f(dis), "disagreement_flag": flag, "n_nuclei": q[s][0],
                            "quality_before_um": _f(q[s][1]), "quality_after_um": _f(q[s][2]),
                            "crop_yx": [int(v) for v in boxes[s]], "error": seg_err[s] or None})
            if per_gap and seeding:
                b = boxes[s]
                mov = [int(b[0]) + sy, int(b[1]) + sy, int(b[2]) + sx, int(b[3]) + sx]
                seed_rows[(k, s)] = {"seed_dz_um": sd[0] if sd is not None else np.nan,
                                     "seed_dy_um": sd[1] if sd is not None else np.nan,
                                     "seed_dx_um": sd[2] if sd is not None else np.nan, "seed_shift_y_vox": sy,
                                     "seed_shift_x_vox": sx, "mov_crop_y0": mov[0], "mov_crop_y1": mov[1],
                                     "mov_crop_x0": mov[2], "mov_crop_x1": mov[3], "seed_source": seed_src,
                                     "cap_basis": "departure" if sd is not None else "step"}
                seg_rep[-1].update(seed_um=_l(sd) if sd is not None else None,
                                   departure_um=_l(deps[s]) if np.all(np.isfinite(deps[s])) else None, mov_crop_yx=mov)
        steps_rep.append({"t_index": k, "reference_t_index": ref_k, "segments": seg_rep})
        if per_gap:
            flags: list[str] = []
            if stage is not None:
                found_flags = {f for s in good for f in PR.gap_flags(raw_steps[s], gp, cap, rigid=False)}
                if sd is not None and PR.xy_mismatch(float(np.hypot(sd[1], sd[2])), gp, cap):
                    found_flags.add("seed_step_mismatch")       # the rigid seed does not fit the recorded stage step
                flags = [f for f in PR.FLAG_ORDER if f in found_flags]
                if sd is None and gp["recorded"] and gp["stage_xy_um"] > margin:
                    why = ("--seed rigid has no rigid step for this gap: register the rigid table over it" if seeding
                           else "use --seed rigid")
                    w = (f"t{ref_k}-t{k}: the stage moved {gp['stage_xy_um']:.1f} um, more than the {margin:g} um crop "
                         f"margin; the piecewise crops cannot follow it ({why})")
                    warnings.append(w)
                    prior_warnings.append(w)
                    log.warning("%s", w)
            gap_rows[k] = PR.gap_row(cap_k, cap_src, gp, flags, with_prior=stage is not None)
            steps_rep[-1].update(cap_um=_f(gap_rows[k]["cap_um"]), cap_source=cap_src)
            if seeding:
                steps_rep[-1]["seed_source"] = seed_src
            if stage is not None:
                steps_rep[-1]["prior"] = PR.gap_summary(gp, flags)
                gap_list.append({**gp, "flags": flags, "cap_um": _f(gap_rows[k]["cap_um"]), "cap_source": cap_src})
        for key in [key for key in cache if key[0] != i]:
            del cache[key]
        last = i

    df = pd.DataFrame(rows, columns=list(schema.SERIES_TABLES[TABLE]))
    for c in ("reference_t_index", "origin_t_index", "segment", "n_segments", "n_nuclei_t0", "n_nuclei",
              "crop_y0", "crop_y1", "crop_x0", "crop_x1"):
        df[c] = pd.array([None if pd.isna(v) else int(v) for v in df[c]], dtype="Int64")
    df["t_seconds"] = pd.to_numeric(df["t_seconds"], errors="coerce")
    for c in ("capped", "filled", "disagreement_flag"):
        df[c] = df[c].astype(bool)
    if per_gap:
        PR.add_gap_columns(df, gap_rows, with_prior=stage is not None)
        if seeding:
            _add_seed_columns(df, seed_rows)
    if guard is not None:           # the axis used, and PC1's numbers for load_piecewise to rebuild it exactly
        df["axis_source"] = [axis_source] * len(df)
        if axis_source == "pc1":
            for key, col in PC1_COLUMNS.items():
                df[col] = np.full(len(df), float(guard["pc1"][key]))
    table_path = write_series_table(df, TABLE, sdir, stem, _series_meta(sdir, stem, record))

    by_seg = []
    for s in range(K):
        b = [x["quality_before_um"] for st_ in steps_rep for x in st_["segments"] if x["segment"] == s
             and x["quality_before_um"] is not None]
        a = [x["quality_after_um"] for st_ in steps_rep for x in st_["segments"] if x["segment"] == s
             and x["quality_after_um"] is not None]
        d = [x["disagreement_um"] for st_ in steps_rep for x in st_["segments"] if x["segment"] == s
             and x["disagreement_um"] is not None]
        by_seg.append({"segment": s, "n_steps": len(b), "median_before_um": float(np.median(b)) if b else None,
                       "median_after_um": float(np.median(a)) if a else None,
                       "median_disagreement_um": float(np.median(d)) if d else None})
    all_b = [x["quality_before_um"] for st_ in steps_rep for x in st_["segments"] if x["quality_before_um"] is not None]
    all_a = [x["quality_after_um"] for st_ in steps_rep for x in st_["segments"] if x["quality_after_um"] is not None]
    report = {
        "series_id": stem, "series_dir": str(sdir), "mode": "piecewise", "method": method,
        "parameters": {"segments": K, "margin_um": margin, "taper": tap, "reference": "previous",
                       "upsample_factor": up, "normalization": norm, "max_step_um": cap if cap > 0 else None,
                       "disagreement_um": thr, "channel_role": primary, "cross_check_role": cross, "stride": stride,
                       "spacing_um": _l(spacing), "origin_t_index": k_origin, "axis_source": axis_source,
                       "interpolation": "linear between segment centres, constant beyond the outermost centres"},
        "segments": [{"segment": s, "axis_start_um": float(edges[s]), "axis_end_um": float(edges[s + 1]),
                      "axis_centre_um": float(centres[s]), "n_nuclei_t0": int(n0[s])} for s in range(K)],
        "n_timepoints": n, "n_registered": len(steps_rep), "capped": capped_list,
        "brightfield_disagreement": disagree_list, "filled": filled_list, "errors": errors, "warnings": warnings,
        "quality": {"by_segment": by_seg, "median_before_um": float(np.median(all_b)) if all_b else None,
                    "median_after_um": float(np.median(all_a)) if all_a else None},
        "cumulative_last_um": [_l(v) for v in cum[last]], "last_t_index": int(tps[last]["t_index"]),
        "steps": steps_rep,
        "pipeline_version": __version__, "git_sha": provenance.git_sha(), "run_timestamp": provenance.run_timestamp(),
    }
    if per_gap:
        PR.add_report_options(report, msa, pri, stage, gap_list, z_range=z_range, warnings=prior_warnings)
        if seeding:
            report["parameters"]["seed"] = "rigid"
            if seed_table is not None:
                report["parameters"]["seed_table"] = str(seed_table)
            n_seeded = sum(1 for st_ in steps_rep if str(st_.get("seed_source") or "").startswith("rigid"))
            report["seed"] = {**seed_info, "n_seeded_gaps": n_seeded, "unseeded_gaps": unseeded_gaps,
                              "warnings": seed_warnings,
                              "z": "recorded, not applied (every plane and the circular z estimate are kept)",
                              "cap_basis": "the departure from the seed (estimate minus (0, seed_y, seed_x))"}
    if guard is not None:
        report["parameters"]["axis_guard"] = True
        report["axis_guard"] = guard
    report_path = sdir / f"{stem}{REPORT_SUFFIX}"
    with open(long_path(report_path), "w", encoding="utf-8") as fh:
        json.dump(_clean(report), fh, indent=1)
    log.info("%s: piecewise registration, %d segment(s), %d timepoint(s), %d step(s) registered -> %s", stem, K, n,
             len(steps_rep), table_path.name)
    return {"series_id": stem, "series_dir": str(sdir), "table": df, "table_path": str(table_path),
            "report": report, "report_path": str(report_path)}


def _row(tp: dict, s: int, K: int, edges, centres, n0, origin_k: int, role: str, stride: int, up: int, method: str,
         taper: float, margin: float, error: str = "") -> dict:
    nan = float("nan")
    r = {c: nan for c in schema.SERIES_TABLES[TABLE]}
    r.update(t_index=int(tp["t_index"]), t_seconds=tp.get("t_seconds"), image_id=tp["image_id"],
             reference_t_index=None, origin_t_index=int(origin_k), channel_role=role, segment=int(s), n_segments=int(K),
             axis_start_um=float(edges[s]), axis_end_um=float(edges[s + 1]), axis_centre_um=float(centres[s]),
             n_nuclei_t0=int(n0[s]), capped=False, filled=False, disagreement_flag=False, n_nuclei=None,
             crop_y0=None, crop_y1=None, crop_x0=None, crop_x1=None, error=error or "", xy_stride=int(stride),
             upsample_factor=int(up), method=method, taper=float(taper), margin_um=float(margin))
    return r


def _add_seed_columns(df: pd.DataFrame, seed_rows: dict) -> None:
    """Append SEED_COLUMNS to the table in place (rows of the origin and error rows: NA and empty text)."""
    keys = list(zip(df["t_index"].tolist(), df["segment"].tolist()))
    for c in SEED_COLUMNS:
        vals = [seed_rows.get((int(t), int(s)), {}).get(c) for t, s in keys]
        if c in ("seed_source", "cap_basis"):
            df[c] = ["" if v is None else str(v) for v in vals]
        elif c.startswith(("seed_shift", "mov_crop")):
            df[c] = pd.array([None if v is None else int(v) for v in vals], dtype="Int64")
        else:
            df[c] = np.asarray([np.nan if v is None else float(v) for v in vals], dtype=float)


def _fill(row: dict, *, step, cum, capped=False, filled=False, bf_step=None, bf_cum=None,
          disagreement=float("nan"), flag=False, n_nuc=None, q_before=float("nan"), q_after=float("nan")) -> None:
    for ax, s_, c_ in zip("zyx", step, cum):
        row[f"step_d{ax}_um"] = float(s_)
        row[f"d{ax}_um"] = float(c_)
    if bf_step is not None:
        for ax, v in zip("zyx", bf_step):
            row[f"brightfield_step_d{ax}_um"] = float(v)
    if bf_cum is not None:
        for ax, v in zip("zyx", bf_cum):
            row[f"brightfield_d{ax}_um"] = float(v)
    row["brightfield_disagreement_um"] = float(disagreement)
    row["disagreement_flag"] = bool(flag)
    row["capped"] = bool(capped)
    row["filled"] = bool(filled)
    row["n_nuclei"] = n_nuc
    row["quality_before_um"] = float(q_before)
    row["quality_after_um"] = float(q_after)


# ---- reading it back (germquant track-link --registration piecewise) ------------------------------------------

class PiecewiseRegistration:
    """The piecewise registration of a series as track-link uses it: ``displacement(t, raw)`` gives the
    displacement (N x 3, um) of raw-frame points at timepoint t, ``register(t, raw)`` the registered points and
    their axis positions. A timepoint the table has no estimate for takes the cumulative displacement of the
    nearest earlier timepoint that has one (zero before the origin), with a warning."""

    def __init__(self, *, stem: str, origin_t: int, edges, centres, cum: dict[int, np.ndarray], axis,
                 axis_source: str, table_path: str | None = None, parameters: dict | None = None):
        self.stem = stem
        self.origin_t = int(origin_t)
        self.edges = np.asarray(edges, dtype=float)
        self.centres = np.asarray(centres, dtype=float)
        self.cum = {int(t): np.asarray(v, dtype=float).reshape(-1, 3) for t, v in cum.items()}
        self.axis = axis
        self.axis_source = axis_source
        self.table_path = table_path
        self.parameters = dict(parameters or {})
        self.substituted: dict[int, int | None] = {}

    @property
    def n_segments(self) -> int:
        return len(self.centres)

    def cumulative(self, t: int) -> np.ndarray:
        t = int(t)
        if t in self.cum:
            return self.cum[t]
        earlier = [k for k in self.cum if k < t]
        src = max(earlier) if earlier else None
        if t not in self.substituted:
            self.substituted[t] = src
            log.warning("%s: no piecewise estimate for t%03d; %s", self.stem, t,
                        f"the displacement of t{src:03d} is used" if src is not None else "zero displacement used")
        return self.cum[src] if src is not None else np.zeros((self.n_segments, 3))

    def register(self, t: int, raw) -> tuple[np.ndarray, np.ndarray]:
        return register_points(raw, self.cumulative(t), self.centres, self.axis)

    def displacement(self, t: int, raw) -> np.ndarray:
        raw = np.asarray(raw, dtype=float).reshape(-1, 3)
        reg, _a = self.register(t, raw)
        return raw - reg

    def summary(self) -> dict:
        return {"table": self.table_path, "n_segments": self.n_segments, "origin_t_index": self.origin_t,
                "segment_edges_um": [round(float(v), 4) for v in self.edges],
                "segment_centres_um": [round(float(v), 4) for v in self.centres], "axis_source": self.axis_source,
                "timepoints_with_estimate": sorted(self.cum),
                "timepoints_substituted": {str(k): v for k, v in sorted(self.substituted.items())},
                "parameters": self.parameters}


def load_piecewise(layout) -> PiecewiseRegistration:
    """The piecewise registration of a series folder (layout.SeriesLayout) from ``<stem>__register_piecewise.csv``,
    with the fixed-axis model rebuilt from the origin timepoint's nuclei table exactly as the registration built
    it; a table whose ``axis_source`` column says ``pc1`` (the folded-axis guard) carries its axis in its
    ``axis_pc1_*`` columns, and that axis is built from them instead (no refit). PiecewiseRegistrationMissing
    without the table; ValueError for a table this code cannot read."""
    p = layout.series_table_path(TABLE)
    if not p.is_file():
        raise PiecewiseRegistrationMissing(
            f"{p.name} not found in {layout.out_dir}: run `germquant register {layout.out_dir} --mode piecewise` first")
    df = pd.read_csv(long_path(p), low_memory=False)
    need = {"t_index", "origin_t_index", "segment", "n_segments", "axis_start_um", "axis_end_um", "axis_centre_um",
            "dz_um", "dy_um", "dx_um"}
    missing = sorted(need - set(df.columns))
    if missing:
        raise ValueError(f"{p.name}: missing column(s) {', '.join(missing)}")
    origin = int(pd.to_numeric(df["origin_t_index"], errors="coerce").dropna().iloc[0])
    K = int(pd.to_numeric(df["n_segments"], errors="coerce").dropna().iloc[0])
    seg_rows = df.drop_duplicates("segment").sort_values("segment")
    if len(seg_rows) != K:
        raise ValueError(f"{p.name}: {len(seg_rows)} segment(s) listed but n_segments is {K}")
    edges = np.concatenate([seg_rows["axis_start_um"].to_numpy(float), seg_rows["axis_end_um"].to_numpy(float)[-1:]])
    centres = seg_rows["axis_centre_um"].to_numpy(float)
    cum: dict[int, np.ndarray] = {}
    for t, g in df.groupby("t_index"):
        g = g.sort_values("segment")
        v = g[["dz_um", "dy_um", "dx_um"]].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        if len(g) == K and np.isfinite(v).all():
            cum[int(t)] = v
    if origin not in cum:
        cum[origin] = np.zeros((K, 3))
    if "axis_source" in df.columns and str(next(iter(df["axis_source"].dropna()), "")) == "pc1":
        axis, source = pc1_model(_pc1_from_table(p)), "pc1"
    else:
        nuc = read_nuclei_table(layout.sub(origin), layout.image_id(origin))
        pts0 = germline_points(nuc, origin)
        if pts0 is None:
            raise ValueError(f"{layout.image_id(origin)}: the origin's nuclei table is missing; the piecewise "
                             f"registration cannot rebuild its axis")
        axis, source = origin_axis(pts0)
    params = {}
    for c in ("channel_role", "method", "taper", "margin_um", "upsample_factor", "xy_stride"):
        if c in df.columns and len(df[c].dropna()):
            v = df[c].dropna().iloc[0]
            params[c] = v.item() if hasattr(v, "item") else v
    return PiecewiseRegistration(stem=layout.stem, origin_t=origin, edges=edges, centres=centres, cum=cum, axis=axis,
                                 axis_source=source, table_path=str(p), parameters=params)
