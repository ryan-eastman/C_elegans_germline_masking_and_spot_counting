"""RAD-51 (or other) spot quantification: ``detect_spots`` dispatches to a registered spot detection
backend (spots/backends.py, docs/BACKENDS.md); the default, ``spotmax``, is the validated SpotMAX chain
implemented here as ``detect_spots_spotmax`` (it replaced the v1 blob_log foci detector).

Per nucleus, detect spots in the spot channel *inside the nucleus mask* using SpotMAX's validated
chain: preprocess -> semantic segmentation (where spots can be) -> peak detection -> per-spot
features incl. effect-size (spot vs background). Returns tidy per-spot and per-nucleus tables.

SpotMAX runs HEADLESS here (no GUI, no interactive prompts). Detection/threshold/effect-size
parameters live in germquant config (`spots:`) and are calibrated against Imaris counts; the GUI
is for *viewing*, not counting. Validated chain: scripts/smoke_spotmax.py.
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from .backends import DEFAULT_BACKEND, PER_NUC_COLS, PER_SPOT_COLS, conform_result, get_backend

log = logging.getLogger(__name__)

_EFFECT = "spot_vs_backgr_effect_size_glass"
# the SpotMAX parameters detect_spots resolves and overlays on the section it hands to a backend
SPOTMAX_KEYS = (
    "spot_radius_um", "gauss_sigma_um", "thresholding_method", "effect_size_metric", "effect_size_min",
    "merge_z_columns", "z_merge_gap_um", "z_merge_valley_frac", "max_spot_candidates",
)
# the SpotMAX defaults detect_spots falls back to when neither the caller nor the config section names a value
SPOTMAX_DEFAULTS = dict(
    spot_radius_um=0.3, gauss_sigma_um=0.08, thresholding_method="threshold_otsu", effect_size_metric=_EFFECT,
    effect_size_min=0.0, merge_z_columns=True, z_merge_gap_um=0.8, z_merge_valley_frac=0.8, max_spot_candidates=30000,
)
# optional SpotMAX keys a pack recipe may set (docs/DESIGN_system_packs.md 4.7, docs/BACKENDS.md): the spotmax
# backend forwards one to detect_spots_spotmax ONLY when the raw spots: section holds it. SPOTMAX_DEFAULTS is
# not extended, so a default run's backend section and its SpotMAX calls are exactly today's.
#   spot_radius_z_um   the z radius: radii = (rz, r, r) / spacing instead of r / spacing
#   remove_hot_pixels  passed to SpotMAX's preprocess_image as do_remove_hot_pixels
#   sharpen            SpotMAX core's DoG sharpening: detection runs on the sharpened copy of the preprocessed
#                      image, and the feature step gets that copy as sharp_spots_image (so the effect size is
#                      on the sharpened scale, as in SpotMAX's own runs) while its image stays the unsharpened one
SPOTMAX_OPTIONAL_KEYS = ("spot_radius_z_um", "remove_hot_pixels", "sharpen")


def detect_spots(
    spot_img: np.ndarray,
    labels: np.ndarray,
    spacing: tuple[float, float, float],
    *,
    marker: str = "RAD-51",
    spot_radius_um: float | None = None,
    gauss_sigma_um: float | None = None,
    thresholding_method: str | None = None,
    effect_size_metric: str | None = None,
    effect_size_min: float | None = None,
    merge_z_columns: bool | None = None,
    z_merge_gap_um: float | None = None,
    z_merge_valley_frac: float | None = None,
    max_spot_candidates: int | None = None,
    backend: str = DEFAULT_BACKEND,
    cfg_section: dict | None = None,
    channels: dict | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (per_spot_df, per_nucleus_df) from the spot detection backend ``backend``: by default
    ``spotmax``, the validated SpotMAX chain of ``detect_spots_spotmax`` below, called with exactly
    these arguments (spots/backends.py lists the others; docs/BACKENDS.md documents them).

    The keyword arguments are the SpotMAX parameters. ``cfg_section`` is the raw ``spots:`` config
    mapping; the backend receives a copy of it with these parameters resolved (an explicit argument wins,
    then the section's own value, then the SpotMAX default), so a backend finds its own keys (``log.diameter_um``) next to the shared ones
    (``max_spot_candidates``). Both tables always have the columns PER_SPOT_COLS / PER_NUC_COLS and
    the ``detector`` column carries the backend name; an unknown name raises ValueError.

    ``channels`` maps further channel roles to (Z, Y, X) arrays for a backend that declares
    ``input_roles`` (net3d reads ``dna``); it reaches only such a backend, so spotmax and log are
    called exactly as before whatever is passed."""
    # precedence: an explicit keyword (not None) wins; else the section's own value; else the SpotMAX default.
    # A caller that passes the shared keys only in cfg_section is therefore not silently overridden.
    section = dict(cfg_section or {})
    explicit = dict(
        spot_radius_um=spot_radius_um, gauss_sigma_um=gauss_sigma_um, thresholding_method=thresholding_method,
        effect_size_metric=effect_size_metric, effect_size_min=effect_size_min, merge_z_columns=merge_z_columns,
        z_merge_gap_um=z_merge_gap_um, z_merge_valley_frac=z_merge_valley_frac,
        max_spot_candidates=max_spot_candidates,
    )
    for key, default in SPOTMAX_DEFAULTS.items():
        if explicit[key] is not None:
            section[key] = explicit[key]
        elif key not in section:
            section[key] = default
    fn = get_backend(backend)
    if getattr(fn, "input_roles", None):
        return conform_result(backend, fn(spot_img, labels, spacing, cfg_section=section, marker=marker,
                                          channels=dict(channels or {})))
    return conform_result(backend, fn(spot_img, labels, spacing, cfg_section=section, marker=marker))


def detect_spots_spotmax(
    spot_img: np.ndarray,
    labels: np.ndarray,
    spacing: tuple[float, float, float],
    *,
    marker: str = "RAD-51",
    spot_radius_um: float = 0.3,
    gauss_sigma_um: float = 0.08,
    thresholding_method: str = "threshold_otsu",
    effect_size_metric: str = _EFFECT,
    effect_size_min: float = 0.0,
    merge_z_columns: bool = True,
    z_merge_gap_um: float = 0.8,
    z_merge_valley_frac: float = 0.8,
    max_spot_candidates: int = 30000,
    spot_radius_z_um: float | None = None,
    remove_hot_pixels: bool | None = None,
    sharpen: bool | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The SpotMAX reference implementation, registered as the ``spotmax`` backend (``detect_spots``
    dispatches here by default; the body is the one the stage always ran, so its outputs are unchanged).

    ``spot_radius_z_um``, ``remove_hot_pixels`` and ``sharpen`` are the optional keys of a pack recipe
    (SPOTMAX_OPTIONAL_KEYS); None means absent, and only when one is given does the call take the branch at the
    top, ``_spotmax_chain_with_options``. Without them the body below runs exactly as it always has.
    Return (per_spot_df, per_nucleus_df). Spots are assigned to the nucleus whose mask they fall
    in; `effect_size_min` (>0) drops spots below that spot-vs-background effect size.

    `merge_z_columns` (default on) collapses the z-axis spot-splitting artifact: confocal axial PSF
    (~0.6-0.8 um) is far wider than the z spot-radius, so SpotMAX detects one focus as 2-3 stacked
    peaks at the SAME (x,y) voxel. Two stacked peaks are merged ONLY if they are within
    `z_merge_gap_um` in z AND there is no real intensity valley between them along z (the dip stays
    above `z_merge_valley_frac` * the dimmer peak) — i.e. one PSF-blurred focus, not two distinct
    foci. A genuine valley keeps them separate, so dense pachytene nuclei are not over-merged.
    Validated vs Imaris on 20251105_N2_HERM_001 (2033 -> ~Imaris 1222). effect_size_min is NOT the right
    knob for the over-count (the extras are real bright peaks of one focus, not dim noise).

    Speed: SpotMAX thins its peak candidates with a minimum-distance filter whose cost grows with the
    square of the candidate count (hours on a flooded gonad, before the candidate cap below is
    reached). The peak step therefore runs with germquant's exact replica of that filter installed
    (spots/peak_filter.py, docs/PERFORMANCE.md); the candidates and every output are unchanged."""
    if spot_radius_z_um is not None or remove_hot_pixels is not None or sharpen is not None:
        return _spotmax_chain_with_options(
            spot_img, labels, spacing, marker=marker, spot_radius_um=spot_radius_um, gauss_sigma_um=gauss_sigma_um,
            thresholding_method=thresholding_method, effect_size_metric=effect_size_metric,
            effect_size_min=effect_size_min, merge_z_columns=merge_z_columns, z_merge_gap_um=z_merge_gap_um,
            z_merge_valley_frac=z_merge_valley_frac, max_spot_candidates=max_spot_candidates,
            spot_radius_z_um=spot_radius_z_um, remove_hot_pixels=remove_hot_pixels, sharpen=sharpen)
    import spotmax.pipe as P

    from .peak_filter import fast_min_distance_filter

    sp = np.asarray(spacing, dtype=float)
    img = np.asarray(spot_img).astype("float32")
    radii_px = spot_radius_um / sp
    gauss = gauss_sigma_um / float(sp.min())

    pre = P.preprocess_image(img, lab=labels, gauss_sigma=gauss, use_gpu=False,
                             spots_zyx_radii_pxl=radii_px)
    segm = P.spots_semantic_segmentation(
        pre, lab=labels, spots_zyx_radii_pxl=radii_px, do_try_all_thresholds=False,
        return_only_segm=True, thresholding_method=thresholding_method, use_gpu=False)
    with fast_min_distance_filter():
        det = P.spot_detection(pre, spots_segmantic_segm=segm, spots_zyx_radii_pxl=radii_px,
                               lab=labels, return_df=True)
    df = det[0] if isinstance(det, tuple) else det
    if df is None or len(df) == 0:
        return _empty(labels, marker, "spotmax")

    zi = np.clip(df["z"].astype(int), 0, labels.shape[0] - 1)
    yi = np.clip(df["y"].astype(int), 0, labels.shape[1] - 1)
    xi = np.clip(df["x"].astype(int), 0, labels.shape[2] - 1)
    df = df.assign(nucleus_id=np.asarray(labels)[zi, yi, xi])
    df = df[df["nucleus_id"] > 0].reset_index(drop=True)

    # FLOOD GUARD: a high-signal or artefact image can propose hundreds of thousands of candidate
    # peaks (n2rad51_noHS_male_002 with its June labels: 479,901 into SpotMAX's minimum-distance
    # filter, 54,487 inside nuclei after it; docs/PERFORMANCE.md), and the per-spot feature step below
    # costs about 5 ms per spot on one core (on one syp-2 mutant gonad the stage once ran for 8 h
    # without finishing). Keep the brightest `max_spot_candidates` before the features so the stage
    # always finishes. A real RAD-51 count is far below this, so hitting the cap means the image is
    # flooded (artefact, bleed-through or over-bright) and its count is a floor, not a measurement:
    # _cap_candidates logs a warning with the candidate count. It is a log line only; the default
    # path adds no flag or column for it.
    df = _cap_candidates(df, img, max_spot_candidates)

    # Per-spot features (effect size + intensity). SpotMAX needs the FULL detection frame
    # (incl. the *_local columns) Cell_ID-indexed, and it regroups by Cell_ID — so use the
    # returned feature frame itself as the per-spot source rather than aligning positionally.
    df["effect_size"] = np.nan
    df["spot_mean_intensity"] = np.nan
    try:
        dfc = df.drop(columns=["nucleus_id"]).copy()
        dfc.index = pd.Index(df["nucleus_id"].to_numpy(), name="Cell_ID")
        out = P.spots_calc_features_and_filter(
            pre, radii_px, dfc, lab=labels, raw_image=img, zyx_voxel_size=tuple(sp),
            gop_filtering_thresholds=None)
        feat = pd.concat(out[1]) if isinstance(out, tuple) and isinstance(out[1], list) and out[1] else None
        if feat is not None and effect_size_metric in feat.columns and {"z", "y", "x"}.issubset(feat.columns):
            feat = feat.reset_index(drop=True)
            fz = np.clip(feat["z"].astype(int), 0, labels.shape[0] - 1)
            fy = np.clip(feat["y"].astype(int), 0, labels.shape[1] - 1)
            fx = np.clip(feat["x"].astype(int), 0, labels.shape[2] - 1)
            feat["nucleus_id"] = np.asarray(labels)[fz, fy, fx]
            feat = feat[feat["nucleus_id"] > 0].reset_index(drop=True)
            feat["effect_size"] = feat[effect_size_metric].astype(float)
            ic = next((c for c in feat.columns if c.startswith("spot_raw_mean")), None)
            feat["spot_mean_intensity"] = feat[ic].astype(float) if ic else np.nan
            df = feat  # has nucleus_id, z, y, x, effect_size, spot_mean_intensity
    except Exception as e:  # noqa: BLE001 - features are an enrichment; detection counts still stand
        log.warning("spot feature computation failed (%s: %s); effect_size unavailable this run.",
                    type(e).__name__, e)

    # collapse z-axis spot-splitting BEFORE filtering/aggregation (see docstring)
    if merge_z_columns:
        df = _merge_z_columns(df, sp, z_merge_gap_um, pre, z_merge_valley_frac)

    if effect_size_min > 0:
        if df["effect_size"].notna().any():
            df = df[df["effect_size"].fillna(np.inf) >= effect_size_min].reset_index(drop=True)
        else:
            log.warning("effect_size_min=%.2f requested but effect sizes unavailable; NOT filtering.",
                        effect_size_min)

    return _finish_tables(df, labels, sp, marker, "spotmax")


def _spotmax_chain_with_options(
    spot_img, labels, spacing, *, marker, spot_radius_um, gauss_sigma_um, thresholding_method, effect_size_metric,
    effect_size_min, merge_z_columns, z_merge_gap_um, z_merge_valley_frac, max_spot_candidates,
    spot_radius_z_um=None, remove_hot_pixels=None, sharpen=None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """detect_spots_spotmax with the optional keys of a pack recipe (SPOTMAX_OPTIONAL_KEYS; None = absent).
    The steps are those of the default body (preprocess, semantic segmentation, peak detection with the
    replica filter, nucleus assignment, flood guard, features, z-split merge, effect-size filter, tables);
    only these differ:

    * ``spot_radius_z_um``: radii_px = (rz, r, r) / spacing instead of r / spacing, for every SpotMAX call;
    * ``remove_hot_pixels``: passed to preprocess_image as do_remove_hot_pixels;
    * ``sharpen``: as SpotMAX core does it (core.py, sharpen_spots and spots_detection). The preprocessed image
      (hot pixels and smoothing) stays ``pre``; its sharpened copy (spotmax.filters.DoG_spots with the radii
      and the object labels) feeds the semantic segmentation and the peak detection and goes to
      spots_calc_features_and_filter as sharp_spots_image, while ``pre`` stays that call's image. SpotMAX
      then computes the spot-versus-background effect size on the sharpened image (spotmax/pipe.py in 1.3.1:
      "When comparing to the background we use the sharpened image"), so effect_size and effect_size_min
      are on the DoG scale, as in the authors' runs; the other intensity features and the z-merge valley
      test stay on ``pre``.

    With every option at its neutral value (spot_radius_z_um equal to spot_radius_um, remove_hot_pixels and
    sharpen false) the tables equal the default body's, which tests/test_pack_knobs.py checks."""
    import spotmax.pipe as P

    from .peak_filter import fast_min_distance_filter

    sp = np.asarray(spacing, dtype=float)
    img = np.asarray(spot_img).astype("float32")
    if spot_radius_z_um is not None:
        radii_px = np.array([float(spot_radius_z_um), spot_radius_um, spot_radius_um], dtype=float) / sp
    else:
        radii_px = spot_radius_um / sp
    gauss = gauss_sigma_um / float(sp.min())
    pre_kw = {} if remove_hot_pixels is None else {"do_remove_hot_pixels": bool(remove_hot_pixels)}

    pre = P.preprocess_image(img, lab=labels, gauss_sigma=gauss, use_gpu=False,
                             spots_zyx_radii_pxl=radii_px, **pre_kw)
    detect_img, feat_kw = pre, {}
    if sharpen:
        from spotmax import filters as smf

        detect_img = smf.DoG_spots(pre, radii_px, use_gpu=False, lab=labels)
        feat_kw["sharp_spots_image"] = detect_img
    segm = P.spots_semantic_segmentation(
        detect_img, lab=labels, spots_zyx_radii_pxl=radii_px, do_try_all_thresholds=False,
        return_only_segm=True, thresholding_method=thresholding_method, use_gpu=False)
    with fast_min_distance_filter():
        det = P.spot_detection(detect_img, spots_segmantic_segm=segm, spots_zyx_radii_pxl=radii_px,
                               lab=labels, return_df=True)
    df = det[0] if isinstance(det, tuple) else det
    if df is None or len(df) == 0:
        return _empty(labels, marker, "spotmax")

    zi = np.clip(df["z"].astype(int), 0, labels.shape[0] - 1)
    yi = np.clip(df["y"].astype(int), 0, labels.shape[1] - 1)
    xi = np.clip(df["x"].astype(int), 0, labels.shape[2] - 1)
    df = df.assign(nucleus_id=np.asarray(labels)[zi, yi, xi])
    df = df[df["nucleus_id"] > 0].reset_index(drop=True)
    df = _cap_candidates(df, img, max_spot_candidates)

    df["effect_size"] = np.nan
    df["spot_mean_intensity"] = np.nan
    try:
        dfc = df.drop(columns=["nucleus_id"]).copy()
        dfc.index = pd.Index(df["nucleus_id"].to_numpy(), name="Cell_ID")
        out = P.spots_calc_features_and_filter(
            pre, radii_px, dfc, lab=labels, raw_image=img, zyx_voxel_size=tuple(sp),
            gop_filtering_thresholds=None, **feat_kw)
        df = _feature_frame(out, df, labels, effect_size_metric)
    except Exception as e:  # noqa: BLE001 - features are an enrichment; detection counts still stand
        log.warning("spot feature computation failed (%s: %s); effect_size unavailable this run.",
                    type(e).__name__, e)

    if merge_z_columns:
        df = _merge_z_columns(df, sp, z_merge_gap_um, pre, z_merge_valley_frac)
    if effect_size_min > 0:
        if df["effect_size"].notna().any():
            df = df[df["effect_size"].fillna(np.inf) >= effect_size_min].reset_index(drop=True)
        else:
            log.warning("effect_size_min=%.2f requested but effect sizes unavailable; NOT filtering.",
                        effect_size_min)
    return _finish_tables(df, labels, sp, marker, "spotmax")


def _feature_frame(out, df: pd.DataFrame, labels, effect_size_metric: str) -> pd.DataFrame:
    """The per-spot frame from the result of spots_calc_features_and_filter, built as the default body of
    detect_spots_spotmax builds it (nucleus_id from the feature frame's own positions, effect_size from the
    metric column, spot_mean_intensity from the first spot_raw_mean column), or ``df`` unchanged when the
    result carries no usable feature frame."""
    feat = pd.concat(out[1]) if isinstance(out, tuple) and isinstance(out[1], list) and out[1] else None
    if feat is None or effect_size_metric not in feat.columns or not {"z", "y", "x"}.issubset(feat.columns):
        return df
    feat = feat.reset_index(drop=True)
    fz = np.clip(feat["z"].astype(int), 0, labels.shape[0] - 1)
    fy = np.clip(feat["y"].astype(int), 0, labels.shape[1] - 1)
    fx = np.clip(feat["x"].astype(int), 0, labels.shape[2] - 1)
    feat["nucleus_id"] = np.asarray(labels)[fz, fy, fx]
    feat = feat[feat["nucleus_id"] > 0].reset_index(drop=True)
    feat["effect_size"] = feat[effect_size_metric].astype(float)
    ic = next((c for c in feat.columns if c.startswith("spot_raw_mean")), None)
    feat["spot_mean_intensity"] = feat[ic].astype(float) if ic else np.nan
    return feat


def _cap_candidates(df: pd.DataFrame, image: np.ndarray, cap: int) -> pd.DataFrame:
    """Flood guard: if more than `cap` candidate spots, keep the brightest `cap` (by raw intensity at
    the peak voxel), in original order. Bounds the per-spot feature step so an artefact or
    over-bright image can not wedge the run for hours. A real RAD-51 count is far below `cap`, so
    triggering this means the image is flooded and its count is a floor, not a measurement. The
    warning carries the candidate count; it is a log line only (no flag, no column)."""
    if len(df) <= cap:
        return df
    z = np.clip(df["z"].astype(int), 0, image.shape[0] - 1)
    y = np.clip(df["y"].astype(int), 0, image.shape[1] - 1)
    x = np.clip(df["x"].astype(int), 0, image.shape[2] - 1)
    bright = np.asarray(image)[z, y, x]
    keep = np.sort(np.argsort(bright)[::-1][:cap])
    log.warning("FLOODED: %d spot candidates inside nuclei exceed the cap max_spot_candidates=%d; "
                "keeping the brightest %d and dropping %d. The image is likely over-bright or carries an "
                "artefact, so its spot count is a FLOOR, not a measurement.",
                len(df), cap, cap, len(df) - cap)
    return df.iloc[keep].reset_index(drop=True)


def _merge_z_columns(df: pd.DataFrame, spacing: np.ndarray, gap_um: float,
                     image: np.ndarray | None = None, valley_frac: float = 0.8) -> pd.DataFrame:
    """Collapse z-axis spot-splitting: spots sharing an (x,y) voxel are merged into one ONLY when a
    consecutive z-pair is within `gap_um` AND has no real intensity valley between them (the dip along
    z stays above `valley_frac` * the dimmer peak) — one PSF-blurred focus. A genuine valley (two
    bright blobs with a dip between) or a z-gap > gap_um starts a new focus, so distinct stacked foci
    in dense nuclei are preserved. `image` is the (smoothed) spot channel used for the valley test;
    without it, falls back to a gap-only rule. Keeps the brightest spot per merged cluster."""
    if df.empty:
        return df
    bright = "spot_mean_intensity" if df["spot_mean_intensity"].notna().any() else "effect_size"
    gap_vox = gap_um / float(spacing[0])
    nz = image.shape[0] if image is not None else None
    work = df.assign(_xp=np.round(df["x"].to_numpy()).astype(int),
                     _yp=np.round(df["y"].to_numpy()).astype(int))
    keep = []
    for (xp, yp), g in work.groupby(["_xp", "_yp"], sort=False):
        if len(g) == 1:
            keep.append(g.index[0])
            continue
        g = g.sort_values("z")
        zc = g["z"].to_numpy()
        zi = np.clip(np.round(zc).astype(int), 0, (nz - 1) if nz else 2**30)
        prof = None
        if image is not None and 0 <= yp < image.shape[1] and 0 <= xp < image.shape[2]:
            prof = image[:, yp, xp]
        cluster_id = np.zeros(len(g), dtype=int)
        for k in range(1, len(g)):
            same = (zc[k] - zc[k - 1]) <= gap_vox
            if same and prof is not None and zi[k] > zi[k - 1]:
                seg = prof[zi[k - 1]:zi[k] + 1]
                same = seg.min() >= valley_frac * min(prof[zi[k - 1]], prof[zi[k]])  # no real valley
            cluster_id[k] = cluster_id[k - 1] if same else cluster_id[k - 1] + 1
        for _, gc in g.groupby(cluster_id):
            keep.append(gc[bright].fillna(-np.inf).idxmax())
    return df.loc[keep].reset_index(drop=True)


def _finish_tables(df: pd.DataFrame, labels, sp: np.ndarray, marker: str, detector: str):
    """The two output tables from a frame of accepted spots with columns nucleus_id, z, y, x (voxel
    coordinates, a fresh RangeIndex), effect_size and spot_mean_intensity. Shared by every built-in
    backend; the construction is the one the SpotMAX path always used, with the detector name as the
    only variable."""
    per_spot = pd.DataFrame({
        "spot_id": np.arange(len(df)),
        "nucleus_id": df["nucleus_id"].astype(int),
        "marker": marker,
        "z_um": df["z"].to_numpy() * sp[0],
        "y_um": df["y"].to_numpy() * sp[1],
        "x_um": df["x"].to_numpy() * sp[2],
        "effect_size": df["effect_size"].to_numpy(),
        "spot_mean_intensity": df["spot_mean_intensity"].to_numpy(),
        "detector": detector,
    }, columns=PER_SPOT_COLS)

    # one row per segmented nucleus (0 spots is a real zero, not a dropped row)
    counts = per_spot.groupby("nucleus_id").size()
    means = per_spot.groupby("nucleus_id")["spot_mean_intensity"].mean()
    nuc_ids = [int(v) for v in np.unique(labels) if v != 0]
    per_nuc = pd.DataFrame({
        "nucleus_id": nuc_ids,
        "marker": marker,
        "n_spots": [int(counts.get(n, 0)) for n in nuc_ids],
        "mean_spot_intensity": [float(means.get(n, np.nan)) for n in nuc_ids],
        "detector": detector,
    }, columns=PER_NUC_COLS)
    return per_spot, per_nuc


def _empty(labels, marker, detector: str = "spotmax"):
    nuc_ids = [int(v) for v in np.unique(labels) if v != 0]
    per_nuc = pd.DataFrame({
        "nucleus_id": nuc_ids, "marker": marker, "n_spots": 0,
        "mean_spot_intensity": np.nan, "detector": detector,
    }, columns=PER_NUC_COLS)
    return pd.DataFrame(columns=PER_SPOT_COLS), per_nuc
