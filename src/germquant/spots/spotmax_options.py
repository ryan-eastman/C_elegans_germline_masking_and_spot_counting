"""The SpotMAX chain with the optional keys of a pack recipe (docs/DESIGN_system_packs.md section 4.7,
docs/BACKENDS.md; tranche 5 package P2a).

detect.detect_spots_spotmax branches here, at its top, only when one of detect.SPOTMAX_OPTIONAL_KEYS
(spot_radius_z_um, remove_hot_pixels, sharpen) is given. Its own body, the validated default chain, is not
edited, and this code lives in a module of its own (design rule R3: new code in new modules). The chain
repeats the default body's steps and calls the same helpers of the detect module (the flood guard
_cap_candidates, the z-split merge _merge_z_columns, _empty and _finish_tables), looked up on that module at
call time, so a fix there reaches both chains. tests/test_pack_knobs.py spies on every SpotMAX call and on
those helpers, and checks with the real SpotMAX that neutral option values give the default body's tables.
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from . import detect as _detect

log = logging.getLogger(__name__)


def spotmax_chain_with_options(
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
      test stay on ``pre``, and the flood guard ranks candidates by the raw image, as in the default body.

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
        return _detect._empty(labels, marker, "spotmax")

    zi = np.clip(df["z"].astype(int), 0, labels.shape[0] - 1)
    yi = np.clip(df["y"].astype(int), 0, labels.shape[1] - 1)
    xi = np.clip(df["x"].astype(int), 0, labels.shape[2] - 1)
    df = df.assign(nucleus_id=np.asarray(labels)[zi, yi, xi])
    df = df[df["nucleus_id"] > 0].reset_index(drop=True)
    df = _detect._cap_candidates(df, img, max_spot_candidates)

    df["effect_size"] = np.nan
    df["spot_mean_intensity"] = np.nan
    try:
        dfc = df.drop(columns=["nucleus_id"]).copy()
        dfc.index = pd.Index(df["nucleus_id"].to_numpy(), name="Cell_ID")
        out = P.spots_calc_features_and_filter(
            pre, radii_px, dfc, lab=labels, raw_image=img, zyx_voxel_size=tuple(sp),
            gop_filtering_thresholds=None, **feat_kw)
        df = feature_frame(out, df, labels, effect_size_metric)
    except Exception as e:  # noqa: BLE001 - features are an enrichment; detection counts still stand
        log.warning("spot feature computation failed (%s: %s); effect_size unavailable this run.",
                    type(e).__name__, e)

    if merge_z_columns:
        df = _detect._merge_z_columns(df, sp, z_merge_gap_um, pre, z_merge_valley_frac)
    if effect_size_min > 0:
        if df["effect_size"].notna().any():
            df = df[df["effect_size"].fillna(np.inf) >= effect_size_min].reset_index(drop=True)
        else:
            log.warning("effect_size_min=%.2f requested but effect sizes unavailable; NOT filtering.",
                        effect_size_min)
    return _detect._finish_tables(df, labels, sp, marker, "spotmax")


def feature_frame(out, df: pd.DataFrame, labels, effect_size_metric: str) -> pd.DataFrame:
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
