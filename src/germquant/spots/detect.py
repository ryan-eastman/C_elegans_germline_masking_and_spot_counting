"""RAD-51 (or other) spot quantification via SpotMAX — replaces the v1 blob_log foci detector.

Per nucleus, detect spots in the spot channel *inside the nucleus mask* using SpotMAX's validated
chain: preprocess -> semantic segmentation (where spots can be) -> peak detection -> per-spot
features incl. effect-size (spot vs background). Returns tidy per-spot and per-nucleus tables.

SpotMAX runs HEADLESS here (no GUI, no interactive prompts). Detection/threshold/effect-size
parameters live in germquant config (`spots:`) and are calibrated against Imaris counts — the GUI
is for *viewing*, not counting. Validated chain: scripts/smoke_spotmax.py.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

PER_SPOT_COLS = [
    "spot_id", "nucleus_id", "marker", "z_um", "y_um", "x_um",
    "effect_size", "spot_mean_intensity", "detector",
]
PER_NUC_COLS = ["nucleus_id", "marker", "n_spots", "mean_spot_intensity", "detector"]
_EFFECT = "spot_vs_backgr_effect_size_glass"


def detect_spots(
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
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (per_spot_df, per_nucleus_df). Spots are assigned to the nucleus whose mask they fall
    in; `effect_size_min` (>0) drops spots below that spot-vs-background effect size."""
    import spotmax.pipe as P

    sp = np.asarray(spacing, dtype=float)
    img = np.asarray(spot_img).astype("float32")
    radii_px = spot_radius_um / sp
    gauss = gauss_sigma_um / float(sp.min())

    pre = P.preprocess_image(img, lab=labels, gauss_sigma=gauss, use_gpu=False,
                             spots_zyx_radii_pxl=radii_px)
    segm = P.spots_semantic_segmentation(
        pre, lab=labels, spots_zyx_radii_pxl=radii_px, do_try_all_thresholds=False,
        return_only_segm=True, thresholding_method=thresholding_method, use_gpu=False)
    det = P.spot_detection(pre, spots_segmantic_segm=segm, spots_zyx_radii_pxl=radii_px,
                           lab=labels, return_df=True)
    df = det[0] if isinstance(det, tuple) else det
    if df is None or len(df) == 0:
        return _empty(labels, marker)

    zi = np.clip(df["z"].astype(int), 0, labels.shape[0] - 1)
    yi = np.clip(df["y"].astype(int), 0, labels.shape[1] - 1)
    xi = np.clip(df["x"].astype(int), 0, labels.shape[2] - 1)
    df = df.assign(nucleus_id=np.asarray(labels)[zi, yi, xi])
    df = df[df["nucleus_id"] > 0].reset_index(drop=True)

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
    except Exception:  # noqa: BLE001 - features are an enrichment; detection counts still stand
        pass

    if effect_size_min > 0 and df["effect_size"].notna().any():
        df = df[df["effect_size"].fillna(np.inf) >= effect_size_min].reset_index(drop=True)

    per_spot = pd.DataFrame({
        "spot_id": np.arange(len(df)),
        "nucleus_id": df["nucleus_id"].astype(int),
        "marker": marker,
        "z_um": df["z"].to_numpy() * sp[0],
        "y_um": df["y"].to_numpy() * sp[1],
        "x_um": df["x"].to_numpy() * sp[2],
        "effect_size": df["effect_size"].to_numpy(),
        "spot_mean_intensity": df["spot_mean_intensity"].to_numpy(),
        "detector": "spotmax",
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
        "detector": "spotmax",
    }, columns=PER_NUC_COLS)
    return per_spot, per_nuc


def _empty(labels, marker):
    nuc_ids = [int(v) for v in np.unique(labels) if v != 0]
    per_nuc = pd.DataFrame({
        "nucleus_id": nuc_ids, "marker": marker, "n_spots": 0,
        "mean_spot_intensity": np.nan, "detector": "spotmax",
    }, columns=PER_NUC_COLS)
    return pd.DataFrame(columns=PER_SPOT_COLS), per_nuc
