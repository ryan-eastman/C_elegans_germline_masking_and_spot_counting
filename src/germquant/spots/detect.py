"""RAD-51 (or other) spot quantification via SpotMAX — replaces the v1 blob_log foci detector.

Per nucleus, detect spots in the spot channel *inside the nucleus mask* using SpotMAX's validated
chain: preprocess -> semantic segmentation (where spots can be) -> peak detection -> per-spot
features incl. effect-size (spot vs background). Returns tidy per-spot and per-nucleus tables.

SpotMAX runs HEADLESS here (no GUI, no interactive prompts). Detection/threshold/effect-size
parameters live in germquant config (`spots:`) and are calibrated against Imaris counts — the GUI
is for *viewing*, not counting. Validated chain: scripts/smoke_spotmax.py.
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

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
    merge_z_columns: bool = True,
    z_merge_gap_um: float = 0.8,
    z_merge_valley_frac: float = 0.8,
    max_spot_candidates: int = 30000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (per_spot_df, per_nucleus_df). Spots are assigned to the nucleus whose mask they fall
    in; `effect_size_min` (>0) drops spots below that spot-vs-background effect size.

    `merge_z_columns` (default on) collapses the z-axis spot-splitting artifact: confocal axial PSF
    (~0.6-0.8 um) is far wider than the z spot-radius, so SpotMAX detects one focus as 2-3 stacked
    peaks at the SAME (x,y) voxel. Two stacked peaks are merged ONLY if they are within
    `z_merge_gap_um` in z AND there is no real intensity valley between them along z (the dip stays
    above `z_merge_valley_frac` * the dimmer peak) — i.e. one PSF-blurred focus, not two distinct
    foci. A genuine valley keeps them separate, so dense pachytene nuclei are not over-merged.
    Validated vs Imaris on 20251105_N2_HERM_001 (2033 -> ~Imaris 1222). effect_size_min is NOT the right
    knob for the over-count (the extras are real bright peaks of one focus, not dim noise)."""
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

    # FLOOD GUARD: a high-signal or artefact image can propose tens-to-hundreds of thousands of
    # candidate peaks, and the per-spot feature step below is O(n) and CPU-bound — on one syp-2 mutant
    # gonad it ran for 8 h without finishing. Cap to the brightest `max_spot_candidates` before
    # features so the pipeline can never wedge. A real RAD-51 count is far below this, so hitting the
    # cap means the image is flooded (artefact / bleed-through / over-bright) and its count is a floor,
    # not a measurement — surfaced via a loud warning.
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


def _cap_candidates(df: pd.DataFrame, image: np.ndarray, cap: int) -> pd.DataFrame:
    """Flood guard: if more than `cap` candidate spots, keep the brightest `cap` (by raw intensity at
    the peak voxel), in original order. Bounds the O(n) per-spot feature step so an artefact /
    over-bright image can't wedge the run for hours. A real RAD-51 count is far below `cap`, so
    triggering this means the image is flooded and its count is a floor, not a measurement."""
    if len(df) <= cap:
        return df
    z = np.clip(df["z"].astype(int), 0, image.shape[0] - 1)
    y = np.clip(df["y"].astype(int), 0, image.shape[1] - 1)
    x = np.clip(df["x"].astype(int), 0, image.shape[2] - 1)
    bright = np.asarray(image)[z, y, x]
    keep = np.sort(np.argsort(bright)[::-1][:cap])
    log.warning("FLOODED: %d spot candidates exceed cap %d — likely artefact/over-bright; keeping the "
                "brightest %d, count is a FLOOR not a measurement.", len(df), cap, cap)
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


def _empty(labels, marker):
    nuc_ids = [int(v) for v in np.unique(labels) if v != 0]
    per_nuc = pd.DataFrame({
        "nucleus_id": nuc_ids, "marker": marker, "n_spots": 0,
        "mean_spot_intensity": np.nan, "detector": "spotmax",
    }, columns=PER_NUC_COLS)
    return pd.DataFrame(columns=PER_SPOT_COLS), per_nuc
