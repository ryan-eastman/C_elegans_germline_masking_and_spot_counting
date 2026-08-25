"""Validate the PGL-1 granule + SYP<->PGL-1 coloc output against Imaris ground truth.

Two comparisons, both source-agnostic (they take tidy DataFrames, so they work whether the Imaris
data came from an .ims via `read_imaris_ims` or an .xlsx Statistics export):

  * compare_granules — match our surfaced PGL-1 granules to Imaris' PGL-1 Surfaces by nearest
    centroid, then report detection recall/precision and per-granule VOLUME agreement (Lin's CCC,
    Bland-Altman) via the shared `agreement_stats`.
  * compare_coloc_metrics — put our Manders M1/M2 + Pearson next to the numbers Imaris' Colocalization
    module reports for the same SYP x PGL-1 pair, with absolute differences. (Our Manders are
    mask-restricted and Imaris' are Costes-thresholded, so expect the same ballpark, not identity —
    this is a sanity cross-check, not a pass/fail.)
"""
from __future__ import annotations

import pandas as pd

from .compare import agreement_stats

_XYZ = ["z_um", "y_um", "x_um"]


def _match(our_df: pd.DataFrame, ref_df: pd.DataFrame, max_match_um: float):
    """Nearest-centroid match our objects -> reference objects. Returns (dist, ref_idx, matched_mask)."""
    from scipy.spatial import cKDTree

    ours = our_df[_XYZ].to_numpy(float)
    refs = ref_df[_XYZ].to_numpy(float)
    dist, idx = cKDTree(refs).query(ours)
    return dist, idx, dist <= max_match_um


def compare_granules(our_df: pd.DataFrame, ref_df: pd.DataFrame, *, max_match_um: float = 1.5) -> dict:
    """Detection + volume agreement of our PGL-1 granules vs an Imaris PGL-1 Surfaces table.

    Both frames need z_um/y_um/x_um; volume agreement additionally needs volume_um3. `max_match_um`
    is the centroid cutoff for calling two objects the same granule."""
    out = {"n_our": int(len(our_df)), "n_ref": int(len(ref_df)), "n_matched": 0,
           "recall_of_ref": float("nan"), "precision": float("nan"),
           "mean_match_dist_um": float("nan"), "volume_stats": {}}
    if our_df.empty or ref_df.empty:
        return out
    dist, idx, matched = _match(our_df, ref_df, max_match_um)
    n_match = int(matched.sum())
    out["n_matched"] = n_match
    out["recall_of_ref"] = n_match / len(ref_df)
    out["precision"] = n_match / len(our_df)
    out["mean_match_dist_um"] = float(dist[matched].mean()) if n_match else float("nan")
    if n_match and "volume_um3" in our_df.columns and "volume_um3" in ref_df.columns:
        ov = our_df["volume_um3"].to_numpy(float)[matched]
        rv = ref_df["volume_um3"].to_numpy(float)[idx[matched]]
        out["volume_stats"] = agreement_stats(ov, rv)
    return out


def compare_coloc_metrics(our_coloc_df: pd.DataFrame, imaris_metrics: dict | None, *,
                          operand: str = "shell_voxel") -> pd.DataFrame:
    """Side-by-side of our vs Imaris coloc metrics (Manders M1/M2, Pearson) with absolute diffs.

    `imaris_metrics` is a dict with any of {manders_m1, manders_m2, pearson_r} from the Imaris Coloc
    module; missing keys read NaN. `operand` picks which of our coloc rows to compare — default
    `shell_voxel` (the validated headline: SYP<->PGL-1 in the lamin-defined perinuclear shell)."""
    imaris_metrics = imaris_metrics or {}
    sub = our_coloc_df[our_coloc_df["sc_operand"] == operand] if not our_coloc_df.empty else our_coloc_df
    ours = sub.iloc[0] if len(sub) else {}
    rows = []
    for key in ("manders_m1", "manders_m2", "pearson_r"):
        o = float(ours[key]) if key in getattr(ours, "index", []) and pd.notna(ours[key]) else float("nan")
        im = float(imaris_metrics.get(key, float("nan")))
        rows.append({"metric": key, "ours": o, "imaris": im, "abs_diff": abs(o - im)})
    return pd.DataFrame(rows, columns=["metric", "ours", "imaris", "abs_diff"])
