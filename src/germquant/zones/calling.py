"""Call meiotic zones (transition zone vs pachytene) along the germline axis.

Two methods (auto-selected):
  * dapi_crescent  — when DNA/DAPI present. Per-nucleus chromatin POLARIZATION (offset of
    the DAPI intensity-weighted centroid from the geometric centroid, in nuclear radii)
    detects the crescent/clustered morphology of leptotene/zygotene. Transition-zone
    boundary = distal-most axis row with >=2 crescent nuclei (count_2plus) OR >=60%
    crescent (pct_60) — NEVER 80% (see ARCHITECTURE.md §3.4).
  * synapsis_state — no-DAPI fallback (axis vs central-element). Stub; novel contribution,
    validate against DAPI gonads before use. Wild-type only.

  *** v1 heuristic — calibrate the polarization threshold against hand-scored gonads. ***
"""
from __future__ import annotations

import numpy as np
import pandas as pd

ZONE_COLS = [
    "zone", "length_um", "length_rows", "fraction_of_germline", "n_nuclei",
    "distal_boundary_position_norm", "proximal_boundary_position_norm", "boundary_method",
]


def chromatin_polarization(
    labels: np.ndarray, dna: np.ndarray, spacing: tuple[float, float, float]
) -> dict[int, float]:
    """Per-nucleus polarization = |weighted_centroid - geom_centroid| / nuclear_radius."""
    from skimage.measure import regionprops

    sp = np.asarray(spacing, float)
    vox_vol = float(np.prod(sp))
    out: dict[int, float] = {}
    for rp in regionprops(labels, intensity_image=dna):
        geom = np.asarray(rp.centroid) * sp
        wt = np.asarray(rp.centroid_weighted) * sp
        offset = float(np.linalg.norm(wt - geom))
        radius = (3 * rp.area * vox_vol / (4 * np.pi)) ** (1 / 3)
        out[int(rp.label)] = offset / radius if radius > 0 else 0.0
    return out


def call_zones(
    nuclei: pd.DataFrame,
    *,
    dna: np.ndarray | None = None,
    labels: np.ndarray | None = None,
    spacing: tuple[float, float, float] = (1, 1, 1),
    method: str = "auto",
    crescent_boundary: str = "count_2plus",
    axis_bin_um: float = 5.0,
    polarization_threshold: float = 0.18,
    pachytene_thirds: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Return (nuclei_annotated, zones_table, flags)."""
    flags: list[str] = []
    df = nuclei.copy()
    use_dapi = (method in ("auto", "dapi_crescent")) and dna is not None and labels is not None
    if method == "synapsis_state" or (method == "auto" and not use_dapi):
        flags.append("zones:synapsis_state_no_dapi_NOVEL_validate")
        df["zone_call"] = "unknown"
        df["zone_call_method"] = "synapsis_state_stub"
        df["pachytene_subzone"] = "NA"
        df["zone_confidence"] = float("nan")
        df["chromatin_polarized"] = False
        return df, pd.DataFrame(columns=ZONE_COLS), flags

    pol = chromatin_polarization(labels, dna, spacing)
    df["chromatin_polarization"] = df["nucleus_id"].map(pol).fillna(0.0)
    df["chromatin_polarized"] = df["chromatin_polarization"] >= polarization_threshold
    df["zone_call_method"] = "dapi_crescent"

    if "axis_position_um" not in df or df["axis_position_um"].isna().all():
        flags.append("zones:no_axis_position")
        df["zone_call"] = "unknown"
        df["pachytene_subzone"] = "NA"
        df["zone_confidence"] = float("nan")
        return df, pd.DataFrame(columns=ZONE_COLS), flags

    # bin nuclei into rows along the axis
    pos = df["axis_position_um"].to_numpy()
    total_len = float(np.nanmax(pos)) if np.isfinite(pos).any() else 0.0
    n_rows = max(1, int(np.ceil(total_len / axis_bin_um)))
    row_idx = np.clip((pos / axis_bin_um).astype(int), 0, n_rows - 1)
    df["axis_row"] = row_idx

    # per-row crescent statistics
    tz_start_row, tz_end_row = _transition_zone_rows(df, n_rows, crescent_boundary)

    def row_to_zone(r: int) -> str:
        if tz_start_row is None:
            return "pachytene"
        if r < tz_start_row:
            return "premeiotic"
        if r <= tz_end_row:
            return "transition_zone"
        return "pachytene"

    df["zone_call"] = df["axis_row"].map(row_to_zone)
    df["zone_confidence"] = float(np.clip((df["chromatin_polarization"] - polarization_threshold).abs().mean() * 5, 0, 1))

    # pachytene thirds
    df["pachytene_subzone"] = "NA"
    if pachytene_thirds:
        pach = df["zone_call"] == "pachytene"
        if pach.any():
            pp = df.loc[pach, "axis_position_norm"]
            lo, hi = pp.min(), pp.max()
            if hi > lo:
                third = (pp - lo) / (hi - lo)
                df.loc[pach, "pachytene_subzone"] = np.where(
                    third < 1 / 3, "early", np.where(third < 2 / 3, "mid", "late")
                )

    zones_tbl = _zones_table(df, total_len, n_rows, crescent_boundary)
    return df, zones_tbl, flags


def _transition_zone_rows(df, n_rows, boundary):
    rows_cresc = df.groupby("axis_row")["chromatin_polarized"].agg(["sum", "mean", "count"])
    def is_tz(r):
        if r not in rows_cresc.index:
            return False
        s = rows_cresc.loc[r]
        return (s["sum"] >= 2) if boundary == "count_2plus" else (s["mean"] >= 0.60)
    flags = [is_tz(r) for r in range(n_rows)]
    if not any(flags):
        return None, None
    # The transition zone is a COMPACT DISTAL band, not every scattered crescent row. Take the
    # distal-most contiguous run of crescent rows (tolerating single-row gaps) and stop at the
    # first sustained (>1 row) gap — otherwise a lone polarized nucleus deep in pachytene would
    # stretch the TZ across most of the gonad (observed on real N2 data: TZ called at 60%).
    start = next(r for r in range(n_rows) if flags[r])
    end, gap = start, 0
    for r in range(start + 1, n_rows):
        if flags[r]:
            end, gap = r, 0
        else:
            gap += 1
            if gap > 1:
                break
    return start, end


def _zones_table(df, total_len, n_rows, boundary):
    rows = []
    for zone in ("transition_zone", "pachytene"):
        sub = df[df["zone_call"] == zone]
        if sub.empty:
            continue
        span = float(sub["axis_position_um"].max() - sub["axis_position_um"].min())
        rows.append(
            {
                "zone": zone,
                "length_um": span,
                "length_rows": int(sub["axis_row"].nunique()),
                "fraction_of_germline": (span / total_len) if total_len else float("nan"),
                "n_nuclei": int(len(sub)),
                "distal_boundary_position_norm": float(sub["axis_position_norm"].min()),
                "proximal_boundary_position_norm": float(sub["axis_position_norm"].max()),
                "boundary_method": boundary,
            }
        )
    return pd.DataFrame(rows, columns=ZONE_COLS)
