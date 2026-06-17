"""Read Imaris Statistics .xlsx exports (the lab's RAD-51 ground truth) and compute per-nucleus
RAD-51 counts the SAME way the lab's R analysis does — by matching Spots to Surfaces on the
coloc-channel max intensity.

Confirmed export structure (~66 sheets):
  * 'Position' sheet (skip 1 header row): Position X/Y/Z, Unit, Category (Spot|Surface|
    MeasurementPoint), 'Surpass Object' ('Spots 1' = RAD-51 foci, 'Surfaces 1' = nuclei,
    'Measurement Points 1' = the hand-drawn gonad centerline), ID.
  * 'Intensity Max Ch=N Img=1' sheets: per-object max intensity. The LAST one is the coloc
    channel the lab uses; a spot and the surface it sits in share that channel's max value, so
    grouping spots by Coloc_Ch_Max and joining to surfaces recovers spots-per-nucleus.

The lab's channel order on this dataset: 488 = SYP-2, 555 = RAD-51 ('Spots 1').
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)
_POS_REN = {"Position X": "X", "Position Y": "Y", "Position Z": "Z", "Surpass Object": "Object"}


def read_imaris_xlsx(path: str | Path):
    """Return (position_df, coloc_df, intensity_max_sheet_names). position_df has X/Y/Z, Category,
    Object, ID; coloc_df has Coloc_Ch_Max (last Intensity-Max channel), Category, Object, ID."""
    path = str(path)
    xl = pd.ExcelFile(path)
    pos = pd.read_excel(path, sheet_name="Position", skiprows=1).rename(columns=_POS_REN)
    imax_sheets = [s for s in xl.sheet_names if s.startswith("Intensity Max")]
    if not imax_sheets:
        raise ValueError(f"no 'Intensity Max' sheet in {path}")
    coloc = (pd.read_excel(path, sheet_name=imax_sheets[-1], skiprows=1)
             .rename(columns={"Intensity Max": "Coloc_Ch_Max", "Surpass Object": "Object"}))
    return pos, coloc, imax_sheets


def per_nucleus_rad51(path: str | Path, spots_object: str = "Spots 1"):
    """Per-nucleus RAD-51 counts via the lab's coloc-max matching. Returns (surfaces_df, spots_df):
    surfaces_df has one row per nucleus (Surface) with X/Y/Z, ID and an integer ``RAD51`` count."""
    pos, coloc, _ = read_imaris_xlsx(path)
    pos = pos[pos["Category"].isin(["Spot", "Surface"])].copy()
    coloc = coloc[coloc["Category"].isin(["Spot", "Surface"])][["Coloc_Ch_Max", "Category", "Object", "ID"]]
    d = pos.merge(coloc, on=["Category", "Object", "ID"], how="left")

    spots = d[(d["Category"] == "Spot") & (d["Object"] == spots_object)].copy()
    surfaces = d[d["Category"] == "Surface"].copy()

    counts = spots.groupby("Coloc_Ch_Max").size().rename("_n").reset_index()
    # Surfaces should each have a UNIQUE coloc label (Imaris fills every nucleus with a distinct random
    # value). When two surfaces collide on a value (~3% of nuclei on some gonads), the lab's coloc join
    # gives EACH the full spot count -> double-counts those spots and biases the GT up. Split the count
    # evenly across colliding surfaces so the total is preserved, and warn so it can't pass silently.
    n_share = surfaces.groupby("Coloc_Ch_Max")["ID"].transform("size").to_numpy()
    n_collide = int((n_share > 1).sum())
    if n_collide:
        log.warning("%s: %d surfaces share a coloc label (Imaris mask-value collision); "
                    "splitting RAD-51 counts evenly across them", Path(path).name, n_collide)
    surfaces = surfaces.merge(counts, on="Coloc_Ch_Max", how="left")
    surfaces["RAD51"] = (surfaces["_n"].fillna(0).to_numpy() / n_share).round().astype(int)
    surfaces = surfaces.drop(columns="_n")
    return surfaces, spots


def summarize(path: str | Path) -> dict:
    """One-line GT summary for a gonad: nuclei, total RAD-51 spots, and spots/nucleus stats."""
    surfaces, spots = per_nucleus_rad51(path)
    assigned = int(surfaces["RAD51"].sum())
    r = surfaces["RAD51"]
    return {
        "file": Path(path).name,
        "n_nuclei": int(len(surfaces)),
        "n_spots_total": int(len(spots)),
        "n_spots_assigned": assigned,
        "rad51_per_nucleus_mean": float(r.mean()) if len(r) else float("nan"),
        "rad51_per_nucleus_median": float(r.median()) if len(r) else float("nan"),
        "rad51_per_nucleus_max": int(r.max()) if len(r) else 0,
        "frac_nuclei_with_rad51": float((r > 0).mean()) if len(r) else float("nan"),
    }
