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


# ---------------------------------------------------------------------------------------------------------
# Per-object statistic sheets (protocol changelog entry 9). The full Imaris export has 66 sheets: a title row,
# a header row, then one row per object (Category Spot, Surface or MeasurementPoint, 'Surpass Object' the
# object name, ID the Imaris object ID). Read with openpyxl in read-only mode, one sheet at a time.
# ---------------------------------------------------------------------------------------------------------

_XYZ = {"Position X": "x", "Position Y": "y", "Position Z": "z"}


def _workbook(path):
    import openpyxl

    from .benchmark_split import refuse_sealed

    refuse_sealed(path)
    return openpyxl.load_workbook(str(path), read_only=True, data_only=True)


def read_sheet(path_or_wb, name: str) -> pd.DataFrame:
    """One sheet of an Imaris statistics export as a table (the second row is the header). KeyError when
    the export has no such sheet."""
    wb = path_or_wb if hasattr(path_or_wb, "sheetnames") else _workbook(path_or_wb)
    if name not in wb.sheetnames:
        raise KeyError(f"no sheet {name!r} (the export has {len(wb.sheetnames)} sheets)")
    rows = wb[name].iter_rows(values_only=True)
    next(rows, None)                                   # the title row
    header = [str(h) if h is not None else "" for h in (next(rows, None) or [])]
    data = [r for r in rows if r is not None and any(v not in (None, "") for v in r)]
    return pd.DataFrame(data, columns=header) if data else pd.DataFrame(columns=header)


def overall_totals(path_or_wb, spots_object: str = "Spots 1", surfaces_object: str = "Surfaces 1") -> dict:
    """The totals of the Overall sheet: ``n_spots`` (Total Number of Spots of ``spots_object``), ``n_surfaces``
    (Total Number of Surfaces of ``surfaces_object``) and ``n_voxels_image`` (Number of Voxels of the image);
    None where the sheet has no such row."""
    ov = read_sheet(path_or_wb, "Overall")

    def one(variable, obj=None):
        m = ov["Variable"] == variable
        if obj is not None and "Surpass Object" in ov.columns:
            m &= ov["Surpass Object"] == obj
        v = ov.loc[m, "Value"]
        return int(v.iloc[0]) if len(v) else None

    return {"n_spots": one("Total Number of Spots", spots_object),
            "n_surfaces": one("Total Number of Surfaces", surfaces_object),
            "n_voxels_image": one("Number of Voxels")}


def _object_rows(df: pd.DataFrame, category: str, obj: str) -> pd.DataFrame:
    m = df["Category"] == category
    if "Surpass Object" in df.columns:
        m &= df["Surpass Object"] == obj
    return df[m]


def _stat(wb, sheet: str, value_col: str, category: str, obj: str) -> pd.Series:
    df = _object_rows(read_sheet(wb, sheet), category, obj)
    return pd.Series(pd.to_numeric(df[value_col], errors="coerce").to_numpy(float),
                     index=df["ID"].astype(np.int64).to_numpy())


def read_export_objects(path, *, spot_channel: int, mask_channel: int, spots_object: str = "Spots 1",
                        surfaces_object: str = "Surfaces 1") -> dict:
    """The spots and surfaces of an Imaris statistics export with the statistics entry 9 (e) checks.
    ``spot_channel`` and ``mask_channel`` are 0-based DataSet channel indices (the sheets count from 1).

    ``spots``: ID, x, y, z (Imaris stage um), diameter_x, diameter_y, diameter_z, ic_spot and ic_mask (Intensity
    Center of the spot and mask channels) and imax_mask (Intensity Max of the mask channel: the lab's coloc
    value). ``surfaces``: ID, x, y, z (the Position sheet), volume and imax_mask. ``overall``: the totals of
    ``overall_totals``. ``sheets``: the sheet names."""
    wb = _workbook(path)
    try:
        pos = read_sheet(wb, "Position").rename(columns=_XYZ)
        sp = _object_rows(pos, "Spot", spots_object)
        su = _object_rows(pos, "Surface", surfaces_object)
        spots = pd.DataFrame({"ID": sp["ID"].astype(np.int64).to_numpy(),
                              **{c: pd.to_numeric(sp[c]).to_numpy(float) for c in ("x", "y", "z")}})
        surfaces = pd.DataFrame({"ID": su["ID"].astype(np.int64).to_numpy(),
                                 **{c: pd.to_numeric(su[c]).to_numpy(float) for c in ("x", "y", "z")}})
        dia = _object_rows(read_sheet(wb, "Diameter"), "Spot", spots_object).set_index("ID")
        for c, k in (("Diameter X", "diameter_x"), ("Diameter Y", "diameter_y"), ("Diameter Z", "diameter_z")):
            spots[k] = pd.to_numeric(dia[c], errors="coerce").reindex(spots["ID"]).to_numpy(float)
        sc, mc = int(spot_channel) + 1, int(mask_channel) + 1
        spots["ic_spot"] = _stat(wb, f"Intensity Center Ch={sc} Img=1", "Intensity Center", "Spot",
                                 spots_object).reindex(spots["ID"]).to_numpy()
        spots["ic_mask"] = _stat(wb, f"Intensity Center Ch={mc} Img=1", "Intensity Center", "Spot",
                                 spots_object).reindex(spots["ID"]).to_numpy()
        imax = f"Intensity Max Ch={mc} Img=1"
        spots["imax_mask"] = _stat(wb, imax, "Intensity Max", "Spot", spots_object).reindex(spots["ID"]).to_numpy()
        surfaces["imax_mask"] = _stat(wb, imax, "Intensity Max", "Surface",
                                      surfaces_object).reindex(surfaces["ID"]).to_numpy()
        surfaces["volume"] = _stat(wb, "Volume", "Volume", "Surface", surfaces_object).reindex(surfaces["ID"]).to_numpy()
        return {"spots": spots, "surfaces": surfaces,
                "overall": overall_totals(wb, spots_object, surfaces_object), "sheets": list(wb.sheetnames)}
    finally:
        wb.close()
