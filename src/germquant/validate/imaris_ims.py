"""Read objects directly from an Imaris .ims (HDF5) project — Spots (RAD-51 foci) and Surfaces
(nuclei) — in IMAGE-RELATIVE microns, so they can be matched against our pipeline output for the
exact same image (our centroids are voxel*spacing, origin at the image corner).

Imaris stores stage coordinates offset by the image extent minimum (DataSetInfo/Image ExtMin{0,1,2}
for X/Y/Z). Subtracting ExtMin puts objects in the same origin-at-corner frame as our pipeline.
Surface VOXEL masks live in Imaris's proprietary MegaSurfaces block tree (not extracted here); we
read per-surface centroid + volume from SurfaceModelInfo, which is enough for count/centroid checks.

The second half of this module is the Imaris TRUTH READER (``export_imaris_truth``, the CLI
``germquant imaris-truth``). It writes the truth file contract for one .ims, the files the benchmark
(package B3) consumes: ``<stem>__imaris_spots.csv``, ``<stem>__imaris_nuclei.csv``,
``<stem>__imaris_labels.tif`` and ``<stem>__imaris_params.json``. Two frames matter:

* Imaris frame: micrometres from the stage origin, voxel i of an axis covering
  ``ExtMin + i * v_imaris .. ExtMin + (i + 1) * v_imaris`` with Imaris's ROUNDED voxel size
  ``v_imaris = (ExtMax - ExtMin) / size`` (0.108 um here, the nd2 says 0.108333). An object position
  is converted to a voxel index with Imaris's own voxel size, ``floor((p - ExtMin) / v_imaris)``.
* Our frame: voxel index times the nd2 voxel size (what the pipeline tables carry). Matching raw
  micrometres between the frames drifts by almost a micrometre across a 2600-pixel field.

Nucleus masks are not read from the block tree either: the lab's export fills every surface with a
distinct value in an extra DataSet channel ("Masked 405"; the coloc-channel trick). The reader finds
that channel, matches each surface to the blob whose centroid coincides with the surface's centre of
mass, splits shared values and drops values no surface claims. docs/IMARIS_TRUTH.md has the details
and what was found in the real files.
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

from ..fsutil import long_path

log = logging.getLogger(__name__)

_SPOTS = "Scene8/Content/Points0/Spot"
_SURF = "Scene8/Content/MegaSurfaces0/SurfaceModelInfo"


def _fattr(group, key: str) -> float:
    v = group.attrs[key]
    if isinstance(v, np.ndarray) and v.dtype.kind == "S":
        return float(v.tobytes().decode("latin1").strip("\x00"))
    if isinstance(v, np.ndarray):
        return float(v.ravel()[0])
    return float(v)


def read_imaris_ims(path: str | Path):
    """Return (spots_df, surfaces_df, geom). spots_df: x_um/y_um/z_um/radius_um (image-relative).
    surfaces_df: x_um/y_um/z_um/volume_um3 (image-relative). geom: extent + shape + voxel size."""
    import h5py

    with h5py.File(str(path), "r") as f:
        info = f["DataSetInfo"]["Image"]
        ext_min = np.array([_fattr(info, f"ExtMin{i}") for i in range(3)])  # x, y, z (µm)
        ext_max = np.array([_fattr(info, f"ExtMax{i}") for i in range(3)])
        shape_xyz = np.array([int(_fattr(info, k)) for k in ("X", "Y", "Z")])
        spot = f[_SPOTS][:] if _SPOTS in f else None
        smi = f[_SURF][:] if _SURF in f else None

    spots = pd.DataFrame(columns=["x_um", "y_um", "z_um", "radius_um"])
    if spot is not None and len(spot):
        spots = pd.DataFrame({
            "x_um": spot["PositionX"] - ext_min[0],
            "y_um": spot["PositionY"] - ext_min[1],
            "z_um": spot["PositionZ"] - ext_min[2],
            "radius_um": spot["Radius"],
        })
    surfaces = pd.DataFrame(columns=["x_um", "y_um", "z_um", "volume_um3"])
    if smi is not None and len(smi):
        surfaces = pd.DataFrame({
            "x_um": smi["CenterOfMassX"] - ext_min[0],
            "y_um": smi["CenterOfMassY"] - ext_min[1],
            "z_um": smi["CenterOfMassZ"] - ext_min[2],
            "volume_um3": smi["Volume"],
        })
    voxel_xyz = (ext_max - ext_min) / np.maximum(shape_xyz, 1)
    geom = {"ext_min_xyz": ext_min, "ext_max_xyz": ext_max,
            "shape_xyz": shape_xyz, "voxel_xyz_um": voxel_xyz}
    return spots, surfaces, geom


def spots_per_surface(spots: pd.DataFrame, surfaces: pd.DataFrame, max_dist_um: float = 3.0,
                      use_surface_radius: bool = True) -> np.ndarray:
    """Assign each spot to the nearest surface centroid and count per surface. By default the cutoff
    is the surface's own sphere-equivalent radius from `volume_um3` (capped at max_dist_um) rather
    than a flat max_dist_um — a flat 3 um exceeds the ~2 um nucleus radius and over-counts. This is
    still an approximation (centroid-based, no voxel masks) — the lab's coloc method (imaris_xlsx) is
    canonical for per-nucleus counts."""
    if surfaces.empty:
        return np.array([], dtype=int)
    sc = surfaces[["x_um", "y_um", "z_um"]].to_numpy()
    if use_surface_radius and "volume_um3" in surfaces.columns:
        cutoff = np.minimum((3 * surfaces["volume_um3"].to_numpy() / (4 * np.pi)) ** (1 / 3), max_dist_um)
    else:
        cutoff = np.full(len(surfaces), max_dist_um)
    counts = np.zeros(len(surfaces), dtype=int)
    if spots.empty:
        return counts
    sp = spots[["x_um", "y_um", "z_um"]].to_numpy()
    for p in sp:
        d = np.linalg.norm(sc - p, axis=1)
        j = int(d.argmin())
        if d[j] <= cutoff[j]:
            counts[j] += 1
    return counts


# ======================================================================================================
# Imaris truth reader: the truth file contract (package B1 writes it, package B3 reads it)
# ======================================================================================================

TRUTH_SUFFIXES = {"spots": "__imaris_spots.csv", "nuclei": "__imaris_nuclei.csv",
                  "labels": "__imaris_labels.tif", "params": "__imaris_params.json"}
#: the contract columns, in order; the extra columns follow them (a consumer selects by name)
SPOT_COLUMNS = ["spot_id", "x_um", "y_um", "z_um", "xi", "yi", "zi", "radius_um", "quality", "spot_class",
                "imaris_x_um", "imaris_y_um", "imaris_z_um"]
SPOT_EXTRA_COLUMNS = ["imaris_mask_label", "spots_object"]
NUCLEI_COLUMNS = ["surface_id", "label", "x_um", "y_um", "z_um", "xi", "yi", "zi", "volume_um3",
                  "n_voxels_mask", "imaris_x_um", "imaris_y_um", "imaris_z_um"]
NUCLEI_EXTRA_COLUMNS = ["mask_centroid_dist_um", "surface_class", "surfaces_object"]

_CLASS_SEP = ";"          # an object carrying several labels: "Class A;Class B"
_TOL_UM = 1.0             # surface centre of mass to mask-blob centroid (real files: at most 0.04 um)
_PROBE_HALF = 3           # pixel probe: a (2*3+1)^2 patch of the centroid plane around each surface centroid
_STATS_CHUNK = 4_000_000  # StatisticsValue rows per read (32 bytes each; 772k-spot exports hold 22M rows)
_RELABEL_PLANES = 2       # z planes relabelled per step (fancy indexing makes an 8-byte index copy)
_BOX_MARGIN = 2           # voxels added around Imaris's surface bounding box when looking for its mask
_BOX_COLS = ["_bmin_x", "_bmin_y", "_bmin_z", "_bmax_x", "_bmax_y", "_bmax_z"]   # hidden, Imaris um


def truth_paths(out_dir: str | Path, stem: str) -> dict[str, Path]:
    """The four truth files of an .ims named ``<stem>.ims`` inside ``out_dir``."""
    out_dir = Path(out_dir)
    return {k: out_dir / f"{stem}{suffix}" for k, suffix in TRUTH_SUFFIXES.items()}


# ---- low-level HDF5 helpers --------------------------------------------------------------------------

def _decode(v) -> str:
    """Imaris stores text attributes as arrays of single bytes and text datasets as fixed-length byte
    strings, both NUL padded; return plain text either way."""
    if isinstance(v, np.ndarray):
        if v.dtype.kind == "S":
            return v.tobytes().decode("latin1", "replace").rstrip("\x00")
        if v.dtype.kind == "U":
            return "".join(v.ravel().tolist()).rstrip("\x00")
        v = v.ravel()[0] if v.size else ""
    if isinstance(v, (bytes, np.bytes_)):
        return bytes(v).decode("latin1", "replace").rstrip("\x00")
    return str(v).rstrip("\x00")


def _attr_text(group, key: str, default: str = "") -> str:
    return _decode(group.attrs[key]) if key in group.attrs else default


def _content_group(f):
    """Scene8/Content (Imaris 8+) or Scene/Content; None when the file holds no scene."""
    for scene in ("Scene8", "Scene"):
        if scene in f and "Content" in f[scene]:
            return f[scene]["Content"]
    return None


def _object_groups(content, prefix: str) -> list[tuple[str, object]]:
    """The Points*/MegaSurfaces* groups of a scene in index order, as (name, group)."""
    if content is None:
        return []
    names = [k for k in content.keys() if k.startswith(prefix) and k[len(prefix):].isdigit()]
    return [(k, content[k]) for k in sorted(names, key=lambda k: int(k[len(prefix):]))]


def read_geometry(f) -> dict:
    """Image geometry of an open .ims: extent (um), size, Imaris voxel size, channel names and the
    padded on-disk shape of the DataSet arrays (Imaris pads every axis to a chunk multiple; the image
    occupies ``[:Z, :Y, :X]``)."""
    info = f["DataSetInfo"]["Image"]
    ext_min = np.array([_fattr(info, f"ExtMin{i}") for i in range(3)], float)
    ext_max = np.array([_fattr(info, f"ExtMax{i}") for i in range(3)], float)
    shape_xyz = np.array([int(_fattr(info, k)) for k in ("X", "Y", "Z")], int)
    voxel = (ext_max - ext_min) / np.maximum(shape_xyz, 1)
    tp = f["DataSet"]["ResolutionLevel 0"]["TimePoint 0"]
    chans = sorted((k for k in tp.keys() if k.startswith("Channel ")), key=lambda k: int(k.split()[-1]))
    names = [(_attr_text(f["DataSetInfo"][k], "Name", k) if k in f["DataSetInfo"] else k) for k in chans]
    first = tp[chans[0]]["Data"] if chans else None
    version = _attr_text(f["DataSetInfo"]["Imaris"], "Version") if "Imaris" in f["DataSetInfo"] else ""
    return {"ext_min_xyz": ext_min, "ext_max_xyz": ext_max, "shape_xyz": shape_xyz, "voxel_xyz_um": voxel,
            "channel_names": names, "n_channels": len(chans),
            "data_shape_zyx": tuple(int(s) for s in first.shape) if first is not None else None,
            "data_dtype": str(first.dtype) if first is not None else None,
            "imaris_version": version, "image_name": _attr_text(info, "Name"),
            "source_file": _attr_text(info, "Filename")}


def _channel_dataset(f, channel: int):
    return f["DataSet"]["ResolutionLevel 0"]["TimePoint 0"][f"Channel {channel}"]["Data"]


def _nullable_int(values) -> pd.arrays.IntegerArray:
    """Floats with NaN to a pandas nullable integer column (written as empty cells, not 'nan')."""
    a = np.asarray(values, float)
    ok = np.isfinite(a)
    out = pd.array(np.where(ok, np.round(a), 0).astype(np.int64), dtype="Int64")
    out[~ok] = pd.NA
    return out


# ---- frame conversion --------------------------------------------------------------------------------

def imaris_to_voxel_index(pos_xyz_um, geom: dict) -> np.ndarray:
    """Imaris positions (N, 3) in stage micrometres to integer voxel indices (N, 3) with Imaris's own
    voxel size: ``floor((p - ExtMin) / v_imaris)``. Not clipped: a spot just outside the image keeps
    its -1 or size index so a consumer can see it."""
    rel = np.asarray(pos_xyz_um, dtype=float).reshape(-1, 3) - geom["ext_min_xyz"]
    return np.floor(rel / geom["voxel_xyz_um"]).astype(np.int64)


def voxel_index_to_our_um(idx_xyz, our_voxel_xyz_um) -> np.ndarray:
    """Voxel indices to our frame: index times the nd2 voxel size (the pipeline convention)."""
    return np.asarray(idx_xyz, dtype=float).reshape(-1, 3) * np.asarray(our_voxel_xyz_um, dtype=float)


def resolve_our_voxel(ims: Path, nd2: str | Path | None, nd2_voxel_xyz_um, geom: dict) -> tuple[np.ndarray, str]:
    """The voxel size of OUR frame, (dx, dy, dz) um, and where it came from: an explicit tuple, the
    header of the given nd2, a sibling ``<stem>.nd2`` next to the .ims, else Imaris's voxel size (the
    frames then coincide, which the params JSON records as a warning)."""
    if nd2_voxel_xyz_um is not None:
        return np.asarray(nd2_voxel_xyz_um, dtype=float).reshape(3), "given"
    if nd2 is None:
        sib = Path(ims).with_suffix(".nd2")
        if sib.exists():
            nd2 = sib
    if nd2 is not None:
        from .. import io as _io

        dz, dy, dx = _io.read_image_metadata(str(nd2))["spacing"]
        return np.array([dx, dy, dz], float), f"nd2:{nd2}"
    return np.array(geom["voxel_xyz_um"], float), "imaris_fallback"


# ---- creation parameters, ROI, classes, statistics ---------------------------------------------------

_ROI_RE = re.compile(r"<bpRegionOfInterest\s+([^>]*?)/?>")
_XML_ATTR_RE = re.compile(r'(\w+)="([^"]*)"')


def read_creation_parameters(group) -> str:
    """The CreationParameters XML of a Points/Surfaces object. The dataset holds the full text; the
    group attribute of the same name is a shorter copy, used only when the dataset is missing."""
    if "CreationParameters" in group:
        v = group["CreationParameters"][()]
        if isinstance(v, np.ndarray):
            v = v.ravel()[0] if v.size else b""
        return _decode(v)
    return _attr_text(group, "CreationParameters")


def parse_roi_boxes(text: str) -> list[dict]:
    """Region-of-interest boxes recorded in a CreationParameters text, in Imaris voxel indices
    (inclusive bounds, x y z order), as Imaris wrote them."""
    boxes = []
    for m in _ROI_RE.finditer(text or ""):
        attrs = dict(_XML_ATTR_RE.findall(m.group(1)))

        def ints(key):
            return [int(float(s)) for s in attrs[key].split()] if key in attrs else None

        boxes.append({"name": attrs.get("mName"), "min_xyz_voxel": ints("mMinXYZ"), "max_xyz_voxel": ints("mMaxXYZ"),
                      "min_t": ints("mMinT"), "max_t": ints("mMaxT"),
                      "active": attrs.get("mActive", "").lower() == "true"})
    return boxes


def xml_flag(text: str, name: str):
    """True/False for an ``mName="true"`` style attribute of the creation parameters, None when absent."""
    m = re.search(rf'\b{re.escape(name)}="(\w+)"', text or "")
    return None if m is None else m.group(1).lower() == "true"


def read_object_classes(group) -> tuple[dict[int, list[str]], dict[str, list[str]]]:
    """Per-object label names (Imaris "Class" labels) of a Points/Surfaces object.

    Layout: LabelValues lists the label names, LabelGroupNames the groups they belong to (cumulative end
    index), LabelSets groups the objects by their label combination: row k applies label ids
    ``LabelSetLabelIDs[prev_end:EndLabelIDs]`` to objects ``LabelSetObjectIDs[prev_end:EndObjectIDs]``.
    Returns (object id -> label names, group name -> label names). Checked against Imaris's own
    "Number of labelled Spots per Class" statistic on the real files."""
    need = ("LabelValues", "LabelSets", "LabelSetLabelIDs", "LabelSetObjectIDs")
    if any(k not in group for k in need):
        return {}, {}
    values = [_decode(r["LabelValue"]) for r in group["LabelValues"][:]]
    groups: dict[str, list[str]] = {}
    if "LabelGroupNames" in group:
        start = 0
        for r in group["LabelGroupNames"][:]:
            end = int(r["EndLabelValue"])
            groups[_decode(r["LabelGroupName"])] = values[start:end]
            start = end
    label_ids = group["LabelSetLabelIDs"][:]["IDLabel"]
    object_ids = group["LabelSetObjectIDs"][:]["IDObject"]
    out: dict[int, list[str]] = {}
    l0 = o0 = 0
    for row in group["LabelSets"][:]:
        l1, o1 = int(row["EndLabelIDs"]), int(row["EndObjectIDs"])
        names = [values[int(i)] for i in label_ids[l0:l1] if 0 <= int(i) < len(values)]
        for oid in object_ids[o0:o1]:
            out.setdefault(int(oid), []).extend(names)
        l0, o0 = l1, o1
    return out, groups


def _statistic_types(group) -> list[dict]:
    """StatisticsType rows with their factors resolved (e.g. {"Channel": "4", "Image": "Image 1"})
    and their category name ("Spot"/"Surface" per object, "Overall" aggregate)."""
    if "StatisticsType" not in group:
        return []
    factors: dict[int, dict[str, str]] = {}
    if "Factor" in group:
        for r in group["Factor"][:]:
            factors.setdefault(int(r["ID_List"]), {})[_decode(r["Name"])] = _decode(r["Level"])
    cats = {int(r["ID"]): _decode(r["CategoryName"]) for r in group["Category"][:]} if "Category" in group else {}
    out = []
    for r in group["StatisticsType"][:]:
        out.append({"id": int(r["ID"]), "name": _decode(r["Name"]), "unit": _decode(r["Unit"]),
                    "category": cats.get(int(r["ID_Category"]), str(int(r["ID_Category"]))),
                    "factors": factors.get(int(r["ID_FactorList"]), {})})
    return out


def statistic_names(group) -> list[str]:
    return sorted({t["name"] for t in _statistic_types(group) if t["category"] != "Overall"})


def read_object_statistics(group, wanted: dict[str, tuple[str, int | None]]) -> dict[str, pd.Series]:
    """Per-object statistics as Series indexed by object id. ``wanted`` maps an output name to
    (statistic name, 0-based channel index or None), e.g. {"quality": ("Quality", None),
    "mask": ("Intensity Max", 3)}. StatisticsValue is scanned in chunks (the file keeps one row per
    object and statistic, 22M rows for a 772k-spot export). Names absent from the file are left out."""
    sel: dict[int, str] = {}
    types = _statistic_types(group)
    for out, (name, ch) in wanted.items():
        for t in types:
            if t["name"] != name or t["category"] == "Overall":
                continue
            if ch is not None and t["factors"].get("Channel") != str(int(ch) + 1):
                continue
            sel[t["id"]] = out
            break
    if not sel or "StatisticsValue" not in group:
        return {}
    sv = group["StatisticsValue"]
    n = sv.shape[0]
    ids = np.array(sorted(sel), dtype=np.int64)
    parts: dict[str, tuple[list, list]] = {out: ([], []) for out in sel.values()}
    for a in range(0, n, _STATS_CHUNK):
        block = sv[a:min(a + _STATS_CHUNK, n)]
        keep = np.isin(block["ID_StatisticsType"], ids) & (block["ID_Object"] >= 0)
        block = block[keep]
        for tid, out in sel.items():
            m = block["ID_StatisticsType"] == tid
            parts[out][0].append(block["ID_Object"][m].astype(np.int64))
            parts[out][1].append(block["Value"][m].astype(float))
    series = {}
    for out, (objs, vals) in parts.items():
        o = np.concatenate(objs) if objs else np.array([], np.int64)
        v = np.concatenate(vals) if vals else np.array([], float)
        s = pd.Series(v, index=o)
        series[out] = s[~s.index.duplicated(keep="first")]
    return series


def _class_strings(ids, classes: dict[int, list[str]]) -> list[str]:
    return [_CLASS_SEP.join(classes.get(int(i), [])) for i in ids]


def _object_report(gname: str, g, kind: str, n: int, ids, classes, groups) -> dict:
    text = read_creation_parameters(g)
    counts = pd.Series(_class_strings(ids, classes)).value_counts().to_dict() if n else {}
    return {"group": gname, "name": _attr_text(g, "Name", gname), "kind": kind, "n": int(n),
            "creation_parameters": text, "roi_boxes": parse_roi_boxes(text),
            "region_of_interest_enabled": xml_flag(text, "mEnableRegionOfInterest"),
            "region_growing": xml_flag(text, "mEnableRegionGrowing"),
            "label_groups": groups, "class_counts": {str(k): int(v) for k, v in counts.items()},
            "statistics": statistic_names(g)}


# ---- objects -----------------------------------------------------------------------------------------

def read_surfaces_table(f, geom: dict, our_voxel_xyz_um) -> tuple[pd.DataFrame, list[dict]]:
    """One row per surface of every MegaSurfaces object (centre of mass and volume from
    SurfaceModelInfo, id from SurfaceModel), in both frames, plus one report per object (name,
    creation parameters, ROI boxes, classes). The mask columns are filled in by the export."""
    frames, reports = [], []
    for gname, g in _object_groups(_content_group(f), "MegaSurfaces"):
        if "SurfaceModelInfo" not in g:
            continue
        smi = g["SurfaceModelInfo"][:]
        n = len(smi)
        ids = (g["SurfaceModel"][:]["ID"].astype(np.int64) if "SurfaceModel" in g and len(g["SurfaceModel"]) == n
               else np.arange(n, dtype=np.int64))
        classes, groups = read_object_classes(g)
        reports.append(_object_report(gname, g, "surfaces", n, ids, classes, groups))
        if n == 0:
            continue
        com = np.c_[smi["CenterOfMassX"], smi["CenterOfMassY"], smi["CenterOfMassZ"]].astype(float)
        idx = imaris_to_voxel_index(com, geom)
        ours = voxel_index_to_our_um(idx, our_voxel_xyz_um)
        # Imaris's own bounding box of each surface (stage um), the search region for its mask blob
        sm = g["SurfaceModel"][:] if "SurfaceModel" in g and len(g["SurfaceModel"]) == n else None
        box = {}
        for key, field in zip(_BOX_COLS, ("ExtendMinX", "ExtendMinY", "ExtendMinZ", "ExtendMaxX", "ExtendMaxY", "ExtendMaxZ")):
            box[key] = sm[field].astype(float) if sm is not None and field in sm.dtype.names else np.full(n, np.nan)
        frames.append(pd.DataFrame({
            "surface_id": ids, "label": np.zeros(n, np.int64),
            "x_um": ours[:, 0], "y_um": ours[:, 1], "z_um": ours[:, 2],
            "xi": idx[:, 0], "yi": idx[:, 1], "zi": idx[:, 2],
            "volume_um3": smi["Volume"].astype(float), "n_voxels_mask": np.full(n, np.nan),
            "imaris_x_um": com[:, 0], "imaris_y_um": com[:, 1], "imaris_z_um": com[:, 2],
            "mask_centroid_dist_um": np.full(n, np.nan), "surface_class": _class_strings(ids, classes),
            "surfaces_object": reports[-1]["name"], "_group": gname, **box}))
    cols = NUCLEI_COLUMNS + NUCLEI_EXTRA_COLUMNS + ["_group"] + _BOX_COLS
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)
    return df[cols].copy(), reports


def read_spots_table(f, geom: dict, our_voxel_xyz_um, mask_channel: int | None = None) -> tuple[pd.DataFrame, list[dict]]:
    """One row per spot of every Points object in both frames, with the Quality statistic (NaN when the
    export did not store it) and the Class labels, plus Imaris's own per-spot maximum of the mask channel
    (``imaris_mask_label``, the lab's spot-to-nucleus assignment) when ``mask_channel`` is known."""
    frames, reports = [], []
    for gname, g in _object_groups(_content_group(f), "Points"):
        if "Spot" not in g:
            continue
        sp = g["Spot"][:]
        n = len(sp)
        ids = sp["ID"].astype(np.int64)
        classes, groups = read_object_classes(g)
        wanted = {"quality": ("Quality", None)}
        if mask_channel is not None:
            wanted["imaris_mask_label"] = ("Intensity Max", int(mask_channel))
        stats = read_object_statistics(g, wanted)
        rep = _object_report(gname, g, "spots", n, ids, classes, groups)
        rep["quality_present"] = "quality" in stats
        rep["mask_label_statistic_present"] = "imaris_mask_label" in stats
        if n == 0:
            rep["n_spots_outside_image"] = 0
            reports.append(rep)
            continue
        pos = np.c_[sp["PositionX"], sp["PositionY"], sp["PositionZ"]].astype(float)
        idx = imaris_to_voxel_index(pos, geom)
        ours = voxel_index_to_our_um(idx, our_voxel_xyz_um)
        quality = stats["quality"].reindex(ids).to_numpy() if "quality" in stats else np.full(n, np.nan)
        mask_lab = (stats["imaris_mask_label"].reindex(ids).to_numpy() if "imaris_mask_label" in stats
                    else np.full(n, np.nan))
        inside = ((idx >= 0) & (idx < geom["shape_xyz"])).all(axis=1)
        rep["n_spots_outside_image"] = int((~inside).sum())
        reports.append(rep)
        frames.append(pd.DataFrame({
            "spot_id": ids, "x_um": ours[:, 0], "y_um": ours[:, 1], "z_um": ours[:, 2],
            "xi": idx[:, 0], "yi": idx[:, 1], "zi": idx[:, 2],
            "radius_um": sp["Radius"].astype(float), "quality": quality,
            "spot_class": _class_strings(ids, classes),
            "imaris_x_um": pos[:, 0], "imaris_y_um": pos[:, 1], "imaris_z_um": pos[:, 2],
            "imaris_mask_label": _nullable_int(mask_lab), "spots_object": rep["name"]}))
    cols = SPOT_COLUMNS + SPOT_EXTRA_COLUMNS
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)
    return df[cols].copy(), reports


# ---- mask channel ------------------------------------------------------------------------------------

def _mask_hint(name: str) -> bool:
    return "mask" in (name or "").lower()


def _pick_channel(cands: dict, report: dict, n_total: int, what: str):
    accepted = [c for c, r in cands.items() if r["fraction_constant"] >= 0.5]
    if not accepted:
        report["result"] = f"no channel is constant {what} at least half of the surfaces"
        return None
    best = max(accepted, key=lambda c: (cands[c]["fraction_constant"], cands[c]["uniqueness"],
                                        _mask_hint(cands[c]["name"])))
    r = cands[best]
    report["result"] = (f"channel {best} ('{r['name']}'): one nonzero constant {what} "
                        f"{r['n_constant']}/{n_total} surfaces, {r['n_distinct']} distinct values")
    report["values_unique"] = r["uniqueness"] >= 0.5
    return best


def surface_mask_statistic(f, geom: dict, channel: int, surfaces: pd.DataFrame):
    """Imaris's own per-surface "Intensity Max" in ``channel``, aligned to ``surfaces`` (NaN where
    missing), or None when the export holds no such statistic. For the mask channel this is the
    surface's fill value, the lab's coloc trick; it serves as a cross-check of the recovered labels and
    as the label source when the masks are not read. (Intensity Min cannot serve: Imaris counts
    boundary voxels outside the fill, so Min is 0 for every surface on the real exports.)"""
    out = np.full(len(surfaces), np.nan)
    found = False
    for gname, g in _object_groups(_content_group(f), "MegaSurfaces"):
        pos = np.flatnonzero((surfaces["_group"] == gname).to_numpy())
        if len(pos) == 0:
            continue
        st = read_object_statistics(g, {"max": ("Intensity Max", int(channel))})
        if "max" not in st:
            continue
        found = True
        out[pos] = st["max"].reindex(surfaces["surface_id"].to_numpy()[pos]).to_numpy()
    return out if found else None


def identify_mask_channel_by_probe(f, geom: dict, surfaces: pd.DataFrame) -> tuple[int | None, dict]:
    """Find the DataSet channel that holds the surface masks: read the plane of each surface centroid and
    look at the patch around the centroid voxel. A mask channel shows exactly one distinct nonzero value
    in the patch of most surfaces (a cavity at the centroid still counts), and those values differ
    between surfaces; a fluorescence channel varies voxel to voxel. Cheap: one plane per distinct
    centroid z per channel, nothing else is read."""
    X, Y, Z = (int(v) for v in geom["shape_xyz"])
    idx = surfaces[["xi", "yi", "zi"]].to_numpy()
    idx = idx[((idx >= 0) & (idx < geom["shape_xyz"])).all(axis=1)]
    report = {"method": "centroid_pixel_probe", "candidates": {}}
    if len(idx) == 0:
        report["result"] = "no surface centroid inside the image"
        return None, report
    h = _PROBE_HALF
    for c in range(geom["n_channels"]):
        ds = _channel_dataset(f, c)
        hits = []
        zeros = total = 0
        for z in np.unique(idx[:, 2]):
            plane = ds[int(z), :Y, :X]
            zeros += int((plane == 0).sum())
            total += plane.size
            for x, y, _ in idx[idx[:, 2] == z]:
                patch = plane[max(y - h, 0):y + h + 1, max(x - h, 0):x + h + 1]
                vals = np.unique(patch[patch > 0])
                hits.append(int(vals[0]) if len(vals) == 1 else 0)
        hits = np.array(hits)
        n_hit = int((hits > 0).sum())
        n_unique = int(len(np.unique(hits[hits > 0])))
        report["candidates"][c] = {
            "name": geom["channel_names"][c], "n_probed": int(len(hits)), "n_constant": n_hit,
            "fraction_constant": float(n_hit / len(hits)), "n_distinct": n_unique,
            "uniqueness": float(n_unique / n_hit) if n_hit else 0.0,
            "zero_fraction_of_probed_planes": float(zeros / max(total, 1))}
    best = _pick_channel(report["candidates"], report, len(idx), "in the centroid patch of")
    return best, report


def read_channel_as_labels(f, geom: dict, channel: int) -> np.ndarray:
    """The image part of a DataSet channel as an int32 ZYX array, read in z-chunks (never the whole
    padded array at once, never as float)."""
    X, Y, Z = (int(v) for v in geom["shape_xyz"])
    ds = _channel_dataset(f, channel)
    labels = np.zeros((Z, Y, X), np.int32)
    step = int(ds.chunks[0]) if ds.chunks else 8
    for z0 in range(0, Z, step):
        z1 = min(z0 + step, Z)
        labels[z0:z1] = ds[z0:z1, :Y, :X]
    return labels


def _label_counts(labels: np.ndarray, minlength: int = 1) -> np.ndarray:
    """Voxels per label value, plane by plane (bincount on the whole volume would first copy it as
    8-byte integers)."""
    cnt = np.zeros(max(int(labels.max()) + 1, minlength), np.int64)
    for z in range(labels.shape[0]):
        c = np.bincount(labels[z].ravel(), minlength=len(cnt))
        cnt[:len(c)] += c[:len(cnt)]
    return cnt


def _voxel_centres_um(zyx: np.ndarray, geom: dict) -> np.ndarray:
    xyz = np.asarray(zyx, float)[:, ::-1]
    return geom["ext_min_xyz"] + (xyz + 0.5) * geom["voxel_xyz_um"]


def _surface_search_boxes(surfaces: pd.DataFrame, geom: dict) -> np.ndarray:
    """Per surface, the voxel box (z0, z1, y0, y1, x0, x1, inclusive, clipped) in which its mask blob is
    looked for: Imaris's own bounding box of the surface widened by ``_BOX_MARGIN`` voxels, or, when the
    export has no bounding boxes, a cube of 1.5 sphere-equivalent radii around the centre of mass."""
    n = len(surfaces)
    com = surfaces[["imaris_x_um", "imaris_y_um", "imaris_z_um"]].to_numpy(float)
    bmin = surfaces[_BOX_COLS[:3]].to_numpy(float, copy=True)
    bmax = surfaces[_BOX_COLS[3:]].to_numpy(float, copy=True)
    r = 1.5 * np.cbrt(3.0 * np.clip(surfaces["volume_um3"].to_numpy(float), 0, None) / (4.0 * np.pi))
    missing = ~(np.isfinite(bmin).all(axis=1) & np.isfinite(bmax).all(axis=1))
    bmin[missing] = com[missing] - r[missing, None]
    bmax[missing] = com[missing] + r[missing, None]
    vox, ext = geom["voxel_xyz_um"], geom["ext_min_xyz"]
    lo = np.floor((bmin - ext) / vox).astype(np.int64) - _BOX_MARGIN
    hi = np.ceil((bmax - ext) / vox).astype(np.int64) - 1 + _BOX_MARGIN
    size = np.asarray(geom["shape_xyz"], int)
    lo = np.clip(lo, 0, size - 1)
    hi = np.clip(hi, 0, size - 1)
    boxes = np.zeros((n, 6), np.int64)
    boxes[:, 0], boxes[:, 1] = lo[:, 2], hi[:, 2]
    boxes[:, 2], boxes[:, 3] = lo[:, 1], hi[:, 1]
    boxes[:, 4], boxes[:, 5] = lo[:, 0], hi[:, 0]
    return boxes


def _match_surface(labels: np.ndarray, box, com, geom: dict, tol_um: float) -> tuple[int, float]:
    """The mask value inside a surface's search box whose centroid (of the voxels inside the box) lies
    nearest to the surface's centre of mass, within ``tol_um``; (0, nan) when none does. A neighbour's
    blob reaching into the box has its centroid near the box edge and is not chosen."""
    z0, z1, y0, y1, x0, x1 = (int(v) for v in box)
    sub = labels[z0:z1 + 1, y0:y1 + 1, x0:x1 + 1]
    best, best_d = 0, np.nan
    for v in np.unique(sub[sub > 0]):
        vox = np.argwhere(sub == v) + np.array([z0, y0, x0])
        d = float(np.linalg.norm(_voxel_centres_um(vox, geom).mean(axis=0) - com))
        if d <= tol_um and (best == 0 or d < best_d):
            best, best_d = int(v), d
    return best, best_d


def _split_shared_blob(labels, v, claimants, com, geom, tol_um, next_label, assigned, dist, boxes):
    """One mask value claimed by several surfaces: connected components first (one value on separate
    nuclei), then nearest-centroid assignment inside a component holding several centroids. The first
    claimant keeps the value, the others get fresh labels. Returns (next_label, n_components, n_voronoi)."""
    from scipy import ndimage

    sl = boxes[v - 1]
    sub = labels[sl]
    comp, n = ndimage.label(sub == v)
    offs = np.array([s.start for s in sl], float)
    ccent = np.array([_voxel_centres_um(np.argwhere(comp == k) + offs, geom).mean(axis=0)
                      for k in range(1, n + 1)]).reshape(-1, 3)
    d2 = np.linalg.norm(com[claimants][:, None, :] - ccent[None, :, :], axis=2)
    comp_of, d_of = d2.argmin(axis=1) + 1, d2.min(axis=1)
    per_comp: dict[int, list[int]] = {}
    for pos, si in enumerate(claimants):
        if d_of[pos] <= tol_um:
            per_comp.setdefault(int(comp_of[pos]), []).append(si)
        else:
            assigned[si], dist[si] = 0, np.nan
    n_voronoi = 0
    keeper = claimants[0]

    def fresh(si):
        nonlocal next_label
        if si == keeper:
            return v
        lab = next_label
        next_label += 1
        return lab

    for k in range(1, n + 1):
        members = per_comp.get(k, [])
        where = comp == k
        if not members:
            sub[where] = 0                              # a component no surface claims
            continue
        if len(members) == 1:
            si = members[0]
            lab = fresh(si)
            sub[where] = lab
            assigned[si], dist[si] = lab, d_of[claimants.index(si)]
            continue
        n_voronoi += 1
        vox = np.argwhere(where)
        near = np.linalg.norm(_voxel_centres_um(vox + offs, geom)[:, None, :] - com[members][None, :, :],
                              axis=2).argmin(axis=1)
        for m_i, si in enumerate(members):
            lab = fresh(si)
            sel = vox[near == m_i]
            sub[sel[:, 0], sel[:, 1], sel[:, 2]] = lab
            assigned[si], dist[si] = lab, d_of[claimants.index(si)]
    return next_label, int(n > 1), n_voronoi


def recover_surface_masks(f, geom: dict, channel: int, surfaces: pd.DataFrame, *, values_unique: bool = True,
                          tol_um: float = _TOL_UM) -> tuple[np.ndarray, pd.DataFrame, dict]:
    """Turn the mask channel into a label image whose labels are the traced surfaces.

    1. Read the channel as int32 labels. When its values are not unique per surface (a binary mask), the
       connected components of the nonzero channel become the blobs instead.
    2. Match every surface, inside Imaris's own bounding box of that surface, to the mask value whose
       voxels there have their centroid nearest to the surface's centre of mass (within ``tol_um``). On
       the real exports the distance is at most 0.04 um and the next blob is more than 2 um away.
    3. A value claimed by several surfaces is split: first into its connected components (one value on
       two separate nuclei), then, for a component holding several centroids, by nearest surface
       centroid. The extra surfaces get new labels above the largest value.
    4. Values no surface claims (nuclei curated away after masking, small fragments) are set to 0.
    Returns (labels ZYX int32, per-surface frame with label / n_voxels_mask / mask_centroid_dist_um, report).
    """
    from scipy import ndimage

    labels = read_channel_as_labels(f, geom, channel)
    strategy = "value_per_surface"
    if not values_unique:
        comp, _n = ndimage.label(labels > 0)
        labels = comp.astype(np.int32, copy=False)
        del comp
        strategy = "connected_components_of_nonzero"
    cnt = _label_counts(labels)
    values = np.nonzero(cnt)[0]
    values = values[values > 0]
    n_surf = len(surfaces)
    com = surfaces[["imaris_x_um", "imaris_y_um", "imaris_z_um"]].to_numpy(float)
    assigned = np.zeros(n_surf, np.int64)
    dist = np.full(n_surf, np.nan)
    for i, box in enumerate(_surface_search_boxes(surfaces, geom)):
        assigned[i], dist[i] = _match_surface(labels, box, com[i], geom, tol_um)
    next_label = int(labels.max()) + 1
    claims: dict[int, list[int]] = {}
    for i, v in enumerate(assigned):
        if v:
            claims.setdefault(int(v), []).append(i)
    shared = {v: idx for v, idx in claims.items() if len(idx) > 1}
    n_split_components = n_split_voronoi = 0
    if shared:
        boxes = ndimage.find_objects(labels, max_label=max(shared))
        for v, idx in shared.items():
            next_label, n_c, n_v = _split_shared_blob(labels, v, idx, com, geom, tol_um, next_label,
                                                      assigned, dist, boxes)
            n_split_components += n_c
            n_split_voronoi += n_v
    claimed_values = {int(v) for v in assigned if v}
    orphans = [int(v) for v in values if int(v) not in claimed_values]
    n_orphan_voxels = int(cnt[orphans].sum()) if orphans else 0
    if orphans:
        lut = np.arange(max(next_label, int(labels.max()) + 1), dtype=np.int32)
        lut[orphans] = 0
        for z0 in range(0, labels.shape[0], _RELABEL_PLANES):
            labels[z0:z0 + _RELABEL_PLANES] = lut[labels[z0:z0 + _RELABEL_PLANES]]
    counts = _label_counts(labels)
    n_vox = np.where(assigned > 0, counts[np.minimum(assigned, len(counts) - 1)], np.nan).astype(float)
    # the traced region: the bounding box of every kept label, plane by plane (the recorded ROI box of the
    # creation wizard does not contain the surfaces on the real exports, this does)
    any_xy = np.zeros(labels.shape[1:], bool)
    zs = []
    for z in range(labels.shape[0]):
        pl = labels[z] > 0
        if pl.any():
            zs.append(z)
            any_xy |= pl
    ys, xs = np.flatnonzero(any_xy.any(axis=1)), np.flatnonzero(any_xy.any(axis=0))
    bbox = ({"min_xyz": [int(xs[0]), int(ys[0]), int(zs[0])], "max_xyz": [int(xs[-1]), int(ys[-1]), int(zs[-1])]}
            if zs else None)
    # sanity: the label under each surface's centroid voxel (own label / 0 = cavity / another label)
    idx = surfaces[["xi", "yi", "zi"]].to_numpy()
    inside = ((idx >= 0) & (idx < geom["shape_xyz"])).all(axis=1)
    under = np.zeros(n_surf, np.int64)
    under[inside] = labels[idx[inside, 2], idx[inside, 1], idx[inside, 0]]
    has = assigned > 0
    per = pd.DataFrame({"label": assigned, "n_voxels_mask": n_vox, "mask_centroid_dist_um": dist},
                       index=surfaces.index)
    report = {
        "strategy": strategy, "match_tolerance_um": float(tol_um),
        "n_surfaces": int(n_surf), "n_surfaces_with_mask": int(has.sum()),
        "max_match_dist_um": float(np.nanmax(dist)) if has.any() else None,
        "median_match_dist_um": float(np.nanmedian(dist)) if has.any() else None,
        "n_values_in_channel": int(len(values)), "n_shared_values": int(len(shared)),
        "n_shared_split_by_components": int(n_split_components),
        "n_shared_split_by_nearest_centroid": int(n_split_voronoi),
        "n_orphan_values": int(len(orphans)), "n_orphan_voxels": n_orphan_voxels, "orphan_values": orphans[:50],
        "n_centroid_inside_own_label": int((has & (under == assigned)).sum()),
        "n_centroid_in_cavity": int((has & (under == 0)).sum()),
        "n_centroid_in_other_label": int((has & (under != 0) & (under != assigned)).sum()),
        "labels_shape_zyx": [int(v) for v in labels.shape], "labels_dtype": "int32",
        "n_labels_written": int(len(claimed_values)), "labels_bbox_xyz_voxel": bbox,
    }
    return labels, per, report


# ---- the frame check (protocol changelog entry 9 (d)(3)) --------------------------------------------
#
# Imaris stores the image extents as TEXT with six significant digits. On the 2600-pixel fields of the
# registered gonads 2600 x 0.108 = 280.8 is exact at that precision, so both xy axes read 0.108. On a field
# above 10000 um of stage position (the 2026-01-05 males: ExtMin0 16249.9) the last stored digit is 0.1 um,
# ExtMin is known only to +-0.05 um and the two xy axes read 0.10801 and 0.108. Spot positions are sub-voxel,
# so the frame is tested, not assumed: Imaris's own per-spot "Intensity Center" of the spot channel (the pixel
# value at the voxel holding the spot centre, an integer) and of the mask channel must equal the image at
# floor((p - ExtMin) / v). ``frame_check`` runs the rule declared in entry 9.

FRAME_PASS_FRACTION = 0.995      # entry 9 (d)(3): overall and in each third of the field in x and in y
FRAME_SEARCH_STEP_UM = 0.001     # the ExtMin search grid inside the text's rounding interval


def text_rounding_halfwidth(value: float, significant: int = 6) -> float:
    """Half the last digit of ``value`` printed with ``significant`` significant digits (C's %g, which
    Imaris uses for the extents): 0.05 for 16249.9, 0.005 for 8150.1 (printed 8150.10 before the zero is
    dropped)."""
    v = abs(float(value))
    if v == 0 or not np.isfinite(v):
        return 0.0
    return 0.5 * 10.0 ** (int(np.floor(np.log10(v))) - (significant - 1))


def read_nd2_channel_zyx(path, channel: int) -> np.ndarray:
    """One channel of an .nd2 as a (Z, Y, X) array in the file's dtype (T and P loops take index 0)."""
    import nd2

    with nd2.ND2File(str(path)) as f:
        sizes = dict(f.sizes)
        darr = f.to_dask()
    dims = list(sizes)
    sl = tuple(0 if d in ("T", "P") else (int(channel) if d == "C" else slice(None)) for d in dims)
    keep = [d for d in dims if d not in ("T", "P", "C")]
    arr = darr[sl]
    if "Z" not in keep:
        arr = arr[None]
        keep = ["Z"] + keep
    return np.asarray(arr.transpose([keep.index(d) for d in ("Z", "Y", "X")]))


FLOAT32_ULPS = 2                 # a spot this many float32 steps from a voxel boundary may read either side


def _values_at(img, idx, ok) -> np.ndarray:
    got = np.full(len(idx), np.nan)
    got[ok] = np.asarray(img[idx[ok, 2], idx[ok, 1], idx[ok, 0]], float)
    return got


def frame_agreement(pos_xyz, ext_min_xyz, voxel_xyz, shape_xyz, expected: dict, images: dict) -> dict:
    """For one candidate frame (ExtMin, voxel): the fraction of spots whose expected value (Imaris's per-spot
    Intensity Center, one array per key of ``expected``) equals ``images[key]`` (ZYX) at the voxel
    floor((p - ExtMin) / v), overall and in each third of the field in x and in y (by the spot's voxel under
    this frame). Imaris stores positions as float32, so a spot within ``FLOAT32_ULPS`` float32 steps of a voxel
    boundary also agrees when the voxel across that boundary holds the value (``agreement``; the plain floor
    comparison is ``agreement_strict``). On the registered gonads every strict disagreement lies within 0.003
    voxel of a boundary. A spot outside the image or without an expected value counts as a disagreement."""
    pos = np.asarray(pos_xyz, float).reshape(-1, 3)
    vox = np.asarray(voxel_xyz, float)
    rel = (pos - np.asarray(ext_min_xyz, float)) / vox
    idx = np.floor(rel).astype(np.int64)
    frac = rel - idx
    tol = FLOAT32_ULPS * np.spacing(np.abs(pos).astype(np.float32)).astype(float) / vox
    shape = np.asarray(shape_xyz, int)
    inside = ((idx >= 0) & (idx < shape)).all(axis=1)
    # the alternative voxels of a spot at a boundary: a step of -1 or +1 on each axis where it is that close
    steps = np.zeros_like(idx)
    steps[frac < tol] = -1
    steps[(1.0 - frac) < tol] = 1
    alts = []
    for mask_bits in range(1, 8):
        use = np.array([(mask_bits >> a) & 1 for a in range(3)], bool)
        alt = idx + np.where(use, steps, 0)
        valid = (steps[:, use] != 0).all(axis=1) & ((alt >= 0) & (alt < shape)).all(axis=1)
        if valid.any():
            alts.append((alt, valid))
    out = {"n_spots": int(len(pos)), "n_outside_image": int((~inside).sum()),
           "n_within_float32_of_a_boundary": int((steps != 0).any(axis=1).sum()), "channels": {}}
    thirds = {ax: np.clip(idx[:, i] * 3 // max(int(shape[i]), 1), 0, 2) for i, ax in ((0, "x"), (1, "y"))}
    for key, exp in expected.items():
        exp = np.asarray(exp, float)
        img = images[key]
        strict = np.isfinite(exp) & (_values_at(img, idx, inside) == exp)
        same = strict.copy()
        for alt, valid in alts:
            same |= valid & (_values_at(img, alt, valid) == exp)
        n = len(pos)
        rec = {"agreement": float(same.mean()) if n else float("nan"), "n_agree": int(same.sum()),
               "agreement_strict": float(strict.mean()) if n else float("nan"), "n_agree_strict": int(strict.sum())}
        for ax, t in thirds.items():
            rec[f"thirds_{ax}"] = [float(same[t == k].mean()) if (t == k).any() else None for k in range(3)]
            rec[f"thirds_{ax}_n"] = [int((t == k).sum()) for k in range(3)]
        out["channels"][key] = rec
    return out


def frame_passes(rec: dict, spot_key: str, mask_key: str | None, min_fraction: float = FRAME_PASS_FRACTION) -> bool:
    """The pass rule of entry 9 (d)(3): the spot channel's agreement at least ``min_fraction`` overall and in
    every third of x and of y that holds spots, and the mask channel's at least ``min_fraction`` overall."""
    s = rec["channels"][spot_key]
    ok = s["agreement"] >= min_fraction and all(v is None or v >= min_fraction
                                                for ax in ("x", "y") for v in s[f"thirds_{ax}"])
    if mask_key is not None and mask_key in rec["channels"]:
        ok = ok and rec["channels"][mask_key]["agreement"] >= min_fraction
    return bool(ok)


def frame_check(pos_xyz, geom: dict, expected: dict, images: dict, *, spot_key: str = "spot_channel",
                mask_key: str | None = "mask_channel", min_fraction: float = FRAME_PASS_FRACTION,
                step_um: float = FRAME_SEARCH_STEP_UM) -> dict:
    """The frame rule of protocol changelog entry 9 (d)(3), in order:

    1. ``stored``: ExtMin as stored and the stored-extent voxel ((ExtMax - ExtMin) / size), the reader's rule;
    2. ``rounded``: ExtMin as stored and that voxel rounded to 4 decimals (0.10801 -> 0.108);
    3. ``search``: per axis, ExtMin moved inside the rounding interval of its six-significant-digit text in
       steps of ``step_um``, for each of the two voxels, coordinate-wise (x, y, z, twice) to the highest
       spot-channel agreement;
    the first frame that passes (``frame_passes``) is ``chosen``; None when none does (a truth-side defect:
    the spot truth of the gonad is then not usable). Every candidate's agreement is returned."""
    ext = np.asarray(geom["ext_min_xyz"], float)
    shape = np.asarray(geom["shape_xyz"], int)
    v_stored = np.asarray(geom["voxel_xyz_um"], float)
    v_round = np.round(v_stored, 4)
    cands = [{"name": "stored", "ext_min_xyz": ext.tolist(), "voxel_xyz_um": v_stored.tolist()},
             {"name": "rounded", "ext_min_xyz": ext.tolist(), "voxel_xyz_um": v_round.tolist()}]
    results = []
    chosen = None
    for c in cands:
        rec = frame_agreement(pos_xyz, c["ext_min_xyz"], c["voxel_xyz_um"], shape, expected, images)
        rec.update(c, passes=frame_passes(rec, spot_key, mask_key, min_fraction))
        results.append(rec)
        if chosen is None and rec["passes"]:
            chosen = rec
    search = None
    if chosen is None:
        half = np.array([text_rounding_halfwidth(v) for v in ext])
        best_overall = None
        for vname, vox in (("stored", v_stored), ("rounded", v_round)):
            e = ext.copy()
            for _ in range(2):
                for ax in range(3):
                    n = int(round(half[ax] / step_um))
                    trial = []
                    for k in range(-n, n + 1):
                        e2 = e.copy()
                        e2[ax] = ext[ax] + k * step_um
                        r = frame_agreement(pos_xyz, e2, vox, shape, {spot_key: expected[spot_key]},
                                            {spot_key: images[spot_key]})["channels"][spot_key]
                        # ties on the float32-tolerant agreement go to the higher strict agreement: the edge of
                        # the agreeing plateau puts a spot on the wrong side of a voxel boundary (entry 10 (f))
                        trial.append((r["agreement"], r["agreement_strict"], -abs(k), k))
                    a, _strict, _neg, k = max(trial)
                    e[ax] = ext[ax] + k * step_um
            rec = frame_agreement(pos_xyz, e, vox, shape, expected, images)
            rec.update(name=f"search_{vname}_voxel", ext_min_xyz=e.tolist(), voxel_xyz_um=vox.tolist(),
                       ext_min_shift_um=(e - ext).round(6).tolist(), passes=frame_passes(rec, spot_key, mask_key,
                                                                                        min_fraction))
            results.append(rec)
            if best_overall is None or rec["channels"][spot_key]["agreement"] > \
                    best_overall["channels"][spot_key]["agreement"]:
                best_overall = rec
            if chosen is None and rec["passes"]:
                chosen = rec
        search = {"rounding_halfwidth_um": half.tolist(), "step_um": float(step_um),
                  "best": best_overall["name"] if best_overall else None}
    return {"rule": "entry 9 (d)(3): stored, then rounded voxel, then an ExtMin search inside the text's "
                    "rounding interval; the first frame passing on the spot and mask channels is chosen",
            "min_fraction": float(min_fraction), "candidates": results, "search": search,
            "chosen": None if chosen is None else {"name": chosen["name"], "ext_min_xyz": chosen["ext_min_xyz"],
                                                    "voxel_xyz_um": chosen["voxel_xyz_um"]}}



# ---- the frame check and the xlsx spot source of an export (protocol changelog entry 9) ---------------

#: entry 9 (e) V7 and (d)(5): the section 4.4 agreement must reach the lowest of the registered train and val
#: gonads (a minimum, entry 10 (b); benchmark_split.SECTION_4_4_MIN applies it to every addition)
TRUTH_ASSIGNMENT_MIN = 0.9935
XLSX_POSITION_TOL_UM = 0.005     # entry 9 (e) V3: surface position, xlsx against the .ims
XLSX_VOLUME_RTOL = 1e-5          # entry 9 (e) V3: surface volume, relative
XLSX_RADIUS_UM = 0.2             # entry 9 (e) V6: every spot radius (Diameter X and Y / 2)


def read_channel_zyx(f, geom: dict, channel: int) -> np.ndarray:
    """The image part of a DataSet channel in its own dtype, (Z, Y, X), read in z-chunks."""
    X, Y, Z = (int(v) for v in geom["shape_xyz"])
    ds = _channel_dataset(f, channel)
    out = np.zeros((Z, Y, X), ds.dtype)
    step = int(ds.chunks[0]) if ds.chunks else 8
    for z0 in range(0, Z, step):
        out[z0:min(z0 + step, Z)] = ds[z0:min(z0 + step, Z), :Y, :X]
    return out


def spot_source_channel(f) -> int | None:
    """The DataSet channel the first Spots object was detected on (``mSourceChannelIndex``)."""
    for _gname, g in _object_groups(_content_group(f), "Points"):
        m = re.search(r'mSourceChannelIndex="(\d+)"', read_creation_parameters(g))
        if m:
            return int(m.group(1))
    return None


def _slim_frame(res: dict) -> dict:
    """A frame-check result without the per-third spot counts (they stay in the candidates)."""
    return {k: v for k, v in res.items() if k != "candidates"} | {"candidates": res["candidates"]}


def _frame_check_ims(f, geom: dict, mask_channel: int | None, nd2) -> tuple[dict, dict]:
    """Entry 9 (d)(3) on the .ims's own Spots: their stored Intensity Center of the spot source channel against
    the nd2 pixel (the .ims channel when no nd2 is given) and of the mask channel against the .ims mask channel.
    Returns (the ``frame_check`` result with the images' source, the images for reuse)."""
    src = spot_source_channel(f)
    groups = [(n, g) for n, g in _object_groups(_content_group(f), "Points") if "Spot" in g]
    if src is None or not groups:
        return {"chosen": None, "untestable": "no Spots object or no source channel in its creation parameters"}, {}
    g = groups[0][1]
    sp = g["Spot"][:]
    ids = sp["ID"].astype(np.int64)
    pos = np.c_[sp["PositionX"], sp["PositionY"], sp["PositionZ"]].astype(float)
    wanted = {"spot_channel": ("Intensity Center", src)}
    if mask_channel is not None:
        wanted["mask_channel"] = ("Intensity Center", int(mask_channel))
    st = read_object_statistics(g, wanted)
    if "spot_channel" not in st:
        return {"chosen": None, "untestable": "the Spots object stores no Intensity Center of its source channel"}, {}
    expected = {k: v.reindex(ids).to_numpy() for k, v in st.items()}
    images = {"spot_channel": (read_nd2_channel_zyx(nd2, src) if nd2 is not None else read_channel_zyx(f, geom, src))}
    if tuple(images["spot_channel"].shape) != tuple(int(v) for v in geom["shape_xyz"][::-1]):
        raise ValueError(f"the spot channel image is {images['spot_channel'].shape}, the .ims image "
                         f"{tuple(int(v) for v in geom['shape_xyz'][::-1])} (ZYX)")
    if "mask_channel" in expected:
        images["mask_channel"] = read_channel_zyx(f, geom, int(mask_channel))
    res = frame_check(pos, geom, expected, images, mask_key="mask_channel" if "mask_channel" in expected else None)
    res.update(image_source=(f"nd2:{nd2}" if nd2 is not None else "ims"), spot_channel=int(src),
               mask_channel=None if mask_channel is None else int(mask_channel), spots_object=groups[0][0],
               n_spots=int(len(pos)))
    return res, images


def _surface_info(f) -> pd.DataFrame:
    """surface_id, Position (px, py, pz), centre of mass (cx, cy, cz) and volume of every surface, in the
    order of ``read_surfaces_table``."""
    rows = []
    for _gname, g in _object_groups(_content_group(f), "MegaSurfaces"):
        if "SurfaceModelInfo" not in g:
            continue
        smi = g["SurfaceModelInfo"][:]
        n = len(smi)
        ids = (g["SurfaceModel"][:]["ID"].astype(np.int64) if "SurfaceModel" in g and len(g["SurfaceModel"]) == n
               else np.arange(n, dtype=np.int64))
        names = smi.dtype.names
        pxyz = [smi[k].astype(float) if k in names else np.full(n, np.nan) for k in ("PositionX", "PositionY", "PositionZ")]
        rows.append(pd.DataFrame({"surface_id": ids, "px": pxyz[0], "py": pxyz[1], "pz": pxyz[2],
                                  "cx": smi["CenterOfMassX"].astype(float), "cy": smi["CenterOfMassY"].astype(float),
                                  "cz": smi["CenterOfMassZ"].astype(float), "volume": smi["Volume"].astype(float)}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["surface_id", "px", "py", "pz", "cx", "cy", "cz", "volume"])


def compare_xlsx_surfaces(xl: pd.DataFrame, info: pd.DataFrame, ims_mask_max) -> dict:
    """Entry 9 (e) V3: the export's surfaces against the .ims Surfaces, one to one by ID (by nearest position
    when the ID sets differ): position within ``XLSX_POSITION_TOL_UM`` (the .ims Position, else its centre of
    mass), volume within ``XLSX_VOLUME_RTOL`` relative, and the mask channel's Intensity Max equal."""
    from scipy.spatial import cKDTree

    info = info.reset_index(drop=True).copy()
    info["imax"] = np.asarray(ims_mask_max, float) if ims_mask_max is not None else np.nan
    pos_cols = ["px", "py", "pz"] if info[["px", "py", "pz"]].notna().all().all() else ["cx", "cy", "cz"]
    out = {"n_xlsx": int(len(xl)), "n_ims": int(len(info)), "position_fields": "Position" if pos_cols[0] == "px"
           else "CenterOfMass"}
    if set(xl["ID"].astype(int)) == set(info["surface_id"].astype(int)):
        out["matched_by"] = "ID"
        m = info.set_index("surface_id").loc[xl["ID"].astype(int).to_numpy()]
    else:
        out["matched_by"] = "position (the ID sets differ)"
        out["ids_only_in_xlsx"] = sorted(set(xl["ID"].astype(int)) - set(info["surface_id"].astype(int)))[:20]
        out["ids_only_in_ims"] = sorted(set(info["surface_id"].astype(int)) - set(xl["ID"].astype(int)))[:20]
        if len(info) == 0 or len(xl) == 0:
            out.update(passes=False, n_matched=0)
            return out
        d, j = cKDTree(info[pos_cols].to_numpy(float)).query(xl[["x", "y", "z"]].to_numpy(float))
        if len(set(j.tolist())) != len(j):
            out.update(passes=False, n_matched=int(len(set(j.tolist()))), reason="not one to one")
            return out
        m = info.iloc[j]
    dpos = np.linalg.norm(xl[["x", "y", "z"]].to_numpy(float) - m[pos_cols].to_numpy(float), axis=1)
    dvol = np.abs(xl["volume"].to_numpy(float) - m["volume"].to_numpy(float)) / np.maximum(np.abs(m["volume"].to_numpy(float)), 1e-12)
    same_imax = np.round(xl["imax_mask"].to_numpy(float)) == np.round(m["imax"].to_numpy(float))
    out.update(n_matched=int(len(m)), max_position_diff_um=float(dpos.max()) if len(dpos) else None,
               max_volume_rel_diff=float(dvol.max()) if len(dvol) else None,
               n_mask_intensity_max_equal=int(same_imax.sum()))
    out["passes"] = bool(len(xl) == len(info) == len(m) and len(m) > 0 and dpos.max() <= XLSX_POSITION_TOL_UM
                         and dvol.max() <= XLSX_VOLUME_RTOL and same_imax.all())
    return out


def truth_assignment_agreement(labels: np.ndarray, zi, yi, xi, coloc) -> dict:
    """Protocol section 4.4 for one gonad: the geometric assignment of the truth spots (the mask label at the
    centre voxel, else the 3 x 3 x 3 majority, fragments below the protocol's size removed) against the lab's
    mask-channel value ``coloc``; the rule of scripts/check_truth_assignment.py."""
    from ..benchmark.spots import TRUTH_MASK_MIN_VOXELS, assign_to_masks, truth_mask_lut

    lut, _keep, _sizes, _frag = truth_mask_lut(labels, TRUTH_MASK_MIN_VOXELS)
    geo = assign_to_masks(labels, np.asarray(zi, np.int64), np.asarray(yi, np.int64), np.asarray(xi, np.int64), lut=lut)
    col = np.nan_to_num(np.asarray(coloc, float), nan=0.0).astype(np.int64)
    agree = geo == col
    return {"n_spots": int(len(col)), "n_agree": int(agree.sum()),
            "agreement": float(agree.mean()) if len(col) else float("nan")}


def xlsx_spots_table(xl: pd.DataFrame, geom: dict, our_vox) -> pd.DataFrame:
    """The truth spots table of an xlsx spot source (entry 9 (e)): placed with the .ims ExtMin and voxel of
    ``geom`` (floor((p - ExtMin) / v), then index times the nd2 voxel), the xlsx positions as the Imaris
    columns, the xlsx Intensity Max of the mask channel as ``imaris_mask_label``, radius half the Diameter X,
    no Quality and no class (the export holds no per-spot class)."""
    pos = xl[["x", "y", "z"]].to_numpy(float)
    idx = imaris_to_voxel_index(pos, geom)
    ours = voxel_index_to_our_um(idx, our_vox)
    n = len(xl)
    df = pd.DataFrame({
        "spot_id": xl["ID"].astype(np.int64).to_numpy(), "x_um": ours[:, 0], "y_um": ours[:, 1], "z_um": ours[:, 2],
        "xi": idx[:, 0], "yi": idx[:, 1], "zi": idx[:, 2], "radius_um": xl["diameter_x"].to_numpy(float) / 2.0,
        "quality": np.full(n, np.nan), "spot_class": [""] * n,
        "imaris_x_um": pos[:, 0], "imaris_y_um": pos[:, 1], "imaris_z_um": pos[:, 2],
        "imaris_mask_label": _nullable_int(xl["imax_mask"].to_numpy(float)), "spots_object": "Spots 1 (xlsx)"})
    return df[SPOT_COLUMNS + SPOT_EXTRA_COLUMNS]


# ---- the export --------------------------------------------------------------------------------------

def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if np.isfinite(v) else None
    if isinstance(o, Path):
        return str(o)
    return o


def file_sha256_path(path, chunk_bytes: int = 1 << 24) -> str:
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_bytes), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_labels(labels: np.ndarray, path: Path, spacing_zyx) -> None:
    """Write the truth label image exactly as the pipeline writes its label images (int32 OME-TIFF with
    PhysicalSize tags, zlib), through the same function; raise when nothing was written."""
    from ..pipeline import _save_labels

    _save_labels(labels, path, spacing_zyx)
    if not Path(long_path(path)).exists():
        raise OSError(f"the label image could not be written: {path}")


def export_imaris_truth(ims: str | Path, out_dir: str | Path, *, nd2: str | Path | None = None,
                        nd2_voxel_xyz_um=None, masks: bool = True, tol_um: float = _TOL_UM,
                        frame_check_on: bool = False, spots_xlsx: str | Path | None = None,
                        xlsx_expect: dict | None = None, allow_test: bool = False) -> dict:
    """Write the truth file contract for one Imaris .ims into ``out_dir`` and return the params dict.

    ``nd2`` (or ``nd2_voxel_xyz_um`` as (dx, dy, dz) um) defines OUR frame; without either a sibling
    ``<stem>.nd2`` is used, else Imaris's own voxel size. ``masks=False`` skips the mask channel work
    (no label image; the nuclei ``label`` then comes from the per-surface intensity statistic of the
    mask channel when that statistic exists).

    ``frame_check_on`` runs the frame test of protocol changelog entry 9 (d)(3) (``frame_check``) on the .ims's
    own spots, with the spot channel from ``nd2`` when it is given; the frame it chooses converts every Imaris
    position (``ext_min_xyz_um`` and ``voxel_imaris_xyz_um`` of the params are then that frame; the stored ones
    are kept under ``*_stored``), and ``frame_check`` in the params records every candidate. When no frame
    passes, the stored frame is kept and ``frame_check.chosen`` is None (the manifest then marks the spot truth
    unusable). ``spots_xlsx`` (entry 9 (e)) takes the spots from an Imaris statistics export instead of the
    .ims's saved Spots object, verifies them (criteria V1 to V7, ``spots_source.verification``; ``xlsx_expect``
    gives the recorded sha256 and Overall totals) and implies the frame check.

    With the frame check or the xlsx source, the files are guarded before anything is opened
    (``benchmark_split.guard_truth_read``, protocol changelog entry 10): a registered gonad's project passes the
    role guard, so a test gonad is refused unless ``allow_test`` (the pre-flight of docs/FINAL_SCORING.md), and
    any other sealed file (a 2025-11-05 export, an aggregate table) is refused."""
    import h5py

    ims = Path(ims)
    frame_check_on = bool(frame_check_on or spots_xlsx is not None)
    if frame_check_on:
        from .benchmark_split import guard_truth_read

        guard_truth_read([ims, nd2, spots_xlsx], allow_test=allow_test)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = truth_paths(out_dir, ims.stem)
    warnings: list[str] = []
    frame_rep, xlsx_rep = None, None
    with h5py.File(str(ims), "r") as f:
        geom = read_geometry(f)
        geom_stored = dict(geom)
        our_vox, vox_source = resolve_our_voxel(ims, nd2, nd2_voxel_xyz_um, geom)
        if vox_source == "imaris_fallback":
            warnings.append("no nd2 voxel size: our frame uses Imaris's rounded voxel size (the frames coincide)")
        surfaces, surf_reports = read_surfaces_table(f, geom, our_vox)
        mask_rep: dict = {"requested": bool(masks), "channel": None, "channel_name": None, "identification": None,
                          "intensity_max_statistic_present": False, "labels_file": None, "reason_no_labels": None,
                          "label_source": None}
        channel, values_unique, stat_label, labels = None, True, None, None
        if len(surfaces):
            channel, ident = identify_mask_channel_by_probe(f, geom, surfaces)
            mask_rep["identification"] = ident
            if channel is not None:
                values_unique = bool(ident.get("values_unique", True))
                mask_rep["channel"] = int(channel)
                mask_rep["channel_name"] = geom["channel_names"][channel]
                stat_label = surface_mask_statistic(f, geom, channel, surfaces)
                mask_rep["intensity_max_statistic_present"] = stat_label is not None
            elif masks:
                mask_rep["reason_no_labels"] = "no DataSet channel is a constant-per-surface mask (see identification)"
                warnings.append(mask_rep["reason_no_labels"])
        else:
            mask_rep["reason_no_labels"] = "the file holds no surfaces"
        images = {}
        if frame_check_on:
            frame_rep, images = _frame_check_ims(f, geom, channel, nd2)
            chosen = frame_rep.get("chosen")
            if chosen is None:
                warnings.append("no frame passed the Intensity Center test of protocol changelog entry 9 (d)(3); "
                                "the stored frame is kept and the spot truth is not usable")
            elif chosen["name"] != "stored":
                geom = {**geom, "ext_min_xyz": np.asarray(chosen["ext_min_xyz"], float),
                        "voxel_xyz_um": np.asarray(chosen["voxel_xyz_um"], float)}
                surfaces_new, surf_reports = read_surfaces_table(f, geom, our_vox)
                for c in ("label", "n_voxels_mask", "mask_centroid_dist_um"):
                    surfaces_new[c] = surfaces[c].to_numpy()
                surfaces = surfaces_new
                warnings.append(f"frame {chosen['name']} used (entry 9 (d)(3)): ExtMin {chosen['ext_min_xyz']}, "
                                f"voxel {chosen['voxel_xyz_um']}")
        if spots_xlsx is not None:
            from .imaris_xlsx import read_export_objects

            src_ch = frame_rep.get("spot_channel") if frame_rep else spot_source_channel(f)
            xl = read_export_objects(spots_xlsx, spot_channel=int(src_ch if src_ch is not None else 2),
                                     mask_channel=int(channel if channel is not None else 3))
            xlsx_rep = {"source": "xlsx", "file": str(spots_xlsx), "verification": {}}
            ver = xlsx_rep["verification"]
            exp = dict(xlsx_expect or {})
            sha = file_sha256_path(spots_xlsx)
            ver["V1_sha256"] = {"sha256": sha, "expected": exp.get("sha256"), "passes": bool(exp.get("sha256")) and sha == exp.get("sha256")}
            ov = xl["overall"]
            n_img = int(np.prod([int(v) for v in geom["shape_xyz"]]))
            ver["V2_totals"] = {"overall": ov, "position_sheet_spots": int(len(xl["spots"])),
                                "position_sheet_surfaces": int(len(xl["surfaces"])),
                                "expected_spots": exp.get("n_spots"), "expected_surfaces": exp.get("n_surfaces"),
                                "image_voxels_match_ims": ov.get("n_voxels_image") == n_img,
                                "passes": bool(exp.get("n_spots") is not None and exp.get("n_surfaces") is not None
                                               and ov.get("n_spots") == exp["n_spots"] == len(xl["spots"])
                                               and ov.get("n_surfaces") == exp["n_surfaces"] == len(xl["surfaces"]))}
            ver["V3_surfaces"] = compare_xlsx_surfaces(xl["surfaces"], _surface_info(f), stat_label)
            xs = xl["spots"]
            vchk = frame_check(xs[["x", "y", "z"]].to_numpy(float), geom_stored,
                               {"spot_channel": xs["ic_spot"].to_numpy(float), "mask_channel": xs["ic_mask"].to_numpy(float)},
                               images, mask_key="mask_channel" if "mask_channel" in images else None) if images else \
                {"chosen": None, "untestable": "no images (the .ims frame check could not run)"}
            ver["V5_intensity_center"] = {"chosen": vchk.get("chosen"), "candidates": vchk.get("candidates"),
                                          "search": vchk.get("search"),
                                          "same_frame_as_ims_spots": bool(vchk.get("chosen") and frame_rep and frame_rep.get("chosen")
                                                                          and vchk["chosen"]["name"] == frame_rep["chosen"]["name"]),
                                          "passes": vchk.get("chosen") is not None}
            xgeom = geom if vchk.get("chosen") is None else {**geom_stored,
                                                              "ext_min_xyz": np.asarray(vchk["chosen"]["ext_min_xyz"], float),
                                                              "voxel_xyz_um": np.asarray(vchk["chosen"]["voxel_xyz_um"], float)}
            xidx = imaris_to_voxel_index(xs[["x", "y", "z"]].to_numpy(float), xgeom)
            inside = ((xidx >= 0) & (xidx < np.asarray(geom["shape_xyz"]))).all(axis=1)
            ver["V4_inside_extent"] = {"n_outside": int((~inside).sum()), "passes": bool(inside.all() and len(xs))}
            r = np.concatenate([xs["diameter_x"].to_numpy(float), xs["diameter_y"].to_numpy(float)]) / 2.0
            ver["V6_radius"] = {"min_um": float(np.nanmin(r)) if len(r) else None, "max_um": float(np.nanmax(r)) if len(r) else None,
                                "diameter_z_values_um": sorted({round(float(v), 6) for v in xs["diameter_z"]})[:10],
                                "passes": bool(len(r) and np.all(np.abs(r - XLSX_RADIUS_UM) <= 1e-6))}
            xlsx_rep["frame_used"] = {"name": (vchk.get("chosen") or {}).get("name"),
                                      "ext_min_xyz_um": xgeom["ext_min_xyz"], "voxel_imaris_xyz_um": xgeom["voxel_xyz_um"]}
            xlsx_spots = xlsx_spots_table(xs, xgeom, our_vox)
        del images
        if masks and channel is not None:
            labels, per, rec = recover_surface_masks(f, geom, channel, surfaces, values_unique=values_unique,
                                                     tol_um=tol_um)
            surfaces["label"] = per["label"].to_numpy()
            surfaces["n_voxels_mask"] = per["n_voxels_mask"].to_numpy()
            surfaces["mask_centroid_dist_um"] = per["mask_centroid_dist_um"].to_numpy()
            mask_rep.update(rec)
            mask_rep["label_source"] = "mask channel blobs matched to the surface centres of mass"
            if stat_label is not None:
                lab = surfaces["label"].to_numpy()
                same = (lab > 0) & (np.round(np.nan_to_num(stat_label, nan=-1.0)) == lab)
                mask_rep["n_labels_equal_to_intensity_statistic"] = int(same.sum())
                mask_rep["n_labels_differing_from_intensity_statistic"] = int(((lab > 0) & ~same).sum())
            if rec["n_surfaces_with_mask"] == 0:
                labels = None
                mask_rep["reason_no_labels"] = "no surface matched a mask blob within the tolerance"
                warnings.append(mask_rep["reason_no_labels"])
            elif rec["n_surfaces_with_mask"] < len(surfaces):
                warnings.append(f"{len(surfaces) - rec['n_surfaces_with_mask']} of {len(surfaces)} surfaces have no mask")
        elif channel is not None and stat_label is not None:
            surfaces["label"] = np.round(np.nan_to_num(stat_label, nan=0.0)).astype(np.int64)
            mask_rep["label_source"] = (f"Intensity Max statistic of channel {channel} per surface "
                                        f"(masks not read; unverified)")
        if not masks:
            mask_rep["reason_no_labels"] = "masks not requested (--no-masks)"
        spots, spot_reports = read_spots_table(f, geom, our_vox, mask_channel=channel)
        if xlsx_rep is not None:
            from scipy.spatial import cKDTree

            ver = xlsx_rep["verification"]
            if labels is not None:
                ta = truth_assignment_agreement(labels, xlsx_spots["zi"], xlsx_spots["yi"], xlsx_spots["xi"],
                                                xlsx_spots["imaris_mask_label"].astype("float").to_numpy())
                ver["V7_truth_assignment"] = {**ta, "minimum": TRUTH_ASSIGNMENT_MIN,
                                              "passes": bool(ta["agreement"] >= TRUTH_ASSIGNMENT_MIN)}
            else:
                ver["V7_truth_assignment"] = {"passes": False, "reason": "no label image"}
            ims_pos = spots[["imaris_x_um", "imaris_y_um", "imaris_z_um"]].to_numpy(float)
            if len(ims_pos) and len(xlsx_spots):
                d, _j = cKDTree(ims_pos).query(xlsx_spots[["imaris_x_um", "imaris_y_um", "imaris_z_um"]].to_numpy(float))
                frac = float((d <= XLSX_POSITION_TOL_UM).mean())
            else:
                frac = None
            xlsx_rep["reported_not_a_criterion"] = {
                "n_ims_spots_saved": int(len(spots)),
                "fraction_of_xlsx_spots_within_0.005_um_of_a_saved_ims_spot": frac}
            ver["passed"] = bool(all(v.get("passes") for k, v in ver.items() if k.startswith("V")))
            xlsx_rep["criteria"] = ("protocol changelog entry 9 (e): V1 sha256 as recorded; V2 Overall totals as "
                                    "recorded; V3 surfaces equal the .ims Surfaces one to one; V4 every spot inside "
                                    "the .ims extent; V5 the Intensity Center test of (d)(3) on the export's own "
                                    "values; V6 every radius 0.2 um; V7 the section 4.4 agreement at least "
                                    f"{TRUTH_ASSIGNMENT_MIN}")
            if not ver["passed"]:
                warnings.append("the xlsx spots failed the verification of entry 9 (e): the gonad is nucleus truth only")
            spot_reports = [{**r, "used_as_spot_truth": False} for r in spot_reports]
            spots = xlsx_spots

    # ---- write ----
    spots[SPOT_COLUMNS + SPOT_EXTRA_COLUMNS].to_csv(long_path(paths["spots"]), index=False)
    nuc_out = surfaces[NUCLEI_COLUMNS + NUCLEI_EXTRA_COLUMNS].copy()
    nuc_out["n_voxels_mask"] = _nullable_int(nuc_out["n_voxels_mask"].to_numpy())
    nuc_out.to_csv(long_path(paths["nuclei"]), index=False)
    if labels is not None:
        _write_labels(labels, paths["labels"], spacing_zyx=(our_vox[2], our_vox[1], our_vox[0]))
        mask_rep["labels_file"] = paths["labels"].name
        del labels
    roi_boxes = [{"object": r["name"], **b} for r in surf_reports + spot_reports for b in r["roi_boxes"]]
    roi_check = None
    if roi_boxes and len(surfaces) and roi_boxes[0]["min_xyz_voxel"] and roi_boxes[0]["max_xyz_voxel"]:
        box = roi_boxes[0]
        idx = surfaces[["xi", "yi", "zi"]].to_numpy()
        inside = ((idx >= np.array(box["min_xyz_voxel"])) & (idx <= np.array(box["max_xyz_voxel"]))).all(axis=1)
        roi_check = {"box": box["name"], "object": box["object"], "n_surface_centroids_inside": int(inside.sum()),
                     "n_surfaces": int(len(surfaces))}
        if not inside.all():
            warnings.append(f"only {int(inside.sum())}/{len(surfaces)} surface centroids lie inside the recorded "
                            f"ROI box; treat the box as a record of the wizard, not a spatial filter")
    classes_spots = sorted({c for s in spots["spot_class"] for c in str(s).split(_CLASS_SEP) if c})
    classes_surf = sorted({c for s in surfaces["surface_class"] for c in str(s).split(_CLASS_SEP) if c})
    params = {
        "ims_file": str(ims), "stem": ims.stem, "imaris_version": geom["imaris_version"],
        "image_name": geom["image_name"], "imaris_source_file": geom["source_file"],
        "nd2_file": str(nd2) if nd2 is not None else None, "our_voxel_source": vox_source,
        "voxel_imaris_xyz_um": (xlsx_rep["frame_used"]["voxel_imaris_xyz_um"] if xlsx_rep else geom["voxel_xyz_um"]),
        "voxel_ours_xyz_um": our_vox,
        "ext_min_xyz_um": (xlsx_rep["frame_used"]["ext_min_xyz_um"] if xlsx_rep else geom["ext_min_xyz"]),
        "ext_max_xyz_um": geom["ext_max_xyz"],
        "voxel_imaris_xyz_um_stored": geom_stored["voxel_xyz_um"], "ext_min_xyz_um_stored": geom_stored["ext_min_xyz"],
        "shape_xyz": geom["shape_xyz"], "data_shape_zyx_padded": geom["data_shape_zyx"],
        "data_dtype": geom["data_dtype"], "channel_names": geom["channel_names"], "n_channels": geom["n_channels"],
        "frame_note": "xi,yi,zi = floor((imaris - ext_min) / voxel_imaris); x_um,y_um,z_um = index * voxel_ours",
        "objects": surf_reports + spot_reports,
        "creation_parameters": {r["name"]: r["creation_parameters"] for r in surf_reports + spot_reports},
        "roi_boxes": roi_boxes, "roi_check": roi_check,
        "n_spots": int(len(spots)), "n_surfaces": int(len(surfaces)),
        "spot_classes_present": classes_spots, "surface_classes_present": classes_surf,
        "region_growing": bool(any(r.get("region_growing") for r in spot_reports)),
        "quality_present": bool(any(r.get("quality_present") for r in spot_reports)),
        "spot_statistics_present": sorted({n for r in spot_reports for n in r["statistics"]}),
        "mask": mask_rep,
        "files": {k: p.name for k, p in paths.items() if k != "labels" or mask_rep["labels_file"]},
        "warnings": warnings,
        "spots_source": xlsx_rep if xlsx_rep is not None else {"source": "ims"},
    }
    if frame_rep is not None:
        params["frame_check"] = frame_rep
    try:
        from .benchmark_split import creation_parameter_summary

        params["creation_parameter_summary"] = creation_parameter_summary(
            params["creation_parameters"].get(next((r["name"] for r in spot_reports), ""), None),
            params["creation_parameters"].get(next((r["name"] for r in surf_reports), ""), None))
    except Exception as e:  # noqa: BLE001 - a summary, never a reason to lose the export
        params["creation_parameter_summary"] = {"error": f"{type(e).__name__}: {e}"}
    for r in params["objects"]:
        r.pop("creation_parameters", None)   # kept once, under creation_parameters
    with open(long_path(paths["params"]), "w", encoding="utf-8") as fh:
        json.dump(_jsonable(params), fh, indent=1)
    for w in warnings:
        log.warning("%s: %s", ims.name, w)
    return params
