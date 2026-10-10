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


def _write_labels(labels: np.ndarray, path: Path, spacing_zyx) -> None:
    """Write the truth label image exactly as the pipeline writes its label images (int32 OME-TIFF with
    PhysicalSize tags, zlib), through the same function; raise when nothing was written."""
    from ..pipeline import _save_labels

    _save_labels(labels, path, spacing_zyx)
    if not Path(long_path(path)).exists():
        raise OSError(f"the label image could not be written: {path}")


def export_imaris_truth(ims: str | Path, out_dir: str | Path, *, nd2: str | Path | None = None,
                        nd2_voxel_xyz_um=None, masks: bool = True, tol_um: float = _TOL_UM) -> dict:
    """Write the truth file contract for one Imaris .ims into ``out_dir`` and return the params dict.

    ``nd2`` (or ``nd2_voxel_xyz_um`` as (dx, dy, dz) um) defines OUR frame; without either a sibling
    ``<stem>.nd2`` is used, else Imaris's own voxel size. ``masks=False`` skips the mask channel work
    (no label image; the nuclei ``label`` then comes from the per-surface intensity statistic of the
    mask channel when that statistic exists)."""
    import h5py

    ims = Path(ims)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = truth_paths(out_dir, ims.stem)
    warnings: list[str] = []
    with h5py.File(str(ims), "r") as f:
        geom = read_geometry(f)
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
        "voxel_imaris_xyz_um": geom["voxel_xyz_um"], "voxel_ours_xyz_um": our_vox,
        "ext_min_xyz_um": geom["ext_min_xyz"], "ext_max_xyz_um": geom["ext_max_xyz"],
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
    }
    for r in params["objects"]:
        r.pop("creation_parameters", None)   # kept once, under creation_parameters
    with open(long_path(paths["params"]), "w", encoding="utf-8") as fh:
        json.dump(_jsonable(params), fh, indent=1)
    for w in warnings:
        log.warning("%s: %s", ims.name, w)
    return params
