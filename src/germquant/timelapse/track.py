"""The tracking exchange of a time-lapse series (docs/DESIGN_timelapse.md section 5, docs/TIMELAPSE.md
"Tracking"): the per-timepoint nucleus labels go out to a tracker as files, the tracker's answer comes back
as files, and the result is two series tables.

* ``export_cellacdc`` (``germquant track-export``): the germline nucleus labels of every timepoint, placed on
  the registered grid by the whole-voxel rounding of the shift table, cropped to the union bounding box of the
  germline across timepoints and sub-sampled by a stride, as a Cell-ACDC ``Position_1/Images`` folder
  (``<stem>_segm.npz`` with a (T, Z, Y, X) integer array, ``<stem>_<channel>.tif`` with the central_element
  channel on the same grid, ``<stem>_metadata.csv``) plus ``export.json`` (the crop box, the stride, the
  shifts applied and the label id map of every timepoint).
* ``link_series`` (``germquant track-link``): links the nuclei across timepoints with one of the backends
  ``overlap`` (Cell-ACDC's own overlap tracker, headless), ``overlap_inhouse`` (the same majority-overlap rule
  in NumPy, used when Cell-ACDC is not importable), ``laptrack`` (LAP linking of the registered centroids with
  gap closing, splitting and merging off) or ``none`` (no linking: the points table only), in the frame of the
  rigid registration (one shift per timepoint) or of the piecewise one (timelapse.piecewise: one displacement
  per nucleus, interpolated along the gonad axis; the overlap backends then move each label by its own shift).
* ``import_tracks`` (``germquant track-import``): a TrackMate spots CSV or XML, an Imaris points table or a
  Cell-ACDC tracked ``segm.npz`` becomes the same two tables with ``link_source`` ``imported:<tool>``.
* ``track_strips`` (``germquant track-strips``): one PNG per track with a tile per timepoint, the
  non-programmer's check that an identity did not jump.

The tables (schema.SERIES_TABLES, append-only): ``<stem>__track_points.csv``, one row per nucleus per
timepoint with ``time_track_id``, the raw and the registered centroid, the step from the track's previous
point and the nucleus measurements joined on (image_id, nucleus_id); ``<stem>__tracks.csv``, one row per track
with its length, gaps, speed, signed axis velocity (sign aligned to t0, see ``align_axis_sign``), SYP intensity
slopes raw and median-normalised, and the SC length change. Registered coordinates are raw minus the
cumulative shift of the timepoint; nothing is ever resampled, and a tracker change never reruns segmentation.
"""
from __future__ import annotations

import json
import logging
import math
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .. import provenance, schema
from ..fsutil import long_path
from .layout import (
    SeriesLayout,
    SeriesLayoutError,
    find_series,
    load_shifts,
    object_table,
    voxel_shifts,
)
from .series import write_series_table

log = logging.getLogger(__name__)

# code defaults of every timelapse.track.* key (read with cfg.get(key, default); no config file is edited)
DEFAULTS = {
    "export_stride": 2,          # xy sub-sampling of the exported labels and image
    "export_crop": "germline",   # germline (in_germline labels, union bounding box) | all (every label, whole grid)
    "export_max_gb": 2.0,        # warn above this estimated size of the exported label array
    "max_dist_um": 3.0,          # laptrack: largest registered centroid displacement linked between frames
    "gap_frames": 2,             # laptrack: gap_closing_max_frame_count (2 bridges one missing timepoint)
    "min_overlap": 0.4,          # overlap backends: IoA = intersection / previous object's volume (Cell-ACDC's IoA_thresh)
    "link_stride": 2,            # overlap backends: xy stride of the label grid they run on
    "strip_size_um": 8.0,        # strips: tile width in microns
    "strip_planes": 3,           # strips: planes max-projected around the centroid
    "strip_tile_px": 64,         # strips: tile side in output pixels (nearest-neighbour upscale)
    "registration": "rigid",     # track-link: rigid (<stem>__register_shifts.csv) | piecewise (<stem>__register_piecewise.csv)
}
BACKENDS = ("overlap", "overlap_inhouse", "laptrack", "none")
REGISTRATIONS = ("rigid", "piecewise")
IMPORT_SOURCES = ("trackmate", "imaris", "cellacdc")
LINK_RECORD = "track_link"           # <stem>__track_link.json
EXPORT_RECORD = "export.json"
POSITION_FOLDER = ("Position_1", "Images")


def cfg_get(cfg, key: str, default):
    """``cfg.get(key, default)`` for a Config, a plain (possibly nested) dict or None."""
    if cfg is None:
        return default
    if isinstance(cfg, dict):
        node = cfg
        for part in key.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return default if node is None else node
    try:
        v = cfg.get(key, default)
    except Exception:  # noqa: BLE001 - a Config-like without a dotted get
        return default
    return default if v is None else v


# =====================================================================================================
# the registered, cropped label grid (shared by the export, the overlap backends and the Cell-ACDC import)
# =====================================================================================================

def _keep_ids(layout: SeriesLayout, t: int, crop: str) -> np.ndarray | None:
    """The nucleus ids a frame keeps: the in_germline ones for ``crop="germline"`` (every id when the table has
    no such column), None (= all labels) for ``crop="all"``."""
    if crop == "all":
        return None
    df = layout.read_nuclei(t)
    if "in_germline" not in df.columns or not len(df):
        return None
    from .layout import _truthy

    return df.loc[df["in_germline"].map(_truthy), "nucleus_id"].astype(int).to_numpy()


def _lut(ids: np.ndarray | None, max_label: int) -> np.ndarray | None:
    if ids is None:
        return None
    lut = np.zeros(max(int(max_label), int(ids.max()) if ids.size else 0) + 1, dtype=bool)
    lut[ids[ids > 0]] = True
    return lut


def _mask_ids(lab: np.ndarray, ids: np.ndarray | None) -> np.ndarray:
    """``lab`` with every label not in ``ids`` set to 0 (a copy only when something is dropped)."""
    if ids is None:
        return lab
    lut = _lut(ids, int(lab.max()) if lab.size else 0)
    return np.where(lut[lab], lab, 0)


def _bbox(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    """(start, stop) voxel corners of the True region of a (Z, Y, X) mask, or None when empty."""
    if not mask.any():
        return None
    z = np.flatnonzero(mask.any(axis=(1, 2)))
    y = np.flatnonzero(mask.any(axis=(0, 2)))
    x = np.flatnonzero(mask.any(axis=(0, 1)))
    return np.array([z[0], y[0], x[0]]), np.array([z[-1] + 1, y[-1] + 1, x[-1] + 1])


def crop_box(layout: SeriesLayout, vox_shifts: dict[int, np.ndarray], crop: str = "germline") -> dict:
    """The union bounding box, on the registered grid (t0's raw voxel grid), of the labels every timepoint
    keeps: ``{"start": [z0, y0, x0], "stop": [z1, y1, x1], "shape_raw": [Z, Y, X], "n_labels": {t: n}}``.
    ``crop="all"`` takes the whole grid of every frame (shifted); ``"germline"`` reads each label image once
    and takes the in_germline labels' extent (the export's second pass reads them again; the files are
    read twice so that no frame is held in memory)."""
    start = None
    stop = None
    n_labels: dict[int, int] = {}
    shape = None
    for t in layout.timepoints:
        sv = vox_shifts.get(int(t), np.zeros(3, int))
        if crop == "all":
            shp = np.array(layout.labels_shape(t))
            shape = shp if shape is None else shape
            lo, hi = -sv, shp - sv
            n_labels[int(t)] = int(len(layout.read_nuclei(t)))
        else:
            lab = layout.read_labels(t)
            shape = np.array(lab.shape) if shape is None else shape
            ids = _keep_ids(layout, t, crop)
            keep = _mask_ids(lab, ids)
            bb = _bbox(keep > 0)
            n_labels[int(t)] = int(len(ids)) if ids is not None else int(len(layout.read_nuclei(t)))
            del lab, keep
            if bb is None:
                continue
            lo, hi = bb[0] - sv, bb[1] - sv
        start = lo if start is None else np.minimum(start, lo)
        stop = hi if stop is None else np.maximum(stop, hi)
    if start is None:
        raise SeriesLayoutError(f"{layout.stem}: no label to export (every timepoint is empty for crop={crop!r})")
    return {"start": [int(v) for v in start], "stop": [int(v) for v in stop],
            "shape_raw": [int(v) for v in shape], "n_labels": n_labels}


def place_frame(lab: np.ndarray, sv: np.ndarray, box: dict, stride: int, ids: np.ndarray | None = None,
                dtype=None) -> np.ndarray:
    """``lab`` (a raw-frame (Z, Y, X) label image) placed on the registered grid: registered voxel r holds the
    raw voxel r + sv (sv = the frame's whole-voxel cumulative shift), cropped to ``box`` and sub-sampled by
    ``stride`` in y and x (z is kept whole). Labels not in ``ids`` are dropped; voxels that fall outside the raw
    grid are 0. The same function serves the export, the overlap backends and the Cell-ACDC import, so the
    three always agree on the grid."""
    s = max(1, int(stride))
    z0, y0, x0 = box["start"]
    z1, y1, x1 = box["stop"]
    zi = np.arange(z0, z1) + int(sv[0])
    yi = np.arange(y0, y1, s) + int(sv[1])
    xi = np.arange(x0, x1, s) + int(sv[2])
    Z, Y, X = lab.shape
    vz, vy, vx = (zi >= 0) & (zi < Z), (yi >= 0) & (yi < Y), (xi >= 0) & (xi < X)
    out = np.zeros((len(zi), len(yi), len(xi)), dtype=dtype or lab.dtype)
    if vz.any() and vy.any() and vx.any():
        sub = lab[np.ix_(zi[vz], yi[vy], xi[vx])]
        out[np.ix_(vz, vy, vx)] = _mask_ids(sub, ids)
    return out


def registered_frames(layout: SeriesLayout, vox_shifts: dict[int, np.ndarray] | None, box: dict, stride: int,
                      crop: str = "germline", dtype=None, label_shifts: dict[int, dict[int, np.ndarray]] | None = None):
    """Yield ``(t_index, frame)`` for every timepoint on disk: the label image of the timepoint placed on the
    registered, cropped, strided grid (place_frame; with ``label_shifts``, the per-nucleus shifts of the piecewise
    registration, place_frame_labels). One label image is in memory at a time."""
    for t in layout.timepoints:
        lab = layout.read_labels(t)
        if label_shifts is not None:
            yield int(t), place_frame_labels(lab, label_shifts.get(int(t), {}), box, stride,
                                             _keep_ids(layout, t, crop), dtype=dtype)
        else:
            yield int(t), place_frame(lab, (vox_shifts or {}).get(int(t), np.zeros(3, int)), box, stride,
                                      _keep_ids(layout, t, crop), dtype=dtype)


def place_frame_labels(lab: np.ndarray, shifts_by_label: dict[int, np.ndarray], box: dict, stride: int,
                       ids: np.ndarray | None = None, dtype=None) -> np.ndarray:
    """``lab`` placed on the registered grid with a whole-voxel shift per label (the piecewise registration gives
    every nucleus its own displacement): the voxels of label L land at their raw voxel minus shift(L), on the grid
    of ``box`` sub-sampled by ``stride`` in y and x, as place_frame does with one shift for the frame (the two
    agree when every label has the same shift). A label without an entry keeps a zero shift; where two labels
    land on one registered voxel the later one in raw C order wins."""
    s = max(1, int(stride))
    z0, y0, x0 = box["start"]
    z1, y1, x1 = box["stop"]
    out = np.zeros((z1 - z0, len(range(y0, y1, s)), len(range(x0, x1, s))), dtype=dtype or lab.dtype)
    keep = _mask_ids(lab, ids)
    zz, yy, xx = np.nonzero(keep)
    if zz.size == 0:
        return out
    L = keep[zz, yy, xx].astype(np.int64)
    lut = np.zeros((int(L.max()) + 1, 3), dtype=np.int64)
    for lbl, v in shifts_by_label.items():
        if 0 < int(lbl) < len(lut):
            lut[int(lbl)] = np.asarray(v, dtype=np.int64)
    rz, ry, rx = zz - lut[L, 0], yy - lut[L, 1], xx - lut[L, 2]
    m = ((rz >= z0) & (rz < z1) & (ry >= y0) & (ry < y1) & (rx >= x0) & (rx < x1)
         & ((ry - y0) % s == 0) & ((rx - x0) % s == 0))
    out[rz[m] - z0, (ry[m] - y0) // s, (rx[m] - x0) // s] = L[m]
    return out


def label_voxel_shifts(objects: pd.DataFrame, spacing) -> dict[int, dict[int, np.ndarray]]:
    """``{t_index: {nucleus_id: whole-voxel shift (z, y, x)}}`` from the object table's per-nucleus displacement
    (``disp_*_um``, rounded to voxels; a NaN displacement is no shift)."""
    sp = np.asarray(spacing, dtype=float)
    d = objects[["disp_z_um", "disp_y_um", "disp_x_um"]].to_numpy(float)
    v = np.rint(np.nan_to_num(d, nan=0.0) / sp).astype(int)
    out: dict[int, dict[int, np.ndarray]] = {}
    for t, nid, vv in zip(objects["t_index"].astype(int), objects["nucleus_id"].astype(int), v):
        out.setdefault(int(t), {})[int(nid)] = vv
    return out


def crop_box_labels(layout: SeriesLayout, label_shifts: dict[int, dict[int, np.ndarray]], crop: str = "germline") -> dict:
    """crop_box for per-label shifts: the union bounding box, on the registered grid, of every kept label moved by
    its own whole-voxel shift (scipy.ndimage.find_objects gives each label's extent in one pass per frame)."""
    from scipy import ndimage as ndi

    start = stop = shape = None
    n_labels: dict[int, int] = {}
    for t in layout.timepoints:
        lab = layout.read_labels(t)
        shape = np.array(lab.shape) if shape is None else shape
        ids = _keep_ids(layout, t, crop)
        keep = _mask_ids(lab, ids)
        n_labels[int(t)] = int(len(ids)) if ids is not None else int(len(layout.read_nuclei(t)))
        sh = label_shifts.get(int(t), {})
        for i, slc in enumerate(ndi.find_objects(keep)):
            if slc is None:
                continue
            v = np.asarray(sh.get(i + 1, np.zeros(3, int)), dtype=int)
            lo = np.array([slc[0].start, slc[1].start, slc[2].start]) - v
            hi = np.array([slc[0].stop, slc[1].stop, slc[2].stop]) - v
            start = lo if start is None else np.minimum(start, lo)
            stop = hi if stop is None else np.maximum(stop, hi)
        del lab, keep
    if start is None:
        raise SeriesLayoutError(f"{layout.stem}: no label to place (every timepoint is empty for crop={crop!r})")
    return {"start": [int(v) for v in start], "stop": [int(v) for v in stop],
            "shape_raw": [int(v) for v in shape], "n_labels": n_labels}


def _label_dtype(max_label: int):
    return np.uint16 if int(max_label) < np.iinfo(np.uint16).max else np.uint32


# =====================================================================================================
# track-export: the Cell-ACDC folder
# =====================================================================================================

def export_cellacdc(series_dir: str | Path, out_dir: str | Path, *, crop: str | None = None, stride: int | None = None,
                    image: str = "ref", channel_name: str = "central_element", stem: str | None = None,
                    cfg=None) -> dict:
    """Write the series' germline labels (and the central_element channel) as a Cell-ACDC position folder
    under ``out_dir`` (module docstring). ``image``: ``"ref"`` builds the image stack from the downsampled
    reference volumes the loop keeps (``t<NNN>/<id>__ref_stride<k>.npy``, nearest-neighbour resampled onto the
    export grid, for display), ``"file"`` reads each timepoint of the series file at the export stride (slow
    over the network: one full frame per timepoint), ``"none"`` writes no image. Returns the export record
    (also written as ``out_dir/export.json``)."""
    t0 = time.perf_counter()
    crop = str(crop or cfg_get(cfg, "timelapse.track.export_crop", DEFAULTS["export_crop"])).lower()
    if crop not in ("germline", "all"):
        raise ValueError(f"crop must be 'germline' or 'all', got {crop!r}")
    stride = int(stride or cfg_get(cfg, "timelapse.track.export_stride", DEFAULTS["export_stride"]))
    if stride < 1:
        raise ValueError("stride must be >= 1")
    if image not in ("ref", "file", "none"):
        raise ValueError(f"image must be 'ref', 'file' or 'none', got {image!r}")
    layout = find_series(series_dir, stem)
    spacing = layout.spacing()
    shifts, shift_source = load_shifts(layout)
    vsh = voxel_shifts(shifts, spacing)
    box = crop_box(layout, vsh, crop)
    z0, y0, x0 = box["start"]
    z1, y1, x1 = box["stop"]
    shape_exp = (len(layout.timepoints), z1 - z0, len(range(y0, y1, stride)), len(range(x0, x1, stride)))
    max_label = 0
    for t in layout.timepoints:
        df = layout.read_nuclei(t)
        if len(df):
            max_label = max(max_label, int(df["nucleus_id"].max()))
    dtype = _label_dtype(max_label)
    est_gb = math.prod(shape_exp) * np.dtype(dtype).itemsize / 1e9
    cap = float(cfg_get(cfg, "timelapse.track.export_max_gb", DEFAULTS["export_max_gb"]))
    if est_gb > cap:
        log.warning("%s: the exported label array is about %.1f GB in memory (cap %.1f GB); raise --stride or "
                    "lower timelapse.track.export_max_gb knowingly", layout.stem, est_gb, cap)

    images_dir = Path(out_dir).joinpath(*POSITION_FOLDER)
    os.makedirs(long_path(images_dir), exist_ok=True)
    basename = f"{layout.stem}_"
    segm_path = images_dir / f"{basename}segm.npz"
    image_path = images_dir / f"{basename}{channel_name}.tif"
    meta_path = images_dir / f"{basename}metadata.csv"

    segm = np.zeros(shape_exp, dtype=dtype)
    label_maps: dict[str, dict[str, int]] = {}
    frames: list[int] = []
    for i, (t, frame) in enumerate(registered_frames(layout, vsh, box, stride, crop, dtype=dtype)):
        segm[i] = frame
        frames.append(t)
        ids = np.unique(frame)
        label_maps[str(t)] = {str(int(v)): int(v) for v in ids if v != 0}     # exported value -> nucleus_id (identity, explicit)
    np.savez_compressed(long_path(segm_path), segm)

    image_info: dict = {"source": image, "channel": channel_name, "file": None}
    if image != "none":
        img = _export_image(layout, vsh, box, stride, image, channel_name, shape_exp)
        if img is not None:
            import tifffile

            sp_img = (float(spacing[0]), float(spacing[1]) * stride, float(spacing[2]) * stride)
            tifffile.imwrite(long_path(image_path), img, compression="zlib", ome=True,
                             metadata={"axes": "TZYX",
                                       "PhysicalSizeX": sp_img[2], "PhysicalSizeXUnit": "µm",
                                       "PhysicalSizeY": sp_img[1], "PhysicalSizeYUnit": "µm",
                                       "PhysicalSizeZ": sp_img[0], "PhysicalSizeZUnit": "µm"})
            image_info["file"] = image_path.name
            image_info["dtype"] = str(img.dtype)
        else:
            image_info["source"] = "none"
            image_info["reason"] = "no reference volume on disk" if image == "ref" else "the series file could not be read"

    sp_exp = (float(spacing[0]), float(spacing[1]) * stride, float(spacing[2]) * stride)
    meta_rows = [("SizeT", shape_exp[0]), ("SizeZ", shape_exp[1]), ("SizeY", shape_exp[2]), ("SizeX", shape_exp[3]),
                 ("TimeIncrement", float(layout.period_s) if layout.period_s else 1.0),
                 ("PhysicalSizeZ", sp_exp[0]), ("PhysicalSizeY", sp_exp[1]), ("PhysicalSizeX", sp_exp[2]),
                 ("segm_isSegm3D", True), ("basename", basename), ("channel_0_name", channel_name)]
    pd.DataFrame(meta_rows, columns=["Description", "values"]).to_csv(long_path(meta_path), index=False)

    record = {
        "series_id": layout.stem, "series_dir": str(layout.out_dir), "file_path": layout.file_path,
        "timepoints": frames, "frame_of_t_index": {str(t): i for i, t in enumerate(frames)},
        "image_ids": {str(t): layout.image_id(t) for t in frames},
        "t_seconds": {str(t): layout.t_seconds_of(t) for t in frames},
        "crop": crop, "box": {"start": box["start"], "stop": box["stop"], "shape_raw": box["shape_raw"]},
        "stride": stride, "shape": list(shape_exp), "dtype": str(np.dtype(dtype)),
        "spacing_raw_um": [float(v) for v in spacing], "spacing_export_um": list(sp_exp),
        "shift_source": shift_source,
        "shifts": {str(int(r.t_index)): {"dz_um": float(r.dz_um), "dy_um": float(r.dy_um), "dx_um": float(r.dx_um),
                                          "voxels": [int(v) for v in vsh[int(r.t_index)]]}
                   for r in shifts.itertuples(index=False)},
        "n_labels": {str(k): v for k, v in box["n_labels"].items()},
        "label_maps": label_maps,
        "image": image_info,
        "files": {"segm": str(segm_path.relative_to(Path(out_dir))), "image": image_info["file"] and str(image_path.relative_to(Path(out_dir))),
                  "metadata": str(meta_path.relative_to(Path(out_dir)))},
        "cellacdc_layout": {"position": POSITION_FOLDER[0], "images_folder": POSITION_FOLDER[1], "basename": basename},
        "coordinates": "exported voxel (z, y, x) = registered raw voxel (z0 + z, y0 + stride*y, x0 + stride*x); "
                       "registered = raw - shift of the timepoint; nothing is resampled",
        "elapsed_s": round(time.perf_counter() - t0, 3), "written": provenance.run_timestamp(),
    }
    with open(long_path(Path(out_dir) / EXPORT_RECORD), "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=1, default=str)
    log.info("%s: exported %d timepoint(s) as %s (%s, stride %d, crop %s) -> %s", layout.stem, len(frames),
             segm_path.name, "x".join(str(v) for v in shape_exp), stride, crop, out_dir)
    return record


def _export_image(layout: SeriesLayout, vsh, box, stride, image, channel_name, shape_exp) -> np.ndarray | None:
    """The central_element channel on the export grid as a (T, Z, Y, X) uint16 stack, from the reference
    volumes (``ref``) or the series file (``file``); None when neither is available."""
    out = np.zeros(shape_exp, dtype=np.uint16)
    z0, y0, x0 = box["start"]
    z1, y1, x1 = box["stop"]
    role_idx = layout.channel_index("central_element")
    for i, t in enumerate(layout.timepoints):
        sv = vsh.get(int(t), np.zeros(3, int))
        if image == "ref":
            ref = layout.ref_path(t)
            if ref is None:
                log.warning("%s: no reference volume for t%03d; no image written", layout.stem, t)
                return None
            vol = np.load(long_path(ref[0]), mmap_mode="r")
            r = ref[1]
        else:
            vol = _read_channel_from_file(layout, t, stride, role_idx)
            if vol is None:
                return None
            r = stride
        Z, Yr, Xr = vol.shape
        zi = np.arange(z0, z1) + int(sv[0])
        yi = (np.arange(y0, y1, stride) + int(sv[1])) // r         # nearest lower reference pixel
        xi = (np.arange(x0, x1, stride) + int(sv[2])) // r
        vz, vy, vx = (zi >= 0) & (zi < Z), (yi >= 0) & (yi < Yr), (xi >= 0) & (xi < Xr)
        if vz.any() and vy.any() and vx.any():
            sub = np.asarray(vol[np.ix_(zi[vz], yi[vy], xi[vx])])
            out[i][np.ix_(vz, vy, vx)] = np.clip(sub, 0, np.iinfo(np.uint16).max).astype(np.uint16)
    return out


def _read_channel_from_file(layout: SeriesLayout, t: int, stride: int, role_idx: int | None) -> np.ndarray | None:
    """One timepoint's central_element channel read from the series file at ``stride`` relative to the run's
    label grid: the file is read at the run's own xy stride times ``stride`` (and on the run's z range), so the
    volume sits on the grid label/stride exactly like a reference volume of stride ``stride``. A full frame read:
    25 s per timepoint over the lab network for the October series."""
    if not layout.file_path or not Path(layout.file_path).exists():
        log.warning("%s: the series file %s is not reachable; no image written", layout.stem, layout.file_path)
        return None
    from ..io import read_image

    try:
        st = read_image(layout.file_path, xy_stride=int(stride) * layout.run_xy_stride(),
                        z_range=layout.run_z_range(), t=int(t))
    except Exception as e:  # noqa: BLE001 - the labels are the product; the image is a convenience
        log.warning("%s: reading t%03d from %s failed (%s: %s); no image written", layout.stem, t, layout.file_path,
                    type(e).__name__, e)
        return None
    idx = role_idx if role_idx is not None and role_idx < st.n_channels else 0
    return np.asarray(st.data[idx])


# =====================================================================================================
# track-link: the backends
# =====================================================================================================

def cellacdc_overlap_available() -> tuple[bool, str]:
    """Whether Cell-ACDC's overlap tracker imports in this process without a Qt binding."""
    try:
        from cellacdc.trackers.CellACDC import CellACDC_tracker  # noqa: F401
    except Exception as e:  # noqa: BLE001 - ImportError or anything a GUI dependency raises
        return False, f"{type(e).__name__}: {e}"
    return True, "cellacdc.trackers.CellACDC.CellACDC_tracker"


def overlap_assign(prev_tracked: np.ndarray, cur: np.ndarray, *, min_ioa: float, next_id: int) -> tuple[dict[int, int], int]:
    """Cell-ACDC's overlap rule in NumPy: every current object takes the previous (tracked) object it covers
    most, measured as intersection over the previous object's volume (IoA), when that IoA is at least
    ``min_ioa`` and no other current object covers that previous object better; every other current object
    starts a new track at ``next_id`` upwards. Returns ``({current id: track id}, next free id)``."""
    cur_ids = np.unique(cur)
    cur_ids = cur_ids[cur_ids > 0]
    assign: dict[int, int] = {}
    both = (prev_tracked > 0) & (cur > 0)
    if both.any() and cur_ids.size:
        pc = prev_tracked[both].astype(np.int64)
        cc = cur[both].astype(np.int64)
        pmax = int(pc.max()) + 1
        pairs, inter = np.unique(cc * pmax + pc, return_counts=True)
        c_of, p_of = pairs // pmax, pairs % pmax
        prev_area = np.bincount(prev_tracked.ravel().astype(np.int64))
        ioa = inter / prev_area[p_of]
        best: dict[int, tuple[int, float]] = {}                 # current id -> (previous id, IoA)
        for c, p, v in zip(c_of, p_of, ioa):
            if v >= float(min_ioa) and (int(c) not in best or v > best[int(c)][1]):
                best[int(c)] = (int(p), float(v))
        claims: dict[int, tuple[int, float]] = {}               # previous id -> (winning current id, IoA)
        for c, (p, v) in best.items():
            if p not in claims or v > claims[p][1]:
                claims[p] = (c, v)
        for p, (c, _v) in claims.items():
            assign[c] = p
    nid = int(next_id)
    for c in cur_ids:
        c = int(c)
        if c not in assign:
            assign[c] = nid
            nid += 1
    return assign, nid


def majority_map(ours: np.ndarray, other: np.ndarray) -> dict[int, tuple[int, int, list[tuple[int, int]]]]:
    """For every label of ``ours`` (> 0), the value of ``other`` covering most of its voxels: ``{label:
    (majority value, its voxel count, [(value, count), ...] sorted by count)}``, from one pass over the
    voxels (no per-label scan)."""
    both = ours > 0
    o = ours[both].astype(np.int64)
    v = other[both].astype(np.int64)
    out: dict[int, tuple[int, int, list[tuple[int, int]]]] = {}
    if o.size == 0:
        return out
    vmax = int(v.max()) + 1
    pairs, counts = np.unique(o * vmax + v, return_counts=True)
    per_label: dict[int, list[tuple[int, int]]] = {}
    for a, b, c in zip(pairs // vmax, pairs % vmax, counts):
        per_label.setdefault(int(a), []).append((int(b), int(c)))
    for a, lst in per_label.items():
        lst.sort(key=lambda bc: -bc[1])
        out[a] = (lst[0][0], lst[0][1], lst)
    return out


def _relabel(frame: np.ndarray, assign: dict[int, int], dtype=np.int64) -> np.ndarray:
    lut = np.zeros(int(frame.max()) + 1 if frame.size else 1, dtype=dtype)
    for c, tid in assign.items():
        if c < lut.size:
            lut[c] = tid
    return lut[frame]


def _link_overlap(layout: SeriesLayout, vsh, box, stride, *, min_ioa: float, use_cellacdc: bool,
                  germline_only: bool, label_shifts=None) -> tuple[dict[tuple[int, int], int], dict]:
    """Frame-by-frame overlap linking on the registered label grid. With ``use_cellacdc`` each frame pair goes
    through ``CellACDC_tracker.track_frame`` (the function Cell-ACDC's own tracker loop calls); otherwise
    through overlap_assign. ``label_shifts`` (the piecewise registration) moves every label by its own
    whole-voxel shift instead of the frame's one. Returns ``{(t_index, nucleus_id): time_track_id}`` and the
    backend's details. The first timepoint keeps its nucleus ids as track ids."""
    if use_cellacdc:
        from cellacdc.trackers.CellACDC import CellACDC_tracker as CT
        from skimage.measure import regionprops
    crop = "germline" if germline_only else "all"
    assign: dict[tuple[int, int], int] = {}
    prev_tracked = None
    next_id = 1
    n_new, n_linked = 0, 0
    for t, frame in registered_frames(layout, vsh, box, stride, crop, dtype=np.int32, label_shifts=label_shifts):
        ids = np.unique(frame)
        ids = ids[ids > 0]
        if prev_tracked is None:
            for v in ids:
                assign[(t, int(v))] = int(v)
            prev_tracked = frame.astype(np.int64)
            next_id = int(ids.max()) + 1 if ids.size else 1
            continue
        if use_cellacdc:
            if ids.size and prev_tracked.any():
                tracked = CT.track_frame(prev_tracked, regionprops(prev_tracked), frame, regionprops(frame),
                                         IoA_thresh=float(min_ioa), uniqueID=next_id)
            else:
                # Cell-ACDC skips an empty frame and keeps a frame after an empty one unchanged; new ids here
                tracked = frame.astype(np.int64) + np.where(frame > 0, next_id - 1, 0)
            mapping = {a: m for a, (m, _cnt, _lst) in majority_map(frame, np.asarray(tracked).astype(np.int64)).items()}
        else:
            mapping, _ = overlap_assign(prev_tracked, frame, min_ioa=min_ioa, next_id=next_id)
        for v, tid in mapping.items():
            assign[(t, int(v))] = int(tid)
            if tid >= next_id:
                n_new += 1
            else:
                n_linked += 1
        if mapping:
            next_id = max(next_id, max(mapping.values()) + 1)
        prev_tracked = _relabel(frame, mapping) if mapping else np.zeros_like(frame, dtype=np.int64)
    info = {"min_overlap": float(min_ioa), "stride": int(stride), "box": {"start": box["start"], "stop": box["stop"]},
            "n_linked_steps": n_linked, "n_new_tracks_after_t0": n_new,
            "implementation": "cellacdc.trackers.CellACDC.CellACDC_tracker.track_frame" if use_cellacdc else "germquant overlap_assign"}
    return assign, info


def _link_laptrack(objects: pd.DataFrame, *, max_dist_um: float, gap_frames: int) -> tuple[dict[tuple[int, int], int], dict]:
    """LAP linking (laptrack) of the registered centroids: the squared Euclidean cost with cutoff
    ``max_dist_um ** 2`` for frame-to-frame links and for gap closing over at most ``gap_frames`` frames (a
    track ending at t joins one starting at t + d when 1 <= d <= gap_frames; 0 = no gap closing); splitting
    and merging off. A timepoint with no object is an empty frame, which gap closing bridges like a dropout."""
    try:
        from laptrack import LapTrack
    except ImportError as e:
        raise RuntimeError("the laptrack backend needs the timelapse extra: uv pip install -e \".[timelapse]\" "
                           "(laptrack>=0.17); the overlap backend needs nothing extra") from e

    df = objects[["t_index", "nucleus_id", "zr_um", "yr_um", "xr_um"]].reset_index(drop=True)
    if df.empty:
        return {}, {"n_points": 0}
    tmin = int(df["t_index"].min())
    df["frame__"] = (df["t_index"].astype(int) - tmin)
    cutoff = float(max_dist_um) ** 2
    lt = LapTrack(metric="sqeuclidean", cutoff=cutoff,
                  gap_closing_metric="sqeuclidean",
                  gap_closing_cutoff=cutoff if int(gap_frames) > 0 else False,
                  gap_closing_max_frame_count=max(1, int(gap_frames)),
                  splitting_cutoff=False, merging_cutoff=False)
    track_df, split_df, merge_df = lt.predict_dataframe(df, ["zr_um", "yr_um", "xr_um"], frame_col="frame__",
                                                        only_coordinate_cols=False)
    track_df = track_df.reset_index()
    assign = {(int(r.t_index), int(r.nucleus_id)): int(r.track_id) + 1 for r in track_df.itertuples(index=False)}
    info = {"max_dist_um": float(max_dist_um), "gap_frames": int(gap_frames), "cutoff_sq_um2": cutoff,
            "n_points": int(len(track_df)), "n_splits": int(len(split_df)), "n_merges": int(len(merge_df)),
            "implementation": "laptrack.LapTrack.predict_dataframe"}
    return assign, info


def link_series(series_dir: str | Path, *, backend: str = "overlap", gap: int | None = None,
                max_dist_um: float | None = None, min_overlap: float | None = None, stride: int | None = None,
                germline_only: bool = True, stem: str | None = None, cfg=None, write: bool = True,
                registration: str | None = None) -> dict:
    """``germquant track-link``: link the nuclei of a finished series across its timepoints with ``backend``
    (BACKENDS) and write ``<stem>__track_points.csv`` and ``<stem>__tracks.csv`` next to the timepoint folders,
    plus ``<stem>__track_link.json`` (the backend, its parameters, the registration and shift source, the axis
    sign per timepoint). ``overlap`` falls back to ``overlap_inhouse`` (and says so) when Cell-ACDC cannot be
    imported. ``registration`` (REGISTRATIONS; default the ``timelapse.track.registration`` key, code default
    ``rigid``): ``rigid`` subtracts the timepoint's cumulative shift of ``<stem>__register_shifts.csv`` (zero
    without the table), ``piecewise`` each nucleus's displacement interpolated along the gonad axis from
    ``<stem>__register_piecewise.csv`` (``germquant register --mode piecewise``; required). Returns the points and
    tracks DataFrames, the link record and the file paths."""
    t_start = time.perf_counter()
    backend = str(backend).lower()
    if backend not in BACKENDS:
        raise ValueError(f"backend must be one of {', '.join(BACKENDS)}, got {backend!r}")
    reg_mode = str(registration if registration is not None
                   else cfg_get(cfg, "timelapse.track.registration", DEFAULTS["registration"])).strip().lower()
    if reg_mode not in REGISTRATIONS:
        raise ValueError(f"registration must be one of {', '.join(REGISTRATIONS)}, got {reg_mode!r}")
    layout = find_series(series_dir, stem)
    spacing = layout.spacing()
    shifts, shift_source = load_shifts(layout)
    piecewise = None
    if reg_mode == "piecewise":
        from .piecewise import load_piecewise

        piecewise = load_piecewise(layout)
        shift_source = "register_piecewise"
    objects = object_table(layout, shifts, germline_only=germline_only, registration=piecewise)
    gap = int(gap if gap is not None else cfg_get(cfg, "timelapse.track.gap_frames", DEFAULTS["gap_frames"]))
    max_dist_um = float(max_dist_um if max_dist_um is not None else cfg_get(cfg, "timelapse.track.max_dist_um", DEFAULTS["max_dist_um"]))
    min_overlap = float(min_overlap if min_overlap is not None else cfg_get(cfg, "timelapse.track.min_overlap", DEFAULTS["min_overlap"]))
    stride = int(stride if stride is not None else cfg_get(cfg, "timelapse.track.link_stride", DEFAULTS["link_stride"]))

    used = backend
    note = ""
    if backend == "overlap":
        ok, why = cellacdc_overlap_available()
        if not ok:
            used = "overlap_inhouse"
            note = f"Cell-ACDC's tracker could not be imported ({why}); the in-house overlap linking was used"
            log.warning("%s", note)
    if used in ("overlap", "overlap_inhouse"):
        crop = "germline" if germline_only else "all"
        if piecewise is None:
            vsh = voxel_shifts(shifts, spacing)
            box = crop_box(layout, vsh, crop)
            assign, binfo = _link_overlap(layout, vsh, box, stride, min_ioa=min_overlap,
                                          use_cellacdc=(used == "overlap"), germline_only=germline_only)
        else:
            lsh = label_voxel_shifts(objects, spacing)
            box = crop_box_labels(layout, lsh, crop)
            assign, binfo = _link_overlap(layout, None, box, stride, min_ioa=min_overlap,
                                          use_cellacdc=(used == "overlap"), germline_only=germline_only,
                                          label_shifts=lsh)
            binfo["placement"] = "per-nucleus whole-voxel shift (piecewise registration)"
    elif used == "laptrack":
        assign, binfo = _link_laptrack(objects, max_dist_um=max_dist_um, gap_frames=gap)
    else:
        assign, binfo = {}, {}
    link_source = used if used != "none" else "none"
    points, tracks, tinfo = build_track_tables(layout, objects, assign, link_source=link_source)
    record = {
        "series_id": layout.stem, "backend_requested": backend, "backend": used, "link_source": link_source, "note": note,
        "germline_only": bool(germline_only), "registration": reg_mode, "shift_source": shift_source,
        "registration_info": piecewise.summary() if piecewise is not None else None, "timepoints": layout.timepoints,
        "missing_timepoints": layout.missing, "n_objects": int(len(objects)),
        "n_tracks": int(tracks["time_track_id"].notna().sum()) if len(tracks) else 0,
        "n_full_length": int(tracks["is_full_length"].fillna(False).astype(bool).sum()) if len(tracks) else 0,
        "backend_info": binfo, **tinfo,
        "elapsed_s": round(time.perf_counter() - t_start, 3), "written": provenance.run_timestamp(),
    }
    paths = {}
    if write:
        paths = write_tracking_tables(layout, points, tracks, record)
    log.info("%s: %s linked %d nucleus-timepoints into %d track(s), %d full length -> %s", layout.stem, used,
             record["n_objects"], record["n_tracks"], record["n_full_length"], layout.out_dir)
    return {"points": points, "tracks": tracks, "record": record, "paths": paths, "layout": layout}


# =====================================================================================================
# the two tables
# =====================================================================================================

def build_track_tables(layout: SeriesLayout, objects: pd.DataFrame, assign: dict[tuple[int, int], int], *,
                       link_source: str) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """The track_points and tracks DataFrames from the object table and ``{(t_index, nucleus_id):
    time_track_id}``. A nucleus without an assignment keeps its row with time_track_id NA and link_source
    ``none`` (so every measurement is still there per timepoint). Steps are between consecutive points of a
    track on the registered coordinates; ``is_gap`` marks a point that follows a gap of one or more missing
    timepoints in its track. Returns the two frames and the axis-sign record (align_axis_sign)."""
    pts = objects.copy().reset_index(drop=True)
    keys = list(zip(pts["t_index"].astype(int), pts["nucleus_id"].astype(int)))
    tid = [assign.get(k) for k in keys]
    pts["time_track_id"] = pd.array([None if v is None else int(v) for v in tid], dtype="Int64")
    pts["link_source"] = np.where(pts["time_track_id"].notna(), link_source, "none")
    pts = pts.sort_values(["time_track_id", "t_index", "nucleus_id"], kind="stable", na_position="last").reset_index(drop=True)

    for c in ("step_um", "step_z_um", "step_y_um", "step_x_um"):
        pts[c] = np.nan
    pts["is_gap"] = False
    linked = pts["time_track_id"].notna()
    if linked.any():
        g = pts[linked].groupby("time_track_id", sort=False)
        prev_t = g["t_index"].shift(1)
        for axis, col in (("z", "zr_um"), ("y", "yr_um"), ("x", "xr_um")):
            pts.loc[linked, f"step_{axis}_um"] = pts.loc[linked, col] - g[col].shift(1)
        pts.loc[linked, "step_um"] = np.sqrt(pts.loc[linked, "step_z_um"] ** 2 + pts.loc[linked, "step_y_um"] ** 2
                                             + pts.loc[linked, "step_x_um"] ** 2)
        pts.loc[linked, "is_gap"] = ((pts.loc[linked, "t_index"] - prev_t) > 1).fillna(False).to_numpy()
    pts["reg_z_um"], pts["reg_y_um"], pts["reg_x_um"] = pts["zr_um"], pts["yr_um"], pts["xr_um"]

    # the per-timepoint quantities the track metrics need: the axis sign and the SYP intensity median
    signs, lengths = align_axis_sign(layout, pts)
    if "axis_position_um" in pts.columns:
        p = pd.to_numeric(pts["axis_position_um"], errors="coerce")
        s = pts["t_index"].map(signs).fillna(1.0)
        L = pts["t_index"].map(lengths)
        pts["axis_position_aligned_um"] = np.where(s < 0, L - p, p)
    # one reference axis for every timepoint (the per-frame axis is refitted and its origin moves; see
    # fixed_axis_positions), used for the axis velocity of the tracks
    pts["axis_position_fixed_um"] = fixed_axis_positions(pts)
    # the same axis at the RAW centroids (the microscope frame, whatever the registration): its slope is
    # axis_velocity_raw_um_per_min, the flow along the gonad when the registered frame follows the flow itself
    raw_view = pd.DataFrame({"t_index": pts["t_index"], "zr_um": pts["z_um"], "yr_um": pts["y_um"],
                             "xr_um": pts["x_um"]}, index=pts.index)
    for c in ("in_germline", "axis_position_um"):
        if c in pts.columns:
            raw_view[c] = pts[c]
    pts["axis_position_fixed_raw_um"] = fixed_axis_positions(raw_view)
    syp = "central_element_mean_intensity"
    if syp in pts.columns:
        med = pts.groupby("t_index")[syp].transform("median")
        pts[f"{syp}_norm"] = pts[syp] / med.where(med > 0)
    pts["minutes"] = pts["t_index"].map(layout.minutes_of)

    tracks = _tracks_from_points(pts, layout, link_source)
    pts = pts.drop(columns=["minutes"])
    info = {"axis_sign_by_t": {str(k): int(v) for k, v in signs.items()},
            "axis_length_by_t": {str(k): (None if (v is None or not np.isfinite(v)) else float(v)) for k, v in lengths.items()}}
    return pts, tracks, info


def align_axis_sign(layout: SeriesLayout, pts: pd.DataFrame) -> tuple[dict[int, int], dict[int, float]]:
    """The axis orientation of every timepoint relative to the first linked timepoint (t0): +1 when the
    distal-proximal order of the nuclei linked between t0 and the timepoint agrees (Spearman rank correlation
    of their axis_position_um >= 0), -1 when it is reversed (the axis stage's orientation heuristic flipped).
    A timepoint with fewer than three nuclei linked to t0 inherits the sign of the previous timepoint. Also
    returns the axis length per timepoint used to flip a reversed frame (``p' = L - p``): total_germline_length_um
    from the timepoint's image_summary, else the largest axis_position_um of the frame."""
    signs: dict[int, int] = {}
    lengths: dict[int, float] = {}
    ts = list(layout.timepoints)
    if not ts:
        return signs, lengths
    has_axis = "axis_position_um" in pts.columns
    for t in ts:
        L = np.nan
        s = layout.read_summary(t)
        if s is not None and "total_germline_length_um" in s.index:
            try:
                L = float(s["total_germline_length_um"])
            except (TypeError, ValueError):
                L = np.nan
        if has_axis and not np.isfinite(L):
            v = pd.to_numeric(pts.loc[pts["t_index"] == t, "axis_position_um"], errors="coerce")
            L = float(v.max()) if v.notna().any() else np.nan
        lengths[int(t)] = L
    t0 = ts[0]
    signs[int(t0)] = 1
    if not has_axis:
        for t in ts[1:]:
            signs[int(t)] = 1
        return signs, lengths
    from scipy.stats import spearmanr

    linked = pts[pts["time_track_id"].notna()]
    ref = linked[linked["t_index"] == t0].set_index("time_track_id")["axis_position_um"]
    prev = 1
    for t in ts[1:]:
        cur = linked[linked["t_index"] == t].set_index("time_track_id")["axis_position_um"]
        common = ref.index.intersection(cur.index)
        a = pd.to_numeric(ref.loc[common], errors="coerce")
        b = pd.to_numeric(cur.loc[common], errors="coerce")
        ok = a.notna() & b.notna()
        sign = prev
        if int(ok.sum()) >= 3:
            rho = spearmanr(a[ok], b[ok]).statistic
            if np.isfinite(rho) and rho != 0:
                sign = 1 if rho > 0 else -1
        signs[int(t)] = sign
        prev = sign
    return signs, lengths


def fixed_axis_positions(pts: pd.DataFrame, k_neighbours: int = 16) -> pd.Series:
    """The position of every track point along ONE reference axis: the first timepoint's axis, in microns.

    The axis stage refits the gonad axis at every timepoint and puts its origin at that frame's most distal
    labelled nucleus, so ``axis_position_um`` of a fixed point in the worm moves from frame to frame (on the
    October series by 54 to 83 um at x = 100 um across 25 frames at stride 4) and its slope over time mixes
    flow with that jitter. Here the first timepoint's axis is applied to every frame: the axis position at a
    registered point is a local linear fit of the first timepoint's axis_position_um over its ``k_neighbours``
    nearest germline nuclei (registered coordinates), evaluated at the point. Within the gonad this follows
    the curved axis; a nucleus that moves past the last first-timepoint nucleus is extrapolated with the local
    gradient instead of being clamped. When the first timepoint carries no axis positions, a centerline is
    fitted to its nuclei with the axis stage's curve fit instead. NaN for every point when the first timepoint
    has fewer than 5 usable nuclei."""
    out = pd.Series(np.nan, index=pts.index, dtype=float)
    model = fixed_axis_model(pts, k_neighbours)
    if model is None:
        return out
    coords = pts[["zr_um", "yr_um", "xr_um"]].to_numpy(float)
    ok = np.isfinite(coords).all(axis=1)
    out.loc[pts.index[ok]] = model(coords[ok])
    return out


def fixed_axis_model(pts: pd.DataFrame, k_neighbours: int = 16):
    """The reference axis of fixed_axis_positions as a function: a callable that maps an (N, 3) array of
    registered coordinates (z, y, x um) to positions along the first timepoint's axis, built from the first
    timepoint's rows of ``pts`` (``zr_um``, ``yr_um``, ``xr_um``, ``t_index``, optional ``in_germline`` and
    ``axis_position_um``) by the rule fixed_axis_positions documents. None when the first timepoint has fewer
    than 5 usable nuclei. Each query point's value depends only on that point, so evaluating a subset gives the
    values the whole set would give (the piecewise registration evaluates it repeatedly, timelapse.piecewise)."""
    if pts.empty or not {"zr_um", "yr_um", "xr_um", "t_index"} <= set(pts.columns):
        return None
    coords = pts[["zr_um", "yr_um", "xr_um"]].to_numpy(float)
    ok = np.isfinite(coords).all(axis=1)
    if not ok.any():
        return None
    t_idx = pd.to_numeric(pts["t_index"], errors="coerce").to_numpy(float)
    t0 = float(np.nanmin(t_idx[ok]))
    sel = ok & (t_idx == t0)
    if "in_germline" in pts.columns:
        germ = pts["in_germline"].astype("boolean").fillna(False).to_numpy(bool)
        if int((sel & germ).sum()) >= 5:
            sel = sel & germ
    p0 = (pd.to_numeric(pts["axis_position_um"], errors="coerce").to_numpy(float)
          if "axis_position_um" in pts.columns else np.full(len(pts), np.nan))
    have_p = sel & np.isfinite(p0)
    if int(have_p.sum()) >= 5:
        ref, val = coords[have_p].copy(), p0[have_p].copy()

        def local_linear(x):
            return _local_linear_field(ref, val, np.asarray(x, dtype=float).reshape(-1, 3), k_neighbours)

        return local_linear
    if int(sel.sum()) < 5:
        return None
    return _centerline_model(coords[sel].copy())


def _local_linear_field(ref: np.ndarray, val: np.ndarray, x: np.ndarray, k: int) -> np.ndarray:
    """``val`` (a scalar field sampled at the points ``ref``) evaluated at ``x`` by a local linear least-squares
    fit over the ``k`` nearest samples of each query point. The intercept of the fit centred on the query point
    is the estimate, so a field that is linear in the neighbourhood is reproduced exactly; a tiny ridge on the
    gradient keeps degenerate neighbourhoods (all samples in one plane) solvable without biasing the estimate."""
    from scipy.spatial import cKDTree

    k = int(min(max(k, 4), len(ref)))
    tree = cKDTree(ref)
    res = np.empty(len(x))
    span = float(np.ptp(ref, axis=0).max()) or 1.0
    for i in range(0, len(x), 8192):
        xb = x[i:i + 8192]
        _d, idx = tree.query(xb, k=k)
        idx = np.atleast_2d(idx)
        d = (ref[idx] - xb[:, None, :]) / span                     # (n, k, 3), centred on the query, scaled
        A = np.concatenate([np.ones(d.shape[:2] + (1,)), d], axis=2)  # (n, k, 4)
        AtA = np.einsum("nki,nkj->nij", A, A)
        AtA[:, 1:, 1:] += 1e-9 * np.eye(3)[None]
        Atb = np.einsum("nki,nk->ni", A, val[idx])
        res[i:i + 8192] = np.linalg.solve(AtA, Atb[..., None])[:, 0, 0]
    return res


def _centerline_positions(ref: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Arc length along a centerline fitted to ``ref`` (the axis stage's curve fit), for the points ``x``, with
    a straight extension beyond both ends; a principal straight line when no curve can be fitted. Zeroed at
    the most distal reference point; the orientation is arbitrary (no axis positions to orient it by)."""
    return _centerline_model(ref)(x)


def _centerline_model(ref: np.ndarray):
    """The centerline of _centerline_positions fitted once to ``ref``, as a callable on query points."""
    from ..axis.linearize import _fit_principal_curve

    center = ref.mean(axis=0)
    scale = float(np.std(ref - center)) or 1.0
    ref_s = (ref - center) / scale
    curve, _method = _fit_principal_curve(ref_s)
    if curve is None or len(curve) < 2:
        _u, _sv, vt = np.linalg.svd(ref_s, full_matrices=False)

        def project(q: np.ndarray) -> np.ndarray:
            return q @ vt[0]
    else:
        seg = np.linalg.norm(np.diff(curve, axis=0), axis=1)
        s_curve = np.concatenate([[0.0], np.cumsum(seg)])
        t_start = (curve[1] - curve[0]) / (np.linalg.norm(curve[1] - curve[0]) or 1.0)
        t_end = (curve[-1] - curve[-2]) / (np.linalg.norm(curve[-1] - curve[-2]) or 1.0)

        def project(q: np.ndarray) -> np.ndarray:
            res = np.empty(len(q))
            for i in range(0, len(q), 4096):
                qb = q[i:i + 4096]
                kk = ((qb[:, None, :] - curve[None, :, :]) ** 2).sum(axis=2).argmin(axis=1)
                a = s_curve[kk].copy()
                first, last = kk == 0, kk == len(curve) - 1
                a[first] = (qb[first] - curve[0]) @ t_start
                a[last] = s_curve[-1] + (qb[last] - curve[-1]) @ t_end
                res[i:i + 4096] = a
            return res

    arc_min = float(np.min(project(ref_s)))

    def positions(x) -> np.ndarray:
        x_s = (np.asarray(x, dtype=float).reshape(-1, 3) - center) / scale
        return (project(x_s) - arc_min) * scale

    return positions


def _axis_velocity(minutes: np.ndarray, g: pd.DataFrame) -> float:
    """Velocity along the gonad in um per minute: the slope on the fixed reference axis (fixed_axis_positions)
    when the track has it, else on the per-frame axis after sign alignment."""
    for col in ("axis_position_fixed_um", "axis_position_aligned_um"):
        if col in g.columns:
            v = g[col].to_numpy(float)
            if int(np.isfinite(v).sum()) >= 2:
                return _slope_per_min(minutes, v)
    return float("nan")


def _slope_per_min(minutes: np.ndarray, values: np.ndarray) -> float:
    m = np.asarray(minutes, float)
    v = np.asarray(values, float)
    ok = np.isfinite(m) & np.isfinite(v)
    if ok.sum() < 2 or np.ptp(m[ok]) == 0:
        return float("nan")
    return float(np.polyfit(m[ok], v[ok], 1)[0])


def _tracks_from_points(pts: pd.DataFrame, layout: SeriesLayout, link_source: str) -> pd.DataFrame:
    n_t = len(layout.timepoints)
    t_lo, t_hi = (layout.timepoints[0], layout.timepoints[-1]) if n_t else (0, 0)
    rows = []
    linked = pts[pts["time_track_id"].notna()]
    for tid, g in linked.groupby("time_track_id", sort=True):
        g = g.sort_values("t_index")
        t = g["t_index"].astype(int).to_numpy()
        minutes = g["minutes"].to_numpy(float)
        n = len(g)
        steps = g["step_um"].to_numpy(float)
        path = float(np.nansum(steps)) if n > 1 else 0.0
        first, last = g.iloc[0], g.iloc[-1]
        net = float(math.sqrt((last["zr_um"] - first["zr_um"]) ** 2 + (last["yr_um"] - first["yr_um"]) ** 2
                              + (last["xr_um"] - first["xr_um"]) ** 2)) if n > 1 else 0.0
        dur = float(minutes[-1] - minutes[0]) if n > 1 else 0.0
        row = {
            "time_track_id": int(tid), "n_points": int(n), "n_gaps": int(t[-1] - t[0] + 1 - n),
            "t_start": int(t[0]), "t_end": int(t[-1]), "t_first": int(t[0]), "t_last": int(t[-1]),
            "is_full_length": bool(n == n_t and t[0] == t_lo and t[-1] == t_hi),
            "path_length_um": path, "net_displacement_um": net,
            "mean_step_um": float(np.nanmean(steps[1:])) if n > 1 else float("nan"),
            "speed_um_per_min": (path / dur) if dur > 0 else float("nan"),
            "axis_velocity_um_per_min": _axis_velocity(minutes, g),
            "axis_velocity_refit_um_per_min": _slope_per_min(minutes, g["axis_position_aligned_um"].to_numpy(float))
            if "axis_position_aligned_um" in g.columns else float("nan"),
            "intensity_slope_raw_per_min": _slope_per_min(minutes, g["central_element_mean_intensity"].to_numpy(float))
            if "central_element_mean_intensity" in g.columns else float("nan"),
            "intensity_slope_norm_per_min": _slope_per_min(minutes, g["central_element_mean_intensity_norm"].to_numpy(float))
            if "central_element_mean_intensity_norm" in g.columns else float("nan"),
            "sc_length_change_um": float("nan"), "link_source": link_source,
            "axis_velocity_raw_um_per_min": _slope_per_min(minutes, g["axis_position_fixed_raw_um"].to_numpy(float))
            if "axis_position_fixed_raw_um" in g.columns else float("nan"),
        }
        if "sc_total_length_um" in g.columns:
            sc = pd.to_numeric(g["sc_total_length_um"], errors="coerce").dropna()
            if len(sc) >= 2:
                row["sc_length_change_um"] = float(sc.iloc[-1] - sc.iloc[0])
        row["intensity_slope_raw"] = row["intensity_slope_raw_per_min"]
        row["intensity_slope_norm"] = row["intensity_slope_norm_per_min"]
        rows.append(row)
    cols = list(schema.TRACKS)
    df = pd.DataFrame(rows, columns=cols) if rows else pd.DataFrame(columns=cols)
    return df


def write_tracking_tables(layout: SeriesLayout, points: pd.DataFrame, tracks: pd.DataFrame, record: dict) -> dict:
    """Write the two series tables (schema order: declared columns, then the joined measurements, then
    SERIES_META) and the link record ``<stem>__track_link.json``. Returns the paths."""
    meta = layout.series_meta()
    p_points = write_series_table(points, "track_points", layout.out_dir, layout.stem, meta)
    p_tracks = write_series_table(tracks, "tracks", layout.out_dir, layout.stem, meta)
    p_rec = layout.out_dir / f"{layout.stem}__{LINK_RECORD}.json"
    with open(long_path(p_rec), "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=1, default=str)
    return {"track_points": str(p_points), "tracks": str(p_tracks), "record": str(p_rec)}


def read_track_points(layout: SeriesLayout) -> pd.DataFrame:
    p = layout.series_table_path("track_points")
    if not p.is_file():
        raise SeriesLayoutError(f"{p.name} not found: run `germquant track-link` or `track-import` first")
    df = pd.read_csv(long_path(p), low_memory=False)
    df["time_track_id"] = pd.array(pd.to_numeric(df["time_track_id"], errors="coerce").round().astype("Int64"))
    for c in ("t_index", "nucleus_id"):
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    return df


# =====================================================================================================
# track-import
# =====================================================================================================

def _col(df: pd.DataFrame, *names: str, required: bool = True) -> str | None:
    low = {str(c).strip().lower(): c for c in df.columns}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    if required:
        raise ValueError(f"no column named any of {', '.join(names)} in {list(df.columns)[:12]}")
    return None


def read_trackmate_csv(path: str | Path) -> pd.DataFrame:
    """TrackMate's spots export (ID, TRACK_ID, POSITION_X/Y/Z, FRAME; the extra header rows TrackMate writes
    under the column names are dropped) as ``track_id, t_index, z_um, y_um, x_um``. Spots without a track
    (TRACK_ID empty) are dropped and counted in the log."""
    df = pd.read_csv(long_path(Path(path)), low_memory=False)
    c_id, c_tr, c_fr = _col(df, "ID"), _col(df, "TRACK_ID"), _col(df, "FRAME")
    cx, cy, cz = _col(df, "POSITION_X"), _col(df, "POSITION_Y"), _col(df, "POSITION_Z", required=False)
    out = pd.DataFrame({
        "spot_id": pd.to_numeric(df[c_id], errors="coerce"), "track_id": pd.to_numeric(df[c_tr], errors="coerce"),
        "t_index": pd.to_numeric(df[c_fr], errors="coerce"),
        "x_um": pd.to_numeric(df[cx], errors="coerce"), "y_um": pd.to_numeric(df[cy], errors="coerce"),
        "z_um": pd.to_numeric(df[cz], errors="coerce") if cz else 0.0,
    })
    out = out[out["spot_id"].notna() & out["t_index"].notna()]          # the unit / description header rows go
    n_untracked = int(out["track_id"].isna().sum())
    if n_untracked:
        log.info("%s: %d spot(s) without a TRACK_ID dropped", Path(path).name, n_untracked)
    out = out[out["track_id"].notna()].copy()
    out["track_id"] = out["track_id"].astype(int)
    out["t_index"] = out["t_index"].astype(int)
    return out[["track_id", "t_index", "z_um", "y_um", "x_um"]].reset_index(drop=True)


def read_trackmate_xml(path: str | Path) -> pd.DataFrame:
    """A TrackMate model XML: spots from AllSpots, tracks from AllTracks' edges (FilteredTracks honoured when
    present) as ``track_id, t_index, z_um, y_um, x_um``."""
    import xml.etree.ElementTree as ET

    root = ET.parse(long_path(Path(path))).getroot()
    model = root.find("Model") if root.tag != "Model" else root
    if model is None:
        raise ValueError(f"{Path(path).name}: no <Model> element (is this a TrackMate XML?)")
    spots: dict[int, tuple[int, float, float, float]] = {}
    for sif in model.iter("SpotsInFrame"):
        frame = int(float(sif.get("frame", "0")))
        for sp in sif.iter("Spot"):
            sid = int(float(sp.get("ID")))
            fr = int(float(sp.get("FRAME", frame)))
            spots[sid] = (fr, float(sp.get("POSITION_Z", "0")), float(sp.get("POSITION_Y")), float(sp.get("POSITION_X")))
    filtered = {int(float(e.get("TRACK_ID"))) for e in model.iter("TrackID") if e.get("TRACK_ID") is not None}
    rows = []
    for tr in model.iter("Track"):
        tid = int(float(tr.get("TRACK_ID")))
        if filtered and tid not in filtered:
            continue
        members = set()
        for e in tr.iter("Edge"):
            members.add(int(float(e.get("SPOT_SOURCE_ID"))))
            members.add(int(float(e.get("SPOT_TARGET_ID"))))
        for sid in sorted(members):
            if sid in spots:
                fr, z, y, x = spots[sid]
                rows.append({"track_id": tid, "t_index": fr, "z_um": z, "y_um": y, "x_um": x})
    return pd.DataFrame(rows, columns=["track_id", "t_index", "z_um", "y_um", "x_um"])


def read_points_table(path: str | Path, *, time_base: int = 0) -> pd.DataFrame:
    """An Imaris (or hand) points table with ``id, t, z, y, x`` (aliases: track_id / hand_track_id / TrackID;
    t_index / frame / time; z_um ...) in raw-frame microns. ``time_base`` is the index of the first timepoint
    in the file (Imaris counts time from 1)."""
    df = pd.read_csv(long_path(Path(path)), low_memory=False)
    c_id = _col(df, "id", "track_id", "hand_track_id", "trackid", "track")
    c_t = _col(df, "t", "t_index", "frame", "time", "timepoint")
    c_z, c_y, c_x = _col(df, "z", "z_um", "position z"), _col(df, "y", "y_um", "position y"), _col(df, "x", "x_um", "position x")
    out = pd.DataFrame({"track_id": pd.to_numeric(df[c_id], errors="coerce"),
                        "t_index": pd.to_numeric(df[c_t], errors="coerce") - int(time_base),
                        "z_um": pd.to_numeric(df[c_z], errors="coerce"), "y_um": pd.to_numeric(df[c_y], errors="coerce"),
                        "x_um": pd.to_numeric(df[c_x], errors="coerce")})
    out = out.dropna(subset=["track_id", "t_index", "z_um", "y_um", "x_um"]).copy()
    out["track_id"] = out["track_id"].astype(int)
    out["t_index"] = out["t_index"].round().astype(int)
    return out.reset_index(drop=True)


def lookup_labels(layout: SeriesLayout, points: pd.DataFrame, spacing, *, radius_um: float = 0.0) -> pd.Series:
    """The nucleus label each raw-frame point (``t_index``, ``z_um``, ``y_um``, ``x_um``) lands in, as an Int64
    Series aligned with ``points`` (NA = outside every label, or a timepoint not on disk). With ``radius_um``
    > 0 a point in the background takes the nearest label centroid within that radius instead."""
    sp = np.asarray(spacing, float)
    out = pd.array([None] * len(points), dtype="Int64")
    for t, g in points.groupby("t_index"):
        if int(t) not in layout.timepoints:
            continue
        lab = layout.read_labels(int(t))
        vox = np.rint(g[["z_um", "y_um", "x_um"]].to_numpy(float) / sp).astype(int)
        inside = np.all((vox >= 0) & (vox < np.array(lab.shape)), axis=1)
        vals = np.zeros(len(g), dtype=np.int64)
        if inside.any():
            vi = vox[inside]
            vals[inside] = lab[vi[:, 0], vi[:, 1], vi[:, 2]]
        if radius_um > 0 and (vals == 0).any():
            nuc = layout.read_nuclei(int(t))
            if len(nuc):
                from scipy.spatial import cKDTree

                cents = nuc[["centroid_z_um", "centroid_y_um", "centroid_x_um"]].to_numpy(float)
                tree = cKDTree(cents)
                miss = np.flatnonzero(vals == 0)
                d, idx = tree.query(g[["z_um", "y_um", "x_um"]].to_numpy(float)[miss])
                ok = d <= float(radius_um)
                vals[miss[ok]] = nuc["nucleus_id"].to_numpy(int)[idx[ok]]
        pos = points.index.get_indexer(g.index)
        for p, v in zip(pos, vals):
            out[p] = int(v) if v > 0 else None
    return pd.Series(out, index=points.index)


def _export_to_raw(points: pd.DataFrame, export: dict) -> pd.DataFrame:
    """Points measured on the exported stack (microns of the export grid) back to raw-frame microns: registered
    = export box origin + coordinate; raw = registered + the timepoint's shift (the frame index of the stack
    is mapped back to t_index with the export record)."""
    sp = np.asarray(export["spacing_raw_um"], float)
    origin = np.asarray(export["box"]["start"], float) * sp
    frame_to_t = {int(v): int(k) for k, v in export["frame_of_t_index"].items()}
    out = points.copy()
    out["t_index"] = out["t_index"].map(lambda f: frame_to_t.get(int(f), int(f)))
    reg = out[["z_um", "y_um", "x_um"]].to_numpy(float) + origin
    sh = np.array([[export["shifts"].get(str(int(t)), {}).get(k, 0.0) for k in ("dz_um", "dy_um", "dx_um")]
                   for t in out["t_index"]], dtype=float)
    out[["z_um", "y_um", "x_um"]] = reg + sh
    return out


def _find_export_record(file: Path, given: str | Path | None) -> dict | None:
    cands = [Path(given)] if given else []
    cands += [p / EXPORT_RECORD for p in (file.parent, *file.parents[1:4])]
    for c in cands:
        if c.is_file():
            with open(long_path(c), encoding="utf-8") as fh:
                return json.load(fh)
    return None


def _assign_from_points(layout: SeriesLayout, points: pd.DataFrame, spacing, *, snap_um: float = 0.0) -> tuple[dict, dict]:
    """``{(t_index, nucleus_id): track_id}`` from imported points, with the counts of what did not map."""
    labels = lookup_labels(layout, points, spacing, radius_um=snap_um)
    pts = points.assign(nucleus_id=labels)
    n_outside = int(pts["nucleus_id"].isna().sum())
    n_off_series = int((~pts["t_index"].isin(layout.timepoints)).sum())
    assign: dict[tuple[int, int], int] = {}
    owner: dict[tuple[int, int], int] = {}
    n_conflict = 0
    n_duplicate = 0
    n_multi = 0
    seen_tt: set[tuple[int, int]] = set()
    for r in pts.dropna(subset=["nucleus_id"]).sort_values(["t_index", "track_id"]).itertuples(index=False):
        tt = (int(r.track_id), int(r.t_index))
        if tt in seen_tt:
            n_multi += 1                      # a second point of the same track in one frame (an Imaris merge, a typo): the first keeps it
            continue
        seen_tt.add(tt)
        key = (int(r.t_index), int(r.nucleus_id))
        if key in owner:
            if owner[key] != int(r.track_id):
                n_conflict += 1                   # two imported tracks in one nucleus at one timepoint: the first keeps it
            else:
                n_duplicate += 1
            continue
        owner[key] = int(r.track_id)
        assign[key] = int(r.track_id)
    info = {"n_points": int(len(pts)), "n_outside_any_label": n_outside, "n_off_series_timepoints": n_off_series,
            "n_conflicts_same_nucleus": n_conflict, "n_duplicate_points": n_duplicate,
            "n_multi_point_frames": n_multi,
            "n_imported_tracks": int(pts["track_id"].nunique())}
    return assign, info


def _assign_from_cellacdc_npz(layout: SeriesLayout, npz: Path, export: dict, spacing) -> tuple[dict, dict]:
    """Relabel a Cell-ACDC tracked ``segm.npz`` (the export's grid) by majority voxel overlap with our labels:
    every one of our exported labels takes the tracked value covering most of its voxels; a label whose
    majority is 0 was deleted (no track), a tracked value claimed by two of our labels in one frame is a merge
    (the larger overlap keeps it), a label split between tracked values is logged as a split."""
    with np.load(long_path(npz)) as arc:
        tracked = arc[arc.files[0]]
    tracked = np.asarray(tracked)
    if tracked.ndim == 3:
        tracked = tracked[:, None]
    frames = [int(t) for t in export["timepoints"]]
    if tracked.shape[0] != len(frames):
        raise ValueError(f"{npz.name}: {tracked.shape[0]} frames but the export has {len(frames)} timepoints")
    vsh = {int(k): np.asarray(v["voxels"], int) for k, v in export["shifts"].items()}
    box = export["box"]
    stride = int(export["stride"])
    crop = export.get("crop", "germline")
    assign: dict[tuple[int, int], int] = {}
    n_deleted = n_split = n_merged = n_added = 0
    for i, t in enumerate(frames):
        if t not in layout.timepoints:
            continue
        lab = layout.read_labels(t)
        ours = place_frame(lab, vsh.get(t, np.zeros(3, int)), box, stride, _keep_ids(layout, t, crop), dtype=np.int64)
        tr = tracked[i].astype(np.int64)
        if tr.shape != ours.shape:
            raise ValueError(f"{npz.name}: frame {i} has shape {tr.shape}, the export grid is {ours.shape}")
        idmap = {int(k): int(v) for k, v in export["label_maps"].get(str(t), {}).items()}
        maj_map = majority_map(ours, tr)
        if not maj_map:
            continue
        claims: dict[int, tuple[int, int]] = {}       # tracked id -> (our label, overlap)
        for a, (maj, cnt, lst) in maj_map.items():
            total = sum(c for _b, c in lst)
            if len(lst) > 1 and (total - cnt) / total > 0.05:
                n_split += 1
            if maj == 0:
                n_deleted += 1
                continue
            if maj in claims:
                n_merged += 1
                if cnt <= claims[maj][1]:
                    continue
                prev_label = claims[maj][0]
                assign.pop((t, idmap.get(prev_label, prev_label)), None)
            claims[maj] = (a, cnt)
            assign[(t, idmap.get(a, a))] = int(maj)
        n_added += int(len(set(np.unique(tr[ours == 0]).tolist()) - {0} - set(claims)))
    info = {"n_deleted_labels": n_deleted, "n_split_labels": n_split, "n_merged_labels": n_merged,
            "n_tracked_objects_without_our_label": n_added, "n_imported_tracks": len(set(assign.values()))}
    return assign, info


def import_tracks(series_dir: str | Path, source: str, file: str | Path, *, space: str = "raw", time_base: int = 0,
                  export_json: str | Path | None = None, snap_um: float = 0.0, germline_only: bool = True,
                  stem: str | None = None, write: bool = True) -> dict:
    """``germquant track-import``: a tracker's result becomes the two series tables with link_source
    ``imported:<source>``. ``source``: ``trackmate`` (spots CSV, or the model XML), ``imaris`` (a points table
    id, t, z, y, x in raw-frame microns; ``time_base`` 1 for Imaris's 1-based time index), ``cellacdc`` (a
    tracked ``segm.npz`` on the export grid; the export record is found next to it or given). ``space``:
    ``raw`` when the coordinates are raw-frame microns, ``export`` when they were measured on the exported
    stack. ``snap_um`` lets a point just outside every label take the nearest nucleus within that distance."""
    t_start = time.perf_counter()
    source = str(source).lower()
    if source not in IMPORT_SOURCES:
        raise ValueError(f"source must be one of {', '.join(IMPORT_SOURCES)}, got {source!r}")
    file = Path(file)
    layout = find_series(series_dir, stem)
    spacing = layout.spacing()
    shifts, shift_source = load_shifts(layout)
    objects = object_table(layout, shifts, germline_only=germline_only)
    export = None
    if source == "cellacdc" or space == "export":
        export = _find_export_record(file, export_json)
        if export is None:
            raise FileNotFoundError(f"no {EXPORT_RECORD} found near {file}; pass --export-json")
    if source == "cellacdc":
        assign, info = _assign_from_cellacdc_npz(layout, file, export, spacing)
    else:
        if source == "trackmate":
            pts = read_trackmate_xml(file) if file.suffix.lower() == ".xml" else read_trackmate_csv(file)
        else:
            pts = read_points_table(file, time_base=time_base)
        if space == "export":
            pts = _export_to_raw(pts, export)
        assign, info = _assign_from_points(layout, pts, spacing, snap_um=snap_um)
    link_source = f"imported:{source}"
    points, tracks, tinfo = build_track_tables(layout, objects, assign, link_source=link_source)
    record = {"series_id": layout.stem, "backend": link_source, "link_source": link_source, "source_file": str(file),
              "space": space, "time_base": int(time_base), "shift_source": shift_source, "germline_only": bool(germline_only),
              "timepoints": layout.timepoints, "n_objects": int(len(objects)),
              "n_tracks": int(tracks["time_track_id"].notna().sum()) if len(tracks) else 0,
              "n_full_length": int(tracks["is_full_length"].fillna(False).astype(bool).sum()) if len(tracks) else 0,
              "n_objects_unlinked": int(points["time_track_id"].isna().sum()),
              "backend_info": info, **tinfo,
              "elapsed_s": round(time.perf_counter() - t_start, 3), "written": provenance.run_timestamp()}
    paths = write_tracking_tables(layout, points, tracks, record) if write else {}
    log.info("%s: imported %s from %s: %d track(s), %d nucleus-timepoints unlinked", layout.stem, source, file.name,
             record["n_tracks"], record["n_objects_unlinked"])
    return {"points": points, "tracks": tracks, "record": record, "paths": paths, "layout": layout}


# =====================================================================================================
# track-strips
# =====================================================================================================

def track_strips(series_dir: str | Path, out_dir: str | Path, *, size_um: float | None = None, planes: int | None = None,
                 image: str = "ref", stride: int | None = None, track_ids=None, max_tracks: int | None = None,
                 stem: str | None = None, cfg=None) -> list[Path]:
    """``germquant track-strips``: one PNG per track under ``out_dir`` with a tile per timepoint of the series:
    the central_element channel around the track's nucleus (its raw centroid in that timepoint's frame, so
    the nucleus the track names is always centred), ``size_um`` wide, max-projected over ``planes`` planes
    around the centroid, intensity scaled per track so bleaching shows as dimming. A timepoint the track has
    no point in is a dark tile marked with a dash. ``image``: ``ref`` (the reference volumes the loop keeps) or
    ``file`` (the series file at ``stride``). Returns the PNG paths."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    layout = find_series(series_dir, stem)
    spacing = layout.spacing()
    pts = read_track_points(layout)
    size_um = float(size_um or cfg_get(cfg, "timelapse.track.strip_size_um", DEFAULTS["strip_size_um"]))
    planes = int(planes or cfg_get(cfg, "timelapse.track.strip_planes", DEFAULTS["strip_planes"]))
    tile_px = int(cfg_get(cfg, "timelapse.track.strip_tile_px", DEFAULTS["strip_tile_px"]))
    stride = int(stride or cfg_get(cfg, "timelapse.track.export_stride", DEFAULTS["export_stride"]))
    out_dir = Path(out_dir)
    os.makedirs(long_path(out_dir), exist_ok=True)

    # the image of every timepoint (one volume at a time) and the tiles of every track cut from it
    linked = pts[pts["time_track_id"].notna()]
    ids = sorted(int(v) for v in linked["time_track_id"].unique())
    if track_ids is not None:
        want = {int(v) for v in track_ids}
        ids = [i for i in ids if i in want]
    if max_tracks is not None:
        ids = ids[: int(max_tracks)]
    if not ids:
        log.warning("%s: no track to draw", layout.stem)
        return []
    tiles: dict[int, dict[int, np.ndarray | None]] = {i: {} for i in ids}
    role_idx = layout.channel_index("central_element")
    for t in layout.timepoints:
        if image == "ref":
            ref = layout.ref_path(t)
            if ref is None:
                raise SeriesLayoutError(f"{layout.stem}: no reference volume for t{t:03d}; run with --image file or rerun the "
                                        f"series with timelapse.keep_reference_volume")
            vol, r = np.load(long_path(ref[0]), mmap_mode="r"), ref[1]
        else:
            vol = _read_channel_from_file(layout, t, stride, role_idx)
            if vol is None:
                raise SeriesLayoutError(f"{layout.stem}: t{t:03d} could not be read from the series file")
            r = stride
        half = max(1, int(round(size_um / (2.0 * spacing[1] * r))))
        here = linked[linked["t_index"] == t]
        for row in here.itertuples(index=False):
            tid = int(row.time_track_id)
            if tid not in tiles:
                continue
            tiles[tid][int(t)] = _tile(vol, (float(row.z_um), float(row.y_um), float(row.x_um)),
                                       (spacing[0], spacing[1] * r, spacing[2] * r), half, planes)
    written = []
    T = layout.timepoints
    for tid in ids:
        strip = np.zeros((tile_px, len(T) * (tile_px + 2) - 2), dtype=float)
        have = [v for v in tiles[tid].values() if v is not None and v.size]
        lo, hi = (np.percentile(np.concatenate([v.ravel() for v in have]), [1, 99.5]) if have else (0.0, 1.0))
        hi = hi if hi > lo else lo + 1.0
        for j, t in enumerate(T):
            tile = tiles[tid].get(int(t))
            x0 = j * (tile_px + 2)
            if tile is None or not tile.size:
                strip[:, x0:x0 + tile_px] = 0.0
                strip[tile_px // 2 - 1: tile_px // 2 + 1, x0 + tile_px // 4: x0 + 3 * tile_px // 4] = 0.35   # a dash: no point here
                continue
            strip[:, x0:x0 + tile_px] = np.clip((_upscale(tile, tile_px) - lo) / (hi - lo), 0, 1)
        fig_w = max(2.0, 0.5 * len(T))
        fig, ax = plt.subplots(figsize=(fig_w, 0.9), dpi=100)
        ax.imshow(strip, cmap="gray", vmin=0, vmax=1, interpolation="nearest", aspect="equal")
        ax.set_xticks([j * (tile_px + 2) + tile_px / 2 for j in range(len(T))])
        ax.set_xticklabels([str(t) for t in T], fontsize=6)
        ax.set_yticks([])
        ax.set_title(f"{layout.stem}  track {tid}  ({size_um:g} um tiles, {planes} planes, {image})", fontsize=7, loc="left")
        for s in ax.spines.values():
            s.set_visible(False)
        p = out_dir / f"{layout.stem}__track{tid:05d}_strip.png"
        fig.savefig(long_path(p), bbox_inches="tight")
        plt.close(fig)
        written.append(p)
    log.info("%s: %d strip(s) -> %s", layout.stem, len(written), out_dir)
    return written


def _tile(vol, centre_um, spacing, half: int, planes: int) -> np.ndarray:
    """A (2 half + 1) square max projection over ``planes`` planes around the centre (clipped at the edges and
    zero-padded, so every tile has the same size)."""
    Z, Y, X = vol.shape
    cz, cy, cx = (int(round(c / s)) for c, s in zip(centre_um, spacing))
    hz = max(0, (int(planes) - 1) // 2)
    z0, z1 = max(0, cz - hz), min(Z, cz + hz + 1)
    out = np.zeros((2 * half + 1, 2 * half + 1), dtype=float)
    if z1 <= z0:
        return out
    y0, y1 = max(0, cy - half), min(Y, cy + half + 1)
    x0, x1 = max(0, cx - half), min(X, cx + half + 1)
    if y1 <= y0 or x1 <= x0:
        return out
    sub = np.asarray(vol[z0:z1, y0:y1, x0:x1]).max(axis=0).astype(float)
    out[y0 - (cy - half): y1 - (cy - half), x0 - (cx - half): x1 - (cx - half)] = sub
    return out


def _upscale(tile: np.ndarray, px: int) -> np.ndarray:
    """Nearest-neighbour resize of a square tile to ``px`` x ``px``."""
    n = tile.shape[0]
    idx = np.minimum((np.arange(px) * n) // px, n - 1)
    return tile[np.ix_(idx, idx)]
