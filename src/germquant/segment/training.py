"""The nucleus training loop (docs/TRAINING.md; docs/DESIGN_system_packs.md section 10).

These commands share this module:

* ``germquant export-training RESULTS_DIR --out PAIRS`` (``export_training_pairs``): for every finished run
  under RESULTS_DIR (one results folder, a time-lapse series folder with its ``t<NNN>/`` subfolders, or a
  batch tree) the stack is re-read from the source file at the stride, z range and timepoint recorded in
  that run's ``run_manifest.json``, and the pairs are written as ``image/<id>_<orient>_<index>.tif``
  (uint16) and ``labels/<id>_<orient>_<index>.tif`` (int32) with a ``manifest.csv``. What is exported
  follows the run's geometry (``--orientations auto``, the default): a 3D run's segmented channel is sliced
  in xy, xz and yz together with the run's own nucleus labels, and the xz and yz slices are resampled along
  z to the xy pixel size (linear interpolation for the image, nearest neighbour for the labels), so every
  exported slice is isotropic and a 2D model trained on them applies in 3D with ``anisotropy = dz / dy``,
  which is how the segment stage runs Cellpose; a projected run (segmentation.nuclei.project) gives one pair
  of orientation ``proj``, the 2D image its segmenter saw (the recorded plane, or the max or mean projection)
  with its 2D labels; a single-plane run (io.planar) gives its plane as orientation ``xy``, cut above
  TILE_PX (1024) pixels into the fewest tiles of at most that size, each tile's offset in the manifest
  (tile_y0, tile_x0). A person corrects the labels in napari or the Cellpose GUI and saves them back under
  the same name. ``germquant import-labels`` (germquant.import_labels) writes labels made elsewhere into the
  same layout, through the plane and slice writers here.
* ``germquant train-nuclei --pairs PAIRS --out MODELDIR`` (``train_nuclei``): ``cellpose.train.train_seg``
  in process, as scripts/train_nucleus_model_combined.py ran it (AdamW, learning rate 1e-5, weight decay
  0.1, batch 8), with the pairs split by source image or by series and never by slice, so no slice of a
  validation image, and no neighbouring timepoint of the same series, is seen in training. It writes the
  model file, ``training_log.json`` and the model card ``<name>.card.json`` with the sha256 of every pair,
  the hold-out metrics (germquant.validate.segmentation_metrics, F1 at IoU 0.5, per orientation and per
  zone) and what the data was: dims, pixel size, object kind, system, the pack recipe and the segmented role
  of the runs the pairs came from and, for public data, the dataset (docs/DESIGN_system_packs.md section 3.6).
* ``germquant register-model`` (models_registry.register_local_model) adds the model to the registry, after
  which ``segmentation.nuclei.cellpose_model: <name>`` in a profile runs it.

Nothing here touches the per-image tables or the fixed-image path.
"""
from __future__ import annotations

import hashlib
import inspect
import json
import logging
import os
import re
import time
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from .. import __version__
from ..fsutil import long_path

log = logging.getLogger(__name__)

ORIENTATIONS = ("xy", "xz", "yz")          # the slices of a 3D run; an explicit --orientations list names these
PROJ = "proj"                              # the one 2D image of a projected run (design section 10, step 2)
TILE_PX = 1024                             # a single plane above this size a side is exported as tiles of it
_AXIS = {"xy": 0, "xz": 1, "yz": 2}        # the array axis a slice index runs along (labels are (Z, Y, X))
_T_SUFFIX = re.compile(r"_t(\d{3})$")
_PAIR_NAME = re.compile(r"^(?P<id>.+)_(?P<orient>xy|xz|yz|proj)_(?P<index>\d+)$")
MANIFEST_COLS = ["image_id", "series_id", "t_index", "orientation", "slice_index", "n_labels", "zone",
                 "image", "labels", "source_file", "xy_stride", "z_range", "channel_role", "channel_index",
                 "pixel_um", "z_zoom", "height", "width"]
# columns a manifest gains only when a pair needs them, after MANIFEST_COLS: source (import:<tool>, a pair of
# germquant import-labels; empty for an exported pair) and tile_y0, tile_x0 (the offset of a 2D plane's pair in
# its plane, on the run's grid: a single-plane run's export, or an import of 2D labels). A manifest of 3D and
# projected runs only has exactly MANIFEST_COLS, as before.
MANIFEST_EXTRA_COLS = ["source", "tile_y0", "tile_x0"]
_TILE_COLS = ("tile_y0", "tile_x0")
LICENCE_NOTE = ("Cellpose-SAM weights are BSD-3-Clause; Cellpose's own training data is distributed under a "
                "non-commercial (CC BY-NC) licence, so a model fine-tuned from cpsam inherits that question for "
                "commercial use (README, License).")


# ---- discovering finished runs -------------------------------------------------------------------------

@dataclass
class RunRef:
    """One finished run (a results folder with its nucleus labels), as export-training sees it."""
    folder: str
    image_id: str
    series_id: str
    source_file: str
    xy_stride: int
    z_range: list[int] | None
    t_index: int | None
    channel_role: str
    channel_index: int | None
    labels_path: str
    nuclei_path: str | None


def series_id_of(image_id: str, t_index) -> str:
    """The series an image id belongs to: the id without its ``_t<NNN>`` suffix when it is a timepoint of a
    series (``t_index`` given and matching the suffix), else the id itself."""
    if t_index is None or (isinstance(t_index, float) and np.isnan(t_index)):
        return image_id
    suffix = f"_t{int(t_index):03d}"
    return image_id[:-len(suffix)] if image_id.endswith(suffix) else image_id


def find_runs(root: str | Path) -> list[RunRef]:
    """Every finished run under ``root``: a folder holding ``run_manifest.json``, ``<id>__nuclei_labels.tif``
    and ``<id>__image_summary.csv`` (a single results folder, the ``t<NNN>/`` folders of a series, or a batch
    tree). The source file is the ``file_path`` the run stamped on its tables, the geometry (xy_stride,
    z_range, t_index) and the segmented channel role (segmentation.nuclei.channel, default dna) come from the
    manifest, and the channel index from the run's ``channel_map_json``. Sorted by series, timepoint, id."""
    root = Path(root)
    if not root.is_dir():
        raise FileNotFoundError(f"results folder not found: {root}")
    runs: list[RunRef] = []
    root_lp = long_path(root)
    for dirpath, _dirs, names in os.walk(root_lp):
        rel = os.path.relpath(dirpath, root_lp)
        folder = root if rel == "." else root / rel
        if "run_manifest.json" not in names and not any(n.endswith("__stages.json") for n in names):
            continue
        for n in sorted(names):
            if not n.endswith("__nuclei_labels.tif"):
                continue
            iid = n[:-len("__nuclei_labels.tif")]
            summary = folder / f"{iid}__image_summary.csv"
            if not summary.exists():
                log.warning("%s: labels without an image_summary table; skipped", folder / n)
                continue
            try:
                runs.append(_run_ref(folder, iid, folder / n, summary))
            except Exception as e:  # noqa: BLE001 - one unreadable run must not stop the export
                log.warning("%s: unreadable run (%s: %s); skipped", folder, type(e).__name__, e)
    runs.sort(key=lambda r: (r.series_id, -1 if r.t_index is None else r.t_index, r.image_id))
    return runs


def _run_ref(folder: Path, iid: str, labels_path: Path, summary_path: Path) -> RunRef:
    man_path = folder / "run_manifest.json"
    if man_path.exists():
        with open(long_path(man_path), encoding="utf-8") as fh:
            man = json.load(fh)
    else:
        # germquant batch writes one manifest at the batch root; the run's own geometry (xy_stride, z_range,
        # t_index) is in its stage record, the segmented role defaults to dna
        with open(long_path(folder / f"{iid}__stages.json"), encoding="utf-8") as fh:
            rec = json.load(fh)
        man = {k: rec.get(k) for k in ("xy_stride", "z_range", "t_index") if rec.get(k) is not None}
        for nested in ("run_geometry", "geometry"):
            if isinstance(rec.get(nested), dict):
                man.update({k: v for k, v in rec[nested].items() if k in ("xy_stride", "z_range", "t_index")})
        if isinstance(rec.get("segment_projection"), dict) or rec.get("planar") is True:
            # a projected or single-plane run of a batch: its 2D pair must be the channel the segmenter saw, so
            # the segmented role comes from the batch's own config (3D batch runs keep the dna default, as before)
            batch_man = _batch_manifest(folder, rec.get("config_hash"))
            if batch_man is not None:
                man["config"] = batch_man.get("config")
    s = pd.read_csv(long_path(summary_path), low_memory=False)
    row = s.iloc[0] if len(s) else pd.Series(dtype=object)
    src = row.get("file_path")
    if not isinstance(src, str) or not src:
        raise ValueError("image_summary has no file_path (the source image is unknown)")
    cfg = man.get("config") or {}
    role = str(((cfg.get("segmentation") or {}).get("nuclei") or {}).get("channel") or "dna")
    cm = row.get("channel_map_json")
    idx = None
    if isinstance(cm, str) and cm:
        try:
            idx = json.loads(cm).get(role)
        except ValueError:
            idx = None
    t_index = man.get("t_index")
    zr = man.get("z_range")
    nuclei = folder / f"{iid}__nuclei.csv"
    return RunRef(folder=str(folder), image_id=iid, series_id=series_id_of(iid, t_index), source_file=src,
                  xy_stride=int(man.get("xy_stride", 1) or 1), z_range=[int(v) for v in zr] if zr else None,
                  t_index=None if t_index is None else int(t_index), channel_role=role,
                  channel_index=None if idx is None else int(idx), labels_path=str(labels_path),
                  nuclei_path=str(nuclei) if nuclei.exists() else None)


def _batch_manifest(folder: Path, config_hash) -> dict | None:
    """The run_manifest.json of the batch a run folder belongs to: the nearest one above the folder, when it
    records the run's config_hash; else None."""
    if not config_hash:
        return None
    for parent in Path(os.path.abspath(folder)).parents:      # absolute: a relative "." has no parents
        p = parent / "run_manifest.json"
        if os.path.exists(long_path(p)):
            try:
                with open(long_path(p), encoding="utf-8") as fh:
                    m = json.load(fh)
            except (OSError, ValueError):
                return None
            return m if isinstance(m, dict) and m.get("config_hash") == config_hash else None
    return None


def _read_json(path: str | Path) -> dict:
    with open(long_path(path), encoding="utf-8") as fh:
        return json.load(fh)


def _run_record(folder: str | Path, image_id: str) -> dict:
    """The record a finished run's geometry is read from, the file find_runs reads: the folder's
    ``run_manifest.json``, else the run's ``<id>__stages.json`` (a batch writes its manifest at the batch root
    only). Keys a stage record nests under run_geometry or geometry are lifted to the top when absent there."""
    folder = Path(folder)
    man = folder / "run_manifest.json"
    rec = _read_json(man if os.path.exists(long_path(man)) else folder / f"{image_id}__stages.json")
    for nested in ("run_geometry", "geometry"):
        if isinstance(rec.get(nested), dict):
            for k, v in rec[nested].items():
                rec.setdefault(k, v)
    return rec


def _geometry(record: dict) -> tuple[str, dict | None]:
    """``("planar", None)`` for a single-plane run (its record holds ``"planar": true``, docs/PLANAR.md; checked
    first, since such a run records no segment_projection, decision D7), ``("projected", {mode, plane, z})`` for a
    run whose record holds segment_projection (a projected run), else ``("3d", None)``: every other run exports
    its xy, xz and yz slices."""
    if record.get("planar") is True:
        return "planar", None
    proj = record.get("segment_projection")
    if isinstance(proj, dict):
        return "projected", dict(proj)
    return "3d", None


def _run_spacing(run: RunRef, record: dict, stack) -> tuple[float, float, float]:
    """The (dz, dy, dx) a run's pixels are exported on: the stack's as read, unless the run records spacing_source
    (a single-plane run, or io.voxel_size_um in any run: decision D6), whose spacing is applied after the read and
    so is not the file's. Then the voxel size of the run's image_summary, which process_image writes from the
    spacing it ran with; a single plane's NaN dz becomes planar.PLANAR_DZ, the placeholder its run used."""
    if not record.get("spacing_source"):
        return tuple(float(s) for s in stack.spacing)
    s = pd.read_csv(long_path(Path(run.folder) / f"{run.image_id}__image_summary.csv"), low_memory=False,
                    float_precision="round_trip")
    row = s.iloc[0]
    dz, dy, dx = (float(row[f"voxel_{a}_um"]) for a in ("dz", "dy", "dx"))
    if not np.isfinite(dz):
        from ..planar import PLANAR_DZ

        dz = PLANAR_DZ
    if not (np.isfinite(dy) and np.isfinite(dx) and dy > 0 and dx > 0):
        raise ValueError(f"{run.image_id}: the run's image_summary holds no pixel size (voxel_dy_um {dy}, "
                         f"voxel_dx_um {dx}), and the run's spacing ({record.get('spacing_source')}) is not the file's")
    return dz, dy, dx


def plane_tiles(height: int, width: int, tile_px: int = TILE_PX) -> list[tuple[int, int, int, int]]:
    """The tiles a plane of ``height`` x ``width`` pixels is exported as: ``(y0, x0, h, w)`` row by row, the fewest
    tiles of at most ``tile_px`` pixels a side, of near-equal size (on each axis the sizes differ by at most one
    pixel, so no thin sliver is left at an edge), every pixel in exactly one tile. A plane that fits is one tile,
    ``(0, 0, height, width)``."""
    def cuts(n: int) -> list[tuple[int, int]]:
        n = int(n)
        k = max(1, -(-n // int(tile_px)))
        edges = [(i * n) // k for i in range(k + 1)]
        return [(edges[i], edges[i + 1] - edges[i]) for i in range(k)]

    return [(y0, x0, h, w) for y0, h in cuts(height) for x0, w in cuts(width)]


# ---- slicing ---------------------------------------------------------------------------------------------

def slice_counts(labels: np.ndarray, orientation: str) -> list[int]:
    """The number of distinct nuclei each slice of the orientation cuts."""
    ax = _AXIS[orientation]
    out = []
    for i in range(labels.shape[ax]):
        sl = np.take(labels, i, axis=ax)
        out.append(int(np.count_nonzero(np.unique(sl))))
    return out


def choose_slices(counts: list[int], n: int, min_labels: int, rng: np.random.Generator) -> list[int]:
    """``n`` slice indices drawn at random (seeded) from the slices that cut at least ``min_labels`` nuclei;
    when fewer than ``n`` qualify, from every slice that cuts at least one. Sorted; empty when none does."""
    cands = [i for i, c in enumerate(counts) if c >= max(1, int(min_labels))]
    if len(cands) < n:
        cands = [i for i, c in enumerate(counts) if c >= 1]
    if not cands:
        return []
    k = min(int(n), len(cands))
    return sorted(int(v) for v in rng.choice(cands, size=k, replace=False))


def isotropic_slice(img: np.ndarray, labels: np.ndarray, orientation: str, index: int, spacing) -> tuple[np.ndarray, np.ndarray, float, float]:
    """The 2D image and label slice at ``index`` along the orientation's axis, resampled to the xy pixel
    size (``dy``): an xz or yz slice is stretched along z by ``dz / dy`` (linear for the image, nearest
    neighbour for the labels); an xy slice is left alone unless dx differs from dy. Returns
    ``(image2d, labels2d, pixel_um, z_zoom)``."""
    from scipy import ndimage as ndi

    ax = _AXIS[orientation]
    im = np.take(np.asarray(img), index, axis=ax).astype(np.float32)
    lab = np.take(np.asarray(labels), index, axis=ax)
    dz, dy, dx = (float(s) for s in spacing)
    px = dy
    if orientation == "xy":
        zoom = (1.0, dx / dy)
        z_zoom = 1.0
    else:                                  # rows are z in both xz (Z, X) and yz (Z, Y)
        z_zoom = dz / dy
        zoom = (z_zoom, (dx / dy) if orientation == "xz" else 1.0)
    if any(abs(f - 1.0) > 1e-9 for f in zoom):
        im = ndi.zoom(im, zoom, order=1, mode="nearest")
        lab = ndi.zoom(lab, zoom, order=0, mode="nearest")
    return im, np.asarray(lab).astype(np.int32), px, z_zoom


def _to_uint16(im: np.ndarray) -> np.ndarray:
    return np.clip(np.rint(np.asarray(im, dtype=np.float64)), 0, 65535).astype(np.uint16)


def _mode_zone(zones: dict, label_ids) -> str:
    present = [zones[int(v)] for v in label_ids if int(v) in zones]
    if not present:
        return ""
    counts = Counter(present)
    best = max(counts.values())
    return sorted(z for z, c in counts.items() if c == best)[0]


def _rng_for(seed: int, image_id: str) -> np.random.Generator:
    """A generator seeded by the export seed and the image id, so the slices chosen for one image do not
    depend on which other images are exported alongside it."""
    h = int(hashlib.sha256(image_id.encode("utf-8")).hexdigest()[:8], 16)
    return np.random.default_rng([int(seed), h])


# ---- export ----------------------------------------------------------------------------------------------

def _listed_orientations(orientations) -> tuple[str, ...] | None:
    """None for ``auto`` (None or the text "auto"); else the normalised list (a sequence, or a text with
    commas), which names slices of 3D runs and keeps its order: the order decides the seeded draw."""
    if orientations is None:
        return None
    if isinstance(orientations, str):
        if orientations.strip().lower() == "auto":
            return None
        orientations = orientations.split(",")
    listed = tuple(str(o).strip().lower() for o in orientations if str(o).strip())
    bad = [o for o in listed if o not in ORIENTATIONS]
    if bad or not listed:
        raise ValueError(f"orientations must be 'auto' or a list among {ORIENTATIONS} (the slices of 3D runs; a "
                         f"projected run always gives its {PROJ} pair), got {listed}")
    return listed


def _projected_pair(run: RunRef, projection: dict | None, img: np.ndarray, labels: np.ndarray, spacing,
                    channel_index: int, zones: dict, out_dir: Path, slices_per_image: int):
    """The one pair of a projected run: the 2D image its segmenter saw (planar.project_image on the re-read
    channel, with the recorded plane z in mode plane) and its 2D labels, labels[0] (the labels are copied
    through every plane). Written as ``<id>_proj_<z>`` (mode plane) or ``<id>_proj_0000`` (max, mean), uint16
    like every exported image and unresampled: the image as the model saw it. ``(row, bytes)``, or None when
    the labels hold no object or ``slices_per_image`` is below 1."""
    import tifffile

    from .. import planar

    projection = projection or {}
    mode = str(projection.get("mode") or "").strip().lower()
    if mode not in ("max", "mean", "plane"):
        raise ValueError(f"{run.image_id}: the run's segment_projection {projection!r} names no projection mode "
                         f"(max, mean or plane)")
    z = None
    if mode == "plane":
        z = projection.get("z")
        if z is None:                     # every projected run records z; resolve it as the run did otherwise
            z = planar.resolve_plane(img, projection.get("plane") if projection.get("plane") is not None else "mid")
        z = int(z)
    img2d = planar.project_image(img, mode, z)
    lab2d = np.asarray(labels[0]).astype(np.int32)
    ids = np.unique(lab2d)
    ids = ids[ids > 0]
    if int(slices_per_image) < 1 or not ids.size:
        return None
    index = z if mode == "plane" else 0
    stem = f"{run.image_id}_{PROJ}_{index:04d}.tif"
    p_im, p_lab = out_dir / "image" / stem, out_dir / "labels" / stem
    tifffile.imwrite(long_path(p_im), _to_uint16(img2d), compression="zlib")
    tifffile.imwrite(long_path(p_lab), lab2d, compression="zlib")
    row = {"image_id": run.image_id, "series_id": run.series_id, "t_index": run.t_index, "orientation": PROJ,
           "slice_index": int(index), "n_labels": int(ids.size), "zone": _mode_zone(zones, ids),
           "image": f"image/{stem}", "labels": f"labels/{stem}", "source_file": run.source_file,
           "xy_stride": run.xy_stride, "z_range": json.dumps(run.z_range) if run.z_range else "",
           "channel_role": run.channel_role, "channel_index": int(channel_index), "pixel_um": float(spacing[1]),
           "z_zoom": 1.0, "height": int(lab2d.shape[0]), "width": int(lab2d.shape[1])}
    return row, os.path.getsize(long_path(p_im)) + os.path.getsize(long_path(p_lab))


def plane_pairs(base: dict, img2d: np.ndarray, lab2d: np.ndarray, zones: dict, out_dir: Path,
                tile_px: int = TILE_PX):
    """Write the pairs of one 2D plane, the image a single-plane run segmented (or an imported 2D labels file's
    image, germquant.import_labels): the plane itself, or above ``tile_px`` pixels a side its tiles (plane_tiles),
    as ``<id>_xy_<k>.tif`` with k the tile's number row by row (0 for a whole plane), uint16 like every exported
    image and not resampled, the labels int32 with their ids as given. A tile whose labels hold no object is not
    written (its number is skipped). Yields ``(row, bytes written)`` per pair, writing each when asked, so the
    caller can stop at its size cap; ``base`` holds the row's fields that do not change per tile (image_id,
    series_id, t_index, source_file, xy_stride, z_range, channel_role, channel_index, pixel_um and, for an import,
    source)."""
    import tifffile

    img2d, lab2d = np.asarray(img2d), np.asarray(lab2d)
    if img2d.shape != lab2d.shape or lab2d.ndim != 2:
        raise ValueError(f"{base.get('image_id')}: a plane and its labels must be 2D of one shape, got {img2d.shape} "
                         f"and {lab2d.shape}")
    for k, (y0, x0, h, w) in enumerate(plane_tiles(*lab2d.shape, tile_px)):
        lab = np.ascontiguousarray(lab2d[y0:y0 + h, x0:x0 + w]).astype(np.int32)
        ids = np.unique(lab)
        ids = ids[ids > 0]
        if not ids.size:
            continue
        stem = f"{base['image_id']}_xy_{k:04d}.tif"
        p_im, p_lab = out_dir / "image" / stem, out_dir / "labels" / stem
        tifffile.imwrite(long_path(p_im), _to_uint16(img2d[y0:y0 + h, x0:x0 + w]), compression="zlib")
        tifffile.imwrite(long_path(p_lab), lab, compression="zlib")
        row = {**base, "orientation": "xy", "slice_index": int(k), "n_labels": int(ids.size),
               "zone": _mode_zone(zones, ids), "image": f"image/{stem}", "labels": f"labels/{stem}",
               "z_zoom": 1.0, "height": int(h), "width": int(w), "tile_y0": int(y0), "tile_x0": int(x0)}
        yield row, os.path.getsize(long_path(p_im)) + os.path.getsize(long_path(p_lab))


def slice_pairs(base: dict, img: np.ndarray, labels: np.ndarray, spacing, out_dir: Path, *, orientations=ORIENTATIONS,
                slices_per_image: int = 3, min_labels: int = 3, seed: int = 0, zones: dict | None = None):
    """Write the xy, xz and yz slice pairs of a (Z, Y, X) image and its labels by export-training's rule for a 3D
    run (choose_slices, seeded per image with _rng_for, isotropic_slice), for germquant.import_labels. Yields
    ``(row, bytes written)`` per pair; ``base`` as for plane_pairs. ``pixel_um`` is taken from ``base`` when it
    holds one (None there for an image whose voxel size is unknown), else the slice's."""
    import tifffile

    rng = _rng_for(seed, str(base["image_id"]))
    zones = zones or {}
    for orient in orientations:
        for i in choose_slices(slice_counts(labels, orient), slices_per_image, min_labels, rng):
            im2, lab2, px, z_zoom = isotropic_slice(img, labels, orient, i, spacing)
            stem = f"{base['image_id']}_{orient}_{i:04d}.tif"
            p_im, p_lab = out_dir / "image" / stem, out_dir / "labels" / stem
            tifffile.imwrite(long_path(p_im), _to_uint16(im2), compression="zlib")
            tifffile.imwrite(long_path(p_lab), lab2, compression="zlib")
            ids = np.unique(lab2)
            ids = ids[ids > 0]
            row = {"pixel_um": px, **base, "orientation": orient, "slice_index": int(i), "n_labels": int(ids.size),
                   "zone": _mode_zone(zones, ids), "image": f"image/{stem}", "labels": f"labels/{stem}",
                   "z_zoom": z_zoom, "height": int(lab2.shape[0]), "width": int(lab2.shape[1])}
            yield row, os.path.getsize(long_path(p_im)) + os.path.getsize(long_path(p_lab))


def merge_manifest(out_dir: Path, rows: list[dict]) -> pd.DataFrame:
    """The manifest of a pairs folder after ``rows`` join it: the rows of an existing manifest.csv whose image is
    not among ``rows`` stay (in their order), then ``rows``. The columns are MANIFEST_COLS, then those of
    MANIFEST_EXTRA_COLS that a row of either side holds (none for 3D and projected runs, so such a manifest is
    exactly MANIFEST_COLS); the tile offsets are whole numbers, empty where a pair has none. Nothing is written."""
    extra = [c for c in MANIFEST_EXTRA_COLS if any(c in r for r in rows)]
    new = pd.DataFrame(rows, columns=MANIFEST_COLS + extra)
    manifest_path = out_dir / "manifest.csv"
    if manifest_path.exists():
        old = pd.read_csv(long_path(manifest_path), low_memory=False)
        old = old[~old["image_id"].astype(str).isin(set(new["image_id"].astype(str)))] if "image_id" in old else old.iloc[0:0]
        # an extra column stays while a kept row holds a value in it
        cols = MANIFEST_COLS + [c for c in MANIFEST_EXTRA_COLS
                                if c in extra or (c in old.columns and old[c].notna().any())]
        for c in cols:
            if c not in old.columns:
                old[c] = pd.NA
            if c not in new.columns:
                new[c] = pd.NA
        new = pd.concat([old[cols], new[cols]], ignore_index=True)
    for c in _TILE_COLS:
        if c in new.columns:
            new[c] = pd.array(pd.to_numeric(new[c], errors="coerce"), dtype="Int64")
    return new


def export_training_pairs(results_dir: str | Path, out_dir: str | Path, *, orientations=None,
                          slices_per_image: int = 3, seed: int = 0, zone_column: str | None = "zone",
                          min_labels: int = 3, max_bytes: float = 2e9) -> pd.DataFrame:
    """Write paired 2D training images from every finished run under ``results_dir`` (module docstring).

    What a run gives follows its geometry (``orientations`` None or "auto", the default; docs/TRAINING.md):

    * a 3D run: per orientation, ``slices_per_image`` slices drawn at random (seeded per image) from those
      that cut at least ``min_labels`` nuclei; each is written as ``image/<id>_<orient>_<index>.tif`` (uint16,
      the raw intensities) and ``labels/<id>_<orient>_<index>.tif`` (int32, the run's labels, resampled with
      the image). The orientations are xy, xz and yz, or the ones an explicit ``orientations`` list names,
      in its order;
    * a projected run (its record holds segment_projection): one pair of orientation ``proj``, the 2D image
      the segmenter saw and the run's 2D labels (``_projected_pair``), whatever an explicit list names;
    * a single-plane run (its record holds ``"planar": true``): its plane as orientation ``xy``, whole, or above
      TILE_PX pixels a side the fewest tiles of at most that size (``plane_pairs``), whatever an explicit list
      names; each row records the tile's offset (tile_y0, tile_x0), and its slice_index is the tile's number.

    A run that records spacing_source (a single-plane run, or io.voxel_size_um in any run) is exported on the
    spacing it ran with, from its image_summary (``_run_spacing``): the xz and yz stretch and pixel_um are the
    run's, not the file's. ``zone_column`` names the nuclei-table column (the staging stage's ``zone``) whose
    most frequent value among the nuclei in a pair becomes its zone tag; absent column, empty tag. The export
    stops with a warning once ``max_bytes`` have been written (a 56 GB series must not become a 56 GB training
    folder). ``manifest.csv`` is merged with an existing one in ``out_dir`` (``merge_manifest``): rows of the
    images exported now replace their old rows, other rows stay; its columns are MANIFEST_COLS, plus the tile
    and import columns of MANIFEST_EXTRA_COLS only when a row holds them. ``export.json`` records the settings
    and the runs; its ``orientations`` is the union over the runs (the slice list, then xy for single-plane runs
    when the list lacks it, then ``proj``), and only when a run is not 3D, ``run_orientations`` gives each run's
    geometry and orientations, in the order of ``runs``. Returns the manifest."""
    import tifffile

    from .. import pipeline

    listed = _listed_orientations(orientations)
    orientations = ORIENTATIONS if listed is None else listed       # the slices of the 3D runs
    runs = find_runs(results_dir)
    if not runs:
        raise FileNotFoundError(f"no finished run (run_manifest.json + <id>__nuclei_labels.tif + "
                                f"<id>__image_summary.csv) under {results_dir}")
    records = [_run_record(r.folder, r.image_id) for r in runs]
    geometry = [_geometry(rec) for rec in records]
    out_dir = Path(out_dir)
    for sub in ("image", "labels"):
        os.makedirs(long_path(out_dir / sub), exist_ok=True)
    rows: list[dict] = []
    written = 0
    stopped = False
    for run, record, (kind, projection) in zip(runs, records, geometry):
        _meta, stack = pipeline._read_input(Path(run.source_file), xy_stride=run.xy_stride,
                                            z_range=tuple(run.z_range) if run.z_range else None, t=run.t_index)
        labels = np.asarray(tifffile.imread(long_path(run.labels_path)))
        if kind == "planar" and labels.ndim == 2:
            labels = labels[None]               # a single-plane run writes its labels as YX (planar.write_labels_yx)
        # the spacing the run ran with: the stack's as read, unless the run records spacing_source
        spacing = _run_spacing(run, record, stack)
        idx = run.channel_index
        if idx is None:
            log.warning("%s: the segmented role %r did not resolve in that run; exporting channel 0 as the "
                        "segment stage did", run.image_id, run.channel_role)
            idx = 0
        img = np.asarray(stack.data[int(idx)])
        if img.shape != labels.shape:
            raise ValueError(f"{run.image_id}: the re-read stack {img.shape} does not match the labels "
                             f"{labels.shape}; the run's stride or z range in run_manifest.json does not "
                             f"describe {run.source_file}")
        zones: dict[int, str] = {}
        if zone_column and run.nuclei_path:
            nuc = pd.read_csv(long_path(run.nuclei_path), low_memory=False)
            if zone_column in nuc.columns:
                z = nuc[["nucleus_id", zone_column]].dropna()
                zones = {int(k): str(v) for k, v in zip(z["nucleus_id"], z[zone_column])}
        rng = _rng_for(seed, run.image_id)
        slice_orients = orientations if kind == "3d" else ()        # a projected or single-plane run is not sliced
        if kind == "projected":
            pair = _projected_pair(run, projection, img, labels, spacing, int(idx), zones, out_dir,
                                   slices_per_image)
            if pair is not None:
                rows.append(pair[0])
                written += pair[1]
                stopped = written > float(max_bytes)
        elif kind == "planar" and int(slices_per_image) >= 1:
            # the plane the run segmented, whole or in tiles above TILE_PX (design section 10, step 2)
            base = {"image_id": run.image_id, "series_id": run.series_id, "t_index": run.t_index,
                    "source_file": run.source_file, "xy_stride": run.xy_stride,
                    "z_range": json.dumps(run.z_range) if run.z_range else "", "channel_role": run.channel_role,
                    "channel_index": int(idx), "pixel_um": float(spacing[1])}
            for row, nbytes in plane_pairs(base, img[0], labels[0], zones, out_dir):
                rows.append(row)
                written += nbytes
                if written > float(max_bytes):
                    stopped = True
                    break
        for orient in slice_orients:
            counts = slice_counts(labels, orient)
            for i in choose_slices(counts, slices_per_image, min_labels, rng):
                im2, lab2, px, z_zoom = isotropic_slice(img, labels, orient, i, spacing)
                stem = f"{run.image_id}_{orient}_{i:04d}.tif"
                p_im, p_lab = out_dir / "image" / stem, out_dir / "labels" / stem
                tifffile.imwrite(long_path(p_im), _to_uint16(im2), compression="zlib")
                tifffile.imwrite(long_path(p_lab), lab2, compression="zlib")
                written += p_im.stat().st_size + p_lab.stat().st_size
                ids = np.unique(lab2)
                ids = ids[ids > 0]
                rows.append({
                    "image_id": run.image_id, "series_id": run.series_id, "t_index": run.t_index,
                    "orientation": orient, "slice_index": int(i), "n_labels": int(ids.size),
                    "zone": _mode_zone(zones, ids), "image": f"image/{stem}", "labels": f"labels/{stem}",
                    "source_file": run.source_file, "xy_stride": run.xy_stride,
                    "z_range": json.dumps(run.z_range) if run.z_range else "", "channel_role": run.channel_role,
                    "channel_index": int(idx), "pixel_um": px, "z_zoom": z_zoom,
                    "height": int(lab2.shape[0]), "width": int(lab2.shape[1]),
                })
                if written > float(max_bytes):
                    stopped = True
                    break
            if stopped:
                break
        log.info("%s: %d slice pair(s) exported", run.image_id, sum(r["image_id"] == run.image_id for r in rows))
        if stopped:
            log.warning("export stopped at %.2f GB (max_bytes); %d of %d runs exported", written / 1e9,
                        len({r["image_id"] for r in rows}), len(runs))
            break
    new = merge_manifest(out_dir, rows)
    new.to_csv(long_path(out_dir / "manifest.csv"), index=False)
    kinds = {kind for kind, _ in geometry}
    # the union over the runs: the slice list as given (a 3D-only tree records exactly that), then the xy of
    # single-plane runs when the list lacks it, then proj
    resolved = (list(orientations) if "3d" in kinds else [])
    if "planar" in kinds and "xy" not in resolved:
        resolved.append("xy")
    resolved += [PROJ] if "projected" in kinds else []
    per_run = {}
    if kinds != {"3d"}:
        # each run's geometry and orientations, in the order of "runs"; written only when a run is not 3D
        per_run["run_orientations"] = [
            {"image_id": r.image_id, "geometry": kind,
             "orientations": list(orientations) if kind == "3d" else (["xy"] if kind == "planar" else [PROJ]),
             **({"segment_projection": projection} if projection is not None else {}),
             **({"tile_px": TILE_PX, "spacing_source": rec.get("spacing_source")} if kind == "planar" else {})}
            for r, rec, (kind, projection) in zip(runs, records, geometry)]
    with open(long_path(out_dir / "export.json"), "w", encoding="utf-8") as fh:
        json.dump({"results_dir": str(results_dir), "orientations": resolved, **per_run,
                   "slices_per_image": int(slices_per_image), "seed": int(seed), "zone_column": zone_column,
                   "min_labels": int(min_labels), "max_bytes": float(max_bytes), "bytes_written": int(written),
                   "stopped_at_cap": stopped, "runs": [asdict(r) for r in runs], "germquant_version": __version__,
                   "date": datetime.now(timezone.utc).isoformat()}, fh, indent=1, default=str)
    log.info("%d slice pair(s) from %d run(s) -> %s", len(rows), len({r["image_id"] for r in rows}), out_dir)
    return new


# ---- pairs and the split ---------------------------------------------------------------------------------

def load_pairs(pairs_dir: str | Path) -> list[dict]:
    """The slice pairs in a training folder: from ``manifest.csv`` when present (pairs whose files were
    deleted are skipped with a log line), else every ``image/*.tif`` with a ``labels/`` file of the same
    name, the id, orientation (xy, xz, yz or proj) and index parsed from the name. Each pair: ``image``,
    ``labels`` (absolute paths), ``image_id``, ``series_id``, ``orientation``, ``slice_index``, ``zone``."""
    pairs_dir = Path(pairs_dir)
    manifest = pairs_dir / "manifest.csv"
    pairs: list[dict] = []
    if manifest.exists():
        df = pd.read_csv(long_path(manifest), low_memory=False)
        for r in df.to_dict("records"):
            p_im, p_lab = pairs_dir / str(r["image"]), pairs_dir / str(r["labels"])
            if not (p_im.exists() and p_lab.exists()):
                log.info("pair %s dropped: its files are gone", r["image"])
                continue
            zone = r.get("zone")
            sid = r.get("series_id")
            if sid is None or (isinstance(sid, float) and np.isnan(sid)):
                sid = series_id_of(str(r["image_id"]), r.get("t_index"))
            pairs.append({"image": str(p_im), "labels": str(p_lab), "image_id": str(r["image_id"]),
                          "series_id": str(sid),
                          "orientation": str(r["orientation"]), "slice_index": int(r["slice_index"]),
                          "zone": "" if zone is None or (isinstance(zone, float) and np.isnan(zone)) else str(zone)})
        return pairs
    for p_im in sorted((pairs_dir / "image").glob("*.tif")):
        p_lab = pairs_dir / "labels" / p_im.name
        if not p_lab.exists():
            log.info("%s has no labels file; skipped", p_im.name)
            continue
        m = _PAIR_NAME.match(p_im.stem)
        iid = m.group("id") if m else p_im.stem
        pairs.append({"image": str(p_im), "labels": str(p_lab), "image_id": iid,
                      "series_id": _T_SUFFIX.sub("", iid), "orientation": m.group("orient") if m else "xy",
                      "slice_index": int(m.group("index")) if m else -1, "zone": ""})
    return pairs


def split_pairs(pairs: list[dict], by: str = "series", val_frac: float = 0.2, seed: int = 0) -> tuple[list[dict], list[dict], dict]:
    """Hold out whole groups, never slices: ``by`` is ``series`` (the default, so the timepoints of one
    series stay together) or ``image``. ``val_frac`` of the groups (at least one, at most all but one) go to
    validation; a single group trains alone and the hold-out set is empty. Returns (train, val, record)."""
    if by not in ("series", "image"):
        raise ValueError(f"split must be 'series' or 'image', got {by!r}")
    key = "series_id" if by == "series" else "image_id"
    groups = sorted({str(p[key]) for p in pairs})
    rng = np.random.default_rng(int(seed))
    perm = rng.permutation(len(groups))
    n_val = 0 if len(groups) < 2 else min(len(groups) - 1, max(1, int(round(float(val_frac) * len(groups)))))
    val_groups = sorted(groups[i] for i in perm[:n_val])
    val = [p for p in pairs if str(p[key]) in set(val_groups)]
    train = [p for p in pairs if str(p[key]) not in set(val_groups)]
    record = {"by": by, "key": key, "val_frac": float(val_frac), "seed": int(seed), "n_groups": len(groups),
              "train_groups": [g for g in groups if g not in set(val_groups)], "val_groups": val_groups,
              "n_train": len(train), "n_val": len(val)}
    return train, val, record


# ---- training --------------------------------------------------------------------------------------------

def _sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(long_path(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def _resolve_base(base: str) -> dict:
    """What ``--base`` names: ``none`` (random initialisation, for smoke tests), a registered model name
    (models_registry), a model file, or a Cellpose built-in name that Cellpose resolves itself."""
    from ..models_registry import resolve_registered

    b = str(base).strip()
    if b.lower() in ("none", "random", "scratch"):
        return {"name": "none", "kind": "random_init", "path": None, "sha256": None}
    reg = resolve_registered(b)
    if reg is not None:
        sha = _sha256(reg)
        from ..models_registry import check_model_pin, is_pack_model_name

        if is_pack_model_name(b):
            # a local entry under a pinned pack name (written before the form was reserved) must hold the pinned
            # weights, as on the run path (docs/DESIGN_system_packs.md, section 3.5); ModelPinError otherwise
            check_model_pin(b, sha)
        return {"name": b, "kind": "registered", "path": str(reg), "sha256": sha}
    p = Path(b)
    if p.is_file():
        return {"name": b, "kind": "file", "path": str(p.resolve()), "sha256": _sha256(p)}
    from ..models_registry import CELLPOSE_BUILTIN_NAMES, is_pack_model_name, known_models

    if is_pack_model_name(b):
        # a pinned pack model (cpsam_v2@0f1cc3f7): the verified blob of an installed pack's model spec, checked
        # against its pin (docs/DESIGN_system_packs.md, section 3.5); only reached where the lines below raise
        from ..models_registry import resolve_model_name

        path = resolve_model_name(b)
        return {"name": b, "kind": "pack", "path": str(path), "sha256": _sha256(path)}
    if b in known_models():
        raise FileNotFoundError(f"--base {b!r} is a registered model but its file is missing; fetch or register it first")
    if b.lower() not in CELLPOSE_BUILTIN_NAMES:
        # Cellpose fine-tunes from its default weights for any unknown name without saying so
        raise ValueError(f"--base {b!r} is not 'none', a registered model with its file on disk, a model file or a "
                         f"Cellpose built-in name {CELLPOSE_BUILTIN_NAMES}")
    return {"name": b, "kind": "cellpose_builtin", "path": None, "sha256": None}


def _build_net(base: dict, device: str):
    """The network ``train_seg`` trains: a fresh Cellpose-SAM network for ``none`` (no weights loaded, no
    download), else the net of ``CellposeModel(pretrained_model=...)`` with the explicit device."""
    import torch

    from ..device import torch_device

    dev = torch_device(device)
    if base["kind"] == "random_init":
        from cellpose.vit import CPSAM

        return CPSAM(dtype=torch.float32).to(dev)
    from cellpose import models

    kw = {"gpu": device != "cpu", "device": dev, "use_bfloat16": device == "cuda"}
    if base["kind"] == "cellpose_builtin" and base["name"].lower() in ("cpsam", "sam"):
        # the stock Cellpose-SAM model exactly as the segment stage loads it (nuclei._cellpose): Cellpose's
        # own default weights, whatever this Cellpose version calls them; the card records the file it used
        model = models.CellposeModel(**kw)
    else:
        model = models.CellposeModel(pretrained_model=base["path"] if base["path"] else base["name"], **kw)
    used = getattr(model, "pretrained_model", None)
    if base["path"] is None and isinstance(used, str) and Path(used).is_file():
        base["path"] = str(Path(used).resolve())
        base["sha256"] = _sha256(used)
    return model.net


def _train_seg(net, **kw):
    """``cellpose.train.train_seg`` (a seam the tests replace). Returns (model path, train losses, test losses)."""
    from cellpose import train

    return train.train_seg(net, **kw)


def _predict(model_path: str | Path, images: list[np.ndarray], device: str) -> list[np.ndarray]:
    """2D labels for each validation image with the trained model (a seam the tests replace)."""
    from cellpose import models

    from ..device import torch_device

    m = models.CellposeModel(gpu=device != "cpu", device=torch_device(device), pretrained_model=str(model_path),
                             use_bfloat16=False)
    out = []
    for img in images:
        masks = m.eval(np.asarray(img, dtype=np.float32), channel_axis=None)[0]
        out.append(np.asarray(masks).astype(np.int32))
    return out


def holdout_metrics(val: list[dict], gt: list[np.ndarray], pred: list[np.ndarray], iou_threshold: float = 0.5) -> dict:
    """segmentation_metrics on every validation pair, pooled overall, per orientation and per zone (tp, fp,
    fn summed; precision, recall, F1 and the miss rate fn / n_gt from the sums), plus the per-pair rows."""
    from ..validate import segmentation_metrics

    per_pair = []
    for p, g, pr in zip(val, gt, pred):
        m = segmentation_metrics(np.asarray(pr), np.asarray(g), iou_threshold=iou_threshold)
        per_pair.append({"image": Path(p["image"]).name, "image_id": p["image_id"], "orientation": p["orientation"],
                         "zone": p.get("zone", ""), **m})

    def pooled(rows):
        tp, fp, fn = (int(sum(r[k] for r in rows)) for k in ("tp", "fp", "fn"))
        prec = tp / (tp + fp) if tp + fp else 0.0
        rec = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
        ious = [float(r["mean_iou"]) for r in rows if np.isfinite(r["mean_iou"])]
        return {"n_pairs": len(rows), "n_gt": tp + fn, "n_pred": tp + fp, "tp": tp, "fp": fp, "fn": fn,
                "precision": prec, "recall": rec, "f1": f1, "miss_rate": fn / (tp + fn) if tp + fn else float("nan"),
                "mean_matched_iou": float(np.mean(ious)) if ious else float("nan")}

    return {"iou_threshold": float(iou_threshold), "overall": pooled(per_pair),
            "per_orientation": {o: pooled([r for r in per_pair if r["orientation"] == o])
                                for o in sorted({r["orientation"] for r in per_pair})},
            "per_zone": {z: pooled([r for r in per_pair if r["zone"] == z])
                         for z in sorted({r["zone"] for r in per_pair if r["zone"]})},
            "per_pair": per_pair}


def _versions() -> dict:
    from importlib import metadata

    out = {}
    for name in ("cellpose", "torch"):
        try:
            out[name] = metadata.version(name)
        except Exception:  # noqa: BLE001
            out[name] = "not-installed"
    return out


# ---- what the card says about the data (docs/DESIGN_system_packs.md, sections 3.6 and 10 step 5) ----------

CARD_DIMS = {"xy": "3d", "xz": "3d", "yz": "3d", PROJ: "2d"}    # the dims of a model trained on an orientation
_DATASET_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,63}")


def _filled(value) -> bool:
    """Whether a manifest cell holds a value (not None, NaN or empty text)."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    return not (isinstance(value, str) and not value.strip()) and value is not pd.NA


def _pair_dims(pair: dict, row: dict, sources: dict) -> str:
    """The dims of one pair: "2d" for a 2D plane or a tile of one (its manifest row holds a tile offset: a
    single-plane run's export, or imported 2D labels; without a manifest row, a run export.json calls planar) and
    for a proj pair, else "3d" (a slice of a 3D run, CARD_DIMS)."""
    if _filled(row.get("tile_y0")):
        return "2d"
    src = sources.get(str(pair["image_id"]))
    if not row and src is not None and src.get("geometry") == "planar":
        return "2d"
    return CARD_DIMS.get(str(pair.get("orientation")), "3d")


def _dataset_record(dataset: str) -> dict:
    """The card's ``dataset`` for public data (design 3.6): ``{"id", "licence", "split_sha256"}``. The id names
    the folder benchmarks/ext/<id>/ of the committed benchmark records, looked up under the config roots (the
    checkout, the working directory, $GERMQUANT_CONFIG_ROOT): the licence comes from its dataset.json and the
    split digest is the sha256 of its split.json on LF-normalised bytes. Without that folder both are None
    (a warning says so). ValueError for an id that is not a plain folder name."""
    from ..config import config_roots

    text = "" if dataset is None else str(dataset).strip()
    if not _DATASET_ID.fullmatch(text):
        raise ValueError(f"dataset {dataset!r} is not a dataset id: letters, digits, '_', '.' and '-', at most 64, "
                         f"starting with a letter or digit (the folder name under benchmarks/ext/)")
    rec = {"id": text, "licence": None, "split_sha256": None}
    for root in config_roots():
        folder = Path(root) / "benchmarks" / "ext" / text
        if not os.path.isdir(long_path(folder)):
            continue
        if os.path.isfile(long_path(folder / "dataset.json")):
            try:
                rec["licence"] = _read_json(folder / "dataset.json").get("licence")
            except (OSError, ValueError, AttributeError) as e:
                log.warning("%s: unreadable (%s); the card records no licence", folder / "dataset.json", e)
        if os.path.isfile(long_path(folder / "split.json")):
            with open(long_path(folder / "split.json"), "rb") as fh:
                rec["split_sha256"] = hashlib.sha256(fh.read().replace(b"\r\n", b"\n")).hexdigest()
        return rec
    log.warning("dataset %r: no committed benchmark record benchmarks/ext/%s/ under %s; the card records the id "
                "only", text, text, ", ".join(str(r) for r in config_roots()))
    return rec


def _source_runs(pairs_dir: Path) -> dict[str, dict | None]:
    """``{image_id: {"geometry", "pack"} or None}`` for the runs the folder's export.json lists (the last export
    into it): each run's record read again from its results folder; None when that folder or record cannot be
    read. Empty without an export.json (pairs made by hand, or merged from elsewhere)."""
    path = pairs_dir / "export.json"
    if not os.path.isfile(long_path(path)):
        return {}
    try:
        runs = _read_json(path).get("runs") or []
    except (OSError, ValueError, AttributeError) as e:
        log.warning("%s: unreadable (%s); the source runs of the pairs are unknown", path, e)
        return {}
    out: dict[str, dict | None] = {}
    for r in runs:
        if not isinstance(r, dict) or not r.get("image_id"):
            continue
        iid = str(r["image_id"])
        try:
            if not r.get("folder"):
                raise FileNotFoundError("no folder")
            rec = _run_record(r["folder"], iid)
        except (OSError, ValueError, TypeError):
            out[iid] = None
            continue
        pack = rec.get("pack")
        out[iid] = {"geometry": _geometry(rec)[0], "pack": dict(pack) if isinstance(pack, dict) else None}
    return out


def _manifest_rows(pairs_dir: Path) -> dict[str, dict]:
    """The manifest.csv rows of a pairs folder by the image path load_pairs gives (``str(pairs_dir / image)``)."""
    path = pairs_dir / "manifest.csv"
    if not os.path.isfile(long_path(path)):
        return {}
    df = pd.read_csv(long_path(path), low_memory=False)
    if "image" not in df.columns:
        return {}
    return {str(pairs_dir / str(r["image"])): r for r in df.to_dict("records")}


def _one_or_all(values: list, key: str) -> dict:
    """``{key: the one value}``, ``{key: None, <key>s: [sorted values]}`` when they differ, ``{key: None}`` when
    there is none (pixel_um_values and segmented_roles name several)."""
    uniq = sorted(set(values))
    if len(uniq) == 1:
        return {key: uniq[0]}
    extra = {"pixel_um": "pixel_um_values", "segmented_role": "segmented_roles"}[key]
    return {key: None, **({extra: uniq} if uniq else {})}


def _card_fields(pairs_dir: str | Path, pairs: list[dict], sources: dict | None = None, *, system: str | None = None,
                kind: str | None = None, dataset: str | None = None) -> dict:
    """What a model card says about its training data (design 3.6 and 10 step 5), in card order:

    * ``kind`` "cellpose"; ``dims`` "3d" (slices of 3D runs, or of imported 3D labels), "2d" (proj pairs, and
      the planes or tiles of single-plane runs and of imported 2D labels: _pair_dims) or "2d+3d";
    * ``pixel_um``: the pairs' manifest pixel size when they share one, else None (``pixel_um_values`` lists
      them when they differ);
    * ``object_kind``: ``kind`` when given (nucleus or cell), else the object kind that the packs of all the
      pairs' runs name (two recipes of one pack included), else nucleus;
    * ``system``: the text given, else None;
    * ``pack``: {id, version, recipe} of the pack recipe every source run used (export.json's runs, read again
      from their results folders), plus ``"edited": true`` when that is a project's edited recipe (the record
      project_config writes); None for runs without a pack, and None with ``pack_note`` when the runs differ
      (an edited recipe differs from its pack recipe) or are unknown, or when labels were imported (manifest
      source import:<tool>: they come from no run);
    * ``segmented_role``: the role the runs segmented (the manifest's channel_role), like pixel_um;
    * ``dataset`` (_dataset_record, or a record made before) only when ``dataset`` is given; ``source_url``
      None (a local model).

    ValueError for an unknown ``kind`` or a malformed ``dataset`` id, before anything is trained."""
    from ..packs.format import OBJECT_KINDS

    if kind is not None and str(kind) not in OBJECT_KINDS:
        raise ValueError(f"kind {kind!r} is not an object kind {OBJECT_KINDS} (what one label of the model is)")
    data = dict(dataset) if isinstance(dataset, dict) else (_dataset_record(dataset) if dataset is not None
                                                           else None)
    pairs_dir = Path(pairs_dir)
    sources = _source_runs(pairs_dir) if sources is None else sources
    rows = _manifest_rows(pairs_dir)
    dims = sorted({_pair_dims(p, rows.get(str(p["image"]), {}), sources) for p in pairs})
    pixel, roles = [], []
    for p in pairs:
        row = rows.get(str(p["image"]), {})
        try:
            v = float(row.get("pixel_um"))
        except (TypeError, ValueError):
            v = float("nan")
        if np.isfinite(v) and v > 0:
            pixel.append(round(v, 9))
        role = row.get("channel_role")
        if isinstance(role, str) and role.strip():
            roles.append(role.strip())
    ids = sorted({str(p["image_id"]) for p in pairs})
    # pairs of germquant import-labels (manifest source import:<tool>) come from no run: no pack recipe made them
    imported: dict[str, set] = {}
    for p in pairs:
        src = rows.get(str(p["image"]), {}).get("source")
        if isinstance(src, str) and src.startswith("import:"):
            imported.setdefault(str(p["image_id"]), set()).add(src)
    unknown = [i for i in ids if i not in imported and sources.get(i) is None]
    records = {i: sources[i]["pack"] for i in ids if i not in imported and sources.get(i) is not None}
    # a project run records edited (project_config, section 9.1): an edited recipe is not the pack recipe itself
    keys = {None if rec is None else (rec.get("id"), rec.get("version"), rec.get("recipe"), bool(rec.get("edited")))
            for rec in records.values()}
    pack, note = None, None
    notes = []
    if imported:
        notes.append(f"the labels of {len(imported)} of {len(ids)} image(s) were imported, not made by a run of a "
                     f"pack recipe ({', '.join(sorted(set().union(*imported.values())))}): "
                     f"{', '.join(sorted(imported)[:5])}" + (" ..." if len(imported) > 5 else ""))
    if unknown:
        notes.append(f"the source runs of {len(unknown)} of {len(ids)} image(s) are unknown (not in the pairs folder's "
                     f"export.json, or a results folder that cannot be read): {', '.join(unknown[:5])}"
                     + (" ..." if len(unknown) > 5 else ""))
    if notes:
        note = "; ".join(notes)
    elif len(keys) == 1 and None not in keys:
        pack_id, version, recipe, edited = next(iter(keys))
        pack = {"id": pack_id, "version": version, "recipe": recipe, **({"edited": True} if edited else {})}
    elif len(keys) > 1:
        note = "the pairs come from runs of different packs or recipes: " + "; ".join(
            sorted("no pack" if k is None else f"{k[0]}@{k[1]}/{k[2]}" + (" (edited)" if k[3] else "")
                   for k in keys))
    if kind is None:
        # the object kind every known source run's pack names (two recipes of one pack included), else nucleus
        found = {(rec or {}).get("object_kind") for rec in records.values()}
        one = next(iter(found)) if len(found) == 1 and not unknown and not imported else None
        kind = one if one in OBJECT_KINDS else "nucleus"
    out = {"kind": "cellpose", "dims": dims[0] if len(dims) == 1 else ("2d+3d" if dims else None),
           **_one_or_all(pixel, "pixel_um"), "object_kind": str(kind),
           "system": (str(system).strip() or None) if system is not None else None, "pack": pack}
    if note:
        out["pack_note"] = note
    out.update(_one_or_all(roles, "segmented_role"))
    if data is not None:
        out["dataset"] = data
    out["source_url"] = None
    return out


def _model_name(name, out_dir: Path) -> str:
    name = str(name or out_dir.name).strip()
    if not name or any(ch in name for ch in "/\\"):
        raise ValueError(f"model name {name!r} must be a bare name")
    return name


def train_nuclei(pairs_dir: str | Path, out_dir: str | Path, *, base: str = "cpsam", name: str | None = None,
                 epochs: int = 150, split_by: str = "series", val_frac: float = 0.2, seed: int = 0,
                 cpu: bool = False, learning_rate: float = 1e-5, weight_decay: float = 0.1, batch_size: int = 8,
                 min_train_masks: int = 5, system: str | None = None, kind: str | None = None,
                 dataset: str | None = None) -> dict:
    """Train a Cellpose nucleus model on the slice pairs in ``pairs_dir`` (module docstring) into ``out_dir``:
    the model file ``out_dir/<name>`` (``name`` defaults to the folder's name), ``training_log.json`` (losses
    per epoch, the split, counts, settings) and the model card ``out_dir/<name>.card.json``. ``base`` is a
    registered model name, a model file, a Cellpose built-in name (``cpsam``), a pinned pack model
    (``cpsam_v2@0f1cc3f7``) or ``none`` for a randomly initialised network (smoke tests only). ``cpu`` keeps
    torch off every accelerator; otherwise the device is germquant.device.select_device(). ``system``,
    ``kind`` (the card's object_kind) and ``dataset`` go on the card (_card_fields). Returns the paths, the
    split record, the hold-out metrics and the card."""
    from ..device import select_device

    pairs_dir, out_dir = Path(pairs_dir), Path(out_dir)
    pairs = load_pairs(pairs_dir)
    if not pairs:
        raise FileNotFoundError(f"no slice pairs (image/*.tif with labels/*.tif) under {pairs_dir}")
    name = _model_name(name, out_dir)
    train, val, split = split_pairs(pairs, by=split_by, val_frac=val_frac, seed=seed)
    fields = _card_fields(pairs_dir, pairs, system=system, kind=kind, dataset=dataset)
    device = "cpu" if cpu else select_device()
    return _train_model(pairs_dir, out_dir, name, pairs, train, val, split, base=base, device=device,
                        fields=fields, epochs=epochs, learning_rate=learning_rate, weight_decay=weight_decay,
                        batch_size=batch_size, min_train_masks=min_train_masks, val_frac=val_frac, seed=seed)


def _train_model(pairs_dir: Path, out_dir: Path, name: str, pairs: list[dict], train: list[dict],
                 val: list[dict], split: dict, *, base, device: str, fields: dict, epochs: int,
                 learning_rate: float, weight_decay: float, batch_size: int, min_train_masks: int,
                 val_frac: float, seed: int) -> dict:
    """One training of train_nuclei (and of every train_curve model) on the given train and val pairs:
    ``base`` is the --base text (resolved here, after the slices are read and checked) or a record that
    _resolve_base already made. ``fields`` (_card_fields) go on the card after its ``what``."""
    import tifffile

    def read(p, dtype):
        return np.asarray(tifffile.imread(long_path(p))).astype(dtype)

    tr_img = [read(p["image"], np.float32) for p in train]
    tr_lab = [read(p["labels"], np.int32) for p in train]
    va_img = [read(p["image"], np.float32) for p in val]
    va_lab = [read(p["labels"], np.int32) for p in val]
    n_masks = [int(np.count_nonzero(np.unique(lab))) for lab in tr_lab]
    thin = [Path(p["image"]).name for p, k in zip(train, n_masks) if k < int(min_train_masks)]
    if thin:
        log.warning("%d training slice(s) hold fewer than %d nuclei and will be dropped by Cellpose: %s",
                    len(thin), min_train_masks, ", ".join(thin[:8]) + (" ..." if len(thin) > 8 else ""))
    if len(thin) == len(train):
        raise ValueError(f"every training slice holds fewer than min_train_masks={min_train_masks} nuclei; "
                         f"lower --min-train-masks or export richer slices")
    base_info = base if isinstance(base, dict) else _resolve_base(base)
    os.makedirs(long_path(out_dir), exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    t0 = time.perf_counter()
    log.info("training %s from %s (%s) on %d slice(s), %d held out, %d epoch(s), device %s", name,
             base_info["name"], base_info["kind"], len(train), len(val), epochs, device)
    net = _build_net(base_info, device)
    model_path, train_losses, test_losses = _train_seg(
        net, train_data=tr_img, train_labels=tr_lab,
        test_data=va_img if val else None, test_labels=va_lab if val else None,
        n_epochs=int(epochs), learning_rate=float(learning_rate), weight_decay=float(weight_decay),
        batch_size=int(batch_size), save_path=str(out_dir), model_name=name, min_train_masks=int(min_train_masks))
    elapsed = time.perf_counter() - t0
    # cellpose writes <save_path>/models/<name>; the model file lives at <out_dir>/<name>
    src, dest = Path(str(model_path)), out_dir / name
    if src.resolve() != dest.resolve():
        os.replace(long_path(src), long_path(dest))
        try:
            os.rmdir(long_path(src.parent))
        except OSError:
            pass
    model_sha = _sha256(dest)

    metrics: dict = {"iou_threshold": 0.5, "overall": None, "per_orientation": {}, "per_zone": {}, "per_pair": []}
    if val:
        try:
            pred = _predict(dest, va_img, device)
            metrics = holdout_metrics(val, va_lab, pred, iou_threshold=0.5)
        except Exception as e:  # noqa: BLE001 - the model is trained and saved; the card says the metrics are missing
            log.warning("hold-out prediction failed (%s: %s); the card carries no metrics", type(e).__name__, e)
            metrics["error"] = f"{type(e).__name__}: {e}"
    else:
        log.warning("a single %s group: nothing held out, the card carries no hold-out metrics", split["by"])

    def pair_record(p, split_name):
        return {"image": os.path.relpath(p["image"], pairs_dir).replace(os.sep, "/"),
                "labels": os.path.relpath(p["labels"], pairs_dir).replace(os.sep, "/"),
                "image_sha256": _sha256(p["image"]), "labels_sha256": _sha256(p["labels"]),
                "image_id": p["image_id"], "series_id": p["series_id"], "orientation": p["orientation"],
                "slice_index": p["slice_index"], "zone": p.get("zone", ""), "split": split_name}

    training_set = [pair_record(p, "train") for p in train] + [pair_record(p, "val") for p in val]
    versions = _versions()
    settings = {"epochs": int(epochs), "learning_rate": float(learning_rate), "weight_decay": float(weight_decay),
                "batch_size": int(batch_size), "min_train_masks": int(min_train_masks), "device": device,
                "split_by": split["by"], "val_frac": float(val_frac), "seed": int(seed)}
    finished = datetime.now(timezone.utc).isoformat()
    log_rec = {
        "name": name, "model_file": str(dest), "model_sha256": model_sha, "base_model": base_info,
        "settings": settings, "split": split,
        "counts": {"n_pairs": len(pairs), "n_train": len(train), "n_val": len(val),
                   "per_orientation": dict(Counter(p["orientation"] for p in pairs)),
                   "train_per_orientation": dict(Counter(p["orientation"] for p in train)),
                   "val_per_orientation": dict(Counter(p["orientation"] for p in val)),
                   "n_train_below_min_masks": len(thin)},
        "train_files": [Path(p["image"]).name for p in train], "val_files": [Path(p["image"]).name for p in val],
        "losses": {"train": [float(v) for v in np.asarray(train_losses).ravel()],
                   "test": [float(v) for v in np.asarray(test_losses).ravel()]},
        "started": started, "finished": finished, "elapsed_s": round(elapsed, 1),
        "germquant_version": __version__, **{f"{k}_version": v for k, v in versions.items()},
    }
    with open(long_path(out_dir / "training_log.json"), "w", encoding="utf-8") as fh:
        json.dump(log_rec, fh, indent=1, default=str)
    card = {
        "name": name,
        "what": (f"Cellpose {fields.get('object_kind') or 'nucleus'} model trained with germquant train-nuclei "
                 f"from {base_info['name']} on "
                 f"{len(train)} slice pair(s) ({', '.join(sorted({p['orientation'] for p in train}))}) from "
                 f"{len(split['train_groups'])} {split['by']} group(s); {len(val)} pair(s) held out."),
        **fields,
        "model_file": str(dest), "model_sha256": model_sha, "model_size": dest.stat().st_size,
        "base_model": base_info, "date": finished, "germquant_version": __version__,
        **{f"{k}_version": v for k, v in versions.items()},
        "training": {**settings, "n_train": len(train), "n_val": len(val),
                     "train_groups": split["train_groups"], "val_groups": split["val_groups"]},
        "orientations": sorted({p["orientation"] for p in pairs}),
        "zones": dict(Counter(p.get("zone", "") or "untagged" for p in pairs)),
        "training_set": training_set,
        "holdout_metrics": {k: v for k, v in metrics.items() if k != "per_pair"},
        "holdout_per_pair": metrics.get("per_pair", []),
        "licence_note": LICENCE_NOTE,
    }
    card_path = out_dir / f"{name}.card.json"
    with open(long_path(card_path), "w", encoding="utf-8") as fh:
        json.dump(card, fh, indent=1, default=str)
    overall = metrics.get("overall") or {}
    log.info("%s trained in %.0f s -> %s (sha256 %s...); hold-out F1 at IoU 0.5: %s", name, elapsed, dest,
             model_sha[:12], f"{overall['f1']:.3f}" if overall else "n/a")
    return {"name": name, "model_path": str(dest), "model_sha256": model_sha, "card_path": str(card_path),
            "log_path": str(out_dir / "training_log.json"), "split": split, "metrics": metrics,
            "device": device, "base_model": base_info, "card": card}


# ---- learning curves (docs/TRAINING.md, step 6; design section 10, step 6; decision D22) -----------------

KEEP_MODELS = ("none", "best")
DEFAULT_SCORER = "holdout_metrics"       # the hold-out metrics train-nuclei computes; the score is overall F1
CURVE_TRAIN_KEYS = ("epochs", "split_by", "val_frac", "cpu", "learning_rate", "weight_decay", "batch_size",
                    "min_train_masks", "system", "kind", "dataset")
CURVE_COLS = ["size", "size_label", "repeat", "n_groups", "groups", "n_train", "n_val", "score", "f1", "precision",
              "recall", "miss_rate", "mean_matched_iou", "model_sha256", "model_size", "weights_kept", "model",
              "card", "elapsed_s", "status", "error"]


def parse_sizes(sizes) -> list:
    """The curve sizes, counted in training groups: whole numbers of at least 1 and "all" (every training
    group), from a list or a text with commas such as "1,2,4,8,16,all". Duplicates are dropped, the order is
    kept. ValueError naming the entry that is neither."""
    items = sizes.split(",") if isinstance(sizes, str) else list(sizes)
    out: list = []
    for s in items:
        if isinstance(s, str) and not s.strip():
            continue
        text = s.strip().lower() if isinstance(s, str) else None
        if text == "all":
            v = "all"
        elif text is not None and re.fullmatch(r"[0-9]+", text):
            v = int(text)
        elif text is None and isinstance(s, (int, np.integer)) and not isinstance(s, bool):
            v = int(s)
        else:
            raise ValueError(f"sizes: {s!r} is neither a whole number of training groups nor 'all'")
        if v != "all" and v < 1:
            raise ValueError(f"sizes: {s!r} is below 1 (a size counts training groups)")
        if v not in out:
            out.append(v)
    if not out:
        raise ValueError("sizes: no size given")
    return out


def _scorer_name(scorer) -> str:
    if scorer is None:
        return DEFAULT_SCORER
    name = getattr(scorer, "name", None)
    if isinstance(name, str) and name.strip():
        return name.strip()
    qual = getattr(scorer, "__qualname__", type(scorer).__qualname__)
    return f"{getattr(scorer, '__module__', '') or ''}.{qual}"


def _finite(v) -> float | None:
    """A float, or None for a missing or non-finite value (curve.csv and curve.json hold no NaN)."""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if np.isfinite(f) else None


def _write_curve(out_dir: Path, record: dict) -> None:
    pd.DataFrame(record["models"], columns=CURVE_COLS).to_csv(long_path(out_dir / "curve.csv"), index=False)
    with open(long_path(out_dir / "curve.json"), "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=1, default=str)


def _set_card_curve(card_path: Path, curve: dict | None = None, **changes) -> None:
    """Write the card's ``curve`` record: ``curve`` itself, or the card's own record with ``changes``."""
    card = _read_json(card_path)
    card["curve"] = dict(curve) if curve is not None else {**card["curve"], **changes}
    with open(long_path(card_path), "w", encoding="utf-8") as fh:
        json.dump(card, fh, indent=1, default=str)


def _drop_weights(path: Path) -> None:
    lp = long_path(path)
    if os.path.exists(lp):
        os.remove(lp)


def train_curve(pairs: str | Path, out: str | Path, *, sizes, repeats: int, seed: int, base: str,
                keep_models: str = "none", scorer=None, **train) -> dict:
    """A learning curve (design section 10, step 6; decision D22): train-nuclei on nested, seeded subsets of
    the training groups of the pairs folder ``pairs``, every model scored on the same held-out pairs.

    * The split is split_pairs(by=split_by, val_frac, seed), the split train-nuclei makes with that seed; its
      val groups are held out from every model, so a curve needs at least two groups.
    * The subsets: per repeat r, the training groups (series or images, the split's key) in the order of
      numpy.random.default_rng([seed, r]); the subset of size k is the first k of that order, so within a
      repeat every subset holds the smaller ones. ``sizes`` counts groups (parse_sizes; "all" is every
      training group); a size above the number of training groups is skipped and recorded.
    * Each model trains from ``base`` (resolved once) into ``out/models/n<kkk>_r<r>/`` as train-nuclei
      trains, and is scored on the val pairs right after: by default with the hold-out metrics train-nuclei
      computes (score = overall F1 at IoU 0.5, scorer "holdout_metrics"); a ``scorer`` callable,
      ``scorer(model_path, val_pairs, *, device) -> {"score": float, "metrics": dict}`` with the higher score
      the better, replaces it and its ``name`` attribute (else its qualified name) is recorded. P7's
      end-to-end runner plugs in there.
    * ``keep_models`` "none" (default) deletes each model's weights once it is scored; "best" keeps only the
      model with the highest score so far (ties keep the earlier, smaller one), so at most two weight files
      exist at a time. Every card (with a ``curve`` record), training_log.json, metric and digest stays.
    * ``curve.csv`` (CURVE_COLS, one row per model) and ``curve.json`` (sizes, subsets, val groups, scorer,
      settings, the rows and the best) are rewritten after every model.

    ``train`` passes train-nuclei's settings on (CURVE_TRAIN_KEYS: epochs, split_by, val_frac, cpu,
    learning_rate, weight_decay, batch_size, min_train_masks, system, kind, dataset), with train-nuclei's
    defaults. A subset whose every slice is too thin for Cellpose is recorded as failed, and the curve goes
    on. Returns curve.json's record plus ``out``, ``curve_csv`` and ``curve_json``."""
    from ..device import select_device

    unknown = sorted(set(train) - set(CURVE_TRAIN_KEYS))
    if unknown:
        raise TypeError(f"train_curve got unexpected keyword argument(s) {', '.join(unknown)}; it passes on "
                        f"{', '.join(CURVE_TRAIN_KEYS)}")
    keep = str(keep_models).strip().lower()
    if keep not in KEEP_MODELS:
        raise ValueError(f"keep_models must be one of {KEEP_MODELS}, got {keep_models!r}")
    asked = parse_sizes(sizes)
    if isinstance(repeats, bool) or not isinstance(repeats, (int, np.integer)) or int(repeats) < 1:
        raise ValueError(f"repeats must be a whole number of at least 1, got {repeats!r}")
    repeats, seed = int(repeats), int(seed)
    params = inspect.signature(train_nuclei).parameters
    opts = {**{k: params[k].default for k in CURVE_TRAIN_KEYS}, **train}
    pairs_dir, out_dir = Path(pairs), Path(out)
    all_pairs = load_pairs(pairs_dir)
    if not all_pairs:
        raise FileNotFoundError(f"no slice pairs (image/*.tif with labels/*.tif) under {pairs_dir}")
    pool_pairs, val, split = split_pairs(all_pairs, by=opts["split_by"], val_frac=opts["val_frac"], seed=seed)
    if not val:
        raise ValueError(f"a learning curve scores every model on pairs held out from all of them, but the pairs "
                         f"under {pairs_dir} form a single {split['by']} group: nothing can be held out")
    pool = list(split["train_groups"])
    n = len(pool)
    resolved = sorted({n if s == "all" else s for s in asked if s == "all" or s <= n})
    skipped = sorted(s for s in asked if s != "all" and s > n)
    if not resolved:
        raise ValueError(f"sizes {asked}: none fits the {n} training group(s) under {pairs_dir}")
    sources = _source_runs(pairs_dir)
    # a bad kind or dataset id is refused here, before anything trains; the dataset record is made once
    data = _card_fields(pairs_dir, all_pairs, sources, system=opts["system"], kind=opts["kind"],
                        dataset=opts["dataset"]).get("dataset")
    key = split["key"]
    subsets = {}
    for r in range(repeats):
        order = [pool[i] for i in np.random.default_rng([seed, r]).permutation(n)]
        subsets[f"r{r}"] = {str(k): sorted(order[:k]) for k in resolved}
    device = "cpu" if opts["cpu"] else select_device()
    base_info = _resolve_base(base)
    scorer_name = _scorer_name(scorer)
    os.makedirs(long_path(out_dir / "models"), exist_ok=True)
    record = {"format": 1, "pairs": str(pairs_dir), "sizes_requested": asked, "sizes": resolved,
              "skipped_sizes": skipped, "repeats": repeats, "seed": seed, "keep_models": keep,
              "scorer": scorer_name, "base_model": base_info, "device": device, "split": split,
              "val_groups": split["val_groups"], "train_pool": pool, "subsets": subsets,
              "settings": {k: opts[k] for k in ("epochs", "learning_rate", "weight_decay", "batch_size",
                                                "min_train_masks", "split_by", "val_frac")},
              "models": [], "best": None, "germquant_version": __version__,
              "started": datetime.now(timezone.utc).isoformat(), "finished": None}
    rows, best = record["models"], None
    log.info("learning curve: sizes %s of %d training %s group(s), %d repeat(s), %d pair(s) held out -> %s",
             resolved, n, split["by"], repeats, len(val), out_dir)
    for k in resolved:
        for r in range(repeats):
            groups = subsets[f"r{r}"][str(k)]
            chosen = set(groups)
            train_pairs = [p for p in pool_pairs if str(p[key]) in chosen]
            name = f"n{k:03d}_r{r}"
            row = {c: None for c in CURVE_COLS}
            row.update(size=k, size_label="all" if k == n else str(k), repeat=r, n_groups=len(groups),
                       groups=";".join(groups), n_train=len(train_pairs), n_val=len(val), weights_kept=False,
                       model=f"models/{name}/{name}", card=f"models/{name}/{name}.card.json")
            t0 = time.perf_counter()
            try:
                fields = _card_fields(pairs_dir, train_pairs + val, sources, system=opts["system"],
                                      kind=opts["kind"], dataset=data)
                res = _train_model(pairs_dir, out_dir / "models" / name, name, train_pairs + val, train_pairs, val,
                                   {**split, "train_groups": groups, "n_train": len(train_pairs)}, base=base_info,
                                   device=device, fields=fields, epochs=opts["epochs"],
                                   learning_rate=opts["learning_rate"], weight_decay=opts["weight_decay"],
                                   batch_size=opts["batch_size"], min_train_masks=opts["min_train_masks"],
                                   val_frac=opts["val_frac"], seed=seed)
            except ValueError as e:      # every slice of this subset too thin to train: recorded, the curve goes on
                log.warning("%s: not trained (%s)", name, e)
                row.update(status="failed", error=str(e), elapsed_s=round(time.perf_counter() - t0, 1))
                rows.append(row)
                _write_curve(out_dir, record)
                continue
            overall = res["metrics"].get("overall") or {}
            extra = None
            if scorer is None:
                score = _finite(overall.get("f1"))
            else:
                got = scorer(Path(res["model_path"]), list(val), device=device)
                if not isinstance(got, dict) or "score" not in got:
                    raise TypeError(f"the scorer {scorer_name} must return a dict with a 'score', got {got!r}")
                score, extra = _finite(got["score"]), got.get("metrics")
            better = score is not None and (best is None or score > best["score"])
            keep_this = keep == "best" and better
            if keep_this and best is not None and best["weights_kept"]:
                _drop_weights(out_dir / best["model"])
                best["weights_kept"] = False
                _set_card_curve(out_dir / best["card"], weights_kept=False)
            if not keep_this:
                _drop_weights(Path(res["model_path"]))
            row.update(score=score, f1=_finite(overall.get("f1")), precision=_finite(overall.get("precision")),
                       recall=_finite(overall.get("recall")), miss_rate=_finite(overall.get("miss_rate")),
                       mean_matched_iou=_finite(overall.get("mean_matched_iou")), model_sha256=res["model_sha256"],
                       model_size=int(res["card"]["model_size"]), weights_kept=keep_this,
                       elapsed_s=round(time.perf_counter() - t0, 1),
                       status="scored" if score is not None else "unscored", error=res["metrics"].get("error"))
            curve = {"size": k, "size_label": row["size_label"], "repeat": r, "groups": groups,
                     "scorer": scorer_name, "score": score, "weights_kept": keep_this, "curve": "curve.json"}
            if extra is not None:
                curve["metrics"] = extra
            _set_card_curve(Path(res["card_path"]), curve)
            if better:
                best = row
            rows.append(row)
            record["best"] = best
            _write_curve(out_dir, record)
            log.info("%s: score %s (%s)%s", name, "n/a" if score is None else f"{score:.3f}", scorer_name,
                     "; weights kept" if keep_this else "")
    record["finished"] = datetime.now(timezone.utc).isoformat()
    _write_curve(out_dir, record)
    return {**record, "out": str(out_dir), "curve_csv": str(out_dir / "curve.csv"),
            "curve_json": str(out_dir / "curve.json")}
