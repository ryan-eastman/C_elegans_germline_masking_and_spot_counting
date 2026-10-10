"""The benchmark manifest, the pipeline output folder and the truth files.

Manifest: a JSON list of entries, or an object whose ``entries`` (or ``gonads``) key holds that list.
Each entry: ``tag``, ``split``, ``results_dir``, ``truth_dir``, ``nd2``, ``ims`` (see the package
docstring). Only ``tag`` and ``results_dir`` are required; a missing ``truth_dir`` means the scripts
must build the truth from the .ims themselves (their fallback path). An optional ``stem`` overrides
the truth file stem, which is otherwise the .ims file name without its extension.

Truth file contract, for an .ims named ``<stem>.ims`` (produced by the Imaris reader package):

* ``<stem>__imaris_spots.csv``: spot_id, x_um, y_um, z_um (OUR frame), xi, yi, zi (voxel indices,
  floor(imaris_coord / imaris_voxel)), radius_um, quality (NaN when absent), spot_class (string or
  empty), imaris_x_um, imaris_y_um, imaris_z_um (the Imaris values).
* ``<stem>__imaris_nuclei.csv``: surface_id, label (the mask value), x_um, y_um, z_um (our frame),
  xi, yi, zi, volume_um3, n_voxels_mask (NaN without a mask), imaris_x_um, imaris_y_um, imaris_z_um.
* ``<stem>__imaris_labels.tif``: int32 ZYX label image the size of the image, 0 where no traced
  surface. Absent when no mask channel could be identified.
* ``<stem>__imaris_params.json``: voxel sizes of both frames, ext_min, shape, the CreationParameters
  text of every Spots and Surfaces object, the ROI box if found, the mask channel and how it was
  identified, n_spots, n_surfaces, the spot classes present, the region_growing flag.

Frame conversion (the lesson from the prototype): Imaris reports voxel CENTRES scaled by its rounded
voxel size (0.108 um where the .nd2 says 0.108333), so matching in raw micrometres drifts by almost a
micrometre across a 2600 pixel field. The voxel index floor(imaris_coord / imaris_voxel) is exact, and
index times the .nd2 voxel size is the pipeline's own frame.
"""
from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

TRUTH_SUFFIXES = {
    "spots": "__imaris_spots.csv",
    "nuclei": "__imaris_nuclei.csv",
    "labels": "__imaris_labels.tif",
    "params": "__imaris_params.json",
}
TRUTH_SPOT_COLUMNS = ["spot_id", "x_um", "y_um", "z_um", "xi", "yi", "zi", "radius_um", "quality",
                      "spot_class", "imaris_x_um", "imaris_y_um", "imaris_z_um"]
TRUTH_NUCLEUS_COLUMNS = ["surface_id", "label", "x_um", "y_um", "z_um", "xi", "yi", "zi", "volume_um3",
                         "n_voxels_mask", "imaris_x_um", "imaris_y_um", "imaris_z_um"]
VOXEL_COLUMNS = ("voxel_dz_um", "voxel_dy_um", "voxel_dx_um")


# ------------------------------------------------------------------------------------ manifest
@dataclass
class Entry:
    tag: str
    results_dir: Path
    split: str = ""
    truth_dir: Path | None = None
    nd2: Path | None = None
    ims: Path | None = None
    stem: str | None = None
    extra: dict = field(default_factory=dict)

    @property
    def truth_stem(self) -> str:
        """The stem that names the truth files: ``stem`` when given, else the .ims file name without
        its extension, else the .nd2 stem, else the tag."""
        if self.stem:
            return self.stem
        if self.ims is not None:
            return Path(self.ims).stem
        if self.nd2 is not None:
            return Path(self.nd2).stem
        return self.tag


def _opt_path(v) -> Path | None:
    return Path(str(v)) if v not in (None, "") else None


def load_manifest(path: str | Path) -> list[Entry]:
    """Read a benchmark manifest (a list, or an object with ``entries``/``gonads``)."""
    with open(path, encoding="utf-8") as fh:
        doc = json.load(fh)
    if isinstance(doc, dict):
        items = doc.get("entries", doc.get("gonads"))
        if items is None:
            raise ValueError(f"{path}: a manifest object needs an 'entries' list")
    else:
        items = doc
    if not isinstance(items, list):
        raise ValueError(f"{path}: the manifest must be a list of entries")
    known = {"tag", "split", "results_dir", "truth_dir", "nd2", "ims", "stem"}
    entries = []
    for i, it in enumerate(items):
        if not isinstance(it, dict) or "tag" not in it or "results_dir" not in it:
            raise ValueError(f"{path}: entry {i} needs at least 'tag' and 'results_dir'")
        entries.append(Entry(
            tag=str(it["tag"]), results_dir=Path(str(it["results_dir"])),
            split=str(it.get("split") or ""), truth_dir=_opt_path(it.get("truth_dir")),
            nd2=_opt_path(it.get("nd2")), ims=_opt_path(it.get("ims")),
            stem=(str(it["stem"]) if it.get("stem") else None),
            extra={k: v for k, v in it.items() if k not in known}))
    return entries


# ----------------------------------------------------------------------------- pipeline output
@dataclass
class PipelineOutputs:
    image_id: str
    nuclei: pd.DataFrame
    spots: pd.DataFrame | None
    labels_path: Path | None
    voxel_zyx_um: tuple[float, float, float] | None
    voxel_source: str


def _single(results_dir: Path, pattern: str, required: bool) -> Path | None:
    hits = sorted(results_dir.glob(pattern))
    if not hits:
        if required:
            raise FileNotFoundError(f"no {pattern} in {results_dir}")
        return None
    if len(hits) > 1:
        raise ValueError(f"several {pattern} in {results_dir}: {[h.name for h in hits]}")
    return hits[0]


def _voxel_from_table(df: pd.DataFrame | None):
    if df is None or df.empty or not all(c in df.columns for c in VOXEL_COLUMNS):
        return None
    v = tuple(float(df[c].iloc[0]) for c in VOXEL_COLUMNS)
    return v if all(np.isfinite(s) and s > 0 for s in v) else None


def find_pipeline_outputs(results_dir: str | Path) -> PipelineOutputs:
    """Read ``<image>__nuclei.csv`` and ``<image>__spots.csv`` from a pipeline output folder and locate
    ``<image>__nuclei_labels.tif`` (not loaded here). The voxel size comes from the tables' own
    ``voxel_d{z,y,x}_um`` columns when present."""
    results_dir = Path(results_dir)
    nuc_path = _single(results_dir, "*__nuclei.csv", required=True)
    image_id = nuc_path.name[: -len("__nuclei.csv")]
    nuclei = pd.read_csv(nuc_path)
    spots_path = _single(results_dir, "*__spots.csv", required=False)
    spots = pd.read_csv(spots_path) if spots_path is not None else None
    labels_path = _single(results_dir, "*__nuclei_labels.tif", required=False)
    voxel = _voxel_from_table(nuclei) or _voxel_from_table(spots)
    return PipelineOutputs(image_id, nuclei, spots, labels_path, voxel,
                           "pipeline tables" if voxel else "")


def load_labels(path: str | Path) -> np.ndarray:
    """Read a ZYX label TIFF as the smallest integer dtype that holds its labels (a 2600 x 2600 x 71
    int32 image is 1.9 GB; as uint16 it is half that)."""
    import tifffile

    arr = np.asarray(tifffile.imread(str(path)))
    if arr.ndim == 2:
        arr = arr[None]
    if arr.ndim != 3:
        raise ValueError(f"{path}: expected a ZYX label image, got shape {arr.shape}")
    if arr.dtype.kind in "iu":
        lo, hi = int(arr.min()), int(arr.max())
        if lo >= 0 and hi <= np.iinfo(np.uint16).max:
            return arr.astype(np.uint16, copy=False)
        if lo >= 0 and hi <= np.iinfo(np.uint32).max:
            return arr.astype(np.uint32, copy=False)
        return arr.astype(np.int64, copy=False)
    return arr.astype(np.int64)


def our_voxel_zyx(outputs: PipelineOutputs | None, entry: Entry | None = None,
                  truth: Truth | None = None) -> tuple[tuple[float, float, float], str]:
    """The pipeline frame's voxel size (dz, dy, dx) in um and where it came from: the pipeline tables,
    else the .nd2 header, else the truth params JSON, else Imaris's own voxel size (with a warning,
    since that frame drifts against ours)."""
    if outputs is not None and outputs.voxel_zyx_um is not None:
        return outputs.voxel_zyx_um, outputs.voxel_source
    if entry is not None and entry.nd2 is not None and Path(entry.nd2).exists():
        try:
            from germquant.io.nd2_reader import read_nd2_metadata

            md = read_nd2_metadata(entry.nd2)
            if md.get("spacing_ok"):
                return tuple(float(s) for s in md["spacing"]), "nd2 header"
        except Exception as e:  # noqa: BLE001 - fall through to the next source
            log.warning("could not read the voxel size from %s: %s", entry.nd2, e)
    if truth is not None and truth.params:
        v = truth_voxel_from_params(truth.params, "ours")
        if v is not None:
            return v, "truth params (nd2 frame)"
        v = truth_voxel_from_params(truth.params, "imaris")
        if v is not None:
            log.warning("using Imaris's voxel size for our frame; coordinates may drift")
            return v, "truth params (Imaris frame)"
    raise ValueError("no voxel size: the tables lack voxel columns, no readable .nd2, no truth params")


# ---------------------------------------------------------------------------------------- truth
@dataclass
class Truth:
    stem: str
    spots: pd.DataFrame
    nuclei: pd.DataFrame
    labels_path: Path | None
    params: dict | None
    source: str = "truth files"

    @property
    def has_labels(self) -> bool:
        return self.labels_path is not None and Path(self.labels_path).exists()

    @property
    def has_spot_class(self) -> bool:
        return "spot_class" in self.spots.columns and bool(
            (self.spots["spot_class"].astype(str).str.strip() != "").any())


def truth_file_names(stem: str) -> dict:
    """The four truth file names for an .ims stem."""
    return {k: f"{stem}{suffix}" for k, suffix in TRUTH_SUFFIXES.items()}


def _coerce_truth_spots(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for c in TRUTH_SPOT_COLUMNS:
        if c not in df.columns:
            df[c] = "" if c == "spot_class" else np.nan
    for c in ("xi", "yi", "zi"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["spot_class"] = df["spot_class"].fillna("").astype(str).str.strip()
    df["quality"] = pd.to_numeric(df["quality"], errors="coerce")
    return df


def _coerce_truth_nuclei(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for c in TRUTH_NUCLEUS_COLUMNS:
        if c not in df.columns:
            df[c] = np.nan
    return df


def load_truth(truth_dir: str | Path | None, stem: str) -> Truth | None:
    """Read the truth files for ``stem`` from ``truth_dir``; None when the spots table is missing (the
    caller then falls back to the .ims). A missing nuclei table gives an empty frame; the label image
    and the params JSON are optional."""
    if truth_dir is None:
        return None
    truth_dir = Path(truth_dir)
    names = truth_file_names(stem)
    spots_path = truth_dir / names["spots"]
    if not spots_path.exists():
        return None
    spots = _coerce_truth_spots(pd.read_csv(spots_path, keep_default_na=True))
    nuc_path = truth_dir / names["nuclei"]
    nuclei = (_coerce_truth_nuclei(pd.read_csv(nuc_path)) if nuc_path.exists()
              else _coerce_truth_nuclei(pd.DataFrame(columns=TRUTH_NUCLEUS_COLUMNS)))
    labels_path = truth_dir / names["labels"]
    params = None
    params_path = truth_dir / names["params"]
    if params_path.exists():
        with open(params_path, encoding="utf-8") as fh:
            params = json.load(fh)
    return Truth(stem, spots, nuclei, labels_path if labels_path.exists() else None, params)


def _find_key(d: dict, candidates) -> object:
    low = {str(k).lower(): k for k in d}
    for c in candidates:
        if c.lower() in low:
            return d[low[c.lower()]]
    return None


def truth_voxel_from_params(params: dict, frame: str) -> tuple[float, float, float] | None:
    """Voxel size (dz, dy, dx) of one frame from the params JSON. Accepts a few spellings: for our
    frame ``voxel_zyx_um_ours`` / ``nd2_voxel_zyx_um`` / ``voxel_zyx_um_nd2`` (zyx order) or the
    xyz variants; for the Imaris frame ``voxel_xyz_um_imaris`` / ``imaris_voxel_xyz_um`` or zyx."""
    if not isinstance(params, dict):
        return None
    if frame == "ours":
        zyx = _find_key(params, ["voxel_zyx_um_ours", "nd2_voxel_zyx_um", "voxel_zyx_um_nd2",
                                 "our_voxel_zyx_um", "voxel_zyx_um"])
        xyz = _find_key(params, ["voxel_xyz_um_ours", "nd2_voxel_xyz_um", "voxel_xyz_um_nd2",
                                 "our_voxel_xyz_um"])
    else:
        zyx = _find_key(params, ["voxel_zyx_um_imaris", "imaris_voxel_zyx_um"])
        xyz = _find_key(params, ["voxel_xyz_um_imaris", "imaris_voxel_xyz_um", "voxel_xyz_um"])
    if zyx is not None:
        v = [float(s) for s in zyx]
    elif xyz is not None:
        v = [float(s) for s in xyz][::-1]
    else:
        return None
    if len(v) != 3 or not all(np.isfinite(s) and s > 0 for s in v):
        return None
    return (v[0], v[1], v[2])


_ROI_RE = re.compile(r'mMinXYZ="([-\d.\s]+)"\s+mMaxXYZ="([-\d.\s]+)"')


def _triplet(v):
    try:
        t = [int(round(float(s))) for s in (v.split() if isinstance(v, str) else v)]
    except (TypeError, ValueError):
        return None
    return t if len(t) == 3 else None


def _roi_from_text(obj) -> dict | None:
    """Search any string inside ``obj`` for Imaris's ``mMinXYZ="x y z" mMaxXYZ="x y z"`` ROI."""
    if isinstance(obj, str):
        m = _ROI_RE.search(obj)
        if m:
            lo, hi = _triplet(m.group(1)), _triplet(m.group(2))
            if lo and hi:
                return {"min_xyz": lo, "max_xyz": hi}
        return None
    if isinstance(obj, dict):
        for v in obj.values():
            r = _roi_from_text(v)
            if r:
                return r
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            r = _roi_from_text(v)
            if r:
                return r
    return None


def roi_box_from_params(params: dict | None) -> dict | None:
    """The Imaris ROI box in voxel indices as ``{"min_xyz": [x, y, z], "max_xyz": [x, y, z]}``
    (inclusive), or None. Reads ``roi_box`` / ``roi`` / ``region_of_interest`` written as a dict
    with min/max triplets or as a pair of triplets, and otherwise falls back to the
    ``mMinXYZ="..." mMaxXYZ="..."`` text of the Surfaces CreationParameters stored in the JSON."""
    if not isinstance(params, dict):
        return None
    box = _find_key(params, ["roi_box", "roi", "region_of_interest", "roi_box_voxels", "roi_voxels"])
    if isinstance(box, dict):
        lo = _find_key(box, ["min_xyz", "min_xyz_voxel", "min", "mMinXYZ", "min_voxel_xyz"])
        hi = _find_key(box, ["max_xyz", "max_xyz_voxel", "max", "mMaxXYZ", "max_voxel_xyz"])
        lo, hi = _triplet(lo) if lo is not None else None, _triplet(hi) if hi is not None else None
        if lo and hi:
            return {"min_xyz": lo, "max_xyz": hi}
    elif isinstance(box, (list, tuple)) and len(box) == 2:
        lo, hi = _triplet(box[0]), _triplet(box[1])
        if lo and hi:
            return {"min_xyz": lo, "max_xyz": hi}
    return _roi_from_text(params)


# ------------------------------------------------------------------------------ frame helpers
def imaris_um_to_voxel_index(coords_um, imaris_voxel_um) -> np.ndarray:
    """floor(coord / voxel) per axis; ``coords_um`` is (n, k) and ``imaris_voxel_um`` has k entries in
    the same axis order."""
    c = np.asarray(coords_um, float)
    v = np.asarray(imaris_voxel_um, float)
    return np.floor(c / v).astype(np.int64)


def voxel_index_to_um(indices, voxel_um) -> np.ndarray:
    """index times voxel size: the pipeline's own coordinate convention."""
    return np.asarray(indices, float) * np.asarray(voxel_um, float)


def points_to_voxel_index(coords_um, voxel_um) -> np.ndarray:
    """Our own coordinates back to voxel indices (round, since they were formed as index * voxel)."""
    return np.rint(np.asarray(coords_um, float) / np.asarray(voxel_um, float)).astype(np.int64)
