"""Inputs of ``germquant sc-straighten``: discovery, image ids, run geometry, the channel-subset reader, labels,
spot tables, refusals and the two guards (golden folders, the sealed SC validation split).

Image ids are read from the ``image_id`` column of ``<id>__nuclei.csv`` and every other input is looked up by
its exact name ``f"{image_id}__<suffix>"`` in the same folder, so ids that end in an underscore or hold spaces
(``20251105_N2_nohs_HERM _``) work. A refusal is an ``InputRefused`` with a code and a message; the command
skips that image and runs the others.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from ...fsutil import long_path

log = logging.getLogger(__name__)

NUCLEI_SUFFIX = "__nuclei.csv"


class InputRefused(Exception):
    def __init__(self, code: str, message: str):
        super().__init__(message)
        self.code = code
        self.message = message


@dataclass
class Found:
    image_id: str
    folder: Path
    nuclei_path: Path


@dataclass
class ImageInputs:
    image_id: str
    folder: Path
    nuclei: pd.DataFrame
    meta: dict
    spacing: tuple
    channel_map: dict
    file_path: str
    source_path: Path
    xy_stride: int
    z_range: list | None
    geometry_source: str
    labels_path: Path
    spot_paths: dict = field(default_factory=dict)
    image_summary_path: Path | None = None


# ------------------------------------------------------------------------------------ discovery
def _walk(root: Path, suffix: str) -> list[Path]:
    root_lp = long_path(root)
    hits = []
    for dirpath, _dirs, names in os.walk(root_lp):
        rel = os.path.relpath(dirpath, root_lp)
        base = root if rel == "." else root / rel
        for n in names:
            if n.endswith(suffix) and not n.startswith("batch_"):
                hits.append(base / n)
    return sorted(hits)


def discover(root: str | Path) -> list[Found]:
    """Every finished image under ``root``: one entry per ``<id>__nuclei.csv`` (id from the table)."""
    root = Path(root)
    out = []
    for f in _walk(root, NUCLEI_SUFFIX):
        try:
            head = pd.read_csv(long_path(f), usecols=["image_id"], nrows=1, dtype=str, keep_default_na=False)
            iid = str(head["image_id"].iloc[0]) if len(head) else f.name[: -len(NUCLEI_SUFFIX)]
        except Exception as e:  # noqa: BLE001 - an unreadable table is reported, not fatal
            log.warning("sc-straighten: cannot read %s (%s)", f, e)
            continue
        if f.name != iid + NUCLEI_SUFFIX:
            log.warning("sc-straighten: %s does not match its image_id %r; using the file name", f.name, iid)
            iid = f.name[: -len(NUCLEI_SUFFIX)]
        out.append(Found(image_id=iid, folder=f.parent, nuclei_path=f))
    return out


def select(found: list[Found], wanted: list[str] | None) -> list[Found]:
    """Keep the images named by image id or by result-folder name (all when ``wanted`` is empty)."""
    if not wanted:
        return list(found)
    w = set(wanted)
    return [f for f in found if f.image_id in w or f.folder.name in w]


# ------------------------------------------------------------------------------------ guards
def golden_root(path: str | Path) -> Path | None:
    """The folder holding a golden_manifest.json at ``path`` or above it, else None."""
    p = Path(path).resolve()
    for d in (p, *p.parents):
        if (d / "golden_manifest.json").is_file():
            return d
    return None


def stem_of(file_path: str) -> str:
    from ...io import image_stem

    base = str(file_path).replace("\\", "/").rsplit("/", 1)[-1]
    try:
        return image_stem(base)
    except Exception:  # noqa: BLE001 - an unknown suffix: strip the last extension
        return os.path.splitext(base)[0]


def _fold(s: str) -> str:
    return str(s).strip().casefold()


def sealed_role(image_id: str, file_path: str) -> str | None:
    """'test' or 'reserved' when the image is sealed by tests/golden/sc_split.json (repo checkouts only). Stems
    and ids are compared case-insensitively and without surrounding white space, so a renamed copy of a sealed
    image (``20260105_N2_HS_HERM_1``) is refused as well."""
    try:
        from ...benchmark import sc_split

        if sc_split.default_path() is None:
            return None
        split = sc_split.load_split()
    except Exception as e:  # noqa: BLE001 - a broken split record must not open the seal
        raise InputRefused("sealed_split_unreadable", f"tests/golden/sc_split.json could not be read: {e}") from e
    stems = {_fold(s): s for s in sc_split.sealed_stems(split)}
    for key in (stem_of(file_path), image_id):
        hit = stems.get(_fold(key))
        if hit is not None:
            return sc_split.role_of(hit, split)
    return None


def sealed_import(path: str) -> bool:
    """True when ``path`` is one of the registered live Filament truth files (test truth, check S3 of the
    protocol), by file name, case-insensitively. Silent (False) outside a repository checkout."""
    try:
        from ...benchmark import sc_split

        if sc_split.default_path() is None:
            return False
        ft = sc_split.load_split().get("filament_truth", {})
    except Exception:  # noqa: BLE001 - a broken split record must not open the seal
        return True
    names = {_fold(Path(str(p).replace("\\", "/")).name) for p in
             (ft.get("ims", {}).get("path", ""), ft.get("nas_copy", "")) if p}
    return _fold(Path(str(path).replace("\\", "/")).name) in names


# ------------------------------------------------------------------------------------ one image
def golden_manifest_entry(folder: str | Path, levels: int = 2) -> dict | None:
    """The entry for ``folder`` (by its name, the golden tag) in a ``golden_manifest.json`` one or two levels up
    (scripts/make_goldens.py records each tag's xy_stride and z_range there), else None."""
    folder = Path(folder)
    for d in list(folder.parents)[:levels]:
        f = d / "golden_manifest.json"
        if f.is_file():
            try:
                man = json.loads(f.read_text(encoding="utf-8"))
            except ValueError:
                return None
            e = man.get(folder.name) if isinstance(man, dict) else None
            return e if isinstance(e, dict) else None
    return None


def _source_mapped(file_path: str, source_map: list[tuple[str, str]]) -> str:
    fp = str(file_path)
    for old, new in source_map or []:
        a = fp.replace("\\", "/")
        o = old.replace("\\", "/")
        if a.startswith(o):
            rest = a[len(o):].lstrip("/")
            return new.rstrip("/\\") + ("/" + rest if rest else "")
    return fp


def resolve(found: Found, *, source_map=None, nuclei_mode: str = "germline", nucleus_ids=None,
            need_source: bool = True) -> ImageInputs:
    folder, iid = found.folder, found.image_id
    nuc = pd.read_csv(long_path(found.nuclei_path), low_memory=False)
    if not len(nuc):
        raise InputRefused("empty_nuclei", f"{iid}: nuclei.csv has no rows")
    if (folder / f"{iid}__series.json").exists() or any(p.is_dir() and p.name[:1] == "t" and p.name[1:].isdigit()
                                                         for p in folder.iterdir()):
        raise InputRefused("series_folder", f"{iid}: a time-lapse series folder; series are not supported in v1")
    first = nuc.iloc[0]
    meta = {c: first[c] for c in nuc.columns if c in _shared_meta()}
    # run geometry
    marker = folder / f"{iid}__done.json"
    geom_src = "done_marker"
    if marker.is_file():
        m = json.loads(marker.read_text(encoding="utf-8"))
    else:
        mf = folder / "run_manifest.json"
        m = json.loads(mf.read_text(encoding="utf-8")) if mf.is_file() else {}
        geom_src = "run_manifest"
        if "xy_stride" not in m:
            gm = golden_manifest_entry(folder)
            if gm is None or "xy_stride" not in gm:
                raise InputRefused("no_run_geometry", f"{iid}: no done marker, no xy_stride in run_manifest.json and "
                                                      "no entry for this folder in a parent golden_manifest.json")
            m, geom_src = gm, "golden_manifest"
    stride = int(m.get("xy_stride", 1))
    if stride != 1:
        raise InputRefused("undersampled_xy", f"{iid}: xy_stride {stride}; sc-straighten needs full-resolution runs")
    z_range = m.get("z_range") or None
    # channel roles as the run resolved them
    try:
        cmap = json.loads(str(first.get("channel_map_json", "{}")))
    except ValueError:
        cmap = {}
    if cmap.get("central_element") is None:
        raise InputRefused("no_central_element", f"{iid}: the run has no central_element channel")
    spacing = tuple(float(first[c]) for c in ("voxel_dz_um", "voxel_dy_um", "voxel_dx_um"))
    file_path = str(first.get("file_path", ""))
    src = Path(_source_mapped(file_path, source_map or []))
    if need_source and not Path(long_path(src)).exists():
        raise InputRefused("source_missing",
                           f"{iid}: the source image {src} is not there; map its folder with --source-map OLD=NEW")
    labels = folder / f"{iid}__nuclei_labels.tif"
    if not labels.is_file():
        raise InputRefused("no_labels", f"{iid}: {labels.name} is missing (the run did not write label images)")
    spots = {}
    for p in sorted(folder.glob("*.csv")):
        n = p.name
        if n == f"{iid}__spots.csv":
            spots["spots"] = p
        elif n.startswith(f"{iid}__spots_") and n.endswith(".csv"):
            spots[n[len(iid) + 2:-4]] = p
    summ = folder / f"{iid}__image_summary.csv"
    return ImageInputs(image_id=iid, folder=folder, nuclei=nuc, meta=meta, spacing=spacing, channel_map=cmap,
                       file_path=file_path, source_path=src, xy_stride=stride, z_range=z_range,
                       geometry_source=geom_src, labels_path=labels, spot_paths=spots,
                       image_summary_path=summ if summ.is_file() else None)


def _shared_meta() -> list[str]:
    from ...schema import SHARED_META

    return SHARED_META


def select_nuclei(nuc: pd.DataFrame, mode: str, nucleus_ids=None) -> list[int]:
    """The considered nuclei: ``germline`` (in_germline rows), ``pachytene`` (is_pachytene), ``all``; then
    restricted to ``nucleus_ids`` when given."""
    n = nuc.copy()
    if mode == "germline":
        col = n["in_germline"] if "in_germline" in n else pd.Series(True, index=n.index)
        n = n[col.astype(str).str.lower().isin(["true", "1"])]
    elif mode == "pachytene":
        if "is_pachytene" not in n or n["is_pachytene"].isna().all():
            raise InputRefused("no_pachytene", "--nuclei pachytene needs is_pachytene (run the staging stage first)")
        n = n[n["is_pachytene"].astype(str).str.lower().isin(["true", "1"])]
    elif mode != "all":
        raise ValueError(f"unknown nuclei selection {mode!r}")
    ids = sorted(int(v) for v in n["nucleus_id"])
    if nucleus_ids:
        keep = {int(v) for v in nucleus_ids}
        ids = [i for i in ids if i in keep]
    return ids


# ------------------------------------------------------------------------------------ pixels
def canonical_dask(darr, dims: list[str], *, xy_stride: int = 1, z_range=None):
    """The dask construction of germquant.io.read_stack: T and P loops at index 0, C and Z ensured, C Z Y X order,
    then z_range and the xy stride, all before any compute."""
    import dask.array as da

    sl, keep = [], []
    for d in dims:
        if d in ("T", "P"):
            sl.append(0)
        else:
            sl.append(slice(None))
            keep.append(d)
    darr = darr[tuple(sl)]
    dims = keep
    if "C" not in dims:
        darr = da.expand_dims(darr, 0)
        dims = ["C"] + dims
    if "Z" not in dims:
        ci = dims.index("C")
        darr = da.expand_dims(darr, ci + 1)
        dims = dims[: ci + 1] + ["Z"] + dims[ci + 1:]
    darr = darr.transpose([dims.index(t) for t in ["C", "Z", "Y", "X"]])
    if z_range is not None:
        darr = darr[:, z_range[0]: z_range[1]]
    if xy_stride > 1:
        darr = darr[:, :, ::xy_stride, ::xy_stride]
    return darr


class ND2Planes:
    """``planes(role, z)`` -> one (Y, X) plane of a role's channel of an .nd2, with the dask construction of
    germquant.io.read_stack; ``shape`` is (C, Z, Y, X) after z_range and the xy stride. Each call computes one
    frame, so only one frame is in memory at a time."""

    def __init__(self, path: str | Path, roles: dict, *, xy_stride: int = 1, z_range=None):
        import nd2

        with nd2.ND2File(str(path)) as f:
            sizes = dict(f.sizes)
            darr = f.to_dask()
        self.darr = canonical_dask(darr, list(sizes.keys()), xy_stride=xy_stride, z_range=z_range)
        self.roles = dict(roles)
        self.shape = tuple(int(v) for v in self.darr.shape)

    def __call__(self, role: str, z: int) -> np.ndarray:
        return np.asarray(self.darr[int(self.roles[role]), int(z)])


def nd2_role_arrays(path: str | Path, roles: dict, *, xy_stride: int = 1, z_range=None) -> dict:
    """The roles' channels of an .nd2 as (Z, Y, X) arrays, filled plane by plane: computing a whole channel at
    once keeps every frame's other channels alive until the end (the frames are the dask chunks), which on a
    106 x 2720 x 2720 stack is three channels, 4.7 GB, for one channel of 1.6 GB."""
    planes = ND2Planes(path, roles, xy_stride=xy_stride, z_range=z_range)
    out = {}
    for role in roles:
        first = planes(role, 0)
        arr = np.empty((planes.shape[1],) + first.shape, first.dtype)
        arr[0] = first
        for z in range(1, planes.shape[1]):
            arr[z] = planes(role, z)
        out[role] = arr
    return out


def nd2_role_planes(path: str | Path, roles: dict, *, xy_stride: int = 1, z_range=None) -> ND2Planes:
    """``plane(role, z)`` -> one (Y, X) plane of a role's channel of an .nd2 (one frame per call)."""
    return ND2Planes(path, roles, xy_stride=xy_stride, z_range=z_range)


def role_planes(path: str | Path, roles: dict, *, xy_stride: int = 1, z_range=None):
    """``plane(role, z)`` for any source: an .nd2 frame by frame; other formats are read once (all channels)."""
    if str(path).lower().endswith(".nd2"):
        return nd2_role_planes(path, roles, xy_stride=xy_stride, z_range=z_range)
    full = read_roles(path, roles, xy_stride=xy_stride, z_range=z_range)
    return lambda role, z: full[role][int(z)]


class CropBank:
    """Per-nucleus crops of one channel, keyed by the crop's slices and filled plane by plane, so the whole
    channel is never in memory; ``bank[slices]`` returns the crop (crops.make_crop indexes a role like an
    array). The crops are the same arrays slicing the whole channel would give."""

    def __init__(self, shape):
        self.shape = tuple(int(v) for v in shape)
        self.crops: dict[tuple, np.ndarray] = {}

    @staticmethod
    def key(sl) -> tuple:
        return tuple((int(s.start), int(s.stop)) for s in sl)

    def __getitem__(self, sl):
        return self.crops[self.key(sl)]

    @property
    def nbytes(self) -> int:
        return int(sum(a.nbytes for a in self.crops.values()))


def crop_banks(plane_fn, roles, shape, slices) -> dict:
    """{role: CropBank} holding the given crops; ``plane_fn(role, z)`` returns one (Y, X) plane of the role."""
    Z = int(shape[0])
    uniq: dict[tuple, tuple] = {}
    for sl in slices:
        uniq.setdefault(CropBank.key(sl), sl)
    by_z: list[list[tuple]] = [[] for _ in range(Z)]
    for sl in uniq.values():
        for z in range(sl[0].start, sl[0].stop):
            by_z[z].append(sl)
    banks = {}
    for r in roles:
        b = CropBank(shape)
        for z in range(Z):
            if not by_z[z]:
                continue
            plane = np.asarray(plane_fn(r, z))
            for sl in by_z[z]:
                k = CropBank.key(sl)
                if k not in b.crops:
                    b.crops[k] = np.empty(tuple(s.stop - s.start for s in sl), plane.dtype)
                b.crops[k][z - sl[0].start] = plane[sl[1], sl[2]]
        banks[r] = b
    return banks


def read_roles(path: str | Path, roles: dict, *, xy_stride: int = 1, z_range=None) -> dict:
    """{role: (Z, Y, X) array} for the given {role: channel index}."""
    p = str(path)
    if p.lower().endswith(".nd2"):
        return nd2_role_arrays(p, roles, xy_stride=xy_stride, z_range=z_range)
    from ...io import read_image

    log.warning("sc-straighten: %s is not an .nd2; the whole image is read once (all channels in memory)", p)
    st = read_image(p, xy_stride=xy_stride, z_range=tuple(z_range) if z_range else None)
    return {role: np.asarray(st.data[int(ci)]) for role, ci in roles.items()}


def read_labels(path: str | Path, max_label: int | None = None) -> np.ndarray:
    """The label TIFF page by page, into uint16 when every label fits (else int32)."""
    import tifffile

    with tifffile.TiffFile(long_path(path)) as tif:
        s = tif.series[0]
        shape = tuple(s.shape)
        if len(shape) == 2:
            return s.asarray().astype(np.uint16 if (max_label or 0) <= 65535 else np.int32)
        pages = list(s.pages)
        dtype = np.uint16 if (max_label is not None and max_label <= 65535) else np.int32
        out = np.zeros(shape, dtype)
        if len(pages) == shape[0]:
            for i, pg in enumerate(pages):
                a = pg.asarray()
                if dtype == np.uint16 and a.size and int(a.max()) > 65535:
                    out = out.astype(np.int32)
                    dtype = np.int32
                out[i] = a
        else:
            out = s.asarray()
            out = out.astype(np.uint16) if out.max() <= 65535 else out.astype(np.int32)
    return out


def read_spots(paths: dict) -> dict:
    out = {}
    for name, p in paths.items():
        try:
            df = pd.read_csv(long_path(p), low_memory=False)
        except Exception as e:  # noqa: BLE001
            log.warning("sc-straighten: cannot read %s (%s)", p, e)
            continue
        need = {"z_um", "y_um", "x_um", "nucleus_id"}
        if not need <= set(df.columns):
            continue
        cols = [c for c in ("spot_id", "nucleus_id", "marker", "z_um", "y_um", "x_um") if c in df.columns]
        out[name] = df[cols].copy()
    return out


def file_digest(path: str | Path, head_bytes: int | None = None) -> str:
    h = hashlib.sha256()
    with open(long_path(path), "rb") as fh:
        if head_bytes is None:
            for chunk in iter(lambda: fh.read(1 << 24), b""):
                h.update(chunk)
        else:
            h.update(fh.read(int(head_bytes)))
    return h.hexdigest()


def detector_max(inp: ImageInputs, ce: np.ndarray) -> float:
    """The detector maximum of the central_element channel: 2^bit_depth - 1 from the acquisition record when the
    run kept one, else the array's dtype maximum."""
    if inp.image_summary_path is not None:
        try:
            s = pd.read_csv(long_path(inp.image_summary_path), low_memory=False)
            col = "bit_depth_central_element"
            if col in s.columns and pd.notna(s[col].iloc[0]):
                return float(2 ** int(float(s[col].iloc[0])) - 1)
        except Exception:  # noqa: BLE001 - advisory
            pass
    if np.issubdtype(ce.dtype, np.integer):
        return float(np.iinfo(ce.dtype).max)
    return float(np.nanmax(ce)) if ce.size else 65535.0
