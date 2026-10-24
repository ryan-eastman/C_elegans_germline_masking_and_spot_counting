"""Labels from files: segment nothing, read the objects (docs/LABELS.md; docs/DESIGN_system_packs.md section 4.11,
tranche 5 package P2c).

``segmentation.nuclei.labels_from`` (code default null) is a path template with the fields {image_dir} (the
image's folder, absolute), {image_stem} (its file name without the format suffix) and {image_id} (the image id,
which for a fixed image is the stem). When it is set the segment stage reads the label image the template names
instead of segmenting: a TIFF or OME-TIFF label image, or a Cell-ACDC ``segm.npz``. A relative result resolves
against the recipe's folder (Config.base_dir), as io.channel_map does.

How a file is read, here and only here:

* TIFF and OME-TIFF through tifffile (never bioio), npz files with numpy and allow_pickle=False (the first array,
  as Cell-ACDC saves it). A Cellpose ``_seg.npy`` is a pickle, so the run path never reads a .npy file: `germquant
  import-labels --allow-pickle` converts one once (read_seg_npy is its reader).
* The axes: a TIFF's own axes (T, Z; a C or colour-sample axis of more than one is refused). A segm.npz has none,
  so its layout follows Cell-ACDC: (Y, X), (Z, Y, X) or (T, Z, Y, X). A 3D array without axis names, (A, Y, X), is
  a z stack when A is the image's number of planes and a series of 2D labels (T, Y, X) otherwise.
* ``segmentation.nuclei.labels_frame`` (code default 0, decision D13) picks the frame of a T axis; a file without
  one has frame 0 only.
* The labels sit on the image's own grid: one plane (2D), or every plane of the image. z_range and xy_stride cut
  them as the image readers cut the stack, a 2D file is copied through z, and the ids are relabelled
  consecutively, 1..N in the order of the file's ids. There is no size gate.

process_image calls labels_for_run once per image, before any stage: it finds the file, hashes its bytes (the
sha256 enters run_geometry labels_from {path, sha256} and the segment sub-hash) and reads the labels from those
same bytes. A missing file, a .npy, a file that does not fit the image, a time series (decision D15) and, in a
projected run, a file with more than one plane (D14) are refused there. Nothing here is reached when the key is
unset.
"""
from __future__ import annotations

import hashlib
import io
import os
from pathlib import Path

import numpy as np

from ..fsutil import long_path
from ..recipe_keys import labels_template_problem
from .formats import image_stem

TEMPLATE_KEY = "segmentation.nuclei.labels_from"
FRAME_KEY = "segmentation.nuclei.labels_frame"
TIFF_SUFFIXES = (".ome.tif", ".ome.tiff", ".tif", ".tiff")
LABEL_SUFFIXES = (*TIFF_SUFFIXES, ".npz")
SERIES_REFUSAL = ("segmentation.nuclei.labels_from does not work with a time series in this release: run the "
                  "image without --t (or t), or segment the series instead (docs/LABELS.md)")
SEG_NPY_HINT = ("the run path never reads a .npy file (a Cellpose _seg.npy is a pickle): convert it once with "
                "`germquant import-labels SRC --format cellpose-npy --allow-pickle --image IMG --pairs DIR`, or save "
                "the masks as a TIFF, and point segmentation.nuclei.labels_from at that file (docs/LABELS.md)")


class LabelsError(ValueError):
    """A labels file or setting that cannot run: a .npy or another format, a pickled npz, axes or a shape that do not
    fit the image, a frame that does not exist, a time series, more than one plane in a projected run. Fatal for
    the image, like a failed read."""


# ---- settings ----------------------------------------------------------------------------------------

def labels_template(cfg):
    """segmentation.nuclei.labels_from as the config holds it, None when unset: a labels run is one where this is
    not None (recipe_keys.check refuses a value that is not a template)."""
    return cfg.get(TEMPLATE_KEY)


def labels_frame(cfg) -> int:
    """segmentation.nuclei.labels_frame (code default 0): the frame of a T axis, read in labels runs only."""
    value = cfg.get(FRAME_KEY)
    return 0 if value is None else int(value)


def labels_path(cfg, image_path) -> Path | None:
    """The labels file of one image: the template formatted with its {image_dir}, {image_stem} and {image_id}, a
    relative result resolved against the recipe's folder. None when segmentation.nuclei.labels_from is unset.
    LabelsError for a value that is not a template. The file is not opened."""
    template = labels_template(cfg)
    if template is None:
        return None
    problem = labels_template_problem(TEMPLATE_KEY, template)
    if problem:
        raise LabelsError(problem)
    image = Path(image_path)
    stem = image_stem(image)
    path = Path(template.format(image_dir=os.path.dirname(os.path.abspath(image)), image_stem=stem,
                                image_id=stem.strip()))
    if not path.is_absolute():
        base = getattr(cfg, "base_dir", None)
        path = (Path(base) if base is not None else Path.cwd()) / path
    return path


def labels_format(path) -> str:
    """``"tiff"`` or ``"npz"`` from the file name (case-insensitive). LabelsError for a .npy file, with the
    import-labels hint, and for every other suffix."""
    name = Path(path).name
    low = name.lower()
    if low.endswith(".npy"):
        raise LabelsError(f"{name}: {SEG_NPY_HINT}")
    if any(low.endswith(s) and len(low) > len(s) for s in TIFF_SUFFIXES):
        return "tiff"
    if low.endswith(".npz") and len(low) > 4:
        return "npz"
    raise LabelsError(f"{name}: a labels file is a TIFF, an OME-TIFF or a Cell-ACDC segm.npz; supported suffixes: "
                      f"{', '.join(LABEL_SUFFIXES)}")


def labels_sha256(path) -> str:
    """sha256 of the file's bytes, read in chunks (FileNotFoundError when it is missing)."""
    h = hashlib.sha256()
    with open(long_path(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


# ---- reading -----------------------------------------------------------------------------------------

def read_labels(path, *, nz, xy_stride=1, z_range=None, frame=0) -> np.ndarray:
    """int32 (Z, Y, X) labels, 1..N, from a TIFF, an OME-TIFF or a segm.npz. ``nz`` is the image's number of planes
    (the file's own z extent, before z_range). The file holds one plane, copied through the planes z_range keeps,
    or every plane of the image, cut by z_range; xy_stride keeps every n-th pixel; ``frame`` picks the frame of a T
    axis. LabelsError for what cannot be a label image of that image (see the module docstring)."""
    path = Path(path)
    fmt = labels_format(path)
    with open(long_path(path), "rb") as fh:
        data = fh.read()
    arr = _frame(data, fmt, path.name, nz=nz, frame=frame)
    return _cut(arr, path.name, nz=nz, xy_stride=xy_stride, z_range=z_range)


def read_seg_npy(path, *, allow_pickle) -> np.ndarray:
    """The masks of a Cellpose GUI ``_seg.npy`` (a pickled dict), for `germquant import-labels` only: the run path
    never calls it. Without ``allow_pickle`` the file is refused unread. LabelsError when the file holds no
    masks."""
    path = Path(path)
    if not allow_pickle:
        raise LabelsError(f"{path.name}: a Cellpose _seg.npy is a pickle, read only with allow_pickle "
                          f"(`germquant import-labels --allow-pickle`), never on the run path")
    obj = np.load(long_path(path), allow_pickle=True)
    rec = obj.item() if isinstance(obj, np.ndarray) and obj.dtype == object and obj.shape == () else None
    masks = rec.get("masks") if isinstance(rec, dict) else None
    if masks is None:
        raise LabelsError(f"{path.name}: not a Cellpose _seg.npy (no masks entry in a saved dict)")
    return np.asarray(masks)


def image_planes(meta, stack, z_range) -> int:
    """The image file's number of planes for read_labels: the stack's, or with a z_range (the stack then holds the
    cut) the header's Z (1 when the header has no Z axis)."""
    n = int(stack.data.shape[1])
    if not z_range:
        return n
    sizes = meta.get("sizes") if isinstance(meta, dict) else None
    if not isinstance(sizes, dict):
        return n
    return int(sizes.get("Z", 1) or 1)


def labels_for_run(cfg, image_path, *, nz, shape_zyx, xy_stride=1, z_range=None, t=None,
                   projected=False) -> dict | None:
    """What a labels run needs, before any stage: None unless segmentation.nuclei.labels_from is set. Otherwise
    ``{"path", "sha256", "labels", "planes"}``: the file, the sha256 of its bytes, the labels read from those same
    bytes (int32 of ``shape_zyx``, the shape of the stack as read) and the file's planes in the chosen frame.

    Refused, before anything is read when it can be: a time series (``t``; decision D15), a .npy file or another
    format, a missing file (FileNotFoundError naming the file, the template and the image), labels that do not fit
    the image and, in a projected run, a file with more than one plane (D14)."""
    path = labels_path(cfg, image_path)
    if path is None:
        return None
    if t is not None:
        raise LabelsError(SERIES_REFUSAL)
    fmt = labels_format(path)
    try:
        with open(long_path(path), "rb") as fh:
            data = fh.read()
    except FileNotFoundError:
        raise FileNotFoundError(f"labels file not found: {path} (segmentation.nuclei.labels_from "
                                f"{labels_template(cfg)!r} for image {image_stem(Path(image_path)).strip()}); the "
                                f"image is not processed") from None
    arr = _frame(data, fmt, path.name, nz=nz, frame=labels_frame(cfg))
    if projected and arr.ndim == 3:
        raise LabelsError(f"{path.name} has {arr.shape[0]} planes: a projected run (segmentation.nuclei.project) "
                          f"takes a labels file with one plane, which is copied through z, and refuses one with more "
                          f"than one plane (save one plane of it, or switch the projection off)")
    labels = _cut(arr, path.name, nz=nz, xy_stride=xy_stride, z_range=z_range)
    if tuple(labels.shape) != tuple(shape_zyx):
        raise LabelsError(f"{path.name}: its labels are {tuple(arr.shape)} and the image read at xy_stride "
                          f"{int(xy_stride)} is {tuple(shape_zyx)} (planes, rows, columns); a labels file holds the "
                          f"image's own grid, one plane or every plane, before --z-range and --xy-stride")
    return {"path": path, "sha256": hashlib.sha256(data).hexdigest(), "labels": labels,
            "planes": 1 if arr.ndim == 2 else int(arr.shape[0])}


# ---- internals ---------------------------------------------------------------------------------------

def _frame(data: bytes, fmt: str, name: str, *, nz: int, frame: int) -> np.ndarray:
    """The labels of the chosen frame: (Y, X), or (Z, Y, X) with Z > 1."""
    arr, axes = _npz(data, name) if fmt == "npz" else _tiff(data, name)
    how = ""
    if axes == "QYX":                                   # no axis names: the image's planes decide
        n = int(arr.shape[0])
        axes = "ZYX" if n == int(nz) else "TYX"
        how = (f"its first axis of {n} was read as z, because the image has {int(nz)} planes" if axes == "ZYX" else
               f"its first axis of {n} was read as time, because the image has {int(nz)} "
               f"plane{'' if int(nz) == 1 else 's'}")
    frame = int(frame)
    if axes.startswith("T"):
        n = int(arr.shape[0])
        if not 0 <= frame < n:
            raise LabelsError(f"segmentation.nuclei.labels_frame {frame} is out of range: {name} has {n} "
                              f"frame{'' if n == 1 else 's'} (0 to {n - 1})" + (f"; {how}" if how else ""))
        arr = arr[frame]
    elif frame != 0:
        raise LabelsError(f"segmentation.nuclei.labels_frame {frame} names a frame of a T axis, but {name} has no T "
                          f"axis" + (f" ({how})" if how else "") + "; frame 0 is the file itself")
    if arr.ndim == 3 and arr.shape[0] == 1:
        arr = arr[0]
    return arr


def _cut(arr: np.ndarray, name: str, *, nz: int, xy_stride: int, z_range) -> np.ndarray:
    """z_range and xy_stride applied, ids relabelled 1..N (int32), a 2D file copied through the planes kept."""
    nz = int(nz)
    sl = slice(int(z_range[0]), int(z_range[1])) if z_range else slice(None)
    if arr.ndim == 3:
        if arr.shape[0] != nz:
            raise LabelsError(f"{name} has {arr.shape[0]} planes and the image has {nz}: a labels file holds one "
                              f"plane, or every plane of the image (before --z-range)")
        arr = arr[sl]
    step = max(1, int(xy_stride))
    lab = _relabel(arr[..., ::step, ::step], name)
    if lab.ndim == 2:
        lab = np.repeat(lab[None], len(range(nz)[sl]), axis=0)
    return np.ascontiguousarray(lab)


def _relabel(arr: np.ndarray, name: str) -> np.ndarray:
    """int32 labels 1..N, in the order of the file's ids (the consecutive relabel of the segmenters' size gates)."""
    from skimage.segmentation import relabel_sequential

    a = np.asarray(arr)
    if a.dtype == np.bool_:
        a = a.astype(np.uint8)
    if a.dtype.kind == "f":
        if not np.isfinite(a).all() or (a != np.floor(a)).any():
            raise LabelsError(f"{name}: label ids must be whole numbers, 0 for the background; the file holds "
                              f"fractions ({a.dtype})")
        a = a.astype(np.int64)
    elif a.dtype.kind not in "iu":
        raise LabelsError(f"{name}: holds {a.dtype} values; label ids are integers, 0 for the background")
    if a.size and a.dtype.kind == "i" and int(a.min()) < 0:
        raise LabelsError(f"{name}: label ids must be 0 (the background) or positive; the file holds {int(a.min())}")
    out, _fw, _inv = relabel_sequential(a)
    if out.size and int(out.max()) > np.iinfo(np.int32).max:
        raise LabelsError(f"{name}: holds more objects than an int32 label image can number")
    return out.astype(np.int32)


def _npz(data: bytes, name: str) -> tuple[np.ndarray, str]:
    """The first array of an npz (allow_pickle=False) and its Cell-ACDC axes."""
    import zipfile

    try:
        obj = np.load(io.BytesIO(data), allow_pickle=False)
    except ValueError as e:                  # a pickle, or a .npy of objects, under an .npz name
        raise LabelsError(f"{name}: not an npz archive of arrays that loads without a pickle (allow_pickle=False): "
                          f"{e}") from None
    except (OSError, EOFError, zipfile.BadZipFile) as e:
        raise LabelsError(f"{name}: not an npz archive ({e})") from None
    if not isinstance(obj, np.lib.npyio.NpzFile):
        raise LabelsError(f"{name}: not an npz archive (one array saved under an .npz name); save the labels with "
                          f"numpy.savez_compressed, as Cell-ACDC does")
    with obj:
        if not obj.files:
            raise LabelsError(f"{name}: holds no array")
        try:
            arr = np.asarray(obj[obj.files[0]])
        except ValueError:
            raise LabelsError(f"{name}: its array holds Python objects, which only a pickle can load; labels files "
                              f"are read with allow_pickle=False (save integer labels)") from None
    axes = {2: "YX", 3: "QYX", 4: "TZYX"}.get(arr.ndim)
    if axes is None:
        raise LabelsError(f"{name}: a {arr.ndim}-D array; a segm.npz holds (Y, X), (Z, Y, X) or (T, Z, Y, X) labels "
                          f"(axes)")
    return arr, axes


def _tiff(data: bytes, name: str) -> tuple[np.ndarray, str]:
    """The one image of a TIFF or OME-TIFF and its axes, reduced to T, Z, Q (z or time, decided by the caller), Y
    and X."""
    import tifffile

    try:
        with tifffile.TiffFile(io.BytesIO(data)) as tf:
            series = list(tf.series)
            if not series:
                raise LabelsError(f"{name}: holds no image")
            if len(series) > 1:
                raise LabelsError(f"{name}: holds {len(series)} images; a labels file holds one label image")
            arr = np.asarray(series[0].asarray())
            axes = str(series[0].axes).upper()
    except LabelsError:
        raise
    except (ValueError, OSError, EOFError, KeyError, IndexError) as e:     # TiffFileError is a ValueError
        raise LabelsError(f"{name}: not a readable TIFF ({type(e).__name__}: {e})") from None
    if len(axes) != arr.ndim or not axes.endswith(("YX", "YXS")):
        raise LabelsError(f"{name}: axes {axes!r} for a {arr.ndim}-D image; a label image has rows and columns (YX)")
    keep, kinds = [], []
    for i, a in enumerate(axes):
        size = int(arr.shape[i])
        if a in ("Y", "X", "T", "Z"):
            keep.append(i)
            kinds.append(a)
        elif a == "C" or (a == "S" and i > axes.index("Y")):
            if size > 1:
                what = f"{size} channels (axis C)" if a == "C" else f"a colour image ({size} samples per pixel)"
                raise LabelsError(f"{name}: {what}; a label image holds one integer id per pixel")
        elif size > 1:                       # Q, I and the like, or S before YX (tifffile's planes of 3 or 4)
            keep.append(i)
            kinds.append("Q")
    arr = arr.reshape(tuple(arr.shape[i] for i in keep)) if len(keep) != arr.ndim else arr
    return _resolve_axes(arr, "".join(kinds), name)


def _resolve_axes(arr: np.ndarray, axes: str, name: str) -> tuple[np.ndarray, str]:
    """Axes over T, Z, Y and X in that order, or "QYX" when one unnamed axis remains (z or time)."""
    rest = axes[:-2]
    if len(rest) > 2 or rest.count("T") > 1 or rest.count("Z") > 1 or (len(rest) == 2 and rest.count("Q") == 0
                                                                        and set(rest) != {"T", "Z"}):
        raise LabelsError(f"{name}: axes {axes!r}; a label image is (Y, X), (Z, Y, X), (T, Y, X) or (T, Z, Y, X)")
    if len(rest) == 2:
        if rest == "QQ":
            rest = "TZ"                      # Cell-ACDC's order
        elif "Q" in rest:
            rest = rest.replace("Q", "Z" if "T" in rest else "T")
        if rest == "ZT":
            arr = np.swapaxes(arr, 0, 1)
            rest = "TZ"
    return arr, rest + "YX"
