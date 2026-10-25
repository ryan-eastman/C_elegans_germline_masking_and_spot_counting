"""Labels made elsewhere, written as training pairs: ``germquant import-labels`` (docs/TRAINING.md, step 3;
docs/DESIGN_system_packs.md sections 6.6 and 10 step 3; tranche 5 package P5, decisions D23 and D28 of the Train B
plan).

``import_labels(src, fmt, image, pairs_dir)`` reads a label image made in another program (the Cellpose GUI,
Cell-ACDC, ilastik, Labkit, micro-sam) and writes it, with the image it was drawn on, into a pairs folder in the
layout export-training writes: ``image/<id>_<orientation>_<index>.tif`` (uint16), ``labels/...`` (int32) and
manifest.csv rows whose ``source`` is ``import:<tool>``. The formats:

* ``tiff``: a TIFF or OME-TIFF label image, and ``cellacdc-npz``: a Cell-ACDC segm.npz (read with numpy and
  allow_pickle=False). Both are read by germquant.io.labels, the reader of segmentation.nuclei.labels_from, with
  its axes and frame rules (docs/LABELS.md): (Y, X), (Z, Y, X), with a T axis, ``frame`` picks the frame;
* ``png``: a PNG of one plane (grayscale or indexed colour: an index is the object's id), read with Pillow;
* ``cellpose-npy``: the masks of a Cellpose GUI ``_seg.npy``. Such a file is a Python pickle, which can run code
  when it is loaded, so it is read only with ``allow_pickle=True`` (``--allow-pickle``), with a warning.

The labels sit on the image's own grid at full resolution. A label image of one plane, for an image of one plane,
gives that plane as orientation xy, cut above 1024 pixels a side into tiles (the export of a single-plane run,
segment.training.plane_pairs); labels of every plane give the xy, xz and yz slices export-training draws for a 3D
run (3 per orientation among those that cut at least 3 objects, seeded per image, xz and yz stretched to the xy
pixel size: segment.training.slice_pairs). The ids are relabelled 1 to N in the order of the file's ids. When the
image is a time series, ``frame`` is also its timepoint, and the pairs carry that timepoint's image id
(``<stem>_t<NNN>``). ``channel`` names the image channel the labels were drawn on (an index, or the channel's
name); by default the only channel, or channel 0 with a warning.

Refused, before anything is written: a format the file name does not carry, a ``_seg.npy`` without
allow_pickle, labels of one plane for an image of several, labels on another grid, a label image without objects,
and an image that a committed benchmark split holds as test data: its file sha256 or its pixels_sha256
(germquant.io.pixels) among the test items of any benchmarks/ext/<dataset>/split.json, through that dataset's
converted.json, under the config roots (the checkout, the working directory, $GERMQUANT_CONFIG_ROOT). With no
committed split at all the import warns "no committed splits found" instead of passing silently (a container
holds no benchmarks/ folder).

import_labels imports germquant.pipeline, and with it torch, before it reads the image: a TIFF header read through
bioio loads dask, and dask loads CuPy, before torch otherwise (decision D28). Importing this module loads neither.
"""
from __future__ import annotations

import hashlib
import io
import json
import logging
import os
import re
from pathlib import Path

import numpy as np

from .fsutil import long_path

log = logging.getLogger(__name__)

FORMATS = ("tiff", "png", "cellacdc-npz", "cellpose-npy")
FORMAT_SUFFIXES = {"tiff": (".ome.tif", ".ome.tiff", ".tif", ".tiff"), "png": (".png",), "cellacdc-npz": (".npz",),
                   "cellpose-npy": (".npy",)}
DEFAULT_TOOLS = {"tiff": "tiff", "png": "png", "cellacdc-npz": "cellacdc", "cellpose-npy": "cellpose"}
SLICES_PER_IMAGE, MIN_LABELS, SEED = 3, 3, 0      # export-training's defaults, for the slices of 3D labels
PICKLE_WARNING = "a _seg.npy file is a Python pickle: loading it can run code; allow only files you made"
NO_SPLITS = "no committed splits found"
_TOOL = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,39}")
_PNG_LABEL_MODES = ("1", "L", "P", "I", "I;16", "I;16L", "I;16B", "I;16N", "F")   # one value per pixel


class ImportLabelsError(ValueError):
    """A labels file and an image that cannot become training pairs: a format the file name does not carry,
    labels that do not fit the image, a channel the image does not have, labels without objects."""


class SplitTestImageError(ImportLabelsError):
    """The image is test data of a committed benchmark split: test images are scored once and never trained on
    (docs/DESIGN_system_packs.md, section 6.6)."""


def import_labels(src, fmt, image, pairs_dir, *, allow_pickle=False, frame=0, tool=None, channel=None) -> dict:
    """Write the labels file ``src`` (format ``fmt``, one of FORMATS) and the image it was drawn on into the pairs
    folder ``pairs_dir`` (module docstring). ``tool`` names the program the labels came from (letters, digits, '_',
    '.' and '-'; default: tiff, png, cellacdc or cellpose by the format), recorded as the source ``import:<tool>``;
    ``channel`` the image channel (an index or a name); ``frame`` the frame of a T axis and a time series'
    timepoint. Rows of the same image id already in the folder's manifest are replaced; other rows stay.

    Returns ``{"image_id", "series_id", "t_index", "image", "image_sha256", "pixels_sha256", "labels",
    "labels_sha256", "format", "tool", "source", "dims" ("2d" or "3d"), "channel_index", "channel_name",
    "n_objects", "pairs" (the manifest rows written), "n_pairs", "replaced" (earlier rows of that image),
    "pairs_dir", "manifest", "splits_checked" (the benchmark datasets compared), "warnings"}``. Refusals raise
    ValueError (ImportLabelsError, SplitTestImageError, germquant.io.labels.LabelsError) or FileNotFoundError
    before anything is written (module docstring)."""
    src, image, pairs_dir = Path(src), Path(image), Path(pairs_dir)
    fmt = _check_format(fmt, src)
    tool = _check_tool(tool, fmt)
    if isinstance(frame, bool) or not isinstance(frame, (int, np.integer)) or int(frame) < 0:
        raise ValueError(f"frame must be a whole number of 0 or more, got {frame!r}")
    frame = int(frame)
    if fmt == "cellpose-npy" and not allow_pickle:
        from .io.labels import LabelsError

        raise LabelsError(f"{src.name}: a Cellpose _seg.npy is a Python pickle, which can run code when it is loaded, "
                          f"so it is read only with --allow-pickle (allow_pickle=True), for a file you made; or save "
                          f"its masks as a TIFF and import that (docs/LABELS.md)")
    if not os.path.isfile(long_path(src)):
        raise FileNotFoundError(f"labels file not found: {src}")
    if not os.path.exists(long_path(image)):
        raise FileNotFoundError(f"image not found: {image}")
    warnings: list[str] = []

    from . import pipeline        # torch first, before any image reader or backend (decision D28)
    from .io.formats import image_stem
    from .io.pixels import pixels_sha256

    meta = pipeline.read_input_metadata(image)
    n_t = _timepoints(meta)
    t = None
    if n_t > 1:
        if frame >= n_t:
            raise ImportLabelsError(f"--frame {frame} is out of range: {image.name} is a time series of {n_t} "
                                    f"timepoints (0 to {n_t - 1}), and the frame is also its timepoint")
        t = frame
    _meta, stack = pipeline._read_input(image, xy_stride=1, z_range=None, t=t)
    names = [str(n) for n in (stack.channel_names or meta.get("channel_names") or [])]
    c = _channel_index(channel, int(stack.n_channels), names, image)
    name = names[c] if c < len(names) else ""
    if channel is None and stack.n_channels > 1:
        _warn(warnings, f"{image.name} has {stack.n_channels} channels and none was named: the pairs take channel 0 "
                        f"({name or 'unnamed'}); name the channel the labels were drawn on with --channel")
    image_sha = _file_sha256(image)
    pix_sha = pixels_sha256(stack.data)
    checked = _check_splits(image, image_sha, pix_sha, warnings)

    with open(long_path(src), "rb") as fh:
        data = fh.read()
    nz = int(stack.data.shape[1])
    arr = _labels_array(src, fmt, data, nz=nz, frame=frame, warnings=warnings)
    from .io import labels as labels_io

    if arr.ndim == 2 and nz != 1:
        raise ImportLabelsError(f"{src.name} holds one plane (2D labels) and {image.name} has {nz} planes: the plane "
                                f"or projection the labels were drawn on is not known. Save that plane as an image of "
                                f"one plane and import the labels against it, or import labels of every plane")
    labels = labels_io._cut(arr, src.name, nz=nz, xy_stride=1, z_range=None)    # 1..N; one plane stays one
    img = np.asarray(stack.data[c])
    if labels.shape != img.shape:
        raise ImportLabelsError(f"{src.name}: its labels are {tuple(arr.shape)} and {image.name} is {nz} "
                                f"plane{'' if nz == 1 else 's'} of {tuple(img.shape[1:])} (rows, columns): a labels "
                                f"file holds the image's own grid at full resolution, one plane or every plane")
    n_objects = int(labels.max()) if labels.size else 0
    if n_objects == 0:
        raise ImportLabelsError(f"{src.name}: the labels hold no object (every pixel is 0), so there is nothing to train on")

    stem = image_stem(image).strip()
    image_id = stem if t is None else f"{stem}_t{t:03d}"
    dims = "2d" if arr.ndim == 2 else "3d"
    base = {"image_id": image_id, "series_id": stem, "t_index": t, "source_file": os.path.abspath(image),
            "xy_stride": 1, "z_range": "", "channel_role": "", "channel_index": int(c), "source": f"import:{tool}"}
    from .segment import training as TR

    if dims == "2d":
        yx = getattr(stack, "spacing_yx", None) or (tuple(stack.spacing[1:]) if stack.spacing_ok else None)
        if yx is None:
            _warn(warnings, f"{image.name} carries no pixel size: the pairs' pixel_um is left empty")
        pairs_gen = TR.plane_pairs({**base, "pixel_um": None if yx is None else float(yx[0])}, img[0], labels[0], {},
                                   pairs_dir)
    else:
        if not stack.spacing_ok:
            _warn(warnings, f"{image.name} carries no voxel size: xz and yz are stretched by the reader's 1 um guess, "
                            f"and the pairs' pixel_um is left empty")
        pairs_gen = TR.slice_pairs({**base, **({} if stack.spacing_ok else {"pixel_um": None})}, img, labels,
                                   stack.spacing, pairs_dir, slices_per_image=SLICES_PER_IMAGE,
                                   min_labels=MIN_LABELS, seed=SEED)
    replaced = _rows_of(pairs_dir / "manifest.csv", image_id)
    for sub in ("image", "labels"):
        os.makedirs(long_path(pairs_dir / sub), exist_ok=True)
    rows = [row for row, _nbytes in pairs_gen]
    merged = TR.merge_manifest(pairs_dir, rows)
    manifest = _write_manifest(merged, pairs_dir / "manifest.csv")
    log.info("%s: %d pair(s) imported from %s (%s, %s) -> %s", image_id, len(rows), src.name, dims,
             base["source"], manifest)
    return {"image_id": image_id, "series_id": stem, "t_index": t, "image": str(image), "image_sha256": image_sha,
            "pixels_sha256": pix_sha, "labels": str(src), "labels_sha256": hashlib.sha256(data).hexdigest(),
            "format": fmt, "tool": tool, "source": base["source"], "dims": dims, "channel_index": int(c),
            "channel_name": name, "n_objects": n_objects, "pairs": [_plain(r) for r in rows], "n_pairs": len(rows),
            "replaced": replaced, "pairs_dir": str(pairs_dir), "manifest": str(manifest), "splits_checked": checked,
            "warnings": warnings}


# ---- arguments -------------------------------------------------------------------------------------------

def _check_format(fmt, src: Path) -> str:
    text = str(fmt).strip().lower() if fmt is not None else ""
    if text not in FORMATS:
        raise ValueError(f"format {fmt!r} is not one of {', '.join(FORMATS)}")
    low = src.name.lower()
    if not any(low.endswith(s) and len(low) > len(s) for s in FORMAT_SUFFIXES[text]):
        raise ImportLabelsError(f"{src.name}: --format {text} reads {' or '.join(FORMAT_SUFFIXES[text])} files; name "
                                f"the format of this file (tiff, png, cellacdc-npz for a .npz, cellpose-npy for a "
                                f"_seg.npy)")
    return text


def _check_tool(tool, fmt: str) -> str:
    if tool is None:
        return DEFAULT_TOOLS[fmt]
    text = str(tool).strip()
    if not _TOOL.fullmatch(text):
        raise ValueError(f"tool {tool!r} is not a tool name: letters, digits, '_', '.' and '-', at most 40, starting "
                         f"with a letter or digit (it is recorded as the source import:<tool>)")
    return text


def _timepoints(meta) -> int:
    n = meta.get("n_timepoints") if isinstance(meta, dict) else None
    if n is None and isinstance(meta, dict) and isinstance(meta.get("sizes"), dict):
        n = meta["sizes"].get("T", 1)
    return max(1, int(n or 1))


def _channel_index(channel, n: int, names: list[str], image: Path) -> int:
    """The image channel the labels were drawn on: None (the only channel, else 0), an index, or a channel name
    (exact, else the one name that matches without regard to case)."""
    shown = ", ".join(f"{i} {nm or 'unnamed'}" for i, nm in enumerate(names[:n])) or f"{n} unnamed"
    if channel is None:
        return 0
    if isinstance(channel, bool):
        raise ImportLabelsError(f"channel {channel!r} is neither an index nor a channel name of {image.name}")
    if isinstance(channel, (int, np.integer)):
        if not 0 <= int(channel) < n:
            raise ImportLabelsError(f"channel {int(channel)} is out of range: {image.name} has {n} channel(s) "
                                    f"({shown})")
        return int(channel)
    text = str(channel).strip()
    if text in names[:n]:
        return names.index(text)
    folded = [i for i, nm in enumerate(names[:n]) if nm.lower() == text.lower()]
    if len(folded) == 1:
        return folded[0]
    raise ImportLabelsError(f"channel {text!r} is not a channel of {image.name} (channels {shown}); name one, or give "
                            f"its index")


def _warn(warnings: list[str], text: str) -> None:
    warnings.append(text)
    log.warning("%s", text)


# ---- reading the labels ----------------------------------------------------------------------------------

def _labels_array(src: Path, fmt: str, data: bytes, *, nz: int, frame: int, warnings: list[str]) -> np.ndarray:
    """The labels of the chosen frame as stored, (Y, X) or (Z, Y, X), before the relabel and the grid check."""
    from .io import labels as labels_io

    if fmt in ("tiff", "cellacdc-npz"):
        # the run path's reader (segmentation.nuclei.labels_from): its axes, T axis and frame rules
        return labels_io._frame(data, labels_io.labels_format(src), src.name, nz=nz, frame=frame)
    if frame != 0:
        raise labels_io.LabelsError(f"--frame {frame} names a frame of a T axis, but {src.name} ({fmt}) has none; "
                                    f"frame 0 is the file itself")
    if fmt == "cellpose-npy":
        _warn(warnings, f"{src.name}: {PICKLE_WARNING}")
        arr = np.asarray(labels_io.read_seg_npy(src, allow_pickle=True))
    else:
        arr = _png(data, src.name)
    if arr.ndim not in (2, 3):
        raise labels_io.LabelsError(f"{src.name}: a {arr.ndim}-D label image; labels are (Y, X) or (Z, Y, X)")
    return arr


def _png(data: bytes, name: str) -> np.ndarray:
    """A PNG label image of one plane, with Pillow: grayscale (8, 16 or 32 bits, or 1 bit) or indexed colour, whose
    indices are the ids. A colour or grey-and-alpha PNG, or an animated one, is refused."""
    from PIL import Image

    from .io.labels import LabelsError

    try:
        with Image.open(io.BytesIO(data)) as im:
            if getattr(im, "n_frames", 1) > 1:
                raise LabelsError(f"{name}: an animated PNG of {im.n_frames} frames; a label image has one plane")
            if im.mode not in _PNG_LABEL_MODES:
                raise LabelsError(f"{name}: a {im.mode} PNG holds more than one value per pixel (a colour image); a "
                                  f"label image holds one integer id per pixel (grayscale, or indexed colour)")
            arr = np.asarray(im).copy()
    except LabelsError:
        raise
    except Exception as e:  # noqa: BLE001 - any decode failure of the file is a refusal of the file
        raise LabelsError(f"{name}: not a readable PNG ({type(e).__name__}: {e})") from None
    return arr


# ---- the test-split check (design 6.6, D23) ----------------------------------------------------------------

def _file_sha256(path: Path) -> str | None:
    """sha256 of a file's bytes, None for a folder (an OME-Zarr image)."""
    if not os.path.isfile(long_path(path)):
        return None
    h = hashlib.sha256()
    with open(long_path(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_json(path: Path):
    with open(long_path(path), encoding="utf-8-sig") as fh:
        return json.load(fh)


def committed_splits() -> tuple[list[dict], list[str]]:
    """The committed benchmark splits under the config roots: ``[{"dataset", "split", "file": {sha256: item},
    "pixels": {sha256: item}, "blind": [test items without digests]}]`` from benchmarks/ext/<dataset>/split.json
    and converted.json (format 1, PLAN section 9.3), each dataset folder once; and the problems found (a split or
    converted.json that cannot be read)."""
    from . import config

    out, problems, seen = [], [], set()
    for root in config.config_roots():
        ext = Path(root) / "benchmarks" / "ext"
        if not os.path.isdir(long_path(ext)):
            continue
        for entry in sorted(os.listdir(long_path(ext))):
            folder = ext / entry
            split_path = folder / "split.json"
            if not os.path.isfile(long_path(split_path)):
                continue
            key = os.path.normcase(os.path.realpath(folder))
            if key in seen:
                continue
            seen.add(key)
            try:
                split = _read_json(split_path)
                if not isinstance(split, dict) or split.get("format") != 1 or not isinstance(split.get("test"), list):
                    raise ValueError("not a format 1 split with a test list")
                test = [str(i) for i in split["test"]]
            except (OSError, ValueError) as e:
                problems.append(f"benchmark split {entry}: {split_path} cannot be read ({e}); its test images are "
                                f"not checked")
                continue
            items: dict = {}
            conv_path = folder / "converted.json"
            if os.path.isfile(long_path(conv_path)):
                try:
                    conv = _read_json(conv_path)
                    if not isinstance(conv, dict) or conv.get("format") != 1 or not isinstance(conv.get("items"), dict):
                        raise ValueError("not a format 1 converted.json with items")
                    items = conv["items"]
                except (OSError, ValueError) as e:
                    problems.append(f"benchmark split {entry}: {conv_path} cannot be read ({e})")
            rec = {"dataset": entry, "split": str(split_path), "file": {}, "pixels": {}, "blind": []}
            for item in test:
                info = items.get(item) if isinstance(items.get(item), dict) else {}
                digests = [(k, info.get(f"{k}_sha256")) for k in ("file", "pixels")]
                for k, value in digests:
                    if isinstance(value, str) and value:
                        rec[k][value.lower()] = item
                if not any(isinstance(v, str) and v for _, v in digests):
                    rec["blind"].append(item)
            out.append(rec)
    return out, problems


def _check_splits(image: Path, image_sha: str | None, pix_sha: str, warnings: list[str]) -> list[str]:
    """Refuse an image that a committed split holds as test data (SplitTestImageError); else the datasets checked.
    Warns when there is no committed split at all, and for test items a split cannot recognise."""
    from . import config

    splits, problems = committed_splits()
    for text in problems:
        _warn(warnings, text)
    if not splits and not problems:
        roots = ", ".join(str(Path(r) / "benchmarks" / "ext") for r in config.config_roots())
        _warn(warnings, f"{NO_SPLITS}: no benchmarks/ext/<dataset>/split.json under {roots}; the test-split check "
                        f"(docs/DESIGN_system_packs.md, section 6.6) had nothing to compare {image.name} with")
        return []
    for s in splits:
        for what, digest, table in (("file sha256", image_sha, s["file"]), ("pixels_sha256", pix_sha, s["pixels"])):
            if digest and digest.lower() in table:
                raise SplitTestImageError(
                    f"{image.name} is test data of the committed benchmark split {s['dataset']} (item "
                    f"{table[digest.lower()]}, matched by its {what}; {s['split']}): test images are scored once and "
                    f"never trained on, so their labels are not imported (docs/DESIGN_system_packs.md, section 6.6)")
        if s["blind"]:
            _warn(warnings, f"benchmark split {s['dataset']}: {len(s['blind'])} test item(s) carry no file or pixels "
                            f"digest in converted.json ({', '.join(s['blind'][:5])}"
                            f"{' ...' if len(s['blind']) > 5 else ''}), so an image of theirs cannot be recognised")
    return sorted({s["dataset"] for s in splits})


# ---- writing ---------------------------------------------------------------------------------------------

def _rows_of(manifest: Path, image_id: str) -> int:
    """How many rows of ``image_id`` a pairs folder's manifest.csv holds now (0 without one)."""
    if not os.path.isfile(long_path(manifest)):
        return 0
    import pandas as pd

    old = pd.read_csv(long_path(manifest), low_memory=False)
    return int((old["image_id"].astype(str) == image_id).sum()) if "image_id" in old.columns else 0


def _write_manifest(df, path: Path) -> Path:
    """manifest.csv written as export-training writes it, through a temporary file next to it and os.replace,
    so a failure part way leaves the old manifest as it was."""
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        df.to_csv(long_path(tmp), index=False)
        os.replace(long_path(tmp), long_path(path))
    finally:
        if os.path.exists(long_path(tmp)):
            try:
                os.remove(long_path(tmp))
            except OSError:
                pass
    return path


def _plain(row: dict) -> dict:
    """A manifest row with plain JSON values (numpy numbers as Python numbers)."""
    out = {}
    for k, v in row.items():
        out[k] = v.item() if isinstance(v, np.generic) else v
    return out
