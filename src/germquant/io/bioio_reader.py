"""Read TIFF, OME-TIFF, CZI, LIF and OME-Zarr stacks through bioio into the Stack that read_stack returns.

bioio (github.com/bioio-devs/bioio) is a thin front end with one reader plugin per format. Every format
other than ``.nd2`` lands here and comes out as the same (C, Z, Y, X) Stack the pipeline expects: the
voxel size (dz, dy, dx) in microns as the file reports it (``spacing_ok`` False and a 1 um isotropic guess
when the file carries no usable voxel size, the rule the nd2 reader applies), the channel names as the
file reports them, and ``xy_stride`` / ``z_range`` applied on the lazy array before anything is loaded,
exactly as read_stack does.

The packages are optional and are imported only when a non-nd2 file is read: install the ``readers``
extra (``pip install "germquant[readers]"``), or bioio plus the plugin for the format at hand (the
``PLUGIN_FOR_SUFFIX`` table in io/formats.py names it). The ``.nd2`` path never touches this module.

Scope: one scene (position) per file, and one time point per read. A time-lapse (T > 1) is read one
timepoint at a time with ``t=<index>`` (selected on the lazy array, as the nd2 reader does; docs/TIMELAPSE.md);
without ``t`` it raises UnsupportedImageError, as does a multi-scene file, with a message saying what to do.
"""
from __future__ import annotations

import logging
import math
from pathlib import Path

import numpy as np

from .formats import PLUGIN_FOR_SUFFIX, image_suffix
from .nd2_reader import Stack, check_timepoint

log = logging.getLogger(__name__)

READERS_EXTRA = "readers"

_CORE_DIMS = "TCZYX"


class UnsupportedImageError(ValueError):
    """The file opens, but holds something the pipeline does not take yet: more than one time point,
    more than one scene or position, or a samples (RGB) dimension."""


def _plugin_for(path: Path) -> str:
    return PLUGIN_FOR_SUFFIX.get(image_suffix(path), "the bioio plugin for this format")


def _open(path: Path):
    """A ``bioio.BioImage`` for ``path``. ImportError naming the ``readers`` extra when bioio is not
    installed, and naming the format plugin when bioio has no reader for this file."""
    suffix = image_suffix(path)
    plugin = _plugin_for(path)
    try:
        from bioio import BioImage
    except ImportError as e:
        raise ImportError(
            f"reading {suffix} files needs the optional '{READERS_EXTRA}' extra: from the germquant folder, "
            f'uv pip install -e ".[{READERS_EXTRA}]" (or just bioio plus {plugin}); .nd2 files need nothing extra'
        ) from e
    try:
        from bioio_base.exceptions import UnsupportedFileFormatError
    except ImportError:  # pragma: no cover - bioio always ships bioio_base
        UnsupportedFileFormatError = Exception  # type: ignore[assignment,misc]
    try:
        return BioImage(str(path))
    except UnsupportedFileFormatError as e:
        raise ImportError(
            f"no bioio reader for {suffix} files is installed: uv pip install {plugin} "
            f"(the '{READERS_EXTRA}' extra installs every reader: uv pip install -e \".[{READERS_EXTRA}]\")"
        ) from e


def _n_timepoints(img) -> int:
    return int(getattr(img.dims, "T", 1) or 1)


def _check_single(img, path: Path, t: int | None = None) -> int | None:
    """Refuse what the pipeline cannot take yet, with a message that says what to do. Returns the timepoint
    to read: ``t`` checked against the T axis (ValueError when out of range), or None when no timepoint was
    asked for, which is only allowed for a file with a single time point."""
    dims = img.dims
    n_t = _n_timepoints(img)
    if t is not None:
        t = check_timepoint(t, n_t, path)
    elif n_t > 1:
        raise UnsupportedImageError(
            f"{path.name} is a time-lapse ({n_t} time points); pass t=<index> (0 to {n_t - 1}) to read one "
            f"time point (docs/TIMELAPSE.md)")
    scenes = tuple(getattr(img, "scenes", ()) or ())
    if len(scenes) > 1:
        shown = ", ".join(str(s) for s in scenes[:5]) + (", ..." if len(scenes) > 5 else "")
        raise UnsupportedImageError(
            f"{path.name} holds {len(scenes)} scenes or positions ({shown}); germquant reads a single "
            f"scene, so export or split the file into one scene per image")
    for d in dims.order:
        n = int(getattr(dims, d, 1) or 1)
        if d not in _CORE_DIMS and n > 1:
            raise UnsupportedImageError(
                f"{path.name} has a {d} dimension of size {n} (samples / RGB?); germquant reads one "
                f"grey-level value per channel")
    return t


def _spacing(img) -> tuple[tuple[float, float, float], bool]:
    """(dz, dy, dx) in microns as the file reports them, or a (1, 1, 1) guess with ``spacing_ok``
    False when any of the three is missing, non-finite or not positive (the nd2 reader's rule)."""
    guess = (1.0, 1.0, 1.0)
    try:
        pps = img.physical_pixel_sizes
        raw = (pps.Z, pps.Y, pps.X)
    except Exception:  # noqa: BLE001 - a plugin failing on its metadata must not block the read; spacing_ok says so
        return guess, False
    out = []
    for s in raw:
        try:
            v = float(s)
        except (TypeError, ValueError):
            return guess, False
        if not math.isfinite(v) or v <= 0:
            return guess, False
        out.append(v)
    return (out[0], out[1], out[2]), True


def _channel_names(img) -> list[str]:
    try:
        return [str(c) for c in (img.channel_names or [])]
    except Exception:  # noqa: BLE001 - no names: the channel map falls back to indices and flags it
        return []


def _subsample(arr, xy_stride: int, z_range):
    """``z_range`` planes and every ``xy_stride``-th pixel in y and x, as read_stack slices them."""
    if z_range is not None:
        arr = arr[:, z_range[0] : z_range[1]]
    if xy_stride > 1:
        arr = arr[:, :, ::xy_stride, ::xy_stride]
    return arr


def _load_czyx(img, path: Path, xy_stride: int, z_range, t: int | None = None) -> np.ndarray:
    """The (C, Z, Y, X) pixels of the single scene and one time point (``t`` when given, else index 0; every
    other dimension at index 0), sub-sampled on the lazy dask array before compute, as read_stack does. The
    TIFF plugins stream through tifffile's zarr store, which needs the zarr package; when that import fails
    the whole stack is read at once and sub-sampled in memory instead (the same bytes, more memory for a
    moment)."""
    extra = {d: 0 for d in img.dims.order if d not in "CZYX"}
    if t is not None and "T" in extra:
        extra["T"] = int(t)
    try:
        return np.ascontiguousarray(_subsample(img.get_image_dask_data("CZYX", **extra), xy_stride, z_range))
    except ImportError as e:
        log.info("%s: lazy read unavailable (%s); reading the whole stack", path.name, e)
    return np.ascontiguousarray(_subsample(np.asarray(img.get_image_data("CZYX", **extra)), xy_stride, z_range))


def read_bioio(
    path: str | Path,
    *,
    xy_stride: int = 1,
    z_range: tuple[int, int] | None = None,
    as_float: bool = False,
    t: int | None = None,
) -> Stack:
    """Load a non-nd2 stack as (C, Z, Y, X) through bioio; same contract as nd2_reader.read_stack, ``t``
    included (one timepoint of a time series, selected on the lazy array)."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"image not found: {path}")
    img = _open(path)
    t = _check_single(img, path, t)
    (dz, dy, dx), spacing_ok = _spacing(img)
    names = _channel_names(img)

    data = _load_czyx(img, path, xy_stride, z_range, t)
    if xy_stride > 1:
        dy *= xy_stride
        dx *= xy_stride
    if as_float:
        data = data.astype(np.float32)

    return Stack(
        data=data,
        spacing=(dz, dy, dx),
        channel_names=names,
        path=path,
        downsample_xy=xy_stride,
        spacing_ok=spacing_ok,
        t_index=t,
    )


def read_bioio_metadata(path: str | Path, *, t: int | None = None) -> dict:
    """Header read, no pixels: the keys read_nd2_metadata returns (sizes, spacing, spacing_ok,
    channel_names, dtype, is_2d, n_timepoints, period_s, t_seconds) plus ``n_scenes`` and ``reader`` (the
    bioio plugin that opened it). ``period_s`` and ``t_seconds`` are None here: the interval and the time
    stamps are not read from these formats yet. ``t`` only checks that the timepoint exists."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"image not found: {path}")
    img = _open(path)
    sp, ok = _spacing(img)
    dims = img.dims
    sizes = {d: int(getattr(dims, d, 1) or 1) for d in dims.order}
    n_t = _n_timepoints(img)
    if t is not None:
        check_timepoint(t, n_t, path)
    try:
        dtype = str(img.dtype)
    except (AttributeError, TypeError, ValueError):
        dtype = "unknown"
    reader = type(getattr(img, "reader", img)).__module__.split(".")[0]
    return {
        "sizes": sizes,
        "spacing": sp,
        "spacing_ok": ok,
        "channel_names": _channel_names(img),
        "dtype": dtype,
        "is_2d": sizes.get("Z", 1) == 1,
        "n_scenes": len(tuple(getattr(img, "scenes", ()) or ())),
        "reader": reader,
        "n_timepoints": n_t,
        "period_s": None,
        "t_seconds": None,
    }
