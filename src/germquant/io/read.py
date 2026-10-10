"""One entry point for every input format: ``read_image(path)`` -> Stack.

``.nd2`` goes to nd2_reader.read_stack exactly as before (the same function object: the validated numbers
depend on its bytes, and nothing in that path changed). Every other supported suffix (TIFF, OME-TIFF, CZI,
LIF, OME-Zarr; see io/formats.py) goes to bioio_reader, which needs the optional ``readers`` extra and
returns the same Stack. Dispatch is by file name only, so a missing or unreadable file fails inside the
reader with that reader's error.
"""
from __future__ import annotations

from pathlib import Path

from .bioio_reader import read_bioio, read_bioio_metadata
from .formats import image_format
from .nd2_reader import Stack, read_nd2_metadata, read_stack


def read_image(
    path: str | Path,
    xy_stride: int = 1,
    z_range: tuple[int, int] | None = None,
    as_float: bool = False,
    t: int | None = None,
) -> Stack:
    """Load any supported image as a (C, Z, Y, X) Stack with its voxel size in microns.

    ``xy_stride`` keeps every n-th pixel in y and x (and scales dy, dx by n); ``z_range`` keeps planes
    ``z0:z1``; ``as_float`` converts to float32; ``t`` reads one timepoint of a time series (ValueError when
    it is out of range; docs/TIMELAPSE.md). Raises ValueError for an unknown suffix, ImportError (naming the
    ``readers`` extra) when a non-nd2 file needs a bioio package that is not installed, and
    UnsupportedImageError for a multi-scene file, or a time-lapse read without ``t`` through bioio.
    """
    # t is passed on only when given, so the fixed-image call is exactly the call it has always been
    kw = {} if t is None else {"t": t}
    if image_format(path) == "nd2":
        return read_stack(path, xy_stride=xy_stride, z_range=z_range, as_float=as_float, **kw)
    return read_bioio(path, xy_stride=xy_stride, z_range=z_range, as_float=as_float, **kw)


def read_image_metadata(path: str | Path, t: int | None = None) -> dict:
    """Cheap header read (no pixels) for any supported image: ``sizes``, ``spacing`` (dz, dy, dx um),
    ``spacing_ok``, ``channel_names``, ``dtype``, ``is_2d``, ``n_timepoints``, ``period_s``, ``t_seconds``.
    The same function as read_nd2_metadata for an ``.nd2``; the bioio reader adds ``n_scenes`` and
    ``reader`` for the other formats. ``t`` only checks that the timepoint exists (ValueError otherwise)."""
    kw = {} if t is None else {"t": t}
    if image_format(path) == "nd2":
        return read_nd2_metadata(path, **kw)
    return read_bioio_metadata(path, **kw)
