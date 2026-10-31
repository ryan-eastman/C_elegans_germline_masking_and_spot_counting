"""Light image header reads for projects (docs/PROJECTS.md, "The light entry"; PLAN.md decision D28; rule R5).

``read_header(path)`` returns ``{"sizes", "n_images", "channel_names", "spacing", "spacing_ok", "spacing_yx_um",
"single_plane", "reader"}`` without loading pixels, and without bioio for the two formats a project meets most:

* an ``.nd2`` through io.read_image_metadata, that is read_nd2_metadata: the nd2 package reads the header and
  loads no dask, CuPy or torch. ``n_images`` is the number of positions (the P axis);
* a TIFF or OME-TIFF with tifffile alone. A file that carries OME-XML (any .tif or .tiff name: bioio tries its
  OME-TIFF reader first) gives the sizes, PhysicalSizeZ, Y and X and the channel names of its first OME image,
  the channel name falling back to the channel ID, as bioio-ome-tiff reads them, and ``n_images`` counts its
  OME images (bioio's scenes). Any other TIFF gives what bioio-tifffile reports: the axes of its first series
  (an unknown axis guessed from the shape, as bioio does), the channel names Channel:0:<i>, the ImageJ spacing
  rule of io.ome_header, and one image per series.

``spacing`` and ``spacing_ok`` follow the readers' rule: (1, 1, 1) and False when any of dz, dy and dx is
missing. ``spacing_yx_um`` is (dy, dx) whenever the file states both as lengths, even without a z step (a
calibrated single-plane image, where io.ome_header and bioio give the (1, 1, 1) guess), by the rule of the run's
reader (io.bioio_reader._spacing_yx), which a single-plane run takes its pixel size from: an ImageJ TIFF without a
unit or in pixels states no length (None), and a TIFF that is neither OME nor ImageJ states its pixel size in its
resolution tags, scaled by its ResolutionUnit as bioio-tifffile scales them (None for the unit none, an aspect
ratio). ``single_plane`` comes from the sizes: no Z, or a Z of 1 (decision D11).

Any other format (CZI, LIF, OME-Zarr) is read only after ``import germquant.pipeline``: bioio imports dask.array,
and dask imports CuPy whenever CuPy is installed, before torch, which on the RTX 5090 makes a later Cellpose run in
the same process fail or silently fall back to the watershed (the cuBLAS load order of rule R5).
"""
from __future__ import annotations

import math
import os
import xml.etree.ElementTree as ET
from pathlib import Path

from ..fsutil import long_path
from ..io.formats import image_format
from ..io.ome_header import is_tiff

_CORE = "TCZYX"
_UNKNOWN_AXES = "QI"          # tifffile's unknown axis and image sequence, which bioio-tifffile guesses from the shape
# micrometres per ResolutionUnit (TIFF tag 296: 2 inch, 3 centimetre, 4 millimetre, 5 micrometre), as bioio-tifffile
# scales the resolution tags of a TIFF that is neither OME nor ImageJ; 1 (none) is an aspect ratio, not a length
_RESOLUTION_UNIT_UM = {2: 25400.0, 3: 10000.0, 4: 1000.0, 5: 1.0}
_NO_LENGTH = ("", "pixel", "pixels")   # ImageJ units that are not lengths (io.bioio_reader._tiff_without_length_unit)


def _nd2_metadata(path) -> dict:
    """The .nd2 header (the seam the tests replace with a stand-in of any sizes)."""
    from ..io import read_image_metadata

    return read_image_metadata(path)


def _other_metadata(path) -> dict:
    """The header of a format only bioio reads, after the pipeline (torch first)."""
    import germquant.pipeline  # noqa: F401  (side effect: torch claims cuBLAS before dask can load CuPy)

    from ..io import read_image_metadata

    return read_image_metadata(path)


def _positive(value) -> float | None:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) and v > 0 else None


def _record(sizes: dict, n_images: int, names: list, z, y, x, reader: str) -> dict:
    ok = all(v is not None for v in (z, y, x))
    return {"sizes": sizes, "n_images": int(n_images), "channel_names": [str(n) for n in names],
            "spacing": (float(z), float(y), float(x)) if ok else (1.0, 1.0, 1.0), "spacing_ok": ok,
            "spacing_yx_um": (float(y), float(x)) if y is not None and x is not None else None,
            "single_plane": int(sizes.get("Z", 1) or 1) == 1, "reader": reader}


def _from_reader_metadata(meta: dict, *, n_images: int, reader: str) -> dict:
    sizes = {str(k): int(v) for k, v in dict(meta.get("sizes") or {}).items()}
    ok = bool(meta.get("spacing_ok"))
    z, y, x = (float(v) for v in meta.get("spacing", (1.0, 1.0, 1.0)))
    rec = _record(sizes, n_images, list(meta.get("channel_names") or []), z if ok else None, y if ok else None,
                  x if ok else None, reader)
    if "spacing_yx_um" in meta:          # the bioio reader's own xy calibration of a one-plane file (D11), which
        yx = meta["spacing_yx_um"]       # a single-plane run takes even when the file states no z step
        rec["spacing_yx_um"] = None if yx is None else (float(yx[0]), float(yx[1]))
    return rec


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _ome_header(root: ET.Element) -> dict | None:
    images = [el for el in root if _local(el.tag) == "Image"]
    pixels = next((el for el in (images[0] if images else root).iter() if _local(el.tag) == "Pixels"), None)
    if pixels is None:
        return None

    def size(axis: str) -> int:
        try:
            return max(1, int(pixels.get(f"Size{axis}", 1)))
        except (TypeError, ValueError):
            return 1

    sizes = {d: size(d) for d in _CORE}
    names = []
    for ch in pixels:
        if _local(ch.tag) == "Channel":
            name = ch.get("Name")
            names.append(name if name is not None else ch.get("ID", ""))
    return _record(sizes, max(1, len(images)), names, _positive(pixels.get("PhysicalSizeZ")),
                   _positive(pixels.get("PhysicalSizeY")), _positive(pixels.get("PhysicalSizeX")),
                   "tifffile (OME-XML)")


def _merge_axes(axes: str, shape: tuple) -> str:
    """bioio-tifffile's guess for unknown axes: each unknown axis takes the next axis of the shape-based guess
    (TCZYX trimmed to the number of dimensions) that is not already named."""
    guess = _CORE[max(0, len(_CORE) - len(shape)):]
    out: list[str] = []
    for a in axes:
        if a not in _UNKNOWN_AXES:
            out.append(a)
            continue
        pick = next((g for g in guess if g not in out and g not in axes), a)
        out.append(pick)
    return "".join(out)


def _series_header(tf) -> dict:
    from ..io.ome_header import _imagej_sizes

    s0 = tf.series[0]
    axes, shape = str(s0.axes), tuple(int(n) for n in s0.shape)
    if any(a in _UNKNOWN_AXES for a in axes):
        axes = _merge_axes(axes, shape)
    raw = dict(zip(axes, shape))
    sizes = {d: int(raw.get(d, 1)) for d in _CORE}
    sizes.update({d: n for d, n in raw.items() if d not in sizes})
    names = [f"Channel:0:{i}" for i in range(sizes["C"])]
    length = True                                 # whether the xy calibration is a length (the run reader's rule)
    if tf.is_imagej:
        z, y, x = _imagej_sizes(tf)
        unit = (tf.imagej_metadata or {}).get("unit")
        length = unit is not None and str(unit).strip().lower() not in _NO_LENGTH
        reader = "tifffile (ImageJ)"
    else:
        z = None
        y, x = _resolution_yx(tf)
        reader = "tifffile"
    rec = _record(sizes, len(tf.series), names, z, y, x, reader)
    if not length:
        rec["spacing_yx_um"] = None               # pixels, not micrometres: a single-plane run assumes 1 um and flags it
    return rec


def _resolution_yx(tf) -> tuple[float | None, float | None]:
    """(dy, dx) in micrometres from the resolution tags of the first page of a TIFF that is neither OME nor ImageJ,
    scaled by its ResolutionUnit as bioio-tifffile scales them; (None, None) for the unit none (1, an aspect ratio)
    and for a missing tag, where the run reader gives no xy calibration either (io.bioio_reader._spacing_yx)."""
    try:
        page = tf.series[0].pages[0] if tf.series else tf.pages[0]
        tags = page.tags
        unit = int(tags["ResolutionUnit"].value)
        x_npix, x_units = tags["XResolution"].value
        y_npix, y_units = tags["YResolution"].value
        if unit == 1:
            return None, None
        scale = _RESOLUTION_UNIT_UM.get(unit, 1.0)
        return _positive(scale * y_units / y_npix), _positive(scale * x_units / x_npix)
    except (KeyError, IndexError, TypeError, ValueError, ZeroDivisionError, AttributeError):
        return None, None


def _tiff_header(path: Path) -> dict:
    import tifffile

    if not os.path.exists(long_path(path)):           # the long form, as it is opened (a deep NAS tree)
        raise FileNotFoundError(f"image not found: {path}")
    with tifffile.TiffFile(long_path(path)) as tf:
        xml = tf.ome_metadata
        if xml:
            try:
                rec = _ome_header(ET.fromstring(xml))
            except ET.ParseError:
                rec = None
            if rec is not None:
                return rec
        return _series_header(tf)


def read_header(path) -> dict:
    """The header of one image (module docstring). Raises what the reader raises for a file it cannot open,
    FileNotFoundError for a missing file and ValueError for an unsupported suffix."""
    path = Path(path)
    fmt = image_format(path)
    if fmt == "nd2":
        meta = _nd2_metadata(path)
        sizes = {str(k): int(v) for k, v in dict(meta.get("sizes") or {}).items()}
        return _from_reader_metadata(meta, n_images=max(1, int(sizes.get("P", 1) or 1)), reader="nd2")
    if is_tiff(path):
        return _tiff_header(path)
    meta = _other_metadata(path)
    return _from_reader_metadata(meta, n_images=max(1, int(meta.get("n_scenes", 1) or 1)),
                                 reader=str(meta.get("reader", "bioio")))
