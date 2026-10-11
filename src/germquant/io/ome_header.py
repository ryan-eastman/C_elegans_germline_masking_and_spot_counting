"""Light TIFF header reads: the voxel size and channel names, with tifffile alone, no bioio.

The shallow pack check (``germquant pack check --image``, docs/DESIGN_system_packs.md, section 3.9) warns when
an image's pixel size or z step lies outside what a pack expects. Reading the header through bioio would
import dask.array, and dask imports CuPy whenever CuPy is installed: before torch, the cuBLAS load order that
rule R5 protects (torch must load before CuPy on the RTX 5090, or Cellpose silently falls back to the
watershed), and against gate G6 (a shallow check leaves torch, CuPy, SpotMAX and Cellpose out of
sys.modules). read_tiff_header reads the same values the bioio plugins report for a TIFF:

* a TIFF that carries OME-XML (any .tif or .tiff name; bioio tries its OME-TIFF reader first): PhysicalSizeZ,
  Y and X of the first image, taken as micrometres without unit conversion, as bioio-ome-tiff does, and the
  OME channel names;
* an ImageJ TIFF without OME-XML: the ImageJ ``spacing`` for z and the resolution tags for y and x, scaled by
  the ImageJ ``unit``, as bioio-tifffile does;
* any other TIFF: no z step, so ``spacing_ok`` is False, as bioio-tifffile reports no z size either.

A missing, non-finite or non-positive size gives (1, 1, 1) with ``spacing_ok`` False, the rule of both
readers. Nothing here imports more than tifffile (on first use) and the standard library. The pack golden
runner (germquant.packs_golden) reads its examples' OME-XML with ome_fields too.
"""
from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from pathlib import Path

from .formats import image_suffix

TIFF_SUFFIXES = (".ome.tif", ".ome.tiff", ".tif", ".tiff")

# micrometres per ImageJ unit, as bioio-tifffile reads them (an unknown unit counts as micrometres, as there).
# The micro sign, the Greek mu and ImageJ's escaped spelling (a backslash, u00B5, m) are built with chr() so
# this source stays ASCII.
_UNIT_UM = {
    "pm": 1e-6, "picometer": 1e-6, "nm": 1e-3, "nanometer": 1e-3,
    "micron": 1.0, "um": 1.0, chr(0xB5) + "m": 1.0, chr(0x3BC) + "m": 1.0, chr(0x5C) + "u00B5m": 1.0,
    "mm": 1e3, "millimeter": 1e3, "cm": 1e4, "centimeter": 1e4, "cal": 2.54e4,
}


def is_tiff(path: str | Path) -> bool:
    """Whether a file name is a TIFF or OME-TIFF (by its suffix, as io.formats reads it)."""
    return image_suffix(path) in TIFF_SUFFIXES


def _positive(value) -> float | None:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) and v > 0 else None


def ome_fields(xml_text: str | None) -> tuple[list[str], tuple[float | None, float | None, float | None]]:
    """(channel names, (PhysicalSizeZ, Y, X)) of the first image in OME-XML; namespaces ignored. A size that is
    missing, non-finite or not positive is None."""
    if not xml_text:
        return [], (None, None, None)
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return [], (None, None, None)
    pixels = next((el for el in root.iter() if el.tag.endswith("Pixels")), None)
    if pixels is None:
        return [], (None, None, None)
    names = [str(ch.get("Name", "")) for ch in pixels.iter() if ch.tag.endswith("}Channel") or ch.tag == "Channel"]
    return names, (_positive(pixels.get("PhysicalSizeZ")), _positive(pixels.get("PhysicalSizeY")),
                   _positive(pixels.get("PhysicalSizeX")))


def _imagej_sizes(tf) -> tuple[float | None, float | None, float | None]:
    """(z, y, x) in micrometres of an ImageJ TIFF, the rule of bioio-tifffile: z from the ImageJ ``spacing``,
    y and x from the resolution tags, all scaled by the ImageJ ``unit``; None where it cannot be read."""
    ij = tf.imagej_metadata or {}
    if "unit" not in ij:
        return None, None, None
    scale = _UNIT_UM.get(ij.get("unit"), 1.0)
    try:
        page = tf.series[0].pages[0] if tf.series else tf.pages[0]
        x_npix, x_units = page.tags["XResolution"].value
        y_npix, y_units = page.tags["YResolution"].value
        x = scale * x_units / x_npix
        y = scale * y_units / y_npix
    except (KeyError, IndexError, TypeError, ValueError, ZeroDivisionError, AttributeError):
        return None, None, None
    z = ij.get("spacing")
    return (_positive(z) * scale if _positive(z) is not None else None), _positive(y), _positive(x)


def read_tiff_header(path: str | Path) -> dict:
    """``{"spacing": (dz, dy, dx), "spacing_ok", "channel_names", "reader"}`` of a TIFF or OME-TIFF, read with
    tifffile (module docstring). Raises what tifffile raises for a file it cannot open."""
    import tifffile

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"image not found: {path}")
    names: list[str] = []
    with tifffile.TiffFile(path) as tf:
        root_ok = False
        if tf.ome_metadata:
            try:
                ET.fromstring(tf.ome_metadata)
                root_ok = True
            except ET.ParseError:
                root_ok = False
        if root_ok:
            names, (z, y, x) = ome_fields(tf.ome_metadata)
            reader = "tifffile (OME-XML)"
        elif tf.is_imagej:
            z, y, x = _imagej_sizes(tf)
            reader = "tifffile (ImageJ)"
        else:
            z = y = x = None
            reader = "tifffile"
    ok = all(v is not None for v in (z, y, x))
    return {"spacing": (float(z), float(y), float(x)) if ok else (1.0, 1.0, 1.0), "spacing_ok": ok,
            "channel_names": names, "reader": reader}
