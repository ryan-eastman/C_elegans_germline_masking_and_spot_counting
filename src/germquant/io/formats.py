"""Which file names the pipeline reads, and how a file name maps to an image id.

Pure name helpers (no reader is imported here), shared by the dispatcher (io/read.py), the bioio reader
and the sample-metadata parser. The Nikon ``.nd2`` is the validated input; every other suffix is read
through bioio and needs the optional ``readers`` extra (see io/bioio_reader.py).
"""
from __future__ import annotations

from pathlib import Path

ND2_SUFFIX = ".nd2"

# file suffix -> the bioio plugin (distribution name) that reads it; the OME double suffixes come first
# in the lookups below because ".ome.tif" must not be read as a plain ".tif"
PLUGIN_FOR_SUFFIX: dict[str, str] = {
    ".ome.tif": "bioio-ome-tiff", ".ome.tiff": "bioio-ome-tiff",
    ".tif": "bioio-tifffile", ".tiff": "bioio-tifffile",
    ".czi": "bioio-czi",
    ".lif": "bioio-lif",
    ".ome.zarr": "bioio-ome-zarr", ".zarr": "bioio-ome-zarr",
}

SUPPORTED_SUFFIXES: tuple[str, ...] = (ND2_SUFFIX, *PLUGIN_FOR_SUFFIX)

_BY_LENGTH = sorted(SUPPORTED_SUFFIXES, key=len, reverse=True)


def image_suffix(path: str | Path) -> str:
    """The format suffix of a file name, lower case, keeping an OME double suffix whole
    (``g.ome.tif`` -> ``.ome.tif``). An unknown suffix comes back as ``Path.suffix`` would give it."""
    name = Path(path).name.lower()
    for suf in _BY_LENGTH:
        if name.endswith(suf) and len(name) > len(suf):
            return suf
    return Path(name).suffix


def image_format(path: str | Path) -> str:
    """The reader family for a file name: ``"nd2"`` or ``"bioio"``. ValueError for an unknown suffix."""
    suffix = image_suffix(path)
    if suffix == ND2_SUFFIX:
        return "nd2"
    if suffix in PLUGIN_FOR_SUFFIX:
        return "bioio"
    raise ValueError(f"{Path(path).name}: unsupported image format {suffix!r}; "
                     f"supported: {', '.join(SUPPORTED_SUFFIXES)}")


def image_stem(path: str | Path) -> str:
    """The file name without its format suffix, which is the image id: ``g.nd2`` -> ``g``,
    ``g.ome.tif`` -> ``g``, ``g.zarr`` -> ``g``. For an ``.nd2`` (and any unknown suffix) this is exactly
    ``Path.stem``, so existing image ids do not change."""
    name = Path(path).name
    suffix = image_suffix(name)
    if suffix in SUPPORTED_SUFFIXES:
        return name[: len(name) - len(suffix)]
    return Path(name).stem
