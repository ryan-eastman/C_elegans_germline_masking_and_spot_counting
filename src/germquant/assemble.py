"""Stack image files into one OME-TIFF with channel names and a voxel size (docs/PACKS.md, "Images with one file
per channel"; docs/DESIGN_system_packs.md, section 3.13, package P3).

germquant reads one image per run and binds channels to roles by name. Many labs and tools (Cell-ACDC among
them) keep one TIFF per channel instead (<stem>_Bright.tif, <stem>_Cy3.tif), and public datasets ship
multi-channel TIFFs whose page order a reader would have to guess. `germquant assemble` writes the one
OME-TIFF germquant reads:

    germquant assemble --channels A.tif B.tif --names Bright Cy3 --spacing 0.24 0.072 0.072 --out X.ome.tif

What it is told, never guessed:

* ``axes``: the declared axes of each input file, from the letters C, Z, Y and X, outermost first, written
  as an array would be indexed. The file's pages are taken in file order (tifffile's first series) and the
  page axes (YX, or YX plus a samples axis) must close the declared axes; the axes before them split the
  pages, in file order. ``ZYX`` (the default) is a z-stack of grey pages, or a single plane; ``ZCYX`` is a
  stack whose channels alternate page by page (DeepFoci's "XYCZ order" as ImageJ stores it); ``ZYXC`` is a
  stack whose pages hold the channels as samples (an RGB-like page per plane). A leading C is sized by the
  channel names (divided among the input files), a leading Z by the pages left.
* ``names``: one name per output channel, in order (every channel of every input, in input order). The
  channel map of a pack binds roles to these names.
* ``spacing``: the voxel size dz, dy, dx in micrometres (a z-stack's file often records no z step).
* ``crop``: optional half-open ranges of the input's pixels, ``{"z": [z0, z1], "y": [y0, y1], "x": [x0, x1]}``.

The declared axes are what is read, and nothing is guessed for a file without metadata. A file whose own ImageJ
or OME metadata records both channels and planes is the exception: a declaration that runs its pages the other
way (CZYX for a ZCYX hyperstack), or that leaves out its channels or its planes, would scramble the channels
without a sign, so it is refused unless ``force`` (``--force``) says to read the pages as declared anyway. The
output never replaces one of the inputs.

The output is CZYX in the inputs' data type, written with tifffile (zlib, OME-XML with the channel names and
PhysicalSizeZ, Y and X in micrometres). ``pixels_sha256`` pins the pixels independent of the file around them:
tifffile writes a fresh UUID into every OME-XML header, so two assemblies of the same pixels differ as files
but share pixels_sha256 (docs/DESIGN_system_packs.md, section 6.9). It is the sha256 of a header line
``germquant-pixels 1 CZYX <dtype> <C> <Z> <Y> <X>`` and a line feed, then the pixels as C-order little-endian
bytes.

Only numpy and tifffile are imported (no bioio, torch or CuPy), so the torch-free example fetcher
(``python -m germquant.packs fetch-example``) assembles with it.
"""
from __future__ import annotations

import hashlib
import os
import re
from collections.abc import Sequence
from pathlib import Path

import numpy as np

OUT_AXES = "CZYX"
MICRON = chr(0xB5) + "m"          # the OME unit name (the schema refuses "um"); built with chr() to keep this ASCII
_AXES_RE = re.compile(r"^[CZYX]+$")


class AssembleError(ValueError):
    """Inputs or settings that cannot make the declared image (axes, names, sizes, crop)."""


def pixels_sha256(data: np.ndarray, axes: str = OUT_AXES) -> str:
    """The sha256 of an image's pixels (module docstring): a header line with the axes, data type and shape,
    then the C-order little-endian bytes. ``data`` is normally a CZYX array."""
    a = np.ascontiguousarray(data)
    if a.dtype.byteorder == ">":
        a = a.astype(a.dtype.newbyteorder("<"))
    head = f"germquant-pixels 1 {axes} {a.dtype.str} {' '.join(str(int(n)) for n in a.shape)}\n"
    h = hashlib.sha256(head.encode("ascii"))
    h.update(memoryview(a.reshape(-1)).cast("B"))
    return h.hexdigest()


def _parse_axes(axes: str) -> str:
    ax = str(axes or "").strip().upper()
    if not _AXES_RE.fullmatch(ax) or len(set(ax)) != len(ax) or "Y" not in ax or "X" not in ax:
        raise AssembleError(f"axes {axes!r}: give the letters C, Z, Y and X, each at most once, with Y and X "
                            f"(for example ZYX, ZCYX or ZYXC)")
    return ax


def _read_pages(path: Path) -> tuple[np.ndarray, str, tuple[str, str, dict] | None]:
    """(pages, page axes, metadata axes) of a TIFF's first series: every page in file order, shape (P, *page
    shape), with page axes YX, YXC or CYX (a samples axis read as channels); the metadata axes are (kind,
    series axes, sizes) when the file is an ImageJ or OME TIFF, else None."""
    import tifffile

    if not path.is_file():
        raise FileNotFoundError(f"input image not found: {path}")
    with tifffile.TiffFile(path) as tf:
        series = tf.series[0]
        data = series.asarray()
        page = series.pages[0]
        page_shape = tuple(int(n) for n in page.shape)
        page_axes = str(page.axes).upper().replace("S", "C")
        meta = None
        if tf.is_imagej or tf.is_ome:
            ax = str(series.axes).upper()
            meta = ("ImageJ" if tf.is_imagej else "OME", ax, dict(zip(ax, (int(n) for n in series.shape))))
    if not set(page_axes) <= set("CYX") or "Y" not in page_axes or "X" not in page_axes:
        raise AssembleError(f"{path.name}: pages with axes {page.axes!r} are not supported (grey pages, or pages "
                            f"with samples, are)")
    if data.size % int(np.prod(page_shape)):
        raise AssembleError(f"{path.name}: {data.shape} does not divide into pages of {page_shape}")
    return data.reshape((-1,) + page_shape), page_axes, meta


def _axes_conflict(declared: str, page_axes: str, meta: tuple[str, str, dict] | None) -> str | None:
    """Why the declared axes contradict a file's own ImageJ or OME metadata, or None. Checked only for grey
    pages and a file that records both channels and planes (more than one of each): the declaration must then
    split the pages by both, in the file's order."""
    if meta is None or "C" in page_axes:
        return None
    kind, file_axes, sizes = meta
    n_c, n_z = sizes.get("C", 1), sizes.get("Z", 1)
    if n_c < 2 or n_z < 2:
        return None
    lead = declared[:len(declared) - len(page_axes)]
    where = (f"its {kind} metadata gives the axes {file_axes} ({n_c} channels, {n_z} planes, pages "
             f"{'channel within plane' if file_axes.index('Z') < file_axes.index('C') else 'plane within channel'})")
    if "C" not in lead or "Z" not in lead:
        return f"{where}, but the declared axes {declared} do not split the pages by both"
    if (lead.index("C") < lead.index("Z")) != (file_axes.index("C") < file_axes.index("Z")):
        return f"{where}, but the declared axes {declared} run them the other way"
    return None


def _to_czyx(pages: np.ndarray, page_axes: str, axes: str, n_channels: int | None, name: str) -> np.ndarray:
    """One input's pages as a CZYX array, read with its declared ``axes`` (module docstring). ``n_channels`` sizes
    a leading C (the names this input carries)."""
    if not axes.endswith(page_axes):
        raise AssembleError(f"{name}: the declared axes {axes} must end with the axes of its pages, {page_axes} "
                            f"(the pages are {pages.shape[1:]})")
    lead = axes[:len(axes) - len(page_axes)]
    n_pages = pages.shape[0]
    sizes: dict[str, int] = {}
    if "C" in lead:
        if not n_channels or n_channels < 1:
            raise AssembleError(f"{name}: the declared axes {axes} split its pages by channel, so the channel names "
                                f"must divide evenly among the input files")
        sizes["C"] = int(n_channels)
    if "Z" in lead:
        rest = n_pages // sizes.get("C", 1)
        if rest * sizes.get("C", 1) != n_pages:
            raise AssembleError(f"{name}: {n_pages} pages do not split into {sizes.get('C', 1)} channels")
        sizes["Z"] = rest
    shape = [sizes[a] for a in lead]
    if int(np.prod(shape, dtype=np.int64)) != n_pages:
        raise AssembleError(f"{name}: the declared axes {axes} need {int(np.prod(shape, dtype=np.int64))} pages, "
                            f"the file has {n_pages}")
    arr = pages.reshape(shape + list(pages.shape[1:]))
    present = lead + page_axes
    for a in OUT_AXES:
        if a not in present:
            arr = arr[np.newaxis]
            present = a + present
    return np.transpose(arr, [present.index(a) for a in OUT_AXES])


def _parse_crop(crop, shape_czyx: tuple[int, ...]) -> dict[str, list[int]] | None:
    if not crop:
        return None
    if not isinstance(crop, dict):
        raise AssembleError(f"crop must be a mapping such as {{'y': [0, 512], 'x': [0, 512]}}, got {crop!r}")
    out: dict[str, list[int]] = {}
    sizes = dict(zip(OUT_AXES, shape_czyx))
    for key, rng in crop.items():
        ax = str(key).upper()
        if ax not in "ZYX" or len(ax) != 1:
            raise AssembleError(f"crop axis {key!r}: only z, y and x can be cropped")
        if (not isinstance(rng, (list, tuple)) or len(rng) != 2
                or not all(isinstance(v, (int, np.integer)) and not isinstance(v, bool) for v in rng)):
            raise AssembleError(f"crop {key}: give [start, stop) as two whole numbers, got {rng!r}")
        a, b = int(rng[0]), int(rng[1])
        if not 0 <= a < b <= sizes[ax]:
            raise AssembleError(f"crop {key} [{a}, {b}) lies outside the {sizes[ax]} pixels of the image")
        out[ax.lower()] = [a, b]
    return out


def _spacing(spacing_zyx_um) -> tuple[float, float, float]:
    try:
        vals = tuple(float(v) for v in spacing_zyx_um)
    except (TypeError, ValueError):
        vals = ()
    if len(vals) != 3 or not all(np.isfinite(v) and v > 0 for v in vals):
        raise AssembleError(f"spacing must be three positive numbers dz dy dx in micrometres, got {spacing_zyx_um!r}")
    return vals


def build(channels: Sequence[str | Path], names: Sequence[str], *, axes: str = "ZYX", crop=None,
          force: bool = False) -> tuple[np.ndarray, dict | None]:
    """(CZYX array, the crop applied) from the input files, without writing anything (module docstring).
    ``force`` reads the pages as declared even where a file's own metadata says otherwise."""
    files = [Path(p) for p in channels]
    names = [str(n) for n in names]
    if not files:
        raise AssembleError("give at least one input file")
    if not names or any(not n.strip() for n in names) or len(set(names)) != len(names):
        raise AssembleError(f"channel names must be unique and not empty, got {names!r}")
    ax = _parse_axes(axes)
    # a leading C (one that splits the pages, not a samples axis of each page) is sized by the names each input
    # carries; whether C splits pages shows only once a file's page axes are read (_to_czyx)
    per_file = len(names) // len(files) if len(names) % len(files) == 0 else None
    blocks = []
    for f in files:
        pages, page_axes, meta = _read_pages(f)
        conflict = None if force else _axes_conflict(ax, page_axes, meta)
        if conflict:
            raise AssembleError(f"{f.name}: {conflict}. Declare axes that match the file, or pass force (--force) to "
                                f"read the pages as declared")
        blocks.append(_to_czyx(pages, page_axes, ax, per_file, f.name))
    first = blocks[0]
    for f, b in zip(files[1:], blocks[1:]):
        if b.shape[1:] != first.shape[1:]:
            raise AssembleError(f"{f.name}: ZYX {b.shape[1:]} differs from {files[0].name}'s {first.shape[1:]}")
        if b.dtype != first.dtype:
            raise AssembleError(f"{f.name}: data type {b.dtype} differs from {files[0].name}'s {first.dtype}")
    data = np.concatenate(blocks, axis=0) if len(blocks) > 1 else first
    if data.shape[0] != len(names):
        raise AssembleError(f"the inputs hold {data.shape[0]} channels, {len(names)} names were given ({names})")
    crop_used = _parse_crop(crop, data.shape)
    if crop_used:
        sl = [slice(None)] + [slice(*crop_used[a]) if a in crop_used else slice(None) for a in "zyx"]
        data = data[tuple(sl)]
    return np.ascontiguousarray(data), crop_used


def write_ome_tiff(path: str | Path, data_czyx: np.ndarray, names: Sequence[str],
                   spacing_zyx_um: Sequence[float]) -> Path:
    """Write a CZYX array as an OME-TIFF with channel names and voxel size (zlib), through a temporary file
    that replaces ``path`` only when complete."""
    import tifffile

    dz, dy, dx = _spacing(spacing_zyx_um)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    meta = {"axes": OUT_AXES, "Channel": {"Name": [str(n) for n in names]},
            "PhysicalSizeZ": dz, "PhysicalSizeY": dy, "PhysicalSizeX": dx,
            "PhysicalSizeZUnit": MICRON, "PhysicalSizeYUnit": MICRON, "PhysicalSizeXUnit": MICRON}
    try:
        tifffile.imwrite(tmp, data_czyx, ome=True, photometric="minisblack", compression="zlib", metadata=meta)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)
    return path


def assemble(channels: Sequence[str | Path], names: Sequence[str], spacing_zyx_um: Sequence[float],
             out: str | Path, *, axes: str = "ZYX", crop=None, force: bool = False) -> dict:
    """Stack ``channels`` (input files) into the OME-TIFF ``out`` (module docstring). Returns ``{"path",
    "axes", "shape", "dtype", "channel_names", "spacing_zyx_um", "crop", "pixels_sha256", "inputs"}``.
    AssembleError when ``out`` is one of the inputs (the original would be replaced)."""
    out_r = Path(out).resolve()                      # a symlink to an input resolves to the input
    for p in channels:
        if Path(p).resolve() == out_r:
            raise AssembleError(f"--out {out} is one of the inputs ({p}); write the OME-TIFF to a new file, so the "
                                f"original stays as it is")
    spacing = _spacing(spacing_zyx_um)
    data, crop_used = build(channels, names, axes=axes, crop=crop, force=force)
    path = write_ome_tiff(out, data, names, spacing)
    return {"path": str(path), "axes": OUT_AXES, "shape": [int(n) for n in data.shape], "dtype": str(data.dtype),
            "channel_names": [str(n) for n in names], "spacing_zyx_um": list(spacing), "crop": crop_used,
            "pixels_sha256": pixels_sha256(data), "inputs": [Path(p).name for p in channels],
            "declared_axes": _parse_axes(axes)}


def read_ome_tiff(path: str | Path) -> dict:
    """``{"data" (CZYX), "channel_names", "spacing_zyx_um" (None when not recorded), "pixels_sha256"}`` of an
    OME-TIFF such as assemble writes, read with tifffile alone."""
    import tifffile

    from .io.ome_header import ome_fields

    path = Path(path)
    with tifffile.TiffFile(path) as tf:
        series = tf.series[0]
        axes = str(series.axes).upper()
        data = series.asarray()
        names, (pz, py, px) = ome_fields(tf.ome_metadata)
    for a in list(axes):
        if a not in OUT_AXES:
            if data.shape[axes.index(a)] != 1:
                raise AssembleError(f"{path.name}: axis {a} of size {data.shape[axes.index(a)]} is not supported")
            data = np.take(data, 0, axis=axes.index(a))
            axes = axes.replace(a, "")
    for a in OUT_AXES:
        if a not in axes:
            data = data[np.newaxis]
            axes = a + axes
    data = np.ascontiguousarray(np.transpose(data, [axes.index(a) for a in OUT_AXES]))
    spacing = [pz, py, px] if all(v is not None for v in (pz, py, px)) else None
    return {"data": data, "channel_names": list(names), "spacing_zyx_um": spacing,
            "pixels_sha256": pixels_sha256(data)}


def file_pixels_sha256(path: str | Path) -> str:
    """pixels_sha256 of an OME-TIFF on disk (read_ome_tiff)."""
    return read_ome_tiff(path)["pixels_sha256"]
