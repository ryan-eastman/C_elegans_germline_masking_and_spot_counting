"""Truth points: the annotation files spot tuning scores against (design 6.5, 7.1 and 8.7; plan section 9.3).

A truth-point file is ``annotations/<image_id>__points_<annotator>.csv`` in a project, where the annotator
is 1 to 32 letters, digits, dashes or underscores. Its columns, in this order:

    point_id   0, 1, 2, ... (unique)
    z_um       the point's position in micrometres, on the image file's own frame: plane 0 of the file is
    y_um       z 0, pixel 0 is y and x 0, whatever stride and z range the run that was annotated used
    x_um       (z_um is 0.0 on a single-plane image)
    zi         the voxel of the point on the file's own grid (stride 1, the full z range):
    yi         round(z_um / dz), round(y_um / dy), round(x_um / dx) with the file's voxel size
    xi
    object_id  the run's label at the point (the label of the voxel, or the most frequent label right around
               it, as the RAD-51 harness maps expert spots; 0 outside every object)
    projected  true when the point was placed on a projection, so it is matched in 2D (design 6.7)

The file is ASCII with LF line endings, written atomically. ``write_points`` takes the points in the frame of
the run they were placed on (napari's world coordinates of that run's layers: z counted from the first plane
the run read, y and x scaled by the run's spacing) and fills the columns; ``read_points`` reads a file back
with its types; ``to_run_frame`` puts a file's points on the grid of another run of the same image. The
object ids are the RAD-51 harness's own ``label_at_points`` (germquant.benchmark.spots, one of the modules
the tune seal allows), so the plugin writes points without importing germquant.benchmark itself.
"""
from __future__ import annotations

import math
import os
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

POINT_COLUMNS = ["point_id", "z_um", "y_um", "x_um", "zi", "yi", "xi", "object_id", "projected"]
ANNOTATOR_RE = re.compile(r"[A-Za-z0-9_-]{1,32}")
POINTS_NAME_RE = re.compile(r"(?P<image_id>.+)__points_(?P<annotator>[A-Za-z0-9_-]{1,32})\.csv")
_COORD_DECIMALS = 6          # micrometres to the picometre: drops float noise such as 2.8000000000000003
_TRUE = {"true", "1", "yes"}
_FALSE = {"false", "0", "no", ""}


def points_file_name(image_id: str, annotator: str) -> str:
    """``<image_id>__points_<annotator>.csv``; ValueError for an annotator outside ^[A-Za-z0-9_-]{1,32}$ or an
    image id that cannot name a file. An image id may hold non-ASCII letters, as a project's image ids may (they
    come from file names); the file's content stays ASCII."""
    if not isinstance(annotator, str) or not ANNOTATOR_RE.fullmatch(annotator):
        raise ValueError(f"annotator {annotator!r}: use 1 to 32 letters, digits, '-' or '_'")
    if not isinstance(image_id, str) or not image_id.strip() or any(c in image_id for c in '\\/:*?"<>|') \
            or any(ord(c) < 32 for c in image_id):
        raise ValueError(f"image id {image_id!r} cannot name an annotation file")
    return f"{image_id}__points_{annotator}.csv"


def parse_points_name(name: str) -> tuple[str, str] | None:
    """(image_id, annotator) of a truth-point file name, or None when the name is not in that form."""
    m = POINTS_NAME_RE.fullmatch(Path(name).name)
    return (m.group("image_id"), m.group("annotator")) if m else None


def _labels3d(labels) -> np.ndarray:
    lab = np.asarray(labels)
    if lab.ndim == 2:
        lab = lab[None]
    if lab.ndim != 3:
        raise ValueError(f"labels must be a (Z, Y, X) or (Y, X) label image, got shape {lab.shape}")
    return lab


def _spacing(spacing, single_plane: bool) -> tuple[float, float, float]:
    try:
        dz, dy, dx = (float(s) for s in spacing)
    except (TypeError, ValueError):
        raise ValueError(f"spacing must be the run's voxel size (dz, dy, dx) in micrometres, got {spacing!r}") from None
    if not (math.isfinite(dy) and dy > 0 and math.isfinite(dx) and dx > 0):
        raise ValueError(f"spacing {spacing!r}: dy and dx must be positive")
    if not single_plane and not (math.isfinite(dz) and dz > 0):
        raise ValueError(f"spacing {spacing!r}: dz must be positive for a stack with several planes")
    return dz, dy, dx


def _stride(xy_stride) -> int:
    s = int(xy_stride)
    if s < 1 or s != xy_stride:
        raise ValueError(f"xy_stride must be a whole number of at least 1, got {xy_stride!r}")
    return s


def _z0(z_range) -> int:
    if z_range is None:
        return 0
    z0 = int(z_range[0])
    if z0 < 0:
        raise ValueError(f"z_range {z_range!r}: the first plane must be at least 0")
    return z0


def _flags(projected, n: int) -> np.ndarray:
    if isinstance(projected, (bool, np.bool_)):
        return np.full(n, bool(projected))
    arr = np.asarray(list(projected))
    if arr.shape != (n,) or arr.dtype.kind != "b":
        raise ValueError(f"projected must be true, false or one boolean per point ({n}), got {projected!r}")
    return arr.astype(bool)


def write_points(path, zyx_um, *, labels, spacing, xy_stride: int = 1, z_range=None, projected=False) -> Path:
    """Write the truth-point file ``path`` for points placed on a run, and return its path.

    ``zyx_um`` holds the points as (z, y, x) micrometres in the run's frame (napari's world coordinates of the
    run's layers); ``labels`` is the run's label image (Z, Y, X), or (Y, X) for a single plane; ``spacing`` the
    run's voxel size (dz, dy, dx) as its tables record it, so dy and dx include the run's stride; ``xy_stride``
    and ``z_range`` the run's geometry; ``projected`` true, false or one flag per point. The file holds the
    positions on the image file's own frame and grid, with each point's object (see the module docstring), in
    ASCII with LF line endings, written atomically (an existing file is replaced)."""
    from germquant.benchmark.spots import label_at_points

    lab = _labels3d(labels)
    single = lab.shape[0] == 1
    dz, dy, dx = _spacing(spacing, single)
    s, z0 = _stride(xy_stride), _z0(z_range)
    pts = np.asarray(zyx_um, dtype=float)
    pts = pts.reshape(-1, 3) if pts.size else np.zeros((0, 3))
    if not np.isfinite(pts).all():
        raise ValueError("every point needs finite z, y and x coordinates")
    n = len(pts)
    zr = np.zeros(n, np.int64) if single else np.rint(pts[:, 0] / dz).astype(np.int64)
    yr = np.rint(pts[:, 1] / dy).astype(np.int64)
    xr = np.rint(pts[:, 2] / dx).astype(np.int64)
    obj = label_at_points(lab, zr, yr, xr) if n else np.zeros(0, np.int64)
    if single:                       # one plane: plane z0 of the file (0.0 on a single-plane image)
        z_file = np.full(n, z0 * dz if (math.isfinite(dz) and dz > 0) else 0.0)
    else:
        z_file = pts[:, 0] + z0 * dz
    df = pd.DataFrame({
        "point_id": np.arange(n, dtype=np.int64),
        "z_um": np.round(z_file, _COORD_DECIMALS) + 0.0,
        "y_um": np.round(pts[:, 1], _COORD_DECIMALS) + 0.0,
        "x_um": np.round(pts[:, 2], _COORD_DECIMALS) + 0.0,
        "zi": zr + z0,
        "yi": np.rint(pts[:, 1] / (dy / s)).astype(np.int64),
        "xi": np.rint(pts[:, 2] / (dx / s)).astype(np.int64),
        "object_id": np.asarray(obj, np.int64),
        "projected": _flags(projected, n),
    }, columns=POINT_COLUMNS)
    out = Path(os.fspath(path))
    write_text_atomic(out, points_csv_text(df))
    return out


def points_csv_text(df: pd.DataFrame) -> str:
    """The CSV text of a points table in the file format (the columns of POINT_COLUMNS, then any extra
    columns; projected as true or false; LF line endings)."""
    cols = POINT_COLUMNS + [c for c in df.columns if c not in POINT_COLUMNS]
    out = df[cols].copy()
    out["projected"] = np.where(out["projected"].astype(bool), "true", "false")
    text = out.to_csv(index=False, lineterminator="\n")
    if not text.isascii():
        raise ValueError("a truth-point file is ASCII; an extra column holds other characters")
    return text


def read_points(path) -> pd.DataFrame:
    """The truth-point file ``path`` as a DataFrame: the columns of POINT_COLUMNS first (point_id, zi, yi, xi
    and object_id as integers, the positions as floats, projected as booleans), then any extra columns.
    ValueError naming the problem for a missing column, an empty cell, a repeated point_id or a projected
    value other than true or false; a byte-order mark is accepted."""
    from germquant.fsutil import long_path

    p = Path(os.fspath(path))
    try:
        df = pd.read_csv(long_path(p), encoding="utf-8-sig", dtype={"projected": str}, keep_default_na=False,
                         na_values=[""])
    except pd.errors.EmptyDataError:
        raise ValueError(f"{p}: empty file; a truth-point file has the header {','.join(POINT_COLUMNS)}") from None
    return check_points(df, p)


def check_points(df: pd.DataFrame, p="points") -> pd.DataFrame:
    """A points table checked and typed as ``read_points`` returns it (``p`` names it in the messages); the
    tuners pass a DataFrame item through it."""
    missing = [c for c in POINT_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{p}: missing columns {missing}; the header is {','.join(POINT_COLUMNS)}")
    num = [c for c in POINT_COLUMNS if c != "projected"]
    for c in num:
        bad = df[c].isna() | pd.to_numeric(df[c], errors="coerce").isna()
        if bad.any():
            raise ValueError(f"{p}: column {c} needs a number in every row (row {int(np.flatnonzero(bad)[0]) + 1})")
    out = df.copy()
    for c in ("point_id", "zi", "yi", "xi", "object_id"):
        vals = pd.to_numeric(out[c])
        if (vals != np.round(vals)).any():
            raise ValueError(f"{p}: column {c} holds whole numbers only")
        out[c] = vals.astype(np.int64)
    for c in ("z_um", "y_um", "x_um"):
        out[c] = pd.to_numeric(out[c]).astype(float)
    flags = out["projected"].fillna("").astype(str).str.strip().str.lower()
    unknown = ~flags.isin(_TRUE | _FALSE)
    if unknown.any():
        raise ValueError(f"{p}: projected must be true or false (row {int(np.flatnonzero(unknown)[0]) + 1}: "
                         f"{out['projected'].iloc[int(np.flatnonzero(unknown)[0])]!r})")
    out["projected"] = flags.isin(_TRUE).to_numpy()
    if out["point_id"].duplicated().any():
        raise ValueError(f"{p}: point_id {int(out.loc[out['point_id'].duplicated(), 'point_id'].iloc[0])} appears twice")
    if (out["object_id"] < 0).any():
        raise ValueError(f"{p}: object_id is 0 or a label, never negative")
    return out[POINT_COLUMNS + [c for c in out.columns if c not in POINT_COLUMNS]].reset_index(drop=True)


def to_run_frame(points: pd.DataFrame, *, shape, spacing, xy_stride: int = 1, z_range=None) -> dict:
    """A points table on the grid of a run of the same image: the run's ``shape`` (Z, Y, X) after its stride
    and z range, its ``spacing`` (dz, dy, dx; dy and dx include the stride), ``xy_stride`` and ``z_range``.

    Returns ``z_um``, ``y_um``, ``x_um`` (the run's frame, as its spots table writes positions: z counted
    from the run's first plane; 0.0 on a single plane), ``zr``, ``yr``, ``xr`` (the run's voxel of each
    point) and ``inside`` (the voxel lies in the run's field: a point on a plane the run did not read, or
    outside the image, cannot be scored)."""
    z_n, y_n, x_n = (int(v) for v in shape)
    single = z_n == 1
    dz, dy, dx = _spacing(spacing, single)
    z0 = _z0(z_range)
    _stride(xy_stride)
    zf = points["z_um"].to_numpy(dtype=float)
    y = points["y_um"].to_numpy(dtype=float)
    x = points["x_um"].to_numpy(dtype=float)
    if single:                       # one plane: only the points on plane z0 of the file are in it
        z = np.zeros(len(points))
        zr = points["zi"].to_numpy(dtype=np.int64) - z0
    else:
        z = zf - z0 * dz
        zr = np.rint(z / dz).astype(np.int64)
    yr = np.rint(y / dy).astype(np.int64)
    xr = np.rint(x / dx).astype(np.int64)
    inside = (zr >= 0) & (zr < z_n) & (yr >= 0) & (yr < y_n) & (xr >= 0) & (xr < x_n)
    return {"z_um": z, "y_um": y, "x_um": x, "zr": zr, "yr": yr, "xr": xr, "inside": inside}


def write_text_atomic(path: Path, text: str) -> Path:
    """Write ASCII ``text`` to ``path`` atomically: a temporary file next to it, flushed and fsynced, then
    os.replace through fsutil.long_path, retried five times over two seconds on PermissionError (an antivirus
    scanner or an indexer holding the old file on Windows). ValueError for text that is not ASCII."""
    from germquant.fsutil import long_path

    if not text.isascii():
        raise ValueError(f"{path}: tuning files are ASCII")
    path = Path(path)
    os.makedirs(long_path(path.parent), exist_ok=True)
    tmp = path.with_name(f"{path.name}.{os.getpid()}.tmp")
    try:
        with open(long_path(tmp), "w", encoding="ascii", newline="\n") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        for attempt in range(6):
            try:
                os.replace(long_path(tmp), long_path(path))
                break
            except PermissionError:
                if attempt == 5:
                    raise
                time.sleep(0.4)
    finally:
        if os.path.exists(long_path(tmp)):
            try:
                os.remove(long_path(tmp))
            except OSError:
                pass
    return path
