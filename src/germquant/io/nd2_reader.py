"""Read Nikon .nd2 z-stacks into a canonical (C, Z, Y, X) array + real voxel size.

Uses the `nd2` package (pure-python, fast metadata, lazy dask access). Voxel size is
pulled ONCE here as (dz, dy, dx) microns and carried on the Stack — every downstream
3D op takes it as `spacing`. Supports cheap sub-sampling (xy_stride / z_range) so a
2600x2600x70 stack can be smoke-tested on a laptop without loading ~3 GB.

Time series (docs/TIMELAPSE.md): ``read_stack(path, t=7)`` reads one timepoint of a T series. nd2 hands
out one dask chunk per (T, Z) frame, so only that timepoint's planes are read. ``t=None`` is the old
behaviour (a T axis takes index 0), and a file without a T axis is read exactly as before.
``read_nd2_metadata`` reports ``n_timepoints``, ``period_s`` and ``t_seconds`` for every file (1, None,
None for a fixed image).
"""
from __future__ import annotations

import dataclasses
import operator
from pathlib import Path

import numpy as np


@dataclasses.dataclass
class Stack:
    data: np.ndarray                      # (C, Z, Y, X)
    spacing: tuple[float, float, float]   # (dz, dy, dx) microns
    channel_names: list[str]
    path: Path
    downsample_xy: int = 1
    spacing_ok: bool = True               # False => voxel size unreadable, spacing is a guess
    t_index: int | None = None            # the timepoint read when one was asked for (read_stack(t=...)); None otherwise
    spacing_yx: tuple[float, float] | None = None   # (dy, dx) um the bioio reader finds in a one-plane file; None otherwise

    @property
    def n_channels(self) -> int:
        return self.data.shape[0]

    @property
    def shape_zyx(self) -> tuple[int, int, int]:
        return tuple(self.data.shape[1:])  # type: ignore[return-value]

    @property
    def anisotropy(self) -> float:
        """dz / dy — the z:xy ratio Cellpose and friends need."""
        return self.spacing[0] / self.spacing[1]

    def channel(self, idx: int | None) -> np.ndarray | None:
        if idx is None:
            return None
        return self.data[idx]


def _channel_names(f) -> list[str]:
    names: list[str] = []
    try:
        for c in (f.metadata.channels or []):
            names.append(str(getattr(c.channel, "name", "")))
    except Exception:
        pass
    return names


def _voxel(f) -> tuple[float, float, float] | None:
    """(dz, dy, dx) microns, or None if the .nd2 voxel size is unreadable/degenerate.

    Returning None (rather than silently defaulting to 1 µm isotropic) lets callers flag the
    file — a wrong spacing silently corrupts every length/area/volume (ARCHITECTURE.md §2).
    """
    try:
        vs = f.voxel_size()
        sp = (float(vs.z), float(vs.y), float(vs.x))
    except Exception:
        return None
    if not all(np.isfinite(s) and s > 0 for s in sp):
        return None
    return sp


def check_timepoint(t, n_timepoints: int, path: str | Path) -> int:
    """``t`` as a plain index into the T axis of a file with ``n_timepoints`` timepoints (1 for a file
    without a T axis). ValueError when ``t`` is not an integer, is out of range, or when the file has no T
    axis and ``t`` is not 0: a fixed image has exactly one timepoint, index 0, so ``t=0`` reads it as
    ``t=None`` does. Shared by the nd2 and the bioio reader so both refuse the same way."""
    if isinstance(t, bool):
        raise ValueError(f"{Path(path).name}: t must be an integer timepoint index, got {t!r}")
    try:
        ti = operator.index(t)
    except TypeError:
        raise ValueError(f"{Path(path).name}: t must be an integer timepoint index, got {t!r}") from None
    if 0 <= ti < n_timepoints:
        return ti
    if n_timepoints <= 1:
        raise ValueError(f"{Path(path).name} has no T axis (a single timepoint): t must be None or 0, got {ti}")
    raise ValueError(f"{Path(path).name}: t={ti} is out of range; the series has {n_timepoints} timepoints "
                     f"(0 to {n_timepoints - 1})")


def _time_loop(f):
    """The experiment's time loop (a TimeLoop, or an NETimeLoop for a non-equidistant series), or None
    when the file has none or the experiment record is unreadable."""
    try:
        for loop in f.experiment:
            if getattr(loop, "type", None) in ("TimeLoop", "NETimeLoop"):
                return loop
    except Exception:
        return None
    return None


def _period_s(f) -> float | None:
    """The programmed interval between timepoints in seconds (the TimeLoop ``periodMs``). None for a fixed
    image, for a non-equidistant loop with more than one phase, and for an unreadable or degenerate value."""
    loop = _time_loop(f)
    if loop is None:
        return None
    try:
        if loop.type == "TimeLoop":
            ms = loop.parameters.periodMs
        else:
            periods = list(loop.parameters.periods)
            if len(periods) != 1:
                return None
            ms = periods[0].periodMs
        s = float(ms) / 1000.0
    except Exception:
        return None
    return s if np.isfinite(s) and s > 0 else None


def _t_seconds(f, n_timepoints: int) -> list[float] | None:
    """The recorded time of each timepoint's first plane, in seconds since the acquisition started (the
    per-frame camera time stamp NIS-Elements stores, ``relativeTimeMs``), one value per timepoint. None for a
    file without a T axis and when the per-frame record is unreadable. Reads one small metadata chunk per
    timepoint, no pixels."""
    if "T" not in f.sizes:
        return None
    try:
        first: dict[int, int] = {}          # timepoint -> sequence index of its first frame
        for i, loop_index in enumerate(f.loop_indices):
            ti = loop_index.get("T")
            if ti is not None and int(ti) not in first:
                first[int(ti)] = i
        out = []
        for ti in range(n_timepoints):
            fm = f.frame_metadata(first[ti])
            out.append(float(fm.channels[0].time.relativeTimeMs) / 1000.0)
    except Exception:
        return None
    return out


def read_nd2_metadata(path: str | Path, *, t: int | None = None) -> dict:
    """Cheap header read: no pixels loaded.

    Besides the file's sizes, spacing, channel names and dtype, the dict says how many timepoints the file
    holds (``n_timepoints``, 1 for a fixed image), the programmed interval in seconds (``period_s``, None
    when there is no time loop) and the recorded time stamp of each timepoint's first plane in seconds since
    the acquisition started (``t_seconds``, a list, or None for a fixed image). ``t`` only checks that the
    timepoint exists (ValueError as for read_stack); the dict is the same for every timepoint.
    """
    import nd2

    with nd2.ND2File(str(path)) as f:
        sp = _voxel(f)
        n_timepoints = int(f.sizes.get("T", 1))
        if t is not None:
            check_timepoint(t, n_timepoints, path)
        return {
            "sizes": dict(f.sizes),
            "spacing": sp if sp is not None else (1.0, 1.0, 1.0),  # (dz, dy, dx) microns
            "spacing_ok": sp is not None,
            "channel_names": _channel_names(f),
            "dtype": str(f.dtype),
            "is_2d": "Z" not in f.sizes,
            "n_timepoints": n_timepoints,
            "period_s": _period_s(f),
            "t_seconds": _t_seconds(f, n_timepoints),
        }


def read_stack(
    path: str | Path,
    *,
    xy_stride: int = 1,
    z_range: tuple[int, int] | None = None,
    as_float: bool = False,
    t: int | None = None,
) -> Stack:
    """Load a stack as (C, Z, Y, X). Singleton T/P loops are dropped (index 0).

    ``t`` picks one timepoint of a time series: the T axis takes index ``t`` and, because nd2 chunks the
    lazy array one frame per (T, Z), only that timepoint's planes are read. ``t=None`` keeps the old
    behaviour (index 0 of a T axis, and a file without a T axis is read exactly as before). ValueError
    when ``t`` is out of range, or when the file has no T axis and ``t`` is not 0 (see check_timepoint).
    """
    import dask.array as da
    import nd2

    path = Path(path)
    with nd2.ND2File(str(path)) as f:
        sizes = dict(f.sizes)
        names = _channel_names(f)
        sp = _voxel(f)
        spacing_ok = sp is not None
        dz, dy, dx = sp if sp is not None else (1.0, 1.0, 1.0)
        darr = f.to_dask()
        dims = list(sizes.keys())

    # the timepoint to read: index t of the T axis when one was asked for, else the first (as always)
    t_idx = 0 if t is None else check_timepoint(t, int(sizes.get("T", 1)), path)

    # drop acquisition loops we don't use (positions / time): take the first (or the timepoint asked for)
    keep_dims: list[str] = []
    sl: list = []
    for d in dims:
        if d == "T":
            sl.append(t_idx)
        elif d == "P":
            sl.append(0)
        else:
            sl.append(slice(None))
            keep_dims.append(d)
    darr = darr[tuple(sl)]
    dims = keep_dims

    # ensure C and Z axes exist
    if "C" not in dims:
        darr = da.expand_dims(darr, 0)
        dims = ["C"] + dims
    if "Z" not in dims:
        ci = dims.index("C")
        darr = da.expand_dims(darr, ci + 1)
        dims = dims[: ci + 1] + ["Z"] + dims[ci + 1 :]

    # transpose to canonical C, Z, Y, X
    target = ["C", "Z", "Y", "X"]
    darr = darr.transpose([dims.index(t) for t in target])

    # cheap sub-sampling BEFORE compute
    if z_range is not None:
        darr = darr[:, z_range[0] : z_range[1]]
    if xy_stride > 1:
        darr = darr[:, :, ::xy_stride, ::xy_stride]
        dy *= xy_stride
        dx *= xy_stride

    data = np.asarray(darr)
    if as_float:
        data = data.astype(np.float32)

    return Stack(
        data=data,
        spacing=(dz, dy, dx),
        channel_names=names,
        path=path,
        downsample_xy=xy_stride,
        spacing_ok=spacing_ok,
        t_index=None if t is None else t_idx,
    )
