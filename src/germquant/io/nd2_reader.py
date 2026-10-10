"""Read Nikon .nd2 z-stacks into a canonical (C, Z, Y, X) array + real voxel size.

Uses the `nd2` package (pure-python, fast metadata, lazy dask access). Voxel size is
pulled ONCE here as (dz, dy, dx) microns and carried on the Stack — every downstream
3D op takes it as `spacing`. Supports cheap sub-sampling (xy_stride / z_range) so a
2600x2600x70 stack can be smoke-tested on a laptop without loading ~3 GB.
"""
from __future__ import annotations

import dataclasses
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


def _time_loop_period_s(f) -> float | None:
    """The TimeLoop period (s) of the experiment, or None (a fixed image, or no readable loop).
    PACKAGE A2 STUB: package A1 supplies the real time-lapse reader (period_s, per-frame t_seconds); the
    merge keeps A1's version of this file."""
    try:
        for loop in f.experiment or []:
            if getattr(loop, "type", "") == "TimeLoop":
                ms = float(loop.parameters.periodMs)
                return ms / 1000.0 if ms > 0 else None
    except Exception:
        return None
    return None


def read_nd2_metadata(path: str | Path) -> dict:
    """Cheap header read — no pixels loaded. ``n_timepoints`` (1 for a fixed image), ``period_s`` and
    ``t_seconds`` describe a time-lapse (A2 stub, see _time_loop_period_s; package A1 owns the real fields)."""
    import nd2

    with nd2.ND2File(str(path)) as f:
        sp = _voxel(f)
        sizes = dict(f.sizes)
        n_t = int(sizes.get("T", 1))
        period_s = _time_loop_period_s(f) if "T" in sizes else None
        return {
            "sizes": sizes,
            "spacing": sp if sp is not None else (1.0, 1.0, 1.0),  # (dz, dy, dx) microns
            "spacing_ok": sp is not None,
            "channel_names": _channel_names(f),
            "dtype": str(f.dtype),
            "is_2d": "Z" not in f.sizes,
            "n_timepoints": n_t,
            "period_s": period_s,
            "t_seconds": [k * period_s for k in range(n_t)] if period_s is not None else None,
        }


def read_stack(
    path: str | Path,
    *,
    xy_stride: int = 1,
    z_range: tuple[int, int] | None = None,
    as_float: bool = False,
    t: int | None = None,
) -> Stack:
    """Load a stack as (C, Z, Y, X). Singleton T/P loops are dropped (index 0). ``t`` selects a timepoint of a
    T axis instead of index 0 (only that timepoint's planes are read through the dask array); ValueError when
    it is out of range, or when the file has no T axis and t is neither None nor 0. PACKAGE A2 STUB of the
    reader contract: package A1 supplies the real one and the merge keeps A1's version of this file."""
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

    if t is not None:
        t = int(t)
        n_t = int(sizes.get("T", 1))
        if "T" not in sizes and t != 0:
            raise ValueError(f"{path.name} has no T axis: t={t} is not a timepoint of it (only t=0 or None)")
        if not 0 <= t < n_t:
            raise ValueError(f"t={t} is out of range for {path.name} (T={n_t})")

    # drop acquisition loops we don't use (positions / time): take the first (or the requested timepoint)
    keep_dims: list[str] = []
    sl: list = []
    for d in dims:
        if d == "T":
            sl.append(0 if t is None else t)
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
    )
