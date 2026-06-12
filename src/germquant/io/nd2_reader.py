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


def read_nd2_metadata(path: str | Path) -> dict:
    """Cheap header read — no pixels loaded."""
    import nd2

    with nd2.ND2File(str(path)) as f:
        sp = _voxel(f)
        return {
            "sizes": dict(f.sizes),
            "spacing": sp if sp is not None else (1.0, 1.0, 1.0),  # (dz, dy, dx) microns
            "spacing_ok": sp is not None,
            "channel_names": _channel_names(f),
            "dtype": str(f.dtype),
            "is_2d": "Z" not in f.sizes,
        }


def read_stack(
    path: str | Path,
    *,
    xy_stride: int = 1,
    z_range: tuple[int, int] | None = None,
    as_float: bool = False,
) -> Stack:
    """Load a stack as (C, Z, Y, X). Singleton T/P loops are dropped (index 0)."""
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

    # drop acquisition loops we don't use (positions / time): take the first
    keep_dims: list[str] = []
    sl: list = []
    for d in dims:
        if d in ("T", "P"):
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
