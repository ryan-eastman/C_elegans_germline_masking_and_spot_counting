"""Synthetic single-plane images for the single-plane tests (tranche 5 package P2b; docs/PLANAR.md).

One calibrated 256 x 256 plane with two channels, named DAPI and foci (the names the generic channel map
resolves): nine nuclei planted as flat disks on a flat background, two specks of debris that the area gate of
config/profiles/generic_2d.yaml (min_area_um2 20.0) drops, and two to five Gaussian puncta per nucleus in the
foci channel. The pixel size is 0.8 um, square. At that pixel size the classical watershed (sigma 0.4 um)
recovers every planted disk pixel for pixel, so the area check of the synthetic_2d golden (design 12.3, P2b
acceptance 2) tests the planar measure, not the segmenter's edge.

The plane is written with tifffile in four forms: an ImageJ TIFF with its xy calibration, an OME-TIFF YX (the
DAPI channel), an OME-TIFF CYX and an OME-TIFF CZYX with Z = 1 (that one with a z step); and as a five-plane
CZYX stack for the case of a 3D file cut to one plane. The pixels are uint16 in every form, and the same array
feeds the reader seam of process_image (patch_reader), so a run through bioio and a run through the seam see
the same data.

compare_tables and write_golden are the golden rules of tests/test_golden.py: exact on the workstation, and in
CI the relative tolerance GERMQUANT_GOLDEN_RTOL with the spots in coordinate order and spot_id as a set.

This module is imported by tests that also run on the tranche 5 base commit (the single_plane_legacy golden is
made there, docs/PLANAR.md), so it uses nothing that the base lacks at import time.
"""
from __future__ import annotations

import functools
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PROFILE = ROOT / "config" / "profiles" / "generic_2d.yaml"

SIZE = 256                          # Y = X
DY = DX = 0.8                       # um
DZ = 0.5                            # um: the z step of the CZYX (Z = 1) form and of the 3D stack
CHANNELS = ("DAPI", "foci")         # the generic map resolves dna and foci by these names
STEM = "plane_001"                  # the image id (generic_2d's filename regex matches nothing)
N3D = 5                             # planes of the 3D stack
UM = chr(0x00B5) + "m"              # the OME unit of length (written with chr: no non-ASCII text in source)

_BG_DNA, _NUC_DNA, _NOISE_DNA = 200.0, 1400.0, 8.0
_BG_FOCI, _NUC_FOCI, _NOISE_FOCI = 150.0, 300.0, 12.0
_PUNCTUM_AMP, _PUNCTUM_SIGMA_UM = 1500.0, 0.6
_RADII_UM = (4.7, 5.3, 6.0)        # each disk's pixel count is within 1 percent of the circle's area
_DEBRIS_R_UM = 1.6
_SEED = 23


def _centres() -> list[tuple[int, int]]:
    step = SIZE // 3
    return [(step // 2 + i * step, step // 2 + j * step) for i in range(3) for j in range(3)]


@functools.lru_cache(maxsize=1)
def _make() -> tuple[np.ndarray, dict]:
    rng = np.random.default_rng(_SEED)
    yy, xx = np.indices((SIZE, SIZE)).astype(np.float64)
    dna = np.full((SIZE, SIZE), _BG_DNA)
    foci = np.full((SIZE, SIZE), _BG_FOCI)
    nuclei, puncta, debris = [], [], []
    for k, (cy, cx) in enumerate(_centres()):
        r_px = _RADII_UM[k % len(_RADII_UM)] / DY
        mask = (yy - cy) ** 2 + (xx - cx) ** 2 <= r_px ** 2
        dna[mask] = _NUC_DNA
        foci[mask] = _NUC_FOCI
        nuclei.append({"y": cy, "x": cx, "r_px": r_px, "n_px": int(mask.sum()),
                       "area_um2": float(mask.sum()) * DY * DX, "disk_um2": float(np.pi * (r_px * DY) ** 2)})
        n = 2 + k % 4
        for j in range(n):
            ang = 2 * np.pi * j / n + 0.5 * k
            puncta.append((k, cy + 0.55 * r_px * np.sin(ang), cx + 0.55 * r_px * np.cos(ang)))
    for _, py, px in puncta:
        s = _PUNCTUM_SIGMA_UM / DY
        foci += _PUNCTUM_AMP * np.exp(-((yy - py) ** 2 + (xx - px) ** 2) / (2 * s * s))
    step = SIZE // 3
    for cy, cx in ((step, step), (2 * step, 2 * step)):         # between the nuclei, below the area gate
        r_px = _DEBRIS_R_UM / DY
        mask = (yy - cy) ** 2 + (xx - cx) ** 2 <= r_px ** 2
        dna[mask] = _NUC_DNA
        debris.append({"y": cy, "x": cx, "n_px": int(mask.sum()), "area_um2": float(mask.sum()) * DY * DX})
    dna += rng.normal(0.0, _NOISE_DNA, dna.shape)
    foci += rng.normal(0.0, _NOISE_FOCI, foci.shape)
    data = np.clip(np.round(np.stack([dna, foci])), 0, 65535).astype(np.uint16)
    return data, {"nuclei": nuclei, "puncta": puncta, "debris": debris}


def plane_cyx() -> np.ndarray:
    """The (2, 256, 256) uint16 plane (a fresh copy)."""
    return _make()[0].copy()


def planted() -> dict:
    """What was planted, in pixels of the plane: ``nuclei`` [{y, x, r_px, n_px, area_um2, disk_um2}] (area_um2
    is the planted disk's pixel count times dy * dx, disk_um2 the area of the circle), ``puncta`` [(nucleus
    index, y, x)] and ``debris`` [{y, x, n_px, area_um2}]."""
    info = _make()[1]
    return {k: [dict(v) if isinstance(v, dict) else v for v in info[k]] for k in info}


def stack_3d() -> np.ndarray:
    """(2, 5, 256, 256) uint16: the plane on z = 2, dimmer copies of it on the other planes."""
    gain = (0.6, 0.8, 1.0, 0.8, 0.6)
    p = plane_cyx().astype(np.float64)
    return np.stack([np.round(p * g) for g in gain], axis=1).astype(np.uint16)


# ---- files ---------------------------------------------------------------------------------------------

def _ome(path: Path, data: np.ndarray, axes: str, names, dz: float | None = None) -> Path:
    import tifffile

    md = {"axes": axes, "PhysicalSizeX": DX, "PhysicalSizeXUnit": UM, "PhysicalSizeY": DY, "PhysicalSizeYUnit": UM,
          "Channel": {"Name": list(names)}}
    if dz is not None:
        md.update({"PhysicalSizeZ": dz, "PhysicalSizeZUnit": UM})
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(path, data, ome=True, metadata=md)
    return path


def write_ome_cyx(path) -> Path:
    """OME-TIFF, axes CYX, PhysicalSizeX and Y only (no z step)."""
    return _ome(path, plane_cyx(), "CYX", CHANNELS)


def write_ome_czyx_z1(path, dz: float = DZ) -> Path:
    """OME-TIFF, axes CZYX with Z = 1, with a z step (PhysicalSizeZ ``dz``)."""
    return _ome(path, plane_cyx()[:, None], "CZYX", CHANNELS, dz=dz)


def write_ome_yx(path) -> Path:
    """OME-TIFF, axes YX: the DAPI channel alone."""
    return _ome(path, plane_cyx()[0], "YX", CHANNELS[:1])


def write_tiff(path, calibrated: bool = True) -> Path:
    """An ImageJ TIFF, axes CYX, with the xy calibration in its resolution tags (unit um); with ``calibrated``
    False a plain tifffile TIFF that carries no calibration (resolution unit none)."""
    import tifffile

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if calibrated:
        tifffile.imwrite(path, plane_cyx(), imagej=True, resolution=(1.0 / DX, 1.0 / DY),
                         metadata={"axes": "CYX", "unit": "um"})
    else:
        tifffile.imwrite(path, plane_cyx())
    return path


def write_ome_3d(path, dz: float = DZ) -> Path:
    """OME-TIFF, axes CZYX with Z = 5 (stack_3d), calibrated in z, y and x."""
    return _ome(path, stack_3d(), "CZYX", CHANNELS, dz=dz)


FORMS = {"tiff": write_tiff, "ome_yx": write_ome_yx, "ome_cyx": write_ome_cyx, "ome_czyx_z1": write_ome_czyx_z1}


# ---- the reader seam of process_image ------------------------------------------------------------------

def nd2_meta(data: np.ndarray, spacing, *, spacing_ok: bool = True, n_timepoints: int = 1, two_d: bool = True) -> dict:
    """A header as read_nd2_metadata returns it for ``data`` (C, Z, Y, X): an .nd2 single plane has no Z axis
    (``two_d``), a series a T axis."""
    c, z, y, x = data.shape
    sizes = ({"T": n_timepoints} if n_timepoints > 1 else {}) | {"C": c} | ({} if two_d else {"Z": z}) | {"Y": y, "X": x}
    return {"sizes": sizes, "spacing": tuple(spacing), "spacing_ok": spacing_ok, "channel_names": list(CHANNELS),
            "dtype": str(data.dtype), "is_2d": "Z" not in sizes, "n_timepoints": n_timepoints, "period_s": None,
            "t_seconds": None}


def patch_reader(monkeypatch, data: np.ndarray | None = None, *, dz: float = DZ, spacing_ok: bool = True,
                 two_d: bool | None = None, series: int = 1, calls: list | None = None) -> np.ndarray:
    """Feed process_image through its .nd2 reader names (pipeline.read_nd2_metadata, pipeline.read_stack) with
    ``data`` (C, Z, Y, X; the plane by default) at spacing (dz, DY, DX). The fake read_stack honours xy_stride,
    z_range and t the way nd2_reader.read_stack does (stride scales dy and dx). ``series`` > 1 makes a T axis
    of that many identical timepoints. ``calls`` collects (xy_stride, z_range, t) per read. Returns the data."""
    from germquant import pipeline
    from germquant.io.nd2_reader import Stack

    data = plane_cyx()[:, None] if data is None else np.asarray(data)
    two_d = data.shape[1] == 1 if two_d is None else two_d
    meta = nd2_meta(data, (dz, DY, DX), spacing_ok=spacing_ok, n_timepoints=series, two_d=two_d)

    def read_meta(path, t=None):
        return dict(meta)

    def read_stack(path, *, xy_stride=1, z_range=None, as_float=False, t=None):
        if calls is not None:
            calls.append((xy_stride, z_range, t))
        arr = data if z_range is None else data[:, z_range[0]:z_range[1]]
        arr = arr[:, :, ::xy_stride, ::xy_stride] if xy_stride > 1 else arr
        sp = (dz, DY * xy_stride, DX * xy_stride) if spacing_ok else (1.0, 1.0 * xy_stride, 1.0 * xy_stride)
        return Stack(data=np.ascontiguousarray(arr), spacing=sp, channel_names=list(CHANNELS), path=Path(path),
                     downsample_xy=xy_stride, spacing_ok=spacing_ok, t_index=t)

    monkeypatch.setattr(pipeline, "read_nd2_metadata", read_meta)
    monkeypatch.setattr(pipeline, "read_stack", read_stack)
    return data


def generic_2d(**overrides):
    """config/profiles/generic_2d.yaml with the classical segmenter forced the way tests/test_golden.py forces it
    (outside the hash, so the config hash is the profile's own), then ``overrides`` applied with Config.set."""
    from germquant.config import load_config

    cfg = load_config(PROFILE)
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    for key, value in overrides.items():
        cfg.set(key, value)
    return cfg


# ---- goldens -------------------------------------------------------------------------------------------

VOLATILE = ["git_sha", "run_timestamp", "pipeline_version", "file_path", "config_hash"]
_RTOL_ENV = os.environ.get("GERMQUANT_GOLDEN_RTOL")       # CI: another CPU's summation order (tests/README.md)
_FLOAT_TOLERANCE = ({"check_exact": False, "rtol": float(_RTOL_ENV), "atol": 0.0} if _RTOL_ENV
                    else {"check_exact": True})
_ORDER_KEYS = {"spots": ["z_um", "y_um", "x_um"]}
_ORDINAL_COLUMNS = {"spots": {"spot_id"}}


def read_table(p: Path) -> pd.DataFrame:
    df = pd.read_csv(p, low_memory=False, float_precision="round_trip")
    return df.drop(columns=[c for c in VOLATILE if c in df.columns])


def _canonical(df: pd.DataFrame, table: str) -> pd.DataFrame:
    keys = _ORDER_KEYS.get(table)
    if not _RTOL_ENV or not keys or not set(keys) <= set(df.columns) or df.empty:
        return df
    k = df[keys].astype(float).round(6)
    return df.loc[k.sort_values(keys, kind="mergesort").index].reset_index(drop=True)


def write_golden(out: Path, image_id: str, gold: Path, tables, config_hash: str, *, sentinel: Path | None = None) -> None:
    """Write the run's tables (volatile columns dropped) into ``gold`` and the config hash into
    ``sentinel`` (default ``gold``/config_hash.txt)."""
    gold.mkdir(parents=True, exist_ok=True)
    for t in tables:
        read_table(out / f"{image_id}__{t}.csv").to_csv(gold / f"{t}.csv", index=False)
    (sentinel or gold / "config_hash.txt").write_text(config_hash + "\n", encoding="utf-8")


def compare_tables(out: Path, image_id: str, gold: Path, tables) -> None:
    """Every table of the run equals the golden's, cell for cell (the rules of tests/test_golden.py)."""
    for t in tables:
        got = _canonical(read_table(out / f"{image_id}__{t}.csv"), t)
        ref = _canonical(pd.read_csv(gold / f"{t}.csv", low_memory=False, float_precision="round_trip"), t)
        common = [c for c in got.columns if c in set(ref.columns)]
        assert common == list(ref.columns), f"{t}: existing columns renamed, removed or reordered"
        assert len(got) == len(ref), f"{t}: row count {len(got)} vs {len(ref)}"
        for col in ref.columns:
            a, b = got[col], ref[col]
            if _RTOL_ENV and col in _ORDINAL_COLUMNS.get(t, ()):
                assert sorted(a.astype(float)) == sorted(b.astype(float)), f"{t}.{col}: not the same set of ids"
                continue
            if pd.api.types.is_numeric_dtype(b):
                pd.testing.assert_series_equal(a.astype(float), b.astype(float), check_names=False,
                                               obj=f"{t}.{col}", **_FLOAT_TOLERANCE)
            else:
                assert (a.astype("string").fillna("<NA>") == b.astype("string").fillna("<NA>")).all(), \
                    f"{t}.{col} differs"


def matched_puncta(spots: pd.DataFrame, labels2d: np.ndarray, max_px: float = 1.5) -> tuple[int, int]:
    """(planted puncta found, extra spots): a planted punctum is found when a spot lies within ``max_px``
    pixels of it in the nucleus it was planted in; every spot not matched to a punctum is an extra."""
    info = planted()
    yx = spots[["y_um", "x_um"]].to_numpy(dtype=float) / np.array([DY, DX]) if len(spots) else np.zeros((0, 2))
    used: set[int] = set()
    found = 0
    for k, py, px in info["puncta"]:
        nuc = info["nuclei"][k]
        want = int(labels2d[int(nuc["y"]), int(nuc["x"])])
        d = np.hypot(yx[:, 0] - py, yx[:, 1] - px) if len(yx) else np.zeros(0)
        for j in np.argsort(d):
            if d[j] > max_px:
                break
            if j not in used and int(spots["nucleus_id"].iloc[j]) == want:
                used.add(int(j))
                found += 1
                break
    return found, len(spots) - len(used)
