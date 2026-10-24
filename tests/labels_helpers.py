"""Helpers for the labels-from-files tests (tests/test_labels_io.py, tests/test_labels_run.py; docs/LABELS.md).

* writers for the formats segmentation.nuclei.labels_from reads: TIFF, OME-TIFF and Cell-ACDC segm.npz, plus a
  pickled npz and a Cellpose-style _seg.npy for the refusals;
* the planted truth labels of the synthetic gonad of tests/test_end_to_end.py (one ellipsoid per planted
  nucleus, the twelve germline nuclei and the three blobs of off-gonad debris);
* a reader patch for process_image (the synthetic stacks are never files) and the table comparison of the
  synthetic goldens (tests/test_golden.py's rules: exact on the workstation, GERMQUANT_GOLDEN_RTOL in CI).

A module, not a conftest.py: the suite shares helpers as modules (tests/packs_helpers.py).
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

from germquant import pipeline

# the synthetic gonad of tests/test_end_to_end.py: 20 x 80 x 320 voxels of (0.4, 0.2, 0.2) um, its nuclei at
# z = 10, y = 40 and these x positions (twelve germline nuclei, then three SYP-negative debris blobs)
GONAD_SHAPE = (20, 80, 320)
GONAD_SPACING = (0.4, 0.2, 0.2)
GONAD_X = [*range(20, 210, 17), 280, 295, 310]
TRUTH_RADIUS_UM = 1.4          # the planted blobs have sigma 1.6 um; 1.4 um keeps neighbours (3 um apart) apart
TEMPLATE = "{image_dir}/{image_stem}_truth_labels.tif"
VOLATILE = ["git_sha", "run_timestamp", "pipeline_version", "file_path", "config_hash"]


def truth_labels(shape=GONAD_SHAPE, spacing=GONAD_SPACING, xs=GONAD_X, radius_um=TRUTH_RADIUS_UM,
                 zc: int = 10, yc: int = 40) -> np.ndarray:
    """int32 (Z, Y, X): an ellipsoid of ``radius_um`` (micrometres, so anisotropic in voxels) around each planted
    nucleus, ids 1..N in the order of ``xs``."""
    zz, yy, xx = np.indices(shape, dtype=np.float32)
    lab = np.zeros(shape, np.int32)
    for i, cx in enumerate(xs, 1):
        d2 = (((zz - zc) * spacing[0]) ** 2 + ((yy - yc) * spacing[1]) ** 2 + ((xx - cx) * spacing[2]) ** 2)
        lab[d2 <= radius_um ** 2] = i
    return lab


def disks(shape_yx, centres, radius_px: int, first_id: int = 1) -> np.ndarray:
    """int32 (Y, X): one disk of ``radius_px`` per centre, ids from ``first_id`` in order."""
    yy, xx = np.indices(shape_yx)
    lab = np.zeros(shape_yx, np.int32)
    for i, (cy, cx) in enumerate(centres, first_id):
        lab[(yy - cy) ** 2 + (xx - cx) ** 2 <= radius_px ** 2] = i
    return lab


# ---- writers -------------------------------------------------------------------------------------------

def write_tiff(path, arr) -> Path:
    """A plain TIFF, as tifffile.imwrite writes it with no metadata."""
    import tifffile

    tifffile.imwrite(str(path), np.asarray(arr))
    return Path(path)


def write_ome(path, arr, axes: str) -> Path:
    """An OME-TIFF with the given axes (YX, ZYX, TYX, TZYX, CZYX ...)."""
    import tifffile

    tifffile.imwrite(str(path), np.asarray(arr), ome=True, metadata={"axes": axes})
    return Path(path)


def write_npz(path, arr) -> Path:
    """A Cell-ACDC style segm.npz: one array, saved the way Cell-ACDC saves it (np.savez_compressed)."""
    np.savez_compressed(str(path), np.asarray(arr))
    return Path(path)


def write_pickled_npz(path) -> Path:
    """An npz whose array holds Python objects, which only a pickle can load."""
    np.savez(str(path), np.array([{"masks": np.ones((4, 4), np.int32)}], dtype=object))
    return Path(path)


def write_seg_npy(path, masks) -> Path:
    """A Cellpose GUI style _seg.npy: a pickled dict with the masks under "masks"."""
    np.save(str(path), {"masks": np.asarray(masks), "outlines": np.zeros_like(masks), "filename": "x.tif"},
            allow_pickle=True)
    return Path(path)


# ---- process_image -------------------------------------------------------------------------------------

def patch_reader(monkeypatch, stack) -> None:
    """process_image reads ``stack`` (and a header that matches it) for any path."""
    c, z, y, x = stack.data.shape
    meta = {"sizes": {"C": c, "Z": z, "Y": y, "X": x}, "spacing": stack.spacing, "spacing_ok": True,
            "channel_names": list(stack.channel_names), "dtype": str(stack.data.dtype), "is_2d": z == 1}
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    monkeypatch.setattr(pipeline, "read_stack", lambda *a, **k: stack)


def read_table(path: Path) -> pd.DataFrame:
    """A results table without the columns that change from run to run."""
    df = pd.read_csv(path, low_memory=False, float_precision="round_trip")
    return df.drop(columns=[c for c in VOLATILE if c in df.columns])


# tests/test_golden.py's tolerance rules: exact on the workstation; in CI GERMQUANT_GOLDEN_RTOL allows another
# CPU's summation order, spots are put in coordinate order on both sides and spot_id is compared as a set
_RTOL_ENV = os.environ.get("GERMQUANT_GOLDEN_RTOL")
_FLOAT_TOLERANCE = ({"check_exact": False, "rtol": float(_RTOL_ENV), "atol": 0.0} if _RTOL_ENV
                    else {"check_exact": True})
_ORDER_KEYS = {"spots": ["z_um", "y_um", "x_um"]}
_ORDINAL_COLUMNS = {"spots": {"spot_id"}}


def _canonical(df: pd.DataFrame, table: str) -> pd.DataFrame:
    keys = _ORDER_KEYS.get(table)
    if not _RTOL_ENV or not keys or not set(keys) <= set(df.columns):
        return df
    k = df[keys].astype(float).round(6)
    return df.loc[k.sort_values(keys, kind="mergesort").index].reset_index(drop=True)


def assert_tables_equal(got: pd.DataFrame, ref: pd.DataFrame, table: str, *, append_only: bool = True) -> None:
    """Every column of ``ref`` equal in ``got`` (exactly, or within the CI tolerance). With ``append_only`` the
    golden rule of tests/test_golden.py: ``got`` may add columns, never rename, remove or reorder them."""
    got, ref = _canonical(got, table), _canonical(ref, table)
    if append_only:
        common = [c for c in got.columns if c in set(ref.columns)]
        assert common == list(ref.columns), f"{table}: existing columns renamed, removed or reordered"
    else:
        assert list(got.columns) == list(ref.columns), f"{table}: columns differ"
    assert len(got) == len(ref), f"{table}: row count {len(got)} vs {len(ref)}"
    for col in ref.columns:
        a, b = got[col], ref[col]
        if _RTOL_ENV and col in _ORDINAL_COLUMNS.get(table, ()):
            assert sorted(a.astype(float)) == sorted(b.astype(float)), f"{table}.{col}: not the same set of ids"
            continue
        if pd.api.types.is_numeric_dtype(b):
            pd.testing.assert_series_equal(a.astype(float), b.astype(float), check_names=False,
                                           obj=f"{table}.{col}", **_FLOAT_TOLERANCE)
        else:
            assert (a.astype("string").fillna("<NA>") == b.astype("string").fillna("<NA>")).all(), \
                f"{table}.{col} differs"
