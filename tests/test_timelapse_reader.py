"""The T-aware reader (docs/TIMELAPSE.md, "Reading one timepoint").

read_stack(t=...) and read_image(t=...) select one timepoint of a time series and read only its planes;
read_nd2_metadata reports n_timepoints, period_s and t_seconds; process_image(t=...) suffixes the image id
and records the timepoint in the run geometry and the markers. A fixed image (t=None) is read, named and
recorded exactly as before, and a timepoint run on the same stack writes the same numbers.

The nd2 reader is exercised against a fake ND2File: the surface read_stack and read_nd2_metadata use, over an
in-memory array, with nd2's one-dask-chunk-per-(T, Z)-frame layout and a record of every frame read. The
bioio reader is exercised against a fake BioImage and, when bioio and bioio-ome-tiff are installed, against
a real OME-TIFF with a T axis written by tifffile.
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace as NS

import dask.array as da
import nd2
import numpy as np
import pandas as pd
import pytest
from test_end_to_end import CHANNELS, SPACING, _synthetic_stack
from test_read_image import NAMES, _write_ome, needs_ome_tiff

from germquant import api, cli, pipeline, provenance
from germquant.config import load_config
from germquant.io import Stack, bioio_reader, read_image, read_image_metadata
from germquant.io import read as read_mod
from germquant.io.bioio_reader import UnsupportedImageError
from germquant.io.nd2_reader import check_timepoint, read_nd2_metadata, read_stack

ROOT = Path(__file__).resolve().parents[1]
VOXEL = (0.2, 0.108, 0.108)        # (dz, dy, dx) of the October 2026 GFP::SYP-2 series
PERIOD_MS = 300000.0               # its TimeLoop periodMs (5 minutes)
T0_MS = 358.0                      # the recorded time stamp of its first plane
SERIES_NAMES = ["Brightfield", "EGFP"]
NAME = "20251105_n2_nohs_HERM_001.nd2"
IMAGE_ID = "20251105_n2_nohs_HERM_001"
TABLES = ("nuclei", "spots", "granules", "coloc", "image_summary")


class _FakeND2File:
    """What read_stack and read_nd2_metadata use of an nd2.ND2File, over an in-memory (T, Z, C, Y, X) array
    (or (Z, C, Y, X) for a fixed image), chunked as nd2 chunks it (one dask block per (T, Z) frame holding
    (C, Y, X)); every frame the dask array reads is recorded in ``frames_read`` as (t, z)."""

    def __init__(self, data: np.ndarray, *, has_t: bool, names=SERIES_NAMES, voxel=VOXEL):
        self._data = data
        self.has_t = has_t
        if has_t:
            T, Z, C, Y, X = data.shape
            self.sizes = {"T": T, "Z": Z, "C": C, "Y": Y, "X": X}
        else:
            Z, C, Y, X = data.shape
            self.sizes = {"Z": Z, "C": C, "Y": Y, "X": X}
        self.dtype = data.dtype
        self.metadata = NS(channels=[NS(channel=NS(name=n)) for n in names])
        self._voxel = voxel
        self.frames_read: list[tuple[int, int]] = []

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def voxel_size(self):
        dz, dy, dx = self._voxel
        return NS(x=dx, y=dy, z=dz)

    @property
    def experiment(self):
        z_loop = NS(type="ZStackLoop", count=self.sizes["Z"], parameters=NS(stepUm=self._voxel[0]))
        if self.has_t:
            return [NS(type="TimeLoop", count=self.sizes["T"], parameters=NS(periodMs=PERIOD_MS)), z_loop]
        return [NS(type="XYPosLoop", count=1, parameters=NS()), z_loop]

    @property
    def loop_indices(self):
        Z = self.sizes["Z"]
        if self.has_t:
            return [{"T": t, "Z": z} for t in range(self.sizes["T"]) for z in range(Z)]
        return [{"P": 0, "Z": z} for z in range(Z)]

    def frame_metadata(self, index: int):
        li = self.loop_indices[index]
        ms = T0_MS + li.get("T", 0) * PERIOD_MS + li["Z"] * 50.0      # planes 50 ms apart within a timepoint
        return NS(channels=[NS(time=NS(relativeTimeMs=ms))])

    def to_dask(self):
        if self.has_t:
            T, Z, C, Y, X = self._data.shape
            chunks = ((1,) * T, (1,) * Z, (C,), (Y,), (X,))
        else:
            Z, C, Y, X = self._data.shape
            chunks = ((1,) * Z, (C,), (Y,), (X,))
        meta = np.empty((0,) * self._data.ndim, dtype=self.dtype)
        return da.map_blocks(self._block, chunks=chunks, dtype=self.dtype, meta=meta)

    def _block(self, block_id=None):
        if block_id is None or isinstance(block_id, np.ndarray):
            return None
        if self.has_t:
            t, z = int(block_id[0]), int(block_id[1])
            self.frames_read.append((t, z))
            return self._data[t, z][None, None]
        z = int(block_id[0])
        self.frames_read.append((0, z))
        return self._data[z][None]


def _series(T=4, Z=3, C=2, Y=8, X=10, seed=0) -> np.ndarray:
    return np.random.default_rng(seed).integers(0, 4000, size=(T, Z, C, Y, X), dtype=np.uint16)


def _czyx(tzcyx: np.ndarray, t: int) -> np.ndarray:
    """The (C, Z, Y, X) stack of timepoint t of a (T, Z, C, Y, X) array."""
    return np.ascontiguousarray(np.moveaxis(tzcyx[t], 1, 0))


@pytest.fixture
def fake_nd2(monkeypatch):
    """Install a fake ND2File over ``data`` and return it (``frames_read`` says which planes a read touched)."""
    def install(data, *, has_t=True, **kw):
        fake = _FakeND2File(data, has_t=has_t, **kw)
        monkeypatch.setattr(nd2, "ND2File", lambda path: fake)
        return fake
    return install


# ---- the nd2 reader ----

def test_t_selects_that_timepoint_and_reads_only_its_planes(fake_nd2):
    data = _series()
    fake = fake_nd2(data)
    st = read_stack("series.nd2", t=2)
    assert isinstance(st, Stack) and np.array_equal(st.data, _czyx(data, 2)) and st.data.dtype == np.uint16
    assert st.spacing == VOXEL and st.spacing_ok and st.channel_names == SERIES_NAMES
    assert st.t_index == 2 and st.downsample_xy == 1 and st.path == Path("series.nd2")
    assert sorted(fake.frames_read) == [(2, 0), (2, 1), (2, 2)]          # that timepoint's planes, nothing else

    # the last timepoint with the usual sub-sampling: still only the planes asked for
    fake = fake_nd2(data)
    sub = read_stack("series.nd2", t=3, xy_stride=2, z_range=(1, 3), as_float=True)
    assert np.array_equal(sub.data, _czyx(data, 3)[:, 1:3, ::2, ::2].astype(np.float32))
    assert sub.spacing == (VOXEL[0], VOXEL[1] * 2, VOXEL[2] * 2) and sub.downsample_xy == 2 and sub.t_index == 3
    assert sorted(fake.frames_read) == [(3, 1), (3, 2)]
    assert read_stack("series.nd2", t=np.int64(1)).t_index == 1        # numpy integers are integers


def test_t_none_is_t0_for_a_series(fake_nd2):
    data = _series()
    fake = fake_nd2(data)
    first = read_stack("series.nd2")
    assert np.array_equal(first.data, _czyx(data, 0)) and first.t_index is None
    assert sorted(fake.frames_read) == [(0, 0), (0, 1), (0, 2)]
    fake = fake_nd2(data)
    zero = read_stack("series.nd2", t=0)
    assert np.array_equal(zero.data, first.data) and zero.t_index == 0
    assert sorted(fake.frames_read) == [(0, 0), (0, 1), (0, 2)]
    for field in ("spacing", "spacing_ok", "channel_names", "downsample_xy", "path"):
        assert getattr(zero, field) == getattr(first, field), field


def test_out_of_range_t_is_a_value_error_before_any_plane_is_read(fake_nd2):
    fake = fake_nd2(_series())
    for bad in (4, -1, 99):
        with pytest.raises(ValueError, match="out of range"):
            read_stack("series.nd2", t=bad)
    for bad in (1.5, "2", object()):
        with pytest.raises(ValueError, match="integer"):
            read_stack("series.nd2", t=bad)
    assert fake.frames_read == []
    with pytest.raises(ValueError, match="out of range"):
        read_nd2_metadata("series.nd2", t=4)
    assert read_nd2_metadata("series.nd2", t=3)["n_timepoints"] == 4


def test_series_metadata_reports_timepoints_period_and_stamps(fake_nd2):
    fake = fake_nd2(_series())
    m = read_nd2_metadata("series.nd2")
    assert m["sizes"] == {"T": 4, "Z": 3, "C": 2, "Y": 8, "X": 10} and m["spacing"] == VOXEL and m["spacing_ok"]
    assert m["channel_names"] == SERIES_NAMES and m["dtype"] == "uint16" and m["is_2d"] is False
    assert m["n_timepoints"] == 4 and m["period_s"] == 300.0
    assert m["t_seconds"] == pytest.approx([0.358, 300.358, 600.358, 900.358])   # the first plane of each timepoint
    assert fake.frames_read == []                                                # a header read loads no pixels


def test_a_fixed_image_is_unchanged(fake_nd2):
    data = _series(T=1)[0]                      # (Z, C, Y, X): no T axis, as every fixed gonad stack
    fake = fake_nd2(data, has_t=False, names=["405", "477"])
    st = read_stack("fixed.nd2")
    expected = np.ascontiguousarray(np.moveaxis(data, 1, 0))
    assert np.array_equal(st.data, expected) and st.t_index is None and st.channel_names == ["405", "477"]
    assert sorted(fake.frames_read) == [(0, 0), (0, 1), (0, 2)]
    zero = read_stack("fixed.nd2", t=0)                            # t=0 is allowed and identical
    assert np.array_equal(zero.data, expected) and zero.t_index == 0
    with pytest.raises(ValueError, match="no T axis"):
        read_stack("fixed.nd2", t=1)
    m = read_nd2_metadata("fixed.nd2")
    assert m["sizes"] == {"Z": 3, "C": 2, "Y": 8, "X": 10} and m["is_2d"] is False
    assert m["n_timepoints"] == 1 and m["period_s"] is None and m["t_seconds"] is None
    assert set(m) >= {"sizes", "spacing", "spacing_ok", "channel_names", "dtype", "is_2d"}   # the old keys stay
    with pytest.raises(ValueError, match="no T axis"):
        read_nd2_metadata("fixed.nd2", t=1)


def test_check_timepoint():
    assert check_timepoint(0, 1, "a.nd2") == 0 and check_timepoint(24, 25, "a.nd2") == 24
    with pytest.raises(ValueError, match="0 to 24"):
        check_timepoint(25, 25, "a.nd2")
    with pytest.raises(ValueError, match="no T axis"):
        check_timepoint(1, 1, "a.nd2")
    with pytest.raises(ValueError, match="integer"):
        check_timepoint(2.0, 25, "a.nd2")


# ---- the dispatcher ----

def test_read_image_passes_t_on_only_when_given(monkeypatch):
    calls = []

    def fake_stack(path, **kw):
        calls.append(("stack", Path(path).name, kw))
        return "stack"

    def fake_meta(path, **kw):
        calls.append(("meta", Path(path).name, kw))
        return {"n_timepoints": 25}

    monkeypatch.setattr(read_mod, "read_stack", fake_stack)
    monkeypatch.setattr(read_mod, "read_nd2_metadata", fake_meta)
    monkeypatch.setattr(read_mod, "read_bioio", fake_stack)
    monkeypatch.setattr(read_mod, "read_bioio_metadata", fake_meta)
    read_image("g.nd2", xy_stride=2)
    read_image("g.nd2", t=7)
    read_image_metadata("g.nd2")
    read_image_metadata("g.nd2", t=7)
    read_image("g.ome.tif", t=3)
    read_image_metadata("g.czi", t=1)
    assert calls == [
        ("stack", "g.nd2", {"xy_stride": 2, "z_range": None, "as_float": False}),     # today's call: no t key
        ("stack", "g.nd2", {"xy_stride": 1, "z_range": None, "as_float": False, "t": 7}),
        ("meta", "g.nd2", {}),
        ("meta", "g.nd2", {"t": 7}),
        ("stack", "g.ome.tif", {"xy_stride": 1, "z_range": None, "as_float": False, "t": 3}),
        ("meta", "g.czi", {"t": 1}),
    ]


# ---- the bioio reader ----

class _FakeTimeLapseImage:
    """The BioImage surface bioio_reader uses, over a (T, C, Z, Y, X) array; get_image_dask_data honours T."""

    def __init__(self, data, *, names=NAMES, pps=SPACING):
        self._data = data
        self.scenes = ("Image:0",)
        self.channel_names = list(names)
        T, C, Z, Y, X = data.shape
        self.dims = NS(order="TCZYX", T=T, C=C, Z=Z, Y=Y, X=X)
        self.physical_pixel_sizes = NS(Z=pps[0], Y=pps[1], X=pps[2])
        self.dtype = data.dtype
        self.reader = self
        self.asked: list[dict] = []

    def get_image_dask_data(self, order, **kw):
        assert order == "CZYX" and set(kw) == {"T"}, (order, kw)
        self.asked.append(dict(kw))
        return da.from_array(self._data[kw["T"]], chunks=self._data.shape[1:])


def test_bioio_time_lapse_reads_one_timepoint_with_t(tmp_path, monkeypatch):
    data = np.random.default_rng(1).integers(0, 4000, size=(3, 2, 4, 6, 8), dtype=np.uint16)
    img = _FakeTimeLapseImage(data)
    monkeypatch.setattr(bioio_reader, "_open", lambda path: img)
    p = tmp_path / "g.czi"
    p.write_bytes(b"x")
    with pytest.raises(UnsupportedImageError, match="time-lapse") as ei:
        read_image(p)
    assert "t=<index>" in str(ei.value) and "0 to 2" in str(ei.value) and img.asked == []
    st = read_image(p, t=1)
    assert np.array_equal(st.data, data[1]) and st.t_index == 1 and img.asked == [{"T": 1}]
    assert st.spacing == SPACING and st.channel_names == NAMES
    sub = read_image(p, t=2, xy_stride=2, z_range=(1, 3))
    assert np.array_equal(sub.data, data[2][:, 1:3, ::2, ::2])
    assert sub.spacing == (SPACING[0], SPACING[1] * 2, SPACING[2] * 2) and sub.t_index == 2
    for bad in (3, -1):
        with pytest.raises(ValueError, match="out of range"):
            read_image(p, t=bad)
    m = read_image_metadata(p)
    assert m["sizes"] == {"T": 3, "C": 2, "Z": 4, "Y": 6, "X": 8} and m["n_timepoints"] == 3
    assert m["period_s"] is None and m["t_seconds"] is None and m["n_scenes"] == 1
    with pytest.raises(ValueError, match="out of range"):
        read_image_metadata(p, t=3)
    assert read_image_metadata(p, t=2)["n_timepoints"] == 3

    # a single time point: read as before, t=0 allowed, anything else refused
    fixed = _FakeTimeLapseImage(data[:1])
    monkeypatch.setattr(bioio_reader, "_open", lambda path: fixed)
    assert np.array_equal(read_image(p).data, data[0]) and read_image(p).t_index is None
    assert np.array_equal(read_image(p, t=0).data, data[0]) and read_image(p, t=0).t_index == 0
    with pytest.raises(ValueError, match="no T axis"):
        read_image(p, t=1)
    assert read_image_metadata(p)["n_timepoints"] == 1


@needs_ome_tiff
def test_ome_tiff_with_a_t_axis_reads_one_timepoint(tmp_path):
    data = np.random.default_rng(2).integers(0, 4000, size=(3, 3, 5, 16, 24), dtype=np.uint16)   # (T, C, Z, Y, X)
    p = tmp_path / "t.ome.tif"
    _write_ome(p, data, axes="TCZYX")
    with pytest.raises(UnsupportedImageError, match="time-lapse"):
        read_image(p)
    for t in range(3):
        st = read_image(p, t=t)
        assert np.array_equal(st.data, data[t]) and st.t_index == t and st.spacing == SPACING, t
    sub = read_image(p, t=2, xy_stride=2, z_range=(1, 4))
    assert np.array_equal(sub.data, data[2][:, 1:4, ::2, ::2])
    with pytest.raises(ValueError, match="out of range"):
        read_image(p, t=3)
    m = read_image_metadata(p)
    assert m["sizes"]["T"] == 3 and m["n_timepoints"] == 3 and m["channel_names"] == NAMES
    assert m["period_s"] is None and m["t_seconds"] is None
    # a fixed OME-TIFF reports one timepoint and takes t=0 only
    f = tmp_path / "f.ome.tif"
    _write_ome(f, data[0])
    assert read_image_metadata(f)["n_timepoints"] == 1
    assert np.array_equal(read_image(f, t=0).data, data[0])
    with pytest.raises(ValueError, match="no T axis"):
        read_image(f, t=1)


# ---- process_image and the API ----

@pytest.fixture
def recording_reader(monkeypatch):
    """The synthetic 3-channel gonad through the pipeline's reader seam, as test_end_to_end's patched_reader,
    with a series header (25 timepoints, 300 s apart) and a record of the keywords the pipeline passed."""
    stack = _synthetic_stack()
    meta = {"sizes": {"T": 25, "Z": 20, "C": 3, "Y": 80, "X": 320}, "spacing": SPACING, "spacing_ok": True,
            "channel_names": CHANNELS, "dtype": "float32", "is_2d": False,
            "n_timepoints": 25, "period_s": 300.0, "t_seconds": [0.358 + 300.0 * i for i in range(25)]}
    calls = []

    def fake_meta(path, **kw):
        calls.append(("meta", kw))
        return meta

    def fake_stack(path, **kw):
        calls.append(("stack", kw))
        return stack

    monkeypatch.setattr(pipeline, "read_nd2_metadata", fake_meta)
    monkeypatch.setattr(pipeline, "read_stack", fake_stack)
    return calls


def _classical():
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    return cfg


def test_process_image_with_t_suffixes_the_id_and_records_the_timepoint(tmp_path, recording_reader):
    out = tmp_path / "t7"
    res = pipeline.process_image(NAME, _classical(), out, prov=None, t=7)
    iid = f"{IMAGE_ID}_t007"
    assert res["image_id"] == iid and res["n_nuclei"] >= 1
    assert recording_reader == [("meta", {"t": 7}), ("stack", {"xy_stride": 1, "z_range": None, "t": 7})]
    # every per-image file carries the suffixed id, none the bare stem
    names = sorted(f.name for f in out.iterdir())
    for suffix in ("__nuclei.csv", "__image_summary.csv", "__stages.json", "__done.json", "__montage.png",
                   "__nuclei_labels.tif"):
        assert f"{iid}{suffix}" in names, suffix
    assert not [n for n in names if n.startswith(IMAGE_ID) and not n.startswith(iid)], names
    # the design columns come from the raw stem; no new table column
    nuclei = pd.read_csv(out / f"{iid}__nuclei.csv")
    assert (nuclei["image_id"] == iid).all()
    assert nuclei["sex"].iloc[0] == "herm" and nuclei["treatment"].iloc[0] == "control"
    assert nuclei["genotype"].iloc[0] == "N2" and str(nuclei["acquisition_date"].iloc[0]) == "20251105"
    assert "t_index" not in nuclei.columns and "t_seconds" not in nuclei.columns
    # the timepoint is in the stage record, the manifest and the completion marker
    stages = json.loads((out / f"{iid}__stages.json").read_text(encoding="utf-8"))
    assert stages["image_id"] == iid and stages["t_index"] == 7 and stages["t_seconds"] == pytest.approx(2100.358)
    manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["t_index"] == 7 and manifest["t_seconds"] == pytest.approx(2100.358)
    marker = provenance.read_done_marker(out, iid)
    assert marker["image_id"] == iid and marker["t_index"] == 7 and marker["t_seconds"] == pytest.approx(2100.358)
    assert marker["xy_stride"] == 1 and marker["z_range"] is None and marker["failed_stages"] == []
    assert provenance.read_done_marker(out, IMAGE_ID) is None            # nothing under the bare stem


def test_timepoint_run_equals_the_fixed_run_and_the_fixed_run_records_nothing_new(tmp_path, recording_reader):
    cfg = _classical()
    fixed, series = tmp_path / "fixed", tmp_path / "t3"
    res0 = pipeline.process_image(NAME, cfg, fixed, prov=None)
    assert recording_reader == [("meta", {}), ("stack", {"xy_stride": 1, "z_range": None})]   # today's calls
    res3 = pipeline.process_image(NAME, cfg, series, prov=None, t=3)
    assert res0["image_id"] == IMAGE_ID and res3["image_id"] == f"{IMAGE_ID}_t003"
    assert res0["n_nuclei"] == res3["n_nuclei"] and res0["qc_flags"] == res3["qc_flags"]
    # the stages ran unchanged: every table is identical apart from the image id and the time stamp
    for t in TABLES:
        a = pd.read_csv(fixed / f"{IMAGE_ID}__{t}.csv", float_precision="round_trip")
        b = pd.read_csv(series / f"{IMAGE_ID}_t003__{t}.csv", float_precision="round_trip")
        assert list(a.columns) == list(b.columns), t
        if len(a):
            assert (a["image_id"] == IMAGE_ID).all() and (b["image_id"] == f"{IMAGE_ID}_t003").all()
        drop = [c for c in ("image_id", "run_timestamp") if c in a.columns]
        pd.testing.assert_frame_equal(a.drop(columns=drop), b.drop(columns=drop), check_exact=True, obj=t)
    # the fixed image's stage record, manifest and marker carry no timepoint fields at all
    stages0 = json.loads((fixed / f"{IMAGE_ID}__stages.json").read_text(encoding="utf-8"))
    manifest0 = json.loads((fixed / "run_manifest.json").read_text(encoding="utf-8"))
    marker0 = provenance.read_done_marker(fixed, IMAGE_ID)
    for rec in (stages0, manifest0, marker0):
        assert "t_index" not in rec and "t_seconds" not in rec
    assert list(marker0) == ["image_id", "config_hash", "enabled_stages", "stage_hashes", "xy_stride", "z_range",
                             "failed_stages", "finished"]


def test_api_run_passes_t_on_only_when_given(tmp_path, monkeypatch):
    seen = []

    def fake_process_image(path, cfg, out_dir, *, xy_stride=1, z_range=None, prov=None, **kw):
        seen.append(kw)
        return {"image_id": "x", "n_nuclei": 0, "qc_pass": True, "qc_flags": [], "out_dir": str(out_dir),
                "tables": {}, "stages": {}}

    monkeypatch.setattr(api, "process_image", fake_process_image)
    cfg = load_config(ROOT / "config" / "config.yaml")
    api.run("g.nd2", config=cfg, out=tmp_path / "a")
    api.run("g.nd2", config=cfg, out=tmp_path / "b", t=7)
    api.run_loaded("g.nd2", cfg, tmp_path / "c", t=0)
    assert seen == [{}, {"t": 7}, {"t": 0}]


def test_info_says_how_many_timepoints_a_series_has(monkeypatch, capsys):
    meta = {"sizes": {"T": 25, "Z": 76, "C": 2, "Y": 2720, "X": 2720}, "dtype": "uint16", "is_2d": False,
            "spacing": VOXEL, "spacing_ok": True, "channel_names": SERIES_NAMES,
            "n_timepoints": 25, "period_s": 300.0, "t_seconds": None}
    monkeypatch.setattr(cli, "read_image_metadata", lambda path: meta)
    monkeypatch.setattr(cli, "_print_acquisition", lambda path: None)
    assert cli.main(["info", "series.nd2"]) == 0
    assert "timepoints : 25 every 300 s" in capsys.readouterr().out
    fixed = {**meta, "sizes": {"Z": 71, "C": 3, "Y": 2600, "X": 2600}, "n_timepoints": 1, "period_s": None}
    monkeypatch.setattr(cli, "read_image_metadata", lambda path: fixed)
    assert cli.main(["info", "fixed.nd2"]) == 0
    assert "timepoints" not in capsys.readouterr().out
