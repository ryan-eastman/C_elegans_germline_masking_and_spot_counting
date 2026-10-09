"""read_image: one entry point for every input format (io/formats.py, io/read.py, io/bioio_reader.py).

The .nd2 branch is the existing read_stack (same function object). The bioio branch is exercised two
ways: with a fake BioImage (no bioio needed: stride / z_range / spacing rules and the refusals for a
time-lapse, several scenes and a samples dimension) and, when bioio and its OME-TIFF plugin are installed,
on real files written with tifffile. The last test runs the whole pipeline on an OME-TIFF of the synthetic
gonad and checks the tables equal the ones the patched .nd2 reader produces.
"""
from __future__ import annotations

import importlib.util
import sys
from collections import namedtuple
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_end_to_end import CHANNELS, SPACING, _synthetic_stack, patched_reader  # noqa: F401

from germquant import pipeline
from germquant.config import load_config
from germquant.io import Stack, bioio_reader, nd2_reader, read_image, read_image_metadata
from germquant.io import read as read_mod
from germquant.io.bioio_reader import UnsupportedImageError
from germquant.io.formats import SUPPORTED_SUFFIXES, image_format, image_stem

ROOT = Path(__file__).resolve().parents[1]
HAVE_BIOIO = importlib.util.find_spec("bioio") is not None
HAVE_OME_TIFF = HAVE_BIOIO and importlib.util.find_spec("bioio_ome_tiff") is not None
HAVE_TIFFFILE_PLUGIN = HAVE_BIOIO and importlib.util.find_spec("bioio_tifffile") is not None
HAVE_CZI = importlib.util.find_spec("bioio_czi") is not None
needs_ome_tiff = pytest.mark.skipif(not HAVE_OME_TIFF, reason="bioio + bioio-ome-tiff not installed")

NAMES = ["405", "477", "545"]
PPS = namedtuple("PhysicalPixelSizes", "Z Y X")


def _data():
    rng = np.random.default_rng(0)
    return rng.integers(0, 4000, size=(3, 5, 16, 24), dtype=np.uint16)


def _write_ome(path, data, axes="CZYX", names=NAMES, spacing=SPACING):
    import tifffile

    dz, dy, dx = spacing
    tifffile.imwrite(path, data, ome=True, metadata={
        "axes": axes,
        "PhysicalSizeX": dx, "PhysicalSizeXUnit": "µm",
        "PhysicalSizeY": dy, "PhysicalSizeYUnit": "µm",
        "PhysicalSizeZ": dz, "PhysicalSizeZUnit": "µm",
        "Channel": {"Name": list(names)},
    })


# ---- name helpers and dispatch (no reader needed) ----

def test_formats_and_stems():
    assert image_format("a.nd2") == "nd2" and image_format("A.ND2") == "nd2"
    for suf in (".tif", ".tiff", ".ome.tif", ".ome.tiff", ".czi", ".lif", ".zarr", ".ome.zarr"):
        assert image_format(f"x{suf}") == "bioio", suf
    with pytest.raises(ValueError, match="unsupported"):
        image_format("x.png")
    assert ".nd2" in SUPPORTED_SUFFIXES and ".ome.tif" in SUPPORTED_SUFFIXES
    # .nd2 (and unknown suffixes) keep Path.stem exactly: image ids do not move
    for name in ("20251105_n2_nohs_HERM_001.nd2", "20251105_N2_nohs_HERM _001.nd2", "a.b.nd2", "plain", "x.png"):
        assert image_stem(name) == Path(name).stem, name
    assert image_stem("g.ome.tif") == "g" and image_stem("g.ome.tiff") == "g" and image_stem("g.OME.ZARR") == "g"
    assert image_stem("g.tif") == "g" and image_stem("g.zarr") == "g" and image_stem(Path("d") / "g.czi") == "g"


def test_nd2_goes_to_the_existing_read_stack(monkeypatch):
    assert read_mod.read_stack is nd2_reader.read_stack          # the same object, not a copy
    calls = []

    def fake(path, *, xy_stride=1, z_range=None, as_float=False):
        calls.append((Path(path).name, xy_stride, z_range, as_float))
        return "stack"

    monkeypatch.setattr(read_mod, "read_stack", fake)
    assert read_image("g.nd2", xy_stride=2, z_range=(1, 3)) == "stack"
    assert calls == [("g.nd2", 2, (1, 3), False)]
    monkeypatch.setattr(read_mod, "read_nd2_metadata", lambda p: {"from": Path(p).name})
    assert read_image_metadata("g.nd2") == {"from": "g.nd2"}


def test_unknown_suffix_and_missing_file(tmp_path):
    with pytest.raises(ValueError, match="unsupported"):
        read_image(tmp_path / "g.png")
    with pytest.raises(FileNotFoundError):
        read_image(tmp_path / "nope.tif")            # checked before any bioio import


def test_missing_readers_extra_is_a_clear_import_error(tmp_path, monkeypatch):
    p = tmp_path / "g.ome.tif"
    p.write_bytes(b"II*\x00")
    monkeypatch.setitem(sys.modules, "bioio", None)  # `import bioio` now raises ImportError
    with pytest.raises(ImportError, match="readers") as ei:
        read_image(p)
    assert "bioio-ome-tiff" in str(ei.value) and '".[readers]"' in str(ei.value)
    with pytest.raises(ImportError, match="readers"):
        read_image_metadata(p)


# ---- the bioio branch with a fake BioImage (no bioio needed) ----

class _Dims:
    def __init__(self, sizes: dict):
        self.order = "".join(sizes)
        self.__dict__.update(sizes)


class _FakeImage:
    """Just the surface bioio_reader uses: dims, scenes, physical_pixel_sizes, channel_names, dtype,
    reader and get_image_dask_data."""

    def __init__(self, data, *, t=1, scenes=("Image:0",), pps=SPACING, names=NAMES, extra=None):
        import dask.array as da

        self._dask = da.from_array(data, chunks=data.shape)
        self.scenes = tuple(scenes)
        self.channel_names = list(names)
        sizes = {"T": t, "C": data.shape[0], "Z": data.shape[1], "Y": data.shape[2], "X": data.shape[3]}
        sizes.update(extra or {})
        self.dims = _Dims(sizes)
        self.physical_pixel_sizes = PPS(*pps)
        self.dtype = data.dtype
        self.reader = self

    def get_image_dask_data(self, order, **kw):
        assert order == "CZYX" and set(kw) <= {"T", "S"} and all(v == 0 for v in kw.values()), (order, kw)
        return self._dask


def _fake(monkeypatch, img):
    monkeypatch.setattr(bioio_reader, "_open", lambda path: img)


def test_fake_bioio_stack_stride_and_z_range(tmp_path, monkeypatch):
    data = _data()
    p = tmp_path / "g.czi"
    p.write_bytes(b"x")
    _fake(monkeypatch, _FakeImage(data))
    st = read_image(p)
    assert isinstance(st, Stack) and np.array_equal(st.data, data) and st.data.dtype == np.uint16
    assert st.spacing == SPACING and st.spacing_ok and st.channel_names == NAMES
    assert st.path == p and st.downsample_xy == 1 and st.n_channels == 3 and st.shape_zyx == (5, 16, 24)
    sub = read_image(p, xy_stride=2, z_range=(1, 4), as_float=True)
    assert np.array_equal(sub.data, data[:, 1:4, ::2, ::2].astype(np.float32)) and sub.data.dtype == np.float32
    assert sub.spacing == (SPACING[0], SPACING[1] * 2, SPACING[2] * 2) and sub.downsample_xy == 2
    meta = read_image_metadata(p)
    assert meta["sizes"] == {"T": 1, "C": 3, "Z": 5, "Y": 16, "X": 24} and meta["spacing"] == SPACING
    assert meta["spacing_ok"] and meta["channel_names"] == NAMES and meta["is_2d"] is False
    assert meta["dtype"] == "uint16" and meta["n_scenes"] == 1
    # a singleton samples dimension is fine (read at index 0)
    _fake(monkeypatch, _FakeImage(data, extra={"S": 1}))
    assert np.array_equal(read_image(p).data, data)


def test_fake_bioio_unreadable_voxel_size(tmp_path, monkeypatch):
    p = tmp_path / "g.tif"
    p.write_bytes(b"x")
    for pps in ((None, 0.2, 0.2), (0.4, float("nan"), 0.2), (0.4, 0.2, 0.0)):
        _fake(monkeypatch, _FakeImage(_data(), pps=pps))
        st = read_image(p)
        assert st.spacing == (1.0, 1.0, 1.0) and st.spacing_ok is False, pps


def test_fake_bioio_refuses_time_lapse_scenes_and_samples(tmp_path, monkeypatch):
    p = tmp_path / "g.lif"
    p.write_bytes(b"x")
    _fake(monkeypatch, _FakeImage(_data(), t=3))
    with pytest.raises(UnsupportedImageError, match="time-lapse"):
        read_image(p)
    _fake(monkeypatch, _FakeImage(_data(), scenes=("Image:0", "Image:1")))
    with pytest.raises(UnsupportedImageError, match="single scene"):
        read_image(p)
    _fake(monkeypatch, _FakeImage(_data(), extra={"S": 3}))
    with pytest.raises(UnsupportedImageError, match="S dimension of size 3"):
        read_image(p)
    assert issubclass(UnsupportedImageError, ValueError)


# ---- real files through bioio (skipped when the readers extra is not installed) ----

@needs_ome_tiff
def test_ome_tiff_equals_a_stack_built_directly(tmp_path):
    data = _data()
    p = tmp_path / "g.ome.tif"
    _write_ome(p, data)
    ref = Stack(data=data, spacing=SPACING, channel_names=NAMES, path=p)
    st = read_image(p)
    assert np.array_equal(st.data, ref.data) and st.data.dtype == ref.data.dtype
    assert st.spacing == ref.spacing and st.spacing_ok is True
    assert st.channel_names == ref.channel_names and st.path == ref.path and st.downsample_xy == 1
    sub = read_image(p, xy_stride=2, z_range=(1, 4))
    assert np.array_equal(sub.data, data[:, 1:4, ::2, ::2])
    assert sub.spacing == (SPACING[0], SPACING[1] * 2, SPACING[2] * 2) and sub.downsample_xy == 2
    meta = read_image_metadata(p)
    assert meta["channel_names"] == NAMES and meta["spacing"] == SPACING and meta["spacing_ok"]
    assert meta["sizes"]["Z"] == 5 and meta["sizes"]["C"] == 3 and meta["is_2d"] is False
    assert meta["reader"] == "bioio_ome_tiff" and meta["n_scenes"] == 1


@needs_ome_tiff
def test_ome_tiff_time_lapse_and_series_are_refused(tmp_path):
    import tifffile

    data = _data()
    p = tmp_path / "t.ome.tif"
    _write_ome(p, np.repeat(data[None], 2, axis=0), axes="TCZYX")
    with pytest.raises(UnsupportedImageError, match="time-lapse"):
        read_image(p)
    assert read_image_metadata(p)["sizes"]["T"] == 2            # the header still reads
    s = tmp_path / "s.ome.tif"
    with tifffile.TiffWriter(s, ome=True) as tw:                 # two series = two bioio scenes
        tw.write(data, metadata={"axes": "CZYX"})
        tw.write(data, metadata={"axes": "CZYX"})
    with pytest.raises(UnsupportedImageError, match="2 scenes"):
        read_image(s)


@pytest.mark.skipif(not HAVE_TIFFFILE_PLUGIN, reason="bioio + bioio-tifffile not installed")
def test_plain_tiff_reads_with_unreadable_voxel_size(tmp_path):
    import tifffile

    data = _data()
    p = tmp_path / "plain.tif"
    tifffile.imwrite(p, data)
    st = read_image(p)
    assert np.array_equal(st.data, data) and st.spacing_ok is False and st.spacing == (1.0, 1.0, 1.0)
    assert len(st.channel_names) == 3


@pytest.mark.skipif(not HAVE_BIOIO or HAVE_CZI, reason="needs bioio without the CZI plugin")
def test_missing_plugin_names_the_package(tmp_path):
    p = tmp_path / "g.czi"
    p.write_bytes(b"\x00" * 16)
    with pytest.raises(ImportError, match="bioio-czi"):
        read_image(p)


@needs_ome_tiff
def test_process_image_reads_an_ome_tiff_like_the_nd2(tmp_path, patched_reader):  # noqa: F811
    """The whole pipeline on an OME-TIFF of the synthetic gonad: same image id, same tables as the run
    fed by the patched .nd2 reader (only file_path and the timestamp differ)."""
    stack = _synthetic_stack()
    tif = tmp_path / "20251105_n2_nohs_HERM_001.ome.tif"
    _write_ome(tif, stack.data, names=CHANNELS, spacing=SPACING)

    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    res_nd2 = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / "nd2", prov=None)
    res_tif = pipeline.process_image(tif, cfg, tmp_path / "tif", prov=None)
    assert res_tif["image_id"] == res_nd2["image_id"] == "20251105_n2_nohs_HERM_001"
    assert res_tif["n_nuclei"] == res_nd2["n_nuclei"] >= 1
    assert res_tif["stages"]["read"]["status"] == "ran"
    assert res_tif["stages"]["acquisition"]["status"] == "skipped"
    drop = ["file_path", "run_timestamp"]
    for t in ("nuclei", "spots", "image_summary"):
        a = pd.read_csv(tmp_path / "nd2" / f"{res_nd2['image_id']}__{t}.csv", float_precision="round_trip")
        b = pd.read_csv(tmp_path / "tif" / f"{res_tif['image_id']}__{t}.csv", float_precision="round_trip")
        assert b["file_path"].iloc[0] == str(tif) if len(b) else True
        pd.testing.assert_frame_equal(a.drop(columns=drop), b.drop(columns=drop), check_exact=True, obj=t)
