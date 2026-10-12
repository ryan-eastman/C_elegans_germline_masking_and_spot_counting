"""germquant.assemble and `germquant assemble` (design 3.13 and 6.9; package P3, acceptance 4): per-channel TIFFs and
multi-channel TIFFs with a declared page order are stacked into one OME-TIFF; the round trip gives the declared
axes, the channel names, the voxel size, the crop and a pixels_sha256 that pins the pixels, not the file (tifffile
writes a fresh UUID into every OME-XML header). A declaration that contradicts a file's own ImageJ or OME record of
its channels and planes is refused unless forced, and the output never replaces an input."""
from __future__ import annotations

import hashlib

import numpy as np
import pytest
import tifffile

from germquant import assemble as A
from germquant import cli
from germquant.io.ome_header import read_tiff_header

SP = (0.24, 0.072, 0.072)


def _stack(seed, shape=(5, 40, 48)):
    rng = np.random.default_rng(seed)
    return rng.integers(100, 4000, size=shape, dtype=np.uint16)


def test_per_channel_files_round_trip_with_names_voxel_size_and_crop(tmp_path):
    bright, cy3 = _stack(1), _stack(2)
    tifffile.imwrite(tmp_path / "s_Bright.tif", bright, imagej=True)        # an ImageJ hyperstack (axes ZYX)
    tifffile.imwrite(tmp_path / "s_Cy3.tif", cy3)                           # a plain multipage TIFF
    out = tmp_path / "s.ome.tif"
    rec = A.assemble([tmp_path / "s_Bright.tif", tmp_path / "s_Cy3.tif"], ["Bright", "Cy3"], SP, out,
                     crop={"z": [1, 4], "y": [5, 35], "x": [8, 40]})
    want = np.stack([bright, cy3])[:, 1:4, 5:35, 8:40]
    back = A.read_ome_tiff(out)
    assert back["data"].shape == (2, 3, 30, 32) and np.array_equal(back["data"], want)
    assert back["channel_names"] == ["Bright", "Cy3"] and back["spacing_zyx_um"] == pytest.approx(list(SP))
    assert rec["shape"] == [2, 3, 30, 32] and rec["axes"] == "CZYX" and rec["dtype"] == "uint16"
    assert rec["crop"] == {"z": [1, 4], "y": [5, 35], "x": [8, 40]} and rec["declared_axes"] == "ZYX"
    assert rec["pixels_sha256"] == A.pixels_sha256(want) == back["pixels_sha256"] == A.file_pixels_sha256(out)
    hdr = read_tiff_header(out)                                             # what pack check and bioio read
    assert hdr["spacing_ok"] and hdr["spacing"] == pytest.approx(SP) and hdr["channel_names"] == ["Bright", "Cy3"]
    with tifffile.TiffFile(out) as tf:
        assert tf.is_ome and tf.series[0].axes == "CZYX"


def test_pixels_sha256_pins_pixels_not_files(tmp_path):
    a = _stack(3)
    tifffile.imwrite(tmp_path / "a.tif", a)
    r1 = A.assemble([tmp_path / "a.tif"], ["DAPI"], SP, tmp_path / "one.ome.tif")
    r2 = A.assemble([tmp_path / "a.tif"], ["DAPI"], SP, tmp_path / "two.ome.tif")
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()  # noqa: E731
    assert digest(tmp_path / "one.ome.tif") != digest(tmp_path / "two.ome.tif")     # a fresh UUID per file
    assert r1["pixels_sha256"] == r2["pixels_sha256"]
    # the same bytes with another shape or data type are other pixels
    flat = a[None]
    assert A.pixels_sha256(flat) != A.pixels_sha256(flat.reshape(1, 5, 48, 40))
    assert A.pixels_sha256(flat) != A.pixels_sha256(flat.view(np.int16))
    assert A.pixels_sha256(flat.astype(">u2")) == A.pixels_sha256(flat)              # byte order does not count


def test_a_single_plane_file_is_one_plane(tmp_path):
    plane = _stack(4, (32, 36))
    tifffile.imwrite(tmp_path / "p.tif", plane)
    rec = A.assemble([tmp_path / "p.tif"], ["phase"], SP, tmp_path / "p.ome.tif")
    assert rec["shape"] == [1, 1, 32, 36]
    assert np.array_equal(A.read_ome_tiff(tmp_path / "p.ome.tif")["data"][0, 0], plane)


def test_channels_as_samples_like_deepfoci(tmp_path):
    """DeepFoci's data TIFFs: one page per plane, three uint16 samples per pixel (RGB, contiguous). Declared ZYXC,
    the samples become the channels in their order."""
    zyxs = np.stack([_stack(5), _stack(6), _stack(7)], axis=-1)                 # (Z, Y, X, 3)
    tifffile.imwrite(tmp_path / "data_001.tif", zyxs, photometric="rgb", compression="zlib")
    rec = A.assemble([tmp_path / "data_001.tif"], ["53BP1", "gH2AX", "DAPI"], (0.3, 0.0647, 0.0647),
                     tmp_path / "d.ome.tif", axes="ZYXC", crop={"y": [0, 20]})
    back = A.read_ome_tiff(tmp_path / "d.ome.tif")
    assert back["data"].shape == (3, 5, 20, 48) and rec["channel_names"] == ["53BP1", "gH2AX", "DAPI"]
    for c in range(3):
        assert np.array_equal(back["data"][c], zyxs[:, :20, :, c])
    with pytest.raises(A.AssembleError, match="must end with the axes of its pages, YXC"):
        A.assemble([tmp_path / "data_001.tif"], ["a", "b", "c"], SP, tmp_path / "x.ome.tif", axes="ZCYX")


def test_channels_alternating_page_by_page(tmp_path):
    """An ImageJ-style stack whose pages alternate channels within each plane (pages z0c0, z0c1, z1c0, ...):
    declared ZCYX, read without trusting any metadata (here a plain multipage TIFF of grey pages)."""
    zcyx = np.stack([np.stack([_stack(10 + z, (24, 28)), _stack(20 + z, (24, 28))]) for z in range(4)])
    tifffile.imwrite(tmp_path / "hyper.tif", zcyx.reshape(8, 24, 28))
    A.assemble([tmp_path / "hyper.tif"], ["DAPI", "foci"], SP, tmp_path / "h.ome.tif", axes="ZCYX")
    back = A.read_ome_tiff(tmp_path / "h.ome.tif")["data"]
    assert np.array_equal(back, np.transpose(zcyx, (1, 0, 2, 3)))
    with pytest.raises(A.AssembleError, match="do not split into 3 channels"):
        A.assemble([tmp_path / "hyper.tif"], ["a", "b", "c"], SP, tmp_path / "x.ome.tif", axes="ZCYX")


@pytest.mark.parametrize("kwargs, needle", [
    ({"names": ["only_one"]}, "2 channels, 1 names"),
    ({"names": ["a", "a"]}, "unique"),
    ({"axes": "ZYXQ"}, "letters C, Z, Y and X"),
    ({"crop": {"y": [0, 99]}}, "outside the 40 pixels"),
    ({"crop": {"c": [0, 1]}}, "only z, y and x"),
    ({"spacing": (0.2, 0.0, 0.1)}, "three positive numbers"),
])
def test_bad_settings_are_refused(tmp_path, kwargs, needle):
    for name, seed in (("a.tif", 1), ("b.tif", 2)):
        tifffile.imwrite(tmp_path / name, _stack(seed))
    args = {"names": ["a", "b"], "spacing": SP, "axes": "ZYX", "crop": None, **kwargs}
    with pytest.raises(A.AssembleError, match=needle):
        A.assemble([tmp_path / "a.tif", tmp_path / "b.tif"], args["names"], args["spacing"], tmp_path / "o.ome.tif",
                   axes=args["axes"], crop=args["crop"])
    assert not (tmp_path / "o.ome.tif").exists()


def test_inputs_that_do_not_match_are_refused(tmp_path):
    tifffile.imwrite(tmp_path / "a.tif", _stack(1))
    tifffile.imwrite(tmp_path / "b.tif", _stack(2, (6, 40, 48)))
    tifffile.imwrite(tmp_path / "c.tif", _stack(3).astype(np.float32))
    with pytest.raises(A.AssembleError, match="differs from a.tif"):
        A.assemble([tmp_path / "a.tif", tmp_path / "b.tif"], ["a", "b"], SP, tmp_path / "o.ome.tif")
    with pytest.raises(A.AssembleError, match="data type float32 differs"):
        A.assemble([tmp_path / "a.tif", tmp_path / "c.tif"], ["a", "c"], SP, tmp_path / "o.ome.tif")
    with pytest.raises(FileNotFoundError):
        A.assemble([tmp_path / "missing.tif"], ["a"], SP, tmp_path / "o.ome.tif")


def test_a_declaration_that_contradicts_the_files_own_axes_is_refused(tmp_path):
    """An ImageJ hyperstack records its axes (ZCYX: the channels of a plane are neighbouring pages). Declaring CZYX,
    or ZYX (no channels), would mix the channels into the planes without a sign: refused, naming the file's axes;
    force reads the pages as declared. A file without such a record is read as declared, as always."""
    zcyx = np.stack([np.stack([_stack(30 + z, (24, 28)), _stack(40 + z, (24, 28))]) for z in range(4)])
    tifffile.imwrite(tmp_path / "hyper.tif", zcyx, imagej=True, metadata={"axes": "ZCYX"})
    with tifffile.TiffFile(tmp_path / "hyper.tif") as tf:
        assert tf.is_imagej and tf.series[0].axes == "ZCYX"
    with pytest.raises(A.AssembleError, match="ImageJ metadata gives the axes ZCYX .*run them the other way"):
        A.assemble([tmp_path / "hyper.tif"], ["DAPI", "foci"], SP, tmp_path / "a.ome.tif", axes="CZYX")
    with pytest.raises(A.AssembleError, match="do not split the pages by both"):
        A.assemble([tmp_path / "hyper.tif"], ["DAPI"], SP, tmp_path / "a.ome.tif", axes="ZYX")
    assert not (tmp_path / "a.ome.tif").exists()
    A.assemble([tmp_path / "hyper.tif"], ["DAPI", "foci"], SP, tmp_path / "ok.ome.tif", axes="ZCYX")
    assert np.array_equal(A.read_ome_tiff(tmp_path / "ok.ome.tif")["data"], np.transpose(zcyx, (1, 0, 2, 3)))
    forced = A.assemble([tmp_path / "hyper.tif"], ["DAPI", "foci"], SP, tmp_path / "f.ome.tif", axes="CZYX", force=True)
    assert forced["shape"] == [2, 4, 24, 28]                                 # read as declared, as asked
    assert cli.main(["assemble", "--channels", str(tmp_path / "hyper.tif"), "--names", "DAPI", "foci", "--spacing",
                     "1", "1", "1", "--axes", "CZYX", "--out", str(tmp_path / "c.ome.tif")]) == 1
    assert cli.main(["assemble", "--channels", str(tmp_path / "hyper.tif"), "--names", "DAPI", "foci", "--spacing",
                     "1", "1", "1", "--axes", "CZYX", "--force", "--out", str(tmp_path / "c.ome.tif")]) == 0
    plain = zcyx.reshape(8, 24, 28)                                          # the same pages, no metadata
    tifffile.imwrite(tmp_path / "plain.tif", plain)
    A.assemble([tmp_path / "plain.tif"], ["DAPI", "foci"], SP, tmp_path / "p.ome.tif", axes="CZYX")


def test_the_output_never_replaces_an_input(tmp_path):
    a = _stack(1)
    tifffile.imwrite(tmp_path / "S.tif", a)
    before = (tmp_path / "S.tif").read_bytes()
    with pytest.raises(A.AssembleError, match="is one of the inputs"):
        A.assemble([tmp_path / "S.tif"], ["DAPI"], SP, tmp_path / "S.tif")
    with pytest.raises(A.AssembleError, match="is one of the inputs"):
        A.assemble([tmp_path / "S.tif"], ["DAPI"], SP, tmp_path / "sub" / ".." / "S.tif")
    assert (tmp_path / "S.tif").read_bytes() == before


def test_germquant_assemble_command(tmp_path, capsys):
    tifffile.imwrite(tmp_path / "x_Bright.tif", _stack(1))
    tifffile.imwrite(tmp_path / "x_Cy3.tif", _stack(2))
    out = tmp_path / "x.ome.tif"
    rc = cli.main(["assemble", "--channels", str(tmp_path / "x_Bright.tif"), str(tmp_path / "x_Cy3.tif"),
                   "--names", "Bright", "Cy3", "--spacing", "0.24", "0.072", "0.072", "--out", str(out),
                   "--crop-x", "0", "32"])
    assert rc == 0 and "pixels_sha256" in capsys.readouterr().out
    back = A.read_ome_tiff(out)
    assert back["data"].shape == (2, 5, 40, 32) and back["channel_names"] == ["Bright", "Cy3"]
    rc = cli.main(["assemble", "--channels", str(tmp_path / "x_Bright.tif"), "--names", "a", "b",
                   "--spacing", "1", "1", "1", "--out", str(tmp_path / "bad.ome.tif")])
    assert rc == 1 and "ERROR" in capsys.readouterr().err
