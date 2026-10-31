"""Truth-point files (germquant.tune.points; design 6.5 and 8.7; plan section 9.3; plugin need I9):
annotations/<image_id>__points_<annotator>.csv with point_id, z_um, y_um, x_um, zi, yi, xi, object_id and
projected. zi, yi and xi lie on the image file's own grid (stride 1, the full z range) whatever run the points
were placed on, object_id is the run's label at the point (0 outside every object), z_um is 0.0 on a single
plane, and the file is ASCII with LF line endings."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from germquant.tune import points as P

DZ, DY, DX = 0.5, 0.2, 0.2


def _labels(shape=(10, 40, 40)):
    lab = np.zeros(shape, np.int32)
    lab[2:8, 5:15, 5:15] = 3
    lab[2:8, 25:35, 20:30] = 7
    return lab


def test_write_and_read_on_the_file_grid(tmp_path):
    full = _labels()
    s, z0 = 2, 2
    run_labels = full[z0:9, ::s, ::s]                  # a run with xy_stride 2 and z_range (2, 9)
    spacing = (DZ, DY * s, DX * s)
    # world coordinates of the run's layers: z from the run's first plane, y and x on the physical grid
    zyx = np.array([
        [1.0, 2.0, 2.0],             # nucleus 3 (file voxel z 4, y 10, x 10)
        [2.0, 6.0, 5.0],             # nucleus 7
        [0.0, 7.8, 7.8],             # background, far from both
        [1.5, 3.05, 1.0],            # just outside nucleus 3 in y: the neighbourhood rule assigns it
    ])
    out = P.write_points(tmp_path / P.points_file_name("img_01", "ryan"), zyx, labels=run_labels, spacing=spacing,
                         xy_stride=s, z_range=(z0, 9), projected=[False, False, True, False])
    raw = out.read_bytes()
    assert raw.isascii() and b"\r" not in raw
    assert raw.splitlines()[0] == b"point_id,z_um,y_um,x_um,zi,yi,xi,object_id,projected"
    df = P.read_points(out)
    assert list(df.columns) == P.POINT_COLUMNS
    assert df["point_id"].tolist() == [0, 1, 2, 3]
    assert df["z_um"].tolist() == [2.0, 3.0, 1.0, 2.5]            # the file's frame: plane z0 is z 1.0 um
    assert df["y_um"].tolist() == [2.0, 6.0, 7.8, 3.05] and df["x_um"].tolist() == [2.0, 5.0, 7.8, 1.0]
    assert df["zi"].tolist() == [4, 6, 2, 5]
    assert df["yi"].tolist() == [10, 30, 39, 15] and df["xi"].tolist() == [10, 25, 39, 5]
    assert df["object_id"].tolist() == [3, 7, 0, 3]
    assert df["projected"].tolist() == [False, False, True, False]
    for c in ("point_id", "zi", "yi", "xi", "object_id"):
        assert df[c].dtype == np.int64
    assert df["projected"].dtype == bool
    # every point lies at its (zi, yi, xi) of the full-resolution file
    assert [int(full[z, y, x]) for z, y, x in zip(df["zi"], df["yi"], df["xi"])][:2] == [3, 7]


def test_to_run_frame_moves_points_between_runs(tmp_path):
    full = _labels()
    out = P.write_points(tmp_path / "a__points_x.csv", [[1.0, 2.0, 2.0], [2.0, 6.0, 5.0]],
                         labels=full[2:9, ::2, ::2], spacing=(DZ, DY * 2, DX * 2), xy_stride=2, z_range=(2, 9))
    df = P.read_points(out)
    # onto a full-resolution run of the whole stack: the same physical points, the file's own voxels
    f = P.to_run_frame(df, shape=full.shape, spacing=(DZ, DY, DX))
    assert f["z_um"].tolist() == [2.0, 3.0] and f["zr"].tolist() == [4, 6]
    assert f["yr"].tolist() == [10, 30] and f["xr"].tolist() == [10, 25] and f["inside"].all()
    # onto a run of planes 5 to 9 only: the first point lies on a plane that run did not read
    f = P.to_run_frame(df, shape=(4, 40, 40), spacing=(DZ, DY, DX), z_range=(5, 9))
    assert f["inside"].tolist() == [False, True] and f["zr"].tolist() == [-1, 1]
    assert f["z_um"].tolist() == [-0.5, 0.5]


def test_single_plane_points(tmp_path):
    lab = _labels()[4]                                # (Y, X), as tifffile reads a single-plane label image
    out = P.write_points(tmp_path / "s__points_x.csv", [[0.0, 2.0, 2.0], [0.0, 6.0, 5.0]], labels=lab,
                         spacing=(float("nan"), DY, DX))
    df = P.read_points(out)
    assert df["z_um"].tolist() == [0.0, 0.0] and df["zi"].tolist() == [0, 0]
    assert df["object_id"].tolist() == [3, 7]
    f = P.to_run_frame(df, shape=(1, 40, 40), spacing=(1.0, DY, DX))
    assert f["z_um"].tolist() == [0.0, 0.0] and f["inside"].all()


def test_names_and_annotators():
    assert P.points_file_name("img_01", "ryan-2") == "img_01__points_ryan-2.csv"
    assert P.parse_points_name("img_01__points_ryan-2.csv") == ("img_01", "ryan-2")
    assert P.parse_points_name("x.ome__points_a.csv") == ("x.ome", "a")
    assert P.parse_points_name("img_01_points_ryan.csv") is None
    assert P.parse_points_name("img_01__points_ry an.csv") is None
    for bad in ("", "a" * 33, "ry an", "r" + chr(0xE9)):
        with pytest.raises(ValueError, match="annotator"):
            P.points_file_name("img", bad)
    for bad in ("a/b", "", " ", "a\tb"):
        with pytest.raises(ValueError, match="image id"):
            P.points_file_name(bad, "ryan")
    # an image id may hold non-ASCII letters, as project ids from file names do (review P5-TUNE N7)
    gonad = "gonad_" + chr(0xFC) + chr(0x3B3)
    assert P.points_file_name(gonad, "ryan") == gonad + "__points_ryan.csv"
    assert P.parse_points_name(P.points_file_name(gonad, "ryan")) == (gonad, "ryan")


def test_read_points_refusals(tmp_path):
    good = "point_id,z_um,y_um,x_um,zi,yi,xi,object_id,projected\n0,1.0,2.0,2.0,2,10,10,3,false\n"
    p = tmp_path / "a__points_x.csv"
    p.write_bytes(b"\xef\xbb\xbf" + good.replace("\n", "\r\n").encode("ascii"))     # Excel's BOM and CRLF
    assert P.read_points(p)["object_id"].tolist() == [3]
    cases = {
        "missing columns": good.replace(",projected", "").replace(",false", ""),
        "projected must be true or false": good.replace("false", "maybe"),
        "appears twice": good + "0,1.0,2.0,2.0,2,10,10,3,true\n",
        "needs a number": good + "1,,2.0,2.0,2,10,10,3,true\n",
        "whole numbers": good + "1,1.0,2.0,2.0,2.5,10,10,3,true\n",
        "never negative": good + "1,1.0,2.0,2.0,2,10,10,-1,true\n",
        "empty file": "",
    }
    for message, text in cases.items():
        p.write_text(text, encoding="ascii")
        with pytest.raises(ValueError, match=message):
            P.read_points(p)


def test_write_points_refusals(tmp_path):
    lab = _labels()
    with pytest.raises(ValueError, match="dz must be positive"):
        P.write_points(tmp_path / "a.csv", [[1, 1, 1]], labels=lab, spacing=(float("nan"), DY, DX))
    with pytest.raises(ValueError, match="finite"):
        P.write_points(tmp_path / "a.csv", [[np.nan, 1, 1]], labels=lab, spacing=(DZ, DY, DX))
    with pytest.raises(ValueError, match="xy_stride"):
        P.write_points(tmp_path / "a.csv", [[1, 1, 1]], labels=lab, spacing=(DZ, DY, DX), xy_stride=0)
    with pytest.raises(ValueError, match="one boolean per point"):
        P.write_points(tmp_path / "a.csv", [[1, 1, 1]], labels=lab, spacing=(DZ, DY, DX), projected=[True, False])
    empty = P.write_points(tmp_path / "e__points_x.csv", np.zeros((0, 3)), labels=lab, spacing=(DZ, DY, DX))
    assert P.read_points(empty).empty and pd.read_csv(empty).columns.tolist() == P.POINT_COLUMNS
