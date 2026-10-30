"""Imported SC geometry: a tiny Imaris file written with h5py (Filaments0 with a branched filament, a path and a
single-vertex stub), SWC and CSV readers, the continuous Imaris frame rule, and an imported path measured through
the pipeline (verified class, as-drawn ends)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from germquant.sc.straighten import imports, params, synth
from germquant.sc.straighten.pipeline import straighten_image

h5py = pytest.importorskip("h5py")

EXT_MIN = (100.0, 200.0, 10.0)          # x, y, z (Imaris stage um)
V_IMS = (0.108, 0.108, 0.2)              # x, y, z Imaris voxel
SHAPE_XYZ = (64, 48, 20)
OUR_ZYX = (0.2, 0.108333, 0.108333)


def _s(text: str) -> np.ndarray:
    return np.frombuffer(text.encode(), dtype="S1")


def _rec(rows, fields):
    return np.array(rows, dtype=fields)


def _write_ims(path):
    with h5py.File(path, "w") as f:
        info = f.create_group("DataSetInfo")
        img = info.create_group("Image")
        for i in range(3):
            img.attrs[f"ExtMin{i}"] = _s(f"{EXT_MIN[i]:g}")
            img.attrs[f"ExtMax{i}"] = _s(f"{EXT_MIN[i] + SHAPE_XYZ[i] * V_IMS[i]:g}")
        img.attrs["X"], img.attrs["Y"], img.attrs["Z"] = (_s(str(v)) for v in SHAPE_XYZ)
        info.create_group("Channel 0").attrs["Name"] = _s("mCherry")
        tp = f.create_group("DataSet/ResolutionLevel 0/TimePoint 0/Channel 0")
        tp.create_dataset("Data", data=np.zeros(SHAPE_XYZ[::-1], np.uint16))
        g = f.create_group("Scene8/Content/Filaments0")
        # vertices at voxel centres: index (i, j, k) -> ExtMin + (index + 0.5) * v
        def pos(ix, iy, iz):
            return (EXT_MIN[0] + (ix + 0.5) * V_IMS[0], EXT_MIN[1] + (iy + 0.5) * V_IMS[1],
                    EXT_MIN[2] + (iz + 0.5) * V_IMS[2], 0.1)
        verts = [pos(10, 10, 5), pos(20, 10, 5), pos(30, 10, 5), pos(20, 20, 5),        # Y: centre index 1
                 pos(5, 30, 8), pos(15, 30, 8), pos(25, 31, 9),                             # a path (global edges)
                 pos(40, 40, 10)]                                                           # a stub
        g.create_dataset("Vertex", data=_rec(verts, [("PositionX", "<f4"), ("PositionY", "<f4"),
                                                     ("PositionZ", "<f4"), ("Radius", "<f4")]))
        edges = [(0, 1), (1, 2), (1, 3), (4, 5), (5, 6)]
        g.create_dataset("Edge", data=_rec(edges, [("VertexA", "<i8"), ("VertexB", "<i8")]))
        fil = [(0, 100, 0, 0, 4, 0, 3), (0, 101, 4, 4, 7, 3, 5), (0, 102, 7, 7, 8, 5, 5)]
        g.create_dataset("Filament", data=_rec(fil, [("ID_Time", "<i8"), ("ID", "<i8"), ("IndexVertexRoot", "<i8"),
                                                     ("IndexVertexBegin", "<i8"), ("IndexVertexEnd", "<i8"),
                                                     ("IndexEdgeBegin", "<i8"), ("IndexEdgeEnd", "<i8")]))
        st = [(1, 0, 0, b"Filament Length (sum)", b"um")]
        g.create_dataset("StatisticsType", data=_rec(st, [("ID", "<i8"), ("ID_Category", "<i8"),
                                                          ("ID_FactorList", "<i8"), ("Name", "S256"), ("Unit", "S16")]))
        g.create_dataset("StatisticsValue", data=_rec([(0, 101, 1, 2.2)], [("ID_Time", "<i8"), ("ID_Object", "<i8"),
                                                                            ("ID_StatisticsType", "<i8"),
                                                                            ("Value", "<f4")]))


def test_ims_filaments_paths_frame_and_stubs(tmp_path):
    p = tmp_path / "t.ims"
    _write_ims(p)
    paths = imports.read_ims_filaments(p, OUR_ZYX)
    ids = sorted(x.path_id for x in paths)
    assert len([i for i in ids if ":100:" in i]) == 3          # the Y splits into three arms at its branch point
    assert len([i for i in ids if ":101:" in i]) == 1          # a path with global edge indices
    assert not any(":102:" in i for i in ids)                  # the single-vertex stub is ignored
    p101 = next(x for x in paths if ":101:" in x.path_id)
    # continuous frame: voxel centre (ix, iy, iz) -> (iz * dz, iy * dy, ix * dx)
    want = np.array([[8 * OUR_ZYX[0], 30 * OUR_ZYX[1], 5 * OUR_ZYX[2]], [9 * OUR_ZYX[0], 31 * OUR_ZYX[1], 25 * OUR_ZYX[2]]])
    got = p101.points[[0, -1]]
    if np.linalg.norm(got[0] - want[0]) > np.linalg.norm(got[0] - want[1]):
        got = got[::-1]
    assert np.allclose(got, want, atol=1e-4)
    assert p101.stored_length_um == pytest.approx(2.2, abs=1e-5)


def test_swc_and_csv_readers(tmp_path):
    swc = tmp_path / "a.swc"
    swc.write_text("# id type x y z r parent\n1 0 0 0 0 1 -1\n2 0 1 0 0 1 1\n3 0 2 0 0 1 2\n4 0 1 1 0 1 2\n")
    sp = imports.read_swc(swc)
    assert len(sp) == 3                                      # split at node 2 (degree 3)
    csv = tmp_path / "p.csv"
    pd.DataFrame({"image_id": ["a"] * 3 + ["b"] * 2, "strand_id": [1, 1, 1, 2, 2], "point_index": [2, 0, 1, 0, 1],
                  "z_um": [0, 0, 0, 1, 1], "y_um": [0, 0, 0, 1, 1], "x_um": [2, 0, 1, 0, 1],
                  "verified": [True] * 5}).to_csv(csv, index=False)
    a = imports.read_paths_csv(csv, "a")
    assert len(a) == 1 and np.allclose(a[0].points[:, 2], [0, 1, 2]) and a[0].verified


def test_imported_path_measured_as_drawn():
    sc = synth.iso_scene(31, L=4.0, peak=500, tilt_deg=10)
    rib = sc.ribbons[0]
    drawn = rib.points[::40]                                   # a hand trace every 0.2 um
    res = straighten_image(sc.channels["central_element"], sc.labels, sc.spacing, sc.nuclei, params=params.resolve(),
                           imported={1: [{"points": drawn, "verified": True, "id": "hand:1"}]}, geometry="import",
                           refine_imported=False, ends_mode="as_drawn")
    st = res.strands
    assert len(st) == 1 and st["strand_class"].iloc[0] == "verified"
    assert st["geometry_source"].iloc[0] == "imported" and "imported" in st["qc_flags"].iloc[0]
    assert abs(st["length_um"].iloc[0] - rib.length_um) < 0.05


def test_assign_nucleus():
    lab = np.zeros((10, 20, 20), np.uint16)
    lab[2:8, 2:10, 2:10] = 3
    lab[2:8, 12:18, 12:18] = 4
    sp = (0.2, 0.1, 0.1)
    pts = np.array([[1.0, 0.5, 0.5], [1.0, 0.6, 0.6], [1.0, 1.4, 1.4]])
    assert imports.assign_nucleus(pts, lab, sp) == 3
    near = np.array([[1.0, 1.05, 1.05]])                        # just outside label 3 (index 10.5)
    assert imports.assign_nucleus(near, lab, sp) == 3
