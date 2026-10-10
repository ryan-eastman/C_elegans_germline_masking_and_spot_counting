"""The Imaris truth reader (``germquant imaris-truth``): frame conversion, mask recovery, the truth file
contract and the CLI wiring, on a tiny synthetic HDF5 written here in the .ims layout (padded DataSet
arrays, Spot and SurfaceModelInfo compound tables, statistics, class label sets, creation parameters)."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

h5py = pytest.importorskip("h5py")

from germquant.validate import imaris_ims as I  # noqa: E402

# image geometry: X, Y, Z voxels, Imaris's voxel size and extent; the arrays on disk are padded to PAD (ZYX)
X, Y, Z = 40, 30, 12
VOX = np.array([0.1, 0.1, 0.5])
EXT_MIN = np.array([100.0, 200.0, 50.0])
EXT_MAX = EXT_MIN + np.array([X, Y, Z]) * VOX
PAD = (16, 32, 48)
OUR_VOX = (0.1003, 0.1003, 0.5)        # (dx, dy, dz): the nd2 voxel size that defines our frame
MASK_CH = 1

# surfaces: id -> (box (z0, z1, y0, y1, x0, x1) inclusive, mask value); 8 and 9 SHARE a value on separate nuclei
BOXES = {3: ((2, 6, 5, 12, 5, 12), 7), 8: ((4, 9, 10, 19, 20, 29), 12), 9: ((6, 10, 20, 27, 30, 36), 12)}
NO_MASK_ID = 15                        # a surface whose centroid lies in a zero region (no mask)
NO_MASK_VOXEL = (15, 25, 1)
ORPHAN = ((10, 10, 25, 26, 2, 3), 3)   # a tiny blob no surface claims
SURFACE_IDS = [3, 8, 9, NO_MASK_ID]
# spots: id -> (voxel (xi, yi, zi), fractional position inside the voxel, class, quality, mask value, radius)
SPOTS = {11: ((8, 8, 4), (0.3, 0.6, 0.5), "Class A", 50.5, 7, 0.2),
         12: ((25, 15, 6), (0.5, 0.5, 0.5), "Class A", 61.0, 12, 0.2),
         13: ((33, 23, 8), (0.25, 0.75, 0.5), "Class A", 42.1, 12, 0.2),
         14: ((1, 1, 0), (0.9, 0.9, 0.5), "Class B", 45.0, 0, 0.2),
         15: ((39, 29, 11), (0.5, 0.5, 0.5), "Class B", 99.0, 0, 0.35)}
SPOTS_CP = ('<bpPointsCreationParameters mApplicationType="eApplication3D" mFormatVersion="10.2" '
            'mEnableRegionOfInterest="false" mEnableRegionGrowing="false" mSourceChannelIndex="0">\n'
            '<bpRegionOfInterestContainer>\n</bpRegionOfInterestContainer>\n</bpPointsCreationParameters>\n')
SURF_CP = ('<bpSurfacesCreationParameters mApplicationType="eApplication3D" mFormatVersion="10.2" '
           'mEnableRegionOfInterest="true" mEnableRegionGrowing="true" mSurfaceGrainSize="0.4" '
           'mManualThresholdValue="36.2618">\n<bpRegionOfInterestContainer>\n'
           '<bpRegionOfInterest mMinXYZ="1 2 0" mMaxXYZ="38 28 11" mMinT="0" mMaxT="0" '
           'mName="Region of Interest 1" mActive="true" mSelected="false"/>\n'
           '</bpRegionOfInterestContainer>\n</bpSurfacesCreationParameters>\n')


def _s(text: str) -> np.ndarray:
    """An Imaris text attribute: an array of single bytes."""
    return np.frombuffer(text.encode("latin1"), dtype="S1")


def _rec(rows, dtype):
    return np.array([tuple(r) for r in rows], dtype=dtype)


def _fill(arr, box, value):
    z0, z1, y0, y1, x0, x1 = box
    arr[z0:z1 + 1, y0:y1 + 1, x0:x1 + 1] = value


def _box_centroid_um(box):
    z0, z1, y0, y1, x0, x1 = box
    mean_xyz = np.array([(x0 + x1) / 2, (y0 + y1) / 2, (z0 + z1) / 2])
    return EXT_MIN + (mean_xyz + 0.5) * VOX


def _box_voxels(box):
    z0, z1, y0, y1, x0, x1 = box
    return (z1 - z0 + 1) * (y1 - y0 + 1) * (x1 - x0 + 1)


def _spot_position_um(voxel, frac):
    return EXT_MIN + (np.array(voxel) + np.array(frac)) * VOX


def _label_tables(g, sets, label_ids, object_ids):
    """The Imaris class-label layout: values 'Class A' / 'Class B' in group 'Set 1', objects grouped by
    their label combination."""
    g.create_dataset("LabelGroupNames", data=_rec([("Set 1", 2)], [("LabelGroupName", "S256"), ("EndLabelValue", "<i8")]))
    g.create_dataset("LabelValues", data=_rec([("Class A",), ("Class B",)], [("LabelValue", "S256")]))
    g.create_dataset("LabelSets", data=_rec(sets, [("EndLabelIDs", "<i8"), ("EndObjectIDs", "<i8")]))
    g.create_dataset("LabelSetLabelIDs", data=_rec([(i,) for i in label_ids], [("IDLabel", "<i8")]))
    g.create_dataset("LabelSetObjectIDs", data=_rec([(i,) for i in object_ids], [("IDObject", "<i8")]))


def build_ims(path, *, mask_mode="unique", surface_stats=True, surfaces=True):
    """Write a synthetic .ims. ``mask_mode`` 'unique' gives every nucleus its own value (8 and 9 collide),
    'binary' fills every nucleus with 1; ``surface_stats=False`` leaves out the per-surface intensity
    statistics so the pixel probe has to find the mask channel."""
    rng = np.random.default_rng(0)
    raw = np.zeros(PAD, np.uint16)
    raw[:Z, :Y, :X] = rng.integers(60, 400, size=(Z, Y, X))
    mask = np.zeros(PAD, np.uint16)
    for box, value in BOXES.values():
        _fill(mask, box, 1 if mask_mode == "binary" else value)
    _fill(mask, ORPHAN[0], 1 if mask_mode == "binary" else ORPHAN[1])
    stype = [("ID", "<i8"), ("ID_Category", "<i8"), ("ID_FactorList", "<i8"), ("Name", "S256"), ("Unit", "S256")]
    svalue = [("ID_Time", "<i8"), ("ID_Object", "<i8"), ("ID_StatisticsType", "<i8"), ("Value", "<f8")]
    factor = [("ID_List", "<i8"), ("Name", "S256"), ("Level", "S256")]
    category = [("ID", "<i8"), ("CategoryName", "S256"), ("Name", "S256")]
    with h5py.File(path, "w") as f:
        info = f.create_group("DataSetInfo")
        img = info.create_group("Image")
        for i in range(3):
            img.attrs[f"ExtMin{i}"] = _s(f"{EXT_MIN[i]:g}")
            img.attrs[f"ExtMax{i}"] = _s(f"{EXT_MAX[i]:g}")
        img.attrs["X"], img.attrs["Y"], img.attrs["Z"] = _s(str(X)), _s(str(Y)), _s(str(Z))
        img.attrs["Unit"] = _s("um")
        img.attrs["Name"] = _s("synthetic gonad")
        img.attrs["Filename"] = _s("/Volumes/somewhere/synthetic.nd2")
        info.create_group("Imaris").attrs["Version"] = _s("10.2")
        for c, name in enumerate(["405", "Masked 405"]):
            info.create_group(f"Channel {c}").attrs["Name"] = _s(name)
        tp = f.create_group("DataSet/ResolutionLevel 0/TimePoint 0")
        for c, arr in enumerate([raw, mask]):
            g = tp.create_group(f"Channel {c}")
            g.create_dataset("Data", data=arr, chunks=(4, 16, 16), compression="gzip")
            g.attrs["ImageSizeX"], g.attrs["ImageSizeY"], g.attrs["ImageSizeZ"] = _s(str(X)), _s(str(Y)), _s(str(Z))
        content = f.create_group("Scene8/Content")
        content.attrs["NumberOfPoints"] = np.uint64(1)
        content.attrs["NumberOfMegaSurfaces"] = np.uint64(1 if surfaces else 0)

        # ---- Spots 1 ----
        pts = content.create_group("Points0")
        pts.attrs["Name"] = _s("Spots 1")
        pts.attrs["CreationParameters"] = _s(SPOTS_CP[:80])     # Imaris keeps a shorter copy as attribute
        pts.create_dataset("CreationParameters", data=np.array([SPOTS_CP.encode()], dtype=f"S{len(SPOTS_CP)}"))
        spot_dtype = [("ID", "<i8"), ("PositionX", "<f4"), ("PositionY", "<f4"), ("PositionZ", "<f4"), ("Radius", "<f4")]
        rows = [(sid, *_spot_position_um(v, fr), rad) for sid, (v, fr, _c, _q, _m, rad) in SPOTS.items()]
        pts.create_dataset("Spot", data=_rec(rows, spot_dtype))
        pts.create_dataset("StatisticsType", data=_rec([(1, 0, 0, "Quality", ""), (2, 0, 1, "Intensity Max", ""),
                                                        (3, 1, 0, "Total Number of Spots", "")], stype))
        pts.create_dataset("Category", data=_rec([(0, "Spot", "Spot"), (1, "Overall", "Overall")], category))
        pts.create_dataset("Factor", data=_rec([(1, "Channel", str(MASK_CH + 1)), (1, "Image", "Image 1")], factor))
        vals = [(0, -1, 3, float(len(SPOTS)))]
        for sid, (_v, _fr, _c, q, m, _r) in SPOTS.items():
            vals += [(0, sid, 1, q), (0, sid, 2, float(m if mask_mode == "unique" else min(m, 1)))]
        pts.create_dataset("StatisticsValue", data=_rec(vals, svalue))
        _label_tables(pts, sets=[(1, 3), (2, 5)], label_ids=[0, 1], object_ids=[11, 13, 12, 14, 15])

        if not surfaces:
            return
        # ---- Surfaces 1 ----
        sg = content.create_group("MegaSurfaces0")
        sg.attrs["Name"] = _s("Surfaces 1")
        sg.create_dataset("CreationParameters", data=np.array([SURF_CP.encode()], dtype=f"S{len(SURF_CP)}"))
        smi_dtype = [("Volume", "<f4"), ("Area", "<f4"), ("CenterOfMassX", "<f8"), ("CenterOfMassY", "<f8"),
                     ("CenterOfMassZ", "<f8")]
        sm_dtype = [("ID_Time", "<i8"), ("ID", "<i8"), ("ExtendMinX", "<f4"), ("ExtendMinY", "<f4"), ("ExtendMinZ", "<f4"),
                    ("ExtendMaxX", "<f4"), ("ExtendMaxY", "<f4"), ("ExtendMaxZ", "<f4")]
        smi_rows, sm_rows, st_vals = [], [], []
        vox_vol = float(np.prod(VOX))
        for sid in SURFACE_IDS:
            if sid in BOXES:
                box, value = BOXES[sid]
                com = _box_centroid_um(box)
                vol = _box_voxels(box) * vox_vol
                z0, z1, y0, y1, x0, x1 = box
                bmin = EXT_MIN + np.array([x0, y0, z0]) * VOX
                bmax = EXT_MIN + np.array([x1 + 1, y1 + 1, z1 + 1]) * VOX
                mval = float(value if mask_mode == "unique" else 1)
            else:
                com = EXT_MIN + (np.array(NO_MASK_VOXEL) + 0.5) * VOX
                vol = 10.0
                bmin, bmax = com - 1.0, com + 1.0
                mval = 0.0
            smi_rows.append((vol, 5.0, *com))
            sm_rows.append((0, sid, *bmin, *bmax))
            # channel 1 (raw) varies inside a surface, channel 2 (mask) is constant: the statistics signature
            st_vals += [(0, sid, 2103, 390.0 + sid), (0, sid, 2115, 61.0), (0, sid, 2102, mval), (0, sid, 2114, mval)]
        sg.create_dataset("SurfaceModelInfo", data=_rec(smi_rows, smi_dtype))
        sg.create_dataset("SurfaceModel", data=_rec(sm_rows, sm_dtype))
        if surface_stats:
            sg.create_dataset("StatisticsType", data=_rec([(2102, 0, 14, "Intensity Max", ""), (2114, 0, 14, "Intensity Min", ""),
                                                           (2103, 0, 13, "Intensity Max", ""), (2115, 0, 13, "Intensity Min", ""),
                                                           (245, 1, 0, "Total Number of Surfaces", "")], stype))
            sg.create_dataset("Category", data=_rec([(0, "Surface", "Surface"), (1, "Overall", "Overall")], category))
            sg.create_dataset("Factor", data=_rec([(13, "Channel", "1"), (13, "Image", "Image 1"),
                                                   (14, "Channel", "2"), (14, "Image", "Image 1")], factor))
            sg.create_dataset("StatisticsValue", data=_rec(st_vals + [(0, -1, 245, 4.0)], svalue))
        _label_tables(sg, sets=[(1, 1), (3, 2), (4, 4)], label_ids=[0, 0, 1, 1], object_ids=SURFACE_IDS)


@pytest.fixture(scope="module")
def truth(tmp_path_factory):
    d = tmp_path_factory.mktemp("ims")
    ims = d / "20251105_N2_test_HERM _001.ims"      # the lab's stems end in a space and an underscore
    build_ims(ims)
    out = d / "truth"
    params = I.export_imaris_truth(ims, out, nd2_voxel_xyz_um=OUR_VOX)
    return ims, out, params


def test_truth_paths_follow_the_contract(tmp_path):
    p = I.truth_paths(tmp_path, "20251105_N2_nohs_HERM _")
    assert p["spots"].name == "20251105_N2_nohs_HERM ___imaris_spots.csv"
    assert p["nuclei"].name == "20251105_N2_nohs_HERM ___imaris_nuclei.csv"
    assert p["labels"].name == "20251105_N2_nohs_HERM ___imaris_labels.tif"
    assert p["params"].name == "20251105_N2_nohs_HERM ___imaris_params.json"


def test_spots_table_and_frame_conversion(truth):
    ims, out, params = truth
    spots = pd.read_csv(I.truth_paths(out, ims.stem)["spots"])
    assert list(spots.columns) == I.SPOT_COLUMNS + I.SPOT_EXTRA_COLUMNS
    assert len(spots) == len(SPOTS) and spots["spot_id"].tolist() == list(SPOTS)
    for _, r in spots.iterrows():
        voxel, frac, cls, q, m, rad = SPOTS[int(r["spot_id"])]
        # voxel index from Imaris's own voxel size (floor), our micrometres from the nd2 voxel size
        assert (int(r["xi"]), int(r["yi"]), int(r["zi"])) == voxel
        assert r["x_um"] == pytest.approx(voxel[0] * OUR_VOX[0], abs=1e-9)
        assert r["y_um"] == pytest.approx(voxel[1] * OUR_VOX[1], abs=1e-9)
        assert r["z_um"] == pytest.approx(voxel[2] * OUR_VOX[2], abs=1e-9)
        raw = _spot_position_um(voxel, frac)
        assert r["imaris_x_um"] == pytest.approx(raw[0], abs=1e-3)   # stored as float32 by Imaris
        assert r["imaris_z_um"] == pytest.approx(raw[2], abs=1e-3)
        assert r["radius_um"] == pytest.approx(rad, abs=1e-6)
        assert r["quality"] == pytest.approx(q)
        assert r["spot_class"] == cls
        assert int(r["imaris_mask_label"]) == m
        assert r["spots_object"] == "Spots 1"
    # the frame lesson: index-wise conversion is not "Imaris micrometres minus the extent minimum"
    last = spots[spots["spot_id"] == 15].iloc[0]
    assert last["x_um"] == pytest.approx(39 * OUR_VOX[0]) and abs(last["x_um"] - 39.5 * VOX[0]) > 0.03
    assert params["n_spots"] == 5 and params["quality_present"] is True and params["region_growing"] is False
    assert params["spot_classes_present"] == ["Class A", "Class B"]


def test_nuclei_table_and_mask_recovery(truth):
    import tifffile

    ims, out, params = truth
    paths = I.truth_paths(out, ims.stem)
    nuc = pd.read_csv(paths["nuclei"])
    assert list(nuc.columns) == I.NUCLEI_COLUMNS + I.NUCLEI_EXTRA_COLUMNS
    assert nuc["surface_id"].tolist() == SURFACE_IDS
    by_id = nuc.set_index("surface_id")
    # the two surfaces sharing value 12 on separate nuclei: the first keeps 12, the second gets a fresh label
    assert int(by_id.loc[3, "label"]) == 7 and int(by_id.loc[8, "label"]) == 12 and int(by_id.loc[9, "label"]) == 13
    assert int(by_id.loc[NO_MASK_ID, "label"]) == 0 and pd.isna(by_id.loc[NO_MASK_ID, "n_voxels_mask"])
    for sid, (box, _v) in BOXES.items():
        assert int(by_id.loc[sid, "n_voxels_mask"]) == _box_voxels(box)
        assert by_id.loc[sid, "mask_centroid_dist_um"] < 1e-6
        c = _box_centroid_um(box)
        idx = np.floor((c - EXT_MIN) / VOX).astype(int)
        assert (int(by_id.loc[sid, "xi"]), int(by_id.loc[sid, "yi"]), int(by_id.loc[sid, "zi"])) == tuple(idx)
        assert by_id.loc[sid, "x_um"] == pytest.approx(idx[0] * OUR_VOX[0])
        assert by_id.loc[sid, "volume_um3"] == pytest.approx(_box_voxels(box) * float(np.prod(VOX)), rel=1e-5)
    assert by_id.loc[8, "surface_class"] == "Class A;Class B" and by_id.loc[9, "surface_class"] == "Class B"
    # the label image: int32 ZYX of the IMAGE size (not the padded array), orphan value gone, OME voxel size tags
    with tifffile.TiffFile(paths["labels"]) as t:
        labels = t.asarray()
        assert t.is_ome and "PhysicalSizeX" in (t.ome_metadata or "") and f'PhysicalSizeX="{OUR_VOX[0]}"' in t.ome_metadata
        assert f'PhysicalSizeZ="{OUR_VOX[2]}"' in t.ome_metadata
    assert labels.shape == (Z, Y, X) and labels.dtype == np.int32
    assert sorted(np.unique(labels).tolist()) == [0, 7, 12, 13]
    for sid, (box, _v) in BOXES.items():
        z0, z1, y0, y1, x0, x1 = box
        block = labels[z0:z1 + 1, y0:y1 + 1, x0:x1 + 1]
        assert (block == int(by_id.loc[sid, "label"])).all()
    oz0, oz1, oy0, oy1, ox0, ox1 = ORPHAN[0]
    assert (labels[oz0:oz1 + 1, oy0:oy1 + 1, ox0:ox1 + 1] == 0).all()
    m = params["mask"]
    assert m["channel"] == MASK_CH and m["channel_name"] == "Masked 405" and m["labels_file"] == paths["labels"].name
    assert m["identification"]["method"] == "centroid_pixel_probe" and m["identification"]["values_unique"] is True
    assert m["identification"]["candidates"][MASK_CH]["n_constant"] == 3 and m["identification"]["candidates"][0]["n_constant"] == 0
    assert m["intensity_max_statistic_present"] is True
    assert m["n_surfaces_with_mask"] == 3 and m["n_surfaces"] == 4 and m["strategy"] == "value_per_surface"
    assert m["labels_bbox_xyz_voxel"] == {"min_xyz": [5, 5, 2], "max_xyz": [36, 27, 10]}      # the orphan is gone
    assert m["n_shared_values"] == 1 and m["n_shared_split_by_components"] == 1 and m["n_shared_split_by_nearest_centroid"] == 0
    assert m["n_orphan_values"] == 1 and m["orphan_values"] == [3] and m["n_orphan_voxels"] == 4
    assert m["n_centroid_inside_own_label"] == 3 and m["n_centroid_in_cavity"] == 0 and m["n_centroid_in_other_label"] == 0
    # surface 9 carries the fresh label 13 while Imaris's statistic says 12: the cross-check reports it
    assert m["n_labels_equal_to_intensity_statistic"] == 2 and m["n_labels_differing_from_intensity_statistic"] == 1
    assert m["labels_shape_zyx"] == [Z, Y, X]


def test_params_json(truth):
    ims, out, params = truth
    paths = I.truth_paths(out, ims.stem)
    with open(paths["params"], encoding="utf-8") as fh:
        p = json.load(fh)
    assert p["stem"] == ims.stem and p["imaris_version"] == "10.2"
    assert p["voxel_imaris_xyz_um"] == pytest.approx(VOX.tolist()) and p["voxel_ours_xyz_um"] == pytest.approx(list(OUR_VOX))
    assert p["our_voxel_source"] == "given"
    assert p["ext_min_xyz_um"] == pytest.approx(EXT_MIN.tolist()) and p["shape_xyz"] == [X, Y, Z]
    assert tuple(p["data_shape_zyx_padded"]) == PAD and p["channel_names"] == ["405", "Masked 405"]
    assert p["n_spots"] == 5 and p["n_surfaces"] == 4
    assert p["creation_parameters"]["Spots 1"] == SPOTS_CP and p["creation_parameters"]["Surfaces 1"] == SURF_CP
    assert p["roi_boxes"] == [{"object": "Surfaces 1", "name": "Region of Interest 1", "min_xyz_voxel": [1, 2, 0],
                               "max_xyz_voxel": [38, 28, 11], "min_t": [0], "max_t": [0], "active": True}]
    assert p["roi_check"]["n_surface_centroids_inside"] == 4
    assert p["surface_classes_present"] == ["Class A", "Class B"]
    kinds = {o["name"]: o for o in p["objects"]}
    assert kinds["Spots 1"]["region_growing"] is False and kinds["Surfaces 1"]["region_growing"] is True
    assert kinds["Spots 1"]["class_counts"] == {"Class A": 3, "Class B": 2}
    assert kinds["Surfaces 1"]["label_groups"] == {"Set 1": ["Class A", "Class B"]}
    assert "Quality" in p["spot_statistics_present"]
    assert sorted(p["files"]) == ["labels", "nuclei", "params", "spots"]
    assert p["warnings"] == ["1 of 4 surfaces have no mask"]      # surface 15 was built without one


def test_no_masks_keeps_tables_and_statistic_labels(tmp_path):
    ims = tmp_path / "g.ims"
    build_ims(ims)
    params = I.export_imaris_truth(ims, tmp_path / "out", nd2_voxel_xyz_um=OUR_VOX, masks=False)
    paths = I.truth_paths(tmp_path / "out", "g")
    assert paths["spots"].exists() and paths["nuclei"].exists() and paths["params"].exists()
    assert not paths["labels"].exists() and params["mask"]["labels_file"] is None
    assert params["mask"]["reason_no_labels"] == "masks not requested (--no-masks)"
    assert params["mask"]["channel"] == MASK_CH and "unverified" in params["mask"]["label_source"]
    nuc = pd.read_csv(paths["nuclei"]).set_index("surface_id")
    # labels from the per-surface statistic (the lab's coloc trick), no voxel counts
    assert nuc["label"].to_dict() == {3: 7, 8: 12, 9: 12, NO_MASK_ID: 0}
    assert nuc["n_voxels_mask"].isna().all()
    assert "labels" not in params["files"]


def test_binary_mask_falls_back_to_connected_components(tmp_path):
    ims = tmp_path / "b.ims"
    build_ims(ims, mask_mode="binary")
    params = I.export_imaris_truth(ims, tmp_path / "out", nd2_voxel_xyz_um=OUR_VOX)
    m = params["mask"]
    assert m["channel"] == MASK_CH and m["strategy"] == "connected_components_of_nonzero"
    assert m["identification"]["values_unique"] is False
    nuc = pd.read_csv(I.truth_paths(tmp_path / "out", "b")["nuclei"]).set_index("surface_id")
    labs = {sid: int(nuc.loc[sid, "label"]) for sid in BOXES}
    assert len(set(labs.values())) == 3 and all(v > 0 for v in labs.values())
    for sid, (box, _v) in BOXES.items():
        assert int(nuc.loc[sid, "n_voxels_mask"]) == _box_voxels(box)
    assert int(nuc.loc[NO_MASK_ID, "label"]) == 0
    assert m["n_orphan_values"] == 1 and m["n_orphan_voxels"] == 4     # the tiny blob is its own component


def test_export_without_surface_statistics(tmp_path):
    ims = tmp_path / "p.ims"
    build_ims(ims, surface_stats=False)
    params = I.export_imaris_truth(ims, tmp_path / "out", nd2_voxel_xyz_um=OUR_VOX)
    m = params["mask"]
    assert m["channel"] == MASK_CH and m["intensity_max_statistic_present"] is False
    nuc = pd.read_csv(I.truth_paths(tmp_path / "out", "p")["nuclei"]).set_index("surface_id")
    assert nuc["label"].to_dict() == {3: 7, 8: 12, 9: 13, NO_MASK_ID: 0}
    assert "n_labels_equal_to_intensity_statistic" not in m
    # without the statistic and without reading the masks there is no label source at all
    params = I.export_imaris_truth(ims, tmp_path / "out2", nd2_voxel_xyz_um=OUR_VOX, masks=False)
    nuc = pd.read_csv(I.truth_paths(tmp_path / "out2", "p")["nuclei"]).set_index("surface_id")
    assert (nuc["label"] == 0).all() and params["mask"]["label_source"] is None


def test_imaris_fallback_voxel_and_no_surfaces(tmp_path):
    ims = tmp_path / "s.ims"
    build_ims(ims, surfaces=False)
    params = I.export_imaris_truth(ims, tmp_path / "out")
    assert params["our_voxel_source"] == "imaris_fallback"
    assert params["voxel_ours_xyz_um"] == pytest.approx(VOX.tolist())
    assert params["n_surfaces"] == 0 and params["mask"]["reason_no_labels"] == "the file holds no surfaces"
    spots = pd.read_csv(I.truth_paths(tmp_path / "out", "s")["spots"])
    assert len(spots) == 5 and spots["imaris_mask_label"].isna().all()
    assert spots.loc[spots["spot_id"] == 11, "x_um"].iloc[0] == pytest.approx(8 * VOX[0])
    assert any("Imaris's rounded voxel size" in w for w in params["warnings"])


def test_nd2_header_defines_our_frame(tmp_path, monkeypatch):
    from germquant import io as gio

    ims = tmp_path / "n.ims"
    build_ims(ims, surfaces=False)
    seen = {}

    def fake_meta(path):
        seen["path"] = path
        return {"spacing": (OUR_VOX[2], OUR_VOX[1], OUR_VOX[0])}      # (dz, dy, dx) like read_nd2_metadata

    monkeypatch.setattr(gio, "read_image_metadata", fake_meta)
    params = I.export_imaris_truth(ims, tmp_path / "out", nd2=tmp_path / "n.nd2")
    assert seen["path"].endswith("n.nd2") and params["our_voxel_source"].startswith("nd2:")
    assert params["voxel_ours_xyz_um"] == pytest.approx(list(OUR_VOX))


def test_read_imaris_ims_still_works_on_the_synthetic_file(tmp_path):
    ims = tmp_path / "r.ims"
    build_ims(ims)
    spots, surfaces, geom = I.read_imaris_ims(ims)
    assert len(spots) == 5 and len(surfaces) == 4
    assert geom["voxel_xyz_um"] == pytest.approx(VOX)
    assert spots["x_um"].iloc[0] == pytest.approx((8 + 0.3) * VOX[0], abs=1e-3)


def test_cli_wiring(monkeypatch):
    from germquant import cli
    from germquant.validate import imaris_ims

    called = {}

    def fake_export(ims, out, *, nd2=None, masks=True):
        called.update(ims=ims, out=out, nd2=nd2, masks=masks)
        return {"stem": "x", "n_spots": 3, "n_surfaces": 2, "warnings": ["w1"],
                "mask": {"labels_file": None, "reason_no_labels": "masks not requested (--no-masks)"}}

    monkeypatch.setattr(imaris_ims, "export_imaris_truth", fake_export)
    assert cli.main(["imaris-truth", "a.ims", "--out", "o", "--nd2", "a.nd2", "--no-masks"]) == 0
    assert called == {"ims": "a.ims", "out": "o", "nd2": "a.nd2", "masks": False}
    assert cli.main(["imaris-truth", "b.ims", "--out", "o2"]) == 0
    assert called["nd2"] is None and called["masks"] is True
