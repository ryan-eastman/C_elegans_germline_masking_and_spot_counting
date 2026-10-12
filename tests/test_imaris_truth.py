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


def build_ims(path, *, mask_mode="unique", surface_stats=True, surfaces=True, ic_stats=False):
    """Write a synthetic .ims. ``mask_mode`` 'unique' gives every nucleus its own value (8 and 9 collide),
    'binary' fills every nucleus with 1; ``surface_stats=False`` leaves out the per-surface intensity
    statistics so the pixel probe has to find the mask channel; ``ic_stats`` stores Imaris's per-spot
    Intensity Center of the spot channel (0) and of the mask channel (1), the frame check's input."""
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
        ic_types = [(4, 0, 2, "Intensity Center", ""), (5, 0, 1, "Intensity Center", "")] if ic_stats else []
        pts.create_dataset("StatisticsType", data=_rec([(1, 0, 0, "Quality", ""), (2, 0, 1, "Intensity Max", ""),
                                                        (3, 1, 0, "Total Number of Spots", "")] + ic_types, stype))
        pts.create_dataset("Category", data=_rec([(0, "Spot", "Spot"), (1, "Overall", "Overall")], category))
        ic_factors = [(2, "Channel", "1"), (2, "Image", "Image 1")] if ic_stats else []
        pts.create_dataset("Factor", data=_rec([(1, "Channel", str(MASK_CH + 1)), (1, "Image", "Image 1")] + ic_factors,
                                               factor))
        vals = [(0, -1, 3, float(len(SPOTS)))]
        for sid, (_v, _fr, _c, q, m, _r) in SPOTS.items():
            vals += [(0, sid, 1, q), (0, sid, 2, float(m if mask_mode == "unique" else min(m, 1)))]
            if ic_stats:
                xi, yi, zi = _v
                vals += [(0, sid, 4, float(raw[zi, yi, xi])), (0, sid, 5, float(mask[zi, yi, xi]))]
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

    def fake_export(ims, out, *, nd2=None, masks=True, frame_check_on=False, spots_xlsx=None, xlsx_expect=None,
                    allow_test=False):
        called.update(ims=ims, out=out, nd2=nd2, masks=masks)
        called.update(frame=frame_check_on, xlsx=spots_xlsx, expect=xlsx_expect, allow_test=allow_test)
        return {"stem": "x", "n_spots": 3, "n_surfaces": 2, "warnings": ["w1"],
                "mask": {"labels_file": None, "reason_no_labels": "masks not requested (--no-masks)"}}

    monkeypatch.setattr(imaris_ims, "export_imaris_truth", fake_export)
    assert cli.main(["imaris-truth", "a.ims", "--out", "o", "--nd2", "a.nd2", "--no-masks"]) == 0
    assert called == {"ims": "a.ims", "out": "o", "nd2": "a.nd2", "masks": False, "frame": False, "xlsx": None,
                      "expect": None, "allow_test": False}
    assert cli.main(["imaris-truth", "c.ims", "--out", "o3", "--frame-check", "--spots-xlsx", "c.xlsx",
                     "--xlsx-expect-sha256", "ab", "--xlsx-expect-spots", "5", "--xlsx-expect-surfaces", "4"]) == 0
    assert called["frame"] is True and called["xlsx"] == "c.xlsx" and called["allow_test"] is False
    assert called["expect"] == {"sha256": "ab", "n_spots": 5, "n_surfaces": 4}
    assert cli.main(["imaris-truth", "b.ims", "--out", "o2"]) == 0
    assert called["nd2"] is None and called["masks"] is True
    assert cli.main(["imaris-truth", "d.ims", "--out", "o4", "--frame-check", "--pre-flight"]) == 0
    assert called["frame"] is True and called["allow_test"] is True


# ---------------------------------------------------------------- the frame check (entry 9 (d)(3))
def _frame_case(seed=5, n=3000, shape_xyz=(200, 150, 10)):
    rng = np.random.default_rng(seed)
    v = np.array([0.1, 0.1, 0.2])
    ext = np.array([16249.93, 4614.7, 2742.24])
    img = rng.integers(0, 4000, size=shape_xyz[::-1]).astype(np.uint16)
    msk = rng.integers(0, 50, size=shape_xyz[::-1]).astype(np.uint16)
    idx = np.stack([rng.integers(0, s, n) for s in shape_xyz], axis=1)
    frac = rng.uniform(0.02, 0.98, size=(n, 3))
    pos = (ext + (idx + frac) * v).astype(np.float32).astype(float)       # Imaris stores float32
    exp = {"spot_channel": img[idx[:, 2], idx[:, 1], idx[:, 0]].astype(float),
           "mask_channel": msk[idx[:, 2], idx[:, 1], idx[:, 0]].astype(float)}
    geom = {"ext_min_xyz": ext, "voxel_xyz_um": v, "shape_xyz": np.array(shape_xyz)}
    return pos, exp, {"spot_channel": img, "mask_channel": msk}, geom


def test_frame_check_rule():
    assert I.text_rounding_halfwidth(16249.9) == pytest.approx(0.05)
    assert I.text_rounding_halfwidth(8150.1) == pytest.approx(0.005)
    assert I.text_rounding_halfwidth(-10921.1) == pytest.approx(0.05)
    assert I.text_rounding_halfwidth(2742.24) == pytest.approx(0.005)
    pos, exp, img, geom = _frame_case()
    r = I.frame_check(pos, geom, exp, img)
    assert r["chosen"]["name"] == "stored" and r["search"] is None
    assert r["candidates"][0]["channels"]["spot_channel"]["agreement"] == 1.0
    # a voxel read from six-significant-digit extents (0.10004 where the grid is 0.1): the far third of x
    # drifts, the stored frame fails there and the rounded voxel passes
    drift = dict(geom, voxel_xyz_um=np.array([0.10004, 0.1, 0.2]))
    r = I.frame_check(pos, drift, exp, img)
    stored = r["candidates"][0]["channels"]["spot_channel"]
    assert not r["candidates"][0]["passes"] and stored["thirds_x"][2] < 0.995 and stored["thirds_x"][0] > 0.995
    assert r["chosen"]["name"] == "rounded" and r["chosen"]["voxel_xyz_um"] == pytest.approx([0.1, 0.1, 0.2])
    # an ExtMin whose text lost 0.03 um: both fail, the search inside the +-0.05 um of the text finds it
    lost = dict(geom, ext_min_xyz=np.array([16249.9, 4614.7, 2742.24]))
    r = I.frame_check(pos, lost, exp, img)
    assert [c["passes"] for c in r["candidates"][:2]] == [False, False]
    assert r["chosen"]["name"].startswith("search") and r["chosen"]["ext_min_xyz"][0] == pytest.approx(16249.93, abs=0.004)
    assert r["search"]["rounding_halfwidth_um"][0] == pytest.approx(0.05)
    # values that belong to no frame: nothing passes, nothing is chosen
    rng = np.random.default_rng(1)
    junk = {k: rng.permutation(v) for k, v in exp.items()}
    assert I.frame_check(pos, geom, junk, img)["chosen"] is None
    # a spot within float32 precision of a voxel boundary may read the voxel across it
    one = np.array([[geom["ext_min_xyz"][0] + 50 * 0.1 - 1e-7, 4614.7 + 0.05, 2742.24 + 0.1]])
    vals = {"spot_channel": np.array([float(img["spot_channel"][0, 0, 50])])}
    a = I.frame_agreement(one, geom["ext_min_xyz"], geom["voxel_xyz_um"], geom["shape_xyz"], vals,
                          {"spot_channel": img["spot_channel"]})["channels"]["spot_channel"]
    assert a["agreement"] == 1.0


def test_frame_passes_needs_every_third_and_the_mask_channel():
    """The pass rule of entry 9 (d)(3) on hand-built results: a third of x or of y below 99.5 percent fails even
    when the overall agreement passes, and so does the mask channel below 99.5 percent."""
    good = {"agreement": 0.999, "thirds_x": [1.0, 1.0, 1.0], "thirds_y": [1.0, 1.0, None]}
    assert I.frame_passes({"channels": {"s": good}}, "s", None)
    assert not I.frame_passes({"channels": {"s": {**good, "thirds_x": [1.0, 1.0, 0.9]}}}, "s", None)
    assert not I.frame_passes({"channels": {"s": {**good, "thirds_y": [0.99, 1.0, 1.0]}}}, "s", None)
    assert not I.frame_passes({"channels": {"s": {**good, "agreement": 0.99}}}, "s", None)
    assert I.frame_passes({"channels": {"s": good, "m": {"agreement": 0.996}}}, "s", "m")
    assert not I.frame_passes({"channels": {"s": good, "m": {"agreement": 0.99}}}, "s", "m")


def test_the_rounded_voxel_is_rounded_to_four_decimals():
    """The real case of the 2026-01-05 males: a 0.108 um grid of 2720 voxels whose stored extents give 0.10801 on
    x. Rounded to 4 decimals the voxel is 0.108 and passes; rounded to 2 decimals (0.11) it would not."""
    rng = np.random.default_rng(7)
    shape_xyz = (2720, 20, 4)
    v = np.array([0.108, 0.108, 0.2])
    ext = np.array([2889.2, 4614.7, 2742.24])
    img = rng.integers(0, 4000, size=shape_xyz[::-1]).astype(np.uint16)
    n = 3000
    idx = np.stack([rng.integers(0, s, n) for s in shape_xyz], axis=1)
    pos = (ext + (idx + rng.uniform(0.02, 0.98, size=(n, 3))) * v).astype(np.float32).astype(float)
    exp = {"spot_channel": img[idx[:, 2], idx[:, 1], idx[:, 0]].astype(float)}
    geom = {"ext_min_xyz": ext, "voxel_xyz_um": np.array([0.10801, 0.108, 0.2]), "shape_xyz": np.array(shape_xyz)}
    r = I.frame_check(pos, geom, exp, {"spot_channel": img}, mask_key=None)
    assert not r["candidates"][0]["passes"]
    assert r["chosen"]["name"] == "rounded" and r["chosen"]["voxel_xyz_um"] == pytest.approx([0.108, 0.108, 0.2])


def test_export_with_the_frame_check(tmp_path):
    ims = tmp_path / "f.ims"
    build_ims(ims, ic_stats=True)
    plain = I.export_imaris_truth(ims, tmp_path / "plain", nd2_voxel_xyz_um=OUR_VOX)
    p = I.export_imaris_truth(ims, tmp_path / "fc", nd2_voxel_xyz_um=OUR_VOX, frame_check_on=True)
    fc = p["frame_check"]
    assert fc["chosen"]["name"] == "stored" and fc["image_source"] == "ims" and fc["spot_channel"] == 0
    assert fc["mask_channel"] == MASK_CH and fc["n_spots"] == len(SPOTS)
    assert fc["candidates"][0]["channels"]["spot_channel"]["agreement"] == 1.0
    assert fc["candidates"][0]["channels"]["mask_channel"]["agreement"] == 1.0
    a = pd.read_csv(I.truth_paths(tmp_path / "plain", "f")["spots"])
    b = pd.read_csv(I.truth_paths(tmp_path / "fc", "f")["spots"])
    pd.testing.assert_frame_equal(a, b)                    # the stored frame: the same truth
    assert p["spots_source"] == {"source": "ims"} and "frame_check" not in plain
    assert p["creation_parameter_summary"]["spots"]["source_channel_index"] == 0
    assert p["voxel_imaris_xyz_um_stored"] == pytest.approx(VOX.tolist())
    # an export without the statistic cannot be tested: no frame is chosen, the stored one is kept
    build_ims(tmp_path / "n.ims")
    q = I.export_imaris_truth(tmp_path / "n.ims", tmp_path / "n", nd2_voxel_xyz_um=OUR_VOX, frame_check_on=True)
    assert q["frame_check"]["chosen"] is None and "untestable" in q["frame_check"]


def _shifted_export(tmp_path, monkeypatch, ulps):
    """The synthetic project with its stored x extent 0.03 um low (as a six-significant-digit text above 10000 um
    can be); the positions and Intensity Center statistics stay those of the true frame."""
    ims = tmp_path / "d.ims"
    build_ims(ims, ic_stats=True)
    with h5py.File(ims, "r+") as f:
        img = f["DataSetInfo/Image"]
        img.attrs["ExtMin0"] = _s(f"{EXT_MIN[0] - 0.03:g}")
        img.attrs["ExtMax0"] = _s(f"{EXT_MAX[0] - 0.03:g}")
    monkeypatch.setattr(I, "text_rounding_halfwidth", lambda v, significant=6: 0.05)
    monkeypatch.setattr(I, "FLOAT32_ULPS", ulps)
    return I.export_imaris_truth(ims, tmp_path / "o", nd2_voxel_xyz_um=OUR_VOX, frame_check_on=True)


def _check_chosen_frame(tmp_path, p):
    chosen = p["frame_check"]["chosen"]
    assert chosen["name"] != "stored"
    assert p["ext_min_xyz_um"] == pytest.approx(chosen["ext_min_xyz"])
    assert p["ext_min_xyz_um_stored"][0] == pytest.approx(EXT_MIN[0] - 0.03)
    sp = pd.read_csv(I.truth_paths(tmp_path / "o", "d")["spots"]).set_index("spot_id")
    for sid, (v, *_rest) in SPOTS.items():
        assert tuple(int(sp.loc[sid, c]) for c in ("xi", "yi", "zi")) == v, sid
        assert sp.loc[sid, "x_um"] == pytest.approx(v[0] * OUR_VOX[0])
    nuc = pd.read_csv(I.truth_paths(tmp_path / "o", "d")["nuclei"])
    ext = np.asarray(chosen["ext_min_xyz"], float)
    vox = np.asarray(chosen["voxel_xyz_um"], float)
    idx = np.floor((nuc[["imaris_x_um", "imaris_y_um", "imaris_z_um"]].to_numpy(float) - ext) / vox).astype(int)
    assert (nuc[["xi", "yi", "zi"]].to_numpy() == idx).all()


def test_export_applies_a_non_stored_frame(tmp_path, monkeypatch):
    """The path that produced the four real truth tables of entry 9: every spot and surface is converted with the
    frame the check chose, not the stored one (no float32 allowance, so the search's choice is unambiguous)."""
    _check_chosen_frame(tmp_path, _shifted_export(tmp_path, monkeypatch, 0))


def test_search_does_not_stop_on_a_voxel_boundary(tmp_path, monkeypatch):
    """With the float32 allowance, ties on the tolerant agreement go to the higher strict agreement: the edge of
    the agreeing plateau would accept a spot on a voxel boundary and floor it into the wrong voxel."""
    _check_chosen_frame(tmp_path, _shifted_export(tmp_path, monkeypatch, 2))


def test_the_frame_check_and_xlsx_source_refuse_test_and_sealed_files(tmp_path):
    """Protocol changelog entry 10: with the frame check or the xlsx source, a test gonad's project and any sealed
    file are refused before anything is opened (the files need not exist, and no output folder is made); the
    pre-flight may open a test gonad's project with allow_test."""
    from germquant.benchmark.detector import TestSetGuardError

    out = tmp_path / "never"
    with pytest.raises(TestSetGuardError):
        I.export_imaris_truth(tmp_path / "20251105_n2_hs_MALE _001.ims", out, frame_check_on=True)
    with pytest.raises(TestSetGuardError):
        I.export_imaris_truth(tmp_path / "x.ims", out, nd2=tmp_path / "20251105_N2_nohs_HERM _002.nd2",
                              frame_check_on=True)
    with pytest.raises(PermissionError, match="sealed"):
        I.export_imaris_truth(tmp_path / "x.ims", out, spots_xlsx=tmp_path / "20251105_n2_hs_MALE_1.xlsx")
    with pytest.raises(PermissionError, match="sealed"):
        I.export_imaris_truth(tmp_path / "x.ims", out, nd2=tmp_path / "compiled_dataEP_flat.nd2", frame_check_on=True)
    assert not out.exists()
    test_like = tmp_path / "20251105_n2_hs_MALE _001.ims"
    build_ims(test_like, ic_stats=True)
    p = I.export_imaris_truth(test_like, tmp_path / "pre", nd2_voxel_xyz_um=OUR_VOX, frame_check_on=True,
                              allow_test=True)
    assert p["frame_check"]["chosen"]["name"] == "stored"


# ------------------------------------------------------------- the xlsx spot source (entry 9 (e))
def _write_export(path, ims_path, *, drop_spot=None, move_surface=0.0, radius=0.2, extra_spot=None,
                  imax_override=None):
    """A synthetic Imaris statistics export of the synthetic .ims: Overall, Position, Diameter, the Intensity
    Center of channels 1 and 2, the Intensity Max of channel 2 and Volume, in the export's layout.
    ``imax_override`` ({spot id: value}) gives a spot an Intensity Max of channel 2 other than its Intensity
    Center (a neighbour's fill value reached by the spot's volume)."""
    imax_override = dict(imax_override or {})
    import h5py
    import openpyxl

    with h5py.File(ims_path, "r") as f:
        raw = f["DataSet/ResolutionLevel 0/TimePoint 0/Channel 0/Data"][:Z, :Y, :X]
        mask = f["DataSet/ResolutionLevel 0/TimePoint 0/Channel 1/Data"][:Z, :Y, :X]
        smi = f["Scene8/Content/MegaSurfaces0/SurfaceModelInfo"][:]
    spots = [(sid, *_spot_position_um(v, fr), v) for sid, (v, fr, *_rest) in SPOTS.items() if sid != drop_spot]
    if extra_spot is not None:
        spots.append(extra_spot)
    surf = []
    for k, sid in enumerate(SURFACE_IDS):
        mval = BOXES[sid][1] if sid in BOXES else 0
        surf.append((sid, smi["CenterOfMassX"][k] + (move_surface if k == 0 else 0.0), smi["CenterOfMassY"][k],
                     smi["CenterOfMassZ"][k], float(smi["Volume"][k]), mval))
    wb = openpyxl.Workbook()
    wb.remove(wb.active)

    def sheet(name, header, rows):
        ws = wb.create_sheet(name)
        ws.append([name] + [""] * (len(header) - 1))
        ws.append(header)
        for r in rows:
            ws.append(list(r))

    sheet("Overall", ["Variable", "Value", "Unit", "Channel", "Class", "Image", "Name", "Surpass Object", "Time", "ID"],
          [("Number of Voxels", X * Y * Z, "", "", "", "Image 1", "", "Volume", "", ""),
           ("Total Number of Spots", len(spots), "", "", "", "", "", "Spots 1", "", ""),
           ("Total Number of Surfaces", len(surf), "", "", "", "", "", "Surfaces 1", "", "")])
    pos_h = ["Position X", "Position Y", "Position Z", "Unit", "Category", "Collection", "Name", "Surpass Object",
             "Time", "ID"]
    sheet("Position", pos_h, [(x, y, z, "um", "Spot", "Position", "", "Spots 1", 1, sid) for sid, x, y, z, _v in spots]
          + [(x, y, z, "um", "Surface", "Position", "", "Surfaces 1", 1, sid) for sid, x, y, z, _vol, _m in surf])
    sheet("Diameter", ["Diameter X", "Diameter Y", "Diameter Z", "Unit", "Category", "Collection", "Surpass Object",
                       "Time", "ID"],
          [(2 * radius, 2 * radius, 0.9, "um", "Spot", "Diameter", "Spots 1", 1, sid) for sid, *_r in spots])

    def at(arr, v):
        x, y, z = v
        return int(arr[z, y, x]) if 0 <= x < X and 0 <= y < Y and 0 <= z < Z else 0

    ic_h = ["Intensity Center", "Unit", "Category", "Channel", "Image", "Name", "Surpass Object", "Time", "ID"]
    sheet("Intensity Center Ch=1 Img=1", ic_h, [(at(raw, v), "", "Spot", 1, "Image 1", "", "Spots 1", 1, sid)
                                                for sid, _x, _y, _z, v in spots])
    sheet("Intensity Center Ch=2 Img=1", ic_h, [(at(mask, v), "", "Spot", 2, "Image 1", "", "Spots 1", 1, sid)
                                                for sid, _x, _y, _z, v in spots])
    mx_h = ["Intensity Max", "Unit", "Category", "Channel", "Image", "Surpass Object", "Time", "ID"]
    sheet("Intensity Max Ch=2 Img=1", mx_h, [(m, "", "Surface", 2, "Image 1", "Surfaces 1", 1, sid)
                                             for sid, _x, _y, _z, _vol, m in surf]
          + [(imax_override.get(sid, at(mask, v)), "", "Spot", 2, "Image 1", "Spots 1", 1, sid)
             for sid, _x, _y, _z, v in spots])
    sheet("Volume", ["Volume", "Unit", "Category", "Surpass Object", "Time", "ID"],
          [(vol, "um^3", "Surface", "Surfaces 1", 1, sid) for sid, _x, _y, _z, vol, _m in surf])
    wb.save(path)


def test_xlsx_spot_source_and_its_verification(tmp_path, monkeypatch):
    import hashlib

    from germquant.benchmark import spots as bsp

    pytest.importorskip("openpyxl")
    ims = tmp_path / "h.ims"
    build_ims(ims, ic_stats=True)
    # the synthetic nuclei are far below the protocol's 1000-voxel fragment size and two share a coloc value
    monkeypatch.setattr(bsp, "TRUTH_MASK_MIN_VOXELS", 1)
    monkeypatch.setattr(I, "TRUTH_ASSIGNMENT_MIN", 0.75)

    def run(name, **kw):
        x = tmp_path / f"{name}.xlsx"
        _write_export(x, ims, **kw)
        sha = hashlib.sha256(x.read_bytes()).hexdigest()
        n = len(SPOTS) - (1 if kw.get("drop_spot") else 0) + (1 if kw.get("extra_spot") else 0)
        expect = {"sha256": "0" * 64 if name == "badsha" else sha, "n_spots": n, "n_surfaces": len(SURFACE_IDS)}
        return I.export_imaris_truth(ims, tmp_path / name, nd2_voxel_xyz_um=OUR_VOX, spots_xlsx=x, xlsx_expect=expect)

    p = run("good", radius=0.2)
    ver = p["spots_source"]["verification"]
    assert {k: v["passes"] for k, v in ver.items() if k.startswith("V")} == {
        "V1_sha256": True, "V2_totals": True, "V3_surfaces": True, "V5_intensity_center": True, "V4_inside_extent": True,
        "V6_radius": True, "V7_truth_assignment": True}
    assert ver["passed"] is True and ver["V3_surfaces"]["matched_by"] == "ID"
    assert ver["V7_truth_assignment"]["agreement"] == pytest.approx(0.8)   # spot 13: coloc 12, its split label 13
    assert p["frame_check"]["chosen"]["name"] == "stored" and p["spots_source"]["frame_used"]["name"] == "stored"
    rep = p["spots_source"]["reported_not_a_criterion"]
    assert rep["n_ims_spots_saved"] == len(SPOTS) and rep["fraction_of_xlsx_spots_within_0.005_um_of_a_saved_ims_spot"] == 1.0
    spots = pd.read_csv(I.truth_paths(tmp_path / "good", "h")["spots"])
    assert list(spots.columns) == I.SPOT_COLUMNS + I.SPOT_EXTRA_COLUMNS and len(spots) == len(SPOTS)
    assert spots["spot_class"].isna().all() and spots["quality"].isna().all() and (spots["radius_um"] == 0.2).all()
    assert spots["spots_object"].eq("Spots 1 (xlsx)").all()
    by = spots.set_index("spot_id")
    for sid, (v, _fr, _c, _q, m, _r) in SPOTS.items():
        assert (int(by.loc[sid, "xi"]), int(by.loc[sid, "yi"]), int(by.loc[sid, "zi"])) == v
        assert int(by.loc[sid, "imaris_mask_label"]) == m and by.loc[sid, "x_um"] == pytest.approx(v[0] * OUR_VOX[0])
    # each criterion fails on its own defect, and then the gonad is nucleus truth only
    assert run("badsha")["spots_source"]["verification"]["V1_sha256"]["passes"] is False
    moved = run("moved", move_surface=0.02)["spots_source"]["verification"]
    assert moved["V3_surfaces"]["passes"] is False and moved["passed"] is False
    assert run("radius", radius=0.35)["spots_source"]["verification"]["V6_radius"]["passes"] is False
    out = run("outside", extra_spot=(99, EXT_MIN[0] - 0.5, EXT_MIN[1] + 0.05, EXT_MIN[2] + 0.1, (-5, 0, 0)))
    ver = out["spots_source"]["verification"]
    assert ver["V4_inside_extent"]["n_outside"] == 1 and ver["passed"] is False
    assert any("nucleus truth only" in w for w in out["warnings"])
    # the coloc value is the export's Intensity Max of the mask channel, not its Intensity Center: spot 14 (outside
    # every nucleus, Intensity Center 0) reaches the fill value 7 of a neighbour
    p = run("imax", imax_override={14: 7})
    by = pd.read_csv(I.truth_paths(tmp_path / "imax", "h")["spots"]).set_index("spot_id")
    assert int(by.loc[14, "imaris_mask_label"]) == 7 and int(by.loc[11, "imaris_mask_label"]) == 7
    assert p["spots_source"]["verification"]["V5_intensity_center"]["passes"] is True


def test_each_xlsx_criterion_fails_on_its_own(tmp_path, monkeypatch):
    """V2, the volume part of V3, V5 and V7 (V3 volume and V7 decided 20260105_n2_hs_1_MALE_1) each fail alone on
    a defect of their own."""
    import hashlib

    from germquant.benchmark import spots as bsp
    from germquant.validate import imaris_xlsx as X

    pytest.importorskip("openpyxl")
    monkeypatch.setattr(bsp, "TRUTH_MASK_MIN_VOXELS", 1)
    ims = tmp_path / "e.ims"
    build_ims(ims, ic_stats=True)
    x = tmp_path / "e.xlsx"
    _write_export(x, ims)
    sha = hashlib.sha256(x.read_bytes()).hexdigest()
    real = X.read_export_objects

    def run(name, mutate=None, expect=None, vmin=0.75):
        monkeypatch.setattr(I, "TRUTH_ASSIGNMENT_MIN", vmin)

        def fake(*a, **k):
            xl = real(*a, **k)
            if mutate:
                mutate(xl)
            return xl

        monkeypatch.setattr(X, "read_export_objects", fake)
        e = {"sha256": sha, "n_spots": len(SPOTS), "n_surfaces": len(SURFACE_IDS), **(expect or {})}
        v = I.export_imaris_truth(ims, tmp_path / name, nd2_voxel_xyz_um=OUR_VOX, spots_xlsx=x,
                                  xlsx_expect=e)["spots_source"]["verification"]
        return {k for k, r in v.items() if k.startswith("V") and not r["passes"]}, v["passed"]

    assert run("ok") == (set(), True)
    assert run("v2", expect={"n_spots": len(SPOTS) + 1}) == ({"V2_totals"}, False)

    def vol(xl):
        xl["surfaces"]["volume"] = xl["surfaces"]["volume"] * (1 + 5e-5)

    assert run("v3", mutate=vol) == ({"V3_surfaces"}, False)

    def ic(xl):
        xl["spots"]["ic_spot"] = xl["spots"]["ic_spot"] + 1

    assert run("v5", mutate=ic) == ({"V5_intensity_center"}, False)
    assert run("v7", vmin=0.9935) == ({"V7_truth_assignment"}, False)


def test_addition_image_check_compares_planes(tmp_path):
    """scripts/check_addition_images.py, entry 9 (d)(2): the .ims planes against the .nd2 planes."""
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "check_addition_images", Path(__file__).resolve().parents[1] / "scripts" / "check_addition_images.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    ims = tmp_path / "p.ims"
    build_ims(ims)
    with h5py.File(ims, "r") as f:
        raw = f["DataSet/ResolutionLevel 0/TimePoint 0/Channel 0/Data"][:Z, :Y, :X]
    r = mod.compare_planes(ims, {0: raw.copy()}, [0, 5, 99])
    assert [(p["plane"], p["equal"]) for p in r] == [(0, True), (5, True)]          # plane 99 is outside
    other = raw.copy()
    other[5, 3, 4] += 1
    r = mod.compare_planes(ims, {0: other}, [5])
    assert r == [{"channel": 0, "plane": 5, "equal": False, "n_different": 1}]
    assert mod.compare_planes(ims, {0: raw[:, :10]}, [0])[0]["equal"] is False
