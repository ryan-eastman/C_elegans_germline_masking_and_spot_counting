"""The tracking exchange and the biologist validation (germquant.timelapse.track / track_validate,
docs/TIMELAPSE.md "Tracking" and "Validating the linking").

A synthetic five-timepoint series is written straight into a results folder the way process_series lays one
out (t<NNN>/<stem>_t<NNN>__nuclei.csv, __nuclei_labels.tif, __image_summary.csv, __ref_stride2.npy, the series
record, a register_shifts table): spherical nuclei with known identities, a division-free flow of one voxel per
timepoint along x, a planted drift of (0, 1, 2) voxels per timepoint, one nucleus that drops out of the
segmentation for one timepoint, two off-gonad nuclei, and per-timepoint label ids shuffled so that nucleus_id
never equals the identity. Everything is checked against those identities. The fixed-image path is not
touched (test_golden runs unchanged); the last test runs process_series on the three-timepoint gonad of
test_timelapse_series and links its real file set.
"""
from __future__ import annotations

import json
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tifffile
from test_timelapse_series import ND2, STEM as GONAD_STEM, _cfg, patched_series_reader  # noqa: F401  (fixture)

from germquant import cli, schema
from germquant.timelapse import (
    export_cellacdc,
    find_series,
    import_tracks,
    link_series,
    load_shifts,
    object_table,
    pick_validation_nuclei,
    process_series,
    track_strips,
    validate_tracks,
)
from germquant.timelapse import track as TR
from germquant.timelapse.track_validate import format_verdict, read_hand_points

ROOT = Path(__file__).resolve().parents[1]
STEM = "20261005_dlw93_nohs_herm_001"
SPACING = (0.4, 0.2, 0.2)
SHAPE = (16, 96, 176)                 # (Z, Y, X) voxels
T = 5
PERIOD_S = 300.0
DRIFT_PX = np.array([0, 1, 2])        # (z, y, x) voxels of drift per timepoint (the sample moves by this)
FLOW_PX = 1                           # +x voxels per timepoint for every germline nucleus (division-free flow)
R_UM = 1.0                            # sphere radius
DROPOUT = (6, 2)                      # identity 6 is missing from the segmentation at t=2
AXIS_ANCHOR_PX = 10.0                 # the distal end of the axis, in registered x voxels
AXIS_LENGTH_UM = 40.0                 # total_germline_length_um written to every image_summary
REF_STRIDE = 2


def _germline_bases() -> dict[int, np.ndarray]:
    """identity -> t0 centre (z, y, x) in voxels of the 18 germline nuclei (two rows of 8, two at the edge)."""
    bases = {}
    ident = 1
    for y in (40, 56):
        for x in range(20, 20 + 8 * 14, 14):
            bases[ident] = np.array([8, y, x], float)
            ident += 1
    for x in (60, 100):                                  # 1.6 um from the top border: the edge nuclei
        bases[ident] = np.array([8, 8, x], float)
        ident += 1
    return bases


JUNK_BASES = {101: np.array([8, 84, 140], float), 102: np.array([8, 84, 156], float)}   # off-gonad, no flow


def raw_centre(identity: int, k: int) -> np.ndarray:
    """The raw-frame centre (voxels) of an identity at timepoint k: base + flow + drift (junk: drift only)."""
    if identity in JUNK_BASES:
        return JUNK_BASES[identity] + k * DRIFT_PX
    return _germline_bases()[identity] + np.array([0, 0, k * FLOW_PX]) + k * DRIFT_PX


def registered_centre_um(identity: int, k: int) -> np.ndarray:
    return (raw_centre(identity, k) - k * DRIFT_PX) * np.array(SPACING)


def _sphere(lab: np.ndarray, centre: np.ndarray, value: int, zz, yy, xx) -> None:
    d2 = (((zz - centre[0]) * SPACING[0]) ** 2 + ((yy - centre[1]) * SPACING[1]) ** 2 + ((xx - centre[2]) * SPACING[2]) ** 2)
    lab[d2 <= R_UM ** 2] = value


def write_synthetic_series(root: Path, *, flip_at: int | None = None, shifts: bool = True) -> dict:
    """Write the series folder; return the truth: ``ids[k] = {nucleus_id: identity}`` per timepoint."""
    from skimage.measure import regionprops

    out = root / STEM
    zz, yy, xx = np.indices(SHAPE).astype(np.float32)
    rng = np.random.default_rng(2026)
    ids: dict[int, dict[int, int]] = {}
    series_tps = []
    for k in range(T):
        present = [i for i in _germline_bases() if not (i == DROPOUT[0] and k == DROPOUT[1])] + list(JUNK_BASES)
        order = rng.permutation(len(present))              # nucleus_id never equals the identity
        nid_of = {present[j]: int(order[j]) + 1 for j in range(len(present))}
        lab = np.zeros(SHAPE, np.int32)
        ref = np.zeros(SHAPE, np.float32)
        for ident in present:
            c = raw_centre(ident, k)
            _sphere(lab, c, nid_of[ident], zz, yy, xx)
            _sphere(ref, c, 1000.0 + 10.0 * ident - 50.0 * k, zz, yy, xx)       # bleaching: -50 per timepoint
        sub = out / f"t{k:03d}"
        sub.mkdir(parents=True)
        iid = f"{STEM}_t{k:03d}"
        tifffile.imwrite(sub / f"{iid}__nuclei_labels.tif", lab, compression="zlib")
        np.save(sub / f"{iid}__ref_stride{REF_STRIDE}.npy", ref[:, ::REF_STRIDE, ::REF_STRIDE].astype(np.uint16))
        rows = []
        ident_of = {v: i for i, v in nid_of.items()}
        for rp in regionprops(lab):
            ident = ident_of[int(rp.label)]
            cz, cy, cx = (np.asarray(rp.centroid) * np.array(SPACING)).tolist()
            germ = ident not in JUNK_BASES
            x_reg_vox = raw_centre(ident, k)[2] - k * DRIFT_PX[2]
            pos = (x_reg_vox - AXIS_ANCHOR_PX) * SPACING[2]
            if flip_at is not None and k == flip_at:
                pos = AXIS_LENGTH_UM - pos                      # the axis stage's orientation flipped at this timepoint
            rows.append({
                "nucleus_id": int(rp.label), "volume_um3": float(rp.area) * np.prod(SPACING), "n_voxels": int(rp.area),
                "centroid_z_um": cz, "centroid_y_um": cy, "centroid_x_um": cx,
                "in_germline": germ, "axis_position_norm": (pos / AXIS_LENGTH_UM) if germ else np.nan,
                "axis_position_um": pos if germ else np.nan, "n_spots": 0,
                "central_element_mean_intensity": 1000.0 + 10.0 * ident - 50.0 * k,
                "sc_total_length_um": (20.0 + ident + 0.5 * k) if germ else np.nan,
            })
        nuclei = pd.DataFrame(rows).sort_values("nucleus_id")
        for col in schema.NUCLEI:
            if col not in nuclei.columns:
                nuclei[col] = pd.NA
        nuclei = nuclei[list(schema.NUCLEI) + [c for c in nuclei.columns if c not in schema.NUCLEI]]
        meta = {"image_id": iid, "file_path": "x.nd2", "genotype": "DLW93", "sex": "herm", "germ_cell": "oocyte",
                "treatment": "noHS", "replicate": "001", "acquisition_date": "20261005",
                "voxel_dz_um": SPACING[0], "voxel_dy_um": SPACING[1], "voxel_dx_um": SPACING[2],
                "n_channels": 2, "dapi_present": False, "channel_map_json": json.dumps({"central_element": 1, "brightfield": 0}),
                "pipeline_version": "0.1.0", "git_sha": "test", "config_hash": "cafe", "run_timestamp": "2026-10-09T00:00:00+00:00"}
        for c, v in meta.items():
            nuclei[c] = v
        nuclei.to_csv(sub / f"{iid}__nuclei.csv", index=False)
        pd.DataFrame([{"n_nuclei": len(nuclei), "n_germline_nuclei": int(nuclei["in_germline"].sum()), "mean_spots": 0.0,
                       "total_germline_length_um": AXIS_LENGTH_UM, "qc_pass": True, "qc_flags": "", **meta}]).to_csv(
            sub / f"{iid}__image_summary.csv", index=False)
        (sub / f"{iid}__done.json").write_text(json.dumps({"image_id": iid, "config_hash": "cafe", "t_index": k,
                                                           "t_seconds": k * PERIOD_S, "failed_stages": []}))
        ids[k] = ident_of
        series_tps.append({"t_index": k, "t_seconds": k * PERIOD_S, "image_id": iid, "out_dir": sub.name, "status": "done"})
    record = {"series_id": STEM, "file_path": "x.nd2", "config_hash": "cafe", "pipeline_version": "0.1.0", "git_sha": "test",
              "n_timepoints": T, "period_s": PERIOD_S, "t_seconds": [k * PERIOD_S for k in range(T)],
              "requested": list(range(T)), "done": list(range(T)), "failed": [], "timepoints": series_tps, "complete": True}
    (out / f"{STEM}__series.json").write_text(json.dumps(record, indent=1))
    (out / "run_manifest.json").write_text(json.dumps({"config": {"timelapse": {"track": {"max_dist_um": 2.5}}}}))
    if shifts:
        sh = pd.DataFrame({"t_index": range(T),
                           "dz_um": [k * DRIFT_PX[0] * SPACING[0] for k in range(T)],
                           "dy_um": [k * DRIFT_PX[1] * SPACING[1] for k in range(T)],
                           "dx_um": [k * DRIFT_PX[2] * SPACING[2] for k in range(T)],
                           "step_dz_um": 0.0, "step_dy_um": DRIFT_PX[1] * SPACING[1], "step_dx_um": DRIFT_PX[2] * SPACING[2],
                           "error": 0.0, "quality_before_um": 0.45, "quality_after_um": 0.2, "capped": False,
                           "brightfield_dz_um": np.nan, "brightfield_dy_um": np.nan, "brightfield_dx_um": np.nan,
                           "brightfield_disagreement_um": np.nan, "method": "synthetic"})
        sh.to_csv(out / f"{STEM}__register_shifts.csv", index=False)
    return {"out": out, "ids": ids}


@pytest.fixture(scope="module")
def series(tmp_path_factory):
    return write_synthetic_series(tmp_path_factory.mktemp("series"))


def _identity_of_points(points: pd.DataFrame, ids: dict) -> pd.Series:
    return pd.Series([ids[int(t)][int(n)] for t, n in zip(points["t_index"], points["nucleus_id"])], index=points.index)


def _check_identities(points: pd.DataFrame, ids: dict, *, split_identities=()):
    """Zero switches: every track holds one identity; every identity has one track (two for the split ones)."""
    linked = points[points["time_track_id"].notna()].copy()
    linked["identity"] = _identity_of_points(linked, ids)
    per_track = linked.groupby("time_track_id")["identity"].nunique()
    assert (per_track == 1).all(), per_track[per_track != 1]
    per_identity = linked.groupby("identity")["time_track_id"].nunique()
    for ident, n in per_identity.items():
        assert n == (2 if ident in split_identities else 1), (ident, n)
    return linked


# ---- the series folder -------------------------------------------------------------------------------

def test_layout_shifts_and_registered_coordinates(series):
    layout = find_series(series["out"])
    assert layout.stem == STEM and layout.timepoints == list(range(T)) and layout.missing == []
    assert layout.spacing() == SPACING and layout.period_s == PERIOD_S and layout.t_seconds_of(3) == 3 * PERIOD_S
    assert layout.minutes_of(2) == 10.0 and layout.ref_path(1)[1] == REF_STRIDE and layout.channel_index("central_element") == 1
    shifts, source = load_shifts(layout)
    assert source == "register_shifts" and list(shifts["t_index"]) == list(range(T))
    assert shifts["dx_um"].tolist() == pytest.approx([k * DRIFT_PX[2] * SPACING[2] for k in range(T)])
    obj = object_table(layout, shifts)                              # germline only
    assert len(obj) == 18 * T - 1 and "in_germline" in obj.columns and obj["in_germline"].all()
    assert {"central_element_mean_intensity", "sc_total_length_um", "axis_position_um", "volume_um3"} <= set(obj.columns)
    for r in obj.itertuples(index=False):
        ident = series["ids"][int(r.t_index)][int(r.nucleus_id)]
        np.testing.assert_allclose([r.zr_um, r.yr_um, r.xr_um], registered_centre_um(ident, int(r.t_index)), atol=1e-6)
        np.testing.assert_allclose([r.z_um, r.y_um, r.x_um], raw_centre(ident, int(r.t_index)) * np.array(SPACING), atol=1e-6)
    everything = object_table(layout, shifts, germline_only=False)
    assert len(everything) == len(obj) + 2 * T
    # without a shift table the coordinates are taken as aligned (zero shifts), with a note
    bare = write_synthetic_series(series["out"].parent / "bare", shifts=False)
    sh0, src0 = load_shifts(find_series(bare["out"]))
    assert src0 == "none" and (sh0[["dz_um", "dy_um", "dx_um"]] == 0).all().all()


# ---- track-export ------------------------------------------------------------------------------------

@pytest.mark.parametrize("stride", [1, 2])
def test_export_round_trips_through_the_npz(series, tmp_path, stride):
    from skimage.measure import regionprops

    out = tmp_path / f"export{stride}"
    rec = export_cellacdc(series["out"], out, stride=stride, crop="germline", image="ref")
    images = out / "Position_1" / "Images"
    assert sorted(p.name for p in images.iterdir()) == [f"{STEM}_central_element.tif", f"{STEM}_metadata.csv", f"{STEM}_segm.npz"]
    with np.load(images / f"{STEM}_segm.npz") as arc:
        segm = arc[arc.files[0]]
    assert segm.ndim == 4 and segm.shape[0] == T and segm.dtype == np.uint16 and list(segm.shape) == rec["shape"]
    img = tifffile.imread(images / f"{STEM}_central_element.tif")
    assert img.shape == segm.shape and img.dtype == np.uint16 and img.max() > 0
    meta = pd.read_csv(images / f"{STEM}_metadata.csv").set_index("Description")["values"]
    assert int(meta["SizeT"]) == T and int(meta["SizeZ"]) == segm.shape[1] and float(meta["PhysicalSizeX"]) == pytest.approx(SPACING[2] * stride)
    assert float(meta["TimeIncrement"]) == PERIOD_S and meta["basename"] == f"{STEM}_"
    assert rec["stride"] == stride and rec["shift_source"] == "register_shifts" and rec["crop"] == "germline"
    for k in range(T):
        assert rec["shifts"][str(k)]["voxels"] == [int(v) for v in k * DRIFT_PX]
    z0, y0, x0 = rec["box"]["start"]
    sp = np.array(SPACING)
    for k in range(T):
        frame = segm[k]
        values = {int(v) for v in np.unique(frame) if v}
        germ_ids = {nid for nid, ident in series["ids"][k].items() if ident not in JUNK_BASES}
        assert values == germ_ids                                       # exactly the germline nuclei, by nucleus_id
        assert {int(a) for a in rec["label_maps"][str(k)]} == germ_ids
        for rp in regionprops(frame):
            ident = series["ids"][k][int(rp.label)]
            cz, cy, cx = rp.centroid
            reg_vox = np.array([z0 + cz, y0 + stride * cy, x0 + stride * cx])
            tol = 1e-6 if stride == 1 else stride * 0.6                 # sub-sampling moves a centroid by under a stride
            np.testing.assert_allclose(reg_vox * sp, registered_centre_um(ident, k), atol=tol * sp.max() + 1e-6)
    # the union box holds every germline nucleus with the flow, and nothing of the off-gonad junk
    assert rec["box"]["stop"][2] <= int(raw_centre(16, T - 1)[2] - (T - 1) * DRIFT_PX[2]) + 6
    assert rec["box"]["start"][1] <= 8 - 5
    # crop=all keeps the junk and the whole (shifted) field
    rec_all = export_cellacdc(series["out"], tmp_path / f"all{stride}", stride=stride, crop="all", image="none")
    with np.load(tmp_path / f"all{stride}" / "Position_1" / "Images" / f"{STEM}_segm.npz") as arc:
        segm_all = arc[arc.files[0]]
    assert rec_all["box"]["shape_raw"] == list(SHAPE) and rec_all["image"]["source"] == "none"
    assert {int(v) for v in np.unique(segm_all[0]) if v} == set(series["ids"][0])
    assert not (tmp_path / f"all{stride}" / "Position_1" / "Images" / f"{STEM}_central_element.tif").exists()


# ---- track-link --------------------------------------------------------------------------------------

@pytest.mark.parametrize("backend", ["overlap_inhouse", "overlap"])
def test_overlap_linking_recovers_identities(series, backend):
    if backend == "overlap" and not TR.cellacdc_overlap_available()[0]:
        pytest.skip("Cell-ACDC's tracker is not importable here")
    res = link_series(series["out"], backend=backend, stride=1)
    rec = res["record"]
    assert rec["backend"] == backend and rec["note"] == "" and rec["shift_source"] == "register_shifts"
    points, tracks = res["points"], res["tracks"]
    # zero switches; the dropout nucleus is two tracks (overlap cannot bridge a missing frame)
    linked = _check_identities(points, series["ids"], split_identities=(DROPOUT[0],))
    assert len(linked) == 18 * T - 1 and (linked["link_source"] == backend).all()
    assert rec["n_tracks"] == 19 and rec["n_full_length"] == 17
    assert list(points.columns[:len(schema.TRACK_POINTS)]) == list(schema.TRACK_POINTS) or set(schema.TRACK_POINTS) <= set(points.columns)
    assert (points["zr_um"] == points["reg_z_um"]).all()
    # steps are the flow alone (one voxel in x per timepoint), not the drift
    steps = linked[linked["step_um"].notna()]
    assert steps["step_um"].to_numpy() == pytest.approx(FLOW_PX * SPACING[2], abs=1e-6)
    assert steps["step_x_um"].to_numpy() == pytest.approx(FLOW_PX * SPACING[2], abs=1e-6)
    assert (steps["step_y_um"].abs() < 1e-6).all() and not linked["is_gap"].any()
    # the first timepoint's track ids are its nucleus ids
    t0 = linked[linked["t_index"] == 0]
    assert (t0["time_track_id"] == t0["nucleus_id"]).all()
    # the tracks table: flow speed and axis velocity, bleaching slopes, SC length change
    full = tracks[tracks["is_full_length"]]
    assert len(full) == 17 and (full["n_points"] == T).all() and (full["n_gaps"] == 0).all()
    assert (full["t_first"] == 0).all() and (full["t_last"] == T - 1).all() and (full["t_start"] == full["t_first"]).all()
    assert full["speed_um_per_min"].to_numpy() == pytest.approx(FLOW_PX * SPACING[2] / (PERIOD_S / 60), rel=1e-6)
    assert full["axis_velocity_um_per_min"].to_numpy() == pytest.approx(FLOW_PX * SPACING[2] / (PERIOD_S / 60), rel=1e-6)
    assert full["intensity_slope_raw_per_min"].to_numpy() == pytest.approx(-50.0 / (PERIOD_S / 60), rel=1e-6)
    assert full["intensity_slope_raw"].to_numpy() == pytest.approx(full["intensity_slope_raw_per_min"].to_numpy())
    assert np.abs(full["intensity_slope_norm_per_min"].to_numpy()).max() < 0.01          # median-normalised: bleaching removed
    assert full["sc_length_change_um"].to_numpy() == pytest.approx(0.5 * (T - 1), abs=1e-9)
    assert full["mean_step_um"].to_numpy() == pytest.approx(FLOW_PX * SPACING[2], abs=1e-6)
    assert (full["path_length_um"] - full["net_displacement_um"]).abs().max() < 1e-6
    assert set(rec["axis_sign_by_t"].values()) == {1}
    # the files
    tp = pd.read_csv(series["out"] / f"{STEM}__track_points.csv")
    tt = pd.read_csv(series["out"] / f"{STEM}__tracks.csv")
    assert list(tp.columns[:len(schema.TRACK_POINTS)]) == list(schema.TRACK_POINTS)
    assert list(tp.columns[-len(schema.SERIES_META):]) == list(schema.SERIES_META)
    assert list(tt.columns[:len(schema.TRACKS)]) == list(schema.TRACKS)
    assert (tp["series_id"] == STEM).all() and (tp["sex"] == "herm").all() and (tp["n_timepoints"] == T).all()
    assert {"volume_um3", "in_germline", "axis_position_um", "central_element_mean_intensity", "sc_total_length_um"} <= set(tp.columns)
    assert len(tp) == 18 * T - 1 and len(tt) == 19
    assert json.loads((series["out"] / f"{STEM}__track_link.json").read_text())["backend"] == backend


def test_overlap_assign_rule():
    """The rule itself: IoA against the previous object's volume, best current object wins, others are new."""
    prev = np.zeros((1, 10, 10), np.int64)
    prev[0, 1:4, 1:4] = 5                        # 9 voxels
    prev[0, 6:9, 6:9] = 7
    cur = np.zeros((1, 10, 10), np.int64)
    cur[0, 2:5, 2:5] = 1                         # overlaps 5 by 4 voxels: IoA 4/9 = 0.44
    cur[0, 1:2, 1:2] = 2                         # overlaps 5 by 1 voxel: IoA 0.11 -> new
    cur[0, 6:9, 6:9] = 3                         # 7 exactly
    assign, nxt = TR.overlap_assign(prev, cur, min_ioa=0.4, next_id=8)
    assert assign == {1: 5, 3: 7, 2: 8} and nxt == 9
    assign, _ = TR.overlap_assign(prev, cur, min_ioa=0.5, next_id=8)
    assert assign == {3: 7, 1: 8, 2: 9}


def test_laptrack_bridges_the_dropout(series):
    res = link_series(series["out"], backend="laptrack", gap=2, max_dist_um=3.0)
    points, tracks = res["points"], res["tracks"]
    linked = _check_identities(points, series["ids"])                        # one track per identity, the dropout bridged
    assert res["record"]["n_tracks"] == 18 and res["record"]["backend_info"]["n_splits"] == 0
    linked["identity"] = _identity_of_points(linked, series["ids"])
    drop = linked[linked["identity"] == DROPOUT[0]].sort_values("t_index")
    assert list(drop["t_index"]) == [0, 1, 3, 4] and list(drop["is_gap"]) == [False, False, True, False]
    assert drop["step_um"].iloc[2] == pytest.approx(2 * FLOW_PX * SPACING[2], abs=1e-6)
    row = tracks[tracks["time_track_id"] == drop["time_track_id"].iloc[0]].iloc[0]
    assert row["n_points"] == 4 and row["n_gaps"] == 1 and not row["is_full_length"]
    assert tracks["is_full_length"].sum() == 17
    # without gap closing the dropout nucleus is two tracks again
    res0 = link_series(series["out"], backend="laptrack", gap=0, write=False)
    _check_identities(res0["points"], series["ids"], split_identities=(DROPOUT[0],))
    assert res0["record"]["n_tracks"] == 19
    # the config in the series' run_manifest.json supplies max_dist_um when the flag is absent
    res_cfg = link_series(series["out"], backend="laptrack", write=False,
                          cfg=json.loads((series["out"] / "run_manifest.json").read_text())["config"])
    assert res_cfg["record"]["backend_info"]["max_dist_um"] == 2.5


def test_none_backend_writes_points_without_tracks(series):
    res = link_series(series["out"], backend="none")
    assert res["points"]["time_track_id"].isna().all() and (res["points"]["link_source"] == "none").all()
    assert len(res["tracks"]) == 0 and list(res["tracks"].columns) == list(schema.TRACKS)
    tp = pd.read_csv(series["out"] / f"{STEM}__track_points.csv")
    assert len(tp) == 18 * T - 1 and tp["time_track_id"].isna().all() and tp["central_element_mean_intensity"].notna().all()
    assert len(pd.read_csv(series["out"] / f"{STEM}__tracks.csv")) == 0


def test_axis_sign_is_aligned_to_t0(tmp_path):
    """A timepoint whose axis orientation flipped (p -> L - p) is recognised and un-flipped before velocities."""
    flipped = write_synthetic_series(tmp_path / "flip", flip_at=3)       # t=3, off centre, so a raw fit would feel it
    res = link_series(flipped["out"], backend="overlap_inhouse", stride=1)
    signs = res["record"]["axis_sign_by_t"]
    assert signs == {"0": 1, "1": 1, "2": 1, "3": -1, "4": 1}
    full = res["tracks"][res["tracks"]["is_full_length"]]
    assert full["axis_velocity_um_per_min"].to_numpy() == pytest.approx(FLOW_PX * SPACING[2] / (PERIOD_S / 60), rel=1e-6)
    pts = res["points"]
    at3 = pts[(pts["t_index"] == 3) & pts["time_track_id"].notna()]
    np.testing.assert_allclose(at3["axis_position_aligned_um"].to_numpy(float),
                               AXIS_LENGTH_UM - at3["axis_position_um"].to_numpy(float), atol=1e-9)
    at2 = pts[(pts["t_index"] == 2) & pts["time_track_id"].notna()]
    np.testing.assert_allclose(at2["axis_position_aligned_um"].to_numpy(float), at2["axis_position_um"].to_numpy(float), atol=1e-9)
    # with the raw column the flipped frame would have given a nonsense velocity for the distal tracks
    raw = pts[pts["time_track_id"].notna()].groupby("time_track_id").apply(
        lambda g: np.polyfit(g["t_index"], g["axis_position_um"], 1)[0] if len(g) > 2 else np.nan, include_groups=False)
    assert raw.dropna().abs().max() > 1.0


# ---- track-import ------------------------------------------------------------------------------------

def _truth_points(ids: dict, identities=None) -> pd.DataFrame:
    rows = []
    for k in range(T):
        for nid, ident in ids[k].items():
            if ident in JUNK_BASES or (identities is not None and ident not in identities):
                continue
            z, y, x = raw_centre(ident, k) * np.array(SPACING)
            rows.append({"track_id": ident, "t_index": k, "z_um": z, "y_um": y, "x_um": x, "nucleus_id": nid})
    return pd.DataFrame(rows)


def test_import_trackmate_imaris_and_cellacdc(series, tmp_path):
    truth = _truth_points(series["ids"])
    # TrackMate spots CSV with the three description rows TrackMate writes under the header
    cols = ["LABEL", "ID", "TRACK_ID", "QUALITY", "POSITION_X", "POSITION_Y", "POSITION_Z", "POSITION_T", "FRAME", "RADIUS", "VISIBILITY"]
    tm = pd.DataFrame({"LABEL": [f"ID{i}" for i in range(len(truth))], "ID": range(len(truth)), "TRACK_ID": truth["track_id"],
                       "QUALITY": 1.0, "POSITION_X": truth["x_um"], "POSITION_Y": truth["y_um"], "POSITION_Z": truth["z_um"],
                       "POSITION_T": truth["t_index"] * PERIOD_S, "FRAME": truth["t_index"], "RADIUS": 1.0, "VISIBILITY": 1})
    extra = pd.DataFrame([["Label", "Spot ID", "Track ID", "Quality", "X", "Y", "Z", "T", "Frame", "Radius", "Visibility"],
                          ["Label", "Spot ID", "Track ID", "Quality", "X", "Y", "Z", "T", "Frame", "R", "Visibility"],
                          ["", "", "", "(quality)", "(micron)", "(micron)", "(micron)", "(sec)", "", "(micron)", ""]], columns=cols)
    csv = tmp_path / "spots.csv"
    pd.concat([extra, tm[cols].astype(str)], ignore_index=True).to_csv(csv, index=False)
    res = import_tracks(series["out"], "trackmate", csv)
    linked = _check_identities(res["points"], series["ids"])
    assert (linked["link_source"] == "imported:trackmate").all() and (linked["time_track_id"] == linked.pipe(_identity_of_points, series["ids"])).all()
    assert res["record"]["n_tracks"] == 18 and res["record"]["backend_info"]["n_outside_any_label"] == 0
    assert res["points"]["time_track_id"].notna().all()
    # the same tracks from the TrackMate model XML
    root = ET.Element("TrackMate")
    model = ET.SubElement(root, "Model", spatialunits="micron")
    allspots = ET.SubElement(model, "AllSpots")
    for k in range(T):
        sif = ET.SubElement(allspots, "SpotsInFrame", frame=str(k))
        for i, r in truth[truth["t_index"] == k].iterrows():
            ET.SubElement(sif, "Spot", ID=str(i), FRAME=str(k), POSITION_X=str(r.x_um), POSITION_Y=str(r.y_um), POSITION_Z=str(r.z_um))
    alltracks = ET.SubElement(model, "AllTracks")
    for tid, g in truth.groupby("track_id"):
        tr = ET.SubElement(alltracks, "Track", TRACK_ID=str(tid))
        idx = g.sort_values("t_index").index.tolist()
        for a, b in zip(idx[:-1], idx[1:]):
            ET.SubElement(tr, "Edge", SPOT_SOURCE_ID=str(a), SPOT_TARGET_ID=str(b))
    xml = tmp_path / "model.xml"
    ET.ElementTree(root).write(xml)
    res_xml = import_tracks(series["out"], "trackmate", xml, write=False)
    pd.testing.assert_frame_equal(res_xml["points"][["t_index", "nucleus_id", "time_track_id"]].reset_index(drop=True),
                                  res["points"][["t_index", "nucleus_id", "time_track_id"]].reset_index(drop=True))
    # an Imaris points table with 1-based time, two identities only; the other nuclei stay in the table unlinked
    im = _truth_points(series["ids"], identities=(1, 2))
    imaris = tmp_path / "imaris.csv"
    pd.DataFrame({"id": im["track_id"], "t": im["t_index"] + 1, "z": im["z_um"], "y": im["y_um"], "x": im["x_um"]}).to_csv(imaris, index=False)
    res_im = import_tracks(series["out"], "imaris", imaris, time_base=1, write=False)
    pts = res_im["points"]
    assert pts["time_track_id"].notna().sum() == 2 * T and set(pts["link_source"]) == {"imported:imaris", "none"}
    assert len(pts) == 18 * T - 1
    _check_identities(pts, series["ids"])
    # a Cell-ACDC tracked segm.npz on the export grid: the exported labels relabelled by identity, one object deleted
    out = tmp_path / "exp"
    rec = export_cellacdc(series["out"], out, stride=2, image="none")
    npz = out / "Position_1" / "Images" / f"{STEM}_segm.npz"
    with np.load(npz) as arc:
        segm = arc[arc.files[0]].astype(np.int64)
    tracked = np.zeros_like(segm)
    for i, k in enumerate(rec["timepoints"]):
        for nid, ident in series["ids"][k].items():
            tracked[i][segm[i] == nid] = ident
    deleted_nid = next(nid for nid, ident in series["ids"][3].items() if ident == 10)
    tracked[3][segm[3] == deleted_nid] = 0
    np.savez_compressed(npz, tracked)                       # Cell-ACDC overwrites segm.npz in place
    res_ac = import_tracks(series["out"], "cellacdc", npz)
    info = res_ac["record"]["backend_info"]
    assert info["n_deleted_labels"] == 1 and info["n_merged_labels"] == 0 and info["n_split_labels"] == 0
    pts = res_ac["points"]
    assert pts["time_track_id"].isna().sum() == 1
    gone = pts[pts["time_track_id"].isna()].iloc[0]
    assert gone["t_index"] == 3 and series["ids"][3][int(gone["nucleus_id"])] == 10
    _check_identities(pts, series["ids"])
    assert (pts.loc[pts["time_track_id"].notna(), "link_source"] == "imported:cellacdc").all()
    # points measured on the exported stack come back through the export record
    exp_pts = _truth_points(series["ids"], identities=(3,))
    origin = np.array(rec["box"]["start"]) * np.array(SPACING)
    rows = []
    for r in exp_pts.itertuples(index=False):
        reg = registered_centre_um(3, int(r.t_index)) - origin
        rows.append({"id": 3, "t": rec["frame_of_t_index"][str(int(r.t_index))], "z": reg[0], "y": reg[1], "x": reg[2]})
    pd.DataFrame(rows).to_csv(tmp_path / "on_export.csv", index=False)
    res_sp = import_tracks(series["out"], "imaris", tmp_path / "on_export.csv", space="export", export_json=out / "export.json", write=False)
    got = res_sp["points"][res_sp["points"]["time_track_id"].notna()]
    assert len(got) == T and (got.pipe(_identity_of_points, series["ids"]) == 3).all()


# ---- track-strips ------------------------------------------------------------------------------------

def test_strips_are_written(series, tmp_path):
    link_series(series["out"], backend="laptrack")
    out = tmp_path / "strips"
    written = track_strips(series["out"], out, size_um=6.0, max_tracks=4)
    assert len(written) == 4 and all(p.exists() and p.stat().st_size > 1000 for p in written)
    assert written[0].name == f"{STEM}__track00001_strip.png"
    import matplotlib.image as mpimg

    img = mpimg.imread(written[0])
    assert img.ndim >= 2 and img.shape[1] > img.shape[0]          # a strip: T tiles side by side
    everything = track_strips(series["out"], tmp_path / "all_strips", planes=1)
    assert len(everything) == 18


# ---- the validation pick and track-validate ----------------------------------------------------------

def test_pick_validation_nuclei_is_stratified_and_deterministic(series, tmp_path):
    a = pick_validation_nuclei(series["out"], n=8, n_tz=3, n_pach=3, n_edge=2, edge_um=3.0, seed=1, write=False)
    b = pick_validation_nuclei(series["out"], n=8, n_tz=3, n_pach=3, n_edge=2, edge_um=3.0, seed=1, write=False)
    pd.testing.assert_frame_equal(a, b)
    c = pick_validation_nuclei(series["out"], n=8, n_tz=3, n_pach=3, n_edge=2, edge_um=3.0, seed=2, write=False)
    assert set(a["nucleus_id"]) != set(c["nucleus_id"])
    assert list(a["hand_track_id"]) == list(range(1, 9)) and a["role"].value_counts().to_dict() == {"holdout": 4, "constraint": 4}
    assert a["zone"].value_counts().to_dict() == {"transition_zone": 3, "pachytene": 3, "edge": 2}
    idents = {nid: ident for nid, ident in series["ids"][0].items()}
    assert {idents[n] for n in a.loc[a["zone"] == "edge", "nucleus_id"]} == {17, 18}      # the two 1.6 um from the border
    tz = a[a["zone"] == "transition_zone"]
    pach = a[a["zone"] == "pachytene"]
    assert tz["axis_position_um"].max() < pach["axis_position_um"].min()               # distal third before the middle band
    assert (a["t_index"] == 0).all() and (a["image_id"] == f"{STEM}_t000").all()
    assert not {idents[n] for n in a["nucleus_id"]} & set(JUNK_BASES)
    with pytest.raises(ValueError):
        pick_validation_nuclei(series["out"], n=8, n_tz=3, n_pach=3, n_edge=3, write=False)
    # the worksheet, through the script
    import importlib.util
    import sys

    spec = importlib.util.spec_from_file_location("pick_validation_nuclei", ROOT / "scripts" / "pick_validation_nuclei.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    out = tmp_path / "picks.csv"
    assert mod.main([str(series["out"]), "--n", "8", "--tz", "3", "--pach", "3", "--edge", "2", "--edge-um", "3.0", "--seed", "1",
                     "--out", str(out)]) == 0
    sheet = read_hand_points(out)
    assert len(sheet) == 8 and list(sheet["hand_track_id"]) == list(range(1, 9)) and set(sheet["role"]) == {"holdout", "constraint"}
    assert out.read_text().startswith("# ")
    # the default location is the series folder
    pick_validation_nuclei(series["out"], n=8, n_tz=3, n_pach=3, n_edge=2, edge_um=3.0, seed=1)
    assert (series["out"] / f"{STEM}__validation_picks.csv").exists()


def _hand_set(ids: dict, identities, swap: dict | None = None, miss: tuple | None = None) -> pd.DataFrame:
    """A hand-followed set: one row per timepoint per identity (raw-frame microns), hand_track_id = identity.
    ``swap`` {hand id: identity} from t=2 on exchanges what a hand track follows; ``miss`` (hand id, t) puts that
    point in the background."""
    rows = []
    for hid in identities:
        for k in range(T):
            ident = swap.get(hid, hid) if (swap and k >= 2) else hid
            if ident == DROPOUT[0] and k == DROPOUT[1]:
                continue
            z, y, x = raw_centre(ident, k) * np.array(SPACING)
            if miss and miss == (hid, k):
                y = 90 * SPACING[1]                       # empty background
            rows.append({"hand_track_id": hid, "t_index": k, "z_um": z, "y_um": y, "x_um": x,
                         "role": "holdout" if hid % 2 else "constraint", "zone": "pachytene" if hid <= 2 else "transition_zone"})
    return pd.DataFrame(rows)


def test_track_validate_scores_a_perfect_and_a_swapped_hand_set(series, tmp_path):
    link_series(series["out"], backend="laptrack")
    perfect = tmp_path / "perfect.csv"
    _hand_set(series["ids"], (1, 2, 3, 4, DROPOUT[0])).to_csv(perfect, index=False)
    res = validate_tracks(series["out"], perfect)
    tr, v = res["tracks"], res["verdict"]
    assert len(tr) == 5 and (tr["segmentation_recall"] == 1.0).all() and (tr["linking_recall"] == 1.0).all()
    assert (tr["edge_precision"] == 1.0).all() and (tr["switches"] == 0).all() and (tr["edges_missing"] == 0).all()
    assert tr.set_index("hand_track_id").loc[DROPOUT[0], "n_points"] == T - 1 and tr.set_index("hand_track_id").loc[DROPOUT[0], "n_edges"] == T - 2
    assert v["pass"] and v["criteria"]["switches_holdout"]["value"] == 0 and v["criteria"]["linking_recall"]["value"] == 1.0
    assert v["criteria"]["segmentation_recall"]["set"] == "pachytene" and v["criteria"]["segmentation_recall"]["value"] == 1.0
    assert set(res["points"]["outcome"]) == {"linked"}
    table = pd.read_csv(series["out"] / f"{STEM}__track_validation.csv")
    assert list(table.columns[:len(schema.TRACK_VALIDATION)]) == list(schema.TRACK_VALIDATION)
    assert len(table) == 5 * T - 1 and set(table["role"]) == {"holdout", "constraint"} and (table["n_points"] >= T - 1).all()
    assert table["time_track_id"].notna().all() and table["nucleus_id"].notna().all()
    assert json.loads((series["out"] / f"{STEM}__track_validation.json").read_text())["verdict"]["pass"] is True
    assert "VERDICT: PASS" in format_verdict(res)
    # hand tracks 1 and 2 exchange nuclei from t=2 on, and one point of hand track 3 lies in the background
    swapped = tmp_path / "swapped.csv"
    _hand_set(series["ids"], (1, 2, 3, 4), swap={1: 2, 2: 1}, miss=(3, 2)).to_csv(swapped, index=False)
    res2 = validate_tracks(series["out"], swapped)
    tr2 = res2["tracks"].set_index("hand_track_id")
    assert tr2.loc[1, "switches"] == 1 and tr2.loc[2, "switches"] == 1 and tr2.loc[1, "linking_recall"] == pytest.approx(3 / 4)
    assert tr2.loc[1, "edge_precision"] == pytest.approx(3 / 4)
    assert tr2.loc[3, "segmentation_recall"] == pytest.approx(4 / 5) and tr2.loc[3, "n_in_label"] == 4
    assert tr2.loc[3, "linking_recall"] == 1.0 and tr2.loc[3, "switches"] == 0          # the miss is not a linking error
    assert tr2.loc[4, "linking_recall"] == 1.0 and tr2.loc[4, "switches"] == 0
    v2 = res2["verdict"]
    assert not v2["pass"] and v2["criteria"]["switches_holdout"]["value"] == 1          # hand track 1 is holdout
    assert not v2["criteria"]["linking_recall"]["pass"]
    outcomes = res2["points"].set_index(["hand_track_id", "t_index"])["outcome"]
    assert outcomes.loc[(1, 2)] == "switch" and outcomes.loc[(3, 2)] == "segmentation_miss" and outcomes.loc[(1, 3)] == "linked"
    assert "VERDICT: FAIL" in format_verdict(res2)
    # a fragmented track (overlap cannot bridge the dropout) is a missing edge, not a switch
    link_series(series["out"], backend="overlap_inhouse", stride=1)
    res3 = validate_tracks(series["out"], perfect, write=False)
    tr3 = res3["tracks"].set_index("hand_track_id")
    assert tr3.loc[DROPOUT[0], "edges_missing"] == 1 and tr3.loc[DROPOUT[0], "switches"] == 0
    assert tr3.loc[DROPOUT[0], "linking_recall"] == pytest.approx(2 / 3)
    assert res3["points"].set_index(["hand_track_id", "t_index"])["outcome"].loc[(DROPOUT[0], 3)] == "edge_missing"


# ---- the command line --------------------------------------------------------------------------------

def test_cli_track_commands(series, tmp_path, capsys):
    sd = str(series["out"])
    assert cli.main(["track-export", sd, "--out", str(tmp_path / "exp"), "--stride", "2"]) == 0
    assert (tmp_path / "exp" / "export.json").exists()
    assert cli.main(["track-link", sd, "--backend", "none"]) == 0
    assert cli.main(["track-link", sd, "--backend", "laptrack", "--gap", "2", "--max-dist-um", "3.0"]) == 0
    out = capsys.readouterr().out
    assert "laptrack" in out and "18 track(s)" in out
    assert cli.main(["track-strips", sd, "--out", str(tmp_path / "strips"), "--max-tracks", "2"]) == 0
    assert len(list((tmp_path / "strips").glob("*.png"))) == 2
    hand = tmp_path / "hand.csv"
    _hand_set(series["ids"], (1, 2)).to_csv(hand, index=False)
    assert cli.main(["track-validate", sd, "--points", str(hand)]) == 0
    assert "VERDICT: PASS" in capsys.readouterr().out
    _hand_set(series["ids"], (1, 2), swap={1: 2, 2: 1}).to_csv(hand, index=False)
    assert cli.main(["track-validate", sd, "--points", str(hand)]) == 2
    tm = _truth_points(series["ids"])
    csv = tmp_path / "tm.csv"
    pd.DataFrame({"ID": range(len(tm)), "TRACK_ID": tm["track_id"], "POSITION_X": tm["x_um"], "POSITION_Y": tm["y_um"],
                  "POSITION_Z": tm["z_um"], "FRAME": tm["t_index"]}).to_csv(csv, index=False)
    assert cli.main(["track-import", sd, "--from", "trackmate", str(csv)]) == 0
    assert "imported:trackmate" in pd.read_csv(series["out"] / f"{STEM}__track_points.csv")["link_source"].unique()


def test_tracking_schema_is_append_only():
    base_points = ["time_track_id", "t_index", "t_seconds", "image_id", "nucleus_id", "z_um", "y_um", "x_um",
                   "reg_z_um", "reg_y_um", "reg_x_um", "step_um", "step_z_um", "step_y_um", "step_x_um", "is_gap", "link_source"]
    base_tracks = ["time_track_id", "n_points", "n_gaps", "t_start", "t_end", "is_full_length", "path_length_um",
                   "net_displacement_um", "speed_um_per_min", "axis_velocity_um_per_min", "intensity_slope_raw_per_min",
                   "intensity_slope_norm_per_min", "sc_length_change_um", "link_source"]
    base_val = ["hand_track_id", "validation_set", "zone", "t_index", "t_seconds", "image_id", "z_um", "y_um", "x_um",
                "nucleus_id", "time_track_id", "outcome"]
    assert schema.TRACK_POINTS[:len(base_points)] == base_points and schema.TRACK_POINTS[-3:] == ["zr_um", "yr_um", "xr_um"]
    assert schema.TRACKS[:len(base_tracks)] == base_tracks and "t_first" in schema.TRACKS and "mean_step_um" in schema.TRACKS
    assert schema.TRACK_VALIDATION[:len(base_val)] == base_val and schema.TRACK_VALIDATION[-1] == "role"
    for cols in (schema.TRACK_POINTS, schema.TRACKS, schema.TRACK_VALIDATION):
        assert len(set(cols)) == len(cols)


# ---- a folder process_series wrote --------------------------------------------------------------------

def test_tracking_reads_a_folder_process_series_wrote(tmp_path, patched_series_reader):  # noqa: F811
    """The real file set (no reference volume, no shift table): the points table is built from it, the overlap
    backend links it, the pick works, and the export says it wrote no image."""
    out = tmp_path / "series"
    process_series(ND2, _cfg(), out, timepoints="all")
    res = link_series(out, backend="none")
    pts = res["points"]
    assert len(pts) and set(pts["t_index"]) == {0, 1, 2} and (pts["image_id"].str.startswith(GONAD_STEM + "_t")).all()
    assert pts["in_germline"].all() and res["record"]["shift_source"] == "none"
    res2 = link_series(out, backend="overlap_inhouse", stride=1)
    assert res2["record"]["n_tracks"] >= 1 and (out / f"{GONAD_STEM}__tracks.csv").exists()
    picks = pick_validation_nuclei(out, n=4, n_tz=2, n_pach=1, n_edge=1, write=False)
    assert len(picks) == 4
    rec = export_cellacdc(out, tmp_path / "exp", stride=2, image="ref")
    assert rec["image"]["source"] == "none" and rec["image"]["reason"] == "no reference volume on disk"
    with pytest.raises(Exception):
        track_strips(out, tmp_path / "strips")                  # no reference volume and no readable series file
