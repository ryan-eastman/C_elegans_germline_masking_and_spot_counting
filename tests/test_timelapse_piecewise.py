"""Piecewise registration along the gonad axis (germquant.timelapse.piecewise, `germquant register --mode
piecewise`, `germquant track-link --registration piecewise`; docs/TIMELAPSE.md "Piecewise registration").

The synthetic series mirrors what the first real series showed: the worm (the brightfield) does not move, and
the germline slides inside it non-rigidly, the distal end still and the proximal end sliding. Sixty nuclei in
two jittered rows along x (a straight gonad, axis_position_um = x - 10 um at t0) move by a shear flow: each
nucleus moves u(s) = U s / L per timepoint along x, with s its axis position at t0 (0 at the distal end, L at
the proximal end) and U = 4.6 um, more than half the 4 um spacing along a row. One translation per frame leaves
the two ends about 2.3 um off in opposite directions, past half the spacing, so the rigid registration sends the
end nuclei to their neighbours; the piecewise one follows the shear. The series folder is written the way
process_series lays one out (reference volumes of the central_element and brightfield roles, label images,
nuclei tables, image summaries, the series record), so these tests do not depend on the segmenter.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import tifffile

from germquant import cli, schema
from germquant.timelapse import (
    ReferenceVolumesMissing,
    find_series,
    link_series,
    load_piecewise,
    register_series,
    register_series_piecewise,
)
from germquant.timelapse import piecewise as PW
from germquant.timelapse import track as TR

STEM = "20261005_DLW93_nohs_gfpsyp2_5min_2hr_herm"
SP = (0.4, 0.2, 0.2)                     # (dz, dy, dx) um of the label grid
SHAPE = (16, 120, 800)                   # (Z, Y, X) voxels: 6.4 x 24 x 160 um
REF_STRIDE = 4                           # reference volumes: every 4th pixel, 0.8 um voxels in y and x
T = 5
PERIOD_S = 300.0
U_UM = 4.6                               # proximal flow per timepoint (um, +x); the distal end is still
X0_UM = 10.0                             # the distal end of the axis
R_UM = 1.4                               # label sphere radius


def _bases(seed=7):
    """identity -> t0 centre (z, y, x) um of the 60 nuclei: two rows (y 8 and 16 um), 4 um apart along x with a
    little jitter so the pattern is not a lattice, and their relative brightness."""
    rng = np.random.default_rng(seed)
    out, amp = {}, {}
    ident = 1
    for y in (8.0, 16.0):
        for x in np.arange(X0_UM, X0_UM + 30 * 4.0, 4.0):
            out[ident] = np.array([3.2 + rng.uniform(-0.3, 0.3), y + rng.uniform(-1.0, 1.0), x + rng.uniform(-0.45, 0.45)])
            amp[ident] = rng.uniform(0.6, 1.4)
            ident += 1
    return out, amp


BASES, AMP = _bases()
L_UM = max(c[2] for c in BASES.values()) - X0_UM


def flow_per_step(identity: int) -> float:
    """u(s): the x step per timepoint of an identity, from its t0 axis position s (shear: 0 distal, U proximal)."""
    return U_UM * (BASES[identity][2] - X0_UM) / L_UM


def centre_um(identity: int, k: int) -> np.ndarray:
    return BASES[identity] + np.array([0.0, 0.0, k * flow_per_step(identity)])


def _render_ref(centres, amps, shape, spacing, sigma_um=(0.6, 1.0, 1.0)):
    """Gaussian nuclei sampled on a reference grid (voxel j at j * spacing um), plus a flat background."""
    v = np.full(shape, 100.0, np.float32)
    zz, yy, xx = (np.arange(n) * s for n, s in zip(shape, spacing))
    for c, a in zip(centres, amps):
        gz = np.exp(-0.5 * ((zz - c[0]) / sigma_um[0]) ** 2)
        gy = np.exp(-0.5 * ((yy - c[1]) / sigma_um[1]) ** 2)
        gx = np.exp(-0.5 * ((xx - c[2]) / sigma_um[2]) ** 2)
        v += np.float32(2000.0 * a) * (gz[:, None, None] * gy[None, :, None] * gx[None, None, :]).astype(np.float32)
    return v


def _brightfield(shape, spacing, seed=11):
    """The worm body in transmitted light: a bright field with dark blobs that do not move between frames."""
    rng = np.random.default_rng(seed)
    blobs = np.column_stack([rng.uniform(1.0, 5.4, 160), rng.uniform(1.0, 23.0, 160), rng.uniform(1.0, 159.0, 160)])
    return 3000.0 - _render_ref(blobs, rng.uniform(0.2, 0.5, 160), shape, spacing, (0.8, 1.2, 1.2)) + 100.0


def write_piecewise_series(root: Path, *, n_t: int = T) -> dict:
    """Write the series folder; returns the folder and ``ids[k] = {nucleus_id: identity}`` per timepoint."""
    out = root / STEM
    zz, yy, xx = np.indices(SHAPE).astype(np.float32)
    ref_shape = (SHAPE[0], -(-SHAPE[1] // REF_STRIDE), -(-SHAPE[2] // REF_STRIDE))
    ref_sp = (SP[0], SP[1] * REF_STRIDE, SP[2] * REF_STRIDE)
    bf = np.clip(np.rint(_brightfield(ref_shape, ref_sp)), 0, 65535).astype(np.uint16)
    rng = np.random.default_rng(2026)
    ids: dict[int, dict[int, int]] = {}
    entries = []
    for k in range(n_t):
        idents = sorted(BASES)
        order = rng.permutation(len(idents))
        nid_of = {idents[j]: int(order[j]) + 1 for j in range(len(idents))}     # nucleus_id never equals the identity
        lab = np.zeros(SHAPE, np.int32)
        cents = {i: centre_um(i, k) for i in idents}
        for i in idents:
            c = cents[i]
            d2 = (zz * SP[0] - c[0]) ** 2 + (yy * SP[1] - c[1]) ** 2 + (xx * SP[2] - c[2]) ** 2
            lab[d2 <= R_UM ** 2] = nid_of[i]
        sub = out / f"t{k:03d}"
        sub.mkdir(parents=True)
        iid = f"{STEM}_t{k:03d}"
        tifffile.imwrite(sub / f"{iid}__nuclei_labels.tif", lab, compression="zlib")
        ref = _render_ref([cents[i] for i in idents], [AMP[i] for i in idents], ref_shape, ref_sp)
        np.save(sub / f"{iid}__ref_stride{REF_STRIDE}.npy", np.clip(np.rint(ref), 0, 65535).astype(np.uint16))
        np.save(sub / f"{iid}__ref_brightfield_stride{REF_STRIDE}.npy", bf)
        rows = []
        for i in idents:
            c = cents[i]
            rows.append({"nucleus_id": nid_of[i], "volume_um3": 4.0 / 3.0 * np.pi * R_UM ** 3, "n_voxels": int((lab == nid_of[i]).sum()),
                         "centroid_z_um": c[0], "centroid_y_um": c[1], "centroid_x_um": c[2], "in_germline": True,
                         "axis_position_um": c[2] - X0_UM, "axis_position_norm": (c[2] - X0_UM) / 150.0, "n_spots": 0,
                         "central_element_mean_intensity": 1000.0 * AMP[i]})
        nuclei = pd.DataFrame(rows).sort_values("nucleus_id")
        for col in schema.NUCLEI:
            if col not in nuclei.columns:
                nuclei[col] = pd.NA
        nuclei = nuclei[list(schema.NUCLEI) + [c for c in nuclei.columns if c not in schema.NUCLEI]]
        meta = {"image_id": iid, "file_path": "x.nd2", "genotype": "DLW93", "sex": "herm", "germ_cell": "oocyte",
                "treatment": "control", "replicate": "NA", "acquisition_date": "20261005",
                "voxel_dz_um": SP[0], "voxel_dy_um": SP[1], "voxel_dx_um": SP[2], "n_channels": 2, "dapi_present": False,
                "channel_map_json": json.dumps({"central_element": 1, "brightfield": 0}), "pipeline_version": "0.1.0",
                "git_sha": "test", "config_hash": "cafe", "run_timestamp": "2026-10-10T00:00:00+00:00"}
        for c, v in meta.items():
            nuclei[c] = v
        nuclei.to_csv(sub / f"{iid}__nuclei.csv", index=False)
        pd.DataFrame([{"n_nuclei": len(nuclei), "n_germline_nuclei": len(nuclei), "mean_spots": 0.0,
                       "total_germline_length_um": 150.0, "qc_pass": True, "qc_flags": "", **meta}]).to_csv(
            sub / f"{iid}__image_summary.csv", index=False)
        (sub / f"{iid}__done.json").write_text(json.dumps({"image_id": iid, "config_hash": "cafe", "t_index": k,
                                                           "t_seconds": k * PERIOD_S, "failed_stages": []}))
        ids[k] = {v: i for i, v in nid_of.items()}
        entries.append({"t_index": k, "t_seconds": k * PERIOD_S, "image_id": iid, "out_dir": sub.name, "status": "done",
                        "reference_volumes": {"stride": REF_STRIDE, "role": "central_element", "channel_index": 1,
                                              "role_resolved": True, "spacing_um": list(ref_sp), "shape": list(ref_shape),
                                              "dtype": "uint16",
                                              "files": {"central_element": f"t{k:03d}/{iid}__ref_stride{REF_STRIDE}.npy",
                                                        "brightfield": f"t{k:03d}/{iid}__ref_brightfield_stride{REF_STRIDE}.npy"}}})
    record = {"series_id": STEM, "file_path": "x.nd2", "config_hash": "cafe", "pipeline_version": "0.1.0", "git_sha": "test",
              "n_timepoints": n_t, "period_s": PERIOD_S, "t_seconds": [k * PERIOD_S for k in range(n_t)], "xy_stride": 1,
              "z_range": None, "requested": list(range(n_t)), "done": list(range(n_t)), "failed": [], "timepoints": entries,
              "complete": True}
    (out / f"{STEM}__series.json").write_text(json.dumps(record, indent=1))
    (out / "run_manifest.json").write_text(json.dumps({"config": {"timelapse": {"register": {"piecewise_segments": 4}}}}))
    return {"out": out, "ids": ids}


@pytest.fixture(scope="module")
def series(tmp_path_factory):
    s = write_piecewise_series(tmp_path_factory.mktemp("pw"))
    s["piecewise"] = register_series_piecewise(s["out"], max_step_um=6.0)
    return s


def _identity_errors(points: pd.DataFrame, ids: dict) -> tuple[int, int]:
    """(tracks holding more than one identity, consecutive links that join two different identities)."""
    linked = points[points["time_track_id"].notna()].copy()
    linked["identity"] = [ids[int(t)][int(n)] for t, n in zip(linked["t_index"], linked["nucleus_id"])]
    mixed = int((linked.groupby("time_track_id")["identity"].nunique() > 1).sum())
    wrong = 0
    for _tid, g in linked.sort_values("t_index").groupby("time_track_id"):
        idents = g["identity"].to_numpy()
        wrong += int(np.sum(idents[1:] != idents[:-1]))
    return mixed, wrong


# ---- the pieces ------------------------------------------------------------------------------------------------

def test_segments_interpolation_and_settings():
    a0 = np.arange(40, dtype=float)
    edges, centres, seg = PW.segment_bounds(a0, 4)
    assert np.allclose(edges, [0, 9.75, 19.5, 29.25, 39]) and list(np.bincount(seg)) == [10, 10, 10, 10]
    assert np.allclose(centres, [4.5, 14.5, 24.5, 34.5])
    assert list(PW.assign_segments([-5.0, 9.8, 39.0, 100.0, np.nan], edges)) == [0, 1, 3, 3, -1]
    with pytest.raises(ValueError, match="too few"):
        PW.segment_bounds(np.arange(12.0), 4)
    cum = np.array([[0, 0, 0], [0, 0, 1], [0, 0, 3], [0, 1, 3]], float)
    d = PW.interp_displacement(cum, centres, [0.0, 4.5, 9.5, 19.5, 34.5, 50.0, np.nan])
    assert np.allclose(d[:6, 2], [0.0, 0.0, 0.5, 2.0, 3.0, 3.0])            # linear between centres, constant beyond
    assert np.allclose(d[:6, 1], [0, 0, 0, 0, 1, 1]) and np.isnan(d[6]).all()
    assert np.allclose(PW.interp_displacement([[1, 2, 3]], [5.0], [0.0, 9.0]), [[1, 2, 3], [1, 2, 3]])
    # the registered position solves reg = raw - D(a(reg)) on a straight axis (a = x)
    axis = lambda p: np.asarray(p, float)[:, 2]  # noqa: E731
    raw = np.array([[0.0, 0.0, 12.0], [0.0, 0.0, 30.0], [np.nan, 0.0, 1.0]])
    reg, a = PW.register_points(raw, cum, centres, axis)
    assert np.allclose(reg[:2, 2] + PW.interp_displacement(cum, centres, reg[:2, 2])[:, 2], raw[:2, 2], atol=1e-6)
    assert np.allclose(a[:2], reg[:2, 2]) and np.isnan(reg[2]).all() and np.isnan(a[2])
    assert PW.piecewise_settings(None) == {"segments": 4, "margin_um": 6.0, "taper": 0.25}
    cfg = {"timelapse": {"register": {"piecewise_segments": 3, "piecewise_margin_um": 2.5}}}
    assert PW.piecewise_settings(cfg) == {"segments": 3, "margin_um": 2.5, "taper": 0.25}
    assert PW.piecewise_settings(cfg, segments=6)["segments"] == 6                 # an explicit value wins over the key
    for bad in ({"segments": 0}, {"margin_um": -1.0}, {"taper": 1.5}):
        with pytest.raises(ValueError):
            PW.piecewise_settings(None, **bad)


def test_place_frame_labels_matches_place_frame_for_one_shift():
    rng = np.random.default_rng(0)
    lab = np.zeros((6, 30, 40), np.int32)
    for v in range(1, 9):
        z, y, x = rng.integers(1, 5), rng.integers(3, 27), rng.integers(3, 37)
        lab[z - 1:z + 1, y - 2:y + 2, x - 2:x + 2] = v
    box = {"start": [0, 2, 3], "stop": [6, 28, 38]}
    sv = np.array([1, -2, 3])
    for stride in (1, 2):
        a = TR.place_frame(lab, sv, box, stride)
        b = TR.place_frame_labels(lab, {v: sv for v in range(1, 9)}, box, stride)
        assert np.array_equal(a, b)
    # two labels with their own shifts: each lands where its own shift puts it
    c = TR.place_frame_labels(lab, {1: np.array([0, 0, 0]), 2: np.array([0, 0, 5])}, {"start": [0, 0, 0], "stop": [6, 30, 40]}, 1)
    z, y, x = np.nonzero(lab == 2)
    m = x - 5 >= 0
    assert (c[z[m], y[m], x[m] - 5] == 2).all()
    assert (c[lab == 1] == 1).sum() >= 1


# ---- the registration on the planted shear -----------------------------------------------------------------------

def test_piecewise_registration_recovers_the_segment_shifts(series):
    res = series["piecewise"]
    csv = pd.read_csv(res["table_path"])
    assert Path(res["table_path"]).name == f"{STEM}__register_piecewise.csv"
    assert Path(res["table_path"]).with_suffix(".parquet").exists()
    assert list(csv.columns) == schema.REGISTER_PIECEWISE + schema.SERIES_META
    assert len(csv) == T * 4 and list(csv["segment"].unique()) == [0, 1, 2, 3] and (csv["n_segments"] == 4).all()
    assert (csv["origin_t_index"] == 0).all() and csv["error"].fillna("").eq("").all() and not csv["filled"].any()
    assert (csv["channel_role"] == "central_element").all() and (csv["xy_stride"] == REF_STRIDE).all()
    assert (csv["taper"] == 0.25).all() and (csv["margin_um"] == 6.0).all() and (csv["series_id"] == STEM).all()
    seg = csv[csv["t_index"] == 0].sort_values("segment")
    assert (seg[["dx_um", "dy_um", "dz_um", "step_dx_um"]] == 0).all().all()
    assert list(seg["n_nuclei_t0"]) == [15, 15, 15, 15]
    centres = seg["axis_centre_um"].to_numpy(float)
    assert np.all(np.diff(centres) > 20) and centres[0] < 20 and centres[-1] > 95
    # each segment's step is the shear at its centre: the distal segment barely moves, the proximal one slides
    for r in csv[csv["t_index"] > 0].itertuples(index=False):
        planted = U_UM * float(r.axis_centre_um) / L_UM
        assert abs(r.step_dx_um - planted) < 0.3, (r.t_index, r.segment, r.step_dx_um, planted)
        assert abs(r.step_dy_um) < 0.2 and abs(r.step_dz_um) < 0.1
        assert abs(r.dx_um - r.t_index * planted) < 0.3 * r.t_index + 1e-9
        assert int(r.reference_t_index) == int(r.t_index) - 1
    # the worm (brightfield) does not move: the cross-check reads zero and disagrees where the germline slides
    later = csv[csv["t_index"] > 0]
    assert later[["brightfield_step_dx_um", "brightfield_step_dy_um", "brightfield_step_dz_um"]].abs().max().max() < 0.15
    assert list(later.loc[later["disagreement_flag"], "segment"].unique()) == [1, 2, 3]
    assert not later.loc[later["segment"] == 0, "disagreement_flag"].any()
    # quality: per segment, the registered nuclei sit on the previous frame's; raw, they sit up to U away (the raw
    # number saturates near half the spacing, because a nucleus that moved most of a spacing lands next to its
    # neighbour's old place)
    assert (later["quality_after_um"] < 0.35).all()
    assert (later["quality_after_um"] < 0.5 * later["quality_before_um"]).all()
    assert later["n_nuclei"].between(13, 17).all() and later.groupby("t_index")["n_nuclei"].sum().eq(60).all()
    rep = json.loads(Path(res["report_path"]).read_text(encoding="utf-8"))
    assert Path(res["report_path"]).name == f"{STEM}__register_piecewise_report.json"
    par = rep["parameters"]
    assert par["segments"] == 4 and par["reference"] == "previous" and par["origin_t_index"] == 0
    assert par["axis_source"] == "axis_position_um" and par["cross_check_role"] == "brightfield" and par["max_step_um"] == 6.0
    assert rep["mode"] == "piecewise" and rep["n_registered"] == T - 1 and len(rep["steps"]) == T - 1
    assert [s["n_nuclei_t0"] for s in rep["segments"]] == [15, 15, 15, 15] and rep["errors"] == {} and rep["filled"] == []
    q = rep["quality"]["by_segment"]
    assert all(x["median_after_um"] < x["median_before_um"] for x in q[1:])
    assert rep["cumulative_last_um"][3][2] == pytest.approx((T - 1) * U_UM * centres[3] / L_UM, abs=0.3 * (T - 1))
    # the rigid table was not touched (none existed before)
    assert not (series["out"] / f"{STEM}__register_shifts.csv").exists()


def test_interpolated_displacements_follow_the_planted_flow(series):
    pw = load_piecewise(find_series(series["out"]))
    assert pw.n_segments == 4 and pw.origin_t == 0 and sorted(pw.cum) == list(range(T))
    lo, hi = pw.centres[0], pw.centres[-1]
    for k in range(T):
        idents = sorted(BASES)
        raw = np.array([centre_um(i, k) for i in idents])
        d = pw.displacement(k, raw)
        planted = np.array([k * flow_per_step(i) for i in idents])
        s0 = np.array([BASES[i][2] - X0_UM for i in idents])
        inner = (s0 >= lo) & (s0 <= hi)
        assert np.abs(d[inner, 2] - planted[inner]).max() < 0.3 * k + 1e-9, k            # linear field: exact between centres
        # beyond the outermost centres the displacement is held constant: off by at most the shear over that end
        assert np.abs(d[~inner, 2] - planted[~inner]).max() < k * (U_UM * max(lo, L_UM - hi) / L_UM + 0.3) + 1e-9
        assert np.abs(d[:, :2]).max() < 0.2 * k + 1e-9
        reg, a = pw.register(k, raw)
        assert np.allclose(raw - reg, d)
        assert np.abs(a[inner] - s0[inner]).max() < 0.3 * k + 1e-6                      # the registered axis position is the t0 one
    # a timepoint without an estimate takes the last one that has one
    assert np.allclose(pw.cumulative(99), pw.cum[T - 1]) and pw.summary()["timepoints_substituted"] == {"99": T - 1}


# ---- linking: rigid breaks, piecewise recovers ------------------------------------------------------------------

def test_linking_that_rigid_registration_breaks(series):
    pytest.importorskip("laptrack")    # the timelapse extra; CI installs it
    out = series["out"]
    register_series(out, max_step_um=6.0)                        # one translation per frame
    rigid = link_series(out, backend="laptrack", gap=2, max_dist_um=3.0, write=False)
    assert rigid["record"]["registration"] == "rigid" and rigid["record"]["shift_source"] == "register_shifts"
    mixed_r, wrong_r = _identity_errors(rigid["points"], series["ids"])
    assert wrong_r >= 8, (mixed_r, wrong_r)                       # the end nuclei jump to their neighbours
    pw = link_series(out, backend="laptrack", gap=2, max_dist_um=3.0, registration="piecewise")
    rec = pw["record"]
    assert rec["registration"] == "piecewise" and rec["shift_source"] == "register_piecewise"
    assert rec["registration_info"]["n_segments"] == 4 and rec["registration_info"]["origin_t_index"] == 0
    assert _identity_errors(pw["points"], series["ids"]) == (0, 0)
    assert rec["n_tracks"] == 60 and rec["n_full_length"] == 60
    pts = pw["points"]
    assert (pts["link_source"] == "laptrack").all()
    # the registered steps are what is left after the flow (small); the raw-frame axis velocity is the flow
    steps = pts["step_um"].dropna()
    assert steps.max() < 1.6 and steps.median() < 0.3
    assert {"disp_x_um", "disp_y_um", "disp_z_um", "axis_position_fixed_raw_um"} <= set(pts.columns)
    np.testing.assert_allclose(pts["x_um"] - pts["xr_um"], pts["disp_x_um"], atol=1e-9)
    tracks = pw["tracks"].set_index("time_track_id")
    first = pts[pts["t_index"] == 0].set_index("time_track_id")
    ident = {tid: series["ids"][0][int(n)] for tid, n in first["nucleus_id"].items()}
    v_flow = np.array([flow_per_step(ident[t]) / (PERIOD_S / 60.0) for t in tracks.index])
    np.testing.assert_allclose(tracks["axis_velocity_raw_um_per_min"].to_numpy(float), v_flow, atol=0.01)
    v_reg = np.abs(tracks["axis_velocity_um_per_min"].to_numpy(float))                     # relative to the flow: small
    assert v_reg.max() < 0.2 and np.median(v_reg) < 0.05 and v_flow.max() > 0.9
    tt = pd.read_csv(out / f"{STEM}__tracks.csv")
    assert list(tt.columns[:len(schema.TRACKS)]) == list(schema.TRACKS)
    assert json.loads((out / f"{STEM}__track_link.json").read_text())["registration"] == "piecewise"
    # the overlap rule on per-nucleus shifted labels follows the shear too
    ov = link_series(out, backend="overlap_inhouse", stride=1, registration="piecewise", write=False)
    assert _identity_errors(ov["points"], series["ids"]) == (0, 0)
    assert ov["record"]["n_full_length"] == 60 and "per-nucleus" in ov["record"]["backend_info"]["placement"]
    ov_r = link_series(out, backend="overlap_inhouse", stride=1, write=False)
    assert ov_r["record"]["n_full_length"] < 60


def test_cli_register_piecewise_and_track_link(tmp_path, capsys):
    s = write_piecewise_series(tmp_path, n_t=3)
    sd = str(s["out"])
    # without the table, track-link --registration piecewise says what to run and exits 1
    assert cli.main(["track-link", sd, "--backend", "none", "--registration", "piecewise"]) == 1
    assert "--mode piecewise" in capsys.readouterr().err
    assert cli.main(["register", sd, "--mode", "piecewise", "--reference", "t0"]) == 1
    assert cli.main(["register", sd, "--mode", "piecewise", "--rigid"]) == 1
    assert cli.main(["register", sd, "--segments", "3"]) == 1                 # --segments needs --mode piecewise
    capsys.readouterr()
    # the run_manifest.json config supplies the segment count; the flag wins over it
    (s["out"] / "run_manifest.json").write_text(json.dumps({"config": {"timelapse": {"register": {"piecewise_segments": 2}}}}))
    assert cli.main(["register", sd, "--mode", "piecewise", "--max-step-um", "6"]) == 0
    out = capsys.readouterr().out
    assert "piecewise registration, 2 segment(s)" in out and "segment 1" in out and "segment 2" not in out
    assert cli.main(["register", sd, "--mode", "piecewise", "--segments", "3", "--margin-um", "4", "--max-step-um", "6"]) == 0
    csv = pd.read_csv(s["out"] / f"{STEM}__register_piecewise.csv")
    assert (csv["n_segments"] == 3).all() and (csv["margin_um"] == 4.0).all() and len(csv) == 3 * 3
    assert not (s["out"] / f"{STEM}__register_shifts.csv").exists()
    pytest.importorskip("laptrack")
    assert cli.main(["track-link", sd, "--backend", "laptrack", "--registration", "piecewise"]) == 0
    assert "piecewise registration, shifts from register_piecewise" in capsys.readouterr().out
    rec = json.loads((s["out"] / f"{STEM}__track_link.json").read_text())
    assert rec["registration"] == "piecewise" and rec["registration_info"]["n_segments"] == 3
    # the default stays rigid (no shift table here: the series is taken as aligned)
    assert cli.main(["track-link", sd, "--backend", "laptrack"]) == 0
    rec = json.loads((s["out"] / f"{STEM}__track_link.json").read_text())
    assert rec["registration"] == "rigid" and rec["shift_source"] == "none" and rec["registration_info"] is None


def test_piecewise_needs_reference_volumes_and_germline_nuclei(tmp_path):
    s = write_piecewise_series(tmp_path / "a", n_t=2)
    for p in s["out"].rglob("*.npy"):
        p.unlink()
    with pytest.raises(ReferenceVolumesMissing, match="timelapse.keep_reference_volume"):
        register_series_piecewise(s["out"])
    assert cli.main(["register", str(s["out"]), "--mode", "piecewise"]) == 1
    s2 = write_piecewise_series(tmp_path / "b", n_t=2)
    f = s2["out"] / "t000" / f"{STEM}_t000__nuclei.csv"
    df = pd.read_csv(f)
    df.iloc[:4].to_csv(f, index=False)                                      # four nuclei at the origin: too few
    with pytest.raises(ValueError, match="fewer than 5|too few"):
        register_series_piecewise(s2["out"])
    with pytest.raises(PW.PiecewiseRegistrationMissing, match="--mode piecewise"):
        load_piecewise(find_series(s2["out"]))


def test_series_tables_append_register_piecewise():
    assert list(schema.SERIES_TABLES)[-1] == "register_piecewise"
    assert list(schema.SERIES_TABLES)[:5] == ["register_shifts", "track_points", "tracks", "track_validation", "series_summary"]
    cols = schema.SERIES_TABLES["register_piecewise"]
    assert len(set(cols)) == len(cols) and "track_id" not in cols
    assert schema.TRACKS[-1] == "axis_velocity_raw_um_per_min" and schema.TRACKS[-2] == "axis_velocity_refit_um_per_min"
