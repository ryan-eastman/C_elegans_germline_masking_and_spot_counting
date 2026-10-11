"""Drift registration as a shift table (germquant.timelapse.register, `germquant register`;
docs/TIMELAPSE.md "Registration") and the reference volumes the series loop keeps for it.

The registration is exercised on a synthetic series of moving blobs with a planted sub-voxel drift per step
(0.4 planes, 1.3 and 2.7 voxels: never a whole voxel, so the upsampled estimate is what is tested), a
division-free shear flow of the nuclei (zero mean, so the drift stays identifiable while every nucleus moves
relative to its neighbours) and a brightfield channel that drifts with the sample. The folder a series run
leaves (reference volumes, nuclei and image_summary tables, the series record) is written directly, so these
tests do not depend on the segmenter. The loop itself (reference volumes written next to the fixed-image file
set, read back on --resume, `germquant register` on the result) is checked on the three-timepoint gonad of
test_timelapse_series. Tolerances are in reference-volume voxels: phase correlation works on the strided
volume and its upsampled estimate is quantised to a tenth of one.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_timelapse_series import DRIFT_PX, ND2, STEM, T, _series_stacks, patched_series_reader  # noqa: F401  (fixture)
from test_timelapse_series import SPACING as GONAD_SPACING

from germquant import api, cli, pipeline, schema
from germquant.config import load_config
from germquant.io.nd2_reader import Stack
from germquant.timelapse import process_series, read_series_record
from germquant.timelapse import register as R
from germquant.timelapse.series import BRIGHTFIELD_ROLE, reference_settings, reference_volume_path, write_reference_volumes

ROOT = Path(__file__).resolve().parents[1]
SP = (0.3, 0.108, 0.108)                  # (dz, dy, dx) um of the synthetic series
SHAPE = (24, 256, 256)                    # (Z, Y, X); deep enough that no blob touches a z edge
STEP_VOX = np.array([0.4, 1.3, 2.7])      # planted drift per step, voxels (z planes, y, x): sub-voxel on purpose
NAMES = ["Brightfield", "DAPI"]           # channel order: brightfield first, as in the lab series
ROLES = {"dna": 1, BRIGHTFIELD_ROLE: 0}
STEM_SYN = "20261005_DLW93_nohs_gfpsyp2_5min_2hr_herm"
PERIOD = 300.0
STEP_UM = float(np.linalg.norm(STEP_VOX * np.asarray(SP)))


# ---- the synthetic series ----------------------------------------------------------------------------------

def _render(centres, sigma, amp):
    """Gaussian blobs at continuous (z, y, x) voxel positions on the SHAPE grid (each within 4 sigma)."""
    v = np.zeros(SHAPE, np.float32)
    for c in centres:
        lo = [max(0, int(np.floor(ci - 4 * si))) for ci, si in zip(c, sigma)]
        hi = [min(n, int(np.ceil(ci + 4 * si)) + 1) for ci, si, n in zip(c, sigma, SHAPE)]
        if any(h <= lo_ for lo_, h in zip(lo, hi)):
            continue
        zz, yy, xx = np.meshgrid(*[np.arange(a, b, dtype=np.float32) for a, b in zip(lo, hi)], indexing="ij")
        v[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] += np.exp(
            -0.5 * (((zz - c[0]) / sigma[0]) ** 2 + ((yy - c[1]) / sigma[1]) ** 2 + ((xx - c[2]) / sigma[2]) ** 2))
    return v * np.float32(amp)


def make_series(n_t=5, *, flow_vox=0.25, jump=None, jump_in_brightfield=True, rotation_deg_per_step=0.0, seed=0):
    """A synthetic series: 16 nuclei (gaussian blobs, sigma 2 planes x 10 x 10 voxels, about a 2.5 um nucleus at
    the lab's voxel size) on a jittered 4 x 4 grid, and a brightfield channel (a bright field with 24 dark blobs
    of its own); both move by STEP_VOX per step. ``flow_vox`` adds a division-free shear flow to the nuclei only:
    an x velocity per step proportional to the distance from the field centre (zero mean over the field).
    ``jump`` = (t_index, dz, dy, dx) adds a one-off displacement of the whole sample at that step, to the
    brightfield too unless ``jump_in_brightfield`` is False. ``rotation_deg_per_step`` turns the nuclei about the
    z axis through the field centre. Returns stacks (one (2, Z, Y, X) float32 array per timepoint: brightfield,
    nuclei), centres_vox (the nuclei per timepoint) and drift_vox (the cumulative sample displacement)."""
    rng = np.random.default_rng(seed)
    g = np.linspace(44.0, 212.0, 4)
    gy, gx = np.meshgrid(g, g, indexing="ij")
    base = np.column_stack([rng.uniform(9.5, 12.5, 16), gy.ravel() + rng.uniform(-6, 6, 16),
                            gx.ravel() + rng.uniform(-6, 6, 16)])
    bf_base = np.column_stack([rng.uniform(9.0, 13.0, 24), rng.uniform(28, 228, 24), rng.uniform(28, 228, 24)])
    centre = np.array([128.0, 128.0])
    shear = (base[:, 2] - centre[1]) / 84.0          # -1 .. 1 across the grid
    stacks, centres, drifts = [], [], []
    drift = np.zeros(3)
    for t in range(n_t):
        if t > 0:
            drift = drift + STEP_VOX
            if jump is not None and jump[0] == t:
                drift = drift + np.asarray(jump[1:], float)
        bf_drift = drift
        if jump is not None and not jump_in_brightfield and t >= jump[0]:
            bf_drift = drift - np.asarray(jump[1:], float)
        rel = base[:, 1:] - centre
        if rotation_deg_per_step:
            a = np.radians(rotation_deg_per_step * t)
            rel = np.column_stack([rel[:, 0] * np.cos(a) - rel[:, 1] * np.sin(a),
                                   rel[:, 0] * np.sin(a) + rel[:, 1] * np.cos(a)])
        c = np.column_stack([base[:, 0], rel[:, 0] + centre[0], rel[:, 1] + centre[1]]) + drift
        c[:, 2] += flow_vox * shear * t
        nuc = _render(c, (2.0, 10.0, 10.0), 3000.0)
        bf = np.float32(3000.0) - _render(bf_base + bf_drift, (2.0, 6.0, 6.0), 1200.0)
        stacks.append(np.stack([bf, nuc]))
        centres.append(c)
        drifts.append(drift.copy())
    return {"stacks": stacks, "centres_vox": centres, "drift_vox": drifts}


def write_fake_series(series_dir, series, stride, *, brightfield=True, record=True, volumes=True, jitter_um=0.03, seed=1):
    """Write what a series run leaves behind for `germquant register`: per timepoint the reference volumes
    (through series.write_reference_volumes, from a Stack), the nuclei table (the planted centroids in microns
    plus a small segmentation jitter) and the image_summary row (the voxel columns), then the series record and
    the series_summary table. Returns the image ids."""
    cfg = load_config(ROOT / "config" / "config.yaml")
    rng = np.random.default_rng(seed)
    series_dir = Path(series_dir)
    series_dir.mkdir(parents=True, exist_ok=True)
    roles = dict(ROLES) if brightfield else {"dna": 1}
    entries, ids, rows = [], [], []
    for t, (data, c) in enumerate(zip(series["stacks"], series["centres_vox"])):
        sub = series_dir / f"t{t:03d}"
        sub.mkdir(exist_ok=True)
        iid = f"{STEM_SYN}_t{t:03d}"
        ids.append(iid)
        info = None
        if volumes:
            st = Stack(data=data, spacing=SP, channel_names=NAMES, path="x.nd2", spacing_ok=True)
            info = write_reference_volumes(st, roles, cfg, sub, iid, stride, series_dir)
        cen = c * np.asarray(SP) + rng.normal(0.0, jitter_um, c.shape)
        pd.DataFrame({"image_id": iid, "nucleus_id": np.arange(1, len(c) + 1), "volume_um3": 8.0,
                      "centroid_z_um": cen[:, 0], "centroid_y_um": cen[:, 1], "centroid_x_um": cen[:, 2],
                      "in_germline": True}).to_csv(sub / f"{iid}__nuclei.csv", index=False)
        pd.DataFrame([{"image_id": iid, "n_nuclei": len(c), "voxel_dz_um": SP[0], "voxel_dy_um": SP[1],
                       "voxel_dx_um": SP[2]}]).to_csv(sub / f"{iid}__image_summary.csv", index=False)
        entries.append({"t_index": t, "t_seconds": t * PERIOD, "image_id": iid, "out_dir": sub.name, "status": "done",
                        "resumed": False, "qc_pass": True, "qc_flags": [], "failed_stages": [], "elapsed_s": 1.0,
                        "reference_volumes": info})
        rows.append({"t_index": t, "t_seconds": t * PERIOD, "image_id": iid, "n_nuclei": len(c)})
    if record:
        n = len(ids)
        rec = {"series_id": STEM_SYN, "file_path": f"{STEM_SYN}.nd2", "config_hash": cfg.hash, "pipeline_version": "0.1.0",
               "git_sha": "test", "n_timepoints": n, "period_s": PERIOD, "t_seconds": [t * PERIOD for t in range(n)],
               "xy_stride": 1, "z_range": None, "enabled_stages": [], "keep_reference_volume": volumes,
               "reference_stride": stride, "resume": False, "requested": list(range(n)), "done": list(range(n)),
               "failed": [], "elapsed_s": {}, "timepoints": entries, "started": "x", "finished": "y", "complete": True}
        (series_dir / f"{STEM_SYN}__series.json").write_text(json.dumps(rec, indent=1), encoding="utf-8")
        meta = {"series_id": STEM_SYN, "file_path": f"{STEM_SYN}.nd2", "genotype": "DLW93", "sex": "herm", "germ_cell": None,
                "treatment": "nohs", "replicate": "NA", "acquisition_date": "2026-10-05", "n_timepoints": n,
                "period_s": PERIOD, "pipeline_version": "0.1.0", "git_sha": "test", "config_hash": cfg.hash,
                "run_timestamp": "now"}
        summ = pd.DataFrame(rows)
        for k, v in meta.items():
            summ[k] = v
        summ.to_csv(series_dir / f"{STEM_SYN}__series_summary.csv", index=False)
    return ids


def _scale(stride):
    """Microns per reference-volume voxel (z planes are not strided)."""
    return np.array([SP[0], SP[1] * stride, SP[2] * stride])


def _study_config(folder: Path, *, stride=4, keep=True) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    p = folder / "study.yaml"
    p.write_text(f"extends: {(ROOT / 'config' / 'config.yaml').as_posix()}\nsegmentation:\n  nuclei:\n    method: classical\n"
                 f"timelapse:\n  reference_stride: {stride}\n  keep_reference_volume: {'true' if keep else 'false'}\n")
    return p


# ---- the estimate and the reference volume writer ----------------------------------------------------------

def test_displacement_sign_and_cap():
    series = make_series(n_t=2, flow_vox=0.0)
    ref, mov = series["stacks"][0][1][:, ::4, ::4], series["stacks"][1][1][:, ::4, ::4]
    d = R.displacement(ref, mov, upsample=10)
    assert np.allclose(d, STEP_VOX / np.array([1, 4, 4]), atol=0.051)    # positive = the sample moved by that much; exact to the 1/10 voxel quantisation
    d_ph = R.displacement(ref, mov, upsample=10, normalization="phase")
    assert np.allclose(d_ph, STEP_VOX / np.array([1, 4, 4]), atol=1.0)   # whitened: same sign, less exact on smooth blobs
    with pytest.raises(ValueError, match="shapes differ"):
        R.displacement(ref, mov[:, :-1])
    step, capped = R.cap_step([0.0, 3.0, 4.0], 1.0)
    assert capped and np.allclose(step, [0.0, 0.6, 0.8])                   # the direction is kept
    step, capped = R.cap_step([0.0, 3.0, 4.0], 0.0)
    assert not capped and np.allclose(step, [0.0, 3.0, 4.0])               # 0 = no cap
    assert R.cap_step([0.1, 0.2, 0.2], 5.0)[1] is False
    assert np.isnan(R.nn_median(None, np.zeros((2, 3)))) and np.isnan(R.nn_spacing(np.zeros((1, 3))))
    assert R.nn_median([[0, 0, 0], [0, 0, 1]], [[0, 0, 0.5]]) == pytest.approx(0.5)


def test_write_reference_volumes(tmp_path):
    cfg = load_config(ROOT / "config" / "config.yaml")
    data = make_series(n_t=1)["stacks"][0]
    st = Stack(data=data, spacing=SP, channel_names=NAMES, path="x.nd2", spacing_ok=True)
    info = write_reference_volumes(st, ROLES, cfg, tmp_path / "t003", "s_t003", 8, tmp_path)
    seg, bf = tmp_path / "t003" / "s_t003__ref_stride8.npy", tmp_path / "t003" / "s_t003__ref_brightfield_stride8.npy"
    assert seg.exists() and bf.exists()
    assert reference_volume_path(tmp_path / "t003", "s_t003", 8) == seg
    assert reference_volume_path(tmp_path / "t003", "s_t003", 8, BRIGHTFIELD_ROLE) == bf
    a = np.load(seg)
    assert a.dtype == np.uint16 and a.shape == (SHAPE[0], 32, 32)
    assert np.array_equal(a, np.rint(data[1][:, ::8, ::8]).astype(np.uint16))     # the dna channel, rounded
    assert np.array_equal(np.load(bf), np.rint(data[0][:, ::8, ::8]).astype(np.uint16))
    assert info == {"stride": 8, "role": "dna", "channel_index": 1, "role_resolved": True,
                    "spacing_um": [0.3, 0.864, 0.864], "shape": [SHAPE[0], 32, 32], "dtype": "uint16",
                    "files": {"dna": "t003/s_t003__ref_stride8.npy", "brightfield": "t003/s_t003__ref_brightfield_stride8.npy"}}
    # without a brightfield role only the segmentation volume is written; an unresolved role falls back to channel 0
    info2 = write_reference_volumes(st, {"dna": 1}, cfg, tmp_path / "u", "u_t000", 4, tmp_path)
    assert list(info2["files"]) == ["dna"] and info2["shape"] == [SHAPE[0], 64, 64]
    info3 = write_reference_volumes(st, {"dna": None}, cfg, tmp_path / "v", "v_t000", 4, tmp_path)
    assert info3["channel_index"] == 0 and info3["role_resolved"] is False
    # the config keys and their code defaults
    assert reference_settings(cfg) == (True, 8)
    cfg.set("timelapse.reference_stride", 4)
    cfg.set("timelapse.keep_reference_volume", False)
    assert reference_settings(cfg) == (False, 4)
    cfg.set("timelapse.reference_stride", 0)
    with pytest.raises(ValueError, match="timelapse.reference_stride"):
        reference_settings(cfg)
    # without an explicit stride, the run's pixel size picks k so the reference pixel is at most 1 um
    cfg2 = load_config(ROOT / "config" / "config.yaml")
    assert reference_settings(cfg2, 0.108333) == (True, 9)              # full resolution
    assert reference_settings(cfg2, 0.216667) == (True, 4)              # xy stride 2
    assert reference_settings(cfg2, 0.433333) == (True, 2)              # xy stride 4
    assert reference_settings(cfg2, 1.5) == (True, 1)
    cfg2.set("timelapse.reference_pixel_um", 2.0)
    assert reference_settings(cfg2, 0.433333) == (True, 4)


# ---- the registration on the synthetic series ----------------------------------------------------------------

@pytest.mark.parametrize("stride,tol", [(4, 0.2), (8, 0.5)])
def test_register_recovers_the_planted_drift(tmp_path, stride, tol):
    series = make_series()
    write_fake_series(tmp_path, series, stride)
    res = R.register_series(tmp_path)
    assert res["series_id"] == STEM_SYN
    csv = pd.read_csv(res["table_path"])
    assert list(csv.columns) == schema.REGISTER_SHIFTS + schema.SERIES_META
    assert Path(res["table_path"]).name == f"{STEM_SYN}__register_shifts.csv"
    assert Path(res["table_path"]).with_suffix(".parquet").exists()
    assert list(csv["t_index"]) == [0, 1, 2, 3, 4] and csv["error"].fillna("").eq("").all()
    assert (csv["method"] == "cross_correlation").all() and (csv["xy_stride"] == stride).all()
    assert (csv["upsample_factor"] == 10).all() and (csv["channel_role"] == "dna").all()
    assert list(csv["t_seconds"]) == [0.0, 300.0, 600.0, 900.0, 1200.0]
    assert (csv["series_id"] == STEM_SYN).all() and (csv["sex"] == "herm").all() and (csv["n_timepoints"] == 5).all()
    # the cumulative shift is the planted drift, within the tolerance in reference-volume voxels
    scale = _scale(stride)
    for t in range(5):
        est = csv.loc[t, ["dz_um", "dy_um", "dx_um"]].to_numpy(float) / scale
        planted = series["drift_vox"][t] / np.array([1, stride, stride])
        assert np.all(np.abs(est - planted) <= tol), (stride, t, est, planted)
    # t0 is the origin; every later row references the previous timepoint; steps add up to the cumulative shift
    assert pd.isna(csv.loc[0, "reference_t_index"]) and list(csv.loc[1:, "reference_t_index"]) == [0, 1, 2, 3]
    assert (csv.loc[0, ["dz_um", "dy_um", "dx_um", "step_dz_um", "step_dy_um", "step_dx_um"]] == 0).all()
    for ax in "zyx":
        assert np.allclose(csv[f"step_d{ax}_um"].cumsum(), csv[f"d{ax}_um"])
        assert np.allclose(csv[f"shift_{ax}_um"], csv[f"step_d{ax}_um"])            # the declared aliases
        assert np.allclose(csv[f"cum_shift_{ax}_um"], csv[f"d{ax}_um"])
    # the brightfield drifted with the sample: the cross-check agrees, nothing is flagged or capped
    assert not csv["capped"].any() and not csv["disagreement_flag"].any()
    bf = csv.loc[1:, "brightfield_disagreement_um"]
    assert bf.notna().all() and (bf < 0.3).all()
    assert np.allclose(csv["bf_disagreement_um"], csv["brightfield_disagreement_um"], equal_nan=True)
    for ax in "zyx":
        assert np.allclose(csv[f"brightfield_d{ax}_um"], csv[f"d{ax}_um"], atol=0.3)
        assert np.allclose(csv.loc[1:, f"bf_shift_{ax}_um"], csv.loc[1:, f"step_d{ax}_um"], atol=0.3)
    # quality: the step moves every nucleus by about the planted step; after the shift only jitter and flow remain
    q = csv.loc[1:]
    assert q["quality_before_um"].between(0.6 * STEP_UM, 1.4 * STEP_UM).all()
    assert (q["quality_after_um"] < 0.5 * q["quality_before_um"]).all()
    assert np.allclose(csv["nn_disp_before_um"], csv["quality_before_um"], equal_nan=True)
    assert np.allclose(csv["nn_disp_after_um"], csv["quality_after_um"], equal_nan=True)
    # the report
    rep = json.loads(Path(res["report_path"]).read_text(encoding="utf-8"))
    assert Path(res["report_path"]).name == f"{STEM_SYN}__register_report.json"
    par = rep["parameters"]
    assert par == {"reference": "previous", "upsample_factor": 10, "normalization": "none", "taper": 0.0,
                   "max_step_um": 5.0, "disagreement_um": 1.0, "channel_role": "dna", "cross_check_role": "brightfield",
                   "stride": stride, "spacing_um": [round(float(v), 6) for v in scale], "rigid": False, "origin_t_index": 0}
    assert rep["n_timepoints"] == 5 and rep["n_registered"] == 4 and len(rep["steps"]) == 4
    assert rep["capped"] == [] and rep["brightfield_disagreement"] == [] and rep["errors"] == {} and rep["rigid"] is None
    assert rep["quality"]["n_steps"] == 4 and rep["quality"]["median_after_um"] < rep["quality"]["median_before_um"]
    assert rep["cumulative_last_um"] == pytest.approx(list(csv.loc[4, ["dz_um", "dy_um", "dx_um"]]), abs=1e-5)
    assert rep["steps"][0]["t_index"] == 1 and rep["steps"][0]["reference_t_index"] == 0
    assert rep["experimental"] == [] and rep["method"] == "cross_correlation"
    assert res["table"].shape[0] == 5 and list(res["table"].columns) == schema.REGISTER_SHIFTS
    # the whitened variant is available and named in the table; on these smooth blobs it is the less exact one
    res_ph = R.register_series(tmp_path, normalization="phase")
    assert (res_ph["table"]["method"] == "phase_cross_correlation").all()
    assert res_ph["report"]["method"] == "phase_cross_correlation" and res_ph["report"]["parameters"]["normalization"] == "phase"
    with pytest.raises(ValueError, match="normalization"):
        R.register_series(tmp_path, normalization="whitened")
    # the edge taper (off by default) keeps the recovery within the tolerance on a field nothing cuts
    res_tap = R.register_series(tmp_path, taper=0.25)
    assert res_tap["report"]["parameters"]["taper"] == 0.25
    for t in range(5):
        est = res_tap["table"].loc[t, ["dz_um", "dy_um", "dx_um"]].to_numpy(float) / scale
        assert np.all(np.abs(est - series["drift_vox"][t] / np.array([1, stride, stride])) <= tol), (stride, t, est)
    with pytest.raises(ValueError, match="taper"):
        R.register_series(tmp_path, taper=1.5)


def test_register_against_t0(tmp_path):
    series = make_series()
    write_fake_series(tmp_path, series, 4)
    res = R.register_series(tmp_path, reference="t0")
    csv = pd.read_csv(res["table_path"])
    assert list(csv.loc[1:, "reference_t_index"]) == [0, 0, 0, 0]
    for t in range(5):
        est = csv.loc[t, ["dz_um", "dy_um", "dx_um"]].to_numpy(float) / _scale(4)
        planted = series["drift_vox"][t] / np.array([1, 4, 4])
        assert np.all(np.abs(est - planted) <= 0.2), (t, est, planted)
    for ax in "zyx":
        assert np.allclose(csv[f"step_d{ax}_um"].cumsum(), csv[f"d{ax}_um"])
    q = csv.loc[1:]
    assert (q["quality_after_um"] < 0.5 * q["quality_before_um"]).all()      # against t0 the raw displacement grows with t
    assert q["quality_before_um"].is_monotonic_increasing
    assert res["report"]["parameters"]["reference"] == "t0"
    with pytest.raises(ValueError, match="reference"):
        R.register_series(tmp_path, reference="next")
    with pytest.raises(ValueError, match="upsample"):
        R.register_series(tmp_path, upsample=0)


def test_register_caps_and_flags_a_jump(tmp_path):
    series = make_series(jump=(3, 0.0, 0.0, 15.0))           # 1.62 um in x on top of the 0.35 um drift step
    write_fake_series(tmp_path, series, 4)
    res = R.register_series(tmp_path, max_step_um=1.0)
    csv = pd.read_csv(res["table_path"])
    assert list(csv.loc[csv["capped"], "t_index"]) == [3]
    step3 = csv.loc[3, ["step_dz_um", "step_dy_um", "step_dx_um"]].to_numpy(float)
    assert np.linalg.norm(step3) == pytest.approx(1.0, abs=1e-6)
    assert step3[2] > 0.95                                   # shortened along the estimated direction
    assert not csv["disagreement_flag"].any()                # the brightfield jumped with the sample: the two agree
    rep = res["report"]
    assert rep["capped"] == [3] and rep["n_capped"] == 1 and rep["parameters"]["max_step_um"] == 1.0
    s3 = rep["steps"][2]
    assert s3["t_index"] == 3 and s3["capped"] is True and s3["raw_step_um"][2] == pytest.approx((2.7 + 15) * SP[2], abs=0.1)
    assert rep["steps"][1]["raw_step_um"] is None and rep["max_step_um"] == pytest.approx(1.0, abs=1e-6)
    # the cumulative shift carries the capped step from there on (the cost of chaining; the quality shows it)
    assert csv.loc[3, "dx_um"] == pytest.approx(csv.loc[2, "dx_um"] + step3[2], abs=1e-6)
    assert csv.loc[3, "quality_after_um"] > csv.loc[2, "quality_after_um"] + 0.3
    # without the cap the jump is recovered and nothing is flagged
    csv2 = pd.read_csv(R.register_series(tmp_path, max_step_um=0)["table_path"])
    assert not csv2["capped"].any()
    assert csv2.loc[3, "step_dx_um"] == pytest.approx((2.7 + 15) * SP[2], abs=0.2 * 4 * SP[2])
    assert csv2.loc[3, "quality_after_um"] < 0.5 * csv2.loc[3, "quality_before_um"]
    assert cli.main(["register", str(tmp_path), "--max-step-um", "1.0", "--upsample", "10"]) == 0
    assert json.loads((tmp_path / f"{STEM_SYN}__register_report.json").read_text(encoding="utf-8"))["capped"] == [3]


def test_register_flags_a_brightfield_disagreement(tmp_path):
    series = make_series(jump=(2, 0.0, 0.0, 14.0), jump_in_brightfield=False)   # the nuclei jump 1.5 um, the brightfield does not
    write_fake_series(tmp_path, series, 4)
    res = R.register_series(tmp_path)
    csv = pd.read_csv(res["table_path"])
    assert not csv["capped"].any()
    assert list(csv.loc[csv["disagreement_flag"], "t_index"]) == [2]
    assert csv.loc[2, "brightfield_disagreement_um"] == pytest.approx(14 * SP[2], abs=0.1)
    assert csv.loc[2, "bf_disagreement_um"] == csv.loc[2, "brightfield_disagreement_um"]
    assert csv.loc[1, "brightfield_disagreement_um"] < 0.3 and csv.loc[3, "brightfield_disagreement_um"] < 0.3
    # the brightfield cumulative shift follows the drift without the jump; the primary one carries the jump
    assert csv.loc[4, "brightfield_dx_um"] == pytest.approx(4 * STEP_VOX[2] * SP[2], abs=0.15)
    assert csv.loc[4, "dx_um"] == pytest.approx((4 * STEP_VOX[2] + 14) * SP[2], abs=0.15)
    rep = res["report"]
    assert rep["brightfield_disagreement"] == [2] and rep["n_disagreement"] == 1
    assert rep["steps"][1]["disagreement_flag"] is True and rep["steps"][0]["disagreement_flag"] is False
    # a stricter threshold flags more steps, a looser one none
    assert R.register_series(tmp_path, disagreement_um=10.0)["report"]["brightfield_disagreement"] == []


def test_register_without_a_brightfield_volume_and_the_channel_option(tmp_path):
    series = make_series()
    write_fake_series(tmp_path / "nobf", series, 4, brightfield=False)
    res = R.register_series(tmp_path / "nobf")
    csv = pd.read_csv(res["table_path"])
    assert csv["brightfield_disagreement_um"].isna().all() and csv["bf_disagreement_um"].isna().all()
    assert csv[["brightfield_dz_um", "bf_shift_x_um"]].isna().all().all() and not csv["disagreement_flag"].any()
    assert res["report"]["parameters"]["cross_check_role"] is None
    assert np.abs(csv.loc[4, "dx_um"] / _scale(4)[2] - series["drift_vox"][4][2] / 4) <= 0.2
    with pytest.raises(ValueError, match="brightfield"):
        R.register_series(tmp_path / "nobf", channel="brightfield")
    # with both volumes, --channel brightfield registers on the brightfield and cross-checks against the nuclei
    write_fake_series(tmp_path / "bf", series, 4)
    res2 = R.register_series(tmp_path / "bf", channel="brightfield")
    csv2 = pd.read_csv(res2["table_path"])
    assert (csv2["channel_role"] == "brightfield").all()
    assert res2["report"]["parameters"] | {"channel_role": "brightfield", "cross_check_role": "dna"} == res2["report"]["parameters"]
    for t in range(5):
        est = csv2.loc[t, ["dz_um", "dy_um", "dx_um"]].to_numpy(float) / _scale(4)
        assert np.all(np.abs(est - series["drift_vox"][t] / np.array([1, 4, 4])) <= 0.2), t
    with pytest.raises(ValueError, match="foci") as e:
        R.register_series(tmp_path / "bf", channel="foci")
    assert "brightfield" in str(e.value) and "dna" in str(e.value)


def test_register_needs_reference_volumes(tmp_path):
    series = make_series(n_t=3)
    write_fake_series(tmp_path, series, 8, volumes=False)
    with pytest.raises(R.ReferenceVolumesMissing, match="timelapse.keep_reference_volume") as e:
        R.register_series(tmp_path)
    assert "timelapse.reference_stride" in str(e.value) and "--resume" in str(e.value)
    assert issubclass(R.ReferenceVolumesMissing, FileNotFoundError)
    assert cli.main(["register", str(tmp_path)]) == 1                     # the command prints the reason and exits 1
    assert not (tmp_path / f"{STEM_SYN}__register_shifts.csv").exists()
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError, match="no timepoint folders"):
        R.register_series(empty)
    with pytest.raises(FileNotFoundError, match="not a folder"):
        R.register_series(tmp_path / "nowhere")


def test_register_discovers_a_series_without_a_record_and_skips_a_missing_volume(tmp_path):
    series = make_series()
    write_fake_series(tmp_path, series, 4, record=False)
    assert not list(tmp_path.glob("*__series.json"))
    res = R.register_series(tmp_path)
    assert res["series_id"] == STEM_SYN                                       # the stem from the image ids
    csv = pd.read_csv(res["table_path"])
    assert list(csv["t_index"]) == [0, 1, 2, 3, 4] and (csv["channel_role"] == "segmentation").all()
    assert csv["t_seconds"].isna().all() and csv["error"].fillna("").eq("").all()
    assert np.abs(csv.loc[4, "dy_um"] / _scale(4)[1] - series["drift_vox"][4][1] / 4) <= 0.2
    assert res["report"]["parameters"]["cross_check_role"] == "brightfield"   # found by file name
    # a timepoint whose volumes are gone is a row with an error; the chain continues across it
    write_fake_series(tmp_path / "gap", series, 4)
    for p in (tmp_path / "gap" / "t002").glob("*.npy"):
        p.unlink()
    res = R.register_series(tmp_path / "gap")
    csv = pd.read_csv(res["table_path"]).set_index("t_index")
    assert str(csv.loc[2, "error"]).startswith("no reference volume") and pd.isna(csv.loc[2, "dx_um"])
    assert pd.isna(csv.loc[2, "reference_t_index"]) and bool(csv.loc[2, "capped"]) is False
    assert csv.loc[3, "reference_t_index"] == 1 and csv.loc[4, "reference_t_index"] == 3
    assert abs(csv.loc[3, "dx_um"] / _scale(4)[2] - series["drift_vox"][3][2] / 4) <= 0.2
    assert abs(csv.loc[3, "step_dx_um"] - 2 * STEP_VOX[2] * SP[2]) <= 0.2 * 4 * SP[2]   # two frames in one step
    assert res["report"]["errors"] == {"2": "no reference volume"} and res["report"]["n_registered"] == 3


def test_rigid_fit_is_experimental_and_recovers_a_rotation(tmp_path):
    # a pure translation: the Kabsch fit finds no rotation; the shift table is unchanged by --rigid
    series = make_series()
    write_fake_series(tmp_path / "a", series, 4)
    res = R.register_series(tmp_path / "a", rigid=True)
    rep = res["report"]
    assert rep["experimental"] == ["rigid"] and rep["parameters"]["rigid"] is True and len(rep["rigid"]) == 4
    for r in rep["rigid"]:
        assert r["n_linked"] >= 12 and abs(r["rotation_deg"]) < 0.3 and r["rmsd_after_um"] <= r["rmsd_before_um"] + 1e-9
    assert list(pd.read_csv(res["table_path"]).columns) == schema.REGISTER_SHIFTS + schema.SERIES_META
    plain = pd.read_csv(R.register_series(tmp_path / "a")["table_path"])
    pd.testing.assert_frame_equal(plain.drop(columns=["run_timestamp"]),
                                  pd.read_csv(R.register_series(tmp_path / "a", rigid=True)["table_path"]).drop(columns=["run_timestamp"]))
    # a 3 degree turn per step about z is recovered after the translation
    write_fake_series(tmp_path / "b", make_series(rotation_deg_per_step=3.0), 4)
    rep_r = R.register_series(tmp_path / "b", rigid=True)["report"]
    assert len(rep_r["rigid"]) == 4
    for r in rep_r["rigid"]:
        assert r["rotation_deg"] == pytest.approx(3.0, abs=0.5), r
        assert abs(abs(r["rotation_axis"][0]) - 1.0) < 0.05               # about z (components ordered z, y, x)
        assert r["rmsd_after_um"] < 0.5 * r["rmsd_before_um"]
        assert np.allclose(np.asarray(r["rotation_matrix"]) @ np.asarray(r["rotation_matrix"]).T, np.eye(3), atol=1e-6)
    assert cli.main(["register", str(tmp_path / "b"), "--rigid"]) == 0
    # kabsch itself: an exact rotation about the centroid plus a translation small against the point spacing
    rng = np.random.default_rng(3)
    P = rng.uniform(0, 100, (20, 3))
    a = np.radians(7.0)
    Rz = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    Q = (P - P.mean(0)) @ Rz.T + P.mean(0) + np.array([1.0, 2.0, 3.0])
    Rk, tk = R.kabsch(P, Q)
    assert np.allclose(Rk, Rz, atol=1e-9) and np.allclose(P @ Rk.T + tk, Q, atol=1e-9)
    assert R.rotation_angle_deg(Rz) == pytest.approx(7.0)
    fit = R.rigid_fit(P, Q, link_radius_um=25.0)
    assert fit["n_linked"] == 20 and fit["rotation_deg"] == pytest.approx(7.0, abs=1e-6) and fit["rmsd_after_um"] < 1e-9
    assert fit["rotation_axis"] == pytest.approx([1.0, 0.0, 0.0], abs=1e-6)       # about the first axis (z)
    assert R.rigid_fit(P[:3], Q[:3], link_radius_um=25.0)["rotation_deg"] is None     # too few links


# ---- the loop keeps the volumes; the command runs on its output ---------------------------------------------

def test_series_loop_keeps_reference_volumes_and_register_runs(tmp_path, patched_series_reader):  # noqa: F811
    cfg_path = _study_config(tmp_path / "cfg", stride=2)
    out = tmp_path / "series"
    assert cli.main(["run", ND2, "--config", str(cfg_path), "--out", str(out), "--t", "all"]) == 0
    stacks = _series_stacks()
    for k in range(T):
        iid = f"{STEM}_t{k:03d}"
        p = out / f"t{k:03d}" / f"{iid}__ref_stride2.npy"
        assert p.exists()
        arr = np.load(p)
        assert arr.dtype == np.uint16 and arr.shape == (20, 40, 160)         # Z x Y/2 x X/2 of the 20 x 80 x 320 gonad
        # the segmentation channel (dna, channel 0) of that timepoint's stack, strided and rounded
        assert np.array_equal(arr, np.clip(np.rint(stacks[k].data[0][:, ::2, ::2]), 0, 65535).astype(np.uint16))
        assert not (out / f"t{k:03d}" / f"{iid}__ref_brightfield_stride2.npy").exists()   # the n2 map has no brightfield role
    rec = read_series_record(out, STEM)
    assert rec["keep_reference_volume"] is True and rec["reference_stride"] == 2
    e = rec["timepoints"][1]["reference_volumes"]
    assert e == {"stride": 2, "role": "dna", "channel_index": 0, "role_resolved": True,
                 "spacing_um": [GONAD_SPACING[0], GONAD_SPACING[1] * 2, GONAD_SPACING[2] * 2], "shape": [20, 40, 160],
                 "dtype": "uint16", "files": {"dna": f"t001/{STEM}_t001__ref_stride2.npy"}}
    man = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert man["keep_reference_volume"] is True and man["reference_stride"] == 2
    # a fixed image writes none, and its result dict and marker are what they were
    fixed = tmp_path / "fixed"
    res_fixed = pipeline.process_image(ND2, load_config(cfg_path), fixed, prov=None)
    assert not list(fixed.glob("*.npy")) and set(res_fixed) == {"image_id", "n_nuclei", "qc_pass", "qc_flags", "out_dir", "tables", "stages"}
    # register on the loop's output, through the command line
    assert cli.main(["register", str(out), "--reference", "t0"]) == 0
    csv = pd.read_csv(out / f"{STEM}__register_shifts.csv")
    assert (out / f"{STEM}__register_shifts.parquet").exists() and (out / f"{STEM}__register_report.json").exists()
    assert list(csv.columns) == schema.REGISTER_SHIFTS + schema.SERIES_META and list(csv["t_index"]) == [0, 1, 2]
    assert (csv["series_id"] == STEM).all() and (csv["sex"] == "herm").all()
    assert csv["config_hash"].iloc[0] == load_config(cfg_path).hash and (csv["n_timepoints"] == T).all()
    assert list(csv["t_seconds"]) == [0.0, 300.0, 600.0] and list(csv["image_id"]) == [f"{STEM}_t00{k}" for k in range(T)]
    # the planted whole-pixel drift of (1, 3) pixels per step: y and z within 0.2 of a reference-volume voxel. In x
    # this gonad is a poor target on purpose: a bright off-gonad blob is cut by the field edge, and a feature cut by
    # the edge does not move, so it pulls each estimate towards zero (by about a quarter against t0, by about half
    # when two such steps are chained). The sub-voxel series above carries the accuracy claim, this test the plumbing
    for k in range(T):
        assert csv.loc[k, "dy_um"] == pytest.approx(k * DRIFT_PX[0] * GONAD_SPACING[1], abs=0.2 * 2 * GONAD_SPACING[1]), k
        assert csv.loc[k, "dx_um"] == pytest.approx(k * DRIFT_PX[1] * GONAD_SPACING[2], abs=0.4), k
        assert csv.loc[k, "dx_um"] >= 0.5 * k * DRIFT_PX[1] * GONAD_SPACING[2], k
        assert csv.loc[k, "dz_um"] == pytest.approx(0.0, abs=0.2 * GONAD_SPACING[0]), k
    assert csv["brightfield_disagreement_um"].isna().all() and not csv["disagreement_flag"].any() and not csv["capped"].any()
    assert (csv["channel_role"] == "dna").all() and (csv["xy_stride"] == 2).all()
    assert csv.loc[1:, "quality_before_um"].notna().all() and csv["error"].fillna("").eq("").all()
    assert list(csv.loc[1:, "reference_t_index"]) == [0, 0]
    # the default (previous frame, accumulated) through the API name: the same plumbing, the chain of references
    res = api.register_series(out)
    tab = res["table"]
    assert res["report"]["parameters"]["reference"] == "previous" and tab.shape[0] == T
    assert list(tab.loc[1:, "reference_t_index"]) == [0, 1] and (tab.loc[1:, "dx_um"] > 0).all()
    assert tab["dy_um"].to_numpy() == pytest.approx([k * DRIFT_PX[0] * GONAD_SPACING[1] for k in range(T)], abs=0.3 * 2 * GONAD_SPACING[1])


def test_series_loop_picks_the_reference_stride_from_the_pixel_size(tmp_path, patched_series_reader, monkeypatch):  # noqa: F811
    """Without timelapse.reference_stride the loop takes k from the run's pixel size and the 1 um target: the
    synthetic gonad has 0.2 um pixels, so k is 5 at full resolution (1.0 // 0.2 would give 4) and 2 at xy
    stride 2; a header without calibration (spacing_ok false, placeholder spacing) falls back to k 8."""
    p = tmp_path / "nostride.yaml"
    p.write_text(f"extends: {(ROOT / 'config' / 'config.yaml').as_posix()}\nsegmentation:\n  nuclei:\n"
                 f"    method: classical\n")
    cfg = load_config(p)
    assert cfg.get("timelapse.reference_stride") is None
    r1 = process_series(ND2, cfg, tmp_path / "s1", timepoints="0")
    assert r1["record"]["reference_stride"] == 5
    assert (tmp_path / "s1" / "t000" / f"{STEM}_t000__ref_stride5.npy").exists()
    r2 = process_series(ND2, cfg, tmp_path / "s2", timepoints="0", xy_stride=2)
    assert r2["record"]["reference_stride"] == 2
    real = pipeline.read_nd2_metadata
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: {**real(*a, **k), "spacing_ok": False})
    r3 = process_series(ND2, cfg, tmp_path / "s3", timepoints="0")
    assert r3["record"]["reference_stride"] == 8


def test_series_loop_without_reference_volumes_and_resume_reads_them_back(tmp_path, patched_series_reader, monkeypatch):  # noqa: F811
    out = tmp_path / "series"
    cfg_off = load_config(_study_config(tmp_path / "off", stride=2, keep=False))
    res = process_series(ND2, cfg_off, out, timepoints="all")
    assert not list(out.rglob("*.npy")) and res["record"]["keep_reference_volume"] is False
    assert all(e["reference_volumes"] is None for e in res["record"]["timepoints"])
    with pytest.raises(R.ReferenceVolumesMissing, match="timelapse.keep_reference_volume"):
        R.register_series(out)
    # a series run with the key on, then its volumes deleted (a run from before they existed): --resume reads
    # them back from the file, one read per timepoint, without running a stage
    out2 = tmp_path / "series2"
    cfg_on = load_config(_study_config(tmp_path / "on", stride=2))
    first = process_series(ND2, cfg_on, out2, timepoints="all")
    assert all(e["reference_volumes"]["files"] == {"dna": f"t{e['t_index']:03d}/{e['image_id']}__ref_stride2.npy"}
               for e in first["record"]["timepoints"])
    stages_before = (out2 / "t000" / f"{STEM}_t000__stages.json").read_bytes()
    for p in out2.rglob("*.npy"):
        p.unlink()
    patched_series_reader.clear()
    res2 = process_series(ND2, cfg_on, out2, timepoints="all", resume=True)
    assert patched_series_reader == [0, 1, 2]
    assert (out2 / "t000" / f"{STEM}_t000__stages.json").read_bytes() == stages_before
    assert res2["done"] == [0, 1, 2]
    for e in res2["record"]["timepoints"]:
        assert e["resumed"] is True and e["reference_volumes"]["regenerated"] is True
        assert e["reference_volumes"]["shape"] == [20, 40, 160] and e["reference_volumes"]["role"] == "dna"
    assert sorted(p.name for p in out2.rglob("*.npy")) == [f"{STEM}_t{k:03d}__ref_stride2.npy" for k in range(T)]
    # a second resume finds them on disk, reads nothing and keeps the richer entry of the run that wrote them
    patched_series_reader.clear()
    res3 = process_series(ND2, cfg_on, out2, timepoints="all", resume=True)
    assert patched_series_reader == []
    assert all(e["reference_volumes"].get("regenerated") is True for e in res3["record"]["timepoints"])
    assert cli.main(["register", str(out2)]) == 0
    csv = pd.read_csv(out2 / f"{STEM}__register_shifts.csv")
    assert csv.loc[2, "dy_um"] == pytest.approx(2 * DRIFT_PX[0] * GONAD_SPACING[1], abs=0.3 * 2 * GONAD_SPACING[1])
    assert 0 < csv.loc[2, "dx_um"] <= 2 * DRIFT_PX[1] * GONAD_SPACING[2]       # short in x: see the note in the test above
    # a timepoint whose read failed has no volume (the callback never ran); the registration chains across it
    good = pipeline.read_stack

    def boom(path, *, t=None, **k):
        if t == 1:
            raise RuntimeError("boom at t001")
        return good(path, t=t, **k)

    out3 = tmp_path / "series3"
    monkeypatch.setattr(pipeline, "read_stack", boom)
    res4 = process_series(ND2, cfg_on, out3, timepoints="all")
    monkeypatch.setattr(pipeline, "read_stack", good)
    assert res4["failed"][0]["t_index"] == 1 and res4["record"]["timepoints"][1]["reference_volumes"] is None
    res5 = R.register_series(out3)
    csv5 = pd.read_csv(res5["table_path"]).set_index("t_index")
    assert csv5.loc[1, "error"] == "no reference volume (timepoint failed)" and csv5.loc[2, "reference_t_index"] == 0
    assert csv5.loc[2, "dy_um"] == pytest.approx(2 * DRIFT_PX[0] * GONAD_SPACING[1], abs=0.3 * 2 * GONAD_SPACING[1])
    assert csv5.loc[2, "dx_um"] == pytest.approx(2 * DRIFT_PX[1] * GONAD_SPACING[2], abs=0.4)   # one step across the gap, against t0
