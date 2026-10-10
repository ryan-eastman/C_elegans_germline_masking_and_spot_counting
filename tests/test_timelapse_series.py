"""The time-lapse loop (germquant.timelapse.series, docs/TIMELAPSE.md) and process_image's ``t`` keyword.

A synthetic three-timepoint series is built the way test_golden / test_end_to_end build their gonad (the same
3-channel stack) with a planted drift of (1, 3) pixels per timepoint in (y, x), and fed through the reader
seam those fixtures patch (pipeline.read_nd2_metadata / pipeline.read_stack). The fixed-image path is not
touched: test_golden runs its goldens unchanged, and a process_image call without ``t`` writes the same id
and the same marker keys it always did (checked here as well).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_end_to_end import CHANNELS, SPACING, _synthetic_stack, patched_reader  # noqa: F401  (fixture)

from germquant import batch, cli, pipeline, provenance, schema
from germquant.config import load_config
from germquant.io.nd2_reader import Stack
from germquant.timelapse import parse_timepoints, process_series, read_series_record

ROOT = Path(__file__).resolve().parents[1]
ND2 = "20251105_n2_nohs_HERM_001.nd2"
STEM = "20251105_n2_nohs_HERM_001"
T = 3
PERIOD_S = 300.0
DRIFT_PX = (1, 3)                    # (y, x) pixels of drift per timepoint
VOLATILE_COLS = ["run_timestamp"]    # the only per-run column that differs between two identical runs
VOLATILE_KEYS = {"elapsed_s", "finished", "run_timestamp", "started"}


def _series_stacks() -> list[Stack]:
    base = _synthetic_stack()
    _, _, Y, X = base.data.shape
    out = []
    for k in range(T):
        dy, dx = DRIFT_PX[0] * k, DRIFT_PX[1] * k
        data = np.zeros_like(base.data)
        data[:, :, dy:, dx:] = base.data[:, :, : Y - dy, : X - dx]         # whole-pixel shift, zero fill
        out.append(Stack(data=data, spacing=SPACING, channel_names=CHANNELS, path="x.nd2", spacing_ok=True))
    return out


@pytest.fixture
def patched_series_reader(monkeypatch):
    """The reader seam for a T=3 series: returns the list of ``t`` values read_stack was asked for."""
    stacks = _series_stacks()
    meta = {"sizes": {"T": T, "C": 3, "Z": 20, "Y": 80, "X": 320}, "spacing": SPACING, "spacing_ok": True,
            "channel_names": CHANNELS, "dtype": "float32", "is_2d": False,
            "n_timepoints": T, "period_s": PERIOD_S, "t_seconds": [k * PERIOD_S for k in range(T)]}
    calls: list = []

    def read_stack(path, *, xy_stride=1, z_range=None, as_float=False, t=None):
        calls.append(t)
        k = 0 if t is None else int(t)
        if not 0 <= k < T:
            raise ValueError(f"t={k} out of range")
        return stacks[k]

    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    monkeypatch.setattr(pipeline, "read_stack", read_stack)
    return calls


def _cfg():
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    return cfg


def _study_config(tmp_path) -> Path:
    p = tmp_path / "study.yaml"
    p.write_text(f"extends: {(ROOT / 'config' / 'config.yaml').as_posix()}\nsegmentation:\n  nuclei:\n    method: classical\n")
    return p


# ---- parse_timepoints ---------------------------------------------------------------------------------

def test_parse_timepoints():
    assert parse_timepoints("all", 3) == [0, 1, 2]
    assert parse_timepoints(7, 25) == [7] and parse_timepoints("7", 25) == [7] and parse_timepoints(" 7 ", 25) == [7]
    assert parse_timepoints("0-24", 25) == list(range(25)) and parse_timepoints("3-3", 5) == [3]
    assert parse_timepoints([2, 0, 2], 3) == [0, 2] and parse_timepoints(range(2), 3) == [0, 1]
    for bad, n in (("3", 3), (-1, 3), ("0-5", 3), ("2-1", 3), ("x", 3), ("", 3), (True, 3), ([], 3), ("all", 0)):
        with pytest.raises(ValueError):
            parse_timepoints(bad, n)


# ---- schema -------------------------------------------------------------------------------------------

def test_series_tables_are_separate_and_use_time_track_id():
    assert set(schema.SERIES_TABLES) == {"register_shifts", "track_points", "tracks", "track_validation", "series_summary"}
    assert not set(schema.SERIES_TABLES) & set(schema.TABLES)
    assert schema.SERIES_SUMMARY[:2] == ["t_index", "t_seconds"]
    for name in ("track_points", "tracks", "track_validation"):
        assert "time_track_id" in schema.SERIES_TABLES[name], name
    for name, cols in schema.SERIES_TABLES.items():
        assert "track_id" not in cols, f"{name}: track_id is the SC fragment id"
        assert len(set(cols)) == len(cols), name
    assert schema.SERIES_META[0] == "series_id" and "config_hash" in schema.SERIES_META


# ---- the fixed path: no t, no change ------------------------------------------------------------------

def test_process_image_without_t_is_unchanged(tmp_path, monkeypatch, patched_reader):  # noqa: F811
    seen: list[dict] = []
    orig = pipeline.read_stack

    def spy(path, *a, **k):
        seen.append(dict(k))
        return orig(path, *a, **k)

    monkeypatch.setattr(pipeline, "read_stack", spy)
    res = pipeline.process_image(ND2, _cfg(), tmp_path / "r", prov=None)
    assert res["image_id"] == STEM
    assert seen == [{"xy_stride": 1, "z_range": None}]                  # no t keyword reaches the reader
    marker = provenance.read_done_marker(tmp_path / "r", STEM)
    assert marker is not None and "t_index" not in marker and "t_seconds" not in marker
    assert list(marker)[-1] == "finished"
    rec = json.loads((tmp_path / "r" / f"{STEM}__stages.json").read_text(encoding="utf-8"))
    assert "t_index" not in rec and "t_seconds" not in rec
    man = json.loads((tmp_path / "r" / "run_manifest.json").read_text(encoding="utf-8"))
    assert "t_index" not in man


def test_process_image_with_t(tmp_path, patched_series_reader):
    res = pipeline.process_image(ND2, _cfg(), tmp_path / "r", prov=None, t=2)
    assert res["image_id"] == f"{STEM}_t002" and patched_series_reader == [2]
    marker = provenance.read_done_marker(tmp_path / "r", f"{STEM}_t002")
    assert marker["t_index"] == 2 and marker["t_seconds"] == 2 * PERIOD_S and list(marker)[-1] == "finished"
    rec = json.loads((tmp_path / "r" / f"{STEM}_t002__stages.json").read_text(encoding="utf-8"))
    assert rec["image_id"] == f"{STEM}_t002" and rec["t_index"] == 2 and rec["t_seconds"] == 2 * PERIOD_S
    man = json.loads((tmp_path / "r" / "run_manifest.json").read_text(encoding="utf-8"))
    assert man["t_index"] == 2 and man["t_seconds"] == 2 * PERIOD_S
    nuclei = pd.read_csv(tmp_path / "r" / f"{STEM}_t002__nuclei.csv")
    assert (nuclei["image_id"] == f"{STEM}_t002").all()
    assert nuclei["sex"].iloc[0] == "herm" and nuclei["replicate"].iloc[0] == 1     # design columns from the raw stem
    assert "t_index" not in nuclei.columns                                           # no per-image column change
    with pytest.raises(ValueError):
        pipeline.process_image(ND2, _cfg(), tmp_path / "neg", prov=None, t=-1)
    with pytest.raises(ValueError):
        pipeline.process_image(ND2, _cfg(), tmp_path / "far", prov=None, t=3)       # the reader refuses


def test_is_done_checks_the_timepoint(tmp_path):
    cfg = _cfg()
    provenance.write_done_marker(tmp_path, "s_t001", config_hash=cfg.hash, enabled_stages=batch.enabled_stage_names(cfg),
                                 stage_hashes={"segment": "abc"}, extra={"t_index": 1, "t_seconds": 300.0})
    assert batch.is_done(tmp_path, "s_t001", cfg)[0] is True
    assert batch.is_done(tmp_path, "s_t001", cfg, t=1)[0] is True
    done, why = batch.is_done(tmp_path, "s_t001", cfg, t=2)
    assert done is False and "t_index" in why


# ---- the loop -----------------------------------------------------------------------------------------

def test_process_series_three_timepoints(tmp_path, patched_series_reader):
    cfg = _cfg()
    out = tmp_path / "series"
    res = process_series(ND2, cfg, out, timepoints="all")
    assert res["series_id"] == STEM and res["n_timepoints"] == T
    assert res["requested"] == [0, 1, 2] and res["done"] == [0, 1, 2] and res["failed"] == []
    assert patched_series_reader == [0, 1, 2]                      # each timepoint read once, with its t
    assert sorted(p.name for p in out.iterdir() if p.is_dir()) == ["t000", "t001", "t002"]

    nuclei = {}
    for k in range(T):
        sub, iid = out / f"t{k:03d}", f"{STEM}_t{k:03d}"
        for name in schema.TABLES:                                # the fixed-image file set, suffixed id
            assert (sub / f"{iid}__{name}.csv").exists(), name
        for suffix in ("__montage.png", "__nuclei_labels.tif", "__stages.json", "__done.json"):
            assert (sub / f"{iid}{suffix}").exists(), suffix
        assert (sub / "run_manifest.json").exists()
        marker = provenance.read_done_marker(sub, iid)
        assert marker["t_index"] == k and marker["t_seconds"] == k * PERIOD_S
        nuclei[k] = pd.read_csv(sub / f"{iid}__nuclei.csv")
        assert (nuclei[k]["image_id"] == iid).all()

    # the planted drift is in the numbers: timepoint k is the same gonad shifted by k * (1, 3) pixels. The
    # RAD-51 foci are single bright voxels, so the spots sit at exact voxel positions and the median
    # displacement of nearest-neighbour matched spots is the shift itself (the watershed fragments of this
    # blobby synthetic tube move with the Li threshold, so nucleus centroids are not a clean readout).
    from scipy.spatial import cKDTree

    def spots_xyz(k):
        s = pd.read_csv(out / f"t{k:03d}" / f"{STEM}_t{k:03d}__spots.csv")
        return s[["z_um", "y_um", "x_um"]].to_numpy(float)
    ref = spots_xyz(0)
    assert len(ref) >= 6
    for k in (1, 2):
        pts = spots_xyz(k)
        assert len(pts) >= 6
        _, idx = cKDTree(ref).query(pts)
        step = np.median(pts - ref[idx], axis=0)
        assert step[2] == pytest.approx(k * DRIFT_PX[1] * SPACING[2], abs=0.02), ("x", k, step)
        assert step[1] == pytest.approx(k * DRIFT_PX[0] * SPACING[1], abs=0.02), ("y", k, step)
        assert step[0] == pytest.approx(0.0, abs=0.02), ("z", k, step)

    # the series record
    rec = read_series_record(out, STEM)
    assert rec is not None and rec == res["record"]
    assert rec["requested"] == [0, 1, 2] and rec["done"] == [0, 1, 2] and rec["failed"] == []
    assert rec["period_s"] == PERIOD_S and rec["t_seconds"] == [0.0, 300.0, 600.0] and rec["n_timepoints"] == T
    assert set(rec["elapsed_s"]) == {"0", "1", "2"} and all(v > 0 for v in rec["elapsed_s"].values())
    assert rec["config_hash"] == cfg.hash and rec["pipeline_version"] and rec["complete"] is True and rec["finished"]
    assert [e["image_id"] for e in rec["timepoints"]] == [f"{STEM}_t00{k}" for k in range(T)]
    assert all(e["status"] == "done" and e["resumed"] is False and e["failed_stages"] == [] for e in rec["timepoints"])
    man = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert man["series_id"] == STEM and man["timepoints_requested"] == [0, 1, 2]

    # the series summary: one row per timepoint, declared columns first, then bg / snr per resolved role, then the meta block
    ss = pd.read_csv(out / f"{STEM}__series_summary.csv")
    assert (out / f"{STEM}__series_summary.parquet").exists()
    assert list(ss.columns[:len(schema.SERIES_SUMMARY)]) == schema.SERIES_SUMMARY
    assert list(ss.columns[-len(schema.SERIES_META):]) == schema.SERIES_META
    assert list(ss["t_index"]) == [0, 1, 2] and list(ss["t_seconds"]) == [0.0, 300.0, 600.0]
    assert list(ss["image_id"]) == [f"{STEM}_t00{k}" for k in range(T)]
    sig_cols = [c for c in ss.columns if c.startswith("sig_")]
    assert sig_cols == ["sig_dna_bg", "sig_dna_snr", "sig_central_element_bg", "sig_central_element_snr",
                        "sig_foci_bg", "sig_foci_snr"]                 # the three resolved roles, not granule / lamin
    assert ss["sig_dna_snr"].notna().all()
    for k in range(T):
        isum = pd.read_csv(out / f"t{k:03d}" / f"{STEM}_t{k:03d}__image_summary.csv")
        assert ss.loc[k, "n_nuclei"] == isum["n_nuclei"].iloc[0] and ss.loc[k, "n_germline_nuclei"] == isum["n_germline_nuclei"].iloc[0]
        assert ss.loc[k, "mean_spots"] == pytest.approx(isum["mean_spots"].iloc[0]) and bool(ss.loc[k, "qc_pass"]) == bool(isum["qc_pass"].iloc[0])
        assert ss.loc[k, "elapsed_s"] == rec["elapsed_s"][str(k)]
    assert (ss["series_id"] == STEM).all() and (ss["n_timepoints"] == T).all() and (ss["period_s"] == PERIOD_S).all()
    assert (ss["sex"] == "herm").all() and (ss["config_hash"] == cfg.hash).all()
    assert res["summary"].shape[0] == T and str(res["summary_path"]).endswith("__series_summary.csv")
    assert ss["n_nuclei"].dtype.kind in "iu"                      # whole counts in the CSV


def test_process_series_resume_skips_done_timepoints(tmp_path, patched_series_reader):
    cfg = _cfg()
    out = tmp_path / "series"
    first = process_series(ND2, cfg, out, timepoints="0-1")
    assert first["done"] == [0, 1] and not (out / "t002").exists()
    ss1 = pd.read_csv(out / f"{STEM}__series_summary.csv")
    assert len(ss1) == 2
    patched_series_reader.clear()
    res = process_series(ND2, cfg, out, timepoints="all", resume=True)
    assert patched_series_reader == [2]                           # only the missing timepoint was processed
    assert res["done"] == [0, 1, 2] and res["failed"] == []
    rec = res["record"]
    assert [e["resumed"] for e in rec["timepoints"]] == [True, True, False]
    assert rec["elapsed_s"]["0"] == first["record"]["elapsed_s"]["0"]        # carried over from the first pass
    ss = pd.read_csv(out / f"{STEM}__series_summary.csv")
    assert len(ss) == 3 and list(ss["t_index"]) == [0, 1, 2]
    assert ss.loc[0, "n_nuclei"] == ss1.loc[0, "n_nuclei"] and ss.loc[1, "mean_spots"] == pytest.approx(ss1.loc[1, "mean_spots"])
    # a different geometry does not satisfy the markers: everything is reprocessed
    patched_series_reader.clear()
    process_series(ND2, cfg, out, timepoints="all", resume=True, xy_stride=2)
    assert patched_series_reader == [0, 1, 2]
    # without resume, a finished series is simply rerun
    patched_series_reader.clear()
    process_series(ND2, cfg, out, timepoints="1")
    assert patched_series_reader == [1]


def test_process_series_records_a_failed_timepoint_and_goes_on(tmp_path, patched_series_reader, monkeypatch):
    cfg = _cfg()
    out = tmp_path / "series"
    good = pipeline.read_stack

    def boom(path, *, t=None, **k):
        if t == 1:
            raise RuntimeError("boom at t001")
        return good(path, t=t, **k)

    monkeypatch.setattr(pipeline, "read_stack", boom)
    res = process_series(ND2, cfg, out, timepoints="all")
    assert res["done"] == [0, 2]
    assert res["failed"] == [{"t_index": 1, "flags": ["EXCEPTION:RuntimeError: boom at t001"]}]
    assert not (out / "t001" / f"{STEM}_t001__done.json").exists()
    rec = read_series_record(out, STEM)
    assert rec["failed"] == res["failed"] and rec["timepoints"][1]["status"] == "failed" and rec["complete"] is True
    ss = pd.read_csv(out / f"{STEM}__series_summary.csv")
    assert len(ss) == 3 and list(ss["t_index"]) == [0, 1, 2]
    bad = ss.set_index("t_index").loc[1]
    assert bad["qc_pass"] == False and str(bad["qc_flags"]).startswith("EXCEPTION:RuntimeError")  # noqa: E712
    assert pd.isna(bad["n_nuclei"]) and pd.isna(bad["mean_spots"]) and bad["t_seconds"] == PERIOD_S
    assert ss.set_index("t_index").loc[2, "qc_pass"] == True  # noqa: E712
    # the reader fixed, --resume reruns only the failed timepoint
    monkeypatch.setattr(pipeline, "read_stack", good)
    patched_series_reader.clear()
    res = process_series(ND2, cfg, out, timepoints="all", resume=True)
    assert patched_series_reader == [1] and res["done"] == [0, 1, 2] and res["failed"] == []
    ss = pd.read_csv(out / f"{STEM}__series_summary.csv")
    assert ss["qc_pass"].all() and ss["n_nuclei"].notna().all()


def _strip(obj):
    """A JSON value without the per-run keys (timing, timestamps), recursively."""
    if isinstance(obj, dict):
        return {k: _strip(v) for k, v in obj.items() if k not in VOLATILE_KEYS}
    if isinstance(obj, list):
        return [_strip(v) for v in obj]
    return obj


def test_cli_run_single_timepoint_equals_the_loop(tmp_path, patched_series_reader):
    """The design's self-test: ``run --t 1`` writes the same t001 folder, file for file, as the whole loop."""
    import tifffile

    cfg_path = _study_config(tmp_path)
    out_all, out_one = tmp_path / "all", tmp_path / "one"
    assert cli.main(["run", ND2, "--config", str(cfg_path), "--out", str(out_all), "--t", "all"]) == 0
    assert cli.main(["run", ND2, "--config", str(cfg_path), "--out", str(out_one), "--t", "1"]) == 0
    assert sorted(p.name for p in out_one.iterdir() if p.is_dir()) == ["t001"]
    rec = read_series_record(out_one, STEM)
    assert rec["requested"] == [1] and rec["done"] == [1]
    assert len(pd.read_csv(out_one / f"{STEM}__series_summary.csv")) == 1

    a, b = out_all / "t001", out_one / "t001"
    names = sorted(p.name for p in a.iterdir())
    assert names == sorted(p.name for p in b.iterdir()) and f"{STEM}_t001__nuclei.csv" in names
    for name in names:
        fa, fb = a / name, b / name
        if name.endswith(".csv"):
            da = pd.read_csv(fa, low_memory=False, float_precision="round_trip").drop(columns=VOLATILE_COLS, errors="ignore")
            db = pd.read_csv(fb, low_memory=False, float_precision="round_trip").drop(columns=VOLATILE_COLS, errors="ignore")
            pd.testing.assert_frame_equal(da, db, check_exact=True, obj=name)
        elif name.endswith(".parquet"):
            da = pd.read_parquet(fa).drop(columns=VOLATILE_COLS, errors="ignore")
            db = pd.read_parquet(fb).drop(columns=VOLATILE_COLS, errors="ignore")
            pd.testing.assert_frame_equal(da, db, check_exact=True, obj=name)
        elif name.endswith(".json"):
            ja = json.loads(fa.read_text(encoding="utf-8"))
            jb = json.loads(fb.read_text(encoding="utf-8"))
            assert _strip(ja) == _strip(jb), name
        elif name.endswith(".tif"):
            assert np.array_equal(tifffile.imread(fa), tifffile.imread(fb)), name
        else:
            assert fa.read_bytes() == fb.read_bytes(), name


def test_cli_resume_needs_t(tmp_path):
    cfg_path = _study_config(tmp_path)
    with pytest.raises(SystemExit, match="--t"):
        cli.main(["run", ND2, "--config", str(cfg_path), "--out", str(tmp_path / "o"), "--resume"])


def test_collect_stacks_series_summaries(tmp_path, patched_series_reader):
    cfg = _cfg()
    root = tmp_path / "res"
    process_series(ND2, cfg, root / STEM, timepoints="all")
    counts = batch.collect(root)
    assert counts["image_summary"] == 3 and counts["batch_summary"] == 3 and counts["series_summary"] == 3
    summ = pd.read_csv(root / "batch_summary.csv")
    assert sorted(summ["image_id"]) == [f"{STEM}_t00{k}" for k in range(T)]
    stacked = pd.read_csv(root / "batch_series_summary.csv")
    assert len(stacked) == 3 and list(stacked["t_index"]) == [0, 1, 2] and (stacked["series_id"] == STEM).all()
    assert cli.main(["collect", str(root)]) == 0
    # a results tree without a series gets no series file
    other = tmp_path / "fixed"
    other.mkdir()
    counts = batch.collect(other)
    assert "series_summary" not in counts and not (other / "batch_series_summary.csv").exists()
