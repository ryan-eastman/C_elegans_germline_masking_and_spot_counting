"""Regression tests for the tranche 3 integration fixes.

* the fixed reference axis of the time-lapse tracks (track.fixed_axis_positions): the per-frame axis of the
  axis stage is refitted at every timepoint and its origin moves, so a velocity measured on it mixes flow
  with origin jitter; the fixed axis does not;
* ``--image file`` reads the series file on the run's own grid (run xy stride times the export stride, the
  run's z range);
* a bare model name that is not a file, not registered and not a Cellpose built-in is refused instead of
  silently becoming Cellpose's default weights;
* net3d refuses a patch whose Y and X differ (the y/x transpose augmentation swaps them).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from germquant.timelapse import track
from germquant.timelapse.layout import SeriesLayout


def _curved_series(n_nuclei=300, n_t=4, v_um=3.0, jitter=(0.0, 9.0, -6.0, 12.0), seed=3):
    """Nuclei on a gently curved centreline (a circular arc of radius 120 um) that all move ``v_um`` along it
    per timepoint. axis_position_um is the true arc length at t0 and the arc length plus a per-frame origin
    offset (``jitter``) later, as the axis stage's refitted, re-anchored axis would report it."""
    rng = np.random.default_rng(seed)
    R = 120.0
    s0 = np.sort(rng.uniform(5.0, 150.0, n_nuclei))
    off = rng.normal(0.0, 1.0, size=(n_nuclei, 2))                 # scatter across the tube (y and z)
    rows = []
    for k in range(n_t):
        s = s0 + k * v_um
        theta = s / R
        x = R * np.sin(theta)
        y = R * (1.0 - np.cos(theta)) + off[:, 0]
        z = 10.0 + off[:, 1]
        for i in range(n_nuclei):
            rows.append({"t_index": k, "nucleus_id": i + 1, "time_track_id": i + 1, "in_germline": True,
                         "zr_um": z[i], "yr_um": y[i], "xr_um": x[i],
                         "axis_position_um": s[i] + jitter[k], "true_s": s[i]})
    return pd.DataFrame(rows)


def test_fixed_axis_follows_the_first_timepoint_axis_through_origin_jitter():
    pts = _curved_series()
    fixed = track.fixed_axis_positions(pts)
    assert fixed.notna().all()
    # the true arc length, within the local-linear smoothing error on a curved tube at a realistic density
    # (300 nuclei over 145 um, like the October series' ~350 over ~240 um)
    t0 = pts["t_index"] == 0
    assert np.abs(fixed[t0] - pts.loc[t0, "true_s"]).median() < 0.15
    later = ~t0 & (pts["true_s"] < 150.0)
    assert np.abs(fixed[later] - pts.loc[later, "true_s"]).median() < 0.2
    minutes = np.arange(4) * 5.0
    pts["axis_position_fixed_um"] = fixed
    pts["axis_position_aligned_um"] = pts["axis_position_um"]
    v_true = 3.0 / 5.0
    vf, vr = [], []
    for _tid, g in pts.groupby("time_track_id"):
        g = g.sort_values("t_index")
        vf.append(track._axis_velocity(minutes, g))
        vr.append(track._slope_per_min(minutes, g["axis_position_aligned_um"].to_numpy(float)))
    assert np.median(vf) == pytest.approx(v_true, rel=0.02)
    assert abs(np.median(vr) - v_true) > 0.2 * v_true                # the refitted axis is biased by the jitter


def test_fixed_axis_falls_back_to_a_centerline_without_axis_positions():
    pts = _curved_series().drop(columns=["axis_position_um"])
    fixed = track.fixed_axis_positions(pts)
    assert fixed.notna().all()
    minutes = np.arange(4) * 5.0
    pts["axis_position_fixed_um"] = fixed
    v = [abs(track._axis_velocity(minutes, g.sort_values("t_index"))) for _t, g in pts.groupby("time_track_id")]
    assert np.median(v) == pytest.approx(0.6, rel=0.05)          # orientation is arbitrary without axis positions


def test_fixed_axis_needs_five_reference_nuclei():
    pts = _curved_series(n_nuclei=4)
    assert track.fixed_axis_positions(pts).isna().all()


def test_image_file_reads_the_series_on_the_runs_grid(tmp_path, monkeypatch):
    import germquant.io as gio

    series_file = tmp_path / "series.nd2"
    series_file.write_bytes(b"x")
    layout = SeriesLayout(out_dir=tmp_path, stem="s", record={"xy_stride": 2, "z_range": [3, 9]}, timepoints=[0],
                          file_path=str(series_file))
    seen = {}

    class _Stack:
        n_channels = 2
        data = np.zeros((2, 6, 4, 4), np.uint16)

    def fake_read_image(path, xy_stride=1, z_range=None, as_float=False, t=None):
        seen.update(xy_stride=xy_stride, z_range=z_range, t=t)
        return _Stack()

    monkeypatch.setattr(gio, "read_image", fake_read_image)
    vol = track._read_channel_from_file(layout, 0, 3, 1)
    assert vol.shape == (6, 4, 4)
    assert seen == {"xy_stride": 6, "z_range": (3, 9), "t": 0}     # run stride 2 x export stride 3
    assert layout.run_xy_stride() == 2 and layout.run_z_range() == (3, 9)
    assert SeriesLayout(out_dir=tmp_path, stem="s", record=None, timepoints=[0]).run_xy_stride() == 1


def test_unregistered_bare_model_name_is_refused(tmp_path, monkeypatch):
    from germquant import pipeline
    from germquant.config import load_config

    monkeypatch.setenv("GERMQUANT_MODELS", str(tmp_path / "empty_models"))
    root = Path(__file__).resolve().parents[1]
    prof = tmp_path / "p.yaml"
    prof.write_text(f"extends: {(root / 'config' / 'config.yaml').as_posix()}\n"
                    "segmentation:\n  nuclei:\n    cellpose_model: live_syp2_nuclei_not_registered\n", encoding="utf-8")
    cfg = load_config(prof)
    with pytest.raises(FileNotFoundError, match="not a file, not a model registered"):
        pipeline.resolve_model_path(cfg)
    prof.write_text(f"extends: {(root / 'config' / 'config.yaml').as_posix()}\n"
                    "segmentation:\n  nuclei:\n    cellpose_model: cpsam\n", encoding="utf-8")
    assert pipeline.resolve_model_path(load_config(prof)) == "cpsam"


def test_net3d_refuses_a_non_square_patch(tmp_path):
    pytest.importorskip("torch")
    from germquant.spots.net3d import train

    with pytest.raises(ValueError, match="Y and X must be equal"):
        train.train_spots(tmp_path / "cache", tmp_path / "split.json", tmp_path / "out", patch="8,16,32", epochs=1)
