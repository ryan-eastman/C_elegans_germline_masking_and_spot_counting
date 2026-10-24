"""Registration priors (docs/TIMELAPSE.md, "Per-gap caps and the acquisition record", "Seeding the piecewise
registration", "The folded-axis guard"): the nd2 stage record (germquant.io.nd2_stage), the series-side record and
the gap prior (germquant.timelapse.priors), the per-gap caps and flags of `germquant register`, the seeded piecewise
registration and the folded-axis guard.

The reader is exercised against a fake ND2File (the surface read_stage_record uses: loop_indices, the reader's
_custom_tags(), events(), frame_metadata() and experiment) for a 3 T x 2 P x 4 Z movie with a stage move into t2 and a
refocus pause in the gap t0 to t1; only position 0 may be read. The registrations run on the synthetic series of
test_timelapse_register (a planted jump) and on a wider copy of the piecewise gonad of test_timelapse_piecewise with a
planted whole-field move of (0, +8, +30) um into t2, the brightfield moving with it (a stage re-centring), on top of
the shear flow; the folded axis is that gonad's origin axis turned back at its middle.
"""
from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace as NS

import nd2
import numpy as np
import pandas as pd
import pytest
import tifffile
from test_timelapse_piecewise import BASES as PW_BASES
from test_timelapse_piecewise import L_UM as PW_L
from test_timelapse_piecewise import STEM as PW_STEM
from test_timelapse_piecewise import T as PW_T
from test_timelapse_piecewise import U_UM as PW_U
from test_timelapse_piecewise import X0_UM as PW_X0
from test_timelapse_piecewise import centre_um as pw_centre
from test_timelapse_piecewise import flow_per_step as pw_flow
from test_timelapse_piecewise import write_piecewise_series
from test_timelapse_register import SP, STEM_SYN, make_series, write_fake_series

from germquant import schema
from germquant.io import nd2_stage as S
from germquant.timelapse import find_series, link_series, load_piecewise, register_series_piecewise
from germquant.timelapse import piecewise as PW
from germquant.timelapse import priors as PR
from germquant.timelapse import register as R

# ---- the nd2 stage record ----------------------------------------------------------------------------------------

T0_MS = 358.0
PERIOD_MS = 300000.0
MU = chr(0xFFFD) + "m"                     # the unit suffix decodes as mojibake on the lab's files


class _Rdr:
    def __init__(self, tags):
        self._tags = tags

    def _custom_tags(self):
        return self._tags


class _FakeStageND2:
    """What read_stage_record uses of an nd2.ND2File, for a T x P x Z movie in loop order (T, P, Z). Position 0 sits
    at x 100 (plus 0.01 um per plane, so the median is tested), y -50, with Z Coord base + 1.0 um per plane; position
    1 sits far away (x 900, y 400, z 3000), so a record that reads it is visibly wrong. ``x_move`` moves position 0
    by that much in x into t2, ``z_base`` is the stack base per timepoint, ``zdrive`` the ZDrive per timepoint (None:
    no ZDrive tag)."""

    def __init__(self, *, T=3, P=2, Z=4, x0=100.0, x_move=20.0, z_base=(2000.0, 2000.0, 2001.0),
                 zdrive=(2100.0, 2100.0, 2101.0), events=None, drop=(), short=None, no_tags=False, version=(3, 0)):
        self.T, self.P, self.Z = T, P, Z
        self.sizes = {"T": T, "P": P, "Z": Z, "C": 2, "Y": 8, "X": 8}
        self.version = version
        self._loop = [{"T": t, "P": p, "Z": z} for t in range(T) for p in range(P) for z in range(Z)]
        n = len(self._loop)
        x, y, zc, zd, ex = (np.zeros(n) for _ in range(5))
        for i, d in enumerate(self._loop):
            t, p, z = d["T"], d["P"], d["Z"]
            if p == 0:
                x[i] = x0 + (x_move if t >= 2 else 0.0) + 0.01 * z
                y[i] = -50.0
                zc[i] = z_base[t] + 1.0 * z
                zd[i] = zdrive[t] if zdrive is not None else 0.0
            else:
                x[i], y[i], zc[i], zd[i] = 900.0, 400.0, 3000.0 + z, 3100.0
            ex[i] = 150.0
        tags = {f"X Coord [{MU}]": x, f"Y Coord [{MU}]": y, f"Z Coord [{MU}]": zc, "Exposure Time [ms]": ex,
                "Exposure Time(EGFP) [ms]": np.array([150.0])}
        if zdrive is not None:
            tags[f"Ti2 ZDrive [{MU}]"] = zd
        for prefix in drop:
            tags = {k: v for k, v in tags.items() if not k.startswith(prefix)}
        if short:
            tags = {k: (v[:short] if k.startswith("X Coord") else v) for k, v in tags.items()}
        self._rdr = NS() if no_tags else _Rdr(tags)
        if events is None:
            events = [(0.1, "Experiment Paused for Refocusing"), (100.0, "Experiment Paused for Refocusing"),
                      (110.0, "Experiment Resumed"), (700.0, "Experiment Paused for Refocusing")]
        self._events = events

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    @property
    def loop_indices(self):
        return list(self._loop)

    @property
    def experiment(self):
        return [NS(type="TimeLoop", count=self.T, parameters=NS(periodMs=PERIOD_MS)),
                NS(type="XYPosLoop", count=self.P, parameters=NS()),
                NS(type="ZStackLoop", count=self.Z, parameters=NS(stepUm=1.0, homeIndex=self.Z - 1, bottomToTop=True,
                                                                  deviceName="NIDAQ Piezo Z (name: Piezo Z)"))]

    def frame_metadata(self, index: int):
        d = self._loop[index]
        ms = T0_MS + d["T"] * PERIOD_MS + (d["P"] * self.Z + d["Z"]) * 50.0
        return NS(channels=[NS(time=NS(relativeTimeMs=ms))])

    def events(self, orient="records"):
        rows = [{"Time [s]": 0.37, "Index": 0.0, f"X Coord [{MU}]": 100.0}]       # the one frame row events() keeps
        rows += [{"Time [s]": t, "Events": desc} for t, desc in self._events]
        return rows


@pytest.fixture
def fake_nd2(monkeypatch, tmp_path):
    """Install fake ND2 files by name: ``install({"a.nd2": fake, ...})`` creates the (empty) files under tmp_path and
    makes nd2.ND2File open the fake of that name. Returns the paths."""
    def install(fakes: dict):
        paths = {}
        for name in fakes:
            p = tmp_path / name
            p.write_bytes(b"")
            paths[name] = p
        monkeypatch.setattr(nd2, "ND2File", lambda path: fakes[Path(path).name])
        return paths
    return install


def test_stage_record_reads_position_zero_per_timepoint(fake_nd2):
    p = fake_nd2({"movie.nd2": _FakeStageND2()})["movie.nd2"]
    rec = S.read_stage_record(p)
    tps = rec["timepoints"]
    assert [tp["t_index"] for tp in tps] == [0, 1, 2] and rec["n_frames"] == 3 * 2 * 4
    assert [tp["n_frames"] for tp in tps] == [4, 4, 4]                           # position 0 only
    assert [tp["stage_x_um"] for tp in tps] == pytest.approx([100.015, 100.015, 120.015])    # the median of 4 planes
    assert [tp["stage_y_um"] for tp in tps] == [-50.0, -50.0, -50.0]
    assert [(tp["z_first_um"], tp["z_last_um"]) for tp in tps] == [(2000.0, 2003.0), (2000.0, 2003.0), (2001.0, 2004.0)]
    assert [tp["zdrive_um"] for tp in tps] == [2100.0, 2100.0, 2101.0]
    assert [tp["t_start_s"] for tp in tps] == pytest.approx([0.358, 300.358, 600.358])
    assert rec["stack"] == {"n_planes": 4, "z_step_um": 1.0, "z_device": "NIDAQ Piezo Z (name: Piezo Z)",
                            "bottom_to_top": True, "home_index": 3}
    # the events, in time order, each placed after the last timepoint that started before it
    assert [(e["time_s"], e["description"], e["after_t_index"]) for e in rec["events"]] == [
        (0.1, "Experiment Paused for Refocusing", None), (100.0, "Experiment Paused for Refocusing", 0),
        (110.0, "Experiment Resumed", 0), (700.0, "Experiment Paused for Refocusing", 2)]
    assert rec["events_error"] is None
    src = rec["source"]
    assert src["path"] == str(p) and src["size_bytes"] == 0 and src["nd2_version"] == nd2.__version__
    assert src["reader_version"] == S.READER_VERSION and src["file_version"] == [3, 0]


def test_stage_record_without_zdrive_and_refusals(fake_nd2):
    paths = fake_nd2({"nozd.nd2": _FakeStageND2(zdrive=None), "nox.nd2": _FakeStageND2(drop=("X Coord",)),
                      "short.nd2": _FakeStageND2(short=5), "notags.nd2": _FakeStageND2(no_tags=True)})
    rec = S.read_stage_record(paths["nozd.nd2"])
    assert [tp["zdrive_um"] for tp in rec["timepoints"]] == [None, None, None]
    assert rec["timepoints"][2]["stage_x_um"] == pytest.approx(120.015)
    for name, why in (("nox.nd2", "X Coord"), ("short.nd2", "5 value"), ("notags.nd2", "_custom_tags")):
        with pytest.raises(S.StageRecordUnavailable, match=why) as e:
            S.read_stage_record(paths[name])
        assert nd2.__version__ in str(e.value), name
        assert isinstance(e.value, ValueError)


def test_gap_prior_and_the_cap_rule(fake_nd2):
    rec = S.read_stage_record(fake_nd2({"movie.nd2": _FakeStageND2()})["movie.nd2"])
    g01 = PR.gap_prior(rec, 0, 1)
    assert g01["recorded"] and g01["refocus_pauses"] == 1                 # "Resumed" is not a pause; t0 to t1 only
    assert g01["stage_xy_um"] == pytest.approx(0.0) and g01["stage_dz_um"] == pytest.approx(0.0)
    assert g01["stack_overlap_um"] == pytest.approx(3.0) and g01["stack_depth_um"] == pytest.approx(3.0)
    g12 = PR.gap_prior(rec, 1, 2)
    assert g12["stage_dx_um"] == pytest.approx(20.0) and g12["stage_dy_um"] == 0.0
    assert g12["stage_xy_um"] == pytest.approx(20.0) and g12["stage_dz_um"] == pytest.approx(1.0)
    assert g12["stack_overlap_um"] == pytest.approx(2.0) and g12["refocus_pauses"] == 0
    # a gap across a missing timepoint (an error row): the steps add, the pauses of both gaps count
    g02 = PR.gap_prior(rec, 0, 2)
    assert g02["stage_dx_um"] == pytest.approx(20.0) and g02["refocus_pauses"] == 1
    # the run's z_range (planes 0 to 2 of 4) restricts the shared depth, not the move
    gz = PR.gap_prior(rec, 1, 2, z_range=[0, 3])
    assert gz["stack_overlap_um"] == pytest.approx(1.0) and gz["stage_dz_um"] == pytest.approx(1.0)
    assert gz["stack_depth_um"] == pytest.approx(2.0)
    # a timepoint the record lacks
    g13 = PR.gap_prior(rec, 1, 3)
    assert g13["recorded"] is False and np.isnan(g13["stage_dx_um"])
    assert PR.gap_flags([0.0, 0.0, 1.0], g13, 5.0, rigid=True) == ["no_record"]
    # the cap: the stage xy step counts from 1.0 um, the stack z step from 0.5 um; a seeded step adds no xy
    assert PR.gap_cap(2, 5.0, gp=g12) == (pytest.approx(26.0), "prior:stage")
    assert PR.gap_cap(2, 5.0, gp=g12, seeded=True) == (pytest.approx(6.0), "prior:stage")
    assert PR.gap_cap(1, 5.0, gp=g01) == (5.0, "max_step_um")
    small = dict(g12, stage_xy_um=0.95, stage_dz_um=0.4)
    assert PR.gap_cap(2, 5.0, gp=small) == (5.0, "max_step_um")
    edge = dict(g12, stage_xy_um=1.0, stage_dz_um=-0.5)
    assert PR.gap_cap(2, 5.0, gp=edge) == (pytest.approx(6.5), "prior:stage")
    assert PR.gap_cap(2, 0.0, gp=g12) == (0.0, "max_step_um")                # a base of 0 stays "no cap"
    assert PR.gap_cap(2, 5.0, gp=g12, max_step_at={2: 0.0}) == (0.0, "max_step_at")
    assert PR.gap_cap(2, 5.0, max_step_at={1: 3.0}) == (5.0, "max_step_um")
    base = 5.0
    assert PR.gap_cap(1, base)[0] is base                                    # the same float when nothing raises it
    # the flags
    assert PR.gap_flags([0.0, 0.0, 19.0], g12, 5.0, rigid=True) == []
    assert PR.gap_flags([0.0, 0.0, 0.3], g12, 5.0, rigid=True) == ["stage_step_mismatch"]
    assert PR.gap_flags([0.0, 0.0, 0.3], g12, 5.0, rigid=False) == []          # rigid only
    assert PR.gap_flags([0.0, 0.0, 0.3], g12, 0.0, rigid=True) == ["stage_step_mismatch"]   # 5 um without a base
    assert PR.gap_flags([0.0, 0.0, 15.2], g12, 0.0, rigid=True) == []
    assert PR.gap_flags([3.6, 0.0, 19.5], g12, 5.0, rigid=True) == ["z_beyond_record"]      # 3.6 > 1.0 + 2.5
    assert PR.gap_flags([3.4, 0.0, 19.5], g12, 5.0, rigid=True) == []
    deep = dict(g12, stack_overlap_um=-0.4)
    assert PR.gap_flags([0.0, 0.0, 19.5], deep, 5.0, rigid=True) == ["no_shared_depth"]
    assert PR.gap_flags(None, None, 5.0, rigid=True) == []


def test_parse_and_check_max_step_at():
    assert PR.parse_max_step_at(["6=40", " t9 = 45.5 ", "3=0"]) == {6: 40.0, 9: 45.5, 3: 0.0}
    assert PR.parse_max_step_at(None) == {}
    for bad in ("6", "6=", "a=4", "6=-1", "6=nan", "6=4=5"):
        with pytest.raises(ValueError, match="--max-step-at"):
            PR.parse_max_step_at([bad])
    with pytest.raises(ValueError, match="twice"):
        PR.parse_max_step_at(["6=4", "6=5"])
    assert PR.check_max_step_at({6: 40}, [1, 2, 6]) == {6: 40.0}
    with pytest.raises(ValueError, match="t000 is not a registrable timepoint"):
        PR.check_max_step_at({0: 40}, [1, 2, 6])


def test_stage_record_file_round_trip_and_precedence(fake_nd2, tmp_path):
    paths = fake_nd2({"a.nd2": _FakeStageND2(x0=100.0), "b.nd2": _FakeStageND2(x0=500.0)})
    sdir = tmp_path / "series"
    sdir.mkdir()
    stem = "movie"
    record = {"series_id": stem, "file_path": str(paths["b.nd2"]), "timepoints": [{"t_index": 1}, {"t_index": 2}]}
    # 2. the series record's file_path: read, restricted to the series' own timepoints, and written next to it
    rec = PR.load_stage_record(sdir, stem, record)
    assert rec["loaded_from"] == "nd2 (the series record's file_path)"
    assert [tp["t_index"] for tp in rec["timepoints"]] == [1, 2] and rec["series_t_indices"] == [1, 2]
    assert rec["timepoints"][0]["stage_x_um"] == pytest.approx(500.015)
    assert rec["events"] == []          # the pause after t0 and the one after t2 lie outside t1 to t2
    f = PR.stage_record_path(sdir, stem)
    assert f.name == "movie__stage_record.json" and rec["record_file"] == str(f)
    body = json.loads(f.read_text(encoding="utf-8"))
    assert body["format"] == PR.RECORD_FORMAT and body["series_id"] == stem and "written" in body
    assert [tp["t_index"] for tp in body["timepoints"]] == [1, 2] and "loaded_from" not in body
    # 1. --nd2 wins over the series record's file
    rec_a = PR.load_stage_record(sdir, stem, record, nd2=paths["a.nd2"])
    assert rec_a["loaded_from"] == "nd2 (given)" and rec_a["timepoints"][0]["stage_x_um"] == pytest.approx(100.015)
    # 3. without an nd2, the record file (as written last: a.nd2's)
    moved = dict(record, file_path=str(tmp_path / "gone.nd2"))
    rec_f = PR.load_stage_record(sdir, stem, moved)
    assert rec_f["loaded_from"] == "record file"
    assert rec_f["timepoints"] == rec_a["timepoints"] and rec_f["source"] == rec_a["source"]
    sim = dict(record, file_path="simulated:movie")
    assert PR.load_stage_record(sdir, stem, sim)["loaded_from"] == "record file"
    # a series of t0 to t2 keeps the pause after t0 and drops the one after t2
    rec_all = PR.load_stage_record(sdir, stem, dict(record, timepoints=[{"t_index": k} for k in range(3)]),
                                   nd2=paths["a.nd2"])
    assert [e["after_t_index"] for e in rec_all["events"]] == [0, 0]
    # 4. nothing: the refusal says what to give
    f.unlink()
    with pytest.raises(ValueError, match="--nd2"):
        PR.load_stage_record(sdir, stem, moved)
    with pytest.raises(FileNotFoundError, match="no such file"):
        PR.load_stage_record(sdir, stem, record, nd2=tmp_path / "missing.nd2")
    (sdir / "movie__stage_record.json").write_text(json.dumps({"format": "something else"}), encoding="utf-8")
    with pytest.raises(ValueError, match="not a germquant.stage_record/1 record"):
        PR.load_stage_record(sdir, stem, moved)


# ---- per-gap caps and flags: the rigid registration ---------------------------------------------------------------

# the report keys and step keys of `germquant register` at 0c29949: with no new option the report is main's
MAIN_REPORT_KEYS = ["series_id", "series_dir", "method", "parameters", "n_timepoints", "n_registered", "n_capped",
                    "n_disagreement", "capped", "brightfield_disagreement", "errors", "warnings", "max_step_um",
                    "median_disagreement_um", "quality", "cumulative_last_um", "last_t_index", "steps", "rigid",
                    "experimental", "pipeline_version", "git_sha", "run_timestamp"]
MAIN_STEP_KEYS = ["t_index", "reference_t_index", "step_um", "raw_step_um", "cumulative_um", "capped",
                  "brightfield_step_um", "brightfield_cumulative_um", "disagreement_um", "disagreement_flag",
                  "quality_before_um", "quality_after_um"]
MAIN_PARAMETERS = ["reference", "upsample_factor", "normalization", "taper", "max_step_um", "disagreement_um",
                   "channel_role", "cross_check_role", "stride", "spacing_um", "rigid", "origin_t_index"]
CAP_COLUMNS = ["cap_um", "cap_source"]


def write_record(series_dir: Path, stem: str, n_t: int, *, xy_moves=None, z_moves=None, pauses=(), n_planes=24,
                 z_step=0.3, drop=()) -> Path:
    """A germquant.stage_record/1 file for a synthetic series of ``n_t`` timepoints: the stage at (1000, 500) um,
    moved by ``xy_moves[t] = (dx, dy)`` into t; the stack's base at 2000 um, moved by ``z_moves[t]`` into t;
    ``pauses`` refocus pauses placed after the given timepoints; ``drop`` timepoints the record lacks."""
    xy_moves, z_moves = dict(xy_moves or {}), dict(z_moves or {})
    x, y, z = 1000.0, 500.0, 2000.0
    tps = []
    for t in range(n_t):
        dx, dy = xy_moves.get(t, (0.0, 0.0))
        x, y, z = x + dx, y + dy, z + z_moves.get(t, 0.0)
        if t in drop:
            continue
        tps.append({"t_index": t, "t_start_s": 300.0 * t + 0.4, "stage_x_um": x, "stage_y_um": y, "z_first_um": z,
                    "z_last_um": z + (n_planes - 1) * z_step, "zdrive_um": z + 100.0, "n_frames": n_planes})
    events = [{"time_s": 300.0 * a + 150.0, "description": "Experiment Paused for Refocusing", "after_t_index": a}
              for a in pauses]
    rec = {"format": PR.RECORD_FORMAT, "series_id": stem, "source": {"path": "synthetic", "reader_version": "test"},
           "stack": {"n_planes": n_planes, "z_step_um": z_step, "z_device": "test"}, "n_frames": n_t * n_planes,
           "timepoints": tps, "events": events, "events_error": None}
    p = PR.stage_record_path(series_dir, stem)
    p.write_text(json.dumps(rec, indent=1), encoding="utf-8")
    return p


@pytest.fixture
def jump_series(tmp_path):
    """make_series with a 1.62 um x jump into t3 (as test_register_caps_and_flags_a_jump), and a record of a 1.6 um
    stage move into t3 and a refocus pause in the gap t1 to t2."""
    write_fake_series(tmp_path, make_series(jump=(3, 0.0, 0.0, 15.0)), 4)
    write_record(tmp_path, STEM_SYN, 5, xy_moves={3: (1.6, 0.0)}, pauses=(1,))
    return tmp_path


def test_rigid_without_a_new_option_is_main(jump_series):
    res = R.register_series(jump_series, max_step_um=1.0)
    assert list(res["table"].columns) == schema.REGISTER_SHIFTS
    assert list(pd.read_csv(res["table_path"]).columns) == schema.REGISTER_SHIFTS + schema.SERIES_META
    rep = res["report"]
    assert list(rep) == MAIN_REPORT_KEYS and list(rep["parameters"]) == MAIN_PARAMETERS
    assert all(list(s) == MAIN_STEP_KEYS for s in rep["steps"])
    assert rep["capped"] == [3]                                    # the 5 um default would not, 1.0 does
    same = R.register_series(jump_series, max_step_um=1.0, max_step_at={}, prior=None, nd2=None)
    pd.testing.assert_frame_equal(same["table"], res["table"])
    assert list(same["report"]) == MAIN_REPORT_KEYS


def test_rigid_prior_raises_the_cap_where_the_stage_moved(jump_series):
    plain = R.register_series(jump_series, max_step_um=0)["table"]
    res = R.register_series(jump_series, max_step_um=1.0, prior="stage")
    t = res["table"]
    assert list(t.columns) == schema.REGISTER_SHIFTS + CAP_COLUMNS + list(PR.GAP_COLUMNS)
    assert not t["capped"].any() and res["report"]["capped"] == []
    for ax in "zyx":                                                # nothing capped: every step is the uncapped one
        assert np.array_equal(t[f"step_d{ax}_um"].to_numpy(), plain[f"step_d{ax}_um"].to_numpy())
    assert t.loc[3, "step_dx_um"] == pytest.approx((2.7 + 15) * SP[2], abs=0.2 * 4 * SP[2])
    assert t.loc[3, "cap_um"] == pytest.approx(1.0 + 1.6) and t.loc[3, "cap_source"] == "prior:stage"
    assert list(t.loc[[1, 2, 4], "cap_um"]) == [1.0, 1.0, 1.0]
    assert list(t.loc[[1, 2, 4], "cap_source"]) == ["max_step_um"] * 3
    assert np.isnan(t.loc[0, "cap_um"]) and t.loc[0, "cap_source"] == ""          # the origin has no gap
    assert t.loc[3, "stage_dx_um"] == pytest.approx(1.6) and t.loc[3, "stage_dy_um"] == 0.0
    assert list(t.loc[1:, "refocus_pauses"]) == [0, 1, 0, 0] and (t.loc[1:, "stack_overlap_um"] > 6.8).all()
    assert t["prior_flags"].fillna("").eq("").all()
    rep = res["report"]
    assert rep["parameters"]["prior"] == "stage" and "max_step_at" not in rep["parameters"]
    blk = rep["prior"]
    assert blk["loaded_from"] == "record file" and blk["record_file"].endswith(f"{STEM_SYN}__stage_record.json")
    assert blk["constants"] == {"stage_floor_um": 1.0, "z_floor_um": 0.5, "z_excess_um": 2.5, "mismatch_no_cap_um": 5.0}
    assert [g["gap"] for g in blk["gaps_with_a_move_or_refocus"]] == ["t1-t2", "t2-t3"]
    assert blk["flags"] == {}
    s3 = rep["steps"][2]
    assert s3["t_index"] == 3 and s3["cap_um"] == pytest.approx(2.6) and s3["cap_source"] == "prior:stage"
    assert s3["prior"]["stage_dx_um"] == pytest.approx(1.6) and s3["prior"]["flags"] == []
    # against t0 the gap is the same (from the last registrable frame) and so is the cap
    t0 = R.register_series(jump_series, max_step_um=1.0, prior="stage", reference="t0")["table"]
    assert t0.loc[3, "cap_um"] == pytest.approx(2.6) and not t0["capped"].any()


def test_rigid_max_step_at_and_its_refusals(jump_series):
    res = R.register_series(jump_series, max_step_um=1.0, max_step_at={3: 0.0})
    t = res["table"]
    assert list(t.columns) == schema.REGISTER_SHIFTS + CAP_COLUMNS                   # no prior columns
    assert not t["capped"].any() and np.isnan(t.loc[3, "cap_um"]) and t.loc[3, "cap_source"] == "max_step_at"
    assert list(t.loc[[1, 2, 4], "cap_um"]) == [1.0, 1.0, 1.0]
    assert res["report"]["parameters"]["max_step_at"] == {"3": 0.0} and "prior" not in res["report"]
    t2 = R.register_series(jump_series, max_step_um=0, max_step_at={3: 1.0})["table"]
    assert list(t2.loc[t2["capped"], "t_index"]) == [3] and t2.loc[3, "cap_um"] == 1.0
    assert np.isnan(t2.loc[1, "cap_um"]) and t2.loc[1, "cap_source"] == "max_step_um"      # no base cap
    for bad in ({0: 2.0}, {7: 2.0}):                                  # the origin has no gap; t7 is not in the series
        with pytest.raises(ValueError, match="not a registrable timepoint"):
            R.register_series(jump_series, max_step_at=bad)
    with pytest.raises(ValueError, match="prior"):
        R.register_series(jump_series, prior="focus")
    with pytest.raises(ValueError, match="nd2"):
        R.register_series(jump_series, nd2=jump_series / "x.nd2")
    # with the record gone and no nd2 to read, --prior stage says what to give
    PR.stage_record_path(jump_series, STEM_SYN).unlink()
    with pytest.raises(ValueError, match="--nd2"):
        R.register_series(jump_series, prior="stage")


def test_rigid_flags_from_the_record(tmp_path):
    """Three gaps the record cannot explain: a stage move of 8 um into t1 while the image moved 0.35 um
    (stage_step_mismatch), a 3 um z jump into t2 with no recorded focus step (z_beyond_record), and a recorded focus
    step into t4 deeper than the stack (no_shared_depth). The flags change nothing."""
    write_fake_series(tmp_path, make_series(jump=(2, 10.0, 0.0, 0.0)), 4)
    write_record(tmp_path, STEM_SYN, 5, xy_moves={1: (8.0, 0.0)}, z_moves={4: 8.0})
    plain = R.register_series(tmp_path, max_step_um=0)["table"]
    res = R.register_series(tmp_path, max_step_um=0, prior="stage")
    t = res["table"]
    assert plain.loc[2, "step_dz_um"] > 2.5 + 0.3                    # the planted jump, about 3 um
    assert list(t["prior_flags"].fillna("")) == ["", "stage_step_mismatch", "z_beyond_record", "", "no_shared_depth"]
    assert t["cap_um"].isna().all() and (t.loc[1:, "cap_source"] == "max_step_um").all()     # a base of 0: no cap
    assert t.loc[4, "stack_overlap_um"] == pytest.approx(23 * 0.3 - 8.0)
    for ax in "zyx":
        assert np.array_equal(t[f"step_d{ax}_um"].to_numpy(), plain[f"step_d{ax}_um"].to_numpy())
    assert res["report"]["prior"]["flags"] == {"t0-t1": ["stage_step_mismatch"], "t1-t2": ["z_beyond_record"],
                                               "t3-t4": ["no_shared_depth"]}
    # with a base cap the mismatch tolerance is the base: 8 - 0.35 is above 5 and below 9
    t5 = R.register_series(tmp_path, max_step_um=5.0, prior="stage")["table"]
    assert t5.loc[1, "prior_flags"] == "stage_step_mismatch" and t5.loc[1, "cap_um"] == pytest.approx(13.0)
    t9 = R.register_series(tmp_path, max_step_um=9.0, prior="stage")["table"]
    assert t9.loc[1, "prior_flags"] == ""
    # a timepoint the record lacks keeps the base cap and says so
    write_record(tmp_path, STEM_SYN, 5, xy_moves={1: (8.0, 0.0)}, drop=(3,))
    td = R.register_series(tmp_path, max_step_um=5.0, prior="stage")["table"]
    assert list(td.loc[[3, 4], "prior_flags"]) == ["no_record", "no_record"]
    assert list(td.loc[[3, 4], "cap_um"]) == [5.0, 5.0]


# ---- the seeded piecewise registration ----------------------------------------------------------------------------
# A wider copy of test_timelapse_piecewise's gonad: a 240 x 48 um field, 60 nuclei in two rows (y 12 and 22 um) at x 50
# to 166 um, the shear flow of that test (u = U s / L per step, U = 4.6 um), a whole-field move of (0, +8, +30) um into
# t2 that the brightfield shares (a stage re-centring), and into t4 no flow but a uniform +0.25 um in x.

STEM_W = "20261020_DLW93_heatshock_syp2gfp_5min_2hr_herm"
LAB_SP = (0.4, 0.4, 0.4)                 # the label grid (um); its labels are only looked up, never registered
LAB_SHAPE = (16, 120, 600)               # 6.4 x 48 x 240 um
REF_K = 2
REF_SP = (0.4, 0.8, 0.8)                 # the reference volumes: every 2nd label pixel
REF_SHAPE = (16, 60, 300)
X0_W = 50.0                              # the distal end
U_W = 4.6                                # the proximal shear per step (um)
MOVE = (0.0, 8.0, 30.0)                  # the whole-field move into MOVE_T
MOVE_T = 2
STILL_T = 4                              # into t4: no flow, a uniform shift of the whole field
STILL_DX = 0.25
NT_W = 5
R_NUC = 1.4


def _bases_w(seed=7):
    rng = np.random.default_rng(seed)
    out, amp = {}, {}
    i = 1
    for y in (12.0, 22.0):
        for x in np.arange(X0_W, X0_W + 30 * 4.0, 4.0):
            out[i] = np.array([3.2 + rng.uniform(-0.3, 0.3), y + rng.uniform(-1.0, 1.0), x + rng.uniform(-0.45, 0.45)])
            amp[i] = rng.uniform(0.6, 1.4)
            i += 1
    return out, amp


BASES_W, AMP_W = _bases_w()
L_W = max(c[2] for c in BASES_W.values()) - X0_W


def _field_offset(k: int) -> np.ndarray:
    """The whole-field displacement (sample and worm) from t0 to tk."""
    off = np.zeros(3)
    if k >= MOVE_T:
        off += np.asarray(MOVE, float)
    if k >= STILL_T:
        off += np.array([0.0, 0.0, STILL_DX])
    return off


def _centre_w(identity: int, k: int) -> np.ndarray:
    n_flow = sum(1 for j in range(1, k + 1) if j != STILL_T)
    return (BASES_W[identity] + np.array([0.0, 0.0, n_flow * U_W * (BASES_W[identity][2] - X0_W) / L_W])
            + _field_offset(k))


def planted_step(centre_axis_um: float, k: int) -> np.ndarray:
    """The planted step into tk of a segment whose t0 axis centre is ``centre_axis_um``."""
    flow = 0.0 if k == STILL_T else U_W * centre_axis_um / L_W
    return _field_offset(k) - _field_offset(k - 1) + np.array([0.0, 0.0, flow])


def _render_w(centres, amps, sigma=(0.6, 1.0, 1.0), level=2000.0):
    v = np.full(REF_SHAPE, 100.0, np.float32)
    zz, yy, xx = (np.arange(n) * s for n, s in zip(REF_SHAPE, REF_SP))
    for c, a in zip(centres, amps):
        gz = np.exp(-0.5 * ((zz - c[0]) / sigma[0]) ** 2)
        gy = np.exp(-0.5 * ((yy - c[1]) / sigma[1]) ** 2)
        gx = np.exp(-0.5 * ((xx - c[2]) / sigma[2]) ** 2)
        v += np.float32(level * a) * (gz[:, None, None] * gy[None, :, None] * gx[None, None, :]).astype(np.float32)
    return v


def _brightfield_w(k: int, seed=11):
    """The worm body: dark blobs in a bright field, moved only by the whole-field moves."""
    rng = np.random.default_rng(seed)
    blobs = np.column_stack([rng.uniform(1.0, 5.4, 400), rng.uniform(1.0, 47.0, 400), rng.uniform(1.0, 239.0, 400)])
    return 3000.0 - _render_w(blobs + _field_offset(k), rng.uniform(0.2, 0.5, 400), sigma=(0.8, 1.2, 1.2)) + 100.0


def write_moved_series(root: Path, *, n_t: int = NT_W, axis_fn=None) -> dict:
    """The series folder, laid out as process_series lays one out (labels, reference volumes of the central_element
    and brightfield roles, nuclei tables, image summaries, the series record). ``axis_fn(identity, centre)`` replaces
    the origin's axis_position_um (the folded axis). Returns the folder and ids[k] = {nucleus_id: identity}."""
    out = root / STEM_W
    rng = np.random.default_rng(2026)
    ids, entries = {}, []
    r = int(np.ceil(R_NUC / LAB_SP[2])) + 1
    for k in range(n_t):
        idents = sorted(BASES_W)
        order = rng.permutation(len(idents))
        nid_of = {idents[j]: int(order[j]) + 1 for j in range(len(idents))}      # nucleus_id never the identity
        cents = {i: _centre_w(i, k) for i in idents}
        lab = np.zeros(LAB_SHAPE, np.int32)
        for i in idents:
            c = cents[i]
            lo = [max(0, int(np.floor(c[a] / LAB_SP[a])) - r) for a in range(3)]
            hi = [min(LAB_SHAPE[a], int(np.floor(c[a] / LAB_SP[a])) + r + 1) for a in range(3)]
            zz, yy, xx = np.meshgrid(*[np.arange(lo[a], hi[a]) * LAB_SP[a] for a in range(3)], indexing="ij")
            box = lab[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
            box[(zz - c[0]) ** 2 + (yy - c[1]) ** 2 + (xx - c[2]) ** 2 <= R_NUC ** 2] = nid_of[i]
        sub = out / f"t{k:03d}"
        sub.mkdir(parents=True)
        iid = f"{STEM_W}_t{k:03d}"
        tifffile.imwrite(sub / f"{iid}__nuclei_labels.tif", lab, compression="zlib")
        ref = _render_w([cents[i] for i in idents], [AMP_W[i] for i in idents])
        np.save(sub / f"{iid}__ref_stride{REF_K}.npy", np.clip(np.rint(ref), 0, 65535).astype(np.uint16))
        np.save(sub / f"{iid}__ref_brightfield_stride{REF_K}.npy",
                np.clip(np.rint(_brightfield_w(k)), 0, 65535).astype(np.uint16))
        rows = []
        for i in idents:
            c = cents[i]
            a = axis_fn(i, c) if (axis_fn is not None and k == 0) else c[2] - X0_W
            rows.append({"nucleus_id": nid_of[i], "volume_um3": 4.0 / 3.0 * np.pi * R_NUC ** 3,
                         "n_voxels": int((lab == nid_of[i]).sum()), "centroid_z_um": c[0], "centroid_y_um": c[1],
                         "centroid_x_um": c[2], "in_germline": True, "axis_position_um": a,
                         "axis_position_norm": a / 150.0, "n_spots": 0,
                         "central_element_mean_intensity": 1000.0 * AMP_W[i]})
        nuclei = pd.DataFrame(rows).sort_values("nucleus_id")
        for col in schema.NUCLEI:
            if col not in nuclei.columns:
                nuclei[col] = pd.NA
        nuclei = nuclei[list(schema.NUCLEI) + [c for c in nuclei.columns if c not in schema.NUCLEI]]
        meta = {"image_id": iid, "file_path": "x.nd2", "genotype": "DLW93", "sex": "herm", "germ_cell": "oocyte",
                "treatment": "heatshock", "replicate": "NA", "acquisition_date": "20261020",
                "voxel_dz_um": LAB_SP[0], "voxel_dy_um": LAB_SP[1], "voxel_dx_um": LAB_SP[2], "n_channels": 2,
                "dapi_present": False, "channel_map_json": json.dumps({"central_element": 1, "brightfield": 0}),
                "pipeline_version": "0.1.0", "git_sha": "test", "config_hash": "cafe",
                "run_timestamp": "2026-10-20T00:00:00+00:00"}
        for c, v in meta.items():
            nuclei[c] = v
        nuclei.to_csv(sub / f"{iid}__nuclei.csv", index=False)
        pd.DataFrame([{"n_nuclei": len(nuclei), "n_germline_nuclei": len(nuclei), "mean_spots": 0.0,
                       "total_germline_length_um": 150.0, "qc_pass": True, "qc_flags": "", **meta}]).to_csv(
            sub / f"{iid}__image_summary.csv", index=False)
        (sub / f"{iid}__done.json").write_text(json.dumps({"image_id": iid, "config_hash": "cafe", "t_index": k,
                                                           "t_seconds": k * 300.0, "failed_stages": []}))
        ids[k] = {v: i for i, v in nid_of.items()}
        entries.append({"t_index": k, "t_seconds": k * 300.0, "image_id": iid, "out_dir": sub.name, "status": "done",
                        "reference_volumes": {"stride": REF_K, "role": "central_element", "channel_index": 1,
                                              "role_resolved": True, "spacing_um": list(REF_SP),
                                              "shape": list(REF_SHAPE), "dtype": "uint16",
                                              "files": {"central_element": f"t{k:03d}/{iid}__ref_stride{REF_K}.npy",
                                                        "brightfield": f"t{k:03d}/{iid}__ref_brightfield_stride{REF_K}.npy"}}})
    record = {"series_id": STEM_W, "file_path": "x.nd2", "config_hash": "cafe", "pipeline_version": "0.1.0",
              "git_sha": "test", "n_timepoints": n_t, "period_s": 300.0, "t_seconds": [k * 300.0 for k in range(n_t)],
              "xy_stride": 1, "z_range": None, "requested": list(range(n_t)), "done": list(range(n_t)), "failed": [],
              "timepoints": entries, "complete": True}
    (out / f"{STEM_W}__series.json").write_text(json.dumps(record, indent=1))
    (out / "run_manifest.json").write_text(json.dumps({"config": {"timelapse": {"register": {"piecewise_segments": 4}}}}))
    return {"out": out, "ids": ids}


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


def _step_errors(table: pd.DataFrame) -> pd.DataFrame:
    """Per segment step (t > 0): the distance to the planted step."""
    rows = []
    for r in table[table["t_index"] > 0].itertuples(index=False):
        pl = planted_step(float(r.axis_centre_um), int(r.t_index))
        st = np.array([r.step_dz_um, r.step_dy_um, r.step_dx_um], float)
        rows.append({"t_index": int(r.t_index), "segment": int(r.segment), "err_um": float(np.linalg.norm(st - pl)),
                     "capped": bool(r.capped), "filled": bool(r.filled)})
    return pd.DataFrame(rows)


SEED_COLUMNS = ["seed_dz_um", "seed_dy_um", "seed_dx_um", "seed_shift_y_vox", "seed_shift_x_vox", "mov_crop_y0",
                "mov_crop_y1", "mov_crop_x0", "mov_crop_x1", "seed_source", "cap_basis"]


@pytest.fixture(scope="module")
def moved(tmp_path_factory):
    s = write_moved_series(tmp_path_factory.mktemp("moved"))
    rig = R.register_series(s["out"], max_step_um=0)             # the rigid table: the seed
    s["rigid"] = rig["table"].set_index("t_index")
    s["rigid_path"] = Path(rig["table_path"])
    s["unseeded"] = register_series_piecewise(s["out"], max_step_um=40.0)["table"].copy()
    res = register_series_piecewise(s["out"], max_step_um=10.0, seed="rigid")
    s["seeded"], s["seeded_report"] = res["table"].copy(), res["report"]
    return s


def test_unseeded_piecewise_cannot_cross_the_move(moved):
    err = _step_errors(moved["unseeded"])
    at2 = err[err["t_index"] == MOVE_T]
    # one crop box for both volumes cannot follow a 31 um move: on this fixture the four segment steps at t2 were
    # 19.6, 28.2, 34.2 and 34.4 um off the planted step (cap 40, nothing capped or filled)
    assert (at2["err_um"] > 2.0).all() and not at2["capped"].any()
    assert (err.loc[err["t_index"] == 1, "err_um"] < 0.3).all()          # before the move the crops follow the shear
    assert list(moved["unseeded"].columns) == schema.REGISTER_PIECEWISE    # no seed column without --seed


def test_seeded_piecewise_follows_the_move_and_links_every_nucleus(moved):
    t = moved["seeded"]
    err = _step_errors(t)
    assert (err["err_um"] < 0.3).all(), err.sort_values("err_um").tail(4)
    assert not err["capped"].any() and not err["filled"].any()
    assert list(t.columns) == schema.REGISTER_PIECEWISE + ["cap_um", "cap_source"] + SEED_COLUMNS
    rig = moved["rigid"]
    for r in t[t["t_index"] > 0].itertuples(index=False):
        k, a = int(r.t_index), int(r.reference_t_index)
        seed = rig.loc[k, ["dz_um", "dy_um", "dx_um"]].to_numpy(float) - rig.loc[a, ["dz_um", "dy_um", "dx_um"]].to_numpy(float)
        assert [r.seed_dz_um, r.seed_dy_um, r.seed_dx_um] == pytest.approx(seed.tolist(), abs=1e-9)
        assert r.seed_shift_y_vox == int(np.rint(seed[1] / REF_SP[1]))
        assert r.seed_shift_x_vox == int(np.rint(seed[2] / REF_SP[2]))
        assert (r.mov_crop_y0 - r.crop_y0, r.mov_crop_x0 - r.crop_x0) == (r.seed_shift_y_vox, r.seed_shift_x_vox)
        assert (r.mov_crop_y1 - r.mov_crop_y0, r.mov_crop_x1 - r.mov_crop_x0) == (r.crop_y1 - r.crop_y0,
                                                                                  r.crop_x1 - r.crop_x0)
        assert r.seed_source == "rigid" and r.cap_basis == "departure" and r.cap_um == 10.0
    at2 = t[t["t_index"] == MOVE_T]
    assert (at2["seed_shift_y_vox"] == 10).all() and (at2["seed_shift_x_vox"] >= 38).all()
    assert (t.loc[t["t_index"] == STILL_T, ["seed_shift_y_vox", "seed_shift_x_vox"]] == 0).all().all()
    rep = moved["seeded_report"]
    assert rep["parameters"]["seed"] == "rigid" and "seed_table" not in rep["parameters"]     # only when given
    sb = rep["seed"]
    assert sb["table"] == str(moved["rigid_path"]) and len(sb["sha256"]) == 64
    assert sb["rigid_parameters"]["max_step_um"] is None and sb["rigid_parameters"]["channel_role"] == "central_element"
    assert sb["n_seeded_gaps"] == NT_W - 1 and sb["unseeded_gaps"] == []
    seg = rep["steps"][MOVE_T - 1]["segments"][0]
    seed2 = rig.loc[2, ["dz_um", "dy_um", "dx_um"]].to_numpy(float) - rig.loc[1, ["dz_um", "dy_um", "dx_um"]].to_numpy(float)
    assert seg["seed_um"] == pytest.approx(seed2.tolist(), abs=1e-6) and seed2[1] == pytest.approx(8.0, abs=0.05)
    assert len(seg["mov_crop_yx"]) == 4 and seg["departure_um"] is not None
    # linking in the seeded frame: every nucleus is followed through the move, no link joins two nuclei
    pytest.importorskip("laptrack")
    register_series_piecewise(moved["out"], max_step_um=10.0, seed="rigid")          # the folder's table is the seeded one
    pw = link_series(moved["out"], backend="laptrack", gap=2, max_dist_um=3.0, registration="piecewise", write=False)
    assert _identity_errors(pw["points"], moved["ids"]) == (0, 0)
    assert pw["record"]["n_tracks"] == 60 and pw["record"]["n_full_length"] == 60


def test_sub_voxel_seeds_give_the_unseeded_registration_exactly(moved, tmp_path):
    """Seeds below half a reference voxel move no crop: the same displacement call on the same arrays, so with every
    gap's seed under 0.4 um the seeded registration is the unseeded one, number for number."""
    tab = moved["rigid"].reset_index()[["t_index", "dz_um", "dy_um", "dx_um"]].copy()
    tab["dz_um"], tab["dy_um"], tab["dx_um"] = 0.0, -0.3 * tab["t_index"], 0.35 * tab["t_index"]
    p = tmp_path / "tiny__register_shifts.csv"
    tab.to_csv(p, index=False)
    a = register_series_piecewise(moved["out"], max_step_um=40.0)["table"]
    b = register_series_piecewise(moved["out"], max_step_um=40.0, seed="rigid", seed_table=p)["table"]
    assert (b.loc[b["t_index"] > 0, ["seed_shift_y_vox", "seed_shift_x_vox"]] == 0).all().all()
    pd.testing.assert_frame_equal(a, b[list(a.columns)], check_exact=True)
    assert b.loc[b["t_index"] > 0, "seed_source"].eq("rigid").all()


def test_the_departure_cap_clips_only_the_departure(moved):
    res = register_series_piecewise(moved["out"], max_step_um=1.0, seed="rigid")
    t = res["table"].set_index(["t_index", "segment"])
    n_capped = 0
    for st in res["report"]["steps"]:
        for seg in st["segments"]:
            r = t.loc[(st["t_index"], seg["segment"])]
            step = np.array([r.step_dz_um, r.step_dy_um, r.step_dx_um])
            sd = np.array([0.0, r.seed_dy_um, r.seed_dx_um])                  # z is recorded, never applied
            dep = np.asarray(seg["departure_um"], float)
            if bool(r.capped):
                n_capped += 1
                assert np.linalg.norm(step - sd) == pytest.approx(1.0, abs=1e-6)
                assert step - sd == pytest.approx(dep / np.linalg.norm(dep), abs=1e-5)      # the direction is kept
                assert np.linalg.norm(dep) > 1.0
            else:
                assert np.linalg.norm(dep) <= 1.0 + 1e-9
                assert step == pytest.approx(sd + dep, abs=1e-5)
    assert n_capped >= 4
    at2 = t.xs(MOVE_T, level="t_index")
    assert (at2["step_dy_um"] > 7.0).all() and (at2["step_dx_um"] > 29.0).all()      # the move stays in the step


def test_seed_table_and_missing_rigid_steps(tmp_path):
    s = write_moved_series(tmp_path, n_t=3)
    with pytest.raises(ValueError, match="germquant register"):           # no rigid table: say what to run
        register_series_piecewise(s["out"], seed="rigid")
    with pytest.raises(ValueError, match="seed"):
        register_series_piecewise(s["out"], seed="stage")
    with pytest.raises(ValueError, match="--seed rigid"):
        register_series_piecewise(s["out"], seed_table=tmp_path / "x.csv")
    rig = R.register_series(s["out"], max_step_um=0)["table"]
    # a table without the t2 row: the gap into t2 is unseeded and says why; t0 to t1 is seeded from that table
    alt = rig[rig["t_index"] != 2][["t_index", "dz_um", "dy_um", "dx_um"]].copy()
    alt.loc[alt["t_index"] == 1, "dx_um"] += 0.8                            # one voxel more than the folder's table
    p = tmp_path / "alt__register_shifts.csv"
    alt.to_csv(p, index=False)
    t = register_series_piecewise(s["out"], max_step_um=40.0, seed="rigid", seed_table=p)["table"]
    g1, g2 = t[t["t_index"] == 1], t[t["t_index"] == 2]
    assert np.allclose(g1["seed_dx_um"], float(alt.loc[alt["t_index"] == 1, "dx_um"].iloc[0]), rtol=0, atol=1e-9)
    assert g1["seed_source"].eq("rigid").all() and g1["cap_basis"].eq("departure").all()
    assert g2["seed_source"].eq("none: no rigid step for t002").all() and g2["cap_basis"].eq("step").all()
    assert g2["seed_dx_um"].isna().all() and (g2["seed_shift_x_vox"] == 0).all()
    assert (g2["mov_crop_x0"] == g2["crop_x0"]).all()
    err = _step_errors(t)
    assert (err.loc[err["t_index"] == 1, "err_um"] < 0.3).all()            # the shifted seed changes the crop, not the step
    assert (err.loc[err["t_index"] == 2, "err_um"] > 2.0).any()            # unseeded across the move: wrong again


def test_piecewise_prior_caps_and_warns(moved, tmp_path):
    """With --prior stage an unseeded piecewise registration warns at a stage move beyond the margin (its crops
    cannot follow it) and counts the move in that gap's cap; the seeded one counts only the z step."""
    write_record(moved["out"], STEM_W, NT_W, xy_moves={MOVE_T: (30.0, -8.0)}, n_planes=16, z_step=0.4)
    un = register_series_piecewise(moved["out"], max_step_um=10.0, prior="stage")
    t = un["table"]
    assert list(t.columns) == schema.REGISTER_PIECEWISE + ["cap_um", "cap_source"] + list(PR.GAP_COLUMNS)
    g2 = t[t["t_index"] == MOVE_T]
    assert np.allclose(g2["cap_um"], 10.0 + float(np.hypot(30.0, 8.0)), rtol=0, atol=1e-9)
    assert g2["cap_source"].eq("prior:stage").all() and np.allclose(g2["stage_dx_um"], 30.0, rtol=0, atol=1e-9)
    assert t.loc[t["t_index"] == 1, "cap_um"].eq(10.0).all()
    assert any("t1-t2" in w and "--seed rigid" in w for w in un["report"]["warnings"])
    se = register_series_piecewise(moved["out"], max_step_um=10.0, prior="stage", seed="rigid")
    ts = se["table"]
    assert ts.loc[ts["t_index"] == MOVE_T, "cap_um"].eq(10.0).all()       # the seed carries the xy move
    assert ts.loc[ts["t_index"] == MOVE_T, "cap_source"].eq("max_step_um").all()
    assert not any("--seed rigid" in w for w in se["report"]["warnings"])
    assert se["report"]["prior"]["loaded_from"] == "record file"
    PR.stage_record_path(moved["out"], STEM_W).unlink()


# ---- the folded-axis guard ------------------------------------------------------------------------------------------
# test_timelapse_piecewise's straight gonad (x 10 to 126 um, axis_position_um = x - 10) with the origin's axis turned
# back at its middle (a = x - 10 up to x = 70, then 130 - x): equal-count segments then pair the two limbs.

AXIS_COLUMNS = ["axis_source", "axis_pc1_vy", "axis_pc1_vx", "axis_pc1_cy_um", "axis_pc1_cx_um", "axis_pc1_offset_um"]


def _fold_origin_axis(out: Path) -> None:
    f = out / "t000" / f"{PW_STEM}_t000__nuclei.csv"
    df = pd.read_csv(f)
    x = df["centroid_x_um"].to_numpy(float)
    df["axis_position_um"] = np.where(x <= 70.0, x - 10.0, 130.0 - x)
    df.to_csv(f, index=False)


def _nucleus_errors(out: Path) -> np.ndarray:
    """Per timepoint, the largest error of the interpolated x displacement against the planted flow, over the nuclei
    between the outermost segment centres (beyond them the displacement is held constant by design)."""
    pw = load_piecewise(find_series(out))
    idents = sorted(PW_BASES)
    a0 = pw.axis(np.array([PW_BASES[i] for i in idents]))
    inner = (a0 >= pw.centres.min()) & (a0 <= pw.centres.max())
    worst = []
    for k in range(1, PW_T):
        raw = np.array([pw_centre(i, k) for i in idents])
        planted = np.array([k * pw_flow(i) for i in idents])
        worst.append(float(np.abs(pw.displacement(k, raw)[inner, 2] - planted[inner]).max()))
    return np.array(worst)


@pytest.fixture(scope="module")
def folded(tmp_path_factory):
    s = write_piecewise_series(tmp_path_factory.mktemp("fold"))
    _fold_origin_axis(s["out"])
    return s


def test_guard_off_the_folded_axis_mixes_the_limbs(folded):
    res = register_series_piecewise(folded["out"], max_step_um=6.0)
    assert list(res["table"].columns) == schema.REGISTER_PIECEWISE
    assert "axis_guard" not in res["report"] and "axis_guard" not in res["report"]["parameters"]
    assert res["report"]["parameters"]["axis_source"] == "axis_position_um"
    # the segments pair the two limbs (crops 130, 102, 69 and 41 um wide on this 116 um gonad; spread ratio 7.4), and
    # the interpolated displacement is several um off the planted flow (inner nuclei: 2.9 um at t1 to 10.4 um at t4)
    t1 = res["table"][res["table"]["t_index"] == 1]
    widths = (t1["crop_x1"] - t1["crop_x0"]).to_numpy(float) * 0.8
    assert widths.max() > 100.0 and (widths > 60.0).sum() >= 3
    assert (_nucleus_errors(folded["out"]) > 2.0).all()


def test_guard_on_uses_the_straight_axis_and_track_link_rebuilds_it(folded):
    res = register_series_piecewise(folded["out"], max_step_um=6.0, axis_guard=True)
    t, rep = res["table"], res["report"]
    assert list(t.columns) == schema.REGISTER_PIECEWISE + AXIS_COLUMNS
    assert t["axis_source"].eq("pc1").all()
    for c in AXIS_COLUMNS[1:]:
        assert t[c].nunique() == 1 and np.isfinite(t[c]).all()
    g = rep["axis_guard"]
    assert rep["parameters"]["axis_guard"] is True and rep["parameters"]["axis_source"] == "pc1"
    assert g["tripped"] is True and g["axis_used"] == "pc1" and g["threshold"] == 4.0
    assert 6.0 < g["spread_ratio"] < 9.0 and g["spread_ratio"] == pytest.approx(max(g["spread_ratio_by_segment"]))
    assert g["axis_before"] == "axis_position_um" and len(g["segment_spans_before_um"]) == 4
    pc1 = g["pc1"]
    assert abs(pc1["vx"]) > 0.99 and pc1["sign_rule"].startswith("spearman")        # the gonad runs along x
    assert [pc1["vy"], pc1["vx"], pc1["cy_um"], pc1["cx_um"], pc1["offset_um"]] == [
        t[c].iloc[0] for c in AXIS_COLUMNS[1:]]
    # each segment is a stretch of the straight gonad again: its step is the shear at its nuclei's median x
    pw = load_piecewise(find_series(folded["out"]))
    assert pw.axis_source == "pc1" and pw.summary()["axis_source"] == "pc1"
    idents = sorted(PW_BASES)
    c0 = np.array([PW_BASES[i] for i in idents])
    seg = PW.assign_segments(pw.axis(c0), pw.edges)
    for r in t[t["t_index"] > 0].itertuples(index=False):
        xs = c0[seg == int(r.segment), 2]
        planted = PW_U * (float(np.median(xs)) - PW_X0) / PW_L
        assert abs(r.step_dx_um - planted) < 0.3, (r.t_index, r.segment, r.step_dx_um, planted)
        assert (xs.max() - xs.min()) < 40.0                                     # one limb per segment
    assert (_nucleus_errors(folded["out"]) < 0.3 * np.arange(1, PW_T) + 1e-9).all()
    # track-link's axis is the registration's: its registered coordinates give the table's own quality numbers
    lay = find_series(folded["out"])
    for k in range(1, PW_T):
        cm = PW.germline_points(PW.read_nuclei_table(lay.sub(k), lay.image_id(k)), k)[["zr_um", "yr_um", "xr_um"]].to_numpy(float)
        cr = PW.germline_points(PW.read_nuclei_table(lay.sub(k - 1), lay.image_id(k - 1)), k - 1)[
            ["zr_um", "yr_um", "xr_um"]].to_numpy(float)
        reg_m, a_m = pw.register(k, cm)
        reg_r, _a = pw.register(k - 1, cr)
        sm = PW.assign_segments(a_m, pw.edges)
        for s_ in range(4):
            q = R.nn_median(reg_m[sm == s_], reg_r)
            row = t[(t["t_index"] == k) & (t["segment"] == s_)].iloc[0]
            assert q == pytest.approx(float(row["quality_after_um"]), abs=1e-9)
            assert int(row["n_nuclei"]) == int((sm == s_).sum())


def test_guard_on_a_straight_axis_changes_nothing(tmp_path):
    s = write_piecewise_series(tmp_path)
    off = register_series_piecewise(s["out"], max_step_um=6.0)
    on = register_series_piecewise(s["out"], max_step_um=6.0, axis_guard=True)
    g = on["report"]["axis_guard"]
    assert g["tripped"] is False and g["axis_used"] == "axis_position_um" and g["pc1"] is None
    assert g["spread_ratio"] < 1.5
    assert list(on["table"].columns) == schema.REGISTER_PIECEWISE + ["axis_source"]
    assert on["table"]["axis_source"].eq("axis_position_um").all()
    pd.testing.assert_frame_equal(off["table"], on["table"][list(off["table"].columns)], check_exact=True)
    pw = load_piecewise(find_series(s["out"]))                       # not pc1: the axis is rebuilt as on main
    assert pw.axis_source == "axis_position_um"


def test_guard_keys_and_threshold(folded):
    cfg = {"timelapse": {"register": {"piecewise_axis_guard": True, "piecewise_axis_fold_ratio": 10.0}}}
    rep = register_series_piecewise(folded["out"], max_step_um=6.0, cfg=cfg)["report"]
    assert rep["axis_guard"]["threshold"] == 10.0 and rep["axis_guard"]["tripped"] is False
    rep2 = register_series_piecewise(folded["out"], max_step_um=6.0, cfg=cfg, axis_fold_ratio=4.0)["report"]
    assert rep2["axis_guard"]["tripped"] is True                       # an explicit value wins over the key
    assert PW.piecewise_options(None) == {"seed": "none", "axis_guard": False, "axis_fold_ratio": 4.0}
    assert PW.piecewise_options({"timelapse": {"register": {"piecewise_seed": "rigid"}}})["seed"] == "rigid"
    for bad in ({"seed": "stage"}, {"axis_fold_ratio": 0.0}, {"axis_fold_ratio": "x"}):
        with pytest.raises(ValueError):
            PW.piecewise_options(None, **bad)
