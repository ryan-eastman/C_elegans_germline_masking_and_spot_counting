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
from test_timelapse_register import SP, STEM_SYN, make_series, write_fake_series

from germquant import schema
from germquant.io import nd2_stage as S
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
