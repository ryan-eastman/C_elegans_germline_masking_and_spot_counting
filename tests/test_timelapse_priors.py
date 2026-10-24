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
import pytest

from germquant.io import nd2_stage as S
from germquant.timelapse import priors as PR

# ---- the nd2 stage record ----------------------------------------------------------------------------------------

T0_MS = 358.0
PERIOD_MS = 300000.0
MU = "�m"                     # the unit suffix decodes as mojibake on the lab's files


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
