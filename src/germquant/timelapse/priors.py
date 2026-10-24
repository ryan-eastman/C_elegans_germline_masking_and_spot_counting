"""Registration priors from the acquisition record (docs/TIMELAPSE.md, "Per-gap caps and the acquisition record";
``germquant register --prior stage``, ``--max-step-at``).

The record. ``load_stage_record`` finds the per-timepoint stage and focus record of a series (io.nd2_stage) and
keeps it next to the series record as ``<stem>__stage_record.json`` (``"format": "germquant.stage_record/1"``): the
source, the stack geometry, the series' own timepoints only and the events between its first and last timepoint. It
is found, in this order:

1. ``nd2`` when given (``--nd2 PATH``): read it and write the record file;
2. otherwise the series record's ``file_path``, when it is an .nd2 that exists: read it and write the record file;
3. otherwise an existing ``<stem>__stage_record.json`` (a folder moved without its nd2, the simulator);
4. otherwise a ValueError that says which of these to give.

The prior of a gap (``gap_prior``) between two registrable timepoints, which need not be consecutive: the stage xy
step in stage coordinates (never mapped to image axes; the mapping and its scale are scope facts), the move of the
stack's centre in Z Coord, the depth both stacks share (restricted to the run's z_range), the refocus pauses in the
gap and the stack depth.

The cap of a gap (``gap_cap``), into timepoint T: ``max_step_at[T]`` when given (0 = no cap there); else, with the
prior, the base cap plus the recorded moves the estimate still has to contain: the stage xy step when it is at least
STAGE_FLOOR_UM (1.0 um; between frames without a re-centring the shown movies' stage moved at most 0.95 um) and the
stack z step when it is at least Z_FLOOR_UM (0.5 um; every ZDrive step of that size falls in a paused gap). A seeded
piecewise step does not contain the stage xy move (the rigid seed carries it), so only the z term is added there. A
base of 0 stays "no cap". Else the base cap.

The flags of a gap (``gap_flags``), advisory only: ``no_record``, ``no_shared_depth`` (the stacks share no depth),
``z_beyond_record`` (the estimate's |dz| is more than Z_EXCESS_UM above the record's |dz|) and, for the rigid step,
``stage_step_mismatch`` (at a recorded stage move the estimate's xy length differs from the stage step's by more than
the base cap, 5 um when the base is 0).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .. import provenance
from ..fsutil import long_path

RECORD_SUFFIX = "__stage_record.json"
RECORD_FORMAT = "germquant.stage_record/1"
PRIORS = ("stage",)
STAGE_FLOOR_UM = 1.0        # a stage xy step counts towards the cap from this length
Z_FLOOR_UM = 0.5            # a stack z step counts towards the cap from this size
Z_EXCESS_UM = 2.5           # z_beyond_record: the estimate's |dz| beyond the record's |dz| by more than this
MISMATCH_NO_CAP_UM = 5.0    # stage_step_mismatch tolerance when the base cap is 0 (no cap)
REFOCUS_TOKEN = "Refocus"   # an event whose description holds this is a refocus pause
FLAG_ORDER = ("no_record", "no_shared_depth", "z_beyond_record", "stage_step_mismatch")
GAP_COLUMNS = ("stage_dx_um", "stage_dy_um", "stage_dz_um", "refocus_pauses", "stack_overlap_um", "prior_flags")


def stage_record_path(series_dir: str | Path, stem: str) -> Path:
    return Path(series_dir) / f"{stem}{RECORD_SUFFIX}"


def restrict_record(rec: dict, t_indices) -> dict:
    """The record of ``t_indices`` only: those timepoints, and the events placed after one of the series' timepoints
    and before its last one (first <= after_t_index < last). Every timepoint when ``t_indices`` is empty or None."""
    keep = None if not t_indices else {int(t) for t in t_indices}
    tps = [dict(tp) for tp in rec.get("timepoints", []) if keep is None or int(tp["t_index"]) in keep]
    ks = [int(tp["t_index"]) for tp in tps]
    lo, hi = (min(ks), max(ks)) if ks else (None, None)
    events = [dict(e) for e in rec.get("events", [])
              if lo is not None and e.get("after_t_index") is not None and lo <= int(e["after_t_index"]) < hi]
    out = {k: v for k, v in rec.items() if k not in ("timepoints", "events")}
    out.update(timepoints=tps, events=events, series_t_indices=ks)
    return out


def write_stage_record(series_dir: str | Path, stem: str, rec: dict) -> Path:
    path = stage_record_path(series_dir, stem)
    body = {"format": RECORD_FORMAT, "series_id": stem}
    body.update({k: v for k, v in rec.items() if k not in ("format", "series_id", "loaded_from", "record_file",
                                                            "written")})
    body["written"] = provenance.run_timestamp()
    with open(long_path(path), "w", encoding="utf-8") as fh:
        json.dump(body, fh, indent=1)
    return path


def read_stage_record_file(path: str | Path) -> dict:
    p = Path(path)
    with open(long_path(p), encoding="utf-8") as fh:
        rec = json.load(fh)
    if not isinstance(rec, dict) or rec.get("format") != RECORD_FORMAT or "timepoints" not in rec:
        raise ValueError(f"{p.name} is not a {RECORD_FORMAT} record")
    return rec


def _nd2_of(record: dict | None) -> Path | None:
    fp = (record or {}).get("file_path")
    if not fp:
        return None
    try:
        p = Path(str(fp))
        return p if p.suffix.lower() == ".nd2" and p.is_file() else None
    except (OSError, ValueError):
        return None


def load_stage_record(series_dir: str | Path, stem: str, record: dict | None, nd2: str | Path | None = None,
                      t_indices=None) -> dict:
    """The stage record of a series, restricted to its timepoints (``t_indices``, else the series record's), found
    by the precedence of the module docstring. The returned dict carries ``loaded_from`` (how it was found) and
    ``record_file`` (the ``<stem>__stage_record.json`` path, written when the record was read from an nd2)."""
    from ..io.nd2_stage import read_stage_record

    sdir = Path(series_dir)
    ks = list(t_indices) if t_indices else [int(e["t_index"]) for e in (record or {}).get("timepoints", [])
                                            if "t_index" in e]
    rec_file = stage_record_path(sdir, stem)
    if nd2 is not None:
        p = Path(nd2)
        if not p.is_file():
            raise FileNotFoundError(f"--nd2 {p}: no such file")
        rec, how = read_stage_record(p), "nd2 (given)"
    elif _nd2_of(record) is not None:
        rec, how = read_stage_record(_nd2_of(record)), "nd2 (the series record's file_path)"
    elif rec_file.is_file():
        out = restrict_record(read_stage_record_file(rec_file), ks)
        out.update(loaded_from="record file", record_file=str(rec_file))
        return out
    else:
        fp = (record or {}).get("file_path")
        raise ValueError(
            f"{sdir}: no acquisition record for --prior stage: the series record's file_path ({fp!r}) is not an "
            f"existing .nd2 and there is no {rec_file.name}; give the movie with --nd2 PATH (read for its metadata "
            f"only) or put a {RECORD_FORMAT} record next to the series record")
    out = restrict_record(rec, ks)
    write_stage_record(sdir, stem, out)
    out.update(loaded_from=how, record_file=str(rec_file))
    return out


def _z_extent(tp: dict, stack: dict, z_range) -> tuple[float, float]:
    """(low, high) Z Coord of the planes a run reads (its z_range, python slice semantics) of one stack."""
    z0, z1 = float(tp["z_first_um"]), float(tp["z_last_um"])
    n = stack.get("n_planes")
    if z_range and n and int(n) > 1:
        n = int(n)
        i0 = min(max(0, int(z_range[0])), n - 1)
        i1 = min(max(i0 + 1, int(z_range[1])), n) - 1
        step = (z1 - z0) / (n - 1)
        za, zb = z0 + i0 * step, z0 + i1 * step
    else:
        za, zb = z0, z1
    return (za, zb) if za <= zb else (zb, za)


def gap_prior(rec: dict, t_a: int, t_b: int, z_range=None) -> dict:
    """The prior of the gap from timepoint ``t_a`` to ``t_b`` (module docstring): ``stage_dx_um``, ``stage_dy_um``
    (stage coordinates), ``stage_xy_um`` (their length), ``stage_dz_um`` (the stack centre's move in Z Coord),
    ``stack_overlap_um`` (negative when the stacks share no depth), ``stack_depth_um``, ``refocus_pauses`` and
    ``recorded`` (False when the record lacks either timepoint; the stage numbers are then NaN)."""
    by_t = {int(tp["t_index"]): tp for tp in rec.get("timepoints", [])}
    a, b = by_t.get(int(t_a)), by_t.get(int(t_b))
    pauses = sum(1 for e in rec.get("events", []) if REFOCUS_TOKEN in str(e.get("description", ""))
                 and e.get("after_t_index") is not None and int(t_a) <= int(e["after_t_index"]) < int(t_b))
    nan = float("nan")
    out = {"t_from": int(t_a), "t_to": int(t_b), "recorded": a is not None and b is not None, "stage_dx_um": nan,
           "stage_dy_um": nan, "stage_xy_um": nan, "stage_dz_um": nan, "stack_overlap_um": nan, "stack_depth_um": nan,
           "refocus_pauses": int(pauses)}
    if a is None or b is None:
        return out
    stack = rec.get("stack") or {}
    lo_a, hi_a = _z_extent(a, stack, z_range)
    lo_b, hi_b = _z_extent(b, stack, z_range)
    dx = float(b["stage_x_um"]) - float(a["stage_x_um"])
    dy = float(b["stage_y_um"]) - float(a["stage_y_um"])
    out.update(stage_dx_um=dx, stage_dy_um=dy, stage_xy_um=float(math.hypot(dx, dy)),
               stage_dz_um=(lo_b + hi_b) / 2.0 - (lo_a + hi_a) / 2.0,
               stack_overlap_um=min(hi_a, hi_b) - max(lo_a, lo_b), stack_depth_um=hi_a - lo_a)
    return out


def gap_cap(t_index: int, base_um: float, *, max_step_at: dict | None = None, gp: dict | None = None,
            seeded: bool = False) -> tuple[float, str]:
    """``(cap, source)`` of the gap into ``t_index`` (module docstring); a cap of 0 means no cap. ``base_um`` is the
    run's --max-step-um (0 for none) and comes back unchanged, the same float, when nothing raises it."""
    if max_step_at and int(t_index) in max_step_at:
        return float(max_step_at[int(t_index)]), "max_step_at"
    if gp is not None and gp.get("recorded") and base_um > 0:
        extra = 0.0
        xy = float(gp["stage_xy_um"])
        if not seeded and xy >= STAGE_FLOOR_UM:
            extra += xy
        dz = abs(float(gp["stage_dz_um"]))
        if dz >= Z_FLOOR_UM:
            extra += dz
        if extra > 0:
            return base_um + extra, "prior:stage"
    return base_um, "max_step_um"


def gap_flags(step_um, gp: dict | None, base_um: float, *, rigid: bool) -> list[str]:
    """The advisory flags of a gap (module docstring) for the estimate ``step_um`` (dz, dy, dx; the uncapped
    estimate, or None when there is none), in FLAG_ORDER."""
    if gp is None:
        return []
    if not gp.get("recorded"):
        return ["no_record"]
    flags = []
    if float(gp["stack_overlap_um"]) <= 0:
        flags.append("no_shared_depth")
    if step_um is not None:
        s = np.asarray(step_um, dtype=float)
        if np.all(np.isfinite(s)):
            if abs(float(s[0])) > abs(float(gp["stage_dz_um"])) + Z_EXCESS_UM:
                flags.append("z_beyond_record")
            if rigid and float(gp["stage_xy_um"]) >= STAGE_FLOOR_UM:
                tol = base_um if base_um > 0 else MISMATCH_NO_CAP_UM
                if abs(float(math.hypot(s[1], s[2])) - float(gp["stage_xy_um"])) > tol:
                    flags.append("stage_step_mismatch")
    return [f for f in FLAG_ORDER if f in flags]


def parse_max_step_at(items) -> dict[int, float]:
    """``{T: UM}`` from ``["T=UM", ...]`` (``--max-step-at``, repeatable). ValueError for a malformed item, a
    negative or non-finite cap, or a timepoint given twice."""
    out: dict[int, float] = {}
    for it in items or []:
        text = str(it).strip()
        if text.count("=") != 1:
            raise ValueError(f"--max-step-at {it!r}: expected T=UM, for example 6=40")
        t_txt, um_txt = (v.strip() for v in text.split("="))
        try:
            t = int(t_txt.lstrip("tT"))
            um = float(um_txt)
        except ValueError:
            raise ValueError(f"--max-step-at {it!r}: expected an integer timepoint and a cap in microns, for "
                             f"example 6=40") from None
        if t < 0 or not math.isfinite(um) or um < 0:
            raise ValueError(f"--max-step-at {it!r}: the timepoint must be >= 0 and the cap a finite number >= 0 "
                             f"(0 = no cap at that gap)")
        if t in out:
            raise ValueError(f"--max-step-at: timepoint {t} given twice")
        out[t] = um
    return out


def check_max_step_at(max_step_at, registrable: list[int]) -> dict[int, float]:
    """``max_step_at`` as ``{int: float}``, every key a registrable timepoint (one with a gap into it). ValueError
    naming the registrable timepoints otherwise: a cap for a gap that does not exist is an error, not a no-op."""
    out = {}
    for t, um in dict(max_step_at or {}).items():
        try:
            k, v = int(t), float(um)
        except (TypeError, ValueError):
            raise ValueError(f"max_step_at: {t!r}={um!r} is not a timepoint and a cap in microns") from None
        if not math.isfinite(v) or v < 0:
            raise ValueError(f"max_step_at: the cap at t{k:03d} must be a finite number >= 0, got {um!r}")
        out[k] = v
    bad = sorted(k for k in out if k not in set(registrable))
    if bad:
        shown = ", ".join(str(k) for k in registrable[:12]) + (" ..." if len(registrable) > 12 else "")
        raise ValueError(f"--max-step-at: {', '.join(f't{k:03d}' for k in bad)} is not a registrable timepoint of "
                         f"this series (a timepoint with a reference volume after the first one: {shown or 'none'})")
    return out


def gap_record(gp: dict | None, flags: list[str]) -> dict:
    """The gap prior's table columns (GAP_COLUMNS) for one row."""
    nan = float("nan")
    if gp is None:
        return {"stage_dx_um": nan, "stage_dy_um": nan, "stage_dz_um": nan, "refocus_pauses": nan,
                "stack_overlap_um": nan, "prior_flags": ""}
    return {"stage_dx_um": gp["stage_dx_um"], "stage_dy_um": gp["stage_dy_um"], "stage_dz_um": gp["stage_dz_um"],
            "refocus_pauses": gp["refocus_pauses"], "stack_overlap_um": gp["stack_overlap_um"],
            "prior_flags": ";".join(flags)}


def report_block(rec: dict, gaps: list[dict], *, z_range=None) -> dict:
    """The ``prior`` block of a registration report: the record and where it came from, the stack, the constants,
    the gaps with a stage move or a refocus pause, and the flags per gap. ``gaps``: one {gap prior, flags, cap,
    cap_source} dict per registered gap."""
    stack = dict(rec.get("stack") or {})
    depth = None
    tps = rec.get("timepoints") or []
    if tps:
        lo, hi = _z_extent(tps[0], stack, z_range)
        depth = hi - lo
    moved = [g for g in gaps if g.get("recorded") and (g["stage_xy_um"] >= STAGE_FLOOR_UM
                                                        or abs(g["stage_dz_um"]) >= Z_FLOOR_UM
                                                        or g["refocus_pauses"] > 0)]
    return {"prior": "stage", "record_file": rec.get("record_file"), "loaded_from": rec.get("loaded_from"),
            "source": rec.get("source"), "stack": stack, "stack_depth_um": depth, "z_range": list(z_range) if z_range else None,
            "n_timepoints_recorded": len(tps), "n_events": len(rec.get("events") or []),
            "events_error": rec.get("events_error"),
            "constants": {"stage_floor_um": STAGE_FLOOR_UM, "z_floor_um": Z_FLOOR_UM, "z_excess_um": Z_EXCESS_UM,
                          "mismatch_no_cap_um": MISMATCH_NO_CAP_UM},
            "gaps_with_a_move_or_refocus": [
                {"gap": f"t{g['t_from']}-t{g['t_to']}", "stage_dx_um": g["stage_dx_um"], "stage_dy_um": g["stage_dy_um"],
                 "stage_dz_um": g["stage_dz_um"], "refocus_pauses": g["refocus_pauses"],
                 "stack_overlap_um": g["stack_overlap_um"], "cap_um": g.get("cap_um"), "cap_source": g.get("cap_source")}
                for g in moved],
            "flags": {f"t{g['t_from']}-t{g['t_to']}": g["flags"] for g in gaps if g.get("flags")}}
