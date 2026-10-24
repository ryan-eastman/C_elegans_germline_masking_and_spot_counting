"""The acquisition record of an nd2 time series, read from the file's metadata without a pixel (docs/TIMELAPSE.md,
"Per-gap caps and the acquisition record"; ``germquant register --prior stage``).

``read_stage_record(path)`` returns, over the frames of position P = 0 (the position ``read_stack`` reads):

* ``timepoints``, one dict per timepoint in t order: ``t_index``, ``t_start_s`` (the first frame's
  ``relativeTimeMs``, the time stamp ``read_nd2_metadata`` reports as ``t_seconds``), ``stage_x_um`` and
  ``stage_y_um`` (the medians over the stack's frames), ``z_first_um`` and ``z_last_um`` (the Z Coord of the stack's
  first and last frame), ``zdrive_um`` (the median Ti2 ZDrive, None when the file has no such tag) and ``n_frames``;
* ``stack``: ``n_planes``, ``z_step_um``, ``z_device``, ``bottom_to_top`` and ``home_index`` of the ZStackLoop;
* ``events``: ``time_s``, ``description`` and ``after_t_index``, the last timepoint whose first frame started at or
  before the event (None for an event before the first timepoint; the rule of the Kate work's read_nd2_meta.py),
  in time order; ``events_error`` says why there are none when the event list could not be read;
* ``source``: the path, size, mtime, nd2 package version, file version and this reader's version.

Why a private call. nd2's public API (0.11.3, pinned in constraints/) gives no per-frame stage position:
``ND2File.events()`` zips the per-frame columns down to the shortest custom tag (a one-value exposure tag on the
lab's files, so it returns one frame row plus the events), and ``frame_metadata(i).channels[0].position`` is the same
for every frame. The per-frame X, Y, Z and ZDrive columns come from the reader's ``_custom_tags()``, the one private
call here. It is guarded: StageRecordUnavailable (a ValueError) names the nd2 version when the method is missing, or
when the X, Y or Z column is absent or shorter than the frame count. Columns are found by their description prefix
("X Coord", "Y Coord", "Z Coord", "Ti2 ZDrive"), because the unit suffix decodes as mojibake on these files. The
event list and the frame time stamps come through the public ``events()`` and ``frame_metadata()``.
"""
from __future__ import annotations

import datetime as _dt
import os
from pathlib import Path

import numpy as np

READER_VERSION = "1"
COLUMN_PREFIXES = {"x": "X Coord", "y": "Y Coord", "z": "Z Coord", "zdrive": "Ti2 ZDrive"}
EVENT_KEY = "Events"          # the event description in nd2's events() records
TIME_KEY = "Time [s]"         # the time of an event in nd2's events() records


class StageRecordUnavailable(ValueError):
    """The nd2 file, or the installed nd2 package, does not give the per-frame stage record."""


def _nd2_version() -> str:
    try:
        import nd2

        return str(getattr(nd2, "__version__", "unknown"))
    except Exception:  # noqa: BLE001 - only for the message
        return "unknown"


def _column(tags: dict, prefix: str, n_frames: int, *, required: bool, path: Path) -> np.ndarray | None:
    """The first custom-tag column whose header starts with ``prefix``, as floats, cut to the frame count. None when
    the column is absent and not required; StageRecordUnavailable when a required one is absent or short."""
    for key, values in tags.items():
        if not str(key).startswith(prefix):
            continue
        try:
            arr = np.asarray(values, dtype=float).ravel()
        except (TypeError, ValueError):
            arr = None
        if arr is None or len(arr) < n_frames:
            if not required:
                return None
            n = 0 if arr is None else len(arr)
            raise StageRecordUnavailable(
                f"{path.name}: the {prefix!r} column of the per-frame record holds {n} value(s) for {n_frames} "
                f"frame(s) (nd2 {_nd2_version()}); no stage record can be read from this file")
        return arr[:n_frames]
    if required:
        raise StageRecordUnavailable(
            f"{path.name}: no {prefix!r} column in the per-frame record (nd2 {_nd2_version()}); no stage record can be "
            f"read from this file")
    return None


def _source(path: Path, f) -> dict:
    size = mtime = None
    try:
        st = os.stat(path)
        size = int(st.st_size)
        mtime = _dt.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds")
    except OSError:
        pass
    version = getattr(f, "version", None)
    try:
        file_version = [int(v) for v in version] if version is not None else None
    except (TypeError, ValueError):
        file_version = str(version)
    return {"path": str(path), "size_bytes": size, "mtime": mtime, "nd2_version": _nd2_version(),
            "file_version": file_version, "reader_version": READER_VERSION}


def _stack(f, frames_per_t: int) -> dict:
    loop = None
    try:
        loop = next((lp for lp in f.experiment if getattr(lp, "type", None) == "ZStackLoop"), None)
    except Exception:  # noqa: BLE001 - an unreadable experiment record: the stack geometry from the frame count
        loop = None
    if loop is None:
        return {"n_planes": int(frames_per_t), "z_step_um": None, "z_device": None, "bottom_to_top": None,
                "home_index": None}
    par = getattr(loop, "parameters", None)

    def _get(name, cast):
        v = getattr(par, name, None)
        try:
            return None if v is None else cast(v)
        except (TypeError, ValueError):
            return None

    return {"n_planes": int(getattr(loop, "count", frames_per_t) or frames_per_t), "z_step_um": _get("stepUm", float),
            "z_device": _get("deviceName", str), "bottom_to_top": _get("bottomToTop", bool),
            "home_index": _get("homeIndex", int)}


def read_stage_record(path: str | Path) -> dict:
    """The acquisition record of the nd2 time series at ``path`` (module docstring), metadata only.
    StageRecordUnavailable (a ValueError) when the per-frame stage columns cannot be read."""
    import nd2

    p = Path(path)
    with nd2.ND2File(str(p)) as f:
        rdr = getattr(f, "_rdr", None)
        tags_fn = getattr(rdr, "_custom_tags", None)
        if not callable(tags_fn):
            raise StageRecordUnavailable(
                f"{p.name}: the installed nd2 package ({_nd2_version()}) has no reader._custom_tags(), the per-frame "
                f"record the stage record is read from; write <stem>__stage_record.json another way or pin nd2 0.11.3")
        try:
            tags = dict(tags_fn() or {})
        except Exception as e:  # noqa: BLE001 - a file whose custom data cannot be decoded
            raise StageRecordUnavailable(f"{p.name}: the per-frame record could not be read (nd2 {_nd2_version()}): "
                                         f"{type(e).__name__}: {e}") from e
        loop = [dict(d) for d in f.loop_indices]
        n = len(loop)
        if n == 0:
            raise StageRecordUnavailable(f"{p.name}: no frames (nd2 {_nd2_version()})")
        x = _column(tags, COLUMN_PREFIXES["x"], n, required=True, path=p)
        y = _column(tags, COLUMN_PREFIXES["y"], n, required=True, path=p)
        z = _column(tags, COLUMN_PREFIXES["z"], n, required=True, path=p)
        zd = _column(tags, COLUMN_PREFIXES["zdrive"], n, required=False, path=p)
        frames: dict[int, list[int]] = {}
        for i, d in enumerate(loop):
            if int(d.get("P", 0) or 0) != 0:
                continue
            frames.setdefault(int(d.get("T", 0) or 0), []).append(i)
        timepoints = []
        for t in sorted(frames):
            fr = frames[t]
            try:
                t0 = float(f.frame_metadata(fr[0]).channels[0].time.relativeTimeMs) / 1000.0
            except Exception as e:  # noqa: BLE001 - without the time stamps no event can be placed
                raise StageRecordUnavailable(f"{p.name}: the time stamp of timepoint {t} is unreadable "
                                             f"(nd2 {_nd2_version()}): {type(e).__name__}: {e}") from e
            timepoints.append({"t_index": int(t), "t_start_s": t0, "stage_x_um": float(np.median(x[fr])),
                               "stage_y_um": float(np.median(y[fr])), "z_first_um": float(z[fr[0]]),
                               "z_last_um": float(z[fr[-1]]),
                               "zdrive_um": float(np.median(zd[fr])) if zd is not None else None,
                               "n_frames": len(fr)})
        events: list[dict] = []
        events_error = None
        try:
            for r in f.events(orient="records"):
                desc = r.get(EVENT_KEY) if isinstance(r, dict) else None
                if isinstance(desc, str) and desc.strip():
                    events.append({"time_s": float(r.get(TIME_KEY)), "description": desc})
        except Exception as e:  # noqa: BLE001 - the stage columns are still a record; say why the events are missing
            events, events_error = [], f"{type(e).__name__}: {e}"
        events.sort(key=lambda e: e["time_s"])
        starts = [(tp["t_index"], tp["t_start_s"]) for tp in timepoints]
        for e in events:
            after = None
            for t, s in starts:
                if s <= e["time_s"]:
                    after = t
            e["after_t_index"] = after
        stack = _stack(f, len(next(iter(frames.values()))) if frames else 1)
        source = _source(p, f)
    return {"source": source, "stack": stack, "n_frames": n, "timepoints": timepoints, "events": events,
            "events_error": events_error}
