"""Run logs of project runs, runs/<run_id>.jsonl (docs/PROJECTS.md, "Run logs"; design 7.6; PLAN.md D17, D48).

Each run writes its own log, made with an exclusive create, so no two processes ever write one file: shards and
Snakemake jobs run at once, often on SMB or NFS shares where appends from several hosts are not atomic. Every
line is one JSON object in ASCII (json.dumps with ensure_ascii), written, flushed and fsynced at once, so a crash
never leaves half a line behind a complete one:

    {"event", "run_id", "time" (UTC, ISO 8601 with milliseconds), ...the event's fields}

The engine owns the format (D17). Its events and fields are those the plugin writes for its own runs
(napari_germquant/_runner.py, RunLog and execute_run), plus project and shard in run_start and interrupted in
run_end:

* run_start: RUN_START_FIELDS;
* image_done, image_failed and image_skipped: IMAGE_FIELDS (None where the event has no value, as a failed image
  has no marker);
* run_end: counts {done, failed, skipped} and interrupted.

done_sha256 is the sha256 of the image's done marker with CRLF normalised to LF (D48): the marker is written in
text mode, so its raw bytes differ between Windows and the other systems, its normalised digest does not.

Readers (project_status, the plugin, collect) merge the logs by time (read_logs). A run id is the plugin's form,
``YYYYMMDDTHHMMSSZ_<8 hex>`` (new_run_id), and names the log.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import uuid
from datetime import datetime, timezone
from pathlib import Path

from ..fsutil import long_path
from .manifest import ProjectError

RUN_ID_RE = re.compile(r"^\d{8}T\d{6}Z_[0-9a-f]{8}$")
VIA = ("gui-subprocess", "gui-inprocess", "cli", "snakemake")
RUNS_DIR = "runs"
RUN_START_FIELDS = ("via", "host", "os", "interpreter", "interpreter_source", "engine", "python", "constraints",
                    "device_requested", "device_selected", "config_hash", "config_hash_source", "recipe", "xy_stride",
                    "disabled_stages", "argv", "images", "plugin", "log_location", "project", "shard")
IMAGE_FIELDS = ("image_id", "out_dir", "status", "done_sha256", "config_hash", "stage_hashes", "failed_stages",
                "qc_pass", "segment_device_used", "elapsed_s", "error")
IMAGE_EVENTS = ("image_done", "image_failed", "image_skipped")
RUN_END_FIELDS = ("counts", "interrupted")


def new_run_id() -> str:
    """A run id of the plugin's form: the UTC start time and 8 random hex digits."""
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]


def check_run_id(run_id) -> str:
    """``run_id`` when it has the form YYYYMMDDTHHMMSSZ_<8 hex> (it names a file), else ProjectError."""
    text = str(run_id)
    if not RUN_ID_RE.fullmatch(text):
        raise ProjectError(f"run id {text!r} is not of the form YYYYMMDDTHHMMSSZ_<8 hex digits> (for example "
                           f"{new_run_id()}), which names the run's log runs/<run id>.jsonl")
    return text


def runs_dir(root) -> Path:
    return Path(root) / RUNS_DIR


def log_path(root, run_id: str) -> Path:
    return runs_dir(root) / f"{run_id}.jsonl"


def utc_time() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


def jsonable(value):
    """``value`` with paths, tuples, numpy scalars and other objects made JSON values (str as the fallback)."""
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [jsonable(v) for v in value]
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    item = getattr(value, "item", None)                     # a numpy scalar
    if callable(item):
        try:
            return jsonable(item())
        except (TypeError, ValueError):
            pass
    return str(value)


class RunLog:
    """One run's log, runs/<run_id>.jsonl. ``start`` makes the file with an exclusive create (ProjectError when a
    log of that run id exists) and writes run_start; ``event`` appends one fsynced line. Both return the record as
    it reads back from the file."""

    def __init__(self, root, run_id: str) -> None:
        self.run_id = check_run_id(run_id)
        self.path = log_path(root, self.run_id)

    def _line(self, kind: str, fields: dict) -> tuple[bytes, dict]:
        rec = {"event": kind, "run_id": self.run_id, "time": utc_time(), **jsonable(fields)}
        text = json.dumps(rec, ensure_ascii=True)
        return (text + "\n").encode("ascii"), json.loads(text)

    def start(self, **fields) -> dict:
        os.makedirs(long_path(self.path.parent), exist_ok=True)
        data, rec = self._line("run_start", fields)
        try:
            fd = os.open(long_path(self.path), os.O_CREAT | os.O_EXCL | os.O_WRONLY | getattr(os, "O_BINARY", 0),
                         0o644)
        except FileExistsError:
            raise ProjectError(f"a run log for run id {self.run_id} already exists ({self.path}): every run needs an "
                               f"id of its own") from None
        try:
            os.write(fd, data)
            os.fsync(fd)
        finally:
            os.close(fd)
        return rec

    def event(self, kind: str, **fields) -> dict:
        data, rec = self._line(kind, fields)
        with open(long_path(self.path), "ab") as fh:
            fh.write(data)
            fh.flush()
            os.fsync(fh.fileno())
        return rec


def done_sha256(path) -> str | None:
    """The sha256 of a done marker's bytes with CRLF normalised to LF (D48); None when it cannot be read."""
    try:
        with open(long_path(path), "rb") as fh:
            raw = fh.read()
    except OSError:
        return None
    return hashlib.sha256(raw.replace(b"\r\n", b"\n")).hexdigest()


def read_log(path) -> list[dict]:
    """The events of one log: its complete lines that are JSON objects, in file order. A partial last line (a run
    still writing) and a line that is not a JSON object are left out; a missing file has no events."""
    try:
        with open(long_path(path), "rb") as fh:
            raw = fh.read()
    except FileNotFoundError:
        return []
    lines = raw.split(b"\n")
    out = []
    for line in lines[:-1]:                                  # the last piece has no newline yet: not complete
        try:
            rec = json.loads(line.decode("utf-8"))
        except ValueError:
            continue
        if isinstance(rec, dict) and rec.get("event"):
            out.append(rec)
    return out


def log_files(root) -> list[Path]:
    """The run logs of a project: runs/<run_id>.jsonl, sorted by name (a run id starts with its UTC start time)."""
    folder = runs_dir(root)
    try:
        names = os.listdir(long_path(folder))
    except (FileNotFoundError, NotADirectoryError):
        return []
    return [folder / n for n in sorted(names) if n.endswith(".jsonl") and RUN_ID_RE.fullmatch(n[: -len(".jsonl")])]


def read_logs(root) -> list[dict]:
    """Every event of every run log of the project, merged by time (events with one time keep the order of their
    files and lines)."""
    events = [e for p in log_files(root) for e in read_log(p)]
    return sorted(events, key=lambda e: str(e.get("time") or ""))
