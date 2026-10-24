"""The project lock, germquant_project.lock (docs/PROJECTS.md, "The lock"; design 7.2; PLAN.md decision D18).

Every interactive write (new, add, edit, use-model, rounds --new, and the API calls behind them) holds the lock
for its read-modify-write of the manifest and the recipe files. The lock is a file made with an exclusive
create (os.open with O_CREAT | O_EXCL), which behaves the same on Windows, macOS and Linux and on SMB and NFS
shares, holding JSON ``{host, pid, time, command}``. A second writer is refused with a message that names the
holder; the writer that made the lock removes it when it is done, and only if it is still its own.

A lock older than 12 hours (STALE_AFTER_S) is reported as stale, by the refusal and by ``lock_report`` (what
`germquant project status` shows), but never broken silently: only ``germquant project unlock DIR --host H
--pid P`` (break_lock) removes it, and only when H and P are the holder's.
"""
from __future__ import annotations

import contextlib
import json
import os
import socket
from datetime import datetime, timezone
from pathlib import Path

from ..fsutil import long_path
from .manifest import LOCK_NAME, ProjectError, remove_file, utc_now

STALE_AFTER_S = 12 * 3600
STALE_AFTER_TEXT = "12 h"


class ProjectLockedError(ProjectError):
    """Another writer holds the project's lock; ``holder`` is its record (host, pid, time, command)."""

    def __init__(self, message: str, holder: dict):
        super().__init__(message)
        self.holder = dict(holder)


def lock_path(root) -> Path:
    return Path(root) / LOCK_NAME


def read_lock(root) -> dict | None:
    """The lock record ``{host, pid, time, command}``, None when the project is not locked. A lock file that cannot
    be read as JSON gives the same keys as None plus ``unreadable`` (the reason)."""
    p = lock_path(root)
    try:
        with open(long_path(p), "rb") as fh:
            raw = fh.read()
    except FileNotFoundError:
        return None
    except OSError as e:
        return {"host": None, "pid": None, "time": None, "command": None, "unreadable": f"{type(e).__name__}: {e}"}
    try:
        data = json.loads(raw.decode("utf-8"))
        if not isinstance(data, dict):
            raise ValueError("not a JSON object")
    except ValueError as e:
        return {"host": None, "pid": None, "time": None, "command": None,
                "unreadable": f"not a lock record ({e}); a writer may have stopped while it was making the lock"}
    return {"host": data.get("host"), "pid": data.get("pid"), "time": data.get("time"), "command": data.get("command")}


def _age_s(root, holder: dict) -> float | None:
    """Seconds since the lock was taken: from its time field, else from the lock file's modification time."""
    try:
        t = datetime.fromisoformat(str(holder.get("time")))
        if t.tzinfo is None:
            t = t.replace(tzinfo=timezone.utc)
        return (datetime.now(timezone.utc) - t).total_seconds()
    except (TypeError, ValueError):
        pass
    try:
        return datetime.now(timezone.utc).timestamp() - os.stat(long_path(lock_path(root))).st_mtime
    except OSError:
        return None


def is_stale(root, holder: dict) -> bool:
    age = _age_s(root, holder)
    return age is not None and age > STALE_AFTER_S


def lock_report(root) -> dict:
    """``{"held", "host", "pid", "since", "stale"}``: the lock block of the status record (PLAN.md 9.3)."""
    holder = read_lock(root)
    if holder is None:
        return {"held": False, "host": None, "pid": None, "since": None, "stale": False}
    return {"held": True, "host": holder.get("host"), "pid": holder.get("pid"), "since": holder.get("time"),
            "stale": is_stale(root, holder)}


def describe(holder: dict) -> str:
    if holder.get("unreadable"):
        return f"an unreadable lock file ({holder['unreadable']})"
    cmd = f" ({holder['command']})" if holder.get("command") else ""
    return f"{holder.get('host')} pid {holder.get('pid')} since {holder.get('time')}{cmd}"


def unlock_hint(root, holder: dict) -> str:
    if holder.get("unreadable"):
        return (f"if no germquant command is writing this project, delete {lock_path(root)} by hand")
    return (f"if that process is gone, remove the lock with `germquant project unlock \"{root}\" --host "
            f"{holder.get('host')} --pid {holder.get('pid')}`")


def _refusal(root, holder: dict) -> str:
    msg = (f"the project {root} is locked by {describe(holder)}: another germquant command is writing it; try again "
           f"when it has finished.")
    if holder.get("unreadable") or is_stale(root, holder):
        msg += f" The lock is older than {STALE_AFTER_TEXT}: {unlock_hint(root, holder)}."
    return msg


@contextlib.contextmanager
def project_lock(root, command: str):
    """Hold the project's lock for the duration of the block; ProjectLockedError when another writer holds it."""
    root = Path(root)
    p = lock_path(root)
    me = {"host": socket.gethostname(), "pid": os.getpid(), "time": utc_now(), "command": str(command)}
    try:
        fd = os.open(long_path(p), os.O_CREAT | os.O_EXCL | os.O_WRONLY | getattr(os, "O_BINARY", 0), 0o644)
    except FileExistsError:
        holder = read_lock(root) or {"host": None, "pid": None, "time": None, "command": None}
        raise ProjectLockedError(_refusal(root, holder), holder) from None
    try:
        os.write(fd, json.dumps(me).encode("ascii"))
        os.fsync(fd)
    finally:
        os.close(fd)
    try:
        yield me
    finally:
        cur = read_lock(root)
        if cur and all(cur.get(k) == me[k] for k in ("host", "pid", "time")):
            remove_file(p)


def break_lock(root, *, host, pid) -> bool:
    """Remove the lock when its holder is ``host`` and ``pid`` (True); False when the project is not locked;
    ProjectError naming the holder when it is someone else, or when the lock file cannot be read."""
    holder = read_lock(root)
    if holder is None:
        return False
    if holder.get("unreadable"):
        raise ProjectError(f"the lock of {root} is {describe(holder)}; {unlock_hint(root, holder)}")
    same_host = str(holder.get("host")) == str(host)
    try:
        same_pid = int(holder.get("pid")) == int(pid)
    except (TypeError, ValueError):
        same_pid = False
    if not (same_host and same_pid):
        raise ProjectError(f"the lock of {root} is held by {describe(holder)}, not by {host} pid {pid}: it was not "
                           f"removed. Give the holder's host and pid (germquant project status shows them)")
    remove_file(lock_path(root))
    return True
