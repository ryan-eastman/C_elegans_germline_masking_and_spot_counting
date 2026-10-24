"""Training rounds: training/<round>/ in a project (docs/PROJECTS.md, "Project models"; PLAN.md decision D29).

A round folder holds one round of the training loop (docs/TRAINING.md): its pairs, corrections.json and models
with their cards. Round names match ^[A-Za-z0-9_-]{1,40}$; a new round without a name takes the next free of r01,
r02 and so on. Folders under training/ with other names are not rounds and are left alone.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

from ..fsutil import long_path
from .manifest import ProjectError

ROUND_RE = re.compile(r"^[A-Za-z0-9_-]{1,40}$")
TRAINING_DIR = "training"


def rounds_in(root) -> list[str]:
    d = Path(root) / TRAINING_DIR
    if not d.is_dir():
        return []
    return sorted(p.name for p in d.iterdir() if p.is_dir() and ROUND_RE.fullmatch(p.name))


def check_round_name(name: str) -> str:
    from .images import windows_device_name

    name = str(name)
    if not ROUND_RE.fullmatch(name):
        raise ProjectError(f"round name {name!r} must match {ROUND_RE.pattern} (letters, digits, '_' and '-', at most "
                           f"40 characters)")
    if windows_device_name(name):
        raise ProjectError(f"round name {name!r} is a Windows device name, which no folder may have")
    return name


def next_round_name(root) -> str:
    d = Path(root) / TRAINING_DIR
    taken = {p.name.lower() for p in d.iterdir()} if d.is_dir() else set()
    n = 1
    while f"r{n:02d}" in taken:
        n += 1
    return f"r{n:02d}"


def make_round(root, name: str) -> Path:
    """Make training/<name>/. Refused when a folder of that name exists there, also one that differs only in letter
    case (R01 next to r01 would be one folder on Windows and macOS)."""
    path = Path(root) / TRAINING_DIR / check_round_name(name)
    os.makedirs(long_path(path.parent), exist_ok=True)
    same = [n for n in os.listdir(long_path(path.parent)) if n.lower() == name.lower()]
    if same and same != [name]:
        raise ProjectError(f"the training round {name} differs only in letter case from {same[0]} "
                           f"({path.parent / same[0]}), which Windows and macOS would make one folder")
    try:
        os.mkdir(long_path(path))
    except FileExistsError:
        raise ProjectError(f"the training round {name} already exists ({path})") from None
    return path
