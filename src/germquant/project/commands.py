"""The commands that run a project headless (docs/PROJECTS.md, "Export commands"; design 7.8 and 9).

``export_commands(project, *, shell="sh", shards=None)`` returns a script that runs the project as `germquant
project run` does and then collects its tables: for a POSIX shell (sh: bash, zsh, a SLURM job) or for Windows
PowerShell (ps1). With ``shards`` N it holds one run line per shard, ``--shard K/N``, which can run one after the
other, side by side, or one per computer or SLURM array task (each writes only its own run log and its images'
result folders, never the manifest); the header says how to run them as a SLURM array.

The engine is the interpreter the project records for its recipe (run.interpreters, written by `germquant
project edit --interpreter RECIPE=PYTHON`), else the interpreter writing the commands, run as ``python -m
germquant``, so the commands run the engine the project was run with, whatever is first on PATH. Everything else
a run uses (the recipe, the stage switches, the device, stride and z range) is in the project folder, so the run
lines name only the folder. Paths are quoted for the shell. The text is not hashed and carries no time stamp.
"""
from __future__ import annotations

import shlex

from .manifest import ProjectError

SHELLS = ("sh", "ps1")


def _ps_quote(text: str) -> str:
    return "'" + str(text).replace("'", "''") + "'"


def _engine_python(manifest: dict) -> str:
    import sys

    rec = (manifest["run"].get("interpreters") or {}).get(manifest["pack"]["recipe"])
    if isinstance(rec, dict) and rec.get("python"):
        return str(rec["python"])
    return sys.executable


def export_commands(project, *, shell="sh", shards=None) -> str:
    """The commands that run the project headless (module docstring), as the text of a script."""
    from .. import __version__
    from ..packs import use_project
    from ..provenance import git_sha
    from . import Project, _root, project_config
    from .manifest import read_manifest

    if shell not in SHELLS:
        raise ProjectError(f"shell {shell!r} is not one of {', '.join(SHELLS)}")
    if shards is not None:
        if isinstance(shards, bool) or not isinstance(shards, int) or shards < 1:
            raise ProjectError(f"shards {shards!r} must be a whole number of 1 or more")
    root = _root(project)
    m = read_manifest(root)
    with use_project(root):
        cfg = project_config(Project(root=root, manifest=m))
    python = _engine_python(m)
    p = m["pack"]
    n = int(shards) if shards else 1
    edited = " (edited)" if m["recipe"].get("file") else ""
    head = [f"germquant project {m.get('name')}: {root}",
            f"pack {p['id']}@{p['version']} recipe {p['recipe']}{edited}, config_hash {cfg.hash}",
            f"written by germquant {__version__} (git {git_sha()}); the engine: {python}"]
    if shell == "sh":
        exe, q = f"{shlex.quote(python)} -m germquant", shlex.quote(str(root))
        lines = ["#!/bin/sh"] + [f"# {h}" for h in head]
    else:
        exe, q = f"& {_ps_quote(python)} -m germquant", _ps_quote(str(root))
        lines = [f"# {h}" for h in head]
    lines.append("# a run skips the images already done; `germquant project status` shows each image afterwards")
    if n > 1:
        lines.append(f"# {n} shards: run these lines one after the other, side by side, or one per computer, then the "
                     f"collect line once all of them have finished")
        if shell == "sh":
            lines.append(f"# as a SLURM array (sbatch --array=1-{n}) each task runs: {exe} project run {q} --shard "
                         f"${{SLURM_ARRAY_TASK_ID}}/{n}")
    lines += ([f"{exe} project run {q}"] if n == 1 else
              [f"{exe} project run {q} --shard {k}/{n}" for k in range(1, n + 1)])
    lines.append(f"{exe} project collect {q}")
    return "\n".join(lines) + "\n"
