"""Provenance manifest — every run writes git SHA + config hash + tool versions + voxel
size into the output folder, so every CSV row traces back to exact code + parameters.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

from .fsutil import long_path

_TOOLS = ["germquant", "nd2", "numpy", "scipy", "scikit-image", "pandas", "cellpose", "spotmax", "torch",
          "cellacdc", "cupy", "skan", "tifffile", "snakemake", "bioio"]

_FILE_SHA_CACHE: dict[tuple[str, float, int], str] = {}


def file_sha256(path: str | Path | None) -> str:
    """sha256 of a file (the Cellpose model), cached on (path, mtime, size) so a batch hashes a 1 GB
    model once. 'absent' when the path is missing or a directory / a bare model name such as 'cpsam'."""
    if not path:
        return "absent"
    p = Path(path)
    if not p.is_file():
        return "absent"
    st = p.stat()
    key = (str(p.resolve()), st.st_mtime, st.st_size)
    if key not in _FILE_SHA_CACHE:
        h = hashlib.sha256()
        with open(long_path(p), "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 24), b""):
                h.update(chunk)
        _FILE_SHA_CACHE[key] = h.hexdigest()
    return _FILE_SHA_CACHE[key]


def _repo_root() -> Path | None:
    """Walk up from this source file looking for the repo root (a .git dir or pixi.lock)."""
    here = Path(__file__).resolve()
    for d in (here, *here.parents):
        if (d / ".git").exists() or (d / "pixi.lock").exists() or (d / "pyproject.toml").exists():
            return d
    return None


def _repo_dir(repo: str | Path | None) -> str:
    """Default to the germquant source repo, not the process cwd."""
    if repo is not None:
        return str(repo)
    root = _repo_root()
    return str(root) if root is not None else str(Path(__file__).resolve().parent)


def lockfile_sha256() -> str:
    """sha256 of pixi.lock if present — the ARCHITECTURE §4 reproducibility anchor.

    Returns 'absent' when no lockfile is committed (e.g. a pip/uv install), so the manifest
    states honestly whether the run was pinned to a resolved environment.
    """
    root = _repo_root()
    if root is None:
        return "absent"
    lock = root / "pixi.lock"
    if not lock.exists():
        return "absent"
    return hashlib.sha256(lock.read_bytes()).hexdigest()


def git_sha(repo: str | Path | None = None) -> str:
    cwd = _repo_dir(repo)
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=cwd,
            capture_output=True, text=True, check=True,
        )
        dirty = subprocess.run(
            ["git", "status", "--porcelain"], cwd=cwd, capture_output=True, text=True,
        ).stdout.strip()
        return out.stdout.strip() + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def tool_versions() -> dict[str, str]:
    versions = {}
    for name in _TOOLS:
        try:
            versions[name] = metadata.version(name)
        except Exception:
            versions[name] = "not-installed"
    return versions


def run_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_manifest(out_dir: str | Path, *, config_hash: str, config: dict, extra: dict | None = None) -> dict:
    out_dir = Path(out_dir)
    os.makedirs(long_path(out_dir), exist_ok=True)
    manifest = {
        "pipeline_version": _safe_version("germquant"),
        "git_sha": git_sha(),
        "pixi_lock_sha256": lockfile_sha256(),
        "config_hash": config_hash,
        "run_timestamp": run_timestamp(),
        "tool_versions": tool_versions(),
        "config": config,
    }
    if extra:
        manifest.update(extra)
    with open(long_path(out_dir / "run_manifest.json"), "w", encoding="utf-8") as fh:
        fh.write(json.dumps(manifest, indent=2, default=str))
    return manifest


def write_done_marker(out_dir: str | Path, image_id: str, *, config_hash: str, enabled_stages: list[str],
                      stage_hashes: dict[str, str], xy_stride: int = 1, z_range=None,
                      failed_stages: list[str] | None = None, extra: dict | None = None) -> Path:
    """``<image_id>__done.json``: written as the very last act of a successful process_image, so its
    presence (with a matching config_hash, enabled stage set and run geometry) means every output of
    that run is on disk. Batch resume keys on it; a folder without it is reprocessed, and so is one whose
    marker lists failed (non-fatal) stages. ``extra`` adds fields (a time-lapse timepoint's t_index and
    t_seconds) before ``finished``; with none given the marker is the one it always was."""
    out_dir = Path(out_dir)
    p = out_dir / f"{image_id}__done.json"
    payload = {"image_id": image_id, "config_hash": config_hash, "enabled_stages": enabled_stages,
               "stage_hashes": stage_hashes, "xy_stride": int(xy_stride),
               "z_range": list(z_range) if z_range else None,
               "failed_stages": list(failed_stages or [])}
    if extra:
        payload.update(extra)
    payload["finished"] = run_timestamp()
    with open(long_path(p), "w", encoding="utf-8") as fh:
        fh.write(json.dumps(payload, indent=1))
    return p


def clear_done_marker(out_dir: str | Path, image_id: str) -> None:
    """Remove a stale marker before a rerun starts writing, so an interrupted rerun can never be
    mistaken for a complete one."""
    p = Path(out_dir) / f"{image_id}__done.json"
    try:
        os.remove(long_path(p))
    except FileNotFoundError:
        pass


def read_done_marker(out_dir: str | Path, image_id: str) -> dict | None:
    p = Path(out_dir) / f"{image_id}__done.json"
    try:
        with open(long_path(p), encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def _safe_version(name: str) -> str:
    try:
        return metadata.version(name)
    except Exception:
        return "0.1.0"
