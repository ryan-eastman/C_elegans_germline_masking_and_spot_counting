"""Provenance manifest — every run writes git SHA + config hash + tool versions + voxel
size into the output folder, so every CSV row traces back to exact code + parameters.
"""
from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

_TOOLS = ["germquant", "nd2", "numpy", "scipy", "scikit-image", "pandas", "skan", "cellpose", "torch"]


def _repo_dir(repo: str | Path | None) -> str:
    """Default to the germquant source repo, not the process cwd."""
    if repo is not None:
        return str(repo)
    return str(Path(__file__).resolve().parent)


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
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "pipeline_version": _safe_version("germquant"),
        "git_sha": git_sha(),
        "config_hash": config_hash,
        "run_timestamp": run_timestamp(),
        "tool_versions": tool_versions(),
        "config": config,
    }
    if extra:
        manifest.update(extra)
    (out_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    return manifest


def _safe_version(name: str) -> str:
    try:
        return metadata.version(name)
    except Exception:
        return "0.1.0"
