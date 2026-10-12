"""Shared helpers for the system-pack tests (docs/PACKS.md, docs/DESIGN_system_packs.md, package P1).

``pack_env`` isolates every pack root and store a test could touch: no $GERMQUANT_PACKS, no
$GERMQUANT_CONFIG_ROOT, a private per-user folder, an empty model store ($GERMQUANT_MODELS) and an empty
Cellpose cache (CELLPOSE_LOCAL_MODELS_PATH), with the discovery cache cleared before and after.
"""
from __future__ import annotations

import hashlib
import os
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

from germquant import packs
from germquant.io.nd2_reader import Stack

ROOT = Path(__file__).resolve().parents[1]
FIXTURE_PACKS = ROOT / "tests" / "fixtures" / "packs"
TEST_PACK = FIXTURE_PACKS / "test_pack"
CELEGANS_PACK = ROOT / "packs" / "celegans_germline"
TINYNET = "tinynet@f5055f8a"
TINYNET_BYTES = b"germquant test pack: tinynet stand-in weights, not a real model\n"
TINYNET_SHA = hashlib.sha256(TINYNET_BYTES).hexdigest()
PYTHON = sys.executable
HEAVY = ("torch", "cupy", "spotmax", "cellpose")


@pytest.fixture
def pack_env(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    for var in ("GERMQUANT_PACKS", "GERMQUANT_CONFIG_ROOT", "GERMQUANT_DEVICE"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("LOCALAPPDATA", str(home / "localappdata"))
    monkeypatch.setenv("XDG_DATA_HOME", str(home / "xdg"))
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    models = tmp_path / "models_root"
    models.mkdir()
    monkeypatch.setenv("GERMQUANT_MODELS", str(models))
    cache = tmp_path / "cellpose_cache"
    cache.mkdir()
    monkeypatch.setenv("CELLPOSE_LOCAL_MODELS_PATH", str(cache))
    packs.refresh()
    yield SimpleNamespace(tmp=tmp_path, models=models, cache=cache, home=home, blobs=models / "models" / "blobs")
    packs.refresh()


def use_fixture_packs(monkeypatch) -> None:
    """$GERMQUANT_PACKS = tests/fixtures/packs, where the test pack sits as a pack under development."""
    monkeypatch.setenv("GERMQUANT_PACKS", str(FIXTURE_PACKS))
    packs.refresh()


def copy_pack(src: Path, dest: Path, *, edit=None, rebuild: bool = True) -> Path:
    """Copy a pack folder, apply ``edit(manifest_dict) -> None`` to its pack.yaml, and rebuild SHA256SUMS."""
    shutil.copytree(src, dest, ignore=shutil.ignore_patterns("__pycache__", ".DS_Store"))
    if edit is not None:
        data = yaml.safe_load((dest / "pack.yaml").read_text(encoding="ascii"))
        edit(data)
        (dest / "pack.yaml").write_bytes(yaml.safe_dump(data, sort_keys=False, allow_unicode=False).encode("ascii"))
    if rebuild:
        packs.write_sums(dest)
    return dest


def subprocess_env(**extra) -> dict:
    env = {k: v for k, v in os.environ.items()}
    env["PYTHONPATH"] = str(ROOT / "src") + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env.setdefault("CUDA_VISIBLE_DEVICES", "-1")
    env["PYTHONIOENCODING"] = "utf-8"
    env.update({k: str(v) for k, v in extra.items()})
    return env


def small_stack(names=("DAPI", "foci")) -> Stack:
    """A small two-channel field like the test pack's example (five nuclei with two puncta each)."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("make_test_pack", ROOT / "scripts" / "make_test_pack.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    data = mod.example_pixels().astype(np.float32)
    return Stack(data=data, spacing=mod.SPACING, channel_names=list(names), path=Path("x.nd2"), spacing_ok=True)


def patch_nd2_reader(monkeypatch, stack: Stack) -> None:
    """Feed ``stack`` to process_image through the .nd2 reader seam (as test_end_to_end does)."""
    from germquant import pipeline

    c, z, y, x = stack.data.shape
    meta = {"sizes": {"C": c, "Z": z, "Y": y, "X": x}, "spacing": stack.spacing, "spacing_ok": True,
            "channel_names": list(stack.channel_names), "dtype": "float32", "is_2d": False}
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    monkeypatch.setattr(pipeline, "read_stack", lambda *a, **k: stack)
