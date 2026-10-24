"""Shared helpers for the post-stage tests (docs/POST_STAGES.md; plan D35).

* ``fixture_post_stage``: the fixture post-stage of germquant.post_stages.testing registered in this process for
  one test, and removed afterwards;
* ``make_result_folder``: a hand-made result folder in the shape process_image leaves behind (tables, a label
  image, the stage record, the run manifest and the completion marker), made without running the pipeline;
* ``snapshot``: every file of a folder with its bytes, to show that nothing that was there changed;
* ``heavy_after`` and ``run_child``: a fresh interpreter with this checkout's src (and any extra folder, such as
  an on-disk fixture distribution) on PYTHONPATH. A sys.modules assertion runs there, never in the suite's own
  process.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON = sys.executable
HEAVY = ("torch", "cupy", "spotmax", "cellpose")
# reading an image header through bioio imports dask.array (and dask imports CuPy where CuPy is installed)
READERS = ("bioio", "dask.array")
GENERIC = ROOT / "config" / "profiles" / "generic.yaml"
SPACING = (0.5, 0.2, 0.2)
READ_SUB_HASH = "0123456789ab"

_PROBE = ("import sys, atexit\n"
          "atexit.register(lambda: sys.stderr.write('HEAVY=' + ','.join(m for m in {heavy!r} if m in sys.modules)"
          " + '\\n'))\n")


@pytest.fixture
def fixture_post_stage():
    """The fixture post-stage registered in process (germquant.post_stages.testing.register_fixture)."""
    from germquant import post_stages as PS
    from germquant.post_stages import testing as T

    name = T.register_fixture(replace=True)
    try:
        yield name
    finally:
        PS.unregister_post_stage(name)


def child_env(pythonpath=(), **extra) -> dict:
    """The environment of a child interpreter: this checkout's src first, then ``pythonpath``, then whatever
    PYTHONPATH the suite already had."""
    env = dict(os.environ)
    parts = [str(ROOT / "src"), *[str(p) for p in pythonpath]]
    if env.get("PYTHONPATH"):
        parts.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(parts)
    env.setdefault("CUDA_VISIBLE_DEVICES", "-1")
    env["PYTHONIOENCODING"] = "utf-8"
    env.update({k: str(v) for k, v in extra.items()})
    return env


def run_child(args: list[str], *, pythonpath=(), timeout: int = 300, **env) -> subprocess.CompletedProcess:
    """``python <args>`` in a fresh interpreter (``args`` may start with -c or -m)."""
    return subprocess.run([PYTHON, *[str(a) for a in args]], capture_output=True, text=True, timeout=timeout,
                          env=child_env(pythonpath, **env), cwd=ROOT)


def heavy_after(code: str, *, heavy=HEAVY, pythonpath=(), **env) -> tuple[str, subprocess.CompletedProcess]:
    """Run ``code`` in a fresh interpreter; returns (the comma-separated modules of ``heavy`` loaded at exit,
    the completed process)."""
    r = run_child(["-c", _PROBE.format(heavy=tuple(heavy)) + code], pythonpath=pythonpath, **env)
    m = re.search(r"HEAVY=(.*)", r.stderr)
    assert m, f"no probe line; rc={r.returncode} stderr={r.stderr[-3000:]}"
    return m.group(1).strip(), r


def snapshot(folder: Path) -> dict[str, bytes]:
    """{file name: bytes} for every file directly in ``folder``."""
    return {p.name: p.read_bytes() for p in sorted(Path(folder).iterdir()) if p.is_file()}


def make_stack(nz: int = 4, ny: int = 24, nx: int = 32):
    """A two-channel stack (DAPI, foci) with three box nuclei, and its (Z, Y, X) labels."""
    from germquant.io.nd2_reader import Stack

    labels = np.zeros((nz, ny, nx), dtype=np.int32)
    labels[:, 2:10, 2:10] = 1
    labels[1:3, 12:20, 4:14] = 2
    labels[:, 4:12, 18:30] = 3
    rng = np.random.default_rng(7)
    dna = (labels > 0).astype(np.float32) * 500.0 + rng.uniform(0, 5, labels.shape).astype(np.float32)
    foci = rng.uniform(0, 100, labels.shape).astype(np.float32)
    data = np.stack([dna, foci]).astype(np.float32)
    stack = Stack(data=data, spacing=SPACING, channel_names=["DAPI", "foci"], path=Path("img.nd2"), spacing_ok=True)
    return stack, labels


def make_result_folder(root: Path, image_id: str = "img_001", *, xy_stride: int = 1, z_range=None,
                       stack=None, labels=None) -> SimpleNamespace:
    """A hand-made result folder: nuclei, spots, image_summary and an extra table (coloc) the fixture does not
    read, the nuclei label image, the stage record, the run manifest, a montage stand-in and the completion
    marker, written with the generic profile's config_hash and enabled stages (batch.is_done is true)."""
    import tifffile

    from germquant import provenance
    from germquant.batch import enabled_stage_names
    from germquant.config import load_config

    cfg = load_config(GENERIC)
    if stack is None:
        stack, labels = make_stack()
    out = Path(root) / image_id
    out.mkdir(parents=True, exist_ok=True)
    roles = {"dna": 0, "foci": 1, "granule": None}
    shared = {"image_id": image_id, "file_path": str(Path(root) / f"{image_id}.nd2"),
              "voxel_dz_um": SPACING[0], "voxel_dy_um": SPACING[1], "voxel_dx_um": SPACING[2],
              "n_channels": 2, "dapi_present": True, "channel_map_json": json.dumps(roles),
              "config_hash": cfg.hash, "run_timestamp": "2026-10-23T00:00:00+00:00"}
    ids = sorted(int(v) for v in np.unique(labels) if v)
    nuclei = pd.DataFrame({"nucleus_id": ids, "n_voxels": [int((labels == i).sum()) for i in ids]})
    spots = pd.DataFrame({"spot_id": [1, 2, 3], "nucleus_id": [1, 1, 3], "marker": ["foci"] * 3,
                          "z_um": [0.5, 1.0, 1.5], "y_um": [1.0, 1.2, 1.4], "x_um": [1.0, 1.1, 4.0]})
    summary = pd.DataFrame([{"n_nuclei": len(ids), "qc_pass": True, "qc_flags": ""}])
    coloc = pd.DataFrame({"granule_id": [1], "overlap": [0.25]})
    for name, df in (("nuclei", nuclei), ("spots", spots), ("image_summary", summary), ("coloc", coloc)):
        df = df.copy()
        for k, v in shared.items():
            df[k] = v
        df.to_csv(out / f"{image_id}__{name}.csv", index=False)
    tifffile.imwrite(out / f"{image_id}__nuclei_labels.tif", labels.astype(np.int32), compression="zlib")
    (out / f"{image_id}__montage.png").write_bytes(b"\x89PNG stand-in, not an image\n")
    (out / "run_manifest.json").write_text(json.dumps({"config_hash": cfg.hash, "xy_stride": xy_stride}, indent=2),
                                           encoding="utf-8")
    (out / f"{image_id}__stages.json").write_text(json.dumps({"image_id": image_id, "complete": True}, indent=1),
                                                  encoding="utf-8")
    hashes = {"read": READ_SUB_HASH, "segment": "aaaaaaaaaaaa", "spots": "bbbbbbbbbbbb"}
    provenance.write_done_marker(out, image_id, config_hash=cfg.hash, enabled_stages=enabled_stage_names(cfg),
                                 stage_hashes=hashes, xy_stride=xy_stride, z_range=z_range)
    image_path = Path(root) / f"{image_id}.nd2"
    image_path.write_bytes(b"not a real nd2 file: the tests hand the stack over or patch the reader\n")
    return SimpleNamespace(out_dir=out, image_id=image_id, cfg=cfg, stack=stack, labels=labels,
                           image_path=image_path, roles=roles)


def edit_marker(out_dir: Path, image_id: str, **changes) -> None:
    """Rewrite fields of the completion marker in place (a test stand-in for another main run)."""
    p = Path(out_dir) / f"{image_id}__done.json"
    m = json.loads(p.read_text(encoding="utf-8"))
    for k, v in changes.items():
        if isinstance(v, dict) and isinstance(m.get(k), dict):
            m[k].update(v)
        else:
            m[k] = v
    p.write_text(json.dumps(m, indent=1), encoding="utf-8")


def post_record(out_dir: Path, image_id: str) -> dict:
    return json.loads((Path(out_dir) / f"{image_id}__post.json").read_text(encoding="ascii"))
