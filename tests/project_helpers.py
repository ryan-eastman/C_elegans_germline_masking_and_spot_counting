"""Shared helpers for the project tests (docs/PROJECTS.md; tranche 5 Train B unit P5-PRJ).

* ``project_env``: pack_env (every pack root and store isolated) plus $GERMQUANT_PACKS on tests/fixtures/packs,
  where the test pack sits as a pack under development.
* ``new_test_project`` makes a project on test_pack 0.1.0 (recipe basic unless given).
* ``write_ome_tiff`` writes a small OME-TIFF with tifffile: channel names, a z step or none (a calibrated YX
  image), one or more OME images.
* ``fake_nd2_headers`` replaces the .nd2 header seam of germquant.project.headers, so a test can add an .nd2
  stand-in with any sizes (P > 1 included) without a real Nikon file.
* ``snapshot`` maps every file under a folder to its bytes, for "writes nothing" checks.
* ``heavy_after(code)`` runs ``code`` in a fresh interpreter (sys.executable, this checkout's src first on
  PYTHONPATH) with the probe of tests/test_packs_isolation.py and returns the heavy modules it left in
  sys.modules: HEAVY plus READERS (torch, cupy, spotmax, cellpose, bioio, dask.array). A sys.modules check
  never runs in the suite's own process, which has loaded torch long before.
"""
from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import numpy as np
import pytest
from packs_helpers import FIXTURE_PACKS, HEAVY, PYTHON, ROOT, TEST_PACK, pack_env, subprocess_env  # noqa: F401
from packs_helpers import use_fixture_packs

READERS = ("bioio", "dask.array")
LIGHT_FORBIDDEN = HEAVY + READERS
EXAMPLE = TEST_PACK / "example" / "test_pack_synthetic.ome.tif"
TEST_PACK_REF = "test_pack@0.1.0"
BASIC_HASH = "02bb83cf74f9"              # test_pack 0.1.0 recipe basic, its expected_config_hash
GAMMA = chr(0x3B3)                       # a Greek gamma, as in gammaH2AX (built with chr: this source stays ASCII)

_PROBE = ("import sys, atexit\n"
          "atexit.register(lambda: sys.stderr.write('HEAVY=' + ','.join(m for m in {heavy!r} if m in sys.modules)"
          " + '\\n'))\n")


@pytest.fixture
def project_env(pack_env, monkeypatch):  # noqa: F811
    use_fixture_packs(monkeypatch)
    return pack_env


def new_test_project(folder: Path, recipe: str | None = None, **kw):
    from germquant.project import new_project

    return new_project(folder, TEST_PACK_REF, recipe, **kw)


def write_ome_tiff(path: Path, *, names=("DAPI", "foci"), zyx=(4, 16, 16), spacing=(0.5, 0.2, 0.2),
                   n_images: int = 1, seed: int = 0) -> Path:
    """An OME-TIFF with ``n_images`` OME images of (C, Z, Y, X) uint16 pixels and the OME channel names. A z of 0
    writes a calibrated YX image (axes CYX, no PhysicalSizeZ)."""
    import tifffile

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    z, y, x = zyx
    meta = {"Channel": {"Name": list(names)}, "PhysicalSizeY": float(spacing[1]), "PhysicalSizeX": float(spacing[2])}
    if z:
        shape, meta["axes"], meta["PhysicalSizeZ"] = (len(names), z, y, x), "CZYX", float(spacing[0])
    else:
        shape, meta["axes"] = (len(names), y, x), "CYX"
    with tifffile.TiffWriter(path, ome=True) as tw:
        for _ in range(n_images):
            tw.write(rng.integers(0, 4000, size=shape, dtype=np.uint16), metadata=dict(meta))
    return path


def fake_nd2_headers(monkeypatch, headers: dict[str, dict]) -> None:
    """Make germquant.project.headers read the .nd2 files named in ``headers`` (by file name) as headers with
    these sizes; ``{"sizes": {...}}`` entries may also set channel_names and spacing."""
    from germquant.project import headers as H

    def fake(path):
        h = headers[Path(path).name]
        sizes = dict(h["sizes"])
        return {"sizes": sizes, "spacing": tuple(h.get("spacing", (0.5, 0.1, 0.1))), "spacing_ok": True,
                "channel_names": list(h.get("channel_names", ["DAPI", "foci"])), "dtype": "uint16",
                "is_2d": "Z" not in sizes, "n_timepoints": int(sizes.get("T", 1)), "period_s": None,
                "t_seconds": None}

    monkeypatch.setattr(H, "_nd2_metadata", fake)


def snapshot(folder: Path) -> dict[str, bytes]:
    folder = Path(folder)
    return {p.relative_to(folder).as_posix(): p.read_bytes() for p in sorted(folder.rglob("*")) if p.is_file()}


def read_manifest_file(root: Path) -> dict:
    import yaml

    return yaml.safe_load((Path(root) / "germquant_project.yaml").read_bytes().decode("ascii"))


def last_json(text: str) -> dict:
    """The JSON object a --json command printed (the whole of stdout)."""
    return json.loads(text)


def heavy_after(code: str, heavy=LIGHT_FORBIDDEN, cwd=None, **env) -> tuple[str, subprocess.CompletedProcess]:
    """(the heavy modules loaded by the end, the finished process) for ``code`` run in a fresh interpreter."""
    full = _PROBE.format(heavy=tuple(heavy)) + code
    r = subprocess.run([PYTHON, "-c", full], capture_output=True, text=True, timeout=300,
                       env=subprocess_env(**env), cwd=cwd or ROOT, encoding="utf-8", errors="replace")
    m = re.search(r"HEAVY=(.*)", r.stderr)
    assert m, f"no probe line; rc={r.returncode} stderr={r.stderr[-2000:]}"
    return m.group(1).strip(), r


def module_main_code(module: str, argv: list[str]) -> str:
    """Python source that runs ``python -m <module> <argv>`` through runpy (so the probe sees the run) and prints
    the exit code to stderr as ``EXIT <code>``."""
    return ("import runpy, sys\n"
            f"sys.argv = [{module!r}] + {list(argv)!r}\n"
            "try:\n"
            f"    runpy.run_module({module!r}, run_name='__main__', alter_sys=True)\n"
            "except SystemExit as e:\n"
            "    print('EXIT', e.code, file=sys.stderr)\n"
            "else:\n"
            "    print('EXIT', 0, file=sys.stderr)\n")
