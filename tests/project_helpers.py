"""Shared helpers for the project tests (docs/PROJECTS.md; tranche 5 Train B unit P5-PRJ).

* ``heavy_after(code)`` runs ``code`` in a fresh interpreter (sys.executable, this checkout's src first on
  PYTHONPATH) with the probe of tests/test_packs_isolation.py and returns the heavy modules it left in
  sys.modules: HEAVY plus READERS (torch, cupy, spotmax, cellpose, bioio, dask.array). A sys.modules check
  never runs in the suite's own process, which has loaded torch long before.
"""
from __future__ import annotations

import re
import subprocess

from packs_helpers import HEAVY, PYTHON, ROOT, subprocess_env

READERS = ("bioio", "dask.array")
LIGHT_FORBIDDEN = HEAVY + READERS

_PROBE = ("import sys, atexit\n"
          "atexit.register(lambda: sys.stderr.write('HEAVY=' + ','.join(m for m in {heavy!r} if m in sys.modules)"
          " + '\\n'))\n")


def heavy_after(code: str, heavy=LIGHT_FORBIDDEN, cwd=None, **env) -> tuple[str, subprocess.CompletedProcess]:
    """(the heavy modules loaded by the end, the finished process) for ``code`` run in a fresh interpreter."""
    full = _PROBE.format(heavy=tuple(heavy)) + code
    r = subprocess.run([PYTHON, "-c", full], capture_output=True, text=True, timeout=300,
                       env=subprocess_env(**env), cwd=cwd or ROOT, encoding="utf-8", errors="replace")
    m = re.search(r"HEAVY=(.*)", r.stderr)
    assert m, f"no probe line; rc={r.returncode} stderr={r.stderr[-2000:]}"
    return m.group(1).strip(), r
