"""Numeric-drift guards for the modular-pipeline migration (docs/ROADMAP_modular_pipeline.md, step 0).

Two layers:

1. ``test_synthetic_golden`` (always runs, CPU, seconds): processes the synthetic gonads from
   test_end_to_end (the 3-channel DAPI/SYP/RAD-51 one and the 4-channel one with a PGL-1 channel, which
   exercises the granule and coloc stages) with the classical segmenter and compares every output table
   cell-for-cell with the checked-in references under tests/golden/synthetic/ and
   tests/golden/synthetic_4ch/. A ``config_hash.txt`` sentinel next to each reference pins the config
   text: an edit to config/config.yaml or the channel map fails the test even if no number moved.
   Regenerate the references ONLY for a deliberately number-changing step, with
   ``GERMQUANT_UPDATE_GOLDEN=1 pytest tests/test_golden.py`` run from the code that is to be certified
   (for the migration: a worktree of the frozen commit, see scripts/make_goldens.py).

2. ``test_real_golden`` (opt-in): set ``GERMQUANT_GOLDEN=<golden tree>`` and ``GERMQUANT_CANDIDATE=<run
   tree>`` (both produced by ``germquant run`` / scripts/make_goldens.py); the test runs
   scripts/regression_diff.py on them and fails on any difference, and on an empty comparison.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from test_end_to_end import patched_reader, patched_reader_4ch  # noqa: F401  (fixtures)

from germquant import pipeline
from germquant.config import load_config

ROOT = Path(__file__).resolve().parents[1]
GOLD_ROOT = Path(__file__).parent / "golden"
TABLES = ("nuclei", "spots", "granules", "coloc", "image_summary")
VOLATILE = ["git_sha", "run_timestamp", "pipeline_version", "file_path", "config_hash"]
CASES = {
    # name: (fixture name, image id, config tweaks)
    "synthetic": ("patched_reader", "20251105_n2_nohs_HERM_001.nd2", {}),
    "synthetic_4ch": ("patched_reader_4ch", "20251105_n2_nohs_HERM_001.nd2", {"coloc.n_random": 20}),
}


def _run_case(tmp_path, name):
    _, nd2, tweaks = CASES[name]
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    for k, v in tweaks.items():
        cfg.set(k, v)
    out = tmp_path / "results"
    res = pipeline.process_image(nd2, cfg, out, prov=None)
    files = {t: next(out.rglob(f"{res['image_id']}__{t}.csv")) for t in TABLES}
    return files, cfg


def _read(p: Path) -> pd.DataFrame:
    # round_trip: pandas' default float parser can land 1 ulp off the exact value of a 17-digit
    # string, which would make a re-written golden disagree with a freshly parsed run.
    df = pd.read_csv(p, low_memory=False, float_precision="round_trip")
    return df.drop(columns=[c for c in VOLATILE if c in df.columns])


# The goldens were made on the lab workstation and are compared exactly there. Other CPUs sum in a different
# vector order, so a correlation or a mean can differ in its last bits (coloc.pearson_r differed by one part
# in 1e14 between the workstation and the GitHub runners). CI sets GERMQUANT_GOLDEN_RTOL (1e-12) to allow
# exactly that and nothing more; leave it unset on the workstation so a real change of any size fails.
_RTOL_ENV = os.environ.get("GERMQUANT_GOLDEN_RTOL")
_FLOAT_TOLERANCE = ({"check_exact": False, "rtol": float(_RTOL_ENV), "atol": 0.0} if _RTOL_ENV
                    else {"check_exact": True})


@pytest.mark.parametrize("name", list(CASES))
def test_synthetic_golden(tmp_path, request, name):
    request.getfixturevalue(CASES[name][0])          # activate the matching synthetic reader
    files, cfg = _run_case(tmp_path, name)
    gold = GOLD_ROOT / name
    if os.environ.get("GERMQUANT_UPDATE_GOLDEN"):
        gold.mkdir(parents=True, exist_ok=True)
        for t, f in files.items():
            _read(f).to_csv(gold / f"{t}.csv", index=False)
        (gold / "config_hash.txt").write_text(cfg.hash + "\n", encoding="utf-8")
        pytest.skip(f"{name} golden regenerated")
    assert gold.exists(), f"no {name} golden checked in; run with GERMQUANT_UPDATE_GOLDEN=1 once"
    sentinel = (gold / "config_hash.txt").read_text(encoding="utf-8").strip()
    assert cfg.hash == sentinel, ("config/config.yaml or its channel map changed: the config_hash of "
                                  "existing runs would no longer reproduce (roadmap rule)")
    for t, f in files.items():
        got, ref = _read(f), pd.read_csv(gold / f"{t}.csv", low_memory=False, float_precision="round_trip")
        # append-only schema rule: every golden column is still present, in the same relative order
        # (new declared columns may land before the per-role extras; nothing is renamed or removed)
        common = [c for c in got.columns if c in set(ref.columns)]
        assert common == list(ref.columns), \
            f"{t}: existing columns renamed, removed or reordered (only adding columns is allowed)"
        assert len(got) == len(ref), f"{t}: row count {len(got)} vs {len(ref)}"
        for col in ref.columns:
            a, b = got[col], ref[col]
            if pd.api.types.is_numeric_dtype(b):
                pd.testing.assert_series_equal(a.astype(float), b.astype(float), check_names=False,
                                               obj=f"{t}.{col}", **_FLOAT_TOLERANCE)
            else:
                assert (a.astype("string").fillna("<NA>") == b.astype("string").fillna("<NA>")).all(), \
                    f"{t}.{col} differs"


@pytest.mark.skipif(not (os.environ.get("GERMQUANT_GOLDEN") and os.environ.get("GERMQUANT_CANDIDATE")),
                    reason="set GERMQUANT_GOLDEN and GERMQUANT_CANDIDATE to compare real runs")
def test_real_golden(tmp_path):
    report = tmp_path / "report.json"
    rc = subprocess.run([sys.executable, str(ROOT / "scripts" / "regression_diff.py"),
                         os.environ["GERMQUANT_GOLDEN"], os.environ["GERMQUANT_CANDIDATE"],
                         "--json", str(report)], check=False).returncode
    assert rc == 0, "real-data regression diff reported differences (see output above)"
    assert json.loads(report.read_text(encoding="utf-8")), "no golden images were compared"
