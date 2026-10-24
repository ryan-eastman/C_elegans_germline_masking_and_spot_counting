"""Single-plane time series are refused in this release (tranche 5 package P2b; design 4.13 and 12.3, P2b acceptance 8;
decision D9 of the Train B plan; docs/PLANAR.md, "Time series").

A recipe with io.planar true, or auto on an image whose header has one plane, cannot run timepoints:

* `germquant run FILE --t 0-1` exits with code 2 before any timepoint, printing planar.PLANAR_SERIES_REFUSED;
* api.run(..., t=0) raises PlanarConfigError with the same message, before any file is written;
* process_series (timelapse/series.py, which this package does not edit) records every timepoint as failed with
  that message, because process_image refuses a timepoint of a single-plane run.

A series with io.planar false still runs as before; auto on a 3D series is not refused.
"""
from __future__ import annotations

from pathlib import Path

import planar_helpers as H
import pytest

from germquant import api, cli, pipeline, planar
from germquant.timelapse import process_series

ROOT = Path(__file__).resolve().parents[1]
ND2 = f"{H.STEM}.nd2"


def _study(tmp_path: Path, planar_value: str = "auto") -> Path:
    p = tmp_path / "series_2d.yaml"
    p.write_text(f"extends: {(ROOT / 'config' / 'profiles' / 'generic_2d.yaml').as_posix()}\n"
                 f"io:\n  planar: {planar_value}\n"
                 "segmentation:\n  nuclei:\n    method: classical\nspots:\n  enabled: false\n", encoding="utf-8")
    return p


def test_planar_series_refused(tmp_path, monkeypatch, capsys):
    H.patch_reader(monkeypatch, series=2)
    msg = planar.PLANAR_SERIES_REFUSED
    assert "io.planar" in msg and "docs/PLANAR.md" in msg
    seen = []
    real = pipeline.process_image
    monkeypatch.setattr(pipeline, "process_image", lambda *a, **k: seen.append(k.get("t")) or real(*a, **k))
    # the command line: refused before the loop, for auto on a one-plane header and for true
    for value in ("auto", "true"):
        out = tmp_path / f"cli_{value}"
        rc = cli.main(["run", ND2, "--config", str(_study(tmp_path, value)), "--out", str(out), "--t", "0-1"])
        assert rc == 2
        assert msg in capsys.readouterr().err
        assert not out.exists() and not seen
    # the API, given t: process_image refuses before any file is written
    cfg = H.generic_2d(**{"spots.enabled": False})
    with pytest.raises(planar.PlanarConfigError) as err:
        api.run(ND2, config=cfg, out=tmp_path / "api", t=0)
    assert str(err.value) == msg and not list((tmp_path / "api").iterdir())
    # process_series (not edited): every timepoint fails with the message
    res = process_series(ND2, cfg, tmp_path / "series", timepoints="all")
    assert res["done"] == [] and [f["t_index"] for f in res["failed"]] == [0, 1]
    for f in res["failed"]:
        assert f["flags"] == [f"EXCEPTION:PlanarConfigError: {msg}"]
    assert seen == [0, 1]                       # the two timepoints reached process_image, and each was refused


def test_other_series_still_run(tmp_path, monkeypatch, capsys):
    """io.planar false runs a one-plane series as before (each timepoint through the 3D path); under auto a 3D
    series is not refused."""
    H.patch_reader(monkeypatch, series=2)
    out = tmp_path / "off"
    rc = cli.main(["run", ND2, "--config", str(_study(tmp_path, "false")), "--out", str(out), "--t", "0"])
    assert rc == 0 and (out / "t000" / f"{H.STEM}_t000__done.json").exists()
    assert "planar" not in (out / "t000" / f"{H.STEM}_t000__done.json").read_text(encoding="utf-8")
    H.patch_reader(monkeypatch, H.stack_3d(), two_d=False, series=2)
    assert cli._series_refusal(_load(_study(tmp_path, "auto")), ND2) is None


def _load(path):
    from germquant.config import load_config

    return load_config(path)
