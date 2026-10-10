"""Envelope, audit and acquisition stages on a synthetic 4-colour gonad with a lamin ring channel."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import ndimage as ndi

from germquant import pipeline
from germquant.config import load_config
from germquant.envelope import (
    cytoplasm_shell,
    envelope_labels,
    no_envelope_ids,
    ring_scores,
    run_envelope,
    territories,
)
from germquant.io.acquisition import acquisition_fields, parse_description
from germquant.io.nd2_reader import Stack

SPACING = (0.4, 0.2, 0.2)
CHANNELS = ["405", "477", "545", "640"]     # ccw77 map: DAPI / PGL-1 / SYP-3 / LMN-1


def _stack():
    Z, Y, X = 20, 80, 360
    zz, yy, xx = np.indices((Z, Y, X)).astype(np.float32)
    dapi = np.zeros((Z, Y, X), np.float32)
    lamin = np.zeros((Z, Y, X), np.float32)
    syp = np.zeros((Z, Y, X), np.float32)
    pgl = np.zeros((Z, Y, X), np.float32)
    for cx in range(20, 210, 34):                              # 6 well-separated germline nuclei along x
        r = np.sqrt(((zz - 10) * SPACING[0]) ** 2 + ((yy - 40) * SPACING[1]) ** 2 + ((xx - cx) * SPACING[2]) ** 2)
        dapi += np.exp(-r ** 2 / (2 * 1.6 ** 2)) * 1500.0
        lamin += np.exp(-((r - 3.0) ** 2) / (2 * 0.2 ** 2)) * 3000.0     # a ring outside the segmented chromatin
        syp[10, 40, max(0, cx - 6):min(X, cx + 6)] = 2000.0
        pgl[10, 60, min(cx + 2, X - 1)] = 4000.0                          # a perinuclear granule, 1 um outside the ring
    for jx in (290, 330):                                                 # off-germline DAPI junk, no ring
        r = np.sqrt(((zz - 10) * SPACING[0]) ** 2 + ((yy - 40) * SPACING[1]) ** 2 + ((xx - jx) * SPACING[2]) ** 2)
        dapi += np.exp(-r ** 2 / (2 * 1.6 ** 2)) * 1500.0
    rng = np.random.default_rng(0)
    lamin += rng.normal(0, 20, lamin.shape).astype(np.float32).clip(0)
    data = np.stack([dapi, pgl, syp, lamin], axis=0)
    return Stack(data=data, spacing=SPACING, channel_names=CHANNELS, path="x.nd2", spacing_ok=True)


@pytest.fixture
def patched_reader_lamin(monkeypatch):
    stack = _stack()
    meta = {"sizes": {"C": 4, "Z": 20, "Y": 80, "X": 360}, "spacing": SPACING, "spacing_ok": True,
            "channel_names": CHANNELS, "dtype": "float32", "is_2d": False}
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    monkeypatch.setattr(pipeline, "read_stack", lambda *a, **k: stack)


def _cfg(**switches):
    cfg = load_config("config/config_ccw77.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    cfg.set("coloc.enabled", False)
    for k, v in switches.items():
        cfg.set(k, v)
    return cfg


def test_envelope_functions_on_a_ring():
    st = _stack()
    labels, n = ndi.label(st.data[0] > 400)           # chromatin blobs (r about 2.6 um), one label per nucleus
    labels = labels.astype(np.int32)
    ids = list(range(1, n + 1))
    assert n == 8                                       # 6 germline nuclei + 2 junk blobs
    env, per = envelope_labels(labels, ids, st.data[3], SPACING)
    assert set(np.unique(env)) - {0} <= set(ids)
    good = per[~per["envelope_fallback"]]
    assert len(good) >= 5 and (good["envelope_vol_ratio"].between(1.0, 3.0)).all()
    scores, thr = ring_scores(labels, ids, st.data[3], SPACING)
    assert np.isfinite(thr) and thr > 0
    junk = no_envelope_ids(scores)
    # the ring-less junk blobs (x > 280) are flagged, the ringed nuclei are not
    cx = pd.Series({i: c[2] for i, c in zip(ids, ndi.center_of_mass(labels > 0, labels, ids))})
    assert all(cx[j] > 280 for j in junk) and len(junk) == 2
    terr, n_terr, tmap = territories(labels, ids, SPACING)
    assert n_terr >= 2 and set(terr["territory_id"]) and len(terr) == len(ids) and tmap.shape == labels.shape[1:]
    cyto, _dt = cytoplasm_shell(env > 0, SPACING, 2.5)
    assert cyto.any() and not (cyto & (env > 0)).any()
    out = run_envelope(labels, ids, st.data[3], SPACING)
    assert out["envelope_labels"].shape == labels.shape and out["summary"]["n_no_envelope"] == len(junk)


def test_parse_description_and_fields(tmp_path, monkeypatch):
    desc = ("Plane #1:\n Name: 405\n Exposure: 50 ms\n ExW:405; Power: 12.5\n"
            "Plane #2:\n Name: 640\n Exposure: 200 ms\n ExW:640; Power: 30\n")
    p = parse_description(desc)
    # the two original keys, unchanged; the camera / filter keys added later are None in a minimal block
    assert p["405"]["exposure_ms"] == 50.0 and p["405"]["power"] == 12.5 and p["640"]["exposure_ms"] == 200.0
    assert p["405"]["em_filter"] is None and p["405"]["camera_mode"] is None and p["405"]["bit_depth"] is None
    monkeypatch.setattr("germquant.io.acquisition.read_description", lambda path: desc)
    f = acquisition_fields("x.nd2", ["405", "477", "545", "640"], {"dna": 0, "granule": 1, "lamin": 3})
    assert f["exp_ms_dna"] == 50.0 and f["laser_pow_lamin"] == 30.0 and np.isnan(f["exp_ms_granule"])


def test_stages_off_by_default_write_empty_audit(tmp_path, patched_reader_lamin):
    res = pipeline.process_image("20260708_ccw77_IF_pgl1_syp3_lmn1_HS_male_008.nd2", _cfg(), tmp_path / "r", prov=None)
    for s in ("envelope", "audit", "acquisition"):
        assert res["stages"][s]["status"] == "skipped", (s, res["stages"][s])
    assert res["tables"]["mask_audit"].empty
    nuc = pd.read_csv(next((tmp_path / "r").rglob("*__nuclei.csv")))
    assert "has_envelope" in nuc.columns and nuc["has_envelope"].isna().all()


def test_envelope_audit_acquisition_stages_run(tmp_path, patched_reader_lamin, monkeypatch):
    monkeypatch.setattr("germquant.io.acquisition.read_description",
                        lambda path: "Plane #1:\n Name: 640\n Exposure: 90 ms\n ExW:640; Power: 5\n")
    cfg = _cfg(**{"envelope.enabled": True, "audit.enabled": True, "acquisition.read_exposures": True})
    res = pipeline.process_image("20260708_ccw77_IF_pgl1_syp3_lmn1_HS_male_008.nd2", cfg, tmp_path / "r", prov=None)
    for s in ("envelope", "audit", "acquisition"):
        assert res["stages"][s]["status"] == "ran", (s, res["stages"][s])
    out = tmp_path / "r"
    nuc = pd.read_csv(next(out.rglob("*__nuclei.csv")))
    germ = nuc[nuc["in_germline"].astype(bool)]
    assert germ["envelope_volume_um3"].notna().all() and germ["has_envelope"].notna().all()
    assert germ["territory_id"].notna().all()
    assert nuc.loc[~nuc["in_germline"].astype(bool), "envelope_volume_um3"].isna().all()
    isum = pd.read_csv(next(out.rglob("*__image_summary.csv"))).iloc[0]
    assert isum["n_envelope_fallback"] >= 0 and np.isfinite(isum["ring_thr"]) and isum["n_territories"] >= 1
    assert isum["exp_ms_lamin"] == 90.0 and isum["laser_pow_lamin"] == 5.0 and np.isnan(isum["exp_ms_dna"])
    audit = pd.read_csv(next(out.rglob("*__mask_audit.csv")))
    assert set(audit["kind"]) <= {"lamin_candidate", "germ_label"} and (audit["kind"] == "germ_label").sum() == len(germ)
    gl = audit[audit["kind"] == "germ_label"]
    assert gl["ring_shell"].notna().all() and gl["ring_inside"].notna().all()       # the script's ring evidence
    assert np.allclose(gl["ring_ratio"], gl["ring_shell"] / np.maximum(gl["ring_inside"], 1.0))
    assert isum["envelope_vol_ratio"] > 0                                            # total envelope / total DAPI volume
    # label TIFs are real files with the voxel size, not stubs, and read back as int32
    import tifffile

    lab_tif = next(out.rglob("*__nuclei_labels.tif"))
    assert lab_tif.stat().st_size > 1000
    arr = tifffile.imread(lab_tif)
    assert arr.dtype == np.int32 and arr.ndim == 3
    with tifffile.TiffFile(lab_tif) as t:
        assert t.is_ome and "PhysicalSizeX" in (t.ome_metadata or "")
    assert next(out.rglob("*__envelope_labels.tif")).exists()
    assert next(out.rglob("*__mask_audit.png")).exists()
    assert any(f.startswith("envelope:fallback_n=") for f in res["qc_flags"])
    assert any(f.startswith("audit:missed_nuclei=") for f in res["qc_flags"])
