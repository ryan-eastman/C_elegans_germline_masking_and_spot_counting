"""Full single-image pipeline via process_image, with the .nd2 reader monkeypatched to feed
synthetic data. Exercises the orchestration, schema-conformed table writing, montage + label
image, and provenance manifest end to end on CPU — the path no other test covered.
"""
import numpy as np
import pandas as pd
import pytest

from germquant import pipeline, schema
from germquant.config import load_config
from germquant.io.nd2_reader import Stack

SPACING = (0.4, 0.2, 0.2)
CHANNELS = ["405", "477", "545"]  # match the n2 channel map (DAPI / SYP-3 / RAD-51)


def _synthetic_stack():
    Z, Y, X = 20, 80, 320
    dapi = np.zeros((Z, Y, X), np.float32)
    syp = np.zeros((Z, Y, X), np.float32)
    rad = np.zeros((Z, Y, X), np.float32)
    zz, yy, xx = np.indices((Z, Y, X)).astype(np.float32)
    # germline: a connected chain of SYP-positive nuclei (each with an SC filament + a RAD-51 focus)
    for cx in range(20, 210, 17):                             # ~12 nuclei along x at y=40
        r2 = (((zz - 10) * SPACING[0]) ** 2 + ((yy - 40) * SPACING[1]) ** 2
              + ((xx - cx) * SPACING[2]) ** 2)
        dapi += np.exp(-r2 / (2 * 1.6 ** 2)) * 1500.0
        syp[10, 40, max(0, cx - 6):min(X, cx + 6)] = 2000.0   # a short SC filament
        rad[10, 42, min(cx + 1, X - 1)] = 6000.0              # one RAD-51 focus
    # off-germline junk: a SEPARATED DAPI cluster (gut/debris), NO SYP, >12 µm from the gonad ->
    # the seeded connected-component isolation must drop it (it carries no synapsed seed).
    for jx in (280, 295, 310):
        r2 = (((zz - 10) * SPACING[0]) ** 2 + ((yy - 40) * SPACING[1]) ** 2
              + ((xx - jx) * SPACING[2]) ** 2)
        dapi += np.exp(-r2 / (2 * 1.6 ** 2)) * 1500.0
    data = np.stack([dapi, syp, rad], axis=0)                 # (C, Z, Y, X)
    return Stack(data=data, spacing=SPACING, channel_names=CHANNELS,
                 path="x.nd2", spacing_ok=True)


@pytest.fixture
def patched_reader(monkeypatch):
    stack = _synthetic_stack()
    meta = {"sizes": {"C": 3, "Z": 20, "Y": 64, "X": 160}, "spacing": SPACING,
            "spacing_ok": True, "channel_names": CHANNELS, "dtype": "float32", "is_2d": False}
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    monkeypatch.setattr(pipeline, "read_stack", lambda *a, **k: stack)


def test_process_image_end_to_end(tmp_path, patched_reader):
    cfg = load_config("config/config.yaml")
    # force the classical segmenter (no GPU/cellpose in CI) for determinism
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    out = tmp_path / "results"
    nd2 = "20251105_n2_nohs_HERM_001.nd2"

    res = pipeline.process_image(nd2, cfg, out, prov=None)

    assert res["n_nuclei"] >= 1
    # metadata parsed from the filename
    assert (out / f"{res['image_id']}__nuclei.csv").exists()
    assert (out / f"{res['image_id']}__montage.png").exists()
    assert (out / f"{res['image_id']}__nuclei_labels.tif").exists()
    assert (out / "run_manifest.json").exists()

    # every table is written for every declared schema with the shared metadata block
    for name, cols in schema.TABLES.items():
        f = out / f"{res['image_id']}__{name}.csv"
        assert f.exists(), f"missing {name} table"
        df = pd.read_csv(f)
        missing = set(cols) - set(df.columns)
        assert not missing, f"{name} table missing schema columns {missing}"
        for meta_col in ("image_id", "genotype", "sex", "voxel_dz_um", "git_sha"):
            assert meta_col in df.columns

    nuclei = pd.read_csv(out / f"{res['image_id']}__nuclei.csv")
    assert (nuclei["voxel_dz_um"] == SPACING[0]).all()
    assert nuclei["sex"].iloc[0] == "herm"
    assert "axis_position_um" in nuclei and nuclei["axis_position_um"].notna().any()

    # germline isolation: the SYP-negative junk blobs are flagged out, the SYP+ germline kept
    assert "in_germline" in nuclei.columns
    summary = pd.read_csv(out / f"{res['image_id']}__image_summary.csv")
    n_germ = int(summary["n_germline_nuclei"].iloc[0])
    assert 5 <= n_germ < int(summary["n_nuclei"].iloc[0])   # some dropped, enough kept for an axis
    # the axis was fit on germline nuclei only (excluded rows carry NaN axis)
    assert nuclei.loc[nuclei["in_germline"] != True, "axis_position_um"].isna().all()  # noqa: E712

    # a 3-channel run has no PGL-1 (granule) channel -> the coloc stage is SKIPPED silently:
    # empty granules/coloc tables, no failure flag.
    assert pd.read_csv(out / f"{res['image_id']}__granules.csv").empty
    assert pd.read_csv(out / f"{res['image_id']}__coloc.csv").empty
    assert not any("coloc:FAILED" in f for f in res["qc_flags"])
    assert not (out / f"{res['image_id']}__granules_labels.tif").exists()


CHANNELS4 = ["405", "477", "545", "647"]  # DAPI / SYP / RAD-51 / PGL-1


def _synthetic_stack_4ch():
    Z, Y, X = 20, 80, 320
    dapi = np.zeros((Z, Y, X), np.float32)
    syp = np.zeros((Z, Y, X), np.float32)
    rad = np.zeros((Z, Y, X), np.float32)
    pgl = np.zeros((Z, Y, X), np.float32)
    zz, yy, xx = np.indices((Z, Y, X)).astype(np.float32)
    for cx in range(20, 210, 17):                             # ~12 germline nuclei along x at y=40
        r2 = (((zz - 10) * SPACING[0]) ** 2 + ((yy - 40) * SPACING[1]) ** 2
              + ((xx - cx) * SPACING[2]) ** 2)
        dapi += np.exp(-r2 / (2 * 1.6 ** 2)) * 1500.0
        syp[10, 40, max(0, cx - 6):min(X, cx + 6)] = 2000.0   # intranuclear SC ribbon
        rad[10, 42, min(cx + 1, X - 1)] = 6000.0              # RAD-51 focus
        # PGL-1 p-granule just OUTSIDE the nucleus (perinuclear shell, y~44), with a cytoplasmic
        # SYP aggregate planted at the SAME spot so syp_aggregate<->PGL-1 overlap is non-trivial.
        pgl[9:12, 43:46, max(0, cx - 1):min(X, cx + 2)] = 4000.0
        syp[9:12, 43:46, max(0, cx - 1):min(X, cx + 2)] = 3000.0
    data = np.stack([dapi, syp, rad, pgl], axis=0)
    return Stack(data=data, spacing=SPACING, channel_names=CHANNELS4,
                 path="x.nd2", spacing_ok=True)


@pytest.fixture
def patched_reader_4ch(monkeypatch):
    stack = _synthetic_stack_4ch()
    meta = {"sizes": {"C": 4, "Z": 20, "Y": 80, "X": 320}, "spacing": SPACING,
            "spacing_ok": True, "channel_names": CHANNELS4, "dtype": "float32", "is_2d": False}
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    monkeypatch.setattr(pipeline, "read_stack", lambda *a, **k: stack)


def test_process_image_4channel_coloc(tmp_path, patched_reader_4ch):
    """4-channel run: PGL-1 granules are surfaced and colocalized with the two SYP operands."""
    cfg = load_config("config/config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    cfg.set("coloc.n_random", 20)                 # keep the null cheap in CI
    out = tmp_path / "results4"

    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, out, prov=None)

    assert not any("coloc:FAILED" in f for f in res["qc_flags"])
    granules = pd.read_csv(out / f"{res['image_id']}__granules.csv")
    assert len(granules) > 0
    assert (granules["marker"] == "PGL-1").all()
    assert {"overlaps_syp_aggregate", "nearest_syp_aggregate_um"} <= set(granules.columns)

    coloc = pd.read_csv(out / f"{res['image_id']}__coloc.csv")
    assert set(coloc["sc_operand"]) == {"syp_aggregate", "sc_ribbon"}

    summary = pd.read_csv(out / f"{res['image_id']}__image_summary.csv")
    assert int(summary["n_granules"].iloc[0]) > 0
    assert "manders_m1_syp_aggregate" in summary.columns

    nuclei = pd.read_csv(out / f"{res['image_id']}__nuclei.csv")
    assert "n_granules" in nuclei.columns
    assert nuclei["n_granules"].sum() > 0

    # the Imaris label image for granules is written; overlap is reported for the aggregate operand
    assert (out / f"{res['image_id']}__granules_labels.tif").exists()
    agg = coloc[coloc["sc_operand"] == "syp_aggregate"].iloc[0]
    assert agg["n_granules"] > 0
    assert agg["frac_granules_overlapping_sc"] >= 0


def test_process_image_no_spots_segmentation_only(tmp_path, patched_reader):
    """--no-spots path: segmentation still runs, the spot stage is skipped (never touches SpotMAX),
    and it's recorded as a deliberate choice — not flagged as a detection failure."""
    cfg = load_config("config/config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    h_before = cfg.hash
    cfg.set("spots.enabled", False)
    assert cfg.get("spots.enabled") is False
    assert cfg.hash != h_before                          # provenance-distinguishable from a normal run

    out = tmp_path / "seg_only"
    res = pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, out, prov=None)

    assert res["n_nuclei"] >= 1
    assert (out / f"{res['image_id']}__nuclei_labels.tif").exists()   # segmentation produced
    assert not (out / f"{res['image_id']}__spots.tif").exists()       # no spots image (nothing to draw)
    spots = pd.read_csv(out / f"{res['image_id']}__spots.csv")
    assert len(spots) == 0                                            # spot stage skipped entirely
    flags = str(pd.read_csv(out / f"{res['image_id']}__image_summary.csv")["qc_flags"].iloc[0])
    assert "disabled_segmentation_only" in flags and "no_spots_detected" not in flags
