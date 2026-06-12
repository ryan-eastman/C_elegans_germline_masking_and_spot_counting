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
    Z, Y, X = 20, 64, 160
    dapi = np.zeros((Z, Y, X), np.float32)
    syp = np.zeros((Z, Y, X), np.float32)
    rad = np.zeros((Z, Y, X), np.float32)
    zz, yy, xx = np.indices((Z, Y, X)).astype(np.float32)
    for cx in (24, 52, 80, 108, 136, 150 - 2):  # 6 nuclei along x
        r2 = (((zz - 10) * SPACING[0]) ** 2 + ((yy - 32) * SPACING[1]) ** 2
              + ((xx - cx) * SPACING[2]) ** 2)
        dapi += np.exp(-r2 / (2 * 1.6 ** 2)) * 1500.0
        syp[10, 32, max(0, cx - 6):min(X, cx + 6)] = 2000.0   # a short SC filament
        rad[10, 34, min(cx + 1, X - 1)] = 6000.0              # one RAD-51 focus
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
