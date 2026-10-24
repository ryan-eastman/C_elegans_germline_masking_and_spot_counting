"""Single-plane OME-TIFFs read through bioio (tranche 5 package P2b; design 4.10 and 12.3, P2b acceptance 4;
docs/PLANAR.md).

These tests read real files written with tifffile, so they skip without bioio and its TIFF readers (CI installs
both, Q3 of the Train B answers). The calibrated 256 x 256 plane of tests/planar_helpers.py is written as an
OME-TIFF CYX (no z step) and as an OME-TIFF CZYX with Z = 1 (with a z step):

* test_calibrated_plane_runs_planar: under config/profiles/generic_2d.yaml (classical forced) both run as
  single planes with the file's own dy and dx (spacing_source file_yx), raise read:planar and no
  voxel:UNREADABLE flag, and give the tables of the synthetic_2d golden, which the reader seam made from the
  same pixels (tests/test_planar_run.py).
* test_planar_false_is_the_legacy_run: with io.planar false the same files give tables cell-identical to
  tests/golden/single_plane_legacy/<form>/, which was made on the tranche 5 base commit's code (8eb73de), before
  any P2b change, so this compares the new code with the old one and not with itself. The legacy behaviour the
  fence pins for one plane holds as well: the CYX file keeps the 1 um guess and its flag, and every run keeps
  the fence's key sets.
* test_the_reader_fills_spacing_yx_for_one_plane: bioio_reader fills Stack.spacing_yx and the header key
  spacing_yx_um when Z is 1, with xy_stride applied, and leaves _spacing's (1, 1, 1) rule as it is.

Making the legacy golden (docs/PLANAR.md, "The legacy golden"): with a detached worktree of the base commit at
BASE, from this checkout, ``GERMQUANT_UPDATE_GOLDEN=1 PYTHONPATH=BASE/src python -m pytest
tests/test_planar_bioio.py -k legacy``. The update refuses to run on code that has planar.decide.
"""
from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import planar_helpers as H
import pytest

from germquant import pipeline, planar, schema

ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "tests" / "golden" / "single_plane_legacy"
GOLD_2D = ROOT / "tests" / "golden" / "synthetic_2d"
PINS = json.loads((ROOT / "tests" / "golden" / "default_invariants.json").read_text(encoding="utf-8"))
TABLES = tuple(schema.TABLES)
HAVE_OME_TIFF = (importlib.util.find_spec("bioio") is not None
                 and importlib.util.find_spec("bioio_ome_tiff") is not None)
HAVE_TIFFFILE_PLUGIN = (importlib.util.find_spec("bioio") is not None
                        and importlib.util.find_spec("bioio_tifffile") is not None)
pytestmark = pytest.mark.skipif(not HAVE_OME_TIFF, reason="bioio + bioio-ome-tiff not installed")
FORMS = ("ome_cyx", "ome_czyx_z1")
LEGACY_FLAG = "voxel:UNREADABLE_assumed_isotropic_1um"


def _run(tmp_path: Path, form: str, cfg) -> tuple[dict, Path]:
    image = H.FORMS[form](tmp_path / "in" / f"{H.STEM}.ome.tif")
    out = tmp_path / "out"
    return pipeline.process_image(image, cfg, out, prov=None), out


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


@pytest.mark.parametrize("form", FORMS)
def test_planar_false_is_the_legacy_run(tmp_path, form):
    cfg = H.generic_2d(**{"io.planar": False})
    res, out = _run(tmp_path, form, cfg)
    iid = res["image_id"]
    if os.environ.get("GERMQUANT_UPDATE_GOLDEN"):
        if hasattr(planar, "decide"):
            pytest.fail("the single_plane_legacy golden is made on the tranche 5 base commit's code, before any P2b "
                        "change (docs/PLANAR.md), never on code that has planar.decide")
        H.write_golden(out, iid, LEGACY / form, TABLES, cfg.hash, sentinel=LEGACY / "config_hash.txt")
        pytest.skip(f"single_plane_legacy/{form} regenerated on the base code")
    sentinel = (LEGACY / "config_hash.txt").read_text(encoding="utf-8").strip()
    assert cfg.hash == sentinel == PINS["golden_sentinels"]["single_plane_legacy"], (
        "generic_2d.yaml, a base it extends or its channel map changed: the legacy golden's hash no longer reproduces")
    assert iid == H.STEM
    H.compare_tables(out, iid, LEGACY / form, TABLES)
    # the legacy single-plane behaviour the fence pins (tests/test_default_invariants.py)
    pin = PINS["runs"]["single_plane"]
    assert sorted(_json(out / f"{iid}__done.json")) == pin["done_marker_keys"]
    assert sorted(_json(out / "run_manifest.json")) == pin["manifest_keys"]
    assert sorted(_json(out / f"{iid}__stages.json")) == pin["stage_record_keys"]
    flags = res["qc_flags"]
    assert not [f for f in flags if f.startswith(("read:planar", "voxel:override", "voxel:UNREADABLE_xy"))]
    nuclei = pd.read_csv(out / f"{iid}__nuclei.csv")
    assert len(nuclei) > 0 and "area_um2" not in nuclei.columns
    if form == "ome_cyx":              # no z step: today's (1, 1, 1) guess, its flag, volumes as voxel counts
        assert PINS["single_plane"]["flag"] == LEGACY_FLAG and LEGACY_FLAG in flags
        assert (nuclei["volume_um3"] == nuclei["n_voxels"]).all()
        assert (nuclei[["voxel_dz_um", "voxel_dy_um", "voxel_dx_um"]] == 1.0).all().all()
    else:                              # a z step: a calibrated one-plane 3D stack
        assert LEGACY_FLAG not in flags
        np.testing.assert_allclose(nuclei["volume_um3"], nuclei["n_voxels"] * H.DZ * H.DY * H.DX)
        assert (nuclei["voxel_dz_um"] == H.DZ).all() and (nuclei["voxel_dy_um"] == H.DY).all()


@pytest.mark.parametrize("form", FORMS)
def test_calibrated_plane_runs_planar(tmp_path, form):
    cfg = H.generic_2d()
    res, out = _run(tmp_path, form, cfg)
    iid = res["image_id"]
    flags = res["qc_flags"]
    assert "read:planar" in flags and not [f for f in flags if f.startswith("voxel:")], flags
    manifest = _json(out / "run_manifest.json")
    assert manifest["planar"] is True and manifest["spacing_source"] == "file_yx"
    assert _json(out / f"{iid}__done.json")["planar"] is True
    nuclei = pd.read_csv(out / f"{iid}__nuclei.csv")
    assert (nuclei["voxel_dy_um"] == H.DY).all() and (nuclei["voxel_dx_um"] == H.DX).all()
    assert nuclei["voxel_dz_um"].isna().all() and nuclei["volume_um3"].isna().all()
    np.testing.assert_allclose(nuclei["area_um2"], nuclei["n_voxels"] * H.DY * H.DX)
    # the same pixels and the same dy and dx as the synthetic_2d golden, read through bioio instead of the seam
    assert (GOLD_2D / "config_hash.txt").read_text(encoding="utf-8").strip() == cfg.hash
    H.compare_tables(out, iid, GOLD_2D, TABLES)


def test_the_reader_fills_spacing_yx_for_one_plane(tmp_path):
    from germquant.io import read_image, read_image_metadata

    forms = ["ome_yx", "ome_cyx", "ome_czyx_z1"] + (["tiff"] if HAVE_TIFFFILE_PLUGIN else [])
    for form in forms:
        suffix = ".tif" if form == "tiff" else ".ome.tif"
        p = H.FORMS[form](tmp_path / form / f"{H.STEM}{suffix}")
        st, meta = read_image(p), read_image_metadata(p)
        assert st.data.shape[1] == 1 and meta["sizes"]["Z"] == 1, form
        assert st.spacing_yx == (H.DY, H.DX) and meta["spacing_yx_um"] == (H.DY, H.DX), form
        if form == "ome_czyx_z1":       # _spacing is untouched: the z step makes the file calibrated in 3D
            assert st.spacing == (H.DZ, H.DY, H.DX) and st.spacing_ok
        else:                           # no z step: today's (1, 1, 1) guess with spacing_ok False
            assert st.spacing == (1.0, 1.0, 1.0) and st.spacing_ok is False, form
        half = read_image(p, xy_stride=2)
        assert half.spacing_yx == (2 * H.DY, 2 * H.DX) and half.data.shape[2:] == (H.SIZE // 2, H.SIZE // 2)
    p3 = H.write_ome_3d(tmp_path / "3d" / f"{H.STEM}.ome.tif")
    st, meta = read_image(p3), read_image_metadata(p3)
    assert st.spacing_yx is None and "spacing_yx_um" not in meta and st.spacing == (H.DZ, H.DY, H.DX)
    one = read_image(p3, z_range=(2, 3))           # a 3D file cut to one plane is still a 3D file: no spacing_yx
    assert one.data.shape[1] == 1 and one.spacing_yx is None
    if HAVE_TIFFFILE_PLUGIN:
        # a plain TIFF without calibration (resolution unit none): bioio-tifffile reports 1.0, which is not a length
        plain = H.write_tiff(tmp_path / "plain" / f"{H.STEM}.tif", calibrated=False)
        st, meta = read_image(plain), read_image_metadata(plain)
        assert st.spacing_yx is None and meta["spacing_yx_um"] is None and not st.spacing_ok
