"""Projected runs through process_image (docs/DESIGN_system_packs.md sections 4.1 to 4.8, docs/PLANAR.md; package
P2a), on a synthetic yeast-like stack: C = 2 (a transmitted-like channel whose focus peaks on plane 3 of 11, and a
foci channel with 3D puncta), eight cells of 12.6 um2 and two specks of debris that the area gate drops.

The recipe is tests/fixtures/planar/projected_yeast_like.yaml (the shape of the yeast recipe of design 3.12): the
classical watershed on the sharpest transmitted plane, an area gate, the labels copied through z, SpotMAX with the
optional keys. One run is checked against the synthetic projected golden (tests/golden/synthetic_projected/, with
its config_hash sentinel, which tests/golden/default_invariants.json pins); the others check the contract: labels
identical in every plane, the area gate, volume_um3 NA and area_um2 present, every spot inside a cell, the flag and
the run geometry, the montage on the segmented plane, and the refusals before any pixel is processed.

Regenerate the golden only for a deliberate change of this path: ``GERMQUANT_UPDATE_GOLDEN=1 pytest
tests/test_projected_run.py``, as for tests/test_golden.py.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_pack_knobs import fake_cellpose  # noqa: F401  (fixture: Cellpose mocked)
from test_segment_backends import registry  # noqa: F401  (fixture: a private segmentation backend registry)

from germquant import pipeline, planar, render, schema
from germquant.config import load_config
from germquant.io.nd2_reader import Stack
from germquant.segment import backends as B

ROOT = Path(__file__).resolve().parents[1]
RECIPE = ROOT / "tests" / "fixtures" / "planar" / "projected_yeast_like.yaml"
GOLD = ROOT / "tests" / "golden" / "synthetic_projected"
GOLD_TABLES = ("nuclei", "spots", "image_summary")
VOLATILE = ["git_sha", "run_timestamp", "pipeline_version", "file_path", "config_hash"]
IMAGE = "20261011_yeastlike_001.nd2"
IMAGE_ID = "20261011_yeastlike_001"

SP = (0.3, 0.1, 0.1)                 # dz, dy, dx um
NZ, NY, NX = 11, 110, 200
Z_FOCUS = 3                          # the planted sharpest plane (mid would be 5)
CELL_R = 20                          # px: 2 um, so a cell is 12.6 um2
CELLS = [(28, 25), (28, 75), (28, 125), (28, 175), (82, 25), (82, 75), (82, 125), (82, 175)]
DEBRIS = [(55, 50), (55, 150)]       # bright specks between the cells: about 5.6 um2 once segmented
DEBRIS_R = 9                         # px; smaller specks get no watershed seed and vanish without any gate
MIN_AREA_UM2 = 8.0                   # the recipe's segmentation.nuclei.min_area_um2: drops the debris, keeps cells
CHANNELS = ["Bright", "Cy3"]


def yeast_like_stack() -> tuple[Stack, list[tuple[int, float, float, float]]]:
    """(stack, planted spots as (cell index, z, y, x) in voxels). The transmitted channel: bright textured disks
    on a flat field, in focus on plane Z_FOCUS and blurred with every plane away from it; the foci channel: a
    flat field, a brighter cytoplasm and 2 to 4 anisotropic Gaussian puncta per cell at several depths."""
    from scipy import ndimage as ndi

    rng = np.random.default_rng(11)
    yy, xx = np.indices((NY, NX)).astype(np.float32)
    focus = np.full((NY, NX), 1000.0, np.float32)
    for cy, cx in CELLS:
        inside = (yy - cy) ** 2 + (xx - cx) ** 2 <= CELL_R ** 2
        focus[inside] = 1600.0 + 120.0 * np.sin(xx[inside] * 1.3) * np.cos(yy[inside] * 1.1)
    for cy, cx in DEBRIS:
        focus[(yy - cy) ** 2 + (xx - cx) ** 2 <= DEBRIS_R ** 2] = 2600.0
    trans = np.zeros((NZ, NY, NX), np.float32)
    for z in range(NZ):
        s = 1.5 * abs(z - Z_FOCUS)
        trans[z] = ndi.gaussian_filter(focus, s) if s else focus
    trans += rng.normal(0.0, 5.0, trans.shape).astype(np.float32)

    zz, yy3, xx3 = np.indices((NZ, NY, NX)).astype(np.float32)
    foci = np.full((NZ, NY, NX), 100.0, np.float32)
    planted = []
    for k, (cy, cx) in enumerate(CELLS):
        foci[(yy3 - cy) ** 2 + (xx3 - cx) ** 2 <= CELL_R ** 2] = 300.0
        n = 2 + k % 3
        for j in range(n):
            ang = 2 * np.pi * j / n + 0.4 * k
            py, px = cy + 9 * np.sin(ang), cx + 9 * np.cos(ang)
            pz = 2 + (3 * j + k) % 7
            d2 = (((zz - pz) * SP[0]) ** 2 / 0.4 ** 2 + ((yy3 - py) * SP[1]) ** 2 / 0.12 ** 2
                  + ((xx3 - px) * SP[2]) ** 2 / 0.12 ** 2)
            foci += 2000.0 * np.exp(-d2 / 2)
            planted.append((k, float(pz), float(py), float(px)))
    foci += rng.normal(0.0, 10.0, foci.shape).astype(np.float32)
    stack = Stack(data=np.stack([trans, foci]), spacing=SP, channel_names=list(CHANNELS), path="x.nd2",
                  spacing_ok=True)
    return stack, planted


def patch_reader(monkeypatch, stack: Stack) -> None:
    c, z, y, x = stack.data.shape
    meta = {"sizes": {"C": c, "Z": z, "Y": y, "X": x}, "spacing": stack.spacing, "spacing_ok": True,
            "channel_names": list(stack.channel_names), "dtype": "float32", "is_2d": False}
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    monkeypatch.setattr(pipeline, "read_stack", lambda *a, **k: stack)


def recipe(**overrides):
    cfg = load_config(RECIPE)
    for key, value in overrides.items():
        cfg.set(key, value)
    return cfg


@pytest.fixture(scope="module")
def golden_run(tmp_path_factory):
    """The recipe run once on the yeast-like stack (shared by the golden and the contract tests)."""
    stack, planted = yeast_like_stack()
    out = tmp_path_factory.mktemp("projected") / "results"
    with pytest.MonkeyPatch.context() as mp:
        patch_reader(mp, stack)
        cfg = recipe()
        res = pipeline.process_image(IMAGE, cfg, out, prov=None)
    return {"res": res, "out": out, "cfg": cfg, "stack": stack, "planted": planted}


def _read(p: Path) -> pd.DataFrame:
    df = pd.read_csv(p, low_memory=False, float_precision="round_trip")
    return df.drop(columns=[c for c in VOLATILE if c in df.columns])


_RTOL_ENV = os.environ.get("GERMQUANT_GOLDEN_RTOL")      # CI: 1e-12 for another CPU's summation order
_FLOAT_TOLERANCE = ({"check_exact": False, "rtol": float(_RTOL_ENV), "atol": 0.0} if _RTOL_ENV
                    else {"check_exact": True})


def test_projected_golden(golden_run):
    out, cfg = golden_run["out"], golden_run["cfg"]
    files = {t: out / f"{IMAGE_ID}__{t}.csv" for t in GOLD_TABLES}
    if os.environ.get("GERMQUANT_UPDATE_GOLDEN"):
        GOLD.mkdir(parents=True, exist_ok=True)
        for t, f in files.items():
            _read(f).to_csv(GOLD / f"{t}.csv", index=False)
        (GOLD / "config_hash.txt").write_text(cfg.hash + "\n", encoding="utf-8")
        pytest.skip("synthetic_projected golden regenerated")
    sentinel = (GOLD / "config_hash.txt").read_text(encoding="utf-8").strip()
    assert cfg.hash == sentinel, ("the projected test recipe, its channel map or a base it extends changed: the "
                                  "golden's config_hash no longer reproduces")
    for t, f in files.items():
        got, ref = _read(f), pd.read_csv(GOLD / f"{t}.csv", low_memory=False, float_precision="round_trip")
        # the append-only schema rule, as in tests/test_golden.py: every golden column is still there, in order
        common = [c for c in got.columns if c in set(ref.columns)]
        assert common == list(ref.columns), f"{t}: existing columns renamed, removed or reordered"
        assert len(got) == len(ref), f"{t}: row count {len(got)} vs {len(ref)}"
        for col in ref.columns:
            a, b = got[col], ref[col]
            if pd.api.types.is_numeric_dtype(b):
                pd.testing.assert_series_equal(a.astype(float), b.astype(float), check_names=False,
                                               obj=f"{t}.{col}", **_FLOAT_TOLERANCE)
            else:
                assert (a.astype("string").fillna("<NA>") == b.astype("string").fillna("<NA>")).all(), \
                    f"{t}.{col} differs"


def test_projected_run_contract(golden_run):
    import tifffile

    res, out, planted = golden_run["res"], golden_run["out"], golden_run["planted"]
    # labels identical in every plane, cells only (the debris fell to the area gate)
    labels = tifffile.imread(out / f"{IMAGE_ID}__nuclei_labels.tif")
    assert labels.shape == (NZ, NY, NX)
    assert all(np.array_equal(labels[z], labels[0]) for z in range(NZ))
    assert res["n_nuclei"] == len(CELLS) == int(labels.max())
    cell_ids = {int(labels[0, cy, cx]) for cy, cx in CELLS}
    assert len(cell_ids) == len(CELLS) and 0 not in cell_ids
    assert all(labels[0, cy, cx] == 0 for cy, cx in DEBRIS)
    # the nuclei table: volume NA, area and equivalent diameter from the 2D labels, after the declared columns
    nuclei = pd.read_csv(out / f"{IMAGE_ID}__nuclei.csv")
    cols = list(nuclei.columns)
    assert cols[:len(schema.NUCLEI)] == schema.NUCLEI
    assert cols[len(schema.NUCLEI):len(schema.NUCLEI) + 2] == ["area_um2", "equivalent_diameter_um"]
    assert nuclei["volume_um3"].isna().all() and nuclei["area_um2"].notna().all()
    planted_area = np.pi * (CELL_R * SP[1]) ** 2
    assert np.allclose(nuclei["area_um2"], planted_area, rtol=0.15)
    np.testing.assert_allclose(nuclei["equivalent_diameter_um"], 2 * np.sqrt(nuclei["area_um2"] / np.pi))
    counts = np.bincount(labels[0].ravel())
    np.testing.assert_allclose(nuclei.set_index("nucleus_id")["area_um2"].sort_index(),
                               counts[1:] * SP[1] * SP[2])
    np.testing.assert_allclose(nuclei["centroid_z_um"], (NZ - 1) / 2 * SP[0])    # the stack middle, never NaN
    assert (nuclei["n_voxels"] == nuclei["nucleus_id"].map(lambda i: counts[i] * NZ)).all()
    assert nuclei["in_germline"].all()                                 # germline off: every object counts
    # every spot lies inside a cell, and each planted punctum is found in its own cell
    spots = pd.read_csv(out / f"{IMAGE_ID}__spots.csv")
    yi = np.round(spots["y_um"] / SP[1]).astype(int)
    xi = np.round(spots["x_um"] / SP[2]).astype(int)
    assert (labels[0][yi, xi] == spots["nucleus_id"]).all() and (spots["nucleus_id"] > 0).all()
    zyx = spots[["z_um", "y_um", "x_um"]].to_numpy() / np.asarray(SP)
    found = 0
    for k, pz, py, px in planted:
        d = np.sqrt((((zyx - [pz, py, px]) * SP) ** 2).sum(axis=1))      # um
        j = int(np.argmin(d))
        if d[j] <= 0.4 and spots["nucleus_id"].iloc[j] == labels[0, CELLS[k][0], CELLS[k][1]]:
            found += 1
    assert found == len(planted) and len(spots) <= len(planted) + 2
    assert (nuclei.set_index("nucleus_id")["n_spots"].sum()) == len(spots)
    # flags, run geometry and the stage record
    flags = res["qc_flags"]
    assert "segment:projected_plane" in flags
    assert not any("few_nuclei" in f for f in flags)                   # qc.min_nuclei 3
    assert not any(f.startswith("segment:no_") for f in flags)
    record = {"mode": "plane", "plane": "sharpest", "z": Z_FOCUS}
    manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    stage_rec = json.loads((out / f"{IMAGE_ID}__stages.json").read_text(encoding="utf-8"))
    assert manifest["segment_projection"] == record and stage_rec["segment_projection"] == record
    assert "segment:projected_plane" in stage_rec["stages"]["segment"]["flags"]
    done = json.loads((out / f"{IMAGE_ID}__done.json").read_text(encoding="utf-8"))
    assert sorted(done) == ["config_hash", "enabled_stages", "failed_stages", "finished", "image_id",
                            "stage_hashes", "xy_stride", "z_range"]       # no new marker key for a projected run
    status = {k: v["status"] for k, v in res["stages"].items()}
    assert status["segment"] == status["measure"] == status["spots"] == status["signal_qc"] == "ran"
    assert res["stages"]["germline"]["reason"] == "disabled" and res["stages"]["axis"]["reason"] == "disabled"
    assert status["granule"] == "skipped" and "missing channel role" in res["stages"]["granule"]["reason"]
    assert status["coloc"] == "skipped" and not [s for s, v in status.items() if v == "failed"]
    isum = pd.read_csv(out / f"{IMAGE_ID}__image_summary.csv")
    assert isum["segmentation_method"].iloc[0] == "classical" and bool(isum["qc_pass"].iloc[0])
    assert (out / f"{IMAGE_ID}__montage.png").exists() and (out / f"{IMAGE_ID}__spots.tif").exists()
    # 3D tables keep their columns: the planar columns appear in the nuclei table only
    for t in schema.TABLES:
        if t != "nuclei":
            assert "area_um2" not in pd.read_csv(out / f"{IMAGE_ID}__{t}.csv", nrows=0).columns


def test_the_area_gate_and_the_projection_modes(tmp_path, monkeypatch):
    stack, _ = yeast_like_stack()
    patch_reader(monkeypatch, stack)
    no_gate = pipeline.process_image(IMAGE, recipe(**{"segmentation.nuclei.min_area_um2": 0.0,
                                                      "spots.enabled": False}), tmp_path / "nogate", prov=None)
    assert no_gate["n_nuclei"] == len(CELLS) + len(DEBRIS)              # the debris is segmented without the gate
    areas = no_gate["tables"]["nuclei"]["area_um2"].sort_values().to_numpy()
    assert (areas[:2] < MIN_AREA_UM2).all() and (areas[2:] > MIN_AREA_UM2).all()
    assert load_config(RECIPE).get("segmentation.nuclei.min_area_um2") == MIN_AREA_UM2
    import tifffile

    for mode in ("max", "mean"):
        # project_plane is ignored by the projections: even a selector that plane mode would refuse
        res = pipeline.process_image(IMAGE, recipe(**{"segmentation.nuclei.project": mode, "spots.enabled": False,
                                                      "segmentation.nuclei.project_plane": "bogus"}),
                                     tmp_path / mode, prov=None)
        assert f"segment:projected_{mode}" in res["qc_flags"] and "segment:projected_plane" not in res["qc_flags"]
        rec = json.loads((tmp_path / mode / f"{IMAGE_ID}__stages.json").read_text(encoding="utf-8"))
        assert rec["segment_projection"] == {"mode": mode, "plane": None, "z": None}
        labels = tifffile.imread(tmp_path / mode / f"{IMAGE_ID}__nuclei_labels.tif")
        assert all(np.array_equal(labels[z], labels[0]) for z in range(NZ))
        assert len({int(labels[0, cy, cx]) for cy, cx in CELLS} - {0}) == len(CELLS)   # every cell, one label each
    res = pipeline.process_image(IMAGE, recipe(**{"segmentation.nuclei.project_plane": "mid", "spots.enabled": False}),
                                 tmp_path / "mid", prov=None)
    rec = json.loads((tmp_path / "mid" / f"{IMAGE_ID}__stages.json").read_text(encoding="utf-8"))
    assert rec["segment_projection"] == {"mode": "plane", "plane": "mid", "z": NZ // 2}


def test_the_plane_is_resolved_on_the_segmentation_channel(tmp_path, monkeypatch):
    """The sharpest plane is picked on the channel the segmentation role names, wherever it sits in the file: with
    the foci channel first, z is still the transmitted channel's focus plane, not the foci channel's sharpest."""
    stack, _ = yeast_like_stack()
    assert planar.resolve_plane(stack.data[1], "sharpest") != Z_FOCUS   # channel 0 would give another plane
    swapped = Stack(data=stack.data[::-1].copy(), spacing=SP, channel_names=["Cy3", "Bright"], path="x.nd2",
                    spacing_ok=True)
    patch_reader(monkeypatch, swapped)
    pipeline.process_image(IMAGE, recipe(**{"spots.enabled": False}), tmp_path / "swapped", prov=None)
    rec = json.loads((tmp_path / "swapped" / f"{IMAGE_ID}__stages.json").read_text(encoding="utf-8"))
    assert rec["segment_projection"] == {"mode": "plane", "plane": "sharpest", "z": Z_FOCUS}


def test_base_plane_shows_the_segmented_plane(tmp_path, monkeypatch):
    """Projected runs draw the outlines on the 2D image the segmenter saw; the default montage call passes no
    base_plane (its pixels are the fence's montage pin)."""
    seen = []
    real = pipeline.make_montage

    def spy(*args, **kw):
        seen.append(kw)
        return real(*args, **kw)

    monkeypatch.setattr(pipeline, "make_montage", spy)
    stack, _ = yeast_like_stack()
    patch_reader(monkeypatch, stack)
    pipeline.process_image(IMAGE, recipe(**{"spots.enabled": False}), tmp_path / "p", prov=None)
    np.testing.assert_array_equal(seen[-1]["base_plane"], stack.data[0][Z_FOCUS])
    assert seen[-1]["base_role"] == "transmitted"
    # the default path: no base_plane keyword
    from test_end_to_end import _synthetic_stack

    gonad = _synthetic_stack()
    patch_reader(monkeypatch, gonad)
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg.set("segmentation.nuclei.method", "classical")
    cfg.set("spots.enabled", False)
    pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, tmp_path / "d", prov=None)
    assert "base_plane" not in seen[-1]

    # make_montage itself: the outline panel's image is base_plane when given, the projection otherwise
    import matplotlib.axes

    shown = []
    real_imshow = matplotlib.axes.Axes.imshow

    def imshow(self, x, *a, **k):
        shown.append(np.asarray(x))
        return real_imshow(self, x, *a, **k)

    monkeypatch.setattr(matplotlib.axes.Axes, "imshow", imshow)
    roles = {"transmitted": 0, "foci": 1}
    labels = np.zeros(stack.shape_zyx, np.int32)
    labels[:, 20:40, 20:40] = 1
    plane = stack.data[0][Z_FOCUS]
    render.make_montage(stack, labels, roles, tmp_path / "m1.png", base_role="transmitted", base_plane=plane)
    assert np.array_equal(shown[2], render._norm(plane))               # panels: transmitted, foci, outlines
    shown.clear()
    render.make_montage(stack, labels, roles, tmp_path / "m2.png", base_role="transmitted")
    assert np.array_equal(shown[2], render._norm(render._mip(stack.data[0])))


def test_refusals_come_before_any_pixel_is_processed(tmp_path, monkeypatch, registry):  # noqa: F811
    calls = []

    def plain3d(dna, spacing, *, model, diameter_um, device, cfg_section):
        calls.append(dna.shape)
        return np.ones(dna.shape, np.int32)

    B.register_backend("plain3d", plain3d)
    stack, _ = yeast_like_stack()
    patch_reader(monkeypatch, stack)
    cases = {
        "axis": ({"axis.enabled": True}, "the axis stage is 3D-only"),
        "backend": ({"segmentation.nuclei.method": "auto", "segmentation.nuclei.backend": "plain3d"},
                    "'plain3d' does not declare planar_ok"),
        "nuclei3d": ({"segmentation.nuclei.method": "cellpose", "segmentation.nuclei.backend": "nuclei3d"},
                     "'nuclei3d' does not declare planar_ok"),
        "index": ({"segmentation.nuclei.project_plane": 11}, "out of range: the stack has 11 plane"),
        "mode": ({"segmentation.nuclei.project": "median"}, "segmentation.nuclei.project must be one of"),
    }
    for name, (overrides, message) in cases.items():
        out = tmp_path / name
        with pytest.raises(planar.PlanarConfigError, match=message):
            pipeline.process_image(IMAGE, recipe(**overrides), out, prov=None)
        assert out.is_dir() and not list(out.iterdir()), f"{name}: files written before the refusal"
    assert not calls
    # a stack without the transmitted channel: require_channel fails the image instead of segmenting channel 0
    no_trans = Stack(data=stack.data[::-1].copy(), spacing=SP, channel_names=["Cy3", "Cy5"], path="x.nd2",
                     spacing_ok=True)
    patch_reader(monkeypatch, no_trans)
    with pytest.raises(ValueError, match="segmentation role 'transmitted'"):
        pipeline.process_image(IMAGE, recipe(), tmp_path / "no_trans", prov=None)
    assert not (tmp_path / "no_trans" / f"{IMAGE_ID}__nuclei_labels.tif").exists()


def test_a_projected_cellpose_run_makes_the_2d_call(tmp_path, monkeypatch, fake_cellpose):  # noqa: F811
    stack, _ = yeast_like_stack()
    patch_reader(monkeypatch, stack)
    cfg = recipe(**{"segmentation.nuclei.method": "cellpose", "spots.enabled": False,
                    "segmentation.nuclei.min_area_um2": 0.0})
    res = pipeline.process_image(IMAGE, cfg, tmp_path / "cp", prov=None)
    (call,) = fake_cellpose.evals
    np.testing.assert_array_equal(call["x"], stack.data[0][Z_FOCUS])
    assert call["kw"] == {"channel_axis": None}
    assert fake_cellpose.inits == [{"device": "device:cpu", "gpu": False, "use_bfloat16": False}]   # cpsam
    assert res["tables"]["image_summary"]["segmentation_method"].iloc[0] == "cellpose"
    assert res["n_nuclei"] >= 1 and "segment:projected_plane" in res["qc_flags"]
    cfg = recipe(**{"segmentation.nuclei.method": "cellpose", "spots.enabled": False,
                    "segmentation.nuclei.cellpose.object_diameter_um": 4.0})
    pipeline.process_image(IMAGE, cfg, tmp_path / "cp_d", prov=None)
    assert fake_cellpose.evals[-1]["kw"] == {"channel_axis": None, "diameter": 4.0 / SP[1]}
