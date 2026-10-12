"""germquant.planar, the 2D path of tranche 5 (docs/DESIGN_system_packs.md section 4, docs/PLANAR.md; package P2a):
the projection settings and plane selectors (a planted sharpest plane, mid, an index, max and mean, the index
out of range), the 2D size gate, segment_planar's control flow (it mirrors segment_nuclei: the backend under
auto and cellpose, the classical fallback, classical never calling the backend, the output contract, planar=True
only for a backend that declares planar_ok), the measurements of 2D objects, check_planar's refusals and the
planar_ok declarations of the stages and backends. No test here calls the real Cellpose."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_segment_backends import registry  # noqa: F401  (fixture: a private segmentation backend registry)

from germquant import planar, schema, stages
from germquant.config import load_config
from germquant.segment import backends as B
from germquant.segment import nuclei as N
from germquant.spots import backends as SB

ROOT = Path(__file__).resolve().parents[1]
SP = (0.3, 0.1, 0.1)


class _Cfg:
    """A stand-in config: get(dotted, default) over a flat dict."""

    def __init__(self, values: dict):
        self._values = dict(values)

    def get(self, key, default=None):
        return self._values.get(key, default)


def _focus_stack(nz=9, z_focus=6, shape=(48, 64), seed=3):
    """A (Z, Y, X) float32 stack: a textured pattern in focus at z_focus and blurred more with every plane
    away from it, plus a little noise."""
    from scipy import ndimage as ndi

    rng = np.random.default_rng(seed)
    yy, xx = np.indices(shape).astype(np.float32)
    sharp = 500.0 + 200.0 * ((np.sin(xx * 0.9) * np.cos(yy * 0.7)) > 0)
    out = np.zeros((nz, *shape), np.float32)
    for z in range(nz):
        s = 1.2 * abs(z - z_focus)
        out[z] = ndi.gaussian_filter(sharp, s) if s else sharp
    return out + rng.normal(0.0, 2.0, out.shape).astype(np.float32)


def _disks(shape=(1, 80, 80), centres=((20, 20), (20, 60), (60, 20), (60, 60)), radius=12, value=1000.0):
    img = np.full(shape, 100.0, np.float32)
    yy, xx = np.indices(shape[1:])
    for cy, cx in centres:
        img[0][(yy - cy) ** 2 + (xx - cx) ** 2 <= radius ** 2] = value
    return img


# ---- settings and plane selectors ------------------------------------------------------------------------

def test_projection_settings_are_inert_by_default_and_validated():
    assert planar.projection_settings(_Cfg({})) is None                      # the code default: none
    for off in ("none", "None", None, False):
        assert planar.projection_settings(_Cfg({"segmentation.nuclei.project": off})) is None
    assert planar.projection_settings(load_config(ROOT / "config" / "config.yaml")) is None
    assert planar.projection_settings(_Cfg({"segmentation.nuclei.project": "max"})) == {"mode": "max", "plane": None}
    assert planar.projection_settings(_Cfg({"segmentation.nuclei.project": "MEAN"})) == {"mode": "mean", "plane": None}
    assert planar.projection_settings(_Cfg({"segmentation.nuclei.project": "plane"})) == {"mode": "plane", "plane": "mid"}
    got = planar.projection_settings(_Cfg({"segmentation.nuclei.project": "plane",
                                           "segmentation.nuclei.project_plane": "sharpest"}))
    assert got == {"mode": "plane", "plane": "sharpest"}
    got = planar.projection_settings(_Cfg({"segmentation.nuclei.project": "plane",
                                           "segmentation.nuclei.project_plane": 4}))
    assert got == {"mode": "plane", "plane": 4}
    with pytest.raises(planar.PlanarConfigError, match="none, max, mean, plane"):
        planar.projection_settings(_Cfg({"segmentation.nuclei.project": "maximum"}))
    with pytest.raises(planar.PlanarConfigError, match="mid, sharpest or a plane index"):
        planar.projection_settings(_Cfg({"segmentation.nuclei.project": "plane",
                                         "segmentation.nuclei.project_plane": "sharp"}))
    # max and mean use every plane and do not read project_plane, so a stale value cannot fail them
    for mode in ("max", "mean"):
        for stale in ("sharp", True, [1, 2]):
            got = planar.projection_settings(_Cfg({"segmentation.nuclei.project": mode,
                                                   "segmentation.nuclei.project_plane": stale}))
            assert got == {"mode": mode, "plane": None}
    assert issubclass(planar.PlanarConfigError, ValueError)


def test_plane_picker_returns_the_planted_sharpest_plane_mid_max_and_mean():
    arr = _focus_stack(nz=9, z_focus=6)
    scores = planar.plane_scores(arr)
    assert int(np.argmax(scores)) == 6 and scores[6] > 2 * np.sort(scores)[-2]
    assert planar.resolve_plane(arr, "sharpest") == 6
    assert planar.resolve_plane(arr, "mid") == 4                     # Z // 2
    assert planar.resolve_plane(_focus_stack(nz=10), "mid") == 5
    assert planar.resolve_plane(arr, 2) == 2 and planar.resolve_plane(arr, "7") == 7
    np.testing.assert_array_equal(planar.project_image(arr, "plane", "sharpest"), arr[6])
    np.testing.assert_array_equal(planar.project_image(arr, "plane", 2), arr[2])
    np.testing.assert_array_equal(planar.project_image(arr, "plane"), arr[4])
    mx = planar.project_image(arr, "max")
    np.testing.assert_array_equal(mx, arr.max(axis=0))
    assert mx.dtype == arr.dtype
    mean = planar.project_image(arr, "mean")
    assert mean.dtype == np.float32 and mean.shape == arr.shape[1:]
    np.testing.assert_allclose(mean, arr.astype(np.float64).mean(axis=0), rtol=1e-6)
    u16 = (arr * 10).astype(np.uint16)
    assert planar.project_image(u16, "max").dtype == np.uint16
    plane = planar.project_image(arr, "plane", 6)
    plane[:] = 0                                                      # a copy: the stack is untouched
    assert arr[6].any()


def test_sharpest_is_the_variance_of_the_laplacian_and_ties_go_to_the_lowest_z():
    from scipy import ndimage as ndi

    arr = _focus_stack(nz=5, z_focus=1)
    p = arr[3].astype(np.float32)
    want = float(np.var(ndi.laplace(ndi.gaussian_filter(p, 1.0)), dtype=np.float64))
    assert planar.sharpness(arr[3]) == want
    flat = np.full((6, 16, 16), 7.0, np.float32)
    assert planar.resolve_plane(flat, "sharpest") == 0                # every plane ties: the lowest z
    twin = np.repeat(flat[:1], 6, axis=0)
    twin[2] = twin[4] = _focus_stack(nz=1, z_focus=0, shape=(16, 16))[0]
    assert planar.resolve_plane(twin, "sharpest") == 2                # two equal sharpest planes: the lower one
    bad = arr.copy()
    bad[0, 0, 0] = np.nan                                             # a plane with a non-finite score never wins
    assert planar.plane_scores(bad)[0] == -np.inf


def test_an_index_out_of_range_raises():
    arr = _focus_stack(nz=9)
    for bad in (9, 12, -1, "-2"):
        with pytest.raises(planar.PlanarConfigError, match="out of range: the stack has 9 plane"):
            planar.resolve_plane(arr, bad)
    for bad in (True, "top", 2.5, None):
        with pytest.raises(planar.PlanarConfigError):
            planar.resolve_plane(arr, bad)
    with pytest.raises(ValueError, match="Z, Y, X"):
        planar.resolve_plane(arr[0], "mid")
    with pytest.raises(planar.PlanarConfigError, match="max, mean or plane"):
        planar.project_image(arr, "sum")


# ---- the size gate, the copy through z and the measurements ---------------------------------------------

def test_drop_small_area_gates_on_physical_area_and_relabels():
    lab = np.zeros((1, 20, 20), np.int32)
    lab[0, 0:1, 0:3] = 4          # 3 px
    lab[0, 5:7, 0:5] = 7          # 10 px
    lab[0, 10:15, 10:20] = 9      # 50 px
    sp = (5.0, 0.5, 0.5)          # 0.25 um2 per pixel; dz plays no part
    got = planar.drop_small_area(lab, sp, 2.0)                        # 8 px minimum
    assert got.dtype == np.int32 and got.shape == lab.shape
    assert sorted(np.unique(got)) == [0, 1, 2]
    assert (got[lab == 4] == 0).all() and (got[lab == 7] == 1).all() and (got[lab == 9] == 2).all()
    keep_all = planar.drop_small_area(lab, sp, 0.0)                   # no gate: every object, relabelled
    assert sorted(np.unique(keep_all)) == [0, 1, 2, 3]
    np.testing.assert_array_equal(planar.drop_small_area(lab[0], (0.5, 0.5), 2.0), got[0])   # (dy, dx) and (Y, X)
    # the rule of nuclei._drop_small: fewer than int(area / pixel area) pixels goes
    assert planar.drop_small_area(lab, sp, 2.5)[lab == 7].max() == 1   # 10 px >= int(10.0)
    assert planar.drop_small_area(lab, sp, 2.75)[lab == 7].max() == 0  # 10 px < int(11.0)


def test_repeat_through_z_and_the_projected_measurements():
    lab2d = np.zeros((1, 12, 12), np.int32)
    lab2d[0, 1:4, 1:5] = 1        # 12 px
    lab2d[0, 6:11, 6:11] = 2      # 25 px
    col = planar.repeat_through_z(lab2d, 5)
    assert col.shape == (5, 12, 12) and col.dtype == np.int32
    assert all(np.array_equal(col[z], lab2d[0]) for z in range(5))
    np.testing.assert_array_equal(planar.repeat_through_z(lab2d[0], 2), col[:2])
    with pytest.raises(ValueError):
        planar.repeat_through_z(np.zeros((2, 3, 3)), 4)

    from germquant.measure import measure_objects

    sp = (0.5, 0.2, 0.25)
    img = np.arange(5 * 12 * 12, dtype=np.float32).reshape(5, 12, 12)
    df = measure_objects(col, {"foci": img}, sp).rename(columns={"label": "nucleus_id"})
    got = planar.projected_measurements(df, lab2d[0], sp)
    assert got["volume_um3"].isna().all()                             # a prism is not a cell volume
    np.testing.assert_allclose(got["area_um2"], [12 * 0.05, 25 * 0.05])
    np.testing.assert_allclose(got["equivalent_diameter_um"], 2 * np.sqrt(np.array([0.6, 1.25]) / np.pi))
    np.testing.assert_allclose(got["centroid_z_um"], 2 * 0.5)          # the stack middle, never NaN
    np.testing.assert_allclose(got["foci_mean_intensity"], df["foci_mean_intensity"])   # the 3D column mean
    assert (got["n_voxels"] == df["n_voxels"]).all() and df["volume_um3"].notna().all()  # the input is untouched
    empty = pd.DataFrame(columns=["nucleus_id"])
    assert planar.projected_measurements(empty, lab2d, sp) is empty


def test_planar_columns_are_declared_outside_tables():
    from germquant.pipeline import _conform_schema

    assert schema.PLANAR_COLUMNS == {"nuclei": ["area_um2", "equivalent_diameter_um"]}
    assert not set(schema.PLANAR_COLUMNS) - set(schema.TABLES)
    assert not ({"area_um2", "equivalent_diameter_um"} & {c for cols in schema.TABLES.values() for c in cols})
    df = pd.DataFrame({"dna_mean_intensity": [1.0], "nucleus_id": [1], "area_um2": [2.0]})
    plain = _conform_schema(df, "nuclei")
    assert list(plain.columns[:len(schema.NUCLEI)]) == schema.NUCLEI
    assert list(plain.columns[len(schema.NUCLEI):]) == ["dna_mean_intensity", "area_um2"]
    two_d = _conform_schema(df, "nuclei", schema.PLANAR_COLUMNS)
    assert list(two_d.columns) == [*schema.NUCLEI, "area_um2", "equivalent_diameter_um", "dna_mean_intensity"]
    assert two_d["equivalent_diameter_um"].isna().all()
    assert list(_conform_schema(df, "spots", schema.PLANAR_COLUMNS).columns[:len(schema.SPOTS)]) == schema.SPOTS


# ---- segment_planar: the control flow of segment_nuclei, on one plane -------------------------------------

def test_classical_on_a_plane_finds_every_planted_disk(registry):  # noqa: F811
    centres = [(14 + 26 * i, 14 + 26 * j) for i in range(4) for j in range(4)]
    img = _disks(shape=(1, 110, 110), centres=centres, radius=8)
    labels, used = planar.segment_planar(img, SP, method="classical", diameter_um=1.6)
    assert used == "classical" and labels.shape == img.shape and labels.dtype == np.int32
    assert labels.max() == 16
    for cy, cx in centres:                                            # one label per disk, at its centre
        assert labels[0, cy, cx] > 0
    assert len({int(labels[0, cy, cx]) for cy, cx in centres}) == 16
    # the watershed runs unchanged on the (1, Y, X) plane, then the area gate and the relabel
    ref = planar.drop_small_area(N._classical(img, SP, 1.6), SP, 0.0)
    np.testing.assert_array_equal(labels, ref)
    flat, _ = planar.segment_planar(img[0], SP, method="classical", diameter_um=1.6)     # a (Y, X) plane works too
    np.testing.assert_array_equal(flat, labels)
    with pytest.raises(ValueError, match="one plane"):
        planar.segment_planar(np.zeros((2, 8, 8), np.float32), SP, method="classical")


def test_segment_planar_mirrors_the_backend_rules(registry, caplog):  # noqa: F811
    img = _disks()
    seen = []

    def flat2d(dna, spacing, *, model, diameter_um, device, cfg_section, planar=False):
        seen.append({"shape": dna.shape, "planar": planar, "model": model, "diameter_um": diameter_um,
                     "device": device, "cfg_section": cfg_section})
        out = np.zeros(dna.shape, np.int16)
        out[0, :40, :] = 3
        out[0, 40:, :] = 5
        return out

    flat2d.planar_ok = True
    B.register_backend("flat2d", flat2d)
    section = {"backend": "flat2d", "flat2d": {"k": 1}}
    labels, used = planar.segment_planar(img, SP, method="auto", backend="flat2d", cellpose_model="m",
                                         diameter_um=2.5, device="cpu", cfg_section=section)
    assert used == "flat2d" and labels.dtype == np.int32 and sorted(np.unique(labels).tolist()) == [1, 2]
    (c,) = seen
    assert c["planar"] is True and c["shape"] == (1, 80, 80)          # planar_ok: planar=True and one plane
    assert c["model"] == "m" and c["diameter_um"] == 2.5 and c["device"] == "cpu"
    assert c["cfg_section"] == section and c["cfg_section"] is not section          # a copy
    # classical never calls the backend, whatever its name
    seen.clear()
    _, used = planar.segment_planar(img, SP, method="classical", backend="flat2d")
    assert used == "classical" and not seen

    # a backend without planar_ok is called as for a one-plane volume: no planar keyword
    calls = []

    def plain3d(dna, spacing, *, model, diameter_um, device, cfg_section):
        calls.append(dna.shape)
        return (dna > 500).astype(np.int32)

    B.register_backend("plain3d", plain3d)
    labels, used = planar.segment_planar(img, SP, method="cellpose", backend="plain3d")
    assert used == "plain3d" and calls == [(1, 80, 80)] and labels.max() == 1

    # auto falls back to the classical watershed when the backend fails; cellpose re-raises
    def boom(dna, spacing, **kw):
        raise RuntimeError("no GPU today")

    boom.planar_ok = True
    B.register_backend("boom", boom)
    import logging

    with caplog.at_level(logging.WARNING, logger="germquant.planar"):
        labels, used = planar.segment_planar(img, SP, method="auto", backend="boom", diameter_um=2.0)
    assert used == "classical" and labels.max() == 4 and "classical watershed fallback" in caplog.text
    with pytest.raises(RuntimeError, match="no GPU today"):
        planar.segment_planar(img, SP, method="cellpose", backend="boom")

    # an unknown name is a config error under auto and cellpose, never a fallback
    for method in ("auto", "cellpose"):
        with pytest.raises(B.UnknownBackendError, match="nope"):
            planar.segment_planar(img, SP, method=method, backend="nope")
    assert planar.segment_planar(img, SP, method="classical", backend="nope")[1] == "classical"

    # the output contract: an integer (1, Y, X) label image
    def wrong(dna, spacing, **kw):
        return np.zeros(dna.shape[1:], np.int32)

    wrong.planar_ok = True
    B.register_backend("wrong", wrong)
    with pytest.raises(ValueError, match="'wrong'.*shape"):
        planar.segment_planar(img, SP, method="cellpose", backend="wrong")
    assert planar.segment_planar(img, SP, method="auto", backend="wrong")[1] == "classical"

    # the size gate runs after the backend
    labels, _ = planar.segment_planar(img, SP, method="cellpose", backend="flat2d", min_area_um2=1e9)
    assert labels.max() == 0


def test_the_default_backend_name_is_cellpose_and_gets_planar(registry):  # noqa: F811
    seen = []

    def fake_cellpose(dna, spacing, *, model, diameter_um, device, cfg_section, planar=False):
        seen.append(planar)
        return np.ones(dna.shape, np.int32)

    fake_cellpose.planar_ok = True
    B.register_backend("cellpose", fake_cellpose, replace=True)
    for name in (None, "", "cellpose"):
        _, used = planar.segment_planar(_disks(), SP, method="cellpose", backend=name)
        assert used == "cellpose"
    assert seen == [True, True, True]


# ---- check_planar: refuse, never skip silently -----------------------------------------------------------

GONAD_ROLES = {"dna": 0, "central_element": 1, "foci": 2, "granule": None, "lamin": None}
GENERIC_ROLES = {"dna": 0, "foci": 1, "granule": None}


def test_stages_declare_planar_ok_as_the_design_lists():
    field = stages.Stage.__dataclass_fields__["planar_ok"]
    assert list(stages.Stage.__dataclass_fields__)[-1] == "planar_ok" and field.default is False
    ok = {s.name for s in stages.STAGES if s.planar_ok}
    assert ok == {"read", "acquisition", "segment", "measure", "signal_qc", "spots", "qc", "render", "write"}
    assert not stages.BY_NAME["germline"].planar_ok and not stages.BY_NAME["axis"].planar_ok


def test_backends_declare_planar_ok():
    assert planar.planar_ok(B.cellpose_backend) and planar.planar_ok(B.classical_backend)
    assert not planar.planar_ok(B.nuclei3d_backend)
    assert not planar.planar_ok(lambda *a, **k: None)


def test_check_planar_refuses_3d_only_stages_with_their_roles_bound():
    cfg = load_config(ROOT / "config" / "config.yaml")
    with pytest.raises(planar.PlanarConfigError, match="projected run refused") as ei:
        planar.check_planar(cfg, GONAD_ROLES)
    msg = str(ei.value)
    assert "the germline stage" in msg and "germline.enabled: false or --no-germline" in msg
    assert "the axis stage" in msg and "axis.enabled: false or --no-axis" in msg
    assert "granule" not in msg and "coloc" not in msg                # their roles are not bound: they skip
    # the generic profile switches the worm stages off; granule and coloc skip for lack of their roles
    generic = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    planar.check_planar(generic, GENERIC_ROLES)
    generic.set("staging.enabled", True)                              # no roles needed: on means refused
    with pytest.raises(planar.PlanarConfigError, match="the staging stage"):
        planar.check_planar(generic, GENERIC_ROLES)
    generic = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    generic.set("envelope.enabled", True)
    planar.check_planar(generic, GENERIC_ROLES)                       # no lamin channel: it skips itself
    with pytest.raises(planar.PlanarConfigError, match="the envelope stage"):
        planar.check_planar(generic, {**GENERIC_ROLES, "lamin": 2})
    generic = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    with pytest.raises(planar.PlanarConfigError, match="the granule stage.*the coloc stage"):
        planar.check_planar(generic, {**GENERIC_ROLES, "granule": 2, "central_element": 3})


def test_check_planar_refuses_a_segmentation_backend_without_planar_ok(registry):  # noqa: F811
    generic = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    for method in ("auto", "cellpose"):
        generic.set("segmentation.nuclei.method", method)
        for name in ("cellpose", "classical"):
            generic.set("segmentation.nuclei.backend", name)
            planar.check_planar(generic, GENERIC_ROLES)
        generic.set("segmentation.nuclei.backend", "nuclei3d")
        with pytest.raises(planar.PlanarConfigError, match="'nuclei3d' does not declare planar_ok"):
            planar.check_planar(generic, GENERIC_ROLES)
    generic.set("segmentation.nuclei.method", "classical")            # the backend is never called
    planar.check_planar(generic, GENERIC_ROLES)
    generic.set("segmentation.nuclei.method", "auto")
    generic.set("segmentation.nuclei.backend", "nope")
    with pytest.raises(B.UnknownBackendError):
        planar.check_planar(generic, GENERIC_ROLES)


def test_check_planar_spot_backends_only_in_single_plane_runs(monkeypatch):
    def flat(image, labels, spacing, *, cfg_section, marker):
        raise AssertionError("never called")

    def flat_ok(image, labels, spacing, *, cfg_section, marker):
        raise AssertionError("never called")

    flat_ok.planar_ok = True
    monkeypatch.setitem(SB._REGISTRY, "flat_spots", flat)
    monkeypatch.setitem(SB._REGISTRY, "flat_spots_ok", flat_ok)
    generic = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    generic.set("spots.backend", "flat_spots")
    planar.check_planar(generic, GENERIC_ROLES)                       # projected runs count spots in 3D
    with pytest.raises(planar.PlanarConfigError, match="single-plane run refused.*'flat_spots'"):
        planar.check_planar(generic, GENERIC_ROLES, planar=True)
    planar.check_planar(generic, {**GENERIC_ROLES, "foci": None}, planar=True)   # no foci channel: it skips
    generic.set("spots.backend", "flat_spots_ok")
    planar.check_planar(generic, GENERIC_ROLES, planar=True)
    generic.set("spots.instances", [{"name": "second", "role": "granule", "backend": "flat_spots"}])
    planar.check_planar(generic, GENERIC_ROLES, planar=True)          # its role is not bound
    with pytest.raises(planar.PlanarConfigError, match="spots_second"):
        planar.check_planar(generic, {**GENERIC_ROLES, "granule": 2}, planar=True)
    generic.set("spots.enabled", False)
    planar.check_planar(generic, {**GENERIC_ROLES, "granule": 2}, planar=True)
