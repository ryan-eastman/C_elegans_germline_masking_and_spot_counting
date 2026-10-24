"""Spot detection on single-plane images (tranche 5 package P2b, unit P2b-SPOTS; docs/DESIGN_system_packs.md
section 4.7, docs/PLANAR.md):

* without planar, detect_spots makes exactly today's SpotMAX calls and today's z-split merge call (a spy on
  every spotmax.pipe call and on _merge_z_columns, as in tests/test_pack_knobs.py);
* the planar branch of detect_spots_spotmax: radii (PLANAR_RZ_PX, r / dy, r / dx), the Gaussian sigma over
  min(dy, dx), the features spacing (1.0, dy, dx) whatever dz the caller passes, no z-split merge, and the stage
  flag z_merge_off_planar in the attrs of both tables (the spots stage raises it as spots:z_merge_off_planar);
  the optional keys of a pack recipe still apply, except the z radius, which has no meaning on one plane;
* planar=True reaches only a backend that declares planar_ok: spotmax does, log and net3d do not, and
  detect_spots refuses them (a second guard: planar.check_planar refuses them first, before any pixel);
* the PLANAR_RZ_PX table: one synthetic plane through the real SpotMAX at 0.5, 1.0 and 2.0, printed for
  docs/PLANAR.md, with 1.0 the shipped value; the planar tables do not depend on the dz the caller passes.
"""
from __future__ import annotations

import logging
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from germquant import planar
from germquant.config import load_config
from germquant.spots import backends as SB
from germquant.spots import detect as SD
from germquant.spots import detect_spots, get_backend

ROOT = Path(__file__).resolve().parents[1]

# the spots stage's call: pipeline.process_image with the spots: section of config/config.yaml
PIPELINE_SPOT_ARGS = dict(spot_radius_um=0.3, gauss_sigma_um=0.08, thresholding_method="threshold_triangle",
                          effect_size_metric="spot_vs_backgr_effect_size_glass", effect_size_min=3.0,
                          merge_z_columns=True, z_merge_gap_um=0.8, z_merge_valley_frac=0.8,
                          max_spot_candidates=30000)


def _spots_section(**extra) -> dict:
    return {**load_config(ROOT / "config" / "config.yaml").get("spots"), **extra}


def _two_nuclei(nz: int = 4):
    """Two box nuclei in a (nz, 10, 10) float32 image of noise; with nz 1 a single plane."""
    labels = np.zeros((nz, 10, 10), np.int32)
    zs = slice(1, 3) if nz > 1 else slice(0, 1)
    labels[zs, 1:4, 1:4] = 1
    labels[zs, 5:8, 5:8] = 2
    img = np.random.default_rng(5).normal(50.0, 5.0, labels.shape).astype(np.float32)
    return img, labels


class _FakeSpotmax:
    """spotmax.pipe and spotmax.filters as call recorders. Each step returns a new array, so a test can tell
    which image went where; detection finds one spot in each nucleus, at the planes in ``zs``."""

    def __init__(self):
        self.zs = [1, 2]
        self.calls: list[tuple[str, tuple, dict]] = []
        self.out: dict[str, object] = {}
        self.pipe = types.SimpleNamespace(
            preprocess_image=self._rec("preprocess_image", lambda image, **kw: np.asarray(image) + 1.0),
            spots_semantic_segmentation=self._rec("spots_semantic_segmentation",
                                                  lambda image, **kw: np.ones(np.shape(image), bool)),
            spot_detection=self._rec("spot_detection", self._detect),
            spots_calc_features_and_filter=self._rec("spots_calc_features_and_filter", self._features),
        )
        self.filters = types.SimpleNamespace(
            DoG_spots=self._rec("DoG_spots", lambda image, radii, **kw: np.asarray(image) * 2.0))

    def _rec(self, name, fn):
        def wrapped(*args, **kw):
            self.calls.append((name, args, kw))
            out = fn(*args, **kw)
            self.out[name] = out
            return out
        return wrapped

    def _detect(self, image, **kw):
        return pd.DataFrame({"z": list(self.zs), "y": [2, 6], "x": [2, 6]}), None

    @staticmethod
    def _features(image, radii, df, **kw):
        feat = pd.DataFrame({"z": df["z"].to_numpy(), "y": df["y"].to_numpy(), "x": df["x"].to_numpy(),
                             SD._EFFECT: [4.0, 5.0], "spot_raw_mean_in_spot": [70.0, 80.0]})
        return None, [feat]

    def names(self) -> list[str]:
        return [c[0] for c in self.calls]

    def call(self, name):
        hits = [c for c in self.calls if c[0] == name]
        assert len(hits) == 1, f"{name} called {len(hits)} times"
        return hits[0]

    def clear(self) -> None:
        self.calls.clear()
        self.out.clear()


@pytest.fixture
def fake_spotmax(monkeypatch):
    fake = _FakeSpotmax()
    monkeypatch.setitem(sys.modules, "spotmax", types.SimpleNamespace(pipe=fake.pipe, filters=fake.filters))
    monkeypatch.setitem(sys.modules, "spotmax.pipe", fake.pipe)
    monkeypatch.setitem(sys.modules, "spotmax.filters", fake.filters)
    return fake


@pytest.fixture
def merge_spy(monkeypatch):
    """Every call of the z-split merge, which the SpotMAX chains look up on the detect module at call time."""
    calls: list[tuple[tuple, dict]] = []
    real = SD._merge_z_columns

    def spy(*args, **kw):
        calls.append((args, kw))
        return real(*args, **kw)

    monkeypatch.setattr(SD, "_merge_z_columns", spy)
    return calls


def _same(a, b) -> bool:
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        return np.array_equal(np.asarray(a), np.asarray(b)) and np.asarray(a).dtype == np.asarray(b).dtype
    return a == b


def _kw_equal(got: dict, want: dict) -> None:
    assert sorted(got) == sorted(want), f"keywords {sorted(got)} != {sorted(want)}"
    for k, v in want.items():
        assert _same(got[k], v), f"keyword {k}: {got[k]!r} != {v!r}"


def _check_spotmax_calls(fake, img, labels, *, radii, gauss, features_spacing, pre_extra=None, sharpened=False):
    """The four SpotMAX calls of a chain, keyword by keyword (and the DoG call when sharpened)."""
    want = ["preprocess_image", "spots_semantic_segmentation", "spot_detection", "spots_calc_features_and_filter"]
    if sharpened:
        want.insert(1, "DoG_spots")
    assert fake.names() == want
    _, args, kw = fake.call("preprocess_image")
    assert len(args) == 1 and _same(args[0], img.astype("float32"))
    _kw_equal(kw, {"lab": labels, "gauss_sigma": gauss, "use_gpu": False, "spots_zyx_radii_pxl": radii,
                   **(pre_extra or {})})
    pre = fake.out["preprocess_image"]
    image = pre
    if sharpened:
        _, args, kw = fake.call("DoG_spots")
        assert len(args) == 2 and args[0] is pre and _same(args[1], radii)
        _kw_equal(kw, {"use_gpu": False, "lab": labels})
        image = fake.out["DoG_spots"]
    _, args, kw = fake.call("spots_semantic_segmentation")
    assert len(args) == 1 and args[0] is image
    _kw_equal(kw, {"lab": labels, "spots_zyx_radii_pxl": radii, "do_try_all_thresholds": False,
                   "return_only_segm": True, "thresholding_method": "threshold_triangle", "use_gpu": False})
    _, args, kw = fake.call("spot_detection")
    assert len(args) == 1 and args[0] is image
    _kw_equal(kw, {"spots_segmantic_segm": fake.out["spots_semantic_segmentation"], "spots_zyx_radii_pxl": radii,
                   "lab": labels, "return_df": True})
    _, args, kw = fake.call("spots_calc_features_and_filter")
    assert len(args) == 3 and args[0] is pre and _same(args[1], radii)
    assert list(args[2].index) == [1, 2] and args[2].index.name == "Cell_ID"
    want_kw = {"lab": labels, "raw_image": img.astype("float32"), "zyx_voxel_size": features_spacing,
               "gop_filtering_thresholds": None}
    if sharpened:
        want_kw["sharp_spots_image"] = image
    _kw_equal(kw, want_kw)
    return pre


# ---- the default call ------------------------------------------------------------------------------------

def test_default_spotmax_call_unchanged(fake_spotmax, merge_spy):
    """Without planar (and with planar=False) detect_spots, through the spotmax backend, makes exactly today's
    SpotMAX calls and today's z-split merge call, and its tables carry no stage flag."""
    img, labels = _two_nuclei()
    sp = (0.4, 0.2, 0.25)
    spa = np.array([0.4, 0.2, 0.25])
    radii = 0.3 / spa                                     # today's: r / spacing on every axis
    for extra in ({}, {"planar": False}):
        fake_spotmax.clear()
        merge_spy.clear()
        per_spot, per_nuc = detect_spots(img, labels, sp, marker="RAD-51", cfg_section=_spots_section(),
                                         **PIPELINE_SPOT_ARGS, **extra)
        pre = _check_spotmax_calls(fake_spotmax, img, labels, radii=radii, gauss=0.08 / 0.2,
                                   features_spacing=(0.4, 0.2, 0.25))
        assert len(merge_spy) == 1                        # the merge runs, with today's five arguments
        args, kw = merge_spy[0]
        assert kw == {} and len(args) == 5
        assert args[0][["z", "y", "x"]].to_numpy().tolist() == [[1, 2, 2], [2, 6, 6]]
        assert _same(args[1], spa) and args[2] == 0.8 and args[3] is pre and args[4] == 0.8
        assert per_spot["z_um"].tolist() == [0.4, 0.8] and per_nuc["n_spots"].tolist() == [1, 1]
        assert "stage_flags" not in per_spot.attrs and "stage_flags" not in per_nuc.attrs


# ---- the planar branch -----------------------------------------------------------------------------------

def test_planar_branch_arguments(fake_spotmax, merge_spy, monkeypatch):
    """detect_spots(planar=True) on one plane: radii (PLANAR_RZ_PX, r / dy, r / dx), the Gaussian sigma over
    min(dy, dx), the features spacing (1.0, dy, dx) whatever dz the caller passes, the z-split merge off although
    merge_z_columns is true, and the stage flag z_merge_off_planar on both tables."""
    assert SD.PLANAR_RZ_PX == 1.0
    fake_spotmax.zs = [0, 0]
    img, labels = _two_nuclei(nz=1)
    radii = np.array([1.0, 0.3 / 0.1, 0.3 / 0.125])      # (PLANAR_RZ_PX, r / dy, r / dx)
    np.testing.assert_allclose(radii, [1.0, 3.0, 2.4])
    tables = []
    for dz in (0.2, 5.0, 0.0, float("nan")):              # the dz passed never reaches SpotMAX
        fake_spotmax.clear()
        per_spot, per_nuc = detect_spots(img, labels, (dz, 0.1, 0.125), marker="RAD-51",
                                         cfg_section=_spots_section(), planar=True, **PIPELINE_SPOT_ARGS)
        _check_spotmax_calls(fake_spotmax, img, labels, radii=radii, gauss=0.08 / 0.1,
                             features_spacing=(1.0, 0.1, 0.125))
        assert merge_spy == []                            # merge_z_columns true, and still no merge
        for table in (per_spot, per_nuc):
            assert table.attrs["stage_flags"] == ["z_merge_off_planar"]
        assert [f"spots:{f}" for f in per_spot.attrs["stage_flags"]] == ["spots:z_merge_off_planar"]
        assert per_spot["z_um"].tolist() == [0.0, 0.0]    # every planar z_um is 0.0
        np.testing.assert_allclose(per_spot["y_um"], [0.2, 0.6])
        np.testing.assert_allclose(per_spot["x_um"], [0.25, 0.75])
        assert per_nuc["n_spots"].tolist() == [1, 1] and (per_spot["detector"] == "spotmax").all()
        tables.append((per_spot, per_nuc))
    for per_spot, per_nuc in tables[1:]:
        pd.testing.assert_frame_equal(per_spot, tables[0][0])
        pd.testing.assert_frame_equal(per_nuc, tables[0][1])

    # the pipeline passes (PLANAR_DZ = 1.0, dy, dx): the sigma is over min(dy, dx), never over the placeholder
    # dz, so a coarse plane (pixels above 1 um) is smoothed by its own pixel size
    fake_spotmax.clear()
    detect_spots(img, labels, (1.0, 1.3, 1.2), marker="RAD-51", cfg_section=_spots_section(), planar=True,
                 **PIPELINE_SPOT_ARGS)
    _check_spotmax_calls(fake_spotmax, img, labels, radii=np.array([1.0, 0.3 / 1.3, 0.3 / 1.2]),
                         gauss=0.08 / 1.2, features_spacing=(1.0, 1.3, 1.2))

    # PLANAR_RZ_PX is read when the call runs (test_planar_rz_px_table patches it)
    monkeypatch.setattr(SD, "PLANAR_RZ_PX", 2.0)
    fake_spotmax.clear()
    detect_spots(img, labels, (1.0, 0.1, 0.125), marker="RAD-51", cfg_section=_spots_section(), planar=True,
                 **PIPELINE_SPOT_ARGS)
    _check_spotmax_calls(fake_spotmax, img, labels, radii=np.array([2.0, 0.3 / 0.1, 0.3 / 0.125]),
                         gauss=0.08 / 0.1, features_spacing=(1.0, 0.1, 0.125))
    assert merge_spy == []


def test_planar_branch_keeps_the_optional_keys(fake_spotmax, merge_spy, caplog):
    """A pack recipe's remove_hot_pixels and sharpen apply on one plane as in 3D (with the planar radii);
    spot_radius_z_um has no meaning on one plane, so it is logged and the z radius stays PLANAR_RZ_PX."""
    fake_spotmax.zs = [0, 0]
    img, labels = _two_nuclei(nz=1)
    section = _spots_section(remove_hot_pixels=True, sharpen=True, spot_radius_z_um=0.8)
    with caplog.at_level(logging.WARNING, logger="germquant.spots"):
        per_spot, per_nuc = detect_spots(img, labels, (1.0, 0.1, 0.125), marker="RAD-51", cfg_section=section,
                                         planar=True, **PIPELINE_SPOT_ARGS)
    _check_spotmax_calls(fake_spotmax, img, labels, radii=np.array([1.0, 0.3 / 0.1, 0.3 / 0.125]),
                         gauss=0.08 / 0.1, features_spacing=(1.0, 0.1, 0.125),
                         pre_extra={"do_remove_hot_pixels": True}, sharpened=True)
    assert merge_spy == []
    assert "spot_radius_z_um" in caplog.text and "single plane" in caplog.text
    assert per_spot.attrs["stage_flags"] == ["z_merge_off_planar"] == per_nuc.attrs["stage_flags"]


# ---- which backends receive planar ---------------------------------------------------------------------

def test_planar_forwarded_only_to_planar_ok_backends(monkeypatch):
    """spotmax declares planar_ok, log and net3d do not. planar.check_planar refuses log and net3d in a
    single-plane run before any pixel; detect_spots(planar=True) refuses them too (the second guard).
    planar=True reaches only a backend that declares planar_ok, and planar=False never passes the keyword."""
    assert planar.planar_ok(get_backend("spotmax"))
    assert not planar.planar_ok(get_backend("log")) and not planar.planar_ok(get_backend("net3d"))

    generic = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    roles = {"dna": 0, "foci": 1, "granule": None}
    planar.check_planar(generic, roles, planar=True)      # spotmax, the code default, may run on one plane
    for name in ("log", "net3d"):
        generic.set("spots.backend", name)
        with pytest.raises(planar.PlanarConfigError,
                           match=f"the spot backend '{name}' does not declare planar_ok"):
            planar.check_planar(generic, roles, planar=True)

    img, labels = _two_nuclei(nz=1)
    sp = (1.0, 0.1, 0.1)
    for name in ("log", "net3d"):
        with pytest.raises(planar.PlanarConfigError, match=f"spot backend '{name}' does not declare planar_ok"):
            detect_spots(img, labels, sp, backend=name, planar=True)

    seen: list[tuple[str, dict]] = []

    def spy(name):
        def fn(image, labels, spacing, **kw):
            seen.append((name, kw))
            return SD._empty(labels, kw["marker"], name)
        return fn

    flat, flat_ok, roles_ok = spy("flat_test"), spy("flat_ok_test"), spy("roles_ok_test")
    flat_ok.planar_ok = True
    roles_ok.planar_ok = True
    roles_ok.input_roles = ("dna",)
    for name, fn in (("flat_test", flat), ("flat_ok_test", flat_ok), ("roles_ok_test", roles_ok)):
        monkeypatch.setitem(SB._REGISTRY, name, fn)
    dna = np.zeros_like(img)
    with pytest.raises(planar.PlanarConfigError, match="spot backend 'flat_test' does not declare planar_ok"):
        detect_spots(img, labels, sp, backend="flat_test", planar=True)
    assert seen == []                                     # refused before the backend runs
    for name in ("flat_test", "flat_ok_test", "roles_ok_test"):
        detect_spots(img, labels, sp, backend=name, channels={"dna": dna})
        detect_spots(img, labels, sp, backend=name, channels={"dna": dna}, planar=False)
    assert [(n, sorted(kw)) for n, kw in seen] == [
        ("flat_test", ["cfg_section", "marker"]), ("flat_test", ["cfg_section", "marker"]),
        ("flat_ok_test", ["cfg_section", "marker"]), ("flat_ok_test", ["cfg_section", "marker"]),
        ("roles_ok_test", ["cfg_section", "channels", "marker"]),
        ("roles_ok_test", ["cfg_section", "channels", "marker"])]
    seen.clear()
    detect_spots(img, labels, sp, backend="flat_ok_test", channels={"dna": dna}, planar=True)
    detect_spots(img, labels, sp, backend="roles_ok_test", channels={"dna": dna}, planar=True)
    assert [(n, sorted(kw)) for n, kw in seen] == [
        ("flat_ok_test", ["cfg_section", "marker", "planar"]),
        ("roles_ok_test", ["cfg_section", "channels", "marker", "planar"])]
    assert all(kw["planar"] is True for _, kw in seen) and seen[1][1]["channels"]["dna"] is dna

    # a planar call takes one plane: image and labels (1, Y, X)
    vol, vlab = _two_nuclei(nz=4)
    with pytest.raises(ValueError, match="one plane"):
        detect_spots(vol, vlab, (0.2, 0.1, 0.1), backend="flat_ok_test", planar=True)
    with pytest.raises(ValueError, match="one plane"):
        detect_spots(img, vlab, sp, backend="flat_ok_test", planar=True)

    # the spotmax backend forwards planar=True, and only then, to detect_spots_spotmax
    calls: list[dict] = []

    def spotmax_spy(image, labels, spacing, **kw):
        calls.append(kw)
        return SD._empty(labels, kw["marker"])

    monkeypatch.setattr(SD, "detect_spots_spotmax", spotmax_spy)
    for extra in ({"planar": True}, {"planar": False}, {}):
        detect_spots(img, labels, sp, marker="RAD-51", cfg_section=_spots_section(), **PIPELINE_SPOT_ARGS, **extra)
    assert calls == [{"marker": "RAD-51", **PIPELINE_SPOT_ARGS, "planar": True},
                     {"marker": "RAD-51", **PIPELINE_SPOT_ARGS}, {"marker": "RAD-51", **PIPELINE_SPOT_ARGS}]


# ---- the PLANAR_RZ_PX table (real SpotMAX) ---------------------------------------------------------------

def _planted_plane(seed: int = 2, dy: float = 0.1, dx: float = 0.125, n_per_nucleus: int = 6):
    """One (1, Y, X) plane of 160 x 160 pixels (16 x 20 um): four disk nuclei (radius 2.8 um, nucleoplasm 120
    over a background of 30) with Gaussian puncta (sigma 0.15 um, amplitude 300 to 800) planted at known
    positions at least 0.8 um apart, Gaussian noise of sd 10, uint16. Returns the image, the labels and the
    planted (nucleus_id, y_um, x_um)."""
    rng = np.random.default_rng(seed)
    shape = (160, 160)
    yy, xx = np.indices(shape).astype(np.float64)
    yum, xum = yy * dy, xx * dx
    labels = np.zeros(shape, np.int32)
    img = np.full(shape, 30.0)
    planted = []
    for k, (cy, cx) in enumerate(((4.0, 4.0), (4.0, 12.0), (12.0, 4.0), (12.0, 12.0)), start=1):
        inside = (yum - cy) ** 2 + (xum - cx) ** 2 <= 2.8 ** 2
        labels[inside] = k
        img[inside] = 120.0
        pts: list[tuple[float, float]] = []
        while len(pts) < n_per_nucleus:
            v = rng.uniform(-2.0, 2.0, size=2)
            if np.hypot(*v) > 2.0 or any(np.hypot(cy + v[0] - q[0], cx + v[1] - q[1]) < 0.8 for q in pts):
                continue
            pts.append((cy + v[0], cx + v[1]))
        for py, px in pts:
            img += rng.uniform(300, 800) * np.exp(-((yum - py) ** 2 + (xum - px) ** 2) / (2 * 0.15 ** 2))
            planted.append((k, py, px))
    img += rng.normal(0, 10, size=shape)
    img = np.clip(np.round(img), 0, 65535).astype(np.uint16)
    return img[None], labels[None], pd.DataFrame(planted, columns=["nucleus_id", "y_um", "x_um"])


def _match(per_spot, planted, tol_um: float = 0.3) -> list[tuple[int, int]]:
    """Greedy nearest matching in yx: each planted punctum takes the closest unused spot within tol_um."""
    det = per_spot[["y_um", "x_um"]].to_numpy(dtype=float)
    used = np.zeros(len(det), bool)
    pairs = []
    for i, p in enumerate(planted[["y_um", "x_um"]].to_numpy(dtype=float)):
        if not len(det):
            break
        d = np.hypot(*(det - p).T)
        d[used] = np.inf
        j = int(np.argmin(d))
        if d[j] <= tol_um:
            used[j] = True
            pairs.append((i, j))
    return pairs


def test_planar_rz_px_table(monkeypatch):
    """One synthetic plane through the real SpotMAX chain with PLANAR_RZ_PX 0.5, 1.0 and 2.0: the spot counts
    and the matched F1 against the planted puncta are printed (docs/PLANAR.md carries them), and 1.0 is the
    shipped value (design 4.7: "If they differ, 1.0 stays and the docs say so"). The shipped value recovers the
    planted puncta in their own nuclei, and the planar tables do not depend on the dz the caller passes."""
    pytest.importorskip("spotmax.pipe")
    assert SD.PLANAR_RZ_PX == 1.0                         # the shipped value, before any patch
    img, labels, planted = _planted_plane()
    sp = (1.0, 0.1, 0.125)                                # what the pipeline passes: (PLANAR_DZ, dy, dx)
    tables, rows = {}, []
    for rz in (0.5, 1.0, 2.0):
        with monkeypatch.context() as m:
            m.setattr(SD, "PLANAR_RZ_PX", rz)
            tables[rz] = detect_spots(img, labels, sp, marker="RAD-51", planar=True, **PIPELINE_SPOT_ARGS)
        per_spot = tables[rz][0]
        n_match = len(_match(per_spot, planted))
        f1 = 2.0 * n_match / (len(per_spot) + len(planted))
        rows.append((rz, len(per_spot), n_match, f1))
    assert SD.PLANAR_RZ_PX == 1.0

    def same(a, b) -> bool:
        try:
            pd.testing.assert_frame_equal(a[0], b[0])
            pd.testing.assert_frame_equal(a[1], b[1])
        except AssertionError:
            return False
        return True

    print(f"\nPLANAR_RZ_PX table: one synthetic plane, {len(planted)} planted puncta in 4 nuclei, "
          f"dy 0.1 um, dx 0.125 um, the spots stage's default SpotMAX parameters")
    print("rz_px  spots  matched  F1     tables equal to rz 1.0")
    for rz, n, n_match, f1 in rows:
        equal = "-" if rz == 1.0 else ("yes" if same(tables[rz], tables[1.0]) else "no")
        print(f"{rz:<6} {n:<6} {n_match:<8} {f1:.3f}  {equal}")

    per_spot, per_nuc = tables[1.0]
    pairs = _match(per_spot, planted)
    assert len(pairs) >= 0.9 * len(planted), f"recovered {len(pairs)} of {len(planted)} planted puncta"
    assert len(per_spot) - len(pairs) <= 0.1 * len(planted), f"{len(per_spot) - len(pairs)} extra spots"
    assert all(int(per_spot.at[j, "nucleus_id"]) == int(planted.at[i, "nucleus_id"]) for i, j in pairs)
    assert (per_spot["z_um"] == 0.0).all() and list(per_nuc["nucleus_id"]) == [1, 2, 3, 4]
    assert per_spot.attrs["stage_flags"] == ["z_merge_off_planar"]
    for dz in (0.2, 5.0):                                 # G10 at the detector: the file's dz plays no part
        got = detect_spots(img, labels, (dz, 0.1, 0.125), marker="RAD-51", planar=True, **PIPELINE_SPOT_ARGS)
        assert same(got, tables[1.0]), f"the planar tables moved with dz {dz}"
