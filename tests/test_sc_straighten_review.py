"""Cases from the correctness and protocol reviews of sc_straighten v1 (development seeds and constructed geometry
only): contact ends whose partner is gone, linked tips that overlap, the deeper end refit, tilted gaps, foci near
a strand that only points at them, the strip and foci orientation, a steep ribbon's xy length, the per-nucleus
crop banks, and the scoring path (registered frames, frozen thresholds, sealed runs, guard bands)."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from germquant.benchmark import sc_ribbons as B
from germquant.benchmark import sc_split
from germquant.sc.straighten import foci, imports, inputs, link, params, refine, synth
from germquant.sc.straighten import frames as F
from germquant.sc.straighten.account import resample
from germquant.sc.straighten.follow import RawStrand, free_orphan_contacts
from germquant.sc.straighten.pipeline import straighten_image

P = params.resolve()
ROOT = Path(__file__).resolve().parents[1]


def _seg(a, b, n=60, ends=("free", "free")):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return RawStrand(points=a + np.linspace(0, 1, n)[:, None] * (b - a), end_types=list(ends))


def _run(sc, roles=None, spots=None, p=P):
    return straighten_image(sc.channels["central_element"], sc.labels, sc.spacing, sc.nuclei, params=p,
                            roles=roles, spot_tables=spots, image_id="t")


def _script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------------------------ tracing and linking
def test_contact_end_without_a_partner_is_freed():
    a = _seg((2, 2, 2), (2, 2, 6), ends=("contact", "free"))         # its partner was removed
    b = _seg((2, 5, 2), (2, 5, 6), ends=("free", "contact"))         # stops against c
    c = _seg((2, 4.0, 6.1), (2, 8.0, 6.1))
    assert free_orphan_contacts([a, b, c], P) == 1
    assert a.end_types == ["free", "free"] and b.end_types == ["free", "contact"]


def test_axial_shadow_leaves_no_contact_end():
    # development iso seed 2008: a steep shadow at one end used to trim the ribbon and leave a contact end
    sc = synth.make_scene("iso", 2008)
    res = _run(sc)
    main = res.strands.sort_values("length_um").iloc[-1]
    assert "contact" not in (main["end0_class"], main["end1_class"])
    m = B.score_scene(sc, res).matches[0]
    assert m["matched"] and abs(m["err_um"]) <= 0.2


def test_overlapping_tips_join_without_folding():
    x1, x2 = np.linspace(0.0, 3.0, 31), np.linspace(2.0, 5.0, 31)
    a = RawStrand(np.stack([np.full_like(x1, 2.0), np.full_like(x1, 3.0), x1], 1), ["free", "free"])
    b = RawStrand(np.stack([np.full_like(x2, 2.0), np.full_like(x2, 3.2), x2], 1), ["free", "free"])
    out = link.link_strands([a, b], P)
    assert len(out) == 1
    pts = out[0].points
    assert np.all(np.diff(pts[:, 2]) > 0)                       # monotone along x: no backward step
    # through the production spline (no image): the tangent never reverses, the length stays about 5 um
    P0, _, q = resample(pts, float(P["refine.ds_um"]))
    sig = np.asarray(P["refine.sigma_loc_um"], float)
    W = np.tile(1.0 / sig ** 2, (len(P0), 1))
    spl = refine.pspline(q, P0, W, float(P["refine.knot_um"]), float(P["refine.lambda"]) * float(np.mean(1 / sig ** 2)))
    line = refine.Centreline(spl, q[0], q[-1])
    T = line.tangent(np.linspace(0, line.L, 200))
    assert T[:, 2].min() > 0
    assert abs(line.L - 5.0) <= 0.25


def test_end_far_beyond_the_ribbon_is_refitted_deeper():
    # development iso seed 2000: the traced end runs 1.1 um past the ribbon; the first fit sits at the inner
    # bound of its window, the refit from a deeper anchor finds the edge
    sc = synth.make_scene("iso", 2000)
    res = _run(sc)
    m = B.score_scene(sc, res).matches[0]
    assert m["matched"] and abs(m["err_um"]) <= 0.3


# ------------------------------------------------------------------------------------ gaps
@pytest.mark.parametrize("tilt,tol", [(0.0, 0.10), (60.0, 0.15)])
def test_gap_length_at_tilt(tilt, tol):
    pts = synth.arc(5.0, 3.0, tilt, 1.0)
    sc = synth.custom_scene([{"points": pts, "peak": 450.0, "gaps": [(2.1, 2.9)]}], seed=2600)
    res = _run(sc)
    assert len(res.gaps) == 1
    assert res.gaps["length_um"].iloc[0] == pytest.approx(0.8, abs=tol)
    assert res.gaps["bridged"].isin([True, False]).all() and np.isfinite(res.gaps["plateau_rel"]).all()


# ------------------------------------------------------------------------------------ foci
def _line(points):
    P_ = np.asarray(points, float)
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P_, axis=0), axis=1))])
    W = np.tile(1e6, (len(P_), 3))
    spl = refine.pspline(s, P_, W, 0.3, 1e-6)
    return refine.Centreline(spl, s[0], s[-1])


@pytest.mark.parametrize("gap_um", [0.3, 0.6, 0.8, 1.0, 1.3])
def test_strand_that_only_points_at_a_focus_is_no_second_candidate(gap_um):
    x = np.linspace(0.0, 4.0, 81)
    a = _line(np.stack([np.ones_like(x), np.zeros_like(x), x], 1))
    y = np.linspace(-3.0, -gap_um, 61)
    b = _line(np.stack([np.ones_like(y), y, np.full_like(y, 2.0)], 1))   # ends gap_um short of the focus
    spots = pd.DataFrame({"spot_id": [0], "nucleus_id": [1], "z_um": [1.0], "y_um": [0.04], "x_um": [2.0]})
    rows = foci.assign_foci({"spots": spots}, [F.StrandFrame(a), F.StrandFrame(b)], [a, b], [[], []],
                            (0.365, 0.15), P, 1)
    assert rows[0]["assignment"] == "unique" and rows[0]["strand_key"] == 0


def test_focus_between_two_close_strands_is_ambiguous():
    x = np.linspace(0.0, 4.0, 81)
    a = _line(np.stack([np.ones_like(x), np.zeros_like(x), x], 1))
    b = _line(np.stack([np.ones_like(x), np.full_like(x, 0.09), x], 1))
    spots = pd.DataFrame({"spot_id": [0], "nucleus_id": [1], "z_um": [1.0], "y_um": [0.04], "x_um": [2.0]})
    rows = foci.assign_foci({"spots": spots}, [F.StrandFrame(a), F.StrandFrame(b)], [a, b], [[], []],
                            (0.365, 0.15), P, 1)
    assert rows[0]["assignment"] == "ambiguous"


def test_strip_rows_and_foci_offsets_follow_the_frame():
    """A flat ribbon with a foci-channel blob 0.3 um above it (+z) and one 0.25 um beside it: v and u of the
    foci and the peak rows of the axial and flat strips (rows run from -1.0 um in steps of 0.1 um for v, from
    -0.6 um in steps of 0.05 um for u)."""
    pts = synth.arc(5.0, 6.0, 0.0, 0.3)
    sc = synth.custom_scene([{"points": pts, "peak": 450.0}], seed=2650)
    rib = sc.ribbons[0]
    up = rib.point_at(1.5) + np.array([0.3, 0.0, 0.0])
    t = rib.point_at(3.05) - rib.point_at(2.95)
    t /= np.linalg.norm(t)
    side = np.cross(t, np.array([1.0, 0.0, 0.0]))
    side /= np.linalg.norm(side)
    lat = rib.point_at(3.0) + 0.25 * side
    fch = np.full(sc.labels.shape, synth.OFFSET, float)
    for q in (up, lat):
        synth._render_points(fch, q[None, :], sc.spacing, synth.PSF_SIGMA, weight=400.0)
    spots = pd.DataFrame({"spot_id": [0, 1], "nucleus_id": [1, 1], "z_um": [up[0], lat[0]], "y_um": [up[1], lat[1]],
                          "x_um": [up[2], lat[2]]})
    res = _run(sc, roles={"foci": np.rint(fch).astype(np.uint16)}, spots={"spots": spots})
    st = res.strands.sort_values("length_um").iloc[-1]
    f = res.foci.set_index("spot_id")
    assert f.loc[0, "v_um"] == pytest.approx(0.30, abs=0.05) and abs(f.loc[0, "u_um"]) <= 0.05
    assert abs(f.loc[1, "u_um"]) == pytest.approx(0.25, abs=0.05) and abs(f.loc[1, "v_um"]) <= 0.1
    k = res.strip_roles.index("foci")
    ds = float(P["strips.ds_um"])
    col0 = int(st["strip_offset"] + round((f.loc[0, "s_um"] + 0.5) / ds))
    col1 = int(st["strip_offset"] + round((f.loc[1, "s_um"] + 0.5) / ds))
    assert abs(int(np.nanargmax(res.strips_axial[k][:, col0])) - 13) <= 1
    want = 17 if f.loc[1, "u_um"] > 0 else 7
    assert abs(int(np.nanargmax(res.strips_flat[k][:, col1])) - want) <= 1


def test_xy_length_of_a_steep_ribbon_is_short():
    pts = synth.arc(3.0, 20.0, 80.0, 0.5)
    sc = synth.custom_scene([{"points": pts, "peak": 600.0}], seed=2660)
    st = _run(sc).strands.sort_values("length_um").iloc[-1]
    assert st["length_um"] > 2.0
    assert st["length_xy_um"] <= 0.35 * st["length_um"]


# ------------------------------------------------------------------------------------ inputs
def test_crop_banks_hold_the_same_crops_as_slicing():
    rng = np.random.default_rng(0)
    vol = rng.integers(0, 4000, (9, 30, 40)).astype(np.uint16)
    slices = [(slice(0, 5), slice(2, 12), slice(5, 20)), (slice(3, 9), slice(10, 30), slice(0, 40)),
              (slice(0, 5), slice(2, 12), slice(5, 20))]
    banks = inputs.crop_banks(lambda r, z: vol[z], ["foci"], vol.shape, slices)
    for sl in slices:
        assert np.array_equal(banks["foci"][sl], vol[sl]) and banks["foci"][sl].dtype == vol.dtype


def test_seal_is_case_insensitive_and_covers_the_live_filament_file():
    assert inputs.sealed_role("20260105_N2_HS_HERM_1", "") == "test"
    assert inputs.sealed_role("x", "D:/copies/20260105_N2_HS_HERM_1.ND2") == "test"
    assert inputs.sealed_role("20251105_n2_hs_HERM _001", "") is None              # development
    name = "E:/2025_12_05_DLW193_mch_syp3_male1_nohs image 1_twonucleiSCtraced.IMS"
    with pytest.raises(inputs.InputRefused) as e:
        imports.read_spec("ims", name, image_id="x", spacing=(0.2, 0.1, 0.1))
    assert e.value.code == "sealed_import"


# ------------------------------------------------------------------------------------ the scoring path
def _zone_truth(f, u, lt, e, n=8):
    pts = np.stack([np.arange(n) * 20.0, np.zeros(n)], 1)
    return B.ZoneTruth(gonad="g", points_um=pts, names=[B._letter(k) for k in range(n)],
                       idx={"SC_fragstart": f, "SC_fullstart": u, "SC_fraglate": lt, "Pachy_end": e},
                       s_points=np.arange(n) * 20.0)


def test_boundary_where_fragstart_equals_fullstart_gets_a_guard_band():
    zt = _zone_truth(3, 3, 6, 8)
    nuc = pd.DataFrame({"nucleus_id": np.arange(70), "centroid_x_um": np.arange(70) * 2.0 + 1.0,
                        "centroid_y_um": 0.0})
    lab = B.label_nuclei(zt, nuc)
    kinds = {b_["kind"]: b_ for b_ in lab.attrs["boundaries"]}
    assert set(kinds) == {"None|SCfull", "SCfull|LPfrag"}
    assert (kinds["None|SCfull"]["lo"], kinds["None|SCfull"]["hi"]) == (30.0, 50.0)
    lab2 = B.label_nuclei(zt, nuc, guard_frac=0.25)
    b2 = {b_["kind"]: b_ for b_ in lab2.attrs["boundaries"]}["None|SCfull"]
    assert (b2["lo"], b2["hi"]) == (35.0, 45.0)
    assert lab2["guard"].sum() < lab["guard"].sum()


def test_frame_rules_are_compared_on_the_image_corners():
    a = {"ext_min_xy": [100.0, 200.0], "v_imaris_xy": [0.108, 0.108]}
    assert B.frame_disagreement(a, dict(a), (0.108333, 0.108333), (281.7, 281.7)) == pytest.approx(0.0, abs=1e-9)
    b = {"ext_min_xy": [100.0, 200.0], "v_imaris_xy": [0.1083, 0.1083]}            # same ExtMin, other voxel
    assert B.frame_disagreement(a, b, (0.108333, 0.108333), (281.7, 281.7)) > 0.5


def test_registered_frame_rule_and_unopened_ims(tmp_path, monkeypatch):
    bs = _script("benchmark_sc")
    ims, nd2f = tmp_path / "g.ims", tmp_path / "g.nd2"
    ims.write_bytes(b"x")
    nd2f.write_bytes(b"x")
    opened = []
    monkeypatch.setattr(bs.B, "ims_frame", lambda p: opened.append(p) or {"source": "ims header"})
    monkeypatch.setattr(bs.B, "nd2_frame", lambda p: {"source": "nd2 stage rule"})
    fr = bs._frames({"frame_source": "nd2 stage rule", "ims": str(ims), "nd2": str(nd2f)})
    assert fr["registered"] == "nd2" and "ims" not in fr and not opened
    fr = bs._frames({"frame_source": "ims header", "ims": str(ims), "nd2": str(nd2f)})
    assert fr["registered"] == "ims" and fr["ims"]["source"] == "ims header" and fr["nd2"]


def test_final_scoring_needs_the_frozen_values_and_logs_nothing_when_refused(tmp_path, monkeypatch):
    bs = _script("benchmark_sc")
    log = tmp_path / "final_runs.jsonl"
    monkeypatch.setattr(bs, "FINAL_LOG", log)
    for argv, msg in (
            (["zones", "--final", "--runs", str(tmp_path), "--out", str(tmp_path)], "--reason"),
            (["zones", "--final", "--reason", "r", "--runs", str(tmp_path), "--out", str(tmp_path)], "--theta-full"),
            (["zones", "--final", "--reason", "r", "--runs", str(tmp_path), "--out", str(tmp_path),
              "--theta-full", "2.5", "--theta-present", "7.5", "--results-root", str(tmp_path),
              "--expect-params-hash", "x", "--guard-frac", "0.25"], "guard"),
            (["synthetic", "--final", "--out", str(tmp_path)], "--reason")):
        with pytest.raises(SystemExit) as e:
            bs.main(argv)
        assert msg in str(e.value)
    assert not log.exists()


def test_calibration_refuses_sealed_seeds_before_any_scene(tmp_path):
    cal = _script("calibrate_sc_straighten")
    with pytest.raises(synth.SealedSeedError):
        cal.main(["--out", str(tmp_path), "--families", "iso", "--n", "1", "--seed0", "7100000", "--procs", "1"])
    with pytest.raises(SystemExit):
        cal.main(["--out", str(tmp_path), "--families", "iso", "--n", "1", "--seed0", "7100000", "--final"])
    assert synth.check_seed(7_100_000, final=True) == 7_100_000
    assert not (tmp_path / "calibration.jsonl").exists()


def test_registered_generator_draws_dense_ribbons_without_the_clearance():
    a = synth.make_scene("dense", 11)
    b = synth.make_scene("dense", 11, generator="registered")
    assert b.info["min_separation_um"] == 0.0 and a.info["min_separation_um"] == synth.MIN_RIBBON_SEPARATION_UM
    with pytest.raises(ValueError):
        synth.make_scene("dense", 11, generator="nonsense")


def test_split_record_still_registered():
    assert sc_split.load_split()["sealed"] is True
    rec = json.loads((ROOT / "tests" / "golden" / "sc_split.json").read_bytes().replace(b"\r\n", b"\n"))
    assert rec["counts"]["N2"] == {"dev": 18, "test": 12, "reserved": 9}
