"""The SC benchmark harness (germquant.benchmark.sc_ribbons): the Imaris frame rule, zone labels and guard bands by
the protocol's rule (section 4.2), per-gonad AUCs with ties, the pre-specified boundary estimator (5.2), the
threshold choice, the seal on test gonads, and the synthetic matching."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from germquant.benchmark import sc_ribbons as B
from germquant.benchmark import sc_split


def _zt(n=9, step=20.0, f=2, u=3, lt=6, e=8):
    pts = np.stack([np.arange(n) * step, np.zeros(n)], 1)
    return B.ZoneTruth(gonad="g", points_um=pts, names=[B._letter(k) for k in range(n)],
                       idx={"SC_fragstart": f, "SC_fullstart": u, "SC_fraglate": lt, "Pachy_end": e},
                       s_points=np.arange(n) * step)


def test_frame_rule_maps_voxel_centres_to_index_times_nd2_voxel():
    frame = {"ext_min_xy": [100.0, 200.0], "v_imaris_xy": [0.108, 0.108]}
    pos = np.array([[100.0 + 10.5 * 0.108, 200.0 + 3.5 * 0.108]])
    out = B.to_pipeline_xy(pos, frame, (0.108333, 0.108333))
    assert np.allclose(out, [[10 * 0.108333, 3 * 0.108333]])


def test_letters():
    assert [B._letter(k) for k in (0, 1, 25, 26, 27)] == ["A", "B", "Z", "AA", "AB"]
    B.check_letters(["A", "B", "C"])
    with pytest.raises(B.TruthDefect):
        B.check_letters(["A", "C"])


def test_zone_labels_guard_bands_and_cut():
    zt = _zt()
    # points at s = 0, 20, ..., 160; f=2 (B, s 20), u=3 (C, 40), l=6 (F, 100), e=8 (H, 140)
    s = np.array([5, 15, 25, 35, 45, 70, 99, 101, 120, 139, 140, 150, 50, 60])
    r = np.array([1.0] * 12 + [30.0, 1.0])
    nuc = pd.DataFrame({"nucleus_id": np.arange(1, len(s) + 1), "centroid_x_um": s, "centroid_y_um": r})
    lab = B.label_nuclei(zt, nuc)
    z = dict(zip(lab["nucleus_id"], lab["zone"]))
    assert z[1] == "None" and z[2] == "None"          # segment 1 (A-B) is before SC_fragstart (B)
    assert z[3] == "SCfrag" and z[4] == "SCfrag"      # segment 2 (B-C)
    assert z[5] == "SCfull" and z[6] == "SCfull" and z[7] == "SCfull"
    assert z[8] == "LPfrag" and z[9] == "LPfrag" and z[10] == "LPfrag"
    ana = dict(zip(lab["nucleus_id"], lab["in_analysis"]))
    assert ana[11] is True or ana[11] == True        # noqa: E712  s = s(Pachy_end) stays
    assert not ana[12]                                 # beyond Pachy_end
    assert not ana[13]                                 # off axis (r 30 > 20)
    g = dict(zip(lab["nucleus_id"], lab["guard"]))
    # boundary at C (s 40): guard (30, 50); at F (s 100): guard (90, 110); at B (s 20): (10, 30)
    assert g[4] and g[5] and g[7] and g[2] and g[3]
    assert not g[6] and not g[9]


def test_auc_with_ties_and_gonad_aucs():
    assert B.auc([2, 3], [1, 1]) == 1.0
    assert B.auc([1, 2], [1, 2]) == 0.5
    assert B.auc([1], [1]) == 0.5
    zt = _zt(n=12, step=20.0, f=2, u=3, lt=8, e=11)
    rng = np.random.default_rng(0)
    s = rng.uniform(0, 220, 400)
    nuc = pd.DataFrame({"nucleus_id": np.arange(400), "centroid_x_um": s, "centroid_y_um": rng.uniform(0, 3, 400)})
    lab = B.label_nuclei(zt, nuc)
    val = pd.Series(np.where(lab["zone"] == "SCfull", 6.0, 1.0), index=lab["nucleus_id"])
    out = B.gonad_aucs(lab, val)
    assert out["auc_full_vs_frag"] == 1.0


def test_boundary_estimator_on_a_step_profile():
    zt = _zt(n=12, step=20.0, f=2, u=4, lt=9, e=11)
    s = np.arange(0, 200, 1.0)
    nuc = pd.DataFrame({"nucleus_id": np.arange(len(s)), "centroid_x_um": s, "centroid_y_um": 0.5})
    lab = B.label_nuclei(zt, nuc)
    full = (s >= 60) & (s < 160)
    present = s >= 20
    sc = pd.DataFrame({"nucleus_id": np.arange(len(s)), "piece_lw_mean_um": np.where(full, 5.0, 0.5),
                       "total_length_um": np.where(present, 20.0, 0.0)})
    est = B.estimate_boundaries(lab, sc, theta_full=3.0, theta_present=5.0)
    assert abs(est["SC_fullstart"] - 60) <= 5.0
    assert abs(est["SC_fragstart"] - 20) <= 5.0
    assert abs(est["SC_fraglate"] - 160) <= 5.0
    hit, err = B.boundary_hit(zt, "SC_fullstart", est["SC_fullstart"])
    assert hit and err <= 5.0
    th = B.choose_thresholds([(lab, sc)])
    assert th["theta_full"] is not None and th["theta_present"] is not None


def test_truth_loader_refuses_test_gonads():
    split = sc_split.load_split()
    rec = next(g for g in split["gonads"] if g["role"] == "test")
    with pytest.raises(sc_split.SealedError):
        B.zone_truth(rec, "nonexistent.xlsx", (0.108333, 0.108333), frame={"ext_min_xy": [0, 0],
                                                                           "v_imaris_xy": [0.108, 0.108]})


def test_synthetic_matching_of_a_perfect_result():
    from germquant.sc.straighten import synth

    sc = synth.iso_scene(3, L=4.0, peak=500.0, tilt_deg=0.0)
    T = sc.ribbons[0].points[::10]
    samples = pd.DataFrame({"strand_id": 1, "in_strand": True, "s_um": np.arange(len(T)) * 0.05,
                            "z_um": T[:, 0], "y_um": T[:, 1], "x_um": T[:, 2]})
    strands = pd.DataFrame({"strand_id": [1], "length_um": [sc.ribbons[0].length_um], "strand_class": ["complete"]})
    res = {"strands": strands, "samples": samples, "gaps": pd.DataFrame(), "foci": pd.DataFrame()}
    s = B.score_scene(sc, res)
    assert s.coverage == pytest.approx(1.0, abs=1e-3) and s.off_truth_um == pytest.approx(0.0, abs=1e-6)
    assert s.matches[0]["matched"] and s.matches[0]["matchable"] and abs(s.matches[0]["err_um"]) < 1e-6
    summ = B.summarise_scores([s])
    assert summ["complete_precision"] == 1.0
