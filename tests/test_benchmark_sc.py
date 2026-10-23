"""The SC benchmark harness (germquant.benchmark.sc_ribbons): the Imaris frame rule, zone labels and guard bands by
the protocol's rule (section 4.2), per-gonad AUCs with ties, the pre-specified boundary estimator (5.2), the
threshold choice, the seal on test gonads, and the synthetic matching."""
from __future__ import annotations

from pathlib import Path

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


def test_filament_check_on_a_synthetic_live_file(tmp_path):
    """Check S3's code path on a phantom written as an Imaris file plus a TIFF (the registered live file is test
    truth and is refused without final=True)."""
    h5py = pytest.importorskip("h5py")
    tifffile = pytest.importorskip("tifffile")
    from germquant.sc.straighten import synth

    sc = synth.dense_scene(5, germ_cell="spermatocyte")
    sp = sc.spacing
    tif = tmp_path / "t.tif"
    tifffile.imwrite(tif, sc.channels["central_element"])
    v_xyz = np.array([sp[2], sp[1], sp[0]])
    ims = tmp_path / "live.ims"
    Z, Y, X = sc.shape

    def s(text):
        return np.frombuffer(str(text).encode(), dtype="S1")

    verts, edges, fils, stats = [], [], [], []
    for k, r in enumerate(sc.ribbons):
        P = np.vstack([r.points[::40], r.points[-1:]])
        vb, eb = len(verts), len(edges)
        for p in P:
            xyz = p[::-1] + 0.5 * v_xyz              # Imaris positions are voxel centres
            verts.append((xyz[0], xyz[1], xyz[2], 0.1))
        edges.extend((vb + i, vb + i + 1) for i in range(len(P) - 1))
        fils.append((0, 1000 + k, vb, vb, len(verts), eb, len(edges)))
        stats.append((0, 1000 + k, 1, r.length_um))
    with h5py.File(ims, "w") as f:
        img = f.create_group("DataSetInfo/Image")
        for i, (n, v) in enumerate(zip((X, Y, Z), v_xyz)):
            img.attrs[f"ExtMin{i}"] = s(0)
            img.attrs[f"ExtMax{i}"] = s(f"{n * v:.9g}")
        img.attrs["X"], img.attrs["Y"], img.attrs["Z"] = s(X), s(Y), s(Z)
        f.create_group("DataSet/ResolutionLevel 0/TimePoint 0/Channel 0").create_dataset(
            "Data", data=np.zeros((2, 2, 2), np.uint16))
        g = f.create_group("Scene8/Content/Filaments0")
        g.create_dataset("Vertex", data=np.array(verts, dtype=[("PositionX", "<f8"), ("PositionY", "<f8"),
                                                               ("PositionZ", "<f8"), ("Radius", "<f4")]))
        g.create_dataset("Edge", data=np.array(edges, dtype=[("VertexA", "<i8"), ("VertexB", "<i8")]))
        g.create_dataset("Filament", data=np.array(fils, dtype=[("ID_Time", "<i8"), ("ID", "<i8"),
                                                                ("IndexVertexRoot", "<i8"), ("IndexVertexBegin", "<i8"),
                                                                ("IndexVertexEnd", "<i8"), ("IndexEdgeBegin", "<i8"),
                                                                ("IndexEdgeEnd", "<i8")]))
        g.create_dataset("StatisticsType", data=np.array([(1, 0, 0, b"Filament Length (sum)", b"um")],
                                                         dtype=[("ID", "<i8"), ("ID_Category", "<i8"),
                                                                ("ID_FactorList", "<i8"), ("Name", "S64"),
                                                                ("Unit", "S8")]))
        g.create_dataset("StatisticsValue", data=np.array(stats, dtype=[("ID_Time", "<i8"), ("ID_Object", "<i8"),
                                                                        ("ID_StatisticsType", "<i8"),
                                                                        ("Value", "<f8")]))
    out = B.filament_check(ims, tif, sp, time_index=0)
    assert out["n_truth_paths"] == len(sc.ribbons) and out["n_nuclei"] == 1
    for rib in out["ribbons"]:
        assert abs(rib["drawn_um"] - rib["stored_um"]) / rib["stored_um"] < 0.02       # S3a on exact truth
    assert out["nuclei"][0]["recall"] > 0.7
    # the registered live file name is refused without final=True
    reg = tmp_path / Path(sc_split.load_split()["filament_truth"]["ims"]["path"]).name
    reg.write_bytes(ims.read_bytes())
    with pytest.raises(sc_split.SealedError):
        B.filament_check(reg, tif, sp, time_index=0)
