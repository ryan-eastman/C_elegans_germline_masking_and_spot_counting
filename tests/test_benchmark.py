"""The benchmark harness (germquant.benchmark) on synthetic volumes and spot sets with known answers,
plus both scripts on a temporary manifest. Nothing here needs microscopy data."""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from germquant.benchmark import matching, stats
from germquant.benchmark import nuclei3d as n3
from germquant.benchmark import spots as bs
from germquant.benchmark.manifest import (
    PipelineOutputs,
    Truth,
    find_pipeline_outputs,
    imaris_um_to_voxel_index,
    load_manifest,
    load_truth,
    our_voxel_zyx,
    roi_box_from_params,
    voxel_index_to_um,
)

ROOT = Path(__file__).resolve().parents[1]
VOX = (0.2, 0.1, 0.1)   # dz, dy, dx in um
SHAPE = (8, 60, 60)


# ================================================================================= stats
def test_agreement_known_values():
    x = np.array([1, 2, 3, 4, 5], float)
    y = x + 1
    assert abs(stats.pearson(x, y) - 1.0) < 1e-12
    assert abs(stats.bias(y, x) - 1.0) < 1e-12
    # population variance of 1..5 is 2: ccc = 2*cov / (var + var + bias^2) = 4 / 5
    assert abs(stats.lin_ccc(x, y) - 0.8) < 1e-12
    assert abs(stats.lin_ccc(x, x) - 1.0) < 1e-12
    m = stats.agreement_metrics(x, y)
    assert m["n"] == 5 and abs(m["mean_ours"] - 4.0) < 1e-12 and abs(m["bias"] - 1.0) < 1e-12


def test_agreement_degenerate():
    assert np.isnan(stats.pearson([1, 1, 1], [1, 2, 3]))
    assert np.isnan(stats.lin_ccc([2, 2], [2, 2]))
    assert np.isnan(stats.pearson([1], [2]))


def test_precision_recall_f1():
    s = stats.precision_recall_f1(3, 4, 6)
    assert abs(s["precision"] - 0.75) < 1e-12 and abs(s["recall"] - 0.5) < 1e-12
    assert abs(s["f1"] - 0.6) < 1e-12
    s0 = stats.precision_recall_f1(0, 0, 0)
    assert np.isnan(s0["precision"]) and np.isnan(s0["recall"]) and np.isnan(s0["f1"])
    assert stats.precision_recall_f1(0, 5, 5)["f1"] == 0.0


def test_cluster_bootstrap_contains_point_and_is_deterministic():
    rng = np.random.default_rng(3)
    groups = []
    for g in range(6):
        e = rng.poisson(5, size=40).astype(float)
        groups.append((e, e + rng.normal(0.5 + 0.2 * g, 1.0, size=40)))
    b1 = stats.cluster_bootstrap(groups, stats.pooled_agreement, n_boot=200, seed=0)
    b2 = stats.cluster_bootstrap(groups, stats.pooled_agreement, n_boot=200, seed=0)
    assert b1 == b2
    for k in ("ccc", "pearson_r", "bias"):
        s = b1["stats"][k]
        assert s["lo"] <= s["point"] <= s["hi"], (k, s)
        assert s["n_nan"] == 0
    assert b1["n_groups"] == 6 and b1["n_boot"] == 200
    counts = [(30, 40, 35), (20, 25, 30), (50, 55, 60), (10, 20, 12)]
    bf = stats.cluster_bootstrap(counts, stats.pooled_detection, n_boot=200, seed=0)["stats"]
    for k in ("precision", "recall", "f1"):
        assert bf[k]["lo"] <= bf[k]["point"] <= bf[k]["hi"]
    assert abs(bf["precision"]["point"] - 110 / 140) < 1e-12
    scalar = stats.cluster_bootstrap([np.array([1.0, 2.0]), np.array([3.0])], stats.pooled_mean,
                                     n_boot=50, seed=1)["stats"]["value"]
    assert scalar["lo"] <= scalar["point"] <= scalar["hi"] and abs(scalar["point"] - 2.0) < 1e-12


def test_paired_comparison():
    b = np.array([0.5, 0.6, 0.7, 0.55, 0.65, 0.8])
    a = b + 0.1
    r = stats.paired_comparison(a, b, n_boot=200, seed=0, label_a="new", label_b="old")
    assert r["n_gonads"] == 6 and abs(r["mean_diff"] - 0.1) < 1e-12
    assert r["n_a_better"] == 6 and r["n_b_better"] == 0 and r["n_ties"] == 0
    assert r["diff_lo"] <= r["mean_diff"] <= r["diff_hi"]
    assert np.isfinite(r["wilcoxon_p"]) and r["wilcoxon_p"] < 0.1
    with pytest.raises(ValueError):
        stats.paired_comparison([1, 2], [1])


# ============================================================================== matching
def _random_pairs(seed=1, n=300, n_extra=50):
    rng = np.random.default_rng(seed)
    a = rng.uniform(0, 20, size=(n, 3))
    b = np.vstack([a + rng.normal(0, 0.1, size=a.shape), rng.uniform(0, 20, size=(n_extra, 3))])
    i, j, cost, _, _ = bs.candidate_pairs(a, b, 0.3, 0.4)
    return a, b, i, j, cost


def test_chunked_hungarian_equals_dense():
    a, b, i, j, cost = _random_pairs()
    m = np.full((len(a), len(b)), np.inf)
    m[i, j] = cost
    feasible = np.isfinite(m)
    r, c = matching.assign_dense(np.where(feasible, m, 0.0), feasible)
    dense_n, dense_cost = len(r), float(m[r, c].sum())
    # tiny dense budget: every component with more than a 2 x 2 matrix would go greedy, but the
    # components here are small, so the result must still be exact
    for budget in (4, matching.DEFAULT_MAX_DENSE_CELLS):
        s = matching.assign_sparse(i, j, cost, len(a), len(b), max_dense_cells=budget)
        assert s.n_matched == dense_n
        assert abs(s.total_cost - dense_cost) < 1e-9
        assert s.n_components >= 1
    s2 = matching.assign_sparse(i, j, cost, len(a), len(b))
    assert np.array_equal(s.rows, s2.rows) and np.array_equal(s.cols, s2.cols)   # deterministic


def test_hungarian_beats_greedy_on_the_classic_case():
    # a0 can take b0 (0.5) or b1 (0.9); a1 can only take b0 (0.6). Greedy takes (a0, b0) and leaves
    # a1 alone; the optimum is (a0, b1), (a1, b0).
    rows = np.array([0, 0, 1])
    cols = np.array([0, 1, 0])
    cost = np.array([0.5, 0.9, 0.6])
    g = matching.assign_greedy(rows, cols, cost, 2, 2)
    assert len(g) == 1 and rows[g[0]] == 0 and cols[g[0]] == 0
    h = matching.assign_sparse(rows, cols, cost, 2, 2)
    assert h.n_matched == 2 and h.method == "hungarian"
    assert set(zip(h.rows.tolist(), h.cols.tolist())) == {(0, 1), (1, 0)}
    r, c = matching.assign_dense(np.array([[0.5, 0.9], [0.6, 0.0]]), np.array([[True, True], [True, False]]))
    assert set(zip(r.tolist(), c.tolist())) == {(0, 1), (1, 0)}


def test_greedy_fallbacks_are_flagged():
    rows = np.array([0, 0, 1])
    cols = np.array([0, 1, 0])
    cost = np.array([0.5, 0.9, 0.6])
    over_budget = matching.assign_sparse(rows, cols, cost, 2, 2, max_candidates=2)
    assert over_budget.method == "greedy" and over_budget.n_matched == 1 and over_budget.notes
    comp = matching.assign_sparse(rows, cols, cost, 2, 2, max_dense_cells=1)
    assert comp.method == "hungarian+greedy" and comp.n_greedy_components == 1 and comp.n_matched == 1


def test_assign_sparse_validation_and_dedupe():
    with pytest.raises(ValueError):
        matching.assign_sparse([0], [0], [-1.0], 1, 1)
    with pytest.raises(ValueError):
        matching.assign_sparse([0], [0], [np.nan], 1, 1)
    with pytest.raises(ValueError):
        matching.assign_sparse([0], [3], [1.0], 1, 2)
    d = matching.assign_sparse([0, 0], [0, 0], [0.7, 0.2], 1, 1)
    assert d.n_candidates == 1 and d.n_matched == 1 and abs(d.total_cost - 0.2) < 1e-12
    e = matching.assign_sparse([], [], [], 3, 3)
    assert e.n_matched == 0 and e.method == "none"


# ================================================================================= spots
def test_candidate_pairs_cylinder():
    a = np.array([[1, 1, 1], [5, 5, 5], [8, 8, 8], [12, 12, 12], [2, 2, 20]], float)
    b = np.array([[1.29, 1, 1], [5.31, 5, 5], [8, 8, 8.39], [12, 12, 12.41], [2, 2, 20]], float)
    i, j, cost, dxy, dz = bs.candidate_pairs(a, b, 0.3, 0.4)
    pairs = dict(zip(i.tolist(), j.tolist()))
    assert pairs == {0: 0, 2: 2, 4: 4}
    assert cost[i.tolist().index(4)] == 0.0   # a coincident pair survives with zero cost
    assert abs(dxy[i.tolist().index(0)] - 0.29) < 1e-9 and abs(dz[i.tolist().index(2)] - 0.39) < 1e-9
    empty = bs.candidate_pairs(np.zeros((0, 3)), b, 0.3, 0.4)
    assert empty[0].size == 0


def test_match_spots_scores():
    ours = np.array([[1, 1, 1], [2, 2, 2], [9, 9, 9]], float)
    truth = np.array([[1, 1, 1], [2.1, 2, 2.2], [5, 5, 5], [7, 7, 7]], float)
    a = bs.match_spots(ours, truth)
    assert a.n_matched == 2 and a.method == "hungarian"
    s = stats.precision_recall_f1(a.n_matched, len(ours), len(truth))
    assert abs(s["precision"] - 2 / 3) < 1e-12 and abs(s["recall"] - 0.5) < 1e-12


def test_label_at_points_fallback():
    labels = np.zeros((6, 20, 20), np.int32)
    labels[1:4, 5:10, 5:10] = 4
    zi = [2, 2, 2, 2, 9]
    yi = [7, 4, 2, 7, 7]          # inside; one voxel outside (fallback); three outside (none); in-bounds
    xi = [7, 7, 7, 7, 7]
    zi[3] = 5                     # z too far away (fallback reaches one plane only)
    lab = bs.label_at_points(labels, zi, yi, xi)
    assert lab.tolist() == [4, 4, 0, 0, 0]


def test_is_class_b():
    assert bs.is_class_b(["A", "B", "Class B", "class a", "", None]).tolist() == [False, True, True, False, False, False]


def _synthetic_case():
    """Three of our nuclei (two germline), two traced truth nuclei, three of our spots, four expert spots."""
    labels = np.zeros(SHAPE, np.int32)
    labels[1:5, 5:15, 5:15] = 1
    labels[2:6, 30:40, 30:40] = 2
    labels[3:7, 45:55, 5:15] = 3
    truth_labels = np.zeros(SHAPE, np.int32)
    truth_labels[1:5, 5:15, 5:15] = 11
    truth_labels[2:6, 30:40, 31:41] = 12
    dz, dy, dx = VOX
    nuclei = pd.DataFrame({
        "nucleus_id": [1, 2, 3], "in_germline": [True, True, False], "n_spots": [2, 1, 0],
        "centroid_z_um": [2.5 * dz, 3.5 * dz, 4.5 * dz], "centroid_y_um": [9.5 * dy, 34.5 * dy, 49.5 * dy],
        "centroid_x_um": [9.5 * dx, 34.5 * dx, 9.5 * dx], "voxel_dz_um": dz, "voxel_dy_um": dy, "voxel_dx_um": dx})
    our_vox = [(2, 8, 8), (3, 12, 12), (4, 35, 35)]      # (z, y, x) voxel indices
    spots = pd.DataFrame({
        "spot_id": [0, 1, 2], "nucleus_id": [1, 1, 2],
        "z_um": [z * dz for z, _, _ in our_vox], "y_um": [y * dy for _, y, _ in our_vox],
        "x_um": [x * dx for _, _, x in our_vox], "voxel_dz_um": dz, "voxel_dy_um": dy, "voxel_dx_um": dx})
    exp_vox = [(2, 8, 8), (3, 12, 12), (4, 35, 30), (7, 55, 55)]
    tspots = pd.DataFrame({
        "spot_id": [100, 101, 102, 103],
        "x_um": [x * dx for _, _, x in exp_vox], "y_um": [y * dy for _, y, _ in exp_vox],
        "z_um": [z * dz for z, _, _ in exp_vox],
        "xi": [x for _, _, x in exp_vox], "yi": [y for _, y, _ in exp_vox], "zi": [z for z, _, _ in exp_vox],
        "radius_um": 0.2, "quality": [10.0, 12.0, np.nan, 3.0], "spot_class": ["A", "A", "B", "A"],
        "imaris_x_um": [x * 0.108 for _, _, x in exp_vox], "imaris_y_um": [y * 0.108 for _, y, _ in exp_vox],
        "imaris_z_um": [z * 0.2 for z, _, _ in exp_vox]})
    tnuclei = pd.DataFrame({
        "surface_id": [1, 2], "label": [11, 12],
        "x_um": [9.5 * dx, 35.5 * dx], "y_um": [9.5 * dy, 34.5 * dy], "z_um": [2.5 * dz, 3.5 * dz],
        "xi": [9, 35], "yi": [9, 34], "zi": [2, 3], "volume_um3": [400 * dz * dy * dx] * 2,
        "n_voxels_mask": [400, 400], "imaris_x_um": [1.0, 3.8], "imaris_y_um": [1.0, 3.7], "imaris_z_um": [0.5, 0.7]})
    params = {"voxel_xyz_um_imaris": [0.108, 0.108, 0.2], "voxel_zyx_um_ours": list(VOX),
              "shape_xyz": [SHAPE[2], SHAPE[1], SHAPE[0]], "ext_min_xyz": [0, 0, 0],
              "roi_box": {"min_xyz": [0, 0, 0], "max_xyz": [42, 42, 7], "unit": "voxel"},
              "n_spots": 4, "n_surfaces": 2, "spot_classes": ["A", "B"], "region_growing": False}
    return labels, truth_labels, nuclei, spots, tspots, tnuclei, params


def test_per_nucleus_table_rules():
    labels, truth_labels, nuclei, spots, tspots, tnuclei, _ = _synthetic_case()
    sl = bs.label_at_points(labels, tspots["zi"], tspots["yi"], tspots["xi"])
    assert sl.tolist() == [1, 1, 2, 0]
    ours_n = bs.ours_counts(nuclei, spots)
    tbl, rule = bs.per_nucleus_table(nuclei, sl, ours_n, truth_nuclei=tnuclei, labels=labels,
                                     truth_labels=truth_labels)
    assert rule == "overlap" and tbl["nucleus_id"].tolist() == [1, 2]
    assert tbl["expert_n"].tolist() == [2, 1] and tbl["ours_n"].tolist() == [2, 1]
    assert tbl["traced"].tolist() == [True, True]
    tbl_c, rule_c = bs.per_nucleus_table(nuclei, sl, ours_n, truth_nuclei=tnuclei)
    assert rule_c == "centroid" and tbl_c["traced"].tolist() == [True, True]
    far = tnuclei.copy()
    far["x_um"] += 3.0
    tbl_f, _ = bs.per_nucleus_table(nuclei, sl, ours_n, truth_nuclei=far, surface_tol_um=2.5)
    assert tbl_f["traced"].tolist() == [False, False]
    assert bs.ours_counts(nuclei.drop(columns="n_spots"), spots).tolist() == [2, 1, 0]


def test_benchmark_gonad_spots_known_scores():
    labels, truth_labels, nuclei, spots, tspots, tnuclei, params = _synthetic_case()
    outputs = PipelineOutputs("img", nuclei, spots, None, VOX, "pipeline tables")
    truth = Truth("img", tspots, tnuclei, None, params)
    row, per_nuc, matches = bs.benchmark_gonad(outputs, truth, labels, VOX, tag="g1", split="s",
                                               truth_labels=truth_labels)
    assert row["n_matched"] == 2 and abs(row["precision"] - 2 / 3) < 1e-12 and abs(row["recall"] - 0.5) < 1e-12
    assert row["n_expert_spots_in_our_labels"] == 3 and row["n_expert_spots_in_germline_labels"] == 3
    assert row["n_expert_spots_class_b"] == 1
    assert row["no_b_n_expert"] == 3 and row["no_b_n_ours"] == 3 and row["no_b_n_matched"] == 2
    assert abs(row["no_b_recall"] - 2 / 3) < 1e-12
    assert row["det_n_expert"] == 3 and row["det_n_ours"] == 3 and row["det_n_matched"] == 2
    assert abs(row["ccc_all"] - 1.0) < 1e-12 and abs(row["pearson_r_all"] - 1.0) < 1e-12 and row["bias_all"] == 0
    assert row["n_traced"] == 2 and row["traced_rule"] == "overlap" and row["matching_method"] == "hungarian"
    assert per_nuc["tag"].tolist() == ["g1", "g1"]
    assert sorted(matches["expert_spot_id"].tolist()) == [100, 101]
    # without the truth label image: centroid rule, no detector-only score
    row2, _, _ = bs.benchmark_gonad(outputs, truth, labels, VOX, tag="g1")
    assert row2["traced_rule"] == "centroid" and np.isnan(row2["det_f1"]) and row2["n_traced"] == 2
    pooled = bs.pool(pd.DataFrame([row, row2]), pd.concat([per_nuc, per_nuc]), n_boot=50, seed=0)
    assert pooled["n_gonads"] == 2 and pooled["spots_all"]["n_matched"] == 4
    assert pooled["spots_all"]["f1"]["lo"] <= pooled["spots_all"]["f1"]["point"] <= pooled["spots_all"]["f1"]["hi"]
    assert pooled["spots_detector_only"]["n_gonads"] == 1


# =============================================================================== nuclei3d
def _cube(vol, label, z, y, x):
    vol[z[0]:z[1], y[0]:y[1], x[0]:x[1]] = label


def test_overlap_table_matches_brute_force():
    rng = np.random.default_rng(0)
    pred = rng.integers(0, 4, size=(3, 8, 8))
    truth = rng.integers(0, 3, size=(3, 8, 8))
    ot = n3.overlap_table(pred, truth)
    for pi, ti, inter in zip(ot.pair_pred, ot.pair_truth, ot.pair_inter):
        p, t = ot.pred_ids[pi], ot.truth_ids[ti]
        assert inter == int(((pred == p) & (truth == t)).sum())
    assert len(ot.pair_inter) == len({(int(p), int(t)) for p, t in zip(pred[(pred > 0) & (truth > 0)],
                                                                      truth[(pred > 0) & (truth > 0)])})
    for k, p in enumerate(ot.pred_ids):
        assert ot.pred_sizes[k] == int((pred == p).sum())
        zs = np.flatnonzero((pred == p).any(axis=(1, 2)))
        assert ot.pred_zmin[k] == zs.min() and ot.pred_zmax[k] == zs.max()
        c = np.argwhere(pred == p).mean(axis=0)
        assert np.allclose(ot.pred_centroid_zyx[k], c)


def test_nuclei3d_identical():
    truth = np.zeros((6, 30, 30), np.int32)
    _cube(truth, 5, (1, 4), (2, 10), (2, 10))
    _cube(truth, 9, (2, 6), (15, 25), (15, 25))
    row, per_truth = n3.benchmark_gonad(truth.copy(), truth, VOX, tag="same")
    assert row["tp_iou50"] == 2 and row["fp_iou50"] == 0 and row["fn_iou50"] == 0
    assert row["f1_iou50"] == 1.0 and row["ap_iou50_95"] == 1.0
    assert row["n_splits"] == 0 and row["n_merges"] == 0 and row["z_extent_error_mean_um"] == 0.0
    assert row["n_pred_ignored_outside_coverage"] == 0 and row["mean_matched_iou_iou50"] == 1.0
    assert per_truth["matched_iou50"].all() and per_truth["iou"].tolist() == [1.0, 1.0]


def test_nuclei3d_split_case():
    truth = np.zeros((6, 30, 30), np.int32)
    _cube(truth, 1, (1, 5), (2, 12), (2, 12))          # 4 x 10 x 10 = 400 voxels
    pred = np.zeros_like(truth)
    _cube(pred, 7, (1, 5), (2, 8), (2, 12))            # 240 voxels: IoU 0.6
    _cube(pred, 8, (1, 5), (8, 12), (2, 12))           # 160 voxels: IoU 0.4
    row, per_truth = n3.benchmark_gonad(pred, truth, VOX, tag="split")
    assert row["n_splits"] == 1 and row["n_merges"] == 0
    assert row["tp_iou50"] == 1 and row["fp_iou50"] == 1 and row["fn_iou50"] == 0
    assert per_truth["is_split"].tolist() == [True] and per_truth["n_pred_majority_inside"].tolist() == [2]
    assert per_truth["pred_label"].tolist() == [7] and abs(per_truth["iou"].iloc[0] - 0.6) < 1e-12


def test_nuclei3d_merge_case():
    truth = np.zeros((6, 30, 30), np.int32)
    _cube(truth, 1, (1, 5), (2, 12), (2, 12))          # 400 voxels
    _cube(truth, 2, (1, 5), (12, 18), (2, 12))         # 240 voxels, touching
    pred = np.zeros_like(truth)
    _cube(pred, 3, (1, 5), (2, 18), (2, 12))           # 640 voxels covering both
    row, per_truth = n3.benchmark_gonad(pred, truth, VOX, tag="merge")
    assert row["n_merges"] == 1 and row["n_splits"] == 0
    # IoU with truth 1 is 400/640 = 0.625 (matched), with truth 2 is 240/640 = 0.375
    assert row["tp_iou50"] == 1 and row["fp_iou50"] == 0 and row["fn_iou50"] == 1
    assert per_truth["matched_iou50"].tolist() == [True, False]


def test_nuclei3d_partial_annotation_rule():
    truth = np.zeros(SHAPE, np.int32)
    _cube(truth, 1, (1, 5), (2, 12), (2, 12))
    pred = np.zeros_like(truth)
    _cube(pred, 1, (1, 5), (2, 12), (2, 12))           # the match
    _cube(pred, 2, (1, 5), (20, 26), (2, 8))           # no overlap, 0.8 um from the truth: inside the 2 um band
    _cube(pred, 3, (1, 5), (50, 58), (50, 58))         # no overlap, far outside the dilated bounding box
    row, _ = n3.benchmark_gonad(pred, truth, VOX, tag="partial")
    assert row["coverage_rule"].startswith("truth_bbox_dilated")
    assert row["n_pred_ignored_outside_coverage"] == 1 and row["n_pred_counted"] == 2
    assert row["tp_iou50"] == 1 and row["fp_iou50"] == 1 and row["fn_iou50"] == 0
    assert row["n_pred_no_overlap_in_coverage"] == 1
    # the ROI box as the region: a box that holds only the truth corner ignores both extra predictions
    roi = {"min_xyz": [0, 0, 0], "max_xyz": [15, 15, 7]}
    row_roi, _ = n3.benchmark_gonad(pred, truth, VOX, tag="roi", roi_box=roi)
    assert row_roi["coverage_rule"] == "roi_box" and row_roi["n_pred_ignored_outside_coverage"] == 2
    assert row_roi["fp_iou50"] == 0 and row_roi["f1_iou50"] == 1.0
    # no rule: every unmatched prediction is a false positive
    row_none, _ = n3.benchmark_gonad(pred, truth, VOX, tag="none", coverage="none")
    assert row_none["n_pred_ignored_outside_coverage"] == 0 and row_none["fp_iou50"] == 2
    with pytest.raises(ValueError):
        n3.benchmark_gonad(pred, truth, VOX, coverage="roi")


def test_nuclei3d_z_extent_error_and_ap():
    truth = np.zeros((10, 20, 20), np.int32)
    _cube(truth, 1, (2, 6), (2, 10), (2, 10))          # 4 planes: 0.8 um
    pred = np.zeros_like(truth)
    _cube(pred, 1, (1, 7), (2, 10), (2, 10))           # 6 planes: 1.2 um, IoU 4/6
    row, per_truth = n3.benchmark_gonad(pred, truth, VOX, tag="z")
    assert abs(row["z_extent_error_mean_um"] - 0.4) < 1e-9
    assert abs(per_truth["z_extent_error_um"].iloc[0] - 0.4) < 1e-9
    assert abs(per_truth["truth_z_span_um"].iloc[0] - 0.8) < 1e-9
    # IoU 0.667: matched at 0.5, 0.55, 0.6, 0.65 (TP/(TP+FP+FN) = 1) and not at the six higher
    # thresholds (0 / (0 + 1 + 1) = 0): AP = 4 / 10
    assert abs(row["ap_iou50_95"] - 0.4) < 1e-9
    assert row["f1_iou65"] == 1.0 and row["f1_iou70"] == 0.0


def test_nuclei3d_pool():
    truth = np.zeros((6, 30, 30), np.int32)
    _cube(truth, 5, (1, 4), (2, 10), (2, 10))
    _cube(truth, 9, (2, 6), (15, 25), (15, 25))
    r1, p1 = n3.benchmark_gonad(truth.copy(), truth, VOX, tag="a", split="x")
    pred = truth.copy()
    pred[pred == 9] = 0
    r2, p2 = n3.benchmark_gonad(pred, truth, VOX, tag="b", split="y")
    pooled = n3.pool(pd.DataFrame([r1, r2]), pd.concat([p1, p2]), n_boot=50, seed=0)
    assert pooled["n_gonads"] == 2 and pooled["tp_iou50"] == 3 and pooled["fn_iou50"] == 1
    f1 = pooled["f1_iou50"]
    assert f1["lo"] <= f1["point"] <= f1["hi"] and abs(f1["point"] - 2 * 3 / (3 + 4)) < 1e-12
    assert set(pooled["by_split"]) == {"x", "y"}
    assert pooled["z_extent_error_um"]["point"] == 0.0


# ==================================================================== manifest and scripts
def test_frame_conversion():
    imaris_vox = np.array([0.108, 0.108, 0.2])
    idx = np.array([[1200, 37, 12], [0, 2599, 70]])
    centres = (idx + 0.5) * imaris_vox                  # Imaris reports voxel centres
    assert np.array_equal(imaris_um_to_voxel_index(centres, imaris_vox), idx)
    ours = voxel_index_to_um(idx, np.array([0.108333, 0.108333, 0.2]))
    assert abs(ours[0, 0] - 1200 * 0.108333) < 1e-9
    # matching in raw micrometres would be off by 1200 * (0.108333 - 0.108) = 0.4 um at x = 1200
    assert abs(ours[0, 0] - centres[0, 0]) > 0.3


def test_roi_box_from_params():
    assert roi_box_from_params({"roi_box": {"min_xyz": [1, 2, 3], "max_xyz": [4, 5, 6]}}) == {
        "min_xyz": [1, 2, 3], "max_xyz": [4, 5, 6]}
    assert roi_box_from_params({"roi": [[1, 2, 3], [4, 5, 6]]}) == {"min_xyz": [1, 2, 3], "max_xyz": [4, 5, 6]}
    text = {"surfaces": {"creation_parameters": '<bp mMinXYZ="526 982 0" mMaxXYZ="1482 1442 70"/>'}}
    assert roi_box_from_params(text) == {"min_xyz": [526, 982, 0], "max_xyz": [1482, 1442, 70]}
    assert roi_box_from_params({"n_spots": 3}) is None and roi_box_from_params(None) is None


def _write_case(tmp_path: Path):
    import tifffile

    labels, truth_labels, nuclei, spots, tspots, tnuclei, params = _synthetic_case()
    results = tmp_path / "results" / "gonad_a"
    truth_dir = tmp_path / "truth" / "gonad_a"
    results.mkdir(parents=True)
    truth_dir.mkdir(parents=True)
    image_id = "20250101_N2_nohs_HERM _001"
    nuclei.to_csv(results / f"{image_id}__nuclei.csv", index=False)
    spots.to_csv(results / f"{image_id}__spots.csv", index=False)
    tifffile.imwrite(results / f"{image_id}__nuclei_labels.tif", labels)
    tspots.to_csv(truth_dir / f"{image_id}__imaris_spots.csv", index=False)
    tnuclei.to_csv(truth_dir / f"{image_id}__imaris_nuclei.csv", index=False)
    tifffile.imwrite(truth_dir / f"{image_id}__imaris_labels.tif", truth_labels)
    (truth_dir / f"{image_id}__imaris_params.json").write_text(json.dumps(params), encoding="utf-8")
    entry = {"tag": "case_a", "split": "noHS", "results_dir": str(results), "truth_dir": str(truth_dir),
             "nd2": str(tmp_path / f"{image_id}.nd2"), "ims": str(tmp_path / f"{image_id}.ims")}
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"entries": [entry]}), encoding="utf-8")
    return manifest, entry, image_id


def test_manifest_and_truth_loading(tmp_path):
    manifest, entry, image_id = _write_case(tmp_path)
    entries = load_manifest(manifest)
    assert len(entries) == 1 and entries[0].tag == "case_a" and entries[0].truth_stem == image_id
    as_list = tmp_path / "list.json"
    as_list.write_text(json.dumps([entry]), encoding="utf-8")
    assert load_manifest(as_list)[0].split == "noHS"
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps([{"split": "x"}]), encoding="utf-8")
    with pytest.raises(ValueError):
        load_manifest(bad)
    outputs = find_pipeline_outputs(entry["results_dir"])
    assert outputs.image_id == image_id and outputs.voxel_zyx_um == VOX and outputs.labels_path is not None
    truth = load_truth(entry["truth_dir"], image_id)
    assert truth is not None and truth.has_labels and truth.has_spot_class and len(truth.spots) == 4
    assert truth.spots["spot_class"].tolist() == ["A", "A", "B", "A"]
    assert load_truth(entry["truth_dir"], "nope") is None
    assert our_voxel_zyx(outputs, entries[0], truth) == (VOX, "pipeline tables")
    bare = PipelineOutputs("x", pd.DataFrame(), None, None, None, "")
    assert our_voxel_zyx(bare, None, truth) == (VOX, "truth params (nd2 frame)")


def _run(script: str, *args: str, cwd: Path):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    env["CUDA_VISIBLE_DEVICES"] = "-1"
    return subprocess.run([sys.executable, str(ROOT / "scripts" / script), *args], cwd=cwd, env=env,
                          capture_output=True, text=True)


def test_benchmark_spots_script(tmp_path):
    manifest, _, _ = _write_case(tmp_path)
    out = tmp_path / "out_spots"
    rc = _run("benchmark_spots.py", "--manifest", str(manifest), "--out", str(out), "--n-boot", "50",
              "--write-matches", cwd=tmp_path)
    assert rc.returncode == 0, rc.stdout + rc.stderr
    per_gonad = pd.read_csv(out / "per_gonad.csv")
    assert per_gonad["tag"].tolist() == ["case_a"] and per_gonad["truth_source"].iloc[0] == "truth files"
    assert per_gonad["n_matched"].iloc[0] == 2 and abs(per_gonad["recall"].iloc[0] - 0.5) < 1e-12
    assert per_gonad["det_n_matched"].iloc[0] == 2 and per_gonad["traced_rule"].iloc[0] == "overlap"
    per_nuc = pd.read_csv(out / "per_nucleus.csv")
    assert set(["tag", "split", "nucleus_id", "expert_n", "ours_n", "traced"]) <= set(per_nuc.columns)
    assert per_nuc["expert_n"].tolist() == [2, 1]
    pooled = json.loads((out / "pooled.json").read_text(encoding="utf-8"))
    for k in ("per_nucleus_all", "per_nucleus_traced", "spots_all", "spots_without_class_b",
              "spots_detector_only", "bootstrap", "skipped", "settings", "matching_methods", "truth_sources"):
        assert k in pooled, k
    assert pooled["spots_all"]["n_matched"] == 2 and pooled["skipped"] == []
    assert pooled["per_nucleus_all"]["ccc"]["point"] == 1.0
    matches = pd.read_csv(out / "per_spot_matches.csv")
    assert sorted(matches["expert_spot_id"].tolist()) == [100, 101]
    # a manifest entry without outputs is skipped and reported, not fatal
    man2 = tmp_path / "manifest2.json"
    doc = json.loads(manifest.read_text(encoding="utf-8"))
    doc["entries"].append({"tag": "missing", "split": "HS", "results_dir": str(tmp_path / "nowhere")})
    man2.write_text(json.dumps(doc), encoding="utf-8")
    rc2 = _run("benchmark_spots.py", "--manifest", str(man2), "--out", str(tmp_path / "out2"), "--n-boot", "20",
               cwd=tmp_path)
    assert rc2.returncode == 0, rc2.stdout + rc2.stderr
    pooled2 = json.loads((tmp_path / "out2" / "pooled.json").read_text(encoding="utf-8"))
    assert [s["tag"] for s in pooled2["skipped"]] == ["missing"]


def test_benchmark_nuclei3d_script(tmp_path):
    manifest, _, _ = _write_case(tmp_path)
    out = tmp_path / "out_nuc"
    rc = _run("benchmark_nuclei3d.py", "--manifest", str(manifest), "--out", str(out), "--n-boot", "50", cwd=tmp_path)
    assert rc.returncode == 0, rc.stdout + rc.stderr
    per_gonad = pd.read_csv(out / "per_gonad.csv")
    r = per_gonad.iloc[0]
    assert r["n_pred"] == 3 and r["n_truth"] == 2 and r["coverage_rule"] == "roi_box"
    assert r["n_pred_ignored_outside_coverage"] == 1 and r["tp_iou50"] == 2 and r["fp_iou50"] == 0
    assert r["f1_iou50"] == 1.0
    # truth 12 is shifted one voxel in x against our nucleus 2: IoU 360/440 = 0.818, so the pair is
    # lost at the 0.85, 0.9 and 0.95 thresholds where TP/(TP+FP+FN) = 1/3: AP = (7 + 3 * 1/3) / 10
    assert abs(r["ap_iou50_95"] - 0.8) < 1e-9
    per_nuc = pd.read_csv(out / "per_nucleus.csv")
    assert per_nuc["truth_label"].tolist() == [11, 12] and per_nuc["pred_label"].tolist() == [1, 2]
    assert abs(per_nuc["iou"].iloc[1] - 360 / 440) < 1e-9 and per_nuc["z_extent_error_um"].tolist() == [0.0, 0.0]
    pooled = json.loads((out / "pooled.json").read_text(encoding="utf-8"))
    for k in ("f1_iou50", "ap_iou50_95", "n_splits", "n_merges", "z_extent_error_um", "coverage_rules",
              "bootstrap", "skipped", "settings"):
        assert k in pooled, k
    assert pooled["f1_iou50"]["point"] == 1.0
    # with the bounding-box rule the third nucleus (2 um band reaches it) becomes a false positive
    rc2 = _run("benchmark_nuclei3d.py", "--manifest", str(manifest), "--out", str(tmp_path / "out_bbox"),
               "--coverage", "bbox", "--n-boot", "20", cwd=tmp_path)
    assert rc2.returncode == 0, rc2.stdout + rc2.stderr
    r2 = pd.read_csv(tmp_path / "out_bbox" / "per_gonad.csv").iloc[0]
    assert r2["n_pred_ignored_outside_coverage"] == 0 and r2["fp_iou50"] == 1


def test_spots_script_ims_fallback(tmp_path, monkeypatch):
    """The fallback adapter builds the truth tables from read_imaris_ims with the frame conversion."""
    spec = importlib.util.spec_from_file_location("benchmark_spots_script", ROOT / "scripts" / "benchmark_spots.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    import germquant.validate.imaris_ims as ims_mod

    imaris_vox = np.array([0.108, 0.108, 0.2])
    idx = np.array([[8, 8, 2], [12, 12, 3], [30, 35, 4]])          # x, y, z voxel indices
    centres = (idx + 0.5) * imaris_vox
    fake_spots = pd.DataFrame({"x_um": centres[:, 0], "y_um": centres[:, 1], "z_um": centres[:, 2], "radius_um": 0.2})
    fake_surfaces = pd.DataFrame({"x_um": [1.0], "y_um": [1.0], "z_um": [0.5], "volume_um3": [4.0]})
    geom = {"ext_min_xyz": np.zeros(3), "ext_max_xyz": np.array([60 * 0.108, 60 * 0.108, 8 * 0.2]),
            "shape_xyz": np.array([60, 60, 8]), "voxel_xyz_um": imaris_vox}
    monkeypatch.setattr(ims_mod, "read_imaris_ims", lambda path: (fake_spots, fake_surfaces, geom))
    truth = mod.truth_from_ims(tmp_path / "fake.ims", VOX)
    assert truth.source == "ims_fallback" and truth.stem == "fake" and not truth.has_labels
    assert truth.spots["xi"].tolist() == [8, 12, 30] and truth.spots["zi"].tolist() == [2, 3, 4]
    assert np.allclose(truth.spots["x_um"], idx[:, 0] * VOX[2])
    assert truth.spots["spot_class"].tolist() == ["", "", ""] and truth.spots["quality"].isna().all()
    assert truth.nuclei["label"].isna().all() and len(truth.nuclei) == 1
    assert truth.params["truth_source"] == "ims_fallback"
    labels, _, nuclei, spots, _, _, _ = _synthetic_case()
    outputs = PipelineOutputs("img", nuclei, spots, None, VOX, "pipeline tables")
    row, per_nuc, _ = bs.benchmark_gonad(outputs, truth, labels, VOX, tag="fb")
    assert row["truth_source"] == "ims_fallback" and row["n_matched"] == 2 and row["n_expert_spots"] == 3
    assert row["traced_rule"] == "centroid" and np.isnan(row["no_b_f1"]) and np.isnan(row["det_f1"])
