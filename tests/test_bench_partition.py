"""The protocol's section 5.1 primary partition (docs/BENCHMARK_PROTOCOL.md): the Imaris truth masks are
the nuclei, every spot is assigned by the section 5 rule (the mask label at its centre voxel, else the most
frequent non-zero label of the 3 x 3 x 3 neighbourhood, else unassigned), truth mask fragments below 1000
voxels are removed, and the truth is counted from Class A+B, Class A and Class B spots. Known answers on
tiny label images, the equivalence of the assignment rule with the spot region R (the masks dilated by a
3 x 3 x 3 cube), the independence from the method's own segmentation, the detector-only count, the pooled
blocks, the paired cluster bootstrap of protocol section 6.3, and scripts/bench_detector.py end to end on
the planted-foci gonads of test_bench_detector. Synthetic data only.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import ndimage
from test_bench_detector import planted_bench  # noqa: F401  (the fixture)

from germquant.benchmark import detector as bd
from germquant.benchmark import spots as bs
from germquant.benchmark import stats

ROOT = Path(__file__).resolve().parents[1]
SP = (0.2, 0.108, 0.108)


def _load_script(name: str):
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), ROOT / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------- the assignment rule
def test_assign_to_masks_known_answers():
    lab = np.zeros((5, 7, 7), np.int32)
    lab[2, 2:5, 1:3] = 3          # label 3 on the left
    lab[2, 2:5, 4:6] = 7          # label 7 on the right
    # centre voxel inside a mask: that label, whatever the neighbours
    assert bs.assign_to_masks(lab, [2], [3], [4]).tolist() == [7]
    # centre on background between the two: 3 x 3 x 3 neighbourhood holds three voxels of each, the tie goes
    # to the smaller label
    assert bs.assign_to_masks(lab, [2], [3], [3]).tolist() == [3]
    # one plane above a mask: the neighbourhood reaches it
    assert bs.assign_to_masks(lab, [1], [2], [1]).tolist() == [3]
    # two planes above: unassigned
    assert bs.assign_to_masks(lab, [0], [3], [1]).tolist() == [0]
    # a majority wins over the smaller label
    lab2 = lab.copy()
    lab2[1, 2:5, 4:6] = 7
    assert bs.assign_to_masks(lab2, [2], [3], [3]).tolist() == [7]
    # outside the image counts as background (no error), and an empty query is fine
    assert bs.assign_to_masks(lab, [-1, 2], [3, 3], [1, 99]).tolist() == [0, 0]
    assert bs.assign_to_masks(lab, [], [], []).size == 0


def test_lut_removes_fragments_everywhere():
    lab = np.zeros((3, 40, 40), np.int32)
    lab[:, 5:20, 5:20] = 1                 # 675 voxels
    lab[:, 22:38, 5:38] = 2                # 1584 voxels
    lab[1, 21, 10] = 5                     # a one-voxel fragment between them
    lut, keep, sizes, frag = bs.truth_mask_lut(lab, min_voxels=1000)
    assert keep.tolist() == [2] and sizes.tolist() == [3 * 16 * 33] and frag.tolist() == [1, 5]
    assert lut[1] == 0 and lut[2] == 2 and lut[5] == 0
    lut2, keep2, _, frag2 = bs.truth_mask_lut(lab, min_voxels=1)
    assert keep2.tolist() == [1, 2, 5] and frag2.size == 0
    # a spot on the fragment: with the fragment removed it falls back to the kept neighbour (label 2 at y 22)
    assert bs.assign_to_masks(lab, [1], [21], [10], lut=lut2).tolist() == [5]
    assert bs.assign_to_masks(lab, [1], [21], [10], lut=lut).tolist() == [2]
    # the int32 truth-mask image of the detector-only run holds the kept labels only
    masks = bs.truth_mask_labels(lab, min_voxels=1000)
    assert masks.dtype == np.int32 and set(np.unique(masks).tolist()) == {0, 2}
    assert np.array_equal(masks == 2, lab == 2)


def test_assignment_equals_the_spot_region():
    """Section 5: only spots whose centre voxel lies in R (the masks dilated by one voxel in every axis)
    are counted. With the 3 x 3 x 3 cube that is exactly the set of voxels the rule assigns."""
    rng = np.random.default_rng(5)
    lab = np.zeros((9, 30, 30), np.int32)
    for k in range(1, 6):
        z, y, x = rng.integers(1, 8), rng.integers(3, 27), rng.integers(3, 27)
        lab[max(z - 1, 0):z + 2, y - 2:y + 3, x - 2:x + 3] = k
    region = ndimage.binary_dilation(lab > 0, structure=np.ones((3, 3, 3), bool))
    zz, yy, xx = (a.ravel() for a in np.indices(lab.shape))
    got = bs.assign_to_masks(lab, zz, yy, xx)
    assert np.array_equal(got > 0, region.ravel())
    inside = lab.ravel() > 0
    assert np.array_equal(got[inside], lab.ravel()[inside])        # a voxel of a mask keeps its own label


# ------------------------------------------------------------------------------- the table
def _toy():
    """Two truth nuclei (labels 4 and 9, each 3 x 12 x 12 = 432 voxels) and one fragment (label 2, 4 voxels)."""
    lab = np.zeros((6, 40, 40), np.int32)
    lab[1:4, 4:16, 4:16] = 4
    lab[1:4, 20:32, 20:32] = 9
    lab[5, 38:40, 0:2] = 2
    spots = pd.DataFrame({
        "spot_id": np.arange(7),
        "zi": [2, 2, 2, 2, 0, 5, 3], "yi": [6, 8, 24, 26, 10, 39, 18],
        "xi": [6, 8, 24, 26, 10, 0, 18],
        "spot_class": ["Class A", "Class B", "Class A", "Class A", "Class A", "Class A", "Class B"]})
    spots["x_um"], spots["y_um"], spots["z_um"] = spots["xi"] * SP[2], spots["yi"] * SP[1], spots["zi"] * SP[0]
    return lab, spots


def test_truth_mask_table_counts_and_classes():
    lab, spots = _toy()
    # method spots: two in nucleus 4, one one plane below nucleus 9 (reached by the neighbourhood), one far away
    ours = np.array([[6, 6, 2], [9, 9, 2], [25, 25, 4], [35, 3, 2]], float) * np.array([SP[2], SP[1], SP[0]])
    det = ours[:1]
    table, facts = bs.truth_mask_table(lab, spots, ours, SP, min_voxels=100, detector_only_xyz=det)
    assert table["truth_label"].tolist() == [4, 9] and table["n_voxels"].tolist() == [432, 432]
    # truth: nucleus 4 gets spots 0 (A) and 1 (B); spot 4 sits one plane above it (z 0 next to z 1): A;
    # nucleus 9 gets spots 2 and 3 (A); spot 5 is on the removed fragment, spot 6 lies between the nuclei
    assert table["expert_n"].tolist() == [3, 2]
    assert table["expert_n_class_a"].tolist() == [2, 2] and table["expert_n_class_b"].tolist() == [1, 0]
    assert table["ours_n"].tolist() == [2, 1] and table["ours_n_detector_only"].tolist() == [1, 0]
    assert facts["n_nuclei"] == 2 and facts["n_fragments_removed"] == 1 and facts["fragment_labels"] == [2]
    assert facts["n_expert_spots_assigned"] == 5 and facts["n_expert_spots_unassigned"] == 2
    assert facts["n_our_spots_assigned"] == 3 and facts["n_our_spots_unassigned"] == 1
    assert facts["n_our_spots_detector_only_assigned"] == 1
    # with a lower fragment floor the fragment is a nucleus and spot 5 lands in it
    t2, f2 = bs.truth_mask_table(lab, spots, ours, SP, min_voxels=1)
    assert t2["truth_label"].tolist() == [2, 4, 9] and t2["expert_n"].tolist() == [1, 3, 2]
    assert t2["ours_n_detector_only"].isna().all() and not f2["detector_only"]
    # without a class column the class counts are NaN
    t3, f3 = bs.truth_mask_table(lab, spots.drop(columns="spot_class"), ours, SP, min_voxels=100)
    assert t3["expert_n_class_a"].isna().all() and not f3["has_class"]


def test_truth_mask_row_statistics():
    lab, spots = _toy()
    ours = np.array([[6, 6, 2], [9, 9, 2], [25, 25, 4]], float) * np.array([SP[2], SP[1], SP[0]])
    table, facts = bs.truth_mask_table(lab, spots, ours, SP, min_voxels=100)
    row = bs.truth_mask_row(table, facts)
    m = stats.agreement_metrics(table["expert_n"], table["ours_n"])
    for k in ("mean_expert", "mean_ours", "pearson_r", "ccc", "bias", "mae"):
        assert row[f"tm_{k}"] == pytest.approx(m[k], nan_ok=True), k
    # expert (3, 2), ours (2, 1): bias -1, MAE 1
    assert row["tm_bias"] == pytest.approx(-1.0) and row["tm_mae"] == pytest.approx(1.0)
    ma = stats.agreement_metrics(table["expert_n_class_a"], table["ours_n"])
    assert row["tm_class_a_mae"] == pytest.approx(ma["mae"]) and row["tm_class_a_bias"] == pytest.approx(ma["bias"])
    assert row["tm_n_nuclei"] == 2 and row["tm_available"] and "tm_det_ccc" not in row
    empty = bs.truth_mask_row(None, None)
    assert not empty["tm_available"] and np.isnan(empty["tm_ccc"]) and np.isnan(empty["tm_n_nuclei"])


def test_mean_abs_error_and_pooled_mae():
    assert stats.mean_abs_error([1, 2, 6], [1, 4, 3]) == pytest.approx(5 / 3)
    assert np.isnan(stats.mean_abs_error([], []))
    p = stats.pooled_agreement([(np.array([1.0, 4.0]), np.array([1.0, 2.0])), (np.array([3.0]), np.array([6.0]))])
    assert p["mae"] == pytest.approx(5 / 3)


# ---------------------------------------------------------------- the paired bootstrap
def test_paired_cluster_bootstrap():
    rng = np.random.default_rng(1)
    groups_b = []
    for _ in range(6):
        e = rng.poisson(4, size=50).astype(float)
        groups_b.append((e, e + rng.normal(0, 1.5, size=50)))
    groups_a = [(e, e + rng.normal(0, 0.3, size=len(e))) for e, _ in groups_b]      # a agrees better
    res = stats.paired_cluster_bootstrap(groups_a, groups_b, stats.pooled_agreement, n_boot=300, seed=7)
    s = res["stats"]["ccc"]
    assert s["diff"]["point"] == pytest.approx(s["a"] - s["b"]) and s["diff"]["lo"] > 0 and s["excludes_zero"]
    ca = stats.cluster_bootstrap(groups_a, stats.pooled_agreement, n_boot=300, seed=7)["stats"]["ccc"]
    assert ca["point"] == pytest.approx(s["a"])
    # a method against itself: every difference is zero
    same = stats.paired_cluster_bootstrap(groups_a, groups_a, stats.pooled_agreement, n_boot=100, seed=7)["stats"]
    assert same["ccc"]["diff"]["lo"] == same["ccc"]["diff"]["hi"] == 0.0 and not same["ccc"]["excludes_zero"]
    with pytest.raises(ValueError, match="same gonads"):
        stats.paired_cluster_bootstrap(groups_a, groups_b[:3], stats.pooled_agreement)
    # the per-draw pairing: the difference of the two per-method resample statistics, draw by draw
    d = stats.paired_cluster_bootstrap([1.0, 2.0, 3.0], [0.5, 1.0, 4.0], stats.pooled_mean, n_boot=50, seed=3)
    assert d["stats"]["value"]["diff"]["point"] == pytest.approx(2.0 - 5.5 / 3)


# --------------------------------------------------------- the harness on planted gonads
def test_partition_is_independent_of_our_segmentation(planted_bench):  # noqa: F811
    b = planted_bench
    tag = "synth_noHS_a"
    per_spot, per_nuc = bd.run_backend(tag, "log", {}, cache_root=b["root"] / "cache")
    cache = bd.load_cache(tag, b["root"] / "cache", arrays=("labels",), mmap=False)
    row, _pn, _m, tm = bd.score(tag, per_spot, per_nuc, manifest=b["manifest"], cache=cache, return_truth_masks=True)
    planted = b["gonads"][0]["planted"]
    # the truth masks are the four planted nuclei and the truth count is the planted count
    assert tm["truth_label"].tolist() == [1, 2, 3, 4] and row["tm_n_nuclei"] == 4
    assert tm["expert_n"].tolist() == planted.groupby("nucleus_id").size().tolist()
    assert (tm["expert_n"] - tm["expert_n_class_a"]).tolist() == [1, 1, 1, 1]      # one Class B per nucleus
    assert row["tm_n_expert_spots_unassigned"] == 0 and {"tag", "split", "role"} <= set(tm.columns)
    # a different segmentation of OUR nuclei (germline nucleus 1 missing from our label image): the
    # germline numbers move, the section 5.1 partition does not
    other = bd.load_cache(tag, b["root"] / "cache", arrays=("labels",), mmap=False)
    other.labels = np.where(other.labels == 1, 0, other.labels).astype(np.int32)
    row2, _pn2, _m2, tm2 = bd.score(tag, per_spot, per_nuc, manifest=b["manifest"], cache=other,
                                    return_truth_masks=True)
    pd.testing.assert_frame_equal(tm, tm2)
    for k in [c for c in row if c.startswith("tm_")]:
        assert row[k] == row2[k] or (np.isnan(row[k]) and np.isnan(row2[k])), k
    assert row["mean_expert_all"] != row2["mean_expert_all"]
    # three items unless asked for the table, as before
    assert len(bd.score(tag, per_spot, per_nuc, manifest=b["manifest"], cache=cache)) == 3


def test_bench_detector_reports_the_partition(planted_bench, capsys):  # noqa: F811
    b = planted_bench
    mod = _load_script("bench_detector.py")
    out = b["root"] / "tm_train"
    args = ["--manifest", str(b["manifest"]), "--cache-root", str(b["root"] / "cache"), "--trials-log",
            str(b["root"] / "trials.jsonl"), "--n-boot", "50"]
    rc = mod.main(["--backend", "log", "--split", "train", "--detector-only", "--out", str(out),
                   "--note", json.dumps({"logo": {"fold": 1}}), *args])
    text = capsys.readouterr().out
    assert rc == 0, text
    assert "truth masks (protocol 5.1 primary)" in text and "detector-only count" in text
    per_gonad = pd.read_csv(out / "per_gonad.csv")
    tm = pd.read_csv(out / "per_truth_nucleus.csv")
    assert set(tm["tag"]) == {"synth_noHS_a", "synth_HS_b"} and len(tm) == 8
    assert {"truth_label", "n_voxels", "expert_n", "expert_n_class_a", "expert_n_class_b", "ours_n",
            "ours_n_detector_only", "role", "split"} <= set(tm.columns)
    # in these gonads the truth masks are the cached nuclei, so the detector-only run sees the same label
    # image and must give the same spots, hence the same counts
    assert (tm["ours_n"] == tm["ours_n_detector_only"]).all()
    for t in ("synth_noHS_a", "synth_HS_b"):
        a = pd.read_csv(out / "detections" / f"{t}__spots.csv")
        d = pd.read_csv(out / "detections_detector_only" / f"{t}__spots.csv")
        pd.testing.assert_frame_equal(a, d)
    for _, r in per_gonad.iterrows():
        assert r["tm_n_nuclei"] == 4 and r["tm_ccc"] == pytest.approx(r["tm_det_ccc"])
        sub = tm[tm["tag"] == r["tag"]]
        assert r["tm_ccc"] == pytest.approx(stats.lin_ccc(sub["expert_n"], sub["ours_n"]))
    pooled = json.loads((out / "pooled.json").read_text(encoding="utf-8"))
    for k in bs.TRUTH_MASK_BLOCKS:
        assert k in pooled, k
    blk = pooled["per_nucleus_truth_masks"]
    assert blk["n_nuclei"] == 8 and blk["n_gonads"] == 2
    assert blk["ccc"]["point"] == pytest.approx(stats.lin_ccc(tm["expert_n"], tm["ours_n"]))
    assert blk["mae"]["point"] == pytest.approx(stats.mean_abs_error(tm["ours_n"], tm["expert_n"]))
    assert pooled["per_nucleus_truth_masks_class_b"]["mean_expert"] == pytest.approx(1.0)
    assert "per_nucleus_truth_masks" in pooled["by_split"]["HS"] and pooled["detector_only"] is True
    rec = bd.read_jsonl(b["root"] / "trials.jsonl")[-1]
    assert "per_nucleus_truth_masks" in rec["pooled"] and "per_nucleus_truth_masks_detector_only" in rec["pooled"]
    assert "tm_ccc" in rec["per_gonad"]["synth_HS_b"]
    assert rec["note"] == {"logo": {"fold": 1}} and pooled["note"] == {"logo": {"fold": 1}}
    # without --detector-only: no second run, no detector-only block
    out2 = b["root"] / "tm_val"
    assert mod.main(["--backend", "log", "--split", "val", "--no-trial-log", "--out", str(out2), *args]) == 0
    p2 = json.loads((out2 / "pooled.json").read_text(encoding="utf-8"))
    assert "per_nucleus_truth_masks" in p2 and "per_nucleus_truth_masks_detector_only" not in p2
    assert not (out2 / "detections_detector_only").exists()
