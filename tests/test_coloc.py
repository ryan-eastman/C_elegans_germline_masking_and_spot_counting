"""Colocalization metrics on tiny planted volumes with KNOWN overlap."""
import numpy as np

from germquant.coloc import colocalize

SP = (0.4, 0.2, 0.2)


def _one_granule(shape, sl):
    labels = np.zeros(shape, np.int32)
    labels[sl] = 1
    return labels


def test_identical_masks_full_overlap():
    shape = (6, 20, 20)
    region = np.ones(shape, bool)
    m = np.zeros(shape, bool)
    m[2:4, 5:10, 5:10] = True
    labels = m.astype(np.int32)
    img = m.astype(np.float32) * 100.0

    row, per_g = colocalize(m, m, labels, img, img, region, SP, n_random=20)

    assert row["dice"] == 1.0
    assert row["jaccard"] == 1.0
    assert row["frac_granules_overlapping_sc"] == 1.0
    assert row["n_granules"] == 1
    assert row["manders_m1"] > 0.999 and row["manders_m2"] > 0.999
    assert row["overlap_volume_um3"] > 0
    assert per_g["overlaps"].all()
    assert (per_g["nearest_um"] == 0).all()          # a granule inside SC has zero distance


def test_disjoint_masks_zero_overlap():
    shape = (6, 20, 30)
    region = np.ones(shape, bool)
    sc = np.zeros(shape, bool)
    sc[2:4, 5:10, 2:7] = True
    gr = np.zeros(shape, bool)
    gr[2:4, 5:10, 22:27] = True
    labels = gr.astype(np.int32)
    sc_img = sc.astype(np.float32) * 100.0
    gr_img = gr.astype(np.float32) * 100.0

    row, per_g = colocalize(sc, gr, labels, sc_img, gr_img, region, SP, n_random=20)

    assert row["dice"] == 0.0
    assert row["overlap_volume_um3"] == 0.0
    assert row["frac_granules_overlapping_sc"] == 0.0
    assert not per_g["overlaps"].any()
    assert row["mean_granule_to_sc_um"] > 0        # separated -> positive nearest-SC distance


def test_partial_overlap_fraction():
    shape = (6, 20, 30)
    region = np.ones(shape, bool)
    gr = np.zeros(shape, bool)
    gr[2:4, 5:10, 5:15] = True        # 2*5*10 = 100 voxels
    sc = np.zeros(shape, bool)
    sc[2:4, 5:10, 5:10] = True        # covers x 5:10 -> half the granule (50 voxels)
    labels = gr.astype(np.int32)
    sc_img = sc.astype(np.float32) * 100.0
    gr_img = gr.astype(np.float32) * 100.0

    row, per_g = colocalize(sc, gr, labels, sc_img, gr_img, region, SP, n_random=50)

    assert abs(row["frac_granule_in_sc"] - 0.5) < 1e-6
    assert abs(per_g["overlap_frac"].iloc[0] - 0.5) < 1e-6
    assert row["frac_granules_overlapping_sc"] == 1.0   # min_frac 0 -> any shared voxel counts
    # SC is fully inside the granule -> all SC intensity coincides with the granule mask
    assert abs(row["manders_m1"] - 1.0) < 1e-6
    assert abs(row["manders_m2"] - 0.5) < 1e-6
    assert 0.0 <= row["overlap_pvalue"] <= 1.0


def test_min_frac_gate():
    shape = (4, 10, 30)
    region = np.ones(shape, bool)
    gr = np.zeros(shape, bool)
    gr[1:3, 2:7, 5:15] = True         # 100 voxels
    sc = np.zeros(shape, bool)
    sc[1:3, 2:7, 5:7] = True          # 20 voxels -> 20% of the granule
    labels = gr.astype(np.int32)
    img = np.zeros(shape, np.float32)

    _, per_lo = colocalize(sc, gr, labels, img, img, region, SP, n_random=0,
                           object_overlap_min_frac=0.1)
    _, per_hi = colocalize(sc, gr, labels, img, img, region, SP, n_random=0,
                           object_overlap_min_frac=0.5)
    assert per_lo["overlaps"].iloc[0]      # 0.2 >= 0.1
    assert not per_hi["overlaps"].iloc[0]  # 0.2 <  0.5
