import numpy as np

from germquant.validate import agreement_stats, segmentation_metrics


def test_agreement_perfect():
    x = np.array([0, 1, 2, 3, 4, 5], float)
    s = agreement_stats(x, x.copy())
    assert s["n"] == 6
    assert abs(s["bias_mean_diff"]) < 1e-9
    assert s["mae"] == 0
    assert abs(s["ccc"] - 1.0) < 1e-9


def test_agreement_biased():
    truth = np.array([1, 2, 3, 4, 5], float)
    pred = truth + 2  # constant over-count
    s = agreement_stats(pred, truth)
    assert abs(s["bias_mean_diff"] - 2.0) < 1e-9
    assert s["pearson_r"] > 0.99   # perfectly correlated...
    assert s["ccc"] < 0.9          # ...but CCC penalises the bias


def test_segmentation_identical():
    gt = np.zeros((10, 10, 10), int)
    gt[1:4, 1:4, 1:4] = 1
    gt[6:9, 6:9, 6:9] = 2
    m = segmentation_metrics(gt.copy(), gt, iou_threshold=0.5)
    assert m["tp"] == 2 and m["fp"] == 0 and m["fn"] == 0
    assert m["f1"] == 1.0
    assert m["mean_iou"] == 1.0


def test_segmentation_misses_one():
    gt = np.zeros((10, 10, 10), int)
    gt[1:4, 1:4, 1:4] = 1
    gt[6:9, 6:9, 6:9] = 2
    pred = np.zeros_like(gt)
    pred[1:4, 1:4, 1:4] = 7  # only the first object, different id
    m = segmentation_metrics(pred, gt, iou_threshold=0.5)
    assert m["tp"] == 1 and m["fn"] == 1
    assert abs(m["recall"] - 0.5) < 1e-9
