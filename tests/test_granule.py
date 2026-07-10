"""PGL-1 granule segmentation on planted blobs."""
import numpy as np

from germquant.granule import segment_granules

SP = (0.4, 0.2, 0.2)


def _blobs(shape, centers, half=(1, 2, 2), value=500.0):
    img = np.zeros(shape, np.float32)
    for z, y, x in centers:
        img[z - half[0]:z + half[0] + 1, y - half[1]:y + half[1] + 1, x - half[2]:x + half[2] + 1] = value
    return img


def test_recovers_planted_blobs():
    shape = (8, 40, 40)
    region = np.ones(shape, bool)
    centers = [(4, 10, 10), (4, 10, 30), (4, 30, 20)]
    img = _blobs(shape, centers)

    labels, df = segment_granules(img, region, SP, thresholding_method="otsu",
                                  min_volume_um3=0.05, max_volume_um3=50.0)

    assert labels.max() == 3
    assert len(df) == 3
    assert {"granule_id", "z_um", "y_um", "x_um", "volume_um3", "granule_mean_intensity"} <= set(df.columns)
    assert (df["volume_um3"] > 0).all()


def test_size_gate_drops_out_of_range():
    shape = (8, 40, 40)
    region = np.ones(shape, bool)
    img = np.zeros(shape, np.float32)
    img[1:5, 2:38, 2:38] = 500.0      # big sheet ~ 4*36*36 voxels -> ~83 µm³
    img[6:8, 5:8, 5:8] = 500.0        # small blob ~ 2*3*3 = 18 voxels -> ~0.29 µm³

    labels, df = segment_granules(img, region, SP, thresholding_method="otsu",
                                  min_volume_um3=0.05, max_volume_um3=10.0)

    assert labels.max() == 1          # the big sheet exceeds max_volume and is dropped
    assert len(df) == 1


def test_empty_region_yields_no_granules():
    shape = (6, 20, 20)
    img = _blobs(shape, [(3, 10, 10)])
    labels, df = segment_granules(img, np.zeros(shape, bool), SP)
    assert labels.max() == 0
    assert df.empty


def test_region_restricts_detection():
    shape = (6, 20, 40)
    img = _blobs(shape, [(3, 10, 10), (3, 10, 30)])
    region = np.zeros(shape, bool)
    region[:, :, :20] = True          # only the left blob is inside the region
    labels, df = segment_granules(img, region, SP, thresholding_method="otsu",
                                  min_volume_um3=0.05, max_volume_um3=50.0)
    assert labels.max() == 1
    assert df["x_um"].iloc[0] < 20 * SP[2]
