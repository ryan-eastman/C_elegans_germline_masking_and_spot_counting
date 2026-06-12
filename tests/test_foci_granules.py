"""Foci neighborhood-intensity sampling + assignment, and the granule (generic 3D-object) stage."""
import numpy as np

from germquant.foci import detect_foci
from germquant.granules import detect_granules

SPACING = (0.4, 0.15, 0.15)


def test_foci_neighborhood_intensity_and_assignment():
    Z, Y, X = 16, 40, 40
    img = np.zeros((Z, Y, X), np.float32)
    zz, yy, xx = np.indices((Z, Y, X)).astype(np.float32)
    # one peaked gaussian spot
    r2 = (((zz - 8) * SPACING[0]) ** 2 + ((yy - 20) * SPACING[1]) ** 2 + ((xx - 20) * SPACING[2]) ** 2)
    img += np.exp(-r2 / (2 * 0.25 ** 2)) * 5000.0
    labels = np.zeros((Z, Y, X), np.int32)
    labels[4:12, 12:28, 12:28] = 7  # a nucleus around the spot

    foci = detect_foci(img, SPACING, labels=labels, min_sigma_um=0.1, max_sigma_um=0.5, threshold_rel=0.2)
    assert len(foci) >= 1
    f = foci.iloc[0]
    # neighborhood mean is pulled below the single-voxel peak (B10: not mean==max-at-center)
    assert f["intensity_max"] >= f["intensity_mean"]
    assert f["intensity_mean"] < img.max()
    assert int(f["nucleus_id"]) == 7
    assert f["detection_radius_um"] > 0


def test_detect_granules_shape_and_assignment():
    Z, Y, X = 16, 40, 40
    img = np.zeros((Z, Y, X), np.float32)
    img[6:10, 8:14, 8:14] = 3000.0     # granule A
    img[6:10, 26:32, 26:32] = 3000.0   # granule B
    labels = np.zeros((Z, Y, X), np.int32)
    labels[4:12, 4:18, 4:18] = 5       # nucleus over granule A only

    g = detect_granules(img, SPACING, labels=labels, min_volume_um3=0.01)
    assert len(g) == 2
    assert (g["volume_um3"] > 0).all()
    assert g["surface_area_um2"].notna().all()
    assert g["sphericity"].notna().all()
    # the granule inside the nucleus is assigned to it; the other is unassigned (0)
    assert set(g["nucleus_id"]) == {0, 5}
