"""Guard the cardinal rule: anisotropic voxel spacing must flow into every measurement.
A regression here means every length/volume in a paper is silently wrong.
"""
import numpy as np

from germquant.measure import measure_objects, surface_area_um2


def test_volume_uses_anisotropic_voxel():
    labels = np.zeros((10, 10, 10), dtype=int)
    labels[2:7, 2:7, 2:7] = 1  # 5x5x5 = 125 voxels
    spacing = (0.20, 0.108, 0.108)  # the real N2 voxel
    df = measure_objects(labels, {"dna": labels.astype(float)}, spacing)
    expected = 125 * 0.20 * 0.108 * 0.108
    assert np.isclose(df.loc[0, "volume_um3"], expected)
    # isotropic assumption would give 125 * 0.108^3 — must NOT match
    assert not np.isclose(df.loc[0, "volume_um3"], 125 * 0.108**3)


def test_centroid_scaled_to_microns():
    labels = np.zeros((10, 10, 10), dtype=int)
    labels[4:6, 4:6, 4:6] = 1
    spacing = (0.20, 0.108, 0.108)
    df = measure_objects(labels, {"dna": labels.astype(float)}, spacing)
    # voxel centroid ~4.5 on each axis -> physical = 4.5 * spacing
    assert np.isclose(df.loc[0, "centroid_z_um"], 4.5 * 0.20)
    assert np.isclose(df.loc[0, "centroid_y_um"], 4.5 * 0.108)


def test_surface_area_scales_with_spacing():
    mask = np.zeros((12, 12, 12), dtype=bool)
    mask[3:9, 3:9, 3:9] = True
    a_iso = surface_area_um2(mask, (0.108, 0.108, 0.108))
    a_aniso = surface_area_um2(mask, (0.20, 0.108, 0.108))
    assert a_iso > 0 and a_aniso > 0
    assert not np.isclose(a_iso, a_aniso)  # spacing genuinely changes the result
