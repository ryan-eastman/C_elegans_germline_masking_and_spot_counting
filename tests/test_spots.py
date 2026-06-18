"""Unit tests for the spot-detection safety helpers (no SpotMAX / GPU needed):
- _cap_candidates: the flood guard that bounds the per-spot feature step.
- _merge_z_columns: the z-axis spot-split merge (PSF duplicates collapse; valley-separated foci stay).
"""
import numpy as np
import pandas as pd

from germquant.spots import spots_to_image
from germquant.spots.detect import _cap_candidates, _merge_z_columns

SP = np.array([0.2, 0.108, 0.108])  # dz, dy, dx (anisotropic, like the real data)


def test_spots_to_image_blobs_at_voxel_locations():
    # spots given in um -> blobs land at the matching voxel, on the requested grid (Imaris-ready)
    df = pd.DataFrame({"z_um": [2 * 0.2, 5 * 0.2], "y_um": [20 * 0.108, 40 * 0.108],
                       "x_um": [20 * 0.108, 50 * 0.108]})
    img = spots_to_image(df, (8, 70, 70), SP, radius_um=0.3)
    assert img.shape == (8, 70, 70) and img.dtype == np.uint16
    assert img[2, 20, 20] == 65535 and img[5, 40, 50] == 65535   # blob centers are bright
    assert img.min() == 0                                        # background stays empty


def test_spots_to_image_empty_is_blank():
    img = spots_to_image(pd.DataFrame(columns=["z_um", "y_um", "x_um"]), (5, 10, 10), SP)
    assert img.shape == (5, 10, 10) and img.max() == 0


def test_spots_to_image_label_mode_unique_ids():
    # far-apart spots -> distinct integer ids (for Imaris label import)
    df = pd.DataFrame({"z_um": [1 * 0.2, 6 * 0.2], "y_um": [5 * 0.108, 50 * 0.108],
                       "x_um": [5 * 0.108, 50 * 0.108]})
    img = spots_to_image(df, (8, 70, 70), SP, radius_um=0.2, label=True)
    assert set(np.unique(img)) - {0} == {1, 2}


def test_cap_candidates_keeps_brightest():
    img = np.zeros((1, 1, 10), np.float32)
    img[0, 0, :] = np.arange(10)                       # brightness increases with x
    df = pd.DataFrame({"z": [0] * 10, "y": [0] * 10, "x": list(range(10)), "nucleus_id": [1] * 10})
    out = _cap_candidates(df, img, cap=3)
    assert len(out) == 3
    assert set(out["x"]) == {7, 8, 9}                  # the 3 brightest survived


def test_cap_candidates_noop_under_cap():
    img = np.zeros((1, 1, 5), np.float32)
    df = pd.DataFrame({"z": [0] * 5, "y": [0] * 5, "x": list(range(5)), "nucleus_id": [1] * 5})
    assert len(_cap_candidates(df, img, cap=10)) == 5  # under the cap -> untouched


def test_merge_z_columns_collapses_psf_stack():
    # 3 peaks at the SAME (x,y), adjacent in z, with a SMOOTH (no-valley) intensity column -> 1 focus
    df = pd.DataFrame({"z": [1, 2, 3], "y": [2, 2, 2], "x": [2, 2, 2],
                       "spot_mean_intensity": [5.0, 6.0, 5.0], "effect_size": [4.0, 4.0, 4.0],
                       "nucleus_id": [1, 1, 1]})
    img = np.zeros((5, 5, 5), np.float32)
    img[1:4, 2, 2] = [5, 6, 5]
    out = _merge_z_columns(df, SP, gap_um=0.8, image=img, valley_frac=0.8)
    assert len(out) == 1


def test_merge_z_columns_keeps_valley_separated_foci():
    # 2 peaks at the same (x,y) but with a deep intensity DIP between them in z -> two distinct foci
    df = pd.DataFrame({"z": [1, 3], "y": [2, 2], "x": [2, 2],
                       "spot_mean_intensity": [10.0, 10.0], "effect_size": [4.0, 4.0],
                       "nucleus_id": [1, 1]})
    img = np.zeros((5, 5, 5), np.float32)
    img[1, 2, 2] = 10; img[2, 2, 2] = 1; img[3, 2, 2] = 10   # valley at z=2
    out = _merge_z_columns(df, SP, gap_um=0.8, image=img, valley_frac=0.8)
    assert len(out) == 2
