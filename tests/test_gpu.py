"""GPU-only: the real Cellpose-SAM 3D segmentation path. Skipped unless a CUDA GPU + cellpose
are present (CI / laptops use the classical fallback). Guards against Cellpose API drift — e.g.
Cellpose 4.x requires an explicit z_axis for a 3-D array, which broke the first real GPU run.
"""
import numpy as np
import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("cellpose")
pytestmark = pytest.mark.skipif(not torch.cuda.is_available(), reason="no CUDA GPU available")

from germquant.segment.nuclei import segment_nuclei  # noqa: E402

SPACING = (0.4, 0.2, 0.2)


def _dapi():
    Z, Y, X = 20, 96, 160
    zz, yy, xx = np.indices((Z, Y, X)).astype("float32")
    dapi = np.zeros((Z, Y, X), "float32")
    for cx in (24, 52, 80, 108, 136):
        r2 = (((zz - 10) * SPACING[0]) ** 2 + ((yy - 48) * SPACING[1]) ** 2
              + ((xx - cx) * SPACING[2]) ** 2)
        dapi += np.exp(-r2 / (2 * 1.6 ** 2)) * 1500
    return dapi


def test_cellpose_sam_3d_runs_on_gpu():
    labels, method = segment_nuclei(_dapi(), SPACING, method="cellpose", cellpose_model="cpsam",
                                    diameter_um=3.0, min_volume_um3=1.0)
    assert method == "cellpose"
    assert labels.shape == (20, 96, 160)
    assert labels.max() >= 3, "Cellpose-SAM should segment most of the 5 synthetic nuclei"
