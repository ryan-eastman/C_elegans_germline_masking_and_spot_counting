"""GPU-only: the real Cellpose-SAM 3D segmentation path. Skipped unless a CUDA GPU + cellpose
are present (CI / laptops use the classical fallback). Guards against Cellpose API drift — e.g.
Cellpose 4.x requires an explicit z_axis for a 3-D array, which broke the first real GPU run.

Markers (registered in pyproject.toml, explained in tests/README.md): ``gpu`` because it needs an
accelerator, ``model`` because it downloads the Cellpose-SAM weights on first use. CI deselects the
``model`` marker. On a machine with a GPU this test DOES run and occupies the card for a few
seconds; set ``CUDA_VISIBLE_DEVICES=-1`` to keep a local run on the CPU (for example while a batch
job owns the GPU), which makes the skipif below trigger (an empty value does not work on Windows,
which drops empty environment variables). ``pytest.importorskip`` turns a missing torch or cellpose
into a skip instead of a collection error.
"""
import numpy as np
import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("cellpose")
_mps = bool(getattr(torch.backends, "mps", None) and torch.backends.mps.is_available())
pytestmark = [
    pytest.mark.gpu,
    pytest.mark.model,
    pytest.mark.skipif(not (torch.cuda.is_available() or _mps),
                       reason="no CUDA or Apple Metal (MPS) accelerator available"),
]

from germquant.segment.nuclei import segment_nuclei

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
