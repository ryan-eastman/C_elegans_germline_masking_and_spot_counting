"""3D nucleus segmentation.

Default = Cellpose-SAM on GPU (anisotropy from real voxel spacing). When Cellpose/torch
isn't importable (e.g. laptop smoke test), falls back to a classical, fully spacing-aware
watershed so the pipeline still runs end-to-end. The classical path is a baseline, NOT a
substitute for a fine-tuned Cellpose model on crowded pachytene nuclei (see ARCHITECTURE.md).
"""
from __future__ import annotations

import logging

import numpy as np

log = logging.getLogger(__name__)


def segment_nuclei(
    dna: np.ndarray,
    spacing: tuple[float, float, float],
    *,
    method: str = "auto",
    cellpose_model: str = "cpsam",
    diameter_um: float = 3.0,
    min_volume_um3: float = 4.0,
) -> tuple[np.ndarray, str]:
    """Return (label_img (Z,Y,X) int32, method_used)."""
    if method in ("cellpose", "auto"):
        try:
            labels = _cellpose(dna, spacing, cellpose_model, diameter_um)
            return _drop_small(labels, spacing, min_volume_um3), "cellpose"
        except Exception as e:  # noqa: BLE001
            if method == "cellpose":
                raise
            log.warning("Cellpose unavailable (%s); using classical watershed fallback.", e)
    labels = _classical(dna, spacing, diameter_um)
    return _drop_small(labels, spacing, min_volume_um3), "classical"


def _cellpose(dna, spacing, model_name, diameter_um) -> np.ndarray:
    """Run Cellpose 3D. NOTE: the Cellpose API has drifted across 3.x -> SAM(4.x). Verify
    against the version pinned for your GPU before the first real run (this path is untested
    on CPU/CI). Two known sensitivities, handled below:
      * Cellpose-SAM ('cpsam') is diameter-agnostic — `diameter` is ignored for it, so we
        only pass diameter for a non-SAM/fine-tuned model.
      * the built-in SAM model loads via the default CellposeModel; a fine-tuned germline
        model loads via `pretrained_model=<path>`.
    method='auto' swallows any failure here and falls back to the classical watershed; only
    method='cellpose' re-raises (see segment_nuclei).
    """
    from cellpose import models

    anisotropy = spacing[0] / spacing[1]          # dz / dy (the key 3D correctness knob)
    is_sam = str(model_name).lower() in ("cpsam", "sam", "")
    if is_sam:
        model = models.CellposeModel(gpu=True)
        kw = {}                                   # SAM ignores diameter
    else:
        model = models.CellposeModel(gpu=True, pretrained_model=model_name)
        kw = {"diameter": diameter_um / spacing[1]}
    # dna is a single-channel (Z, Y, X) volume: Cellpose 4.x (SAM) requires an explicit z_axis
    # (and no channel axis) for a 3-D array, else it raises "z_axis must be specified".
    out = model.eval(dna, do_3D=True, z_axis=0, channel_axis=None, anisotropy=anisotropy, **kw)
    masks = out[0] if isinstance(out, (list, tuple)) else out
    return np.asarray(masks).astype(np.int32)


def _classical(dna, spacing, diameter_um) -> np.ndarray:
    """Gaussian smooth -> Li threshold -> distance-watershed, all anisotropy-aware."""
    from scipy import ndimage as ndi
    from skimage.feature import peak_local_max
    from skimage.filters import threshold_li
    from skimage.segmentation import watershed

    img = dna.astype(np.float32)
    sigma_um = 0.4
    sigma_vox = tuple(sigma_um / s for s in spacing)
    sm = ndi.gaussian_filter(img, sigma=sigma_vox)

    thr = threshold_li(sm)
    fg = sm > thr
    fg = ndi.binary_fill_holes(fg)  # version-robust (vs skimage remove_small_holes)

    # distance transform in physical units (sampling = spacing)
    dist = ndi.distance_transform_edt(fg, sampling=spacing)

    # seed spacing ~ one nucleus radius, expressed in voxels
    min_dist_vox = max(1, int((diameter_um / 2) / spacing[1]))
    peaks = peak_local_max(
        dist, min_distance=min_dist_vox, labels=fg, exclude_border=False
    )
    markers = np.zeros(fg.shape, dtype=np.int32)
    for i, p in enumerate(peaks, start=1):
        markers[tuple(p)] = i
    if markers.max() == 0:
        markers, _ = ndi.label(fg)

    return watershed(-dist, markers, mask=fg).astype(np.int32)


def _drop_small(labels, spacing, min_volume_um3) -> np.ndarray:
    """Drop objects below a physical volume, then relabel consecutively.

    Implemented directly (bincount) rather than via skimage.remove_small_objects,
    whose `min_size` semantics changed in 0.26 — this is version-robust and never
    merges touching instances.
    """
    from skimage.segmentation import relabel_sequential

    lab = labels.astype(np.int32)
    vox_vol = spacing[0] * spacing[1] * spacing[2]
    min_vox = max(1, int(min_volume_um3 / vox_vol))
    counts = np.bincount(lab.ravel())
    small = np.where(counts < min_vox)[0]
    small = small[small != 0]
    if small.size:
        lab[np.isin(lab, small)] = 0
    lab, _, _ = relabel_sequential(lab)
    return lab.astype(np.int32)
