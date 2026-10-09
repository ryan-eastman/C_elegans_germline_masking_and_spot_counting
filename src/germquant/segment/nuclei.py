"""3D nucleus segmentation.

Default = Cellpose-SAM on GPU (anisotropy from real voxel spacing). When Cellpose/torch
isn't importable (e.g. laptop smoke test), falls back to a classical, fully spacing-aware
watershed so the pipeline still runs end-to-end. The classical path is a baseline, NOT a
substitute for a fine-tuned Cellpose model on crowded pachytene nuclei (see ARCHITECTURE.md).

The segmenter that runs is a pluggable backend (segment/backends.py, docs/BACKENDS.md): the built-in
'cellpose' and 'classical', or any backend an installed package registers under the entry-point group
'germquant.segment_backends'. The config key ``segmentation.nuclei.backend`` (code default 'cellpose')
picks it; the ``segmentation_method`` column reports the one that produced the labels.
"""
from __future__ import annotations

import logging

import numpy as np

from .backends import DEFAULT_BACKEND, get_backend

log = logging.getLogger(__name__)


def segment_nuclei(
    dna: np.ndarray,
    spacing: tuple[float, float, float],
    *,
    method: str = "auto",
    cellpose_model: str = "cpsam",
    diameter_um: float = 3.0,
    min_volume_um3: float = 4.0,
    device: str | None = None,
    backend: str | None = DEFAULT_BACKEND,
    cfg_section: dict | None = None,
) -> tuple[np.ndarray, str]:
    """Return (label_img (Z,Y,X) int32, method_used). `device`: cuda | mps | cpu | auto (None = auto,
    overridden by $GERMQUANT_DEVICE). On an accelerator failure (an operator Metal cannot run, out of
    memory) the same model is retried once on the CPU before the classical fallback is considered.

    `backend` names the segmenter that `method` 'cellpose' and 'auto' run (docs/BACKENDS.md): the
    default 'cellpose' is the Cellpose-SAM path above, 'classical' the watershed, and any other name
    must be registered (an installed package's entry point); an unknown name raises, it never silently
    falls back. method='auto' falls back to the classical watershed when the backend fails,
    method='cellpose' re-raises, method='classical' never calls the backend. `method_used` is the
    backend's name, or 'classical' after a fallback, so the segmentation_method column says which
    segmenter made the labels. `cfg_section` is the segmentation.nuclei config section, handed to the
    backend (as a copy) for keys of its own; the built-ins ignore it.
    """
    if method in ("cellpose", "auto"):
        name = (str(backend).strip() if backend is not None else "") or DEFAULT_BACKEND
        run = get_backend(name)          # a misspelt backend is a config error: raise, never fall back
        what = "Cellpose" if name == DEFAULT_BACKEND else f"segmentation backend {name!r}"
        section = dict(cfg_section or {})   # a copy, so a backend cannot alter the config it was given
        try:
            labels = run(dna, spacing, model=cellpose_model, diameter_um=diameter_um, device=device,
                         cfg_section=section)
            labels = _check_labels(labels, dna.shape, name)
            return _drop_small(labels, spacing, min_volume_um3), name
        except Exception as e:
            if method == "cellpose":
                raise
            log.warning("%s unavailable (%s); using classical watershed fallback.", what, e)
    labels = _classical(dna, spacing, diameter_um)
    return _drop_small(labels, spacing, min_volume_um3), "classical"


def _check_labels(out, shape, name: str) -> np.ndarray:
    """Hold a backend to the contract (docs/BACKENDS.md): an integer (Z, Y, X) label image the size of
    the input. Anything else is the backend's bug and is reported as such, not segmented around."""
    lab = np.asarray(out)
    if lab.shape != tuple(shape):
        raise ValueError(f"segmentation backend {name!r} returned labels of shape {lab.shape}; expected "
                         f"the input shape {tuple(shape)}")
    if lab.dtype.kind not in "iu":
        raise TypeError(f"segmentation backend {name!r} returned {lab.dtype} labels; expected an integer "
                        f"label image (0 = background, one positive label per nucleus)")
    return lab


def _cellpose(dna, spacing, model_name, diameter_um, device_pref: str | None = None) -> np.ndarray:
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

    from ..device import select_device, torch_device

    anisotropy = spacing[0] / spacing[1]          # dz / dy (the key 3D correctness knob)
    is_builtin_sam = str(model_name).lower() in ("cpsam", "sam", "")
    dev = select_device(device_pref)
    # An explicit device covers CUDA, Apple Metal (mps) and CPU alike; Cellpose's own gpu=True would
    # also find CUDA then MPS, but making it explicit keeps the choice in one place (germquant.device),
    # honours $GERMQUANT_DEVICE, and lets MPS run in float32 (bfloat16 is incomplete on Metal).
    kw_model = {"device": torch_device(dev), "gpu": dev != "cpu", "use_bfloat16": dev == "cuda"}
    if is_builtin_sam:
        model = models.CellposeModel(**kw_model)
    else:
        model = models.CellposeModel(pretrained_model=model_name, **kw_model)
    log.info("Cellpose device: %s", dev)
    # cpsam AND models fine-tuned FROM cpsam (our germline_nuclei_* models) are diameter-AGNOSTIC, so
    # run them at NATIVE scale — the regime they were validated in (F1 0.98). Passing `diameter` resizes
    # the image before inference: empirically marginally WORSE on real GT (0.978 vs 0.980) and ~15%
    # different object count. `diameter_um` is kept in the signature/config only for a hypothetical
    # CLASSIC (non-SAM) Cellpose model, which would need it; none of our models do. See diameter_test.py.
    _ = diameter_um
    kw = {}
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
