"""The ``nuclei3d`` segmentation backend (EXPERIMENTAL): a trained nuclei3d model's three maps, the h-maxima
seeds of its distance map and a seeded watershed inside its foreground, returned as instance labels under
the contract of segment/backends.py. ``segment_nuclei`` then applies the pipeline's own size filter
(``segmentation.nuclei.min_volume_um3``) and the consecutive relabel, as for every backend.

Selected only by ``segmentation.nuclei.backend: nuclei3d``; never the default and not validated beyond
the development scores on the train and val gonads (docs/BENCHMARK.md, "nuclei3d"). Keys under
``segmentation.nuclei.nuclei3d``, every one with a code default:

| key | code default | meaning |
|---|---|---|
| ``model`` | none (required) | the model file ``germquant train-nuclei3d`` wrote (``model.pt``), absolute or relative to the working directory |
| ``device`` | the stage's ``segmentation.nuclei.device`` | ``auto`` (germquant.device.select_device, honours ``GERMQUANT_DEVICE``), ``cpu`` or ``cuda`` |
| ``tile_zyx`` | [96, 384, 384] | inference tile; the overlap is the network's receptive radius, so the result does not depend on it beyond float rounding |
| ``fg_threshold`` | the model's own, else 0.5 | foreground probability threshold |
| ``h_um`` | the model's own, else 0.4 | h of the h-maxima seeds, in um of the distance map |
| ``neighbour_weight`` | the model's own, else 1.0 | how far the neighbour closeness presses the distance map down at contacts |

The input is the channel the segment stage hands every backend (DAPI unless the profile names another
role). The model's ``model`` argument (``segmentation.nuclei.cellpose_model``) and ``diameter_um`` are not
used. Provenance: the model file's sha256 enters the segment stage sub-hash and the run manifest
(``segment_backend``), through the ``provenance`` attribute of the registered callable.
"""
from __future__ import annotations

import logging

import numpy as np

log = logging.getLogger(__name__)

NAME = "nuclei3d"
NUCLEI3D_DEFAULTS: dict = {
    "model": None,
    "device": None,
    "tile_zyx": [96, 384, 384],
    "fg_threshold": None,
    "h_um": None,
    "neighbour_weight": None,
}
_warned = False


def nuclei3d_kwargs(cfg_section: dict | None) -> dict:
    """The backend's keys from the ``segmentation.nuclei`` section, with their code defaults."""
    sub = (cfg_section or {}).get(NAME, None)
    if sub is None:
        sub = {}
    if not isinstance(sub, dict):
        raise TypeError(f"segmentation.nuclei.{NAME} must be a mapping of the backend's parameters, "
                        f"got {type(sub).__name__}")
    unknown = sorted(set(sub) - set(NUCLEI3D_DEFAULTS))
    if unknown:
        log.warning("segmentation.nuclei.%s: unknown keys %s ignored (known: %s)", NAME, unknown,
                    sorted(NUCLEI3D_DEFAULTS))
    return {k: sub.get(k, d) for k, d in NUCLEI3D_DEFAULTS.items()}


def nuclei3d_provenance(cfg_section: dict | None) -> dict:
    """What a run records before any pixel is read: the model path and its sha256 ("absent" for a
    missing file; never raises)."""
    from germquant import provenance

    sub = (cfg_section or {}).get(NAME) or {}
    model = sub.get("model") if isinstance(sub, dict) else None
    return {"backend": NAME, "model": str(model) if model else None,
            "model_sha256": provenance.file_sha256(model) if model else "absent"}


def resolve_device(key_device, stage_device) -> str:
    """``cuda`` or ``cpu``: the backend's own key, else the stage's device, else auto. Apple's MPS is
    untested for this model, so it runs on the CPU there."""
    want = str(key_device or stage_device or "auto").strip().lower()
    if want == "cpu":
        return "cpu"
    from germquant.device import select_device

    dev = select_device(None if want == "auto" else want)
    return dev if dev in ("cuda", "cpu") else "cpu"


def segment_nuclei3d(dna, spacing, *, device=None, cfg_section: dict | None = None) -> np.ndarray:
    """int32 instance labels (Z, Y, X) of a DAPI volume with the configured model (module docstring)."""
    _warn_experimental()
    from .infer import segment_volume

    kw = nuclei3d_kwargs(cfg_section)
    if not kw["model"]:
        raise ValueError(f"segmentation.nuclei.{NAME}.model must name a trained nuclei3d model file "
                         f"(germquant train-nuclei3d writes MODELDIR/model.pt)")
    vol = np.asarray(dna)
    if vol.ndim != 3:
        raise ValueError(f"the nuclei3d backend needs a (Z, Y, X) volume, got shape {vol.shape}")
    dev = resolve_device(kw["device"], device)
    tile = tuple(None if t is None else int(t) for t in kw["tile_zyx"]) if kw["tile_zyx"] is not None else None
    post = {k: kw[k] for k in ("fg_threshold", "h_um", "neighbour_weight")}
    labels, info = segment_volume(vol, tuple(float(s) for s in spacing), kw["model"], device=dev, tile_zyx=tile,
                                  post=post)
    log.info("nuclei3d: model %s (sha256 %s...), device %s, post %s: %d instances before the size filter "
             "(%.0f s inference, %.0f s instances)", kw["model"], info["model_sha256"][:12], dev, info["post"],
             info["n_instances"], info["infer_s"], info["instances_s"])
    return labels


def _warn_experimental() -> None:
    global _warned
    if not _warned:
        _warned = True
        log.warning("the 'nuclei3d' segmentation backend is EXPERIMENTAL: a 3D U-Net trained on the Imaris nucleus "
                    "surfaces of six N2 gonads; its masks follow the Imaris convention (tighter than Cellpose's) "
                    "and its volumes are not comparable with the cellpose defaults.")
