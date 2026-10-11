"""The ``net3d`` spot backend (EXPERIMENTAL): a trained net3d model's heatmap, peaks and nucleus
assignment, written into the spot tables of the backend contract (spots/backends.py).

Selected only by ``spots.backend: net3d``; never the default and not validated against Imaris beyond
the development scores in docs/BENCHMARK.md. Keys under ``spots.net3d``, every one with a code default:

| key | code default | meaning |
|---|---|---|
| ``model`` | none (required) | the model file ``germquant train-spots`` wrote (``model.pt``), absolute or relative to the working directory |
| ``threshold`` | the model's own | peak probability threshold; by default the one chosen on the validation gonads and stored in the model file (0.5 when the file has none) |
| ``device`` | ``auto`` | ``auto`` (germquant.device.select_device, honours ``GERMQUANT_DEVICE``), ``cpu`` or ``cuda`` |
| ``tile_zyx`` | [96, 384, 384] | inference tile; the overlap is the network's receptive radius, so the result does not depend on it beyond float rounding |
| ``separation_um`` | the model's own, else [0.6, 0.3] | (z, xy) radii of the suppression ellipsoid |
| ``assign_max_um`` | the model's own, else 0.3 | reach of the nearest-label fallback |

Shared key read from the top of ``spots:``: ``max_spot_candidates`` (30000), the flood guard; the
most probable peaks are kept when it is exceeded, with the SpotMAX path's loud warning. The z-split
merge keys do not apply (the anisotropic suppression replaces them).

Inputs: the foci channel the contract hands every backend and, for a model trained with DAPI, the DNA
channel through the contract's optional ``channels`` mapping (the backend declares ``input_roles =
("dna",)``; the pipeline and the bench harness pass that channel when the image has it). A model that
needs DAPI on an image without it raises a clear error, which the pipeline records as a failed spots
stage.

Columns: ``spot_mean_intensity`` is the mean raw intensity of the core voxels (radius 0.25 um) and
``effect_size`` the log backend's Glass's delta against the local annulus (0.5 to 1.0 um); after the
declared columns come ``prob`` (the heatmap peak), ``assignment`` (``label`` or ``nearest``) and
``contested`` (two or more nuclei within ``assign_max_um``). Provenance: the model file's sha256 enters
the spots stage sub-hash and the run manifest (``spot_backend``).
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from . import peaks as pk

log = logging.getLogger(__name__)

DETECTOR = "net3d"
NET3D_DEFAULTS: dict = {
    "model": None,
    "threshold": None,
    "device": "auto",
    "tile_zyx": [96, 384, 384],
    "separation_um": None,
    "assign_max_um": None,
}
SHARED_DEFAULTS: dict = {"max_spot_candidates": 30000}
MEASURE_DIAMETER_UM = 0.5          # core and annulus of the effect size, as the log backend's defaults
_warned = False


def net3d_kwargs(cfg_section: dict) -> dict:
    """Keyword arguments of ``detect_spots_net3d`` from the raw ``spots:`` section."""
    sub = cfg_section.get("net3d", None)
    if sub is None:
        sub = {}
    if not isinstance(sub, dict):
        raise TypeError(f"spots.net3d must be a mapping of the net3d backend's parameters, got {type(sub).__name__}")
    unknown = sorted(set(sub) - set(NET3D_DEFAULTS))
    if unknown:
        log.warning("spots.net3d: unknown keys %s ignored (known: %s)", unknown, sorted(NET3D_DEFAULTS))
    kw = {k: sub.get(k, d) for k, d in NET3D_DEFAULTS.items()}
    kw.update({k: cfg_section.get(k, d) for k, d in SHARED_DEFAULTS.items()})
    return kw


def net3d_provenance(cfg_section: dict) -> dict:
    """What a run records for this backend before any pixel is read: the model path and its sha256
    (the file digest the spots stage sub-hash includes). Never raises for a missing file ("absent")."""
    from germquant import provenance

    sub = cfg_section.get("net3d") or {}
    model = sub.get("model") if isinstance(sub, dict) else None
    return {"backend": DETECTOR, "model": str(model) if model else None,
            "model_sha256": provenance.file_sha256(model) if model else "absent"}


def resolve_device(device: str | None) -> str:
    want = str(device or "auto").strip().lower()
    if want == "cpu":
        return "cpu"
    from germquant.device import select_device

    dev = select_device(None if want == "auto" else want)
    return dev if dev in ("cuda", "cpu") else "cpu"      # MPS is untested for this model: CPU there


def model_inputs(meta: dict, image, channels: dict | None) -> list:
    """The raw input volumes in the model's order (``input_roles`` of the model file)."""
    roles = list(meta.get("input_roles") or ["foci", "dna"])
    chans = dict(channels or {})
    out = []
    for r in roles:
        if r == "foci":
            out.append(image)
        elif r in chans and chans[r] is not None:
            out.append(chans[r])
        else:
            raise ValueError(f"the net3d model was trained on the {roles} channels but the run provides no "
                             f"{r!r} channel (the image has no such role, or the caller did not pass it)")
    return out


def spots_from_peaks(found: pd.DataFrame, image, labels, spacing, marker: str):
    """The two tables of the contract from ``peaks.find_spots`` output: unassigned peaks dropped, the
    core intensity and effect size measured on the raw foci channel, the extras appended."""
    from ..detect import PER_SPOT_COLS, _empty, _finish_tables
    from ..log_detect import _ellipsoid_offsets, _measure, _resolve, _shell_offsets

    lab = np.asarray(labels)
    sp = np.asarray(spacing, float)
    f = found[found["nucleus_id"] > 0].reset_index(drop=True)
    if f.empty:
        per_spot, per_nuc = _empty(lab, marker, DETECTOR)
        for c, dt in (("prob", float), ("assignment", object), ("contested", bool)):
            per_spot[c] = pd.Series(dtype=dt)
        return per_spot, per_nuc
    p = _resolve(MEASURE_DIAMETER_UM, 0.0, None, None, None)
    cand = pd.DataFrame({"nucleus_id": f["nucleus_id"].astype(int).to_numpy(), "z": f["z"].to_numpy(int),
                         "y": f["y"].to_numpy(int), "x": f["x"].to_numpy(int)})
    meas = _measure(cand, np.asarray(image), lab, sp, _ellipsoid_offsets(p["core_um"], sp),
                    _shell_offsets(p["annulus_inner_um"], p["annulus_outer_um"], sp), p["core_um"],
                    p["annulus_outer_um"])
    out = pd.DataFrame({
        "nucleus_id": f["nucleus_id"].astype(int).to_numpy(),
        "z": f["zr"].to_numpy(float), "y": f["yr"].to_numpy(float), "x": f["xr"].to_numpy(float),
        "effect_size": meas["effect_size"].to_numpy(float),
        "spot_mean_intensity": meas["spot_mean_intensity"].to_numpy(float),
        "prob": f["prob"].to_numpy(float), "assignment": f["rule"].astype(str).to_numpy(),
        "contested": f["contested"].astype(bool).to_numpy(),
    })
    out = out.sort_values(["nucleus_id", "z", "y", "x"], kind="stable").reset_index(drop=True)
    per_spot, per_nuc = _finish_tables(out, lab, sp, marker, DETECTOR)
    per_spot = per_spot[PER_SPOT_COLS].copy()
    for c in ("prob", "assignment", "contested"):
        per_spot[c] = out[c].to_numpy()
    return per_spot, per_nuc


def detect_spots_net3d(spot_img, labels, spacing, *, marker: str = "RAD-51", channels: dict | None = None,
                       model=None, threshold=None, device: str = "auto", tile_zyx=(96, 384, 384),
                       separation_um=None, assign_max_um=None, max_spot_candidates: int = 30000):
    """Return (per_spot, per_nucleus) with the contract's columns (plus the extras of the module
    docstring) and detector ``net3d``."""
    _warn_experimental()
    from .infer import channel_ranges, load_model, predict_heatmap

    if not model:
        raise ValueError("spots.net3d.model must name a trained net3d model file (germquant train-spots writes "
                         "MODELDIR/model.pt)")
    img = np.asarray(spot_img)
    lab = np.asarray(labels)
    if img.ndim != 3 or lab.shape != img.shape:
        raise ValueError("the net3d backend needs a (Z, Y, X) foci image and a label image of the same shape, "
                         f"got {img.shape} and {lab.shape}")
    sp = tuple(float(s) for s in spacing)
    dev = resolve_device(device)
    net, meta, sha = load_model(model, dev)
    pk_meta = meta.get("peaks") or {}
    thr = float(threshold if threshold is not None else pk_meta.get("threshold", 0.5))
    sep = tuple(separation_um if separation_um is not None else pk_meta.get("separation_um", pk.DEFAULT_SEPARATION_UM))
    reach = float(assign_max_um if assign_max_um is not None else pk_meta.get("assign_max_um", pk.DEFAULT_ASSIGN_MAX_UM))
    inputs = model_inputs(meta, img, channels)
    ranges = channel_ranges(inputs, meta.get("norm"))
    tile = tuple(None if t is None else int(t) for t in tile_zyx) if tile_zyx is not None else None
    heat = predict_heatmap(net, inputs, ranges, device=dev, tile_zyx=tile, labels=lab)
    found = pk.find_spots(heat, lab, sp, thr, separation_um=sep, assign_max_um=reach)
    del heat
    n_near = len(found)
    found = found[found["nucleus_id"] > 0].reset_index(drop=True)
    n_cap = int(max_spot_candidates)
    if len(found) > n_cap:
        log.warning("FLOODED: %d net3d spots exceed cap %d; keeping the %d most probable, the count is a FLOOR "
                    "not a measurement.", len(found), n_cap, n_cap)
        found = found.iloc[:n_cap]                    # find_spots returns the most probable first
    log.info("net3d: model %s (sha256 %s...), device %s, threshold %.3f: %d peaks near a nucleus, %d assigned",
             model, sha[:12], dev, thr, n_near, len(found))
    return spots_from_peaks(found, img, lab, sp, marker)


def _warn_experimental() -> None:
    global _warned
    if not _warned:
        _warned = True
        log.warning("the 'net3d' spot backend is EXPERIMENTAL: a learned heatmap detector trained on the lab's "
                    "Imaris spots of a few N2 gonads (its model card lists them); its counts are not validated and "
                    "not comparable with the spotmax defaults.")
