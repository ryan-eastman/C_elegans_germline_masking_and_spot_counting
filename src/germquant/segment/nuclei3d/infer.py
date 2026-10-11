"""Deterministic tiled inference of a nuclei3d model, the model file, and the instance step (h-maxima seeds
of the distance map, then a seeded watershed inside the foreground).

Normalisation. The DAPI channel is mapped per image to ``(raw - lo) / (hi - lo)`` with ``lo`` and ``hi``
the 1st and 99.9th percentiles of the whole volume read on every 4th voxel in y and x (the net3d
function ``robust_range``; the model file records the settings). The whole volume, never a mask, so the
input does not depend on any segmentation.

Tiling. As for net3d (germquant.spots.net3d.infer, whose tile plan is reused): the volume is zero padded
at its far end to a multiple of the network's size divisor and cut into tiles whose origins are
multiples of that divisor and which overlap by the receptive radius; each tile contributes only its core,
so tiled and untiled predictions agree to float rounding. No tile is skipped: the model sees the whole
field, independent of any other segmentation. float32 throughout, cuDNN deterministic, TF32 off.

The instance step (``instances``):

1. foreground: probability above ``fg_threshold`` (0.5);
2. seed map S = cell distance x (1 - ``neighbour_weight`` x neighbour closeness), clipped to [0, cell
   distance], and 0 outside the foreground: the distance map with the contacts between nuclei pressed
   down, so two touching nuclei keep two domes;
3. seeds: the h-maxima of S with ``h_um`` (0.4 um; skimage.morphology.h_maxima, full connectivity),
   each connected group of maxima one seed;
4. a seeded watershed of -S inside the foreground (skimage.segmentation.watershed, face connectivity).

Steps 3 and 4 run per connected component of the foreground, inside its bounding box, which is exact
(neither the h-maxima of S nor the flooding can cross the background, where S is 0) and keeps the
memory at the size of the largest component. A component whose S never reaches ``h_um`` has no
h-maximum and stays background (it is thinner than 2 x ``h_um``, far below the size filter). Labels
are numbered by component in raster order of the component labelling, then by seed, so the result is
deterministic. The pipeline's own ``_drop_small`` (segmentation.nuclei.min_volume_um3) runs afterwards,
in ``segment_nuclei``.
"""
from __future__ import annotations

import logging
import math
from pathlib import Path

import numpy as np

from germquant.spots.net3d.infer import (  # noqa: F401  (reused: one implementation of each)
    NORM_DEFAULTS,
    deterministic_float32,
    file_sha256,
    normalize_block,
    robust_range,
    tile_plan,
)

log = logging.getLogger(__name__)

MODEL_FORMAT = "germquant-nuclei3d"
MODEL_FORMAT_VERSION = 1
DEFAULT_TILE_ZYX = (96, 384, 384)
POST_DEFAULTS = {"fg_threshold": 0.5, "h_um": 0.4, "neighbour_weight": 1.0}
_MODEL_CACHE: dict = {}


# ------------------------------------------------------------------------------ model file
def _plain(obj):
    """Python scalars, lists and dicts only (numpy values converted), for weights_only loading."""
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Path):
        return str(obj)
    return obj


def save_model_file(path, model, meta: dict) -> str:
    """Write a nuclei3d model file (a torch.save dict of plain values and tensors, loadable with
    ``weights_only=True``) under a temporary name, rename it into place, return its sha256."""
    import torch

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    doc = {"format": MODEL_FORMAT, "format_version": MODEL_FORMAT_VERSION, **_plain(meta),
           "arch": model.arch(), "state_dict": {k: v.detach().cpu() for k, v in model.state_dict().items()}}
    tmp = path.with_name(path.name + ".partial")
    torch.save(doc, tmp)
    tmp.replace(path)
    return file_sha256(path)


def read_model_file(path) -> dict:
    import torch

    doc = torch.load(str(path), map_location="cpu", weights_only=True)
    if not isinstance(doc, dict) or doc.get("format") != MODEL_FORMAT:
        got = doc.get("format") if isinstance(doc, dict) else "?"
        raise ValueError(f"{path} is not a nuclei3d model file (format {got!r})")
    if int(doc.get("format_version", 0)) > MODEL_FORMAT_VERSION:
        raise ValueError(f"{path} is nuclei3d model format {doc.get('format_version')}; this code reads up to "
                         f"{MODEL_FORMAT_VERSION}")
    return doc


def load_model(path, device: str = "cpu"):
    """(model in eval mode on ``device``, the file's metadata without the weights, sha256), cached per
    (path, size, mtime, device) so a batch loads the model once."""
    from .model import build_model

    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"nuclei3d model file not found: {p}")
    st = p.stat()
    key = (str(p.resolve()), st.st_size, st.st_mtime, str(device))
    if key not in _MODEL_CACHE:
        doc = read_model_file(p)
        model = build_model(doc["arch"])
        model.load_state_dict(doc["state_dict"])
        model.to(device).eval()
        meta = {k: v for k, v in doc.items() if k != "state_dict"}
        _MODEL_CACHE.clear()                 # one model at a time is all a run needs
        _MODEL_CACHE[key] = (model, meta, file_sha256(p))
    return _MODEL_CACHE[key]


def post_params(meta: dict | None = None, overrides: dict | None = None) -> dict:
    """The instance-step parameters: the code defaults, then the model file's ``post`` block, then the
    non-None ``overrides`` (the backend's config keys)."""
    out = dict(POST_DEFAULTS)
    out.update({k: v for k, v in ((meta or {}).get("post") or {}).items() if k in POST_DEFAULTS and v is not None})
    out.update({k: v for k, v in (overrides or {}).items() if k in POST_DEFAULTS and v is not None})
    return {k: float(v) for k, v in out.items()}


# --------------------------------------------------------------------------------- inference
def _tiles(model, dna, lo_hi, device: str, tile_zyx):
    """Yield (core slices into the volume, (3, cz, cy, cx) float32 raw outputs of the core)."""
    import torch

    vol = np.asarray(dna) if not hasattr(dna, "shape") else dna
    if len(vol.shape) != 3:
        raise ValueError(f"expected a (Z, Y, X) DAPI volume, got shape {vol.shape}")
    if int(model.in_channels) != 1:
        raise ValueError(f"the nuclei3d model takes {model.in_channels} input channels; only DAPI is supported")
    shape = tuple(int(s) for s in vol.shape)
    lo, hi = lo_hi
    dev = torch.device(device)
    was_training = model.training
    model.eval()
    try:
        with deterministic_float32(), torch.inference_mode():
            for core, win in tile_plan(shape, model.levels, tile_zyx):
                cz, cy, cx = (slice(s.start, min(s.stop, n)) for s, n in zip(core, shape))
                if cz.start >= cz.stop or cy.start >= cy.stop or cx.start >= cx.stop:
                    continue                                 # a core made only of padding
                wshape = tuple(s.stop - s.start for s in win)
                x = np.zeros((1, 1) + wshape, np.float32)
                real = tuple(slice(s.start, min(s.stop, n)) for s, n in zip(win, shape))
                inner = tuple(slice(0, r.stop - r.start) for r in real)
                x[(0, 0) + inner] = normalize_block(vol[real], lo, hi)
                out = model(torch.from_numpy(x).to(dev)).float()[0].cpu().numpy()
                off = tuple(r.start - w.start for r, w in zip((cz, cy, cx), win))
                yield (cz, cy, cx), out[:, off[0]:off[0] + (cz.stop - cz.start), off[1]:off[1] + (cy.stop - cy.start),
                                        off[2]:off[2] + (cx.stop - cx.start)]
    finally:
        if was_training:
            model.train()


def _activate(raw: np.ndarray):
    """(foreground probability, cell distance >= 0 um, neighbour closeness in [0, 1]) from raw outputs."""
    fg = 1.0 / (1.0 + np.exp(-raw[0].astype(np.float32)))
    dist = np.maximum(raw[1], 0.0).astype(np.float32)
    nb = np.clip(raw[2], 0.0, 1.0).astype(np.float32)
    return fg.astype(np.float32), dist, nb


def predict_maps(model, dna, lo_hi, *, device: str = "cpu", tile_zyx=DEFAULT_TILE_ZYX) -> dict:
    """The three maps of a (Z, Y, X) DAPI volume: ``fg`` (foreground probability), ``dist`` (cell distance,
    um, >= 0) and ``nb`` (neighbour closeness, 0 to 1), float32 arrays of the volume's shape."""
    shape = tuple(int(s) for s in dna.shape)
    maps = {k: np.zeros(shape, np.float32) for k in ("fg", "dist", "nb")}
    for core, raw in _tiles(model, dna, lo_hi, device, tile_zyx):
        fg, dist, nb = _activate(raw)
        maps["fg"][core], maps["dist"][core], maps["nb"][core] = fg, dist, nb
    return maps


def seed_map_block(fg_prob, dist, nb, fg_threshold: float, neighbour_weight: float):
    """(foreground mask, seed map S) of a block of maps (module docstring, steps 1 and 2)."""
    fg = np.asarray(fg_prob) > np.float32(fg_threshold)
    damp = np.clip(1.0 - np.float32(neighbour_weight) * np.asarray(nb, np.float32), 0.0, 1.0)
    s = np.where(fg, np.asarray(dist, np.float32) * damp, np.float32(0.0)).astype(np.float32)
    return fg, s


def predict_seed_map(model, dna, lo_hi, *, device: str = "cpu", tile_zyx=DEFAULT_TILE_ZYX,
                     fg_threshold: float = POST_DEFAULTS["fg_threshold"],
                     neighbour_weight: float = POST_DEFAULTS["neighbour_weight"]):
    """(foreground mask, seed map S) of a whole volume without keeping the three float maps: what
    ``instances_from_seed_map`` needs, at a third of the memory of ``predict_maps``."""
    shape = tuple(int(s) for s in dna.shape)
    fg = np.zeros(shape, bool)
    s = np.zeros(shape, np.float32)
    for core, raw in _tiles(model, dna, lo_hi, device, tile_zyx):
        f, d, n = _activate(raw)
        fg[core], s[core] = seed_map_block(f, d, n, fg_threshold, neighbour_weight)
    return fg, s


# ------------------------------------------------------------------------------- instances
def instances_from_seed_map(fg: np.ndarray, s: np.ndarray, *, h_um: float = POST_DEFAULTS["h_um"],
                            keep: np.ndarray | None = None) -> np.ndarray:
    """int32 instance labels from a foreground mask and a seed map (module docstring, steps 3 and 4).

    ``keep`` (a boolean volume, optional) restricts the work to the foreground components that touch it;
    the others stay background. Components are processed independently and completely, so every
    instance of a kept component is exactly what the unrestricted call gives (the validation uses this
    with the scoring region: a component that never touches it can only produce predictions the protocol
    ignores)."""
    from scipy import ndimage as ndi
    from skimage.morphology import h_maxima
    from skimage.segmentation import watershed

    fg = np.asarray(fg, bool)
    s = np.asarray(s, np.float32)
    if fg.shape != s.shape or fg.ndim != 3:
        raise ValueError(f"foreground {fg.shape} and seed map {s.shape} must be (Z, Y, X) volumes of one shape")
    out = np.zeros(fg.shape, np.int32)
    if not fg.any():
        return out
    face = ndi.generate_binary_structure(3, 1)
    full = np.ones((3, 3, 3), bool)
    comp, n_comp = ndi.label(fg, structure=face)
    wanted = None
    if keep is not None:
        keep = np.asarray(keep, bool)
        if keep.shape != fg.shape:
            raise ValueError(f"keep {keep.shape} and foreground {fg.shape} differ in shape")
        wanted = np.zeros(n_comp + 1, bool)
        wanted[np.unique(comp[keep])] = True
        wanted[0] = False
    nxt = 0
    h = float(h_um)
    for k, sl in enumerate(ndi.find_objects(comp), start=1):
        if sl is None or (wanted is not None and not wanted[k]):
            continue
        box = tuple(slice(max(0, a.start - 1), min(n, a.stop + 1)) for a, n in zip(sl, fg.shape))
        m = comp[box] == k
        sc = np.where(m, s[box], np.float32(0.0))
        if not float(sc.max()) >= h:
            continue                      # no dome of height h: no h-maximum (module docstring)
        seeds = (h_maxima(sc, h) > 0) & m
        markers, n_seeds = ndi.label(seeds, structure=full)
        if n_seeds == 0:
            continue
        ws = watershed(-sc, markers, mask=m, connectivity=1)
        sel = ws > 0
        out[box][sel] = ws[sel] + nxt
        nxt += int(n_seeds)
    return out


def instances(maps: dict, *, fg_threshold: float = POST_DEFAULTS["fg_threshold"],
              h_um: float = POST_DEFAULTS["h_um"],
              neighbour_weight: float = POST_DEFAULTS["neighbour_weight"], keep: np.ndarray | None = None) -> np.ndarray:
    """int32 instance labels from the three maps of ``predict_maps`` (``keep`` as in
    ``instances_from_seed_map``)."""
    fg = np.zeros(maps["fg"].shape, bool)
    s = np.zeros(maps["fg"].shape, np.float32)
    for z in range(fg.shape[0]):
        fg[z], s[z] = seed_map_block(maps["fg"][z], maps["dist"][z], maps["nb"][z], fg_threshold, neighbour_weight)
    return instances_from_seed_map(fg, s, h_um=h_um, keep=keep)


def segment_volume(dna, spacing, model_path, *, device: str = "cpu", tile_zyx=DEFAULT_TILE_ZYX,
                   post: dict | None = None) -> tuple[np.ndarray, dict]:
    """(int32 labels before the size filter, info) for a raw (Z, Y, X) DAPI volume: normalisation,
    tiled inference, the instance step. ``post`` overrides the model file's instance parameters."""
    import time

    net, meta, sha = load_model(model_path, device)
    params = post_params(meta, post)
    norm = {**NORM_DEFAULTS, **(meta.get("norm") or {})}
    t0 = time.time()
    lo_hi = robust_range(dna, norm["low_pct"], norm["high_pct"], norm["stride_yx"])
    trained = meta.get("spacing_zyx_um")
    if trained and any(abs(float(a) - float(b)) > 0.1 * float(b) for a, b in zip(spacing, trained)):
        log.warning("nuclei3d: the image's voxel size %s differs by more than 10 percent from the training data's "
                    "%s; the model sees nuclei at another scale", tuple(spacing), tuple(trained))
    fg, s = predict_seed_map(net, dna, lo_hi, device=device, tile_zyx=tile_zyx,
                             fg_threshold=params["fg_threshold"], neighbour_weight=params["neighbour_weight"])
    t1 = time.time()
    labels = instances_from_seed_map(fg, s, h_um=params["h_um"])
    t2 = time.time()
    info = {"model": str(model_path), "model_sha256": sha, "device": device, "post": params, "norm_range": list(lo_hi),
            "foreground_fraction": float(fg.mean()), "n_instances": int(labels.max()),
            "infer_s": round(t1 - t0, 1), "instances_s": round(t2 - t1, 1)}
    return labels, info


def parse_zyx(text) -> tuple[int, int, int]:
    """'Z,Y,X' (or a sequence) into three positive ints."""
    if isinstance(text, (list, tuple)):
        vals = [int(v) for v in text]
    else:
        vals = [int(v) for v in str(text).replace("x", ",").split(",") if v.strip()]
    if len(vals) != 3 or any(v <= 0 for v in vals):
        raise ValueError(f"expected three positive integers Z,Y,X, got {text!r}")
    return tuple(vals)


def finite_or_none(x):
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None
