"""Deterministic tiled inference of a net3d model in float32, and the model file.

Normalisation. Every input channel is mapped per gonad (per image) to ``(raw - lo) / (hi - lo)`` with
``lo`` and ``hi`` the ``low_pct`` and ``high_pct`` percentiles (code defaults 1 and 99.9) of the whole
volume, read on every ``stride_yx``-th voxel in y and x of every plane (4 by default; about 30
million voxels on a 71 x 2600 x 2600 stack). The whole volume, not the nucleus masks, so the input of
the network never depends on the segmentation. The same function serves training and inference, and
the model file records the settings.

Tiling. The volume is zero padded (in normalised units) at its far end to a multiple of the network's
size divisor, then cut into tiles whose origins are multiples of that divisor (so every pooling window
of a tile is a pooling window of the whole volume) and which overlap by the network's receptive
radius (``model.tile_margin``). Each tile contributes only its core, the part at least that margin
away from any cut, so a voxel's prediction never sees a tile edge: tiled and untiled predictions agree
to float rounding (the test suite checks this). With ``labels`` given, tiles whose core holds no
nucleus voxel (after a small dilation) are skipped and stay 0; the backend drops every spot without a
nucleus within ``assign_max_um`` anyway, so the skip changes no result.

Determinism: float32 throughout, no autocast, cuDNN deterministic algorithms with benchmarking and
TF32 off for the duration of the call (restored afterwards). A CPU and a CUDA run agree to float
rounding, not bit for bit.
"""
from __future__ import annotations

import contextlib
import hashlib
import logging
import math
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)

MODEL_FORMAT = "germquant-net3d"
MODEL_FORMAT_VERSION = 1
NORM_DEFAULTS = {"low_pct": 1.0, "high_pct": 99.9, "stride_yx": 4}
DEFAULT_TILE_ZYX = (96, 384, 384)
_MODEL_CACHE: dict = {}


# ---------------------------------------------------------------------------- normalisation
def robust_range(arr, low_pct: float = NORM_DEFAULTS["low_pct"], high_pct: float = NORM_DEFAULTS["high_pct"],
                 stride_yx: int = NORM_DEFAULTS["stride_yx"]) -> tuple[float, float]:
    """(lo, hi): the two percentiles of a (Z, Y, X) volume on a fixed stride (module docstring)."""
    a = np.asarray(arr)
    if a.ndim != 3:
        raise ValueError(f"expected a (Z, Y, X) volume, got shape {a.shape}")
    s = max(1, int(stride_yx))
    sub = np.asarray(a[:, ::s, ::s], dtype=np.float32)
    lo, hi = np.percentile(sub, [float(low_pct), float(high_pct)])
    lo, hi = float(lo), float(hi)
    if not math.isfinite(lo) or not math.isfinite(hi):
        raise ValueError("the volume holds non-finite values")
    if hi <= lo:
        hi = lo + 1.0                       # a flat image: keep the scale finite
    return lo, hi


def normalize_block(block, lo: float, hi: float) -> np.ndarray:
    return (np.asarray(block, dtype=np.float32) - np.float32(lo)) / np.float32(hi - lo)


def channel_ranges(channels, norm: dict | None = None) -> list[tuple[float, float]]:
    n = {**NORM_DEFAULTS, **(norm or {})}
    return [robust_range(c, n["low_pct"], n["high_pct"], n["stride_yx"]) for c in channels]


# ------------------------------------------------------------------------------ model file
def file_sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def save_model_file(path, model, meta: dict) -> str:
    """Write a net3d model file (a torch.save dict of plain values and tensors, loadable with
    ``weights_only=True``) under a temporary name, rename it into place and return its sha256."""
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
        raise ValueError(f"{path} is not a net3d model file (format {doc.get('format') if isinstance(doc, dict) else '?'!r})")
    if int(doc.get("format_version", 0)) > MODEL_FORMAT_VERSION:
        raise ValueError(f"{path} is net3d model format {doc.get('format_version')}; this code reads up to "
                         f"{MODEL_FORMAT_VERSION}")
    return doc


def load_model(path, device: str = "cpu"):
    """(model in eval mode on ``device``, the file's metadata without the weights, sha256). Cached
    per (path, size, mtime, device), so a batch loads a model once."""
    from .model import build_model

    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"net3d model file not found: {p}")
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


# --------------------------------------------------------------------------------- tiling
def _axis_tiles(n_padded: int, tile: int, margin: int, div: int) -> list[tuple[int, int, int, int]]:
    """(core_start, core_stop, window_start, window_stop) along one padded axis."""
    tile = max(div, (int(tile) // div) * div)
    if tile >= n_padded:
        return [(0, n_padded, 0, n_padded)]
    core = tile - 2 * margin
    if core < div:
        raise ValueError(f"tile extent {tile} leaves no core inside the {margin}-voxel overlap; use a tile of at "
                         f"least {2 * margin + div}")
    out = []
    for c0 in range(0, n_padded, core):
        c1 = min(c0 + core, n_padded)
        out.append((c0, c1, max(0, c0 - margin), min(n_padded, c1 + margin)))
    return out


def tile_plan(shape_zyx, levels: int, tile_zyx=DEFAULT_TILE_ZYX) -> list[tuple[tuple[slice, ...], tuple[slice, ...]]]:
    """Every tile as (core slices, window slices) in padded coordinates."""
    from .model import size_divisor, tile_margin

    div = size_divisor(levels)
    margin = tile_margin(levels)
    padded = [int(math.ceil(s / d) * d) for s, d in zip(shape_zyx, div)]
    if tile_zyx is None:
        tile_zyx = padded
    tz, ty, tx = (padded[a] if tile_zyx[a] is None else int(tile_zyx[a]) for a in range(3))
    axes = [_axis_tiles(padded[0], tz, margin[0], div[0]), _axis_tiles(padded[1], ty, margin[1], div[1]),
            _axis_tiles(padded[2], tx, margin[2], div[2])]
    plan = []
    for z in axes[0]:
        for y in axes[1]:
            for x in axes[2]:
                plan.append(((slice(z[0], z[1]), slice(y[0], y[1]), slice(x[0], x[1])),
                             (slice(z[2], z[3]), slice(y[2], y[3]), slice(x[2], x[3]))))
    return plan


@contextlib.contextmanager
def deterministic_float32():
    """cuDNN deterministic, no benchmarking, no TF32, for the duration of the block."""
    import torch

    saved = (torch.backends.cudnn.deterministic, torch.backends.cudnn.benchmark,
             torch.backends.cudnn.allow_tf32, torch.backends.cuda.matmul.allow_tf32)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.allow_tf32 = False
    torch.backends.cuda.matmul.allow_tf32 = False
    try:
        yield
    finally:
        (torch.backends.cudnn.deterministic, torch.backends.cudnn.benchmark,
         torch.backends.cudnn.allow_tf32, torch.backends.cuda.matmul.allow_tf32) = saved


def occupancy(labels, dilate_px: int = 4) -> tuple[np.ndarray, np.ndarray]:
    """(z planes with any label, (Y, X) projection of the labels dilated by ``dilate_px``): the tile skip."""
    from scipy import ndimage as ndi

    lab = np.asarray(labels)
    zany = np.zeros(lab.shape[0], bool)
    proj = np.zeros(lab.shape[1:], bool)
    for z in range(lab.shape[0]):
        m = lab[z] > 0
        if m.any():
            zany[z] = True
            proj |= m
    if dilate_px > 0 and proj.any():
        proj = ndi.binary_dilation(proj, structure=np.ones((2 * dilate_px + 1,) * 2, bool))
        zany = ndi.binary_dilation(zany, structure=np.ones(3, bool), iterations=max(1, dilate_px // 2))
    return zany, proj


def predict_heatmap(model, channels, ranges, *, device: str = "cpu", tile_zyx=DEFAULT_TILE_ZYX,
                    labels=None, occ=None, dilate_px: int = 4) -> np.ndarray:
    """Per-voxel focus probability (float32, the shape of the channels) from raw (Z, Y, X) channels
    and their (lo, hi) normalisation ranges (module docstring). ``model`` is a Net3D in eval mode;
    ``labels`` (or a precomputed ``occ`` from ``occupancy``) switches on the tile skip."""
    import torch

    chans = [np.asarray(c) for c in channels]
    if len(chans) != int(model.in_channels):
        raise ValueError(f"the model takes {model.in_channels} input channel(s), got {len(chans)}")
    shape = tuple(chans[0].shape)
    if any(tuple(c.shape) != shape for c in chans) or len(shape) != 3:
        raise ValueError(f"input channels must be (Z, Y, X) volumes of one shape, got {[c.shape for c in chans]}")
    if occ is None and labels is not None:
        if tuple(np.asarray(labels).shape) != shape:
            raise ValueError(f"labels {np.asarray(labels).shape} and image {shape} differ in shape")
        occ = occupancy(labels, dilate_px)
    heat = np.zeros(shape, np.float32)
    dev = torch.device(device)
    was_training = model.training
    model.eval()
    with deterministic_float32(), torch.inference_mode():
        for core, win in tile_plan(shape, model.levels, tile_zyx):
            cz, cy, cx = (slice(s.start, min(s.stop, n)) for s, n in zip(core, shape))
            if cz.start >= cz.stop or cy.start >= cy.stop or cx.start >= cx.stop:
                continue                                     # a core made only of padding
            if occ is not None and not (occ[0][cz].any() and occ[1][cy, cx].any()):
                continue
            wshape = tuple(s.stop - s.start for s in win)
            x = np.zeros((1, len(chans)) + wshape, np.float32)
            real = tuple(slice(s.start, min(s.stop, n)) for s, n in zip(win, shape))
            inner = tuple(slice(0, r.stop - r.start) for r in real)
            for k, (c, (lo, hi)) in enumerate(zip(chans, ranges)):
                x[(0, k) + inner] = normalize_block(c[real], lo, hi)
            prob = torch.sigmoid(model(torch.from_numpy(x).to(dev)).float())[0, 0].cpu().numpy()
            off = tuple(r.start - w.start for r, w in zip((cz, cy, cx), win))
            heat[cz, cy, cx] = prob[off[0]:off[0] + (cz.stop - cz.start), off[1]:off[1] + (cy.stop - cy.start),
                                    off[2]:off[2] + (cx.stop - cx.start)]
    if was_training:
        model.train()
    return heat
