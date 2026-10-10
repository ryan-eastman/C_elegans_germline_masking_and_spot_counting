"""Training data for net3d: patches from the bench cache, heatmap targets at the expert's Imaris
spots and the loss weights (numpy and scipy only; train.py turns the arrays into tensors).

What a patch holds. The input channels (RAD-51 and, unless the model is the foci-only ablation, DAPI)
normalised per gonad with ``infer.robust_range``; the target, a Gaussian of peak 1 at every Class A
truth spot (``sigma_um`` per axis, (z, xy) = (0.25, 0.12) um by default: about 1.25 voxels in z and
1.1 in xy, the size of the confocal image of a focus; where two Gaussians overlap the target is their
maximum); and the per-voxel loss weight.

Where the truth is known, and the weights. The expert's spots were detected by Imaris over the whole
field with one Quality threshold per gonad, set by eye (docs/BENCHMARK_PROTOCOL.md, section 4.1), and
the counts that matter are per nucleus. The loss is restricted to the germline as the pipeline sees
it: the germline nucleus labels of the
golden run (the cached ``labels.npy`` with the ``germline_label_ids`` of ``meta.json``) dilated by
``known_um`` (0.3 um, the nearest-label reach of the assignment). Inside that region:

* weight 1 inside the traced truth coverage, the protocol's region C (the expert's nucleus masks
  dilated by 1.0 um in xy and one plane in z): there the annotator looked at the nuclei the spots
  belong to;
* weight ``w_untraced`` (0.5 by default) elsewhere in the germline: the same automatic detection with the
  same threshold, but nucleus by nucleus nobody looked, so a missed or an extra spot there is weaker
  evidence;
* weight 0 inside an ellipsoid of ``ignore_um`` (0.5 um in z, 0.3 um in xy) around every Class B spot.
  Class B is Imaris's "isolated" class (average distance to the 3 nearest spots above a per-gonad
  cut), not a false positive, so the network is neither taught to find nor to reject those spots;
* every weight is multiplied by ``1 + pos_weight * target`` (``pos_weight`` 9: ten times the weight at a
  spot centre), which spends the gradient where the foci are without moving the per-voxel optimum
  (each voxel's weighted cross-entropy is still minimised by predicting its own target);
* outside the volume (a patch may hang over the stack's edge by a quarter of its extent) the input is
  0 and the weight is 0.

Sampling. Each sample picks a gonad uniformly, then with probability ``p_spot`` (0.5) a Class A spot
inside the weighted region plus a uniform offset of up to a third of the patch, else the centroid of a
random germline nucleus plus up to half a patch. Augmentation: flips in y and x and a y/x transpose (the
xy grid is isotropic), and a log-uniform gain per channel (RAD-51 0.75 to 1.33, DAPI 0.8 to 1.25). Every
sample of an epoch draws from its own generator seeded with (seed, epoch, index), so a run is
reproducible whatever the batching.

The region map (uint8: 0 outside, 1 germline untraced, 2 germline traced) is computed once per gonad
by exact ellipsoidal and cylindrical dilations (per-plane Euclidean distance transforms, restricted to
the bounding box of the mask) and kept as an .npy file in the prep folder, with a JSON of what it was
made from; a mismatch recomputes it.
"""
from __future__ import annotations

import json
import logging
import math
import threading
from dataclasses import dataclass, field
from pathlib import Path
from queue import Queue

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from .infer import channel_ranges, normalize_block

log = logging.getLogger(__name__)

REGION_VERSION = 1
DATA_DEFAULTS = {
    "sigma_um": (0.25, 0.12),        # target Gaussian sigma (z, xy)
    "known_um": 0.3,                 # germline labels dilated by this (the loss region)
    "coverage_xy_um": 1.0,           # truth masks dilated by this in xy ...
    "coverage_z_planes": 1,          # ... and this many planes in z: the protocol's region C
    "w_untraced": 0.5,
    "ignore_um": (0.5, 0.3),         # zero-weight ellipsoid around a Class B spot (z, xy)
    "pos_weight": 9.0,
    "p_spot": 0.5,
    "max_overhang": 0.25,            # a patch may hang over the stack's edge by this fraction
    "gain_foci": (0.75, 1.3333333),
    "gain_dna": (0.8, 1.25),
}


# ------------------------------------------------------------------------------ dilations
def _bbox_yx(mask_any_yx: np.ndarray, pad: int, shape_yx) -> tuple[slice, slice] | None:
    ys, xs = np.nonzero(mask_any_yx)
    if ys.size == 0:
        return None
    return (slice(max(0, ys.min() - pad), min(shape_yx[0], ys.max() + pad + 1)),
            slice(max(0, xs.min() - pad), min(shape_yx[1], xs.max() + pad + 1)))


def _plane_distance(plane: np.ndarray, sampling) -> np.ndarray:
    """Distance (um) of every pixel to the nearest True pixel of the plane (inf when there is none)."""
    if not plane.any():
        return np.full(plane.shape, np.inf, np.float32)
    return ndi.distance_transform_edt(~plane, sampling=sampling).astype(np.float32)


def dilate_planes(mask: np.ndarray, spacing, radius_xy_of_dz) -> np.ndarray:
    """Dilate a (Z, Y, X) boolean mask: voxel (z, y, x) is set when some mask voxel in plane z + k lies
    within ``radius_xy_of_dz[k]`` um in xy, for every plane offset k the mapping holds. Exact for any
    shape that is a stack of discs, which is what the ellipsoid and the cylinder are."""
    m = np.asarray(mask, bool)
    Z, Y, X = m.shape
    out = np.zeros_like(m)
    if not m.any():
        return out
    rmax = max(radius_xy_of_dz.values())
    pad = int(math.ceil(rmax / min(spacing[1], spacing[2]))) + 1
    any_yx = m.any(axis=0)
    box = _bbox_yx(any_yx, pad, (Y, X))
    sub = m[:, box[0], box[1]]
    nz = max(abs(k) for k in radius_xy_of_dz)
    dists: dict[int, np.ndarray] = {}        # a sliding window of per-plane distance maps
    for z in range(Z):
        acc = np.zeros(sub.shape[1:], bool)
        for k, r in radius_xy_of_dz.items():
            zz = z + k
            if zz < 0 or zz >= Z:
                continue
            if zz not in dists:
                dists[zz] = _plane_distance(sub[zz], (spacing[1], spacing[2]))
            acc |= dists[zz] <= r + 1e-6
        out[z, box[0], box[1]] = acc
        for old in [p for p in dists if p < z + 1 - nz]:
            del dists[old]
    return out


def dilate_ellipsoid(mask: np.ndarray, radius_um: float, spacing) -> np.ndarray:
    """Exact dilation by the ellipsoid (a sphere in microns) of ``radius_um``."""
    dz = float(spacing[0])
    nz = int(math.floor(radius_um / dz + 1e-9))
    radii = {k: math.sqrt(max(radius_um ** 2 - (k * dz) ** 2, 0.0)) for k in range(-nz, nz + 1)}
    return dilate_planes(mask, spacing, radii)


def dilate_cylinder(mask: np.ndarray, radius_xy_um: float, n_planes: int, spacing) -> np.ndarray:
    """Dilation by ``radius_xy_um`` in xy and ``n_planes`` planes in z (the protocol's region C)."""
    radii = {k: float(radius_xy_um) for k in range(-int(n_planes), int(n_planes) + 1)}
    return dilate_planes(mask, spacing, radii)


def germline_mask(labels: np.ndarray, germline_ids) -> np.ndarray:
    """Boolean mask of the voxels whose label is a germline nucleus (lookup table, plane by plane)."""
    lab = np.asarray(labels)
    top = int(max(int(lab.max()) if lab.size else 0, max([int(i) for i in germline_ids], default=0)))
    lut = np.zeros(top + 1, bool)
    ids = np.asarray([int(i) for i in germline_ids if int(i) > 0], np.int64)
    lut[ids] = True
    out = np.zeros(lab.shape, bool)
    for z in range(lab.shape[0]):
        out[z] = lut[np.asarray(lab[z])]
    return out


def build_region(labels, germline_ids, truth_labels, spacing, cfg: dict | None = None) -> np.ndarray:
    """The uint8 region map: 0 outside the loss region, 1 germline (dilated) outside the traced
    coverage, 2 germline inside it (module docstring)."""
    c = {**DATA_DEFAULTS, **(cfg or {})}
    known = dilate_ellipsoid(germline_mask(labels, germline_ids), float(c["known_um"]), spacing)
    region = known.astype(np.uint8)
    if truth_labels is not None:
        traced = np.zeros(known.shape, bool)
        tl = np.asarray(truth_labels)
        for z in range(tl.shape[0]):
            traced[z] = tl[z] > 0
        traced = dilate_cylinder(traced, float(c["coverage_xy_um"]), int(c["coverage_z_planes"]), spacing)
        region[known & traced] = 2
    return region


# ----------------------------------------------------------------------------- truth spots
def truth_positions(spots: pd.DataFrame, params: dict | None) -> np.ndarray:
    """(n, 3) continuous (z, y, x) positions in voxel units of the image grid. Imaris reports voxel
    centres, so ``(imaris_coord - ext_min) / imaris_voxel - 0.5`` puts a spot on the grid where voxel i
    is centred on i; its floor-based integer index (the truth files' zi, yi, xi) is that position
    rounded. Without the Imaris columns or the params, the integer indices are used as they are."""
    cols = ["imaris_z_um", "imaris_y_um", "imaris_x_um"]
    p = params or {}
    ext = p.get("ext_min_xyz_um")
    vox = p.get("voxel_imaris_xyz_um") or p.get("imaris_voxel_xyz_um")
    if len(spots) and ext is not None and vox is not None and all(c in spots.columns for c in cols) \
            and spots[cols].notna().all().all():
        ext_zyx = np.asarray(ext, float)[::-1]
        vox_zyx = np.asarray(vox, float)[::-1]
        return (spots[cols].to_numpy(float) - ext_zyx) / vox_zyx - 0.5
    return spots[["zi", "yi", "xi"]].to_numpy(float).reshape(-1, 3)


def class_b_mask(spot_class) -> np.ndarray:
    from germquant.benchmark.spots import is_class_b

    return is_class_b(spot_class)


# --------------------------------------------------------------------------------- a gonad
@dataclass
class TrainGonad:
    tag: str
    spacing: tuple[float, float, float]
    channels: list                    # raw (Z, Y, X) arrays (memory maps), foci first
    ranges: list                      # (lo, hi) per channel
    region: np.ndarray                # uint8 (Z, Y, X)
    spots_a: np.ndarray               # (n, 3) continuous voxel positions of the Class A spots
    spots_b: np.ndarray               # (m, 3) of the Class B spots
    centers: np.ndarray               # (k, 3) germline nucleus centroids in voxels
    condition: str = ""
    sex: str = ""
    info: dict = field(default_factory=dict)

    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(int(s) for s in self.region.shape)

    def __post_init__(self):
        pts = np.rint(self.spots_a).astype(np.int64) if len(self.spots_a) else np.zeros((0, 3), np.int64)
        shape = np.asarray(self.region.shape)
        ok = np.all((pts >= 0) & (pts < shape), axis=1) if len(pts) else np.zeros(0, bool)
        inside = np.zeros(len(pts), bool)
        if ok.any():
            inside[ok] = self.region[pts[ok, 0], pts[ok, 1], pts[ok, 2]] > 0
        self.spot_pool = np.flatnonzero(inside)      # Class A spots that can centre a patch


def region_key(meta: dict, truth_labels_path, cfg: dict) -> dict:
    """What a cached region map was made from (a mismatch recomputes it)."""
    tl = Path(truth_labels_path) if truth_labels_path else None
    return {"version": REGION_VERSION, "tag": meta.get("tag"),
            "labels_sha256": ((meta.get("arrays") or {}).get("labels") or {}).get("data_sha256"),
            "germline_ids_n": len(meta.get("germline_label_ids", [])),
            "germline_ids_sum": int(sum(int(i) for i in meta.get("germline_label_ids", []))),
            "truth_labels": str(tl) if tl else None,
            "truth_labels_size": tl.stat().st_size if tl and tl.is_file() else None,
            "known_um": float(cfg["known_um"]), "coverage_xy_um": float(cfg["coverage_xy_um"]),
            "coverage_z_planes": int(cfg["coverage_z_planes"])}


def load_region(prep_dir, tag: str, meta: dict, labels, truth_labels_path, spacing, cfg: dict) -> np.ndarray:
    """The region map of a gonad from the prep folder, computed and saved when absent or stale."""
    prep = Path(prep_dir)
    prep.mkdir(parents=True, exist_ok=True)
    npy, js = prep / f"{tag}__region.npy", prep / f"{tag}__region.json"
    key = region_key(meta, truth_labels_path, cfg)
    if npy.is_file() and js.is_file():
        try:
            if json.loads(js.read_text(encoding="utf-8")) == key:
                return np.load(npy, mmap_mode="r")
        except (OSError, ValueError):
            pass
    truth_labels = None
    if truth_labels_path:
        from germquant.benchmark.manifest import load_labels

        truth_labels = load_labels(truth_labels_path)
    region = build_region(labels, meta.get("germline_label_ids", []), truth_labels, spacing, cfg)
    del truth_labels
    tmp = npy.with_name(npy.stem + ".partial.npy")
    np.save(tmp, region)
    tmp.replace(npy)
    js.write_text(json.dumps(key, indent=1), encoding="utf-8")
    return np.load(npy, mmap_mode="r")


# ----------------------------------------------------------------------------- one patch
def _window(origin, size, shape):
    """Slices into the volume and into the patch for a patch at ``origin`` (may hang over the edge)."""
    vol, pat = [], []
    for o, s, n in zip(origin, size, shape):
        a, b = max(0, o), min(n, o + s)
        if b <= a:
            return None, None
        vol.append(slice(a, b))
        pat.append(slice(a - o, b - o))
    return tuple(vol), tuple(pat)


def render_target(size, origin, spots: np.ndarray, sigma_vox) -> np.ndarray:
    """Max of unit-peak Gaussians at the spots (global voxel positions) on a patch grid."""
    t = np.zeros(size, np.float32)
    if len(spots) == 0:
        return t
    s = np.asarray(sigma_vox, float)
    reach = np.ceil(4.0 * s).astype(int)
    local = spots - np.asarray(origin, float)
    near = np.all((local >= -reach) & (local < np.asarray(size) + reach), axis=1)
    for p in local[near]:
        lo = np.maximum(np.floor(p - reach).astype(int), 0)
        hi = np.minimum(np.ceil(p + reach).astype(int) + 1, size)
        if np.any(hi <= lo):
            continue
        g = [np.exp(-0.5 * ((np.arange(lo[a], hi[a]) - p[a]) / s[a]) ** 2) for a in range(3)]
        blob = (g[0][:, None, None] * g[1][None, :, None] * g[2][None, None, :]).astype(np.float32)
        sl = tuple(slice(lo[a], hi[a]) for a in range(3))
        np.maximum(t[sl], blob, out=t[sl])
    return t


def ignore_mask(size, origin, spots: np.ndarray, radius_vox) -> np.ndarray:
    """True inside the ellipsoid of ``radius_vox`` (z, y, x) around every spot."""
    m = np.zeros(size, bool)
    if len(spots) == 0:
        return m
    r = np.asarray(radius_vox, float)
    reach = np.ceil(r).astype(int)
    local = spots - np.asarray(origin, float)
    near = np.all((local >= -reach - 1) & (local < np.asarray(size) + reach + 1), axis=1)
    for p in local[near]:
        lo = np.maximum(np.floor(p - r).astype(int), 0)
        hi = np.minimum(np.ceil(p + r).astype(int) + 1, size)
        if np.any(hi <= lo):
            continue
        zz, yy, xx = np.meshgrid(*(np.arange(lo[a], hi[a]) for a in range(3)), indexing="ij")
        inside = ((zz - p[0]) / r[0]) ** 2 + ((yy - p[1]) / r[1]) ** 2 + ((xx - p[2]) / r[2]) ** 2 <= 1.0
        sl = tuple(slice(lo[a], hi[a]) for a in range(3))
        m[sl] |= inside
    return m


def make_patch(g: TrainGonad, origin, size, cfg: dict, rng: np.random.Generator | None = None,
               augment: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(input (C, z, y, x), target (1, z, y, x), weight (1, z, y, x)) float32 for a patch at ``origin``."""
    c = {**DATA_DEFAULTS, **cfg}
    size = tuple(int(s) for s in size)
    sp = np.asarray(g.spacing, float)
    x = np.zeros((len(g.channels),) + size, np.float32)
    region = np.zeros(size, np.uint8)
    vol, pat = _window(origin, size, g.shape)
    if vol is not None:
        for k, (ch, (lo, hi)) in enumerate(zip(g.channels, g.ranges)):
            x[(k,) + pat] = normalize_block(ch[vol], lo, hi)
        region[pat] = g.region[vol]
    sz, sxy = c["sigma_um"]
    target = render_target(size, origin, g.spots_a, (sz / sp[0], sxy / sp[1], sxy / sp[2]))
    iz, ixy = c["ignore_um"]
    ign = ignore_mask(size, origin, g.spots_b, (iz / sp[0], ixy / sp[1], ixy / sp[2]))
    lut = np.array([0.0, float(c["w_untraced"]), 1.0], np.float32)
    w = lut[region] * (1.0 + np.float32(c["pos_weight"]) * target)
    w[ign] = 0.0
    if augment and rng is not None:
        gains = [math.exp(rng.uniform(math.log(c["gain_foci"][0]), math.log(c["gain_foci"][1])))]
        gains += [math.exp(rng.uniform(math.log(c["gain_dna"][0]), math.log(c["gain_dna"][1])))
                  for _ in range(len(g.channels) - 1)]
        for k, gk in enumerate(gains):
            x[k] *= np.float32(gk)
        flips = rng.integers(0, 2, size=3)
        x, target, w = _flip(x, flips), _flip(target[None], flips)[0], _flip(w[None], flips)[0]
    return x, target[None], w[None].astype(np.float32)


def _flip(a: np.ndarray, flips) -> np.ndarray:
    """Flip y, flip x, transpose y/x of a (C, Z, Y, X) array as ``flips`` (three 0/1 flags) says."""
    if flips[0]:
        a = a[:, :, ::-1, :]
    if flips[1]:
        a = a[:, :, :, ::-1]
    if flips[2]:
        a = a.transpose(0, 1, 3, 2)
    return np.ascontiguousarray(a)


def sample_origin(g: TrainGonad, size, cfg: dict, rng: np.random.Generator) -> tuple[int, int, int]:
    """A patch origin: around a Class A spot of the weighted region, else around a germline nucleus."""
    c = {**DATA_DEFAULTS, **cfg}
    size = np.asarray(size, int)
    if len(g.spot_pool) and (rng.random() < float(c["p_spot"]) or not len(g.centers)):
        centre = g.spots_a[g.spot_pool[rng.integers(len(g.spot_pool))]]
        centre = centre + rng.uniform(-size / 3.0, size / 3.0)
    elif len(g.centers):
        centre = g.centers[rng.integers(len(g.centers))] + rng.uniform(-size / 2.0, size / 2.0)
    else:
        centre = rng.uniform(0, np.asarray(g.shape, float))
    origin = np.rint(centre - size / 2.0).astype(int)
    hang = np.floor(size * float(c["max_overhang"])).astype(int)
    lo = -hang
    hi = np.asarray(g.shape) - size + hang
    hi = np.maximum(hi, lo)
    return tuple(int(v) for v in np.clip(origin, lo, hi))


def sample(gonads: list[TrainGonad], size, cfg: dict, seed: int, epoch: int, index: int):
    """One training sample (input, target, weight), reproducible from (seed, epoch, index)."""
    rng = np.random.default_rng([int(seed), int(epoch), int(index)])
    g = gonads[int(rng.integers(len(gonads)))]
    origin = sample_origin(g, size, cfg, rng)
    return make_patch(g, origin, size, cfg, rng, augment=True)


def batches(gonads, size, cfg: dict, seed: int, epoch: int, n_steps: int, batch_size: int):
    """The batches of one epoch, in order: (inputs (B, C, ...), targets (B, 1, ...), weights)."""
    for step in range(int(n_steps)):
        items = [sample(gonads, size, cfg, seed, epoch, step * batch_size + b) for b in range(batch_size)]
        yield (np.stack([i[0] for i in items]), np.stack([i[1] for i in items]), np.stack([i[2] for i in items]))


class Prefetch:
    """Run a generator in a background thread, ``depth`` items ahead (numpy releases the GIL for the
    copies and arithmetic, so patch assembly overlaps the GPU step). Order is preserved."""

    _END = object()

    def __init__(self, gen, depth: int = 3):
        self.q: Queue = Queue(maxsize=max(1, int(depth)))
        self.err: BaseException | None = None
        self.t = threading.Thread(target=self._run, args=(gen,), daemon=True)
        self.t.start()

    def _run(self, gen):
        try:
            for item in gen:
                self.q.put(item)
        except BaseException as e:  # noqa: BLE001 - re-raised in the consumer
            self.err = e
        finally:
            self.q.put(self._END)

    def __iter__(self):
        while True:
            item = self.q.get()
            if item is self._END:
                if self.err is not None:
                    raise self.err
                return
            yield item


def gonad_ranges(channels, norm: dict | None = None) -> list[tuple[float, float]]:
    return channel_ranges(channels, norm)
