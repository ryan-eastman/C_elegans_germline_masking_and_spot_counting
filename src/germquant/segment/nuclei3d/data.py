"""Training data for nuclei3d: targets and per-voxel weights from the recovered Imaris nucleus masks and the
deployed Cellpose labels, and patches of the DAPI channel from the bench cache (numpy and scipy only;
train.py turns the arrays into tensors).

The truth. ``<stem>__imaris_labels.tif`` of each train gonad (germquant imaris-truth): the masks of the
nuclei the annotator traced, 64 to 146 per gonad, every surface-matched label kept whatever its size
(a size cut would delete real surfaces; docs/DESIGN_inhouse_components.md). Only a subset of the
nuclei of a gonad is traced, so the truth says nothing about most nuclei. The deployed Cellpose-SAM
labels of the golden run (the bench cache's ``labels.npy``) say where the other nuclei probably are,
and that is all they are used for: they decide which voxels are left out of the loss, never a target
value inside a traced nucleus.

Which Cellpose labels belong to traced nuclei. A Cellpose label with at least ``traced_frac`` (0.5) of
its voxels inside traced masks is the traced nucleus's own label (a matched pair puts about 0.7 there,
since the Cellpose masks are about 1.4 times the Imaris volume); every other Cellpose label is an
untraced nucleus. On the train and val gonads the fraction is bimodal (below 0.05 for a neighbour that
only touches a traced mask, 0.4 to 1.0 for the nucleus's own label), so the cut is not delicate; a label
that merges a traced nucleus with an untraced one falls below it and is treated as untraced, the safe side.

Classes of the voxels (the uint8 class map) and their loss weights:

* traced (2): inside a traced mask. Weight 1, target foreground with its distances.
* shell (3): inside a traced nucleus's own Cellpose label but outside its traced mask. The annotator
  traced that nucleus and its surface ends at the mask, so by the truth this is background: weight
  ``w_traced_shell`` (1 by default; 0 treats it as unknown, like an untraced label), target background.
* ignore (0): inside an untraced Cellpose label or within ``margin_um`` (0.5 um, an exact ellipsoidal
  dilation) of one, and not inside a traced mask: weight 0 (partial annotation). This takes precedence
  over the shell.
* background (1): everything else, where neither a traced mask nor a deployed Cellpose label is present.
  Weight 1, target background. A nucleus that Cellpose missed is taught as background here, a known
  limit of using the deployed labels as the gate.

Targets inside a traced mask (0 in every other voxel):

* foreground 1;
* cell distance (um): the Euclidean distance (with the voxel size) from the voxel to the nearest voxel
  outside its own nucleus, computed per traced label in its bounding box with the edge of the stack
  counting as outside; it drops to about one voxel at a contact with another traced nucleus;
* neighbour closeness: (1 - d / ``neighbour_radius_um``) ** ``neighbour_power``, clipped to [0, 1], where
  d is the distance (um) to the nearest voxel of another nucleus: another traced mask, or a Cellpose
  label that is not the nucleus's own (one that holds less than ``own_frac`` of the traced mask), the
  latter standing in for untraced neighbours. 1 at a contact, 0 beyond the radius (1.0 um, power 3: the
  inverted and steepened neighbour distance of Scherr et al. 2020).

The maps are computed once per gonad into the prep folder (``<tag>__n3d_cls.npy`` uint8,
``__n3d_dist.npy`` and ``__n3d_nb.npy`` float16, ``__n3d_targets.json`` with what they were made from;
a mismatch recomputes them) and memory-mapped for training.

Patches. Each sample picks a train gonad uniformly, then a centre: with probability ``p_traced`` (0.6) a
traced nucleus (its centroid plus up to a third of the patch), with ``p_germline`` (0.3) a germline
nucleus of the golden run (plus up to half the patch), else a uniform position in the stack. A patch may
hang over the stack's edge by a quarter of its extent (zero input, zero weight there). Input: the DAPI
channel normalised per gonad with ``infer.robust_range`` (1st and 99.9th percentiles). Augmentation:
flips in y and x, a y/x transpose (the xy grid is isotropic; the distances are unchanged by these) and
a log-uniform gain (0.8 to 1.25). Every sample draws from its own generator seeded with (seed, epoch,
index), so a run is reproducible whatever the batching.
"""
from __future__ import annotations

import json
import logging
import math
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

from germquant.spots.net3d.data import Prefetch, dilate_ellipsoid  # noqa: F401  (Prefetch re-exported)

from .infer import normalize_block

log = logging.getLogger(__name__)

TARGET_VERSION = 1
CLASS_IGNORE, CLASS_BACKGROUND, CLASS_TRACED, CLASS_SHELL = 0, 1, 2, 3
DATA_DEFAULTS = {
    "traced_frac": 0.5,             # Cellpose label with this fraction inside traced masks = the traced nucleus's own
    "own_frac": 0.5,                # Cellpose label holding this fraction of a traced mask is not its neighbour
    "margin_um": 0.5,               # untraced Cellpose labels dilated by this: weight 0
    "w_traced_shell": 1.0,          # weight of a traced nucleus's own Cellpose label outside its traced mask
    "neighbour_radius_um": 1.0,
    "neighbour_power": 3.0,
    "p_traced": 0.6,
    "p_germline": 0.3,
    "max_overhang": 0.25,
    "gain": (0.8, 1.25),
}
TARGET_KEYS = ("traced_frac", "own_frac", "margin_um", "neighbour_radius_um", "neighbour_power")


# ------------------------------------------------------------------------------- helpers
def label_centroids(labels) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(ids, sizes, (n, 3) centroids in voxels) of the labels of a (Z, Y, X) image, plane by plane."""
    lab = labels
    top = 0
    for z in range(lab.shape[0]):
        top = max(top, int(np.asarray(lab[z]).max()))
    size = np.zeros(top + 1, np.int64)
    acc = np.zeros((top + 1, 3), float)
    Y, X = lab.shape[1:]
    yy, xx = np.meshgrid(np.arange(Y, dtype=float), np.arange(X, dtype=float), indexing="ij")
    yy, xx = yy.ravel(), xx.ravel()
    for z in range(lab.shape[0]):
        f = np.asarray(lab[z]).ravel().astype(np.int64, copy=False)
        bc = np.bincount(f, minlength=top + 1)
        size += bc
        acc[:, 0] += bc * float(z)
        acc[:, 1] += np.bincount(f, weights=yy, minlength=top + 1)
        acc[:, 2] += np.bincount(f, weights=xx, minlength=top + 1)
    ids = np.flatnonzero(size[1:]) + 1
    return ids, size[ids], acc[ids] / np.maximum(size[ids], 1)[:, None]


def _grow(sl: tuple, pad, shape) -> tuple:
    return tuple(slice(max(0, s.start - p), min(n, s.stop + p)) for s, p, n in zip(sl, pad, shape))


def cellpose_roles(cellpose, truth, cfg: dict | None = None) -> dict:
    """Which Cellpose labels are traced nuclei's own labels, and the overlap of every (Cellpose label,
    traced label) pair: ``traced_ids`` (Cellpose ids), ``pairs`` {(cellpose id, truth id): voxels},
    ``truth_sizes`` {truth id: voxels}, ``inside_frac`` {cellpose id: fraction inside traced masks}."""
    from germquant.benchmark.nuclei3d import overlap_table

    c = {**DATA_DEFAULTS, **(cfg or {})}
    ot = overlap_table(cellpose, truth)
    inside = np.bincount(ot.pair_pred, weights=ot.pair_inter, minlength=ot.n_pred)
    frac = inside / np.maximum(ot.pred_sizes, 1)
    traced = ot.pred_ids[frac >= float(c["traced_frac"]) - 1e-12]
    pairs = {(int(ot.pred_ids[p]), int(ot.truth_ids[t])): int(v)
             for p, t, v in zip(ot.pair_pred, ot.pair_truth, ot.pair_inter)}
    return {"traced_ids": np.asarray(traced, np.int64), "pairs": pairs,
            "truth_sizes": {int(i): int(s) for i, s in zip(ot.truth_ids, ot.truth_sizes)},
            "inside_frac": {int(i): float(f) for i, f, v in zip(ot.pred_ids, frac, inside) if v > 0},
            "n_cellpose": int(ot.n_pred), "n_truth": int(ot.n_truth)}


# ------------------------------------------------------------------------------- targets
def build_targets(truth, cellpose, spacing, cfg: dict | None = None) -> dict:
    """The class map (uint8), the cell distance and the neighbour closeness (float32) of a gonad, and a
    record of how they were made (module docstring). ``truth`` and ``cellpose`` are (Z, Y, X) label
    images of one shape (arrays or memory maps)."""
    c = {**DATA_DEFAULTS, **(cfg or {})}
    sp = tuple(float(s) for s in spacing)
    shape = tuple(int(s) for s in truth.shape)
    if tuple(cellpose.shape) != shape:
        raise ValueError(f"truth {shape} and Cellpose labels {tuple(cellpose.shape)} differ in shape")
    t0 = time.time()
    roles = cellpose_roles(cellpose, truth, c)
    cp_top = max([0] + [int(np.asarray(cellpose[z]).max()) for z in range(shape[0])])
    lut_traced = np.zeros(cp_top + 1, bool)
    lut_traced[roles["traced_ids"]] = True
    untraced = np.zeros(shape, bool)
    for z in range(shape[0]):
        cz = np.asarray(cellpose[z])
        untraced[z] = (cz > 0) & ~lut_traced[cz]
    ignore = dilate_ellipsoid(untraced, float(c["margin_um"]), sp) if untraced.any() else untraced
    del untraced
    cls = np.full(shape, CLASS_BACKGROUND, np.uint8)
    for z in range(shape[0]):
        plane = cls[z]
        plane[np.asarray(cellpose[z]) > 0] = CLASS_SHELL      # untraced labels are overwritten next
        plane[ignore[z]] = CLASS_IGNORE
        plane[np.asarray(truth[z]) > 0] = CLASS_TRACED
    del ignore

    tvol = np.asarray(truth)
    dist = np.zeros(shape, np.float32)
    nb = np.zeros(shape, np.float32)
    R = float(c["neighbour_radius_um"])
    power = float(c["neighbour_power"])
    pad_nb = tuple(int(math.ceil(R / s)) + 1 for s in sp)
    if tvol.dtype.kind not in "iu":
        raise TypeError(f"the truth label image must hold integers, got {tvol.dtype}")
    objs = ndi.find_objects(tvol if tvol.itemsize <= 4 else tvol.astype(np.int64))
    n_with_neighbour = 0
    for t, sl in enumerate(objs, start=1):
        if sl is None:
            continue
        # cell distance: inside the label, to the nearest voxel outside it (the stack edge counts as outside)
        box = _grow(sl, (1, 1, 1), shape)
        m = tvol[box] == t
        d = ndi.distance_transform_edt(np.pad(m, 1), sampling=sp)[1:-1, 1:-1, 1:-1]
        dist[box][m] = d[m]
        # neighbour closeness: distance to the nearest voxel of another nucleus within the radius
        box = _grow(sl, pad_nb, shape)
        tc = tvol[box]
        cc = np.asarray(cellpose[box])
        m = tc == t
        size_t = roles["truth_sizes"].get(t, int(m.sum()))
        own = [cid for (cid, tid), v in roles["pairs"].items() if tid == t and v >= float(c["own_frac"]) * size_t]
        others = (tc > 0) & ~m
        cp_other = (cc > 0) & ~m
        if own:
            cp_other &= ~np.isin(cc, np.asarray(own, cc.dtype))
        others |= cp_other
        if others.any():
            dn = ndi.distance_transform_edt(~others, sampling=sp)
            val = np.clip(1.0 - dn / R, 0.0, 1.0) ** power
            nb[box][m] = val[m]
            n_with_neighbour += int((val[m] > 0).any())
    counts = np.bincount(cls[:, ::4, ::4].ravel(), minlength=4)
    info = {
        "version": TARGET_VERSION, "n_traced": int(roles["n_truth"]), "n_cellpose": int(roles["n_cellpose"]),
        "n_cellpose_traced": int(len(roles["traced_ids"])),
        "n_cellpose_touching_traced": int(len(roles["inside_frac"])),
        "n_traced_with_neighbour_within_radius": int(n_with_neighbour),
        "class_fraction_sampled": {name: float(counts[k] / max(counts.sum(), 1)) for k, name in
                                   ((CLASS_IGNORE, "ignore"), (CLASS_BACKGROUND, "background"),
                                    (CLASS_TRACED, "traced"), (CLASS_SHELL, "shell"))},
        "max_cell_distance_um": float(dist.max()) if dist.size else 0.0,
        "build_s": round(time.time() - t0, 1),
        "settings": {k: c[k] for k in TARGET_KEYS},
    }
    return {"cls": cls, "dist": dist, "nb": nb, "info": info}


def targets_key(tag: str, meta: dict, truth_labels_path, cfg: dict) -> dict:
    """What a cached set of target maps was made from (a mismatch recomputes them)."""
    c = {**DATA_DEFAULTS, **(cfg or {})}
    tl = Path(truth_labels_path)
    st = tl.stat() if tl.is_file() else None
    return {"version": TARGET_VERSION, "tag": tag,
            "labels_sha256": ((meta.get("arrays") or {}).get("labels") or {}).get("data_sha256"),
            "truth_labels": str(tl), "truth_labels_size": st.st_size if st else None,
            "truth_labels_mtime": st.st_mtime if st else None,
            **{k: (list(c[k]) if isinstance(c[k], tuple) else c[k]) for k in TARGET_KEYS}}


def load_targets(prep_dir, tag: str, meta: dict, cellpose, truth_labels_path, spacing, cfg: dict) -> dict:
    """The target maps of a gonad from the prep folder (memory maps), computed and saved when absent or
    stale, plus the record and the traced nuclei's centroids (voxels)."""
    from germquant.benchmark.manifest import load_labels

    prep = Path(prep_dir)
    prep.mkdir(parents=True, exist_ok=True)
    files = {k: prep / f"{tag}__n3d_{k}.npy" for k in ("cls", "dist", "nb")}
    js = prep / f"{tag}__n3d_targets.json"
    key = targets_key(tag, meta, truth_labels_path, cfg)
    if js.is_file() and all(f.is_file() for f in files.values()):
        try:
            doc = json.loads(js.read_text(encoding="utf-8"))
            if doc.get("key") == key:
                return {**{k: np.load(f, mmap_mode="r") for k, f in files.items()}, "info": doc["info"],
                        "traced_centers": np.asarray(doc["traced_centers"], float).reshape(-1, 3)}
        except (OSError, ValueError, KeyError):
            pass
    truth = load_labels(truth_labels_path)
    built = build_targets(truth, cellpose, spacing, cfg)
    _ids, _sizes, centers = label_centroids(truth)
    del truth
    for k, f in files.items():
        arr = built[k] if k == "cls" else built[k].astype(np.float16)
        tmp = f.with_name(f.stem + ".partial.npy")
        np.save(tmp, arr)
        tmp.replace(f)
    doc = {"key": key, "info": built["info"], "traced_centers": centers.round(3).tolist()}
    tmp = js.with_name(js.name + ".partial")
    tmp.write_text(json.dumps(doc, indent=1), encoding="utf-8")
    tmp.replace(js)
    del built
    return {**{k: np.load(f, mmap_mode="r") for k, f in files.items()}, "info": doc["info"],
            "traced_centers": centers}


# --------------------------------------------------------------------------------- a gonad
@dataclass
class TrainGonad:
    tag: str
    spacing: tuple[float, float, float]
    dna: np.ndarray                   # raw (Z, Y, X), a memory map
    lo_hi: tuple[float, float]
    cls: np.ndarray                   # uint8 class map
    dist: np.ndarray                  # cell distance (um), float16
    nb: np.ndarray                    # neighbour closeness, float16
    traced_centers: np.ndarray        # (n, 3) voxels
    germline_centers: np.ndarray      # (m, 3) voxels
    condition: str = ""
    sex: str = ""
    info: dict = field(default_factory=dict)

    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(int(s) for s in self.cls.shape)


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


def class_weights(cfg: dict | None = None) -> np.ndarray:
    """The loss weight of each class (ignore, background, traced, shell)."""
    c = {**DATA_DEFAULTS, **(cfg or {})}
    return np.array([0.0, 1.0, 1.0, float(c["w_traced_shell"])], np.float32)


def _flip(a: np.ndarray, flips) -> np.ndarray:
    """Flip y, flip x, transpose y/x of a (C, Z, Y, X) array as ``flips`` (three 0/1 flags) says."""
    if flips[0]:
        a = a[:, :, ::-1, :]
    if flips[1]:
        a = a[:, :, :, ::-1]
    if flips[2]:
        a = a.transpose(0, 1, 3, 2)
    return np.ascontiguousarray(a)


def make_patch(g: TrainGonad, origin, size, cfg: dict | None = None, rng: np.random.Generator | None = None,
               augment: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(input (1, z, y, x), targets (3, z, y, x): foreground, cell distance, neighbour closeness, weight
    (1, z, y, x)) float32 for a patch at ``origin``."""
    c = {**DATA_DEFAULTS, **(cfg or {})}
    size = tuple(int(s) for s in size)
    x = np.zeros((1,) + size, np.float32)
    cls = np.zeros(size, np.uint8)                      # outside the stack: class ignore, weight 0
    tgt = np.zeros((3,) + size, np.float32)
    vol, pat = _window(origin, size, g.shape)
    if vol is not None:
        x[(0,) + pat] = normalize_block(g.dna[vol], g.lo_hi[0], g.lo_hi[1])
        cls[pat] = g.cls[vol]
        tgt[(1,) + pat] = g.dist[vol]
        tgt[(2,) + pat] = g.nb[vol]
    tgt[0] = cls == CLASS_TRACED
    w = class_weights(c)[cls][None]
    if augment and rng is not None:
        lo, hi = c["gain"]
        x *= np.float32(math.exp(rng.uniform(math.log(lo), math.log(hi))))
        flips = rng.integers(0, 2, size=3)
        x, tgt, w = _flip(x, flips), _flip(tgt, flips), _flip(w, flips)
    return x, tgt, w.astype(np.float32)


def sample_origin(g: TrainGonad, size, cfg: dict, rng: np.random.Generator) -> tuple[int, int, int]:
    """A patch origin: around a traced nucleus, a germline nucleus or anywhere (module docstring)."""
    c = {**DATA_DEFAULTS, **cfg}
    size = np.asarray(size, int)
    u = rng.random()
    if len(g.traced_centers) and u < float(c["p_traced"]):
        centre = g.traced_centers[rng.integers(len(g.traced_centers))] + rng.uniform(-size / 3.0, size / 3.0)
    elif len(g.germline_centers) and u < float(c["p_traced"]) + float(c["p_germline"]):
        centre = g.germline_centers[rng.integers(len(g.germline_centers))] + rng.uniform(-size / 2.0, size / 2.0)
    else:
        centre = rng.uniform(0, np.asarray(g.shape, float))
    origin = np.rint(centre - size / 2.0).astype(int)
    hang = np.floor(size * float(c["max_overhang"])).astype(int)
    lo = -hang
    hi = np.maximum(np.asarray(g.shape) - size + hang, lo)
    return tuple(int(v) for v in np.clip(origin, lo, hi))


def sample(gonads: list[TrainGonad], size, cfg: dict, seed: int, epoch: int, index: int):
    """One training sample (input, targets, weight), reproducible from (seed, epoch, index)."""
    rng = np.random.default_rng([int(seed), int(epoch), int(index)])
    g = gonads[int(rng.integers(len(gonads)))]
    origin = sample_origin(g, size, cfg, rng)
    return make_patch(g, origin, size, cfg, rng, augment=True)


def batches(gonads, size, cfg: dict, seed: int, epoch: int, n_steps: int, batch_size: int):
    """The batches of one epoch, in order: (inputs (B, 1, ...), targets (B, 3, ...), weights (B, 1, ...))."""
    for step in range(int(n_steps)):
        items = [sample(gonads, size, cfg, seed, epoch, step * batch_size + b) for b in range(batch_size)]
        yield (np.stack([i[0] for i in items]), np.stack([i[1] for i in items]), np.stack([i[2] for i in items]))
