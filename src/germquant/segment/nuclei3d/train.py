"""Training of the nuclei3d model (``germquant train-nuclei3d``).

Which gonads. The roles come from the split file (``tests/golden/rad51_split.json``), checked against the
pre-registration (``germquant.benchmark.detector.gonad_role``: a split file that moves a registered gonad
is refused) and the manifest. Training uses every train gonad with a truth label image: for the N2 RAD-51
split all six, n2rad51_HS_herm_000, HS_male_000, HS_male_002, noHS_herm_000, noHS_male_000 and
noHS_male_001 (the two heat-shock males were exported with Imaris spot region growing, which ruins their
spot truth but leaves their nucleus surfaces intact). Early stopping and model selection use the val
gonads, n2rad51_HS_herm_001 and noHS_herm_001. The test gonads are named in the log and never opened:
every gonad this module loads passes ``guard_test`` first.

What is learned (data.py): the foreground, the cell distance (um) and the neighbour closeness inside the
traced Imaris masks, background where neither a traced mask nor a deployed Cellpose label is present and
in a traced nucleus's own Cellpose label outside its mask, nothing inside untraced Cellpose labels plus a
0.5 um margin. The loss is, per voxel, the binary cross-entropy of the foreground logit plus
``lambda_dist`` times the squared error of the cell distance plus ``lambda_nb`` times the squared error of
the neighbour closeness, times the voxel's weight, summed and divided by the summed weight of the batch.

How. AdamW, a linear warm-up then a cosine decay of the learning rate, bf16 autocast on CUDA (float32 on
the CPU), seeds for python, numpy and torch, ``torch.use_deterministic_algorithms(True)`` with cuDNN
benchmarking off and ``CUBLAS_WORKSPACE_CONFIG`` set (the net3d settings, reused), so two runs with the
same seed on the same machine give the same weights.

Validation (every ``val_every`` epochs) and the early stopping. Each val gonad goes through the backend's
own path: tiled float32 inference of the whole stack (infer.py), the instance step with the fixed
parameters ``post`` (foreground 0.5, h 0.4 um, neighbour weight 1), the pipeline's size filter
(``segment.nuclei._drop_small`` with ``min_volume_um3``, 4 um3 as in config/config.yaml), then
``germquant.benchmark.nuclei3d.benchmark_gonad`` with the pre-registered rule (``coverage="protocol"``,
docs/BENCHMARK_PROTOCOL.md section 5.3): truth fragments under 1000 voxels removed, predictions in scope
when at least half of their volume lies within 1 um (xy) and one plane (z) of a truth mask. The two val
gonads are pooled by summing TP, FP and FN at every IoU threshold. Only the foreground components that
touch region C go through the instance step (``instances_from_seed_map(keep=...)``): the others can only
produce predictions the protocol ignores, so every metric is what the whole stack would give, while the
counts of all predictions (``n_pred``, ``n_pred_ignored_outside_coverage``) cover the kept components only
(scripts/benchmark_nuclei3d.py segments the whole stack). Selection: the best pooled F1 at IoU
0.5, then the best pooled AP (IoU 0.5 to 0.95), then the earlier epoch. Early stopping: ``patience``
validations in a row without a new best. Every epoch's weights are kept (``checkpoints/epoch_NNN.pt``,
sha256 in the log); ``model.pt`` holds the selected epoch's weights.

The optional post grid (``post_grid``): for the selected epoch, the instance-step parameters are scored on
the same val gonads over a small grid (the maps are predicted once per gonad), chosen by the same rule,
and written into ``model.pt``. Every configuration scored on val (each validated epoch, each grid point)
is listed in the log; together they are the configurations the protocol's tuning budget counts.
"""
from __future__ import annotations

import itertools
import json
import logging
import math
import shutil
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from germquant.spots.net3d.train import _lr_at, _utc, _versions, _write_json, set_determinism, wait_for_gpu

from . import data as d3
from .infer import (DEFAULT_TILE_ZYX, NORM_DEFAULTS, POST_DEFAULTS, instances, instances_from_seed_map,
                    predict_maps, predict_seed_map, read_model_file, robust_range, save_model_file)

log = logging.getLogger(__name__)

TRAIN_DEFAULTS = {
    "epochs": 40,
    "patch": (48, 160, 160),
    "seed": 0,
    "batch_size": 4,
    "steps_per_epoch": 250,
    "lr": 1e-3,
    "weight_decay": 1e-4,
    "warmup_steps": 250,
    "min_lr_frac": 0.05,
    "patience": 8,
    "val_every": 1,
    "base_channels": 16,
    "levels": 4,
    "lambda_dist": 1.0,
    "lambda_nb": 1.0,
    "min_volume_um3": 4.0,
}
SUMMARY_KEYS = ("f1_iou50", "precision_iou50", "recall_iou50", "ap_iou50_95", "tp_iou50", "fp_iou50", "fn_iou50",
                "n_pred", "n_pred_counted", "n_pred_ignored_outside_coverage", "n_truth", "n_splits", "n_merges",
                "mean_matched_iou_iou50", "z_extent_error_mean_um", "z_extent_error_mean_abs_um", "n_z_pairs")


def _finite(x):
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _num(x, default: float = -math.inf) -> float:
    f = _finite(x)
    return default if f is None else f


# ------------------------------------------------------------------------- gonads and roles
def resolve_gonads(split_file, manifest) -> dict:
    """The train and val gonads with a truth label image, what was left out and why, and the test
    gonads (named, never opened)."""
    from germquant.benchmark import detector as bd
    from germquant.benchmark.manifest import load_manifest, truth_file_names

    split = json.loads(Path(split_file).read_text(encoding="utf-8"))
    roles_doc = split.get("roles")
    if not isinstance(roles_doc, dict):
        raise ValueError(f"{split_file}: the split file needs a 'roles' mapping (train, val, test)")
    entries = {e.tag: e for e in load_manifest(manifest)}
    roles = {}
    for role, tags in roles_doc.items():
        for t in tags:
            roles[str(t)] = bd.gonad_role(str(t), role)     # refuses a split that moves a registered gonad
    out = {"train": [], "val": [], "test": [], "excluded": [], "roles": roles, "entries": entries}
    for t, role in roles.items():
        if role == "test":
            out["test"].append(t)
            continue
        e = entries.get(t)
        if e is None:
            out["excluded"].append({"tag": t, "role": role, "reason": "not in the manifest"})
            continue
        man_role = str(e.extra.get("role") or e.split or "").strip()
        if man_role in ("train", "val", "test") and man_role != role:
            raise ValueError(f"{t}: the split file says {role!r}, the manifest says {man_role!r}")
        labels = Path(e.truth_dir or "") / truth_file_names(e.truth_stem)["labels"]
        if e.truth_dir is None or not labels.is_file():
            out["excluded"].append({"tag": t, "role": role, "reason": f"no truth label image ({labels})"})
            continue
        if role in ("train", "val"):
            out[role].append(t)
    bd.guard_test(out["train"] + out["val"], allow_test=False, roles=roles)
    return out


def _gonad_facts(entry, cache_meta: dict) -> dict:
    return {"condition": str(entry.extra.get("condition") or cache_meta.get("condition") or ""),
            "sex": str(entry.extra.get("sex") or cache_meta.get("sex") or "")}


def load_train_gonad(tag: str, entry, cache_root, prep_dir, *, data_cfg: dict, norm: dict) -> d3.TrainGonad:
    from germquant.benchmark import detector as bd
    from germquant.benchmark.manifest import load_truth

    bd.guard_test([tag], allow_test=False, roles={tag: bd.entry_role(entry)})
    cache = bd.load_cache(tag, cache_root, arrays=("dna", "labels"), mmap=True)
    if cache.dna is None:
        raise FileNotFoundError(f"{tag}: the bench cache holds no DNA channel")
    truth = load_truth(entry.truth_dir, entry.truth_stem)
    if truth is None or not truth.has_labels:
        raise FileNotFoundError(f"{tag}: no truth label image for stem {entry.truth_stem!r} in {entry.truth_dir}")
    t0 = time.time()
    tg = d3.load_targets(prep_dir, tag, cache.meta, cache.labels, truth.labels_path, cache.spacing, data_cfg)
    sp = np.asarray(cache.spacing, float)
    nuc = cache.nuclei
    germ = nuc[nuc["in_germline"]] if "in_germline" in nuc.columns else nuc
    centers = (germ[["centroid_z_um", "centroid_y_um", "centroid_x_um"]].to_numpy(float) / sp
               if len(germ) else np.zeros((0, 3)))
    lo_hi = robust_range(cache.dna, norm["low_pct"], norm["high_pct"], norm["stride_yx"])
    g = d3.TrainGonad(tag=tag, spacing=cache.spacing, dna=cache.dna, lo_hi=lo_hi, cls=tg["cls"], dist=tg["dist"],
                      nb=tg["nb"], traced_centers=tg["traced_centers"], germline_centers=centers,
                      **_gonad_facts(entry, cache.meta))
    g.info = {**tg["info"], "n_germline_nuclei": int(len(centers)), "norm_range": list(lo_hi),
              "prep_s": round(time.time() - t0, 1),
              "labels_sha256": ((cache.meta.get("arrays") or {}).get("labels") or {}).get("data_sha256"),
              "dna_sha256": ((cache.meta.get("arrays") or {}).get("dna") or {}).get("data_sha256"),
              "truth_labels": str(truth.labels_path)}
    return g


@dataclass
class ValGonad:
    tag: str
    split: str
    dna: object
    lo_hi: tuple
    truth: np.ndarray                 # the truth labels with the protocol's fragments removed
    spacing: tuple
    region: np.ndarray | None = None  # the protocol's region C of that truth
    condition: str = ""
    sex: str = ""
    n_truth: int = 0
    n_fragments: int = 0


def load_val_gonad(tag: str, entry, cache_root, *, norm: dict) -> ValGonad:
    from germquant.benchmark import detector as bd
    from germquant.benchmark.manifest import load_labels, load_truth

    bd.guard_test([tag], allow_test=False, roles={tag: bd.entry_role(entry)})
    cache = bd.load_cache(tag, cache_root, arrays=("dna",), mmap=True)
    if cache.dna is None:
        raise FileNotFoundError(f"{tag}: the bench cache holds no DNA channel")
    truth = load_truth(entry.truth_dir, entry.truth_stem)
    if truth is None or not truth.has_labels:
        raise FileNotFoundError(f"{tag}: no truth label image for stem {entry.truth_stem!r} in {entry.truth_dir}")
    from germquant.benchmark import nuclei3d as n3

    tl = load_labels(truth.labels_path)
    if tuple(tl.shape) != tuple(cache.dna.shape):
        raise ValueError(f"{tag}: truth labels {tl.shape} and DAPI {cache.dna.shape} differ in shape")
    tl, fragments = n3.remove_truth_fragments(tl)
    region = n3.protocol_coverage_mask(tl, cache.spacing)
    lo_hi = robust_range(cache.dna, norm["low_pct"], norm["high_pct"], norm["stride_yx"])
    facts = _gonad_facts(entry, cache.meta)
    return ValGonad(tag=tag, split=facts["condition"], dna=cache.dna, lo_hi=lo_hi, truth=tl, spacing=cache.spacing,
                    region=region, n_truth=int(len(d3.label_centroids(tl)[0])), n_fragments=len(fragments), **facts)


# ------------------------------------------------------------------------------- validation
def score_labels(vg: ValGonad, labels: np.ndarray) -> tuple[dict, object]:
    """The protocol's per-gonad row and per-truth table of one label image (the size filter applied)."""
    from germquant.benchmark import nuclei3d as n3

    return n3.benchmark_gonad(labels, vg.truth, vg.spacing, tag=vg.tag, split=vg.split, coverage="protocol",
                              protocol_region=vg.region)


def pool_rows(rows: list[dict]) -> dict:
    """Pooled val scores from summed counts: F1, precision and recall at IoU 0.5, AP over 0.5 to 0.95,
    splits, merges, and the z-extent error averaged over the matched pairs."""
    from germquant.benchmark import nuclei3d as n3
    from germquant.benchmark.stats import precision_recall_f1

    tp = sum(int(r["tp_iou50"]) for r in rows)
    counted = sum(int(r["n_pred_counted"]) for r in rows)
    n_truth = sum(int(r["n_truth"]) for r in rows)
    s = precision_recall_f1(tp, counted, n_truth)
    terms = []
    for t in n3.IOU_THRESHOLDS:
        k = n3._tag(t)
        a = sum(int(r[f"tp_{k}"]) for r in rows)
        b = sum(int(r[f"fp_{k}"]) for r in rows)
        c = sum(int(r[f"fn_{k}"]) for r in rows)
        terms.append(a / (a + b + c) if a + b + c else float("nan"))
    nz = sum(int(r["n_z_pairs"]) for r in rows)
    zsum = sum(float(r["z_extent_error_mean_um"]) * int(r["n_z_pairs"]) for r in rows if int(r["n_z_pairs"]))
    zabs = sum(float(r["z_extent_error_mean_abs_um"]) * int(r["n_z_pairs"]) for r in rows if int(r["n_z_pairs"]))
    miou = [float(r["mean_matched_iou_iou50"]) * int(r["tp_iou50"]) for r in rows if int(r["tp_iou50"])]
    return {"f1_iou50": s["f1"], "precision_iou50": s["precision"], "recall_iou50": s["recall"],
            "ap_iou50_95": float(np.nanmean(terms)) if np.isfinite(terms).any() else float("nan"),
            "tp_iou50": tp, "fp_iou50": counted - tp, "fn_iou50": n_truth - tp, "n_truth": n_truth,
            "n_pred": sum(int(r["n_pred"]) for r in rows), "n_pred_counted": counted,
            "n_pred_ignored_outside_coverage": sum(int(r["n_pred_ignored_outside_coverage"]) for r in rows),
            "n_splits": sum(int(r["n_splits"]) for r in rows), "n_merges": sum(int(r["n_merges"]) for r in rows),
            "mean_matched_iou_iou50": (sum(miou) / tp) if tp else float("nan"),
            "z_extent_error_mean_um": (zsum / nz) if nz else float("nan"),
            "z_extent_error_mean_abs_um": (zabs / nz) if nz else float("nan"), "n_z_pairs": nz}


def selection_key(entry: dict):
    """Higher is better: pooled F1 at IoU 0.5, then pooled AP, then the earlier epoch."""
    return (_num(entry.get("f1_iou50")), _num(entry.get("ap_iou50_95")), -int(entry.get("epoch") or 0))


def _summ(row: dict) -> dict:
    return {k: (_finite(row.get(k)) if isinstance(row.get(k), (float, np.floating)) else row.get(k))
            for k in SUMMARY_KEYS}


def evaluate(model, vals: list[ValGonad], *, device: str = "cpu", tile_zyx=DEFAULT_TILE_ZYX,
             post: dict | None = None, min_volume_um3: float = TRAIN_DEFAULTS["min_volume_um3"]) -> dict:
    """Score the model on the val gonads with one set of instance-step parameters (module docstring)."""
    from germquant.segment.nuclei import _drop_small

    p = {**POST_DEFAULTS, **(post or {})}
    rows, timing = [], {}
    for vg in vals:
        t0 = time.time()
        fg, s = predict_seed_map(model, vg.dna, vg.lo_hi, device=device, tile_zyx=tile_zyx,
                                 fg_threshold=p["fg_threshold"], neighbour_weight=p["neighbour_weight"])
        t1 = time.time()
        fg_frac = float(fg.mean())
        labels = instances_from_seed_map(fg, s, h_um=p["h_um"], keep=vg.region)
        del fg, s
        labels = _drop_small(labels, vg.spacing, min_volume_um3)
        t2 = time.time()
        row, _pt = score_labels(vg, labels)
        del labels
        rows.append(row)
        timing[vg.tag] = {"infer_s": round(t1 - t0, 1), "instances_s": round(t2 - t1, 1),
                          "score_s": round(time.time() - t2, 1), "foreground_fraction": round(fg_frac, 5)}
    pooled = pool_rows(rows)
    return {"post": p, "pooled": pooled, "per_gonad": {r["tag"]: _summ(r) for r in rows}, "timing": timing}


def post_grid(model, vals: list[ValGonad], grid: dict, *, device: str = "cpu", tile_zyx=DEFAULT_TILE_ZYX,
              min_volume_um3: float = TRAIN_DEFAULTS["min_volume_um3"]) -> list[dict]:
    """Every combination of the instance-step parameters in ``grid`` ({name: [values]}, names from
    ``POST_DEFAULTS``) scored on the val gonads; the maps are predicted once per gonad."""
    from germquant.segment.nuclei import _drop_small

    names = [k for k in POST_DEFAULTS if k in grid]
    unknown = sorted(set(grid) - set(POST_DEFAULTS))
    if unknown:
        raise ValueError(f"unknown instance-step parameters {unknown} (known: {sorted(POST_DEFAULTS)})")
    combos = [dict(zip(names, vals_)) for vals_ in itertools.product(*[list(grid[k]) for k in names])]
    combos = [{**POST_DEFAULTS, **c} for c in combos]
    rows = {i: [] for i in range(len(combos))}
    for vg in vals:
        maps = predict_maps(model, vg.dna, vg.lo_hi, device=device, tile_zyx=tile_zyx)
        for i, c in enumerate(combos):
            labels = instances(maps, fg_threshold=c["fg_threshold"], h_um=c["h_um"],
                               neighbour_weight=c["neighbour_weight"], keep=vg.region)
            labels = _drop_small(labels, vg.spacing, min_volume_um3)
            row, _pt = score_labels(vg, labels)
            rows[i].append(row)
            del labels
        del maps
    return [{"post": c, "pooled": pool_rows(rows[i]), "per_gonad": {r["tag"]: _summ(r) for r in rows[i]}}
            for i, c in enumerate(combos)]


# ------------------------------------------------------------------------------------ train
def nuclei3d_loss(outputs, targets, weights, *, lambda_dist: float = TRAIN_DEFAULTS["lambda_dist"],
                  lambda_nb: float = TRAIN_DEFAULTS["lambda_nb"]):
    """(loss, {"bce", "dist", "nb": the weighted mean of each term, detached}) for raw network outputs
    (B, 3, Z, Y, X), targets (B, 3, ...) (foreground 0/1, cell distance um, neighbour closeness) and
    per-voxel weights (B, Z, Y, X) or (B, 1, Z, Y, X). Per voxel: the binary cross-entropy of the
    foreground logit plus ``lambda_dist`` times the squared error of the cell distance plus ``lambda_nb``
    times that of the neighbour closeness, times the weight; summed and divided by the summed weight (at
    least 1). A voxel of weight 0 contributes nothing, whatever the network says there."""
    import torch
    import torch.nn.functional as F

    out = outputs.float()
    w = weights[:, 0] if weights.dim() == out.dim() else weights
    w = w.float()
    bce = F.binary_cross_entropy_with_logits(out[:, 0], targets[:, 0], reduction="none")
    l_d = (out[:, 1] - targets[:, 1]) ** 2
    l_n = (out[:, 2] - targets[:, 2]) ** 2
    wsum = w.sum().clamp_min(1.0)
    loss = ((bce + float(lambda_dist) * l_d + float(lambda_nb) * l_n) * w).sum() / wsum
    with torch.no_grad():
        terms = {"bce": (bce * w).sum() / wsum, "dist": (l_d * w).sum() / wsum, "nb": (l_n * w).sum() / wsum}
    return loss, {k: v.detach() for k, v in terms.items()}


def train_nuclei3d(cache_root, split_file, out_dir, **kwargs) -> dict:
    """Train nuclei3d on the bench cache (module docstring) into ``out_dir``: ``model.pt`` (the weights of
    the selected epoch, with the instance-step parameters), ``checkpoints/epoch_NNN.pt``,
    ``training_log.json`` and ``model_card.json``. Returns the paths and the selected validation. The
    keyword arguments are those of ``_train``. torch's deterministic-algorithm and cuDNN settings are
    restored afterwards."""
    import torch

    saved = (torch.are_deterministic_algorithms_enabled(), torch.backends.cudnn.deterministic,
             torch.backends.cudnn.benchmark)
    try:
        return _train(cache_root, split_file, out_dir, **kwargs)
    finally:
        torch.use_deterministic_algorithms(saved[0])
        torch.backends.cudnn.deterministic, torch.backends.cudnn.benchmark = saved[1], saved[2]


def _train(cache_root, split_file, out_dir, *, manifest=None, epochs: int = TRAIN_DEFAULTS["epochs"],
           patch=TRAIN_DEFAULTS["patch"], seed: int = TRAIN_DEFAULTS["seed"],
           batch_size: int = TRAIN_DEFAULTS["batch_size"], steps_per_epoch: int = TRAIN_DEFAULTS["steps_per_epoch"],
           lr: float = TRAIN_DEFAULTS["lr"], weight_decay: float = TRAIN_DEFAULTS["weight_decay"],
           warmup_steps: int = TRAIN_DEFAULTS["warmup_steps"], patience: int = TRAIN_DEFAULTS["patience"],
           val_every: int = TRAIN_DEFAULTS["val_every"], device: str = "auto",
           base_channels: int = TRAIN_DEFAULTS["base_channels"], levels: int = TRAIN_DEFAULTS["levels"],
           lambda_dist: float = TRAIN_DEFAULTS["lambda_dist"], lambda_nb: float = TRAIN_DEFAULTS["lambda_nb"],
           tile_zyx=DEFAULT_TILE_ZYX, data_cfg: dict | None = None, norm: dict | None = None,
           post: dict | None = None, grid: dict | None = None,
           min_volume_um3: float = TRAIN_DEFAULTS["min_volume_um3"], prep_dir=None, min_free_gb: float = 0.0,
           wait_minutes: float = 120.0, name: str | None = None) -> dict:
    import torch

    from germquant import provenance
    from germquant.benchmark import detector as bd

    from .model import HEADS, build_model, n_parameters, size_divisor

    t_start = time.time()
    started = _utc()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "checkpoints").mkdir(exist_ok=True)
    manifest = Path(manifest) if manifest else bd.DEFAULT_MANIFEST
    patch = tuple(int(v) for v in patch)
    div = size_divisor(levels)
    if len(patch) != 3 or any(p <= 0 or p % d for p, d in zip(patch, div)):
        raise ValueError(f"patch {patch} must be three positive multiples of {div} (the network's size divisor)")
    if patch[1] != patch[2]:
        raise ValueError(f"patch {patch}: Y and X must be equal (the y/x transpose augmentation swaps them)")
    dcfg = {**d3.DATA_DEFAULTS, **(data_cfg or {})}
    ncfg = {**NORM_DEFAULTS, **(norm or {})}
    pcfg = {**POST_DEFAULTS, **(post or {})}

    # the device and determinism come first: CUBLAS_WORKSPACE_CONFIG must be set before CUDA starts
    want = str(device or "auto").lower()
    dev = "cuda" if (want in ("auto", "cuda") and torch.cuda.is_available()) else "cpu"
    if want == "cuda" and dev != "cuda":
        log.warning("CUDA requested but not available; training on the CPU")
    set_determinism(int(seed))
    if dev == "cuda" and min_free_gb > 0:
        log.info("GPU: %.1f GB free", wait_for_gpu(float(min_free_gb), max_minutes=float(wait_minutes)))

    gon = resolve_gonads(split_file, manifest)
    if not gon["train"]:
        raise ValueError("no train gonad with a truth label image")
    if not gon["val"]:
        raise ValueError("no val gonad: early stopping needs one")
    prep = Path(prep_dir) if prep_dir else out / "prep"
    log.info("train gonads %s; val %s; left out %s; test (never opened) %s", gon["train"], gon["val"],
             [e["tag"] for e in gon["excluded"]], gon["test"])
    t0 = time.time()
    train_g = [load_train_gonad(t, gon["entries"][t], cache_root, prep, data_cfg=dcfg, norm=ncfg) for t in gon["train"]]
    val_g = [load_val_gonad(t, gon["entries"][t], cache_root, norm=ncfg) for t in gon["val"]]
    prep_s = time.time() - t0
    log.info("data ready in %.0f s", prep_s)

    model = build_model({"in_channels": 1, "base_channels": base_channels, "levels": levels}).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    total_steps = int(epochs) * int(steps_per_epoch)
    git_sha = provenance.git_sha()
    spacing = [float(s) for s in train_g[0].spacing]
    settings = {"epochs": int(epochs), "patch": list(patch), "seed": int(seed), "batch_size": int(batch_size),
                "steps_per_epoch": int(steps_per_epoch), "lr": float(lr), "weight_decay": float(weight_decay),
                "warmup_steps": int(warmup_steps), "min_lr_frac": TRAIN_DEFAULTS["min_lr_frac"],
                "patience": int(patience), "val_every": int(val_every), "device": dev,
                "amp": "bf16 autocast" if dev == "cuda" else "float32", "deterministic": True,
                "lambda_dist": float(lambda_dist), "lambda_nb": float(lambda_nb), "tile_zyx": list(tile_zyx),
                "data": {k: (list(v) if isinstance(v, tuple) else v) for k, v in dcfg.items()}, "norm": ncfg,
                "post": pcfg, "grid": grid or None, "min_volume_um3": float(min_volume_um3),
                "scoring": "germquant.benchmark.nuclei3d.benchmark_gonad, coverage='protocol' (protocol section 5.3)",
                "heads": list(HEADS), "arch": model.arch(), "n_parameters": n_parameters(model),
                "spacing_zyx_um": spacing}
    gonad_rec = {"train": [{"tag": g.tag, "condition": g.condition, "sex": g.sex, **g.info} for g in train_g],
                 "val": [{"tag": v.tag, "condition": v.condition, "sex": v.sex, "n_truth": v.n_truth,
                          "n_truth_fragments_removed": v.n_fragments, "norm_range": list(v.lo_hi),
                          "region_c_fraction": float(v.region.mean()) if v.region is not None else None}
                         for v in val_g],
                 "excluded": gon["excluded"], "test_never_opened": gon["test"]}
    meta_base = {"input_roles": ["dna"], "norm": ncfg, "heads": list(HEADS), "post": pcfg,
                 "spacing_zyx_um": spacing, "min_volume_um3_at_validation": float(min_volume_um3)}

    epochs_log, checkpoints, entries = [], [], []
    best = None
    since_best = 0
    stop_reason = "max epochs reached"
    step = 0
    log.info("training on %s: %d parameters, %d epochs x %d steps x %d patches %s", dev, n_parameters(model),
             epochs, steps_per_epoch, batch_size, patch)
    for epoch in range(int(epochs)):
        model.train()
        t_ep = time.time()
        parts = {"loss": [], "bce": [], "dist": [], "nb": []}
        if dev == "cuda":
            torch.cuda.reset_peak_memory_stats()
        for xb, yb, wb in d3.Prefetch(d3.batches(train_g, patch, dcfg, int(seed), epoch, steps_per_epoch,
                                                 batch_size), depth=4):
            for gparam in opt.param_groups:
                gparam["lr"] = _lr_at(step, total_steps, lr, warmup_steps, TRAIN_DEFAULTS["min_lr_frac"])
            x = torch.from_numpy(xb).to(dev, non_blocking=True)
            y = torch.from_numpy(yb).to(dev, non_blocking=True)
            w = torch.from_numpy(wb).to(dev, non_blocking=True)[:, 0]
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=(dev == "cuda")):
                outp = model(x)
            loss, terms = nuclei3d_loss(outp, y, w, lambda_dist=lambda_dist, lambda_nb=lambda_nb)
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()
            parts["loss"].append(float(loss.detach()))
            for k, v in terms.items():
                parts[k].append(float(v))
            step += 1
        rec = {"epoch": epoch + 1, "lr_end": opt.param_groups[0]["lr"], "train_s": round(time.time() - t_ep, 1),
               **{f"train_{k}": float(np.mean(v)) for k, v in parts.items()}}
        if dev == "cuda":
            rec["gpu_peak_gb"] = round(torch.cuda.max_memory_allocated() / 1e9, 2)
            rec["gpu_reserved_gb"] = round(torch.cuda.max_memory_reserved() / 1e9, 2)
        ckpt = out / "checkpoints" / f"epoch_{epoch + 1:03d}.pt"
        sha_ep = save_model_file(ckpt, model, {**meta_base, "epoch": epoch + 1, "git_sha": git_sha})
        rec["checkpoint"], rec["checkpoint_sha256"] = f"checkpoints/{ckpt.name}", sha_ep
        checkpoints.append({"epoch": epoch + 1, "file": f"checkpoints/{ckpt.name}", "sha256": sha_ep})
        if (epoch + 1) % max(1, int(val_every)) == 0 or epoch + 1 == int(epochs):
            t_v = time.time()
            if dev == "cuda":
                torch.cuda.empty_cache()     # the card is shared: hand the training blocks back before inference
            ev = evaluate(model, val_g, device=dev, tile_zyx=tile_zyx, post=pcfg, min_volume_um3=min_volume_um3)
            if dev == "cuda":
                rec["gpu_val_reserved_gb"] = round(torch.cuda.max_memory_reserved() / 1e9, 2)
                torch.cuda.empty_cache()
                torch.cuda.reset_peak_memory_stats()
            rec["val"] = ev
            rec["val_s"] = round(time.time() - t_v, 1)
            entry = {"epoch": epoch + 1, **{k: ev["pooled"][k] for k in ("f1_iou50", "ap_iou50_95", "recall_iou50",
                                                                          "precision_iou50")}}
            entries.append(entry)
            if best is None or selection_key(entry) > selection_key(best):
                best = entry
                since_best = 0
            else:
                since_best += 1
            rec["selection"] = {"epoch": best["epoch"], "f1_iou50": best["f1_iou50"],
                                "ap_iou50_95": best["ap_iou50_95"], "validations_without_new_best": since_best}
            pv = ev["pooled"]
            log.info("epoch %d: loss %.4f (bce %.4f, dist %.4f, nb %.4f) | val F1@0.5 %.3f AP %.3f recall %.3f "
                     "(TP %d FP %d FN %d, splits %d, merges %d) | best epoch %d | %.0f s train, %.0f s val",
                     epoch + 1, rec["train_loss"], rec["train_bce"], rec["train_dist"], rec["train_nb"],
                     _num(pv["f1_iou50"], math.nan), _num(pv["ap_iou50_95"], math.nan),
                     _num(pv["recall_iou50"], math.nan), pv["tp_iou50"], pv["fp_iou50"], pv["fn_iou50"],
                     pv["n_splits"], pv["n_merges"], best["epoch"], rec["train_s"], rec["val_s"])
        else:
            log.info("epoch %d: loss %.4f (%.0f s)", epoch + 1, rec["train_loss"], rec["train_s"])
        epochs_log.append(rec)
        _write_json(out / "training_log.json", {"status": "running", "settings": settings, "gonads": gonad_rec,
                                                 "epochs": epochs_log, "best_epoch": best["epoch"] if best else None,
                                                 "git_sha": git_sha})
        if best is not None and since_best >= int(patience):
            stop_reason = f"early stop: no new best in {since_best} validations"
            break

    # the selected epoch's weights, then the optional post grid on val
    best_model = build_model(model.arch())
    best_model.load_state_dict(read_model_file(out / "checkpoints" / f"epoch_{best['epoch']:03d}.pt")["state_dict"])
    best_model.to(dev).eval()
    grid_rec = None
    post_final = dict(pcfg)
    if grid:
        if dev == "cuda":
            torch.cuda.empty_cache()
        t_g = time.time()
        res = post_grid(best_model, val_g, grid, device=dev, tile_zyx=tile_zyx, min_volume_um3=min_volume_um3)
        ranked = sorted(res, key=lambda r: (_num(r["pooled"]["f1_iou50"]), _num(r["pooled"]["ap_iou50_95"])),
                        reverse=True)
        post_final = dict(ranked[0]["post"])
        grid_rec = {"grid": {k: list(v) for k, v in grid.items()}, "n_configurations": len(res), "results": res,
                    "chosen": post_final, "rule": "pooled F1 at IoU 0.5, then pooled AP (first in grid order on a tie)",
                    "elapsed_s": round(time.time() - t_g, 1)}
        log.info("post grid: %d configurations on val in %.0f s; chosen %s (F1 %.3f, AP %.3f)", len(res),
                 grid_rec["elapsed_s"], post_final, _num(ranked[0]["pooled"]["f1_iou50"], math.nan),
                 _num(ranked[0]["pooled"]["ap_iou50_95"], math.nan))
    sel_val = next(r["val"] for r in epochs_log if r["epoch"] == best["epoch"])
    final_val = sel_val if grid_rec is None else next(r for r in grid_rec["results"] if r["post"] == post_final)
    n_configs = len(entries) + (grid_rec["n_configurations"] if grid_rec else 0)
    selection = {"rule": "pooled protocol F1 at IoU 0.5 over the val gonads, then pooled AP, then the earlier epoch",
                 "epoch": best["epoch"], "post": post_final, "n_epochs_validated": len(entries),
                 "n_configurations_scored_on_val": n_configs,
                 "val_pooled": final_val["pooled"], "val_per_gonad": final_val["per_gonad"]}
    final_meta = {**meta_base, "post": post_final, "epoch": best["epoch"], "git_sha": git_sha,
                  "training": {"train_gonads": [g.tag for g in train_g], "val_gonads": [v.tag for v in val_g],
                               "selected_epoch": best["epoch"], "epochs_run": len(epochs_log),
                               "val_selected": {k: final_val["pooled"].get(k) for k in SUMMARY_KEYS},
                               "created": _utc()}}
    model_sha = save_model_file(out / "model.pt", best_model, final_meta)
    finished = _utc()
    elapsed = time.time() - t_start
    versions = _versions()
    log_doc = {"status": "finished", "name": name or out.name, "started": started, "finished": finished,
               "elapsed_s": round(elapsed, 1), "prep_s": round(prep_s, 1), "settings": settings,
               "gonads": gonad_rec, "epochs": epochs_log, "best_epoch": best["epoch"], "selection": selection,
               "post_grid": grid_rec, "stop_reason": stop_reason, "checkpoints": checkpoints,
               "model_file": str(out / "model.pt"), "model_sha256": model_sha, "git_sha": git_sha,
               "versions": versions, "split_file": str(split_file), "manifest": str(manifest),
               "cache_root": str(cache_root)}
    _write_json(out / "training_log.json", log_doc)
    card = model_card(name or out.name, out, model_sha, settings, gonad_rec, selection, stop_reason,
                      len(epochs_log), git_sha, versions, split_file, manifest, cache_root, finished)
    _write_json(out / "model_card.json", card)
    shutil.copyfile(out / "model_card.json", out / "checkpoints" / "model_card.json")
    log.info("nuclei3d trained in %.0f s: selected epoch %s, post %s, model %s (sha256 %s...)", elapsed,
             best["epoch"], post_final, out / "model.pt", model_sha[:12])
    return {"model_path": str(out / "model.pt"), "model_sha256": model_sha, "card_path": str(out / "model_card.json"),
            "log_path": str(out / "training_log.json"), "best_epoch": best["epoch"], "post": post_final,
            "val_best": final_val["pooled"], "selection": selection, "device": dev, "epochs_run": len(epochs_log),
            "stop_reason": stop_reason, "gonads": gonad_rec}


def model_card(name, out, model_sha, settings, gonad_rec, selection, stop_reason, epochs_run, git_sha, versions,
               split_file, manifest, cache_root, finished) -> dict:
    return {
        "name": name,
        "kind": "nuclei3d 3D nucleus segmentation (germquant segmentation backend 'nuclei3d', EXPERIMENTAL)",
        "what": ("A compact residual 3D U-Net (xy-only first pooling) with three heads (foreground, cell distance "
                 "in um, neighbour closeness) trained on the recovered Imaris nucleus masks of the train gonads, "
                 "with zero weight inside untraced deployed Cellpose labels plus a 0.5 um margin; instances from "
                 "the h-maxima of the distance map and a seeded watershed inside the foreground; the epoch chosen "
                 "on the val gonads by the protocol's F1 at IoU 0.5."),
        "model_file": str(Path(out) / "model.pt"), "model_sha256": model_sha,
        "model_size": (Path(out) / "model.pt").stat().st_size, "date": finished, "git_sha": git_sha,
        "versions": versions, "input_roles": ["dna"], "arch": settings["arch"], "n_parameters": settings["n_parameters"],
        "training_gonads": gonad_rec["train"], "val_gonads": gonad_rec["val"], "excluded_gonads": gonad_rec["excluded"],
        "test_gonads_never_opened": gonad_rec["test_never_opened"],
        "epochs": {"max": settings["epochs"], "run": epochs_run, "selected": selection["epoch"],
                   "stop_reason": stop_reason, "steps_per_epoch": settings["steps_per_epoch"],
                   "batch_size": settings["batch_size"], "patch_zyx": settings["patch"]},
        "post": selection["post"],
        "selection": {k: v for k, v in selection.items() if k not in ("val_pooled", "val_per_gonad")},
        "val_metrics_at_selection": selection["val_pooled"],
        "val_metrics_per_gonad": selection["val_per_gonad"],
        "settings": settings,
        "limitations": [
            "Trained on six N2 gonads of one imaging session (2025-11-05), 582 traced nuclei in all; the truth is a "
            "threshold-curated subset of the nuclei (64 to 146 per gonad) traced with per-gonad Imaris thresholds.",
            "The loss ignores untraced deployed Cellpose labels plus a 0.5 um margin, so a nucleus Cellpose missed is "
            "taught as background, and the deployed labels decide where the truth is silent.",
            "Val holds two hermaphrodites only; the epoch (and any post grid) was chosen on them, so the val numbers "
            "are optimistic.",
            "Development numbers only: no test gonad has been scored (docs/BENCHMARK_PROTOCOL.md).",
        ],
        "split_file": str(split_file), "manifest": str(manifest), "cache_root": str(cache_root),
    }


def parse_grid(text: str | None) -> dict | None:
    """'h_um=0.3,0.4;fg_threshold=0.4,0.5' into {name: [values]} (None for an empty text)."""
    if not text:
        return None
    out = {}
    for part in str(text).split(";"):
        if not part.strip():
            continue
        if "=" not in part:
            raise ValueError(f"grid part {part!r} must read name=v1,v2,...")
        k, v = part.split("=", 1)
        vals = [float(x) for x in v.split(",") if x.strip()]
        if not vals:
            raise ValueError(f"grid part {part!r} has no values")
        out[k.strip()] = vals
    unknown = sorted(set(out) - set(POST_DEFAULTS))
    if unknown:
        raise ValueError(f"unknown instance-step parameters {unknown} (known: {sorted(POST_DEFAULTS)})")
    return out
