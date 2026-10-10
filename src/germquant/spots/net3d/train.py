"""Training of the net3d foci detector (``germquant train-spots``).

Which gonads. The roles come from the split file (``tests/golden/rad51_split.json``), checked against the
pre-registration (``germquant.benchmark.detector.gonad_role``: a split file that moves a registered
gonad is refused). Training uses the train gonads whose spot truth is usable (the manifest's
``spot_truth_usable``, and no Imaris region growing in the split record): for the N2 RAD-51 split that
is n2rad51_HS_herm_000, n2rad51_noHS_herm_000, n2rad51_noHS_male_000 and n2rad51_noHS_male_001; the two
heat-shock males of train (000 and 002) were exported with region growing and are left out. Early
stopping, model selection and the threshold use the val gonads (n2rad51_HS_herm_001 and
n2rad51_noHS_herm_001, both hermaphrodites). No heat-shock male is in train or val. The test gonads are
never opened: every gonad this module loads passes ``guard_test`` first.

What is learned. data.py: Gaussian heatmap targets at the Class A spots, Class B ignored, the loss
masked to the germline and down-weighted outside the traced coverage. The loss is the weighted binary
cross-entropy of the logits against the targets, normalised by the summed weight of the batch.

How. AdamW, a linear warm-up then a cosine decay of the learning rate, bf16 autocast on CUDA (float32
on the CPU), seeds for python, numpy and torch, ``torch.use_deterministic_algorithms(True)`` with
cuDNN benchmarking off and ``CUBLAS_WORKSPACE_CONFIG`` set, so two runs with the same seed on the
same machine give the same weights.

Validation (every ``val_every`` epochs, and the early stopping). Each val gonad goes through the
backend's own path: tiled float32 inference (infer.py), peak finding (peaks.py, one pass at the lowest
threshold of the grid, then filtered; see peaks.py for why that is exact) and the backend's nucleus
assignment, then the bench harness's scoring: expert spots mapped into the cached labels, our counts per
nucleus, Lin's CCC over the traced subset (primary) and over all germline nuclei, and Hungarian spot
matching in the 0.3 um by 0.4 um cylinder (F1, secondary). The scoring reuses
``germquant.benchmark.spots`` with its expensive parts (the traced-subset overlap, the expert spot
labels) computed once per gonad; at the first validation the result is checked against
``germquant.benchmark.detector.score`` itself and the run stops if the two differ. The two val gonads are
pooled (nuclei concatenated for the CCC, counts summed for the F1); for every threshold of the grid
the pooled scores are logged, and the epoch's score is that of its best threshold (CCC first, F1 to
break a tie). The best epoch's weights and threshold become ``model.pt``.
"""
from __future__ import annotations

import json
import logging
import math
import os
import random
import shutil
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from . import data as nd
from . import peaks as pk
from .infer import (DEFAULT_TILE_ZYX, NORM_DEFAULTS, channel_ranges, occupancy, predict_heatmap,
                    read_model_file, save_model_file)

log = logging.getLogger(__name__)

DEFAULT_THRESHOLDS = tuple(round(0.05 + 0.05 * i, 2) for i in range(18))       # 0.05 to 0.90
TRAIN_DEFAULTS = {
    "epochs": 40,
    "patch": (32, 128, 128),
    "seed": 0,
    "batch_size": 8,
    "steps_per_epoch": 200,
    "lr": 1e-3,
    "weight_decay": 1e-4,
    "warmup_steps": 200,
    "min_lr_frac": 0.05,
    "patience": 8,
    "val_every": 1,
    "min_delta": 1e-4,
    "base_channels": 16,
    "levels": 4,
    "use_dna": True,
}
SUMMARY_KEYS = ("ccc_traced", "pearson_traced", "bias_traced", "mean_expert_traced", "mean_ours_traced",
                "n_traced", "ccc_all", "bias_all", "n_germline", "n_matched", "n_ours", "n_expert",
                "precision", "recall", "f1")


# ---------------------------------------------------------------------------------- helpers
def parse_patch(text) -> tuple[int, int, int]:
    """'Z,Y,X' (or a sequence) into three positive ints."""
    if isinstance(text, (list, tuple)):
        vals = [int(v) for v in text]
    else:
        vals = [int(v) for v in str(text).replace("x", ",").split(",") if v.strip()]
    if len(vals) != 3 or any(v <= 0 for v in vals):
        raise ValueError(f"--patch needs three positive integers Z,Y,X, got {text!r}")
    return tuple(vals)


def _utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _finite(x):
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _num(x, default: float = -math.inf) -> float:
    """A float, ``default`` for NaN or a missing value (0.0 stays 0.0)."""
    f = _finite(x)
    return default if f is None else f


def _versions() -> dict:
    from importlib import metadata

    out = {}
    for name in ("germquant", "torch", "numpy", "scipy", "pandas"):
        try:
            out[name] = metadata.version(name)
        except Exception:  # noqa: BLE001
            out[name] = "not-installed"
    return out


def gpu_free_gb() -> float:
    """Free memory of GPU 0 in GB. nvidia-smi first: on Windows (WDDM) ``torch.cuda.mem_get_info``
    reports the whole card as free while other processes hold most of it; the torch figure is the
    fallback when nvidia-smi is not on the path."""
    import subprocess

    try:
        out = subprocess.run(["nvidia-smi", "--query-gpu=memory.free", "--format=csv,noheader,nounits"],
                             capture_output=True, text=True, timeout=30, check=True).stdout
        return float(out.strip().splitlines()[0]) * 1024 ** 2 / 1e9          # MiB -> GB
    except Exception:  # noqa: BLE001 - no nvidia-smi: ask torch
        import torch

        return torch.cuda.mem_get_info()[0] / 1e9


def wait_for_gpu(min_free_gb: float, *, poll_s: float = 60.0, max_minutes: float = 120.0) -> float:
    """Block until GPU 0 has at least ``min_free_gb`` GB free (``gpu_free_gb``, polling every ``poll_s``
    seconds, at most ``max_minutes``); returns the free GB. RuntimeError on timeout."""
    t0 = time.time()
    while True:
        free_gb = gpu_free_gb()
        if free_gb >= min_free_gb:
            return free_gb
        if time.time() - t0 > max_minutes * 60:
            raise RuntimeError(f"the GPU never had {min_free_gb} GB free within {max_minutes} min "
                               f"({free_gb:.1f} GB free at the end)")
        log.info("waiting for %.1f GB of free GPU memory (%.1f GB free now)", min_free_gb, free_gb)
        time.sleep(poll_s)


def set_determinism(seed: int) -> None:
    import torch

    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")      # before the first cuBLAS call
    random.seed(seed)
    np.random.seed(seed % (2 ** 32))
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True)


# ------------------------------------------------------------------------- gonads and roles
def resolve_gonads(split_file, manifest) -> dict:
    """The train gonads with usable spot truth, the val gonads, what was left out and why, and the
    test gonads (named, never opened)."""
    from germquant.benchmark import detector as bd
    from germquant.benchmark.manifest import load_manifest

    split = json.loads(Path(split_file).read_text(encoding="utf-8"))
    roles_doc = split.get("roles")
    if not isinstance(roles_doc, dict):
        raise ValueError(f"{split_file}: the split file needs a 'roles' mapping (train, val, test)")
    entries = {e.tag: e for e in load_manifest(manifest)}
    rg = {g.get("tag"): bool((g.get("ims") or {}).get("region_growing")) for g in split.get("gonads", [])
          if isinstance(g, dict)}
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
        if not bd.spot_truth_usable(e) or rg.get(t):
            out["excluded"].append({"tag": t, "role": role,
                                    "reason": "spot truth unusable (Imaris region growing)"})
            continue
        if role in ("train", "val"):
            out[role].append(t)
    bd.guard_test(out["train"] + out["val"], allow_test=False, roles=roles)
    return out


def _gonad_facts(entry, cache_meta: dict) -> dict:
    return {"condition": str(entry.extra.get("condition") or cache_meta.get("condition") or ""),
            "sex": str(entry.extra.get("sex") or cache_meta.get("sex") or "")}


def load_train_gonad(tag: str, entry, cache_root, prep_dir, *, use_dna: bool, data_cfg: dict,
                     norm: dict) -> nd.TrainGonad:
    from germquant.benchmark import detector as bd
    from germquant.benchmark.manifest import load_truth

    bd.guard_test([tag], allow_test=False, roles={tag: bd.entry_role(entry)})
    cache = bd.load_cache(tag, cache_root, arrays=("foci", "dna", "labels"), mmap=True)
    truth = load_truth(entry.truth_dir, entry.truth_stem)
    if truth is None:
        raise FileNotFoundError(f"{tag}: no truth files for stem {entry.truth_stem!r} in {entry.truth_dir}")
    channels = [cache.foci] + ([cache.dna] if use_dna else [])
    if use_dna and cache.dna is None:
        raise FileNotFoundError(f"{tag}: the bench cache holds no DNA channel")
    t0 = time.time()
    region = nd.load_region(prep_dir, tag, cache.meta, cache.labels, truth.labels_path, cache.spacing, data_cfg)
    pos = nd.truth_positions(truth.spots, truth.params)
    is_b = nd.class_b_mask(truth.spots["spot_class"]) if truth.has_spot_class else np.zeros(len(pos), bool)
    sp = np.asarray(cache.spacing, float)
    nuc = cache.nuclei
    germ = nuc[nuc["in_germline"]] if "in_germline" in nuc.columns else nuc
    centers = (germ[["centroid_z_um", "centroid_y_um", "centroid_x_um"]].to_numpy(float) / sp
               if len(germ) else np.zeros((0, 3)))
    g = nd.TrainGonad(tag=tag, spacing=cache.spacing, channels=channels,
                      ranges=channel_ranges(channels, norm), region=region, spots_a=pos[~is_b],
                      spots_b=pos[is_b], centers=centers, **_gonad_facts(entry, cache.meta))
    counts = np.bincount(np.asarray(region[:, ::4, ::4]).ravel(), minlength=3)
    g.info = {"n_spots": int(len(pos)), "n_class_a": int((~is_b).sum()), "n_class_b": int(is_b.sum()),
              "n_class_a_in_region": int(len(g.spot_pool)), "n_germline_nuclei": int(len(centers)),
              "region_fraction_sampled": {"outside": float(counts[0] / counts.sum()),
                                          "germline_untraced": float(counts[1] / counts.sum()),
                                          "germline_traced": float(counts[2] / counts.sum())},
              "ranges": [list(r) for r in g.ranges], "prep_s": round(time.time() - t0, 1),
              "labels_sha256": ((cache.meta.get("arrays") or {}).get("labels") or {}).get("data_sha256"),
              "foci_sha256": ((cache.meta.get("arrays") or {}).get("foci") or {}).get("data_sha256")}
    return g


# ------------------------------------------------------------------------------- validation
class ValScorer:
    """The bench harness's per-gonad scores (germquant.benchmark.spots.benchmark_gonad) for many spot
    sets on one gonad: the expert spots' labels, the germline per-nucleus table and the traced flags
    are computed once, the per-spot set changes. Only the numbers the early stopping reads are
    produced; ``check_against_harness`` compares them with the harness itself."""

    def __init__(self, tag: str, cache, truth, truth_labels, *, tol_xy_um: float, tol_z_um: float):
        from germquant.benchmark import spots as bs

        self.tag = tag
        self.spacing = cache.spacing
        self.tol = (float(tol_xy_um), float(tol_z_um))
        labels = cache.labels
        tz, ty, tx = bs._spot_voxels(truth.spots, cache.spacing)
        spot_labels = bs.label_at_points(labels, tz, ty, tx)
        nuclei = cache.nuclei.copy()
        nuclei["n_spots"] = 0
        per_nuc, rule = bs.per_nucleus_table(nuclei, spot_labels, np.zeros(len(nuclei), int),
                                             truth_nuclei=truth.nuclei, labels=labels, truth_labels=truth_labels)
        self.rule = rule
        self.nucleus_ids = per_nuc["nucleus_id"].to_numpy(np.int64)
        self.expert_n = per_nuc["expert_n"].to_numpy(float)
        self.traced = per_nuc["traced"].to_numpy(bool)
        self.truth_xyz = (truth.spots[["x_um", "y_um", "z_um"]].to_numpy(float) if len(truth.spots)
                          else np.zeros((0, 3)))

    def score(self, nucleus_id: np.ndarray, xyz_um: np.ndarray) -> dict:
        from germquant.benchmark import spots as bs
        from germquant.benchmark.stats import agreement_metrics

        nid = np.asarray(nucleus_id, np.int64)
        counts = pd.Series(nid).value_counts() if len(nid) else pd.Series(dtype=np.int64)
        ours = pd.Series(self.nucleus_ids).map(counts).fillna(0).to_numpy(float)
        a = bs.match_spots(np.asarray(xyz_um, float).reshape(-1, 3), self.truth_xyz,
                           tol_xy_um=self.tol[0], tol_z_um=self.tol[1])
        m_all = agreement_metrics(self.expert_n, ours)
        m_tr = agreement_metrics(self.expert_n[self.traced], ours[self.traced])
        n_ours, n_exp = int(len(nid)), int(len(self.truth_xyz))
        prf = _prf(a.n_matched, n_ours, n_exp)
        return {"ccc_traced": m_tr["ccc"], "pearson_traced": m_tr["pearson_r"], "bias_traced": m_tr["bias"],
                "mean_expert_traced": m_tr["mean_expert"], "mean_ours_traced": m_tr["mean_ours"],
                "n_traced": int(self.traced.sum()), "ccc_all": m_all["ccc"], "bias_all": m_all["bias"],
                "n_germline": int(len(self.nucleus_ids)), "n_matched": int(a.n_matched), "n_ours": n_ours,
                "n_expert": n_exp, **prf, "_expert": self.expert_n, "_ours": ours, "_traced": self.traced}


def _prf(n_matched: int, n_ours: int, n_expert: int) -> dict:
    from germquant.benchmark.stats import precision_recall_f1

    s = precision_recall_f1(n_matched, n_ours, n_expert)
    return {"precision": s["precision"], "recall": s["recall"], "f1": s["f1"]}


def pool_scores(rows: dict) -> dict:
    """Pooled val scores: the per-nucleus agreement over the concatenated nuclei, F1 from summed counts."""
    from germquant.benchmark.stats import agreement_metrics

    e = np.concatenate([r["_expert"] for r in rows.values()])
    o = np.concatenate([r["_ours"] for r in rows.values()])
    tr = np.concatenate([r["_traced"] for r in rows.values()])
    m_all = agreement_metrics(e, o)
    m_tr = agreement_metrics(e[tr], o[tr])
    nm = sum(r["n_matched"] for r in rows.values())
    no = sum(r["n_ours"] for r in rows.values())
    ne = sum(r["n_expert"] for r in rows.values())
    return {"ccc_traced": m_tr["ccc"], "pearson_traced": m_tr["pearson_r"], "bias_traced": m_tr["bias"],
            "mean_expert_traced": m_tr["mean_expert"], "mean_ours_traced": m_tr["mean_ours"],
            "n_traced": int(tr.sum()), "ccc_all": m_all["ccc"], "bias_all": m_all["bias"],
            "n_germline": int(len(e)), "n_matched": int(nm), "n_ours": int(no), "n_expert": int(ne),
            **_prf(nm, no, ne)}


@dataclass
class ValGonad:
    tag: str
    cache: object
    channels: list
    ranges: list
    occ: tuple
    scorer: ValScorer
    marker: str
    condition: str = ""
    sex: str = ""


def load_val_gonad(tag: str, entry, cache_root, *, use_dna: bool, norm: dict, tol_xy_um: float,
                   tol_z_um: float) -> ValGonad:
    from germquant.benchmark import detector as bd
    from germquant.benchmark.manifest import load_labels, load_truth

    bd.guard_test([tag], allow_test=False, roles={tag: bd.entry_role(entry)})
    cache = bd.load_cache(tag, cache_root, arrays=("foci", "dna", "labels"), mmap=True)
    truth = load_truth(entry.truth_dir, entry.truth_stem)
    if truth is None:
        raise FileNotFoundError(f"{tag}: no truth files for stem {entry.truth_stem!r} in {entry.truth_dir}")
    truth_labels = load_labels(truth.labels_path) if truth.has_labels else None
    scorer = ValScorer(tag, cache, truth, truth_labels, tol_xy_um=tol_xy_um, tol_z_um=tol_z_um)
    del truth_labels
    channels = [cache.foci] + ([cache.dna] if use_dna else [])
    return ValGonad(tag=tag, cache=cache, channels=channels, ranges=channel_ranges(channels, norm),
                    occ=occupancy(cache.labels), scorer=scorer, marker=cache.marker,
                    **_gonad_facts(entry, cache.meta))


def _xyz_um(found: pd.DataFrame, spacing) -> np.ndarray:
    dz, dy, dx = (float(s) for s in spacing)
    return np.stack([found["xr"].to_numpy(float) * dx, found["yr"].to_numpy(float) * dy,
                     found["zr"].to_numpy(float) * dz], axis=1) if len(found) else np.zeros((0, 3))


def evaluate(model, vals: list[ValGonad], *, thresholds=DEFAULT_THRESHOLDS, device: str = "cpu",
             tile_zyx=DEFAULT_TILE_ZYX, separation_um=pk.DEFAULT_SEPARATION_UM,
             assign_max_um: float = pk.DEFAULT_ASSIGN_MAX_UM) -> dict:
    """Score the model on the val gonads at every threshold of the grid (module docstring)."""
    thresholds = sorted(float(t) for t in thresholds)
    found, timing = {}, {}
    for vg in vals:
        t0 = time.time()
        heat = predict_heatmap(model, vg.channels, vg.ranges, device=device, tile_zyx=tile_zyx, occ=vg.occ)
        t1 = time.time()
        f = pk.find_spots(heat, vg.cache.labels, vg.cache.spacing, thresholds[0], separation_um=separation_um,
                          assign_max_um=assign_max_um)
        del heat
        found[vg.tag] = f[f["nucleus_id"] > 0].reset_index(drop=True)
        timing[vg.tag] = {"infer_s": round(t1 - t0, 1), "peaks_s": round(time.time() - t1, 1),
                          "n_peaks_at_floor": int(len(f))}
    sweep = []
    for t in thresholds:
        rows = {}
        for vg in vals:
            f = found[vg.tag]
            f = f[f["prob"] >= t]
            rows[vg.tag] = vg.scorer.score(f["nucleus_id"].to_numpy(), _xyz_um(f, vg.cache.spacing))
        pooled = pool_scores(rows)
        sweep.append({"threshold": t, **pooled,
                      "per_gonad": {k: {s: r[s] for s in SUMMARY_KEYS} for k, r in rows.items()}})

    def key(s):
        c = _finite(s["ccc_traced"])
        f1 = _finite(s["f1"])
        return (-math.inf if c is None else round(c, 6), -math.inf if f1 is None else f1)

    best = max(sweep, key=key)
    return {"best_threshold": best["threshold"], "best": best, "sweep": sweep, "found": found, "timing": timing}


def check_against_harness(vg: ValGonad, found: pd.DataFrame, threshold: float, manifest) -> dict:
    """Score one spot set with germquant.benchmark.detector.score itself and with the ValScorer; the
    differences of the numbers the early stopping reads (all zero expected)."""
    from germquant.benchmark import detector as bd
    from germquant.spots.backends import PER_NUC_COLS, PER_SPOT_COLS

    f = found[found["prob"] >= threshold].reset_index(drop=True)
    xyz = _xyz_um(f, vg.cache.spacing)
    per_spot = pd.DataFrame({"spot_id": np.arange(len(f)), "nucleus_id": f["nucleus_id"].astype(int).to_numpy(),
                             "marker": vg.marker, "z_um": xyz[:, 2], "y_um": xyz[:, 1], "x_um": xyz[:, 0],
                             "effect_size": np.nan, "spot_mean_intensity": np.nan, "detector": "net3d"},
                            columns=PER_SPOT_COLS)
    ids = vg.cache.nuclei["nucleus_id"].astype(int).to_numpy()
    counts = per_spot["nucleus_id"].value_counts()
    per_nuc = pd.DataFrame({"nucleus_id": ids, "marker": vg.marker,
                            "n_spots": [int(counts.get(i, 0)) for i in ids], "mean_spot_intensity": np.nan,
                            "detector": "net3d"}, columns=PER_NUC_COLS)
    row, _pn, _m = bd.score(vg.tag, per_spot, per_nuc, manifest=manifest, cache=vg.cache,
                            tol_xy_um=vg.scorer.tol[0], tol_z_um=vg.scorer.tol[1])
    mine = vg.scorer.score(f["nucleus_id"].to_numpy(), xyz)
    pairs = {"ccc_traced": "ccc_traced", "ccc_all": "ccc_all", "bias_traced": "bias_traced", "f1": "f1",
             "n_matched": "n_matched", "n_traced": "n_traced"}
    diffs = {}
    for k_mine, k_h in pairs.items():
        a, b = _finite(mine[k_mine]), _finite(row[k_h])
        diffs[k_mine] = 0.0 if (a is None and b is None) else (math.inf if (a is None or b is None) else abs(a - b))
    return {"tag": vg.tag, "threshold": float(threshold), "max_abs_diff": max(diffs.values()) if diffs else 0.0,
            "diffs": diffs, "harness": {k: _finite(row[k]) for k in pairs.values()}}


# ------------------------------------------------------------------------------------ train
def _lr_at(step: int, total: int, base: float, warmup: int, min_frac: float) -> float:
    if step < warmup:
        return base * (step + 1) / max(1, warmup)
    p = (step - warmup) / max(1, total - warmup)
    return base * (min_frac + (1 - min_frac) * 0.5 * (1 + math.cos(math.pi * min(1.0, p))))


def _strip(ev: dict) -> dict:
    """An evaluation without the per-spot tables and the per-nucleus arrays (for the JSON log)."""
    def clean(d):
        return {k: (_finite(v) if isinstance(v, (float, np.floating)) else v) for k, v in d.items()
                if not k.startswith("_")}

    return {"best_threshold": ev["best_threshold"], "best": clean({k: v for k, v in ev["best"].items()
                                                                 if k != "per_gonad"}),
            "best_per_gonad": {t: clean(r) for t, r in ev["best"]["per_gonad"].items()},
            "sweep": [{"threshold": s["threshold"], **clean({k: v for k, v in s.items() if k not in
                                                             ("threshold", "per_gonad")}),
                       "per_gonad": {t: clean(r) for t, r in s["per_gonad"].items()}} for s in ev["sweep"]],
            "timing": ev["timing"]}


def train_spots(cache_root, split_file, out_dir, **kwargs) -> dict:
    """Train net3d on the bench cache (module docstring) into ``out_dir``: ``model.pt`` (the best
    epoch's weights with the threshold chosen on val), ``checkpoints/best.pt`` and ``last.pt``,
    ``training_log.json`` and ``model_card.json``. Returns the paths and the best validation. The
    keyword arguments are those of ``_train``. torch's deterministic-algorithm and cuDNN settings are
    restored afterwards, so a caller in the same process gets its own back."""
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
           batch_size: int = TRAIN_DEFAULTS["batch_size"],
           steps_per_epoch: int = TRAIN_DEFAULTS["steps_per_epoch"], lr: float = TRAIN_DEFAULTS["lr"],
           weight_decay: float = TRAIN_DEFAULTS["weight_decay"],
           warmup_steps: int = TRAIN_DEFAULTS["warmup_steps"], patience: int = TRAIN_DEFAULTS["patience"],
           val_every: int = TRAIN_DEFAULTS["val_every"], use_dna: bool = True, device: str = "auto",
           base_channels: int = TRAIN_DEFAULTS["base_channels"], levels: int = TRAIN_DEFAULTS["levels"],
           thresholds=DEFAULT_THRESHOLDS, tile_zyx=DEFAULT_TILE_ZYX, data_cfg: dict | None = None,
           norm: dict | None = None, separation_um=pk.DEFAULT_SEPARATION_UM,
           assign_max_um: float = pk.DEFAULT_ASSIGN_MAX_UM, prep_dir=None, min_free_gb: float = 0.0,
           wait_minutes: float = 120.0, name: str | None = None, check_harness: bool = True,
           tol_xy_um: float = 0.3, tol_z_um: float = 0.4) -> dict:
    import torch
    import torch.nn.functional as F

    from germquant import provenance
    from germquant.benchmark import detector as bd

    from .model import build_model, n_parameters, size_divisor

    t_start = time.time()
    started = _utc()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "checkpoints").mkdir(exist_ok=True)
    manifest = Path(manifest) if manifest else bd.DEFAULT_MANIFEST
    patch = parse_patch(patch)
    div = size_divisor(levels)
    if any(p % d for p, d in zip(patch, div)):
        raise ValueError(f"--patch {patch} must be a multiple of {div} (the network's size divisor)")
    dcfg = {**nd.DATA_DEFAULTS, **(data_cfg or {})}
    ncfg = {**NORM_DEFAULTS, **(norm or {})}
    thresholds = sorted(float(t) for t in thresholds)

    # the device and determinism come first: CUBLAS_WORKSPACE_CONFIG must be set before CUDA starts
    want = str(device or "auto").lower()
    dev = "cuda" if (want in ("auto", "cuda") and torch.cuda.is_available()) else "cpu"
    if want == "cuda" and dev != "cuda":
        log.warning("CUDA requested but not available; training on the CPU")
    set_determinism(int(seed))
    if dev == "cuda" and min_free_gb > 0:
        free = wait_for_gpu(float(min_free_gb), max_minutes=float(wait_minutes))
        log.info("GPU: %.1f GB free", free)

    gon = resolve_gonads(split_file, manifest)
    if not gon["train"]:
        raise ValueError("no train gonad with usable spot truth")
    if not gon["val"]:
        raise ValueError("no val gonad: early stopping and the threshold need one")
    prep = Path(prep_dir) if prep_dir else out / "prep"
    log.info("train gonads %s; val %s; left out %s; test (never opened) %s", gon["train"], gon["val"],
             [e["tag"] for e in gon["excluded"]], gon["test"])
    t0 = time.time()
    train_g = [load_train_gonad(t, gon["entries"][t], cache_root, prep, use_dna=use_dna, data_cfg=dcfg, norm=ncfg)
               for t in gon["train"]]
    val_g = [load_val_gonad(t, gon["entries"][t], cache_root, use_dna=use_dna, norm=ncfg, tol_xy_um=tol_xy_um,
                            tol_z_um=tol_z_um) for t in gon["val"]]
    prep_s = time.time() - t0
    log.info("data ready in %.0f s", prep_s)

    in_ch = 2 if use_dna else 1
    model = build_model({"in_channels": in_ch, "base_channels": base_channels, "levels": levels}).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    total_steps = int(epochs) * int(steps_per_epoch)
    roles_in = ["foci", "dna"] if use_dna else ["foci"]
    git_sha = provenance.git_sha()
    settings = {"epochs": int(epochs), "patch": list(patch), "seed": int(seed), "batch_size": int(batch_size),
                "steps_per_epoch": int(steps_per_epoch), "lr": float(lr), "weight_decay": float(weight_decay),
                "warmup_steps": int(warmup_steps), "min_lr_frac": TRAIN_DEFAULTS["min_lr_frac"],
                "patience": int(patience), "val_every": int(val_every), "min_delta": TRAIN_DEFAULTS["min_delta"],
                "device": dev, "amp": "bf16 autocast" if dev == "cuda" else "float32",
                "deterministic": True, "thresholds": thresholds, "tile_zyx": list(tile_zyx),
                "data": {k: (list(v) if isinstance(v, tuple) else v) for k, v in dcfg.items()}, "norm": ncfg,
                "separation_um": list(separation_um), "assign_max_um": float(assign_max_um),
                "tol_xy_um": float(tol_xy_um), "tol_z_um": float(tol_z_um), "input_roles": roles_in,
                "arch": model.arch(), "n_parameters": n_parameters(model)}
    gonad_rec = {"train": [{"tag": g.tag, "condition": g.condition, "sex": g.sex, **g.info} for g in train_g],
                 "val": [{"tag": v.tag, "condition": v.condition, "sex": v.sex,
                          "n_germline_nuclei": int(len(v.scorer.nucleus_ids)), "n_traced": int(v.scorer.traced.sum()),
                          "n_expert_spots": int(len(v.scorer.truth_xyz)), "traced_rule": v.scorer.rule,
                          "ranges": [list(r) for r in v.ranges]} for v in val_g],
                 "excluded": gon["excluded"], "test_never_opened": gon["test"]}
    meta_base = {"input_roles": roles_in, "norm": ncfg,
                 "peaks": {"threshold": 0.5, "separation_um": list(separation_um), "assign_max_um": float(assign_max_um)},
                 "target": {"sigma_um": list(dcfg["sigma_um"]), "classes": ["Class A"]}}

    epochs_log, checkpoints = [], []
    best = {"epoch": None, "score": (-math.inf, -math.inf), "threshold": None, "eval": None, "sha256": None}
    since_best = 0
    stop_reason = "max epochs reached"
    harness_check = None
    step = 0
    log.info("training on %s: %d parameters, %d epochs x %d steps x %d patches %s", dev, n_parameters(model),
             epochs, steps_per_epoch, batch_size, patch)
    for epoch in range(int(epochs)):
        model.train()
        t_ep = time.time()
        losses = []
        if dev == "cuda":
            torch.cuda.reset_peak_memory_stats()
        for xb, yb, wb in nd.Prefetch(nd.batches(train_g, patch, dcfg, int(seed), epoch, steps_per_epoch,
                                                 batch_size), depth=4):
            for gparam in opt.param_groups:
                gparam["lr"] = _lr_at(step, total_steps, lr, warmup_steps, TRAIN_DEFAULTS["min_lr_frac"])
            x = torch.from_numpy(xb).to(dev, non_blocking=True)
            y = torch.from_numpy(yb).to(dev, non_blocking=True)
            w = torch.from_numpy(wb).to(dev, non_blocking=True)
            with torch.autocast(device_type="cuda", dtype=torch.bfloat16, enabled=(dev == "cuda")):
                logits = model(x)
            bce = F.binary_cross_entropy_with_logits(logits.float(), y, reduction="none")
            loss = (bce * w).sum() / w.sum().clamp_min(1.0)
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()
            losses.append(float(loss.detach()))
            step += 1
        rec = {"epoch": epoch + 1, "lr_end": opt.param_groups[0]["lr"], "train_loss": float(np.mean(losses)),
               "train_s": round(time.time() - t_ep, 1)}
        if dev == "cuda":
            rec["gpu_peak_gb"] = round(torch.cuda.max_memory_allocated() / 1e9, 2)
            rec["gpu_reserved_gb"] = round(torch.cuda.max_memory_reserved() / 1e9, 2)
        last_path = out / "checkpoints" / "last.pt"
        sha_last = save_model_file(last_path, model, {**meta_base, "epoch": epoch + 1, "git_sha": git_sha})
        rec["last_sha256"] = sha_last
        if (epoch + 1) % max(1, int(val_every)) == 0 or epoch + 1 == int(epochs):
            t_v = time.time()
            if dev == "cuda":
                # hand the training batches' cached blocks back before the inference tiles (and the tiles'
                # blocks back after): the card is shared, and without this the caching allocator holds the
                # sum of both working sets (about 16 GB measured) instead of the larger one
                torch.cuda.empty_cache()
            ev = evaluate(model, val_g, thresholds=thresholds, device=dev, tile_zyx=tile_zyx,
                          separation_um=separation_um, assign_max_um=assign_max_um)
            if dev == "cuda":
                rec["gpu_val_reserved_gb"] = round(torch.cuda.max_memory_reserved() / 1e9, 2)
                torch.cuda.empty_cache()
                torch.cuda.reset_peak_memory_stats()
            if check_harness and harness_check is None:
                harness_check = [check_against_harness(v, ev["found"][v.tag], ev["best_threshold"], manifest)
                                 for v in val_g]
                worst = max(c["max_abs_diff"] for c in harness_check)
                if not worst <= 1e-9:
                    raise RuntimeError(f"the validation scorer disagrees with the bench harness by {worst}: "
                                       f"{harness_check}")
                log.info("validation scorer agrees with germquant.benchmark.detector.score (max diff %.2e)", worst)
            rec["val"] = _strip(ev)
            rec["val_s"] = round(time.time() - t_v, 1)
            b = ev["best"]
            score = (round(_num(b["ccc_traced"]), 6), _num(b["f1"]))
            improved = (score[0] > best["score"][0] + TRAIN_DEFAULTS["min_delta"]) or (
                abs(score[0] - best["score"][0]) <= TRAIN_DEFAULTS["min_delta"] and score[1] > best["score"][1])
            if best["epoch"] is None or improved:
                sha_best = save_model_file(out / "checkpoints" / "best.pt", model,
                                           {**meta_base, "epoch": epoch + 1, "git_sha": git_sha,
                                            "peaks": {**meta_base["peaks"], "threshold": ev["best_threshold"]}})
                best = {"epoch": epoch + 1, "score": score, "threshold": ev["best_threshold"], "eval": rec["val"],
                        "sha256": sha_best}
                since_best = 0
                rec["best_sha256"] = sha_best
                checkpoints.append({"epoch": epoch + 1, "file": "checkpoints/best.pt", "sha256": sha_best})
            else:
                since_best += 1
            log.info("epoch %d: loss %.4f | val traced CCC %.3f F1 %.3f at threshold %.2f (germline CCC %.3f) | "
                     "best epoch %s | %.0f s train, %.0f s val", epoch + 1, rec["train_loss"],
                     _num(b["ccc_traced"], math.nan), _num(b["f1"], math.nan),
                     ev["best_threshold"], _num(b["ccc_all"], math.nan), best["epoch"], rec["train_s"],
                     rec["val_s"])
        else:
            log.info("epoch %d: loss %.4f (%.0f s)", epoch + 1, rec["train_loss"], rec["train_s"])
        epochs_log.append(rec)
        _write_json(out / "training_log.json", {"status": "running", "settings": settings, "gonads": gonad_rec,
                                                 "epochs": epochs_log, "best_epoch": best["epoch"],
                                                 "harness_check": harness_check, "git_sha": git_sha})
        if best["epoch"] is not None and since_best >= int(patience):
            stop_reason = f"early stop: no val improvement in {since_best} evaluations"
            break

    # the final model: the best epoch's weights and threshold, with its record
    final_meta = {**meta_base, "peaks": {**meta_base["peaks"], "threshold": best["threshold"]},
                  "epoch": best["epoch"], "git_sha": git_sha,
                  "training": {"train_gonads": [g.tag for g in train_g], "val_gonads": [v.tag for v in val_g],
                               "truth_classes": "Class A spots as targets; Class B ignored (zero weight in a ball)",
                               "best_epoch": best["epoch"], "epochs_run": len(epochs_log),
                               "val_best": best["eval"]["best"] if best["eval"] else None,
                               "created": _utc()}}
    best_model = build_model(model.arch())
    best_model.load_state_dict(read_model_file(out / "checkpoints" / "best.pt")["state_dict"])
    model_sha = save_model_file(out / "model.pt", best_model, final_meta)
    finished = _utc()
    elapsed = time.time() - t_start
    versions = _versions()
    log_doc = {"status": "finished", "name": name or out.name, "started": started, "finished": finished,
               "elapsed_s": round(elapsed, 1), "prep_s": round(prep_s, 1), "settings": settings,
               "gonads": gonad_rec, "epochs": epochs_log, "best_epoch": best["epoch"],
               "best_threshold": best["threshold"], "stop_reason": stop_reason, "harness_check": harness_check,
               "checkpoints": checkpoints, "last_sha256": epochs_log[-1].get("last_sha256") if epochs_log else None,
               "model_file": str(out / "model.pt"), "model_sha256": model_sha, "git_sha": git_sha,
               "versions": versions, "split_file": str(split_file), "manifest": str(manifest),
               "cache_root": str(cache_root)}
    _write_json(out / "training_log.json", log_doc)
    card = model_card(name or out.name, out, model_sha, settings, gonad_rec, best, stop_reason, len(epochs_log),
                      git_sha, versions, split_file, manifest, cache_root, finished)
    _write_json(out / "model_card.json", card)
    shutil.copyfile(out / "model_card.json", out / "checkpoints" / "model_card.json")
    log.info("net3d trained in %.0f s: best epoch %s, threshold %s, model %s (sha256 %s...)", elapsed,
             best["epoch"], best["threshold"], out / "model.pt", model_sha[:12])
    return {"model_path": str(out / "model.pt"), "model_sha256": model_sha, "card_path": str(out / "model_card.json"),
            "log_path": str(out / "training_log.json"), "best_epoch": best["epoch"],
            "best_threshold": best["threshold"], "val_best": best["eval"]["best"] if best["eval"] else None,
            "device": dev, "epochs_run": len(epochs_log), "stop_reason": stop_reason, "gonads": gonad_rec}


def model_card(name, out, model_sha, settings, gonad_rec, best, stop_reason, epochs_run, git_sha, versions,
               split_file, manifest, cache_root, finished) -> dict:
    vb = best["eval"]["best"] if best["eval"] else {}
    return {
        "name": name,
        "kind": "net3d RAD-51 foci heatmap detector (germquant spots backend 'net3d', EXPERIMENTAL)",
        "what": ("A compact residual 3D U-Net (xy-only first pooling) trained on the lab's Imaris spots: Gaussian "
                 "heatmap targets at the Class A spots of the train gonads with usable spot truth, Class B "
                 "ignored, loss restricted to the germline (dilated 0.3 um) and down-weighted outside the "
                 "traced truth coverage; epoch and peak threshold chosen on the val gonads with the bench "
                 "harness (per-nucleus CCC over the traced subset first, spot F1 second)."),
        "model_file": str(Path(out) / "model.pt"), "model_sha256": model_sha,
        "model_size": (Path(out) / "model.pt").stat().st_size,
        "date": finished, "git_sha": git_sha, "versions": versions,
        "input_roles": settings["input_roles"], "arch": settings["arch"], "n_parameters": settings["n_parameters"],
        "training_gonads": gonad_rec["train"], "val_gonads": gonad_rec["val"], "excluded_gonads": gonad_rec["excluded"],
        "test_gonads_never_opened": gonad_rec["test_never_opened"],
        "truth_classes_used": {"targets": ["Class A"], "ignored": ["Class B"],
                               "note": "Class B (Imaris's isolated-spot class) gets zero weight inside a "
                                       "0.5 um (z) by 0.3 um (xy) ellipsoid; it is neither taught nor penalised."},
        "epochs": {"max": settings["epochs"], "run": epochs_run, "best": best["epoch"], "stop_reason": stop_reason,
                   "steps_per_epoch": settings["steps_per_epoch"], "batch_size": settings["batch_size"],
                   "patch_zyx": settings["patch"]},
        "threshold": best["threshold"],
        "val_metrics_at_best_epoch": {k: vb.get(k) for k in SUMMARY_KEYS} if vb else None,
        "val_metrics_per_gonad": best["eval"]["best_per_gonad"] if best["eval"] else None,
        "settings": settings,
        "limitations": [
            "Trained on four N2 gonads of one imaging session (2025-11-05): three without heat shock "
            "(one hermaphrodite, two males) and one heat-shocked hermaphrodite.",
            "No heat-shock male is in train or val: both train heat-shock males were exported with Imaris "
            "region growing and are excluded, so the heat-shock undercount of the deployed recipe is not "
            "represented in training or model selection.",
            "Val holds two hermaphrodites only; the epoch and threshold were selected without a male.",
            "The truth is threshold-curated automatic Imaris Spots, not hand-placed points.",
            "Development numbers only: no test gonad has been scored (docs/BENCHMARK_PROTOCOL.md).",
        ],
        "split_file": str(split_file), "manifest": str(manifest), "cache_root": str(cache_root),
    }


def _write_json(path: Path, doc: dict) -> None:
    from germquant.benchmark.report import jsonable

    tmp = path.with_name(path.name + ".partial")
    tmp.write_text(json.dumps(jsonable(doc), indent=1) + "\n", encoding="utf-8")
    os.replace(tmp, path)
