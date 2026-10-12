"""Which base weights a fine-tuned Cellpose model started from (docs/DESIGN_system_packs.md, 3.11 and 15d).

germline_nuclei_combined was trained by scripts/train_nucleus_model_combined.py with CellposeModel(gpu=True),
whose default weights depend on the Cellpose version installed at training time (cpsam_v2 from 4.2.0 on).
File dates point to cpsam_v2 but do not prove it. This module compares the fine-tuned state dict with each
candidate base tensor by tensor:

* for every tensor name the model and a base share with the same shape: whether it is bit-identical, and its
  relative L2 distance ``||model - base|| / ||base||``;
* per base: the share of the model's tensors it covers, the number of identical tensors, and the overall
  relative distance ``sqrt(sum ||model - base||^2) / sqrt(sum ||base||^2)`` over the shared tensors.

Fine-tuning moves the weights a little away from their start, so the true base is much closer than any other.
The verdict names the closest base only when it covers at least ``min_coverage`` of the model's tensors and the
next base is at least ``separation`` times farther away; otherwise it says "not determined", and the card
must say so. torch is imported only to read the files (scripts/compare_base_weights.py runs it on the
workstation; the tests use small fakes).
"""
from __future__ import annotations

import math
from collections.abc import Mapping
from pathlib import Path

import numpy as np

METHOD = ("tensor-by-tensor comparison of the fine-tuned state dict with each candidate base state dict "
          "(germquant.packs.basecheck): bit-identical tensors and the relative L2 distance ||model - base|| / "
          "||base|| over the tensors the two share with the same shape; the closest base is named only when it "
          "covers at least {cov:.0%} of the model's tensors and the next base is at least {sep:g} times farther")


def _as_array(t) -> np.ndarray:
    if isinstance(t, np.ndarray):
        return t
    detach = getattr(t, "detach", None)
    if detach is not None:                         # a torch tensor (bfloat16 has no numpy dtype: go via float32)
        t = detach().cpu()
        if str(getattr(t, "dtype", "")) == "torch.bfloat16":
            t = t.float()
        return t.numpy()
    return np.asarray(t)


def _against(model: Mapping, base: Mapping) -> dict:
    """The record of one model against one base."""
    keys = list(model)
    shared = mismatch = equal = 0
    num = den = 0.0
    rel = []
    for k in keys:
        if k not in base:
            continue
        a, b = _as_array(model[k]), _as_array(base[k])
        if a.shape != b.shape:
            mismatch += 1
            continue
        shared += 1
        if a.dtype == b.dtype and np.array_equal(a, b):
            equal += 1
        a64, b64 = a.astype(np.float64, copy=False), b.astype(np.float64, copy=False)
        d2 = float(np.sum((a64 - b64) ** 2))
        n2 = float(np.sum(b64 ** 2))
        num += d2
        den += n2
        rel.append(math.sqrt(d2) / math.sqrt(n2) if n2 > 0 else (0.0 if d2 == 0 else math.inf))
    return {"n_model_tensors": len(keys), "n_shared": shared, "n_shape_mismatch": mismatch, "n_identical": equal,
            "coverage": shared / len(keys) if keys else 0.0,
            "relative_distance": (math.sqrt(num) / math.sqrt(den)) if den > 0 else (0.0 if num == 0 else math.inf),
            "median_tensor_distance": float(np.median(rel)) if rel else None}


def verdict(per_base: Mapping[str, dict], *, separation: float = 2.0, min_coverage: float = 0.9) -> dict:
    """The comparison record from per-base results: the closest well-covered base, or "not determined"."""
    covered = {n: r for n, r in per_base.items() if r["coverage"] >= min_coverage}
    ranked = sorted(covered.items(), key=lambda kv: kv[1]["relative_distance"])
    result, ratio = "not determined", None
    if len(ranked) == 1:
        result = ranked[0][0]
    elif len(ranked) >= 2:
        d1, d2 = ranked[0][1]["relative_distance"], ranked[1][1]["relative_distance"]
        ratio = math.inf if d1 == 0 else d2 / d1
        if ratio >= separation:
            result = ranked[0][0]
    return {"method": METHOD.format(cov=min_coverage, sep=separation), "bases": dict(per_base), "result": result,
            "separation_ratio": ratio, "separation_required": separation, "min_coverage": min_coverage}


def compare_state_dicts(model: Mapping, bases: Mapping[str, Mapping], *, separation: float = 2.0,
                        min_coverage: float = 0.9) -> dict:
    """The comparison record for one model against named bases (module docstring)."""
    return verdict({name: _against(model, base) for name, base in bases.items()}, separation=separation,
                   min_coverage=min_coverage)


def load_state_dict(path: str | Path) -> dict:
    """A Cellpose model file's state dict on the CPU (torch.load with weights_only)."""
    import torch

    sd = torch.load(str(path), map_location="cpu", weights_only=True)
    if isinstance(sd, Mapping) and "state_dict" in sd and isinstance(sd["state_dict"], Mapping):
        sd = sd["state_dict"]
    return dict(sd)


def compare_files(model_path: str | Path, base_paths: Mapping[str, str | Path], *, separation: float = 2.0,
                  min_coverage: float = 0.9) -> dict:
    """compare_state_dicts on model files; each base is loaded and released in turn to bound the memory."""
    model = load_state_dict(model_path)
    per_base = {}
    for name, path in base_paths.items():
        per_base[name] = _against(model, load_state_dict(path))
    return verdict(per_base, separation=separation, min_coverage=min_coverage)
