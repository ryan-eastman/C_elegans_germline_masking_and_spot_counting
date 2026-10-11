"""Parameter search on the detector bench harness: many parameter sets of one spot detection backend
run on the cached gonads of one split (train or val), each one scored and logged exactly like a
``scripts/bench_detector.py`` run, then ranked by one fixed objective. ``scripts/bench_tune.py`` is
the command line; docs/BENCHMARK.md ("Tuning on train and val") has the searches made with it.

Why a separate runner. ``bench_detector.py`` runs one parameter set per call and reloads every gonad's
cache (about 3 GB) each time. A search over tens of parameter sets keeps each gonad in one worker
process and runs all of its parameter sets there. The detection and the scoring are the harness's own
(``detector.run_backend`` and ``detector.score``), the pooling is ``detector.pool_scores`` and every
parameter set appends one ``detector.trial_record`` line to trials.jsonl, with a ``tuning`` block that
names the search and the stage it belongs to, so the protocol's tuning budget is counted from the same
log. A tuning run never touches a test gonad: there is no way to pass ``allow_test`` here, and
``scripts/bench_detector.py --final`` stays the only path to a test score.

The SpotMAX effect-size filter. ``detect_spots_spotmax`` applies ``effect_size_min`` as its very last
step, after detection, the features and the z-split merge, as a plain row filter: a spot whose effect
size is below the minimum is dropped, a missing effect size counts as passing, and nothing is filtered
when no spot has an effect size or the minimum is 0. Parameter sets that differ only in
``effect_size_min`` therefore share everything before that filter. ``run_configs`` runs SpotMAX once
per group of such sets with ``effect_size_min`` 0 and applies each set's minimum afterwards with
``filter_effect_size``, which rebuilds the two tables exactly as ``detect_spots_spotmax`` builds them
(spot ids renumbered from 0, per-nucleus counts and mean intensities recomputed, one row per nucleus).
The equality is tested on synthetic data (tests/test_bench_tuning.py) and was checked against the
golden spots tables of the real gonads.

The objective. The pooled per-nucleus Lin's CCC over the germline nuclei of the split (Class A+B truth,
end-to-end: the detector's own spots and counts on the deployed nucleus masks; ``per_nucleus_all`` in
pooled.json), rounded to ``CCC_DECIMALS`` decimals, with the pooled spot F1 (``spots_all``) breaking
ties. A missing value ranks last.
"""
from __future__ import annotations

import copy
import math
import re
import time

import numpy as np
import pandas as pd

from . import detector as bd
from .manifest import load_labels, load_truth

CCC_DECIMALS = 3
OBJECTIVE = ("per_nucleus_all", "ccc")       # pooled germline per-nucleus CCC, Class A+B, end-to-end
TIE_BREAK = ("spots_all", "f1")              # pooled spot F1, Class A+B
OBJECTIVE_TEXT = (f"pooled per_nucleus_all CCC rounded to {CCC_DECIMALS} decimals, then pooled spots_all F1 "
                  "(higher is better; a missing value ranks last)")
_NAME_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.+=-]*")


# ------------------------------------------------------------------------------ configurations
def check_configs(configs) -> list[dict]:
    """Validate a list of parameter sets: each a mapping with a ``name`` (letters, digits and
    ``_ . + = -``, unique; it names the output folder) and ``params`` (a mapping of spots: keys without
    ``backend``). Returns deep copies."""
    if not isinstance(configs, list) or not configs:
        raise ValueError("configs must be a non-empty list of {name, params}")
    out, seen = [], set()
    for c in configs:
        if not isinstance(c, dict) or "name" not in c:
            raise ValueError(f"every config needs a name: {c!r}")
        name = str(c["name"])
        if not _NAME_RE.fullmatch(name):
            raise ValueError(f"config name {name!r} is not usable as a folder name")
        if name in seen:
            raise ValueError(f"config name {name!r} appears twice")
        seen.add(name)
        params = copy.deepcopy(c.get("params") or {})
        if not isinstance(params, dict):
            raise ValueError(f"{name}: params must be a mapping of spots: keys")
        if "backend" in params:
            raise ValueError(f"{name}: pass the backend with --backend, not inside params")
        out.append({"name": name, "params": params})
    return out


def effect_size_group_key(spots_config: dict | None, params: dict) -> str:
    """The identity of a SpotMAX run apart from its effect-size filter: the effective spots: section
    (golden section with ``params`` laid over it) without ``effect_size_min``, as canonical JSON."""
    section = bd.effective_section(spots_config, params)
    section.pop("effect_size_min", None)
    return bd.canonical_json(section)


def effective_effect_size_min(spots_config: dict | None, params: dict) -> float:
    """The ``effect_size_min`` a run with ``params`` uses, with the pipeline's code default."""
    return bd.pipeline_spot_kwargs(bd.effective_section(spots_config, params))["effect_size_min"]


def filter_effect_size(per_spot: pd.DataFrame, per_nuc: pd.DataFrame, effect_size_min: float):
    """Apply SpotMAX's effect-size filter to the tables of a run made with ``effect_size_min`` 0, exactly
    as ``detect_spots_spotmax`` applies it before building its tables: rows whose effect size (a missing
    one counting as infinite) is below the minimum are dropped, in order, when the minimum is above 0
    and at least one spot has an effect size; then the spot ids are renumbered from 0 and the
    per-nucleus counts and mean spot intensities are recomputed over every nucleus of ``per_nuc`` with
    the expressions of ``detect._finish_tables``. Returns (per_spot, per_nucleus)."""
    from germquant.spots.backends import PER_NUC_COLS

    e = float(effect_size_min)
    ps = per_spot
    if e > 0 and len(ps) and ps["effect_size"].notna().any():
        keep = ps["effect_size"].astype(float).fillna(np.inf).to_numpy() >= e
        ps = ps.loc[keep]
    ps = ps.reset_index(drop=True)
    if "spot_id" in ps.columns:
        ps = ps.assign(spot_id=np.arange(len(ps)))
    if not len(per_nuc):
        return ps, per_nuc.copy()
    marker = per_nuc["marker"].iloc[0]
    detector = per_nuc["detector"].iloc[0]
    nuc_ids = [int(v) for v in per_nuc["nucleus_id"]]
    counts = ps.groupby("nucleus_id").size()
    means = ps.groupby("nucleus_id")["spot_mean_intensity"].mean()
    pn = pd.DataFrame({
        "nucleus_id": nuc_ids,
        "marker": marker,
        "n_spots": [int(counts.get(n, 0)) for n in nuc_ids],
        "mean_spot_intensity": [float(means.get(n, np.nan)) for n in nuc_ids],
        "detector": detector,
    }, columns=PER_NUC_COLS)
    return ps, pn


def chunk_configs(configs: list[dict], n_chunks: int, *, backend: str = "",
                  spots_config: dict | None = None) -> list[list[dict]]:
    """Split ``configs`` into at most ``n_chunks`` lists of about equal size, in order. For SpotMAX the
    sets that share a detection (``effect_size_group_key``) stay in one list, so the shared run is
    made once."""
    n = max(1, int(n_chunks))
    if backend == "spotmax":
        groups: dict[str, list[dict]] = {}
        for c in configs:
            groups.setdefault(effect_size_group_key(spots_config, c["params"]), []).append(c)
        units = list(groups.values())
    else:
        units = [[c] for c in configs]
    n = min(n, len(units))
    sizes = [len(u) for u in units]
    total = sum(sizes)
    chunks: list[list[dict]] = [[] for _ in range(n)]
    k, filled = 0, 0
    for u in units:
        if chunks[k] and filled >= total * (k + 1) / n and k < n - 1:
            k += 1
        chunks[k].extend(u)
        filled += len(u)
    return [c for c in chunks if c]


# ---------------------------------------------------------------------------------- running
def run_configs(tag: str, backend: str, configs: list[dict], *, cache_root=None, manifest=None,
                score_kwargs: dict | None = None, effect_size_shortcut: bool = True,
                verbose: bool = True) -> list[dict]:
    """Run and score every parameter set of ``configs`` (``{"name", "params"}``) on one gonad in this
    process. The cache is loaded once; every set is detected first (SpotMAX sets that differ only in
    ``effect_size_min`` share one detection when ``effect_size_shortcut``, see ``filter_effect_size``),
    then the foci channel is released, the truth is loaded once and every set is scored with
    ``detector.score``. Returns one record per set, in the order given: ``name``, ``params``, ``row``
    (the per-gonad row of ``bench_detector.py``), ``per_nucleus``, ``matches``, ``per_spot``,
    ``shared_detection`` (the parameter set that was actually run, when the set was derived from a
    shared SpotMAX run, else None) and ``per_truth_nucleus`` (the protocol's section 5.1 partition). A
    test gonad is refused."""
    bd.guard_test([tag])               # the pre-registered test gonads, before anything is read
    configs = check_configs(configs)
    t_start = time.time()
    entries = bd._entries_by_tag(manifest)
    if tag not in entries:
        raise KeyError(f"{tag} is not in the manifest")
    entry = entries[tag]
    bd.guard_test([tag], roles={tag: bd.entry_role(entry)})          # no allow_test: tuning never scores test
    if not bd.spot_truth_usable(entry):
        raise ValueError(f"{tag}: the spot truth is unusable (exported with Imaris region growing)")
    cache = bd.load_cache(tag, cache_root, arrays=("foci", "labels"), mmap=False)
    golden = cache.spots_config

    def say(msg: str) -> None:
        if verbose:
            print(f"  [{tag}] {msg}", flush=True)

    # groups of sets that share one detection (SpotMAX only), in first-appearance order
    groups: dict[str, list[dict]] = {}
    for c in configs:
        key = effect_size_group_key(golden, c["params"]) if (backend == "spotmax" and effect_size_shortcut) \
            else c["name"]
        groups.setdefault(key, []).append(c)

    detections: dict[str, dict] = {}
    for group in groups.values():
        if len(group) > 1:
            base = {k: v for k, v in group[0]["params"].items() if k != "effect_size_min"}
            base["effect_size_min"] = 0.0
            t = time.time()
            ps0, pn0 = bd.run_backend(tag, backend, base, cache=cache)
            dt = time.time() - t
            say(f"shared {backend} run {bd.canonical_json(base)}: {len(ps0)} spots before the effect-size "
                f"filter ({dt:.0f} s)")
            for c in group:
                e = effective_effect_size_min(golden, c["params"])
                ps, pn = filter_effect_size(ps0, pn0, e)
                detections[c["name"]] = {"per_spot": ps, "per_nucleus": pn, "detect_s": dt, "shared": base}
            del ps0, pn0
        else:
            c = group[0]
            t = time.time()
            ps, pn = bd.run_backend(tag, backend, c["params"], cache=cache)
            dt = time.time() - t
            detections[c["name"]] = {"per_spot": ps, "per_nucleus": pn, "detect_s": dt, "shared": None}
            say(f"{c['name']}: {len(ps)} spots ({dt:.0f} s)")
    cache.foci = None

    truth = load_truth(entry.truth_dir, entry.truth_stem)
    if truth is None:
        raise FileNotFoundError(f"{tag}: no truth files for stem {entry.truth_stem!r} in {entry.truth_dir}")
    truth_labels = load_labels(truth.labels_path) if truth.has_labels else None
    entry_list = list(entries.values())
    out = []
    for c in configs:
        d = detections[c["name"]]
        t = time.time()
        row, pn, matches, per_tm = bd.score(tag, d["per_spot"], d["per_nucleus"], manifest=entry_list, cache=cache,
                                            truth=truth, truth_labels=truth_labels, return_truth_masks=True,
                                            **(score_kwargs or {}))
        row["backend"] = backend
        row["detect_s"] = round(d["detect_s"], 1)
        row["score_s"] = round(time.time() - t, 1)
        n_nuc = d["per_nucleus"]
        row["n_spots_reconcile"] = (bool(int(n_nuc["n_spots"].sum()) == len(d["per_spot"])) if len(n_nuc)
                                    else len(d["per_spot"]) == 0)
        out.append({"tag": tag, "name": c["name"], "params": c["params"], "row": row, "per_nucleus": pn,
                    "matches": matches, "per_spot": d["per_spot"], "shared_detection": d["shared"],
                    "per_truth_nucleus": per_tm})
        say(f"{c['name']}: ours {row['n_our_spots']}, germline CCC {row['ccc_all']:.3f}, F1 {row['f1']:.3f}")
    say(f"done, {len(configs)} parameter sets in {time.time() - t_start:.0f} s")
    return out


# ---------------------------------------------------------------------------------- ranking
def _point(pooled: dict | None, block: str, stat: str) -> float | None:
    b = (pooled or {}).get(block) or {}
    s = b.get(stat)
    v = s.get("point") if isinstance(s, dict) else s
    return bd.finite_or_none(v)


def objective(pooled: dict | None) -> tuple[float, float]:
    """The ranking key of a pooled result: (CCC rounded to CCC_DECIMALS, F1), missing values as -inf."""
    ccc = _point(pooled, *OBJECTIVE)
    f1 = _point(pooled, *TIE_BREAK)
    return (round(ccc, CCC_DECIMALS) if ccc is not None else -math.inf,
            f1 if f1 is not None else -math.inf)


def rank(items: list, pooled_of=lambda it: it["pooled"]) -> list:
    """``items`` sorted best first by ``objective`` (stable: equal keys keep their order)."""
    return sorted(items, key=lambda it: objective(pooled_of(it)), reverse=True)


SUMMARY_FIELDS = [
    ("ccc_all", "per_nucleus_all", "ccc"), ("ccc_all_lo", "per_nucleus_all", "ccc", "lo"),
    ("ccc_all_hi", "per_nucleus_all", "ccc", "hi"), ("bias_all", "per_nucleus_all", "bias"),
    ("r_all", "per_nucleus_all", "pearson_r"),
    ("ccc_traced", "per_nucleus_traced", "ccc"), ("bias_traced", "per_nucleus_traced", "bias"),
    ("ccc_all_class_a", "per_nucleus_all_class_a", "ccc"), ("ccc_traced_class_a", "per_nucleus_traced_class_a", "ccc"),
    ("f1", "spots_all", "f1"), ("precision", "spots_all", "precision"), ("recall", "spots_all", "recall"),
    ("no_b_f1", "spots_without_class_b", "f1"), ("det_f1", "spots_detector_only", "f1"),
    ("tm_ccc", "per_nucleus_truth_masks", "ccc"), ("tm_bias", "per_nucleus_truth_masks", "bias"),
    ("tm_mae", "per_nucleus_truth_masks", "mae"), ("tm_ccc_class_a", "per_nucleus_truth_masks_class_a", "ccc"),
]


def summary_row(name: str, params: dict, pooled: dict) -> dict:
    """One flat row of headline numbers for a pooled result (summary.csv of bench_tune.py)."""
    row = {"name": name, "params": bd.canonical_json(params), "params_sha256": bd.params_sha256(params)}
    for spec in SUMMARY_FIELDS:
        col, block, stat = spec[:3]
        which = spec[3] if len(spec) > 3 else "point"
        b = (pooled.get(block) or {}).get(stat)
        row[col] = bd.finite_or_none(b.get(which)) if isinstance(b, dict) else None
    row["n_ours"] = (pooled.get("spots_all") or {}).get("n_ours")
    row["n_expert"] = (pooled.get("spots_all") or {}).get("n_expert")
    row["n_gonads"] = pooled.get("n_gonads")
    key = objective(pooled)
    row["objective_ccc"] = key[0] if math.isfinite(key[0]) else None
    return row


def budget_used(records: list[dict], search: str) -> set[str]:
    """The distinct parameter sets (params_sha256) that trial records of the search ``search`` hold."""
    return {str(r.get("params_sha256")) for r in records
            if isinstance(r.get("tuning"), dict) and r["tuning"].get("search") == search}
