"""Leave-one-gonad-out estimates for the spot detectors (docs/BENCHMARK_PROTOCOL.md, section 6.2).

The folds. Over the train and val gonads with usable spot truth (six for the N2 RAD-51 split: the two
heat-shock males of train were exported with Imaris region growing), every gonad is held out once. A fold
tunes or trains on the other five and scores the held-out gonad with the section 5 metrics. Selection
inside a fold uses only its five gonads: the designated val gonads when they are among them, otherwise one
gonad of the same condition chosen by tag order (the first in sorted tag order among the fold's train
gonads of that condition) in place of the held-out val gonad; the remaining gonads are the fold's train
gonads. ``logo_folds`` builds them.

The tuned classical detectors. The D3 searches (docs/BENCHMARK.md, "Tuning on train and val") were staged
and data dependent: a fixed grid (stage A for ``log``, S1 for SpotMAX), a neighbourhood of the grid's best
(stage B; S2), for ``log`` a finer neighbourhood of the stage-B best (stage C), then the best few sets on
train (5 for ``log``, 3 for SpotMAX) scored on val and the best of them on val selected; the objective is
``tuning.objective`` (pooled germline per-nucleus CCC rounded to 3 decimals, then pooled spot F1). A fold
replays that whole procedure on its own gonads: the grid ranked on the fold's train gonads, the
neighbourhood generated around the fold's own best (``log_stage_b``, ``log_stage_c``,
``spotmax_stage_s2``, the rules of tuning_d3/search_plan.json), the fold's finalists scored on the fold's val
gonads. When the fold's grid best is the D3 one, its neighbourhood is exactly the D3 configurations; when it
is not, the fold needs configurations D3 never ran, and ``replay_log_search`` / ``replay_spotmax_search``
name them (``missing``) so they can be run before the selection is final.

Scores come from the run folders of ``bench_tune.py`` and ``bench_detector.py`` (``index_runs``): for every
parameter set and gonad, the germline per-nucleus counts and the spot matching counts, which is all the
objective needs. The held-out gonad is then scored afresh by ``bench_detector.py`` with the selected
parameters (so its record carries every metric, the section 5.1 partition included).

Pooling (``pool_heldout``, ``paired``): the held-out per-gonad tables of a detector are pooled over the
six folds with the cluster bootstrap by gonad of protocol section 6.3 (1000 resamples, seed 20251105), and
two detectors are compared on the same resamples (``stats.paired_cluster_bootstrap``).
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import detector as bd
from . import tuning as bt
from .stats import cluster_bootstrap, lin_ccc, paired_cluster_bootstrap, pooled_agreement, pooled_detection

LOG_STAGE_A_DIAMETERS = (0.35, 0.42, 0.5, 0.6, 0.7)
LOG_STAGE_A_SNRS = (1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
LOG_DEFAULTS = {"min_snr": 2.0, "diameter_um": 0.5, "pre_gate_sigmas": 3.0, "merge_z_columns": True,
                "separation_ratio": 0.5}
LOG_TOP = 5
SPOTMAX_TOP = 3
SPOTMAX_S2_SHAPES = ((0.25, 0.08), (0.35, 0.08), (0.3, 0.05), (0.3, 0.12))
SPOTMAX_ABBR = {"threshold_triangle": "tri", "threshold_otsu": "otsu", "threshold_li": "li",
                "threshold_yen": "yen"}


# ------------------------------------------------------------------------------------ folds
def logo_folds(entries, *, cache_root=None, require_cache: bool = True) -> list[dict]:
    """The leave-one-gonad-out folds of protocol section 6.2 for the spot detectors (module docstring).
    ``entries`` is the manifest (a list of ``Entry``); roles come from the pre-registration
    (``detector.entry_role``). Only train and val gonads with usable spot truth (and, with
    ``require_cache``, a bench cache) take part; a test gonad never does. Returns one fold per gonad in
    sorted tag order: ``held_out``, its ``role`` and ``condition``, the fold's ``train`` and ``val``
    gonads, and ``substitute`` (the train gonad that stands in for a held-out val gonad, or None)."""
    root = Path(cache_root or bd.DEFAULT_CACHE_ROOT)
    roles = {e.tag: bd.entry_role(e) for e in entries}
    usable = []
    for e in entries:
        if roles[e.tag] not in ("train", "val") or not bd.spot_truth_usable(e):
            continue
        if require_cache and not (root / e.tag / bd.CACHE_FILES["meta"]).is_file():
            continue
        usable.append(e)
    bd.guard_test([e.tag for e in usable], roles=roles)
    cond = {e.tag: str(e.extra.get("condition") or "") for e in usable}
    tags = sorted(e.tag for e in usable)
    designated = [t for t in tags if roles[t] == "val"]
    folds = []
    for h in tags:
        rest = [t for t in tags if t != h]
        val = [t for t in designated if t != h]
        substitute = None
        if roles[h] == "val":
            cands = sorted(t for t in rest if roles[t] == "train" and cond[t] == cond[h])
            if cands:
                substitute = cands[0]
                val.append(substitute)
        train = [t for t in rest if t not in val]
        folds.append({"held_out": h, "role": roles[h], "condition": cond[h], "train": train, "val": val,
                      "substitute": substitute})
    return folds


# ------------------------------------------------------------------------- the run index
@dataclass
class GonadScore:
    """What the D3 objective needs from one gonad of one run: the germline per-nucleus counts and the
    spot matching counts (Class A+B)."""

    expert: np.ndarray
    ours: np.ndarray
    n_matched: int
    n_ours: int
    n_expert: int
    source: str = ""


@dataclass
class RunIndex:
    """Scores by parameter set (``params_sha256``) and gonad, with the parameter sets and their names."""

    scores: dict = field(default_factory=dict)       # sha -> {tag: GonadScore}
    params: dict = field(default_factory=dict)       # sha -> params
    names: dict = field(default_factory=dict)        # sha -> first name seen

    def has(self, sha: str, tag: str) -> bool:
        return tag in self.scores.get(sha, {})


def _read_run_dir(run_dir: Path) -> tuple[dict | None, str | None, dict]:
    pooled = json.loads((run_dir / "pooled.json").read_text(encoding="utf-8"))
    params = pooled.get("params")
    if params is None:
        return None, None, {}
    pg = pd.read_csv(run_dir / "per_gonad.csv")
    pn = pd.read_csv(run_dir / "per_nucleus.csv")
    out = {}
    for _, r in pg.iterrows():
        tag = str(r["tag"])
        sub = pn[pn["tag"] == tag]
        out[tag] = GonadScore(expert=sub["expert_n"].to_numpy(float), ours=sub["ours_n"].to_numpy(float),
                              n_matched=int(r["n_matched"]), n_ours=int(r["n_ours"]), n_expert=int(r["n_expert"]),
                              source=str(run_dir).replace("\\", "/"))
    return params, run_dir.name, out


def index_runs(roots, *, backend: str | None = None) -> RunIndex:
    """Index the run folders under ``roots`` (each a folder of ``bench_tune.py`` set folders or one
    ``bench_detector.py`` output): every folder holding pooled.json, per_gonad.csv and per_nucleus.csv
    adds its per-gonad scores under its parameter set. ``backend`` keeps only runs of that backend. A
    gonad scored twice under the same parameter set must give the same counts (ValueError otherwise)."""
    idx = RunIndex()
    for root in roots:
        root = Path(root)
        dirs = [root] if (root / "pooled.json").is_file() else sorted(p for p in root.iterdir() if p.is_dir())
        for d in dirs:
            if not ((d / "pooled.json").is_file() and (d / "per_gonad.csv").is_file()
                    and (d / "per_nucleus.csv").is_file()):
                continue
            if backend is not None:
                pooled = json.loads((d / "pooled.json").read_text(encoding="utf-8"))
                if pooled.get("backend") not in (None, backend):
                    continue
            params, name, per_tag = _read_run_dir(d)
            if params is None:
                continue
            sha = bd.params_sha256(params)
            idx.params.setdefault(sha, params)
            idx.names.setdefault(sha, name)
            slot = idx.scores.setdefault(sha, {})
            for tag, s in per_tag.items():
                if tag in slot:
                    o = slot[tag]
                    same = (o.n_matched, o.n_ours, o.n_expert) == (s.n_matched, s.n_ours, s.n_expert) and \
                        np.array_equal(o.expert, s.expert) and np.array_equal(o.ours, s.ours)
                    if not same:
                        raise ValueError(f"{tag} under {name}: {o.source} and {s.source} disagree")
                    continue
                slot[tag] = s
    return idx


# ------------------------------------------------------------------------------ objective
def objective_on(idx: RunIndex, sha: str, tags) -> tuple[float, float]:
    """The D3 objective (``tuning.objective``) of one parameter set pooled over ``tags``: the germline
    per-nucleus CCC over the concatenated nuclei, rounded to ``tuning.CCC_DECIMALS``, then the F1 of the
    summed matching counts. KeyError when a gonad has no score for the set."""
    s = [idx.scores[sha][t] for t in tags]
    e = np.concatenate([g.expert for g in s]) if s else np.zeros(0)
    o = np.concatenate([g.ours for g in s]) if s else np.zeros(0)
    ccc = lin_ccc(e, o)
    det = pooled_detection([(g.n_matched, g.n_ours, g.n_expert) for g in s])
    pooled = {bt.OBJECTIVE[0]: {bt.OBJECTIVE[1]: {"point": ccc}}, bt.TIE_BREAK[0]: {bt.TIE_BREAK[1]: {"point": det["f1"]}}}
    return bt.objective(pooled)


def rank_on(idx: RunIndex, configs: list[dict], tags) -> list[dict]:
    """``configs`` sorted best first by the objective pooled over ``tags`` (stable: ties keep their order),
    each with its ``sha`` and ``objective``. Every set needs a score on every gonad (``missing_scores``
    lists the gaps; KeyError otherwise)."""
    items = []
    for c in configs:
        sha = bd.params_sha256(c["params"])
        items.append({**c, "sha": sha, "objective": objective_on(idx, sha, tags)})
    return sorted(items, key=lambda it: it["objective"], reverse=True)


def missing_scores(idx: RunIndex, configs: list[dict], tags) -> list[dict]:
    """The (set, gonad) pairs of ``configs`` x ``tags`` with no score in the index."""
    out = []
    for c in configs:
        sha = bd.params_sha256(c["params"])
        lack = [t for t in tags if not idx.has(sha, t)]
        if lack:
            out.append({"name": c["name"], "params": c["params"], "tags": lack})
    return out


# ----------------------------------------------------------------- the D3 stages, replayed
def _fmt(v: float) -> str:
    return f"{float(v):g}"


def log_point(params: dict) -> dict:
    """(diameter, SNR, separation ratio, pre-gate, z merge) of a ``log`` parameter set, with the code
    defaults (separation = half the diameter, pre-gate 3, z merge on)."""
    lg = dict(params.get("log") or {})
    d = float(lg.get("diameter_um", LOG_DEFAULTS["diameter_um"]))
    sep = lg.get("min_separation_um")
    return {"d": d, "s": float(lg.get("min_snr", LOG_DEFAULTS["min_snr"])),
            "k": round(float(sep) / d, 6) if sep is not None else LOG_DEFAULTS["separation_ratio"],
            "g": float(lg.get("pre_gate_sigmas", LOG_DEFAULTS["pre_gate_sigmas"])),
            "m": bool(params.get("merge_z_columns", LOG_DEFAULTS["merge_z_columns"]))}


def log_config(stage: str, d: float, s: float, k: float = 0.5, g: float = 3.0, m: bool = True) -> dict:
    """One ``log`` parameter set as the D3 configuration files write it: keys at their code default are
    left out (separation at half the diameter, pre-gate 3, z merge on)."""
    d, s = round(float(d), 4), round(float(s), 4)
    lg = {"diameter_um": d, "min_snr": s}
    if float(k) != 0.5:
        lg["min_separation_um"] = round(float(k) * d, 6)
    if float(g) != 3.0:
        lg["pre_gate_sigmas"] = float(g)
    params: dict = {"log": lg}
    if not m:
        params["merge_z_columns"] = False
    name = f"log{stage}_d{_fmt(d)}_s{_fmt(s)}" if stage == "A" else \
        f"log{stage}_d{_fmt(d)}_s{_fmt(s)}_k{_fmt(k)}_g{_fmt(g)}_m{int(bool(m))}"
    return {"name": name, "params": params}


def log_stage_a() -> list[dict]:
    """Stage A: diameter x SNR at the code defaults otherwise (30 sets)."""
    return [log_config("A", d, s) for d in LOG_STAGE_A_DIAMETERS for s in LOG_STAGE_A_SNRS]


def log_stage_b(d: float, s: float) -> list[dict]:
    """Stage B at the stage-A best (d, s): separation {0.5, 0.75, 1.0} x d, pre-gate {2, 3, 4}, z merge on and
    off, less the stage-A point (17 sets)."""
    out = []
    for k in (0.5, 0.75, 1.0):
        for g in (2.0, 3.0, 4.0):
            for m in (True, False):
                if k == 0.5 and g == 3.0 and m:
                    continue
                out.append(log_config("B", d, s, k, g, m))
    return out


def log_stage_c(point: dict) -> list[dict]:
    """Stage C at the stage-B best: diameter {d - 0.04, d, d + 0.04} x SNR {s - 0.25, s, s + 0.25} inside the
    stage-A ranges, with the best's separation ratio, pre-gate and z merge, less the centre (at most 8)."""
    out = []
    for dd in (point["d"] - 0.04, point["d"], point["d"] + 0.04):
        for ss in (point["s"] - 0.25, point["s"], point["s"] + 0.25):
            dd_r, ss_r = round(dd, 4), round(ss, 4)
            if dd_r == round(point["d"], 4) and ss_r == round(point["s"], 4):
                continue
            if not (min(LOG_STAGE_A_DIAMETERS) <= dd_r <= max(LOG_STAGE_A_DIAMETERS)
                    and min(LOG_STAGE_A_SNRS) <= ss_r <= max(LOG_STAGE_A_SNRS)):
                continue
            out.append(log_config("C", dd_r, ss_r, point["k"], point["g"], point["m"]))
    return out


def spotmax_config(stage: str, method: str, r: float, g: float, e: float) -> dict:
    params = {"thresholding_method": method, "spot_radius_um": float(r), "gauss_sigma_um": float(g),
              "effect_size_min": float(e)}
    name = f"sm{stage}_{SPOTMAX_ABBR.get(method, method)}_r{_fmt(r)}_g{_fmt(g)}_e{_fmt(e)}"
    return {"name": name, "params": params}


def spotmax_stage_s1() -> list[dict]:
    """S1: thresholding method x effect size at the deployed radius 0.3 um and smoothing 0.08 um (12 sets)."""
    return [spotmax_config("S1", m, 0.3, 0.08, e) for m in ("threshold_triangle", "threshold_otsu", "threshold_li")
            for e in (2.0, 3.0, 4.0, 5.0)]


def spotmax_stage_s2(method: str, e: float) -> list[dict]:
    """S2 at the S1 best (method, effect size e): (radius, smoothing) in ``SPOTMAX_S2_SHAPES`` x effect size
    {e - 0.5, e, e + 0.5} (12 sets)."""
    return [spotmax_config("S2", method, r, g, ee) for r, g in SPOTMAX_S2_SHAPES for ee in (e - 0.5, e, e + 0.5)]


def _same_sets(a: list[dict], b: list[dict]) -> bool:
    return {bd.params_sha256(c["params"]) for c in a} == {bd.params_sha256(c["params"]) for c in b}


def replay_log_search(idx: RunIndex, train: list[str], val: list[str], *, top: int = LOG_TOP) -> dict:
    """The D3 ``log`` search replayed on one fold's gonads (module docstring). Returns the stage bests, the
    fold's own stage-B and stage-C sets, whether they are the D3 ones, the finalists, the selected set,
    and ``missing`` (set and gonad pairs that must be run first; the selection is None until they are)."""
    a = log_stage_a()
    miss = missing_scores(idx, a, train)
    if miss:
        return {"complete": False, "missing": miss, "stage": "A"}
    ra = rank_on(idx, a, train)
    pa = log_point(ra[0]["params"])
    b = log_stage_b(pa["d"], pa["s"])
    miss = missing_scores(idx, b, train)
    if miss:
        return {"complete": False, "missing": miss, "stage": "B", "stage_a_best": ra[0]["name"]}
    rb = rank_on(idx, b + [ra[0]], train)
    pb = log_point(rb[0]["params"])
    c = log_stage_c(pb)
    miss = missing_scores(idx, c, train)
    if miss:
        return {"complete": False, "missing": miss, "stage": "C", "stage_a_best": ra[0]["name"],
                "stage_b_best": rb[0]["name"]}
    allsets = a + b + c
    rall = rank_on(idx, allsets, train)
    finalists = rall[:top]
    miss = missing_scores(idx, finalists, val)
    out = {"stage_a_best": ra[0]["name"], "stage_b_best": rb[0]["name"],
           "stage_b": [x["name"] for x in b], "stage_c": [x["name"] for x in c],
           "n_sets": len(allsets), "finalists": [x["name"] for x in finalists],
           "train_objective": {x["name"]: list(x["objective"]) for x in finalists}}
    if miss:
        return {**out, "complete": False, "missing": miss, "stage": "val"}
    rv = rank_on(idx, [{"name": x["name"], "params": x["params"]} for x in finalists], val)
    return {**out, "complete": True, "missing": [], "selected": rv[0]["name"], "selected_params": rv[0]["params"],
            "val_objective": {x["name"]: list(x["objective"]) for x in rv}}


def replay_spotmax_search(idx: RunIndex, train: list[str], val: list[str], *, top: int = SPOTMAX_TOP) -> dict:
    """The D3 SpotMAX search replayed on one fold's gonads (module docstring); the same returns as
    ``replay_log_search``."""
    s1 = spotmax_stage_s1()
    miss = missing_scores(idx, s1, train)
    if miss:
        return {"complete": False, "missing": miss, "stage": "S1"}
    r1 = rank_on(idx, s1, train)
    best = r1[0]["params"]
    s2 = spotmax_stage_s2(best["thresholding_method"], float(best["effect_size_min"]))
    miss = missing_scores(idx, s2, train)
    if miss:
        return {"complete": False, "missing": miss, "stage": "S2", "stage_s1_best": r1[0]["name"]}
    allsets = s1 + s2
    rall = rank_on(idx, allsets, train)
    finalists = rall[:top]
    miss = missing_scores(idx, finalists, val)
    out = {"stage_s1_best": r1[0]["name"], "stage_s2": [x["name"] for x in s2], "n_sets": len(allsets),
           "finalists": [x["name"] for x in finalists],
           "train_objective": {x["name"]: list(x["objective"]) for x in finalists}}
    if miss:
        return {**out, "complete": False, "missing": miss, "stage": "val"}
    rv = rank_on(idx, [{"name": x["name"], "params": x["params"]} for x in finalists], val)
    return {**out, "complete": True, "missing": [], "selected": rv[0]["name"], "selected_params": rv[0]["params"],
            "val_objective": {x["name"]: list(x["objective"]) for x in rv}}


# ------------------------------------------------------------------------------- pooling
POOL_BLOCKS = {
    "truth_masks": ("per_truth_nucleus", "expert_n", "ours_n"),
    "truth_masks_class_a": ("per_truth_nucleus", "expert_n_class_a", "ours_n"),
    "truth_masks_detector_only": ("per_truth_nucleus", "expert_n", "ours_n_detector_only"),
    "germline": ("per_nucleus", "expert_n", "ours_n"),
    "traced": ("per_nucleus_traced", "expert_n", "ours_n"),
}
DETECTION_BLOCKS = {"spots_all": "", "spots_without_class_b": "no_b_", "spots_detector_only": "det_"}


@dataclass
class Heldout:
    """One detector's held-out scores: one per-gonad row per fold and the per-nucleus tables."""

    name: str
    per_gonad: pd.DataFrame
    per_nucleus: pd.DataFrame
    per_truth_nucleus: pd.DataFrame
    sources: dict = field(default_factory=dict)

    def payloads(self, block: str, tags) -> list:
        if block in POOL_BLOCKS:
            table, e_col, o_col = POOL_BLOCKS[block]
            if table == "per_truth_nucleus":
                df = self.per_truth_nucleus
            else:
                df = self.per_nucleus
                if table == "per_nucleus_traced":
                    df = df[df["traced"].astype(str).str.lower().isin(["true", "1"])]
            if e_col not in df.columns or o_col not in df.columns:
                return []
            return [(df.loc[df["tag"] == t, e_col].to_numpy(float), df.loc[df["tag"] == t, o_col].to_numpy(float))
                    for t in tags]
        prefix = DETECTION_BLOCKS[block]
        rows = self.per_gonad.set_index("tag")
        need = [f"{prefix}n_matched", f"{prefix}n_ours", f"{prefix}n_expert"]
        if any(c not in rows.columns for c in need):
            return []
        return [tuple(int(rows.loc[t, c]) for c in need) for t in tags]


def load_heldout(name: str, run_dirs) -> Heldout:
    """A detector's held-out scores from ``bench_detector.py`` output folders (one per fold, or one folder
    holding several gonads); each gonad must appear once."""
    pg, pn, ptn, src = [], [], [], {}
    for d in run_dirs:
        d = Path(d)
        g = pd.read_csv(d / "per_gonad.csv")
        n = pd.read_csv(d / "per_nucleus.csv")
        t = pd.read_csv(d / "per_truth_nucleus.csv") if (d / "per_truth_nucleus.csv").is_file() else pd.DataFrame()
        for tag in g["tag"]:
            if tag in src:
                raise ValueError(f"{name}: {tag} appears in {src[tag]} and in {d}")
            src[str(tag)] = str(d).replace("\\", "/")
        pg.append(g)
        pn.append(n)
        ptn.append(t)
    return Heldout(name, pd.concat(pg, ignore_index=True), pd.concat(pn, ignore_index=True),
                   pd.concat(ptn, ignore_index=True), src)


def _stat_fn(block: str):
    return pooled_agreement if block in POOL_BLOCKS else pooled_detection


def pool_heldout(h: Heldout, tags, *, n_boot: int = bd.N_BOOT, seed: int = bd.BOOT_SEED) -> dict:
    """The pooled leave-one-gonad-out estimate of one detector over the held-out gonads ``tags``: for every
    block, the pooled statistics with their cluster-bootstrap intervals by gonad (protocol section 6.3)."""
    out = {}
    for block in list(POOL_BLOCKS) + list(DETECTION_BLOCKS):
        pay = h.payloads(block, tags)
        if not pay or (block in POOL_BLOCKS and not any(len(p[0]) for p in pay)):
            continue
        boot = cluster_bootstrap(pay, _stat_fn(block), n_boot=n_boot, seed=seed)["stats"]
        rec = {"n_gonads": len(tags), "stats": boot}
        if block in POOL_BLOCKS:
            rec["n_nuclei"] = int(sum(len(p[0]) for p in pay))
        else:
            rec["n_matched"], rec["n_ours"], rec["n_expert"] = (int(sum(p[i] for p in pay)) for i in range(3))
        out[block] = rec
    return out


def paired(a: Heldout, b: Heldout, tags, *, n_boot: int = bd.N_BOOT, seed: int = bd.BOOT_SEED) -> dict:
    """Detector a against detector b over the same held-out gonads, every block, both statistics on the same
    resample (``stats.paired_cluster_bootstrap``)."""
    out = {}
    for block in list(POOL_BLOCKS) + list(DETECTION_BLOCKS):
        pa, pb = a.payloads(block, tags), b.payloads(block, tags)
        if not pa or not pb:
            continue
        out[block] = paired_cluster_bootstrap(pa, pb, _stat_fn(block), n_boot=n_boot, seed=seed)["stats"]
    return out


def per_gonad_table(h: Heldout) -> pd.DataFrame:
    """The headline held-out numbers per gonad: truth-mask (section 5.1) CCC, bias and MAE (end-to-end, and
    the detector-only CCC when measured), germline and traced CCC, spot F1 with and without Class B."""
    cols = ["tag", "split", "role", "tm_n_nuclei", "tm_ccc", "tm_bias", "tm_mae", "tm_class_a_ccc", "tm_det_ccc",
            "ccc_all", "bias_all", "ccc_traced", "bias_traced", "f1", "precision", "recall", "no_b_f1", "det_f1",
            "n_our_spots", "n_expert_spots"]
    keep = [c for c in cols if c in h.per_gonad.columns]
    return h.per_gonad[keep].copy()


def finite(x) -> float | None:
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None
