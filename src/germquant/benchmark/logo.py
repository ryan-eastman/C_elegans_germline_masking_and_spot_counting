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

The pre-registered SpotMAX grid (protocol section 8 item 2), which D3 did not run. ``prereg_spotmax_grid``
lists its 60 sets: thresholding method (``PREREG_SPOTMAX_METHODS``) by ``effect_size_min``
(``PREREG_SPOTMAX_EFFECT_SIZES``) by ``z_merge_gap_um`` (``PREREG_SPOTMAX_Z_GAPS``). Every other spots key
stays at the golden run's deployed value, so pg_tri_e3.0_z0.8 is the deployed recipe. The rule is that of
sections 6.1 and 6.4, not D3's: every set is scored on the selection gonads and the best on the primary
metric is selected. There is no train stage, so every set is a finalist. ``primary_objective_on`` is that
metric, pooled as ``bench_detector.py`` pools it (the point estimate over the pooled nuclei of the given
gonads). It is the per-nucleus CCC on the section 5.1 partition (``per_nucleus_truth_masks``: the truth
masks, Class A+B truth, end-to-end counts) rounded to 3 decimals, then the section 5.2 spot F1 inside R
(``spots_truth_masks``). A missing or NaN value ranks last. ``select_grid`` ranks the sets on given gonads
(ties keep the grid's order) and names the (set, gonad) pairs still to run. A score without the truth-mask
counts is one of those pairs, so a set is never passed over only because its run folder lacks them.
``replay_grid`` applies the rule to one fold's selection gonads (its ``val``, from ``logo_folds``). It
returns what ``replay_spotmax_search`` returns, so ``bench_logo.py heldout --selection`` reads its output
unchanged. ``designated_split`` gives the split's own train and val gonads; the frozen grid set is chosen
on those val gonads. The index carries the truth-mask counts only when it is built with
``index_runs(..., truth_masks=True)``. The default index, and so the D3 replay, reads nothing more than
before.

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

# protocol section 8 item 2: the pre-registered SpotMAX grid (4 x 5 x 3 = 60 sets), in the protocol's order
PREREG_SPOTMAX_METHODS = ("threshold_otsu", "threshold_li", "threshold_triangle", "threshold_yen")
PREREG_SPOTMAX_EFFECT_SIZES = (2.0, 2.5, 3.0, 3.5, 4.0)
PREREG_SPOTMAX_Z_GAPS = (0.6, 0.8, 1.0)
PREREG_SPOTMAX_SEARCH = "prereg_spotmax_grid"          # the bench_tune.py --search label of its runs
# the primary metric (protocol section 5.1, Class A+B truth, end-to-end; changelog entry 4; the summary column
# tm_ccc of bench_tune.py) and its tie-break, the section 5.2 spot F1 inside R
PRIMARY_OBJECTIVE = ("per_nucleus_truth_masks", "ccc")
PRIMARY_TIE_BREAK = ("spots_truth_masks", "f1")
PRIMARY_OBJECTIVE_TEXT = (f"pooled per_nucleus_truth_masks CCC (protocol section 5.1 partition, Class A+B truth, "
                          f"end-to-end counts) rounded to {bt.CCC_DECIMALS} decimals, then pooled spots_truth_masks "
                          "F1 (section 5.2, matching inside R), then the grid's order (higher is better; a missing "
                          "value ranks last)")
TM_SPOTS_COLUMNS = ("tm_spots_n_matched", "tm_spots_n_ours", "tm_spots_n_expert")

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
    spot matching counts (Class A+B). For the primary objective, read only by ``index_runs(...,
    truth_masks=True)``: the per-truth-nucleus counts of the section 5.1 partition (``tm_expert``,
    ``tm_ours``, from per_truth_nucleus.csv) and the section 5.2 matching counts inside R (``tm_spots``:
    matched, ours, expert, from the tm_spots_ columns of per_gonad.csv). None where the run folder does
    not have them."""

    expert: np.ndarray
    ours: np.ndarray
    n_matched: int
    n_ours: int
    n_expert: int
    source: str = ""
    tm_expert: np.ndarray | None = None
    tm_ours: np.ndarray | None = None
    tm_spots: tuple[int, int, int] | None = None

    def has_truth_masks(self) -> bool:
        """True when the score has everything ``primary_objective_on`` pools: the truth-mask counts and the
        matching counts inside R."""
        return self.tm_expert is not None and self.tm_ours is not None and self.tm_spots is not None


@dataclass
class RunIndex:
    """Scores by parameter set (``params_sha256``) and gonad, with the parameter sets and their names."""

    scores: dict = field(default_factory=dict)       # sha -> {tag: GonadScore}
    params: dict = field(default_factory=dict)       # sha -> params
    names: dict = field(default_factory=dict)        # sha -> first name seen

    def has(self, sha: str, tag: str) -> bool:
        return tag in self.scores.get(sha, {})


def _tm_spot_counts(r) -> tuple[int, int, int] | None:
    """The matching counts inside R (matched, ours, expert) of one per_gonad.csv row, pooled as
    ``spots._detection_block`` pools them. That block leaves out a row whose n_expert + n_ours is not above 0
    (a NaN count included), so such a row is (0, 0, 0) here and adds nothing. None for a row with spots in R
    but no readable matched count, which the scorer never writes; the pair then counts as missing."""
    m, o, e = (finite(r[c]) for c in TM_SPOTS_COLUMNS)
    if o is None or e is None or not o + e > 0:
        return (0, 0, 0)
    return None if m is None else (int(m), int(o), int(e))


def _read_run_dir(run_dir: Path, truth_masks: bool = False) -> tuple[dict | None, str | None, dict]:
    pooled = json.loads((run_dir / "pooled.json").read_text(encoding="utf-8"))
    params = pooled.get("params")
    if params is None:
        return None, None, {}
    pg = pd.read_csv(run_dir / "per_gonad.csv")
    pn = pd.read_csv(run_dir / "per_nucleus.csv")
    ptn = None
    if truth_masks and (run_dir / "per_truth_nucleus.csv").is_file():
        ptn = pd.read_csv(run_dir / "per_truth_nucleus.csv")
        if not {"tag", "expert_n", "ours_n"} <= set(ptn.columns):
            ptn = None
    has_tm_spots = all(c in pg.columns for c in TM_SPOTS_COLUMNS)
    out = {}
    for _, r in pg.iterrows():
        tag = str(r["tag"])
        sub = pn[pn["tag"] == tag]
        out[tag] = GonadScore(expert=sub["expert_n"].to_numpy(float), ours=sub["ours_n"].to_numpy(float),
                              n_matched=int(r["n_matched"]), n_ours=int(r["n_ours"]), n_expert=int(r["n_expert"]),
                              source=str(run_dir).replace("\\", "/"))
        if truth_masks:
            if ptn is not None:          # a gonad without truth-mask rows has none to pool, as in bench_detector.py
                tm = ptn[ptn["tag"].astype(str) == tag]
                out[tag].tm_expert, out[tag].tm_ours = tm["expert_n"].to_numpy(float), tm["ours_n"].to_numpy(float)
            if has_tm_spots:             # a folder written before the section 5.2 columns has none
                out[tag].tm_spots = _tm_spot_counts(r)
    return params, run_dir.name, out


def _merge_truth_masks(o: GonadScore, s: GonadScore, what: str) -> None:
    """A gonad scored twice under one parameter set (``index_runs`` with ``truth_masks``): its truth-mask
    counts must agree where both runs have them; where only the second has them, they are taken from it."""
    for attrs in (("tm_expert", "tm_ours"), ("tm_spots",)):
        mine, theirs = [getattr(o, a) for a in attrs], [getattr(s, a) for a in attrs]
        if any(v is None for v in theirs):
            continue
        if any(v is None for v in mine):
            for a, v in zip(attrs, theirs):
                setattr(o, a, v)
        elif not all(np.array_equal(np.asarray(a), np.asarray(b)) for a, b in zip(mine, theirs)):
            raise ValueError(f"{what}: {o.source} and {s.source} disagree on the truth-mask counts")


def index_runs(roots, *, backend: str | None = None, truth_masks: bool = False) -> RunIndex:
    """Index the run folders under ``roots`` (each a folder of ``bench_tune.py`` set folders or one
    ``bench_detector.py`` output): every folder holding pooled.json, per_gonad.csv and per_nucleus.csv
    adds its per-gonad scores under its parameter set. ``backend`` keeps only runs of that backend. A
    gonad scored twice under the same parameter set must give the same counts (ValueError otherwise).
    ``truth_masks`` also reads what ``primary_objective_on`` pools (per_truth_nucleus.csv and the tm_spots_
    columns, ``GonadScore``), with the same rule for a gonad scored twice; without it nothing more is read
    than the D3 objective needs."""
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
            params, name, per_tag = _read_run_dir(d, truth_masks)
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
                    if truth_masks:
                        _merge_truth_masks(o, s, f"{tag} under {name}")
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


def primary_objective_on(idx: RunIndex, sha: str, tags) -> tuple[float, float]:
    """The primary metric of one parameter set pooled over ``tags``, as ``bench_detector.py`` pools it (the
    point estimates of its ``per_nucleus_truth_masks`` and ``spots_truth_masks`` blocks over exactly these
    gonads): the per-nucleus CCC over the concatenated truth nuclei (protocol section 5.1, Class A+B truth,
    end-to-end counts), rounded to ``tuning.CCC_DECIMALS``, then the F1 of the summed matching counts inside
    R (section 5.2). Ranked with ``tuning.objective``, so a missing or NaN value is -inf and ranks last. A
    gonad without a score for the set makes the value missing, and so does a score without the truth-mask
    counts (an index built without ``truth_masks``, a run folder without them). A gonad with no truth nuclei,
    or with no spot in R, adds nothing, as in ``bench_detector.py``."""
    s = [idx.scores.get(sha, {}).get(t) for t in tags]
    pooled = {}
    if s and all(g is not None for g in s):
        if all(g.tm_expert is not None and g.tm_ours is not None for g in s):
            ccc = lin_ccc(np.concatenate([g.tm_expert for g in s]), np.concatenate([g.tm_ours for g in s]))
            pooled[PRIMARY_OBJECTIVE[0]] = {PRIMARY_OBJECTIVE[1]: {"point": ccc}}
        if all(g.tm_spots is not None for g in s):
            det = pooled_detection([g.tm_spots for g in s])
            pooled[PRIMARY_TIE_BREAK[0]] = {PRIMARY_TIE_BREAK[1]: {"point": det["f1"]}}
    return bt.objective(pooled, keys=(PRIMARY_OBJECTIVE, PRIMARY_TIE_BREAK))


def rank_on(idx: RunIndex, configs: list[dict], tags) -> list[dict]:
    """``configs`` sorted best first by the objective pooled over ``tags`` (stable: ties keep their order),
    each with its ``sha`` and ``objective``. Every set needs a score on every gonad (``missing_scores``
    lists the gaps; KeyError otherwise)."""
    items = []
    for c in configs:
        sha = bd.params_sha256(c["params"])
        items.append({**c, "sha": sha, "objective": objective_on(idx, sha, tags)})
    return sorted(items, key=lambda it: it["objective"], reverse=True)


def missing_scores(idx: RunIndex, configs: list[dict], tags, *, truth_masks: bool = False) -> list[dict]:
    """The (set, gonad) pairs of ``configs`` x ``tags`` with no score in the index. With ``truth_masks`` (the
    primary objective's rule), a score without the truth-mask counts is missing too
    (``GonadScore.has_truth_masks``: a run folder without per_truth_nucleus.csv or the tm_spots_ columns)."""
    out = []
    for c in configs:
        sha = bd.params_sha256(c["params"])
        have = idx.scores.get(sha, {})
        lack = [t for t in tags if t not in have or (truth_masks and not have[t].has_truth_masks())]
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


# ----------------------------------------------------------- the pre-registered SpotMAX grid
def prereg_spotmax_grid() -> list[dict]:
    """The 60 parameter sets of protocol section 8 item 2 as ``{name, params}``: thresholding method
    outermost, in the protocol's order (otsu, li, triangle, yen), then ``effect_size_min`` and then
    ``z_merge_gap_um``, both ascending. ``params`` holds those three keys only, so every other spots key
    stays at the golden run's deployed value (``spot_radius_um`` 0.3, ``gauss_sigma_um`` 0.08,
    ``merge_z_columns`` on, ``z_merge_valley_frac`` 0.8, the Glass effect size, ``max_spot_candidates``
    30000) and pg_tri_e3.0_z0.8 is the deployed recipe. Names: pg_<otsu|li|tri|yen>_e<effect size>_z<z gap>,
    one decimal each."""
    out = []
    for m in PREREG_SPOTMAX_METHODS:
        for e in PREREG_SPOTMAX_EFFECT_SIZES:
            for z in PREREG_SPOTMAX_Z_GAPS:
                out.append({"name": f"pg_{SPOTMAX_ABBR[m]}_e{float(e):.1f}_z{float(z):.1f}",
                            "params": {"thresholding_method": m, "effect_size_min": float(e),
                                       "z_merge_gap_um": float(z)}})
    return out


def designated_split(folds: list[dict]) -> dict:
    """The designated split behind ``folds`` (``logo_folds``): its train and val gonads, that is every
    held-out gonad of the folds by its role (usable spot truth and a bench cache, never a test gonad), in
    tag order. The grid set a competitor takes to test is chosen on these val gonads (protocol section 6.1)."""
    return {"train": [f["held_out"] for f in folds if f["role"] == "train"],
            "val": [f["held_out"] for f in folds if f["role"] == "val"]}


def select_grid(idx: RunIndex, configs: list[dict], tags) -> dict:
    """The grid rule of protocol sections 6.1 and 6.4 on the selection gonads ``tags``: every set of
    ``configs`` ranked by ``primary_objective_on`` over ``tags``, the higher rounded CCC first, then the
    higher spot F1, then the order of ``configs`` (a stable sort). Returns ``selected`` (the best set as
    ``{name, params}``, or None while a pair is missing), ``ranking`` (every set with its ``objective``, best
    first; a set with a pair missing ranks last) and ``missing`` (the pairs still to run, by set). A pair is
    missing when it has no score or its score lacks the truth-mask counts (``missing_scores(...,
    truth_masks=True)``). A run folder without per_truth_nucleus.csv or the tm_spots_ columns therefore holds
    up the selection; it does not let another set be selected. A NaN CCC computed from the counts still ranks
    last. A test gonad in ``tags`` is refused (``detector.guard_test``). ValueError without selection gonads,
    and when nothing is missing but no set has a finite CCC."""
    tags = list(tags)
    bd.guard_test(tags)
    if not configs:
        raise ValueError("select_grid needs at least one parameter set")
    if not tags:
        raise ValueError("select_grid needs at least one selection gonad (a fold or split without val gonads?)")
    ranking = sorted(({"name": c["name"], "params": c["params"],
                       "objective": primary_objective_on(idx, bd.params_sha256(c["params"]), tags)}
                      for c in configs), key=lambda it: it["objective"], reverse=True)
    missing = missing_scores(idx, configs, tags, truth_masks=True)
    selected = None
    if not missing:
        if not math.isfinite(ranking[0]["objective"][0]):
            raise ValueError(f"no parameter set has a finite truth-mask CCC on {tags}")
        selected = {"name": ranking[0]["name"], "params": ranking[0]["params"]}
    return {"selected": selected, "ranking": ranking, "missing": missing}


def replay_grid(idx: RunIndex, configs: list[dict], train: list[str], val: list[str]) -> dict:
    """The grid rule on one fold, or on ``designated_split``: ``select_grid`` over the fold's selection
    gonads ``val``. ``train`` takes no part in the rule (every set has a val score and the best on val is
    selected, protocol section 6.4); it is only checked, with ``val``, to hold no test gonad. Returns what
    ``replay_spotmax_search`` returns and ``bench_logo.py`` reads: ``complete``, ``missing`` (with ``stage``
    "grid" while incomplete), ``n_sets``, ``finalists`` (every set competes on val, so the whole grid in its
    order) and, when complete, ``selected``, ``selected_params`` and ``val_objective`` (every set's
    objective on ``val``, best first)."""
    bd.guard_test(list(train) + list(val))
    sel = select_grid(idx, configs, val)
    out = {"n_sets": len(configs), "finalists": [c["name"] for c in configs]}
    if sel["missing"]:
        return {**out, "complete": False, "missing": sel["missing"], "stage": "grid"}
    return {**out, "complete": True, "missing": [], "selected": sel["selected"]["name"],
            "selected_params": sel["selected"]["params"],
            "val_objective": {x["name"]: list(x["objective"]) for x in sel["ranking"]}}


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
