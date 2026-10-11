"""Agreement statistics and the cluster bootstrap by gonad.

The unit of replication in every benchmark is the gonad: nuclei and spots inside one gonad share its
staining, its focus and its expert, so a confidence interval that resamples nuclei would be far too
narrow. ``cluster_bootstrap`` resamples whole gonads with replacement and recomputes the pooled
statistic each time (1000 resamples, seed 0 by default); the 2.5 and 97.5 percentiles are the
interval. Everything is deterministic for a given seed.

Lin's concordance correlation coefficient (CCC) measures agreement with the identity line, so unlike
Pearson's r it drops when the counts are correlated but biased. ``bias`` is mean(ours minus expert):
positive when we over-count.
"""
from __future__ import annotations

import math
from collections.abc import Callable, Sequence

import numpy as np


def _clean_pair(x, y):
    x = np.asarray(x, float).ravel()
    y = np.asarray(y, float).ravel()
    if len(x) != len(y):
        raise ValueError("paired arrays must have the same length")
    ok = np.isfinite(x) & np.isfinite(y)
    return x[ok], y[ok]


def pearson(x, y) -> float:
    """Pearson's r; NaN with fewer than two pairs or a constant series."""
    x, y = _clean_pair(x, y)
    if len(x) < 2 or x.std() == 0 or y.std() == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def lin_ccc(x, y) -> float:
    """Lin's concordance correlation coefficient with population moments (the same convention as
    germquant.validate.compare.agreement_stats); NaN with fewer than two pairs or when both series are
    the same constant."""
    x, y = _clean_pair(x, y)
    if len(x) < 2:
        return float("nan")
    mx, my = x.mean(), y.mean()
    cov = float(np.mean((x - mx) * (y - my)))
    denom = float(x.var() + y.var() + (mx - my) ** 2)
    if denom == 0:
        return float("nan")
    return float(2 * cov / denom)


def bias(ours, expert) -> float:
    """Mean of (ours minus expert); positive when we over-count."""
    o, e = _clean_pair(ours, expert)
    return float(np.mean(o - e)) if len(o) else float("nan")


def mean_abs_error(ours, expert) -> float:
    """Mean of |ours minus expert| (protocol section 5.1); NaN without pairs."""
    o, e = _clean_pair(ours, expert)
    return float(np.mean(np.abs(o - e))) if len(o) else float("nan")


def agreement_metrics(expert, ours) -> dict:
    """n, means, Pearson's r, Lin's CCC, bias and mean absolute error for paired per-nucleus counts."""
    e, o = _clean_pair(expert, ours)
    return {
        "n": int(len(e)),
        "mean_expert": float(e.mean()) if len(e) else float("nan"),
        "mean_ours": float(o.mean()) if len(o) else float("nan"),
        "pearson_r": pearson(e, o),
        "ccc": lin_ccc(e, o),
        "bias": bias(o, e),
        "mae": mean_abs_error(o, e),
    }


def precision_recall_f1(n_matched, n_pred, n_truth) -> dict:
    """Detection scores from counts. Precision is NaN without predictions, recall is NaN without truth
    objects, and F1 is 0 when both are defined but nothing matched."""
    n_matched, n_pred, n_truth = int(n_matched), int(n_pred), int(n_truth)
    p = n_matched / n_pred if n_pred else float("nan")
    r = n_matched / n_truth if n_truth else float("nan")
    if math.isnan(p) or math.isnan(r):
        f = float("nan")
    elif p + r == 0:
        f = 0.0
    else:
        f = 2 * p * r / (p + r)
    return {"n_matched": n_matched, "n_pred": n_pred, "n_truth": n_truth,
            "precision": p, "recall": r, "f1": f}


def cluster_bootstrap(groups: Sequence, statistic: Callable, *, n_boot: int = 1000, seed: int = 0,
                      alpha: float = 0.05) -> dict:
    """Cluster bootstrap over ``groups`` (one payload per gonad).

    ``statistic(list_of_payloads)`` returns a float or a dict of floats computed on the pooled payloads.
    Gonads are resampled with replacement (the same number as observed) ``n_boot`` times with
    ``numpy.random.default_rng(seed)``. Returns ``{"n_boot", "n_groups", "seed", "alpha", "stats"}``
    where ``stats`` maps each statistic name (``"value"`` for a scalar statistic) to
    ``{"point", "lo", "hi", "n_nan"}``; ``lo``/``hi`` are the alpha/2 and 1 - alpha/2 percentiles of
    the resampled values, NaN resamples left out and counted in ``n_nan``.
    """
    groups = list(groups)
    n = len(groups)
    point = statistic(groups)
    scalar = not isinstance(point, dict)
    point_d = {"value": point} if scalar else dict(point)
    names = list(point_d)
    draws_by_name: dict[str, list[float]] = {k: [] for k in names}
    if n:
        rng = np.random.default_rng(seed)
        draws = rng.integers(0, n, size=(int(n_boot), n))
        for d in draws:
            v = statistic([groups[i] for i in d])
            v = {"value": v} if scalar else v
            for k in names:
                draws_by_name[k].append(float(v.get(k, float("nan"))))
    stats = {}
    for k in names:
        arr = np.asarray(draws_by_name[k], float)
        finite = arr[np.isfinite(arr)]
        lo = float(np.percentile(finite, 100 * alpha / 2)) if finite.size else float("nan")
        hi = float(np.percentile(finite, 100 * (1 - alpha / 2))) if finite.size else float("nan")
        stats[k] = {"point": float(point_d[k]), "lo": lo, "hi": hi,
                    "n_nan": int(arr.size - finite.size)}
    return {"n_boot": int(n_boot), "n_groups": n, "seed": int(seed), "alpha": float(alpha),
            "stats": stats}


def pooled_agreement(payloads) -> dict:
    """Statistic for ``cluster_bootstrap``: each payload is (expert_counts, ours_counts); the pooled
    CCC, Pearson's r, bias and mean absolute error over every nucleus of the given gonads."""
    if not payloads:
        return {"ccc": float("nan"), "pearson_r": float("nan"), "bias": float("nan"), "mae": float("nan")}
    e = np.concatenate([np.asarray(p[0], float).ravel() for p in payloads])
    o = np.concatenate([np.asarray(p[1], float).ravel() for p in payloads])
    return {"ccc": lin_ccc(e, o), "pearson_r": pearson(e, o), "bias": bias(o, e), "mae": mean_abs_error(o, e)}


def pooled_detection(payloads) -> dict:
    """Statistic for ``cluster_bootstrap``: each payload is (n_matched, n_pred, n_truth); the pooled
    precision, recall and F1 from the summed counts."""
    if not payloads:
        return {"precision": float("nan"), "recall": float("nan"), "f1": float("nan")}
    m = sum(int(p[0]) for p in payloads)
    npred = sum(int(p[1]) for p in payloads)
    ntruth = sum(int(p[2]) for p in payloads)
    s = precision_recall_f1(m, npred, ntruth)
    return {"precision": s["precision"], "recall": s["recall"], "f1": s["f1"]}


def pooled_mean(payloads) -> float:
    """Statistic for ``cluster_bootstrap``: the mean over every value of every payload."""
    vals = [np.asarray(p, float).ravel() for p in payloads]
    vals = [v for v in vals if v.size]
    if not vals:
        return float("nan")
    return float(np.nanmean(np.concatenate(vals)))


def paired_cluster_bootstrap(groups_a: Sequence, groups_b: Sequence, statistic: Callable, *,
                             n_boot: int = 1000, seed: int = 0, alpha: float = 0.05) -> dict:
    """The paired comparison of protocol section 6.3: two methods scored on the same gonads (one payload
    per gonad and method, in the same gonad order). Every draw resamples the gonads once, with
    replacement, and computes both methods' pooled statistic on that same resample; the interval of the
    difference (a minus b) is the alpha/2 and 1 - alpha/2 percentiles of the per-draw differences. The
    draws are those of ``cluster_bootstrap`` for the same seed and number of gonads, so the per-method
    intervals of that function come from the same resamples.

    Returns ``{"n_boot", "n_groups", "seed", "alpha", "stats"}``; ``stats`` maps each statistic name to
    ``{"a", "b", "diff": {"point", "lo", "hi", "n_nan"}, "excludes_zero"}`` where ``a`` and ``b`` are the
    point estimates and ``excludes_zero`` says whether the interval of the difference lies wholly above
    or wholly below zero."""
    groups_a, groups_b = list(groups_a), list(groups_b)
    if len(groups_a) != len(groups_b):
        raise ValueError("a paired comparison needs both methods scored on the same gonads, in the same order")
    n = len(groups_a)
    pa, pb = statistic(groups_a), statistic(groups_b)
    scalar = not isinstance(pa, dict)
    pa_d = {"value": pa} if scalar else dict(pa)
    pb_d = {"value": pb} if scalar else dict(pb)
    names = list(pa_d)
    diffs: dict[str, list[float]] = {k: [] for k in names}
    if n:
        rng = np.random.default_rng(seed)
        for d in rng.integers(0, n, size=(int(n_boot), n)):
            va = statistic([groups_a[i] for i in d])
            vb = statistic([groups_b[i] for i in d])
            va = {"value": va} if scalar else va
            vb = {"value": vb} if scalar else vb
            for k in names:
                diffs[k].append(float(va.get(k, float("nan"))) - float(vb.get(k, float("nan"))))
    stats = {}
    for k in names:
        arr = np.asarray(diffs[k], float)
        finite = arr[np.isfinite(arr)]
        lo = float(np.percentile(finite, 100 * alpha / 2)) if finite.size else float("nan")
        hi = float(np.percentile(finite, 100 * (1 - alpha / 2))) if finite.size else float("nan")
        a_pt, b_pt = float(pa_d[k]), float(pb_d.get(k, float("nan")))
        stats[k] = {"a": a_pt, "b": b_pt,
                    "diff": {"point": a_pt - b_pt, "lo": lo, "hi": hi, "n_nan": int(arr.size - finite.size)},
                    "excludes_zero": bool(finite.size and (lo > 0 or hi < 0))}
    return {"n_boot": int(n_boot), "n_groups": n, "seed": int(seed), "alpha": float(alpha), "stats": stats}


def paired_comparison(a, b, *, n_boot: int = 1000, seed: int = 0, label_a: str = "a",
                      label_b: str = "b") -> dict:
    """Two detectors scored on the same gonads: per-gonad differences of a metric (a minus b), the
    mean difference with a cluster-bootstrap interval, how many gonads favour each side, and the
    Wilcoxon signed-rank p-value when it can be computed (at least 5 non-zero differences)."""
    a = np.asarray(a, float).ravel()
    b = np.asarray(b, float).ravel()
    if len(a) != len(b):
        raise ValueError("the two detectors must be scored on the same gonads, in the same order")
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    diff = a - b
    boot = cluster_bootstrap(list(diff), pooled_mean, n_boot=n_boot, seed=seed)["stats"]["value"]
    p_value = float("nan")
    nz = diff[diff != 0]
    if len(nz) >= 5:
        try:
            from scipy.stats import wilcoxon

            p_value = float(wilcoxon(a, b).pvalue)
        except Exception:  # noqa: BLE001 - all-equal or degenerate input, leave NaN
            p_value = float("nan")
    return {
        "label_a": label_a, "label_b": label_b, "n_gonads": int(len(diff)),
        "mean_a": float(a.mean()) if len(a) else float("nan"),
        "mean_b": float(b.mean()) if len(b) else float("nan"),
        "mean_diff": boot["point"], "diff_lo": boot["lo"], "diff_hi": boot["hi"],
        "n_a_better": int((diff > 0).sum()), "n_b_better": int((diff < 0).sum()),
        "n_ties": int((diff == 0).sum()), "wilcoxon_p": p_value,
        "per_gonad_diff": [float(d) for d in diff],
    }
