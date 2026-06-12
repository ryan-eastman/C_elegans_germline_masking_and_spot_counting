"""Validate pipeline output against hand-scored ground truth.

Two kinds of comparison:
  * agreement_stats / compare_table — for per-object counts/lengths (SC fragments,
    RAD-51 foci/nucleus, TZ length): Lin's CCC, Pearson/Spearman, bias, MAE/RMSE, and
    Bland-Altman limits of agreement (the plot reviewers ask for).
  * segmentation_metrics — instance F1 by IoU matching + mean matched IoU/Dice
    (Metrics Reloaded style) for nucleus segmentation vs a labelled ground-truth mask.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


# ---------------------------------------------------------------- counts / lengths
def agreement_stats(pred: np.ndarray, truth: np.ndarray) -> dict:
    pred = np.asarray(pred, float)
    truth = np.asarray(truth, float)
    ok = np.isfinite(pred) & np.isfinite(truth)
    pred, truth = pred[ok], truth[ok]
    n = len(pred)
    out = {"n": n}
    if n == 0:
        return out
    diff = pred - truth
    out["bias_mean_diff"] = float(diff.mean())
    out["mae"] = float(np.abs(diff).mean())
    out["rmse"] = float(np.sqrt((diff**2).mean()))
    if n >= 2:
        sd = float(diff.std(ddof=1))
        out["sd_diff"] = sd
        out["loa_lower"] = float(diff.mean() - 1.96 * sd)
        out["loa_upper"] = float(diff.mean() + 1.96 * sd)
        # Lin's concordance correlation coefficient
        cov = float(np.cov(pred, truth, ddof=0)[0, 1])
        vp, vt = float(pred.var()), float(truth.var())
        denom = vp + vt + (pred.mean() - truth.mean()) ** 2
        out["ccc"] = (2 * cov / denom) if denom > 0 else float("nan")
        from scipy.stats import pearsonr, spearmanr

        try:
            out["pearson_r"] = float(pearsonr(pred, truth)[0])
            out["spearman_rho"] = float(spearmanr(pred, truth)[0])
        except Exception:
            out["pearson_r"] = out["spearman_rho"] = float("nan")
    return out


def compare_table(
    pred_df: pd.DataFrame,
    truth_df: pd.DataFrame,
    key_cols: list[str],
    pred_col: str,
    truth_col: str,
) -> tuple[pd.DataFrame, dict]:
    """Join pred & truth on key_cols, return (paired_df, agreement_stats)."""
    merged = pred_df.merge(truth_df, on=key_cols, how="inner", suffixes=("_pred", "_truth"))
    pc = pred_col if pred_col in merged else f"{pred_col}_pred"
    tc = truth_col if truth_col in merged else f"{truth_col}_truth"
    if pc not in merged or tc not in merged:
        raise KeyError(f"columns not found after join: {pred_col!r}/{truth_col!r}; have {list(merged.columns)}")
    paired = merged[key_cols + [pc, tc]].rename(columns={pc: "pred", tc: "truth"})
    return paired, agreement_stats(paired["pred"].to_numpy(), paired["truth"].to_numpy())


def bland_altman_plot(pred, truth, out_path, *, title="", units=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pred = np.asarray(pred, float)
    truth = np.asarray(truth, float)
    mean = (pred + truth) / 2
    diff = pred - truth
    bias = diff.mean()
    sd = diff.std(ddof=1) if len(diff) > 1 else 0.0

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(mean, diff, s=14, alpha=0.6)
    ax.axhline(bias, color="C3", label=f"bias {bias:.2f}")
    ax.axhline(bias + 1.96 * sd, color="C0", ls="--", label="±1.96 SD")
    ax.axhline(bias - 1.96 * sd, color="C0", ls="--")
    ax.set_xlabel(f"mean of pred & truth {units}")
    ax.set_ylabel(f"pred − truth {units}")
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------- segmentation
def segmentation_metrics(pred_labels: np.ndarray, gt_labels: np.ndarray, iou_threshold: float = 0.5) -> dict:
    """Instance detection F1 by greedy IoU matching + mean matched IoU/Dice."""
    iou, pred_ids, gt_ids = _iou_matrix(pred_labels, gt_labels)
    if iou.size == 0:
        return {"n_pred": int(len(pred_ids)), "n_gt": int(len(gt_ids)), "tp": 0,
                "fp": int(len(pred_ids)), "fn": int(len(gt_ids)), "f1": 0.0,
                "precision": 0.0, "recall": 0.0, "mean_iou": float("nan"),
                "mean_dice": float("nan"), "iou_threshold": iou_threshold}

    # candidate pairs above threshold, matched greedily by descending IoU
    pi, gi = np.where(iou >= iou_threshold)
    order = np.argsort(-iou[pi, gi])
    used_p, used_g, matched_ious = set(), set(), []
    for k in order:
        p, g = int(pi[k]), int(gi[k])
        if p in used_p or g in used_g:
            continue
        used_p.add(p)
        used_g.add(g)
        matched_ious.append(float(iou[p, g]))

    tp = len(matched_ious)
    fp = len(pred_ids) - tp
    fn = len(gt_ids) - tp
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    mean_iou = float(np.mean(matched_ious)) if matched_ious else float("nan")
    mean_dice = float(np.mean([2 * i / (1 + i) for i in matched_ious])) if matched_ious else float("nan")
    return {"n_pred": int(len(pred_ids)), "n_gt": int(len(gt_ids)), "tp": tp, "fp": fp, "fn": fn,
            "precision": precision, "recall": recall, "f1": f1, "mean_iou": mean_iou,
            "mean_dice": mean_dice, "iou_threshold": iou_threshold}


def _iou_matrix(pred: np.ndarray, gt: np.ndarray):
    """IoU between every (pred, gt) instance pair, excluding background (label 0)."""
    p = pred.ravel()
    g = gt.ravel()
    fg = (p > 0) | (g > 0)
    p, g = p[fg], g[fg]
    pu = np.unique(p[p > 0])
    gu = np.unique(g[g > 0])
    if len(pu) == 0 or len(gu) == 0:
        return np.zeros((0, 0)), pu, gu
    pidx = {v: i for i, v in enumerate(pu)}
    gidx = {v: i for i, v in enumerate(gu)}
    inter = np.zeros((len(pu), len(gu)), dtype=np.int64)
    both = (p > 0) & (g > 0)
    pi = np.array([pidx[v] for v in p[both]])
    gi = np.array([gidx[v] for v in g[both]])
    np.add.at(inter, (pi, gi), 1)
    p_area = np.array([int((pred == v).sum()) for v in pu])
    g_area = np.array([int((gt == v).sum()) for v in gu])
    union = p_area[:, None] + g_area[None, :] - inter
    iou = inter / np.maximum(union, 1)
    return iou, pu, gu
