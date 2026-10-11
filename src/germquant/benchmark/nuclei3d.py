"""3D nucleus instance benchmark: our label image against the expert's partial label image.

The expert traced a subset of the nuclei (90 to 150 per gonad, inside a hand-drawn ROI box), so the
truth is partial and a plain F1 would punish every real nucleus the expert did not trace. The rule
here: a prediction that overlaps no truth object AND lies outside the truth coverage region is
ignored, not counted as a false positive. The coverage region is the ROI box from the truth params
JSON when it is known, else the bounding box of all truth labels dilated by 2 um. A prediction "lies
outside" when its centroid is outside the region. The counts of ignored predictions are reported so a
reader can judge the rule. Predictions inside the region that overlap nothing are false positives:
with a threshold-curated truth that is a lower bound on precision, and recall is the number to trust.

Overlaps are counted from a joint histogram of (prediction label, truth label) pairs, plane by plane,
so no dense n_pred by n_truth matrix is formed and a 2600 x 2600 x 71 image is never duplicated.
Matching is the Hungarian algorithm on the sparse IoU graph (``matching.assign_sparse``) with pairs
below the IoU threshold excluded; it maximises the number of matches first and the total IoU second.

Metrics: F1 at IoU 0.5 (plus precision and recall), AP averaged over IoU 0.5 to 0.95 in steps of
0.05 where AP at one threshold is TP / (TP + FP + FN) (the DSB-2018 / StarDist definition, there are
no confidence scores to rank by), splits (a truth object that holds the majority of the volume of two
or more predictions), merges (a prediction that holds the majority of the volume of two or more truth
objects) and the z-extent error of matched pairs (predicted z span minus truth z span, in um, with a
span of (zmax - zmin + 1) planes). Everything is deterministic.

The pre-registered rule (``coverage="protocol"``, docs/BENCHMARK_PROTOCOL.md section 5.3). The rules
above came first and stay the default of this module and of scripts/benchmark_nuclei3d.py, so earlier
outputs keep reproducing. The protocol, registered later, defines the scope differently, and a score
reported under the protocol uses this mode:

* truth fragments below ``min_truth_voxels`` (1000) voxels are removed from the truth first (section
  4.2); they are neither matched nor counted;
* the truth coverage region C is the union of the remaining truth masks dilated by 1.0 um in xy and one
  plane in z (an exact cylindrical dilation);
* a prediction is in scope when at least half of its voxels lie inside C; predictions out of scope are
  ignored (``n_pred_ignored_outside_coverage``), in-scope predictions without a match are false
  positives;
* splits: truth nuclei overlapped by two or more in-scope predictions that each cover at least 20
  percent of the truth volume; merges: in-scope predictions that cover at least 20 percent of the
  volume of two or more truth nuclei.

The ROI box of the Imaris wizard is never used by this mode (section 4.2: it is not a record of the
traced extent).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .matching import DEFAULT_MAX_DENSE_CELLS, Assignment, assign_sparse
from .stats import cluster_bootstrap, pooled_detection, pooled_mean, precision_recall_f1

IOU_THRESHOLDS = tuple(round(0.5 + 0.05 * k, 2) for k in range(10))
DEFAULT_DILATE_UM = 2.0
DEFAULT_MAJORITY_FRAC = 0.5
# docs/BENCHMARK_PROTOCOL.md, sections 4.2 and 5.3
PROTOCOL_RULES = {"min_truth_voxels": 1000, "coverage_xy_um": 1.0, "coverage_z_planes": 1,
                  "scope_frac": 0.5, "split_frac": 0.2}
COVERAGE_CHOICES = ("auto", "roi", "bbox", "none", "protocol")


def _tag(t: float) -> str:
    return f"iou{int(round(t * 100)):02d}"


# ------------------------------------------------------------------------------ joint histogram
@dataclass
class OverlapTable:
    """Per-object sizes, z extents and centroids plus the sparse overlap counts of every
    (prediction, truth) pair that shares at least one voxel."""

    shape: tuple
    pred_ids: np.ndarray
    pred_sizes: np.ndarray
    pred_zmin: np.ndarray
    pred_zmax: np.ndarray
    pred_centroid_zyx: np.ndarray
    truth_ids: np.ndarray
    truth_sizes: np.ndarray
    truth_zmin: np.ndarray
    truth_zmax: np.ndarray
    truth_bbox_min_zyx: np.ndarray | None
    truth_bbox_max_zyx: np.ndarray | None
    pair_pred: np.ndarray      # index into pred_ids
    pair_truth: np.ndarray     # index into truth_ids
    pair_inter: np.ndarray     # shared voxels

    @property
    def n_pred(self) -> int:
        return int(len(self.pred_ids))

    @property
    def n_truth(self) -> int:
        return int(len(self.truth_ids))

    def iou(self) -> np.ndarray:
        union = (self.pred_sizes[self.pair_pred] + self.truth_sizes[self.pair_truth]
                 - self.pair_inter)
        return self.pair_inter / np.maximum(union, 1)

    def pred_has_overlap(self) -> np.ndarray:
        has = np.zeros(self.n_pred, bool)
        has[self.pair_pred] = True
        return has


def _as_volume(a):
    if hasattr(a, "ndim") and a.ndim == 2:
        a = np.asarray(a)[None]
    if a.ndim != 3:
        raise ValueError(f"expected a ZYX label image, got {a.ndim} dimensions")
    return a


def _max_label(vol) -> int:
    m = 0
    for z in range(vol.shape[0]):
        plane = np.asarray(vol[z])
        if plane.size:
            m = max(m, int(plane.max()))
    return m


def overlap_table(pred, truth) -> OverlapTable:
    """Plane-by-plane joint histogram of a predicted and a truth label image of the same shape.
    ``pred`` and ``truth`` may be numpy arrays or any object with ``shape`` and plane indexing."""
    pred = _as_volume(pred)
    truth = _as_volume(truth)
    if tuple(pred.shape) != tuple(truth.shape):
        raise ValueError(f"label images differ in shape: {pred.shape} vs {truth.shape}")
    Z, Y, X = pred.shape
    pmax, tmax = _max_label(pred), _max_label(truth)
    yy, xx = np.meshgrid(np.arange(Y, dtype=float), np.arange(X, dtype=float), indexing="ij")
    yy, xx = yy.ravel(), xx.ravel()

    p_size = np.zeros(pmax + 1, np.int64)
    p_sum = np.zeros((pmax + 1, 3), float)
    p_zmin = np.full(pmax + 1, Z, np.int64)
    p_zmax = np.full(pmax + 1, -1, np.int64)
    t_size = np.zeros(tmax + 1, np.int64)
    t_zmin = np.full(tmax + 1, Z, np.int64)
    t_zmax = np.full(tmax + 1, -1, np.int64)
    bbox_min = np.array([Z, Y, X], np.int64)
    bbox_max = np.array([-1, -1, -1], np.int64)
    keys, counts = [], []
    stride = tmax + 1

    for z in range(Z):
        p = np.asarray(pred[z]).astype(np.int64, copy=False)
        t = np.asarray(truth[z]).astype(np.int64, copy=False)
        if p.min() < 0 or t.min() < 0:
            raise ValueError("label images must not hold negative labels")
        pf = p.ravel()
        bc = np.bincount(pf, minlength=pmax + 1)
        p_size += bc
        present = np.flatnonzero(bc[1:]) + 1
        p_zmin[present] = np.minimum(p_zmin[present], z)
        p_zmax[present] = z
        p_sum[:, 0] += bc * float(z)
        p_sum[:, 1] += np.bincount(pf, weights=yy, minlength=pmax + 1)
        p_sum[:, 2] += np.bincount(pf, weights=xx, minlength=pmax + 1)

        tf = t.ravel()
        tb = np.bincount(tf, minlength=tmax + 1)
        t_size += tb
        tpresent = np.flatnonzero(tb[1:]) + 1
        t_zmin[tpresent] = np.minimum(t_zmin[tpresent], z)
        t_zmax[tpresent] = z
        tm = t > 0
        if tpresent.size:
            rows = np.flatnonzero(tm.any(axis=1))
            cols = np.flatnonzero(tm.any(axis=0))
            bbox_min = np.minimum(bbox_min, [z, rows.min(), cols.min()])
            bbox_max = np.maximum(bbox_max, [z, rows.max(), cols.max()])
            both = tm & (p > 0)
            if both.any():
                k = p[both] * stride + t[both]
                u, c = np.unique(k, return_counts=True)
                keys.append(u)
                counts.append(c)

    pred_ids = np.flatnonzero(p_size[1:]) + 1
    truth_ids = np.flatnonzero(t_size[1:]) + 1
    if keys:
        allk = np.concatenate(keys)
        allc = np.concatenate(counts)
        u, inv = np.unique(allk, return_inverse=True)
        inter = np.bincount(inv, weights=allc).astype(np.int64)
        pl, tl = u // stride, u % stride
        pair_pred = np.searchsorted(pred_ids, pl)
        pair_truth = np.searchsorted(truth_ids, tl)
    else:
        pair_pred = pair_truth = inter = np.zeros(0, np.int64)
    centroid = p_sum[pred_ids] / np.maximum(p_size[pred_ids], 1)[:, None]
    has_truth = truth_ids.size > 0
    return OverlapTable(
        shape=(Z, Y, X), pred_ids=pred_ids, pred_sizes=p_size[pred_ids], pred_zmin=p_zmin[pred_ids],
        pred_zmax=p_zmax[pred_ids], pred_centroid_zyx=centroid, truth_ids=truth_ids,
        truth_sizes=t_size[truth_ids], truth_zmin=t_zmin[truth_ids], truth_zmax=t_zmax[truth_ids],
        truth_bbox_min_zyx=bbox_min if has_truth else None,
        truth_bbox_max_zyx=bbox_max if has_truth else None,
        pair_pred=pair_pred, pair_truth=pair_truth, pair_inter=inter)


# ------------------------------------------------------------------------------------ matching
def match_instances(ot: OverlapTable, iou_threshold: float, *,
                    max_dense_cells: int = DEFAULT_MAX_DENSE_CELLS, pred_ok: np.ndarray | None = None) -> Assignment:
    """One-to-one matching of predictions (rows) to truth objects (cols) among pairs with
    IoU >= threshold, maximising the number of matches, then the summed IoU. ``pred_ok`` (a boolean per
    prediction) restricts the matching to those predictions (the protocol's in-scope predictions)."""
    iou = ot.iou()
    ok = iou >= iou_threshold - 1e-12
    if pred_ok is not None:
        ok &= np.asarray(pred_ok, bool)[ot.pair_pred]
    if not ok.any():
        z = np.zeros(0, np.int64)
        return Assignment(z, z.copy(), np.zeros(0, float), 0, 0, 0, "none", [])
    return assign_sparse(ot.pair_pred[ok], ot.pair_truth[ok], 1.0 - iou[ok], ot.n_pred, ot.n_truth,
                         max_dense_cells=max_dense_cells, max_candidates=None)


def coverage_box(ot: OverlapTable, voxel_zyx, *, dilate_um: float = DEFAULT_DILATE_UM,
                 roi_box: dict | None = None, rule: str = "auto"):
    """The truth coverage region as inclusive voxel bounds (min_zyx, max_zyx, rule_name), or
    (None, None, "none"). ``rule``: "auto" (the ROI box when given, else the dilated truth bounding
    box), "roi" (ROI box only), "bbox" (dilated bounding box only), "none" (no region: nothing is
    ignored, every unmatched prediction is a false positive)."""
    Z, Y, X = ot.shape
    hi_clip = np.array([Z - 1, Y - 1, X - 1], np.int64)
    if rule == "none":
        return None, None, "none"
    if rule in ("auto", "roi") and roi_box:
        lo = np.array(roi_box["min_xyz"], np.int64)[::-1]
        hi = np.array(roi_box["max_xyz"], np.int64)[::-1]
        return np.clip(lo, 0, hi_clip), np.clip(hi, 0, hi_clip), "roi_box"
    if rule == "roi":
        raise ValueError("coverage rule 'roi' needs an ROI box from the truth params")
    if ot.truth_bbox_min_zyx is None:
        return None, None, "no_truth"
    pad = np.ceil(float(dilate_um) / np.asarray(voxel_zyx, float)).astype(np.int64)
    lo = np.clip(ot.truth_bbox_min_zyx - pad, 0, hi_clip)
    hi = np.clip(ot.truth_bbox_max_zyx + pad, 0, hi_clip)
    return lo, hi, f"truth_bbox_dilated_{dilate_um:g}um"


def ignored_predictions(ot: OverlapTable, box_min, box_max) -> np.ndarray:
    """Predictions that overlap no truth object and whose centroid lies outside the coverage box.
    Without a box (``box_min`` None) nothing is ignored; without any truth object everything that
    overlaps nothing is ignored."""
    no_overlap = ~ot.pred_has_overlap()
    if box_min is None:
        return np.zeros(ot.n_pred, bool) if ot.truth_ids.size else no_overlap
    c = ot.pred_centroid_zyx
    inside = np.all((c >= np.asarray(box_min) - 0.5) & (c <= np.asarray(box_max) + 0.5), axis=1)
    return no_overlap & ~inside


def splits_and_merges(ot: OverlapTable, frac: float = DEFAULT_MAJORITY_FRAC):
    """Majority-overlap splits and merges. A prediction belongs to a truth object when more than
    ``frac`` of the prediction's voxels lie in it; a truth object with two or more such predictions
    is split. A truth object belongs to a prediction when more than ``frac`` of the truth voxels lie
    in it; a prediction with two or more such truth objects is a merge. Returns
    (n_pred_inside_per_truth, n_truth_inside_per_pred)."""
    pred_in = ot.pair_inter > frac * ot.pred_sizes[ot.pair_pred]
    n_pred_per_truth = np.bincount(ot.pair_truth[pred_in], minlength=ot.n_truth)
    truth_in = ot.pair_inter > frac * ot.truth_sizes[ot.pair_truth]
    n_truth_per_pred = np.bincount(ot.pair_pred[truth_in], minlength=ot.n_pred)
    return n_pred_per_truth, n_truth_per_pred


# ------------------------------------------------------------------------- the protocol's scope
def remove_truth_fragments(truth, min_voxels: int = PROTOCOL_RULES["min_truth_voxels"]):
    """(truth without the labels of fewer than ``min_voxels`` voxels, the removed label ids). The input
    is returned unchanged (not copied) when nothing is removed; otherwise a new array of its dtype is
    built plane by plane through a lookup table."""
    vol = _as_volume(truth)
    top = _max_label(vol)
    counts = np.zeros(top + 1, np.int64)
    for z in range(vol.shape[0]):
        counts += np.bincount(np.asarray(vol[z]).ravel().astype(np.int64, copy=False), minlength=top + 1)
    small = np.flatnonzero((counts > 0) & (counts < int(min_voxels)))
    small = small[small > 0]
    if small.size == 0:
        return vol, []
    lut = np.arange(top + 1, dtype=np.int64)
    lut[small] = 0
    out = np.empty(vol.shape, dtype=np.asarray(vol[0]).dtype)
    for z in range(vol.shape[0]):
        out[z] = lut[np.asarray(vol[z]).astype(np.int64, copy=False)]
    return out, [int(v) for v in small]


def protocol_coverage_mask(truth, voxel_zyx, *, xy_um: float = PROTOCOL_RULES["coverage_xy_um"],
                           z_planes: int = PROTOCOL_RULES["coverage_z_planes"]) -> np.ndarray:
    """Region C of the protocol: the truth masks dilated by ``xy_um`` in xy and ``z_planes`` planes in
    z (exact: per-plane Euclidean distances in um, restricted to the bounding box of the masks)."""
    from germquant.spots.net3d.data import dilate_cylinder

    vol = _as_volume(truth)
    m = np.zeros(vol.shape, bool)
    for z in range(vol.shape[0]):
        m[z] = np.asarray(vol[z]) > 0
    return dilate_cylinder(m, float(xy_um), int(z_planes), tuple(float(s) for s in voxel_zyx))


def voxels_inside(ot: OverlapTable, pred, region: np.ndarray) -> np.ndarray:
    """Voxels of every prediction (in ``ot.pred_ids`` order) that lie inside the boolean ``region``."""
    pred = _as_volume(pred)
    top = int(ot.pred_ids.max()) if ot.n_pred else 0
    counts = np.zeros(top + 1, np.int64)
    for z in range(pred.shape[0]):
        r = region[z]
        if r.any():
            counts += np.bincount(np.asarray(pred[z])[r].astype(np.int64, copy=False), minlength=top + 1)
    return counts[ot.pred_ids] if ot.n_pred else np.zeros(0, np.int64)


def protocol_splits_and_merges(ot: OverlapTable, in_scope: np.ndarray,
                               frac: float = PROTOCOL_RULES["split_frac"]):
    """Section 5.3 splits and merges: a pair counts when the in-scope prediction covers at least
    ``frac`` of the truth nucleus's volume. Returns (n_pred_covering_per_truth, n_truth_covered_per_pred);
    a truth nucleus with two or more is split, a prediction with two or more is a merge."""
    ok = in_scope[ot.pair_pred] & (ot.pair_inter >= frac * ot.truth_sizes[ot.pair_truth] - 1e-9)
    n_pred_per_truth = np.bincount(ot.pair_truth[ok], minlength=ot.n_truth)
    n_truth_per_pred = np.bincount(ot.pair_pred[ok], minlength=ot.n_pred)
    return n_pred_per_truth, n_truth_per_pred


# ----------------------------------------------------------------------------------- one gonad
def benchmark_gonad(pred, truth, voxel_zyx, *, tag: str = "", split: str = "",
                    roi_box: dict | None = None, coverage: str = "auto",
                    dilate_um: float = DEFAULT_DILATE_UM, majority_frac: float = DEFAULT_MAJORITY_FRAC,
                    iou_thresholds=IOU_THRESHOLDS,
                    max_dense_cells: int = DEFAULT_MAX_DENSE_CELLS,
                    min_truth_voxels: int = PROTOCOL_RULES["min_truth_voxels"],
                    coverage_xy_um: float = PROTOCOL_RULES["coverage_xy_um"],
                    coverage_z_planes: int = PROTOCOL_RULES["coverage_z_planes"],
                    scope_frac: float = PROTOCOL_RULES["scope_frac"],
                    split_frac: float = PROTOCOL_RULES["split_frac"],
                    protocol_region: np.ndarray | None = None) -> tuple[dict, pd.DataFrame]:
    """Score one gonad. Returns (per-gonad row, per-truth-nucleus table). ``coverage`` is one of
    ``COVERAGE_CHOICES``: the box rules of the module docstring ("auto", "roi", "bbox", "none") or
    "protocol", the pre-registered rule of section 5.3 with its parameters ``min_truth_voxels``,
    ``coverage_xy_um``, ``coverage_z_planes``, ``scope_frac`` and ``split_frac`` (ignored by the box
    rules). ``protocol_region`` is region C computed beforehand with ``protocol_coverage_mask`` from the
    same truth with its fragments removed (a caller that scores many label images against one truth,
    such as the nuclei3d validation, saves the dilation each time); None computes it here."""
    if coverage not in COVERAGE_CHOICES:
        raise ValueError(f"coverage must be one of {COVERAGE_CHOICES}, got {coverage!r}")
    voxel_zyx = tuple(float(s) for s in voxel_zyx)
    dz = voxel_zyx[0]
    vox_um3 = float(np.prod(voxel_zyx))
    proto = coverage == "protocol"
    fragments: list[int] = []
    if proto:
        truth, fragments = remove_truth_fragments(truth, min_truth_voxels)
    ot = overlap_table(pred, truth)
    in_scope = None
    if proto:
        if protocol_region is not None:
            region = np.asarray(protocol_region, bool)
            if region.shape != tuple(ot.shape):
                raise ValueError(f"protocol_region {region.shape} and the label images {ot.shape} differ in shape")
        else:
            region = protocol_coverage_mask(truth, voxel_zyx, xy_um=coverage_xy_um, z_planes=coverage_z_planes)
        inside = voxels_inside(ot, pred, region)
        del region
        in_scope = inside >= float(scope_frac) * ot.pred_sizes - 1e-9
        ignored = ~in_scope
        box_min = box_max = None
        cov_rule = (f"protocol_truth_dilated_{coverage_xy_um:g}um_xy_{int(coverage_z_planes)}plane_z_"
                    f"scope_{scope_frac:g}_of_volume")
        n_pred_rule, n_truth_rule = protocol_splits_and_merges(ot, in_scope, split_frac)
        split_rule = f"in_scope_cover_{split_frac:g}_of_truth"
    else:
        box_min, box_max, cov_rule = coverage_box(ot, voxel_zyx, dilate_um=dilate_um, roi_box=roi_box,
                                                  rule=coverage)
        ignored = ignored_predictions(ot, box_min, box_max)
        split_rule = f"majority_{majority_frac:g}"
    n_ignored = int(ignored.sum())
    n_counted = ot.n_pred - n_ignored
    n_pred_in, n_truth_in = splits_and_merges(ot, majority_frac)
    if not proto:
        n_pred_rule, n_truth_rule = n_pred_in, n_truth_in

    row: dict = {
        "tag": tag, "split": split, "n_pred": ot.n_pred, "n_truth": ot.n_truth,
        "n_pred_ignored_outside_coverage": n_ignored, "n_pred_counted": int(n_counted),
        "n_pred_no_overlap_in_coverage": int(((~ot.pred_has_overlap()) & ~ignored).sum()),
        "coverage_rule": cov_rule,
        "coverage_min_zyx": "" if box_min is None else " ".join(str(int(v)) for v in box_min),
        "coverage_max_zyx": "" if box_max is None else " ".join(str(int(v)) for v in box_max),
        "n_splits": int((n_pred_rule >= 2).sum()), "n_merges": int((n_truth_rule >= 2).sum()),
    }
    ap_terms = []
    methods = set()
    a50 = None
    for t in iou_thresholds:
        a = match_instances(ot, t, max_dense_cells=max_dense_cells, pred_ok=in_scope)
        methods.add(a.method)
        tp = a.n_matched
        fp = int(n_counted - tp)
        fn = int(ot.n_truth - tp)
        s = precision_recall_f1(tp, n_counted, ot.n_truth)
        k = _tag(t)
        row[f"tp_{k}"], row[f"fp_{k}"], row[f"fn_{k}"] = tp, fp, fn
        row[f"f1_{k}"] = s["f1"]
        denom = tp + fp + fn
        ap_terms.append(tp / denom if denom else float("nan"))
        if abs(t - 0.5) < 1e-9:
            a50 = a
            row["precision_iou50"], row["recall_iou50"] = s["precision"], s["recall"]
            row["mean_matched_iou_iou50"] = (float(np.mean(1.0 - a.cost)) if a.n_matched
                                             else float("nan"))
    row["ap_iou50_95"] = float(np.nanmean(ap_terms)) if np.isfinite(ap_terms).any() else float("nan")
    row["matching_method"] = "+".join(sorted(m for m in methods if m != "none")) or "none"

    # per truth nucleus (matched at IoU 0.5), with the z-extent error of matched pairs
    n_t = ot.n_truth
    matched_pred = np.full(n_t, -1, np.int64)
    matched_iou = np.full(n_t, np.nan)
    if a50 is not None and a50.n_matched:
        matched_pred[a50.cols] = a50.rows
        matched_iou[a50.cols] = 1.0 - a50.cost
    truth_span = (ot.truth_zmax - ot.truth_zmin + 1) * dz
    has = matched_pred >= 0
    pred_span = np.full(n_t, np.nan)
    pred_vox = np.full(n_t, np.nan)
    pred_label = np.zeros(n_t, np.int64)
    if has.any():
        mp = matched_pred[has]
        pred_span[has] = (ot.pred_zmax[mp] - ot.pred_zmin[mp] + 1) * dz
        pred_vox[has] = ot.pred_sizes[mp]
        pred_label[has] = ot.pred_ids[mp]
    z_err = pred_span - truth_span
    per_truth = pd.DataFrame({
        "tag": tag, "split": split,
        "truth_label": ot.truth_ids, "truth_voxels": ot.truth_sizes,
        "truth_volume_um3": ot.truth_sizes * vox_um3, "truth_z_span_um": truth_span,
        "matched_iou50": has, "pred_label": pred_label, "iou": matched_iou,
        "pred_voxels": pred_vox, "pred_z_span_um": pred_span, "z_extent_error_um": z_err,
        "n_pred_majority_inside": n_pred_in, "is_split": n_pred_rule >= 2,
        "n_pred_split_rule": n_pred_rule, "z_extent_error_slices": np.rint(z_err / dz),
    })
    ze = z_err[np.isfinite(z_err)]
    row["n_z_pairs"] = int(ze.size)
    row["z_extent_error_mean_um"] = float(ze.mean()) if ze.size else float("nan")
    row["z_extent_error_median_um"] = float(np.median(ze)) if ze.size else float("nan")
    row["z_extent_error_mean_abs_um"] = float(np.abs(ze).mean()) if ze.size else float("nan")
    row["voxel_dz_um"], row["voxel_dy_um"], row["voxel_dx_um"] = voxel_zyx
    # appended columns (the earlier ones keep their order and meaning)
    row["z_extent_error_mean_slices"] = float(ze.mean() / dz) if ze.size else float("nan")
    row["z_extent_error_mean_abs_slices"] = float(np.abs(ze).mean() / dz) if ze.size else float("nan")
    row["split_rule"] = split_rule
    row["n_truth_fragments_removed"] = len(fragments)
    row["truth_fragments_removed"] = " ".join(str(v) for v in fragments)
    return row, per_truth


# ------------------------------------------------------------------------------------ pooling
def _pool_rows(rows: pd.DataFrame, per_truth: pd.DataFrame, *, n_boot: int, seed: int,
               iou_thresholds) -> dict:
    tags = list(rows["tag"])
    det = [(int(r["tp_iou50"]), int(r["n_pred_counted"]), int(r["n_truth"])) for _, r in rows.iterrows()]
    boot_det = cluster_bootstrap(det, pooled_detection, n_boot=n_boot, seed=seed)["stats"]

    def ap_stat(payloads):
        terms = []
        for t in iou_thresholds:
            k = _tag(t)
            tp = sum(p[k][0] for p in payloads)
            fp = sum(p[k][1] for p in payloads)
            fn = sum(p[k][2] for p in payloads)
            d = tp + fp + fn
            terms.append(tp / d if d else float("nan"))
        return float(np.nanmean(terms)) if np.isfinite(terms).any() else float("nan")

    ap_payloads = [{_tag(t): (int(r[f"tp_{_tag(t)}"]), int(r[f"fp_{_tag(t)}"]), int(r[f"fn_{_tag(t)}"]))
                    for t in iou_thresholds} for _, r in rows.iterrows()]
    boot_ap = cluster_bootstrap(ap_payloads, ap_stat, n_boot=n_boot, seed=seed)["stats"]["value"]
    z_payloads = [per_truth.loc[per_truth["tag"] == t, "z_extent_error_um"].dropna().to_numpy()
                  for t in tags]
    boot_z = cluster_bootstrap(z_payloads, pooled_mean, n_boot=n_boot, seed=seed)["stats"]["value"]
    tp = int(rows["tp_iou50"].sum())
    n_counted = int(rows["n_pred_counted"].sum())
    n_truth = int(rows["n_truth"].sum())
    return {
        "n_gonads": int(len(rows)), "tags": tags,
        "n_pred": int(rows["n_pred"].sum()), "n_pred_counted": n_counted,
        "n_pred_ignored_outside_coverage": int(rows["n_pred_ignored_outside_coverage"].sum()),
        "n_truth": n_truth, "tp_iou50": tp, "fp_iou50": int(n_counted - tp), "fn_iou50": int(n_truth - tp),
        "precision_iou50": boot_det["precision"], "recall_iou50": boot_det["recall"],
        "f1_iou50": boot_det["f1"], "ap_iou50_95": boot_ap,
        "mean_per_gonad_f1_iou50": float(rows["f1_iou50"].mean()) if len(rows) else float("nan"),
        "mean_per_gonad_ap": float(rows["ap_iou50_95"].mean()) if len(rows) else float("nan"),
        "n_splits": int(rows["n_splits"].sum()), "n_merges": int(rows["n_merges"].sum()),
        "z_extent_error_um": boot_z, "n_z_pairs": int(rows["n_z_pairs"].sum()),
        "coverage_rules": {str(r["tag"]): str(r["coverage_rule"]) for _, r in rows.iterrows()},
        "matching_methods": {str(r["tag"]): str(r["matching_method"]) for _, r in rows.iterrows()},
    }


def pool(rows: pd.DataFrame, per_truth: pd.DataFrame, *, n_boot: int = 1000, seed: int = 0,
         iou_thresholds=IOU_THRESHOLDS) -> dict:
    """Pooled scores over gonads with cluster-bootstrap intervals, overall and by split."""
    out = _pool_rows(rows, per_truth, n_boot=n_boot, seed=seed, iou_thresholds=iou_thresholds)
    out["bootstrap"] = {"n_boot": int(n_boot), "seed": int(seed), "unit": "gonad",
                        "interval": "2.5 and 97.5 percentiles"}
    splits = [s for s in rows["split"].astype(str).unique() if s]
    if len(splits) > 1:
        out["by_split"] = {s: _pool_rows(rows[rows["split"].astype(str) == s],
                                         per_truth[per_truth["split"].astype(str) == s],
                                         n_boot=n_boot, seed=seed, iou_thresholds=iou_thresholds)
                           for s in splits}
    return out
