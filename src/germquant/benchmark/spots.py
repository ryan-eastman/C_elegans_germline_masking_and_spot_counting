"""RAD-51 spot benchmark: the expert's Imaris spots against our detector, per nucleus and per spot.

Per nucleus: every expert spot is mapped into OUR nucleus label image (the label at the spot's voxel,
with a fallback to the most frequent label in a 1-voxel-in-z, 2-voxel-in-xy neighbourhood, as the
prototype did), which gives an expert count per nucleus to set against our own ``n_spots``. Pearson's
r, Lin's CCC and the bias (ours minus expert) are reported over all germline nuclei and over the
traced subset: our nuclei whose label overlaps a truth label by more than half their volume when the
truth label image exists, else our nuclei whose centroid lies within 2.5 um of a truth surface
centroid.

Per spot: one-to-one matching of our spots to the expert's inside an anisotropic cylinder (0.3 um in
xy, 0.4 um in z by default; about three xy pixels and two z planes on the 0.108 x 0.108 x 0.2 um grid)
by the Hungarian algorithm on the connected components of the candidate graph (``matching``). The
cost of a candidate pair is sqrt((dxy / tol_xy)^2 + (dz / tol_z)^2), so a pair at the edge of the
cylinder costs about 1. Precision, recall and F1 are reported for all expert spots; without the
expert's Class B spots when the class column exists (Class B spots and the spots we matched to them
are left out of both counts, so finding them is neither rewarded nor punished); and, when the truth
label image exists, "detector-only": expert spots inside the traced nuclei against our spots inside
the same traced nuclei, which removes our nucleus segmentation from the comparison.

Both sets of coordinates are in the pipeline frame (voxel index times the .nd2 voxel size): the
truth files already carry that conversion (see ``manifest``).

The protocol's primary partition (docs/BENCHMARK_PROTOCOL.md, section 5.1, "traced subset"). The
"traced" columns above are this harness's older subset: OUR germline nuclei with more than half of their
volume inside one truth mask, counted in OUR labels. That is not the section 5.1 partition. Section 5.1
takes the truth masks themselves as the nuclei, the same for every method: the truth label image with
fragments below 1000 voxels removed (section 4.2); every truth spot and every one of the method's spots is
assigned by the section 5 rule (the mask label at its centre voxel; if that voxel is background, the most
frequent non-zero label in the 3 x 3 x 3 neighbourhood; otherwise unassigned); truth count and method
count are the spots assigned to each truth nucleus. ``truth_mask_table`` computes it and the ``tm_``
columns of ``benchmark_gonad`` summarise it (Pearson r, Lin's CCC, bias and mean absolute error, for the
truth counted from Class A+B, Class A and Class B spots). The spot region R of section 5 (the masks dilated
by one voxel in every axis, a 3 x 3 x 3 cube) is exactly the set of centre voxels the rule assigns, so "only
spots in R take part in the counts" needs no separate step. With ``detector_only_spots`` (the method's
spot stage re-run with the truth masks in place of its own nuclei) the detector-only count is reported next
to the end-to-end one (``tm_det_`` columns).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .manifest import PipelineOutputs, Truth, points_to_voxel_index
from .matching import DEFAULT_MAX_CANDIDATES, DEFAULT_MAX_DENSE_CELLS, Assignment, assign_sparse
from .stats import agreement_metrics, cluster_bootstrap, pooled_agreement, pooled_detection, precision_recall_f1

DEFAULT_TOL_XY_UM = 0.3
DEFAULT_TOL_Z_UM = 0.4
DEFAULT_SURFACE_TOL_UM = 2.5
DEFAULT_OVERLAP_FRAC = 0.5
CENTROID_COLS = ["centroid_x_um", "centroid_y_um", "centroid_z_um"]
# protocol section 4.2: the scoring code ignores truth mask fragments below 1000 voxels
TRUTH_MASK_MIN_VOXELS = 1000
# the per-truth-nucleus table of the section 5.1 partition
TRUTH_MASK_COLUMNS = ["truth_label", "n_voxels", "expert_n", "expert_n_class_a", "expert_n_class_b", "ours_n",
                      "ours_n_detector_only"]
_TM_STATS = ("mean_expert", "mean_ours", "pearson_r", "ccc", "bias", "mae")
_TM_CLASS_STATS = ("mean_expert", "pearson_r", "ccc", "bias", "mae")


# ------------------------------------------------------------------------------ label lookups
def label_at_points(labels, zi, yi, xi, *, fallback_z: int = 1, fallback_xy: int = 2) -> np.ndarray:
    """Label under each (zi, yi, xi) voxel; 0 outside the image. Where the voxel itself is background
    the most frequent non-zero label in the (2*fallback_z+1) x (2*fallback_xy+1)^2 neighbourhood is
    taken (ties go to the smallest label), so a spot sitting on the edge of a nucleus still counts."""
    zi = np.asarray(zi, np.int64).ravel()
    yi = np.asarray(yi, np.int64).ravel()
    xi = np.asarray(xi, np.int64).ravel()
    Z, Y, X = labels.shape
    out = np.zeros(len(zi), np.int64)
    inb = (zi >= 0) & (zi < Z) & (yi >= 0) & (yi < Y) & (xi >= 0) & (xi < X)
    idx = np.flatnonzero(inb)
    if idx.size:
        out[idx] = np.asarray(labels[zi[idx], yi[idx], xi[idx]]).astype(np.int64)
    for k in idx[out[idx] == 0]:
        z, y, x = int(zi[k]), int(yi[k]), int(xi[k])
        block = labels[max(z - fallback_z, 0): z + fallback_z + 1,
                       max(y - fallback_xy, 0): y + fallback_xy + 1,
                       max(x - fallback_xy, 0): x + fallback_xy + 1]
        block = np.asarray(block)
        vals = block[block > 0]
        if vals.size:
            u, c = np.unique(vals, return_counts=True)
            out[k] = int(u[np.argmax(c)])
    return out


def germline_mask(nuclei: pd.DataFrame) -> np.ndarray:
    """``in_germline`` as booleans, or every nucleus when the column is absent."""
    if "in_germline" not in nuclei.columns:
        return np.ones(len(nuclei), bool)
    col = nuclei["in_germline"]
    if col.dtype == bool:
        return col.to_numpy()
    return col.astype(str).str.strip().str.lower().isin(["true", "1", "yes"]).to_numpy()


def ours_counts(nuclei: pd.DataFrame, spots: pd.DataFrame | None) -> np.ndarray:
    """Our spot count per nucleus row: the ``n_spots`` column, else counted from the spots table."""
    if "n_spots" in nuclei.columns:
        return nuclei["n_spots"].fillna(0).astype(int).to_numpy()
    if spots is not None and "nucleus_id" in spots.columns:
        counts = spots["nucleus_id"].dropna().astype(int).value_counts()
        return nuclei["nucleus_id"].map(counts).fillna(0).astype(int).to_numpy()
    return np.zeros(len(nuclei), int)


# --------------------------------------------------------------------------------- the traced
def traced_by_centroid(our_xyz, truth_xyz, tol_um: float = DEFAULT_SURFACE_TOL_UM) -> np.ndarray:
    """Our nuclei with a truth surface centroid within ``tol_um``."""
    our_xyz = np.asarray(our_xyz, float).reshape(-1, 3)
    truth_xyz = np.asarray(truth_xyz, float).reshape(-1, 3)
    if len(our_xyz) == 0 or len(truth_xyz) == 0:
        return np.zeros(len(our_xyz), bool)
    d, _ = cKDTree(truth_xyz).query(our_xyz, distance_upper_bound=tol_um)
    return np.isfinite(d)


def traced_by_overlap(labels, truth_labels, nucleus_ids, frac: float = DEFAULT_OVERLAP_FRAC) -> np.ndarray:
    """Our nuclei (by label) with more than ``frac`` of their voxels inside one truth label."""
    from .nuclei3d import overlap_table

    ot = overlap_table(labels, truth_labels)
    best = np.zeros(ot.n_pred, float)
    if ot.pair_inter.size:
        np.maximum.at(best, ot.pair_pred, ot.pair_inter / np.maximum(ot.pred_sizes[ot.pair_pred], 1))
    traced_by_label = dict(zip(ot.pred_ids.tolist(), (best > frac).tolist()))
    return np.array([bool(traced_by_label.get(int(n), False)) for n in nucleus_ids], bool)


def per_nucleus_table(nuclei: pd.DataFrame, spot_labels: np.ndarray, ours_n: np.ndarray, *,
                      truth_nuclei: pd.DataFrame | None = None, labels=None, truth_labels=None,
                      surface_tol_um: float = DEFAULT_SURFACE_TOL_UM,
                      overlap_frac: float = DEFAULT_OVERLAP_FRAC,
                      traced_rule: str = "auto") -> tuple[pd.DataFrame, str]:
    """One row per germline nucleus with ``expert_n`` (expert spots mapped into its label),
    ``ours_n`` and ``traced``. ``traced_rule``: "auto" (overlap when ``truth_labels`` is given, else
    centroid), "overlap", "centroid" or "none". Returns (table, rule used)."""
    germ = germline_mask(nuclei)
    tbl = nuclei.loc[germ, ["nucleus_id"] + [c for c in CENTROID_COLS if c in nuclei.columns]].copy()
    tbl = tbl.reset_index(drop=True)
    expert_per_label = pd.Series(spot_labels[spot_labels > 0]).value_counts()
    tbl["expert_n"] = tbl["nucleus_id"].map(expert_per_label).fillna(0).astype(int)
    tbl["ours_n"] = np.asarray(ours_n)[germ].astype(int)

    rule = traced_rule
    if rule == "auto":
        rule = "overlap" if (truth_labels is not None and labels is not None) else (
            "centroid" if truth_nuclei is not None and len(truth_nuclei) else "none")
    if rule == "overlap":
        if truth_labels is None or labels is None:
            raise ValueError("traced_rule 'overlap' needs both label images")
        tbl["traced"] = traced_by_overlap(labels, truth_labels, tbl["nucleus_id"].to_numpy(), overlap_frac)
    elif rule == "centroid":
        if truth_nuclei is None or not all(c in tbl.columns for c in CENTROID_COLS):
            tbl["traced"] = False
            rule = "none"
        else:
            tbl["traced"] = traced_by_centroid(tbl[CENTROID_COLS].to_numpy(),
                                               truth_nuclei[["x_um", "y_um", "z_um"]].to_numpy(),
                                               surface_tol_um)
    else:
        tbl["traced"] = False
        rule = "none"
    return tbl, rule


# ------------------------------------------------------- the protocol's section 5.1 partition
def label_sizes(labels) -> np.ndarray:
    """Voxels per label value (index = label value), counted plane by plane so that a full-size label
    image is never copied as a whole."""
    top = 0
    for z in range(labels.shape[0]):
        plane = np.asarray(labels[z])
        if plane.size:
            if int(plane.min()) < 0:
                raise ValueError("label images must not hold negative labels")
            top = max(top, int(plane.max()))
    counts = np.zeros(top + 1, np.int64)
    for z in range(labels.shape[0]):
        counts += np.bincount(np.asarray(labels[z]).ravel(), minlength=top + 1)
    return counts


def truth_mask_lut(truth_labels, min_voxels: int = TRUTH_MASK_MIN_VOXELS):
    """The truth masks of protocol section 5: the recovered Imaris masks with fragments below
    ``min_voxels`` voxels removed (section 4.2). Returns (lut, kept labels, their voxel counts, removed
    fragment labels); ``lut[v]`` is ``v`` for a kept label and 0 for the background and for a fragment."""
    sizes = label_sizes(truth_labels)
    ids = np.flatnonzero(sizes[1:]) + 1 if len(sizes) > 1 else np.zeros(0, np.int64)
    keep = ids[sizes[ids] >= int(min_voxels)]
    frag = ids[sizes[ids] < int(min_voxels)]
    lut = np.zeros(len(sizes), np.int64)
    lut[keep] = keep
    return lut, keep, sizes[keep], frag


def truth_mask_labels(truth_labels, min_voxels: int = TRUTH_MASK_MIN_VOXELS) -> np.ndarray:
    """The truth masks as an int32 label image (fragments set to 0), built plane by plane: the label image
    a detector-only run hands the method's spot stage in place of its own nuclei."""
    lut = truth_mask_lut(truth_labels, min_voxels)[0].astype(np.int32)
    out = np.zeros(truth_labels.shape, np.int32)
    for z in range(truth_labels.shape[0]):
        out[z] = lut[np.asarray(truth_labels[z])]
    return out


_OFFSETS_333 = np.array([(dz, dy, dx) for dz in (-1, 0, 1) for dy in (-1, 0, 1) for dx in (-1, 0, 1)], np.int64)
_CENTRE_333 = 13


def assign_to_masks(labels, zi, yi, xi, *, lut=None) -> np.ndarray:
    """The spot to nucleus assignment of protocol section 5 for centre voxels (zi, yi, xi): the mask label
    at the centre voxel; if that voxel is background, the most frequent non-zero label in its 3 x 3 x 3
    neighbourhood (ties go to the smallest label); otherwise 0, unassigned. Voxels outside the image are
    background. ``lut`` (``truth_mask_lut``) maps raw label values to the masks kept, so a removed
    fragment is background everywhere, at the centre and in the neighbourhood alike."""
    zi = np.asarray(zi, np.int64).ravel()
    yi = np.asarray(yi, np.int64).ravel()
    xi = np.asarray(xi, np.int64).ravel()
    n = len(zi)
    out = np.zeros(n, np.int64)
    if n == 0:
        return out
    Z, Y, X = labels.shape
    vals = np.zeros((n, len(_OFFSETS_333)), np.int64)
    for k, (dz, dy, dx) in enumerate(_OFFSETS_333):
        z, y, x = zi + dz, yi + dy, xi + dx
        ok = np.flatnonzero((z >= 0) & (z < Z) & (y >= 0) & (y < Y) & (x >= 0) & (x < X))
        if ok.size:
            vals[ok, k] = np.asarray(labels[z[ok], y[ok], x[ok]]).astype(np.int64)
    if lut is not None:
        lut = np.asarray(lut, np.int64)
        inside = vals < len(lut)
        vals = np.where(inside, lut[np.where(inside, vals, 0)], 0)
    out = vals[:, _CENTRE_333].copy()
    for i in np.flatnonzero((out == 0) & (vals > 0).any(axis=1)):
        nz = vals[i][vals[i] > 0]
        u, c = np.unique(nz, return_counts=True)
        out[i] = int(u[np.argmax(c)])
    return out


def _assign_xyz(labels, xyz, voxel_zyx, lut) -> np.ndarray:
    """``assign_to_masks`` for points in the pipeline frame (x, y, z um); centre voxel = round(coordinate
    / voxel), the convention by which the pipeline forms its coordinates (protocol section 4.3)."""
    xyz = np.asarray(xyz, float).reshape(-1, 3)
    if not len(xyz):
        return np.zeros(0, np.int64)
    idx = points_to_voxel_index(xyz, np.asarray(voxel_zyx, float)[::-1])
    return assign_to_masks(labels, idx[:, 2], idx[:, 1], idx[:, 0], lut=lut)


def truth_mask_table(truth_labels, truth_spots: pd.DataFrame, ours_xyz, voxel_zyx, *,
                     min_voxels: int = TRUTH_MASK_MIN_VOXELS, has_class: bool | None = None,
                     detector_only_xyz=None) -> tuple[pd.DataFrame, dict]:
    """The primary partition of protocol section 5.1 for one gonad: the truth masks are the nuclei.

    ``truth_labels`` is the truth label image (ZYX; fragments below ``min_voxels`` are removed here),
    ``truth_spots`` the truth spot table (centre voxels from its ``zi, yi, xi`` columns, else rounded from
    ``x_um, y_um, z_um``), ``ours_xyz`` the method's end-to-end spots (x, y, z um, the pipeline frame) and
    ``detector_only_xyz`` optionally the spots of its detector-only run (the spot stage re-run with the
    truth masks in place of its own nuclei). Every spot is assigned by ``assign_to_masks``. Returns one row
    per truth nucleus with ``TRUTH_MASK_COLUMNS`` (the truth counts from Class A+B, Class A and Class B
    spots, NaN for the classes without a class column; the method's count; the detector-only count, NaN
    without it) and the per-gonad facts (nuclei, fragments removed, assigned and unassigned spots)."""
    voxel_zyx = tuple(float(s) for s in voxel_zyx)
    lut, keep, sizes, frag = truth_mask_lut(truth_labels, min_voxels)
    nbins = len(lut)
    if len(truth_spots):
        tz, ty, tx = _spot_voxels(truth_spots, voxel_zyx)
    else:
        tz = ty = tx = np.zeros(0, np.int64)
    t_lab = assign_to_masks(truth_labels, tz, ty, tx, lut=lut)
    if has_class is None:
        has_class = "spot_class" in truth_spots.columns and bool(
            (truth_spots["spot_class"].astype(str).str.strip() != "").any())
    b = is_class_b(truth_spots["spot_class"]) if (has_class and len(truth_spots)) else np.zeros(len(t_lab), bool)

    def count(lab: np.ndarray, sel=None) -> np.ndarray:
        ok = lab > 0 if sel is None else (lab > 0) & sel
        return np.bincount(lab[ok], minlength=nbins)[keep].astype(np.int64)

    o_lab = _assign_xyz(truth_labels, ours_xyz, voxel_zyx, lut)
    nan_col = np.full(len(keep), np.nan)
    table = pd.DataFrame({
        "truth_label": keep.astype(np.int64), "n_voxels": np.asarray(sizes, np.int64),
        "expert_n": count(t_lab),
        "expert_n_class_a": count(t_lab, ~b) if has_class else nan_col,
        "expert_n_class_b": count(t_lab, b) if has_class else nan_col.copy(),
        "ours_n": count(o_lab),
        "ours_n_detector_only": nan_col.copy(),
    }, columns=TRUTH_MASK_COLUMNS)
    facts = {"min_voxels": int(min_voxels), "n_nuclei": int(len(keep)), "n_fragments_removed": int(len(frag)),
             "fragment_labels": [int(v) for v in frag], "has_class": bool(has_class),
             "n_expert_spots_assigned": int((t_lab > 0).sum()), "n_expert_spots_unassigned": int((t_lab == 0).sum()),
             "n_our_spots_assigned": int((o_lab > 0).sum()), "n_our_spots_unassigned": int((o_lab == 0).sum()),
             "detector_only": detector_only_xyz is not None}
    if detector_only_xyz is not None:
        d_lab = _assign_xyz(truth_labels, detector_only_xyz, voxel_zyx, lut)
        table["ours_n_detector_only"] = count(d_lab)
        facts["n_our_spots_detector_only_assigned"] = int((d_lab > 0).sum())
        facts["n_our_spots_detector_only_unassigned"] = int((d_lab == 0).sum())
    return table, facts


def truth_mask_row(table: pd.DataFrame | None, facts: dict | None) -> dict:
    """The per-gonad ``tm_`` columns of the section 5.1 partition: the number of truth nuclei, the
    fragments removed, the truth and method spots assigned and unassigned, and the agreement (``_TM_STATS``:
    means, Pearson r, Lin's CCC, bias = mean of method minus truth, mean absolute error) of the method's
    end-to-end count with the truth counted from Class A+B (``tm_``), Class A (``tm_class_a_``) and Class B
    (``tm_class_b_``) spots. With a detector-only count, the same under ``tm_det_``. NaN without a truth
    label image."""
    nan = float("nan")
    have = table is not None and facts is not None
    row: dict = {"tm_available": bool(have), "tm_min_voxels": facts["min_voxels"] if have else nan}
    for k in ("n_nuclei", "n_fragments_removed", "n_expert_spots_assigned", "n_expert_spots_unassigned",
              "n_our_spots_assigned", "n_our_spots_unassigned"):
        row[f"tm_{k}"] = facts[k] if have else nan
    with_det = have and facts.get("detector_only")
    if with_det:
        row["tm_n_our_spots_detector_only_assigned"] = facts["n_our_spots_detector_only_assigned"]
        row["tm_n_our_spots_detector_only_unassigned"] = facts["n_our_spots_detector_only_unassigned"]
    for prefix, ours_col in (("tm_", "ours_n"), ("tm_det_", "ours_n_detector_only")):
        if prefix == "tm_det_" and not with_det:
            continue
        for cls, truth_col, keys in (("", "expert_n", _TM_STATS), ("class_a_", "expert_n_class_a", _TM_CLASS_STATS),
                                     ("class_b_", "expert_n_class_b", _TM_CLASS_STATS)):
            usable = have and table[truth_col].notna().all() and len(table)
            m = agreement_metrics(table[truth_col], table[ours_col]) if usable else {}
            for k in keys:
                row[f"{prefix}{cls}{k}"] = m.get(k, nan)
    return row


# ------------------------------------------------------------------------------ spot matching
def candidate_pairs(a_xyz, b_xyz, tol_xy_um: float = DEFAULT_TOL_XY_UM, tol_z_um: float = DEFAULT_TOL_Z_UM):
    """All (i, j) with points a[i] and b[j] inside the cylinder (dxy <= tol_xy and |dz| <= tol_z),
    with the scaled cost sqrt((dxy/tol_xy)^2 + (dz/tol_z)^2), dxy and dz. Found with a KD-tree on
    the coordinates scaled by the tolerances (the cylinder sits inside the ball of radius sqrt 2)."""
    a = np.asarray(a_xyz, float).reshape(-1, 3)
    b = np.asarray(b_xyz, float).reshape(-1, 3)
    empty = (np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0), np.zeros(0), np.zeros(0))
    if len(a) == 0 or len(b) == 0:
        return empty
    scale = np.array([1.0 / tol_xy_um, 1.0 / tol_xy_um, 1.0 / tol_z_um])
    rec = cKDTree(a * scale).sparse_distance_matrix(cKDTree(b * scale), np.sqrt(2.0) + 1e-9,
                                                    output_type="ndarray")
    i = np.asarray(rec["i"], np.int64)
    j = np.asarray(rec["j"], np.int64)
    if i.size == 0:
        return empty
    d = a[i] - b[j]
    dxy = np.hypot(d[:, 0], d[:, 1])
    dz = np.abs(d[:, 2])
    ok = (dxy <= tol_xy_um + 1e-9) & (dz <= tol_z_um + 1e-9)
    i, j, dxy, dz = i[ok], j[ok], dxy[ok], dz[ok]
    cost = np.sqrt((dxy / tol_xy_um) ** 2 + (dz / tol_z_um) ** 2)
    return i, j, cost, dxy, dz


def match_spots(ours_xyz, truth_xyz, *, tol_xy_um: float = DEFAULT_TOL_XY_UM,
                tol_z_um: float = DEFAULT_TOL_Z_UM, max_candidates: int | None = DEFAULT_MAX_CANDIDATES,
                max_dense_cells: int = DEFAULT_MAX_DENSE_CELLS) -> Assignment:
    """One-to-one matching of our spots (rows) to expert spots (cols) inside the cylinder."""
    i, j, cost, _, _ = candidate_pairs(ours_xyz, truth_xyz, tol_xy_um, tol_z_um)
    n_a = len(np.asarray(ours_xyz).reshape(-1, 3))
    n_b = len(np.asarray(truth_xyz).reshape(-1, 3))
    if i.size == 0:
        z = np.zeros(0, np.int64)
        return Assignment(z, z.copy(), np.zeros(0, float), 0, 0, 0, "none", [])
    return assign_sparse(i, j, cost, n_a, n_b, max_dense_cells=max_dense_cells,
                         max_candidates=max_candidates)


def is_class_b(spot_class) -> np.ndarray:
    """True for Imaris "Class B" spots (the column holds "B" or "Class B", any case)."""
    s = pd.Series(spot_class).fillna("").astype(str).str.strip().str.upper()
    s = s.str.replace("CLASS", "", regex=False).str.strip()
    return (s == "B").to_numpy()


def _scores(n_matched: int, n_ours: int, n_truth: int, prefix: str) -> dict:
    s = precision_recall_f1(n_matched, n_ours, n_truth)
    return {f"{prefix}n_matched": s["n_matched"], f"{prefix}n_ours": s["n_pred"],
            f"{prefix}n_expert": s["n_truth"], f"{prefix}precision": s["precision"],
            f"{prefix}recall": s["recall"], f"{prefix}f1": s["f1"]}


def _nan_scores(prefix: str) -> dict:
    nan = float("nan")
    return {f"{prefix}n_matched": 0, f"{prefix}n_ours": 0, f"{prefix}n_expert": 0,
            f"{prefix}precision": nan, f"{prefix}recall": nan, f"{prefix}f1": nan}


def _spot_voxels(truth_spots: pd.DataFrame, voxel_zyx) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(zi, yi, xi) of the expert spots: the contract's integer columns, else rounded from x_um."""
    cols = ["zi", "yi", "xi"]
    if all(c in truth_spots.columns for c in cols) and truth_spots[cols].notna().all().all():
        return tuple(truth_spots[c].to_numpy().astype(np.int64) for c in cols)
    idx = points_to_voxel_index(truth_spots[["x_um", "y_um", "z_um"]].to_numpy(),
                                np.asarray(voxel_zyx, float)[::-1])
    return idx[:, 2], idx[:, 1], idx[:, 0]


# ----------------------------------------------------------------------------------- one gonad
def benchmark_gonad(outputs: PipelineOutputs, truth: Truth, labels, voxel_zyx, *, tag: str = "",
                    split: str = "", truth_labels=None, tol_xy_um: float = DEFAULT_TOL_XY_UM,
                    tol_z_um: float = DEFAULT_TOL_Z_UM, surface_tol_um: float = DEFAULT_SURFACE_TOL_UM,
                    overlap_frac: float = DEFAULT_OVERLAP_FRAC, traced_rule: str = "auto",
                    max_candidates: int | None = DEFAULT_MAX_CANDIDATES,
                    max_dense_cells: int = DEFAULT_MAX_DENSE_CELLS,
                    truth_mask_min_voxels: int = TRUTH_MASK_MIN_VOXELS,
                    detector_only_spots: pd.DataFrame | None = None, return_truth_masks: bool = False):
    """Score one gonad. ``labels`` is our nucleus label image (ZYX), ``truth_labels`` the expert's
    (optional). Returns (per-gonad row, per-nucleus table, matched spot pairs), and with
    ``return_truth_masks`` a fourth item: the per-truth-nucleus table of the protocol's section 5.1
    partition (``truth_mask_table`` with ``tag`` and ``split`` columns; empty without ``truth_labels``).
    The row carries that partition's ``tm_`` columns (``truth_mask_row``) whenever ``truth_labels`` is
    given; ``detector_only_spots`` (the spot stage re-run on the truth masks) adds the ``tm_det_`` ones."""
    voxel_zyx = tuple(float(s) for s in voxel_zyx)
    nuclei = outputs.nuclei
    our_spots = outputs.spots if outputs.spots is not None else pd.DataFrame(
        columns=["spot_id", "x_um", "y_um", "z_um"])
    tspots = truth.spots
    tz, ty, tx = _spot_voxels(tspots, voxel_zyx)
    spot_labels = label_at_points(labels, tz, ty, tx)
    ours_n = ours_counts(nuclei, outputs.spots)
    per_nuc, rule = per_nucleus_table(nuclei, spot_labels, ours_n, truth_nuclei=truth.nuclei,
                                      labels=labels, truth_labels=truth_labels,
                                      surface_tol_um=surface_tol_um, overlap_frac=overlap_frac,
                                      traced_rule=traced_rule)
    per_nuc.insert(0, "split", split)
    per_nuc.insert(0, "tag", tag)
    germ_ids = set(per_nuc["nucleus_id"].astype(int).tolist())

    ours_xyz = our_spots[["x_um", "y_um", "z_um"]].to_numpy(float) if len(our_spots) else np.zeros((0, 3))
    truth_xyz = tspots[["x_um", "y_um", "z_um"]].to_numpy(float) if len(tspots) else np.zeros((0, 3))
    a = match_spots(ours_xyz, truth_xyz, tol_xy_um=tol_xy_um, tol_z_um=tol_z_um,
                    max_candidates=max_candidates, max_dense_cells=max_dense_cells)

    row: dict = {"tag": tag, "split": split, "truth_source": truth.source,
                 "voxel_dz_um": voxel_zyx[0], "voxel_dy_um": voxel_zyx[1], "voxel_dx_um": voxel_zyx[2],
                 "n_nuclei_total": int(len(nuclei)), "n_germline_nuclei": int(len(per_nuc)),
                 "traced_rule": rule, "n_traced": int(per_nuc["traced"].sum()),
                 "n_expert_spots": int(len(tspots)),
                 "n_expert_spots_in_our_labels": int((spot_labels > 0).sum()),
                 "n_expert_spots_in_germline_labels": int(np.isin(spot_labels, list(germ_ids)).sum()),
                 "n_our_spots": int(len(our_spots))}
    if "nucleus_id" in our_spots.columns and len(our_spots):
        row["n_our_spots_in_germline"] = int(our_spots["nucleus_id"].isin(germ_ids).sum())
    else:
        row["n_our_spots_in_germline"] = 0

    m_all = agreement_metrics(per_nuc["expert_n"], per_nuc["ours_n"])
    tr = per_nuc[per_nuc["traced"]]
    m_tr = agreement_metrics(tr["expert_n"], tr["ours_n"])
    for k in ("mean_expert", "mean_ours", "pearson_r", "ccc", "bias"):
        row[f"{k}_all"] = m_all[k]
        row[f"{k}_traced"] = m_tr[k]

    row.update(_scores(a.n_matched, len(ours_xyz), len(truth_xyz), ""))
    row["matching_method"] = a.method
    row["n_candidate_pairs"] = int(a.n_candidates)
    row["n_greedy_components"] = int(a.n_greedy_components)
    row["matching_notes"] = "; ".join(a.notes)

    # without Class B: B spots and the spots we matched to them leave both counts
    has_class = truth.has_spot_class
    row["has_spot_class"] = bool(has_class)
    if has_class and len(tspots):
        b = is_class_b(tspots["spot_class"])
        n_b = int(b.sum())
        matched_to_b = int(b[a.cols].sum()) if a.n_matched else 0
        row["n_expert_spots_class_b"] = n_b
        row.update(_scores(a.n_matched - matched_to_b, len(ours_xyz) - matched_to_b,
                           len(truth_xyz) - n_b, "no_b_"))
    else:
        row["n_expert_spots_class_b"] = 0
        row.update(_nan_scores("no_b_"))

    # detector-only: both sets restricted to the traced nuclei of the truth label image
    if truth_labels is not None:
        t_in = label_at_points(truth_labels, tz, ty, tx) > 0
        if len(ours_xyz):
            oi = points_to_voxel_index(ours_xyz, np.asarray(voxel_zyx)[::-1])
            o_in = label_at_points(truth_labels, oi[:, 2], oi[:, 1], oi[:, 0]) > 0
        else:
            o_in = np.zeros(0, bool)
        a_det = match_spots(ours_xyz[o_in], truth_xyz[t_in], tol_xy_um=tol_xy_um, tol_z_um=tol_z_um,
                            max_candidates=max_candidates, max_dense_cells=max_dense_cells)
        row.update(_scores(a_det.n_matched, int(o_in.sum()), int(t_in.sum()), "det_"))
    else:
        row.update(_nan_scores("det_"))
    row["tol_xy_um"], row["tol_z_um"] = tol_xy_um, tol_z_um
    row["surface_tol_um"], row["overlap_frac"] = surface_tol_um, overlap_frac

    # the protocol's section 5.1 primary partition: the truth masks are the nuclei
    tm_table = pd.DataFrame(columns=["tag", "split", *TRUTH_MASK_COLUMNS])
    tm_facts = None
    if truth_labels is not None:
        det_xyz = None
        if detector_only_spots is not None:
            det_xyz = (detector_only_spots[["x_um", "y_um", "z_um"]].to_numpy(float) if len(detector_only_spots)
                       else np.zeros((0, 3)))
        tm_table, tm_facts = truth_mask_table(truth_labels, tspots, ours_xyz, voxel_zyx,
                                              min_voxels=truth_mask_min_voxels, has_class=has_class,
                                              detector_only_xyz=det_xyz)
        tm_table.insert(0, "split", split)
        tm_table.insert(0, "tag", tag)
    row.update(truth_mask_row(tm_table if tm_facts is not None else None, tm_facts))
    if detector_only_spots is not None:
        row["n_our_spots_detector_only"] = int(len(detector_only_spots))

    if a.n_matched:
        d = ours_xyz[a.rows] - truth_xyz[a.cols]
        our_ids = (our_spots["spot_id"].to_numpy() if "spot_id" in our_spots.columns
                   else np.arange(len(our_spots)))
        t_ids = (tspots["spot_id"].to_numpy() if "spot_id" in tspots.columns
                 else np.arange(len(tspots)))
        matches = pd.DataFrame({"tag": tag, "our_spot_id": our_ids[a.rows],
                                "expert_spot_id": t_ids[a.cols], "cost": a.cost,
                                "dxy_um": np.hypot(d[:, 0], d[:, 1]), "dz_um": np.abs(d[:, 2])})
    else:
        matches = pd.DataFrame(columns=["tag", "our_spot_id", "expert_spot_id", "cost", "dxy_um", "dz_um"])
    if return_truth_masks:
        return row, per_nuc, matches, tm_table
    return row, per_nuc, matches


# ------------------------------------------------------------------------------------ pooling
def _detection_block(rows: pd.DataFrame, prefix: str, *, n_boot: int, seed: int) -> dict:
    sub = rows[rows[f"{prefix}n_expert"] + rows[f"{prefix}n_ours"] > 0]
    if sub.empty:
        return {"n_gonads": 0, "note": "no gonad has this variant"}
    payloads = [(int(r[f"{prefix}n_matched"]), int(r[f"{prefix}n_ours"]), int(r[f"{prefix}n_expert"]))
                for _, r in sub.iterrows()]
    boot = cluster_bootstrap(payloads, pooled_detection, n_boot=n_boot, seed=seed)["stats"]
    return {"n_gonads": int(len(sub)), "n_matched": int(sub[f"{prefix}n_matched"].sum()),
            "n_ours": int(sub[f"{prefix}n_ours"].sum()), "n_expert": int(sub[f"{prefix}n_expert"].sum()),
            "precision": boot["precision"], "recall": boot["recall"], "f1": boot["f1"],
            "mean_per_gonad_f1": float(sub[f"{prefix}f1"].mean())}


def _agreement_block(per_nuc: pd.DataFrame, tags, *, n_boot: int, seed: int) -> dict:
    payloads = [(per_nuc.loc[per_nuc["tag"] == t, "expert_n"].to_numpy(float),
                 per_nuc.loc[per_nuc["tag"] == t, "ours_n"].to_numpy(float)) for t in tags]
    payloads = [p for p in payloads if len(p[0])]
    boot = cluster_bootstrap(payloads, pooled_agreement, n_boot=n_boot, seed=seed)["stats"]
    m = agreement_metrics(per_nuc["expert_n"], per_nuc["ours_n"])
    return {"n_nuclei": m["n"], "n_gonads": int(len(payloads)), "mean_expert": m["mean_expert"],
            "mean_ours": m["mean_ours"], "ccc": boot["ccc"], "pearson_r": boot["pearson_r"],
            "bias": boot["bias"], "mae": boot["mae"]}


TRUTH_MASK_BLOCKS = ("per_nucleus_truth_masks", "per_nucleus_truth_masks_class_a",
                     "per_nucleus_truth_masks_class_b", "per_nucleus_truth_masks_detector_only")


def truth_mask_blocks(per_tm: pd.DataFrame | None, tags, *, n_boot: int, seed: int) -> dict:
    """The pooled blocks of the section 5.1 partition from per-truth-nucleus tables (``truth_mask_table``
    rows with a ``tag`` column): the method's count against the truth counted from Class A+B spots
    (``per_nucleus_truth_masks``), from Class A and from Class B spots, and the detector-only count against
    the Class A+B truth when it was measured. Each block is that of the other per-nucleus partitions (n,
    means, and CCC, Pearson r, bias and mean absolute error with their cluster-bootstrap intervals by
    gonad). Empty when there is no table."""
    if per_tm is None or not len(per_tm):
        return {}
    present = set(per_tm["tag"])
    tags = [t for t in tags if t in present]
    out = {"per_nucleus_truth_masks": _agreement_block(per_tm, tags, n_boot=n_boot, seed=seed)}
    for cls in ("class_a", "class_b"):
        col = f"expert_n_{cls}"
        if col in per_tm.columns and per_tm[col].notna().any():
            sub = per_tm[per_tm[col].notna()].copy()
            sub["expert_n"] = sub[col]
            out[f"per_nucleus_truth_masks_{cls}"] = _agreement_block(
                sub, [t for t in tags if t in set(sub["tag"])], n_boot=n_boot, seed=seed)
    col = "ours_n_detector_only"
    if col in per_tm.columns and per_tm[col].notna().any():
        sub = per_tm[per_tm[col].notna()].copy()
        sub["ours_n"] = sub[col]
        out["per_nucleus_truth_masks_detector_only"] = _agreement_block(
            sub, [t for t in tags if t in set(sub["tag"])], n_boot=n_boot, seed=seed)
    return out


def _pool_rows(rows: pd.DataFrame, per_nuc: pd.DataFrame, *, n_boot: int, seed: int,
               per_tm: pd.DataFrame | None = None) -> dict:
    tags = list(rows["tag"])
    traced = per_nuc[per_nuc["traced"]]
    out = {
        "n_gonads": int(len(rows)), "tags": tags,
        "n_expert_spots": int(rows["n_expert_spots"].sum()), "n_our_spots": int(rows["n_our_spots"].sum()),
        "per_nucleus_all": _agreement_block(per_nuc, tags, n_boot=n_boot, seed=seed),
        "per_nucleus_traced": _agreement_block(traced, tags, n_boot=n_boot, seed=seed),
        "spots_all": _detection_block(rows, "", n_boot=n_boot, seed=seed),
        "spots_without_class_b": _detection_block(rows, "no_b_", n_boot=n_boot, seed=seed),
        "spots_detector_only": _detection_block(rows, "det_", n_boot=n_boot, seed=seed),
        "matching_methods": {str(r["tag"]): str(r["matching_method"]) for _, r in rows.iterrows()},
        "truth_sources": {str(r["tag"]): str(r["truth_source"]) for _, r in rows.iterrows()},
    }
    out.update(truth_mask_blocks(per_tm, tags, n_boot=n_boot, seed=seed))
    return out


def pool(rows: pd.DataFrame, per_nuc: pd.DataFrame, *, n_boot: int = 1000, seed: int = 0,
         per_truth_nucleus: pd.DataFrame | None = None) -> dict:
    """Pooled scores over gonads with cluster-bootstrap intervals, overall and by split. With
    ``per_truth_nucleus`` (the concatenated fourth item of ``benchmark_gonad(..., return_truth_masks=True)``)
    the blocks of the section 5.1 partition are added (``truth_mask_blocks``)."""
    per_tm = per_truth_nucleus if per_truth_nucleus is not None and len(per_truth_nucleus) else None
    out = _pool_rows(rows, per_nuc, n_boot=n_boot, seed=seed, per_tm=per_tm)
    out["settings"] = {c: float(rows[c].iloc[0]) for c in ("tol_xy_um", "tol_z_um", "surface_tol_um",
                                                            "overlap_frac") if c in rows.columns and len(rows)}
    out["bootstrap"] = {"n_boot": int(n_boot), "seed": int(seed), "unit": "gonad",
                        "interval": "2.5 and 97.5 percentiles"}
    splits = [s for s in rows["split"].astype(str).unique() if s]
    if len(splits) > 1:
        out["by_split"] = {s: _pool_rows(rows[rows["split"].astype(str) == s],
                                         per_nuc[per_nuc["split"].astype(str) == s],
                                         n_boot=n_boot, seed=seed,
                                         per_tm=(per_tm[per_tm["split"].astype(str) == s] if per_tm is not None
                                                 else None))
                           for s in splits}
    return out
