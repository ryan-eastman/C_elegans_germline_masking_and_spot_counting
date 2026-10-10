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
                    ) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """Score one gonad. ``labels`` is our nucleus label image (ZYX), ``truth_labels`` the expert's
    (optional). Returns (per-gonad row, per-nucleus table, matched spot pairs)."""
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
            "bias": boot["bias"]}


def _pool_rows(rows: pd.DataFrame, per_nuc: pd.DataFrame, *, n_boot: int, seed: int) -> dict:
    tags = list(rows["tag"])
    traced = per_nuc[per_nuc["traced"]]
    return {
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


def pool(rows: pd.DataFrame, per_nuc: pd.DataFrame, *, n_boot: int = 1000, seed: int = 0) -> dict:
    """Pooled scores over gonads with cluster-bootstrap intervals, overall and by split."""
    out = _pool_rows(rows, per_nuc, n_boot=n_boot, seed=seed)
    out["settings"] = {c: float(rows[c].iloc[0]) for c in ("tol_xy_um", "tol_z_um", "surface_tol_um",
                                                            "overlap_frac") if c in rows.columns and len(rows)}
    out["bootstrap"] = {"n_boot": int(n_boot), "seed": int(seed), "unit": "gonad",
                        "interval": "2.5 and 97.5 percentiles"}
    splits = [s for s in rows["split"].astype(str).unique() if s]
    if len(splits) > 1:
        out["by_split"] = {s: _pool_rows(rows[rows["split"].astype(str) == s],
                                         per_nuc[per_nuc["split"].astype(str) == s],
                                         n_boot=n_boot, seed=seed) for s in splits}
    return out
