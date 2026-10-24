"""Scoring sc_straighten against SC truth (docs/SC_VALIDATION_PROTOCOL.md sections 4 to 7).

* Synthetic (T3, check S1): one-to-one matching of our strands to the truth ribbons of a phantom, coverage,
  off-truth length, gaps, foci and the precision of the ``complete`` class (``score_scene``).
* Zones (T1, check S2): the per-gonad zone polyline from the lab's Imaris Measurement Points export, the frame
  rule into the pipeline frame, zone labels and guard bands per nucleus, per-gonad AUCs and the pre-specified
  boundary estimator (``zone_truth``, ``label_nuclei``, ``auc``, ``estimate_boundaries``).

Truth loaders refuse test and reserved gonads unless ``final=True`` (germquant.benchmark.sc_split.guard).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from . import sc_split

TOL_XY, TOL_Z = 0.3, 0.6              # tolerance ellipsoid (um) between truth and our samples
MATCH_MAX_UM = 0.5                    # admissible mean symmetric centreline distance for a matched pair
MATCHABLE_FRAC = 0.8
BOOT_SEED = 20261023


# ------------------------------------------------------------------------------------ polylines
def resample(P: np.ndarray, ds: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    P = np.asarray(P, float)
    if len(P) < 2:
        return P.copy(), np.zeros(len(P))
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
    n = max(2, int(np.ceil(s[-1] / ds)) + 1)
    q = np.linspace(0, s[-1], n)
    return np.stack([np.interp(q, s, P[:, k]) for k in range(3)], 1), q


def _scaled(P):
    P = np.asarray(P, float)
    return np.stack([P[:, 0] / TOL_Z, P[:, 1] / TOL_XY, P[:, 2] / TOL_XY], 1)


def within_tol(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """For each point of A: is a point of B inside the tolerance ellipsoid?"""
    if len(A) == 0:
        return np.zeros(0, bool)
    if len(B) == 0:
        return np.zeros(len(A), bool)
    d, _ = cKDTree(_scaled(B)).query(_scaled(A), distance_upper_bound=1.0)
    return np.isfinite(d)


def seg_lengths(P: np.ndarray) -> np.ndarray:
    """Length attributed to each point (half of each neighbouring segment)."""
    if len(P) < 2:
        return np.zeros(len(P))
    d = np.linalg.norm(np.diff(P, axis=0), axis=1)
    w = np.zeros(len(P))
    w[:-1] += 0.5 * d
    w[1:] += 0.5 * d
    return w


def mean_sym_distance(A: np.ndarray, B: np.ndarray) -> float:
    da, _ = cKDTree(B).query(A)
    db, _ = cKDTree(A).query(B)
    return float(0.5 * (da.mean() + db.mean()))


# ------------------------------------------------------------------------------------ synthetic
@dataclass
class SceneScore:
    family: str
    seed: int
    truth_length: float
    our_length: float
    coverage: float
    off_truth_um: float
    n_strands: int
    n_strands_1um: int
    matches: list = field(default_factory=list)       # per truth ribbon
    strands: list = field(default_factory=list)       # per our strand
    gaps: dict = field(default_factory=dict)
    foci: list = field(default_factory=list)
    extra: dict = field(default_factory=dict)


def strand_paths(samples: pd.DataFrame) -> dict[int, np.ndarray]:
    """In-strand centreline points per strand_id from a samples table."""
    out = {}
    if samples is None or not len(samples):
        return out
    s = samples[samples["in_strand"].astype(bool)]
    for sid, g in s.groupby("strand_id"):
        g = g.sort_values("s_um")
        out[int(sid)] = g[["z_um", "y_um", "x_um"]].to_numpy(float)
    return out


def score_scene(scene, result, *, min_strand_um: float = 1.0) -> SceneScore:
    """Match ``result`` (an ImageResult, or a dict with strands/samples/gaps/foci tables) to the phantom's truth."""
    strands = result.strands if hasattr(result, "strands") else result["strands"]
    samples = result.samples if hasattr(result, "samples") else result["samples"]
    gaps = result.gaps if hasattr(result, "gaps") else result["gaps"]
    foci = result.foci if hasattr(result, "foci") else result["foci"]
    paths = {sid: resample(P)[0] for sid, P in strand_paths(samples).items()}
    srow = strands.set_index("strand_id") if len(strands) else pd.DataFrame()
    truths = [resample(r.points)[0] for r in scene.ribbons]
    tlen = [r.length_um for r in scene.ribbons]
    T_all = np.vstack(truths) if truths else np.zeros((0, 3))
    O_all = np.vstack(list(paths.values())) if paths else np.zeros((0, 3))
    cov = float(sum((seg_lengths(T)[within_tol(T, O_all)]).sum() for T in truths) / max(sum(tlen), 1e-9)) \
        if truths else float("nan")
    off = float(sum((seg_lengths(P)[~within_tol(P, T_all)]).sum() for P in paths.values()))
    our_len = float(strands["length_um"].sum()) if len(strands) else 0.0
    sc = SceneScore(family=scene.family, seed=scene.seed, truth_length=float(sum(tlen)), our_length=our_len,
                    coverage=cov, off_truth_um=off, n_strands=int(len(strands)),
                    n_strands_1um=int((strands["length_um"] >= min_strand_um).sum()) if len(strands) else 0)
    # one-to-one matching (Hungarian) between truth ribbons and our strands of at least 1 um
    from scipy.optimize import linear_sum_assignment

    sids = [sid for sid in paths if sid in srow.index and srow.loc[sid, "length_um"] >= min_strand_um]
    nt, ns = len(truths), len(sids)
    big = 1e6
    cost = np.full((nt, ns), big)
    for i, T in enumerate(truths):
        for j, sid in enumerate(sids):
            P = paths[sid]
            if len(P) < 2:
                continue
            # cheap reject: no overlap at all
            if not within_tol(P[:: max(1, len(P) // 20)], T).any():
                continue
            d = mean_sym_distance(T, P)
            if d <= MATCH_MAX_UM:
                cost[i, j] = d
    pairs = {}
    if nt and ns:
        ri, ci = linear_sum_assignment(cost)
        for i, j in zip(ri, ci):
            if cost[i, j] < big:
                pairs[i] = j
    for i, T in enumerate(truths):
        rec = {"ribbon_id": scene.ribbons[i].ribbon_id, "truth_um": tlen[i], "matched": i in pairs,
               "strand_id": None, "our_um": float("nan"), "err_um": float("nan"), "dist_um": float("nan"),
               "cover": float("nan"), "purity": float("nan"), "matchable": False,
               "n_overlapping_1um": int(sum(1 for sid in sids if within_tol(paths[sid], T).mean() >= 0.5)),
               "tilt_deg": getattr(scene.ribbons[i], "tilt_deg", float("nan"))}
        if i in pairs:
            sid = sids[pairs[i]]
            P = paths[sid]
            cover = float(seg_lengths(T)[within_tol(T, P)].sum() / max(tlen[i], 1e-9))
            pur = float(seg_lengths(P)[within_tol(P, T)].sum() / max(seg_lengths(P).sum(), 1e-9))
            ours = float(srow.loc[sid, "length_um"])
            d, _ = cKDTree(T).query(P)
            rec.update(strand_id=int(sid), our_um=ours, err_um=ours - tlen[i], dist_um=float(np.median(d)),
                       cover=cover, purity=pur, matchable=bool(cover >= MATCHABLE_FRAC and pur >= MATCHABLE_FRAC))
        sc.matches.append(rec)
    matched_strand = {m["strand_id"]: m for m in sc.matches if m["strand_id"] is not None}
    for sid in (srow.index if len(srow) else []):
        r = srow.loc[sid]
        m = matched_strand.get(int(sid))
        sc.strands.append({"strand_id": int(sid), "length_um": float(r["length_um"]),
                           "strand_class": r["strand_class"], "matchable": bool(m and m["matchable"]),
                           "abs_err_um": abs(m["err_um"]) if m else float("nan")})
    # gaps: a detected gap matches a truth gap when their s intervals overlap on the truth ribbon
    tg = [(i, g0, g1) for i, r in enumerate(scene.ribbons) for g0, g1 in r.gaps]
    det = []
    if len(gaps):
        for _, g in gaps.iterrows():
            mid = np.array([g["mid_z_um"], g["mid_y_um"], g["mid_x_um"]], float)
            best = None
            for i, r in enumerate(scene.ribbons):
                T = r.points
                d, k = cKDTree(T).query(mid)
                if best is None or d < best[0]:
                    best = (d, i, float(r.s[k]))
            if best is None:
                continue
            d, i, s_mid = best
            det.append({"ribbon": i, "s0": s_mid - 0.5 * g["length_um"], "s1": s_mid + 0.5 * g["length_um"],
                        "length_um": float(g["length_um"]), "dist_um": d})
    used = set()
    tp, errs = 0, []
    for ti, (i, g0, g1) in enumerate(tg):
        best, bo = None, 0.0
        for di, dg in enumerate(det):
            if di in used or dg["ribbon"] != i or dg["dist_um"] > 0.5:
                continue
            ov = min(g1, dg["s1"]) - max(g0, dg["s0"])
            if ov > bo:
                best, bo = di, ov
        if best is not None:
            used.add(best)
            tp += 1
            errs.append(abs(det[best]["length_um"] - (g1 - g0)))
    sc.gaps = {"n_truth": len(tg), "n_detected": len(det), "tp": tp, "abs_len_err": errs}
    # foci: on the strand matched to the truth ribbon
    if len(scene.truth_foci) and foci is not None and len(foci):
        fz = foci.set_index("spot_id")
        for _, tf in scene.truth_foci.iterrows():
            rid = int(tf["ribbon_id"])
            ti = next(k for k, r in enumerate(scene.ribbons) if r.ribbon_id == rid)
            m = sc.matches[ti]
            row = fz.loc[tf["spot_id"]] if tf["spot_id"] in fz.index else None
            rec = {"spot_id": int(tf["spot_id"]), "assignment": None if row is None else row["assignment"],
                   "correct_ribbon": False, "s_err_um": float("nan")}
            if row is not None and m["strand_id"] is not None and not pd.isna(row["strand_id"]) and \
                    int(row["strand_id"]) == m["strand_id"] and row["assignment"] in ("unique", "ambiguous"):
                rec["correct_ribbon"] = True
                # our s runs from the end with the smaller rounded (z, y, x); align with the truth direction
                T = scene.ribbons[ti].points
                P = paths[m["strand_id"]]
                flip = np.linalg.norm(P[0] - T[0]) > np.linalg.norm(P[0] - T[-1])
                L_ours = float(srow.loc[m["strand_id"], "length_um"])
                s_ours = float(row["s_um"])
                s_al = (L_ours - s_ours) if flip else s_ours
                # truth s measured from the truth start; compare positions through the nearest truth point
                # a focus may sit up to 0.3 um beyond an end (s < 0 or s > L): clamp, never wrap around
                p_ours = P[int(np.clip(int(round(s_ours / 0.05)), 0, len(P) - 1))] if len(P) else None
                if p_ours is not None:
                    _d, k = cKDTree(T).query(p_ours)
                    s_truth_at = float(scene.ribbons[ti].s[k])
                    rec["s_err_um"] = s_truth_at - float(tf["s_um"])
                rec["s_aligned_um"] = s_al
            sc.foci.append(rec)
    return sc


def summarise_scores(scores: list[SceneScore]) -> dict:
    """The S1 metrics over a list of scenes of one family (point estimates)."""
    out: dict = {"n_scenes": len(scores)}
    if not scores:
        return out
    fam = scores[0].family
    m_all = [m for s in scores for m in s.matches]
    if fam in ("iso", "gaps", "foci"):
        # the length error by tilt (descriptive; the end placement over-measures tilted ribbons at real SNR)
        out["tilt_bins"] = {}
        for lo, hi in ((0.0, 30.0), (30.0, 60.0), (60.0, 75.0)):
            e = np.array([m["err_um"] for m in m_all if m["matched"] and np.isfinite(m["tilt_deg"])
                          and (lo < m["tilt_deg"] <= hi or (lo == 0.0 and m["tilt_deg"] == 0.0))], float)
            out["tilt_bins"][f"{int(lo)}-{int(hi)}"] = {
                "n": int(len(e)), "len_err_mean": float(np.mean(e)) if len(e) else float("nan"),
                "len_err_sd": float(np.std(e, ddof=1)) if len(e) > 1 else float("nan")}
        flat = [m for m in m_all if not (np.isfinite(m["tilt_deg"]) and m["tilt_deg"] > 60)]
        steep = [m for m in m_all if np.isfinite(m["tilt_deg"]) and m["tilt_deg"] > 60]
        for name, ms in (("tilt_le60", flat), ("tilt_gt60", steep)):
            e = np.array([m["err_um"] for m in ms if m["matched"]], float)
            out[name] = {
                "n": len(ms), "n_matched": int(len(e)),
                "len_err_mean": float(np.mean(e)) if len(e) else float("nan"),
                "len_err_sd": float(np.std(e, ddof=1)) if len(e) > 1 else float("nan"),
                "loa": [float(np.mean(e) - 1.96 * np.std(e, ddof=1)), float(np.mean(e) + 1.96 * np.std(e, ddof=1))]
                if len(e) > 1 else [float("nan")] * 2,
                "median_centreline_dist_um": float(np.nanmedian([m["dist_um"] for m in ms if m["matched"]]))
                if len(e) else float("nan"),
                "one_strand_frac": float(np.mean([m["matched"] and m["n_overlapping_1um"] == 1 for m in ms]))
                if ms else float("nan")}
    if fam in ("dense", "snr", "phantom_gonad"):
        tl = np.array([s.truth_length for s in scores])
        ol = np.array([s.our_length for s in scores])
        rel = (ol - tl) / np.maximum(tl, 1e-9)
        from .stats import lin_ccc

        out["total_rel_err_mean"] = float(np.mean(rel))
        out["total_rel_err_sd"] = float(np.std(rel, ddof=1)) if len(rel) > 1 else float("nan")
        out["total_ccc"] = float(lin_ccc(tl, ol)) if len(rel) > 2 else float("nan")
        out["coverage_mean"] = float(np.mean([s.coverage for s in scores]))
        out["off_truth_frac_mean"] = float(np.mean([s.off_truth_um / max(s.truth_length, 1e-9) for s in scores]))
        out["strands_per_nucleus"] = float(np.mean([s.n_strands for s in scores]))
        out["matchable_frac"] = float(np.mean([m["matchable"] for m in m_all])) if m_all else float("nan")
    if fam == "gaps":
        tp = sum(s.gaps["tp"] for s in scores)
        nt = sum(s.gaps["n_truth"] for s in scores)
        nd = sum(s.gaps["n_detected"] for s in scores)
        e = [x for s in scores for x in s.gaps["abs_len_err"]]
        out["gap_recall"] = tp / nt if nt else float("nan")
        out["gap_precision"] = tp / nd if nd else float("nan")
        out["gap_len_mae"] = float(np.mean(e)) if e else float("nan")
    if fam == "foci":
        f = [x for s in scores for x in s.foci]
        e = np.array([abs(x["s_err_um"]) for x in f if x["correct_ribbon"] and np.isfinite(x["s_err_um"])])
        out["foci_n"] = len(f)
        out["foci_correct_ribbon"] = float(np.mean([x["correct_ribbon"] for x in f])) if f else float("nan")
        out["foci_abs_s_err_median"] = float(np.median(e)) if len(e) else float("nan")
        out["foci_abs_s_err_p90"] = float(np.percentile(e, 90)) if len(e) else float("nan")
    if fam == "noise":
        out["frac_total_under_1um"] = float(np.mean([s.our_length < 1.0 for s in scores]))
    st = [x for s in scores for x in s.strands]
    comp = [x for x in st if x["strand_class"] == "complete"]
    out["n_complete"] = len(comp)
    out["complete_precision"] = float(np.mean([x["matchable"] and x["abs_err_um"] <= 0.4 for x in comp])) \
        if comp else float("nan")
    n_truth = len(m_all)
    out["complete_recall"] = len([x for x in comp if x["matchable"] and x["abs_err_um"] <= 0.4]) / n_truth \
        if n_truth else float("nan")
    return out


# ------------------------------------------------------------------------------------ zones (T1)
ZONES = ("None", "SCfrag", "SCfull", "LPfrag")
LETTERS_KEYS = ("SC_fragstart", "SC_fullstart", "SC_fraglate", "Pachy_end")


class TruthDefect(ValueError):
    """A truth-side defect (rule 4 of the protocol): recorded in the changelog before any metric."""


def read_position_sheet(path: str | Path) -> pd.DataFrame:
    """The ``Position`` sheet of a per-gonad Imaris export: one title row skipped, then Position X, Position Y,
    Position Z and Name, in row order."""
    df = pd.read_excel(path, sheet_name="Position", skiprows=1)
    cols = {c: str(c).strip() for c in df.columns}
    df = df.rename(columns=cols)
    need = ["Position X", "Position Y", "Name"]
    for c in need:
        if c not in df.columns:
            raise TruthDefect(f"{path}: Position sheet has no column {c!r}")
    df = df[df["Name"].notna()].reset_index(drop=True)
    return df


def check_letters(names: list[str]) -> None:
    for k, n in enumerate(names):
        want = _letter(k)
        if str(n).strip() != want:
            raise TruthDefect(f"point {k + 1} is named {n!r}, expected {want!r}")


def _letter(k: int) -> str:
    s = ""
    k += 1
    while k:
        k, r = divmod(k - 1, 26)
        s = chr(65 + r) + s
    return s


def master_row(master_path: str | Path, gonad: str) -> dict:
    """The master table row of ONE gonad (rows are filtered to it before anything else is done)."""
    df = pd.read_excel(master_path, sheet_name="measurement_points")
    df.columns = [str(c).strip() for c in df.columns]
    sel = df[df["Germline"].astype(str).str.strip() == gonad]
    if len(sel) != 1:
        raise TruthDefect(f"master table: {len(sel)} rows for {gonad}")
    r = sel.iloc[0]
    return {k: (str(r[k]).strip() if k in r and pd.notna(r[k]) else None) for k in ("gonad_start", *LETTERS_KEYS)}


@dataclass
class ZoneTruth:
    gonad: str
    points_um: np.ndarray       # (n, 2) pipeline frame (x, y) um
    names: list
    idx: dict                   # SC_fragstart / SC_fullstart / SC_fraglate / Pachy_end -> 1-based point index
    s_points: np.ndarray        # arc length of each point along the xy polyline
    frame: dict = field(default_factory=dict)


def ims_frame(ims_path: str | Path) -> dict:
    """ExtMin and the Imaris voxel (x, y) from the .ims header attributes only (no dataset read)."""
    import h5py

    from germquant.validate.imaris_ims import read_geometry

    with h5py.File(str(ims_path), "r") as f:
        g = read_geometry(f)
    return {"ext_min_xy": [float(g["ext_min_xyz"][0]), float(g["ext_min_xyz"][1])],
            "v_imaris_xy": [float(g["voxel_xyz_um"][0]), float(g["voxel_xyz_um"][1])], "source": "ims header"}


def nd2_frame(nd2_path: str | Path) -> dict:
    """ExtMin xy = the nd2 stage xy of the first frame; Imaris xy voxel = the nd2 voxel rounded to 0.001 um."""
    import nd2

    with nd2.ND2File(str(nd2_path)) as f:
        vs = f.voxel_size()
        pos = f.frame_metadata(0).channels[0].position.stagePositionUm
    return {"ext_min_xy": [float(pos.x), float(pos.y)], "v_imaris_xy": [round(float(vs.x), 3), round(float(vs.y), 3)],
            "source": "nd2 stage rule"}


def to_pipeline_xy(pos_xy: np.ndarray, frame: dict, v_nd2_xy) -> np.ndarray:
    """Continuous index (coord - ExtMin) / v_imaris - 0.5 (Imaris positions are voxel centres), times the nd2 voxel."""
    P = np.asarray(pos_xy, float)
    em = np.asarray(frame["ext_min_xy"], float)
    vi = np.asarray(frame["v_imaris_xy"], float)
    return ((P - em) / vi - 0.5) * np.asarray(v_nd2_xy, float)


def frame_disagreement(frame_a: dict, frame_b: dict, v_nd2_xy, extent_xy) -> float:
    """The largest |difference| (um, x or y) between two frame rules over the image corners: a stage position
    that frame_a maps to each corner of the run's xy extent is mapped again with frame_b (ExtMin and the voxel
    both enter, so a voxel mismatch shows at the far corner)."""
    v = np.asarray(v_nd2_xy, float)
    em_a, vi_a = np.asarray(frame_a["ext_min_xy"], float), np.asarray(frame_a["v_imaris_xy"], float)
    corners = np.array([[0.0, 0.0], [extent_xy[0], 0.0], [0.0, extent_xy[1]], [extent_xy[0], extent_xy[1]]])
    stage = (corners / v + 0.5) * vi_a + em_a
    other = to_pipeline_xy(stage, frame_b, v_nd2_xy)
    return float(np.max(np.abs(other - corners)))


def zone_truth(record: dict, master_path, v_nd2_xy, *, frame: dict | None = None, final: bool = False,
               split: dict | None = None) -> ZoneTruth:
    """Guarded loader: the zone polyline of one gonad in the pipeline frame and its named zone points."""
    sc_split.guard(record["gonad"], final=final, split=split)
    exp = record["zone_export_copy"]["path"]
    df = read_position_sheet(exp)
    names = [str(n).strip() for n in df["Name"]]
    check_letters(names)
    row = master_row(master_path, record["gonad"])
    idx = {}
    for k in LETTERS_KEYS:
        letter = row.get(k)
        if letter is None or letter not in names:
            raise TruthDefect(f"{record['gonad']}: {k} letter {letter!r} not among the points")
        idx[k] = names.index(letter) + 1
    f, u, lt, e = (idx[k] for k in LETTERS_KEYS)
    if not (f <= u <= lt <= e):
        raise TruthDefect(f"{record['gonad']}: zone letters out of order")
    if frame is None:
        raise ValueError("frame required (ims_frame or nd2_frame)")
    pts = to_pipeline_xy(df[["Position X", "Position Y"]].to_numpy(float), frame, v_nd2_xy)
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    return ZoneTruth(gonad=record["gonad"], points_um=pts, names=names, idx=idx, s_points=s, frame=frame)


GUARD_FRAC_REGISTERED = 0.5          # protocol 4.2: half the neighbouring point spacing on each side


def label_nuclei(zt: ZoneTruth, nuclei: pd.DataFrame, *, x_col="centroid_x_um", y_col="centroid_y_um",
                 off_axis_um: float = 20.0, guard_frac: float = GUARD_FRAC_REGISTERED) -> pd.DataFrame:
    """Zone label, segment, guard band and off-axis cut per nucleus (protocol 4.2). ``guard_frac`` is the guard
    band's share of the neighbouring point spacing on each side of a boundary (the registered rule is 0.5; another
    value is a protocol change that needs a changelog entry before it is scored)."""
    from germquant.staging.zones import project_to_polyline

    s, r, L = project_to_polyline(nuclei[[x_col, y_col]].to_numpy(float), zt.points_um)
    cut = float(off_axis_um)
    inside0 = (r <= off_axis_um) & (s >= 0) & (s <= L)
    if inside0.sum() >= 20:
        cut = float(min(off_axis_um, 2.5 * np.median(r[inside0])))
    sp = zt.s_points
    n = len(sp)
    f, u, lt, e = (zt.idx[k] for k in LETTERS_KEYS)
    s_end = sp[e - 1]
    keep = (r <= cut) & (s >= 0) & (s <= s_end)
    b = np.searchsorted(sp, s, side="right")          # 1-based segment index: s(p_b) <= s < s(p_b+1)
    b = np.where(s >= sp[-1], n - 1, b)
    zone = np.full(len(s), None, dtype=object)
    zone[b < f] = "None"
    zone[(b >= f) & (b < u)] = "SCfrag"
    zone[(b >= u) & (b < lt)] = "SCfull"
    zone[(b >= lt) & (b <= e)] = "LPfrag"
    guard = np.zeros(len(s), bool)
    bounds = []

    def seg_zone(bb: int) -> str | None:
        # the zone of segment bb (from point bb to point bb + 1), by the rule above; None outside the zoned
        # segments 1 .. e - 1 (the zones end at Pachy_end)
        if bb < 1 or bb > min(e - 1, n - 1):
            return None
        return "None" if bb < f else ("SCfrag" if bb < u else ("SCfull" if bb < lt else "LPfrag"))

    gf = float(guard_frac)
    for k in sorted({f, u, lt}):
        if not (1 < k < n):
            continue
        # a named point is a boundary when the zones of the segments on both sides of it exist and differ
        # (so SC_fragstart == SC_fullstart gives a None|SCfull boundary)
        before, after = seg_zone(k - 1), seg_zone(k)
        if before is None or after is None or before == after:
            continue
        sk = sp[k - 1]
        lo = sk - gf * (sk - sp[k - 2])
        hi = sk + gf * (sp[k] - sk)
        guard |= (s > lo) & (s < hi)
        bounds.append({"k": k, "s_um": float(sk), "lo": float(lo), "hi": float(hi), "kind": f"{before}|{after}"})
    out = nuclei[["nucleus_id"]].copy()
    out["s_um"], out["r_um"], out["zone"], out["in_analysis"], out["guard"] = s, r, zone, keep, guard
    out["off_axis_cut_um"] = cut
    out.attrs["boundaries"] = bounds
    out.attrs["s_end"] = float(s_end)
    return out


def auc(pos: np.ndarray, neg: np.ndarray) -> float:
    """Mann-Whitney AUC, ties count one half."""
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if not len(pos) or not len(neg):
        return float("nan")
    allv = np.concatenate([pos, neg])
    from scipy.stats import rankdata

    rk = rankdata(allv)
    rp = rk[:len(pos)].sum()
    return float((rp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def gonad_aucs(lab: pd.DataFrame, values: pd.Series, *, min_n: int = 10) -> dict:
    """S2a (SCfull vs SCfrag + LPfrag) and S2b (any SC zone vs None) AUCs outside guard bands."""
    v = values.reindex(lab["nucleus_id"]).to_numpy(float)
    ok = lab["in_analysis"].to_numpy(bool) & ~lab["guard"].to_numpy(bool) & np.isfinite(v)
    z = lab["zone"].to_numpy(object)
    full = ok & (z == "SCfull")
    frag = ok & ((z == "SCfrag") | (z == "LPfrag"))
    anysc = ok & np.isin(z, ["SCfrag", "SCfull", "LPfrag"])
    none = ok & (z == "None")
    out = {"n_full": int(full.sum()), "n_frag": int(frag.sum()), "n_sc": int(anysc.sum()), "n_none": int(none.sum())}
    out["auc_full_vs_frag"] = auc(v[full], v[frag]) if full.sum() >= min_n and frag.sum() >= min_n else float("nan")
    out["auc_sc_vs_none"] = auc(v[anysc], v[none]) if anysc.sum() >= min_n and none.sum() >= min_n else float("nan")
    out["auc_full_vs_scfrag"] = auc(v[full], v[ok & (z == "SCfrag")]) \
        if full.sum() >= min_n and (ok & (z == "SCfrag")).sum() >= min_n else float("nan")
    out["auc_full_vs_lpfrag"] = auc(v[full], v[ok & (z == "LPfrag")]) \
        if full.sum() >= min_n and (ok & (z == "LPfrag")).sum() >= min_n else float("nan")
    return out


def balanced_accuracy(call: np.ndarray, truth: np.ndarray) -> float:
    call = np.asarray(call, bool)
    truth = np.asarray(truth, bool)
    if truth.all() or (~truth).all():
        return float("nan")
    return float(0.5 * (np.mean(call[truth]) + np.mean(~call[~truth])))


def local_fraction(s: np.ndarray, call: np.ndarray, window: float = 5.0, min_n: int = 3):
    order = np.argsort(s)
    ss = s[order]
    cc = call[order].astype(float)
    f = np.full(len(ss), np.nan)
    for j in range(len(ss)):
        m = np.abs(ss - ss[j]) <= window
        if m.sum() >= min_n:
            f[j] = cc[m].mean()
    return ss, f


def _first_run(ss, f, start_ok, cond, span=20.0, min_def=3, s_limit=None):
    for j in range(len(ss)):
        if not start_ok(ss[j]) or not np.isfinite(f[j]) or not cond(f[j]):
            continue
        hi = ss[j] + span if s_limit is None else min(ss[j] + span, s_limit)
        m = (ss >= ss[j]) & (ss <= hi) & np.isfinite(f)
        if m.sum() < min_def:
            continue
        if all(cond(x) for x in f[m]):
            return float(ss[j])
    return None


def estimate_boundaries(lab: pd.DataFrame, nuc: pd.DataFrame, theta_full: float, theta_present: float) -> dict:
    """The pre-specified boundary estimator (protocol 5.2) over the analysis nuclei (guard bands included)."""
    sel = lab["in_analysis"].to_numpy(bool)
    ids = lab["nucleus_id"].to_numpy()[sel]
    s = lab["s_um"].to_numpy(float)[sel]
    n2 = nuc.set_index("nucleus_id")
    full = n2.reindex(ids)["piece_lw_mean_um"].to_numpy(float) >= theta_full
    present = n2.reindex(ids)["total_length_um"].to_numpy(float) >= theta_present
    s_end = float(lab.attrs.get("s_end", np.nanmax(s) if len(s) else 0.0))
    out = {}
    ss, ff = local_fraction(s, full)
    _, fp = local_fraction(s, present)
    fs = _first_run(ss, ff, lambda x: True, lambda v: v >= 0.5)
    out["SC_fullstart"] = fs
    out["SC_fragstart"] = _first_run(ss, fp, lambda x: True, lambda v: v >= 0.5)
    if fs is not None:
        fl = _first_run(ss, ff, lambda x: x > fs + 20.0, lambda v: v < 0.5, s_limit=s_end)
        out["SC_fraglate"] = fl if fl is not None else s_end
    else:
        out["SC_fraglate"] = None
    return out


def boundary_hit(zt: ZoneTruth, key: str, est: float | None) -> tuple[bool, float]:
    k = zt.idx[key]
    sp = zt.s_points
    if est is None:
        return False, float("nan")
    lo = sp[max(k - 2, 0)]
    hi = sp[min(k, len(sp) - 1)]
    return bool(lo <= est <= hi), float(abs(est - sp[k - 1]))


def choose_thresholds(per_gonad: list[tuple[pd.DataFrame, pd.DataFrame]]) -> dict:
    """theta_full and theta_present on development gonads (protocol 5.2): maximise the mean per-gonad balanced
    accuracy of each call outside guard bands; ties go to the smaller value."""
    res = {}
    for name, grid_, col, truth_fn in (
            ("theta_full", np.round(np.arange(1.0, 6.0001, 0.5), 3), "piece_lw_mean_um",
             lambda z: (z == "SCfull"), ),
            ("theta_present", np.round(np.arange(0.5, 10.0001, 0.5), 3), "total_length_um",
             lambda z: np.isin(z, ["SCfrag", "SCfull", "LPfrag"]))):
        best, bv = None, -1.0
        table = []
        for t in grid_:
            accs = []
            for lab, nuc in per_gonad:
                ok = lab["in_analysis"].to_numpy(bool) & ~lab["guard"].to_numpy(bool)
                z = lab["zone"].to_numpy(object)
                if name == "theta_full":
                    ok &= np.isin(z, ["SCfull", "SCfrag", "LPfrag"])
                else:
                    ok &= np.isin(z, ["None", "SCfrag", "SCfull", "LPfrag"])
                v = nuc.set_index("nucleus_id").reindex(lab["nucleus_id"])[col].to_numpy(float)
                ok &= np.isfinite(v)
                if ok.sum() < 2:
                    continue
                accs.append(balanced_accuracy(v[ok] >= t, truth_fn(z[ok])))
            m = float(np.nanmean(accs)) if accs and np.isfinite(accs).any() else float("nan")
            table.append((float(t), m))
            if np.isfinite(m) and m > bv + 1e-12:
                best, bv = float(t), m
        res[name] = best
        res[name + "_table"] = table
    return res


def nan_mean(xs) -> float:
    xs = [x for x in xs if x is not None and math.isfinite(x)]
    return float(np.mean(xs)) if xs else float("nan")


# ------------------------------------------------------------------------------------ live filaments (T2, S3)
def _is_registered_filament_file(path) -> bool:
    try:
        ft = sc_split.load_split()["filament_truth"]
    except Exception:  # noqa: BLE001
        return True                                  # when unsure, treat it as the sealed file
    names = {Path(ft["ims"]["path"]).name, Path(str(ft.get("nas_copy", ""))).name}
    return Path(str(path)).name in names


def _cluster(cents: np.ndarray, link_um: float = 3.0) -> np.ndarray:
    """Single-linkage groups of path centroids (labels 0..k-1, ordered by first appearance)."""
    n = len(cents)
    lab = -np.ones(n, int)
    k = 0
    for i in range(n):
        if lab[i] >= 0:
            continue
        stack = [i]
        lab[i] = k
        while stack:
            j = stack.pop()
            d = np.linalg.norm(cents - cents[j], axis=1)
            for m in np.nonzero((d <= link_um) & (lab < 0))[0]:
                lab[m] = k
                stack.append(int(m))
        k += 1
    return lab


def filament_check(ims_path, tif_path, voxel_zyx, *, time_index: int | None = None, params=None,
                   final: bool = False) -> dict:
    """Check S3 of the protocol (4.3, 6.3): our strands against the live Filament truth of one traced time point.
    The registered live file is test truth: refused unless ``final=True``."""
    import tifffile

    from germquant.sc.straighten import imports
    from germquant.sc.straighten import params as PRM
    from germquant.sc.straighten.pipeline import straighten_image

    if _is_registered_filament_file(ims_path) and not final:
        raise sc_split.SealedError("the live Filament file is test truth (protocol 2.2, 3.1 rule 6); scored only "
                                   "with final=True by the orchestrator")
    P = params or PRM.resolve()
    sp = tuple(float(v) for v in voxel_zyx)
    paths = [p for p in imports.read_ims_filaments(ims_path, sp, time_index=time_index) if len(p.points) >= 2]
    if not paths:
        raise TruthDefect(f"{ims_path}: no filament with two or more vertices at time {time_index}")
    img = np.asarray(tifffile.imread(str(tif_path)))
    img = img.reshape((-1,) + img.shape[-2:])
    cents = np.array([p.points.mean(0) for p in paths])
    grp = _cluster(cents)
    shape = img.shape
    zz, yy, xx = np.meshgrid(*[np.arange(n) * s for n, s in zip(shape, sp)], indexing="ij")
    labels = np.zeros(shape, np.uint16)
    centres, radii = [], []
    for g in range(grp.max() + 1):
        V = np.vstack([p.points for p, gg in zip(paths, grp) if gg == g])
        c = V.mean(0)
        centres.append(c)
        radii.append(float(np.linalg.norm(V - c, axis=1).max() + 0.5))
    best = np.full(shape, np.inf)
    for g, (c, r) in enumerate(zip(centres, radii), start=1):
        d = np.sqrt((zz - c[0]) ** 2 + (yy - c[1]) ** 2 + (xx - c[2]) ** 2)
        take = (d <= r) & (d < best)
        labels[take] = g
        best = np.where(take, d, best)
    nuc = pd.DataFrame({"nucleus_id": np.arange(1, len(centres) + 1), "germ_cell": "spermatocyte",
                        "in_germline": True,
                        "volume_um3": [float((labels == g).sum() * np.prod(sp)) for g in range(1, len(centres) + 1)]})
    ours = straighten_image(img, labels, sp, nuc, params=P, image_id="s3")
    imp = {}
    for p, g in zip(paths, grp):
        imp.setdefault(int(g) + 1, []).append({"points": p.points, "verified": True, "id": p.path_id})
    drawn = straighten_image(img, labels, sp, nuc, params=P, imported=imp, geometry="import", refine_imported=False,
                             ends_mode="as_drawn", image_id="s3_drawn")
    refined = straighten_image(img, labels, sp, nuc, params=P, imported=imp, geometry="import",
                               refine_imported=True, ends_mode="fit", image_id="s3_refined")
    out = {"n_truth_paths": len(paths), "n_nuclei": len(centres), "nuclei": [], "ribbons": []}
    # S3a: as drawn against the stored Imaris lengths
    dl = drawn.strands.set_index("import_id")["length_um"] if len(drawn.strands) else pd.Series(dtype=float)
    rl = refined.strands.set_index("import_id")["length_um"] if len(refined.strands) else pd.Series(dtype=float)
    ours_paths = {sid: resample(P_)[0] for sid, P_ in strand_paths(ours.samples).items()}
    srow = ours.strands.set_index("strand_id") if len(ours.strands) else pd.DataFrame()
    for g in range(1, len(centres) + 1):
        tr = [(p, resample(p.points)[0]) for p, gg in zip(paths, grp) if gg + 1 == g]
        mine = [sid for sid in ours_paths if sid in srow.index and int(srow.loc[sid, "nucleus_id"]) == g]
        Q_ours = np.vstack([ours_paths[s] for s in mine]) if mine else np.zeros((0, 3))
        T_all = np.vstack([t for _p, t in tr])
        rec_len = sum(seg_lengths(t)[within_tol(t, Q_ours)].sum() for _p, t in tr)
        tot_truth = sum(seg_lengths(t).sum() for _p, t in tr)
        prec_len = sum(seg_lengths(ours_paths[s])[within_tol(ours_paths[s], T_all)].sum() for s in mine)
        tot_ours = sum(seg_lengths(ours_paths[s]).sum() for s in mine)
        refined_total = float(sum(rl.get(p.path_id, np.nan) for p, _t in tr))
        out["nuclei"].append({"nucleus": g, "recall": rec_len / max(tot_truth, 1e-9),
                              "precision": prec_len / max(tot_ours, 1e-9), "our_total_um": float(
                                  srow.loc[mine, "length_um"].sum()) if mine else 0.0,
                              "refined_truth_total_um": refined_total,
                              "total_rel_err": (float(srow.loc[mine, "length_um"].sum()) - refined_total) /
                              max(refined_total, 1e-9) if mine else -1.0})
        for p, t in tr:
            best, bd = None, np.inf
            for s in mine:
                if srow.loc[s, "length_um"] < 1.0:
                    continue
                d = mean_sym_distance(t, ours_paths[s])
                if d < bd:
                    best, bd = s, d
            row = {"path_id": p.path_id, "nucleus": g, "stored_um": p.stored_length_um,
                   "drawn_um": float(dl.get(p.path_id, np.nan)), "refined_um": float(rl.get(p.path_id, np.nan)),
                   "matched_strand": best if bd <= MATCH_MAX_UM else None, "matchable": False,
                   "our_um": np.nan}
            if row["matched_strand"] is not None:
                P_ = ours_paths[best]
                cover = float(seg_lengths(t)[within_tol(t, P_)].sum() / max(seg_lengths(t).sum(), 1e-9))
                pur = float(seg_lengths(P_)[within_tol(P_, t)].sum() / max(seg_lengths(P_).sum(), 1e-9))
                row.update(matchable=bool(cover >= MATCHABLE_FRAC and pur >= MATCHABLE_FRAC),
                           our_um=float(srow.loc[best, "length_um"]))
            out["ribbons"].append(row)
    return out


# ------------------------------------------------------------------------------------ S1 criteria (protocol 8.1)
def evaluate_s1(fams: dict) -> list[dict]:
    """The synthetic criteria of protocol 8.1 on the per-family summaries (summarise_scores). Each row: id,
    criterion, value, pass (None when the family is missing)."""
    rows = []

    def add(cid, text, value, ok):
        rows.append({"id": cid, "criterion": text, "value": value, "pass": None if ok is None else bool(ok)})

    iso = fams.get("iso", {}).get("tilt_le60")
    for b_, v_ in (fams.get("iso", {}).get("tilt_bins") or {}).items():
        add("S1a", f"iso mean length error, tilt {b_} deg (descriptive, n {v_['n']})", v_["len_err_mean"], None)
    if iso:
        loa = iso.get("loa", [float("nan")] * 2)
        add("S1a", "iso mean length error within +/-0.10 um", iso["len_err_mean"], abs(iso["len_err_mean"]) <= 0.10)
        add("S1a", "iso 95% limits within +/-0.40 um", loa, max(abs(loa[0]), abs(loa[1])) <= 0.40)
        add("S1a", "iso median centreline distance <= 0.05 um", iso["median_centreline_dist_um"],
            iso["median_centreline_dist_um"] <= 0.05)
        add("S1a", "iso one strand per ribbon >= 0.95", iso["one_strand_frac"], iso["one_strand_frac"] >= 0.95)
    d = fams.get("dense")
    if d:
        add("S1b", "dense mean total relative error within +/-5%", d["total_rel_err_mean"],
            abs(d["total_rel_err_mean"]) <= 0.05)
        add("S1b", "dense SD of the relative error <= 7%", d["total_rel_err_sd"], d["total_rel_err_sd"] <= 0.07)
        add("S1b", "dense CCC of totals >= 0.80", d["total_ccc"], d["total_ccc"] >= 0.80)
        add("S1b", "dense coverage >= 0.95", d["coverage_mean"], d["coverage_mean"] >= 0.95)
        add("S1b", "dense off-truth <= 0.10", d["off_truth_frac_mean"], d["off_truth_frac_mean"] <= 0.10)
    g = fams.get("gaps")
    if g:
        add("S1c", "gap precision >= 0.80", g["gap_precision"], g["gap_precision"] >= 0.80)
        add("S1c", "gap recall >= 0.80", g["gap_recall"], g["gap_recall"] >= 0.80)
        add("S1c", "gap length error <= 0.15 um", g["gap_len_mae"], g["gap_len_mae"] <= 0.15)
    f = fams.get("foci")
    if f:
        add("S1d", "foci 90th percentile |s error| <= 0.15 um", f["foci_abs_s_err_p90"], f["foci_abs_s_err_p90"] <= 0.15)
        add("S1d", "foci on the correct ribbon >= 0.95", f["foci_correct_ribbon"], f["foci_correct_ribbon"] >= 0.95)
    num = den = 0.0
    for name, v in fams.items():
        if isinstance(v, dict) and v.get("n_complete") and v.get("complete_precision") == v.get("complete_precision"):
            num += v["complete_precision"] * v["n_complete"]
            den += v["n_complete"]
    if den:
        add("S1e", "complete precision (pooled over families) >= 0.90", num / den, num / den >= 0.90)
    for name, v in sorted(fams.items()):
        if isinstance(v, dict) and v.get("n_complete"):
            add("S1e", f"{name} family: complete precision (descriptive, {v['n_complete']} complete strands)",
                v.get("complete_precision"), None)
    n = fams.get("noise")
    if n:
        add("S1f", "noise-only nuclei with a total under 1 um >= 0.95", n["frac_total_under_1um"],
            n["frac_total_under_1um"] >= 0.95)
    if "snr_floor_counts" in fams:
        add("S1g", "SNR floor (lowest scene peak with the dense mean total within +/-5%; descriptive)",
            fams["snr_floor_counts"], None)
    return rows
