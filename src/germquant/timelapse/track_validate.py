"""Biologist validation of the linking (docs/DESIGN_timelapse.md section 7, docs/TIMELAPSE.md "Validating the
linking").

``pick_validation_nuclei`` (scripts/pick_validation_nuclei.py) chooses, seeded and stratified, the nuclei of
the first timepoint the biologist follows by hand in Imaris: by default 16 per series, 6 in the transition
zone, 7 in pachytene and 3 at the edge of the field, half marked ``holdout`` (never looked at while tuning)
and half ``constraint`` (used to tune). The worksheet ``<stem>__validation_picks.csv`` holds the t0 point of
every pick (hand_track_id, t_index, z_um, y_um, x_um in raw-frame microns); the biologist appends one row per
later timepoint.

``validate_tracks`` (``germquant track-validate``) looks every hand point up in that timepoint's label image.
A point outside every label is a segmentation miss, not a linking error. The in-label points of a hand track
are mapped to ``time_track_id`` through the track_points table; consecutive in-label points form the hand
edges. Per hand track: ``segmentation_recall`` (points in a label / points), ``linking_recall`` (hand edges the
pipeline kept in one track / hand edges, so a bridged gap counts and a fragmented track does not),
``edge_precision`` (pipeline edges leaving a hand point that land on the hand nucleus / pipeline edges whose
target timepoint has a hand point), ``switches`` (hand edges where the pipeline's identity jumped to or from
another nucleus) and ``edges_missing`` (hand edges the pipeline simply broke). The verdict compares the holdout
switches, the linking recall and the (pachytene, when zones are known) segmentation recall with the
thresholds ``timelapse.validate.max_switches_holdout`` (0), ``min_linking_recall`` (0.95) and
``min_segmentation_recall`` (0.85), every one read with its code default. ``subset="constraint"``
(``germquant track-validate --set constraint``) scores the constraint tracks alone, which is how linking settings
are compared without looking at the holdout tracks; switches and linking recall are then judged on them.
"""
from __future__ import annotations

import json
import logging
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

from .. import provenance, schema
from ..fsutil import long_path
from .layout import SeriesLayout, SeriesLayoutError, _truthy, find_series
from .series import write_series_table
from .track import cfg_get, lookup_labels, read_track_points

log = logging.getLogger(__name__)

DEFAULTS = {
    "n_nuclei": 16, "n_transition": 6, "n_pachytene": 7, "n_edge": 3, "seed": 0,
    "edge_um": 5.0,              # a nucleus within this distance of the xy border of the field is an edge candidate
    "tz_quantile": 0.35,         # axis_position_um quantile below which a nucleus counts as transition zone (no staging)
    "pach_quantile": 0.90,       # ... and up to which it counts as pachytene
    "max_switches_holdout": 0,
    "min_linking_recall": 0.95,
    "min_segmentation_recall": 0.85,
}
PICKS_TABLE = "validation_picks"            # <stem>__validation_picks.csv (the biologist's worksheet)
VALIDATION_RECORD = "track_validation"      # <stem>__track_validation.json next to the table
ROLES = ("holdout", "constraint")
STAGING_TZ = ("pre",)                       # staging zones counted as transition zone when the staging stage ran
STAGING_PACH = ("early", "mid", "late")


# =====================================================================================================
# the stratified pick
# =====================================================================================================

def pick_validation_nuclei(series_dir: str | Path, *, n: int | None = None, n_tz: int | None = None,
                           n_pach: int | None = None, n_edge: int | None = None, seed: int | None = None,
                           edge_um: float | None = None, stem: str | None = None, cfg=None,
                           write: bool = True, out: str | Path | None = None) -> pd.DataFrame:
    """The worksheet (module docstring): one row per pick at the first timepoint with ``hand_track_id``,
    ``role``, ``zone`` (transition_zone | pachytene | edge | fill), ``t_index``, ``z_um``, ``y_um``, ``x_um``
    (raw frame), ``nucleus_id``, ``image_id`` and ``axis_position_um``. Zones come from the staging stage's
    ``zone`` column when it ran, else from axis_position_um quantiles (centroid_x_um when there is no axis);
    edge nuclei lie within ``edge_um`` of the xy border of the field (the most proximal nuclei fill in when
    there are too few). The same seed gives the same picks."""
    n = int(n if n is not None else cfg_get(cfg, "timelapse.validate.n_nuclei", DEFAULTS["n_nuclei"]))
    n_tz = int(n_tz if n_tz is not None else cfg_get(cfg, "timelapse.validate.n_transition", DEFAULTS["n_transition"]))
    n_pach = int(n_pach if n_pach is not None else cfg_get(cfg, "timelapse.validate.n_pachytene", DEFAULTS["n_pachytene"]))
    n_edge = int(n_edge if n_edge is not None else cfg_get(cfg, "timelapse.validate.n_edge", DEFAULTS["n_edge"]))
    seed = int(seed if seed is not None else cfg_get(cfg, "timelapse.validate.seed", DEFAULTS["seed"]))
    edge_um = float(edge_um if edge_um is not None else cfg_get(cfg, "timelapse.validate.edge_um", DEFAULTS["edge_um"]))
    if n_tz + n_pach + n_edge != n:
        raise ValueError(f"the strata must add up to n: {n_tz} + {n_pach} + {n_edge} != {n}")
    layout = find_series(series_dir, stem)
    t0 = layout.timepoints[0]
    df = layout.read_nuclei(t0)
    if "in_germline" in df.columns:
        df = df[df["in_germline"].map(_truthy)]
    if df.empty:
        raise SeriesLayoutError(f"{layout.stem}: no germline nucleus at t{t0:03d} to pick from")
    df = df.sort_values("nucleus_id").reset_index(drop=True)
    pos = pd.to_numeric(df.get("axis_position_um", pd.Series(np.nan, index=df.index)), errors="coerce")
    if pos.notna().sum() < max(3, len(df) // 2):
        log.warning("%s: axis_position_um is missing at t%03d; centroid_x_um orders the nuclei instead", layout.stem, t0)
        pos = pd.to_numeric(df["centroid_x_um"], errors="coerce")
    rank = pos.rank(method="average", pct=True)
    tzq = float(cfg_get(cfg, "timelapse.validate.tz_quantile", DEFAULTS["tz_quantile"]))
    pq = float(cfg_get(cfg, "timelapse.validate.pach_quantile", DEFAULTS["pach_quantile"]))
    zones = df.get("zone")
    staged = zones is not None and zones.astype(str).str.lower().isin(STAGING_TZ + STAGING_PACH).any()
    if staged:
        z = zones.astype(str).str.lower()
        is_tz, is_pach = z.isin(STAGING_TZ), z.isin(STAGING_PACH)
        zone_source = "staging zone column"
    else:
        is_tz, is_pach = rank <= tzq, (rank > tzq) & (rank <= pq)
        zone_source = f"axis_position_um quantiles (transition zone <= {tzq:g}, pachytene <= {pq:g})"
    spacing = layout.spacing()
    shape = layout.labels_shape(t0)
    yx = df[["centroid_y_um", "centroid_x_um"]].to_numpy(float)
    extent = np.array([shape[1] * spacing[1], shape[2] * spacing[2]])
    dist_border = np.minimum(yx, extent - yx).min(axis=1)
    is_edge = pd.Series(dist_border <= edge_um, index=df.index)

    rng = np.random.default_rng(seed)
    taken: set[int] = set()
    picks: list[tuple[int, str]] = []

    def take(mask: pd.Series, k: int, zone: str, order=None):
        cands = df.loc[mask & ~df["nucleus_id"].isin(taken), "nucleus_id"].astype(int).tolist()
        if order is not None:
            cands = [c for c in order if c in set(cands)]
            chosen = cands[:k]
        else:
            chosen = sorted(rng.choice(cands, size=min(k, len(cands)), replace=False).tolist()) if cands else []
        for c in chosen:
            taken.add(int(c))
            picks.append((int(c), zone))
        return k - len(chosen)

    short = take(is_edge, n_edge, "edge")
    if short:                       # too few nuclei at the border: the most proximal ones stand in
        by_pos = df.assign(_r=rank).sort_values("_r", ascending=False)["nucleus_id"].astype(int).tolist()
        short = take(pd.Series(True, index=df.index), short, "edge", order=by_pos)
    short += take(is_tz, n_tz, "transition_zone")
    short += take(is_pach, n_pach, "pachytene")
    if short:
        short = take(pd.Series(True, index=df.index), short, "fill")
        if short:
            log.warning("%s: only %d germline nuclei at t%03d; %d pick(s) short", layout.stem, len(df), t0, short)
    roles = [ROLES[i % 2] for i in range(len(picks))]
    rows = []
    by_id = df.set_index("nucleus_id")
    pos_by_id = pos.set_axis(df["nucleus_id"].astype(int))
    for i, ((nid, zone), role) in enumerate(zip(picks, roles), 1):
        r = by_id.loc[nid]
        ap = float(pos_by_id.loc[nid])
        rows.append({"hand_track_id": i, "role": role, "zone": zone, "t_index": int(t0),
                     "z_um": float(r["centroid_z_um"]), "y_um": float(r["centroid_y_um"]), "x_um": float(r["centroid_x_um"]),
                     "nucleus_id": int(nid), "image_id": layout.image_id(t0),
                     "axis_position_um": ap if np.isfinite(ap) else float("nan")})
    out_df = pd.DataFrame(rows, columns=["hand_track_id", "role", "zone", "t_index", "z_um", "y_um", "x_um",
                                         "nucleus_id", "image_id", "axis_position_um"])
    out_df.attrs.update({"seed": seed, "zone_source": zone_source, "edge_um": edge_um, "series_id": layout.stem})
    if write:
        p = Path(out) if out else layout.series_table_path(PICKS_TABLE)
        os.makedirs(long_path(p.parent), exist_ok=True)
        with open(long_path(p), "w", encoding="utf-8", newline="") as fh:
            fh.write(f"# {layout.stem}: {len(out_df)} nuclei picked at t{t0:03d} with seed {seed}; zones from {zone_source}; "
                     f"edge = within {edge_um:g} um of the xy border. Follow each nucleus in Imaris and append one row per "
                     f"timepoint (hand_track_id, t_index, z_um, y_um, x_um in raw-frame microns); leave the other columns empty.\n")
            out_df.to_csv(fh, index=False)
        log.info("%s: %d validation nuclei picked at t%03d (%s) -> %s", layout.stem, len(out_df), t0, zone_source, p)
    return out_df


def read_hand_points(path: str | Path) -> pd.DataFrame:
    """The biologist's completed worksheet (or any table with hand_track_id, t_index, z_um, y_um, x_um; a
    leading ``#`` comment line is skipped; ``role`` / ``validation_set`` and ``zone`` ride along when present).
    Rows without coordinates (the blank rows of the worksheet) are dropped."""
    df = pd.read_csv(long_path(Path(path)), comment="#", low_memory=False)
    low = {str(c).strip().lower(): c for c in df.columns}

    def col(*names):
        for n in names:
            if n in low:
                return low[n]
        return None

    c_id = col("hand_track_id", "id", "track_id")
    c_t = col("t_index", "t", "frame", "time")
    c_z, c_y, c_x = col("z_um", "z"), col("y_um", "y"), col("x_um", "x")
    if None in (c_id, c_t, c_z, c_y, c_x):
        raise ValueError(f"{Path(path).name}: needs hand_track_id, t_index, z_um, y_um, x_um columns (got {list(df.columns)})")
    out = pd.DataFrame({"hand_track_id": pd.to_numeric(df[c_id], errors="coerce"),
                        "t_index": pd.to_numeric(df[c_t], errors="coerce"),
                        "z_um": pd.to_numeric(df[c_z], errors="coerce"), "y_um": pd.to_numeric(df[c_y], errors="coerce"),
                        "x_um": pd.to_numeric(df[c_x], errors="coerce")})
    c_role, c_zone = col("role", "validation_set"), col("zone")
    out["role"] = df[c_role].astype(str).str.strip().str.lower().where(df[c_role].notna(), None) if c_role else None
    out["zone"] = df[c_zone].astype(str).str.strip().where(df[c_zone].notna(), None) if c_zone else None
    out = out.dropna(subset=["hand_track_id", "t_index", "z_um", "y_um", "x_um"]).copy()
    out["hand_track_id"] = out["hand_track_id"].astype(int)
    out["t_index"] = out["t_index"].round().astype(int)
    return out.sort_values(["hand_track_id", "t_index"]).reset_index(drop=True)


# =====================================================================================================
# track-validate
# =====================================================================================================

def validate_tracks(series_dir: str | Path, hand_csv: str | Path, *, stem: str | None = None, cfg=None,
                    snap_um: float = 0.0, write: bool = True, subset: str | None = None) -> dict:
    """Score the series' track_points against the biologist's hand-followed nuclei (module docstring). Returns
    ``points`` (one row per hand point: the TRACK_VALIDATION columns with the per-track scores repeated),
    ``tracks`` (one row per hand track), ``verdict`` (the thresholds, the measured values and PASS / FAIL per
    criterion) and the written ``paths``.

    ``subset`` ``constraint`` or ``holdout`` scores only the hand tracks of that role (``None`` or ``all``: every
    track). Use ``constraint`` to compare linking settings, so the holdout tracks stay unseen until the chosen
    setting is judged on them once. A subset run writes only ``<stem>__track_validation_<subset>.json``; the
    track_validation table is written by a run over every track."""
    sub = None if subset in (None, "", "all") else str(subset).strip().lower()
    if sub is not None and sub not in ROLES:
        raise ValueError(f"subset must be one of all, {', '.join(ROLES)}; got {subset!r}")
    layout = find_series(series_dir, stem)
    spacing = layout.spacing()
    hand = read_hand_points(hand_csv)
    if hand.empty:
        raise ValueError(f"{Path(hand_csv).name}: no hand point with coordinates")
    # the role of every hand track: the worksheet's column, else the picks file, else holdout (the strict reading)
    roles = _roles(hand, layout)
    zones = _zones(hand, layout)
    if sub is not None:
        hand = hand[hand["hand_track_id"].map(roles) == sub].reset_index(drop=True)
        if hand.empty:
            raise ValueError(f"{Path(hand_csv).name}: no {sub} hand track")
    off = ~hand["t_index"].isin(layout.timepoints)
    if off.any():
        log.warning("%s: %d hand point(s) at timepoints not in the series dropped: %s", layout.stem, int(off.sum()),
                    sorted(hand.loc[off, "t_index"].unique().tolist()))
        hand = hand[~off].reset_index(drop=True)
    hand["nucleus_id"] = lookup_labels(layout, hand, spacing, radius_um=snap_um)

    tp = read_track_points(layout)
    tp = tp[tp["time_track_id"].notna()]
    tid_of = {(int(r.t_index), int(r.nucleus_id)): int(r.time_track_id) for r in tp.itertuples(index=False)}
    track_at: dict[int, dict[int, int]] = {}                 # time_track_id -> {t_index: nucleus_id}
    for r in tp.itertuples(index=False):
        track_at.setdefault(int(r.time_track_id), {})[int(r.t_index)] = int(r.nucleus_id)

    hand["time_track_id"] = pd.array([tid_of.get((int(t), int(n))) if pd.notna(n) else None
                                      for t, n in zip(hand["t_index"], hand["nucleus_id"])], dtype="Int64")
    hand["image_id"] = hand["t_index"].map(layout.image_id)
    hand["t_seconds"] = hand["t_index"].map(layout.t_seconds_of)
    hand["outcome"] = ""
    per_track = []
    for hid, g in hand.groupby("hand_track_id", sort=True):
        g = g.sort_values("t_index")
        idx = g.index.tolist()
        inlab = [i for i in idx if pd.notna(hand.at[i, "nucleus_id"])]
        for i in idx:
            if i not in inlab:
                hand.at[i, "outcome"] = "segmentation_miss"
        n_edges = n_repro = switches = missing = 0
        n_pe = n_pe_ok = 0
        for k, i in enumerate(inlab):
            t, nid = int(hand.at[i, "t_index"]), int(hand.at[i, "nucleus_id"])
            tid = hand.at[i, "time_track_id"]
            tid = None if pd.isna(tid) else int(tid)
            if k == 0:
                hand.at[i, "outcome"] = "linked" if tid is not None else "unlinked"
            else:
                j = inlab[k - 1]
                t_prev, nid_prev = int(hand.at[j, "t_index"]), int(hand.at[j, "nucleus_id"])
                tid_prev = hand.at[j, "time_track_id"]
                tid_prev = None if pd.isna(tid_prev) else int(tid_prev)
                n_edges += 1
                if tid is not None and tid_prev is not None and tid == tid_prev:
                    n_repro += 1
                    hand.at[i, "outcome"] = "linked"
                else:
                    jumped = ((tid_prev is not None and track_at.get(tid_prev, {}).get(t) not in (None, nid))
                              or (tid is not None and track_at.get(tid, {}).get(t_prev) not in (None, nid_prev)))
                    if jumped:
                        switches += 1
                        hand.at[i, "outcome"] = "switch"
                    else:
                        missing += 1
                        hand.at[i, "outcome"] = "edge_missing" if tid is not None else "unlinked"
            # the pipeline's edge leaving this hand point: where does the track go next, and does the hand agree?
            if tid is not None:
                later = sorted(tt for tt in track_at.get(tid, {}) if tt > t)
                if later:
                    t_next = later[0]
                    row_next = g[g["t_index"] == t_next]
                    if len(row_next) and pd.notna(row_next["nucleus_id"].iloc[0]):
                        n_pe += 1
                        n_pe_ok += int(track_at[tid][t_next] == int(row_next["nucleus_id"].iloc[0]))
        tids = [int(hand.at[i, "time_track_id"]) for i in inlab if pd.notna(hand.at[i, "time_track_id"])]
        dominant = max(sorted(set(tids)), key=tids.count) if tids else None
        per_track.append({
            "hand_track_id": int(hid), "role": roles.get(int(hid), "holdout"), "zone": zones.get(int(hid)),
            "dominant_time_track_id": dominant, "n_points": int(len(idx)), "n_in_label": int(len(inlab)),
            "segmentation_recall": len(inlab) / len(idx) if idx else float("nan"),
            "n_edges": n_edges, "n_edges_reproduced": n_repro,
            "linking_recall": (n_repro / n_edges) if n_edges else float("nan"),
            "n_pipeline_edges": n_pe, "n_pipeline_edges_correct": n_pe_ok,
            "edge_precision": (n_pe_ok / n_pe) if n_pe else float("nan"),
            "switches": switches, "edges_missing": missing,
        })
    tracks = pd.DataFrame(per_track)
    per_point = hand.drop(columns=["role", "zone"]).merge(tracks.drop(columns=["zone"]), on="hand_track_id", how="left")
    per_point["role"] = per_point["hand_track_id"].map(roles)
    per_point["validation_set"] = per_point["role"]
    per_point["zone"] = per_point["hand_track_id"].map(zones)
    per_point = per_point.sort_values(["hand_track_id", "t_index"]).reset_index(drop=True)

    verdict = _verdict(tracks, hand, cfg)
    verdict["subset"] = sub or "all"
    paths = {}
    if write and sub is not None:
        p_rec = layout.out_dir / f"{layout.stem}__{VALIDATION_RECORD}_{sub}.json"
        with open(long_path(p_rec), "w", encoding="utf-8") as fh:
            json.dump({"series_id": layout.stem, "hand_file": str(hand_csv), "subset": sub, "verdict": verdict,
                       "tracks": json.loads(tracks.to_json(orient="records")), "written": provenance.run_timestamp()},
                      fh, indent=1, default=str)
        paths = {"record": str(p_rec)}
    elif write:
        meta = layout.series_meta()
        cols = list(schema.TRACK_VALIDATION) + [c for c in ("edges_missing", "dominant_time_track_id", "n_edges",
                                                          "n_edges_reproduced", "n_pipeline_edges", "n_pipeline_edges_correct")]
        table = per_point[[c for c in cols if c in per_point.columns]]
        p_tab = write_series_table(table, "track_validation", layout.out_dir, layout.stem, meta)
        p_rec = layout.out_dir / f"{layout.stem}__{VALIDATION_RECORD}.json"
        with open(long_path(p_rec), "w", encoding="utf-8") as fh:
            json.dump({"series_id": layout.stem, "hand_file": str(hand_csv), "verdict": verdict,
                       "tracks": json.loads(tracks.to_json(orient="records")), "written": provenance.run_timestamp()},
                      fh, indent=1, default=str)
        paths = {"track_validation": str(p_tab), "record": str(p_rec)}
    return {"points": per_point, "tracks": tracks, "verdict": verdict, "paths": paths, "layout": layout}


def _roles(hand: pd.DataFrame, layout: SeriesLayout) -> dict[int, str]:
    roles: dict[int, str] = {}
    if hand["role"].notna().any():
        for hid, g in hand.groupby("hand_track_id"):
            v = g["role"].dropna()
            if len(v):
                roles[int(hid)] = str(v.iloc[0]).lower()
    picks = layout.series_table_path(PICKS_TABLE)
    if picks.is_file():
        try:
            pk = pd.read_csv(long_path(picks), comment="#", low_memory=False)
            for r in pk.itertuples(index=False):
                roles.setdefault(int(r.hand_track_id), str(r.role).lower())
        except Exception as e:  # noqa: BLE001 - the worksheet is optional
            log.warning("%s: could not read %s (%s)", layout.stem, picks.name, e)
    for hid in hand["hand_track_id"].unique():
        role = roles.get(int(hid))
        if role not in ROLES:
            if role is not None:
                log.warning("%s: hand track %d has role %r; treated as holdout", layout.stem, int(hid), role)
            roles[int(hid)] = "holdout"
    return roles


def _zones(hand: pd.DataFrame, layout: SeriesLayout) -> dict[int, str | None]:
    zones: dict[int, str | None] = {}
    if hand["zone"].notna().any():
        for hid, g in hand.groupby("hand_track_id"):
            v = g["zone"].dropna()
            if len(v):
                zones[int(hid)] = str(v.iloc[0])
    picks = layout.series_table_path(PICKS_TABLE)
    if picks.is_file():
        try:
            pk = pd.read_csv(long_path(picks), comment="#", low_memory=False)
            for r in pk.itertuples(index=False):
                zones.setdefault(int(r.hand_track_id), None if pd.isna(r.zone) else str(r.zone))
        except Exception:  # noqa: BLE001
            pass
    return {int(h): zones.get(int(h)) for h in hand["hand_track_id"].unique()}


def _ratio(num: float, den: float) -> float:
    return (num / den) if den else float("nan")


def _verdict(tracks: pd.DataFrame, hand: pd.DataFrame, cfg) -> dict:
    max_sw = int(cfg_get(cfg, "timelapse.validate.max_switches_holdout", DEFAULTS["max_switches_holdout"]))
    min_lr = float(cfg_get(cfg, "timelapse.validate.min_linking_recall", DEFAULTS["min_linking_recall"]))
    min_sr = float(cfg_get(cfg, "timelapse.validate.min_segmentation_recall", DEFAULTS["min_segmentation_recall"]))
    hold = tracks[tracks["role"] == "holdout"]
    judged = hold if len(hold) else tracks
    if len(hold):
        judged_set = "holdout"
    elif len(tracks) and (tracks["role"] == "constraint").all():
        judged_set = "constraint (no holdout track scored)"
    else:
        judged_set = "all (no holdout track)"
    # the switches of the holdout tracks; without any (a constraint-only run), of every track scored
    switches_holdout = int(judged["switches"].sum()) if len(judged) else 0
    lr = _ratio(float(judged["n_edges_reproduced"].sum()), float(judged["n_edges"].sum()))
    lr_all = _ratio(float(tracks["n_edges_reproduced"].sum()), float(tracks["n_edges"].sum()))
    sr_all = _ratio(float(tracks["n_in_label"].sum()), float(tracks["n_points"].sum()))
    pach = tracks[tracks["zone"].astype(str).str.lower().isin(["pachytene", *STAGING_PACH])]
    sr_pach = _ratio(float(pach["n_in_label"].sum()), float(pach["n_points"].sum())) if len(pach) else float("nan")
    sr = sr_pach if np.isfinite(sr_pach) else sr_all
    sr_set = "pachytene" if np.isfinite(sr_pach) else "all zones"
    by_zone = {str(z): _ratio(float(g["n_in_label"].sum()), float(g["n_points"].sum()))
               for z, g in tracks.groupby(tracks["zone"].fillna("unknown"))}
    crit = {
        "switches_holdout": {"value": switches_holdout, "threshold": max_sw, "rule": "<=",
                             "pass": switches_holdout <= max_sw, "set": judged_set},
        "linking_recall": {"value": lr, "threshold": min_lr, "rule": ">=", "pass": bool(np.isfinite(lr) and lr >= min_lr),
                           "set": judged_set},
        "segmentation_recall": {"value": sr, "threshold": min_sr, "rule": ">=", "pass": bool(np.isfinite(sr) and sr >= min_sr),
                                "set": sr_set},
    }
    return {"pass": all(c["pass"] for c in crit.values()), "criteria": crit,
            "n_hand_tracks": int(len(tracks)), "n_holdout": int(len(hold)), "n_hand_points": int(len(hand)),
            "linking_recall_all": lr_all, "segmentation_recall_all": sr_all, "segmentation_recall_by_zone": by_zone,
            "switches_all": int(tracks["switches"].sum()), "edges_missing_all": int(tracks["edges_missing"].sum()),
            "edge_precision_all": _ratio(float(tracks["n_pipeline_edges_correct"].sum()), float(tracks["n_pipeline_edges"].sum()))}


def format_verdict(res: dict) -> str:
    """The printed report of validate_tracks: the per-track table and the PASS / FAIL lines."""
    tracks, v = res["tracks"], res["verdict"]
    scored = "" if v.get("subset", "all") == "all" else f" (only the {v['subset']} hand tracks are scored)"
    lines = [f"{res['layout'].stem}: {v['n_hand_tracks']} hand track(s), {v['n_hand_points']} point(s), "
             f"{v['n_holdout']} holdout{scored}"]
    cols = ["hand_track_id", "role", "zone", "dominant_time_track_id", "n_points", "n_in_label", "segmentation_recall",
            "linking_recall", "edge_precision", "switches", "edges_missing"]
    with pd.option_context("display.width", 200, "display.max_columns", 20, "display.float_format", "{:.3f}".format):
        lines.append(tracks[cols].to_string(index=False))
    for name, c in v["criteria"].items():
        val = c["value"]
        shown = f"{val:.3f}" if isinstance(val, float) and math.isfinite(val) else str(val)
        lines.append(f"  {'PASS' if c['pass'] else 'FAIL'}  {name} = {shown} ({c['rule']} {c['threshold']}; {c['set']})")
    lines.append(f"  linking recall over every hand track {v['linking_recall_all']:.3f}; segmentation recall by zone: "
                 + ", ".join(f"{z} {r:.3f}" for z, r in v["segmentation_recall_by_zone"].items()))
    lines.append(f"  edge precision {v['edge_precision_all']:.3f}; switches {v['switches_all']}; missing edges {v['edges_missing_all']}")
    lines.append(f"VERDICT: {'PASS' if v['pass'] else 'FAIL'}" + ("" if v.get("subset", "all") == "all"
                                                                   else f" (on the {v['subset']} hand tracks only)"))
    return "\n".join(lines)
