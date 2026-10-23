"""Flow-prior linking, the ``laptrack_flow`` backend of ``germquant track-link`` (docs/TIMELAPSE.md, "Flow-prior
linking").

Why. The germline moves as a tissue: neighbouring nuclei move together. On the first real series the proximal half
slid by up to 3.6 um per 5 minutes while the distal end barely moved (C:/Users/ryane/live_results/FIRST_LOOK.md,
sections 7 and 8). The piecewise registration follows that motion with one translation per axial segment, but a
slide that starts inside a segment, a step longer than the registration's step cap, or the rigid registration alone
leaves the nuclei of a region displaced by a large fraction of their spacing, and plain LAP linking then joins
neighbours. Here every nucleus is matched at the position the motion of its NEIGHBOURS predicts for it.

What, for every step between two consecutive timepoints on disk (t, t'), in the registered frame of the chosen
registration (rigid or piecewise; the flow estimated here is the motion that registration left):

1. Local consensus (``search_um`` > 0). For every nucleus at t, the candidate displacements are the vectors to the
   nuclei at t' within ``search_um``. Its ``k`` nearest neighbours at t vote: the candidate displacement that the
   largest number of them share (another candidate of theirs within VOTE_BANDWIDTH_UM) is the nucleus's initial
   predicted displacement, when at least VOTE_MIN_SUPPORT of the neighbourhood agrees and no rival displacement
   (more than two bandwidths away) has VOTE_AMBIGUITY of that support; otherwise it is zero (the registration).
   Quasi-random packing makes the true common displacement stand out; a wrong partner of each nucleus points in its
   own direction.
2. Confident links. Mutual nearest neighbours between the predicted positions at t' and the nuclei at t', at most
   ``link_radius_um`` apart, that pass the ratio test both ways (nearest distance at most ``ratio`` times the
   second-nearest one).
3. The flow field. Each confident link is a sample of the displacement from t to t'. A link whose displacement
   deviates from the median of its ``k`` nearest other links by more than max(OUTLIER_MIN_UM, OUTLIER_FACTOR x the
   median deviation) is dropped. The field at any point is the component-wise median of the displacements of the
   ``k`` links nearest to it (FlowField): forward with the links placed at their t ends, backward (the negated
   displacement) at their t' ends. A step with fewer than MIN_LINKS confident links keeps the registration alone
   (zero field) and says so in the record.
4. Steps 2 and 3 are repeated ``iterations`` times, each round predicting with the field of the previous one.
5. Predictions. Every nucleus is carried forward along the step fields to each of the next ``gap`` timepoints
   (through the timepoints on disk in between; a missing timepoint is spanned by the step between its neighbours),
   and, with ``fb``, backward the same way.
6. Linking. laptrack as the ``laptrack`` backend runs it (squared Euclidean cost, the cutoff, gap closing over
   ``gap`` frames, splitting and merging off), except that the cost of joining a nucleus at t to one at t + d is the
   squared distance between the first's predicted position at t + d and the second's registered position. With
   ``fb`` (forward-backward consistency) it is the mean of that and the squared distance between the second's
   backward prediction at t and the first's position, so both directions must agree. The cutoff is ``cutoff_um``
   when given, else ``max_dist_um`` (the laptrack backend's cutoff): on the simulator a tighter cutoff on the
   predicted positions removed a few identity switches but cut more true links, because the residual that is left is
   dominated by centroid jumps of the segmentation (partial pieces, splits, nuclei cut by the stack edge), not by
   motion (docs/TIMELAPSE.md). One more setting differs from the laptrack backend, for the same reason the costs are
   small: the gap-closing segment start and end costs are ALTERNATIVE_FACTOR x cutoff squared (laptrack's frame
   linking uses that factor times its largest candidate cost) instead of laptrack's automatic rule, 1.05 x the 90th
   percentile of the gap-closing candidate costs, which assumes the candidates are mostly wrong pairs; on predicted
   positions most are right ones and the rule would refuse gap closings that miss by more than a micron or two
   (GAP_COST_RULE "laptrack" restores it, for comparison).

Nothing is resampled and no table changes frame: the track_points and tracks tables are in the registered frame of
the chosen registration, as for the laptrack backend; the predicted step of every nucleus to the next timepoint
on disk is appended as ``flow_step_z_um``, ``flow_step_y_um``, ``flow_step_x_um``. Every operation is deterministic
(k-d tree queries, medians, a stable sort for the vote, laptrack's LAP): the same series gives the same tracks.
"""
from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

BACKEND = "laptrack_flow"
# code defaults of the timelapse.track.flow_* keys (read with track.cfg_get; no config file names them)
DEFAULTS = {
    "search_um": 6.0,         # local consensus: largest candidate displacement, um (0: start from the registration)
    "link_radius_um": 3.0,    # confident links: largest distance from the predicted position, um
    "ratio": 0.5,             # confident links: nearest over second-nearest distance, both ways
    "k": 8,                   # neighbourhood: nuclei that vote, links the field takes its median over
    "iterations": 2,          # rounds of confident links, field and prediction
    "cutoff_um": None,        # LAP cutoff on the predicted positions, um; None: max_dist_um (timelapse.track.max_dist_um)
    "fb": False,              # forward-backward consistency: the cost averages the forward and the backward miss
}
KEYS = {name: f"timelapse.track.flow_{name}" for name in DEFAULTS}
MIN_LINKS = 8                 # a step with fewer confident links keeps the registration alone (zero field)
VOTE_BANDWIDTH_UM = 1.0       # consensus: candidate displacements this close support each other
VOTE_MIN_SUPPORT = 0.5        # consensus: at least this fraction of the neighbourhood must share the displacement
VOTE_AMBIGUITY = 0.8          # consensus: a rival displacement with this fraction of the best support: no consensus
OUTLIER_FACTOR = 3.0          # field: a link deviating from its neighbours' median by more than
OUTLIER_MIN_UM = 1.0          #   max(OUTLIER_MIN_UM, OUTLIER_FACTOR x the median deviation) is not used
ALTERNATIVE_FACTOR = 1.05     # gap closing: segment start / end cost = this x cutoff squared (step 6)
GAP_COST_RULE = "cutoff"      # "cutoff": the line above; "laptrack": laptrack's automatic percentile rule (comparison only)
_COORDS = ("zr_um", "yr_um", "xr_um")


def flow_settings(cfg=None, **given) -> dict:
    """The flow-prior settings: an explicit value (not None) wins, else the ``timelapse.track.flow_*`` key of ``cfg``
    (a Config, a plain dict or None), else the code default (DEFAULTS). ValueError for an unknown name or a value out
    of range."""
    from .track import cfg_get

    unknown = sorted(set(given) - set(DEFAULTS))
    if unknown:
        raise ValueError(f"unknown flow setting(s) {', '.join(unknown)}; known: {', '.join(DEFAULTS)}")
    raw = {}
    for name, default in DEFAULTS.items():
        v = given.get(name)
        raw[name] = v if v is not None else cfg_get(cfg, KEYS[name], default)
    out = {
        "search_um": _number(raw["search_um"], KEYS["search_um"], minimum=0.0),
        "link_radius_um": _number(raw["link_radius_um"], KEYS["link_radius_um"], minimum=0.0, strict=True),
        "ratio": _number(raw["ratio"], KEYS["ratio"], minimum=0.0, strict=True),
        "k": _integer(raw["k"], KEYS["k"]),
        "iterations": _integer(raw["iterations"], KEYS["iterations"]),
        "cutoff_um": (None if raw["cutoff_um"] is None
                      else _number(raw["cutoff_um"], KEYS["cutoff_um"], minimum=0.0, strict=True)),
        "fb": _flag(raw["fb"], KEYS["fb"]),
    }
    if out["ratio"] > 1.0:
        raise ValueError(f"{KEYS['ratio']} must be in (0, 1], got {raw['ratio']!r}")
    return out


def _number(v, key: str, *, minimum: float, strict: bool = False) -> float:
    if isinstance(v, bool):
        raise ValueError(f"{key} must be a number, got {v!r}")
    try:
        x = float(v)
    except (TypeError, ValueError) as e:
        raise ValueError(f"{key} must be a number, got {v!r}") from e
    if not math.isfinite(x) or x < minimum or (strict and x == minimum):
        raise ValueError(f"{key} must be {'>' if strict else '>='} {minimum:g}, got {v!r}")
    return x


def _integer(v, key: str) -> int:
    if isinstance(v, bool) or not isinstance(v, (int, np.integer)) and not (isinstance(v, float) and v.is_integer()):
        raise ValueError(f"{key} must be a positive integer, got {v!r}")
    if int(v) < 1:
        raise ValueError(f"{key} must be a positive integer, got {v!r}")
    return int(v)


def _flag(v, key: str) -> bool:
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, str) and v.strip().lower() in ("true", "false", "yes", "no", "1", "0"):
        return v.strip().lower() in ("true", "yes", "1")
    if isinstance(v, (int, np.integer)) and int(v) in (0, 1):
        return bool(v)
    raise ValueError(f"{key} must be true or false, got {v!r}")


# =====================================================================================================
# the pieces: the field, the consensus, the confident links
# =====================================================================================================

class FlowField:
    """A displacement field from link samples: at a point, the component-wise median of the displacements of the
    ``k`` samples nearest to it (all samples when there are fewer); zero everywhere without a sample. NaN for a
    non-finite query point."""

    def __init__(self, positions, displacements, k: int):
        self.positions = np.asarray(positions, dtype=float).reshape(-1, 3)
        self.displacements = np.asarray(displacements, dtype=float).reshape(-1, 3)
        self.k = int(k)
        self._tree = None
        if len(self.positions):
            from scipy.spatial import cKDTree

            self._tree = cKDTree(self.positions)

    @classmethod
    def zero(cls) -> FlowField:
        return cls(np.zeros((0, 3)), np.zeros((0, 3)), 1)

    @property
    def n_samples(self) -> int:
        return len(self.positions)

    def __call__(self, x) -> np.ndarray:
        x = np.asarray(x, dtype=float).reshape(-1, 3)
        out = np.zeros_like(x)
        ok = np.isfinite(x).all(axis=1)
        out[~ok] = np.nan
        if self._tree is None or not ok.any():
            return out
        k = min(self.k, self.n_samples)
        _d, idx = self._tree.query(x[ok], k=k)
        idx = np.asarray(idx).reshape(int(ok.sum()), k)
        out[ok] = np.median(self.displacements[idx], axis=1)
        return out


def consensus_displacement(A, B, *, search_um: float, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Step 1 of the module docstring: ``(u, support)`` for the points ``A`` (t) against ``B`` (t'): u (n x 3) the
    displacement most of each point's ``k`` nearest neighbours in A share with some point of B within ``search_um``,
    support (n,) the fraction of the neighbourhood that shares it; u = 0 and support = 0 where there is no clear
    consensus."""
    from scipy.spatial import cKDTree

    A = np.asarray(A, dtype=float).reshape(-1, 3)
    B = np.asarray(B, dtype=float).reshape(-1, 3)
    n = len(A)
    u = np.zeros((n, 3))
    support = np.zeros(n)
    if n == 0 or len(B) == 0 or search_um <= 0:
        return u, support
    cand = cKDTree(B).query_ball_point(A, r=float(search_um))
    vecs = [B[np.asarray(sorted(c), dtype=int)] - A[i] if len(c) else np.zeros((0, 3)) for i, c in enumerate(cand)]
    kk = min(int(k), n)
    _d, nb = cKDTree(A).query(A, k=kk)
    nb = np.asarray(nb).reshape(n, kk)
    bw2 = VOTE_BANDWIDTH_UM ** 2
    for i in range(n):
        parts = [(m, vecs[m]) for m in nb[i] if len(vecs[m])]
        if not parts:
            continue
        V = np.concatenate([v for _m, v in parts])
        owner = np.concatenate([np.full(len(v), j) for j, (_m, v) in enumerate(parts)])
        near = ((V[:, None, :] - V[None, :, :]) ** 2).sum(axis=2) <= bw2           # (c, c)
        has = np.zeros((len(V), len(parts)), dtype=bool)                              # candidate x neighbour
        r, c = np.nonzero(near)
        has[r, owner[c]] = True
        votes = has.sum(axis=1)
        order = np.lexsort((np.arange(len(V)), np.einsum("ij,ij->i", V, V), -votes))  # most votes, then shortest
        best = int(order[0])
        top = int(votes[best])
        rival = (((V - V[best]) ** 2).sum(axis=1) > 4.0 * bw2)
        second = int(votes[rival].max()) if rival.any() else 0
        if top < VOTE_MIN_SUPPORT * kk or second >= VOTE_AMBIGUITY * top:
            continue
        # the displacement: per agreeing neighbour, its candidate closest to the winner; their median
        d2 = ((V - V[best]) ** 2).sum(axis=1)
        picks = []
        for j in range(len(parts)):
            sel = np.flatnonzero((owner == j) & near[best])
            if len(sel):
                picks.append(V[sel[np.argmin(d2[sel])]])
        u[i] = np.median(np.asarray(picks), axis=0)
        support[i] = top / float(kk)
    return u, support


def mutual_links(P, B, *, radius_um: float, ratio: float) -> tuple[np.ndarray, np.ndarray]:
    """Step 2 of the module docstring: the index pairs (i, j) such that B[j] is the nearest point of B to P[i] and
    P[i] the nearest point of P to B[j], at most ``radius_um`` apart, with the nearest distance at most ``ratio``
    times the second-nearest one on both sides (a lone point passes the ratio test)."""
    from scipy.spatial import cKDTree

    P = np.asarray(P, dtype=float).reshape(-1, 3)
    B = np.asarray(B, dtype=float).reshape(-1, 3)
    empty = (np.zeros(0, dtype=int), np.zeros(0, dtype=int))
    okp = np.flatnonzero(np.isfinite(P).all(axis=1))
    okb = np.flatnonzero(np.isfinite(B).all(axis=1))
    if not len(okp) or not len(okb):
        return empty
    Pf, Bf = P[okp], B[okb]

    def two_nearest(tree_pts, query):
        k = min(2, len(tree_pts))
        d, idx = cKDTree(tree_pts).query(query, k=k)
        d = np.asarray(d, dtype=float).reshape(len(query), k)
        idx = np.asarray(idx).reshape(len(query), k)
        if k == 1:
            d = np.concatenate([d, np.full_like(d, np.inf)], axis=1)
        return d, idx

    dB, iB = two_nearest(Bf, Pf)
    dP, iP = two_nearest(Pf, Bf)
    i = np.arange(len(Pf))
    j = iB[:, 0]
    keep = (iP[j, 0] == i) & (dB[:, 0] <= radius_um)
    with np.errstate(divide="ignore", invalid="ignore"):
        keep &= (dB[:, 0] <= ratio * dB[:, 1]) & (dP[j, 0] <= ratio * dP[j, 1])
    return okp[i[keep]], okb[j[keep]]


def _deviation_from_neighbours(S: np.ndarray, W: np.ndarray, k: int) -> np.ndarray:
    """|W_l - median of the W of the k samples nearest to S_l, itself excluded| for every sample l."""
    from scipy.spatial import cKDTree

    n = len(S)
    if n < 2:
        return np.zeros(n)
    kk = min(int(k) + 1, n)
    _d, idx = cKDTree(S).query(S, k=kk)
    idx = np.asarray(idx).reshape(n, kk)
    others = np.empty((n, kk - 1), dtype=int)
    for r in range(n):                      # drop the sample itself wherever the query put it (ties)
        row = idx[r][idx[r] != r]
        others[r] = row[: kk - 1]
    return np.linalg.norm(W - np.median(W[others], axis=1), axis=1)


# =====================================================================================================
# one step, the series
# =====================================================================================================

@dataclass
class StepFlow:
    """The flow between two consecutive timepoints on disk: ``forward`` at points of t_from (displacement to t_to),
    ``backward`` at points of t_to (displacement to t_from), and the step's record entry."""
    t_from: int
    t_to: int
    forward: FlowField
    backward: FlowField
    info: dict = field(default_factory=dict)


def estimate_step(A, B, settings: dict, *, t_from: int = 0, t_to: int = 1) -> StepFlow:
    """Steps 1 to 4 of the module docstring for the registered positions ``A`` (t_from) and ``B`` (t_to)."""
    A = np.asarray(A, dtype=float).reshape(-1, 3)
    B = np.asarray(B, dtype=float).reshape(-1, 3)
    k = int(settings["k"])
    info: dict = {"t_from": int(t_from), "t_to": int(t_to), "n_from": int(len(A)), "n_to": int(len(B)),
                  "n_consensus": 0, "n_confident": 0, "n_used": 0, "fallback": None}
    zero = StepFlow(int(t_from), int(t_to), FlowField.zero(), FlowField.zero(), info)
    fa = np.isfinite(A).all(axis=1)
    if not fa.any() or not np.isfinite(B).all(axis=1).any():
        info["fallback"] = "a frame without nuclei: the registration alone"
        return zero
    u = np.zeros_like(A)
    if settings["search_um"] > 0:
        u[fa], support = consensus_displacement(A[fa], B[np.isfinite(B).all(axis=1)], search_um=settings["search_um"], k=k)
        info["n_consensus"] = int((support > 0).sum())
    flow = None
    for it in range(int(settings["iterations"])):
        ii, jj = mutual_links(A + u, B, radius_um=settings["link_radius_um"], ratio=settings["ratio"])
        if len(ii) < MIN_LINKS:
            if flow is None:
                info.update(n_confident=int(len(ii)),
                            fallback=f"{len(ii)} confident link(s), fewer than {MIN_LINKS}: the registration alone")
                return zero
            break                            # keep the field of the previous round
        W = B[jj] - A[ii]
        dev = _deviation_from_neighbours(A[ii], W, k)
        keep = dev <= max(OUTLIER_MIN_UM, OUTLIER_FACTOR * float(np.median(dev)))
        flow = StepFlow(int(t_from), int(t_to), FlowField(A[ii][keep], W[keep], k), FlowField(B[jj][keep], -W[keep], k), info)
        u = flow.forward(A)
        info.update(n_confident=int(len(ii)), n_used=int(keep.sum()), iterations_run=it + 1,
                    residual_median_um=float(np.median(dev[keep])), residual_p95_um=float(np.percentile(dev[keep], 95)))
    mag = np.linalg.norm(u[fa], axis=1)
    info.update(field_median_um=float(np.median(mag)), field_max_um=float(np.max(mag)))
    info["_residuals"] = dev[keep]
    return flow


def estimate_flow(objects: pd.DataFrame, settings: dict) -> dict[int, StepFlow]:
    """``{t_from: StepFlow}`` for every pair of consecutive timepoints that have nuclei in ``objects`` (``t_index``
    and the registered ``zr_um``, ``yr_um``, ``xr_um``)."""
    ts = sorted(int(t) for t in pd.unique(objects["t_index"]))
    pos = {t: objects.loc[objects["t_index"] == t, list(_COORDS)].to_numpy(float) for t in ts}
    return {a: estimate_step(pos[a], pos[b], settings, t_from=a, t_to=b) for a, b in zip(ts[:-1], ts[1:])}


def predict(objects: pd.DataFrame, flows: dict[int, StepFlow], depth: int, *, backward: bool = False) -> np.ndarray:
    """Step 5 of the module docstring: an (N, depth, 3) array, entry [n, d - 1] the predicted registered position of
    row n of ``objects`` at its timepoint plus d (minus d with ``backward``), carried along the step fields through
    the timepoints in between; NaN where that timepoint has no nuclei or lies beyond ``depth``."""
    t_all = objects["t_index"].astype(int).to_numpy()
    xyz = objects[list(_COORDS)].to_numpy(float)
    out = np.full((len(objects), int(depth), 3), np.nan)
    ts = sorted(set(t_all.tolist()))
    back = {f.t_to: f for f in flows.values()}
    for t in ts:
        rows = np.flatnonzero(t_all == t)
        x = xyz[rows].copy()
        cur = t
        while True:
            step = back.get(cur) if backward else flows.get(cur)
            if step is None:
                break
            nxt = step.t_from if backward else step.t_to
            d = abs(nxt - t)
            if d > depth:
                break
            x = x + (step.backward(x) if backward else step.forward(x))
            out[rows, d - 1] = x
            cur = nxt
    return out


# =====================================================================================================
# the backend
# =====================================================================================================

def link_flow(objects: pd.DataFrame, *, max_dist_um: float, gap_frames: int, settings: dict
              ) -> tuple[dict[tuple[int, int], int], dict, pd.DataFrame]:
    """The ``laptrack_flow`` backend (module docstring): ``({(t_index, nucleus_id): time_track_id}, the backend
    record, the flow_step_* columns aligned with ``objects``)``. ``objects`` is the object table of link_series (the
    registered centroids ``zr_um``, ``yr_um``, ``xr_um``); ``max_dist_um`` and ``gap_frames`` mean what they mean for
    the laptrack backend; ``settings`` comes from flow_settings."""
    try:
        from laptrack import LapTrack
    except ImportError as e:
        raise RuntimeError("the laptrack_flow backend needs the timelapse extra: uv pip install -e \".[timelapse]\" "
                           "(laptrack>=0.17)") from e

    extra = pd.DataFrame({"flow_step_z_um": np.nan, "flow_step_y_um": np.nan, "flow_step_x_um": np.nan},
                         index=objects.index)
    cutoff = float(settings["cutoff_um"]) if settings.get("cutoff_um") is not None else float(max_dist_um)
    info = {"max_dist_um": float(max_dist_um), "gap_frames": int(gap_frames), "cutoff_um": cutoff,
            "cutoff_sq_um2": cutoff ** 2,
            "cutoff_source": "flow_cutoff_um" if settings.get("cutoff_um") is not None else "max_dist_um",
            "flow": {**settings, "min_links": MIN_LINKS, "vote_bandwidth_um": VOTE_BANDWIDTH_UM,
                     "vote_min_support": VOTE_MIN_SUPPORT, "vote_ambiguity": VOTE_AMBIGUITY,
                     "outlier_factor": OUTLIER_FACTOR, "outlier_min_um": OUTLIER_MIN_UM},
            "implementation": "germquant.timelapse.flowlink (local consensus, mutual nearest neighbours with a ratio "
                              "test, k-nearest median field) + laptrack.LapTrack.predict_dataframe on the predicted "
                              "positions"}
    df = objects[["t_index", "nucleus_id", *_COORDS]].reset_index(drop=True)
    if df.empty:
        info.update(n_points=0, n_steps=0, steps=[])
        return {}, info, extra
    flows = estimate_flow(df, settings)
    D = max(1, int(gap_frames))
    fwd = predict(df, flows, D)
    cols = {"fl_frame": (df["t_index"].astype(int) - int(df["t_index"].min())).astype(float)}
    for c in _COORDS:
        cols[c] = df[c].to_numpy(float)
    for d in range(D):
        for a, c in enumerate(_COORDS):
            cols[f"fwd{d + 1}_{c}"] = fwd[:, d, a]
    fb = bool(settings["fb"])
    if fb:
        bwd = predict(df, flows, D, backward=True)
        for d in range(D):
            for a, c in enumerate(_COORDS):
                cols[f"bwd{d + 1}_{c}"] = bwd[:, d, a]
    coord_cols = list(cols)
    L = pd.DataFrame(cols)
    L["frame__"] = (df["t_index"].astype(int) - int(df["t_index"].min())).to_numpy()
    L["row__"] = np.arange(len(df))            # the object row, mapped back below (a narrow frame keeps laptrack fast)
    off_f, off_b = 4, 4 + 3 * D

    def cost(u, v, **_kw):
        """Squared miss of the earlier point's prediction at the later point's frame (with fb: averaged with the
        later point's backward miss); inf beyond the prediction depth or without a prediction."""
        d = int(v[0] - u[0] + 0.5)
        if d < 1 or d > D:
            return math.inf
        o = off_f + 3 * (d - 1)
        a0, a1, a2 = u[o] - v[1], u[o + 1] - v[2], u[o + 2] - v[3]
        c = a0 * a0 + a1 * a1 + a2 * a2
        if fb:
            o = off_b + 3 * (d - 1)
            b0, b1, b2 = v[o] - u[1], v[o + 1] - u[2], v[o + 2] - u[3]
            c = 0.5 * (c + b0 * b0 + b1 * b1 + b2 * b2)
        return c if c == c else math.inf

    def cost_next(u, v, **_kw):
        """cost() for consecutive frames (laptrack's frame linking), without the frame arithmetic."""
        a0, a1, a2 = u[4] - v[1], u[5] - v[2], u[6] - v[3]
        c = a0 * a0 + a1 * a1 + a2 * a2
        if fb:
            b0, b1, b2 = v[off_b] - u[1], v[off_b + 1] - u[2], v[off_b + 2] - u[3]
            c = 0.5 * (c + b0 * b0 + b1 * b1 + b2 * b2)
        return c if c == c else math.inf

    cut2 = cutoff ** 2
    alt = {}
    if GAP_COST_RULE == "cutoff":
        # laptrack's automatic segment costs (1.05 x the 90th percentile of the gap-closing candidate costs) assume
        # the candidates are mostly wrong pairs; on predicted positions most are the right ones and the rule would
        # refuse gap closings that miss by more than a micron or two (half a micron on a clean synthetic series).
        # With 1.05 x cutoff squared, as in laptrack's frame linking, a gap is closed whenever the miss is within the
        # cutoff and the assignment has no better use for the two ends.
        alt = {"segment_start_cost": ALTERNATIVE_FACTOR * cut2, "segment_end_cost": ALTERNATIVE_FACTOR * cut2,
               "no_splitting_cost": ALTERNATIVE_FACTOR * cut2, "no_merging_cost": ALTERNATIVE_FACTOR * cut2}
    lt = LapTrack(metric=cost_next, cutoff=cut2, gap_closing_metric=cost,
                  gap_closing_cutoff=cut2 if int(gap_frames) > 0 else False,
                  gap_closing_max_frame_count=max(1, int(gap_frames)),
                  splitting_cutoff=False, merging_cutoff=False, **alt)
    track_df, split_df, merge_df = lt.predict_dataframe(L, coord_cols, frame_col="frame__", only_coordinate_cols=False)
    track_df = track_df.reset_index()
    t_of = df["t_index"].astype(int).to_numpy()
    n_of = df["nucleus_id"].astype(int).to_numpy()
    assign = {(int(t_of[r]), int(n_of[r])): int(tid) + 1
              for r, tid in zip(track_df["row__"].astype(int), track_df["track_id"].astype(int))}

    # the predicted step of every nucleus to the next timepoint on disk, and the record of every step
    steps = []
    residuals = []
    t_all = df["t_index"].astype(int).to_numpy()
    xyz = df[list(_COORDS)].to_numpy(float)
    step_cols = np.full((len(df), 3), np.nan)
    for t_from in sorted(flows):
        f = flows[t_from]
        rows = np.flatnonzero(t_all == t_from)
        step_cols[rows] = f.forward(xyz[rows])
        entry = {k: v for k, v in f.info.items() if not k.startswith("_")}
        steps.append(entry)
        if "_residuals" in f.info:
            residuals.append(np.asarray(f.info["_residuals"], dtype=float))
    extra.loc[:, ["flow_step_z_um", "flow_step_y_um", "flow_step_x_um"]] = step_cols
    res = np.concatenate(residuals) if residuals else np.zeros(0)
    info.update({
        "n_points": int(len(track_df)), "n_splits": int(len(split_df)), "n_merges": int(len(merge_df)),
        "n_steps": len(steps), "n_steps_registration_only": int(sum(1 for s in steps if s["fallback"])),
        "residual_median_um": float(np.median(res)) if len(res) else None,
        "residual_p95_um": float(np.percentile(res, 95)) if len(res) else None,
        "field_median_um": float(np.nanmedian(np.linalg.norm(step_cols, axis=1))) if np.isfinite(step_cols).any() else None,
        "steps": steps,
    })
    return assign, info, extra
