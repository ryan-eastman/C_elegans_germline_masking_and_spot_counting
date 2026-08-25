"""Axis v2: trace an actual PATH along the gonad instead of using a geodesic distance field.

WHY. The old `geodesic_axis` returned `DA`, the geodesic distance from one graph-diameter endpoint, and
called that the axis coordinate. On a straight tube that is fine, but on a branched or two-arm cloud the
level sets of a distance field are iso-distance CONTOURS, so two different arms at the same distance
collapse into the SAME row. Measured consequence on N2 noHS_male_001: 53% of the nuclei in a given row
belonged to a second, fused germline. `largest_component()` cannot catch that, because two touching
gonads form one connected component.

WHAT THIS DOES INSTEAD.
  1. pre-filter obvious non-germline objects by size (spermatid dots / over-segmentation fragments and
     giant somatic-carcass blobs), robustly, before any geometry;
  2. build a kNN graph and take the SHORTEST PATH between the two graph-diameter endpoints - an actual
     ordered polyline through the tube, not a distance field;
  3. smooth that polyline, resample it evenly, and treat it as the tube centreline;
  4. project every nucleus onto the centreline -> arc length s (position along the gonad) and
     perpendicular distance r (how far off the traced tube it sits);
  5. REJECT nuclei whose r exceeds a robust cutoff = the tube radius estimated from the on-path
     population. Nuclei on a second, fused arm sit far off the traced path and are dropped by
     construction, which is exactly what the distance-field version could not do.

Returns a dataframe of the kept, ordered nuclei plus diagnostics (how many were dropped and why), so a
bad gonad announces itself instead of silently producing scrambled rows."""
import numpy as np
from scipy import stats
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, shortest_path
from scipy.spatial import cKDTree

KNN = 10
FRAG_FRAC = 0.30      # objects below this fraction of median volume = spermatid dots / fragments
GIANT_FRAC = 3.0      # objects above this multiple of median volume = somatic / carcass / merged blobs
PERP_MULT = 2.0       # keep nuclei within PERP_MULT x (robust tube radius) of the traced centreline
SMOOTH_PTS = 9        # moving-average window (in path nodes) used to smooth the centreline


def _prefilter(df):
    v = df.vol_env.to_numpy()
    med = float(np.median(v))
    keep = (v >= FRAG_FRAC * med) & (v <= GIANT_FRAC * med)
    return df[keep].reset_index(drop=True), {"dropped_small": int((v < FRAG_FRAC * med).sum()),
                                             "dropped_giant": int((v > GIANT_FRAC * med).sum()),
                                             "median_vol_um3": round(med, 2)}


def _knn_graph(P, k=KNN):
    n = len(P)
    d, j = cKDTree(P).query(P, k=min(k + 1, n))
    src = np.repeat(np.arange(n), d.shape[1] - 1)
    G = coo_matrix((d[:, 1:].ravel(), (src, j[:, 1:].ravel())), shape=(n, n)).tocsr()
    return G.maximum(G.T)


def _largest_cc(P, G):
    ncomp, lab = connected_components(G, directed=False)
    if ncomp == 1:
        return np.ones(len(P), bool), ncomp
    sizes = np.bincount(lab)
    return lab == int(np.argmax(sizes)), ncomp


def _trace_path(P, G):
    """ordered node indices of the shortest path between the two graph-diameter endpoints."""
    D0, pred0 = shortest_path(G, directed=False, indices=0, return_predecessors=True)
    D0[~np.isfinite(D0)] = -1
    A = int(np.argmax(D0))
    DA, predA = shortest_path(G, directed=False, indices=A, return_predecessors=True)
    finite = np.isfinite(DA)
    B = int(np.argmax(np.where(finite, DA, -1)))
    path = []
    cur = B
    while cur != A and cur >= 0:
        path.append(cur)
        cur = predA[cur]
    path.append(A)
    return np.array(path[::-1]), float(DA[B])


def _centreline(P, path, smooth=SMOOTH_PTS):
    C = P[path]
    if len(C) >= smooth:
        ker = np.ones(smooth) / smooth
        C = np.stack([np.convolve(C[:, i], ker, mode="same") for i in range(3)], 1)
        C[:smooth] = P[path][:smooth]      # convolution edge effects: keep the raw ends
        C[-smooth:] = P[path][-smooth:]
    seg = np.linalg.norm(np.diff(C, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(seg)])
    return C, arc


def _project(P, C, arc):
    """nearest centreline node for each nucleus -> (arc length s, perpendicular distance r)."""
    d, idx = cKDTree(C).query(P)
    return arc[idx], d


def orient_by_volume(df, s):
    """germ nuclei grow distal->proximal; Theil-Sen slope over the middle of the axis decides polarity."""
    v = df.vol_env.to_numpy()
    L = float(np.nanmax(s))
    nb, ctr, med = 20, [], []
    edges = np.linspace(0, L, nb + 1)
    for i in range(nb):
        m = (s >= edges[i]) & (s < edges[i + 1])
        if m.sum() >= 5:
            ctr.append(0.5 * (edges[i] + edges[i + 1])); med.append(float(np.median(v[m])))
    ctr, med = np.array(ctr), np.array(med)
    inner = (ctr > 0.1 * L) & (ctr < 0.9 * L)
    slope = float(stats.theilslopes(med[inner], ctr[inner])[0]) if inner.sum() >= 4 else np.nan
    return bool(np.isfinite(slope) and slope < 0), (None if not np.isfinite(slope) else round(slope, 4))


def build_axis(df, perp_mult=PERP_MULT):
    """df needs cz, cy, cx (um) and vol_env. Returns (kept_df_with_s_um_and_r_um, info)."""
    n_in = len(df)
    df, pf = _prefilter(df)
    P = df[["cz", "cy", "cx"]].to_numpy()
    G = _knn_graph(P)
    cc, ncomp = _largest_cc(P, G)
    df, P = df[cc].reset_index(drop=True), P[cc]
    G = _knn_graph(P)
    path, L = _trace_path(P, G)
    C, arc = _centreline(P, path)
    s, r = _project(P, C, arc)
    # robust tube radius from the on-path population, then reject off-path (second-arm) nuclei
    rad = float(np.percentile(r, 75))
    cutoff = max(perp_mult * rad, 3.0)
    on = r <= cutoff
    df = df[on].reset_index(drop=True)
    s, r = s[on], r[on]
    flip, slope = orient_by_volume(df, s)
    if flip:
        s = float(np.nanmax(s)) - s
    df["s_um"], df["r_um"] = s, r
    df = df.sort_values("s_um").reset_index(drop=True)
    info = {"n_in": n_in, "n_out": len(df), **pf, "n_components": int(ncomp),
            "path_len_um": round(float(L), 1), "tube_radius_um": round(rad, 2),
            "perp_cutoff_um": round(cutoff, 2), "dropped_off_path": int((~on).sum()),
            "off_path_frac": round(float((~on).mean()), 3), "flipped": flip, "vol_slope": slope}
    return df, info
