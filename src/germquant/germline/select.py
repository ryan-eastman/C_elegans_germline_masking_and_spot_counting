"""Germline isolation: flag which segmented nuclei belong to the gonad vs gut autofluorescence,
debris, or nuclei outside the germline (the "nuclei called outside the germline" montage problem).

Uses the SYP / central-element signal as the germline marker plus spatial connectivity — and NEVER
the spot count, so the selection is independent of the RAD-51 readout it protects. That independence
lets us validate it with spot-recall (spot-bearing cells are germline) without circularity.

Methods (config `germline.method`):
  syp_seeded_cc — DEFAULT. Seed on SYP-bright (synapsed) nuclei, then keep the WHOLE spatially-connected
                nuclear component(s) containing those seeds. Crucially this keeps the SYP-NEGATIVE distal
                tip (mitotic + transition zone, before SYP-3 loads) because it is contiguous with the
                pachytene core — the other methods, which threshold each nucleus on SYP independently,
                wrongly drop the distal third of the gonad (verified on real data). 14-gonad validation:
                mean spot-recall ~0.99, ~7% dropped (the separated gut/debris clusters).
  multi_cc    — keep the largest SYP-positive component + any other >= 10% of its size. Recovers split
                gonads but, like all SYP-threshold methods, drops the SYP-negative distal tip.
  largest_cc  — single largest SYP-positive component only. Simplest; truncates split gonads.
  seed_mass   — components carrying enough strong-SYP seeds; over-drops (~47%).

All add a boolean `in_germline` column and return (df, flags). Geometry is in microns. The selection
NEVER reads the spot count, so spot-recall is an independent (non-circular) validation — see
scripts/validate_germline_select.py.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

_XYZ = ["centroid_z_um", "centroid_y_um", "centroid_x_um"]
_SYP = "central_element_mean_intensity"


def _components(xyz: np.ndarray, radius: float):
    """Connected components of nuclei whose centroids are within `radius` µm (single-linkage)."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    n = len(xyz)
    if n == 0:
        return np.array([], dtype=int), 0
    pairs = cKDTree(xyz).query_pairs(radius, output_type="ndarray")
    if len(pairs) == 0:
        return np.arange(n), n
    rows = np.r_[pairs[:, 0], pairs[:, 1]]
    cols = np.r_[pairs[:, 1], pairs[:, 0]]
    g = csr_matrix((np.ones(rows.size), (rows, cols)), shape=(n, n))
    ncomp, lab = connected_components(g, directed=False)
    return lab, ncomp


def _largest_cc(xyz, syp, *, syp_percentile=25.0, link_radius_um=12.0):
    out = np.zeros(len(syp), dtype=bool)
    keep = syp >= np.percentile(syp, syp_percentile)
    idx = np.where(keep)[0]
    if idx.size == 0:
        return out
    lab, _ = _components(xyz[idx], link_radius_um)
    out[idx[lab == np.bincount(lab).argmax()]] = True
    return out


def _syp_seeded_cc(xyz, syp, *, syp_percentile=25.0, link_radius_um=12.0, min_seed_frac=0.10):
    """Seed on the SYP-bright synapsed core, then keep the WHOLE spatially-connected nuclear component
    containing it — so the SYP-NEGATIVE distal tip (mitotic + transition zone, before SYP loads) is
    kept because it is contiguous with the pachytene core, while separate gut/debris clusters (which
    carry no synapsed seeds) are dropped. Fixes thresholding each nucleus on SYP independently, which
    wrongly cuts off the distal third of the gonad."""
    out = np.zeros(len(syp), dtype=bool)
    seed = syp >= np.percentile(syp, syp_percentile)
    total_seed = int(seed.sum())
    if total_seed == 0:
        return out
    lab, ncomp = _components(xyz, link_radius_um)     # connect ALL nuclei, not just SYP+
    for c in range(ncomp):
        members = lab == c
        if seed[members].sum() >= max(min_seed_frac * total_seed, 5):
            out[members] = True
    return out


def _multi_cc(xyz, syp, *, syp_percentile=25.0, link_radius_um=12.0, size_frac=0.10):
    """Keep the largest SYP-positive component PLUS any other SYP-positive component that is
    substantial (>= size_frac x the largest's size). Members already passed the SYP>=p threshold,
    so no extra brightness gate is applied — that lets a DIMMER second arm (e.g. a folded/split male
    gonad, whose spermatocytes carry less SYP) be recovered, which a brightness-ratio gate drops."""
    out = np.zeros(len(syp), dtype=bool)
    keep = syp >= np.percentile(syp, syp_percentile)
    idx = np.where(keep)[0]
    if idx.size == 0:
        return out
    lab, ncomp = _components(xyz[idx], link_radius_um)
    counts = np.bincount(lab)
    big = int(counts.argmax())
    keep_comps = {c for c in range(ncomp) if counts[c] >= size_frac * counts[big]}
    out[idx[np.isin(lab, list(keep_comps))]] = True
    return out


def _seed_mass(xyz, syp, *, link_radius_um=12.0, smooth_radius_um=8.0, cand_percentile=35.0,
               min_component=15):
    """Robust-z SYP, smooth over neighbours, keep components carrying enough strong-SYP seeds."""
    from scipy.spatial import cKDTree

    n = len(syp)
    out = np.zeros(n, dtype=bool)
    med = np.median(syp)
    mad = np.median(np.abs(syp - med)) * 1.4826 + 1e-9
    z = (syp - med) / mad
    tree = cKDTree(xyz)
    zsm = np.array([np.median(z[ix]) for ix in tree.query_ball_point(xyz, smooth_radius_um)])
    cand = zsm >= np.percentile(zsm, cand_percentile)
    seed = (z >= 1.0) & (zsm >= 0.5)
    idx = np.where(cand)[0]
    if idx.size < min_component:
        return out
    lab, ncomp = _components(xyz[idx], link_radius_um)
    seed_sub = seed[idx]
    comp_seed = np.array([seed_sub[lab == c].sum() for c in range(ncomp)])
    comp_size = np.array([(lab == c).sum() for c in range(ncomp)])
    if comp_seed.max() <= 0:
        return out
    keep_comp = (comp_seed >= max(0.25 * comp_seed.max(), 5)) & (comp_size >= min_component)
    out[idx[keep_comp[lab]]] = True
    return out


def select_germline(
    nuclei: pd.DataFrame,
    *,
    method: str = "syp_seeded_cc",
    syp_percentile: float = 25.0,
    link_radius_um: float = 12.0,
    min_seed_frac: float = 0.10,
    size_frac: float = 0.10,
    smooth_radius_um: float = 8.0,
    cand_percentile: float = 35.0,
    min_component: int = 15,
) -> tuple[pd.DataFrame, list[str]]:
    """Add a boolean `in_germline` column. Returns (df, flags).

    If the SYP/central-element column or the centroids are missing (no SC channel), every nucleus is
    kept (`in_germline=True`) and a flag is raised — never silently drop everything.
    """
    df = nuclei.copy()
    flags: list[str] = []
    if df.empty:
        df["in_germline"] = pd.Series([], dtype=bool)
        return df, flags
    if _SYP not in df.columns or not set(_XYZ).issubset(df.columns):
        df["in_germline"] = True
        flags.append("germline:no_syp_channel_kept_all")
        return df, flags

    xyz = df[_XYZ].to_numpy(dtype=float)
    syp = df[_SYP].to_numpy(dtype=float)
    if method == "seed_mass":
        sel = _seed_mass(xyz, syp, link_radius_um=link_radius_um, smooth_radius_um=smooth_radius_um,
                         cand_percentile=cand_percentile, min_component=min_component)
    elif method == "multi_cc":
        sel = _multi_cc(xyz, syp, syp_percentile=syp_percentile, link_radius_um=link_radius_um,
                        size_frac=size_frac)
    elif method == "largest_cc":
        sel = _largest_cc(xyz, syp, syp_percentile=syp_percentile, link_radius_um=link_radius_um)
    else:   # syp_seeded_cc (default)
        sel = _syp_seeded_cc(xyz, syp, syp_percentile=syp_percentile, link_radius_um=link_radius_um,
                             min_seed_frac=min_seed_frac)

    if not sel.any():   # degenerate (e.g. flat SYP) — keep all rather than emit an empty germline
        df["in_germline"] = True
        flags.append(f"germline:{method}_selected_none_kept_all")
        return df, flags

    df["in_germline"] = sel
    n_drop = int((~sel).sum())
    if n_drop:
        flags.append(f"germline:dropped_{n_drop}_of_{len(df)}_nuclei")
    return df, flags
