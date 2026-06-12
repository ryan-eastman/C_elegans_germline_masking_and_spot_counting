"""Per-nucleus synaptonemal-complex tracing -> length & FRAGMENTATION.

For each nucleus we crop the SC channel, enhance the filament with a Sato (tubeness)
ridge filter at physical scales, threshold, skeletonize in 3D, then measure with `skan`
using the real voxel spacing. The heat phenotype readout is **fragmentation**: the number
of disconnected skeleton components (fragments) per nucleus and their length distribution.

  *** v1 — VALIDATE BEFORE PUBLISHING. ***
Dense pachytene SCs can be mis-merged or fragmented by automatic skeletonization. Per
ARCHITECTURE.md Risk 1, validate these lengths against an Imaris-FilamentTracer / Fiji-SNT
ground-truth subset (~15-25 nuclei) and report Bland-Altman agreement. `skan` is an optional
dependency (pulls numba); if absent, SC tracing is skipped and a QC flag is raised.
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

TRACK_COLS = [
    "track_id", "nucleus_id", "channel_role", "marker",
    "length_um", "n_branches", "n_junctions", "tortuosity", "trace_method",
]
PER_NUC_COLS = [
    "nucleus_id", "marker", "n_fragments", "sc_total_length_um",
    "sc_mean_fragment_um", "sc_median_fragment_um", "sc_longest_fragment_um",
    "sc_mean_intensity", "expected_n_tracks",
]


def skan_available() -> bool:
    try:
        import skan  # noqa: F401

        return True
    except Exception:
        return False


def trace_sc(
    sc_img: np.ndarray,
    labels: np.ndarray,
    spacing: tuple[float, float, float],
    *,
    marker: str = "SYP-3",
    ridge_sigmas_um=(0.15, 0.25, 0.40),
    intensity_percentile: float = 99.0,
    min_fragment_length_um: float = 0.5,
    expected_n_tracks: dict[int, int] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (sc_tracks_df, sc_per_nucleus_df)."""
    if not skan_available():
        log.warning("skan not installed; skipping SC tracing (install germquant[sc]).")
        return pd.DataFrame(columns=TRACK_COLS), pd.DataFrame(columns=PER_NUC_COLS)

    from scipy import ndimage as ndi
    from skimage.filters import sato
    from skimage.morphology import skeletonize

    sp = np.asarray(spacing, dtype=float)
    sigmas_vox = [float(np.mean(s / sp)) for s in ridge_sigmas_um]
    min_len = float(min_fragment_length_um)

    track_rows, nuc_rows = [], []
    track_id = 0
    objects = ndi.find_objects(labels)
    for nid_minus1, sl in enumerate(objects):
        nid = nid_minus1 + 1
        if sl is None:
            continue
        sub_lab = labels[sl] == nid
        sub_sc = sc_img[sl].astype(np.float32) * sub_lab

        ridge = sato(sub_sc, sigmas=sigmas_vox, black_ridges=False)
        inside = ridge[sub_lab]
        if inside.size == 0 or inside.max() <= 0:
            continue
        thr = np.percentile(inside[inside > 0], intensity_percentile) if (inside > 0).any() else 0
        mask = (ridge >= thr) & sub_lab
        if mask.sum() < 2:
            continue

        skel = skeletonize(mask)
        if skel.sum() < 2:
            continue

        lengths, n_br, n_jn = _skan_lengths(skel, spacing)
        lengths = [length for length in lengths if length >= min_len]
        if not lengths:
            continue

        for length in lengths:
            track_rows.append(
                {
                    "track_id": track_id, "nucleus_id": nid, "channel_role": "central_element",
                    "marker": marker, "length_um": length, "n_branches": n_br,
                    "n_junctions": n_jn, "tortuosity": float("nan"), "trace_method": "skan",
                }
            )
            track_id += 1

        arr = np.asarray(lengths)
        exp = (expected_n_tracks or {}).get(nid)
        nuc_rows.append(
            {
                "nucleus_id": nid, "marker": marker, "n_fragments": int(arr.size),
                "sc_total_length_um": float(arr.sum()),
                "sc_mean_fragment_um": float(arr.mean()),
                "sc_median_fragment_um": float(np.median(arr)),
                "sc_longest_fragment_um": float(arr.max()),
                "sc_mean_intensity": float(sub_sc[sub_lab].mean()),
                "expected_n_tracks": exp if exp is not None else float("nan"),
            }
        )

    return (
        pd.DataFrame(track_rows, columns=TRACK_COLS),
        pd.DataFrame(nuc_rows, columns=PER_NUC_COLS),
    )


def _skan_lengths(skel: np.ndarray, spacing) -> tuple[list[float], int, int]:
    """Per-connected-component skeleton length (µm) via skan, plus branch/junction counts."""
    import skan

    try:
        sk = skan.Skeleton(skel.astype(bool), spacing=spacing)
    except ValueError:
        return [], 0, 0
    summary = skan.summarize(sk, separator="_")
    # each disconnected skeleton (a "fragment") has a unique skeleton_id
    id_col = "skeleton_id" if "skeleton_id" in summary else "skeleton-id"
    len_col = "branch_distance" if "branch_distance" in summary else "branch-distance"
    lengths = summary.groupby(id_col)[len_col].sum().tolist()
    n_branches = int(len(summary))
    n_junctions = int((summary.get("branch_type", 0) == 2).sum()) if "branch_type" in summary else 0
    return [float(length) for length in lengths], n_branches, n_junctions
