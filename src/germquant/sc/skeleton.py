"""Per-nucleus synaptonemal-complex tracing -> length & FRAGMENTATION.

For each nucleus we crop the SC channel, resample it to **isotropic** voxels (so the ridge
filter and skeleton aren't biased by z-anisotropy), enhance the filament with a Sato (tubeness)
ridge filter at physical scales, threshold *inside the nucleus*, skeletonize in 3D, then measure
each disconnected component with `skan` using the isotropic voxel size. The heat-phenotype
readout is **fragmentation**: the number of disconnected skeleton components (fragments) per
nucleus and their length distribution.

Two correctness details that matter for the fragment count:
  * the ridge filter runs on the *unmasked* crop and is only confined to the nucleus afterwards
    (eroded by 1 voxel) — hard-zeroing the crop first would create a sharp edge the tubeness
    filter latches onto, inflating fragment counts;
  * every attempted nucleus emits a per-nucleus row, with ``n_fragments=0`` when nothing traces,
    so a fully-desynapsed nucleus is a real zero in the distribution, not a dropped row.

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
    channel_role: str = "central_element",
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
    target = float(sp.min())                       # isotropic voxel size for resampling (µm)
    zoom = sp / target                             # per-axis upsample factors -> isotropic
    iso_sp = (target, target, target)
    sigmas_iso = [float(s / target) for s in ridge_sigmas_um]
    min_len = float(min_fragment_length_um)
    exp_map = expected_n_tracks or {}

    track_rows, nuc_rows = [], []
    track_id = 0
    objects = ndi.find_objects(labels)
    for nid_minus1, sl in enumerate(objects):
        nid = nid_minus1 + 1
        if sl is None:
            continue
        sub_lab = labels[sl] == nid
        sub_sc = sc_img[sl].astype(np.float32)
        mean_int = float(sub_sc[sub_lab].mean()) if sub_lab.any() else 0.0
        exp = exp_map.get(nid)

        lengths, perfrag = _trace_one(
            sub_sc, sub_lab, zoom, iso_sp, sigmas_iso, intensity_percentile, min_len,
            ndi, sato, skeletonize,
        )
        if not lengths:                            # attempted but nothing traced -> real zero
            nuc_rows.append(_empty_nuc_row(nid, marker, mean_int, exp))
            continue

        for length, (n_br, n_jn, tort) in zip(lengths, perfrag):
            track_rows.append({
                "track_id": track_id, "nucleus_id": nid, "channel_role": channel_role,
                "marker": marker, "length_um": length, "n_branches": n_br,
                "n_junctions": n_jn, "tortuosity": tort, "trace_method": "skan_iso",
            })
            track_id += 1

        arr = np.asarray(lengths)
        nuc_rows.append({
            "nucleus_id": nid, "marker": marker, "n_fragments": int(arr.size),
            "sc_total_length_um": float(arr.sum()),
            "sc_mean_fragment_um": float(arr.mean()),
            "sc_median_fragment_um": float(np.median(arr)),
            "sc_longest_fragment_um": float(arr.max()),
            "sc_mean_intensity": mean_int,
            "expected_n_tracks": exp if exp is not None else float("nan"),
        })

    return (
        pd.DataFrame(track_rows, columns=TRACK_COLS),
        pd.DataFrame(nuc_rows, columns=PER_NUC_COLS),
    )


def _empty_nuc_row(nid, marker, mean_int, exp):
    return {
        "nucleus_id": nid, "marker": marker, "n_fragments": 0,
        "sc_total_length_um": 0.0, "sc_mean_fragment_um": float("nan"),
        "sc_median_fragment_um": float("nan"), "sc_longest_fragment_um": 0.0,
        "sc_mean_intensity": mean_int,
        "expected_n_tracks": exp if exp is not None else float("nan"),
    }


def _trace_one(sub_sc, sub_lab, zoom, iso_sp, sigmas_iso, pct, min_len, ndi, sato, skeletonize):
    """Trace one nucleus crop. Returns (fragment_lengths, per_fragment_stats)."""
    try:
        # resample crop + mask to isotropic voxels (intensity: linear; mask: nearest)
        sc_iso = ndi.zoom(sub_sc, zoom, order=1)
        lab_iso = ndi.zoom(sub_lab.astype(np.float32), zoom, order=0) > 0.5
        if lab_iso.sum() < 2:
            return [], []

        # ridge filter on the UNMASKED crop (no artificial edge), then confine to the
        # nucleus eroded by 1 voxel so we keep only interior filament, not the boundary.
        ridge = sato(sc_iso, sigmas=sigmas_iso, black_ridges=False)
        interior = ndi.binary_erosion(lab_iso)
        if interior.sum() < 2:
            interior = lab_iso
        vals = ridge[interior]
        pos = vals[vals > 0]
        if pos.size == 0:
            return [], []
        thr = np.percentile(pos, pct)
        mask = (ridge >= thr) & interior
        if mask.sum() < 2:
            return [], []

        skel = skeletonize(mask)
        if skel.sum() < 2:
            return [], []
        return _skan_fragments(skel, iso_sp, min_len)
    except Exception as e:  # noqa: BLE001 - one bad nucleus must not kill the batch
        log.debug("SC trace failed for a nucleus (%s)", e)
        return [], []


def _skan_fragments(skel: np.ndarray, spacing, min_len: float):
    """Per-connected-component (fragment) length (µm) + per-fragment branch/junction/tortuosity.

    tortuosity = path length / summed straight-line branch chords (>=1; 1.0 = straight); for a
    single-branch fragment this is exactly branch_distance / euclidean_distance.
    """
    import skan

    try:
        sk = skan.Skeleton(skel.astype(bool), spacing=spacing)
    except ValueError:
        return [], []
    summary = skan.summarize(sk, separator="_")
    id_col = "skeleton_id" if "skeleton_id" in summary else "skeleton-id"
    len_col = "branch_distance" if "branch_distance" in summary else "branch-distance"
    euc_col = "euclidean_distance" if "euclidean_distance" in summary else "euclidean-distance"
    src_col = "node_id_src" if "node_id_src" in summary else "node-id-src"
    dst_col = "node_id_dst" if "node_id_dst" in summary else "node-id-dst"

    lengths, stats = [], []
    for _, grp in summary.groupby(id_col):
        total = float(grp[len_col].sum())
        if total < min_len:
            continue
        n_branches = int(len(grp))
        # junction nodes = skeleton nodes of degree >= 3 (incident to >=3 branch endpoints)
        if src_col in grp and dst_col in grp:
            deg = pd.concat([grp[src_col], grp[dst_col]]).value_counts()
            n_junctions = int((deg >= 3).sum())
        else:
            n_junctions = 0
        euc = float(grp[euc_col].sum()) if euc_col in grp else 0.0
        tort = float(total / euc) if euc > 0 else float("nan")
        lengths.append(total)
        stats.append((n_branches, n_junctions, tort))
    return lengths, stats
