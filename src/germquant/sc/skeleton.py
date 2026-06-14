"""Per-nucleus synaptonemal-complex tracing -> length & a FRAGMENTATION index.

For each nucleus we crop the SC channel, resample to **isotropic** voxels (so the ridge filter
and skeleton aren't biased by z-anisotropy), enhance the filament with a Sato (tubeness) ridge
filter at physical scales, build a **continuous** mask with hysteresis thresholding *inside the
nucleus* (a single high percentile shreds the thin SC into disconnected blobs — verified on
synthetic data), skeletonize in 3D, and measure with `skan`.

  *** WHAT IS AND ISN'T RECOVERABLE (validated on synthetic SC ground truth, scripts/validate_sc_tracer.py) ***
  * sc_total_length_um — RELIABLE. Continuous-mask skeleton length tracks the true SC length
    (corr ~0.67, ~35 vs 35 µm). This is the primary, trustworthy readout.
  * n_fragments (disconnected components) — a LOWER BOUND, NOT the biological fragment count. The
    ~6 SCs in a pachytene nucleus physically overlap in 3D, so even the clean signal yields ~2-3
    components, not 6 — and lateral merging bridges each strand's heat-gaps. Per-nucleus fragment
    COUNT is not recoverable from light microscopy at this density (ARCHITECTURE.md Risk 1); use
    Imaris/SNT manual tracing for absolute counts.
  * sc_fragmentation_index — SYP intensity coefficient-of-variation within the nucleus. Heat
    fragmentation raises it (gaps add dark/bright contrast). PER-NUCLEUS it is noisy (corr ~0.21
    with true frag count) but at the POPULATION level it separates control from heat (Cohen's
    d ~0.57). Use it for the control-vs-heat comparison, and CALIBRATE against Imaris before
    quoting absolute fragmentation. It is reported alongside n_fragments, never instead of it.

`skan` is an optional dependency (pulls numba); if absent, SC tracing is skipped + a QC flag set.
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
    "sc_mean_intensity", "sc_fragmentation_index", "expected_n_tracks",
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
    intensity_percentile: float = 99.0,   # deprecated (old single-threshold method); kept for back-compat
    ridge_hyst_low_pct: float = 45.0,     # hysteresis grow level (continuous strand mask)
    ridge_hyst_high_pct: float = 80.0,    # hysteresis seed level
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
        in_vals = sub_sc[sub_lab] if sub_lab.any() else np.zeros(1, np.float32)
        mean_int = float(in_vals.mean())
        # fragmentation INDEX = SYP intensity coefficient-of-variation within the nucleus. Heat
        # fragmentation raises it (gaps add dark/bright contrast). Population-level proxy (per-nucleus
        # noisy) — the heat-vs-control readout, since exact fragment COUNT is unrecoverable (docstring).
        frag_index = float(in_vals.std() / mean_int) if mean_int > 0 else float("nan")
        exp = exp_map.get(nid)

        lengths, perfrag = _trace_one(
            sub_sc, sub_lab, zoom, iso_sp, sigmas_iso, ridge_hyst_low_pct, ridge_hyst_high_pct, min_len,
            ndi, sato, skeletonize,
        )
        if not lengths:                            # attempted but nothing traced -> real zero
            nuc_rows.append(_empty_nuc_row(nid, marker, mean_int, frag_index, exp))
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
            "sc_fragmentation_index": frag_index,
            "expected_n_tracks": exp if exp is not None else float("nan"),
        })

    return (
        pd.DataFrame(track_rows, columns=TRACK_COLS),
        pd.DataFrame(nuc_rows, columns=PER_NUC_COLS),
    )


def _empty_nuc_row(nid, marker, mean_int, frag_index, exp):
    return {
        "nucleus_id": nid, "marker": marker, "n_fragments": 0,
        "sc_total_length_um": 0.0, "sc_mean_fragment_um": float("nan"),
        "sc_median_fragment_um": float("nan"), "sc_longest_fragment_um": 0.0,
        "sc_mean_intensity": mean_int, "sc_fragmentation_index": frag_index,
        "expected_n_tracks": exp if exp is not None else float("nan"),
    }


def _trace_one(sub_sc, sub_lab, zoom, iso_sp, sigmas_iso, low_pct, high_pct, min_len,
               ndi, sato, skeletonize):
    """Trace one nucleus crop. Returns (fragment_lengths, per_fragment_stats).

    Uses HYSTERESIS thresholding (seed at high_pct, grow down to low_pct) to keep the thin SC as a
    CONTINUOUS strand rather than the disconnected blobs a single high percentile produces — this is
    what recovers the SC length (the old single-threshold under-measured it ~2x; verified on synthetic
    SC ground truth, scripts/validate_sc_tracer.py).
    """
    try:
        from skimage.filters import apply_hysteresis_threshold

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
        pos = ridge[interior]
        pos = pos[pos > 0]
        if pos.size == 0:
            return [], []
        lo, hi = np.percentile(pos, low_pct), np.percentile(pos, high_pct)
        if hi <= lo:                                  # near-flat ridge -> single-level fallback
            mask = (ridge >= hi) & interior
        else:
            mask = apply_hysteresis_threshold(ridge, lo, hi) & interior
        if mask.sum() < 2:
            return [], []

        # Mask smoothing before skeletonize -> far fewer spurs/loops/speckles on dense nuclei (the
        # user's "skeleton looks bad" case). The Sato ridge renders a thin SC as a double-walled tube
        # that skeletonizes into two parallel rails + rungs; in-plane closing fuses it to ONE
        # centerline, fill-holes removes spurious loops, small-object removal drops orphan speckle.
        # Strand-preserving: NO opening (it erodes the 1-2 voxel SC strands). Picked over spur-pruning
        # and ridge-pre-smoothing by a validated bake-off (the latter bridges nearby strands into
        # false rings). Effect: detected SC length ~7% shorter (calibrate vs Imaris); per-nucleus
        # length correlation unchanged; the fragmentation index is on raw SYP so it is untouched.
        from skimage.morphology import disk

        mask = ndi.binary_closing(mask, disk(1)[None, :, :])         # fuse the double-rail tube
        # fill small interior holes (spurious skeleton loops) and drop orphan speckle blobs, via
        # ndi/bincount rather than skimage.remove_small_* (whose kwarg semantics changed in 0.26 —
        # same version-robust convention as segment/nuclei._drop_small).
        holes = ndi.binary_fill_holes(mask) & ~mask
        hl, hn = ndi.label(holes)
        if hn:
            hsz = np.bincount(hl.ravel())
            small = np.where(hsz < 27)[0]
            mask = mask | np.isin(hl, small[small != 0])
        ol, on = ndi.label(mask)
        if on:
            osz = np.bincount(ol.ravel())
            big = np.where(osz >= 30)[0]
            mask = np.isin(ol, big[big != 0]) & interior
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
