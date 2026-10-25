"""Drift registration of a time-lapse series as a shift table (docs/DESIGN_timelapse.md section 4,
docs/TIMELAPSE.md "Registration"; ``germquant register``).

Nothing the stages see is ever resampled. ``register_series`` reads the downsampled reference volumes the
series loop kept next to each timepoint (``t<NNN>/<image_id>__ref_stride<k>.npy``, the segmentation channel
at xy stride k, and ``<image_id>__ref_brightfield_stride<k>.npy`` when a brightfield role resolved;
timelapse.series), estimates the translation between each timepoint and its reference frame with skimage's
``phase_cross_correlation`` (``upsample_factor`` 10 by default, so an estimate is quantised to a tenth of a
reference-volume voxel), accumulates the steps to the first timepoint and writes

* ``<stem>__register_shifts.csv`` and ``.parquet`` (schema.SERIES_TABLES["register_shifts"]): one row per
  timepoint with the per-frame step, the cumulative shift, the brightfield cross-check, the cap flag and the
  quality numbers;
* ``<stem>__register_report.json``: the parameters, a per-step summary, the capped steps, the brightfield
  disagreements and, with ``rigid=True``, the experimental Kabsch rotation per step.

Sign convention (the contract shared by the time-lapse packages): a positive shift means the sample moved by
that much in the raw frame, so a feature at position p at t0 sits at p + shift at t, and a registered
coordinate is the raw centroid minus the cumulative shift of its timepoint. ``phase_cross_correlation``
reports the opposite (the shift that brings the moving volume back onto the reference), so its result is
negated here. The tracking package computes registered coordinates from this table; this module writes
shifts only.

Normalization: by default the cross-power spectrum is NOT whitened (``normalization="none"``, skimage's
``normalization=None``, plain cross-correlation with the same upsampled peak refinement; the method column
says ``cross_correlation``). Nuclei at stride 4 to 8 are smooth blobs with little texture, and whitening
("phase", the method column's ``phase_cross_correlation``) then amplifies the empty high frequencies: on
the synthetic series of tests/test_timelapse_register.py it biased a step by up to 0.6 reference voxel
where the plain correlation was exact to the quantisation. ``phase`` stays available for textured data.

Quality: for every step, the median nearest-neighbour displacement of the nucleus centroids (the
per-timepoint nuclei tables; germline nuclei when the in_germline column says which) between the timepoint
and its reference frame, before and after the step is subtracted. It is a sanity number, not a loss: after
a good registration what remains is the real flow of the nuclei plus segmentation jitter.

Per-gap caps and the acquisition record (opt-in; timelapse.priors, docs/TIMELAPSE.md "Per-gap caps and the
acquisition record"). ``max_step_at`` sets the cap of single gaps and ``prior="stage"`` raises each gap's cap by
the stage and focus moves the nd2 recorded in that gap (the focus move as far as the estimate's own z uses it), and
flags gaps the record cannot explain. Without either the registration runs main's statements in main's order: the
same table columns and the same report keys.
"""
from __future__ import annotations

import json
import logging
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from .. import __version__, provenance, schema
from ..fsutil import long_path
from .series import BRIGHTFIELD_ROLE, SERIES_RECORD_SUFFIX, read_series_record, write_series_table

log = logging.getLogger(__name__)

TABLE = "register_shifts"
REPORT_SUFFIX = "__register_report.json"
NORMALIZATIONS = ("none", "phase")
METHODS = {"none": "cross_correlation", "phase": "phase_cross_correlation"}    # the method column per normalization
REFERENCE_MODES = ("previous", "t0")
DEFAULTS = {"reference": "previous", "upsample": 10, "max_step_um": 5.0, "disagreement_um": 1.0,
            "normalization": "none", "taper": 0.0, "link_radius_factor": 0.7, "min_links": 6}
CENTROID_COLS = ("centroid_z_um", "centroid_y_um", "centroid_x_um")
SEGMENTATION_ROLE_UNKNOWN = "segmentation"      # the main volume's role when no series record names it
_T_DIR = re.compile(r"t(\d{3})")
_REF_FILE = re.compile(r"__ref(?:_(?P<role>[A-Za-z0-9_]+?))?_stride(?P<k>\d+)\.npy$")


class ReferenceVolumesMissing(FileNotFoundError):
    """The series kept no reference volumes: timelapse.keep_reference_volume was off, or the series was
    processed before the volumes existed (then one ``germquant run --t all --resume`` pass reads them back)."""


# ---- the estimate ---------------------------------------------------------------------------------------

def displacement(ref: np.ndarray, mov: np.ndarray, upsample: int = DEFAULTS["upsample"],
                 normalization: str | None = DEFAULTS["normalization"], taper: float = DEFAULTS["taper"]) -> np.ndarray:
    """The displacement of the sample between ``ref`` and ``mov`` in voxels (z, y, x): positive when a
    feature sits at a higher index in ``mov`` than in ``ref``. skimage's phase_cross_correlation returns the
    shift that brings ``mov`` back onto ``ref``, so its result is negated. ``normalization`` "none" (the
    default; plain cross-correlation, see the module docstring) or "phase". ``taper`` > 0 subtracts each
    volume's mean and multiplies it by a Tukey window with that fraction of every axis tapered (0.25: the
    outer eighth on each side), which quiets content the field edge cuts, a stationary feature that pulls
    the estimate towards zero. ValueError when the two volumes differ in shape or the estimate is not
    finite."""
    from skimage.registration import phase_cross_correlation

    a = np.asarray(ref, dtype=np.float32)
    b = np.asarray(mov, dtype=np.float32)
    if a.shape != b.shape:
        raise ValueError(f"reference volume shapes differ: {a.shape} vs {b.shape} (was timelapse.reference_stride "
                         f"changed between runs of this series?)")
    if a.ndim != 3:
        raise ValueError(f"a reference volume must be (Z, Y, X), got shape {a.shape}")
    if taper and float(taper) > 0:
        a, b = _tapered(a, float(taper)), _tapered(b, float(taper))
    norm = None if normalization in (None, "none") else "phase"
    out = phase_cross_correlation(a, b, upsample_factor=max(1, int(upsample)), normalization=norm)
    shift = np.asarray(out[0] if isinstance(out, tuple) else out, dtype=float)
    if shift.shape != (3,) or not np.all(np.isfinite(shift)):
        raise ValueError("phase correlation returned no finite shift (a flat or empty volume?)")
    return -shift


def _tapered(a: np.ndarray, alpha: float) -> np.ndarray:
    """``a`` minus its mean, multiplied by a separable Tukey window (fraction ``alpha`` of each axis tapered)."""
    from scipy.signal.windows import tukey

    out = a - a.mean(dtype=np.float64).astype(np.float32)
    for ax, n in enumerate(out.shape):
        shape = [1] * out.ndim
        shape[ax] = n
        out = out * tukey(n, alpha=min(1.0, alpha)).astype(np.float32).reshape(shape)
    return out


def cap_step(step_um, max_step_um) -> tuple[np.ndarray, bool]:
    """The step shortened to ``max_step_um`` (its direction kept) when it is longer, with the capped flag;
    ``max_step_um`` None or <= 0 means no cap."""
    step = np.asarray(step_um, dtype=float)
    norm = float(np.linalg.norm(step))
    if max_step_um is not None and float(max_step_um) > 0 and norm > float(max_step_um):
        return step * (float(max_step_um) / norm), True
    return step, False


def nn_median(points, reference) -> float:
    """Median over ``points`` of the distance to the nearest point of ``reference`` (um); NaN when either is
    empty."""
    from scipy.spatial import cKDTree

    a, b = _points(points), _points(reference)
    if a is None or b is None:
        return float("nan")
    d, _ = cKDTree(b).query(a)
    return float(np.median(d))


def nn_spacing(points) -> float:
    """Median nearest-neighbour distance within one set of points (um); NaN with fewer than two points."""
    from scipy.spatial import cKDTree

    a = _points(points)
    if a is None or len(a) < 2:
        return float("nan")
    d, _ = cKDTree(a).query(a, k=2)
    return float(np.median(d[:, 1]))


def kabsch(P, Q) -> tuple[np.ndarray, np.ndarray]:
    """The proper rotation R and translation t minimising |R p + t - q| over the paired rows of P and Q."""
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    cp, cq = P.mean(axis=0), Q.mean(axis=0)
    H = (P - cp).T @ (Q - cq)
    U, _s, Vt = np.linalg.svd(H)
    d = 1.0 if np.linalg.det(Vt.T @ U.T) >= 0 else -1.0
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, cq - R @ cp


def rotation_angle_deg(R) -> float:
    return math.degrees(math.acos(max(-1.0, min(1.0, (float(np.trace(np.asarray(R))) - 1.0) / 2.0))))


def rigid_fit(moving_um, reference_um, *, link_radius_um: float, min_links: int = DEFAULTS["min_links"]) -> dict:
    """EXPERIMENTAL level 2 (docs/DESIGN_timelapse.md section 4): the Kabsch rigid fit between the translated
    moving centroids and the reference centroids, linked as mutual nearest neighbours within
    ``link_radius_um``. Returns n_linked, link_radius_um, rotation_deg, rotation_axis (unit vector, (z, y, x)
    order), rotation_matrix, translation_um (the residual translation after the rotation) and the RMS
    distance of the linked pairs before and after the fit; the fit fields are None with fewer than
    ``min_links`` links."""
    from scipy.spatial import cKDTree

    out = {"n_linked": 0, "link_radius_um": _f(link_radius_um), "rotation_deg": None, "rotation_axis": None,
           "rotation_matrix": None, "translation_um": None, "rmsd_before_um": None, "rmsd_after_um": None}
    A, B = _points(moving_um), _points(reference_um)
    if A is None or B is None or not (float(link_radius_um) > 0):
        return out
    dab, iab = cKDTree(B).query(A)
    _dba, iba = cKDTree(A).query(B)
    pairs = [(i, int(j)) for i, j in enumerate(iab) if iba[j] == i and dab[i] <= float(link_radius_um)]
    out["n_linked"] = len(pairs)
    if len(pairs) < max(3, int(min_links)):
        return out
    P = A[[i for i, _ in pairs]]
    Q = B[[j for _, j in pairs]]
    R, t = kabsch(P, Q)
    fit = P @ R.T + t
    w = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    nw = float(np.linalg.norm(w))
    out.update(rotation_deg=rotation_angle_deg(R), rotation_axis=(w / nw).round(6).tolist() if nw > 1e-12 else None,
               rotation_matrix=np.round(R, 9).tolist(), translation_um=np.round(t, 6).tolist(),
               rmsd_before_um=_rms(P - Q), rmsd_after_um=_rms(fit - Q))
    return out


# ---- what a series folder holds --------------------------------------------------------------------------

def discover_series(series_dir: str | Path) -> dict:
    """What a series folder holds: ``series_dir``, ``stem``, ``record`` (the ``<stem>__series.json`` content,
    None without one) and ``timepoints`` (one dict per timepoint in t_index order: t_index, image_id,
    t_seconds, folder, status, entry). The series record is used when present; ``t<NNN>/`` folders it does
    not list (or every folder, without a record) are read from their file set."""
    sdir = Path(series_dir)
    if not sdir.is_dir():
        raise FileNotFoundError(f"{sdir} is not a folder")
    records = sorted(sdir.glob(f"*{SERIES_RECORD_SUFFIX}"))
    if len(records) > 1:
        raise ValueError(f"{sdir} holds {len(records)} series records ({', '.join(p.name for p in records)}); "
                         f"one series per folder")
    record, stem = None, None
    if records:
        stem = records[0].name[: -len(SERIES_RECORD_SUFFIX)]
        record = read_series_record(sdir, stem)
        if record is None:
            raise ValueError(f"{records[0]} is not a readable series record")
        stem = str(record.get("series_id") or stem)
    tps: dict[int, dict] = {}
    for e in (record or {}).get("timepoints", []):
        try:
            k = int(e["t_index"])
        except (KeyError, TypeError, ValueError):
            continue
        iid = str(e.get("image_id") or f"{stem}_t{k:03d}")
        tps[k] = {"t_index": k, "image_id": iid, "t_seconds": e.get("t_seconds"),
                  "folder": sdir / str(e.get("out_dir") or f"t{k:03d}"), "status": e.get("status"), "entry": e}
    for p in sorted(sdir.iterdir()):
        m = _T_DIR.fullmatch(p.name)
        if not m or not p.is_dir():
            continue
        k = int(m.group(1))
        if k in tps:
            continue
        iid = _image_id_of(p)
        if iid is None:
            continue
        marker = provenance.read_done_marker(p, iid) or {}
        tps[k] = {"t_index": k, "image_id": iid, "t_seconds": marker.get("t_seconds"), "folder": p,
                  "status": "done" if marker else "unknown", "entry": None}
        if stem is None:
            stem = re.sub(r"_t\d{3}$", "", iid)
    if not tps:
        raise FileNotFoundError(f"{sdir}: no timepoint folders (t000/, t001/, ...) and no series record; run "
                                f"`germquant run SERIES.nd2 --t all --out {sdir}` first (docs/TIMELAPSE.md)")
    return {"series_dir": sdir, "stem": stem, "record": record, "timepoints": [tps[k] for k in sorted(tps)]}


def collect_reference_volumes(series_dir: Path, tps: list[dict]) -> dict:
    """The reference volumes of every timepoint at one common stride: ``stride`` (the stride most timepoints
    have, the smallest on a tie), ``role`` (the segmentation channel's role as the series record names it,
    'segmentation' without a record), ``files`` (one {role: Path} per timepoint, empty when that timepoint
    has none at the stride), ``spacing_um`` (dz, dy, dx of the volumes) and ``roles`` (every role seen).
    ReferenceVolumesMissing when no timepoint has any."""
    recs = [_record_volumes(series_dir, tp) for tp in tps]
    disks = [_disk_volumes(tp) for tp in tps]
    strides: Counter = Counter()
    for (rk, rfiles, _info), disk in zip(recs, disks):
        ks = set(disk)
        if rfiles and rk is not None:
            ks.add(rk)
        strides.update(ks)
    if not strides:
        raise ReferenceVolumesMissing(
            f"{series_dir}: no reference volumes (t<NNN>/<image_id>__ref_stride<k>.npy) next to the {len(tps)} "
            f"timepoint folder(s). `germquant run --t` keeps them when timelapse.keep_reference_volume is true (the "
            f"code default; k is timelapse.reference_stride when set, else chosen from timelapse.reference_pixel_um), "
            f"so this series was run with the key off "
            f"or before the volumes existed. Rerun `germquant run SERIES.nd2 --t all --out {series_dir} --resume` "
            f"with the key on: finished timepoints are not reprocessed, only their planes are read back.")
    stride = max(strides, key=lambda k: (strides[k], -k))
    role = next((str(info["role"]) for _rk, _f, info in recs if info.get("role")), SEGMENTATION_ROLE_UNKNOWN)
    files: list[dict] = []
    spacing = None
    roles: set[str] = set()
    for tp, (rk, rfiles, info), disk in zip(tps, recs, disks):
        if rfiles and rk == stride:
            f = dict(rfiles)
        else:
            f = {(role if r is None else r): p for r, p in disk.get(stride, {}).items()}
        files.append(f)
        roles.update(f)
        sp = info.get("spacing_um") if (rk == stride and f) else None
        if sp is None and f:
            sp = _spacing_from_summary(tp, stride)
        if sp is not None:
            sp = [float(v) for v in sp]
            if spacing is None:
                spacing = sp
            elif not np.allclose(sp, spacing, rtol=1e-6, atol=0.0):
                log.warning("t%03d: reference volume voxel size %s differs from %s", tp["t_index"], sp, spacing)
    if spacing is None:
        raise ValueError(f"{series_dir}: cannot tell the voxel size of the reference volumes (no reference_volumes "
                         f"entry in the series record and no <image_id>__image_summary.csv with voxel_* columns)")
    return {"stride": stride, "role": role, "files": files, "spacing_um": spacing, "roles": roles}


def load_centroids(folder: str | Path, image_id: str) -> np.ndarray | None:
    """The (N, 3) raw-frame centroids (z, y, x um) of ``<image_id>__nuclei.csv``: the germline nuclei when the
    in_germline column marks any, else every nucleus; None without a table or without a finite centroid."""
    f = Path(folder) / f"{image_id}__nuclei.csv"
    if not f.is_file():
        return None
    try:
        df = pd.read_csv(long_path(f), low_memory=False)
    except Exception as e:  # noqa: BLE001 - a damaged table is a missing quality number, not a failed registration
        log.warning("%s: unreadable nuclei table (%s)", image_id, e)
        return None
    if df.empty or any(c not in df.columns for c in CENTROID_COLS):
        return None
    if "in_germline" in df.columns:
        m = np.asarray([str(v).strip().lower() == "true" for v in df["in_germline"].tolist()], dtype=bool)
        if m.any():
            df = df[m]
    pts = df[list(CENTROID_COLS)].to_numpy(dtype=float)
    pts = pts[np.all(np.isfinite(pts), axis=1)]
    return pts if len(pts) else None


# ---- the command ----------------------------------------------------------------------------------------

def register_series(
    series_dir: str | Path,
    *,
    reference: str = DEFAULTS["reference"],
    upsample: int = DEFAULTS["upsample"],
    max_step_um: float | None = DEFAULTS["max_step_um"],
    channel: str | None = None,
    rigid: bool = False,
    disagreement_um: float = DEFAULTS["disagreement_um"],
    normalization: str = DEFAULTS["normalization"],
    taper: float = DEFAULTS["taper"],
    max_step_at: dict | None = None,
    prior: str | None = None,
    nd2: str | Path | None = None,
) -> dict:
    """Register a finished series folder (the ``--out`` of ``germquant run SERIES.nd2 --t ...``) and write the
    register_shifts table and the registration report next to its series record (module docstring).

    ``reference`` "previous" registers every timepoint against the previous registrable one and accumulates
    the steps; "t0" registers each directly against the first timepoint with a volume (the step is then the
    difference of consecutive estimates). ``upsample`` is phase_cross_correlation's upsample_factor. A step
    longer than ``max_step_um`` is shortened to it, its direction kept, and flagged ``capped`` (None or 0: no
    cap). ``channel`` picks the role whose volume is registered (default: the segmentation channel's);
    the other kept role is the cross-check, whose per-step estimate is compared with the primary one and
    flagged when they differ by more than ``disagreement_um``. ``rigid`` adds the experimental Kabsch
    rotation per step to the report (the table is unchanged). ``normalization`` "none" (plain
    cross-correlation, the default; module docstring) or "phase" (whitened). ``taper`` (0, off, by default)
    is the Tukey fraction `displacement` applies to both volumes, for a field whose edge cuts bright tissue.

    Opt-in (module docstring): ``max_step_at`` ``{T: um}`` caps the gap into timepoint T at that length instead
    (0: no cap there; T must be a registrable timepoint after the first). ``prior="stage"`` reads the series'
    acquisition record (priors.load_stage_record: ``nd2``, else the series record's nd2, else
    ``<stem>__stage_record.json``), raises each gap's cap by the recorded stage xy move and by as much of the stack
    z move as the step's own z uses (priors.gap_cap) and flags the gaps the record cannot explain (priors.gap_flags;
    advisory). With either, the
    table gains ``cap_um`` and ``cap_source`` (and, with the prior, the gap prior's columns) after the declared
    ones, each ``steps`` entry its cap and gap prior, and the report the given options and a ``prior`` block.

    Returns ``series_id``, ``series_dir``, ``table`` (the DataFrame, declared columns), ``table_path`` (the
    CSV; the parquet copy sits next to it), ``report`` and ``report_path``. ReferenceVolumesMissing (a
    FileNotFoundError) when the series kept no reference volumes; ValueError for a bad option."""
    mode = str(reference).strip().lower()
    if mode not in REFERENCE_MODES:
        raise ValueError(f"reference must be one of {REFERENCE_MODES}, got {reference!r}")
    norm = "none" if normalization is None else str(normalization).strip().lower()
    if norm not in NORMALIZATIONS:
        raise ValueError(f"normalization must be one of {NORMALIZATIONS}, got {normalization!r}")
    method = METHODS[norm]
    tap = float(taper or 0.0)
    if not 0.0 <= tap <= 1.0:
        raise ValueError(f"taper must be a fraction between 0 and 1, got {taper!r}")
    up = int(upsample)
    if up < 1:
        raise ValueError(f"upsample must be a positive integer, got {upsample!r}")
    cap = float(max_step_um) if max_step_um is not None else 0.0
    thr = float(disagreement_um)
    pri = _prior_option(prior, nd2)
    per_gap = bool(max_step_at) or pri is not None

    found = discover_series(series_dir)
    sdir, stem, record, tps = found["series_dir"], found["stem"], found["record"], found["timepoints"]
    vols = collect_reference_volumes(sdir, tps)
    stride, main_role, files, roles = vols["stride"], vols["role"], vols["files"], vols["roles"]
    spacing = np.asarray(vols["spacing_um"], dtype=float)
    primary = main_role if channel is None else str(channel).strip()
    if primary not in roles:
        raise ValueError(f"no reference volume for channel role {primary!r}; this series kept {sorted(roles)} "
                         f"(the segmentation channel is <image_id>__ref_stride<k>.npy, a brightfield channel "
                         f"<image_id>__ref_brightfield_stride<k>.npy)")
    if primary != BRIGHTFIELD_ROLE and BRIGHTFIELD_ROLE in roles:
        cross: str | None = BRIGHTFIELD_ROLE
    elif primary == BRIGHTFIELD_ROLE and main_role in roles:
        cross = main_role
    else:
        cross = None
    if per_gap:
        from . import priors as PR

        msa = PR.check_max_step_at(max_step_at, [int(tp["t_index"]) for j, tp in enumerate(tps)
                                                 if primary in files[j]][1:])
        stage = (PR.load_stage_record(sdir, stem, record, nd2=nd2, t_indices=[int(tp["t_index"]) for tp in tps])
                 if pri == "stage" else None)
        z_range = (record or {}).get("z_range")
        gap_rows: dict[int, dict] = {}       # t_index -> the per-gap columns of that row
        gap_list: list[dict] = []            # the registered gaps, for the report's prior block

    cache: dict[tuple[int, str], np.ndarray] = {}

    def load(i: int, role: str) -> np.ndarray:
        key = (i, role)
        if key not in cache:
            cache[key] = np.load(long_path(files[i][role]))
        return cache[key]

    cents: dict[int, np.ndarray | None] = {}

    def centroids(i: int) -> np.ndarray | None:
        if i not in cents:
            cents[i] = load_centroids(tps[i]["folder"], tps[i]["image_id"])
        return cents[i]

    n = len(tps)
    cum: list = [None] * n          # cumulative shift (um; z, y, x) per timepoint, None without an estimate
    est: list = [None] * n          # the uncapped cumulative estimate (t0 mode), so a capped step is not carried into the next
    cum_bf: list = [None] * n
    origin = last = None            # the timepoint every shift is accumulated to; the last registrable one
    rows: list[dict] = []
    steps: list[dict] = []
    rigid_rows: list[dict] = []
    warnings: list[str] = []
    capped_list: list[int] = []
    disagree_list: list[int] = []
    errors: dict[str, str] = {}

    for i, tp in enumerate(tps):
        k = int(tp["t_index"])
        row = _row(tp, primary, stride, up, method)
        if primary not in files[i]:
            why = "no reference volume"
            if tp.get("status") not in (None, "done"):
                why += f" (timepoint {tp['status']})"
            row["error"] = why
            errors[str(k)] = why
            rows.append(row)
            continue
        if origin is None:
            origin = last = i
            cum[i] = np.zeros(3)
            est[i] = np.zeros(3)
            has_bf = bool(cross) and cross in files[i]
            cum_bf[i] = np.zeros(3) if has_bf else None
            _fill(row, ref_k=None, step=np.zeros(3), cum=np.zeros(3),
                  bf_step=np.zeros(3) if has_bf else None, bf_cum=cum_bf[i])
            rows.append(row)
            continue
        ref_i = last if mode == "previous" else origin
        ref_k = int(tps[ref_i]["t_index"])
        try:
            disp_um = displacement(load(ref_i, primary), load(i, primary), up, norm, tap) * spacing
        except Exception as e:  # noqa: BLE001 - one bad volume is a row with an error, not a failed table
            why = f"{type(e).__name__}: {e}"
            row["error"] = why
            errors[str(k)] = why
            log.warning("t%03d: %s", k, why)
            rows.append(row)
            continue
        raw_step = disp_um if mode == "previous" else disp_um - est[last]
        cap_k = cap
        if per_gap:              # the gap from the last registrable frame: its own cap (priors.gap_cap)
            gp = PR.gap_prior(stage, int(tps[last]["t_index"]), k, z_range) if stage is not None else None
            cap_k, cap_src = PR.gap_cap(k, cap, max_step_at=msa, gp=gp, est_dz_um=float(raw_step[0]))
        step, capped = cap_step(raw_step, cap_k)
        base = cum[ref_i] if mode == "previous" else cum[last]
        cum_i = base + step
        cum[i] = cum_i
        est[i] = disp_um if mode != "previous" else cum_i
        applied = cum_i - cum[ref_i]          # the displacement between this timepoint and its reference frame
        # the cross-check: the same estimate on the other channel's volumes, never capped
        bf_step = bf_cum = None
        disagreement = float("nan")
        flag = False
        if cross and cross in files[i] and cross in files[ref_i]:
            try:
                bf_disp = displacement(load(ref_i, cross), load(i, cross), up, norm, tap) * spacing
                if mode == "previous":
                    bf_base = cum_bf[ref_i] if cum_bf[ref_i] is not None else cum[ref_i]
                    bf_step, bf_cum = bf_disp, bf_base + bf_disp
                else:
                    prev_bf = cum_bf[last] if cum_bf[last] is not None else cum[last]
                    bf_step, bf_cum = bf_disp - prev_bf, bf_disp
                disagreement = float(np.linalg.norm(bf_step - raw_step))
                flag = disagreement > thr
            except Exception as e:  # noqa: BLE001 - the cross-check is advisory
                warnings.append(f"t{k:03d}: {cross} cross-check failed: {type(e).__name__}: {e}")
        cum_bf[i] = bf_cum
        # quality: the centroids of this timepoint against the reference frame's, raw and shifted
        c_mov, c_ref = centroids(i), centroids(ref_i)
        q_before = nn_median(c_mov, c_ref)
        q_after = nn_median(c_mov - applied, c_ref) if c_mov is not None else float("nan")
        if rigid and c_mov is not None and c_ref is not None:
            fit = rigid_fit(c_mov - applied, c_ref, link_radius_um=DEFAULTS["link_radius_factor"] * nn_spacing(c_ref),
                            min_links=DEFAULTS["min_links"])
            rigid_rows.append({"t_index": k, "reference_t_index": ref_k, **fit})
        _fill(row, ref_k=ref_k, step=step, cum=cum_i, capped=capped, bf_step=bf_step, bf_cum=bf_cum,
              disagreement=disagreement, flag=flag, q_before=q_before, q_after=q_after)
        steps.append({"t_index": k, "reference_t_index": ref_k, "step_um": _l(step),
                      "raw_step_um": _l(raw_step) if capped else None, "cumulative_um": _l(cum_i), "capped": capped,
                      "brightfield_step_um": _l(bf_step), "brightfield_cumulative_um": _l(bf_cum),
                      "disagreement_um": _f(disagreement), "disagreement_flag": flag,
                      "quality_before_um": _f(q_before), "quality_after_um": _f(q_after)})
        if per_gap:
            flags = PR.gap_flags(raw_step, gp, cap, rigid=True)
            gap_rows[k] = PR.gap_row(cap_k, cap_src, gp, flags, with_prior=stage is not None)
            steps[-1].update(cap_um=_f(gap_rows[k]["cap_um"]), cap_source=cap_src)
            if stage is not None:
                steps[-1]["prior"] = PR.gap_summary(gp, flags)
                gap_list.append({**gp, "flags": flags, "cap_um": _f(gap_rows[k]["cap_um"]), "cap_source": cap_src})
        if capped:
            capped_list.append(k)
        if flag:
            disagree_list.append(k)
        for key in [key for key in cache if key[0] not in (i, origin)]:
            del cache[key]
        last = i
        rows.append(row)

    df = pd.DataFrame(rows, columns=list(schema.SERIES_TABLES[TABLE]))
    df["reference_t_index"] = pd.array([None if pd.isna(v) else int(v) for v in df["reference_t_index"]], dtype="Int64")
    df["t_seconds"] = pd.to_numeric(df["t_seconds"], errors="coerce")
    for c in ("capped", "disagreement_flag"):
        df[c] = df[c].astype(bool)
    if per_gap:
        PR.add_gap_columns(df, gap_rows, with_prior=stage is not None)
    table_path = write_series_table(df, TABLE, sdir, stem, _series_meta(sdir, stem, record))

    before = [s["quality_before_um"] for s in steps if s["quality_before_um"] is not None]
    after = [s["quality_after_um"] for s in steps if s["quality_after_um"] is not None]
    dis = [s["disagreement_um"] for s in steps if s["disagreement_um"] is not None]
    norms = [float(np.linalg.norm(s["step_um"])) for s in steps]
    report = {
        "series_id": stem, "series_dir": str(sdir), "method": method,
        "parameters": {"reference": mode, "upsample_factor": up, "normalization": norm, "taper": tap,
                       "max_step_um": cap if cap > 0 else None,
                       "disagreement_um": thr, "channel_role": primary, "cross_check_role": cross, "stride": stride,
                       "spacing_um": _l(spacing), "rigid": bool(rigid),
                       "origin_t_index": int(tps[origin]["t_index"]) if origin is not None else None},
        "n_timepoints": n, "n_registered": len(steps), "n_capped": len(capped_list), "n_disagreement": len(disagree_list),
        "capped": capped_list, "brightfield_disagreement": disagree_list, "errors": errors, "warnings": warnings,
        "max_step_um": max(norms) if norms else None,
        "median_disagreement_um": float(np.median(dis)) if dis else None,
        "quality": {"n_steps": len(before), "median_before_um": float(np.median(before)) if before else None,
                    "median_after_um": float(np.median(after)) if after else None,
                    "max_after_um": float(np.max(after)) if after else None},
        "cumulative_last_um": _l(cum[last]) if last is not None else None,
        "last_t_index": int(tps[last]["t_index"]) if last is not None else None,
        "steps": steps,
        "rigid": rigid_rows if rigid else None,
        "experimental": ["rigid"] if rigid else [],
        "pipeline_version": __version__, "git_sha": provenance.git_sha(), "run_timestamp": provenance.run_timestamp(),
    }
    if per_gap:
        PR.add_report_options(report, msa, pri, stage, gap_list, z_range=z_range)
    report_path = sdir / f"{stem}{REPORT_SUFFIX}"
    with open(long_path(report_path), "w", encoding="utf-8") as fh:
        json.dump(_clean(report), fh, indent=1)
    log.info("%s: %d timepoint(s), %d step(s) registered, %d capped, %d brightfield disagreement(s) -> %s",
             stem, n, len(steps), len(capped_list), len(disagree_list), table_path.name)
    return {"series_id": stem, "series_dir": str(sdir), "table": df, "table_path": str(table_path),
            "report": report, "report_path": str(report_path)}


# ---- helpers ------------------------------------------------------------------------------------------

def _prior_option(prior, nd2) -> str | None:
    """The ``prior`` option checked: None or one of priors.PRIORS; ``nd2`` only goes with a prior."""
    if prior is None:
        if nd2 is not None:
            raise ValueError("nd2 is read for its stage record only: give it with prior='stage' (--nd2 with "
                             "--prior stage)")
        return None
    from .priors import PRIORS

    pri = str(prior).strip().lower()
    if pri not in PRIORS:
        raise ValueError(f"prior must be one of {PRIORS} (or None), got {prior!r}")
    return pri


def _row(tp: dict, role: str, stride: int, up: int, method: str) -> dict:
    nan = float("nan")
    r = {c: nan for c in schema.SERIES_TABLES[TABLE]}
    r.update(t_index=int(tp["t_index"]), t_seconds=tp.get("t_seconds"), image_id=tp["image_id"],
             reference_t_index=None, channel_role=role, disagreement_flag=False, capped=False,
             xy_stride=int(stride), upsample_factor=int(up), method=method, error="")
    return r


def _fill(row: dict, *, ref_k, step, cum, capped=False, bf_step=None, bf_cum=None, disagreement=float("nan"),
          flag=False, q_before=float("nan"), q_after=float("nan")) -> None:
    row["reference_t_index"] = ref_k
    for ax, s, c in zip("zyx", step, cum):
        row[f"shift_{ax}_um"] = row[f"step_d{ax}_um"] = float(s)
        row[f"cum_shift_{ax}_um"] = row[f"d{ax}_um"] = float(c)
    if bf_step is not None:
        for ax, s in zip("zyx", bf_step):
            row[f"bf_shift_{ax}_um"] = float(s)
    if bf_cum is not None:
        for ax, c in zip("zyx", bf_cum):
            row[f"brightfield_d{ax}_um"] = float(c)
    row["bf_disagreement_um"] = row["brightfield_disagreement_um"] = float(disagreement)
    row["disagreement_flag"] = bool(flag)
    row["capped"] = bool(capped)
    row["nn_disp_before_um"] = row["quality_before_um"] = float(q_before)
    row["nn_disp_after_um"] = row["quality_after_um"] = float(q_after)


def _image_id_of(folder: Path) -> str | None:
    for suffix in ("__image_summary.csv", "__done.json", "__stages.json"):
        hits = sorted(folder.glob(f"*{suffix}"))
        if hits:
            return hits[0].name[: -len(suffix)]
    hits = sorted(folder.glob("*__ref_stride*.npy"))
    if hits:
        return hits[0].name.split("__ref_stride")[0]
    return None


def _record_volumes(series_dir: Path, tp: dict) -> tuple[int | None, dict[str, Path], dict]:
    """(stride, {role: Path present on disk}, the entry) from the series record's reference_volumes entry."""
    entry = tp.get("entry") or {}
    info = entry.get("reference_volumes") if isinstance(entry, dict) else None
    if not isinstance(info, dict) or not isinstance(info.get("files"), dict):
        return None, {}, {}
    files: dict[str, Path] = {}
    for role, rel in info["files"].items():
        p = series_dir / str(rel)
        if p.is_file():
            files[str(role)] = p
        else:
            log.warning("%s: the reference volume %s listed in the series record is missing", tp["image_id"], rel)
    try:
        k = int(info.get("stride"))
    except (TypeError, ValueError):
        k = None
    return k, files, info


def _disk_volumes(tp: dict) -> dict[int, dict[str | None, Path]]:
    """{stride: {role or None (the segmentation channel): Path}} from the file names in the timepoint folder."""
    out: dict[int, dict[str | None, Path]] = {}
    folder = tp["folder"]
    if not folder.is_dir():
        return out
    for p in sorted(folder.glob(f"{tp['image_id']}__ref*_stride*.npy")):
        m = _REF_FILE.search(p.name)
        if m:
            out.setdefault(int(m.group("k")), {})[m.group("role")] = p
    return out


def _spacing_from_summary(tp: dict, stride: int) -> list[float] | None:
    f = tp["folder"] / f"{tp['image_id']}__image_summary.csv"
    if not f.is_file():
        return None
    try:
        s = pd.read_csv(long_path(f), nrows=1, low_memory=False)
        dz, dy, dx = (float(s[c].iloc[0]) for c in ("voxel_dz_um", "voxel_dy_um", "voxel_dx_um"))
    except Exception:  # noqa: BLE001
        return None
    return [dz, dy * int(stride), dx * int(stride)]


def _series_meta(sdir: Path, stem: str, record: dict | None) -> dict:
    """The SERIES_META block of the table: the series record and the series_summary table supply the design
    columns and the run's config hash; version, git sha and timestamp are this registration's."""
    meta: dict = {k: None for k in schema.SERIES_META}
    meta["series_id"] = stem
    rec = record or {}
    for k in ("file_path", "n_timepoints", "period_s", "config_hash"):
        meta[k] = rec.get(k)
    f = sdir / f"{stem}__series_summary.csv"
    if f.is_file():
        try:
            s = pd.read_csv(long_path(f), nrows=1, low_memory=False)
            if len(s):
                for k in schema.SERIES_META:
                    if k in ("pipeline_version", "git_sha", "run_timestamp") or meta.get(k) is not None:
                        continue
                    if k in s.columns and pd.notna(s[k].iloc[0]):
                        v = s[k].iloc[0]
                        meta[k] = v.item() if hasattr(v, "item") else v
        except Exception as e:  # noqa: BLE001 - the design columns are a convenience
            log.warning("%s: series summary unreadable (%s)", stem, e)
    meta.update(pipeline_version=__version__, git_sha=provenance.git_sha(), run_timestamp=provenance.run_timestamp())
    return meta


def _points(p) -> np.ndarray | None:
    if p is None:
        return None
    a = np.asarray(p, dtype=float)
    if a.size == 0:
        return None
    a = a.reshape(-1, 3)
    return a if len(a) else None


def _rms(d: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum(np.asarray(d, dtype=float) ** 2, axis=1))))


def _f(v) -> float | None:
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def _l(v) -> list[float] | None:
    if v is None:
        return None
    return [round(float(x), 6) for x in np.asarray(v, dtype=float).ravel()]


def _clean(obj):
    """A JSON-ready copy: numpy scalars and arrays to Python, NaN and infinities to null."""
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [_clean(v) for v in obj.tolist()]
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        x = float(obj)
        return x if math.isfinite(x) else None
    return obj
