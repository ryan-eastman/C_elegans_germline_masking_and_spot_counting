"""Single-image pipeline: .nd2 -> tidy tables + montage + label mask + manifest.

Stages, in the order of `germquant.stages.STAGES`: read -> segment nuclei (Cellpose) -> measure ->
signal QC -> isolate germline -> linearize axis -> spots (SpotMAX) -> granules (PGL-1) -> coloc -> qc ->
render -> write. Optional stages are switched by one config key each (read with a code default, so config files
never need editing) and run through `stages.run_stage`, which flags a failure (`<stage>:FAILED_<Exc>`)
and continues instead of killing the batch. Voxel spacing from the .nd2 is threaded into every 3D op.
A per-stage outcome record is written next to the tables as ``<image_id>__stages.json``.
"""
from __future__ import annotations

# --- cuBLAS load-order guard (Windows / RTX 5090) --------------------------------------------
# torch (ships cuBLAS 12.8) and CuPy (pulls cuBLAS 12.9 via the nvidia-*-cu12 wheels) each carry
# their own cublas64_12.dll. Whichever is imported first claims that DLL name for the whole
# process. If CuPy wins, torch's bf16/fp32 *batched* GEMMs fail on Blackwell (sm_120) with
# CUBLAS_STATUS_INVALID_VALUE, and segment_nuclei silently drops from Cellpose-SAM to the
# classical watershed fallback (different, worse masks). Importing torch HERE — before any
# SpotMAX/CuPy import happens during a run — makes torch's cuBLAS load first and win.
try:
    import torch  # noqa: F401  (imported for side effect: claim cublas64_12.dll before CuPy)
except Exception:  # CPU-only / torch-not-installed smoke test: watershed fallback still runs
    pass

import json
import logging
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import provenance, qc, schema
from .axis import linearize_germline
from .config import Config
from .fsutil import long_path
from .io import parse_sample, read_image, read_image_metadata, read_nd2_metadata, read_stack
from .io.acquisition import expected_flags, expected_settings, read_acquisition
from .io.formats import image_format
from .measure import measure_objects
from .render import make_montage
from .segment import segment_nuclei
from .stages import STAGES, Outcome, missing_roles, run_stage, stage_enabled
from .stages import stage_hashes as stage_hashes_fn

log = logging.getLogger(__name__)

# envelope stage parameters: config key suffix -> code default (the August analysis constants)
_ENVELOPE_KEYS = {
    "band_um": 2.0, "smooth_um": 0.15, "volume_ratio_min": 1.0, "volume_ratio_max": 3.0,
    "ring_smooth_um": 0.25, "ring_shell_um": 0.4, "ring_ratio_max": 0.97, "shell_over_thr_max": 0.75,
    "territory_dilate_um": 4.0, "cyto_um": 2.5,
}


def spot_instances(cfg) -> list[dict]:
    """The extra spot instances a config declares (spots.instances), validated once: each needs a
    non-empty identifier-like `name` (its table is spots_<name>, its nuclei column n_spots_<name>) that
    is unique and collides with no core table or column, and a channel `role`. Raises ValueError for a
    bad profile before any image is read."""
    raw = cfg.get("spots.instances", []) or []
    out, seen = [], set()
    for i, inst in enumerate(raw):
        if not isinstance(inst, dict):
            raise ValueError(f"spots.instances[{i}] must be a mapping with name and role")
        name = str(inst.get("name", "")).strip()
        if not name or not name.replace("_", "a").isalnum() or name in seen:
            raise ValueError(f"spots.instances[{i}]: name {name!r} must be a unique identifier")
        if f"spots_{name}" in schema.TABLES or f"n_spots_{name}" in schema.NUCLEI:
            raise ValueError(f"spots.instances[{i}]: name {name!r} collides with a core table or column")
        role = str(inst.get("role", "")).strip()
        if not role:
            raise ValueError(f"spots.instances[{i}] ({name}): a channel role is required")
        seen.add(name)
        out.append({**inst, "name": name, "role": role})
    return out


def spot_backend_section(cfg, inst: dict | None = None) -> dict:
    """The raw ``spots:`` mapping a spot detection backend receives (docs/BACKENDS.md), read with a code
    default of an empty mapping, with a spot instance's own keys overlaid; a mapping the instance carries
    (``log:``) merges one level deep over the shared one, so an instance can change one key of it."""
    section = dict(cfg.get("spots", {}) or {})
    for k, v in (inst or {}).items():
        base = section.get(k)
        section[k] = {**base, **v} if isinstance(v, dict) and isinstance(base, dict) else v
    return section


def resolve_model_path(cfg):
    """The Cellpose model the config names, resolved against the config's base_dir when it is a relative
    path to an existing file (so a profile outside the repo still finds models/...), or against
    ``$GERMQUANT_MODELS`` (the HPC / container binding: the config's ``models/...`` path is looked up
    under it, with or without the leading ``models/``); bare built-in names such as 'cpsam' pass through."""
    model_path = cfg.get("segmentation.nuclei.cellpose_model")
    if not model_path:
        return model_path
    p = Path(str(model_path))
    if p.is_absolute():
        return model_path
    root = os.environ.get("GERMQUANT_MODELS")
    if root:
        for cand in (Path(root) / p, Path(root) / Path(*p.parts[1:]) if len(p.parts) > 1 and p.parts[0] == "models" else None):
            if cand is not None and cand.is_file():
                return cand
    if (cfg.base_dir / p).is_file():
        return cfg.base_dir / p
    return model_path


def read_input_metadata(path: str | Path) -> dict:
    """The header of the input image through this module's reader names (``read_nd2_metadata`` for an
    ``.nd2``, io.read_image_metadata otherwise): the seam the synthetic-stack test fixtures patch, so the
    time-lapse loop (timelapse.series) learns the timepoint count from the same reader process_image uses."""
    if image_format(path) == "nd2":
        return read_nd2_metadata(path)
    return read_image_metadata(path)


def _read_input(path: Path, *, xy_stride: int, z_range, t=None):
    """Header and pixels of the input image: ``(metadata dict, Stack)``.

    An ``.nd2`` goes through this module's ``read_nd2_metadata`` and ``read_stack`` names: the functions
    every golden was made with, and the seam the synthetic-stack test fixtures patch. Every other format
    goes through io.read_image_metadata / io.read_image (bioio, optional ``readers`` extra), which return
    the same Stack and honour xy_stride and z_range the same way. ``t`` (one timepoint of a time series)
    is passed on only when given, so a fixed image makes exactly the calls it always has."""
    kw = {} if t is None else {"t": t}
    if image_format(path) == "nd2":
        return read_nd2_metadata(path, **kw), read_stack(path, xy_stride=xy_stride, z_range=z_range, **kw)
    return read_image_metadata(path, **kw), read_image(path, xy_stride=xy_stride, z_range=z_range, **kw)


def _t_seconds_of(meta: dict, t: int) -> float | None:
    """The recorded time stamp of timepoint ``t`` from the header's ``t_seconds`` list, or None when the
    file carries none (bioio formats, or an unreadable per-frame record)."""
    stamps = meta.get("t_seconds")
    if stamps is None or not 0 <= int(t) < len(stamps):
        return None
    return float(stamps[int(t)])


def t_seconds_of(meta: dict, t: int) -> float | None:
    """The acquisition time (s) of timepoint ``t`` from the reader's ``t_seconds`` list, or None when the
    header carries none (a fixed image, or a reader without the field). Public name of ``_t_seconds_of``
    for the time-lapse loop; tolerant of a header that is not a dict."""
    if not isinstance(meta, dict):
        return None
    try:
        return _t_seconds_of(meta, t)
    except (IndexError, TypeError, ValueError):
        return None


def process_image(
    nd2_path: str | Path,
    cfg: Config,
    out_dir: str | Path,
    *,
    xy_stride: int = 1,
    z_range: tuple[int, int] | None = None,
    prov: dict | None = None,
    t: int | None = None,
) -> dict:
    """Run every enabled stage on one image and write its results folder. ``nd2_path`` (the name is kept
    for compatibility) is any supported image: an ``.nd2``, or a TIFF / OME-TIFF, CZI, LIF or OME-Zarr
    file when the ``readers`` extra is installed (io/read.py). Returns image_id, n_nuclei, qc_pass,
    qc_flags, out_dir, the tables (DataFrames by name) and the per-stage record.

    ``t`` processes one timepoint of a time series (docs/TIMELAPSE.md): the reader loads that timepoint's
    stack, every stage runs on it unchanged, the image id (and with it every per-image file name) becomes
    ``<stem>_t<NNN>``, and ``t_index`` / ``t_seconds`` are recorded in the run geometry (run_manifest.json,
    ``__stages.json``) and the completion marker. With ``t=None``, the default, nothing changes: a fixed
    image's id, tables, stage record and marker are exactly what they were."""
    nd2_path = Path(nd2_path)
    out_dir = Path(out_dir)
    os.makedirs(long_path(out_dir), exist_ok=True)
    flags: list[str] = []
    outcomes: dict[str, Outcome] = {}      # per-stage record -> <image_id>__stages.json

    # ---- read + resolve channels ----
    t_read = time.perf_counter()
    acq_expected = expected_settings(cfg)   # a bad acquisition.expected section fails here, before any pixel is read ({} when absent)
    meta, stack = _read_input(nd2_path, xy_stride=xy_stride, z_range=z_range, t=t)
    role_to_idx, ch_flags = cfg.channel_map.resolve(meta["channel_names"])
    flags += ch_flags
    spacing = stack.spacing
    if not stack.spacing_ok:
        # voxel size unreadable -> spacing is a 1 µm isotropic guess; every physical
        # length/area/volume is unreliable. Surface it loudly rather than silently.
        flags.append("voxel:UNREADABLE_assumed_isotropic_1um")
    dapi_present = role_to_idx.get("dna") is not None

    sample = parse_sample(nd2_path, cfg.get("metadata.filename_regex"), cfg.get("metadata.defaults"))
    if t is not None:
        # one timepoint of a series: the design columns were parsed from the raw stem above; the image id,
        # and so every per-image file name, carries the timepoint (three digits, one underscore)
        sample["image_id"] = f"{sample['image_id']}_t{int(t):03d}"
    # a rerun must never leave an older run's completion marker behind (batch --resume trusts it)
    provenance.clear_done_marker(out_dir, sample["image_id"])
    # per-stage provenance (sub-hashes, model digest, run geometry): additive fields on the manifest and
    # the stage record, never table columns, so no CSV header changes.
    spot_instances(cfg)                    # a bad spots.instances declaration fails here, before any pixel is read
    model_file = resolve_model_path(cfg)
    model_sha = (provenance.file_sha256(model_file)
                 if str(cfg.get("segmentation.nuclei.method", "auto")) != "classical" else "unused")
    stage_hashes = stage_hashes_fn(cfg, role_to_idx, provenance.tool_versions(), model_sha)
    run_geometry = {"xy_stride": int(xy_stride), "z_range": list(z_range) if z_range else None,
                    "stage_hashes": stage_hashes, "model_sha256": model_sha,
                    "enabled_stages": [s.name for s in STAGES if stage_enabled(cfg, s.name)]}
    if t is not None:
        # the timepoint enters the run geometry (manifest, stage record) and the marker ONLY for a series
        run_geometry["t_index"] = int(t)
        run_geometry["t_seconds"] = _t_seconds_of(meta, t)
    prov = prov or provenance.write_manifest(out_dir, config_hash=cfg.hash, config=cfg.as_dict(),
                                             extra=run_geometry)

    shared = {
        "image_id": sample["image_id"], "file_path": str(nd2_path),
        "genotype": sample["genotype"], "sex": sample["sex"], "germ_cell": sample["germ_cell"],
        "treatment": sample["treatment"], "replicate": sample["replicate"],
        "acquisition_date": sample["acquisition_date"],
        "voxel_dz_um": spacing[0], "voxel_dy_um": spacing[1], "voxel_dx_um": spacing[2],
        "n_channels": stack.n_channels, "dapi_present": dapi_present,
        "channel_map_json": json.dumps({r: role_to_idx[r] for r in role_to_idx}),
        "pipeline_version": prov["pipeline_version"], "git_sha": prov["git_sha"],
        "config_hash": cfg.hash, "run_timestamp": prov["run_timestamp"],
    }
    outcomes["read"] = Outcome("ran", elapsed_s=time.perf_counter() - t_read, flags=list(flags))

    # ---- acquisition metadata (exposure / laser power / filters / camera per role, from the nd2 description) ----
    # With acquisition.expected set (docs/QC.md), the record is also checked against the study's expected
    # settings and the advisory acq:<role>_* flags are raised; the section is absent from every shipped
    # config, so the default path writes the same flags it always has.
    acq_fields: dict = {}

    def _acquisition():
        rec = read_acquisition(nd2_path, list(meta["channel_names"]), role_to_idx)
        if acq_expected:
            flags.extend(expected_flags(rec.per_role, acq_expected))
        return rec.fields

    acq_on = stage_enabled(cfg, "acquisition")
    is_nd2 = image_format(nd2_path) == "nd2"      # the description text lives in the .nd2 header only
    if acq_expected and not (acq_on and is_nd2):
        log.warning("acquisition.expected is set but the acquisition stage is not running (%s); the "
                    "expected-settings check is skipped",
                    "acquisition.read_exposures is off" if not acq_on else "not an .nd2")
    acq_res = run_stage("acquisition", _acquisition, flags=flags, outcomes=outcomes, enabled=acq_on and is_nd2,
                        skip_reason="disabled" if not acq_on else "not an .nd2 (no acquisition description)")
    if acq_res is not None:
        acq_fields = acq_res

    # ---- segment nuclei (always on; a failure here is fatal, as it always was) ----
    def _segment():
        dna = stack.channel(role_to_idx.get("dna"))
        if dna is None:
            flags.append("segment:no_dna_channel_using_first_channel")
            dna = stack.data[0]
        seg = cfg.segmentation.nuclei
        return segment_nuclei(
            dna, spacing,
            method=seg.get("method", "auto"), cellpose_model=str(model_file) if model_file else "cpsam",
            diameter_um=float(seg.get("diameter_um", 3.0)), min_volume_um3=float(seg.get("min_volume_um3", 4.0)),
            device=seg.get("device"),          # cuda | mps | cpu | auto; $GERMQUANT_DEVICE overrides
            # the segmenter: built-in cellpose / classical or an installed plugin (docs/BACKENDS.md), read
            # with a code default so the config files stay untouched; the section goes along for the
            # backend's own keys
            backend=cfg.get("segmentation.nuclei.backend", "cellpose"),
            cfg_section=cfg.get("segmentation.nuclei", {}),
        )

    labels, seg_method = run_stage("segment", _segment, flags=flags, outcomes=outcomes, fatal=True)
    n_nuclei = int(labels.max())

    # ---- measure nuclei ----
    def _measure():
        intensity = {r: stack.data[i] for r, i in role_to_idx.items() if i is not None}
        df = measure_objects(labels, intensity, spacing, compute_surface=False)
        return df.rename(columns={"label": "nucleus_id"}) if not df.empty else pd.DataFrame(columns=["nucleus_id"])

    nuclei = run_stage("measure", _measure, flags=flags, outcomes=outcomes, fatal=True)

    # ---- signal QC: is each channel measurable, and does one leak into another? (germquant.signal_qc) ----
    # Runs on an xy sub-sample right after measure, before any stage trusts a channel. It adds only NEW
    # image_summary columns (sig_<role>_*) and its own table (channel_crosstalk), and its flags are off
    # unless a profile enables them, so the qc_flags text of the validated default path (and its goldens)
    # stays byte-identical. The acquisition stage's exposure / power fields normalize the levels.
    signal_summary_fields: dict = {}
    crosstalk_table = pd.DataFrame(columns=schema.CHANNEL_CROSSTALK)

    def _signal_qc():
        from .signal_qc import DEFAULTS as SQ
        from .signal_qc import run_signal_qc, signal_flags

        res = run_signal_qc(
            stack.data, labels, spacing, role_to_idx, acq_fields=acq_fields, image_id=sample["image_id"],
            stride=int(cfg.get("signal_qc.stride", SQ["stride"])),
            tissue_dilation_um=float(cfg.get("signal_qc.tissue_dilation_um", SQ["tissue_dilation_um"])),
            ring_um=float(cfg.get("signal_qc.ring_um", SQ["ring_um"])),
            min_bg_voxels=int(cfg.get("signal_qc.min_bg_voxels", SQ["min_bg_voxels"])),
            bit_depth=cfg.get("signal_qc.bit_depth"), source_dtype=meta.get("dtype"),
        )
        if bool(cfg.get("signal_qc.flags.enabled", False)):
            flags.extend(signal_flags(
                res["summary"], res["crosstalk"], res["roles"],
                min_snr=float(cfg.get("signal_qc.min_snr", SQ["min_snr"])),
                max_sat_frac=float(cfg.get("signal_qc.max_sat_frac", SQ["max_sat_frac"])),
                max_r=float(cfg.get("signal_qc.crosstalk.max_r", SQ["crosstalk_max_r"]))))
        return res

    sig_on = stage_enabled(cfg, "signal_qc")
    sig_roles = [r for r, i in role_to_idx.items() if i is not None]
    sig_res = run_stage("signal_qc", _signal_qc, flags=flags, outcomes=outcomes,
                        enabled=sig_on and bool(sig_roles),
                        skip_reason="disabled" if not sig_on else "no resolved channel role")
    if sig_res is not None:
        signal_summary_fields = dict(sig_res["summary"])
        crosstalk_table = sig_res["crosstalk"]

    # ---- isolate germline (drop nuclei segmented OUTSIDE the gonad: gut, debris, off-gonad) ----
    # Uses SYP (central_element) intensity + spatial connectivity, NEVER the spot count. Excluded
    # nuclei are re-attached (flagged in_germline=False) before writing, so nothing is hidden;
    # downstream means (axis, spots/nucleus) operate on the germline subset.
    excluded = nuclei.iloc[0:0].copy()

    def _germline():
        from .germline import select_germline

        df, germ_flags = select_germline(
            nuclei,
            method=cfg.get("germline.method", "syp_seeded_cc"),
            syp_percentile=float(cfg.get("germline.syp_percentile", 25.0)),
            link_radius_um=float(cfg.get("germline.link_radius_um", 12.0)),
            min_seed_frac=float(cfg.get("germline.min_seed_frac", 0.10)),
            size_frac=float(cfg.get("germline.size_frac", 0.10)),
        )
        flags.extend(germ_flags)
        return df[~df["in_germline"]].copy(), df[df["in_germline"]].copy()

    germ_res = run_stage("germline", _germline, flags=flags, outcomes=outcomes, fatal=True,
                         enabled=stage_enabled(cfg, "germline") and not nuclei.empty,
                         skip_reason="disabled" if not stage_enabled(cfg, "germline") else "no nuclei")
    if germ_res is not None:
        excluded, nuclei = germ_res
    elif not stage_enabled(cfg, "germline") and not nuclei.empty:
        # germline isolation switched off (config/profiles/generic.yaml, --no-germline): every segmented
        # nucleus is the analysis set, so in_germline is a real True on every row rather than the NA an
        # unrun stage would leave (schema.py). The default path (germline on) is untouched.
        nuclei = nuclei.assign(in_germline=True)
    n_germline_nuclei = int(len(nuclei))

    # ---- lamin envelope masks (seeded watershed on LMN-1), ring test, territories ----
    env_ctx = None
    env_summary_fields: dict = {}
    lamin_idx = role_to_idx.get("lamin")

    def _envelope():
        nonlocal nuclei
        from .envelope import run_envelope

        ids = [int(v) for v in nuclei["nucleus_id"].tolist()]
        params = {k: cfg.get(f"envelope.{k}", v) for k, v in _ENVELOPE_KEYS.items()}
        ctx = run_envelope(labels, ids, stack.data[lamin_idx], spacing, params)
        nuclei = nuclei.merge(ctx["per_nucleus"], on="nucleus_id", how="left")
        flags.append(f"envelope:fallback_n={ctx['summary']['n_envelope_fallback']}")
        if ctx["summary"]["n_no_envelope"]:
            flags.append(f"envelope:no_envelope_n={ctx['summary']['n_no_envelope']}")
        return ctx

    env_on = stage_enabled(cfg, "envelope")
    env_ctx = run_stage("envelope", _envelope, flags=flags, outcomes=outcomes,
                        enabled=env_on and lamin_idx is not None and n_germline_nuclei > 0,
                        skip_reason=("disabled" if not env_on else
                                     "no lamin channel" if lamin_idx is None else "no germline nuclei"))
    if env_ctx is not None:
        env_summary_fields = dict(env_ctx["summary"])

    # ---- mask audit (lamin-only nuclei missed by the labels; labels with no envelope) ----
    audit_table = pd.DataFrame(columns=schema.MASK_AUDIT)
    audit_summary_fields: dict = {}
    audit_ctx = None

    def _audit():
        from .qc_audit import run_audit

        ids = [int(v) for v in nuclei["nucleus_id"].tolist()]
        return run_audit(labels, ids, stack.data[lamin_idx], spacing, env_ctx,
                         ring_smooth_um=float(cfg.get("envelope.ring_smooth_um", 0.25)),
                         missed_cov_max=float(cfg.get("audit.missed_cov_max", 0.15)))

    audit_on = stage_enabled(cfg, "audit")
    audit_ctx = run_stage("audit", _audit, flags=flags, outcomes=outcomes,
                          enabled=audit_on and env_ctx is not None,
                          skip_reason="disabled" if not audit_on else "envelope stage did not run")
    if audit_ctx is not None:
        audit_table = audit_ctx["table"]
        audit_summary_fields = dict(audit_ctx["summary"])
        flags.append(f"audit:missed_nuclei={audit_ctx['summary']['n_missed_nuclei']}")

    # ---- linearize axis (principal-curve centerline -> per-nucleus distal->proximal position) ----
    axis_conf = float("nan")

    def _axis():
        df, conf, axis_flags = linearize_germline(
            nuclei, confidence_min=float(cfg.get("axis.qc_confidence_min", 0.6)))
        flags.extend(axis_flags)
        return df, conf

    axis_res = run_stage("axis", _axis, flags=flags, outcomes=outcomes, fatal=True,
                         enabled=stage_enabled(cfg, "axis") and not nuclei.empty,
                         skip_reason="disabled" if not stage_enabled(cfg, "axis") else "no germline nuclei")
    if axis_res is not None:
        nuclei, axis_conf = axis_res

    # ---- pachytene staging from the hand-traced axis (staging.traces_file) ----
    zones_table = pd.DataFrame(columns=schema.ZONES)
    staging_summary_fields: dict = {}
    stag_ctx = None
    stag_on = stage_enabled(cfg, "staging")

    def _staging():
        nonlocal nuclei
        from .staging import load_traces, run_staging

        # the traces file is read inside the stage so a missing or malformed file flags staging:FAILED
        tf = cfg.get("staging.traces_file", "staging/pachytene_traces.json")
        tf_path = Path(tf) if Path(tf).is_absolute() else cfg.base_dir / tf
        trace = load_traces(tf_path).get(sample["image_id"])
        if not isinstance(trace, dict) or "points_um" not in trace:
            flags.append("staging:no_trace")
            return None
        if trace.get("status") != "traced":
            flags.append(f"staging:trace_status_{trace.get('status')}")
            return None
        df = nuclei
        use_env = str(cfg.get("staging.centroid", "envelope")) == "envelope" and "envelope_centroid_x_um" in df.columns \
            and df["envelope_centroid_x_um"].notna().any()
        work = df.copy()
        if use_env:                       # project the envelope centroids, as the August analysis did
            for ax_ in ("x", "y"):
                work[f"centroid_{ax_}_um"] = work[f"envelope_centroid_{ax_}_um"].fillna(work[f"centroid_{ax_}_um"])
        res = run_staging(work, trace, spacing, env_ctx=env_ctx,
                          off_axis_um=cfg.get("staging.off_axis_um"), adaptive=bool(cfg.get("staging.adaptive_cutoff", True)))
        nuclei = df.merge(res["zones"], on="nucleus_id", how="left")
        flags.append(f"staging:zoned_n={res['summary']['n_zoned_nuclei']},centroid={'envelope' if use_env else 'dapi'}")
        return res

    stag_ctx = run_stage("staging", _staging, flags=flags, outcomes=outcomes,
                         enabled=stag_on and not nuclei.empty,
                         skip_reason="disabled" if not stag_on else "no germline nuclei")
    if stag_ctx is not None:
        zones_table = stag_ctx["zones"]
        staging_summary_fields = dict(stag_ctx["summary"])
    elif stag_on and outcomes.get("staging") is not None and outcomes["staging"].status == "ran":
        # ran but found no usable trace for this image: record it as a skip with the reason
        why = next((f for f in flags if f.startswith(("staging:no_trace", "staging:trace_status_"))), "no trace")
        outcomes["staging"] = Outcome("skipped", reason=why.replace("staging:", "").replace("_", " "))

    # ---- spots: RAD-51 (or other) foci per nucleus, through the spot detection backend that
    # spots.backend names (code default spotmax, the validated SpotMAX chain; docs/BACKENDS.md).
    # SpotMAX detects peaks ABOVE local background inside each nucleus mask, merges z-axis spot-splits,
    # and tags every spot with its effect size. Detection params are cross-validated vs Imaris (config). ----
    spots = pd.DataFrame(columns=schema.SPOTS)
    per_nuc_spots = None  # kept so off-gonad (excluded) nuclei also get n_spots at re-attach
    spots_idx = role_to_idx.get("foci")

    def _spots():
        from .spots import detect_spots

        backend = cfg.get("spots.backend", "spotmax")
        per_spot, per_nuc = detect_spots(
            stack.data[spots_idx], labels, spacing,
            marker=cfg.channel_map.marker("foci"),
            spot_radius_um=float(cfg.get("spots.spot_radius_um", 0.3)),
            gauss_sigma_um=float(cfg.get("spots.gauss_sigma_um", 0.08)),
            thresholding_method=cfg.get("spots.thresholding_method", "threshold_triangle"),
            effect_size_metric=cfg.get("spots.effect_size_metric", "spot_vs_backgr_effect_size_glass"),
            effect_size_min=float(cfg.get("spots.effect_size_min", 3.0)),
            merge_z_columns=bool(cfg.get("spots.merge_z_columns", True)),
            z_merge_gap_um=float(cfg.get("spots.z_merge_gap_um", 0.8)),
            z_merge_valley_frac=float(cfg.get("spots.z_merge_valley_frac", 0.8)),
            max_spot_candidates=int(cfg.get("spots.max_spot_candidates", 30000)),
            backend=backend, cfg_section=spot_backend_section(cfg),
        )
        df = nuclei
        if not df.empty and not per_nuc.empty:
            df = df.merge(per_nuc[["nucleus_id", "n_spots"]], on="nucleus_id", how="left")
            df["n_spots"] = df["n_spots"].fillna(0).astype(int)
        flags.append(f"spots:{backend}_n={len(per_spot)}")          # the default run still writes spots:spotmax_n=
        if backend != "spotmax":
            flags.append(f"spots:unvalidated_backend={backend}")  # only spotmax is cross-validated vs Imaris
        return per_spot, per_nuc, df

    spots_on = stage_enabled(cfg, "spots")
    spots_res = run_stage("spots", _spots, flags=flags, outcomes=outcomes,
                          enabled=spots_on and spots_idx is not None and n_nuclei > 0,
                          skip_reason=("disabled" if not spots_on else
                                       "no foci channel" if spots_idx is None else "no nuclei"))
    if spots_res is not None:
        spots, per_nuc_spots, nuclei = spots_res

    # ---- further spot instances (COSA-1 crossover foci, ...): the same detector on another channel role.
    # `spots.instances` is read with a code default of none, so the RAD-51 path above is untouched. Each
    # instance gets its own table (spots_<name>), nuclei column (n_spots_<name>), flag and summary; an
    # optional late-pachytene restriction (needs the staging stage) and an expected count per nucleus
    # (6 oocyte / 5 spermatocyte bivalents) turn it into a crossover-designation readout. ----
    extra_spot_tables: dict[str, pd.DataFrame] = {}
    extra_spot_cols: list[str] = []
    extra_spot_per_nuc: dict[str, pd.DataFrame] = {}
    spot_instance_summary: dict = {}
    inst_enabled: list[bool] = []
    for inst in spot_instances(cfg):
        name, role = inst["name"], inst["role"]
        table_name, col = f"spots_{name}", f"n_spots_{name}"
        extra_spot_tables[table_name] = pd.DataFrame(columns=schema.SPOTS)
        extra_spot_cols.append(col)
        r_idx = role_to_idx.get(role)

        def _inst(inst=inst, name=name, role=role, table_name=table_name, col=col, r_idx=r_idx):
            nonlocal nuclei
            from .io.sample_metadata import expected_sc_count
            from .spots import detect_spots

            def p(key, default):
                return inst.get(key, cfg.get(f"spots.{key}", default))

            backend = p("backend", "spotmax")
            per_spot, per_nuc = detect_spots(
                stack.data[r_idx], labels, spacing, marker=cfg.channel_map.marker(role),
                spot_radius_um=float(p("spot_radius_um", 0.3)), gauss_sigma_um=float(p("gauss_sigma_um", 0.08)),
                thresholding_method=p("thresholding_method", "threshold_triangle"),
                effect_size_metric=p("effect_size_metric", "spot_vs_backgr_effect_size_glass"),
                effect_size_min=float(p("effect_size_min", 3.0)), merge_z_columns=bool(p("merge_z_columns", True)),
                z_merge_gap_um=float(p("z_merge_gap_um", 0.8)), z_merge_valley_frac=float(p("z_merge_valley_frac", 0.8)),
                max_spot_candidates=int(p("max_spot_candidates", 30000)),
                backend=backend, cfg_section=spot_backend_section(cfg, inst),
            )
            df = nuclei
            if not df.empty:
                counts = per_nuc.set_index("nucleus_id")["n_spots"] if not per_nuc.empty else pd.Series(dtype=float)
                df = df.copy()
                df[col] = df["nucleus_id"].map(counts).fillna(0).astype(int)
            summ = {f"mean_{col}": float(df[col].mean()) if (not df.empty and col in df) else float("nan")}
            zone = inst.get("restrict_to_zone")
            if zone and "zone" in df.columns:
                sel = df[df["zone"] == zone]
                summ[f"mean_{col}_{zone}"] = float(sel[col].mean()) if len(sel) else float("nan")
                summ[f"n_nuclei_{zone}"] = int(len(sel))
            if inst.get("expected_from_germ_cell", True):
                summ[f"{name}_expected_per_nucleus"] = expected_sc_count(sample["germ_cell"]) or float("nan")
            nuclei = df
            flags.append(f"spots_{name}:{backend}_n={len(per_spot)}")
            if backend != "spotmax":
                flags.append(f"spots_{name}:unvalidated_backend={backend}")
            return per_spot, per_nuc, summ

        on_i = spots_on and r_idx is not None and n_nuclei > 0
        inst_enabled.append(on_i)
        res_i = run_stage(f"spots_{name}", _inst, flags=flags, outcomes=outcomes, enabled=on_i,
                          skip_reason=("disabled" if not spots_on else
                                       f"no {role} channel" if r_idx is None else "no nuclei"))
        if res_i is not None:
            extra_spot_tables[table_name], extra_spot_per_nuc[col], summ = res_i
            spot_instance_summary.update(summ)

    # ---- SC tracing (per-nucleus SC length, fragment lower bound, fragmentation index) ----
    # The June 2026 tracer restored verbatim (germquant.sc.skeleton.trace_sc), run on the GERMLINE
    # nuclei only (the label image is masked to them; the tracer traces every label it is given).
    # Off by default (sc.trace.enabled); independent of coloc and of any granule channel, so a
    # 3-channel DAPI/SYP/RAD-51 image yields spots and SC readouts in one run. docs/SC_TRACING.md.
    sc_tracks = pd.DataFrame(columns=schema.SC_TRACKS)
    sc_per_nuc = pd.DataFrame(columns=schema.SC_PER_NUCLEUS)
    sc_summary_fields: dict = {}
    ce_idx = role_to_idx.get("central_element")

    def _sc_trace():
        from .io.sample_metadata import expected_sc_count
        from .sc import skan_available, trace_sc

        if not skan_available():
            flags.append("sc_trace:skan_missing")
            return sc_tracks, sc_per_nuc, nuclei, {}
        germ_ids = [int(v) for v in nuclei["nucleus_id"].tolist()]
        germ_labels = np.where(np.isin(labels, germ_ids), labels, 0)
        exp = expected_sc_count(sample["germ_cell"])
        tracks, per_nuc = trace_sc(
            stack.data[ce_idx], germ_labels, spacing,
            marker=cfg.channel_map.marker("central_element"),
            ridge_sigmas_um=tuple(cfg.get("sc.ridge_sigmas_um", [0.15, 0.25, 0.40])),
            ridge_hyst_low_pct=float(cfg.get("sc.ridge_hyst_low_pct", 45.0)),
            ridge_hyst_high_pct=float(cfg.get("sc.ridge_hyst_high_pct", 80.0)),
            intensity_percentile=float(cfg.get("sc.intensity_percentile", 90.0)),
            min_fragment_length_um=float(cfg.get("sc.min_fragment_length_um", 0.5)),
            expected_n_tracks={nid: exp for nid in germ_ids} if exp else None,
        )
        df = nuclei
        if not per_nuc.empty:
            cols = per_nuc[["nucleus_id", "sc_total_length_um", "n_fragments", "sc_fragmentation_index",
                            "expected_n_tracks"]].rename(columns={"n_fragments": "sc_n_fragments_lb",
                                                                  "expected_n_tracks": "sc_expected_n_tracks"})
            df = df.merge(cols, on="nucleus_id", how="left")
        summary = {
            "mean_sc_total_length_um": float(per_nuc["sc_total_length_um"].mean()) if len(per_nuc) else float("nan"),
            "mean_sc_fragmentation_index": float(per_nuc["sc_fragmentation_index"].mean()) if len(per_nuc) else float("nan"),
            "mean_sc_n_fragments_lb": float(per_nuc["n_fragments"].mean()) if len(per_nuc) else float("nan"),
        }
        flags.append(f"sc_trace:n_traced={int((per_nuc['n_fragments'] > 0).sum()) if len(per_nuc) else 0}")
        flags.append("sc:uncalibrated")           # absolute lengths await the Imaris Filament calibration
        return tracks, per_nuc, df, summary

    sc_on = stage_enabled(cfg, "sc_trace")
    sc_res = run_stage("sc_trace", _sc_trace, flags=flags, outcomes=outcomes,
                       enabled=sc_on and ce_idx is not None and n_germline_nuclei > 0,
                       skip_reason=("disabled" if not sc_on else
                                    "no central_element channel" if ce_idx is None else "no germline nuclei"))
    if sc_res is not None:
        sc_tracks, sc_per_nuc, nuclei, sc_summary_fields = sc_res

    # ---- p-granule (PGL-1) surfacing, then SYP<->PGL-1 colocalization ----
    # Granules are surfaced as 3D objects inside the germline dilated by a perinuclear shell (P-granules
    # sit just OUTSIDE the nuclear envelope, so the region MUST include the perinuclear cytoplasm or the
    # overlap reads ~0 by construction). The coloc stage then measures their DIRECT voxel/object overlap
    # with two SYP operands: `syp_aggregate` = cytoplasmic SYP blobs in the shell (the headline for
    # P-granule coincidence) and `sc_ribbon` = the intranuclear SC ribbon (control), plus the shell
    # voxel Pearson / partition coefficient headline. Both stages degrade gracefully.
    granules = pd.DataFrame(columns=schema.GRANULES)
    coloc = pd.DataFrame(columns=schema.COLOC)
    masks: dict = {}
    granule_labels = None
    per_nuc_granules = None
    coloc_summary_fields: dict = {}
    gran_ctx = None

    def _granule():
        nonlocal nuclei
        ctx = _run_granule(stack, labels, role_to_idx, nuclei, spacing, cfg, env_ctx=env_ctx)
        if ctx["method"] != "threshold_triangle" or ctx["region_kind"] != "coloc_region":
            # only a non-default recipe is flagged, so the qc_flags text of the validated default path
            # (and its goldens) stays byte-identical
            flags.append(f"granule:method={ctx['method']},region={ctx['region_kind']}")
        per_nuc = ctx["per_nuc"]
        if per_nuc is not None and not nuclei.empty:
            df = nuclei
            if not per_nuc.empty:
                df = df.merge(per_nuc, on="nucleus_id", how="left")
            if "n_granules" not in df.columns:
                df["n_granules"] = 0
                df["granule_volume_um3"] = 0.0
            df["n_granules"] = df["n_granules"].fillna(0).astype(int)
            df["granule_volume_um3"] = df["granule_volume_um3"].fillna(0.0)
            nuclei = df
        return ctx

    gran_on = stage_enabled(cfg, "granule")
    gran_missing = missing_roles(role_to_idx, "granule")
    gran_ctx = run_stage("granule", _granule, flags=flags, outcomes=outcomes,
                         enabled=gran_on and not gran_missing and n_germline_nuclei > 0,
                         skip_reason=("disabled" if not gran_on else
                                      f"missing channel role(s) {gran_missing}" if gran_missing
                                      else "no germline nuclei"))
    if gran_ctx is not None:
        granules, granule_labels, per_nuc_granules = gran_ctx["granules"], gran_ctx["granule_labels"], gran_ctx["per_nuc"]

    def _coloc():
        c, m = _run_coloc(stack, labels, role_to_idx, spacing, cfg, gran_ctx)
        flags.append(f"coloc:granules_n={len(granules)}")
        return c, m

    coloc_on = stage_enabled(cfg, "coloc")
    g_out = outcomes.get("granule")
    coloc_res = run_stage("coloc", _coloc, flags=flags, outcomes=outcomes,
                          enabled=coloc_on and gran_ctx is not None,
                          skip_reason=("disabled" if not coloc_on else
                                       "granule stage failed" if g_out is not None and g_out.status == "failed" else
                                       f"granule stage skipped: {g_out.reason}" if g_out is not None and g_out.status == "skipped"
                                       else "granule stage did not run"))
    if coloc_res is not None:
        coloc, masks = coloc_res
        granules = gran_ctx["granules"]          # now carries the per-operand overlap columns
        # Contract: a coloc FAILURE keeps the granule stage's products (granules table, label image,
        # n_granules on nuclei) but adds no coloc fields to image_summary, so its shape matches a run
        # without coloc; consumers key off __stages.json / the coloc:FAILED_ flag, not on n_granules.
        coloc_summary_fields = _coloc_summary(coloc, granules)

    # ---- partition coefficient (SYP-3 into P granules) and the per-granule lit fraction ----
    # Both work in the padded germline crop of the August analysis; the envelope stage supplies the
    # envelope, else the DAPI nuclei stand in (flagged). The granule mask is the granule stage's.
    partition_table = pd.DataFrame(columns=schema.PARTITION)
    partition_summary: dict = {}
    tail_table = pd.DataFrame(columns=schema.GRANULE_TAIL)
    tail_summary: dict = {}

    def _pc_frame():
        """(crop slice, envelope labels crop, envelope mask crop, source name)."""
        from .envelope import germline_crop

        if env_ctx is not None:
            sl = env_ctx["crop"]
            env_lab = env_ctx["envelope_labels"][sl]
            return sl, env_lab, env_lab > 0, "envelope"
        ids = [int(v) for v in nuclei["nucleus_id"].tolist()]
        sl = germline_crop(labels, ids)
        lab_c = np.where(np.isin(labels[sl], ids), labels[sl], 0).astype(np.int32)
        return sl, lab_c, lab_c > 0, "dapi"

    def _partition():
        from .partition import run_partition

        sl, env_lab, env_mask, src = _pc_frame()
        pre = gran_ctx.get("cyto_ctx") if gran_ctx is not None else None
        if pre is not None and pre["crop"] != sl:
            pre = None
        res = run_partition(gran_ctx["syp"][sl], gran_ctx["granule_labels"][sl] > 0, env_mask, spacing,
                            cyto_um=float(cfg.get("envelope.cyto_um", 2.5)),
                            bg_percentile=float(cfg.get("partition.bg_percentile", 3.0)),
                            env_labels_c=env_lab, zones=zones_table if len(zones_table) else None,
                            dz=int(cfg.get("partition.zshift_planes", 6)), precomputed=pre,
                            whole_min_cyto_voxels=int(cfg.get("partition.whole_min_cyto_voxels", 5000)),
                            whole_min_granule_objects=int(cfg.get("partition.whole_min_granule_objects", 40)))
        flags.append(f"partition:frame={src},granules={gran_ctx['method']}")
        if res.get("gated"):
            flags.append(f"partition:whole_gated_{res['gated']}")
        return res

    part_on = stage_enabled(cfg, "partition")
    part_res = run_stage("partition", _partition, flags=flags, outcomes=outcomes,
                         enabled=part_on and gran_ctx is not None,
                         skip_reason="disabled" if not part_on else "granule stage did not run")
    if part_res is not None:
        partition_table, partition_summary = part_res["table"], dict(part_res["summary"])

    def _tail():
        from .granule.segment import TOPHAT_DEFAULTS
        from .partition import run_granule_tail

        sl, env_lab, _env_mask, src = _pc_frame()
        drop = env_ctx["no_envelope_ids"] if env_ctx is not None else []
        res = run_granule_tail(gran_ctx["syp"][sl], gran_ctx["pgl"][sl], env_lab, spacing, drop_ids=drop,
                               cyto_um=float(cfg.get("envelope.cyto_um", 2.5)),
                               tophat_kw={k: float(cfg.get(f"granule.tophat.{k}", v)) for k, v in TOPHAT_DEFAULTS.items()})
        flags.append(f"granule_tail:frame={src},dropped={len(drop)}")
        return res

    tail_on = stage_enabled(cfg, "granule_tail")
    tail_res = run_stage("granule_tail", _tail, flags=flags, outcomes=outcomes,
                         enabled=tail_on and gran_ctx is not None,
                         skip_reason="disabled" if not tail_on else "granule stage did not run")
    if tail_res is not None:
        tail_table, tail_summary = tail_res["table"], dict(tail_res["summary"])

    # ---- per-image summary fields ----
    # Every stage's image_summary contribution, spread in the order the row has always used (declared
    # columns are reordered to the schema on write; the per-instance extras keep this order). A stage that
    # needs another stage's number reads it here: the partition measurability check below reads the
    # signal_qc stage's SYP signal-to-noise (sig_central_element_snr) from this dict.
    summary_fields: dict = {**coloc_summary_fields, **sc_summary_fields, **env_summary_fields, **audit_summary_fields,
                            **acq_fields, **staging_summary_fields, **partition_summary, **tail_summary,
                            **spot_instance_summary, **signal_summary_fields}

    # ---- partition measurability (partition_syp_snr, partition_measurable; docs/PARTITION_QC.md) ----
    # A SYP channel with no live fluorophore, or recorded through the wrong emission filter, yields a
    # partition coefficient floor above 1 that no rotation null removes. When a profile sets
    # partition.min_syp_snr (code default 0.0: off) an image below it is flagged partition:low_syp_snr;
    # every partition number is still written, the flag is the signal.
    if part_res is not None:
        from .partition import partition_measurability

        meas_fields, meas_flags = partition_measurability(summary_fields,
                                                          min_syp_snr=cfg.get("partition.min_syp_snr", 0.0))
        summary_fields.update(meas_fields)
        flags.extend(meas_flags)
        outcomes["partition"].flags.extend(meas_flags)

    # ---- QC ----
    t_qc = time.perf_counter()
    qc_pass, qc_all = qc.qc_flags(
        n_nuclei=n_nuclei, channel_flags=ch_flags,
        axis_flags=[f for f in flags if f.startswith("axis")],
        spots_found=len(spots) > 0 or any(len(t) > 0 for t in extra_spot_tables.values()),
        spots_enabled=(bool(cfg.get("spots.enabled", True)) and spots_idx is not None) or any(inst_enabled),
    )
    qc_all = sorted(set(flags + qc_all))
    outcomes["qc"] = Outcome("ran", elapsed_s=time.perf_counter() - t_qc,
                             flags=[f for f in qc_all if f.startswith("qc:")])

    image_summary = pd.DataFrame([{
        "n_nuclei": n_nuclei,
        "n_germline_nuclei": n_germline_nuclei,
        "mean_spots": float(nuclei["n_spots"].mean()) if "n_spots" in nuclei else float("nan"),
        "total_germline_length_um": float(nuclei["axis_position_um"].max())
        if "axis_position_um" in nuclei and not nuclei.empty else float("nan"),
        "qc_pass": qc_pass, "qc_flags": ";".join(qc_all),
        "segmentation_method": seg_method, "axis_confidence": axis_conf,
        **summary_fields,
    }])

    # ---- write outputs ----
    # re-attach the off-gonad nuclei (in_germline=False) so the table is a complete, auditable record
    # of every segmented object — analysis filters on in_germline.
    if not excluded.empty:
        # give off-gonad nuclei their spot counts too, so n_spots is complete (spots may land in
        # debris/gut nuclei); without this they'd read NaN and per-spot vs per-nucleus totals diverge.
        if per_nuc_spots is not None and not per_nuc_spots.empty:
            m = per_nuc_spots.set_index("nucleus_id")["n_spots"]
            excluded = excluded.copy()
            excluded["n_spots"] = excluded["nucleus_id"].map(m).fillna(0).astype(int)
        # off-gonad nuclei are outside the coloc region, so their p-granule load is a real zero
        # (never NaN, so per-nucleus vs per-object granule totals reconcile) — mirror the spots re-attach.
        if per_nuc_granules is not None:
            excluded = excluded.copy()
            mg = per_nuc_granules.set_index("nucleus_id") if not per_nuc_granules.empty else None
            excluded["n_granules"] = (
                excluded["nucleus_id"].map(mg["n_granules"]) if mg is not None else 0)
            excluded["n_granules"] = excluded["n_granules"].fillna(0).astype(int)
            excluded["granule_volume_um3"] = (
                excluded["nucleus_id"].map(mg["granule_volume_um3"]) if mg is not None else 0.0)
            excluded["granule_volume_um3"] = excluded["granule_volume_um3"].fillna(0.0)
        # off-gonad nuclei carry their real count for every extra spot instance too (the per-spot table
        # is a complete record, so the nuclei column must reconcile with it), like n_spots above
        for col in extra_spot_cols:
            pn = extra_spot_per_nuc.get(col)
            if pn is not None and not pn.empty:
                excluded[col] = excluded["nucleus_id"].map(pn.set_index("nucleus_id")["n_spots"]).fillna(0).astype(int)
        nuclei = pd.concat([nuclei, excluded], ignore_index=True)
        for col in extra_spot_cols:            # an instance that did not run: a real zero everywhere
            if col in nuclei.columns:
                nuclei[col] = nuclei[col].fillna(0).astype(int)
    tables = {"nuclei": nuclei, "spots": spots, "granules": granules,
              "coloc": coloc, "image_summary": image_summary,
              "sc_tracks": sc_tracks, "sc_per_nucleus": sc_per_nuc, "mask_audit": audit_table,
              "zones": zones_table, "partition": partition_table, "granule_tail": tail_table,
              "channel_crosstalk": crosstalk_table,
              **extra_spot_tables}

    def _write():
        _write_tables(tables, shared, out_dir, sample["image_id"], cfg.get("output.formats", ["csv"]))
        if cfg.get("output.write_label_images", True):
            _save_labels(labels, out_dir / f"{sample['image_id']}__nuclei_labels.tif", spacing)
        if cfg.get("output.write_spots_image", True) and len(spots):
            _save_spots_image(spots, labels.shape, spacing, out_dir / f"{sample['image_id']}__spots.tif",
                              radius_um=float(cfg.get("spots.spot_radius_um", 0.3)))
        # surfaced objects for Imaris: PGL-1 granules as a label image (-> Surfaces), the SC ribbon and
        # cytoplasmic SYP-aggregate masks as calibrated binary TIFs (overlay them on the raw channels).
        if cfg.get("output.write_label_images", True) and granule_labels is not None and granule_labels.max() > 0:
            _save_labels(granule_labels, out_dir / f"{sample['image_id']}__granules_labels.tif", spacing)
        for mname in ("sc_ribbon", "syp_aggregate"):
            m = masks.get(mname)
            if m is not None and m.any():
                _save_mask_image(m, spacing, out_dir / f"{sample['image_id']}__{mname}.tif")
        if env_ctx is not None and cfg.get("output.write_envelope_labels", True):
            _save_labels(env_ctx["envelope_labels"], out_dir / f"{sample['image_id']}__envelope_labels.tif", spacing)

    def _audit_overlay():
        from .qc_audit import write_audit_overlay

        sl = env_ctx["crop"]
        germ = nuclei[nuclei["in_germline"].astype(bool)] if "in_germline" in nuclei.columns else nuclei
        ids = [int(v) for v in germ["nucleus_id"].tolist()]     # the germline set the audit table used
        write_audit_overlay(out_dir / f"{sample['image_id']}__mask_audit.png",
                            stack.data[role_to_idx["dna"]][sl] if role_to_idx.get("dna") is not None else stack.data[0][sl],
                            stack.data[lamin_idx][sl], labels[sl], ids, env_ctx["envelope_labels"][sl],
                            audit_ctx["candidate_labels"], audit_ctx["missed_ids"], env_ctx["no_envelope_ids"],
                            title=sample["image_id"])

    def _render():
        if audit_ctx is not None and env_ctx is not None:
            # an optional stage's PNG must never abort the image: non-fatal, flagged on failure
            run_stage("audit_overlay", _audit_overlay, flags=flags, outcomes=outcomes)
        excl_ids = set(excluded["nucleus_id"]) if not excluded.empty else None
        make_montage(
            stack, labels, role_to_idx, out_dir / f"{sample['image_id']}__montage.png",
            foci_df=spots if len(spots) else None, excluded_ids=excl_ids,
            sc_mask=masks.get("syp_aggregate"),
            granule_mask=(granule_labels > 0) if granule_labels is not None else None,
            scalebar_um=float(cfg.get("render.scalebar_um", 10)),
            title=f"{sample['image_id']}  [{sample['sex']}/{sample['treatment']}]  n={n_nuclei}",
        )

    def _dump_stage_record(complete: bool) -> dict:
        # per-stage outcome record (ran / skipped + reason / failed + error, elapsed, flags); advisory,
        # so it can never fail a finished run, and written even when the montage raises.
        rec = {name: oc.as_dict() for name, oc in outcomes.items()}
        try:
            with open(long_path(out_dir / f"{sample['image_id']}__stages.json"), "w", encoding="utf-8") as fh:
                json.dump({"image_id": sample["image_id"], "config_hash": cfg.hash, "complete": complete,
                           **run_geometry, "stages": rec}, fh, indent=1)
        except Exception as e:  # noqa: BLE001 - the record is advisory
            log.warning("could not write the stage record: %s", e)
        return rec

    # write and render are fatal exactly as before the registry existed (an unwritten table or a broken
    # montage is a real failure of the image); the stage record still lands for post-mortems.
    complete = False
    try:
        run_stage("write", _write, flags=flags, outcomes=outcomes, fatal=True)
        run_stage("render", _render, flags=flags, outcomes=outcomes, fatal=True,
                  enabled=stage_enabled(cfg, "render"), skip_reason="disabled")
        complete = True
    finally:
        stage_record = _dump_stage_record(complete)
    # completion marker: the LAST file written, so `germquant batch --resume` can trust a folder that has
    # it (same config_hash, enabled stages and run geometry, no failed stage) and reprocess one that does not.
    timepoint = {k: run_geometry[k] for k in ("t_index", "t_seconds") if k in run_geometry}   # {} for a fixed image
    provenance.write_done_marker(
        out_dir, sample["image_id"], config_hash=cfg.hash, enabled_stages=run_geometry["enabled_stages"],
        stage_hashes=stage_hashes, xy_stride=run_geometry["xy_stride"], z_range=run_geometry["z_range"],
        failed_stages=[n for n, oc in outcomes.items() if oc.status == "failed"], extra=timepoint or None)

    log.info("%s: %d nuclei, %d germline, spots=%d, qc_pass=%s",
             sample["image_id"], n_nuclei, n_germline_nuclei, len(spots), qc_pass)
    return {"image_id": sample["image_id"], "n_nuclei": n_nuclei, "qc_pass": qc_pass,
            "qc_flags": qc_all, "out_dir": str(out_dir), "tables": tables, "stages": stage_record}


def _write_tables(tables, shared, out_dir, image_id, formats):
    for name, df in tables.items():
        df = _conform_schema(df, name)
        for k, v in shared.items():
            df[k] = v
        base = out_dir / f"{image_id}__{name}"
        if "csv" in formats:
            df.to_csv(long_path(base.with_suffix(".csv")), index=False)
        if "parquet" in formats:
            try:
                df.to_parquet(long_path(base.with_suffix(".parquet")), index=False)
            except Exception as e:  # pyarrow missing
                log.warning("parquet write failed (%s); CSV written.", e)


def _conform_schema(df, name):
    """Guarantee every declared schema column exists (filled NA if a stage didn't produce it)
    and is ordered first, so producer/consumer drift surfaces as an empty column rather than a
    KeyError in R. Extra columns a stage adds (e.g. per-role intensities) are kept after.
    """
    df = df.copy()
    declared = schema.TABLES.get(name, schema.SPOTS if name.startswith("spots_") else [])
    for col in declared:
        if col not in df.columns:
            df[col] = pd.NA
    extras = [c for c in df.columns if c not in declared]
    return df[list(declared) + extras]


def _save_labels(labels, path, spacing=None):
    """int32 label image. With `spacing` (dz, dy, dx um) it is written as OME-TIFF carrying the physical
    voxel size (the ImageJ format rejects int32), so Imaris and Fiji / Bio-Formats load it on the right
    grid; tifffile.imread returns the identical int32 array either way. A failed tagged write falls back
    to the plain write so a truncated file is never left behind."""
    try:
        import tifffile

        arr = labels.astype(np.int32)
        if spacing is not None:
            sp = tuple(float(s) for s in spacing)
            try:
                tifffile.imwrite(long_path(path), arr, compression="zlib", ome=True,
                                 resolution=(1 / sp[2], 1 / sp[1]),
                                 metadata={"axes": "ZYX", "PhysicalSizeX": sp[2], "PhysicalSizeXUnit": "µm",
                                           "PhysicalSizeY": sp[1], "PhysicalSizeYUnit": "µm",
                                           "PhysicalSizeZ": sp[0], "PhysicalSizeZUnit": "µm"})
                return
            except Exception as e:  # noqa: BLE001 - fall back to the untagged file, never a stub
                log.warning("tagged label write failed (%s); writing a plain TIFF", e)
        tifffile.imwrite(long_path(path), arr, compression="zlib")
    except Exception as e:  # noqa: BLE001
        log.warning("could not save label image: %s", e)


def _save_spots_image(spots, shape, spacing, path, radius_um=0.3):
    """3D blob image of detected spots, on the same voxel grid as the label TIF / original image, with
    voxel size baked in — load it in Imaris as a Channel, or run Imaris Spots detection on it."""
    try:
        import tifffile

        from .spots import spots_to_image

        img = spots_to_image(spots, shape, spacing, radius_um=radius_um)
        sp = tuple(float(s) for s in spacing)
        tifffile.imwrite(long_path(path), img, compression="zlib", imagej=True,
                         resolution=(1 / sp[2], 1 / sp[1]),
                         metadata={"spacing": sp[0], "unit": "um", "axes": "ZYX"})
    except Exception as e:  # noqa: BLE001
        log.warning("could not save spots image: %s", e)


def _save_mask_image(mask, spacing, path):
    """Binary 3D mask (SC ribbon / SYP aggregate) as a voxel-calibrated uint8 TIF for Imaris overlay."""
    try:
        import tifffile

        img = (np.asarray(mask) > 0).astype(np.uint8) * 255
        sp = tuple(float(s) for s in spacing)
        tifffile.imwrite(long_path(path), img, compression="zlib", imagej=True,
                         resolution=(1 / sp[2], 1 / sp[1]),
                         metadata={"spacing": sp[0], "unit": "um", "axes": "ZYX"})
    except Exception as e:  # noqa: BLE001
        log.warning("could not save mask image: %s", e)


def _granule_kwargs(cfg) -> dict:
    """Blob-segmentation parameters shared by the PGL-1 granules and the cytoplasmic SYP aggregates."""
    return {
        "thresholding_method": cfg.get("granule.thresholding_method", "threshold_triangle"),
        "gauss_sigma_um": float(cfg.get("granule.gauss_sigma_um", 0.1)),
        "min_volume_um3": float(cfg.get("granule.min_volume_um3", 0.03)),
        "max_volume_um3": float(cfg.get("granule.max_volume_um3", 8.0)),
    }


def _run_granule(stack, labels, role_to_idx, nuclei, spacing, cfg, env_ctx: dict | None = None) -> dict:
    """Granule stage: build the perinuclear region and cytoplasmic shell, surface the PGL-1 granules in
    the region and assign them to nuclei. Returns the context dict the coloc stage consumes
    (region, shell, granule labels/mask, assigned granules table, per-nucleus tallies).

    Two recipes (`granule.method`): the validated default `threshold_triangle` over the coloc region
    (unchanged), and `imaris_tophat` (docs/ROADMAP_modular_pipeline.md step 8), the August 2026
    Imaris-calibrated recipe, normally with `granule.region: cytoplasm_shell` = the cytoplasm within
    `envelope.cyto_um` of the lamin envelope (or of the DAPI nuclei when the envelope stage did not run),
    computed inside the padded germline crop exactly as the analysis scripts did."""
    from .granule import segment_granules, segment_granules_tophat

    germ_ids = {int(v) for v in nuclei["nucleus_id"].tolist()}
    region_name = str(cfg.get("coloc.region", "perinuclear_shell"))
    dilation_um = float(cfg.get("coloc.region_dilation_um", 1.5))
    region, nuc_union, shell = _build_region(labels, germ_ids, spacing, region_name, dilation_um)

    syp = stack.data[role_to_idx["central_element"]]
    pgl = stack.data[role_to_idx["granule"]]

    # Refine the cytoplasmic shell using the lamin (nuclear-envelope) channel when present: anchor it to
    # the real envelope instead of a fixed DAPI dilation. On real ccw77 data the shell voxel coloc
    # (coloc stage), which excludes the bright intranuclear SC ribbon, is what cleanly separates male
    # (high) from herm (low) — see docs/COLOCALIZATION.md.
    lamin_idx = role_to_idx.get("lamin")
    lamin_img = stack.data[lamin_idx] if lamin_idx is not None else None
    use_lamin = bool(cfg.get("coloc.use_lamin_shell", True)) and lamin_img is not None
    shell, shell_region_name = _perinuclear_shell(
        nuc_union, spacing, dilation_um, lamin_img, use_lamin, shell)
    g_kw = _granule_kwargs(cfg)

    method = str(cfg.get("granule.method", "threshold_triangle"))
    region_kind = str(cfg.get("granule.region", "coloc_region"))
    cyto_ctx = None
    if method == "imaris_tophat":
        from .envelope import cytoplasm_shell, germline_crop
        from .granule.segment import TOPHAT_DEFAULTS

        t_kw = {k: float(cfg.get(f"granule.tophat.{k}", v)) for k, v in TOPHAT_DEFAULTS.items()}
        cyto_um = float(cfg.get("envelope.cyto_um", 2.5))
        if region_kind == "cytoplasm_shell":
            # the August frame: the germline bounding box padded by 30 voxels; envelope labels from the
            # envelope stage, else the DAPI nuclei themselves
            sl = env_ctx["crop"] if env_ctx is not None else germline_crop(labels, germ_ids)
            env_mask = (env_ctx["envelope_labels"][sl] > 0) if env_ctx is not None else nuc_union[sl]
            region_kind = "cytoplasm_shell" if env_ctx is not None else "cytoplasm_shell_dapi"
            cyto_c, dt_c = cytoplasm_shell(env_mask, spacing, cyto_um)
            lab_c, _gran_c = segment_granules_tophat(np.asarray(pgl)[sl], cyto_c, spacing,
                                                     marker=cfg.channel_map.marker("granule"), **t_kw)
            granule_labels = np.zeros(labels.shape, np.int32)
            granule_labels[sl] = lab_c
            # re-measure on the full grid so centroids are whole-image microns
            from .granule.segment import _measure

            granules = _measure(granule_labels, np.asarray(pgl).astype(np.float32), np.asarray(spacing, float),
                                cfg.channel_map.marker("granule"))
            if not granules.empty:
                granules["detector"] = "imaris_tophat_cc"
            cyto_ctx = {"crop": sl, "cyto": cyto_c, "dt": dt_c, "env_mask": env_mask, "cyto_um": cyto_um}
        else:
            granule_labels, granules = segment_granules_tophat(
                pgl, region, spacing, marker=cfg.channel_map.marker("granule"), **t_kw)
    else:
        region_kind = "coloc_region"
        # PGL-1 granules across the whole perinuclear region (the validated default; unchanged)
        granule_labels, granules = segment_granules(
            pgl, region, spacing, marker=cfg.channel_map.marker("granule"), **g_kw)
    granules, per_nuc = _assign_granules_to_nuclei(granules, nuclei)
    return {
        "germ_ids": germ_ids, "region_name": region_name, "dilation_um": dilation_um,
        "region": region, "nuc_union": nuc_union, "shell": shell, "shell_region_name": shell_region_name,
        "syp": syp, "pgl": pgl, "granule_labels": granule_labels, "granule_mask": granule_labels > 0,
        "granules": granules, "per_nuc": per_nuc, "method": method, "region_kind": region_kind,
        "cyto_ctx": cyto_ctx,
    }


def _run_coloc(stack, labels, role_to_idx, spacing, cfg, ctx: dict):
    """Coloc stage: surface both SYP operands (intranuclear SC ribbon = control, cytoplasmic aggregate)
    and colocalize each vs the granules from the granule stage, then the shell voxel / partition
    coefficient headline row. Returns (coloc_df, masks). Per-granule overlap columns are merged into
    ``ctx["granules"]`` in place."""
    from .coloc import colocalize, partition_coefficient, shell_voxel_coloc
    from .granule import segment_granules
    from .sc import surface_sc_ribbon

    syp, pgl = ctx["syp"], ctx["pgl"]
    region, shell = ctx["region"], ctx["shell"]
    region_name, dilation_um, shell_region_name = ctx["region_name"], ctx["dilation_um"], ctx["shell_region_name"]
    granule_labels, granule_mask, granules = ctx["granule_labels"], ctx["granule_mask"], ctx["granules"]
    g_kw = _granule_kwargs(cfg)

    # SYP operands: intranuclear SC ribbon (ridge filter) + cytoplasmic aggregate (blob-seg in the shell)
    masks: dict = {}
    if cfg.get("sc.enabled", True):
        masks["sc_ribbon"] = surface_sc_ribbon(
            syp, labels, spacing, keep_nucleus_ids=ctx["germ_ids"],
            ridge_sigmas_um=tuple(cfg.get("sc.ridge_sigmas_um", [0.15, 0.25, 0.40])),
            ridge_hyst_low_pct=float(cfg.get("sc.ridge_hyst_low_pct", 45.0)),
            ridge_hyst_high_pct=float(cfg.get("sc.ridge_hyst_high_pct", 80.0)),
        )
    agg_labels, _ = segment_granules(syp, shell, spacing, marker="SYP-agg", **g_kw)
    masks["syp_aggregate"] = agg_labels > 0

    operand_cfg = str(cfg.get("coloc.sc_operand", "both"))
    operands = ["syp_aggregate", "sc_ribbon"] if operand_cfg == "both" else [operand_cfg]

    coloc_rows = []
    for op in operands:
        mask = masks.get(op)
        if mask is None:
            continue
        row, per_g = colocalize(
            mask, granule_mask, granule_labels, syp, pgl, region, spacing,
            sc_operand=op, region_name=region_name, region_dilation_um=dilation_um,
            n_random=int(cfg.get("coloc.n_random", 100)),
            object_overlap_min_frac=float(cfg.get("coloc.object_overlap_min_frac", 0.0)),
            costes=bool(cfg.get("coloc.costes", False)),
        )
        coloc_rows.append(row)
        if not per_g.empty and not granules.empty:
            per_g = per_g.rename(columns={
                "overlaps": f"overlaps_{op}", "overlap_frac": f"overlap_frac_{op}",
                "nearest_um": f"nearest_{op}_um"})
            granules = granules.merge(per_g, on="granule_id", how="left")

    # HEADLINE: threshold-light SYP<->PGL-1 voxel coloc in the (lamin-defined) cytoplasmic shell,
    # excluding the intranuclear SC ribbon. Robust where the object-overlap metric is not.
    sv = shell_voxel_coloc(syp, pgl, shell, spacing)
    # HEADLINE 2: partition coefficient — SYP-3 enrichment INSIDE the PGL-1 p-granules vs the shell
    # cytoplasm. Exposure-INDEPENDENT (a ratio) + blur-robust (rotation null), so it holds where the
    # shell voxel Pearson is confounded by the dim-mCherry SYP background / variable exposure.
    pc = partition_coefficient(syp, granule_mask, shell, spacing)
    coloc_rows.append({
        "sc_operand": "shell_voxel", "region": shell_region_name, "region_dilation_um": dilation_um,
        "region_voxels": sv["shell_voxels"], "region_volume_um3": sv["shell_volume_um3"],
        "pearson_r": sv["shell_pearson"], "manders_m1": sv["shell_manders_m1"],
        "manders_m2": sv["shell_manders_m2"],
        "partition_coef": pc["partition_coef"], "partition_coef_rot": pc["partition_coef_rot"],
        "pc_gran_voxels": pc["pc_gran_voxels"],
    })

    coloc_df = pd.DataFrame(coloc_rows, columns=list(schema.COLOC))
    ctx["granules"] = granules
    return coloc_df, masks


def _build_region(labels, germ_ids, spacing, region_name, dilation_um):
    """Coloc region R + the perinuclear shell (R minus nucleus interiors). `perinuclear_shell` (default)
    dilates the germline-nucleus union by `dilation_um` (spacing-aware EDT) so P-granules and cytoplasmic
    SYP aggregates just outside the envelope are inside R. Returns (R, nucleus_union, shell)."""
    from scipy import ndimage as ndi

    nuc_union = np.isin(labels, list(germ_ids)) if germ_ids else np.zeros(labels.shape, bool)
    if region_name == "nuclear" or dilation_um <= 0:
        region = nuc_union.copy()
    elif region_name == "germline_bbox":
        region = np.zeros(labels.shape, bool)
        objs = ndi.find_objects(nuc_union.astype(np.int32))
        if objs and objs[0] is not None:
            region[objs[0]] = True
    else:  # perinuclear_shell
        if nuc_union.any():
            dt = ndi.distance_transform_edt(~nuc_union, sampling=tuple(float(s) for s in spacing))
            region = nuc_union | (dt <= dilation_um)
        else:
            region = nuc_union.copy()
    shell = region & ~nuc_union
    return region, nuc_union, shell


def _perinuclear_shell(nuc_union, spacing, dilation_um, lamin_img, use_lamin, dapi_shell):
    """Cytoplasmic shell just OUTSIDE the germline nuclei. When a lamin channel is present and
    `use_lamin`, anchor the shell to the REAL nuclear envelope (thresholded lamin) rather than the
    fixed DAPI dilation — the anatomically correct perinuclear zone where P-granules dock. Falls back
    to the DAPI-dilation shell if lamin is absent/degenerate. Returns (shell_mask, region_label)."""
    if not use_lamin or lamin_img is None:
        return dapi_shell, "dapi_shell"
    from scipy import ndimage as ndi
    from skimage.filters import threshold_triangle

    sp = tuple(float(s) for s in spacing)
    dt_out = ndi.distance_transform_edt(~nuc_union, sampling=sp)
    gonad = nuc_union | (dt_out <= max(3.0, dilation_um * 2))       # germline neighbourhood
    vals = np.asarray(lamin_img)[gonad]
    if vals.size == 0 or float(vals.max()) <= float(vals.min()):
        return dapi_shell, "dapi_shell"
    env = (np.asarray(lamin_img) > threshold_triangle(vals)) & gonad
    if not env.any():
        return dapi_shell, "dapi_shell"
    dt_env = ndi.distance_transform_edt(~env, sampling=sp)
    shell = (dt_env <= dilation_um) & (~nuc_union) & gonad          # near the envelope, outside nucleus
    if int(shell.sum()) < 100:
        return dapi_shell, "dapi_shell"
    return shell, "lamin_shell"


def _assign_granules_to_nuclei(granules, nuclei):
    """Assign each granule to the nearest germline-nucleus centroid (perinuclear granules belong to
    their nucleus) and tally per-nucleus n_granules + summed granule volume. Returns (granules_df with
    nucleus_id, per_nucleus_df)."""
    per_cols = ["nucleus_id", "n_granules", "granule_volume_um3"]
    if granules.empty or nuclei.empty:
        gd = granules.copy()
        gd["nucleus_id"] = pd.NA
        return gd, pd.DataFrame(columns=per_cols)
    from scipy.spatial import cKDTree

    cent = nuclei[["centroid_z_um", "centroid_y_um", "centroid_x_um"]].to_numpy(float)
    ids = nuclei["nucleus_id"].to_numpy()
    _, idx = cKDTree(cent).query(granules[["z_um", "y_um", "x_um"]].to_numpy(float))
    gd = granules.copy()
    gd["nucleus_id"] = ids[idx]
    per_nuc = (gd.groupby("nucleus_id")
               .agg(n_granules=("granule_id", "size"), granule_volume_um3=("volume_um3", "sum"))
               .reset_index())
    return gd, per_nuc


def _coloc_summary(coloc_df, granules):
    """Headline coloc fields for image_summary (the cytoplasmic SYP-aggregate operand is the headline)."""
    def _get(op, col):
        r = coloc_df[coloc_df["sc_operand"] == op] if not coloc_df.empty else coloc_df
        if len(r) and pd.notna(r[col].iloc[0]):
            return float(r[col].iloc[0])
        return float("nan")

    return {
        "n_granules": int(len(granules)),
        # HEADLINE: shell voxel coloc (excludes the SC ribbon; separates male>herm on real data)
        "shell_pearson": _get("shell_voxel", "pearson_r"),
        "shell_manders_m1": _get("shell_voxel", "manders_m1"),
        "shell_manders_m2": _get("shell_voxel", "manders_m2"),
        # HEADLINE: SYP-3 partition coefficient in the p-granules (exposure-independent, blur-robust)
        "partition_coef": _get("shell_voxel", "partition_coef"),
        "partition_coef_rot": _get("shell_voxel", "partition_coef_rot"),
        # secondary object metrics (threshold-fragile — see docs)
        "manders_m1_syp_aggregate": _get("syp_aggregate", "manders_m1"),
        "manders_m2_syp_aggregate": _get("syp_aggregate", "manders_m2"),
        "frac_granules_overlapping_syp_aggregate": _get("syp_aggregate", "frac_granules_overlapping_sc"),
        "overlap_pvalue_syp_aggregate": _get("syp_aggregate", "overlap_pvalue"),
        "frac_granules_overlapping_sc_ribbon": _get("sc_ribbon", "frac_granules_overlapping_sc"),
    }
