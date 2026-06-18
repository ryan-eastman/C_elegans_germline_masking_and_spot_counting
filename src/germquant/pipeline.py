"""Single-image pipeline: .nd2 -> tidy tables + montage + label mask + manifest.

Stages: read -> segment nuclei (Cellpose) -> measure -> isolate germline -> linearize axis ->
count spots (SpotMAX) -> tidy tables. Each stage is wrapped so one failure degrades gracefully
(flags it) instead of killing the batch. Voxel spacing from the .nd2 is threaded into every 3D op.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from . import provenance, qc, schema
from .axis import linearize_germline
from .config import Config
from .io import parse_sample, read_nd2_metadata, read_stack
from .measure import measure_objects
from .render import make_montage
from .segment import segment_nuclei

log = logging.getLogger(__name__)


def process_image(
    nd2_path: str | Path,
    cfg: Config,
    out_dir: str | Path,
    *,
    xy_stride: int = 1,
    z_range: tuple[int, int] | None = None,
    prov: dict | None = None,
) -> dict:
    nd2_path = Path(nd2_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    flags: list[str] = []

    # ---- read + resolve channels ----
    meta = read_nd2_metadata(nd2_path)
    role_to_idx, ch_flags = cfg.channel_map.resolve(meta["channel_names"])
    flags += ch_flags
    stack = read_stack(nd2_path, xy_stride=xy_stride, z_range=z_range)
    spacing = stack.spacing
    if not stack.spacing_ok:
        # voxel size unreadable -> spacing is a 1 µm isotropic guess; every physical
        # length/area/volume is unreliable. Surface it loudly rather than silently.
        flags.append("voxel:UNREADABLE_assumed_isotropic_1um")
    dapi_present = role_to_idx.get("dna") is not None

    sample = parse_sample(nd2_path, cfg.get("metadata.filename_regex"), cfg.get("metadata.defaults"))
    prov = prov or provenance.write_manifest(out_dir, config_hash=cfg.hash, config=cfg.as_dict())

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

    # ---- segment nuclei ----
    dna = stack.channel(role_to_idx.get("dna"))
    if dna is None:
        flags.append("segment:no_dna_channel_using_first_channel")
        dna = stack.data[0]
    seg = cfg.segmentation.nuclei
    labels, seg_method = segment_nuclei(
        dna, spacing,
        method=seg.get("method", "auto"), cellpose_model=seg.get("cellpose_model", "cpsam"),
        diameter_um=float(seg.get("diameter_um", 3.0)), min_volume_um3=float(seg.get("min_volume_um3", 4.0)),
    )
    n_nuclei = int(labels.max())

    # ---- measure nuclei ----
    intensity = {r: stack.data[i] for r, i in role_to_idx.items() if i is not None}
    nuclei = measure_objects(labels, intensity, spacing, compute_surface=False)
    if not nuclei.empty:
        nuclei = nuclei.rename(columns={"label": "nucleus_id"})
    else:
        nuclei = pd.DataFrame(columns=["nucleus_id"])

    # ---- isolate germline (drop nuclei segmented OUTSIDE the gonad: gut, debris, off-gonad) ----
    # Uses SYP (central_element) intensity + spatial connectivity, NEVER the spot count. Excluded
    # nuclei are re-attached (flagged in_germline=False) before writing, so nothing is hidden;
    # downstream means (axis, spots/nucleus) operate on the germline subset.
    excluded = nuclei.iloc[0:0].copy()
    if cfg.get("germline.enabled", True) and not nuclei.empty:
        from .germline import select_germline

        nuclei, germ_flags = select_germline(
            nuclei,
            method=cfg.get("germline.method", "syp_seeded_cc"),
            syp_percentile=float(cfg.get("germline.syp_percentile", 25.0)),
            link_radius_um=float(cfg.get("germline.link_radius_um", 12.0)),
            min_seed_frac=float(cfg.get("germline.min_seed_frac", 0.10)),
            size_frac=float(cfg.get("germline.size_frac", 0.10)),
        )
        flags += germ_flags
        excluded = nuclei[~nuclei["in_germline"]].copy()
        nuclei = nuclei[nuclei["in_germline"]].copy()
    n_germline_nuclei = int(len(nuclei))

    # ---- linearize axis (principal-curve centerline -> per-nucleus distal->proximal position) ----
    axis_conf = float("nan")
    if not nuclei.empty:
        nuclei, axis_conf, axis_flags = linearize_germline(
            nuclei, confidence_min=float(cfg.get("axis.qc_confidence_min", 0.6))
        )
        flags += axis_flags

    # ---- spots (SpotMAX) — RAD-51 (or other) foci per nucleus.
    # Detects peaks ABOVE local background inside each nucleus mask, merges z-axis spot-splits, and
    # tags every spot with its effect size. Detection params are cross-validated vs Imaris (config). ----
    spots = pd.DataFrame(columns=schema.SPOTS)
    per_nuc_spots = None  # kept so off-gonad (excluded) nuclei also get n_spots at re-attach
    spots_idx = role_to_idx.get("foci")
    if cfg.get("spots.enabled", True) and spots_idx is not None and n_nuclei > 0:
        try:
            from .spots import detect_spots

            per_spot, per_nuc_spots = detect_spots(
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
            )
            spots = per_spot
            if not nuclei.empty and not per_nuc_spots.empty:
                nuclei = nuclei.merge(per_nuc_spots[["nucleus_id", "n_spots"]], on="nucleus_id", how="left")
                nuclei["n_spots"] = nuclei["n_spots"].fillna(0).astype(int)
            flags.append(f"spots:spotmax_n={len(spots)}")
        except Exception as e:  # noqa: BLE001 - SpotMAX missing or detection failure shouldn't kill the batch
            log.warning("spot detection failed (%s: %s); continuing without spots.", type(e).__name__, e)
            flags.append(f"spots:FAILED_{type(e).__name__}")
            per_nuc_spots = None

    # ---- QC ----
    qc_pass, qc_all = qc.qc_flags(
        n_nuclei=n_nuclei, channel_flags=ch_flags,
        axis_flags=[f for f in flags if f.startswith("axis")],
        zone_flags=[], sc_traced=False, mean_sc_len=None, foci_found=len(spots) > 0,
    )
    qc_all = sorted(set(flags + qc_all))

    image_summary = pd.DataFrame([{
        "n_nuclei": n_nuclei,
        "n_germline_nuclei": n_germline_nuclei,
        "mean_spots": float(nuclei["n_spots"].mean()) if "n_spots" in nuclei else float("nan"),
        "total_germline_length_um": float(nuclei["axis_position_um"].max())
        if "axis_position_um" in nuclei and not nuclei.empty else float("nan"),
        "qc_pass": qc_pass, "qc_flags": ";".join(qc_all),
        "segmentation_method": seg_method, "axis_confidence": axis_conf,
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
        nuclei = pd.concat([nuclei, excluded], ignore_index=True)
    tables = {"nuclei": nuclei, "spots": spots, "image_summary": image_summary}
    _write_tables(tables, shared, out_dir, sample["image_id"], cfg.get("output.formats", ["csv"]))

    if cfg.get("output.write_label_images", True):
        _save_labels(labels, out_dir / f"{sample['image_id']}__nuclei_labels.tif")
    if cfg.get("render.montage", True):
        excl_ids = set(excluded["nucleus_id"]) if not excluded.empty else None
        make_montage(
            stack, labels, role_to_idx, out_dir / f"{sample['image_id']}__montage.png",
            foci_df=spots if len(spots) else None, excluded_ids=excl_ids,
            scalebar_um=float(cfg.get("render.scalebar_um", 10)),
            title=f"{sample['image_id']}  [{sample['sex']}/{sample['treatment']}]  n={n_nuclei}",
        )

    log.info("%s: %d nuclei, %d germline, spots=%d, qc_pass=%s",
             sample["image_id"], n_nuclei, n_germline_nuclei, len(spots), qc_pass)
    return {"image_id": sample["image_id"], "n_nuclei": n_nuclei, "qc_pass": qc_pass,
            "qc_flags": qc_all, "out_dir": str(out_dir), "tables": tables}


def _write_tables(tables, shared, out_dir, image_id, formats):
    for name, df in tables.items():
        df = _conform_schema(df, name)
        for k, v in shared.items():
            df[k] = v
        base = out_dir / f"{image_id}__{name}"
        if "csv" in formats:
            df.to_csv(base.with_suffix(".csv"), index=False)
        if "parquet" in formats:
            try:
                df.to_parquet(base.with_suffix(".parquet"), index=False)
            except Exception as e:  # pyarrow missing
                log.warning("parquet write failed (%s); CSV written.", e)


def _conform_schema(df, name):
    """Guarantee every declared schema column exists (filled NA if a stage didn't produce it)
    and is ordered first, so producer/consumer drift surfaces as an empty column rather than a
    KeyError in R. Extra columns a stage adds (e.g. per-role intensities) are kept after.
    """
    df = df.copy()
    declared = schema.TABLES.get(name, [])
    for col in declared:
        if col not in df.columns:
            df[col] = pd.NA
    extras = [c for c in df.columns if c not in declared]
    return df[list(declared) + extras]


def _save_labels(labels, path):
    try:
        import tifffile

        tifffile.imwrite(str(path), labels.astype(np.int32), compression="zlib")
    except Exception as e:  # noqa: BLE001
        log.warning("could not save label image: %s", e)
