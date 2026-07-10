"""Single-image pipeline: .nd2 -> tidy tables + montage + label mask + manifest.

Stages: read -> segment nuclei (Cellpose) -> measure -> isolate germline -> linearize axis ->
count spots (SpotMAX) -> tidy tables. Each stage is wrapped so one failure degrades gracefully
(flags it) instead of killing the batch. Voxel spacing from the .nd2 is threaded into every 3D op.
"""
from __future__ import annotations

import json
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd

from . import provenance, qc, schema
from .axis import linearize_germline
from .config import Config
from .fsutil import long_path
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
    os.makedirs(long_path(out_dir), exist_ok=True)
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
                max_spot_candidates=int(cfg.get("spots.max_spot_candidates", 30000)),
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

    # ---- p-granule (PGL-1) surfacing + SYP<->PGL-1 colocalization ----
    # Surface the SC (SYP) signal and the PGL-1 granules as 3D masks and measure their DIRECT
    # voxel/object overlap inside the germline dilated by a perinuclear shell (P-granules sit just
    # OUTSIDE the nuclear envelope, so the region MUST include the perinuclear cytoplasm or the
    # overlap reads ~0 by construction). Two SYP operands are compared: `syp_aggregate` = cytoplasmic
    # SYP blobs in the shell (the headline for P-granule coincidence) and `sc_ribbon` = the
    # intranuclear SC ribbon (control). The stage degrades gracefully — a failure flags, not kills.
    granules = pd.DataFrame(columns=schema.GRANULES)
    coloc = pd.DataFrame(columns=schema.COLOC)
    masks: dict = {}
    granule_labels = None
    per_nuc_granules = None
    coloc_summary_fields: dict = {}
    gran_idx = role_to_idx.get("granule")
    if (cfg.get("coloc.enabled", True) and gran_idx is not None
            and role_to_idx.get("central_element") is not None and n_germline_nuclei > 0):
        try:
            granules, coloc, masks, granule_labels, per_nuc_granules = _run_coloc(
                stack, labels, role_to_idx, nuclei, spacing, cfg)
            if per_nuc_granules is not None and not nuclei.empty:
                if not per_nuc_granules.empty:
                    nuclei = nuclei.merge(per_nuc_granules, on="nucleus_id", how="left")
                if "n_granules" not in nuclei.columns:
                    nuclei["n_granules"] = 0
                    nuclei["granule_volume_um3"] = 0.0
                nuclei["n_granules"] = nuclei["n_granules"].fillna(0).astype(int)
                nuclei["granule_volume_um3"] = nuclei["granule_volume_um3"].fillna(0.0)
            coloc_summary_fields = _coloc_summary(coloc, granules)
            flags.append(f"coloc:granules_n={len(granules)}")
        except Exception as e:  # noqa: BLE001 - surfacing/coloc failure shouldn't kill the batch
            log.warning("coloc/granule stage failed (%s: %s); continuing without it.",
                        type(e).__name__, e)
            flags.append(f"coloc:FAILED_{type(e).__name__}")
            per_nuc_granules = None

    # ---- QC ----
    qc_pass, qc_all = qc.qc_flags(
        n_nuclei=n_nuclei, channel_flags=ch_flags,
        axis_flags=[f for f in flags if f.startswith("axis")],
        spots_found=len(spots) > 0,
        spots_enabled=bool(cfg.get("spots.enabled", True)) and spots_idx is not None,
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
        **coloc_summary_fields,
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
        nuclei = pd.concat([nuclei, excluded], ignore_index=True)
    tables = {"nuclei": nuclei, "spots": spots, "granules": granules,
              "coloc": coloc, "image_summary": image_summary}
    _write_tables(tables, shared, out_dir, sample["image_id"], cfg.get("output.formats", ["csv"]))

    if cfg.get("output.write_label_images", True):
        _save_labels(labels, out_dir / f"{sample['image_id']}__nuclei_labels.tif")
    if cfg.get("output.write_spots_image", True) and len(spots):
        _save_spots_image(spots, labels.shape, spacing, out_dir / f"{sample['image_id']}__spots.tif",
                          radius_um=float(cfg.get("spots.spot_radius_um", 0.3)))
    # surfaced objects for Imaris: PGL-1 granules as a label image (-> Surfaces), the SC ribbon and
    # cytoplasmic SYP-aggregate masks as calibrated binary TIFs (overlay them on the raw channels).
    if cfg.get("output.write_label_images", True) and granule_labels is not None and granule_labels.max() > 0:
        _save_labels(granule_labels, out_dir / f"{sample['image_id']}__granules_labels.tif")
    for mname in ("sc_ribbon", "syp_aggregate"):
        m = masks.get(mname)
        if m is not None and m.any():
            _save_mask_image(m, spacing, out_dir / f"{sample['image_id']}__{mname}.tif")
    if cfg.get("render.montage", True):
        excl_ids = set(excluded["nucleus_id"]) if not excluded.empty else None
        make_montage(
            stack, labels, role_to_idx, out_dir / f"{sample['image_id']}__montage.png",
            foci_df=spots if len(spots) else None, excluded_ids=excl_ids,
            sc_mask=masks.get("syp_aggregate"),
            granule_mask=(granule_labels > 0) if granule_labels is not None else None,
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
    declared = schema.TABLES.get(name, [])
    for col in declared:
        if col not in df.columns:
            df[col] = pd.NA
    extras = [c for c in df.columns if c not in declared]
    return df[list(declared) + extras]


def _save_labels(labels, path):
    try:
        import tifffile

        tifffile.imwrite(long_path(path), labels.astype(np.int32), compression="zlib")
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


def _run_coloc(stack, labels, role_to_idx, nuclei, spacing, cfg):
    """Surface PGL-1 granules + both SYP operands, then colocalize each vs the granules within the
    perinuclear region. Returns (granules_df, coloc_df, masks, granule_labels, per_nucleus_granules_df).
    """
    from .coloc import colocalize
    from .granule import segment_granules
    from .sc import surface_sc_ribbon

    germ_ids = {int(v) for v in nuclei["nucleus_id"].tolist()}
    region_name = str(cfg.get("coloc.region", "perinuclear_shell"))
    dilation_um = float(cfg.get("coloc.region_dilation_um", 1.5))
    region, _nuc_union, shell = _build_region(labels, germ_ids, spacing, region_name, dilation_um)

    syp = stack.data[role_to_idx["central_element"]]
    pgl = stack.data[role_to_idx["granule"]]
    g_kw = dict(
        thresholding_method=cfg.get("granule.thresholding_method", "threshold_triangle"),
        gauss_sigma_um=float(cfg.get("granule.gauss_sigma_um", 0.1)),
        min_volume_um3=float(cfg.get("granule.min_volume_um3", 0.03)),
        max_volume_um3=float(cfg.get("granule.max_volume_um3", 8.0)),
    )

    # PGL-1 granules across the whole perinuclear region
    granule_labels, granules = segment_granules(
        pgl, region, spacing, marker=cfg.channel_map.marker("granule"), **g_kw)
    granule_mask = granule_labels > 0

    # SYP operands: intranuclear SC ribbon (ridge filter) + cytoplasmic aggregate (blob-seg in the shell)
    masks: dict = {}
    if cfg.get("sc.enabled", True):
        masks["sc_ribbon"] = surface_sc_ribbon(
            syp, labels, spacing, keep_nucleus_ids=germ_ids,
            ridge_sigmas_um=tuple(cfg.get("sc.ridge_sigmas_um", [0.15, 0.25, 0.40])),
            ridge_hyst_low_pct=float(cfg.get("sc.ridge_hyst_low_pct", 45.0)),
            ridge_hyst_high_pct=float(cfg.get("sc.ridge_hyst_high_pct", 80.0)),
        )
    agg_labels, _ = segment_granules(syp, shell, spacing, marker="SYP-agg", **g_kw)
    masks["syp_aggregate"] = agg_labels > 0

    operand_cfg = str(cfg.get("coloc.sc_operand", "both"))
    operands = ["syp_aggregate", "sc_ribbon"] if operand_cfg == "both" else [operand_cfg]

    granules, per_nuc = _assign_granules_to_nuclei(granules, nuclei)
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

    coloc_df = pd.DataFrame(coloc_rows, columns=list(schema.COLOC))
    return granules, coloc_df, masks, granule_labels, per_nuc


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
        "manders_m1_syp_aggregate": _get("syp_aggregate", "manders_m1"),
        "manders_m2_syp_aggregate": _get("syp_aggregate", "manders_m2"),
        "frac_granules_overlapping_syp_aggregate": _get("syp_aggregate", "frac_granules_overlapping_sc"),
        "overlap_pvalue_syp_aggregate": _get("syp_aggregate", "overlap_pvalue"),
        "frac_granules_overlapping_sc_ribbon": _get("sc_ribbon", "frac_granules_overlapping_sc"),
    }
