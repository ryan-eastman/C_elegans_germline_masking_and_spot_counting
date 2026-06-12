"""Single-image pipeline: .nd2 -> tidy tables + montage + label mask + manifest.

Each stage is wrapped so one failure degrades gracefully (flags it) instead of killing
the whole batch. Voxel spacing from the .nd2 is threaded into every 3D op.
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
from .foci import detect_foci
from .io import parse_sample, read_nd2_metadata, read_stack
from .io.sample_metadata import expected_sc_count
from .measure import measure_objects
from .render import make_montage
from .sc import trace_sc
from .segment import segment_nuclei
from .zones import call_zones

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

    # ---- axis + zones ----
    axis_conf = float("nan")
    zones_tbl = pd.DataFrame(columns=schema.ZONES)
    if not nuclei.empty:
        nuclei, axis_conf, axis_flags = linearize_germline(
            nuclei, confidence_min=float(cfg.get("axis.qc_confidence_min", 0.6))
        )
        flags += axis_flags
        if cfg.get("zones.enabled", True):
            nuclei, zones_tbl, zone_flags = call_zones(
                nuclei, dna=dna, labels=labels, spacing=spacing,
                method=cfg.get("zones.method", "auto"),
                crescent_boundary=cfg.get("zones.crescent_boundary", "count_2plus"),
                axis_bin_um=float(cfg.get("zones.axis_bin_um", 5.0)),
                pachytene_thirds=bool(cfg.get("zones.pachytene_thirds", True)),
            )
            flags += zone_flags

    # ---- SC tracing + fragmentation ----
    sc_tracks = pd.DataFrame(columns=schema.SC_TRACKS)
    sc_per_nuc = pd.DataFrame(columns=schema.SC_PER_NUCLEUS)
    sc_traced = False
    ce_idx = role_to_idx.get("central_element")
    if cfg.get("sc.enabled", True) and ce_idx is not None and n_nuclei > 0:
        exp = expected_sc_count(sample["germ_cell"])
        exp_map = {int(nid): exp for nid in nuclei["nucleus_id"]} if exp else None
        sc_tracks, sc_per_nuc = trace_sc(
            stack.data[ce_idx], labels, spacing,
            marker=cfg.channel_map.marker("central_element"),
            ridge_sigmas_um=cfg.get("sc.ridge_sigmas_um", [0.15, 0.25, 0.40]),
            intensity_percentile=float(cfg.get("sc.intensity_percentile", 99.0)),
            min_fragment_length_um=float(cfg.get("sc.min_fragment_length_um", 0.5)),
            expected_n_tracks=exp_map,
        )
        sc_traced = not sc_per_nuc.empty
        if sc_traced:
            nuclei = nuclei.merge(
                sc_per_nuc[["nucleus_id", "n_fragments", "sc_total_length_um"]]
                .rename(columns={"n_fragments": "sc_n_fragments"}),
                on="nucleus_id", how="left",
            )

    # ---- RAD-51 foci ----
    foci = pd.DataFrame(columns=schema.FOCI)
    foci_idx = role_to_idx.get("foci")
    if cfg.get("foci.enabled", True) and foci_idx is not None:
        foci = detect_foci(
            stack.data[foci_idx], spacing, marker=cfg.channel_map.marker("foci"),
            min_sigma_um=float(cfg.get("foci.min_sigma_um", 0.10)),
            max_sigma_um=float(cfg.get("foci.max_sigma_um", 0.35)),
            threshold_rel=float(cfg.get("foci.threshold_rel", 0.10)), labels=labels,
        )
        if not nuclei.empty:
            counts = foci[foci["nucleus_id"] > 0].groupby("nucleus_id").size()
            nuclei["n_foci"] = nuclei["nucleus_id"].map(counts).fillna(0).astype(int)

    # ---- granules (optional generic 3D-object module; off the N2 critical path) ----
    granules = pd.DataFrame(columns=schema.GRANULES)
    gran_idx = role_to_idx.get("granule")
    if cfg.get("granules.enabled", False) and gran_idx is not None:
        from .granules import detect_granules

        granules = detect_granules(
            stack.data[gran_idx], spacing, labels=labels,
            marker=cfg.channel_map.marker("granule"),
            threshold_method=cfg.get("granules.threshold_method", "li"),
            min_volume_um3=float(cfg.get("granules.min_volume_um3", 0.05)),
        )

    # ---- QC ----
    mean_sc = float(sc_per_nuc["sc_total_length_um"].mean()) if sc_traced else None
    qc_pass, qc_all = qc.qc_flags(
        n_nuclei=n_nuclei, channel_flags=ch_flags, axis_flags=[f for f in flags if f.startswith("axis")],
        zone_flags=[f for f in flags if f.startswith("zones")], sc_traced=sc_traced,
        mean_sc_len=mean_sc, foci_found=len(foci) > 0,
    )
    qc_all = sorted(set(flags + qc_all))

    image_summary = pd.DataFrame([{
        "n_nuclei": n_nuclei,
        "n_pachytene_nuclei": int((nuclei.get("zone_call") == "pachytene").sum()) if "zone_call" in nuclei else 0,
        "mean_sc_total_length_um": mean_sc if mean_sc is not None else float("nan"),
        "mean_sc_n_fragments": float(sc_per_nuc["n_fragments"].mean()) if sc_traced else float("nan"),
        "mean_foci": float(nuclei["n_foci"].mean()) if "n_foci" in nuclei else float("nan"),
        "total_germline_length_um": float(nuclei["axis_position_um"].max()) if "axis_position_um" in nuclei and not nuclei.empty else float("nan"),
        "qc_pass": qc_pass, "qc_flags": ";".join(qc_all),
        "segmentation_method": seg_method, "axis_confidence": axis_conf,
    }])

    # ---- write outputs ----
    tables = {
        "nuclei": nuclei, "sc_tracks": sc_tracks, "sc_per_nucleus": sc_per_nuc,
        "foci": foci, "zones": zones_tbl, "granules": granules, "image_summary": image_summary,
    }
    _write_tables(tables, shared, out_dir, sample["image_id"], cfg.get("output.formats", ["csv"]))

    if cfg.get("output.write_label_images", True):
        _save_labels(labels, out_dir / f"{sample['image_id']}__nuclei_labels.tif")
    if cfg.get("render.montage", True):
        make_montage(
            stack, labels, role_to_idx, out_dir / f"{sample['image_id']}__montage.png",
            foci_df=foci if len(foci) else None, scalebar_um=float(cfg.get("render.scalebar_um", 10)),
            title=f"{sample['image_id']}  [{sample['sex']}/{sample['treatment']}]  n={n_nuclei}",
        )

    log.info("%s: %d nuclei, sc=%s, foci=%d, qc_pass=%s", sample["image_id"], n_nuclei, sc_traced, len(foci), qc_pass)
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
