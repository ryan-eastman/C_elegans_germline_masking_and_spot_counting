"""Tidy output schema — long format, one row per object. All lengths µm, volumes µm³,
intensities raw a.u. Every table carries the shared metadata block so R/Positron joins on
image_id (+ nucleus_id) and facets by genotype/sex/treatment.
"""

# Stamped onto every row of every table (provenance + design).
SHARED_META = [
    "image_id", "file_path", "genotype", "sex", "germ_cell", "treatment", "replicate",
    "acquisition_date", "voxel_dz_um", "voxel_dy_um", "voxel_dx_um",
    "n_channels", "dapi_present", "channel_map_json",
    "pipeline_version", "git_sha", "config_hash", "run_timestamp",
]

NUCLEI = [
    "nucleus_id", "volume_um3", "n_voxels",
    "centroid_z_um", "centroid_y_um", "centroid_x_um",
    "in_germline", "axis_position_norm", "axis_position_um",
    "n_spots",
    # p-granule (PGL-1) load assigned to the nearest germline nucleus (perinuclear granules)
    "n_granules", "granule_volume_um3",
]

# one row per detected spot (RAD-51 etc.) — SpotMAX detector.
# Carries the spot-vs-background effect size (>0 = above local background) for calibration.
SPOTS = [
    "spot_id", "nucleus_id", "marker",
    "z_um", "y_um", "x_um",
    "effect_size", "spot_mean_intensity", "detector",
]

# one row per surfaced p-granule (PGL-1) object. Overlap columns are filled per SC operand by the
# coloc stage: `*_syp_aggregate` = vs the cytoplasmic SYP-aggregate mask (headline for the
# perinuclear P-granule biology); `*_sc_ribbon` = vs the intranuclear SC ribbon mask (control).
GRANULES = [
    "granule_id", "nucleus_id", "marker",
    "z_um", "y_um", "x_um",
    "volume_um3", "n_voxels",
    "granule_mean_intensity", "granule_max_intensity",
    "overlaps_syp_aggregate", "overlap_frac_syp_aggregate", "nearest_syp_aggregate_um",
    "overlaps_sc_ribbon", "overlap_frac_sc_ribbon", "nearest_sc_ribbon_um",
    "detector",
]

# one row per image PER SC operand — SYP<->PGL-1 direct-overlap colocalization.
# Object/mask overlap (dice/jaccard/frac_granules_overlapping) is the headline; Manders M1/M2 are
# mask-restricted intensity fractions; pearson_r is a region-restricted DIAGNOSTIC only (inflated by
# the shared sparse background); overlap_pvalue/zscore are a translation-null test that the observed
# overlap exceeds chance given each mask's size within the region.
COLOC = [
    "sc_operand",              # "syp_aggregate" | "sc_ribbon"
    "region", "region_dilation_um", "region_voxels", "region_volume_um3",
    "n_granules", "n_granules_overlapping_sc", "frac_granules_overlapping_sc",
    "overlap_volume_um3", "dice", "jaccard",
    "frac_sc_in_granules", "frac_granule_in_sc",
    "manders_m1", "manders_m2", "pearson_r",
    "mean_granule_to_sc_um", "median_granule_to_sc_um",
    "sc_mask_voxels", "granule_mask_voxels",
    "n_random", "overlap_pvalue", "overlap_zscore",
    "costes_threshold_syp", "costes_threshold_pgl",
    "partition_coef", "partition_coef_rot", "pc_gran_voxels",   # condensate PC (exposure-independent)
]

# one row per image (QC summary). The coloc HEADLINE is the shell voxel coloc (SYP<->PGL-1 in the
# lamin-defined perinuclear cytoplasmic shell, SC ribbon excluded) — it separates male>herm on real
# data. The syp_aggregate object metrics are secondary (threshold-fragile).
IMAGE_SUMMARY = [
    "n_nuclei", "n_germline_nuclei", "mean_spots",
    "total_germline_length_um", "qc_pass", "qc_flags",
    "n_granules",
    "shell_pearson", "shell_manders_m1", "shell_manders_m2",
    "partition_coef", "partition_coef_rot",   # HEADLINE: SYP-3 enrichment in p-granules (exposure-independent)
    "manders_m1_syp_aggregate", "manders_m2_syp_aggregate",
    "frac_granules_overlapping_syp_aggregate", "overlap_pvalue_syp_aggregate",
    "frac_granules_overlapping_sc_ribbon",
]

TABLES = {
    "nuclei": NUCLEI,
    "spots": SPOTS,
    "granules": GRANULES,
    "coloc": COLOC,
    "image_summary": IMAGE_SUMMARY,
}
