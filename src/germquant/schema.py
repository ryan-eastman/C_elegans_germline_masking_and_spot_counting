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
]

# one row per detected spot (RAD-51 etc.) — SpotMAX detector.
# Carries the spot-vs-background effect size (>0 = above local background) for calibration.
SPOTS = [
    "spot_id", "nucleus_id", "marker",
    "z_um", "y_um", "x_um",
    "effect_size", "spot_mean_intensity", "detector",
]

# one row per image (QC summary)
IMAGE_SUMMARY = [
    "n_nuclei", "n_germline_nuclei", "mean_spots",
    "total_germline_length_um", "qc_pass", "qc_flags",
]

TABLES = {
    "nuclei": NUCLEI,
    "spots": SPOTS,
    "image_summary": IMAGE_SUMMARY,
}
