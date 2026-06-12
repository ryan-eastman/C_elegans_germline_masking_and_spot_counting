"""Tidy output schema — long format, one row per object. All lengths µm, areas µm²,
volumes µm³, intensities raw a.u. Every table carries the shared metadata block so
R/Positron joins on image_id (+ nucleus_id) and facets by genotype/sex/treatment/zone.
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
    "axis_position_norm", "axis_position_um",
    "zone_call", "zone_call_method", "pachytene_subzone", "zone_confidence",
    "chromatin_polarized", "dna_mean_intensity",
    "sc_total_length_um", "sc_n_fragments", "n_foci",
]

# one row per SC track/fragment (the fragmentation readout lives here)
SC_TRACKS = [
    "track_id", "nucleus_id", "channel_role", "marker",
    "length_um", "n_branches", "n_junctions", "tortuosity", "trace_method",
]

# per-nucleus SC aggregate
SC_PER_NUCLEUS = [
    "nucleus_id", "marker", "n_fragments", "sc_total_length_um",
    "sc_mean_fragment_um", "sc_median_fragment_um", "sc_longest_fragment_um",
    "sc_mean_intensity", "expected_n_tracks",
]

FOCI = [
    "focus_id", "nucleus_id", "marker",
    "z_um", "y_um", "x_um", "sigma_um",
    "intensity_mean", "intensity_max", "detector", "detection_radius_um",
]

GRANULES = [
    "granule_id", "nucleus_id", "marker",
    "volume_um3", "surface_area_um2", "n_voxels", "sphericity",
    "intensity_mean", "intensity_integrated", "intensity_max",
    "centroid_z_um", "centroid_y_um", "centroid_x_um", "threshold_method",
]

# one row per zone per germline
ZONES = [
    "zone", "length_um", "length_rows", "fraction_of_germline", "n_nuclei",
    "distal_boundary_position_norm", "proximal_boundary_position_norm", "boundary_method",
]

# one row per image (QC summary)
IMAGE_SUMMARY = [
    "n_nuclei", "n_pachytene_nuclei",
    "mean_sc_total_length_um", "mean_sc_n_fragments", "mean_foci",
    "total_germline_length_um", "qc_pass", "qc_flags",
]

TABLES = {
    "nuclei": NUCLEI,
    "sc_tracks": SC_TRACKS,
    "sc_per_nucleus": SC_PER_NUCLEUS,
    "foci": FOCI,
    "granules": GRANULES,
    "zones": ZONES,
    "image_summary": IMAGE_SUMMARY,
}
