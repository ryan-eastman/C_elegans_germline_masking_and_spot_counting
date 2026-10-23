"""Declared columns of the sc_straighten tables. They live OUTSIDE germquant.schema.TABLES on purpose: a default
run neither writes nor knows them, so every existing output stays byte-identical. Append-only, like the main
schema. The strands, gaps, foci and nuclei tables carry the run's shared metadata block (germquant.schema
SHARED_META, copied from the run's nuclei.csv) plus ``sc_straighten_version`` and ``sc_straighten_hash``; the
samples table carries only ``image_id`` and its keys (join to the rest through the strands table).
"""
from __future__ import annotations

SUFFIXES = {
    "strands": "__sc_strands.csv",
    "samples": "__sc_strand_samples.csv",
    "gaps": "__sc_strand_gaps.csv",
    "foci": "__sc_strand_foci.csv",
    "nuclei": "__sc_strand_nuclei.csv",
    "strips_flat": "__sc_strand_strips_flat.tif",
    "strips_axial": "__sc_strand_strips_axial.tif",
    "overview": "__sc_strand_overview.png",
    "record": "__sc_straighten.json",
}

STRANDS = [
    "strand_id", "strand_uid", "nucleus_id", "marker", "geometry_source", "verified", "strand_class", "qc_status",
    "qc_flags",
    "length_um", "length_xy_um", "length_traced_um", "e0_um", "e1_um", "length_se_um", "chord_um", "straightness",
    "curvature_mean_per_um", "curvature_max_per_um", "tz_abs_mean", "frac_steep45", "frac_steep60", "z_extent_um",
    "fit_rms_um", "n_links", "link_cost_max", "frac_in_label", "frac_near_other_label", "roll_from_axial_med_deg",
    "end0_class", "end1_class", "end0_z_um", "end0_y_um", "end0_x_um", "end1_z_um", "end1_y_um", "end1_x_um",
    "core_median", "integrated_median", "snr_core_median", "profile_cv_corr", "n_gaps", "gap_um", "gap_fraction",
    "n_pieces", "longest_piece_um", "multiplicity_median", "fwhm_lat_um", "fwhm_ax_um", "width_lat_decon_um",
    "width_ax_decon_um", "frac_contaminated", "n_foci", "strip_offset", "strip_len",
    "e0_se_um", "e1_se_um", "end0_sigma_um", "end1_sigma_um", "n_members", "kink_max_deg", "frac_saturated",
    "near_stack_edge", "merged_run_um", "frac_multiplicity_high", "gap_max_um", "end_types_traced", "import_id",
]

SAMPLES = [
    "image_id", "strand_id", "strand_uid", "nucleus_id", "i", "s_um", "in_strand", "z_um", "y_um", "x_um",
    "t_z", "t_y", "t_x", "nlat_z", "nlat_y", "nlat_x", "nax_z", "nax_y", "nax_x", "core", "core_se", "integrated",
    "foci_core", "dna_core", "contaminated", "in_gap", "in_label", "strength", "lineness", "multiplicity",
    "roll_from_axial_deg",
]

GAPS = [
    "gap_id", "strand_id", "strand_uid", "nucleus_id", "s_start_um", "s_end_um", "length_um", "depth", "z_score",
    "bridged", "mid_z_um", "mid_y_um", "mid_x_um",
]

FOCI = [
    "spot_table", "spot_id", "nucleus_id", "marker", "assignment", "strand_id", "strand_uid", "s_um", "s_frac",
    "s_folded_frac", "u_um", "v_um", "dist_um", "second_strand_id", "in_gap", "rank_on_strand", "next_focus_um",
    "next_focus_frac", "z_um", "y_um", "x_um",
]

NUCLEI = [
    "nucleus_id", "status", "error", "n_strands", "n_main", "n_complete", "n_fragments", "n_merged_suspect",
    "expected", "count_status", "complete_nucleus", "total_length_um", "total_length_xy_um",
    "total_length_complete_um", "n_pieces", "mean_piece_um", "longest_piece_um", "n_gaps", "gaps_per_10um",
    "gap_fraction", "piece_lw_mean_um", "on_strand_frac", "snr_core_median", "coverage_proxy", "multiplicity_frac",
    "noise_sd", "n_seeds", "budget_exceeded", "volume_um3", "in_germline", "axis_position_um", "zone",
    "is_pachytene", "touches_z", "touches_xy", "qc_status", "qc_flags",
]

SUMMARY = [
    "image_id", "status", "sc_straighten_version", "sc_straighten_hash", "params_hash", "n_nuclei", "n_traced",
    "n_no_sc", "n_failed", "n_skipped", "n_strands", "n_complete", "n_partial", "n_merged_suspect", "n_fragment",
    "n_verified", "flags", "elapsed_s", "results_dir",
]

TABLES = {"strands": STRANDS, "samples": SAMPLES, "gaps": GAPS, "foci": FOCI, "nuclei": NUCLEI}
# the per-image tables `germquant collect` stacks (samples are not stacked: they are large and per-image)
COLLECTED = ("strands", "nuclei", "gaps", "foci")
BATCH_NAMES = {"strands": "batch_sc_strands.csv", "nuclei": "batch_sc_strand_nuclei.csv",
               "gaps": "batch_sc_strand_gaps.csv", "foci": "batch_sc_strand_foci.csv",
               "summary": "batch_sc_straighten_summary.csv"}
EXTRA_META = ["sc_straighten_version", "sc_straighten_hash"]
