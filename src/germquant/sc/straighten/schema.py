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
    "error", "punctum_sigma_z_um",
]

SAMPLES = [
    "image_id", "strand_id", "strand_uid", "nucleus_id", "i", "s_um", "in_strand", "z_um", "y_um", "x_um",
    "t_z", "t_y", "t_x", "nlat_z", "nlat_y", "nlat_x", "nax_z", "nax_y", "nax_x", "core", "core_se", "integrated",
    "foci_core", "dna_core", "contaminated", "in_gap", "in_label", "strength", "lineness", "multiplicity",
    "roll_from_axial_deg",
]

GAPS = [
    "gap_id", "strand_id", "strand_uid", "nucleus_id", "s_start_um", "s_end_um", "length_um", "depth", "z_score",
    "bridged", "mid_z_um", "mid_y_um", "mid_x_um", "plateau_rel",
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
    "is_pachytene", "touches_z", "touches_xy", "qc_status", "qc_flags", "error_message", "n_puncta",
    "puncta_signal_frac", "steep_length_frac",
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

# One line per declared column (docs/SC_STRAIGHTEN.md carries the tables docs_columns_table() makes; a test keeps
# them equal and every declared column described). Units: um, counts above the territory background (core),
# counts x um^2 (integrated).
_SHARED = {
    "image_id": "the image id of the run",
    "strand_id": "strand number in the image, assigned by (nucleus_id, length descending, end coordinates)",
    "strand_uid": "first 10 hex of the sha1 of the in-strand centreline rounded to 0.01 um (stable for the same "
                  "geometry)",
    "nucleus_id": "the nucleus label",
}
COLUMN_DOCS: dict[str, dict[str, str]] = {
    "strands": {
        **{k: _SHARED[k] for k in ("strand_id", "strand_uid", "nucleus_id")},
        "marker": "the central_element marker named by the run's channel map",
        "geometry_source": "auto (traced) or imported (--paths-from)",
        "verified": "imported geometry marked verified",
        "strand_class": "complete, partial, merged_suspect, fragment, punctum or verified (see Strand classes)",
        "qc_status": "pass, review or fail: the worst level among qc_flags",
        "qc_flags": "the strand's QC codes, separated by ;",
        "length_um": "3D length between the measured ends (anchored half-maximum end rule)",
        "length_xy_um": "the same length with z ignored (short for steep strands)",
        "length_traced_um": "length of the refined centreline before the end rule",
        "e0_um": "net change of length at end 0 against the traced end (positive: the end lies beyond it)",
        "e1_um": "the same at end 1",
        "length_se_um": "standard error of length_um (both end fits and a curve term)",
        "chord_um": "straight distance between the two ends",
        "straightness": "chord_um / length_um",
        "curvature_mean_per_um": "mean curvature of the centreline (1/um)",
        "curvature_max_per_um": "largest curvature of the centreline (1/um)",
        "tz_abs_mean": "mean absolute z component of the unit tangent",
        "frac_steep45": "share of the strand steeper than 45 degrees to the image plane",
        "frac_steep60": "share of the strand steeper than 60 degrees to the image plane",
        "z_extent_um": "z range of the strand",
        "fit_rms_um": "residual RMS of the centreline spline fit",
        "n_links": "end-matching links joined into the strand",
        "link_cost_max": "largest link cost (NaN without links)",
        "frac_in_label": "share of the strand inside its nucleus label",
        "frac_near_other_label": "share of the strand within 0.3 um of another nucleus label",
        "roll_from_axial_med_deg": "median angle between n_ax and the optical axis projected into the normal plane",
        "end0_class": "end 0: free_fit, free_unfit, contact, territory_edge, stack_edge, cycle_cut, max_length or "
                      "as_drawn",
        "end1_class": "the same for end 1",
        "end0_z_um": "end 0 position (whole-image um, z)",
        "end0_y_um": "end 0 position (y)",
        "end0_x_um": "end 0 position (x)",
        "end1_z_um": "end 1 position (z); s runs from end 0, the end with the smaller rounded (z, y, x)",
        "end1_y_um": "end 1 position (y)",
        "end1_x_um": "end 1 position (x)",
        "core_median": "median core signal over clean inner samples (mean in the ellipse 0.12 x 0.25 um)",
        "integrated_median": "median integrated signal over clean inner samples (ellipse 0.35 x 0.7 um)",
        "snr_core_median": "median of core over its standard error",
        "profile_cv_corr": "noise-corrected coefficient of variation of the integrated profile",
        "n_gaps": "gaps on the strand",
        "gap_um": "summed gap length",
        "gap_fraction": "gap_um / length_um",
        "n_pieces": "pieces of 0.25 um or more after splitting the strand at its gaps",
        "longest_piece_um": "longest piece",
        "multiplicity_median": "median integrated signal over the nucleus's median strand signal",
        "fwhm_lat_um": "lateral FWHM of the mean cross-section (QC)",
        "fwhm_ax_um": "axial FWHM of the mean cross-section (QC)",
        "width_lat_decon_um": "fwhm_lat_um with the PSF FWHM removed in quadrature (QC; the SC is below "
                              "resolution)",
        "width_ax_decon_um": "fwhm_ax_um with the PSF FWHM removed in quadrature (QC)",
        "frac_contaminated": "share of the samples near another strand's centreline (contamination_um)",
        "n_foci": "foci assigned to the strand (unique or ambiguous)",
        "strip_offset": "first column of the strand in the strip TIFFs",
        "strip_len": "number of strip columns of the strand (s from -0.5 um to length_um + 0.5 um every 0.05 um)",
        "e0_se_um": "standard error of the end 0 position",
        "e1_se_um": "standard error of the end 1 position",
        "end0_sigma_um": "fitted erf sigma at end 0",
        "end1_sigma_um": "fitted erf sigma at end 1",
        "n_members": "traced pieces joined into the strand",
        "kink_max_deg": "largest turn of the tangent within qc.kink_um",
        "frac_saturated": "share of core samples at the detector maximum",
        "near_stack_edge": "a sample lies within qc.stack_edge_um of the first or last plane",
        "merged_run_um": "longest run with multiplicity at or above classes.merged_multiplicity",
        "frac_multiplicity_high": "share of samples with multiplicity at or above classes.merged_multiplicity",
        "gap_max_um": "longest gap",
        "end_types_traced": "end types after tracing, end matching and accounting (end 0 and end 1, "
                            "joined by a vertical bar)",
        "import_id": "id of the imported path (imported geometry)",
        "error": "exception type when measuring the strand failed (empty otherwise)",
        "punctum_sigma_z_um": "axial extent of a steep strand at most punctum.max_xy_um long in the image plane: the "
                              "sigma of a Gaussian fitted along s to its core profile times tz_abs_mean (the punctum "
                              "rule; NaN where the rule does not apply)",
    },
    "samples": {
        **_SHARED,
        "i": "sample index along the strand (every 0.1 um)",
        "s_um": "arc position from end 0 (negative, or beyond length_um, in the 0.5 um margins)",
        "in_strand": "s_um lies within [0, length_um]",
        "z_um": "centreline position (whole-image um, z)",
        "y_um": "centreline position (y)",
        "x_um": "centreline position (x)",
        "t_z": "unit tangent (z)",
        "t_y": "unit tangent (y)",
        "t_x": "unit tangent (x)",
        "nlat_z": "lateral normal n_lat = t x n_ax (z); the u axis of the flat strip",
        "nlat_y": "lateral normal (y)",
        "nlat_x": "lateral normal (x)",
        "nax_z": "axial normal of the rotation-minimising frame rolled toward the optical axis (z); the v axis of "
                 "the axial strip",
        "nax_y": "axial normal (y)",
        "nax_x": "axial normal (x)",
        "core": "core signal (counts above the territory background)",
        "core_se": "standard error of core (noise model)",
        "integrated": "integrated signal (counts x um^2)",
        "foci_core": "core signal of the foci role (NA when not read)",
        "dna_core": "core signal of the dna role (NA when not read)",
        "contaminated": "within contamination_um of another strand's centreline",
        "in_gap": "inside a detected gap",
        "in_label": "inside the nucleus label",
        "strength": "ridge strength (noise-null units)",
        "lineness": "ratio of the largest to the middle Hessian eigenvalue magnitude (low: line-like, high: "
                    "a junction or a blob)",
        "multiplicity": "integrated over the nucleus's median strand signal",
        "roll_from_axial_deg": "angle between n_ax and the optical axis projected into the normal plane",
    },
    "gaps": {
        "gap_id": "gap number in the image",
        **{k: _SHARED[k] for k in ("strand_id", "strand_uid", "nucleus_id")},
        "s_start_um": "start of the gap on the strand",
        "s_end_um": "end of the gap on the strand",
        "length_um": "missing signal integrated over the gap window widened by the blur, as a fraction of the "
                     "local plateau",
        "depth": "(plateau - minimum of the smoothed profile) / plateau; above 1 when the local background "
                 "exceeds the signal there",
        "z_score": "(plateau - minimum) / its standard error",
        "bridged": "an end-matching link lies in the gap (the tracer was joined across it)",
        "mid_z_um": "gap midpoint (z)",
        "mid_y_um": "gap midpoint (y)",
        "mid_x_um": "gap midpoint (x)",
        "plateau_rel": "the local plateau of the gap over the strand's median integrated signal",
    },
    "foci": {
        "spot_table": "spots (RAD-51), or the name of a spots_<name> table (COSA-1 and others)",
        "spot_id": "the spot's id in its table",
        "nucleus_id": "the spot's nucleus",
        "marker": "the spot table's marker",
        "assignment": "unique, ambiguous (a second strand about as close), off_strand or no_strand",
        "strand_id": "the closest strand",
        "strand_uid": "its uid",
        "s_um": "position along that strand",
        "s_frac": "s_um / length_um, clipped to [0, 1]",
        "s_folded_frac": "min(s, L - s) / L (does not depend on the orientation)",
        "u_um": "lateral offset from the centreline (along n_lat)",
        "v_um": "axial offset from the centreline (along n_ax)",
        "dist_um": "distance to the closest strand",
        "second_strand_id": "the next closest strand",
        "in_gap": "the spot lies in a gap of its strand",
        "rank_on_strand": "order along the strand among the assigned foci of the table",
        "next_focus_um": "distance along the strand to the next assigned focus",
        "next_focus_frac": "next_focus_um / length_um",
        "z_um": "spot position (z)",
        "y_um": "spot position (y)",
        "x_um": "spot position (x)",
    },
    "nuclei": {
        "nucleus_id": "the nucleus label",
        "status": "traced, no_sc (no SC above the presence rules), skipped or failed",
        "error": "exception type of a failed nucleus",
        "n_strands": "strands after accounting (punctum strands not counted under punctum.exclude)",
        "n_main": "strands of 1 um or more (the male X rule may drop one)",
        "n_complete": "complete or verified strands",
        "n_fragments": "strands under 1 um",
        "n_merged_suspect": "merged_suspect strands",
        "expected": "expected SC count: 6 for oocytes, 5 for spermatocytes",
        "count_status": "low, ok or high: n_main against expected",
        "complete_nucleus": "every main strand complete, n_complete == expected and no nucleus fail flag",
        "total_length_um": "summed strand length (every traced nucleus)",
        "total_length_xy_um": "summed length_xy_um",
        "total_length_complete_um": "total_length_um of complete nuclei only (NaN otherwise)",
        "n_pieces": "pieces of 0.25 um or more after splitting the strands at their gaps",
        "mean_piece_um": "mean of those pieces",
        "longest_piece_um": "longest piece",
        "n_gaps": "gaps on the nucleus's strands",
        "gaps_per_10um": "10 x n_gaps / total_length_um (NaN under 1 um)",
        "gap_fraction": "summed gap length / total_length_um",
        "piece_lw_mean_um": "length-weighted mean piece length, sum(l^2) / sum(l): the piece that holds a random "
                            "point of the nucleus's SC (the continuity score of the validation protocol)",
        "on_strand_frac": "SYP above background inside the strand tubes over SYP above background in the territory",
        "snr_core_median": "median core SNR over the nucleus's strands",
        "coverage_proxy": "share of strong crest voxels within 0.3 um of a strand",
        "multiplicity_frac": "share of the strand length with multiplicity at or above classes.merged_multiplicity",
        "noise_sd": "SD of the transformed noise in the territory (about 1 when the noise model fits)",
        "n_seeds": "crest seeds of the tracer",
        "budget_exceeded": "tracing stopped at the work budget",
        "volume_um3": "label volume (from the run)",
        "in_germline": "the run's germline call",
        "axis_position_um": "position along the gonad axis (from the run)",
        "zone": "the run's zone (when staging ran)",
        "is_pachytene": "the run's pachytene call (when staging ran)",
        "touches_z": "the label touches the first or last plane",
        "touches_xy": "the label touches the image border",
        "qc_status": "pass, review or fail: the worst level among qc_flags",
        "qc_flags": "the nucleus's QC codes, separated by ;",
        "error_message": "the exception message of a failed nucleus",
        "n_puncta": "punctum strands (the axial blur of SYP puncta or aggregates); under punctum.exclude they are left "
                    "out of every other count, length, continuity and intensity column of this table",
        "puncta_signal_frac": "SYP above background inside the punctum tubes (outside the strand tubes) over SYP above "
                              "background in the territory",
        "steep_length_frac": "share of total_length_um steeper than 60 degrees to the image plane",
    },
}


def docs_columns_table(table: str) -> str:
    """The markdown column table docs/SC_STRAIGHTEN.md carries for one table (kept equal by a test)."""
    rows = ["| column | meaning |", "|---|---|"]
    for c in TABLES[table]:
        rows.append(f"| {c} | {COLUMN_DOCS[table][c]} |")
    return "\n".join(rows)
