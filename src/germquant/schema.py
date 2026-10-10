"""Tidy output schema — long format, one row per object. All lengths µm, volumes µm³,
intensities raw a.u. Every table carries the shared metadata block so R/Positron joins on
image_id (+ nucleus_id) and facets by genotype/sex/treatment.

Growth rule (docs/ROADMAP_modular_pipeline.md): the schema is APPEND-ONLY. A new stage appends its
columns at the end of an existing table's declared list or adds a new table in STAGE_TABLES; existing
columns are never renamed, reordered, re-typed or removed. Declared columns are written first and any
per-role extras after them (pipeline._conform_schema), so the columns of a table written by an older
commit always appear in the same relative order in the current one and the golden regression gate
(tests/test_golden.py, scripts/regression_diff.py) can tell "columns added" from "values changed". Every declared table is written on every run, empty when
its stage is off, so consumers never hit a missing file.
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
    # in_germline: the germline stage's verdict (True = in the gonad; False = gut / debris / off-gonad,
    # kept in the table so nothing is hidden). When that stage is switched off (config/profiles/
    # generic.yaml, --no-germline) every segmented nucleus is the analysis set and the column is True on
    # every row, never NA, so downstream filters on in_germline keep working. axis_position_norm /
    # axis_position_um are NA on off-gonad rows and on every row when the axis stage is off. Whether a
    # stage ran is recorded in <image_id>__stages.json, not inferred from missing values.
    "in_germline", "axis_position_norm", "axis_position_um",
    "n_spots",
    # p-granule (PGL-1) load assigned to the nearest germline nucleus (perinuclear granules)
    "n_granules", "granule_volume_um3",
    # SC tracer (sc_trace stage, off by default): total SC length (recoverable), fragment count as a
    # LOWER BOUND (strands overlap in 3D), SYP-intensity-CV fragmentation index, expected SC number
    "sc_total_length_um", "sc_n_fragments_lb", "sc_fragmentation_index", "sc_expected_n_tracks",
    # lamin envelope (envelope stage, off by default): watershed envelope volume and its ratio to the DAPI
    # volume, whether the volume gate fell back to the DAPI label, the ring test (no envelope = sperm /
    # somatic / debris) and the 2D territory the nucleus belongs to
    "envelope_volume_um3", "envelope_vol_ratio", "envelope_fallback", "ring_ratio", "shell_over_thr",
    "has_envelope", "territory_id",
    "envelope_centroid_z_um", "envelope_centroid_y_um", "envelope_centroid_x_um",
    # staging stage (hand-traced pachytene axis): zone early/mid/late/pre/post/off_axis, arc length and
    # perpendicular distance to the trace, pachytene membership, territory crossed by the trace
    "zone", "s_um", "r_um", "is_pachytene", "in_territory", "off_axis_cut_um",
]

# staging stage: one row per germline nucleus (the same columns as appended to nuclei; `germquant
# restage` rewrites only this table)
ZONES = ["nucleus_id", "zone", "s_um", "r_um", "is_pachytene", "in_territory", "off_axis_cut_um"]

# partition stage: one row per region (whole shell; early / mid / late / pach when staged)
PARTITION = [
    "region", "n_granules", "granule_voxels", "outside_voxels", "region_voxels", "bg",
    "partition_coef", "partition_coef_rot", "partition_coef_zshift", "partition_coef_specific",
]

# granule_tail stage: one row per re-segmented granule (v4 envelope set, no-envelope nuclei dropped)
GRANULE_TAIL = ["granule_id", "excess_n", "syp_excess", "pgl_raw", "pgl_n", "vol_um3", "dist_um", "dmin_um"]

# mask audit (audit stage): one row per lamin-only nucleus candidate and per germline label
MASK_AUDIT = [
    "kind", "object_id", "vol_um3", "z_um", "y_um", "x_um", "cov_germ_label", "cov_any_label",
    "cov_envelope", "ring_shell", "ring_inside", "ring_ratio", "shell_over_thr", "no_envelope",
]

# one row per SC track/fragment (sc_trace stage)
SC_TRACKS = [
    "track_id", "nucleus_id", "channel_role", "marker",
    "length_um", "n_branches", "n_junctions", "tortuosity", "trace_method",
]

# per-nucleus SC aggregate (sc_trace stage)
SC_PER_NUCLEUS = [
    "nucleus_id", "marker", "n_fragments", "sc_total_length_um",
    "sc_mean_fragment_um", "sc_median_fragment_um", "sc_longest_fragment_um",
    "sc_mean_intensity", "sc_fragmentation_index", "expected_n_tracks",
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

# signal_qc stage: one block of per-channel signal metrics per channel role, appended to image_summary as
# sig_<role>_<metric> (NaN for a role the image does not carry; a role a custom channel map adds beyond
# these five rides along as an extra column after the declared ones). offset = camera offset (0.5th
# percentile), bg / noise = background median and robust sigma outside the tissue (bg_source says
# "outside_tissue" or the "p5" fallback), p99 = 99th percentile inside the tissue, snr = (p99 - bg) / noise,
# sat_frac = fraction of sampled voxels at the saturation level, and the exposure- and power-normalized
# levels (NaN unless the acquisition stage supplied exposure and laser power). docs/QC.md.
SIGNAL_QC_ROLES = ("dna", "central_element", "foci", "granule", "lamin")
SIGNAL_QC_METRICS = ("offset", "bg", "bg_source", "noise", "p99", "snr", "sat_frac",
                     "bg_per_ms", "p99_per_ms", "bg_per_ms_pct")
SIGNAL_QC = [f"sig_{r}_{m}" for r in SIGNAL_QC_ROLES for m in SIGNAL_QC_METRICS]

# signal_qc stage: one row per unordered pair of resolved channel roles and region (nuclei /
# perinuclear_ring / cytoplasm) on the stage's xy sub-sample: Pearson r of the two channels over the
# region and the least-squares slope of the offset-subtracted role_b on role_a. image_id and the rest of
# the shared block are stamped on by the writer, as for every other table.
CHANNEL_CROSSTALK = ["role_a", "role_b", "region", "n_voxels", "pearson_r", "slope"]

# one row per image (QC summary). The coloc HEADLINE is the shell voxel coloc (SYP<->PGL-1 in the
# lamin-defined perinuclear cytoplasmic shell, SC ribbon excluded) — it separates male>herm on real
# data. The syp_aggregate object metrics are secondary (threshold-fragile).
# n_germline_nuclei equals n_nuclei when the germline stage is off (every nucleus is the analysis
# set); total_germline_length_um is NaN when the axis stage is off; mean_spots is over the analysis set.
IMAGE_SUMMARY = [
    "n_nuclei", "n_germline_nuclei", "mean_spots",
    "total_germline_length_um", "qc_pass", "qc_flags",
    "n_granules",
    "shell_pearson", "shell_manders_m1", "shell_manders_m2",
    "partition_coef", "partition_coef_rot",   # HEADLINE: SYP-3 enrichment in p-granules (exposure-independent)
    "manders_m1_syp_aggregate", "manders_m2_syp_aggregate",
    "frac_granules_overlapping_syp_aggregate", "overlap_pvalue_syp_aggregate",
    "frac_granules_overlapping_sc_ribbon",
    # SC tracer gonad means over germline nuclei (sc_trace stage)
    "mean_sc_total_length_um", "mean_sc_fragmentation_index", "mean_sc_n_fragments_lb",
    # envelope stage
    "n_envelope_fallback", "envelope_vol_ratio", "n_no_envelope", "no_envelope_vol_frac", "n_territories", "ring_thr",
    # audit stage
    "n_lamin_candidates", "n_missed_nuclei", "n_missed_unlabelled",
    # staging stage
    "n_zoned_nuclei", "pachytene_length_um", "off_axis_cut_um", "n_early", "n_mid", "n_late", "n_off_axis",
    "n_territories_on_trace",
    # partition stage (whole shell, and the pooled hand-traced pachytene region when staged)
    "partition_coef_whole", "partition_coef_rot_whole", "partition_coef_zshift_whole", "partition_coef_specific_whole",
    "partition_bg", "partition_coef_early", "partition_coef_specific_early", "partition_coef_mid",
    "partition_coef_specific_mid", "partition_coef_late", "partition_coef_specific_late",
    "partition_coef_pach", "partition_coef_specific_pach",
    # granule_tail stage (lit fraction)
    "tail_n_granules", "tail_n_no_envelope_dropped", "tail_nuclear_syp", "tail_frac_excess_gt_0.25",
    "tail_frac_excess_gt_0.5", "tail_frac_excess_gt_1.0", "tail_excess_p50", "tail_excess_p90", "tail_excess_p99",
    "tail_n_granules_scored",
    # acquisition stage (acquisition.read_exposures, off by default): the complete per-channel record from
    # the .nd2 description as compact JSON, channel names as keys, each with exposure_ms, power, em_filter,
    # ex_filter, dichroic, camera_mode, bit_depth, binning, conversion_gain, readout_speed and modality
    # (null where the description has no such line). The per-role columns of the same stage, exp_ms_<role>
    # and laser_pow_<role> (the original two) followed by em_filter_<role>, camera_mode_<role> and
    # bit_depth_<role>, take their names from the channel map, so they are per-role extras written after
    # the declared columns (docs/QC.md, "Acquisition record").
    "acquisition_json",
    # signal_qc stage (per channel role, in SIGNAL_QC_ROLES x SIGNAL_QC_METRICS order; see above)
    *SIGNAL_QC,
    # partition measurability (partition stage, docs/PARTITION_QC.md): the SYP signal-to-noise the signal_qc
    # stage measured (sig_central_element_snr; NaN when that stage did not run) and whether it clears
    # partition.min_syp_snr (True when the gate is off, the code default). Both NA when the partition stage
    # did not run. An image below the gate is flagged partition:low_syp_snr; its numbers are still written.
    "partition_syp_snr", "partition_measurable",
]

TABLES = {
    "nuclei": NUCLEI,
    "spots": SPOTS,
    "granules": GRANULES,
    "coloc": COLOC,
    "image_summary": IMAGE_SUMMARY,
    "sc_tracks": SC_TRACKS,
    "sc_per_nucleus": SC_PER_NUCLEUS,
    "mask_audit": MASK_AUDIT,
    "zones": ZONES,
    "partition": PARTITION,
    "granule_tail": GRANULE_TAIL,
    "channel_crosstalk": CHANNEL_CROSSTALK,
}

# which stage owns which table (nuclei and image_summary are shared: every stage may append columns).
# A new optional stage registers its table here so it is written (empty) even when the stage is off.
STAGE_TABLES = {
    "spots": ["spots"],
    "granule": ["granules"],
    "coloc": ["coloc"],
    "sc_trace": ["sc_tracks", "sc_per_nucleus"],
    "audit": ["mask_audit"],
    "staging": ["zones"],
    "partition": ["partition"],
    "granule_tail": ["granule_tail"],
    "signal_qc": ["channel_crosstalk"],
}
