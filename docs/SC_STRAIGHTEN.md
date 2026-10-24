# SC ribbon straightening and quantification (`germquant sc-straighten`)

Experimental, version 1.0. Every image it writes carries the flag `sc_straighten:uncalibrated` until the
pre-registered protocol (docs/SC_VALIDATION_PROTOCOL.md) has scored the method on its sealed test set. The design
it implements is docs/SC_STRAIGHTEN_DESIGN.md; the differences are listed at the end of this page.

## What it does

A post-run command on finished fixed-image results (`germquant run` or `batch` output at xy stride 1). For every
germline nucleus it traces each synaptonemal-complex ribbon of the central_element channel (SYP-3 or GFP::SYP-2)
in 3D inside the nucleus territory, refines the centreline, places the ends with an anchored half-maximum rule,
and resamples the ribbon along its centreline into a straight strip (Fiji's Straighten, in 3D, with
rotation-minimising frames). Along each straightened ribbon it measures the length, the intensity profile, the
width, gaps and pieces, and the positions of foci (RAD-51 `spots`, COSA-1 or any other `spots_<name>` table)
along the ribbon. Each strand gets QC flags and an abstaining class (`complete`, `partial`, `merged_suspect`,
`fragment`, `verified`); each nucleus gets a tier that does not need a perfect decomposition (total length,
continuity, the on-strand SYP fraction).

It writes only new files. Nothing in `process_image`, the stage registry, `schema.TABLES`, the config files or
the existing SC tracer changes, so every existing output stays byte-identical (tests/test_sc_straighten_cli.py
checks the sha256 of every pre-existing file of a finished run).

## Running it

```
germquant sc-straighten RESULTS_ROOT [IMAGE_IDS ...] [--out DIR] [--workers N] [--set KEY=VALUE ...]
germquant sc-straighten RESULTS_ROOT n2rad51_noHS_herm_000 --out C:/work/sc_runs --workers 4
germquant sc-view RESULTS_DIR IMAGE_ID [--nucleus N] [--strips] [--stride K]
```

`IMAGE_IDS` are image ids or result-folder names; without them every finished image under the root runs. One
line per image, for example `20251105_N2_nohs_HERM _: 1243 nuclei (traced 663, no_sc 580), 4542 strands
(partial 3041, fragment 713, complete 441, merged_suspect 347), flags: foci_uncalibrated,
label_merged_suspect=217, 1.7 min`. Exit code 0 when at least one image finished or was up to date, 1 when none
did, 2 for usage errors.

Options: `--config` / `--profile` (only their `sc_straighten` section is read), `--set KEY=VALUE` (for example
`follow.seed_min=3.5`; it enters the parameter hash), `--out DIR` (mirror the tree under DIR; required when the
root is inside a golden tree), `--nuclei germline|pachytene|all`, `--nucleus-ids N ...`, `--workers N`
(spawn-safe processes; the outputs are identical for every worker count), `--force` / `--no-resume`,
`--no-strips`, `--no-overview`, `--source-map OLD=NEW` (the source image moved), `--dry-run`, and for imported
geometry `--geometry auto|import|both`, `--paths-from ims:FILE|swc:FILE|csv:FILE ...`,
`--refine-imported on|off`, `--ends fit|as_drawn`. `--allow-sealed --reason TEXT [--audit-log FILE]` runs images
of the sealed SC validation split; only the scoring does that (below).

Resume: an image whose record `<id>__sc_straighten.json` says `done` with the same sub-hash is skipped. The
sub-hash covers the version, the parameters, a digest of this package's source files (`code_sha256` in the
record), the tool versions, the sha256 of nuclei.csv, the labels and the spot tables, the size and first 64 MB of
the source image, the run geometry, the geometry mode and the nucleus selection; a code change therefore reruns
the image.

`sc-view` opens the strands, ends, gaps and foci of one image in napari over the SYP channel (read at xy stride
`--stride`, default 2, for speed; the strands are at full resolution), and with `--nucleus N --strips` a second
viewer with the straightened strips of that nucleus. It needs napari, which lives in the napari-germquant
environment, not in the main one.

## Inputs and refusals

The image id is read from the `image_id` column of `<id>__nuclei.csv` and every other input is found by its exact
name, so ids that end in an underscore or hold spaces work. Inputs: the labels (`<id>__nuclei_labels.tif`, read
page by page into uint16), the run geometry (`<id>__done.json`; else `run_manifest.json`, image flag
`geometry_from_manifest`; else the entry of the result folder in a `golden_manifest.json` one or two levels up,
image flag `geometry_from_golden_manifest`), the channel roles as the run resolved them (`channel_map_json`), the
spot tables, and the source image (`file_path`, through `--source-map`). The central_element channel is read
whole; every other role (the foci channel for profiles and strips) is read plane by plane into per-nucleus crops,
so only one whole channel is in memory.

An image is refused, with a message, and the others run, when: the source image is missing (`source_missing`,
the message names `--source-map`), the run is at xy stride above 1 (`undersampled_xy`), it is a time-lapse series
folder (`series_folder`), the run has no central_element channel (`no_central_element`), there is no run
geometry (`no_run_geometry`), the label image is missing or has another shape (`no_labels`,
`labels_mismatch`), the output folder is inside a golden tree without `--out` (`golden_folder`), or the image is
a test or reserved image of the sealed SC validation split (`sealed`; tests/golden/sc_split.json; image stems are
compared without regard to case). `--allow-sealed` lifts the last refusal only with `--reason`; each sealed image
then writes a start and an end line to the audit log, and its record says `sealed_override: true` with the
reason. Imported geometry from the registered live Filament truth file is refused the same way
(`sealed_import`).

## Outputs (all new files, in the image's results folder or the mirrored folder under `--out`)

| file | one row per | notes |
|---|---|---|
| `<id>__sc_strands.csv` | strand | class, QC, length, ends, shape, profile summaries, gaps, width, foci count, strip offset |
| `<id>__sc_strand_samples.csv` | strand x s every 0.1 um | centreline, frame, core, integrated, foci and dna cores, flags per sample; only `image_id` and keys |
| `<id>__sc_strand_gaps.csv` | gap | s interval, length, depth, z score, whether a traced link lies in it, midpoint, local plateau |
| `<id>__sc_strand_foci.csv` | spot x spot table, in processed nuclei | assignment (unique, ambiguous, off_strand, no_strand), s, folded s, u, v, rank, next-focus distance |
| `<id>__sc_strand_nuclei.csv` | considered nucleus | status (traced, no_sc, skipped, failed), counts, totals, continuity, on-strand fraction, QC |
| `<id>__sc_strand_strips_flat.tif` | | uint16 strips of every strand concatenated along s; axes (role, u, s) |
| `<id>__sc_strand_strips_axial.tif` | | the same with axes (role, v, s) |
| `<id>__sc_strand_overview.png` | | SYP maximum projection with the strands coloured by class |
| `<id>__sc_straighten.json` | image | the record: parameters, inputs and their digests, code digest, noise model, PSF, counts, flags, timing |

Strips: the TIFF description (JSON) lists the roles (`central_element`, then `foci` when the run has it), the
steps and the index `[strand_id, offset, len]`; `strip_offset` and `strip_len` of the strands table give the same.
Columns are s from -0.5 um to length_um + 0.5 um every 0.05 um (0.5 um margins beyond each end). The flat strip's
rows are u (along n_lat) from -0.6 to +0.6 um every 0.05 um (25 rows), averaged over |v| <= 0.3 um; the axial
strip's rows are v (along n_ax, toward +z) from -1.0 to +1.0 um every 0.1 um (21 rows), averaged over
|u| <= 0.15 um. Values are raw counts rounded to uint16 (0 outside the image).

Coordinates are whole-image microns of the run (voxel index times spacing, axis order z, y, x), so the strands,
gaps and foci overlay the run's label image in napari. s runs from the end with the smaller rounded (z, y, x).
The strands, gaps, foci and nuclei tables carry the run's shared metadata block (copied from its nuclei.csv) plus
`sc_straighten_version` and `sc_straighten_hash`. `germquant collect` stacks them into `batch_sc_strands.csv`,
`batch_sc_strand_nuclei.csv`, `batch_sc_strand_gaps.csv`, `batch_sc_strand_foci.csv` and
`batch_sc_straighten_summary.csv` when sc_straighten records exist (restricted to the study's image ids), and
writes nothing new otherwise. A record whose nuclei.csv changed after sc-straighten ran is stale: its tables are
not stacked and its summary row says `stale`.

## Columns

tests/test_sc_straighten_qc.py keeps these tables equal to `germquant.sc.straighten.schema` (declared columns;
the strands, gaps, foci and nuclei tables add the shared metadata block after them).

### `<id>__sc_strands.csv`

| column | meaning |
|---|---|
| strand_id | strand number in the image, assigned by (nucleus_id, length descending, end coordinates) |
| strand_uid | first 10 hex of the sha1 of the in-strand centreline rounded to 0.01 um (stable for the same geometry) |
| nucleus_id | the nucleus label |
| marker | the central_element marker named by the run's channel map |
| geometry_source | auto (traced) or imported (--paths-from) |
| verified | imported geometry marked verified |
| strand_class | complete, partial, merged_suspect, fragment or verified (see Strand classes) |
| qc_status | pass, review or fail: the worst level among qc_flags |
| qc_flags | the strand's QC codes, separated by ; |
| length_um | 3D length between the measured ends (anchored half-maximum end rule) |
| length_xy_um | the same length with z ignored (short for steep strands) |
| length_traced_um | length of the refined centreline before the end rule |
| e0_um | net change of length at end 0 against the traced end (positive: the end lies beyond it) |
| e1_um | the same at end 1 |
| length_se_um | standard error of length_um (both end fits and a curve term) |
| chord_um | straight distance between the two ends |
| straightness | chord_um / length_um |
| curvature_mean_per_um | mean curvature of the centreline (1/um) |
| curvature_max_per_um | largest curvature of the centreline (1/um) |
| tz_abs_mean | mean absolute z component of the unit tangent |
| frac_steep45 | share of the strand steeper than 45 degrees to the image plane |
| frac_steep60 | share of the strand steeper than 60 degrees to the image plane |
| z_extent_um | z range of the strand |
| fit_rms_um | residual RMS of the centreline spline fit |
| n_links | end-matching links joined into the strand |
| link_cost_max | largest link cost (NaN without links) |
| frac_in_label | share of the strand inside its nucleus label |
| frac_near_other_label | share of the strand within 0.3 um of another nucleus label |
| roll_from_axial_med_deg | median angle between n_ax and the optical axis projected into the normal plane |
| end0_class | end 0: free_fit, free_unfit, contact, territory_edge, stack_edge, cycle_cut, max_length or as_drawn |
| end1_class | the same for end 1 |
| end0_z_um | end 0 position (whole-image um, z) |
| end0_y_um | end 0 position (y) |
| end0_x_um | end 0 position (x) |
| end1_z_um | end 1 position (z); s runs from end 0, the end with the smaller rounded (z, y, x) |
| end1_y_um | end 1 position (y) |
| end1_x_um | end 1 position (x) |
| core_median | median core signal over clean inner samples (mean in the ellipse 0.12 x 0.25 um) |
| integrated_median | median integrated signal over clean inner samples (ellipse 0.35 x 0.7 um) |
| snr_core_median | median of core over its standard error |
| profile_cv_corr | noise-corrected coefficient of variation of the integrated profile |
| n_gaps | gaps on the strand |
| gap_um | summed gap length |
| gap_fraction | gap_um / length_um |
| n_pieces | pieces of 0.25 um or more after splitting the strand at its gaps |
| longest_piece_um | longest piece |
| multiplicity_median | median integrated signal over the nucleus's median strand signal |
| fwhm_lat_um | lateral FWHM of the mean cross-section (QC) |
| fwhm_ax_um | axial FWHM of the mean cross-section (QC) |
| width_lat_decon_um | fwhm_lat_um with the PSF FWHM removed in quadrature (QC; the SC is below resolution) |
| width_ax_decon_um | fwhm_ax_um with the PSF FWHM removed in quadrature (QC) |
| frac_contaminated | share of the samples near another strand's centreline (contamination_um) |
| n_foci | foci assigned to the strand (unique or ambiguous) |
| strip_offset | first column of the strand in the strip TIFFs |
| strip_len | number of strip columns of the strand (s from -0.5 um to length_um + 0.5 um every 0.05 um) |
| e0_se_um | standard error of the end 0 position |
| e1_se_um | standard error of the end 1 position |
| end0_sigma_um | fitted erf sigma at end 0 |
| end1_sigma_um | fitted erf sigma at end 1 |
| n_members | traced pieces joined into the strand |
| kink_max_deg | largest turn of the tangent within qc.kink_um |
| frac_saturated | share of core samples at the detector maximum |
| near_stack_edge | a sample lies within qc.stack_edge_um of the first or last plane |
| merged_run_um | longest run with multiplicity at or above classes.merged_multiplicity |
| frac_multiplicity_high | share of samples with multiplicity at or above classes.merged_multiplicity |
| gap_max_um | longest gap |
| end_types_traced | end types after tracing, end matching and accounting (end 0 and end 1, joined by a vertical bar) |
| import_id | id of the imported path (imported geometry) |
| error | exception type when measuring the strand failed (empty otherwise) |

### `<id>__sc_strand_samples.csv`

| column | meaning |
|---|---|
| image_id | the image id of the run |
| strand_id | strand number in the image, assigned by (nucleus_id, length descending, end coordinates) |
| strand_uid | first 10 hex of the sha1 of the in-strand centreline rounded to 0.01 um (stable for the same geometry) |
| nucleus_id | the nucleus label |
| i | sample index along the strand (every 0.1 um) |
| s_um | arc position from end 0 (negative, or beyond length_um, in the 0.5 um margins) |
| in_strand | s_um lies within [0, length_um] |
| z_um | centreline position (whole-image um, z) |
| y_um | centreline position (y) |
| x_um | centreline position (x) |
| t_z | unit tangent (z) |
| t_y | unit tangent (y) |
| t_x | unit tangent (x) |
| nlat_z | lateral normal n_lat = t x n_ax (z); the u axis of the flat strip |
| nlat_y | lateral normal (y) |
| nlat_x | lateral normal (x) |
| nax_z | axial normal of the rotation-minimising frame rolled toward the optical axis (z); the v axis of the axial strip |
| nax_y | axial normal (y) |
| nax_x | axial normal (x) |
| core | core signal (counts above the territory background) |
| core_se | standard error of core (noise model) |
| integrated | integrated signal (counts x um^2) |
| foci_core | core signal of the foci role (NA when not read) |
| dna_core | core signal of the dna role (NA when not read) |
| contaminated | within contamination_um of another strand's centreline |
| in_gap | inside a detected gap |
| in_label | inside the nucleus label |
| strength | ridge strength (noise-null units) |
| lineness | ratio of the largest to the middle Hessian eigenvalue magnitude (low: line-like, high: a junction or a blob) |
| multiplicity | integrated over the nucleus's median strand signal |
| roll_from_axial_deg | angle between n_ax and the optical axis projected into the normal plane |

### `<id>__sc_strand_gaps.csv`

| column | meaning |
|---|---|
| gap_id | gap number in the image |
| strand_id | strand number in the image, assigned by (nucleus_id, length descending, end coordinates) |
| strand_uid | first 10 hex of the sha1 of the in-strand centreline rounded to 0.01 um (stable for the same geometry) |
| nucleus_id | the nucleus label |
| s_start_um | start of the gap on the strand |
| s_end_um | end of the gap on the strand |
| length_um | missing signal integrated over the gap window widened by the blur, as a fraction of the local plateau |
| depth | (plateau - minimum of the smoothed profile) / plateau; above 1 when the local background exceeds the signal there |
| z_score | (plateau - minimum) / its standard error |
| bridged | an end-matching link lies in the gap (the tracer was joined across it) |
| mid_z_um | gap midpoint (z) |
| mid_y_um | gap midpoint (y) |
| mid_x_um | gap midpoint (x) |
| plateau_rel | the local plateau of the gap over the strand's median integrated signal |

### `<id>__sc_strand_foci.csv`

| column | meaning |
|---|---|
| spot_table | spots (RAD-51), or the name of a spots_<name> table (COSA-1 and others) |
| spot_id | the spot's id in its table |
| nucleus_id | the spot's nucleus |
| marker | the spot table's marker |
| assignment | unique, ambiguous (a second strand about as close), off_strand or no_strand |
| strand_id | the closest strand |
| strand_uid | its uid |
| s_um | position along that strand |
| s_frac | s_um / length_um, clipped to [0, 1] |
| s_folded_frac | min(s, L - s) / L (does not depend on the orientation) |
| u_um | lateral offset from the centreline (along n_lat) |
| v_um | axial offset from the centreline (along n_ax) |
| dist_um | distance to the closest strand |
| second_strand_id | the next closest strand |
| in_gap | the spot lies in a gap of its strand |
| rank_on_strand | order along the strand among the assigned foci of the table |
| next_focus_um | distance along the strand to the next assigned focus |
| next_focus_frac | next_focus_um / length_um |
| z_um | spot position (z) |
| y_um | spot position (y) |
| x_um | spot position (x) |

### `<id>__sc_strand_nuclei.csv`

| column | meaning |
|---|---|
| nucleus_id | the nucleus label |
| status | traced, no_sc (no SC above the presence rules), skipped or failed |
| error | exception type of a failed nucleus |
| n_strands | strands after accounting |
| n_main | strands of 1 um or more (the male X rule may drop one) |
| n_complete | complete or verified strands |
| n_fragments | strands under 1 um |
| n_merged_suspect | merged_suspect strands |
| expected | expected SC count: 6 for oocytes, 5 for spermatocytes |
| count_status | low, ok or high: n_main against expected |
| complete_nucleus | every main strand complete, n_complete == expected and no nucleus fail flag |
| total_length_um | summed strand length (every traced nucleus) |
| total_length_xy_um | summed length_xy_um |
| total_length_complete_um | total_length_um of complete nuclei only (NaN otherwise) |
| n_pieces | pieces of 0.25 um or more after splitting the strands at their gaps |
| mean_piece_um | mean of those pieces |
| longest_piece_um | longest piece |
| n_gaps | gaps on the nucleus's strands |
| gaps_per_10um | 10 x n_gaps / total_length_um (NaN under 1 um) |
| gap_fraction | summed gap length / total_length_um |
| piece_lw_mean_um | length-weighted mean piece length, sum(l^2) / sum(l): the piece that holds a random point of the nucleus's SC (the continuity score of the validation protocol) |
| on_strand_frac | SYP above background inside the strand tubes over SYP above background in the territory |
| snr_core_median | median core SNR over the nucleus's strands |
| coverage_proxy | share of strong crest voxels within 0.3 um of a strand |
| multiplicity_frac | share of the strand length with multiplicity at or above classes.merged_multiplicity |
| noise_sd | SD of the transformed noise in the territory (about 1 when the noise model fits) |
| n_seeds | crest seeds of the tracer |
| budget_exceeded | tracing stopped at the work budget |
| volume_um3 | label volume (from the run) |
| in_germline | the run's germline call |
| axis_position_um | position along the gonad axis (from the run) |
| zone | the run's zone (when staging ran) |
| is_pachytene | the run's pachytene call (when staging ran) |
| touches_z | the label touches the first or last plane |
| touches_xy | the label touches the image border |
| qc_status | pass, review or fail: the worst level among qc_flags |
| qc_flags | the nucleus's QC codes, separated by ; |
| error_message | the exception message of a failed nucleus |

## Method in short

Per image: a noise model `var = a + b I` fitted on second differences in the germline labels; a generalised
Anscombe transform to unit noise per nucleus crop; Hessian line fields at about one voxel with a strength
calibrated against unit white noise. Per nucleus: seeded ridge following (predictor-corrector; pure prediction
inside another strand's tube, so crossings pass through), the removal of axial shadows, T contacts (an end that
runs into another strand is trimmed there and becomes a `contact` end), one global end matching (networkx
maximum-weight matching over gaps and crossings; overlapping tips are joined without a backward step), duplicate
removal within and across nuclei, after which a `contact` end whose partner was removed is a free end again;
normal-plane re-centring on raw counts and a penalised cubic spline. Per image, a PSF estimate from flat isolated
bright stretches (fallback 0.365 / 0.15 um sigma with the flag `psf_fallback`). Per strand: the anchored erf end
rule (when the fitted edge sits at the inner bound of its window the end is fitted again from an anchor 0.6 um
deeper), rotation-minimising frames with one roll toward the optical axis, strips, profiles, gaps, widths,
multiplicity, foci, flags and the class.

Gaps: a candidate is a run of the smoothed integrated profile below 0.9 of the strand's median; its window is
widened by one PSF sigma along the strand on each side, and the dip is judged against a straight baseline between
the profile medians just outside the window (0.4 um on each side), because the integrated signal of a ribbon
changes along it with its brightness and its tilt. A gap is a dip below 0.6 of the baseline, deeper than 3
standard errors; its length is the missing fraction of the baseline integrated over the window (blurring along
the strand conserves signal). On development phantoms (5 um arcs, peak 450) gaps of 0.3 to 1.0 um measure within
0.02 um at 0 and 30 degrees of tilt, and about 0.1 um short at 60 degrees, where a 0.3 um gap is not detected.

Details and the reasons for each choice: docs/SC_STRAIGHTEN_DESIGN.md.

## Parameters

A top-level `sc_straighten` section, read only by this command (no shipped config or profile names it, so no
config_hash moves): code defaults < the section of `--config` or `--profile` < `--set`. The keys and defaults
are in `germquant/sc/straighten/params.py` (DEFAULTS). The keys a development sweep may tune are
`follow.seed_min`, `follow.lo`, `link.*`, `refine.knot_um`, `refine.lambda`, `refine.sigma_loc_um`, `gaps.*` and
`classes.*`; everything else is fixed unless a changelog entry of the validation protocol says otherwise.
Development calibration: `scripts/calibrate_sc_straighten.py` (synthetic development seeds only; `--generator
registered` draws the dense and snr phantoms without the 0.25 um clearance of protocol changelog entry 2(c)).

## QC flags

`qc_status` is `fail` when any fail flag is present, else `review` when any review flag is present, else `pass`.
tests/test_sc_straighten_qc.py keeps this table equal to `germquant.sc.straighten.qc.CODES` and gives every code a
positive and a negative case.

| code | scope | level | threshold key | fires when |
|---|---|---|---|---|
| fragment | strand | info | classes.min_complete_um | length_um below the minimum of a whole strand |
| low_snr | strand | review | qc.low_snr_review | median core SNR below the review level |
| low_snr_fail | strand | fail | qc.low_snr_fail | median core SNR below the fail level (replaces low_snr) |
| end_unfit | strand | review | ends.fit_min_amp_se | a free end whose half-maximum fit failed (the traced end is kept) |
| contact_end | strand | info | - | an end where the strand meets another strand |
| truncated_z | strand | fail | qc.stack_edge_um | an end at the stack edge, or a sample within stack_edge_um of the first or last plane |
| costly_link | strand | review | classes.max_link_cost | a link of the end matching costs more than max_link_cost |
| contaminated | strand | review | classes.max_contaminated | more than max_contaminated of the samples lie near another strand |
| merged_suspect | strand | review | classes.merged_multiplicity | integrated signal at or above merged_multiplicity x the nucleus's median strand signal over at least merged_run_um |
| double_width | strand | review | qc.double_width_ratio | lateral or axial FWHM above double_width_ratio x the image's complete-strand median |
| exits_label | strand | review | classes.min_frac_in_label | less than min_frac_in_label of the strand inside its nucleus label |
| neighbour_nucleus | strand | fail | qc.near_other_label_fail | more than near_other_label_fail of the strand near another nucleus label |
| steep_z | strand | review | qc.steep_frac_review | more than steep_frac_review of the strand steeper than 60 degrees to the image plane |
| long | strand | review | qc.long_um | longer than long_um: two chromosomes joined end to end |
| kink | strand | review | qc.kink_deg | a turn above kink_deg within kink_um |
| fit_residual | strand | review | qc.fit_rms_review_um | centreline fit residual RMS above fit_rms_review_um |
| gap_long | strand | review | qc.gap_long_um | a gap of gap_long_um or more |
| saturated | strand | review | qc.saturated_frac | more than saturated_frac of the core samples at the detector maximum |
| imported | strand | info | - | geometry imported from a file (Imaris Filaments, SWC or CSV) |
| verified | strand | info | - | imported geometry marked verified |
| error | strand | fail | - | an exception while measuring the strand (error:<Type>) |
| label_merged_suspect | nucleus | review | qc.merged_label_ratio | label volume above merged_label_ratio x the median of the traced nuclei within merged_label_window_um along the gonad axis (the image median with fewer than 5 of them or no axis position) |
| truncated_z | nucleus | fail | - | the label touches the first or last plane (totals unusable) |
| edge_xy | nucleus | review | - | the label touches the image border |
| count_low | nucleus | info | nucleus.expected | fewer main strands than expected |
| count_high | nucleus | info | nucleus.expected | more main strands than expected |
| noise_mismatch | nucleus | review | noise.nucleus_sd_range | the transformed noise SD of the territory outside nucleus_sd_range |
| budget_exceeded | nucleus | review | budget_s | tracing stopped at the work budget (strands so far kept) |
| error | nucleus | fail | - | an exception in the nucleus (error:<Type>; status failed, no zero row) |
| sc_straighten:uncalibrated | image | info | - | always in v1, until the validation protocol has scored the method |
| noise_fit_poor | image | review | noise.poor_rel_rms | the noise model fits its bins above the darkest decile of the sampled voxels with a relative RMS above poor_rel_rms (or fell back) |
| psf_fallback | image | review | psf.min_stretches | fewer than min_stretches flat isolated bright stretches: the fallback PSF was used |
| sc_channel_suspect | image | review | qc.sc_channel_suspect_frac | among nuclei with SC present, fewer than sc_channel_suspect_frac have a strand of 2 um or more |
| foci_uncalibrated | image | info | - | a spot table was used (spot detection is not calibrated for SC positions) |
| geometry_from_manifest | image | info | - | run geometry (xy_stride, z_range) from run_manifest.json: no done marker |
| geometry_from_golden_manifest | image | info | - | run geometry from the golden_manifest.json above the results folder: no done marker and no xy_stride in run_manifest.json |
| label_merged_suspect | image | info | - | label_merged_suspect=N: the number of such nuclei |
| nucleus_errors | image | info | - | nucleus_errors=N: the number of failed nuclei |

## Strand classes (abstain rather than mislead)

`verified` (imported with verified=true), `fragment` (under 1 um), `merged_suspect`, `complete` (both ends
`free_fit`, link cost at most 60, at most 30 percent contaminated, at least 80 percent in its label, no review or
fail flag) and `partial` (everything else). Nucleus totals use every strand after the accounting step; nuclei
without SC keep their weak strands and their total (a status never zeroes a column).

Per-chromosome statistics are not supported yet. The `complete` class is the only candidate for them, and its
precision is measured only on synthetic data: high for isolated ribbons, about 0.7 in dense phantom nuclei, and
unknown on real nuclei until the fixed-image hand traces of the protocol (S4) exist. On the development gonads
almost no nucleus is a `complete_nucleus`, so `total_length_complete_um` is nearly always empty.

## Accuracy notes (development phantoms)

* Length: isolated ribbons within 30 degrees of the image plane measure within about 0.05 um; tilted ribbons
  (30 to 90 degrees) at real SNR measure about 0.1 to 0.2 um long, mostly in the end placement (the bias
  disappears at ten times the SNR). Arcs curved in a vertical plane measure 0.07 to 0.17 um short over 4 um
  (the axial PSF). Ends kept where the trace stopped (`contact`, `territory_edge`, `stack_edge`, `free_unfit`)
  carry a spline overshoot of about 0.03 to 0.05 um each.
* `length_xy_um` ignores z and runs a little short at steep tilt (-0.07 um at 45 degrees to -0.14 um at 90 degrees
  over a few um).
* Dense nuclei: the totals run about 5 to 7 percent short with a spread of 8 to 9 percent; per-ribbon
  decomposition is rarely complete (see Strand classes).

## napari

`germquant.sc.straighten.layers` returns napari LayerData tuples without importing napari: the strands as a
Shapes path layer (features = the strands table, coloured by class), end points, gap midpoints and foci; a gallery
of one nucleus's straightened strips; an image crop in the run frame. `germquant sc-view` opens them.

## Validation status

Not validated. The pre-registered protocol (docs/SC_VALIDATION_PROTOCOL.md) fixes a development set and a sealed
test set of images, synthetic development and sealed seeds, the metrics and the success criteria; the test set is
scored once, later, by the orchestrator (`scripts/benchmark_sc.py zones|synthetic|filaments --final --reason
TEXT`, with the frozen parameters and thresholds of the freeze entry; every --final run is logged). Development
results are kept under `C:/Users/ryane/sc_ribbon_work/dev_eval/` (summary_<commit>.json); they are not test
results and claim nothing.

Development snapshot at the build, before the review fixes (commit 3e88403; synthetic development seeds 2000 to 2049, 50 scenes per
family; the eight 2025-11-05 development gonads):

* Synthetic criteria of protocol 8.1 met on development seeds: S1a (isolated ribbons up to 60 degrees: mean
  length error +0.05 um, 95 percent limits -0.17 to +0.27 um, median centreline distance 0.02 um, one strand per
  ribbon in all), S1c (gap precision 1.00, recall 0.84, gap length error 0.08 um), S1d (foci: 90th percentile of
  the s error 0.07 um, all on the right ribbon), S1f (no noise-only nucleus reached 1 um). Not met: S1b (dense
  nuclei: mean total -4.7 percent, SD 7.9 percent against 7, CCC 0.57 against 0.80; coverage 0.97 and off-truth
  0.01 pass) and S1e (complete precision 0.89 against 0.90). The dense under-count of 6 to 7 percent does not shrink
  with brightness (scene peak 300, 430, 700), and reaches 15 percent at peak 200: no SNR floor within 5 percent.
  Steep isolated ribbons (60 to 75 degrees) are over-measured (flag steep_z).
* Zones (protocol S2, development gonads, three trials logged in dev_trials.jsonl): SC present against None by
  total_length_um, mean per-gonad AUC 0.983 (bootstrap interval 0.977 to 0.989, eight gonads), level with the raw
  SYP intensity (0.985); SCfull against fragmented by piece_lw_mean_um, mean AUC 0.77 (0.40 to 0.99) over the three
  gonads with at least ten nuclei in each class, below the intensity comparator (0.96); SC_fullstart found within
  one point spacing in 6 of 8 gonads, median error 21 um (trial 2, commit f7dc9a3: 7 of 8, 11 um).

## Differences from the design document

* The package is `germquant.sc.straighten` (next to the SC tracer), not `germquant.sc_straighten`; the CLI module
  is `germquant/cli_sc_straighten.py`. The protocol's generator is therefore `germquant.sc.straighten.synth`.
* Synthetic dense nuclei keep their ribbons at least 0.25 um apart in 3D (two chromosomes' SCs do not pass through
  each other); without that rule 13.5 percent of the truth length lay within half a PSF FWHM of another ribbon.
  `--generator registered` draws them without it, and the scoring reports S1b and S1e on both.
* End matching works from anchors 0.3 um inside each tip (the PSF-blurred tip drifts along the optical axis),
  with the direction from the next 0.9 um and the gap estimated tip to tip; tips that meet or overlap are judged
  by their collinearity and joined without a backward step; two `contact` ends facing each other may link across
  up to 1.5 um (a strand that crosses between them); an unlinked end costs 60 (development calibration: curved
  ribbons broken at sub-micron gaps were left unlinked at 40).
* Accounting also removes axial shadows (short steep strands standing on another strand inside the axial PSF;
  before the T contacts and the end matching, so a shadow never trims or joins the strand it stands on, and again
  after them) and shadowed strands (at least 80 percent of the length inside the contamination reach of a longer
  strand). A `contact` end with no surviving strand within `follow.contact_free_um` (0.6 um axial, 0.4 um lateral)
  is free again.
* A trace that reaches the territory or stack edge after its ridge has already faded ends `free` (trimmed to its
  last strong point, so the end rule applies); before, a bright steep ribbon's axial PSF tail could carry the trace
  1 um or more past the true end.
* The tracer follows by pure prediction while inside another strand's claim tube (crossings pass through).
* Gap length is the missing integrated signal over the widened dip window against a local baseline (see Method);
  the design's 0.35 threshold of the strand median missed most sub-micron gaps of tilted ribbons.
* The P-spline uses knots every 0.5 um and lambda 0.1 (development calibration: a strong second-difference
  penalty bends a clamped spline outward at its ends and adds length); imported paths measured as drawn use a
  light penalty (`refine.lambda_imported` 0.01).
* An end fit that would extend the strand more than `ends.max_extension_um` (0.5 um) beyond its traced end is
  `free_unfit` (no long straight extrapolation); an edge at the inner bound of the profile window is fitted again
  from a deeper anchor.
* Foci: the closest strand and the ambiguity test use the distance to the strand proper (lateral, axial and past
  the end), so a strand that ends short of a focus and only points at it is not a candidate.
* Per-nucleus timings are in the record (`timing`), not in the nuclei table, so the tables are identical for
  every worker count.
* The multiplicity reference is the nucleus's own median strand signal (the image's complete-strand median when
  the nucleus holds under 3 um of SC): SC brightness changes along the gonad, and one image-wide reference called
  most strands of the bright proximal pachytene region merged on n2rad51_noHS_herm_000.
* `label_merged_suspect` compares a label with the median volume of the traced nuclei within 15 um along the gonad
  axis (nuclei grow along the gonad); `noise_fit_poor` judges the noise fit above the darkest decile of the
  sampled voxels (near the camera offset a read-noise floor bends the variance curve, and no ridge is measured
  there).
* Not built in v1 (DESIGN 1.3): the napari review dock and edit replay, time-lapse series, the registry stage and
  its profile, deconvolution, the per-gonad zone estimator as a production output, api.py entry points; and the
  real regression goldens of two development tags.
