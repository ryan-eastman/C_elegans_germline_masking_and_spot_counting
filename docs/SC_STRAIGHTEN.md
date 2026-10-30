# SC ribbon straightening and quantification (`germquant sc-straighten`)

Experimental, version 2.0. Every image it writes carries the flag `sc_straighten:uncalibrated` until the
pre-registered protocol (docs/SC_VALIDATION_PROTOCOL.md) has scored the method on its sealed test set. The design
it implements is docs/SC_STRAIGHTEN_DESIGN.md; the differences are listed at the end of this page.

Version 2.0 is the v2 development round of the protocol (changelog entries 10 to 26, and 27 to 29 after an
independent review of the round). Every change it tried is a switch, and none met its pre-registered adoption
rule, so every switch is off by default and the defaults give the outputs of version 1.0 in every column of 1.0
(tests/test_sc_straighten_v2_compat.py); see "Version 2.0 switches". New in 2.0: the strand column `punctum_sigma_z_um` and the nucleus columns `n_puncta`,
`puncta_signal_frac` (both 0 with the punctum class off) and `steep_length_frac`.

## What it does

A post-run command on finished fixed-image results (`germquant run` or `batch` output at xy stride 1). For every
germline nucleus it traces each synaptonemal-complex ribbon of the central_element channel (SYP-3 or GFP::SYP-2)
in 3D inside the nucleus territory, refines the centreline, places the ends with an anchored half-maximum rule,
and resamples the ribbon along its centreline into a straight strip (Fiji's Straighten, in 3D, with
rotation-minimising frames). Along each straightened ribbon it measures the length, the intensity profile, the
width, gaps and pieces, and the positions of foci (RAD-51 `spots`, COSA-1 or any other `spots_<name>` table)
along the ribbon. Each strand gets QC flags and an abstaining class (`complete`, `partial`, `merged_suspect`,
`fragment`, `verified`, and `punctum` when its switch is on); each nucleus gets a tier that does not need a
perfect decomposition (total length, continuity, the on-strand SYP fraction).

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
line per image, for example `20251105_N2_nohs_HERM _: 1243 nuclei (traced 663, no_sc 580), 4528 strands
(partial 3015, fragment 702, complete 471, merged_suspect 340), flags: foci_uncalibrated, label_merged_suspect=99,
2.8 min`. Exit code 0 when at least one image finished or was up to date, 1 when none
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
compared without regard to case). An image whose file name or folder name is a sealed image stem or RAD-51 tag is
refused before any of its files is read, and discovery does not read its nuclei.csv (protocol changelog entry 27).
`--allow-sealed` lifts the last refusal only with `--reason`; each sealed image
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
| strand_class | complete, partial, merged_suspect, fragment, punctum or verified (see Strand classes) |
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
| punctum_sigma_z_um | axial extent of a steep strand at most punctum.max_xy_um long in the image plane: the sigma of a Gaussian fitted along s to its core profile times tz_abs_mean (the punctum rule; NaN where the rule does not apply) |

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
| n_strands | strands after accounting (punctum strands not counted under punctum.exclude) |
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
| n_puncta | punctum strands (the axial blur of SYP puncta or aggregates; 0 unless punctum.enabled); under punctum.exclude they are left out of every other count, length, continuity and intensity column |
| puncta_signal_frac | SYP above background inside the punctum tubes (outside the strand tubes) over SYP above background in the territory (0 unless punctum.enabled) |
| steep_length_frac | share of total_length_um steeper than 60 degrees to the image plane |

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
`classes.*`; everything else is fixed unless a changelog entry of the validation protocol says otherwise. The
on/off switches of versions 2.0 and 3 (`punctum.one_gaussian`, `gaps.dim_stretch`) take true, false, on, off,
yes or no, and `follow.contact_test` (local, claim_tangent), `link.junction_gap` (anchors, tips) and
`punctum.sigma_source` and `gaps.dim_sigma_source` (psf, foci, foci_depth) their listed values; any other value is
refused (protocol changelog entry 27), so a typo never gives another configuration.
Development calibration: `scripts/calibrate_sc_straighten.py` (synthetic development seeds only; `--generator
registered` draws the dense and snr phantoms without the 0.25 um clearance of protocol changelog entry 2(c)).

## Version 2.0 switches (all off by default)

Each change of the v2 round has a switch whose default gives version 1.0 (`V1_VALUES` in
`germquant/sc/straighten/params.py`; `params.v1_sets()` gives the `--set` lines). Each was designed on the
development seeds 2000 to 2049 and 3000 to 3049, tuned on 10000 to 10099 and, where the tuning met its targets,
confirmed once on 11000 to 11099 (per-image PSF, punctum label, complete guard), by the rule of protocol changelog
entry 11 (b): its targeted metric improves by the amount its entry states, and no criterion of protocol 8.1 gets
worse beyond the tolerances of that entry. None was adopted.

| switch (default) | when on | entries | why it is off |
|---|---|---|---|
| `psf.foci: off` | the image PSF's sigma_z from compact isolated RAD-51 foci (scaled to the central_element emission wavelength; flag `psf_from_foci`), else from flat stretches | 12 to 16 | on the confirmation seeds the median per-image error of sigma_z was +3.4 percent at sigma_z 0.33 um (within 3 required; +0.9, -0.4 and -1.8 percent at 0.38, 0.43 and 0.48 um, against +10.6 to -24.0 percent with the fallback) |
| `punctum.enabled: false` (with `punctum.exclude`, `punctum.before_linking`, `punctum.max_xy_um`, `punctum.sigma_ratio`) | class and flag `punctum` for a steep strand at most 0.5 um long in the image plane whose axial extent fits the axial PSF (the blur of a SYP punctum or aggregate); with `exclude` such strands leave the nucleus tier | 16 to 19 | removing them from the tier cut the excess length of dense nuclei with puncta by 36 and 54 percent (60 required); as a label only, the confirmation labelled one strand of a gapped ribbon (none allowed) |
| `follow.contact_test: local` (`claim_tangent`) | the tracer's contact stop compares the trace with the tangent of the strand that claimed the voxel instead of the local line direction | 20, 21 | on the tuning seeds breaks after the tracer fell by 18 percent (30 required) and strand switches rose by 13.5 percent (10 allowed): most traces released at a crossing do not regain their ribbon behind it (dense totals -6.9 to -1.4 percent, not a targeted metric) |
| `link.junction_dirs: false` (with `link.junction_gap`, `link.max_chain_um`) | end directions measured beyond the points of an end's last 1 um that lie within one PSF FWHM of another strand; such pairs judged by the distance between their tips; optionally a bound on chain length | 21 to 23 | as first specified it made wrong links (gaps estimated from shifted anchors); revision 1 cut wrong links by 31 percent but lost correct ones (179 to 172) and raised the dense SD (8.6 to 8.8 percent) |
| `classes.complete_guard: false` (with `classes.complete_guard_fwhm` 2.0, `classes.max_complete_um` 9.0) | a `complete` strand also needs a length of at most 9 um, no end within the contamination reach of another strand, and no other strand's end along it | 23 to 25 | on the confirmation seeds the pooled precision of `complete` rose from 0.88 to 0.94, but in dense nuclei only from 0.67 to 0.74 (0.80 required) |

Two further changes of the plan were built and run on the design seeds only and not proposed (entries 23 and 24): a
junction-aware end fit (it shortened ends near other strands) and the re-pairing of branches at junction zones (it
did not reduce the turns and end-to-end joins of the final strands). The development phantom families of the
round (`psf_check`, `puncta_pt`, `puncta_agg`, `dense_puncta`, `real_dense` and, recalibrated by protocol changelog
entry 28, `real_dense_r2`; germquant.sc.straighten.synth_dev) are not registered and have no sealed seeds.

## Version 3 rules (both off by default)

The v3 development round (protocol changelog entries 31 to 33) adds two rules, each behind switches whose defaults
give the outputs of version 2.0 (which equal those of 1.0 in every column of 1.0): with the switches off every
column of every table, the strips and the classes equal those of the v2 head (tests/test_sc_straighten_v3_compat.py
against tests/golden/sc_straighten/v2_compat.json). No column is added. `params.V3_PUNCTUM`, `params.V3_GAPS` and
`params.v3_sets()` give the configurations. Neither rule was adopted: both failed targets of their entries on the
tuning seeds and were withdrawn (protocol changelog entry 34; see "Validation status"), so their switches stay off
and the configuration proposed for the freeze is that of version 2.0.

| switches (defaults) | when on | entry |
|---|---|---|
| `punctum.enabled: true`, `punctum.exclude: true`, `punctum.sigma_source: foci_depth` (default `psf`), `punctum.one_gaussian: true` (default false); fixed `punctum.max_rel_rms` 0.2, `psf.ref_depth_min_fits` 30, `psf.ref_depth_min_span_um` 3.0 | the punctum class B1 with the image's measured axial sigma: a steep strand (frac_steep60 above 0.5) at most `punctum.max_xy_um` (0.5 um) long in the image plane whose core profile along s is one Gaussian (its sigma inside the fit bounds, its centre on the strand, residual RMS at most 0.2 of its amplitude) no wider in z than `punctum.sigma_ratio` (1.3) times the reference axial sigma at the strand's depth is class and flag `punctum`; it stays in the strand, sample and strip outputs and leaves every ribbon column of the nucleus tier (`n_puncta` and `puncta_signal_frac` report it). The reference (`psf.axial_reference`) is the sigma_z of the image's compact isolated RAD-51 foci (the foci estimate of entry 12, scaled to the central_element emission wavelength), with a depth term (the Theil-Sen slope of the foci's sigma_z on their depth) when at least 30 kept foci span at least 3 um in z; without a foci estimate it is the image PSF's sigma_z, as in v2. The image PSF itself stays that of v1; the record gives the reference (`punctum_reference`). | 32 |
| `gaps.dim_stretch: true` (default false), `gaps.dim_sigma_source: foci_depth` (default `psf`); `gaps.dim_fill` and `gaps.dim_z` (chosen on the tuning seeds by the rule of entry 33) | dim stretches are not gaps: a gap candidate that passes every v1 test is a dim stretch, and the ribbon continues over it, when (i) no sample of its window lies within `contamination_um` of another strand or within 0.1 um of a link between two traced members, (ii) its relative minimum lies above the minimum that an empty gap of its missing length shows after the blur along the strand (the axial sigma of the measured reference at its depth, the lateral sigma of the image PSF, and the 0.15 um smoothing) by more than `dim_z` standard errors, and (iii) its fill fraction (the share of its flanks' signal that a stretch of that missing length and blurred minimum keeps) exceeds `dim_fill`. An empty gap (fill 0) stays a gap whatever its length and tilt. The strands, ends, links and lengths are those of v1; only the columns that the reading of gaps feeds change. The record counts the dim stretches (`dim_stretches`) and names the reference (`dim_reference`); `ImageResult.dim_stretches` lists them for development scripts. Known ambiguity: the junction of two different ribbons whose tips are dim, with a short empty core that the tracer follows, can look like a filled stretch behind the blur (entry 33 counts such joins in the dense phantoms). | 33 |

Two development families (not registered, no sealed seeds) check the rules against a known PSF: `puncta_psf`
(entry 32) holds dense nuclei with puncta, nuclei of puncta alone and steep isolated arcs, and `gaps_psf` (entry 33)
tilted isolated ribbons with empty gaps of 0.3 to 1.0 um and with dim stretches of 0.4 to 1.5 um at 0.15 to 0.6 of
their brightness; both are rendered with an axial PSF of 0.33 to 0.48 um that changes with depth and hold RAD-51-like
foci with their spot table and wavelengths, so that the measured reference can be checked against the PSF the
scene was rendered with.

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
| punctum | strand | info | punctum.sigma_ratio | a steep strand at most punctum.max_xy_um long in the image plane whose axial extent fits the reference axial sigma (the image PSF's, or with punctum.sigma_source foci or foci_depth the image's compact RAD-51 foci's, version 3 rules) and, with punctum.one_gaussian, whose profile along it is one Gaussian: the axial blur of a SYP punctum or aggregate (class punctum; only with punctum.enabled, off by default) |
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
| sc_straighten:uncalibrated | image | info | - | always, until the validation protocol has scored the method |
| noise_fit_poor | image | review | noise.poor_rel_rms | the noise model fits its bins above the darkest decile of the sampled voxels with a relative RMS above poor_rel_rms (or fell back) |
| psf_fallback | image | review | psf.min_stretches | fewer than min_stretches flat isolated bright stretches and no foci estimate: the fallback PSF was used |
| psf_from_foci | image | info | psf.foci_min | the image PSF's sigma_z is the estimate from compact isolated RAD-51 foci (scaled to the central_element emission wavelength) |
| sc_channel_suspect | image | review | qc.sc_channel_suspect_frac | among nuclei with SC present, fewer than sc_channel_suspect_frac have a strand of 2 um or more |
| foci_uncalibrated | image | info | - | a spot table was used (spot detection is not calibrated for SC positions) |
| geometry_from_manifest | image | info | - | run geometry (xy_stride, z_range) from run_manifest.json: no done marker |
| geometry_from_golden_manifest | image | info | - | run geometry from the golden_manifest.json above the results folder: no done marker and no xy_stride in run_manifest.json |
| label_merged_suspect | image | info | - | label_merged_suspect=N: the number of such nuclei |
| nucleus_errors | image | info | - | nucleus_errors=N: the number of failed nuclei |

## Strand classes (abstain rather than mislead)

`verified` (imported with verified=true), `fragment` (under 1 um), `merged_suspect`, `complete` (both ends
`free_fit`, link cost at most 60, at most 30 percent contaminated, at least 80 percent in its label, no review or
fail flag) and `partial` (everything else); `punctum` only with `punctum.enabled`, and the stricter `complete`
only with `classes.complete_guard` (both off; see Version 2.0 switches). Nucleus totals use every strand after the
accounting step; nuclei without SC keep their weak strands and their total (a status never zeroes a column).

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
* Dense nuclei: the totals run about 7 percent short on fresh seeds (6.7 to 8.3 percent over three blocks of 100
  scenes; 4.8 to 4.9 percent on the design seeds 2000 to 2049 and 3000 to 3049) with a spread of 8 to 10 percent,
  and 14 to 15 percent short when ribbons may pass through each other (the generator as registered); per-ribbon
  decomposition is rarely
  complete (see Strand classes). The v2 round traced the losses to junctions (ends inside another ribbon's blur
  come out short, tracer breaks at crossings, links across junctions); none of its changes fixed them.
* Real data: 39 to 46 percent of the traced strand length on the development gonads runs steeper than 60 degrees
  to the image plane, against 13 percent for randomly oriented ribbons; many short steep strands are probably the
  axial blur of SYP puncta or aggregates. The registered phantoms hold no puncta; the realism phantoms do (below).
* Continuity and brightness (protocol changelog entries 28 and 29; development phantoms calibrated on truth-free
  statistics of the development images): with the brightness of each ribbon varying along it as much as on real
  in-plane strands, the median `piece_lw_mean_um` of dense phantoms falls from 5.0 to 3.6 um (real full SC: 2.6
  um), because dim stretches of continuous ribbons are detected as gaps (gaps per nucleus 1.8 to 4.3) while the
  strands stay whole. On real data part of the low continuity of full SC may therefore be dim stretches read as
  gaps; no change of the method was made for it in version 2.0. With all four realism features (PSF, brightness,
  texture, puncta) the dense totals have an SD of 16 percent and a CCC of 0.18.

## napari

`germquant.sc.straighten.layers` returns napari LayerData tuples without importing napari: the strands as a
Shapes path layer (features = the strands table, coloured by class), end points, gap midpoints and foci; a gallery
of one nucleus's straightened strips; an image crop in the run frame. `germquant sc-view` opens them.

## Validation status

Not validated. The pre-registered protocol (docs/SC_VALIDATION_PROTOCOL.md) fixes a development set and a sealed
test set of images, synthetic development and sealed seeds, the metrics and the success criteria; the test set is
scored once, later, by the orchestrator (`scripts/benchmark_sc.py zones|synthetic|filaments --final --reason
TEXT`, with the frozen parameters and thresholds of the freeze entry; every --final run is logged). Development
results are kept under `C:/Users/ryane/sc_ribbon_work/dev_eval/` (summary_<commit>.json); they are not test results
and claim nothing.

The intensity comparator of S2 (`central_element_mean_intensity`) measures the amount of SYP in a nucleus, not the
continuity of its SC (protocol changelog entry 10): on development phantoms whose ribbons were cut into pieces it
separated continuous from cut ribbons with an AUC of 0.55 to 0.61, against 0.99 to 1.00 for `piece_lw_mean_um`, so
its lead on the zone labels does not show that intensity reads continuity.

Development snapshot of version 2.0 (protocol changelog entry 26; computed once on the fresh seeds 12000 to 12099,
`dev_eval/summary_a711fe5.json`): the configuration proposed for the freeze (every switch off) and version 1.0 give
the same numbers (version 1.0 there is the code of 2.0 with every switch at its v1 value, which
tests/test_sc_straighten_v2_compat.py ties to commit 5cb3cb6, and an independent recompute with the v1 code gave
the same dense S1b; protocol changelog entry 27). Met: S1a (isolated ribbons: mean length error +0.076 um, 95
percent limits -0.07 to +0.22 um, median centreline distance 0.017 um, one strand per ribbon in all), S1c (gap
precision 1.00, recall 0.91, gap length error 0.034 um), S1d (0.064 um, all on the right ribbon) and S1f; not met:
S1b (dense nuclei: mean total -6.7 percent, SD 8.7 percent, CCC 0.49; coverage 0.96 and off-truth 0.01 pass) and
S1e (complete precision 0.890 pooled against 0.90; 0.67 in dense nuclei). With the registered generator the dense totals run 15.3 percent short (SD 9.9 percent,
CCC 0.22). No zone checkpoint was run in v2, since no change that reaches an S2 input was adopted. The hold-out
check H1 (entry 26) had one gonad with a germquant run, 20251021_n2_nohs_HERM_2, so it decides nothing: S2b
0.989; S2a not evaluable under the registered guard and 0.739 under guard_frac 0.25 (intensity comparator 0.670);
SC_fullstart within one point spacing with the thresholds of the twelve (out of sample).

An independent review of the round (protocol changelog entry 27) found no blocker; it disclosed one zone
computation made outside a checkpoint (it compared no configuration), corrected the reports above, and led to
code-only changes of the scoring path and the input guards, and to a recalibrated realism probe (entries 28 and
29), whose brightness feature reproduces the low continuity of real full SC (see Accuracy notes); by entry 11 (e)
nothing in the method changed after the hold-out check.

Development snapshot after the review fixes (commit 6b7be22; synthetic development seeds 2000 to 2049 and the
fresh seeds 3000 to 3049, 50 scenes per family; twelve development gonads, the eight of 2025-11-05 and the four
heat-shock males of 2026-01-05; development trial 4 of dev_trials.jsonl). Values are given as seeds 2000 to 2049 /
seeds 3000 to 3049.

* Synthetic criteria of protocol 8.1 met: S1a (isolated ribbons up to 60 degrees: mean length error +0.078 /
  +0.079 um, 95 percent limits -0.09 to +0.24 / -0.06 to +0.22 um, median centreline distance 0.02 um, one strand
  per ribbon in all; by tilt +0.04 / +0.06 um within 30 degrees, +0.11 / +0.12 um at 30 to 60 degrees and +0.09 /
  +0.13 um at 60 to 75 degrees), S1c (gap precision 1.00 / 1.00, recall 0.88 / 0.89, gap length error 0.03 /
  0.06 um), S1d (foci: 90th percentile of the s error 0.07 / 0.07 um, all on the right ribbon) and S1f (no
  noise-only nucleus reached 1 um). Not met: S1b (dense nuclei: the mean total of -4.9 / -4.8 percent, the
  coverage of 0.97 / 0.96 and the off-truth length of 0.01 pass, the SD of 8.1 / 9.0 percent against 7 and the CCC
  of 0.58 / 0.59 against 0.80 do not; on the fresh blocks of 100 scenes of the v2 round the mean total was -6.9,
  -8.3 and -6.7 percent, outside its band, so the mean passes on these two blocks only, protocol changelog entry
  27) and S1e (complete precision 0.90 / 0.87 pooled against 0.90; 0.71 / 0.58 in
  dense nuclei). With the generator as registered (dense ribbons may pass through each other) the dense totals run
  13.9 percent short (SD 10.3 percent, CCC 0.27). The dense under-count is 6 to 9 percent at scene peaks of 300 to
  700 counts and 15 percent at 200: no SNR floor within 5 percent on seeds 2000 to 2049 (700 counts on seeds 3000
  to 3049).
* Zones (protocol S2 on development gonads, trial 4): SC present against None by total_length_um, mean per-gonad
  AUC 0.987 (bootstrap interval 0.982 to 0.992, twelve gonads), level with the raw SYP intensity (0.988); SCfull
  against fragmented by piece_lw_mean_um, mean AUC 0.83 (0.67 to 0.93) over the seven gonads with at least ten
  nuclei in each class (per gonad 0.45 to 0.98), below the intensity comparator (0.97); SC_fullstart found in all
  twelve gonads and within one point spacing in 11 of 12, median error 5.3 um, with thresholds chosen on the same
  gonads (theta_full 2.0 um; theta_present 10.0 um, the top of its registered grid).
* On the development gonads 46 percent of the traced strand length runs steeper than 60 degrees to the image plane
  (13 percent would for randomly oriented ribbons), and strands of 1 to 2 um are mostly steep. Many are probably the
  axial blur of SYP puncta or aggregates rather than ribbons; the steep_z flag marks half of all strands. Per-strand
  lengths on real data are not interpretable before the hand traces of protocol S4.
* Flag rates (twelve development images, 35,840 strands, 4,740 traced nuclei): strands steep_z 51 percent,
  end_unfit 39, contaminated 27, exits_label 20, double_width 19, fragment 18, costly_link 16, contact_end 11,
  merged_suspect 8, gap_long 8, neighbour_nucleus 2; strand qc_status review 89 percent, pass 8, fail 3; classes
  partial 69 percent, fragment 18, merged_suspect 7, complete 6. Traced nuclei: count_high 49 percent, count_low 45,
  truncated_z 4, label_merged_suspect 3. Images: psf_fallback 10 of 12, noise_fit_poor none.
* Cost, four workers on a shared machine: 1.3 to 3.5 min per gonad. The main process peaked at 4.5 GB on the
  2025-11-05 stacks and 7.7 GB on a 2026-01-05 stack (106 x 2720 x 2720) at commit 6b7be22; with the plane-by-plane
  reader of commit cf43042 (the same numbers) that stack peaks at 3.9 GB.

Version 3 round (protocol changelog entries 31 to 35; development seeds and development gonads only; the scripts
and outputs are under `C:/Users/ryane/sc_ribbon_work/v3/`). Two rules were pre-specified with their targets before
any number, tuned on the seeds 13000 to 13099, and both withdrawn there; neither was confirmed.

* The punctum class with the image's measured axial sigma (entry 32): it met its yield (82 percent of point puncta
  and 73 percent of aggregates labelled with no strand outside the class; 80 to 81 percent of the puncta of
  `puncta_psf` at every PSF level), its reference (the foci estimate within -2.3 to +3.4 percent of the rendered
  sigma_z at the middle depth; the depth term lowered the error at the puncta's depths from 2.5 to 3.9 to 1.6 to
  3.2 percent) and its excess target (the puncta excess of dense nuclei fell by 33.5 percent), but it labelled two
  axial spurs lying on gaps ribbons and 0.61 and 0.70 percent of the ribbon length of `dense_puncta` and
  `real_dense_r2` (0.5 allowed). On the twelve development gonads (entry 35) it lowered S2a by 0.023 under both
  guards (intervals below 0) and the steep share of the traced length only from 0.36 to 0.35 in SCfull nuclei
  (0.71 to 0.62 in the None zone), so most steep length is not puncta.
* Dim stretches that are not gaps (entry 33): it never read a matched truth gap as continuous (0 of 709 at the pair
  it took), but at every threshold of its grid it joined two truth ribbons 2.75 to 5.5 times per 100 dense nuclei
  (2 allowed) and removed at most 22 percent of the gaps of the brightness phantoms (40 percent required). In dense
  phantoms about one gap per nucleus lies at a junction between two ribbons, where the empty-gap model cannot tell a
  dim stretch from two dim ribbon tips with a short empty core. On the development gonads it changed S2a by
  -0.0001. The gaps that brightness variation creates in dense phantoms (entry 29) remain.
* The hold-out check of entry 35 now covers six development gonads with germquant runs: with the thresholds of the
  twelve, S2a 0.833 over six under a quarter-spacing guard (none evaluable under the registered guard; intensity
  comparator 0.840), S2b 0.962, and SC_fullstart within one point spacing in 3 of 6 (11 of 12 in sample).
* Snapshot of the proposed configuration (version 2.0, seeds 15000 to 15099): S1a met (+0.081 um), S1b not met
  (mean -5.9 percent, SD 8.8 percent, CCC 0.54), S1c, S1d and S1f met, S1e not met (0.889 pooled).

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
