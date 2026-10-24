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
germquant sc-view RESULTS_DIR IMAGE_ID [--nucleus N] [--strips]
```

`IMAGE_IDS` are image ids or result-folder names; without them every finished image under the root runs. One
line per image, for example (development gonad n2rad51_noHS_herm_000, six workers) `20251105_N2_nohs_HERM _: 1243
nuclei (traced 663, no_sc 580), 4542 strands (partial 3041, fragment 713, complete 441, merged_suspect 347), flags:
noise_fit_poor, foci_uncalibrated, label_merged_suspect=217, 1.7 min`. Exit code 0 when at least one image
finished or was up to date, 1 when none did, 2 for usage errors.

Cost on the eight 2025-11-05 development gonads (2600 x 2600 x 67 to 71 planes, six workers): 0.7 to 2.2 min per
gonad, at most 5.6 GB in the main process and about 7 GB over all processes; CPU only.

Options: `--config` / `--profile` (only their `sc_straighten` section is read), `--set KEY=VALUE` (for example
`follow.seed_min=3.5`; it enters the parameter hash), `--out DIR` (mirror the tree under DIR; required when the
root is inside a golden tree), `--nuclei germline|pachytene|all`, `--nucleus-ids N ...`, `--workers N`
(spawn-safe processes; the outputs are identical for every worker count), `--force` / `--no-resume`,
`--no-strips`, `--no-overview`, `--source-map OLD=NEW` (the source image moved), `--dry-run`, and for imported
geometry `--geometry auto|import|both`, `--paths-from ims:FILE|swc:FILE|csv:FILE ...`,
`--refine-imported on|off`, `--ends fit|as_drawn`.

Resume: an image whose record `<id>__sc_straighten.json` says `done` with the same sub-hash (version, parameters,
tool versions, the sha256 of nuclei.csv, the labels and the spot tables, the size and first 64 MB of the source
image, the run geometry, the geometry mode and the nucleus selection) is skipped.

## Inputs and refusals

The image id is read from the `image_id` column of `<id>__nuclei.csv` and every other input is found by its exact
name, so ids that end in an underscore or hold spaces work. Inputs: the labels (`<id>__nuclei_labels.tif`, read
page by page into uint16), the run geometry (`<id>__done.json`, else `run_manifest.json` with the image flag
`geometry_from_manifest`), the channel roles as the run resolved them (`channel_map_json`), the spot tables, and
the source image (`file_path`, through `--source-map`). Only the channels it needs are read, one at a time.

An image is refused, with a message, and the others run, when: the source image is missing (`source_missing`,
the message names `--source-map`), the run is at xy stride above 1 (`undersampled_xy`), it is a time-lapse series
folder (`series_folder`), the run has no central_element channel (`no_central_element`), there is no run
geometry (`no_run_geometry`), the label image is missing or has another shape (`no_labels`,
`labels_mismatch`), the output folder is inside a golden tree without `--out` (`golden_folder`), or the image is
a test or reserved image of the sealed SC validation split (`sealed`; tests/golden/sc_split.json; only
`--allow-sealed`, which only `scripts/benchmark_sc.py --final` passes, lifts it, and the record then says
`sealed_override: true`).

## Outputs (all new files, in the image's results folder or the mirrored folder under `--out`)

| file | one row per | notes |
|---|---|---|
| `<id>__sc_strands.csv` | strand | class, QC, length, ends, shape, profile summaries, gaps, width, foci count, strip offset |
| `<id>__sc_strand_samples.csv` | strand x s every 0.1 um | centreline, frame, core, integrated, foci and dna cores, flags per sample; only `image_id` and keys |
| `<id>__sc_strand_gaps.csv` | gap | s interval, length, depth, z score, whether a traced link lies in it, midpoint |
| `<id>__sc_strand_foci.csv` | spot x spot table, in processed nuclei | assignment (unique, ambiguous, off_strand, no_strand), s, folded s, u, v, rank, next-focus distance |
| `<id>__sc_strand_nuclei.csv` | considered nucleus | status (traced, no_sc, skipped, failed), counts, totals, continuity, on-strand fraction, QC |
| `<id>__sc_strand_strips_flat.tif` | | uint16 (role, u, s) strips of every strand concatenated along s; the index is in the TIFF description and in `strip_offset`, `strip_len` |
| `<id>__sc_strand_strips_axial.tif` | | the same along v (the optical-axis direction of the normal plane) |
| `<id>__sc_strand_overview.png` | | SYP maximum projection with the strands coloured by class |
| `<id>__sc_straighten.json` | image | the record: parameters, inputs and their digests, noise model, PSF, counts, flags, timing |

Coordinates are whole-image microns of the run (voxel index times spacing, axis order z, y, x), so the strands,
gaps and foci overlay the run's label image in napari. s runs from the end with the smaller rounded (z, y, x).
The strands, gaps, foci and nuclei tables carry the run's shared metadata block (copied from its nuclei.csv) plus
`sc_straighten_version` and `sc_straighten_hash`. `germquant collect` stacks them into `batch_sc_strands.csv`,
`batch_sc_strand_nuclei.csv`, `batch_sc_strand_gaps.csv`, `batch_sc_strand_foci.csv` and
`batch_sc_straighten_summary.csv` when sc_straighten records exist, and writes nothing new otherwise.

Key strand columns: `length_um` (3D, anchored ends), `length_xy_um` (z ignored), `length_traced_um` (before the
end rule), `e0_um`, `e1_um` (net change at each end against the traced end), `length_se_um`, `end0_class`,
`end1_class` (`free_fit`, `free_unfit`, `contact`, `territory_edge`, `stack_edge`, `cycle_cut`, `max_length`,
`as_drawn`), `snr_core_median`, `integrated_median`, `profile_cv_corr` (noise-corrected heterogeneity),
`n_gaps`, `gap_um`, `n_pieces`, `longest_piece_um`, `fwhm_lat_um`, `fwhm_ax_um` and their PSF-deconvolved widths
(the SC central region is below resolution: widths are QC, not biology), `multiplicity_median` (integrated
signal over the nucleus's median strand signal; the image's complete-strand median for nuclei with under 3 um of
SC), `frac_contaminated`, `n_foci`.

Key nucleus columns: `status`, `n_strands`, `n_main` (strands of 1 um or more), `n_complete`, `expected` (6 for
oocytes, 5 for spermatocytes), `count_status`, `complete_nucleus`, `total_length_um` (every traced nucleus),
`total_length_complete_um` (complete nuclei only), `n_pieces`, `mean_piece_um`, `longest_piece_um`,
`gaps_per_10um`, `gap_fraction`, `piece_lw_mean_um` (the length-weighted mean piece length: the length of the
piece that holds a random point of the nucleus's SC; the continuity score of the validation protocol),
`on_strand_frac` (SYP above background inside the strand tubes over SYP above background in the territory),
`snr_core_median`, `coverage_proxy`, `multiplicity_frac`, `noise_sd`.

## Method in short

Per image: a noise model `var = a + b I` fitted on second differences in the germline labels; a generalised
Anscombe transform to unit noise per nucleus crop; Hessian line fields at about one voxel with a strength
calibrated against unit white noise. Per nucleus: seeded ridge following (predictor-corrector; pure prediction
inside another strand's tube, so crossings pass through), one global end matching (networkx maximum-weight
matching over gaps and crossings), duplicate removal within and across nuclei, normal-plane re-centring on raw
counts and a penalised cubic spline. Per image, a PSF estimate from flat isolated bright stretches (fallback 0.365
/ 0.15 um sigma with the flag `psf_fallback`). Per strand: the anchored erf end rule, rotation-minimising frames
with one roll toward the optical axis, strips, profiles, gaps (the missing integrated signal over a dip, divided
by the median: a gap length that does not depend on the PSF or the tilt), widths, multiplicity, foci, flags and
the class. Details and the reasons for each choice: docs/SC_STRAIGHTEN_DESIGN.md.

## Parameters

A top-level `sc_straighten` section, read only by this command (no shipped config or profile names it, so no
config_hash moves): code defaults < the section of `--config` or `--profile` < `--set`. The keys and defaults
are in `germquant/sc/straighten/params.py` (DEFAULTS). The keys a development sweep may tune are
`follow.seed_min`, `follow.lo`, `link.*`, `refine.knot_um`, `refine.lambda`, `refine.sigma_loc_um`, `gaps.*` and
`classes.*`; everything else is fixed unless a changelog entry of the validation protocol says otherwise.
Development calibration: `scripts/calibrate_sc_straighten.py` (synthetic development seeds only).

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
| label_merged_suspect | nucleus | review | qc.merged_label_ratio | label volume above merged_label_ratio x the median of traced nuclei |
| truncated_z | nucleus | fail | - | the label touches the first or last plane (totals unusable) |
| edge_xy | nucleus | review | - | the label touches the image border |
| count_low | nucleus | info | nucleus.expected | fewer main strands than expected |
| count_high | nucleus | info | nucleus.expected | more main strands than expected |
| noise_mismatch | nucleus | review | noise.nucleus_sd_range | the transformed noise SD of the territory outside nucleus_sd_range |
| budget_exceeded | nucleus | review | budget_s | tracing stopped at the work budget (strands so far kept) |
| error | nucleus | fail | - | an exception in the nucleus (error:<Type>; status failed, no zero row) |
| sc_straighten:uncalibrated | image | info | - | always in v1, until the validation protocol has scored the method |
| noise_fit_poor | image | review | noise.poor_rel_rms | the noise model fits its bins with a relative RMS above poor_rel_rms (or fell back) |
| psf_fallback | image | review | psf.min_stretches | fewer than min_stretches flat isolated bright stretches: the fallback PSF was used |
| sc_channel_suspect | image | review | qc.sc_channel_suspect_frac | among nuclei with SC present, fewer than sc_channel_suspect_frac have a strand of 2 um or more |
| foci_uncalibrated | image | info | - | a spot table was used (spot detection is not calibrated for SC positions) |
| geometry_from_manifest | image | info | - | run geometry (xy_stride, z_range) from run_manifest.json: no done marker |
| label_merged_suspect | image | info | - | label_merged_suspect=N: the number of such nuclei |
| nucleus_errors | image | info | - | nucleus_errors=N: the number of failed nuclei |

## Strand classes (abstain rather than mislead)

`verified` (imported with verified=true), `fragment` (under 1 um), `merged_suspect`, `complete` (both ends
`free_fit`, link cost at most 60, at most 30 percent contaminated, at least 80 percent in its label, no review or
fail flag) and `partial` (everything else). Per-chromosome statistics use only `complete` or `verified` strands.
Nucleus totals use every strand after the accounting step; nuclei without SC keep their weak strands and their
total (a status never zeroes a column).

## napari

`germquant.sc.straighten.layers` returns napari LayerData tuples without importing napari: the strands as a
Shapes path layer (features = the strands table, coloured by class), end points, gap midpoints and foci; a gallery
of one nucleus's straightened strips; an image crop in the run frame. `germquant sc-view` opens them (it needs
napari, which lives in the napari-germquant environment, not in the main one).

## Validation status

Not validated. The pre-registered protocol (docs/SC_VALIDATION_PROTOCOL.md) fixes a development set and a sealed
test set of images, synthetic development and sealed seeds, the metrics and the success criteria; the test set is
scored once, later, by the orchestrator. Development results are kept under
`C:/Users/ryane/sc_ribbon_work/dev_eval/` (summary_<commit>.json); they are not test results and claim nothing.

Development snapshot at the build (commit 3e88403; synthetic development seeds 2000 to 2049, 50 scenes per
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
* End matching works from anchors 0.3 um inside each tip (the PSF-blurred tip drifts along the optical axis),
  with the direction from the next 0.9 um and the gap estimated tip to tip; tips that meet or overlap are judged
  by their collinearity; two `contact` ends facing each other may link across up to 1.5 um (a strand that crosses
  between them); an unlinked end costs 60 (development calibration: curved ribbons broken at sub-micron gaps
  were left unlinked at 40).
* Accounting also removes axial shadows (short steep strands standing on another strand inside the axial PSF;
  before the end matching, so they are never joined to that strand, and again after it) and shadowed strands (at
  least 80 percent of the length inside the contamination reach of a longer strand).
* A trace that reaches the territory or stack edge after its ridge has already faded ends `free` (trimmed to its
  last strong point, so the end rule applies); before, a bright steep ribbon's axial PSF tail could carry the trace
  1 um or more past the true end.
* The tracer follows by pure prediction while inside another strand's claim tube (crossings pass through).
* Gap length is the missing integrated signal over the dip divided by the median (blur-independent), detected
  below 0.6 of the median; the design's 0.35 threshold missed most sub-micron gaps of tilted ribbons.
* The P-spline uses knots every 0.5 um and lambda 0.1 (development calibration: a strong second-difference
  penalty bends a clamped spline outward at its ends and adds length); imported paths measured as drawn use a
  light penalty (`refine.lambda_imported` 0.01).
* An end fit that would extend the strand more than `ends.max_extension_um` (0.5 um) beyond its traced end is
  `free_unfit` (no long straight extrapolation).
* Per-nucleus timings are in the record (`timing`), not in the nuclei table, so the tables are identical for
  every worker count.
* The multiplicity reference is the nucleus's own median strand signal (the image's complete-strand median when
  the nucleus holds under 3 um of SC): SC brightness changes along the gonad, and one image-wide reference called
  most strands of the bright proximal pachytene region merged on n2rad51_noHS_herm_000.
* Not built in v1 (DESIGN 1.3): the napari review dock and edit replay, time-lapse series, the registry stage and
  its profile, deconvolution, the per-gonad zone estimator as a production output, api.py entry points; and the
  real regression goldens of two development tags.
