# Design: time-lapse quantification in germquant

Status: design from a 7-agent pass on 2026-10-09 (codebase map, tracker and registration survey, two proposals, two judges, one synthesis). Not implemented. First dataset: the October 2026 GFP::SYP-2 5-minute, 2-hour series on the lab NAS (25 timepoints, 76 or 66 planes, brightfield plus EGFP, no DNA stain). Planned as tranche 2 after the engine groundwork.

## 1. Recommendation

Build Stance A: a timepoint loop around process_image, a per-timepoint shift table, and linking as a file exchange with existing trackers. Graft from B: laptrack (gap closing), sub-voxel registration with a brightfield cross-check and a costed level 2, orthogonal training slices, a realistic segmentation gate, a synthetic time-lapse golden, and holdout versus constraint validation nuclei. About 60 changed lines in germquant plus a read-only timelapse package; nothing appended to any per-image table and no new STAGES entry, so the 23 real goldens stay byte-identical. About eight focused weeks.

## 2. Reader and loop

read_stack gains t (default None): a T axis takes index t when given, else 0 as today. No golden file has a T key, so the fixed branch is untouched; nd2's per-frame chunking reads only that timepoint's 2.25 GB. read_nd2_metadata adds n_timepoints, period_s, t_seconds. process_image gains t; after parse_sample has parsed the raw stem, image_id becomes <stem>_t012 (single underscore) and the body runs unchanged, keyed on image_id. t_index and t_seconds enter run_geometry and the markers only when t is given, so fixed __stages.json and __done.json stay byte-identical too. New timelapse/series.py: process_series with per-timepoint resume and <stem>__series.json; CLI germquant run STACK.nd2 --t all|7|0-24 [--resume]. Self-test: run --t 7 equals the loop's t007.

## 3. Segmentation without DAPI

New channel map: EGFP is central_element, Brightfield an unconsumed role, no dna role. The segment stage reads segmentation.nuclei.channel (already dna everywhere, read by nothing), so frozen recipes resolve the same pixels; the live profile sets it to central_element. The DAPI model will fragment SC bundles, so this is the first real training loop: germquant export-training writes GFP slices and the run's own labels in xy, xz and yz, t24 included; the biologist corrects about 30 slices with a zone tag; train in-process with train_seg as the combined-model script does (150 epochs), split by series, never by slice; register it with sha256. Gate, per zone and timepoint: F1 at IoU 0.5 of 0.85 pachytene, 0.75 transition zone, single-frame miss rate under 10 percent. volume_um3 here is the SC-signal hull, not the envelope.

## 4. Registration

Never resample what the stages see. germquant register runs skimage phase_cross_correlation on stride-4 GFP volumes kept from the loop, upsample_factor 10 (stride 4 alone quantizes to 0.43 um, the size of the flow), previous frame as reference, accumulated to t0. Brightfield gives cross-check columns and a disagreement flag; capped steps are flagged. Quality: median nearest-neighbour displacement before and after the shift. Level 2 only if that stays above 1 um or hand residuals show a coherent field: a Kabsch rigid fit on linked centroids, one register-link-register pass. Elastix stays out of budget.

## 5. Tracking backend

Write no linking code. The contract is on disk: track-export writes germline labels (gonad crop, stride 2, registered) and the GFP crop as a Cell-ACDC Position_1/Images folder plus an export.json; a tracker returns a relabeled npz or a points table. First: Cell-ACDC's overlap tracker headless (imports without Qt). Second: laptrack on registered centroids, gap closing over 2 frames, splitting and merging off; it closes the single-frame gaps the overlap tracker cannot. Importers: tracked npz by majority voxel overlap (splits and deletions logged); points tables from TrackMate, Imaris or the Cell-ACDC GUI. A tracker change never reruns segmentation. btrack, Trackastra (after the cuBLAS guard) and ultrack (pulls napari; never in .venv) wait for evidence.

## 6. Outputs and schema

Per-timepoint outputs are the fixed-image set under the suffixed id; no per-image column changes. New schema.SERIES_TABLES (register, track_points, tracks, track_validation, series_summary): one file per series with SERIES_META, outside TABLES because test_stages forces every TABLES entry through process_image. track_points is the per-object-per-timepoint table: time_track_id (track_id already means an SC fragment), registered centroids, step, link_source, plus nucleus measurements joined on (image_id, nucleus_id). tracks: n_points, n_gaps, is_full_length, speed, signed axis velocity after t0 sign alignment, intensity slopes raw and median-normalised, SC length change. collect stacks them and scores framing QC once per series.

## 7. Biologist validation

Before tracking, pick 16 nuclei per series, stratified and seeded (6 transition zone, 7 pachytene, 3 at the edge), half holdout, half constraint. The biologist follows them in Imaris and exports id, t, z, y, x in raw-frame microns. track-validate looks each point up in that timepoint's labels, maps it to time_track_id, and reports segmentation recall (a point outside any label is a segmentation miss, not a linking error), linking recall, edge precision and switches. Acceptance: zero switches on holdout nuclei, linking recall at or above 0.95, pachytene segmentation recall at or above 0.85. A strips PNG (25 tiles per track) is the non-programmer's check: swaps show as texture jumps.

## 8. First measurements

After phase 1: n_germline_nuclei per timepoint, the bleaching curve, SC length and fragmentation by timepoint and position. After tracking: per-nucleus SYP-2 intensity trajectories raw and median-normalised (the control baseline for the heat-shock series), drift-corrected flow and signed axis velocity by position, within-track SC length change, and the fraction of transition-zone nuclei that depolarise within 2 hours from an off-by-default polarization stage.

## 9. Phases

0 (0.5 wk): local copies, map, profile, channel-key lever, goldens re-diffed, t0 montage. 1 (1): reader, loop, markers, CLI, synthetic golden, 50 timepoints overnight. 2 (2, overlaps 3): export-training, 30 corrected slices, training, model card, rerun. 3 (1): register. 4 (2): export, overlap tracker, laptrack, importers, series tables, strips. 5 (1.5): hand set, track-validate, tuning, first_look.py, docs/TIMELAPSE.md. 6 (1, conditional): level-2 registration or btrack if acceptance fails. Eight weeks, nine with the conditional.

Note added 2026-10-10: the first real series showed a germline that slides inside a still worm non-rigidly, so
phase 6 became a piecewise registration along the gonad axis (one translation per axial segment, interpolated
for every nucleus; `germquant register --mode piecewise`, `track-link --registration piecewise`), ahead of the
hand set's verdict; see docs/TIMELAPSE.md "Piecewise registration".

## 10. Risks

Segmentation without DAPI is the critical path; population numbers per timepoint are the fallback. Packed rows: golden nearest-neighbour spacing is 3.3 to 4.7 um median, 1.2 um at the 10th percentile; flag cutoffs above 0.7 x nn_median_um. Never a full-resolution 4D array (56 GB); cap exports near 2 GB. Qt stays out of .venv; torch backends import after germquant.pipeline. Axis sign can flip between frames; align to t0, test with a deliberate flip. Build the hand set before tuning or it stops being blind.

## 11. Judges and my call

Both chose A, 41 to 36. Only judge 2 flagged stride-4 registration at upsample 1 as a 0.43 um quantization; I take upsample 10 plus the brightfield cross-check. traccuracy: judge 2 says it installs cleanly, judge 1 found it downgrades click in .venv; I keep it out of the numeric venv. Segmentation gate: judge 1 lowered A's 0.9 to 0.85/0.75, judge 2 added the miss-rate gate; I take both, since the SC-trained Koehler model reached only Jaccard 0.78. Both corrected the epoch count (150) and laptrack's field names; adopted.
