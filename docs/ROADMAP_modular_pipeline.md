# Roadmap: one modular meiosis pipeline

Written 2026-09-10 from Ryan's statement of purpose: germquant should quantify whatever a C. elegans
meiosis question needs (RAD-51 foci, COSA-1 foci, nuclei and germline, pachytene staging, PGL-1 / SYP-3
partitioning, SC tracking for fragmentation) as one pipeline where you switch on the components you
need, instead of installing a separate pipeline per question.

Status (2026-09-10 evening): steps 0 to 12 are implemented on branch feat/modular-stages (local, not
pushed): golden tooling with synthetic 3- and 4-channel goldens from the frozen commit 8666f07, stage
registry and generated CLI switches, granule stage split out of coloc, per-stage provenance and the
completion marker, hash-preserving profiles, batch resume / exclusions / collect, the SC tracer restored
as sc_trace, envelope + audit + acquisition stages, granule.method imaris_tophat, the staging stage with
`germquant trace` / `restage` and the 25 migrated traces, partition + granule_tail stages and the
ccw77_partition profile, spot instances for COSA-1, ImageJ voxel tags on every label TIF. Real goldens
(7 N2 dry-ice gonads plus stride-4 variants) are generated from a worktree of 8666f07 under
C:/Users/ryane/germquant_golden/8666f07/; the RAD-51 cross-validation gonads join them once the NAS is
reachable, and the ccw77 reconciliation of the partition stage needs the nd2 files on the E: drive.
Real-data checks so far (2026-09-10 night): the current head run at full resolution on the dry-ice
gonad noHS_male_001 matches the frozen golden cell-for-cell in every table, label and mask TIF, with
the same config hash (only appended NA columns); the ccw77_partition profile ran every stage on that
gonad without a failure (envelope, audit, staging from its migrated trace, tophat granules, partition
per zone and pooled pachytene, lit fraction; about 55 minutes with the GPU shared). Two adversarial
reviews of tranches A and B were applied (commits cf0cf61 and 22c9d65; the second caught a label-TIF
writer that produced empty files). Still open: step 13 (Imaris SC calibration, needs traces from
Ryan), the rest of step 14 (image build and A100 diff need Docker and Alpine), step 15 (the deliberate
hash bump, Ryan's call; the truthful `spots:no_foci_role` QC wording waits for it too, because it
changes flag text), the ccw77 reconciliation of the partition stage (needs the E: drive nd2 files) and
the RAD-51 goldens (need the NAS).

## Where things stand today (HEAD 8666f07)

* The pipeline runs read, segment (Cellpose), measure, germline isolation, axis, spots (SpotMAX RAD-51),
  granule + coloc (PGL-1 x SYP). Spots and coloc are the only stages with on/off switches.
* The per-nucleus SC tracer (Sato ridge, hysteresis mask, skeleton; SC length, fragment lower bound,
  fragmentation index) was deleted from the fork on 2026-06-17. It survives byte-identical in the frozen
  C_elegans_ml checkout and in this repo's upstream branches (imaris-validation, segmentation-model,
  zone-calibration), with its 4 tests, its synthetic validator and an Imaris comparison script.
  It was validated on synthetic data only (length recoverable, fragment count not recoverable at
  pachytene density). It was never calibrated against Imaris.
* The July "sc_ribbon" mask inside coloc is a control operand, not a tracer. With the corrected channel
  map (N2 dry-ice run, 2026-08-22, 7 gonads) it still touches 26 to 71 percent of granules, so it leaks
  into the perinuclear shell. It costs 5 to 15 minutes per image and nothing downstream reads it.
* The whole August ccw77 analysis (lamin envelope masking, Imaris-calibrated granule recipe, partition
  coefficient with rotation null and z-shift floor, hand pachytene tracing and zones, mask audit, lit
  fraction, acquisition metadata) lives as standalone scripts under analysis/coloc with hard-coded
  paths, a fixed channel order and a 58 GB crop cache. None of it is a pipeline stage.

## Target: the stage list

Order is fixed and explicit (one tuple in pipeline.py). Each stage has one config key, a role check
against the channel map, a try/except that writes a flag and continues, and writes its own schema
table every run (empty when off). Default state is what keeps today's outputs identical.

| stage | from | key | default |
|---|---|---|---|
| read | existing | io.channel_map | on |
| acquisition (exposure, laser power per role) | analysis/coloc build_acquisition_metadata.py | acquisition.read_exposures | off |
| segment (Cellpose nuclei) | existing | segmentation.nuclei.method | on |
| measure | existing | none | on |
| germline | existing | germline.enabled | on |
| envelope (LMN-1 watershed masks, ring test, territories) | pc_lamin_worker.lamin_nuclei, nucleus_filter.py | envelope.enabled | off |
| axis | existing, gains an enable gate | axis.enabled | on |
| staging (zones from a saved hand trace) | trace_pachytene.py geometry, rezone_all.py | staging.enabled | off |
| spots, instance 1: RAD-51 | existing SpotMAX stage | spots.enabled | on |
| spots, instance 2: COSA-1 | same detector, role crossover_foci, table spots_crossover, optional late-pachytene restriction, expected count QC (6 proposed) | spots.instances | off |
| sc_trace (SC length, fragment lower bound, fragmentation index, skeleton TIF) | frozen sc/skeleton.py restored verbatim | sc.trace.enabled | off |
| granule (PGL-1 segmentation; methods threshold_triangle or imaris_tophat) | existing code split out of coloc | granule.enabled | on |
| coloc (shell voxel Pearson / Manders; keeps the sc_ribbon control under the old key) | existing | coloc.enabled | on |
| partition (distance-matched PC, rotation null, z-shift floor, per zone and pooled pachytene) | pc_lamin_worker.py, pc_zone_worker.py | partition.enabled | off |
| granule_tail (per-granule SYP-3 excess, lit fraction) | granule_tail.py | granule_tail.enabled | off |
| audit (lamin-only candidates, ring scores, overlay) | qc_mask_audit.py | audit.enabled | off |
| qc, render, write | existing, extended | existing keys | on |

CLI additions: one generated `--no-<stage>` switch per stage (`--no-spots` and `--no-coloc` kept as
aliases), `--profile NAME` for small overlay configs in config/profiles/ (rad51_foci, sc_fragmentation,
ccw77_partition, n2_sc), `germquant trace` (the polyline tracer as a separate interactive step whose
saved trace the headless run consumes), `germquant restage`, `--resume` keyed on a completion marker.

## What must not move

* RAD-51 counts on the SpotMAX cross-validation gonads must be exactly reproduced after every step.
  Step 0 freezes goldens (14 gonads: 6 N2 noHS, 6 N2 HS, HS_male_008, one clean hermaphrodite) and a
  regression script diffs every table.
* The config hash is a hash of the raw config text, so config.yaml, config_ccw77.yaml and
  config_n2dryice.yaml are not edited during the migration (not even comments). Every new key is read
  with a code default. The one deliberate hash bump (fold the new blocks into config.yaml, rename
  sc.enabled to coloc.sc_ribbon_operand, drop dead keys) is a separate release Ryan schedules.
* Profiles must preserve the hash of a plain config file: the base_dir logic in config.py only works
  when the file's parent directory is literally "config", so profile loading needs its own path fix
  (critic finding).
* The ccw77 partition values (male pooled pachytene 1.088 v 1.300) are the reconciliation target for the
  partition stage: same crop frame (germline bbox padded 30 voxels), background as the 3rd percentile
  in that crop, region gates 50/500 and per-bin minimums 20/50 for the whole gonad, region gates 30/200
  and per-bin minimums 15/40 for zones (the two workers differ; both get their own keys).
* Legacy traces are in crop pixels times the constant 0.1083 um, while the nuclei tables carry
  0.108333 um; the trace migration must convert with the constant, then map through the crop offset.
* granule_tail v4 drops no-envelope labels before rebuilding the cytoplasm mask and re-segments
  granules inside it, so the stage cannot simply consume the granule stage's labels if it is to
  reproduce granule_tail_v4.csv; it re-segments within the filtered cytoplasm like the script does.

## Order of work

Tranche A, zero numeric change, about 6 working days:
0. Goldens and regression script (1 day plus GPU time). Run one gonad twice to measure Cellpose
   run-to-run tolerance first.
1. Stage tuple, `_run_stage` helper, generated CLI switches, granule split out of coloc (1.5 days).
2. Schema additions, append-only (3 hours).
3. Per-stage provenance in run_manifest.json (stage sub-hashes, model sha256, extras versions,
   xy_stride, z_range), completion marker (1 day).
4. Hash-preserving profiles and `--profile` (1 day).
5. Batch hygiene: `--resume`, `--force`, exclusions file, collect step (1 day).

Tranche B, components, about 11 working days:
6. Restore the SC tracer as sc_trace (1.5 days). The frozen trace_sc traces every label it is given,
   so germline restriction means passing a masked label array. Reproduce the synthetic result on the
   current venv before wiring (length 34.7 v 35.2 um before mask smoothing, about 7 percent shorter
   after; fragmentation index d 0.57). Rerun the balanced-8 N2 batch as the first real-data check.
7. Envelope, audit and acquisition stages (2 days). Also port the SYP-3 nuclear/cytoplasm contrast
   from build_session_covariates.py, which figures 3 and 5 depend on (critic finding).
8. granule.method imaris_tophat and region cytoplasm_shell (1 day). Number-changing only for runs
   that select it. Decide the regression target (Imaris 3186 on HS_male_008; the current workers give
   3220 on the envelope set, 3190 on the DAPI set, 3177 after dropping no-envelope labels).
9. Staging stage, `germquant trace`, `germquant restage`, legacy trace migration (2 days). No
   automated stager is restored: every one plateaued near 5 rows of error.
10. Partition and granule_tail stages with a reconciliation report against the published CSVs
    (2.5 days).
11. COSA-1 as a second spots instance (1 day). Unvalidated until a dataset with a COSA-1 channel
    exists; stays out of every profile until then.
12. Imaris-loadable outputs (spacing tags on every label TIF, SC skeleton TIF), truthful QC flags,
    docs (1 day). Keep sc_ribbon_length exported: test_sc_surface.py uses it.

Tranche C, calibration and HPC, about 3.5 days plus Ryan's tracing and Alpine access:
13. Imaris Filament calibration of the SC tracer: reader for Filament Statistics .xlsx, nucleus
    matching by centroid within 3 um, CCC / bias / Bland-Altman for SC length, lower-bound coverage
    for fragment count. Acceptance before absolute lengths are quoted: CCC at or above 0.6 and bias
    within 10 percent, recorded in the config as a calibration record, never as a length multiplier.
    Until then image_summary carries sc:uncalibrated.
14. Container pins from the validated desktop venv, apptainer.def with the sc extra, model path
    override, Alpine profile; golden diff on an A100 against the tolerance from step 0.
15. Optional hash-bump release (config.yaml text change), only when scheduled.

## What stays out

* The sc_ribbon control operand is not deleted or refactored in this plan (it is under the frozen
  config text); it is documented as a control, switched off in every new profile, and renamed in the
  hash-bump release.
* Automated pachytene stagers (upstream zone caller, refine_pachytene, crescent_axis.py, axis_v2.py).
* The blob_log foci detector (superseded by SpotMAX).
* A dependency resolver that could reorder stages; partial reload of intermediates.
* Figure, grant and one-off robustness scripts (pubstyle, fig_pub_*, fig_grant_*, amount_metrics,
  denominator_check, exposure_sim, tail_by_distance bleed correction): analysis layer, not pipeline.
* The crop cache, chload.py fixed channel indices and the hard-coded run/NAS path logic.

## Decisions only Ryan can make

1. SC ground truth: trace SCs with the Imaris Filament Tracer in 15 to 25 pachytene nuclei per
   condition on one N2 control and one heat-shocked gonad from the RAD-51 validation set, and export
   Filament Statistics as .xlsx. Without this the SC stage ships uncalibrated.
2. COSA-1: is there a dataset with a COSA-1 channel, which fluorophore, and is 6 per late-pachytene
   nucleus the right expected count for the QC?
3. Which granule count the imaris_tophat port must reproduce on HS_male_008.
4. Keep config_ccw77.yaml frozen with the August recipe only in a profile (recommended), and whether
   to re-run the July 0708/0622 folders under that profile as the new source of record.
5. When to take the config hash bump.
6. Whether hand tracing with `germquant trace` becomes the lab's standard staging step (it is what
   makes pachytene SC means, the COSA-1 late-pachytene restriction and partition zones available).
7. Where the golden outputs live (E: drive vs the lab NAS), and Alpine account/partition strings.

**Tree-wide golden diff (2026-09-24):** all 11 real goldens (7 full-resolution N2 dry-ice gonads and 4 stride-4 variants, tests/golden/real_manifest.json) rerun with the current code (65709aa through 12a5f75, i.e. main after PR #3 plus PR #4) match the frozen 8666f07 goldens exactly: 0 failing comparisons, only the appended columns differ. The same gonad run under Python 3.14.5 with the same pins is also byte-identical (see the Python note below).

**Python version (2026-09-24):** Python 3.14 validated but not yet adopted: every package in the environment has 3.14 wheels for Windows, Mac and Linux, the README recipe installs unchanged, the test suite is identical, and the golden gonad diff is clean. Switching means bumping 3.11 to 3.14 in the README, docs/MAC_METAL.md, the Dockerfile and both constraints files (pinning the full freeze), and regenerating or removing pixi.toml and environment.yml, which still describe an older environment.
