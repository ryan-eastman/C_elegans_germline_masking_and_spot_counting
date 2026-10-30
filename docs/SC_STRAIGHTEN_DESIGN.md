# SC ribbon straightening and quantification: synthesized design (v1)

Written 2026-10-23 by the synthesis step of the SC ribbon workflow. Status: design only. Nothing in the engine
repo, its worktrees, the goldens or the NAS was changed. Base for the build: origin/main 8eb73de.

Inputs: Design A (accuracy of length and continuity), Design B (robustness and usability), the code brief
(`code_brief.md`), the field review (`field_review.md`) and the truth inventory (`truth_inventory.md`), all in
`C:/Users/ryane/sc_ribbon_work/`. Companion files written in this step:

- `VALIDATION_PROTOCOL.md`: the pre-registered validation protocol, written before any SC output was compared
  with any SC truth.
- `sc_split.json` (sha256 db328c1c...e9eb): the registered split, built by `scripts/make_sc_split.py` from the
  inventory metadata and the bytes of the truth copies (no truth value was read).

Everything below is ASCII so the file can be copied to `docs/SC_STRAIGHTEN_DESIGN.md` unchanged (rule R9).

## 0. Verdict

### 0.1 Judgement

| criterion | Design A | Design B | taken |
|---|---|---|---|
| correctness of the geometry core | seeded predictor-corrector ridge following plus one global end matching; measured on exact-render synthetic scenes: dense nuclei total -0.1 +/- 4.2 %, 6.7 strands for 6 ribbons | non-maximum suppression, hysteresis, skeleton graph, junction pairing; no accuracy measurement, and its own probe found that skeleton length follows the threshold; A's crest-graph runs of the same family gave +17 % or -6 % totals | A |
| correctness of length and ends | re-centred P-spline plus anchored erf half-maximum ends: +0.05 +/- 0.16 um on isolated ribbons at twice the real SNR, +0.28 +/- 0.18 um at real SNR before tuning | smoothing spline with one isotropic residual target, relative half-max end with at most 0.5 um extension; unmeasured | A |
| robustness of the shell | refusals, label QC, a few flags | per-nucleus error isolation (never a zero row), work budget, refusals with messages, source map, flag vocabulary as data, memory plan, lazy channel read | B |
| fit with the existing code | post-run command; nothing in STAGES, schema.TABLES or process_image; but it reuses `results_images` (wrong ids, 0.3 item 1), adds a profile that needs a new pin, and its file sizes do not fit (0.3 item 4) | post-run command wired like cli_packs, collect branch, no profile; but it adds to api.py's stability surface and puts series code next to an area another workflow is changing | A's placement with B's wiring and corrections |
| testability | numeric tolerances per synthetic scene, sha256 fence over every pre-existing file, tiny h5py .ims fixture | one positive and one negative phantom per QC code, determinism across worker counts, import fence | union of both |
| value to the lab | per-ribbon length, gaps and foci where the decomposition is clean; a total for every nucleus; strips | review loop with hand-trace mode, heat-shock readout (on-strand SYP fraction), overview picture, time-lapse | A's numbers and B's readouts now; B's review loop next |

### 0.2 Decision

- **Core from A.** Per-image noise model and variance-stabilising transform; Hessian line fields at one voxel with
  a strength calibrated against a noise null; seeded predictor-corrector following; one global end matching
  across crossings and gaps; length accounting within and across nuclei; normal-plane re-centring on raw counts;
  penalised cubic spline with per-axis weights; anchored half-maximum erf end rule; rotation-minimising frames with
  one constant roll toward the optical axis; strips; profiles; gaps; foci; abstaining strand classes; a nucleus
  tier that does not depend on a perfect decomposition.
- **Grafts from B.** The QC flag vocabulary as data at three levels with a test per code; nucleus status (traced,
  no_sc, skipped, failed) and error isolation; a cooperative per-nucleus work budget; refusals with messages and
  `--source-map`; a channel-subset reader and uint16 labels (memory); one slim sample table; the on-strand SYP
  fraction (an intensity partition that tracing choices barely move, the heat-shock readout); a per-nucleus noise
  check; the `sc_channel_suspect` image flag; the overview PNG; determinism across workers; the collect fence and
  the import fence; a docs flag table kept equal to the code by a test.
- **New in the synthesis.** Image-id discovery that works for ids ending in an underscore; a sealed-split guard and
  a golden-folder guard; a concatenated (ragged) strip layout; the continuity score `piece_lw_mean_um` that the
  protocol scores against zone truth; a benchmark harness in the layout of the RAD-51 harness
  (`germquant.benchmark` plus `scripts/`), with guarded truth loaders.
- **Deferred** (section 1.3): the napari review dock and edits replay (B), time-lapse series and strand linking
  (B), the registry stage and its profile (A, phase 2), deconvolution (B), stride above 1, the production zone
  estimator, api.py entry points.

### 0.3 Corrections found by checking both designs against the code and the data

1. **Image ids that end in an underscore.** `staging.tracer.results_images` and `batch.collect` take
   `name.split("__")[0]`. For `20251105_N2_nohs_HERM ___nuclei.csv` that gives `20251105_N2_nohs_HERM ` (the
   trailing underscore is lost), and `provenance.read_done_marker` then misses `20251105_N2_nohs_HERM ___done.json`
   (checked on the bc62476 golden: found with the true id, not with the parsed one). This hits all four
   replicate-000 gonads of 2025-11-05, which are four of the eight development gonads (herm_000 among them). Design
   A's plan to reuse `results_images` and the done marker would mislabel them and lose their run geometry; Design
   B's statement that every name has exactly one `__` is false for them. v1 reads the id from the table's
   `image_id` column and matches exact suffixes. The same bug makes `germquant collect` drop these images from
   every batch table and from `batch_summary.csv`, with only a logged warning (checked on scratch copies of the
   bc62476 herm_000 and herm_001 tables: herm_001's 617 nuclei were stacked, herm_000's 1274 rows vanished);
   fixing it changes batch outputs, so it goes to Ryan as a separate item (section 13.3) and is not part of this PR.
2. **Input folders.** The 8666f07 goldens have no done marker (only `run_manifest.json`). The bc62476 goldens
   (the T5 baseline: `regression_diff.py` found 0 failing comparisons against 8666f07 on all 23 tags) have done
   markers with `xy_stride` and `z_range`, empty `zones.csv` and empty SC tables. Development reads the bc62476
   folders; `run_manifest.json` is the fallback for run geometry (image flag `geometry_from_manifest`).
3. **Memory.** A 2026-01-05 stack is 106 x 2720 x 2720 x 3 channels: 4.70 GB as `read_stack` returns it, and its
   int32 labels are 3.14 GB. Design A's parent process (all channels plus int32 labels, "peak under 6 GB") would
   exceed the 8 GB per-process budget. v1 reads only the channel roles it needs and reads labels page by page
   into uint16 when the largest label fits (B's two ideas): 1.57 GB per channel plus 1.57 GB of labels there.
4. **Output size.** Design A's gallery, float32 (strand, role, u, s) padded with NaN to the longest strand, is
   about 0.5 GB per gonad for 3000 strands, and its two 0.05 um sample tables carry SHARED_META on every row
   (a few hundred MB each). v1 writes strips concatenated along s in uint16 with an offset index, and one sample
   table every 0.1 um that carries only `image_id` plus keys (B's documented exception).
5. **Evidence gaps in A.** The skew criterion of the end matching appears in no probe (no skew parameter in any
   `logs/designA/track/*.json`); every dense variant was compared on the same 12 synthetic nuclei (seed 11); the
   dense total of -0.1 +/- 4.2 % combines 97.6 % coverage with 3.2 um per nucleus lying off truth, so it is partly
   a cancellation; the gap, foci and `complete` rules are unmeasured; isolated ribbons at real SNR are +0.28 um
   before tuning. Hence sealed synthetic seeds and separate coverage and off-truth gates (protocol S1).
6. **B's tracing core** is the family that both probes found unreliable at pachytene density: B's own probe
   ("any mask-skeleton length is an artefact of the threshold") and A's crest-graph runs, against 0.0 +/- 4.2 %
   for following plus matching.
7. **B's analysis rules.** `in_analysis` admits strands that passed QC but were never reviewed while the pass
   precision is unknown, and totals only for complete nuclei would leave most nuclei without a number (14 % of
   synthetic ribbons come out whole and pure). v1 keeps A's abstaining classes and reports a total for every
   traced nucleus, with coverage and SNR beside it.
8. **Profiles.** A's `n2_sc_straighten.yaml` would add a config_hash pin
   (`test_every_shipped_config_and_profile_is_pinned`) for no v1 benefit. v1 ships no profile.
9. **Truth.** Neither design's validation inputs exist as assumed. Madeleine's "SC lengths" are SC zone
   boundaries along the gonad axis (39 N2 gonads); per-ribbon traces exist for 2 live nuclei only. The protocol
   is built on zone truth, synthetic truth and the live check, and asks for hand traces.
10. **Madeleine's summary code** (`20260121_SC_zone_analysis.R`, read as code, not run): the segment labels start
    each zone at its named point, but the summary takes each zone's START as the cumulative length at the END of
    the zone's first segment (`min(cum_length)` over the zone), one measurement point after the named point (about
    20 um). Zone ends are at the named points. If that is not intended, the derived start positions and the SCfull
    lengths in the lab's tables are shifted by one point spacing. The protocol scores against the named points;
    the question goes to Ryan and Madeleine (section 13.3).

## 1. Scope

### 1.1 v1 (built, tested and run on development data in this workflow)

1. `germquant sc-straighten`: a post-run command on finished fixed-image results at xy stride 1. It traces every
   SC ribbon in 3D, straightens it, measures it, places foci on it and classes it, and writes only new files.
2. Imported geometry (Imaris Filaments in .ims, SWC, a CSV of paths) measured through the same path, with
   refinement and the end rule switchable. Needed for the live check (protocol S3) and the hand traces (S4).
3. `layers.py` (napari layer data without importing napari) and `germquant sc-view`, which imports napari lazily.
4. `germquant collect` stacks the new tables when they exist and writes nothing new when they do not.
5. The benchmark harness (`germquant.benchmark.sc_ribbons`, `germquant.benchmark.sc_split`,
   `scripts/benchmark_sc.py`) with guarded truth loaders, the split record and the synthetic phantoms.
6. Tests and synthetic goldens; development runs on the eight development gonads of 2025-11-05; development
   calibration on development seeds; development zone scoring after the protocol is registered.
7. Should, if time allows (none blocks the PR): the overview PNG; real regression goldens for two development
   tags; development runs on the ten other development gonads (they need germquant runs on the GPU first).

### 1.2 Not in v1 and never in this workflow

- Any sc_straighten output on a test or reserved image, any read of test or reserved truth, any sealed synthetic
  seed. The orchestrator scores the test set once, later (protocol section 10).

### 1.3 Deferred, with the trigger that starts each item

| item | from | trigger |
|---|---|---|
| napari "SC ribbons" dock: review queue, accept/reject, re-trace (snap plus minimum-cost path), split, merge, ends, blind truth mode, append-only edits file with replay by geometry matching | B | engine PR merged; plugin PR in napari-germquant, after Ryan's go |
| time-lapse series (per timepoint in place), strand linking within tracked nuclei, kymographs, labels from stride-k runs | B | the linker branch has landed; the DLW93 blind rule is lifted |
| registry stage `sc_straighten` after `sc_trace` (hash_when_disabled, run_stage only when enabled, tables outside TABLES) and its profile with a pin | A | the protocol's test scoring has passed |
| P5b post-stage adapter (`<id>__post_<name>_<table>.csv` naming) | B | P5b lands; the names live in one constants module, so a rename is one line plus a golden refresh |
| Richardson-Lucy deconvolution for tracing only | B | a development experiment after v1 |
| per-gonad zone estimator as a production output (automated zone lengths) | new | protocol S2c passes |
| api.py entry point | B | the plugin PR needs it |
| parquet copies, SNT .traces reader (SWC covers SNT exports), 3D strip files | A, B | on request |

## 2. Placement and fences

### 2.1 New files

```
src/germquant/sc_straighten/
  __init__.py   VERSION = "1.0"; run_results(), straighten_image() (experimental public functions)
  params.py     the resolved sc_straighten section: code defaults < --config/--profile section < --set; params_hash
  inputs.py     discovery, image ids, run geometry, channel-subset reader, labels as uint16, spot tables,
                refusals, golden-folder guard, sealed-split guard
  noise.py      noise model fit, generalised Anscombe transform (GAT), per-nucleus noise check
  crops.py      per-nucleus crops, territory, label QC
  fields.py     Hessian line fields, noise-null calibration (cached)
  follow.py     seeded ridge following (A's track.py; "track" is taken twice in this repo)
  link.py       global end matching
  account.py    duplicate removal within and across nuclei
  refine.py     re-centring, P-spline, arc length, PSF estimate, end rule
  frames.py     rotation-minimising frames, roll, inverse map (to_strip, from_strip)
  strips.py     cross-section sampling, flat and axial strips, concatenated layout
  measure.py    profiles, gaps, widths, multiplicity, strand measures, classes, nucleus tier
  foci.py       spot assignment to strands
  qc.py         the flag vocabulary as data (code, scope, level, rule, threshold key, description)
  imports.py    Imaris Filaments (.ims via h5py), SWC, CSV paths
  schema.py     declared columns of the new tables (outside germquant.schema.TABLES)
  outputs.py    atomic writers, record, sub_hash, resume, overview PNG, collect_tables()
  layers.py     napari LayerData tuples (no napari import)
  synth.py      exact-render phantoms with truth polylines, gaps and foci
src/germquant/cli_sc_straighten.py    add_parsers(sub), dispatch(args): sc-straighten, sc-view
src/germquant/benchmark/sc_ribbons.py zone truth, filament truth, synthetic scoring, metrics, boundary estimator
src/germquant/benchmark/sc_split.py   split record loader and seal guard
scripts/benchmark_sc.py               manifest | synthetic | zones | filaments | report  (--final is logged)
scripts/calibrate_sc_straighten.py    development sweeps (SNR x tilt x density), development seeds only
scripts/make_sc_straighten_goldens.py real regression goldens for development tags
tests/golden/sc_split.json            byte copy of C:/Users/ryane/sc_ribbon_work/sc_split.json
tests/golden/sc_straighten/           synthetic goldens
tests/test_sc_straighten_*.py, tests/test_benchmark_sc.py, tests/test_sc_split.py
docs/SC_STRAIGHTEN.md                 usage, outputs, columns, flag table (kept equal to qc.CODES by a test)
docs/SC_VALIDATION_PROTOCOL.md        copy of VALIDATION_PROTOCOL.md (from registration on, append-only changelog)
```

### 2.2 Edits to existing files (all additive and inert by default)

- `cli.py`: `from . import cli_sc_straighten`, `cli_sc_straighten.add_parsers(sub)` next to the cli_packs line, and
  two dispatch entries (`"sc-straighten"`, `"sc-view"`). No existing command changes.
- `batch.py`, end of `collect`: `n = sc_straighten.outputs.collect_tables(out_root)`; when it returns None (no
  sc_straighten record under the root) nothing is written and `counts` is unchanged.
- `docs/SC_TRACING.md`: one paragraph pointing to SC_STRAIGHTEN.md and stating the code brief's finding that
  `sc_total_length_um` follows nucleus volume and background on real data.

Untouched: `process_image`, `STAGES`, `stage_hashes`, `_STAGE_TOOLS`, `schema.TABLES`, `sc/skeleton.py`,
`sc/surface.py`, `sc/__init__.py`, `io/*`, `staging/*`, `timelapse/*`, `config/config.yaml`, every profile and
channel map, `api.py`, `default_invariants.json`. The fence therefore holds by construction:
`test_default_invariants.py`, `test_stages.py`, `test_golden.py` and `regression_diff.py` stay green unedited.

### 2.3 Commands

```
germquant sc-straighten RESULTS_ROOT [IMAGE_IDS ...]
    [--config C | --profile P]      only their sc_straighten section is read; default: code defaults
    [--set KEY=VALUE ...]           e.g. --set follow.seed_min=3.5 (enters params_hash)
    [--out DIR]                     mirror the tree under DIR; required when the root is inside a golden tree
    [--geometry auto|import|both] [--paths-from ims:FILE|swc:FILE|csv:FILE ...]
    [--refine-imported on|off] [--ends fit|as_drawn]
    [--nuclei germline|pachytene|all] [--nucleus-ids N ...]
    [--workers N] [--force] [--no-resume] [--no-strips] [--no-overview]
    [--source-map OLD=NEW ...] [--allow-sealed] [--dry-run]
germquant sc-view RESULTS_DIR IMAGE_ID [--nucleus N] [--strips] [--stride K]
```

`IMAGE_IDS` are image ids or result-folder names (for example `n2rad51_noHS_herm_000`). One line per image, for
example
`20251105_N2_nohs_HERM _: 1243 germline nuclei (traced 498, no_sc 742, failed 3), 4105 strands (complete 211,
partial 2810, merged_suspect 96, fragment 988), flags: psf_fallback, 9.6 min`.
Exit code 0 when at least one image finished, 1 when none did, 2 for usage errors.

### 2.4 Guards

- **Golden folders.** The command refuses to write into any folder that has a `golden_manifest.json` in it or in a
  parent; `--out` is then required.
- **Sealed split.** When `tests/golden/sc_split.json` is found next to the package (`provenance._repo_root()`) and
  says `"sealed": true`, the command refuses an image whose source stem equals the `image_stem` of a test or
  reserved gonad, unless `--allow-sealed` is given; the refused image is skipped with a message and the others
  run. With the flag the record carries `sealed_override: true`; only the scoring script passes it, under
  `--final`. Outside the repo (installed elsewhere) the guard is silent. Development runs name their images
  explicitly, because a golden root such as bc62476 also holds the four test tags and the dry-ice tags.
- **Truth loaders** in `benchmark.sc_split` and `benchmark.sc_ribbons` raise `SealedError` for test or reserved
  gonads unless called with `final=True`.

## 3. Inputs, discovery and refusals

1. **Discovery.** Walk `RESULTS_ROOT` for files ending in `__nuclei.csv` (through `fsutil.long_path`). The image id
   is the `image_id` value of the table's first row; the file name must equal `image_id + "__nuclei.csv"`
   (otherwise warn and use the exact-suffix strip). Every other input is looked up by exact name
   `f"{image_id}__<suffix>"` in the same folder. Series folders (`<stem>__series.json` or `t<NNN>/` children) are
   refused in v1 with a message.
2. **Run geometry.** `xy_stride` and `z_range` from `<id>__done.json`; else from `run_manifest.json` in the folder
   (flag `geometry_from_manifest`); else refuse (`no_run_geometry`). `xy_stride` above 1 is refused
   (`undersampled_xy`); live runs come with the time-lapse work.
3. **Channel roles** from the run's `channel_map_json` column (SHARED_META), not from the current config, so a
   session with swapped channels (2025-10-21, SC on 545) is measured as it was run. No central_element role:
   refuse.
4. **Source image.** `file_path` column, rewritten by `--source-map OLD=NEW` (prefix replacement). Missing: refuse
   with a message that names `--source-map`. The record keeps size, mtime and the sha256 of the first 64 MB.
5. **Channel-subset reader** (`inputs.read_roles`). For .nd2: the same dask construction as `io.read_stack` (T and
   P loops dropped, C and Z ensured, transposed to C, Z, Y, X, z_range and stride applied before compute), then
   one channel computed at a time, so the peak is one channel. Other formats: `io.read_image`, then slice (a
   warning names the memory cost). A workstation-only test checks equality with `read_stack` on herm_000. Roles
   read: central_element always; the profile roles in `strips.roles` and `profile.roles` (default foci).
6. **Labels.** `<id>__nuclei_labels.tif`, read page by page (tifffile) into uint16 when the largest nucleus id is
   at most 65535, else int32. Assert shape equals the stack's and the spacing equals `voxel_d*_um` of nuclei.csv
   within 1e-9; refuse on mismatch.
7. **Nuclei.** `in_germline` rows (default); `pachytene` needs `is_pachytene` from `__zones.csv` or nuclei.csv and
   refuses when it is empty (the 2025-11-05 goldens never ran staging); `all` takes every label.
8. **Spot tables.** `<id>__spots.csv` and every `<id>__spots_<name>.csv` (exact prefix), whole-image um, with
   `nucleus_id`. Their use sets `foci_uncalibrated` (RAD-51 detector at 8666f07: spot F1 about 0.6).
9. **Frame.** Whole-image microns of the run: whole-image voxel index x spacing. A crop with origin `org` maps
   crop index c to (c + org) x spacing. Every coordinate written is in this frame, so napari overlays line up with
   the existing layers.

## 4. Algorithm

### 4.1 Flow per image

```
inputs -> refusals -> sub_hash -> (resume: skip if record done with the same sub_hash)
noise model on germline territory voxels (the GAT is applied per crop, never to the whole channel)
null_q = calibrate_null(sigma_vox)                        (cached per (sigma_vox, q, shape, seed))
per nucleus, in parallel over crops:  crop -> GAT -> fields -> follow -> link -> account (within) -> RawStrands
account across nuclei (in nucleus order)
per strand: re-centre -> P-spline
PSF estimate from flat, isolated, bright stretches (fallback with flag)
per strand: end rule -> frames -> samples and strips -> profiles -> gaps -> widths, multiplicity -> foci -> flags, class
per nucleus: tier aggregates;  per image: flags, counts
write tables, strips, overview atomically; write the record last (done: true)
```

### 4.2 Noise model and transform (`noise.py`)

- Fit on second differences along x (variance 6 sigma^2) inside the germline territory, 12 planes, 30 intensity
  bins up to the 85th percentile, robust MAD per bin, weighted least squares `var = a + b * I`.
- GAT: `y = (2 / b) * sqrt(max(b * I + 3/8 * b^2 + a, 0))`, unit noise. Image flag `noise_fit_poor` when the relative
  RMS of the fit exceeds 0.15. On herm_000 Design A measured `var = 4.04 * I - 394`.
- Per-nucleus check (B): robust SD of the GAT second differences in the territory; nucleus flag `noise_mismatch`
  outside [0.7, 1.4].

### 4.3 Crops, territory, label QC (`crops.py`)

- One `ndi.find_objects` pass; boxes padded by (0.6, 0.5, 0.5) um (z, y, x); crops are copies of under 50 MB.
- Territory: voxels within 0.4 um of the nucleus's own label (EDT with `sampling=spacing`) that lie in no other
  label. Nucleus flags: `label_merged_suspect` (volume above 1.6 x the median volume of traced nuclei),
  `truncated_z` (box touches the first or last plane), `edge_xy` (box touches the image border).

### 4.4 Line fields (`fields.py`)

- Gaussian derivatives of the GAT crop at sigma (0.2, 0.108, 0.108) um (one voxel per axis); `eigh` only on
  territory voxels; eigenvalues w0 <= w1 <= w2, eigenvectors v0, v1, v2.
- Strength `S = max(-w1, 0) * sigma^2 / null_q`; Newton offset `p = -sum_k (g . v_k / w_k) v_k` over k in {0, 1};
  lineness `|w2| / |w1|`; crest voxel: `w1 < 0` and `|p|_inf <= 0.6` voxel.
- `null_q`: the 99.9 % quantile of S on crest-qualified voxels of unit white noise on a 40 x 96 x 96 grid, seed 0;
  recorded (0.38 on herm_000).

### 4.5 Seeded following (`follow.py`)

- Seeds: crest voxels with S >= 4.0 in descending S, ties broken by voxel index; a claimed voxel never seeds. Both
  directions are followed from each seed.
- Step: predict 0.1 um along the current direction t; correct with the Newton offset (normal-plane component,
  capped at 0.15 um); update `t = normalise(0.7 t + 0.3 d_local)`, using d_local only where lineness <= 0.5
  (inside junction cores momentum carries the trace through).
- Stops and end types: S < 1.5 for more than 0.6 um (trim to the last strong point; `free`); leaving the
  territory (`territory_edge`); entering another strand's claim tube in parallel, |cos| > 0.85 (`contact`);
  re-entering its own path more than 0.6 um back (`loop`); 30 um (`max_length`); within 0.4 um of the first or
  last plane (`stack_edge`).
- Claim tube: (1, 2, 2) voxels around accepted points. Strands under 0.5 um are dropped.
- Budget (B): a cooperative check every 200 steps against `budget_s` per nucleus; when exceeded, the strands so
  far are kept and the nucleus gets `budget_exceeded`.

### 4.6 End matching (`link.py`)

- Candidate ends: types `free`, `contact`, `loop`. Outward direction from the last 0.15 to 0.9 um.
- Admissible pair: gap <= 1.0 um, both directions within 50 deg of the link, closest approach of the two end
  tangent lines <= 0.3 um (the skew criterion; unmeasured so far, so its effect is reported by the development
  sweep and it can be switched off with `link.max_skew_um: null`).
- Cost `40 * gap_um + angle_a + angle_b + 100 * skew_um` (degrees, um); an unlinked end costs 40;
  `networkx.max_weight_matching` with weight `80 - cost` on a graph built in sorted end order; a strand's two ends
  never pair with each other; chains become strands and a cycle is cut at its most expensive link
  (`cycle_cut` ends). Each link records gap, angles, skew and cost.

### 4.7 Length accounting (`account.py`)

- Within a nucleus: stretches where two strands run within 0.15 um laterally and 0.3 um axially, parallel
  (|cos| > 0.9), for at least 0.4 um stay on the stronger strand and are cut from the other; remainders under
  0.5 um are dropped. Crossings stay on both strands.
- Across nuclei: a strand belongs to the label that holds most of its length; copies traced from a neighbour's
  territory are removed by the same parallel rule.
- A nucleus total is the sum of its strands after this step.

### 4.8 Re-centring and spline (`refine.py`)

- Resample at 0.1 um; initial P-spline; four rounds of normal-plane mean-shift on raw counts minus the territory
  background (Gaussian window 0.15 um lateral, 0.30 um axial, grid +/-0.4 x +/-0.6 um at 0.05 um, shift capped at
  0.3 um; contaminated samples and junction cores weigh 0).
- Final fit: cubic B-spline, knots every 0.3 um, second-difference penalty `lambda * mean(W)`, weights
  `W = 1 / sigma_loc^2` with sigma_loc (0.06, 0.03, 0.03) um (z, y, x). Knots, lambda and sigma_loc are tuned at
  matched SNR on development seeds (`scripts/calibrate_sc_straighten.py`). Arc length by 5-node Gauss-Legendre
  quadrature per knot interval; `fit_rms_um` is the weighted residual RMS.

### 4.9 PSF estimate (`refine.py`)

- From flat (|t_z| <= 0.3), isolated (no other strand within 1 um), bright (core SNR >= 10) stretches of at least
  1 um: Gaussian fits of the mean lateral and axial cross-sections give sigma_xy and sigma_z. With fewer than 20
  stretches the fallback sigma_z 0.365 um, sigma_xy 0.15 um (FWHM 0.86 and 0.35 um, Design A's measurement on
  232 flat strands of herm_000) is used and the image is flagged `psf_fallback`.

### 4.10 End rule (`refine.py`)

- Anchor min(0.6 um, L/3) inside each `free` end. Profile p(s) from -0.6 to +1.4 um in 0.05 um steps: along the
  spline inside the anchor, straight along the anchor tangent outside it (mean-shift re-centring capped at
  0.2 um); value = background-subtracted mean in an ellipse 0.12 um lateral x 0.25 um axial.
- Fit `p = b + a * Phi((e - s) / sigma)` by least squares from the half-maximum crossing; sigma bounded to 0.6 to
  1.6 x sigma_t, `sigma_t = sqrt(sigma_xy^2 (1 - t_z^2) + sigma_z^2 t_z^2)`.
- End classes: `free_fit` (a >= 5 x SE and residual RMS <= 0.25 a), `free_unfit`, `contact`, `territory_edge`,
  `stack_edge`, `cycle_cut`, `max_length`. Non-free ends stay where the trace stopped.
- `length_um` = spline length between the anchors + e0 + e1; `length_se_um = sqrt(se_e0^2 + se_e1^2 + se_curve^2)`
  with se_curve from the development calibration table (by SNR and tilt).

### 4.11 Frames, sampling and strips (`frames.py`, `strips.py`)

- Rotation-minimising frames by double reflection (Wang 2008) on samples every 0.05 um over s from -0.5 to
  L + 0.5 um; one roll per strand, `phi = atan2(sum w (n2 . z_perp), sum w (n1 . z_perp))` with `w = |z_perp|^2`;
  `n_ax = cos(phi) n1 + sin(phi) n2`, `n_lat = t x n_ax`; when `sum w` is about 0 (a strand along z) n_ax = n1.
  `roll_from_axial_deg` is kept per sample. Never Frenet frames.
- Cross-sections: u in [-0.6, 0.6] um at 0.05 um (25 samples), v in [-1.0, 1.0] um at 0.1 um (21 samples);
  points `c(s) + u n_lat + v n_ax`; `map_coordinates` order 1, cval NaN, on float32 raw counts of each role.
- Flat strip: mean over |v| <= 0.3 um, shape (role, u, s). Axial strip: mean over |u| <= 0.15 um, shape
  (role, v, s). Order 3 is used only for pictures.
- Inverse maps: `to_strip(points_um) -> (strand, s, u, v)` (KD-tree on samples plus two Newton steps) and
  `from_strip(strand, s, u, v) -> points_um`; their round trip is a test.
- Neighbours: samples within 0.35 um lateral or 0.7 um axial of another strand's centreline are `contaminated`.

### 4.12 Imported geometry (`imports.py`)

- `read_ims_filaments(path, time_index=None)`: h5py on `/Scene8/Content/Filaments<k>/` (Vertex PositionX/Y/Z,
  Edge VertexA/B, Filament and DendriteSegment ranges, time ids); unbranched paths from the edge graph, split at
  branch points; single-vertex filaments ignored; the stored `Segment Length` statistics kept for the check S3a.
- Frame: `index = (p - ExtMin) / v_imaris - 0.5` (continuous; Imaris positions are voxel centres), then
  um = index x our voxel; ExtMin, ExtMax and size from `germquant.validate.imaris_ims.read_geometry`. Our voxel
  defaults to the nd2 voxel of the run; for a TIFF without z spacing it is given explicitly (the live check: 0.108,
  0.108, 0.223077 um).
- `read_swc(path, spacing)`; `read_paths_csv(path)` with columns image_id, strand_id, point_index, z_um, y_um,
  x_um, verified.
- Assignment to a nucleus: the label holding most of the path's points, else the nearest territory. Measurement
  as for automatic strands; `--refine-imported off` and `--ends as_drawn` separate tracing error from measurement
  error. `verified=true` gives class `verified`.

### 4.13 Determinism, workers, errors

- Workers (`--workers N`, spawn-safe process pool) receive crops; results are assembled in nucleus order;
  `strand_id` is assigned last by (nucleus_id, -length_um, end coordinates); `strand_uid` is the first 10 hex of
  sha1 over the centreline rounded to 0.01 um. Outputs are identical for every worker count (a test).
- An exception in a nucleus gives status `failed` with `error:<Type>` and the message in the record, never a zero
  row; an exception in an image gives a failed record and the next image runs.

## 5. Measurements, flags and classes

### 5.1 Strand measures

`length_um` (3D, anchored ends), `length_xy_um` (z ignored, for xy-only lab measurements), `length_traced_um`,
`e0_um`, `e1_um`, `length_se_um`, `chord_um`, `straightness`, `curvature_mean_per_um`, `curvature_max_per_um`,
`tz_abs_mean`, `frac_steep45` (|t_z| > 0.707), `frac_steep60` (|t_z| > 0.866), `z_extent_um`, `fit_rms_um`,
`n_links`, `link_cost_max`, `frac_in_label`, `frac_near_other_label`, `roll_from_axial_med_deg`.

### 5.2 Profiles and gaps

- Background b: the territory's 15th percentile of raw counts. Every 0.05 um: `core` (mean of I - b in the ellipse
  0.12 x 0.25 um), `core_se` (from the noise model), `integrated` (sum of (I - b) du dv in the ellipse
  0.35 x 0.7 um, counts um^2: the continuity signal), `<role>_core` for each profile role, `contaminated`,
  `strength`, `lineness`, `in_label`.
- Gaps: integrated profile smoothed over 0.1 um; M = its median over clean, non-end samples; a gap is a run below
  0.35 M, at least 0.25 um long, with (M - min) / SE > 3, not touching an end. A matched link counts as a gap only if
  it passes this test. Per strand: `n_gaps`, `gap_um`, `gap_fraction`, pieces (the strand split at its gaps):
  `n_pieces`, `longest_piece_um`.
- `profile_cv_corr = sqrt(max(cv^2 - mean(SE^2) / M^2, 0))` (noise-corrected heterogeneity; replaces raw voxel CV).

### 5.3 Width and multiplicity (QC only)

FWHM lateral and axial from Gaussian fits of the mean cross-section; deconvolved width
`sqrt(max(FWHM^2 - FWHM_psf^2, 0))` (the SC central region is below resolution). `multiplicity(s)` = integrated /
the image median over complete strands; a stretch of at least 0.5 um at 1.7 or above sets `merged_suspect`.

### 5.4 Foci (`foci.py`)

For every spot table (RAD-51 `spots`, COSA-1 or other `spots_<name>`): the nearest point on the strands of the
same nucleus by `to_strip`; assigned when |u| <= 0.3 um and |v| <= 0.6 um; `ambiguous` when a second strand is within
1.25 x the PSF-scaled distance; else `off_strand` (or `no_strand`). Columns `s_um`, `s_frac`,
`s_folded_frac = min(s, L - s) / L`, `u_um`, `v_um`, `dist_um`, `in_gap`, `rank_on_strand`, `next_focus_um`,
`next_focus_frac`. s = 0 is the end with the smaller rounded (z, y, x); folded and inter-focus values do not depend
on orientation. Optional cross-check: peaks of the foci profile (`find_peaks`, prominence 4 SE).

### 5.5 QC flag vocabulary (`qc.py`; levels info, review, fail)

| code | scope | level | fires when (key, default) |
|---|---|---|---|
| fragment | strand | info | length_um < classes.min_complete_um (1.0) |
| low_snr | strand | review; fail below the fail level | snr_core_median < qc.low_snr_review (4.0); < qc.low_snr_fail (2.5) |
| end_unfit | strand | review | an end is free_unfit |
| contact_end | strand | info | an end is contact |
| truncated_z | strand | fail for length | an end is stack_edge, or a sample lies within qc.stack_edge_um (0.4) of the first or last plane |
| costly_link | strand | review | link_cost_max > classes.max_link_cost (60) |
| contaminated | strand | review | frac_contaminated > classes.max_contaminated (0.3) |
| merged_suspect | strand | review | multiplicity >= classes.merged_multiplicity (1.7) over >= classes.merged_run_um (0.5) |
| double_width | strand | review | fwhm_lat or fwhm_ax > qc.double_width_ratio (1.6) x the image's complete-strand median |
| exits_label | strand | review | frac_in_label < classes.min_frac_in_label (0.8) |
| neighbour_nucleus | strand | fail | frac_near_other_label > 0.5 |
| steep_z | strand | review | frac_steep60 > qc.steep_frac_review (0.5) |
| long | strand | review | length_um > qc.long_um (10.0): two chromosomes joined end to end |
| kink | strand | review | turn > qc.kink_deg (60) within qc.kink_um (0.3) |
| fit_residual | strand | review | fit_rms_um > qc.fit_rms_review_um (0.12) |
| gap_long | strand | review | a gap >= qc.gap_long_um (0.6) |
| saturated | strand | review | more than 1 % of core samples at the detector maximum (bit depth from the acquisition record, else the dtype) |
| imported, verified | strand | info | provenance |
| error:Type | strand, nucleus | fail | an exception |
| label_merged_suspect | nucleus | review | volume > qc.merged_label_ratio (1.6) x the median of traced nuclei |
| truncated_z | nucleus | fail for totals | the label touches the first or last plane |
| edge_xy | nucleus | review | the label touches the image border |
| count_low, count_high | nucleus | info | n_main below or above expected |
| noise_mismatch | nucleus | review | per-nucleus GAT noise SD outside [0.7, 1.4] |
| budget_exceeded | nucleus | review | work above budget_s |
| sc_straighten:uncalibrated | image | info | always in v1, until the protocol scores the method |
| noise_fit_poor, psf_fallback | image | review | 4.2, 4.9 |
| sc_channel_suspect | image | review | among nuclei with SC present, fewer than 20 % have a strand of 2 um or more |
| foci_uncalibrated | image | info | a spot table was used |
| geometry_from_manifest | image | info | run geometry from run_manifest.json |
| label_merged_suspect=N, nucleus_errors=N | image | info | counts |

`qc_status` = fail if any fail flag, else review if any review flag, else pass. A test iterates `qc.CODES` with one
positive and one negative phantom per code, and another keeps the table in docs/SC_STRAIGHTEN.md equal to the code.

### 5.6 Strand class (abstain rather than mislead)

- `verified`: imported with verified=true (overrides the rest).
- `fragment`: length under 1.0 um.
- `merged_suspect`: the flag of that name.
- `complete`: both ends `free_fit`, link_cost_max <= 60, frac_contaminated <= 0.3, frac_in_label >= 0.8, no review or
  fail flag.
- `partial`: everything else.

Per-chromosome statistics use only `complete` or `verified` strands. The precision of `complete` is a protocol
metric (S1e, S4).

### 5.7 Nucleus tier

- `status`: traced, no_sc (no seed in the territory, or total under 1.0 um, or median core SNR of its strands
  under 3), skipped (selection), failed. A status never zeroes a column: a no_sc nucleus keeps its weak strands
  (flagged low_snr) and their total.
- Counts: `n_strands`, `n_main` (strands of 1.0 um or more), `n_complete`, `n_fragments`, `n_merged_suspect`,
  `expected` (6 oocyte, 5 spermatocyte, from germ_cell), `count_status` (low, ok, high: n_main against expected),
  `complete_nucleus` (n_complete == expected and every main strand complete and no nucleus fail flag; the Cahoon
  2023 male rule, drop the shortest of six, is an option).
- Length: `total_length_um` (sum after accounting, every traced nucleus), `total_length_xy_um`,
  `total_length_complete_um` (only when complete_nucleus).
- Continuity: pieces are the strands split at their gaps; `n_pieces` (pieces of 0.25 um or more),
  `mean_piece_um`, `longest_piece_um`, `gaps_per_10um` (10 x n_gaps / total_length_um; NA under 1 um), `gap_fraction`,
  and **`piece_lw_mean_um` = sum(l_j^2) / sum(l_j) over all pieces of the nucleus** (0 when there are none): the
  length of the piece that holds a random point of the nucleus's SC. It is the score the protocol uses for SCfull
  against fragmented nuclei (S2a).
- Intensity: `on_strand_frac` (B): the sum of max(I - b, 0) over territory voxels inside the union of strand tubes
  (ellipse 0.35 x 0.7 um around each centreline) divided by the same sum over the territory; `snr_core_median`;
  `coverage_proxy` (share of crest voxels with S >= seed_min within 0.3 um of a strand); `multiplicity_frac`.
- Context: `volume_um3`, `in_germline`, `axis_position_um`, `zone`, `is_pachytene` (copied from the run), label
  flags, `qc_flags`, `elapsed_s`.

## 6. Outputs

All files are new and sit in the image's results folder (or the mirrored folder under `--out`). Tables are CSV
(`float_format="%.6g"`). Per-strand, per-gap, per-focus and per-nucleus tables append SHARED_META copied from the
run's nuclei.csv (it describes the source run) plus `sc_straighten_version` and `sc_straighten_hash` (the
sub_hash); the sample table carries only `image_id` and keys and joins to SHARED_META through the strands table.

| file | one row per | key columns |
|---|---|---|
| `<id>__sc_strands.csv` | strand (duplicates removed in 4.7 are not written) | strand_id, strand_uid, nucleus_id, marker, geometry_source (auto, imported), verified, strand_class, qc_status, qc_flags, the 5.1 measures, end0_class, end1_class, end0_z/y/x_um, end1_z/y/x_um, core_median, integrated_median, snr_core_median, profile_cv_corr, n_gaps, gap_um, gap_fraction, n_pieces, longest_piece_um, multiplicity_median, fwhm_lat_um, fwhm_ax_um, frac_contaminated, n_foci, strip_offset, strip_len |
| `<id>__sc_strand_samples.csv` | strand x s every 0.1 um | image_id, strand_id, strand_uid, nucleus_id, i, s_um, in_strand, z_um, y_um, x_um, t_z, t_y, t_x, nlat_z, nlat_y, nlat_x, nax_z, nax_y, nax_x, core, core_se, integrated, foci_core, dna_core (NA when the role is not read), contaminated, in_gap, in_label, strength, lineness, multiplicity, roll_from_axial_deg |
| `<id>__sc_strand_gaps.csv` | gap | gap_id, strand_id, strand_uid, nucleus_id, s_start_um, s_end_um, length_um, depth, z_score, bridged, mid_z_um, mid_y_um, mid_x_um |
| `<id>__sc_strand_foci.csv` | spot x spot table, in processed nuclei | spot_table, spot_id, nucleus_id, marker, assignment (unique, ambiguous, off_strand, no_strand), strand_id, strand_uid, s_um, s_frac, s_folded_frac, u_um, v_um, dist_um, second_strand_id, in_gap, rank_on_strand, next_focus_um, next_focus_frac, z_um, y_um, x_um |
| `<id>__sc_strand_nuclei.csv` | considered nucleus | the 5.7 columns |
| `<id>__sc_strand_strips_flat.tif` | | uint16 (rounded raw counts; 0 = outside the image), zlib; axes (role, u, s_concat): 25 rows; each strand occupies columns strip_offset to strip_offset + strip_len, s from -0.5 um to L + 0.5 um at 0.05 um; description JSON: format, version, roles, du, dv, ds, s_offset_um, [strand_id, offset, len] |
| `<id>__sc_strand_strips_axial.tif` | | the same with v (21 rows) |
| `<id>__sc_strand_overview.png` | | SYP maximum projection (xy stride 2) with strands coloured by class, failed nuclei outlined, counts box; not part of any golden |
| `<id>__sc_straighten.json` | image | the record (section 9) |

Size estimate for a 2025-11-05 gonad (about 17,000 um of strands): samples about 45 MB, strips about 70 MB before
compression. `germquant collect` writes `batch_sc_strands.csv`, `batch_sc_strand_nuclei.csv`,
`batch_sc_strand_gaps.csv`, `batch_sc_strand_foci.csv` and `batch_sc_straighten_summary.csv` (from the records)
only when sc_straighten records exist, finding images by the record's `image_id` (not by splitting names). Samples
are not stacked. Names avoided: `sc_ribbon` (the coloc mask) and `track` (sc_tracks and time-lapse ids). No new
name ends in `__nuclei.csv`, `__spots.csv` or another existing table suffix, and none contains `__spots_`.

## 7. napari

- `layers.strand_layers(results_dir, image_id, nucleus_ids=None)`: Shapes paths in world um (z, y, x) with the
  strands table as features, coloured by class; Points for end classes; Points for gap midpoints; Points for foci
  coloured by assignment.
- `layers.strip_layers(results_dir, image_id, nucleus_id)`: a gallery Image built on the fly from the concatenated
  strips (strands of the nucleus stacked as rows, padded with NaN; scale du, ds) and the foci as Points at
  (row, u, s).
- `layers.image_layers(...)`: the SYP crop with scale = spacing and translate = origin x spacing.
- `germquant sc-view` imports napari lazily and says, when it is missing, that it runs from the napari-germquant
  environment. The plugin's Results dock will reuse `layers.py` (plugin PR).

## 8. Config

A new top-level section, read with `cfg.get("sc_straighten.<key>", code default)`; no shipped config or profile
names it, so no config_hash changes. Precedence: code defaults < the section of `--config` or `--profile` <
`--set`. Axis order is (z, lateral) wherever two values appear, (z, y, x) for three. If a pack recipe ever sets
these keys they go to `packs/lint.py` ALLOW_KEYS, never to WORM_SWITCHES.

```yaml
sc_straighten:
  enabled: false                  # reserved for the phase-2 stage; the command ignores it
  nuclei: germline                # germline | pachytene | all
  territory_reach_um: 0.4
  profile: {roles: [foci], ds_um: 0.05, core_um: [0.25, 0.12], integrate_um: [0.7, 0.35], background_pct: 15}
  noise: {planes: 12, bins: 30, qmax: 0.85, poor_rel_rms: 0.15, nucleus_sd_range: [0.7, 1.4]}
  fields: {sigma_um: [0.2, 0.108, 0.108], null_quantile: 0.999, null_shape: [40, 96, 96], null_seed: 0, crest_offset_max: 0.6}
  follow: {seed_min: 4.0, lo: 1.5, step_um: 0.1, momentum: 0.7, max_correction_um: 0.15, lineness_max: 0.5,
           weak_run_um: 0.6, claim_vox: [1, 2, 2], contact_cos: 0.85, loop_back_um: 0.6,
           min_length_um: 0.5, max_length_um: 30.0}
  link: {max_gap_um: 1.0, max_angle_deg: 50.0, max_skew_um: 0.3, gap_cost_per_um: 40.0, skew_cost_per_um: 100.0,
         unlinked_cost: 40.0, dir_window_um: [0.15, 0.9]}
  dedupe: {lateral_um: 0.15, axial_um: 0.3, min_run_um: 0.4, parallel_cos: 0.9, min_remainder_um: 0.5}
  refine: {ds_um: 0.1, window_um: [0.30, 0.15], grid_um: [0.6, 0.4], max_shift_um: 0.3, iterations: 4,
           knot_um: 0.3, sigma_loc_um: [0.06, 0.03, 0.03], lambda: 1.0}
  ends: {anchor_um: 0.6, inside_um: 0.6, outside_um: 1.4, step_um: 0.05, max_recentre_um: 0.2,
         sigma_bounds: [0.6, 1.6], fit_min_amp_se: 5.0, fit_max_rel_rms: 0.25}
  psf: {mode: auto, fallback_sigma_um: [0.365, 0.15], min_stretches: 20, flat_tz_max: 0.3, isolation_um: 1.0, min_snr: 10.0}
  gaps: {q: 0.35, min_um: 0.25, z: 3.0, smooth_um: 0.1}
  contamination_um: [0.7, 0.35]
  strips: {write: true, ds_um: 0.05, du_um: 0.05, dv_um: 0.1, half_um: [1.0, 0.6], extend_um: 0.5,
           flat_v_um: 0.3, axial_u_um: 0.15, roles: [central_element, foci]}
  samples: {step_um: 0.1}
  foci: {tables: auto, max_um: [0.6, 0.3], ambiguity_ratio: 1.25}
  classes: {min_complete_um: 1.0, max_link_cost: 60.0, max_contaminated: 0.3, merged_multiplicity: 1.7,
            merged_run_um: 0.5, min_frac_in_label: 0.8}
  qc: {low_snr_review: 4.0, low_snr_fail: 2.5, steep_frac_review: 0.5, long_um: 10.0, kink_deg: 60.0, kink_um: 0.3,
       fit_rms_review_um: 0.12, gap_long_um: 0.6, stack_edge_um: 0.4, merged_label_ratio: 1.6, double_width_ratio: 1.6}
  nucleus: {present_min_um: 1.0, present_min_snr: 3.0, expected: {oocyte: 6, spermatocyte: 5}, male_x_rule: false}
  budget_s: 60
  workers: 1
  overview: true
```

Tunable on development data (synthetic development seeds and development images): follow.seed_min, follow.lo,
link.*, refine.knot_um, refine.lambda, refine.sigma_loc_um, gaps.*, classes.*. Everything else is fixed unless a
changelog entry of the protocol says otherwise.

## 9. Provenance, resume, performance

- Record `<id>__sc_straighten.json`: format `germquant.sc_straighten`, version, image_id, status (done, failed,
  refused), done, sub_hash, params_hash, params (resolved), inputs (results_dir; nuclei.csv, labels TIFF and each
  spot table with path and sha256; source image path, size, mtime, head sha256; xy_stride, z_range, geometry
  source; channel map; roles read; imported files with sha256), tool versions (numpy, scipy, scikit-image,
  networkx, tifffile, h5py, nd2), git sha of this run, noise model (a, b, rel_rms, bin table), null_q, PSF estimate
  and source, counts, flags, nucleus errors, elapsed time, workers, sealed_override, finished.
- `sub_hash` = sha256 over VERSION, params_hash, tool versions, input digests, run geometry, geometry mode and
  imported digests, nucleus selection. Not the git sha and not the worker count. `--resume` (default) skips an
  image whose record is done with the same sub_hash; `--force` reruns.
- Writes are atomic (temporary file in the same folder, flushed, `os.replace`, through `fsutil.long_path`); the
  record is written last.
- Performance: Design A measured 0.1 to 1.1 s per nucleus for the chain without strips; expect 10 to 20 min per
  2025-11-05 gonad on one worker and about 5 min on four. Parent memory: two channels plus uint16 labels (2.9 GB on
  2025-11-05, 4.7 GB on 2026-01-05) plus crops; peak under 6 GB. CPU only. Default `workers: 1` because the machine
  is shared; development runs use `--workers 4`.

## 10. Validation hooks

- `benchmark/sc_split.py`: loads `tests/golden/sc_split.json`, checks it against the registered sha256 (protocol
  changelog entry 1), answers `role_of(stem)`, raises `SealedError` for test and reserved gonads unless `final=True`.
- `benchmark/sc_ribbons.py`:
  - zone truth: per-gonad Position sheet (row order, names A, B, ...), master table rows filtered to the requested
    gonads before anything else is done with them, frame conversion (protocol 4.1), projection with
    `staging.zones.project_to_polyline`, zone labels and guard bands (4.2), per-gonad AUC, the boundary estimator
    (5.2), cluster bootstrap from `benchmark.stats`;
  - synthetic scoring: one-to-one matching of strands to truth ribbons, coverage, off-truth length, gaps, foci,
    class precision (protocol 4.4);
  - filament check: territories from the truth nuclei, `straighten_image` on the TIFF, matching (protocol 4.3).
- `scripts/benchmark_sc.py manifest | synthetic | zones | filaments | report`: development by default; test or
  reserved roles and sealed seeds only with `--final --reason TEXT`, which appends a line to
  `C:/Users/ryane/germquant_golden/benchmarks/sc_ribbon_v1/final_runs.jsonl` (commit, params_hash, split sha256,
  time, reason). Development outputs go to `.../benchmarks/sc_ribbon_v1/dev/`, every configuration tried is a line
  of `dev_trials.jsonl`.

## 11. Tests and goldens

All CPU, small, without heavy parallelism.

- `test_sc_straighten_geometry.py`: P-spline length of an arc (r 1.2 um, L 5 um) within 0.005 um noise-free and
  0.05 um with 0.02 um noise; frames orthonormal and twist-free on a helix; a flat curve gives n_ax = +/-z;
  image -> strip -> image round trip under 1e-6 um; erf end fit bias under 0.02 um at SNR 10; GAT variance within
  5 % of 1; crop-origin frame mapping; end matching identical under permutation of input order; arc length of a
  helix under anisotropic sampling within 0.5 %.
- `test_sc_straighten_synthetic.py` (exact-render scenes from `synth.py`, development seeds): isolated arc at real
  noise: 1 strand, length error <= 0.2 um, centreline distance <= 0.05 um; 60 deg tilt <= 0.4 um; X crossing: 2
  strands, correct continuation; T contact: a `contact` end; 0.5 um gap: 1 strand, 1 gap of 0.5 +/- 0.15 um; foci
  at s 1.0 and 3.0 um within 0.1 um; noise-only nucleus: 0 strands; a fixed dense nucleus: total within 10 %.
- `test_sc_straighten_qc.py`: one positive and one negative phantom per `qc.CODES` entry; the docs table equals
  the code.
- `test_sc_straighten_cli.py` on the end-to-end synthetic gonad of `test_end_to_end.py`: `process_image`, then
  `run_results`: the sha256 of every pre-existing file is unchanged; only the declared new files appear with the
  declared columns; resume does nothing; a parameter change changes sub_hash; `--force`; workers 1 and 2 give
  identical files; refusals (missing source with the --source-map message and exit 1, stride 2, series folder,
  no central_element, golden folder without --out, sealed image without --allow-sealed); an image id ending in
  `_` and containing a space; a folder without a done marker uses run_manifest.json.
- `test_sc_straighten_inputs.py`: the channel-subset dims logic on fake dask arrays; page-wise uint16 labels; a
  workstation-only equivalence with read_stack on herm_000 (skipped when the file is absent).
- `test_sc_straighten_import.py`: a tiny .ims written with h5py (Filaments0 Vertex, Edge, DendriteSegment, a
  branch, a single-vertex stub), SWC and CSV readers, the continuous frame rule.
- `test_sc_straighten_fence.py`: importing germquant.sc_straighten loads no napari, Qt, torch or CuPy; collect on a
  tree without sc_straighten records writes exactly the old file set; layers return valid tuples without napari.
- `test_benchmark_sc.py`: polyline frame conversion; zone labels by the R rule; guard bands; per-gonad AUC with
  ties; the boundary estimator on synthetic step profiles; SealedError for a test gonad; synthetic matching.
- `test_sc_split.py`: the record's counts and rules; the four RAD-51 test gonads are test; no RAD-51 test gonad is
  dev; the sha256 of the truth copies (workstation-only).
- Goldens: `tests/golden/sc_straighten/synthetic_gonad.json` (the end-to-end synthetic gonad: digests exact on the
  workstation, rounded values with GERMQUANT_GOLDEN_RTOL on CI) and `phantom_gonad.json` (a seeded 20-nucleus
  phantom: digests plus accuracy against its truth); real goldens for n2rad51_noHS_herm_000 and n2rad51_HS_male_000
  under `germquant_golden/<sha>/sc_straighten/<tag>/` by `scripts/make_sc_straighten_goldens.py`, compared with
  `regression_diff.py` (workstation only; the overview PNG excluded).

## 12. Build plan for this workflow

Branch `feat/sc-straighten` from origin/main 8eb73de, one PR to main (not stacked), plugin PR later.
Python: the main .venv with `PYTHONPATH=<worktree>/src`; nothing is installed into it (numpy, scipy 1.17.1,
scikit-image 0.26, networkx 3.6.1, h5py, tifffile, openpyxl and matplotlib are already there).

| package | content | depends on | acceptance |
|---|---|---|---|
| P1 core geometry | noise, fields, follow, link, account, refine, frames, synth; port from `scripts/designA/strandlib.py`, `tracklib.py`, `crestlib.py` | data classes agreed first (RawStrand, Strand, ImportedPath, Params) | geometry and synthetic tests pass |
| P2 image pipeline | params, inputs, crops, strips, measure, foci, qc, schema, outputs, `__init__`, cli_sc_straighten, cli.py and batch.py edits | P1 interfaces (stubs until P1 lands) | CLI, QC, fence, inputs tests pass; existing suite unchanged |
| P3 imports and viewing | imports, layers, sc-view | P2 tables | import and layer tests pass |
| P4 benchmark harness | benchmark/sc_split, benchmark/sc_ribbons, scripts/benchmark_sc.py, tests/golden/sc_split.json | P1 synth, P2 tables | benchmark and split tests pass; SealedError proven |
| P5 integration and development | merge P1 to P4; full suite; development runs on the eight 2025-11-05 development gonads (bc62476 inputs, `--out C:/Users/ryane/sc_ribbon_work/runs/<sha>/`); `calibrate_sc_straighten.py` on development seeds; development zone scoring; goldens; docs; PR | P1 to P4 | section 12.1 |

### 12.1 Acceptance of the PR

1. Existing tests pass with no edit; new tests pass on the workstation; CI green on Linux, Windows and macOS.
2. Every pre-existing file of a finished run is byte-identical after `sc-straighten` (fence test), and a default
   `germquant collect` writes the old file set.
3. The eight development runs finish with peak memory under 8 GB per process and under 20 min per gonad on four
   workers; no nucleus error goes unreported.
4. The development synthetic gates of protocol S1 are reported on development seeds; any gate that fails is listed
   with the claim it withholds.
5. Development zone metrics (protocol S2 on development gonads only) are in `dev_trials.jsonl`, with the frozen
   values of section 9 of the protocol prepared but not yet entered.
6. Nothing was run on a test or reserved image, no test or reserved truth was read, no sealed seed was generated.

## 13. Risks, open questions, items for Ryan

### 13.1 Risks

- Per-ribbon decomposition in dense oocyte nuclei (14 % whole and pure on synthetic): per-chromosome numbers stay
  limited to `complete` strands until the review loop exists.
- Recovery below about 300 counts of ribbon peak (-9 % at 300, -13 % at 200 on synthetic): the SNR floor.
- Optical merging (axial FWHM about 0.9 um): `merged_suspect` and `double_width` flag it; it cannot be resolved.
- Steep ribbons: 17 to 33 % of skeleton length runs steeper than 45 to 60 deg on herm_000; length is less certain
  there (+0.39 um at 60 deg at real SNR before tuning).
- Label errors (merged and clipped labels): label flags; the territory reach covers the median 3 % of ridge
  signal outside the label.
- The synthetic geometry is still a model; real decomposition rates need hand traces (protocol S4).
- Heat-shock morphology (broadened or aggregated SYP) is outside the phantom families; `on_strand_frac` and the
  flags describe it, lengths there are not claimed.

### 13.2 Open questions

- Whether `piece_lw_mean_um` separates the zones better than mean SYP intensity (the protocol reports both).
- Whether the skew criterion helps (development sweep).
- Whether the male X rule of Cahoon 2023 should be on by default.

### 13.3 Items for Ryan

1. **Number-changing fix outside this PR:** `batch.collect` and `staging.tracer.results_images` misparse image ids
   that end in an underscore, so every replicate-000 image (4 of the 12 RAD-51 gonads) is silently dropped from
   `batch_summary.csv` and every `batch_*` table, and gets a wrong id in `germquant trace`. A fix changes `batch_*`
   outputs of any study with such names. OK to fix in its own PR?
2. **Madeleine's summary code:** zone starts are taken one measurement point after the named point (0.3 item 10).
   Intended?
3. **Truth requests:** re-save the three empty Filament files; independent 3D hand traces on fixed images
   (protocol section 11): who traces, and how many hours.
4. **Acquisition:** a noise variance slope of about 4 counts per count suggests EM gain or a high-gain camera mode
   (Design A); worth checking the camera settings.
5. **Split:** approve the registered split (18 development, 12 test, 9 reserved N2 gonads; IR and CCW69 reserved).
