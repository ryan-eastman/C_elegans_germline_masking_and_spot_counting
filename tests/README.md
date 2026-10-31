# The test suite

How the suite is laid out, what the markers mean, and how to run it on the workstation, on a Mac
and on GitHub Actions. Everything here runs on the CPU and needs no microscopy data: the synthetic
gonads are built inside the fixtures.

## What is in here

- One module per stage or concern (`test_axis.py`, `test_coloc.py`, `test_granule.py`,
  `test_partition.py`, `test_sc.py`, `test_staging.py`, `test_spots.py`, ...), plus the config and
  profile machinery (`test_profiles.py`), the stage registry (`test_stages.py`), batch bookkeeping
  (`test_batch.py`) and the model download command (`test_models_registry.py`, which fakes the
  network).
- `test_imaris_truth.py`: the Imaris truth reader (`germquant imaris-truth`) on a tiny synthetic
  HDF5 written in the `.ims` layout; skips itself when `h5py` is not installed.
- `test_end_to_end.py`: the whole single-image pipeline on synthetic 3-channel and 4-channel gonads
  with the `.nd2` reader monkeypatched. Its fixtures (`patched_reader`, `patched_reader_4ch`) are
  imported by several other modules.
- `test_timelapse_reader.py`: the T-aware reader (`docs/TIMELAPSE.md`). A fake `ND2File` with nd2's
  one-chunk-per-frame layout records which planes a read touches, so the tests show that `t=<index>`
  reads one timepoint's planes and nothing else, that `t=None` is the old behaviour, that out-of-range
  values are refused before any pixel is read, and that the header reports `n_timepoints`, `period_s`
  and `t_seconds`. The bioio branch runs on a fake `BioImage` and, when `bioio-ome-tiff` is installed,
  on an OME-TIFF with a T axis written by tifffile. `process_image(t=...)` is checked to suffix the
  image id, record the timepoint in the markers only, and write the same numbers as a fixed run on the
  same stack.
- `test_live_profile.py`: the live GFP::SYP-2 profile and channel map (no `dna` role) and the
  `segmentation.nuclei.channel` key: a spy backend shows which channel the segment stage hands over
  (DAPI by default, the named role otherwise, channel 0 with a flag when the role is missing), and a
  synthetic two-channel series with a planted drift runs through `process_series` end to end on the CPU.
- `test_training_loop.py`: `export-training`, `train-nuclei` and `register-model`
  (`docs/TIMELAPSE.md`, "Segmenting without DAPI"). The export is checked file for file against the
  source stack and labels (dtypes, the isotropic xz / yz resampling, the manifest, zone tags, the seeded
  draw, the size cap); training runs with the three Cellpose seams (`_build_net`, `_train_seg`,
  `_predict`) replaced by fakes, so the split rule, the losses, the sha256 of every pair, the per-orientation
  F1 and the card are checked without loading weights or writing a 1.2 GB file; the registry test registers
  a model, verifies it with `get-model`, and shows the segment stage receives the registered file. One test
  is opt-in (below): a real training run of a randomly initialised Cellpose-SAM network.
- `test_golden.py`: the numeric gate. `test_synthetic_golden` runs the same two synthetic gonads
  with the classical nucleus segmenter and SpotMAX and compares every output table cell for cell
  with the references under `golden/synthetic/` and `golden/synthetic_4ch/`. A `config_hash.txt`
  sentinel next to each reference pins the text of `config/config.yaml` and the channel map, so an
  edit to either fails the test even when no number moved. The output schema is append-only: a new
  column is allowed, a renamed, removed or reordered column is not. `test_real_golden` is the opt-in
  real-data comparison (see the environment variables below). The goldens were produced on the
  workstation (Windows, the pinned `constraints/desktop-2026-09.txt` environment) and are checked
  with `check_exact=True`.
- `test_default_invariants.py`: the tranche 5 fence (`docs/DESIGN_system_packs.md`, section 11.2).
  It pins today's default behaviour, read from `golden/default_invariants.json`, which was made once
  at the T5 base commit (bc62476) on the workstation: the config hash of every shipped config and
  profile, the LF-normalised sha256 of every channel map and Snakemake profile, the literal default
  stage hashes with every tool version fixed, and three default runs with every tool version fixed
  (the two synthetic gonads of `test_golden.py`, whose recipe is copied into the fence, and a
  single-plane image through the bioio branch): the stage hashes each run writes, the provenance key
  sets, the stage statuses, the qc flags, the output file names, the column list of all twelve tables,
  the TIF arrays and the montage pixels. It also checks that `batch --resume` accepts the marker each
  run wrote (segment hash computed as `cli._batch` computes it) and a marker as the base writes it,
  and the legacy single-plane behaviour. A failing pin means today's behaviour moved: fix the code, not
  the JSON. A new recipe file arrives with its own pin; existing pins never change, and
  `scripts/t5_scope_check.py` fails a diff that changes or removes one.
- `test_t5_scripts.py`: `scripts/t5_scope_check.py` on throwaway git repositories (a frozen file
  edited, deleted or renamed fails; a new file in a frozen folder passes; the sealed benchmark folder
  takes no new file; a recorded base that has diverged from main charges only the branch's own
  changes; the fence files may not be deleted; the pins file is append-only, where a scope may switch
  between "all" and "workstation") and `scripts/t5_gate.py` (its verdicts with stub steps, the ruff
  version it picks and the file-list check).
- `test_planar.py`, `test_pack_knobs.py` and `test_projected_run.py`: the projected 2D path and the pack
  knobs of tranche 5 package P2a (`docs/PLANAR.md`). `test_planar.py` checks `germquant.planar` on its
  own: the plane selectors (a planted sharpest plane, mid, an index, max and mean, an index out of
  range), the area gate, `segment_planar`'s control flow against `segment_nuclei`'s, the 2D measurements
  and `check_planar`'s refusals. `test_pack_knobs.py` pins every knob's default to today's call: Cellpose
  mocked (the 3D call, the object diameter, the 2D call and its CPU retry), SpotMAX spied through fake
  `spotmax.pipe` and `spotmax.filters` modules (and, where SpotMAX is installed, the effect size on the
  sharpened image), `require_channel` and `qc.min_nuclei`. `test_projected_run.py` runs a synthetic
  yeast-like stack through the test recipe `fixtures/planar/projected_yeast_like.yaml` and compares the
  tables with the golden in `golden/synthetic_projected/`, whose `config_hash.txt` sentinel is pinned
  in `golden/default_invariants.json`; it also checks the run's contract and its refusals before any
  pixel is processed. Regenerate that golden with `GERMQUANT_UPDATE_GOLDEN=1`, as for `test_golden.py`.
- `test_recipe_keys.py`: `germquant.recipe_keys`, the strict checks of the same package's opt-in recipe
  keys (`require_channel`, `qc.min_nuclei`, `cellpose.object_diameter_um` and the optional SpotMAX
  keys): no shipped config or profile sets one; every malformed value, a quoted "false" included, is
  refused with its key named; and `process_image` refuses it before any stage runs or any file is
  written, also under `method: auto`, where the watershed fallback would otherwise hide a bad Cellpose
  setting.
- Tranche 5 Train B (single-plane images, labels from files, projects, the sample sheet and batch formats, the
  training and tuning engine, post-stages; docs/PLANAR.md, docs/LABELS.md, docs/PROJECTS.md, docs/TRAINING.md,
  docs/TUNING.md, docs/POST_STAGES.md):
  - single-plane images: `test_planar_single.py` (planar.decide, the spacing of a single-plane run, the `info`
    line), `test_planar_bioio.py` (single-plane OME-TIFFs through bioio, against the single_plane_legacy and
    synthetic_2d goldens; skips without the readers, which CI installs), `test_planar_spots.py` (today's SpotMAX
    calls without planar, the planar branch with it), `test_planar_run.py` (single-plane runs through
    process_image under config/profiles/generic_2d.yaml, the synthetic_2d golden, a default batch resume deciding
    as before) and `test_planar_series.py` (single-plane time series refused), with the helpers of
    `planar_helpers.py`;
  - labels from files: `test_labels_io.py` (the readers of segmentation.nuclei.labels_from) and `test_labels_run.py`
    (labels runs through process_image and `batch`, the synthetic_labels golden), with `labels_helpers.py`;
  - projects: `test_project_manifest.py` (manifest, lock, new, open, add, the light header reads),
    `test_project_recipe.py` (recipe edits, channel roles, project models), `test_project_light.py` (the frozen
    signatures of PLAN section 9.1 and the torch-free entry, in fresh interpreters), `test_project_run.py`
    (run_project, its run log, models, devices, time series, collect, export-commands, the api re-exports),
    `test_project_status.py` (every status and reason code), `test_project_shards.py` (shards and three parallel
    processes) and `test_project_parity.py` (the project route against `run --pack` and api.run, also with a
    post-stage on, roots that move, the Snakemake project mode), with `project_helpers.py`;
  - the sample sheet and batch formats: `test_sample_sheet.py` (the sheet end to end, the row hash, the
    synthetic_sheet golden) and `test_batch_formats.py` (a default batch discovers exactly what 8eb73de does; other
    formats with io.input_suffixes, their ids and the duplicate_id rule);
  - training: `test_training_export_geometry.py` (export-training by run geometry, the 3D export unchanged),
    `test_import_labels.py` (import-labels from TIFF, PNG and Cell-ACDC files) and `test_train_curve.py` (the card
    fields and the learning curve, with the Cellpose seams replaced);
  - `test_pack_derive.py`: `pack derive`, `--project` on the pack entry and the version features, with new cases in
    `test_packs_isolation.py` (the --project listing and a recipe's needs.post_stages load no torch);
  - tuning: `test_tune_effect_filter.py`, `test_tune_spots.py` (the planted optimum in 3D and 2D, the synthetic_tune
    golden), `test_tune_points.py` and `test_tune_seal.py` (the allow-list and the import order, in fresh
    interpreters), with `tune_helpers.py`;
  - post-stages: `test_post_stages_discovery.py` (discovery from metadata and the declaration file, never importing
    a post-stage, and the torch-free `python -m germquant.post_stages list` and `write`),
    `test_post_stages_runner.py` (the runner on a hand-made results folder: only the declared tables and
    __post.json, the sub-hash, resume, a failure recorded while the image stays done) and
    `test_post_stages_wiring.py` (the wiring into `run`, `batch`, `project run` and `project status`, an error of
    the runner that never fails the image, collect's batch_post tables, pack status, the series refusals, the
    synthetic_post golden), with `post_helpers.py` (it builds the post_pack fixture pack in a temporary folder at
    test time; none is committed);
  - `test_t5_style_check.py`: `scripts/t5_style_check.py`, gate G7, on throwaway repositories, and on this checkout
    on tranche 5 branch names only (decision D55 in docs/DESIGN_system_packs.md, section 14).
  The six goldens of Train B, each with its `config_hash.txt` sentinel pinned in `golden/default_invariants.json`:
  `golden/single_plane_legacy/` (made on the base code before P2b changed anything), `golden/synthetic_2d/`,
  `golden/synthetic_labels/`, `golden/synthetic_sheet/`, `golden/synthetic_tune/` and `golden/synthetic_post/`;
  the pins file also gains the config hash of `config/profiles/generic_2d.yaml`.
- `test_peak_filter.py`: the exact replica of SpotMAX's minimum-distance peak filter
  (`germquant.spots.peak_filter`, `docs/PERFORMANCE.md`) against the installed SpotMAX function on
  random inputs with tied intensities, duplicate points and boundary distances, its fallback domain,
  the context manager that installs it, the SpotMAX stage with and without it on a flooded synthetic
  volume, and the flood guard's warning. The SpotMAX comparisons skip when SpotMAX is not installed.
- `test_benchmark_split.py`: the pre-registered RAD-51 benchmark split (`golden/rad51_split.json`,
  `docs/BENCHMARK_PROTOCOL.md`): the 12 tags and their train/val/test roles, the sha256 helper, the
  Imaris creation-parameter reader on a synthetic HDF5 project, and the manifest builder
  (`scripts/benchmark_manifest.py`) on a temporary tree.
- `test_rad51_additions.py`: the additions of protocol changelog entry 9 (`golden/rad51_additions.json`):
  the registered split file byte-identical and its roles pinned, the record validated with roles fixed in
  code (train only), its refusals, the role guard, the sealed names, the builder on stub files and the
  manifest with additions keeping the registered entries identical; since protocol changelog entry 10 the
  spot-truth rule (the export's own checks, the decisions file `golden/rad51_additions_decisions.json`, the
  section 4.4 minimum), a decision removing an addition's fold, and the truth read guard. The check that the 12
  recorded files exist with their sizes, and the decisions on the real truth files, skip when the 2026-01-05
  copies and truth exports are not on the machine. Since protocol changelog entry 11 the decisions file holds
  the decision that keeps n2rad51_20260105a_HS_male_001.
- `test_rad51_spot_sources.py`: the recorded spot sources of protocol changelog entry 11
  (`golden/rad51_spot_sources.json`, n2rad51_HS_male_000's 2025-12-31 export): the record against the split and
  the pins in code, its refusals (a test or val gonad, an unpinned gonad, another export name, size or sha256),
  the read guard that admits exactly the pinned export and keeps the other 2025-11-05 exports sealed, the
  usability rule, the manifest that changes only that gonad's entry (and the folds that then hold it and not
  n2rad51_HS_male_002), net3d's training gonads and the net3d fold order. The check of the pinned file skips when
  the copy is not on the machine. Since protocol changelog entry 12: the record's split digests are pinned as the
  historical ones it was made against (so the re-export of n2rad51_HS_male_002 does not break it; the additions
  record's test does the same), `load_spot_sources` refuses a malformed record, and a file with the pinned name and
  size but other bytes is hashed and refused by both guards before anything parses it.
- `test_gpu.py`: the real Cellpose-SAM 3D path. See the markers.
- `test_benchmark.py`: the benchmark harness (`germquant.benchmark`, see `docs/BENCHMARK.md`) on
  synthetic label volumes and spot sets with known answers for every metric (a split, a merge, an
  out-of-coverage prediction, the z-extent error, the chunked Hungarian against the dense result, the
  bootstrap intervals), plus both benchmark scripts run on a temporary manifest in a subprocess.
- `test_bench_detector.py`: the detector bench harness (`germquant.benchmark.detector`,
  `scripts/bench_cache.py`, `scripts/bench_detector.py`; `docs/BENCHMARK.md`, "Detector bench harness").
  A real pipeline run on a synthetic uint16 gonad is cached and re-run: the cached label image is the one
  the spots stage received, `run_backend` hands `detect_spots` exactly the pipeline's arguments and
  reproduces its spots table (log backend always, SpotMAX when installed). Tiny planted-foci gonads with
  truth files exercise the runner, the trial log, the test guard and the final-run log (on synthetic
  test-role gonads and temporary logs only), and the process pool against the sequential run.
- `test_net3d.py`: the EXPERIMENTAL net3d detector (`germquant.spots.net3d`; `docs/BENCHMARK.md`,
  "net3d"), on the CPU with tiny networks (skips itself without torch). The U-Net's shapes and the bound
  on its receptive field (a single perturbed voxel never reaches farther), tiled inference equal to
  untiled inference in xy and z and with the tile skip, peak finding on a synthetic heatmap (sub-voxel
  refinement, ellipsoidal suppression, the threshold property the validation sweep relies on, the
  nearest-label fallback and the contested flag), the exact dilations of the loss region and the
  weights and masked loss, the backend contract with a randomly initialised model (including the DNA
  channel through `channels`, the provenance hook, the spots sub-hash and the default path untouched),
  a pipeline run with `spots.backend: net3d`, and `germquant train-spots` for one epoch on tiny planted
  gonads (the test gonad is never opened, the same seed gives the same weights, the validation scorer
  agrees with the bench harness, and `scripts/bench_detector.py --backend net3d` scores the model).

- `test_bench_tuning.py`: the parameter search runner (`germquant.benchmark.tuning`,
  `scripts/bench_tune.py`; `docs/BENCHMARK.md`, "Tuning on train and val") and the two tuned profiles.
  The SpotMAX effect-size filter applied after a run with `effect_size_min` 0 reproduces a run made with
  that minimum, table for table (SpotMAX when installed); the ranking rule; on the planted-foci gonads of
  `test_bench_detector.py`, one trial record per parameter set, the same scores and detections as
  `bench_detector.py`, the process pool, the budget guard, the top-K selection for val and the test
  guard. `config/profiles/bench_log_tuned.yaml` and `bench_spotmax_tuned.yaml` load, extend
  `config/config.yaml` and change only `spots:` keys.

- The system packs (`docs/PACKS.md`, tranche 5 package P1), on the test pack in
  `fixtures/packs/test_pack` (a classical recipe that runs anywhere, and a Cellpose recipe with a 64-byte
  stand-in model) and the built-in `packs/celegans_germline`; `packs_helpers.py` holds the shared
  `pack_env` fixture, which isolates every pack root, the model store and the Cellpose cache:
  - `test_packs_format.py`: the pack.yaml rules (recipes at `config/profiles/<recipe>.yaml`),
    `validate_pack` (tunables the recipe does not set, or that are not scalars), LF-normalised SHA256SUMS
    (CRLF and LF copies agree, junk never hashed), the search order, built-in roots (the checkout and
    `$GERMQUANT_CONFIG_ROOT` only), shadowing, a broken pack listed as unusable, recipe resolution and the
    pack record;
  - `test_packs_extends.py`: the two extends fallbacks, every shipped config resolving to the same file
    as before, and the portable stage hashes of edited pack recipes when the engine root or the pack
    moves, also through a two-level chain (an absolute-path overlay keeps today's hashes);
  - `test_packs_models.py`: the blob store and fetch-models (Cellpose-cache adoption, `--from-file`, a
    local HTTP server with Range requests for download and resume, mismatches discarded, a wrong blob in
    the store removed, the size checked before any hash, nothing done twice), the verification record
    across processes and a store that cannot keep one, the pin check before segmentation (a pinned name no
    pack declares is refused), reserved names, classical runs that never resolve a model,
    `train-nuclei --base <pack model>` with the same pin check, and Cellpose's first call kept as the old
    line (rule R3);
  - `test_packs_status.py`: `pack check` blockers and fixes (modules, undeclared backends, the readers
    extra, models, the engine range, the recipe hash, faulty tunables, engine-file digests), image headers
    read without bioio (.nd2 and TIFF, the TIFF reader checked against bioio; other formats left to the
    deep check), the deep check's device report and accelerator warning, exit code 2 (an unknown pack
    included), and `run --pack` refusing an unusable pack unless `--no-pack-check`;
  - `test_packs_install.py`: install refusals (zip-slip, drive letters, symlinks and Windows junctions,
    bad sums, non-ASCII hashed text anywhere in the pack, recipes outside `config/profiles/`, faulty
    tunables, built-in collisions), a zip's junk left behind, immutability (a damaged installed copy
    included), vendor, every path of every committed SHA256SUMS tracked by git, and the .gitattributes and
    .gitignore rules in a scratch repository;
  - `test_packs_isolation.py`: in a fresh interpreter, `germquant.packs` and `python -m germquant.packs`
    load no torch, CuPy, SpotMAX or Cellpose, with `--image` too (no bioio, no dask); a deep check imports
    the pipeline (torch) before it loads any backend; `python -m germquant --help` imports torch first; the
    version entries; a simulated non-editable install inside another project;
  - `test_packs_run.py`: `--pack celegans_germline` against `--config config/config.yaml` (same hash,
    files, tables and done markers; the manifests differ only by the pack and runtime records), `batch
    --pack` records, the test pack's recipe goldens with an empty model store, and the tifffile example
    reader against the bioio reader;
  - `test_packs_lint.py`: model cards (`validate_card` on a real `train-nuclei` card), every lint rule,
    the `--previous` version rules, the base-weight comparison on small fake state dicts, and the result
    the C. elegans card records (cpsam_v2, from the run on the workstation CPU);
  - `test_packs_cli.py` and `test_datafetch.py`: the command surface and the resumable downloader (a
    server that answers a Range request from another byte is read again from the start).

- The outside packs and example fetching (tranche 5 package P3):
  - `test_outside_packs.py`: `packs/cultured_cells` and `packs/budding_yeast` load, verify and lint clean,
    are provisional with `validated: false`, carry the design's starting values and the QC settings of a field
    that is not a gonad (no deny-listed key, `expected_from_germ_cell: false` on every spots instance), and
    their channel maps bind as specified (no index fallback for the yeast transmitted role); the recipe
    goldens pass with an empty model store, an edit of a `pack_base.yaml` fails them, edits of a scratch
    copy of generic.yaml and config.yaml move no outside pack hash, a misspelt key in a `pack_base.yaml` fails
    lint, and a stack without a transmitted channel fails the yeast recipe with the require_channel message;
    the cards say what is provisional and which scale effect_size is on; the yeast golden sees its own spot
    threshold and smallest-cell gate (scratch copies with effect_size_min 4 or min_area_um2 1 fail on the
    tables); lint asks for require_channel and a current recipe golden, a permissive data licence on every
    real example, and keeps large files and stray TIFFs out of a pack; `pack golden` fails a recipe without
    its golden and exits 1 when none ran; every CITATION.cff has the fields CFF 1.2.0 requires; a frozen base is
    never rewritten silently; the cpsam comparison script refuses a device other than the one asked for;
  - `test_datafetch_zip.py`: one member of a remote zip by HTTP Range requests (classic and ZIP64 archives
    written by hand and by Python's zipfile; stored, deflated and bzip2 members), size, sha256 and CRC-32
    refusals, resume after a broken transfer (stored and deflated members), deflate64, LZMA, zstd and
    encrypted members refused with the method named, a server that ignores Range refused, sizes that cannot be
    right refused before any member data is requested, a deflate or bzip2 bomb stopped one byte past the
    recorded size, a pinned zip_size the server contradicts named, an oversized central directory refused, and
    bytes after the end record accepted;
  - `test_assemble.py`: `germquant assemble` and germquant.assemble (per-channel files, a single plane,
    channels as samples like DeepFoci, channels alternating page by page, crops, refusals) and pixels_sha256;
    a declaration that contradicts a file's own ImageJ axes is refused unless forced, and the output never
    replaces an input;
  - `test_pack_examples.py`: `pack fetch-example` against a local stand-in server (plain files and a zip
    member), the record beside the image (licence address and what was changed), nothing downloaded twice, a
    pack fix of channel names or voxel size reaching an image fetched before, wrong pins refused with nothing
    wrong kept, an unusable pack fetching nothing, an examples folder inside a pack root refused, loopback
    addresses only with `$GERMQUANT_ALLOW_LOOPBACK_EXAMPLES`, the torch-free entry in a fresh interpreter,
    `check_expect`, and the example-spec rules (quoted digests, plain and unique source names).

## Markers

Registered in `pyproject.toml` (`[tool.pytest.ini_options] markers`); CI passes `--strict-markers`,
so a misspelt marker is an error rather than a silently empty selection.

- `gpu`: needs a CUDA or Apple Metal accelerator. `tests/test_gpu.py` skips itself when none is
  present, and `pytest.importorskip("torch")` / `pytest.importorskip("cellpose")` at the top of the
  module turn a missing package into a skip instead of a collection error. On a machine WITH a GPU
  the test runs and occupies the card for a few seconds. To keep a run on the CPU (for example while
  a batch job owns the GPU) set `CUDA_VISIBLE_DEVICES=-1`; an empty value does not work on Windows,
  which drops empty environment variables.
- `model`: needs large model weights, either the 1.2 GB nucleus model that `germquant get-model`
  downloads or the Cellpose-SAM weights that Cellpose fetches on first use. CI deselects it with
  `-m "not model"`. Today only `test_gpu.py` carries it; any future test that loads the fine-tuned
  model must carry it too, so CI never downloads a gigabyte.

Opt-ins driven by environment variables rather than markers:

- `GERMQUANT_GOLDEN=<golden tree>` and `GERMQUANT_CANDIDATE=<run tree>`: `test_real_golden` runs
  `scripts/regression_diff.py` on two `germquant run` output trees and fails on any difference.
- `GERMQUANT_UPDATE_GOLDEN=1`: `test_synthetic_golden` rewrites the synthetic references instead of
  checking them. Only for a deliberately number-changing step, run from the code that is being
  certified; see the docstring of `test_golden.py` and `scripts/make_goldens.py`.
- `GERMQUANT_WORKSTATION_PINS=1`: `test_default_invariants.py` enforces the pins marked
  `"scope": "workstation"` in `golden/default_invariants.json` (today the montage pixels, whose
  rendering may differ by platform). Without it those tests skip, and the skip reason says whether
  this machine agrees with the workstation pin. `scripts/t5_gate.py` sets it on the workstation; that
  script runs ruff (the version CI pins), `pytest -m "not model"`, the scope check and these pins
  (`--no-workstation` skips them on another machine) and, with `GERMQUANT_GOLDEN` and
  `GERMQUANT_CANDIDATE` set, the strict real-golden comparison and the file-list check, plus the
  comparison with an older golden tree when `GERMQUANT_GOLDEN_LEGACY` is set too. It prints one
  PASS/FAIL table for the PR description; a run that skipped ruff, pytest or the scope check ends
  INCOMPLETE, and setting only one of the two golden variables is an error. The scope check runs in CI
  on tranche 5 branches only, so every T5 package and integration branch name must contain `tranche5`
  or a path segment that starts with `t5` (for example `feat/t5-p1`).
- `GERMQUANT_TRAIN_SMOKE=1`: `test_training_loop.py::test_real_random_init_training_on_the_cpu` trains a
  randomly initialised Cellpose-SAM network (no weights downloaded) for two epochs on the CPU through
  `train_nuclei` and checks the model file, log and card. It takes about 90 seconds and writes a 1.2 GB
  model file under the test's temporary folder, which is why it is not in the default suite; set
  `OMP_NUM_THREADS` to leave cores free for other work.
- `GERMQUANT_SNAKEMAKE=<snakemake executable>`: `test_batch_formats.py` and `test_project_parity.py` run their
  Snakemake dry runs (the formats and pack mode, and the project mode) with that executable, else with one on PATH;
  without either they skip. Snakemake is never installed into the .venv: it lives in a throwaway environment (the
  Train B plan's decision D25).
- `GERMQUANT_TRAIN_CARD_BASE=cpsam_v2@0f1cc3f7` (with `GERMQUANT_MODELS` pointing at a store that holds that blob):
  `test_train_curve.py::test_card_with_a_real_pack_base` writes a model card with a real pack base; the Cellpose
  seams stay replaced, so nothing trains.

Other automatic skips: the extended-path check in `test_fsutil.py` runs only on Windows, its
POSIX twin only elsewhere, `test_staging.py` skips its comparison with the August staging tables
when they are not on the machine, and `test_benchmark_split.py` skips the check that the 24
registered `.nd2` and `.ims` files exist with their registered sizes when the N2 RAD-51 data is not
on the machine, as `test_rad51_additions.py` does for the 12 files of the additions.

## Running the suite locally

Workstation (Windows, the pinned `.venv`), from the repository root:

    .venv\Scripts\python.exe -m pytest -q

PowerShell, keeping the GPU out of it:

    $env:CUDA_VISIBLE_DEVICES = "-1"
    .venv\Scripts\python.exe -m pytest -q

Mac or Linux, in the activated environment: `pytest -q`. The whole suite is 1557 tests;
`pytest -q -m "not model"`, what CI and `scripts/t5_gate.py` run, passes 1547, skips 9 and deselects 1 in about
11 minutes on the workstation's CPU. The skips are the environment and opt-in skips of the lists above and
below (the workstation pins, the POSIX twin of the long-path check, the real goldens, the Snakemake dry runs
without `GERMQUANT_SNAKEMAKE`, the card with a real pack base, the training smoke). Useful selections:

    pytest tests/test_golden.py -q        # just the numeric gate
    pytest -q -m "not model"              # what CI runs
    pytest -q -m gpu                      # the accelerator path, on a machine that has one

From a git worktree whose editable install still points at the main checkout, put the worktree's
`src` first so the tests import the code under test:

    $env:PYTHONPATH = "<worktree>\src"
    python -c "import germquant; print(germquant.__file__)"    # must print the worktree path
    python -m pytest tests -q -p no:cacheprovider

Lint, the same command CI runs: `ruff check src tests scripts` (`uvx ruff check src tests scripts`
works without installing ruff). The rule set is written out in `pyproject.toml` under
`[tool.ruff.lint]` so a ruff release cannot move the gate.

## Continuous integration

`.github/workflows/ci.yml` runs on every push to `main` and to `feat/**` branches, on pull requests
and on demand (the "Run workflow" button). Four jobs, with Python 3.11: `ubuntu-latest`, `macos-15` (Apple
Silicon; `macos-14` until 2026-10-30, as GitHub retires the macOS 14 images on 2026-11-02) and `windows-latest`
twice, as shards 0 and 1, because one Windows job took 15.4 minutes against the 15-minute budget of gate G1:

1. `uv` installs torch and torchvision as CPU wheels from the PyTorch CPU index (on macOS, the arm64
   wheels from PyPI, the CPU and Metal build), then germquant with the `gpu`, `spots`, `sc`, `dev` and
   `timelapse` extras (Cellpose, SpotMAX and Cell-ACDC, skan, pytest, ruff and laptrack) against
   `constraints/ci-cpu-2026-10.txt` (macOS: `constraints/ci-mac-arm64-2026-10.txt`, the same pins): the
   validated workstation versions minus the CUDA-only packages. The uv cache directory is saved between
   runs. A test that cannot pass on macOS for an environment reason is marked with the reason.
2. The two bioio TIFF readers (bioio, bioio-base, bioio-tifffile and bioio-ome-tiff) against
   `constraints/readers-2026-10.txt`, the workstation .venv's versions, and a step that fails the job when they
   do not import: the tests that read single-plane and OME-TIFF files through germquant.io skip without them, so a
   broken install would otherwise pass as skips.
3. `ruff check src tests scripts`.
4. `pytest -q -rs -m "not model" --strict-markers`: the full suite including the synthetic golden
   tests; `test_gpu.py` reports itself as skipped. On Windows each shard takes every other file of the sorted
   `tests/test_*.py` list (shard 0 the first, third, ...), so the two run every test file once; a shard fails first
   when pytest would collect a test file outside that list, or when its own list is empty.
5. On tranche 5 branches only (a name that contains `tranche5` or `/t5`, or starts with `t5`), on the Linux job:
   the scope check (`scripts/t5_scope_check.py`: no frozen path modified or deleted, no fence pin changed) and the
   style check (`scripts/t5_style_check.py`, gate G7: no en or em dash in a line the diff adds, pure ASCII in the
   text files it adds under config/, packs/, tests/fixtures/ and tests/golden/ and in every YAML file it adds).
   Both always run, also after a failed Test step, and the step fails when either fails.

`QT_QPA_PLATFORM=offscreen` and `CUDA_VISIBLE_DEVICES=-1` are set for the whole job. So is `GERMQUANT_GOLDEN_RTOL`,
1e-12, and 1e-6 on macOS, whose arm64 runner sums float32 intensity means in another order: the synthetic goldens
are compared exactly on the workstation that made them, but the GitHub runners have other CPUs and sum floats in a
different vector order, which moved `coloc.pearson_r` by one part in 1e14 on the first CI run; the tolerance
allows that and nothing more (leave the variable unset locally so any real change fails). The README badge shows
the result of the latest run on the default branch.

The rehearsal below predates the Windows split and the readers step (both from Train B, 2026-10-23).
Rehearsal on the workstation (2026-10-09, the same commands in a fresh `uv` environment with a cold
cache): torch and torchvision 39 s, the rest of the packages 31 s, ruff clean, pytest about 120 s.
Three and a half minutes in all, well inside the budget; a GitHub runner is slower, but not by
enough to matter. That rehearsal also found that `load_config` used to read the config and
channel-map files with the machine's locale codec, so the `config_hash` (and the golden sentinel)
came out differently in UTF-8 mode, on Linux and on macOS than on the workstation (cp1252). The
codec is now pinned to the workstation's (`fsutil.CONFIG_TEXT_ENCODING`), and the golden tests pass
with and without `PYTHONUTF8=1`. Outputs made before the pin on a Mac, on Linux or in UTF-8 mode carry the
other hash for the same file (config.yaml: d20d4321faad instead of 84bdda5a5230); treat them as a different config version.

What CI does and does not prove: it shows that a clean machine can install the pinned environment
and that the numbers are reproduced exactly on it. The references were made on the workstation, so
a last-digit difference on another CPU or operating system would surface as a golden failure. That
is a finding, not noise: first check whether a pin moved, then decide whether the difference is a
platform effect that deserves a documented tolerance. Until the first green run on `ubuntu-latest`,
Linux is unrehearsed: the install was resolved for Linux with `uv pip compile`, but the suite itself
has only been run on Windows.

## The real golden image sets

Two sets of real gonads are frozen under `C:/Users/ryane/germquant_golden/8666f07/` (tags in
`tests/golden/real_manifest.json`):

* `n2dryice_*` (7 gonads plus 4 stride-4 variants): the 4-channel N2 set of 2026-06-25 (DAPI, SYP-3, PGL-1,
  lamin). That imaging session was a fixation test (dry-ice freeze-crack), not the lab's standard preparation.
  It is a regression input only: it shows whether a code change moves the numbers on a fixed input. Do not use
  it as biological reference data, as training or benchmark data, or as representative images in a paper.
* `n2rad51_*` (12 gonads): the 3-channel N2 RAD-51 cross-validation set of 2025-11-05 (DAPI, SYP, RAD-51),
  no heat shock and heat shock, hermaphrodites and males, imaged with the standard preparation and paired
  with Imaris annotations. These are both regression inputs and the benchmark set for spot detection.
  Their train/val/test split, with the size and sha256 of every `.nd2` and `.ims`, is
  `tests/golden/rad51_split.json` (rebuilt by `scripts/make_rad51_split.py`); the benchmark protocol is
  `docs/BENCHMARK_PROTOCOL.md`.

