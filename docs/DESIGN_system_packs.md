# Tranche 5 design: system packs, a 2D path, projects and training on outside systems

Status: P0, P1, P2a and P3 are built (feat/tranche5 at 0521abf, gate G4 passed on 2026-10-12); the rest is
design. Written 2026-10-11 and revised the same day after a completeness critic
and a verifier.
"T5" below means this tranche.

Basis, read only: the engine at C:/Users/ryane/germquant_wt_grid (branch feat/spotmax-grid, HEAD 3fcb4b5) and
the plugin at C:/Users/ryane/napari-germquant (HEAD aced6e4). This document merges three competing designs. It
keeps the risk-first design that both reviews chose, grafts the parts the reviews asked for, and fixes what they
found wrong. Section 16 lists every change and the reason for it. A completeness critic then checked the draft
against both repositories; section 17 lists what changed in response. A verifier then checked that revision, and
section 17.1 lists its fixes.

## 0. Summary

- A user picks a system pack (C. elegans germline, cultured cells or budding yeast), runs it in napari without
  code or headless with the same numbers, trains on their own data, and cites the exact recipe.
- First milestone (M1, about four weeks): an unvalidated, codeless run of a yeast pack and a cultured-cell pack on
  a public example image. `germquant run --pack` writes the same files.
- Then come projects, the training loop and the post-stage contract (M2), validated packs on public benchmarks
  (M3), and the paper artefacts with the v0.2.0 release (M4).
- The one rule: nothing that exists today changes. Every C. elegans table, label image and config hash stays
  byte-identical. New behaviour is reached only through a new key, file or command, and every new behaviour gets
  its own golden before anything depends on it.

## 1. Purpose and scope

### 1.1 Goal

The PI wants a modular tool that people outside C. elegans meiosis can use. A user picks a system pack (models
and presets trained or tuned for their system), runs it codeless in napari or headless with identical numbers,
and trains on their own data. The first outside systems are cultured cells and yeast, with public datasets as
benchmark data. The target is a Cell Reports Methods article. Its novelty is the framework (packs, provenance,
goldens, training loop, time-lapse) plus the in-house components, validated on several systems.

### 1.2 In scope

1. System packs: format, registry, pinned models, cards, validation records, versions, checksums, install.
2. A 2D path with two opt-in modes. Planar: single-plane images. Projected: objects segmented on one 2D image
   and copied through z, with spots still counted in 3D.
3. A C. elegans pack that points at today's files, plus two outside packs: cultured cells and budding yeast.
4. Public benchmark manifests, sealed splits, one shared scorer and validation records for the outside packs.
5. A project folder that the GUI creates and `germquant project run` runs unchanged.
6. The training loop on outside data: label import, export of the plane the model saw, fine-tuning, learning
   curves, and spot preset tuning from points.
7. Plugin: a pack picker, a runner that executes in the engine environment, a Project dock, a role editor, a
   parameter form, a spot tuning panel and correct 2D display.
8. Paper artefacts: a claim map, parity and device reports, a reproduction guide, a DOI per pack version, DOIs for
   the engine and plugin releases, a deposit of the derived benchmark outputs, a key resources table draft and a
   limitations list.
9. A post-stage contract: measurement steps that an outside lab adds as a separate distribution, run after
   process_image and writing only new tables (section 5.6).
10. CI on macOS for the engine, CI for the plugin, and a container image that carries the packs and the readers
    extra.

### 1.3 Out of scope

Each item can be added later without changing what T5 builds.

- Multi-scene and multi-position files. They are refused with a message.
- Segmentation input_roles (a cytoplasm channel plus a nucleus channel).
- bioimage.io, Spotiflow, StarDist and InstanSeg as backends. They are comparators only.
- Time-lapse for outside systems (yeast and cell tracking). Planar time series are tested, and refused if they
  fail.
- Per-image channel bindings outside the config hash (section 7.3 says why).
- Pack discovery through Python entry points.
- A local modules folder for post-stages: it would run code from a folder, which R5 forbids. Post-stages come only
  from installed distributions (section 5.6).
- Publishing the engine on PyPI or the plugin on napari hub. Version 0.2.0 is installed as an editable checkout
  (section 3.8, open question O21).
- The planned config hash bump and the cpsam drift fix (section 15). Outside packs carry their own frozen base, so
  the bump moves no outside pack hash (section 3.12).

### 1.4 Milestones

| milestone | what a user gets | packages | rough time |
|---|---|---|---|
| M0 | nothing visible; today's behaviour is pinned | P0 | 1.5 days |
| M1 | codeless and headless outside runs, not yet validated; macOS, plugin and container checks | P1, P2a, P3, P4 | about 4 weeks |
| M2 | projects, single-plane images, the training loop, spot tuning, post-stages | P2b, P2c, P5, P5b, P6 | about 4 to 5 weeks |
| M3 | outside packs 1.0.0, validated on public benchmarks; real-Mac smoke test; pilot lab | P7, P8 | about 3.5 weeks plus GPU |
| M4 | paper artefacts, deposits, DOIs, engine v0.2.0, plugin 0.2.0 | P9 | about 2 weeks |

The total is about 14 to 15 focused weeks plus RTX 5090 time and downloads; P7's code is written during M2. M1's
critical path is P1, then P3, then P4: 18 to 19 working days before integration. The revision after the
completeness critic added about 13.5 working days (section 17). If the schedule must shrink, these move to a
follow-up first: the Mac device report, the optional YeaZ backend, the BABY and Spotiflow scoring, and P5b if the
PI agrees (O18). Moving P5b also moves what depends on it: P6 acceptance 10, the post-stage listing in section
8.5, P6's dependency on P5b, item 9 of section 1.2 and the words "and post-stages" in claim C1. The macOS CI job
and the real-Mac smoke test are not on this list.

### 1.5 Paper claims and the artefacts behind them

- C1. One engine runs three systems through versioned, citable recipes. Artefacts: packs with SHA256SUMS, expected
  config hashes, CITATION.cff and one Zenodo DOI per pack version, plus Zenodo DOIs for the engine and plugin
  releases. Outside labs extend the engine without editing it, through backends and post-stages (section 5.6). The
  time-lapse module enters the same framework through the C. elegans pack's live_syp2 recipe; outside time-lapse is
  open question O16.
- C2. Accuracy on public ground truth, scored by one harness. The harness uses sealed splits and protocols
  written before scoring, and scores ceilings and comparators with the same code. Artefacts: benchmarks/ext/ and
  the validation records, deposited with the per-image outputs under a DOI (P9).
- C3. The training loop works outside the worm. Artefacts: import-labels, learning curves, spot tuning from points,
  model cards, one timed codeless correction session by a lab member who did not write the code, and an outside
  pilot lab's report (O20).
- C4. A GUI run and a headless run of the same project write the same files, and every number traces back to a
  project, a config hash, stage hashes and model digests. Artefacts: project folders, parity_report.json and
  device_report.json (agreement across devices, never claimed as identity).
- C5. The in-house components carry over. Artefacts: the z-column merge ablation on DeepFoci and yeast smFISH,
  and the golden and provenance machinery applied to every pack.

### 1.6 Precondition

T5 starts from main, and only after the tranche 4 and SpotMAX-grid work has merged into main. In the local refs on
2026-10-11, feat/tranche4 is 18 commits and feat/spotmax-grid 24 commits ahead of main. The two are already
stacked: `git merge-base --is-ancestor feat/tranche4 feat/spotmax-grid` is true, and spotmax-grid is 6 commits
beyond tranche4. Starting T5 on top of them would stack T5 on unmerged work, which the lab does not do. The order
is therefore:
1. merge feat/tranche4 into main;
2. retarget the spotmax-grid PR to main, rebasing it onto main if tranche4 was squash-merged, and merge it;
3. check that both heads reached main with `git merge-base --is-ancestor <head> main` (after a squash merge, check
   instead that main holds the same tree changes).

Only then does P0 record and pin the base commit (O1).

## 2. Rules that protect today's numbers

- R1 Frozen texts. Never edit config/config.yaml (config_hash 84bdda5a5230), config_ccw77.yaml,
  config_n2dryice.yaml, any existing file under config/profiles/ (generic.yaml, the bench profiles and the
  Snakemake profiles included) or any of the five files in config/channel_maps/, not even a comment. New recipes
  are new files. fsutil.CONFIG_TEXT_ENCODING (cp1252) and the Config.hash payload format do not change.
- R2 Inert defaults. Every new key is read with a code default that reproduces today's behaviour. New provenance
  (manifest keys, stage-hash payload keys, done-marker extras, qc flags, table columns) is written only when the
  new feature is active. spot_backend and segment_backend already follow this pattern.
- R3 New code in new modules. These keep their statements: segment_nuclei, nuclei._drop_small, the Cellpose 3D
  call, measure_objects, _merge_z_columns, the three TIFF writers, spots_to_image, nd2_reader.read_stack and
  read_nd2_metadata, and bioio_reader._spacing. Where a branch is unavoidable (cellpose_backend, detect_spots, the
  top of detect_spots_spotmax, the empty kw that _cellpose builds, the channel fallback of the segment closure, the
  few-nuclei threshold in qc.qc_flags), the default branch is the old lines verbatim, and a mocked-call test pins
  the default call.
- R4 New lookups run only where today's code raises. Extends resolution, model resolution and pack discovery
  never change what an existing file or name resolves to.
- R5 Packs and projects are data. Listing packs and shallow checks parse YAML and JSON with yaml.safe_load and
  never import a backend, torch, Cellpose, SpotMAX or CuPy. This holds for germquant.packs imported on its own and
  for its command entry, `python -m germquant.packs`. Importing germquant.api loads torch, because api imports
  pipeline, and `python -m germquant` enters cli.py, whose first import is torch. No code from a pack or project
  folder ever runs. This also protects the cuBLAS load order on the RTX 5090: torch must load before CuPy, or
  Cellpose silently drops to the watershed.
- R6 Third-party model code and post-stages ship as separate distributions on entry-point groups (the two existing
  backend groups and, from P5b, germquant.post_stages), each with its own constraints file, installed into a cloned
  environment (.venv_packs). The validated .venv changes only after the full golden gate passes in it.
- R7 Every new recipe ships a recipe golden (synthetic image, CPU, CI, config_hash sentinel). Every recipe with a
  model ships a model golden (small public crop, RTX 5090, tagged with the device). Packs at 0.1.x carry recipe
  goldens only; model goldens come with 1.0.0 (P8).
- R8 The RAD-51 seal. T5 never modifies tests/golden/rad51_split.json, src/germquant/benchmark/,
  src/germquant/validate/benchmark_split.py, docs/BENCHMARK_PROTOCOL.md or docs/FINAL_SCORING.md. A diff check
  enforces this (section 11.2) instead of a content pin, because the RAD-51 final scoring work is still active and
  may legitimately update the split record. New modules never import germquant.benchmark.detector, tuning or
  logo, or validate.benchmark_split, and never name the split file.
- R9 Text and line endings. Every hashed file that T5 adds or writes is pure ASCII, because hashed text is
  decoded as cp1252. (Some frozen files, config.yaml among them, already hold UTF-8 dashes; they stay as they
  are.) This checkout has core.autocrlf=true: text files are LF in the index and CRLF on disk. Every digest of a
  text file that T5 pins or ships is taken on LF-normalised bytes, so Windows and Linux agree. Config hashes
  already read with universal newlines.

What moves on purpose: nothing that an existing run writes. Two notes. Raising the version to 0.2.0.dev0 changes
the pipeline_version column, which both golden tools already treat as volatile. Packs pin Cellpose weights by
digest, while the generic profile keeps its current behaviour (section 15a).

## 3. System packs

### 3.1 What a pack is

A pack is a checksummed folder of static files that a user picks by name. Choosing a pack and one of its recipes
resolves to a profile file that load_config reads like any other profile. The engine then runs unchanged, so the
config hash equals that of `--config <the same file>`. A pack names backends and models by name, never by
callable, so the engine can check availability from metadata before a run.

### 3.2 Layout

```
packs/<id>/
  pack.yaml                         manifest (not hashed, ASCII)
  SHA256SUMS                        digest of every other file, text LF-normalised
  README.md                         the pack card, in plain language
  CITATION.cff  CHANGELOG.md
  config/profiles/<recipe>.yaml     recipe text (hashed)
  config/profiles/pack_base.yaml    frozen base of every outside recipe (hashed with it; section 3.12)
  config/channel_maps/<map>.yaml    channel map (hashed with the recipe)
  post_stages/<recipe>.yaml         optional: the post-stages a recipe switches on (not hashed; section 5.6)
  models/<model>.card.json          model cards; weights never live in a pack
  tuning/<recipe>_grid.yaml         the grid that spot tuning searches
  example/<name>.ome.tif            tiny synthetic example, under 2 MB
  example/expected/recipe/<recipe>/       nuclei.csv, spots.csv, image_summary.csv, config_hash.txt
  example/expected/model_cuda/<recipe>/   the same plus env.json, written on the RTX 5090
  validation/<recipe>.json          validation records, written only by the benchmark tool
```

Recipes sit in <pack>/config/profiles/, so config._base_dir_for already anchors their relative channel-map paths
to the pack folder. The loader needs no change for that.

### 3.3 pack.yaml, format 1

Abridged example: the yeast pack as it ships at M1.

```yaml
format: 1
id: budding_yeast                     # ^[a-z][a-z0-9_]{1,40}$
version: 0.1.0                        # semver, section 3.7
title: "Budding yeast: cells and spots"
summary: "Whole cells from one transmitted-light plane, spots counted per cell in 3D."
system: {organism: Saccharomyces cerevisiae, sample: "fixed or live cells", imaging: "widefield or confocal"}
object_kind: cell                     # nucleus | cell: what one row of the nuclei table is
maintainer: "Cahoon lab"
licence: MIT                          # of the pack text; model and data licences sit on their cards
citation: {text: "...", doi: null, cff: CITATION.cff}
requires: {germquant: ">=0.2.0.dev0,<0.3", extras: [gpu, spots, readers]}
root: pack                            # pack (default) | checkout (the C. elegans pack only)
expects: {pixel_um: [0.05, 0.15], z_step_um: [0.15, 0.5]}     # outside the range: a warning, never a block
recipes:
  - id: cells_spots
    title: "Cells and spots, z-stack"
    profile: config/profiles/cells_spots.yaml
    dims: projected                   # 3d | 2d | projected | 3d+t
    default: true
    status: provisional               # provisional | validated
    expected_config_hash: "<12 hex>"
    needs:
      segment_backends: [cellpose]
      spot_backends: [spotmax]
      models: ["cpsam_v2@0f1cc3f7"]
      modules: [cellpose, torch, spotmax]
      accelerator: optional           # optional | recommended
    tunables:
      - {key: spots.effect_size_min, label: "Spot strength threshold", type: float, min: 0.5, max: 10}
      - {key: spots.spot_radius_um, label: "Spot radius (xy)", type: float, unit: um, min: 0.05, max: 1.0}
      - {key: segmentation.nuclei.project_plane, label: "Plane to segment", type: choice,
         choices: [mid, sharpest]}
    validation: null                  # validation/cells_spots.json once a record passes
    tuning_grid: tuning/cells_spots_grid.yaml
models:
  - name: "cpsam_v2@0f1cc3f7"
    kind: cellpose
    file_name: cpsam_v2
    sha256: 0f1cc3f7ecdd8a037a57c6c48d9d8921391be4cbce3fa9f13c3e3a2e1253c667
    size: 1233586851
    sources:
      - {type: cellpose_cache, name: cpsam_v2}
      - {type: url, url: "https://huggingface.co/mouseland/cellpose-sam/resolve/<pinned commit>/cpsam_v2"}
    card: models/cpsam_v2.card.json
    licence: "weights BSD-3-Clause; training data includes CC BY-NC sets"
examples:
  - {id: synthetic, image: example/yeast_synthetic.ome.tif}
  - id: sbiad3750_sut509_p1
    base_url: "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/750/S-BIAD3750/Files/"
    sources:
      - {path: "SpotMAX_gt_dataset/Yeast_smFISH/SUT509_gene/Position_1/Images/<stem>_Bright.tif",
         size: 14701884, sha256: "<recorded at build>"}
      - {path: "SpotMAX_gt_dataset/Yeast_smFISH/SUT509_gene/Position_1/Images/<stem>_Cy3.tif",
         size: 14701884, sha256: "<recorded at build>"}
    assemble: {axes: ZYX, channels: [Bright, Cy3], spacing_zyx_um: [0.24, 0.0720635, 0.0720635]}
    pixels_sha256: "<recorded at build>"
    licence: CC0
    attribution: "Padovani et al. 2026, Science Advances; BioImage Archive S-BIAD3750"
    expect: {n_objects: [lo, hi], spots_per_object: [lo, hi]}
```

Here <stem> is Y7092_WT_SUT509_q570_1xprobe_SCD_08022022-10_s10 (the file list was read from the deposit).

The loader, packs.validate_pack, refuses:
- an unknown format major;
- a pack id, recipe id or model name that breaks its pattern;
- a model without sha256 and size, or whose 8-hex name suffix differs from the start of its sha256;
- a tunable key that the recipe does not set, or whose value is not a scalar;
- status validated, or validated: true, without a passing validation record.

### 3.4 Recipes, hashing and portable extends

- recipes[].profile is a path inside the pack (root: pack) or relative to the engine config root (root:
  checkout).
- pack.yaml, the card, model cards and validation records are never hashed. The recipe chain and its channel map
  are hashed by today's Config.hash rule. So `--pack X --recipe R` and `--config <the recipe file>` give the same
  config_hash.
- config._resolve_extends keeps today's candidates and their order. Only when all of them fail does it try two
  new forms:
  - `pack:<id>@<version>/<recipe>`, resolved through packs.resolve_recipe (imported lazily). The unversioned
    `pack:<id>/<recipe>` picks the newest compatible version and is meant for interactive use only.
  - For a bare name with no directory part: <root>/config/<name>, then <root>/config/profiles/<name>, for each
    built-in root of section 3.8 (the checkout that holds src/germquant, then $GERMQUANT_CONFIG_ROOT; never the
    working directory). config_roots() is not used here, because it also lists the working directory and an
    unchecked _repo_root. config_roots() itself stays as it is, so --profile keeps today's behaviour (R4).
  The error lists every candidate tried and the packs available.
- Effect: a pack anywhere can say `extends: generic.yaml`, and a project recipe can say
  `extends: pack:cultured_cells@1.0.0/nuclei_foci`. The extends line is hashed text with no machine path, so the
  hash is the same on every machine. Today a profile outside the checkout needs an absolute extends path, which the
  plugin's exported overlay writes. The built-in outside packs extend their own frozen base instead (section 3.12).
- Inherited paths and stage hashes. When a file extends a base in another folder, load_config rewrites the
  inherited relative paths (io.channel_map, segmentation.nuclei.cellpose_model, qc.exclusions_file) as absolute
  ones (config.py line 195; tests/test_profiles.py asserts it), and stage_hashes hashes the resolved sections. An
  edited project on celegans_germline would therefore hash the absolute path of
  models/models/germline_nuclei_combined in its segment sub-hash, which is_done compares, so every image would turn
  stale on another machine. The read sub-hash of any edited project would also embed the pack install path.
  load_config now records each anchoring with the original relative value. In a chain that uses one of the two
  new extends forms, stage_hashes hashes that relative value instead of the absolute path. The model digest stays
  in the segment payload, so a different file still changes the hash. Every existing chain, out-of-repo overlays
  with absolute extends paths included, keeps today's payload, so nothing reprocesses.
- The project writer, the plugin and `pack derive` always write the version into a pack: line, and pack lint
  refuses an unversioned one inside a pack. A project therefore never moves when a newer pack is installed.
- test_missing_extends_target_is_an_error and test_bare_extends_name_prefers_the_repo_config_dir pass unchanged.
- Config gains two instance attributes: pack_record (None) and anchors (the inherited paths load_config anchored,
  each with its original relative value, and whether the chain used a new extends form; empty by default). Both
  are set in __init__, so Config.__getattr__ is never consulted for them. Neither is data: as_dict(), the
  manifest's config block and the hash do not change.

### 3.5 Models: pinned names, one blob store, fetching

- Names. A pack model is named `<name>@<first 8 hex of its sha256>`, for example `cpsam_v2@0f1cc3f7`. The recipe
  text carries the name, so the config hash moves when the weights move, and two packs can share one model
  without a name clash. Release models (germline_nuclei_combined) keep their names.
- Store. Weights live once per digest at models_dir()/blobs/<sha256>/<file_name>, with <sha256>.json beside them
  (source, size, fetch date, the packs that use it). models_dir() already honours $GERMQUANT_MODELS.
- Fetch. `germquant pack fetch-models ID [--recipe R] [--from-file NAME=PATH]` tries, in order:
  1. the blob already in the store;
  2. Cellpose's own cache (CELLPOSE_LOCAL_MODELS_PATH or ~/.cellpose/models), after a size and sha256 check,
     hard-linked when possible and copied otherwise;
  3. a file given with --from-file, after the same check;
  4. a download of the pinned URL, resumable, with size and sha256 checks.
  A mismatch is discarded and never registered. A rerun does nothing. On this workstation ~/.cellpose/models
  already holds cpsam_v2 (1,233,586,851 bytes) and cpsam (1,233,587,898 bytes), so no download is needed. Cellpose
  itself downloads from https://huggingface.co/mouseland/cellpose-sam/resolve/main/; packs pin a commit instead of
  main.
- Resolution. pipeline.resolve_model_path keeps every step it has: an absolute path, $GERMQUANT_MODELS, a file
  under base_dir, the installation registry, the release models, the base_dir registry, and the refusal of
  unknown bare names that Cellpose would silently replace with its stock weights. Only where it raises today, a
  value of the form `<name>@<8 hex>` is looked up in the installed packs' model specs and resolves to the verified
  blob. A declared but missing blob raises FileNotFoundError naming `germquant pack fetch-models <pack>`, before
  segmentation starts (process_image has already read the stack; run --pack catches a missing blob earlier through
  pack status). training._resolve_base gets the same last step, so `train-nuclei --base cpsam_v2@0f1cc3f7` works.
- Classical runs. process_image resolves the model before it checks the method (pipeline.py line 289), and today a
  one-part name that is neither registered nor a Cellpose built-in raises (line 183). Under
  segmentation.nuclei.method classical, the new `<name>@<8 hex>` branch therefore returns the value unresolved and
  the pin check is skipped: the classical segmenter never loads the model, and model_sha256 is recorded as "unused"
  anyway. cellpose and auto still raise. Without this, no recipe golden could pass in CI, where no 1.2 GB blob
  exists and a sha256 pin cannot be faked with a small file.
- Pin check. Whatever a `<name>@<8 hex>` value resolves to, a local registry entry included, must have a sha256
  that starts with the suffix and equals the spec's digest, or the run stops before any stage processes pixels
  (classical runs excepted, as above). pack status never looks at the local registry, so only this check catches
  a wrong local entry. register_local_model refuses new names of this form, since the form is reserved for packs.
- Effect on the outputs. _cellpose loads CellposeModel(pretrained_model=<blob path>). Cellpose 4.2.1.1 chooses the
  network from the weights' state dict, not from the file name, so this is the same load as the bare name; P3
  checks it once. The manifest's model_sha256 and the segment sub-hash then carry the real digest, not "absent".
- The cpsam drift. nuclei._cellpose maps the bare name cpsam to CellposeModel() with no weights, and Cellpose's
  default became cpsam_v2 in 4.2.0 (2026-06-11). So generic.yaml runs cpsam_v2 today while its manifests record
  "absent". Packs avoid the ambiguity by digest. generic.yaml is left alone (section 15a).
- Non-Cellpose model keys (segmentation.nuclei.nuclei3d.model, spots.net3d.model) keep resolving as today in T5.
- Verification record. provenance.file_sha256 caches a digest only for the life of a process, and every
  subprocess-runner call is a new process, so each card refresh would re-hash 1.2 GB per model. pack status
  therefore trusts a blob whose size and mtime match its <sha256>.json record (size, mtime, verified_at) and
  re-hashes only when either changes. The run path still hashes the file once per process and applies the pin
  check, so the record only speeds up status.

### 3.6 Cards and validation records

- The pack card (README.md) says what the pack is for and not for; lists recipes with their dims and channel
  naming; describes the models (training data, base model, licence) and how the presets were chosen; summarises
  validation with intervals; states limitations (training-data overlap, prism volumes in projected runs, anisotropy
  limits); and gives the citation. Every outside card carries the object_kind table: in the yeast pack one row of
  the "nuclei" table is a whole cell, because table and column names never change. When a recipe sets
  spots.sharpen, the card says that effect_size and effect_size_min are then on the sharpened (DoG) scale (section
  4.7). The card also tells labs that name channels after genes (MDN1, ACT1) to bind roles with the role editor
  (section 8.5), because no role token matches such names.
- A model card (models/<name>.card.json) has the fields train_nuclei writes today (name, what, model_file,
  model_sha256, model_size, base_model, date, germquant_version and one <tool>_version key per tool, training,
  orientations, zones, training_set, holdout_metrics, holdout_per_pair, licence_note) plus kind, dims, pixel_um,
  object_kind and source_url. Trained models also get system, pack {id, version, recipe} and, for public data,
  dataset {id, licence, split_sha256}. Every Cellpose-derived card states that the Cellpose training data include
  CC BY-NC sets. models_registry.validate_card checks a card against its spec using the keys train_nuclei actually
  writes (segment/training.py): it reads model_size, accepts size only as an alias in pack-written cards, and so
  accepts every existing card. model_file and base_model.path hold absolute machine paths; validate_card treats
  them as non-portable, and pack derive rewrites them (section 10).
- A validation record (validation/<recipe>.json) is written only by `germquant extbench report --final`. It holds
  the pack id, version and recipe, config_hash, model digests, engine version and git sha, constraints file,
  device, dataset id and licence, split file and its sha256, protocol file and its git blob sha, date, metrics
  with 95 percent cluster-bootstrap intervals and n, ceilings, comparators, and the gate rule with its verdict.

### 3.7 Versions, checksums and lint

- Engine: P1 sets 0.2.0.dev0 in pyproject.toml and __init__ (with an editable reinstall). v0.2.0 is tagged at M4.
- Packs use semver. Before 1.0.0, any change to hashed text or a model digest bumps the minor version. From 1.0.0
  on: MAJOR when any hashed text or model digest changes (numbers can move); MINOR for a new recipe, or a widened
  engine range with the goldens re-verified; PATCH for card, docs or metadata only. `germquant pack lint DIR
  --previous ZIP` enforces this against the previous release.
- An installed version is immutable. Install writes <packs_root>/<id>/<version>/ and refuses a different
  SHA256SUMS under the same version.
- SHA256SUMS holds `<sha256>  <relative/path>` lines, sorted, with forward slashes and LF, for every file except
  itself. Text files (.yaml, .yml, .json, .md, .csv, .txt, .cff) are hashed after CRLF to LF normalisation, all
  others as bytes. pack_sha256 is the sha256 of the normalised SHA256SUMS text. A root: checkout pack also records
  the normalised sha256 of each engine file it references.
- Git attributes. P1 adds `packs/** text=auto eol=lf` to .gitattributes, plus `binary` rules for `*.tif`, `*.tiff`,
  `*.png`, `*.npz`, `*.zip` and `*.h5` under packs/ and tests/fixtures/packs/. A plain `text` rule switches off
  git's binary detection: in a scratch repository, `packs/** text eol=lf` turned a 13-byte PNG holding CRLF pairs
  into an 11-byte blob.
- Git ignore rules. Today .gitignore ignores `models/` everywhere (line 52) and `*.tif`, `*.tiff`, `*.ome.tif` and
  `*.ome.tiff` (lines 15 to 18). `git check-ignore -v` reports `packs/<id>/models/*.card.json`,
  `packs/<id>/example/*.ome.tif` and the test pack's example as ignored, so a plain `git add packs/` would leave
  them out of every clone, and pack verify and the recipe goldens would fail there. P1 appends `!/packs/**` and
  `!/tests/fixtures/packs/**` at the end of .gitignore. A file-level negation such as `!packs/*/models/*.card.json`
  does not work, because git never re-includes a file whose parent folder is excluded (checked in the same scratch
  repository). The broad negations also re-include what .gitignore ignores everywhere else: `__pycache__/`,
  `*.py[cod]` and `.DS_Store` (lines 25, 26 and 47). In a scratch repository with the real .gitignore and the two
  negations, `git add -A` staged packs/y/__pycache__/junk.pyc and packs/y/.DS_Store. So P1 appends, after the
  negations, `/packs/**/__pycache__/`, `/packs/**/*.py[cod]`, `/packs/**/.DS_Store`,
  `/tests/fixtures/packs/**/__pycache__/`, `/tests/fixtures/packs/**/*.py[cod]` and
  `/tests/fixtures/packs/**/.DS_Store`. With them, the same test tracks the cards and examples and leaves those
  files out. build-sums skips such files too. A CI test checks that every path listed in every committed
  SHA256SUMS is tracked (`git ls-files --error-unmatch`).
- Lint rules:
  - ASCII only, and no em or en dashes;
  - no config/config.yaml inside a pack, because it would shadow a bare `extends: config.yaml`;
  - outside recipes set every worm stage switch explicitly (germline, envelope, audit, axis, staging, sc.trace,
    granule, coloc, partition, granule_tail);
  - a neutral metadata.filename_regex, '(?!x)x', which matches nothing, unless the pack declares its own;
  - segmentation.nuclei.method: cellpose (or the backend's own equivalent);
  - each recipe's expected_config_hash equals its computed hash, SHA256SUMS is current, and pack: lines are
    versioned;
  - every key an outside recipe chain sets, pack_base.yaml included, is one the engine reads. That is a key that
    config/config.yaml sets (read from the engine checkout, since outside chains no longer reach it), or a key on
    an allow-list kept in packs.lint (the code-default keys that config.yaml does not carry, and the keys T5
    adds). config.yaml also carries keys that nothing reads, so packs.lint keeps a deny-list of them that wins
    over config.yaml. spots.channel is one: config.yaml sets it (line 66), but the spots stage always counts the
    foci role (pipeline.py line 580). The rule refuses typos and deny-listed keys. It does not apply to the
    C. elegans pack, whose recipes are the frozen files themselves;
  - outside recipe chains reach no config/config.yaml: they extend the pack's frozen pack_base.yaml (section 3.12);
  - outside recipe chains switch the batch framing check off through qc.germline_outlier_factor, set qc.min_nuclei,
    and set expected_from_germ_cell: false on every spots instance (section 3.12);
  - no absolute or machine path in a recipe or a card.
- Why the neutral regex: the worm regex is matched case-insensitively, so it reads a file such as
  20261005_hela_hs_m_001 as genotype HELA, treatment heat, sex male, germ_cell spermatocyte (checked with
  parse_sample on this worktree).

### 3.8 Where packs live and how they are installed

packs.pack_roots() searches, in order:
1. <project>/packs/ while a project is open or running (derived packs, section 10);
2. the built-in packs, as <root>/packs/<id>/pack.yaml, where root is the engine checkout (a config root that holds
   src/germquant) or $GERMQUANT_CONFIG_ROOT;
3. each folder of $GERMQUANT_PACKS (an os.pathsep list; the first folder is the install target);
4. the per-user folder: %LOCALAPPDATA%/germquant/packs, ~/Library/Application Support/germquant/packs or
   ~/.local/share/germquant/packs.

Installed packs sit at <root>/<id>/<version>/pack.yaml. The first (id, version) found wins. A built-in id and
version cannot be shadowed: a duplicate elsewhere is listed as shadowed and ignored. Discovery is cached per
process, and packs.refresh() rescans. A broken pack is logged and listed as unusable; it can never break a run
that does not use it. `germquant pack install SRC [--replace]` takes a folder or a zip. It refuses absolute
paths, drive letters, `..` components, symlinks, a missing or wrong SHA256SUMS, and non-ASCII hashed files.

The built-in root. setuptools packages only src/, and provenance._repo_root returns the first parent folder that
holds .git, pixi.lock or pyproject.toml. In a non-editable install inside another project's environment
(napari-germquant/.venv, for example) it therefore returns that project's folder, and the built-in packs and the
generic.yaml fallback would vanish or resolve against the wrong tree. So a root without src/germquant is not taken
as the checkout, for pack discovery and for the bare-name extends fallback alike, and the working directory is
never a built-in root. When no built-in root is found, `pack list` says so and names two fixes: install the engine
as an editable checkout (the supported route for 0.2.0), or set $GERMQUANT_CONFIG_ROOT, as the container does (open
question O21).

Pinned versions outlive an engine update. A built-in pack has no version folder, so a git pull that moves
budding_yeast from 0.1.0 to 1.0.0 removes 0.1.0 from disk, and a project pinned to it would no longer resolve.
packs.vendor(id, version, dest) therefore copies a pack version through the install path (SHA256SUMS verified,
immutable). new_project uses it to copy the chosen version into <project>/packs/<id>/<version>/, and pack derive
does the same for its parent. A pack holds no weights, so the copy is small. A root: checkout pack copies only its
own files; its recipes still point at the frozen engine files, whose recorded digests it verifies. An identical
copy (the same pack_sha256) is not reported as shadowed.

### 3.9 Availability: pack status

packs.pack_status(pack, recipe=None, image=None, deep=False) returns usable, blockers (each with its fix) and
warnings. The shallow check imports nothing heavy:
- the engine version lies inside requires.germquant (packaging.specifiers; P1 declares packaging>=21);
- the recipe loads and its hash equals expected_config_hash;
- each needed backend name is declared (the built-ins plus entry-point names read through importlib.metadata,
  without loading them), and each needed module is found by importlib.util.find_spec (cellpose, torch, spotmax).
  The built-in cellpose name is always declared, so the module check is what greys a pack in an environment
  without Cellpose or torch;
- each model blob is present with the right size; its sha256 is verified when its size or mtime differs from the
  blob's verification record (section 3.5), and trusted otherwise;
- the readers extra is present when the input is not an .nd2;
- with an image header given, the pixel size and z step fall inside expects (a warning only).

`--deep` imports germquant.pipeline first (torch before CuPy), then loads the backends and reports the device
that will run, with a warning when the recipe recommends an accelerator and none is present. `germquant run
--pack` refuses an unusable pack (exit code 2, with the fix printed) unless --no-pack-check is given. Pack
recipes set segmentation.nuclei.method: cellpose, so a missing Cellpose fails loudly instead of falling back to
the watershed.

The shallow check, `pack list`, `pack fetch-models` and `pack fetch-example` are also reachable through
`python -m germquant.packs`, an entry that imports neither cli.py nor pipeline and so never loads torch.
`python -m germquant.packs --version` prints the engine version, git sha and python, which the plugin compares with
its own germquant (section 8.1). `germquant --version` prints the same through the full CLI, which loads torch
first, as today.

### 3.10 Running a pack, and what the outputs record

- api.resolve_config(config=None, profile=None, pack=None, recipe=None) takes exactly one source. With pack it
  loads packs.resolve_recipe(...) through load_config and sets cfg.pack_record. api.run gains pack=None and
  recipe=None.
- CLI: run and batch gain `--pack ID[@VER]` in the existing mutually exclusive group with --config and --profile,
  and `--recipe R` (an error without --pack; the pack's default recipe otherwise).
- process_image gets five guarded hunks:
  1. a read_fn=None keyword passed to _read_input, used only by the pack golden runner;
  2. when cfg.pack_record is set, run_geometry["pack"] = {id, version, recipe, source, pack_sha256, object_kind,
     dims, status} and run_geometry["runtime"] = {device requested, device selected, torch version,
     cuda_available, python, platform};
  3. right after model_sha is computed, the pin check of section 3.5 (skipped under the classical method);
  4. in the segment closure (P2a): when segmentation.nuclei.require_channel is true (code default false) and the
     segmentation role does not resolve, the image fails with a message naming the role, instead of segmenting
     channel 0 with the flag segment:no_<role>_channel_using_first_channel. The default branch is today's lines
     verbatim. Every outside recipe sets the key (section 3.12). It is a recipe key, not a test of
     cfg.pack_record, so `--pack` and `--config <the recipe file>` behave alike. It does not apply when
     segmentation.nuclei.labels_from is set, because nothing is segmented then;
  5. at the qc.qc_flags call (pipeline.py line 893; P2a): qc.min_nuclei (code default 10) sets the threshold of
     the existing qc:few_nuclei_<n> flag. qc.qc_flags gains min_nuclei=10, and process_image passes the key only
     when the config sets it, so the default call and every default flag string stay as they are.
- __stages.json of pack runs also gets post_run.segment_device_used: the device Cellpose actually ran on after any
  CPU retry. A context variable in germquant.device carries it, and cellpose_backend fills it.
- `germquant batch --pack` writes the same pack and runtime records into the batch-root run_manifest.json,
  because batch passes its manifest to process_image and no per-image manifest exists.
- The pack, runtime and device records never enter a table, the done marker or a stage hash: the recipe text and
  the model digest already do.

### 3.11 The C. elegans pack

packs/celegans_germline/pack.yaml uses root: checkout. Its recipes point at the existing files and never copy
them. A copy would change base_dir, so models/models and config/channel_maps would stop resolving even though
identical text hashes the same.
- Recipes: default (config/config.yaml), rad51_foci, segmentation_only, n2_sc, ccw77_partition and live_syp2
  (3d+t). bench_log_tuned, bench_spotmax_tuned and config_n2dryice.yaml (a fixation test, regression use only)
  are not recipes.
- Model: the release asset germline_nuclei_combined (kind release, sha256 e4e32320...8f57, 1,218,646,375 bytes),
  with a card written from the release record: a June 2026 Cellpose-SAM fine-tune and the CC BY-NC note. Before the
  card is written, P1 compares the fine-tuned state dict with both cached bases (cpsam_v2 and cpsam) tensor by
  tensor and records the method and the result. The file dates point to cpsam_v2: it was cached on 2026-06-12, the
  model is dated 2026-06-13, the training script builds CellposeModel(gpu=True), whose default in Cellpose 4.2.x is
  cpsam_v2, and the cpsam (v1) cache file dates from 2026-09-24. Dates do not prove it, so if the distances do not
  separate clearly the card says the base is not determined (section 15d).
- validation/ points at docs/SPOTMAX_VALIDATION.md, docs/BENCHMARK.md and the 23 real goldens. It never computes
  or reads a test-split score.
- Acceptance: `--pack celegans_germline` writes tables cell-identical to `--config config/config.yaml`, with
  config_hash 84bdda5a5230.

### 3.12 The two outside packs

Both packs extend their own frozen base, config/profiles/pack_base.yaml, instead of generic.yaml, which extends
config.yaml. P3 writes it once from the generic profile as resolved at the T5 base (yaml.safe_dump, ASCII, no
comments), with io.channel_map pointing at the pack's own map, no machine path, and the settings both recipes
share: every worm switch off, the neutral regex and the QC settings below. The writer drops the keys on the lint
deny-list, such as spots.channel, which the resolved profile inherits from config.yaml (section 3.7). It is hashed
with every recipe, so the planned config.yaml hash bump moves no outside pack hash, and a validated 1.0.0 pack
keeps its version, DOI and validation record across that release. Recipes name it by its bare name, which today's
extends rule already finds in the recipe's own folder, so it needs no new lookup.

cultured_cells (object_kind nucleus):
- nuclei_foci (default, dims projected): Cellpose-SAM on the DAPI maximum projection, an area gate, nuclei copied
  through z; foci counted in 3D by SpotMAX; a second foci channel (53BP1) as a spots instance. Runs on a CPU.
- nuclei_foci_3d (dims 3d): Cellpose in 3D with anisotropy, exactly as today; GPU recommended.
- nuclei_foci_2d (dims 2d): single-plane images; added in P8, after P2b has merged.

budding_yeast (object_kind cell):
- cells_spots (default, dims projected): Cellpose-SAM on one transmitted-light plane (sharpest or middle, chosen
  by tuning), whole-cell masks copied through z; spots counted per cell in 3D.
- cells_2d (dims 2d): single-plane images; added in P8, after P2b has merged.
- cells_2d_yeaz (dims 2d, optional): the germquant-yeaz backend; greyed until that extra is installed (P8).

Starting values, provisional until tuning (section 6.8):
- Yeast channel map, yeast_transmitted_foci.yaml:
  - transmitted: names phase, contr, PhC, Bright, DIC, trans; no index fallback; required;
  - foci: names smFISH, FISH, Quasar, Cy3, Cy5, mNeon, GFP, YFP, foci, spots; index 1; required;
  - dna: names DAPI, Hoechst; no index fallback; optional.
  These cover the S-BIAD3750 channel names of SUT509 (Bright and Cy3) and mtDNA (phase_contr and mNeon, the channel
  the truth refers to; mKate is a reference channel and binds no role), as read from the deposit. They do not cover
  MDN1: its smFISH channels are named MDN1 and ACT1, which no foci token matches, so ChannelMap.resolve would bind
  foci by index 1, which could be DAPI (a role collision, a QC fail) or the wrong probe. The spotmax_gt adapter
  therefore writes the MDN1 channel as smFISH_MDN1, the channel the truth refers to, and ACT1 under a name no token
  matches (section 6.5). MDN1 Position_4 has no phase_contr channel at all. The token BF is left out on purpose,
  because it would also match BFP. Without an index fallback, a stack with no transmitted channel raises
  channel:transmitted:MISSING_REQUIRED (a QC fail), and require_channel makes the image fail instead of segmenting
  channel 0.
- Yeast spots, seeded from the authors' SUT509 SpotMAX INI (read from the deposit): spot_radius_um 0.244 and
  spot_radius_z_um 1.0 (the INI's (1.0, 0.244, 0.244) um), gauss_sigma_um 0.054 (the INI's 0.75 px at
  0.072 um), threshold_triangle, effect_size_min 2.5, remove_hot_pixels true, sharpen true, merge_z_columns
  false. Segmentation: channel transmitted, project plane, project_plane sharpest, min_area_um2 3.0,
  require_channel true. effect_size_min 2.5 is on the sharpened scale, as in the authors' runs (section 4.7).
- Cultured channel map, cultured_dna_foci.yaml: dna (DAPI, Hoechst, DNA, SiR-DNA, DRAQ5, H2B; index 0;
  required); foci (gH2AX, gammaH2AX, H2AX, RAD51, foci, puncta; index 1); foci2 (53BP1, foci2; index 2;
  optional). Laser-line tokens are left out on purpose: such names bind by position with the existing advisory
  flag.
- Cultured spots: spot_radius_um 0.3, spot_radius_z_um 0.6, gauss_sigma_um 0.08, threshold_triangle,
  effect_size_min 3.0, merge_z_columns true, and spots.instances [{name: foci2, role: foci2,
  expected_from_germ_cell: false}]. Segmentation: channel dna, project max, min_area_um2 20.0, require_channel
  true.
- Both packs: method cellpose, cellpose_model cpsam_v2@0f1cc3f7, the neutral regex, every worm switch off. Both
  also set segmentation.nuclei.diameter_um to the system's typical object size (yeast 4, cultured nuclei 12).
  Only the classical segmenter reads it, and the recipe goldens run the classical segmenter.
- Both packs also set QC values for fields that are not gonads: qc.germline_outlier_factor and qc.min_nuclei in
  pack_base.yaml, and expected_from_germ_cell: false on every spots instance a recipe declares. Each spots instance
  otherwise writes the worm-only column <name>_expected_per_nucleus (pipeline.py line 666). With germline off every
  object counts as germline, so the batch framing check (1.8 times the median) would flag dense fields; a factor of
  1.0e+9 switches it off. batch and batch --resume read the factor (cli.py line 764). The plain
  `germquant collect` has no config and keeps 1.8 (cli.py line 915), so collect_project passes the recipe's
  qc.germline_outlier_factor to batch.collect. docs/PACKS.md tells labs to collect outside results with
  `germquant project collect` or `germquant batch --resume`. qc:few_nuclei_<n> fires below 10 nuclei today, and
  DeepFoci fields hold about 5, so qc.min_nuclei is 1 for cultured cells and 10 for yeast. These flags are
  advisory, not QC fails, but they would mark almost every cultured-cell image.

The yeast recipe, abridged:

```yaml
# packs/budding_yeast/config/profiles/cells_spots.yaml
extends: pack_base.yaml                 # the pack's frozen base (section 3.12), not generic.yaml
io: {channel_map: config/channel_maps/yeast_transmitted_foci.yaml}
metadata: {filename_regex: '(?!x)x'}    # matches nothing; design columns come from a sample sheet
segmentation:
  nuclei: {method: cellpose, channel: transmitted, cellpose_model: cpsam_v2@0f1cc3f7,
           project: plane, project_plane: sharpest, min_area_um2: 3.0, require_channel: true}
spots: {spot_radius_um: 0.244, spot_radius_z_um: 1.0, gauss_sigma_um: 0.054,
        thresholding_method: threshold_triangle, effect_size_min: 2.5, remove_hot_pixels: true,
        sharpen: true, merge_z_columns: false}
germline: {enabled: false}
envelope: {enabled: false}
audit: {enabled: false}
axis: {enabled: false}
staging: {enabled: false}
sc: {trace: {enabled: false}}
granule: {enabled: false}
coloc: {enabled: false}
partition: {enabled: false}
granule_tail: {enabled: false}
```

The recipe sets no spots.channel key: nothing reads one, because the spots stage always counts the foci role
(pipeline.py line 580). The QC values sit in pack_base.yaml.

### 3.13 Commands and API

- `germquant pack list [--all] [--json]`, `show ID[@VER] [--recipe R]`, `check ID[@VER] [--recipe R]
  [--image FILE] [--deep] [--json]`, `fetch-models ID [--recipe R] [--from-file NAME=PATH]`, `fetch-example ID
  [--out DIR]`, `install SRC [--replace]`, `verify ID|DIR`, `build-sums DIR`, `lint DIR [--previous ZIP]`,
  `golden ID [--tier recipe|model] [--update]`, and from P5 `derive` (section 10).
- `germquant run IMAGE --pack ID[@VER] [--recipe R]` and `germquant batch FOLDER --pack ...`.
- `germquant assemble --channels A.tif B.tif --names Bright Cy3 --spacing DZ DY DX --out X.ome.tif` (P3) stacks
  a one-TIFF-per-channel layout, such as Cell-ACDC's, into one OME-TIFF with channel names and voxel size.
- `python -m germquant` works through a new src/germquant/__main__.py that enters cli.main, so the torch-first
  import guard in cli.py still runs first.
- `python -m germquant.packs list|check|fetch-models|fetch-example|--version`: the torch-free entry of section 3.9,
  which the plugin's runner uses for listing and checks. `germquant --version` prints the engine version and git
  sha.
- New api names, appended to __all__: list_packs, find_pack, load_pack, pack_status, fetch_pack_models,
  fetch_example, install_pack, list_segment_backends, list_spot_backends (metadata only), known_models and
  register_local_model; resolve_config and run gain pack and recipe. Later packages add create_project,
  open_project, run_project, project_status, collect_project, import_labels, tune_spots and list_post_stages.
  Importing germquant.api loads torch, because api imports pipeline. Code that must stay light, such as a pack
  card refresh, imports germquant.packs, which offers the same list_packs, pack_status and fetch_pack_models.
- New command groups live in their own modules (cli_packs.py, cli_projects.py, cli_extbench.py), each with
  add_parsers(sub) and dispatch(args), so cli.main gains one line per group and merges stay small.

## 4. The 2D path

### 4.1 Two modes, both opt-in

- Planar: the image has one plane, and every quantity is 2D. The key io.planar turns it on.
- Projected: a 3D stack whose objects are segmented on one 2D image (a projection, or one plane of the
  segmentation channel) and copied through every plane; everything downstream stays 3D. The key
  segmentation.nuclei.project turns it on. This is how the yeast deposit itself was made: Cell-ACDC segmented one
  z-slice of the transmitted channel, and SpotMAX stacks 2D masks through z. Thin monolayer stacks need it too.

Both modes live in a new module, src/germquant/planar.py, and share one 2D segmentation function. Neither is
reachable from a frozen config, and every synthetic fixture stays 3D.

### 4.2 New keys, all inert by default

- io.planar: false (default) | auto | true.
  - auto: planar when the header says the file has one plane (sizes has no Z, or Z is 1) and the stack read has
    one plane. This one rule serves both readers; is_2d is not used, because it means "no Z axis" in the nd2
    reader and "Z is 1" in the bioio reader.
  - true: planar is required; more than one plane is refused before any stage runs, with a message that names
    --z-range and segmentation.nuclei.project.
  - Under auto, a 3D file cut to one plane with --z-range stays 3D.
- io.voxel_size_um: null | a number | [dy, dx] | [dz, dy, dx]. The file's own pixel size before xy_stride, for
  uncalibrated files; stride is applied after it. Recorded as run_geometry spacing_source "override" with the flag
  voxel:override, plus voxel:override_differs_from_file when the file's calibration is more than 1 percent away.
- io.input_suffixes: [".nd2"]. What batch discovers (P5).
- segmentation.nuclei.project: none | max | mean | plane.
- segmentation.nuclei.project_plane: mid | sharpest | an index. sharpest is the plane with the highest variance of
  the Laplacian of a Gaussian-smoothed plane (sigma 1 px, float32); ties go to the lowest z. An index out of range
  raises before segmentation.
- segmentation.nuclei.min_area_um2: 0.0 (no gate). The size gate for 2D objects: min pixels = area / (dy * dx).
- segmentation.nuclei.cellpose.object_diameter_um: null. When set, Cellpose gets diameter = d / dy pixels (dy after
  stride) in both the 2D and the 3D call, through the kw mapping that _cellpose already builds empty. Cellpose
  4.2.1.1 then rescales the image by 30 / diameter. When unset, both calls are exactly today's. It exists for
  objects outside Cellpose-SAM's native size range, such as 100x nuclei of 150 to 300 px. A micrometre key keeps a
  recipe valid across microscopes and xy_stride. The existing segmentation.nuclei.diameter_um keeps its meaning:
  Cellpose-SAM ignores it, and the classical watershed uses it for seed spacing.
- segmentation.nuclei.labels_from: null (P2c, section 4.11).
- spots.spot_radius_z_um, spots.remove_hot_pixels, spots.sharpen: section 4.7.

### 4.3 The planar decision and spacing

planar.decide(cfg, meta, stack) returns None unless io.planar is set, so every existing run takes the existing
path. In a planar run, dy and dx come from io.voxel_size_um, else from the new Stack.spacing_yx, else from
stack.spacing[1:] when spacing_ok, else (1, 1) with the new flag voxel:UNREADABLE_xy_assumed_1um. Array code
receives (PLANAR_DZ = 1.0, dy, dx). The placeholder reaches only code whose result does not depend on dz for one
plane. Planar runs write voxel_dz_um as NaN, raise the flag read:planar, and record planar: true and
spacing_source in run_geometry.

### 4.4 Gating: refuse, never skip silently

Stage gains planar_ok: bool = False as its last field, outside every hash payload. planar_ok is True for read,
acquisition, segment, measure, signal_qc, spots (and its instances), qc, render and write. germline and axis stay
False in T5: they are worm biology, every outside recipe switches them off, and axis weights by volume, which is
NA for 2D objects.

planar.check_planar(cfg, role_to_idx) runs right after the read, in planar and in projected runs. It raises
PlanarConfigError, a ValueError that is fatal like a failed read, when:
- a stage without planar_ok is switched on and its roles are bound;
- the segmentation backend lacks planar_ok (in both modes, since it receives a 2D image);
- in planar runs only, a spot backend of the main instance or of spots.instances lacks planar_ok. Projected runs
  count spots in 3D, so any spot backend may run there.

Stages whose roles are missing skip themselves as today. Recipes switch the 3D-only stages off explicitly.

### 4.5 Segmentation

- planar.segment_planar(img_1yx, spacing, *, method, cellpose_model, min_area_um2, device, backend, cfg_section)
  mirrors the control flow of segment_nuclei without touching it: get_backend, UnknownBackendError, auto falling
  back to classical, cellpose fatal, classical never calling the backend, and _check_labels.
- Backends get planar=True only when they declare planar_ok (the input_roles pattern). cellpose_backend gains
  planar=False. When True it calls a new nuclei._cellpose_2d(img2d, model, device_pref, diameter_px), which builds
  the model exactly as _cellpose does, runs model.eval(img2d, channel_axis=None) with no do_3D, z_axis or
  anisotropy, returns masks[None], and keeps the CPU retry. classical runs unchanged on the plane (its backend
  gains a planar keyword that it ignores); a probe on this worktree segmented 16 of 16 planted disks that way.
  cellpose and classical declare planar_ok; nuclei3d does not.
- The size gate is planar.drop_small_area, followed by a consecutive relabel.
- Projected mode: planar.project_image(img, mode, plane) builds or picks the 2D image, segment_planar runs on it,
  and the (1, Y, X) labels are repeated through Z. Only such runs get the flag segment:projected_<mode> and
  run_geometry segment_projection {mode, plane, z} (z is null for projections).
- render.make_montage gains base_plane=None. Projected runs pass the segmented plane for the segmented role, so the
  overlay sits on the image the model saw rather than on a maximum projection of a bright-field stack.

### 4.6 Measure

- Planar runs: planar.measure_planar calls measure_objects unchanged, then sets volume_um3 to NA and centroid_z_um
  to 0.0, and adds area_um2 (n_voxels x dy x dx) and equivalent_diameter_um.
- Projected runs: volume_um3 is NA, because a prism volume is not a cell volume. area_um2 and
  equivalent_diameter_um come from the 2D labels. Intensities are means over the 3D column. centroid_z_um stays
  the stack middle; it is documented and never NaN, so nothing downstream breaks.

### 4.7 Spots

- detect_spots gains planar=False and forwards planar=True only to backends that declare planar_ok.
  _spotmax_backend forwards it. detect_spots_spotmax(planar=False) gets one branch at the top: radii_px =
  (PLANAR_RZ_PX, r / dy, r / dx), gauss = sigma / min(dy, dx), features spacing (1.0, dy, dx), and the z merge off
  with the flag spots:z_merge_off_planar. The else branch is today's lines verbatim.
- PLANAR_RZ_PX = 1.0. SpotMAX 1.3.1 sets the z resolution to 1 px when SizeZ is 1 (core.py), and its peak finder
  already treats a (1, Y, X) image as 2D: it squeezes the image, takes the maximum of the footprint over z and
  keeps only the yx minimum distance. A test shows the tables for 0.5, 1.0 and 2.0. If they differ, 1.0 stays and
  the docs say so.
- spotmax declares planar_ok. log and net3d do not, and are refused in planar runs. Every planar z_um is 0.0.
- Optional SpotMAX keys for packs. They are read from a new SPOTMAX_OPTIONAL_KEYS tuple and forwarded only when
  present in the raw section. SPOTMAX_DEFAULTS is not extended, so a default run's backend section and its SpotMAX
  calls do not change; a spy test pins this.
  - spots.spot_radius_z_um: radii_px = (rz, r, r) / spacing instead of r / spacing.
  - spots.remove_hot_pixels: passed to preprocess_image as do_remove_hot_pixels.
  - spots.sharpen: mirrors SpotMAX core. The preprocessed image (smoothing and hot-pixel removal) stays `pre`. A
    sharpened copy (SpotMAX's DoG filter of `pre`) feeds spots_semantic_segmentation and spot_detection, and goes
    to spots_calc_features_and_filter as sharp_spots_image, while `pre` stays its image argument. SpotMAX then
    computes the spot-versus-background t-test and effect sizes on the sharpened image (spotmax/pipe.py lines 1236
    to 1248 in 1.3.1: "When comparing to the background we use the sharpened image"), and it uses the image itself
    when no sharpened image is given. So with sharpen on, germquant's effect_size column
    (spot_vs_backgr_effect_size_glass) is a DoG-scale value, exactly as in the authors' runs, and the same value
    that passing the DoG image as the image would give. Passing `pre` keeps the other features (preprocessed and
    local background intensities) on the unsharpened image, as SpotMAX core does. Consequences: effect_size_min is
    tuned together with sharpen, sharpen is part of the run-once group key (section 6.8), the RAD-51 value of 3.0
    is never carried across a sharpen change, and cards and docs/PACKS.md say which scale a recipe's threshold is
    on.
- A known difference from the authors' runs: germquant's chain does not run SpotMAX's spotFIT step. The authors'
  deposited SpotMAX outputs are therefore a comparator, not a parity target.

### 4.8 Writers and tables

- Planar runs call new writers in planar.py: write_labels_yx (int32 OME-TIFF, axes YX, PhysicalSizeX and Y only),
  write_spots_yx (a 2D disk per spot) and write_mask_yx. The 3D writers are not edited. tifffile.imread already
  returns (Y, X) for a (1, Y, X) file, so regression_diff and the plugin's _read_labels see 2D either way.
- schema.PLANAR_COLUMNS = {"nuclei": ["area_um2", "equivalent_diameter_um"]}, outside TABLES. _write_tables and
  _conform_schema gain extra_declared=None, filled only in planar and projected runs, so these columns follow the
  declared ones and never appear in a 3D CSV. There is no new TABLES entry (tests/test_stages.py forces every
  TABLES entry into every results folder).

### 4.9 Provenance and resume

- stage_hashes gains planar=False and labels_sha256=None. planar adds "planar": true to the read, segment,
  measure and spots payloads, in planar runs only. labels_sha256 enters the segment payload in labels runs only.
- The done marker gets the extra {"planar": true} in planar runs only.
- batch.is_done gains segment_hash_planar=None and sample_row_sha256=None. The marker's planar flag decides which
  segment hash to compare. cli._batch and timelapse.series compute both segment hashes once per batch. Today they
  compute one, so without this every planar image would be reprocessed under --resume. In labels runs the segment
  hash is computed per image.
- Default hashes, markers and resume decisions are unchanged, and P0 pins them.

### 4.10 Readers

Stack gains spacing_yx: tuple | None = None as its last field. bioio_reader adds _spacing_yx(img) and fills
spacing_yx (and a metadata key) when Z is 1. _spacing, spacing_ok and the (1, 1, 1) rule pinned in
tests/test_read_image.py stay as they are, so an existing Z=1 TIFF run keeps its numbers unless io.planar is set.
Today a calibrated 2D OME-TIFF loses its xy calibration: a probe on this worktree flagged voxel:UNREADABLE and
wrote volume_um3 equal to the pixel count. In nd2_reader.py, read_stack and read_nd2_metadata are untouched; the
Stack dataclass, which is defined there, gains only this field, with a default. `germquant info` prints "single
plane" when the header has one plane.

### 4.11 Labels from files (optional, P2c)

segmentation.nuclei.labels_from: null, or a path template with {image_dir}, {image_stem} and {image_id}.
- When set, the segment stage reads labels instead of segmenting: TIFF and OME-TIFF label images, and Cell-ACDC
  segm.npz (2D or 3D, frame 0 unless the template names one). All are read with numpy and tifffile; npz files are
  read with allow_pickle=False.
- Cellpose _seg.npy files are pickles, so the run path never reads them. `germquant import-labels --allow-pickle`
  converts them to TIFF once.
- Stride and z_range apply; 2D labels are copied through z; labels are relabelled consecutively; there is no size
  gate.
- segmentation_method reads labels_file; the flag segment:labels_from_file is raised; run_geometry gets
  labels_from {path, sha256}; the labels sha256 enters the segment sub-hash and the per-image resume hash. A
  missing file fails that image.
- Uses: detector-only scoring on public truth masks, and labs that already have masks.

### 4.12 Profile

config/profiles/generic_2d.yaml is new and extends generic.yaml: io.planar auto; io.input_glob "**/*" with
io.input_suffixes [.ome.tif, .ome.tiff, .tif, .tiff, .czi, .nd2]; min_area_um2 20.0; spots.merge_z_columns
false; granule and coloc off; the neutral filename regex. Its hash gets a P0-style pin and a synthetic golden.

### 4.13 Time-lapse

A synthetic Z=1 series goes through process_series and `germquant register` in a test. If the length-1 axis
breaks the phase correlation, planar series are refused with a message in T5 rather than fixed.

### 4.14 Still refused after T5

- Multi-scene and multi-position files. Projects refuse nd2 files with P > 1, because read_stack silently reads
  position 0 and stays untouched.
- log, net3d and nuclei3d in planar runs.
- Any 3D-only stage switched on, with its roles bound, in a run with 2D objects.

## 5. Backends

### 5.1 Contract additions, all optional (docs/BACKENDS.md)

- planar_ok and the planar=True keyword on both registries (P2a).
- declared_backend_names() on both registries: the built-ins plus entry-point names, without loading them (P1).
- models_registry.resolve_model_name(name, expect_sha256=None): the public helper a third-party backend calls for
  its own model key (P1).
- The post-stage contract on a new entry-point group, germquant.post_stages (P5b, section 5.6).
- The pack rule: a pack names backends and models, never callables, and its recipes set method: cellpose (or the
  backend's equivalent), so a missing piece fails loudly.

### 5.2 Cultured cells

The built-in cellpose backend with cpsam_v2@0f1cc3f7. cpsam v1 (sha256 e1440429...abe2, 1,233,587,898 bytes) is
scored on the tune split in P8 and becomes the default only if it wins on tune and val. Comparators, for the
benchmark only: StarDist 2D_versatile_fluo (TensorFlow environment; CPU only on native Windows), InstanSeg
fluorescence_nuclei_and_cells (torch, cloned environment), Spotiflow general and smfish_3d (cloned environment).

### 5.3 Budding yeast

The default is Cellpose-SAM (cpsam_v2@0f1cc3f7) on the transmitted role in projected mode. YeaZ images are in
Cellpose-SAM's training set (sampled at 1 to 2 percent), and the card says so.

### 5.4 Optional germquant-yeaz backend (P8)

plugins/germquant-yeaz/ is a separate distribution with its own pyproject and the entry point
germquant.segment_backends: yeaz.
- It calls yeaz.unet prediction, threshold and segment directly with germquant.device.select_device (Cell-ACDC's
  wrapper picks its own device), and copies Cell-ACDC's preprocessing (img_to_float, then CLAHE).
- It declares planar_ok, refuses non-planar input and imports torch lazily inside the callable.
- provenance() records the yeaz version, the weights' name and sha256 (through resolve_model_name), threshold and
  min_distance.
- Weights yeaz_budding_phc@<sha8> and yeaz_budding_bf@<sha8> are pack model specs that point at the source URLs
  Cell-ACDC lists. They are pinned at build and never redistributed.
- It never imports a Cell-ACDC model module: their installers prompt and pip-install on Enter.
- It is installed only in .venv_packs, with constraints/packs-extras-2026-10.txt, or with
  constraints/packs-extras-mac-arm64-2026-10.txt on a Mac. The plugin runs the cells_2d_yeaz recipe with the
  .venv_packs interpreter chosen for it (section 8.1).
- The recipe sets cellpose_model to the yeaz weights' name, so model_sha256 is the weights' digest.
- Parity with Cell-ACDC's YeaZ_v2 wrapper is checked on 5 images.

### 5.5 Not in T5

Segmentation input_roles; a bioimage.io backend; Cellpose 3 models (they need numpy below 2.1, and Cellpose 4
refuses their weights); Cell-ACDC's StarDist and cellpose_v3 wrappers anywhere near the batch venv; BABY,
YeastMate, deepBlink, cpdino and YeastSAM; a Spotiflow spot backend (only after parity on cv_detect); registry
names for nuclei3d and net3d. Cell-ACDC is never upgraded in the batch venv (1.8.x needs pandas below 3.0).

Some of this departs from the proposal, whose section 11 planned the yeast pack on Cell-ACDC's yeast models and
listed Cell-ACDC's model API and bioimage.io models as segmentation backends. T5 defaults yeast to Cellpose-SAM,
offers YeaZ as a direct optional backend, and leaves the Cell-ACDC model adapter and bioimage.io backends out. O19
asks the PI to accept the change.

### 5.6 Post-stages (P5b)

The proposal promises plugins in two senses: segmentation backends, which exist, and measurement steps that an
outside lab adds without editing the engine (proposal sections 1, 3, 9 and 11). Today the code has only the two
backend groups (segment/backends.py and spots/backends.py). T5 builds a minimal post-stage contract.
- A post-stage is a distribution with an entry point in the group germquant.post_stages. Its metadata, read through
  importlib.metadata without importing it, names the post-stage, its version, the channel roles and tables it
  reads, the tables it writes and its typed parameters, declared like recipe tunables (key, type, unit, range).
  pydantic is not added to the engine.
- A run switches post-stages on through a separate file, post_stages.yaml (a list of {name, params}). It sits in
  the project root next to recipe.yaml, ships in a pack as post_stages/<recipe>.yaml, or is passed with
  --post-stages FILE to run and batch (api.run and api.run_loaded gain post_stages=None). The project's file, or
  one given with --post-stages, replaces the pack's. load_config never reads it, so config_hash, the stage hashes,
  the done marker and is_done (batch.py line 119 compares config_hash first) do not change. Switching one on
  neither reruns the image nor sets validated: false. The file's digest enters the post-stage sub-hash in
  <image_id>__post.json.
- It runs after process_image has finished, on what it wrote: the result folder's tables and label image, plus the
  image read through the same reader with the run's stride and z range. api.run_loaded (which covers germquant
  run, api.run and run_project) and cli._batch call it only when post-stages are set; process_image is not edited.
  In T5 a time-series run refuses post-stages with a message, both in the series loop (process_series calls
  pipeline.process_image directly) and in run_loaded given t.
- It writes only new tables, <image_id>__post_<name>_<table>.csv, with the columns it declared. It never edits an
  existing table, a label image or the done marker. collect gathers them into batch_post_<name>_<table>.csv.
- Its record is <image_id>__post.json: version, parameters, its own sub-hash (its configuration, the digests of its
  inputs and its version) and a done flag. The main done marker and is_done do not change; a changed post-stage
  sub-hash reruns only the post-stage.
- A failure is recorded there and never fails the image.
- Like backends, post-stages install only into .venv_packs or a lab's own environment, never into the validated
  .venv (R6). A pack may need one by name (needs.post_stages), and pack status checks it from metadata.
- docs/POST_STAGES.md describes the contract with a worked example. A fixture post-stage, registered in process as
  the backend tests register backends (register_backend exists in both registries), covers it in CI.

## 6. Public datasets and benchmark manifests

### 6.1 Datasets used in T5

Licences and contents were checked on 2026-10-11 against the archive pages and APIs; "verified" below means that
check. Each dataset.json records the access date, and licences are checked again at download time.

1. SpotMAX ground truth, BioImage Archive S-BIAD3750 (Padovani et al. 2026, Science Advances). Verified, CC0.
   - https://doi.org/10.6019/S-BIAD3750
   - files: https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/750/S-BIAD3750/Files/
   - Yeast smFISH MDN1 (11 positions; 700 x 700 x 25 at 0.072 um and 0.29 um; 2D cell masks) and SUT509 (30
     positions; 15 planes at 0.24 um), two annotators per position. mtDNA nucleoids (189 positions; 43 planes at
     0.35 um; one annotator). C. elegans COSA-1 pachytene (19 positions; three annotators; a stretch goal).
   - Cell-ACDC masks, manual 3D spot centres as CSV, the authors' SpotMAX and RS-FISH outputs and INI files.
     The yeast part is about 6.25 GB.
   - Role: the primary yeast spot benchmark. Yeast segmentation is scored as agreement with curated YeaZ masks.
   - MDN1 Position_4 holds the ACT1, DAPI and MDN1 files but no phase_contr, so it cannot run end to end. It is
     scored detector-only (truth masks through labels_from), and dataset.json and the protocol count 10 end-to-end
     MDN1 positions. Its channel:transmitted:MISSING_REQUIRED flag is expected and recorded, not treated as a
     failure.
2. DeepFoci (Vicar et al. 2021). Verified, CC BY 4.0.
   - https://zenodo.org/records/4067741 (part 1), https://zenodo.org/records/5549971 (part 2),
     https://zenodo.org/records/2576241 (part 3, an older training release).
   - Cultured human cells (NHDF fibroblasts, U-87, head and neck primary cultures) stained for DAPI, gamma-H2AX
     and 53BP1; 1392 x 1040 x 50 voxels at about 0.065 x 0.065 x 0.3 um.
   - Manual 3D nucleus masks (237/10/30 FOVs) and manual foci centroids (239/60/100 FOVs). The test labels were
     clicked on 2D projections. A two-expert validation set gives a median inter-expert Dice of 0.75 at 1.95 um.
   - Part 1 is 122.4 GB. T5 uses the foci validation subset (10.2 GB) and part of the nucleus set, fetched as
     single zip members where possible.
   - Role: the primary cultured-cell benchmark for nuclei and foci.
3. BBBC039, U2OS Hoechst nuclei. Verified, CC0. https://bbbc.broadinstitute.org/BBBC039
   - 2D, 200 FOVs, about 23,000 hand-annotated nuclei with splits; about 81 MB.
   - Role: 2D nuclei (needs P2b). It may sit inside Cellpose's training data; each protocol checks.
4. BBBC033, mouse trophoblast stem cells in 3D. Verified, CC0. https://bbbc.broadinstitute.org/BBBC033
   - One stack with manual 3D instance labels; about 154 MB.
   - Role: a quick 3D sanity check and a model-golden crop.
5. BABY (Pietsch et al.). Verified, CC BY 4.0. https://doi.org/10.7488/ds/3427 (also
   https://datashare.ed.ac.uk/handle/10283/4386)
   - Hand-curated yeast outlines on 5-plane bright-field stacks in microfluidic traps, with train, val and test
     groups and tracks.
   - Role: yeast segmentation independent of YeaZ. The projected cells_spots recipe reads the 5-plane stacks as
     they are and segments one plane, so this score needs no P2b. It scores segmentation only: the stacks have no
     foci channel, so the spots stage skips and channel:foci:MISSING_REQUIRED is expected. cells_2d, a single-plane
     recipe, cannot read 5-plane stacks; it is scored on single planes that the baby adapter writes (needs P2b).
6. Spotiflow benchmarking sets (Dominguez Mantes et al. 2025). Verified; BSD-3-Clause as listed on the record.
   - https://zenodo.org/records/14914082 (spotiflow_datasets.zip, 697 MB)
   - 2D point annotations: telomeres, TERRA, DNA FISH, smFISH, MERFISH, HybISS and synthetic sets.
   - Role: a spot benchmark independent of SpotMAX (needs P2b).

### 6.2 Listed for later

Verified unless noted.
- BBBC038 (CC0), https://bbbc.broadinstitute.org/BBBC038: a heterogeneous 2D nuclei stress test.
- Yeast Image Toolkit (CC BY 4.0), https://yeast-image-toolkit.org/: yeast tracking points with scoring code.
- idr0063 (CC BY 4.0), https://doi.org/10.17867/10000159: yeast meiotic loci in 3D plus time; C.K. Cahoon is an
  author.
- BBBC035 (CC BY 3.0), https://bbbc.broadinstitute.org/BBBC035: synthetic 3D plus time nuclei.
- LMRG FISH (CC BY 4.0), https://zenodo.org/records/6560910: simulated 3D spots with exact truth.
- deepBlink (MIT), https://figshare.com/articles/dataset/Datasets/12958037, and U-FISH (MIT),
  https://huggingface.co/datasets/GangCaoLab/FISH_spots: 2D spot sets that overlap the Spotiflow sets.
- Unlabelled DNA-repair images for an in-house annotation set: S-BIAD3475, RAD51 in mouse ES cells (CC BY 4.0,
  https://doi.org/10.6019/S-BIAD3475); S-BIAD1494, gamma-H2AX and 53BP1 in 3D (CC BY 4.0,
  https://doi.org/10.6019/S-BIAD1494); S-BIAD679, an NBS1-GFP time-lapse (CC0,
  https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD679).
- idr0062 NesSys (CC BY 4.0), https://doi.org/10.17867/10000125: lamin B1 plus DAPI, few cultured samples.
- idr0107, Arabidopsis pachytene 3D-SIM with traced synaptonemal complexes (CC BY 4.0 per the dataset pass; no
  URL was recorded, so not verified here): a possible external test of the SC tracer.
- Kromp et al., S-BIAD634 (CC0), https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD634: mostly tissue.
- idr0047 yeast smFISH (CC BY 4.0), https://idr.openmicroscopy.org/study/idr0047/: algorithmic truth, so for scale
  tests only.
- S-BIAD867 Rad52-YFP screen (CC BY-NC-SA 3.0), https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD867:
  per-strain counts only, and NC-SA limits sharing derived labels.

### 6.3 Excluded

- YeaZ training sets: no licence is stated (https://www.epfl.ch/labs/lpbs/data-and-software/). They are also not
  independent of YeaZ-based backends.
- Cell Tracking Challenge: custom terms; publication outside the challenge needs written permission. Third-party
  Zenodo re-uploads labelled CC BY 4.0 conflict with those terms.
- DynamicNuclearNet and SpotNet (non-commercial), LIVECell (CC BY-NC 4.0), Allen Cell hiPSC (custom
  non-commercial terms), FociRad (no licence file) and the YeastMate data on OSF (no licence).
- LFCT low-frame-rate tracking: not verified; its DOI, https://doi.org/10.5281/zenodo.20186357, returned 404 on
  2026-10-11.
- NuInsSeg: H&E histology, not a fit.
- HARLEY yeast multi-annotator foci: no licence stated; not verified.

### 6.4 What is committed per benchmark

benchmarks/ext/<dataset>/ holds no machine paths:
- dataset.json: id, title, DOI and URL, licence, citation, access date, adapter and its version, truth kinds, the
  files used, the channel names written, and the voxel sizes with their source.
- files.tsv: relative path, size and sha256. BioStudies FTP publishes no checksums, so a maintainer records them
  once and every later fetch is verified against them.
- split.json: tune, val and test ids, assigned by whole position, FOV or source group with a seed, committed with
  its sha256 before any tuning. Published splits are used where they exist (DeepFoci, BBBC039, BABY).
- protocol.md: matching tolerances, primary and secondary metrics, the tuning grid and budget, comparators,
  contamination checks, ablations, and the pre-registered minimums for "validated". It is frozen before the first
  test score.
- final_runs.jsonl: append-only, one line per final test scoring.
- smoke_expected.json: the pinned metrics of the smoke tier.

The data root comes from $GERMQUANT_EXTBENCH_DATA or --data-root, and is recorded only in a gitignored local.json.

### 6.5 Adapters and the truth contract

Each adapter writes, under the data root:
- converted/<item>.ome.tif: CZYX or CYX, with PhysicalSize X, Y and Z and channel names that the pack maps bind by
  name;
- <item>__truth_labels.tif with a _meta.json (complete or partial, source, annotator, overlap rule);
- <item>__truth_points_<annotator>.csv with point_id, z_um, y_um, x_um, zi, yi, xi, object_id and projected;
- converted.json with the output sha256s and pixels_sha256.

Per dataset:
- spotmax_gt stacks the per-channel TIFFs, reads the metadata.csv voxel sizes, converts segm.npz
  (allow_pickle=False; MDN1 metadata says the masks are 2D, mtDNA metadata says 3D), and reads every annotator's
  truth CSV. Their layouts differ. dimitra_MDN1_gt_manual.csv holds frame_i, Cell_ID, z, y, x and id in integer
  voxels. The mtDNA mNeon3Dspot_segm_gt_manual.csv holds spot_id, float z, y and x, is_in_segm_mask,
  PhysicalSizeZ, PhysicalSizeY, PhysicalSizeX and Cell_ID. The leo_manual_coords files (MDN1 and SUT509) and the
  SUT509 yagya_gt_manual.csv files have not been read yet. The adapter maps each layout to the truth-point format.
  Where segmInfo.csv exists (the 11 MDN1 and 30 SUT509 positions; none of the 189 mtDNA positions has one, and
  their masks are 3D), it records the plane that Cell-ACDC segmented, for an oracle-plane ablation. It writes the
  MDN1 channel as smFISH_MDN1, the channel the truth refers to, and ACT1 under a name no role token matches, so
  foci binds MDN1 by name; dataset.json records the renaming. It marks MDN1 Position_4, which has no phase_contr,
  as detector-only.
- deepfoci declares the XYCZ order and reorders with tifffile instead of letting bioio guess, flags the
  projection-clicked test labels as projected, and keeps the two validation experts apart.
- bbbc039, bbbc033, baby and spotiflow complete the set. baby flattens overlapping mother and bud outlines, smaller
  wins, and reports the affected pixel fraction. It names the bright-field channel so that the transmitted role
  binds it by name, writes each stack twice, as it is (5 planes, for cells_spots) and as one plane for cells_2d
  (the middle plane unless the annotations name one), and records the chosen plane in dataset.json and protocol.md.

CI tests every adapter on tiny synthetic stand-ins of its layout. Nothing is downloaded in CI.

### 6.6 Splits, protocols and the final-score gate

- Tuning reads tune, selection reads val, and test is scored once per frozen pack version.
- `germquant extbench score --split test` refuses unless all of these hold: --final is given; the git tree is
  clean; protocol.md is committed (its blob sha goes into the record); the pack verifies against its SHA256SUMS;
  and no earlier final line exists for the same pack id, version and recipe. A second final run needs
  --rescore-reason, which is logged.
- The gate numbers go to the PI before any test score (section 14). Proposed:
  - yeast: cell F1 at IoU 0.5 of at least 0.85 against the curated masks; spot F1 within 0.05 of the
    inter-annotator F1; per-cell count CCC of at least 0.80;
  - cultured cells: nucleus F1 at IoU 0.5 of at least 0.85; foci F1 at 1.95 um of at least 0.64 (DeepFoci's own
    lower bound).
- import-labels refuses images whose file sha256 or pixels_sha256 appears in the test partition of any
  committed split.

### 6.7 Scoring: one code path for every system

- Objects: germquant.benchmark.nuclei3d overlap_table and match_instances on (Z, Y, X) arrays, with 2D items as
  (1, Y, X). F1 at IoU 0.5 is primary; AP over IoU 0.5 to 0.95, matched IoU, splits and merges are secondary.
  Projected runs are scored on one plane of their labels against 2D truth. DeepFoci 3D truth is max-projected for
  the projected recipe and scored in 3D for nuclei_foci_3d.
- Spots: germquant.benchmark.spots match_spots (one to one, Hungarian assignment on the sparse graph, inside an
  anisotropic cylinder) and label_at_points; projected truth is matched in 2D. Tolerances per protocol: DeepFoci's
  published 1.95 um on projections, plus a strict 3D score (0.5 um in xy, 1.0 um in z) on a sealed 3D-labelled
  holdout of its training FOVs; S-BIAD3750 at the SpotMAX paper's own rule if it states one, else 0.3 um in xy and
  0.6 um in z for smFISH and 0.3 and 0.7 um for mtDNA.
- Per-object counts: Lin's CCC (primary), Pearson, bias and MAE from germquant.benchmark.stats, with 95 percent
  intervals from a cluster bootstrap over the item group (1000 resamples, fixed seed).
- Modes: end_to_end, and detector_only (truth masks through segmentation.nuclei.labels_from).
- Ceilings: annotator against annotator under the same scorer (the S-BIAD3750 pairs; the DeepFoci validation
  experts).
- Self-check before any germquant number is quoted: the harness must reproduce DeepFoci's published inter-expert
  median Dice (about 0.75 at 1.95 um) within the protocol's tolerance.
- Comparators, scored by the same code: the authors' SpotMAX and RS-FISH outputs shipped in S-BIAD3750; SpotMAX
  with the authors' INI values; Spotiflow, StarDist and InstanSeg zero-shot; YeaZ for yeast cells.
- Contamination checks in every protocol before a stock score is quoted: Cellpose-SAM's training list (YeaZ
  images, LIVECell, TissueNet and the Cellpose sets; BBBC039 to be checked), the S-BIAD3750 masks being
  hand-corrected YeaZ output, and SpotMAX's AI being trained on annotations from the same groups.
- The seal: extbench and tune may import only germquant.benchmark.stats, matching, nuclei3d and spots (manifest
  comes in through the package) and germquant.validate.compare. An AST test enforces the allow-list and the
  absence of the string rad51_split. The shared scorer ties public numbers to those modules: T5 never edits them,
  and the smoke tier catches drift if the RAD-51 work changes them.

### 6.8 Spot preset tuning

- tune.effect_filter re-implements the effect-size filter. A test proves it equal to
  benchmark.tuning.filter_effect_size on synthetic data; the test imports both, the module does not, and
  benchmark/tuning.py is not edited.
- tune.spots.tune_spotmax(items, grid, objective, split) runs SpotMAX once per parameter group that differs only in
  effect_size_min (at 0), filters per value, scores point F1 and per-object CCC, ranks by a fixed objective that is
  recorded in the output, and writes trials.jsonl, ranked.csv and best.yaml.
- sharpen is part of the group key like every other parameter. Because the effect size is computed on the sharpened
  image when sharpen is on (section 4.7), effect_size_min is searched separately for each sharpen value, and a
  threshold is never carried from one sharpen value to the other.
- The same code serves benchmark tune splits (truth points) and a user's napari points (section 10). The budget
  comes from protocol.md.
- The pre-registered z-merge ablation runs here: merge_z_columns on and off, on DeepFoci and on yeast smFISH.

### 6.9 Golden tiers for packs

- R, recipe tier, in CI. Each recipe runs on its pack's synthetic example with method forced to classical by
  mutating cfg._data, as tests/test_golden.py does, so the sentinel is the recipe's true hash. Tables are compared
  cell for cell with example/expected/recipe/<recipe>/. Classical runs never resolve a pack model (section 3.5), so
  the tier runs with an empty models_dir, and a CI test checks exactly that. CI has no bioio, so
  packs_golden.read_example reads the OME-TIFF with tifffile through process_image's read_fn seam; where bioio
  exists, a test shows that both readers give the same Stack.
- M, model tier, on the RTX 5090. A small public crop fetched by recipe (BBBC039, BBBC033, one S-BIAD3750
  position, one DeepFoci crop) and never committed runs as shipped (Cellpose on CUDA in bfloat16). It must equal
  example/expected/model_cuda/ exactly, twice in a row. env.json records torch, Cellpose, SpotMAX and the device. A
  Mac run is compared with a stated tolerance (object F1 at IoU 0.5 of at least 0.98 against the CUDA labels, spot
  counts within 2 percent) and reported, not gated.
- S, smoke tier, on the workstation (marker model). Three to five images per dataset with metrics pinned to three
  decimals: the drift alarm for models, the environment and the shared scorer.
- Examples are pinned by pixels_sha256, not by file sha256, because tifffile writes a fresh UUID into every
  OME-XML header. No public image is redistributed, and committed real-example outputs carry attribution.

### 6.10 Biases stated next to the numbers

S-BIAD3750 is SpotMAX's own benchmark, SpotMAX's AI was trained on annotations from the same groups, and the yeast
presets start from the authors' INI. germquant is therefore framed as reaching SpotMAX-level accuracy inside a
framework, not as beating SpotMAX. DeepFoci test labels were clicked on projections and matched with a lenient
radius. BABY cells are trapped in microfluidics. The generality claim rests on DeepFoci, BABY and Spotiflow.

## 7. The project file

### 7.1 Layout

```
myproject/
  germquant_project.yaml   manifest (not hashed, ASCII)
  recipe.yaml              hashed; exists only after an edit
  channel_map.yaml         hashed with the recipe; exists only after roles were edited
  post_stages.yaml         not hashed; exists only after a post-stage was switched on (section 5.6)
  samples.csv              optional design sheet
  runs/                    one log per run, runs/<run_id>.jsonl, merged on read
  images/                  optional; images may live anywhere (absolute or UNC paths)
  results/<image_id>/      process_image output, one run_manifest.json per image
  annotations/             <image_id>__points_<annotator>.csv, in the truth-point format
  training/<round>/        pairs, corrections.json, models with their cards
  tuning/<run>/            trials.jsonl, ranked.csv, best.yaml
  packs/                   the pinned pack version (copied at creation) and any derived packs
```

recipe.yaml sits in the project root, so _base_dir_for anchors channel_map.yaml and training/ paths to the
project.

### 7.2 Manifest fields (format 1)

- name; id (uuid4); created and saved stamps {at, germquant version, git sha, plugin version};
- pack {id, version, recipe, sha256};
- recipe {file or null, expected_config_hash, overrides (only --no-<stage> keys)};
- run {xy_stride, z_range, device, runner, interpreters (per recipe: python path, engine version, git sha)};
- models [{key, ref, sha256, card, card_sha256, source}];
- images [{id, path, size, mtime, sha256 (optional), timepoints (optional), excluded, reason}];
- samples (a file name or null);
- validated (computed, with the record that justifies it).

The manifest is written atomically (temporary file, fsync, os.replace) through fsutil.long_path, as ASCII with
escapes (yaml.safe_dump with allow_unicode False). Only the interactive commands write it (new, add, save,
write_recipe and use-model), under a lock file, germquant_project.lock, made with an exclusive create and holding
host, pid and time. A lock file behaves the same on Windows and Linux, unlike the byte-level locks section 16
rejected. A lock older than a set age is reported, never broken silently. `project run`, its shards and Snakemake
jobs never write the manifest.

### 7.3 Recipe edits and the hash

- An unmodified project has no recipe.yaml. project_config loads the pack recipe file itself, so the project's
  config_hash equals `--pack` and the pack's expected_config_hash, and it is comparable across every lab that uses
  that pack version. (A recipe.yaml holding only an extends line would hash differently: the base text, then NUL
  extends NUL, then its own text.)
- Any edit (parameter form, role editor, model choice) writes recipe.yaml with
  `extends: pack:<id>@<version>/<recipe>` plus the edits. That gives a new hash and validated: false.
- Switching a post-stage on is not such an edit. It writes post_stages.yaml, which changes neither the hash nor
  validated (section 5.6).
- Channel roles: the role editor writes channel_map.yaml, derived from the pack map with the image's exact channel
  names as tokens (index fallback and required flags kept), and points recipe.yaml at it. That is a new hash, and
  correctly so: it changes which channel is analysed.
- Per-image bindings outside the hash are not offered. is_done would not see them, because the segment sub-hash
  carries no channel roles, so a rebinding would keep stale results.
- On open, the hash is recomputed. A mismatch with expected_config_hash raises a warning that names the pack
  version change.

### 7.4 Sample sheet

- samples.csv has image_id (or file) plus any of genotype, sex, germ_cell, treatment, replicate and
  acquisition_date. Extra columns are allowed.
- process_image and api.run_loaded gain sample=None: the design fields for that image, which override
  parse_sample's. The sheet is applied per image, so every per-image table carries the design columns.
- The sheet is not hashed, just like the file name it replaces, so adding a row never invalidates finished images.
  run_geometry gets sample_sheet {source, row_sha256} and the done marker gets the extra sample_row_sha256, both in
  sheet runs only. is_done compares the row sha256 when given, so a corrected genotype reruns that image only.
- An image missing from the sheet falls back to the filename regex and raises the flag
  metadata:not_in_sample_sheet, in sheet runs only. Extra columns are joined into batch_<table>.csv at collect.
- CLI: `germquant run IMAGE --samples FILE` and `germquant batch FOLDER --samples FILE`. Projects use samples.csv.

### 7.5 Project models

A model trained in the project lives under training/<round>/models/<name>/, and recipe.yaml names it by relative
path. resolve_model_path already resolves a relative file against cfg.base_dir, so no new lookup is needed. It
tries $GERMQUANT_MODELS first, though, for every relative path (pipeline.py line 157), not only for paths that
start with models/, so a file at training/<round>/models/<name>/<file> under a shared $GERMQUANT_MODELS would win
over the project's own file. The manifest records the sha256 and the card. run_project compares the resolved file's
sha256 with the manifest's models[].sha256 before any pixel is read and refuses on a mismatch, and status reports
the same mismatch. `germquant pack derive` promotes a model into a derived pack with a pinned name (section 10).

### 7.6 Run logs

Each run writes its own log, runs/<run_id>.jsonl, so no two processes ever append to one file. That matters for
`project run --shard` and Snakemake jobs, which run at once, often on SMB or NFS shares where appends from several
hosts are not atomic. Readers (project_status, the plugin and collect) merge the logs by time. Each log holds one
fsynced JSON object per event:
- run_start {run_id, via (gui-subprocess, gui-inprocess, cli, snakemake), host, os, interpreter, engine version and
  git sha, python, constraints file, device requested and selected, config_hash, images};
- image_done, image_failed or image_skipped {run_id, image_id, out_dir, status, sha256 of __done.json,
  stage_hashes, failed_stages, qc_pass, segment_device_used, elapsed_s, error};
- run_end {run_id, counts}.

The plugin's subprocess runner tails the log of its own run for progress.

### 7.7 Status and portability

Status values are done, stale (with the reason is_done gives: config_hash, stages, geometry, planar flag, segment
hash or sample row), failed and missing. Status also reports the installed against the pinned pack version and
whether the project's copy of it is present, whether models are present and verified (a resolved file whose digest
differs from the manifest is reported), and an environment-mix warning when runs name more than one device,
interpreter or engine sha. A project zipped on the workstation and opened on a Mac recomputes the same config_hash,
resolves relative image paths, and keeps UNC paths verbatim; those fail with a clear message without the NAS or
VPN.

### 7.8 Commands, API, batch and Snakemake

- germquant.project: Project, new_project, open_project, save_project, add_images, project_config, write_recipe,
  run_project(project, images=None, resume=True, device=None, shard=None, progress=None), project_status,
  collect_project and use_model, re-exported from api. add_images checks headers and refuses multi-scene files,
  nd2 files with P > 1 and duplicate image ids (the id is the stem and names the results folder). new_project
  copies the pinned pack version into the project's packs/ (section 3.8). collect_project passes the recipe's
  qc.germline_outlier_factor to batch.collect (section 3.12).
- CLI: `germquant project new DIR --pack ID[@VER] [--recipe R] [--images GLOB ...]`, `add`, `run DIR [--images ID
  ...] [--no-resume] [--device D] [--shard K/N]`, `status`, `collect`, `use-model` and `export-commands`. --shard
  gives a SLURM array a deterministic slice. Shards write only their own run logs and result folders.
- run_project calls api.run_loaded per image into results/<image_id>/ (process_series for images with
  timepoints). It sets GERMQUANT_DEVICE from run.device for the call and restores it afterwards.
  project_config sets cfg.pack_record from the manifest's pack entry, noting whether the recipe was edited, so
  project runs carry the pack and runtime records like `--pack` runs.
- Batch over outside formats: io.input_suffixes (default [".nd2"]) goes to discover_files(suffixes=...).
  cli._batch keys out_dir, resume and collect ids on io.formats.image_stem, which equals Path.stem for .nd2 and is
  correct for x.ome.tif. Discovery skips germquant's own outputs: files with its output suffixes, and anything
  under the output root. The README's ".nd2 only" line changes.
- workflow/Snakefile accepts project: DIR (one `germquant project run DIR --images <id>` per job) or
  germquant_pack: ID@VER with recipe:, and uses io.input_suffixes and image_stem.

## 8. Plugin changes (napari-germquant)

### 8.1 Runners: run where the numbers are made

On the workstation the plugin venv has the CPU torch wheel, no SpotMAX, and an editable germquant from
C:/Users/ryane/germquant_wt_dev. The engine .venv has CUDA torch, SpotMAX and an editable germquant from
C:/Users/ryane/C_elegans_germline_masking_and_spot_counting. An in-process GUI run of an outside pack in the
plugin venv would fail its spots stage and segment on the CPU. So:
- A new module, _runner.py, holds EngineRunner with two implementations. SubprocessRunner runs whenever an engine
  interpreter is configured (stored per machine in QSettings). Detection at <checkout>/.venv/Scripts/python.exe or
  <checkout>/.venv/bin/python is only a suggestion: on this workstation the plugin's germquant is an editable
  install from C:/Users/ryane/germquant_wt_dev, which has no .venv, while the engine .venv lives in
  C:/Users/ryane/C_elegans_germline_masking_and_spot_counting, so detection would find nothing. On first run the
  dock therefore asks for the engine interpreter, offering any it detected. InProcessRunner, which calls the api in
  a thread worker, runs only when the user explicitly picks this environment. The dock always says which
  environment and device run.
- The runner checks the chosen interpreter with `<python> -m germquant.packs --version` (engine version, git sha,
  python) and warns when they differ from the plugin's own germquant.
- A recipe may have its own interpreter, chosen from the needs that pack status reports: cells_2d_yeaz runs in
  .venv_packs. The choice is stored per machine and recorded in the project (run.interpreters) and in every
  run_start.
- SubprocessRunner runs commands with list arguments and no shell. Listing and checks go through the torch-free
  entry: `<python> -m germquant.packs list --json`, `check ID --json`, `fetch-models ID` and `fetch-example ID`.
  Runs go through `<python> -m germquant run IMAGE --pack ... --out DIR` and, from M2,
  `project run DIR --images ...`, tailing the run's own log.
- Availability, models and examples therefore reflect the environment that will run. Both environments use one
  model store when $GERMQUANT_MODELS points both at the same folder, which the docs recommend.

### 8.2 Run dock: the pack picker (M1, P4)

- A System pack combo with its recipes sits above the existing profile combo, which moves unchanged into an
  "Advanced: profile or config file" section. The bench profiles are hidden there unless "show all" is ticked.
- A card shows the title, system, object kind, dims, version and engine range; the models with their status and a
  Get models button; the backends, available or missing with the extra to install; the validation status with
  dataset and metric, or "not validated"; the licence and citation; and a Get example button.
- An unusable pack or recipe is greyed with its reason and its fix, and Run is disabled.
- _engine.Recipe gains the (pack, recipe) form: cli_args() gives --pack ID@VER --recipe R, and load() calls
  api.resolve_config(pack=..., recipe=...). Export writes --pack commands and an overlay that extends
  pack:<id>@<version>/<recipe>, in ASCII, with no absolute path.

### 8.3 Results (M1, P4)

- The Labels layer is named "<id> cells" or "<id> nuclei" from the manifest's pack.object_kind.
- A summary line: N cells, M spots, mean spots per cell.
- A Provenance box: pack id and version, recipe, config_hash, model sha256, git sha and the device used.
- A banner when image_summary.segmentation_method differs from the pack's backend.

### 8.4 Project dock (M2, P6)

New, Open and Save; pack and recipe chosen at creation; an image list (files, folders, NAS paths, header checks
and refusals shown); a samples.csv editor; a status board from api.project_status; Run selected and Run all
through the runner; double-click opens results; the config hash and a validated: false badge after any edit; the
environment-mix warning; and a device choice (auto, cuda, mps, cpu) and the interpreter per recipe, both recorded
in the project.

### 8.5 Role editor and parameter form (P6)

The role editor offers one combo per role, listing the image's channels, and writes channel_map.yaml and
recipe.yaml through the engine's write_recipe, in ASCII. The parameter form shows only the recipe's tunables, with
their types, units and ranges, so an undeclared key cannot be written. Every edit shows the new config_hash. The
form also lists the installed post-stages with their declared parameters. Switching one on writes
post_stages.yaml, which leaves config_hash as it is (after P5b).

### 8.6 Training dock (P6)

It defaults to the project's training/<round>/ and adds: Import labels (Cell-ACDC npz, Cellpose _seg.npy with an
explicit "allow pickle" choice, label TIFFs); orientations that follow the run (xy for planar, proj for projected,
xy, xz and yz for 3D); base models from the pack and the project; Use in project, which writes recipe.yaml; and
Derive pack.

### 8.7 Spot tuning panel (P6, core)

A Points layer per image is saved as annotations/<image_id>__points_<annotator>.csv. Tune runs `germquant
tune-spots` through the runner. A table ranks the top five sets, and Apply writes best.yaml into recipe.yaml.

### 8.8 2D display (P6)

- ResultsSet.spacing returns (1.0, dy, dx) when voxel_dz_um is NaN and dy and dx are finite. Today any non-finite
  value drops all three to 1.
- spots_image_layer adds the missing axis to a YX TIFF, as _read_labels already does.
- The Points size comes from the run's spots.spot_radius_um, not from the fixed 0.6 um.
- The Time-lapse dock uses 1.0 for a NaN dz.

### 8.9 Versions

pyproject requires germquant>=0.2.0.dev0,<0.3 during T5, and >=0.2,<0.3 after the tag. PEP 440 sorts 0.2.0.dev0
before 0.2, so a plain >=0.2 would reject the development engine. Pack features are gated on hasattr(api,
"list_packs"), and an older engine gets one clear message. The plugin README's "No system packs, no project files"
line goes.

### 8.10 Plugin CI (P4)

napari-germquant has no CI and no git remote today. P4 creates its GitHub repository (public or private: O22) and a
workflow on ubuntu, windows and macos: Python 3.11, pytest-qt with QT_QPA_PLATFORM=offscreen, and the engine
checked out at a pinned commit and installed editable with the readers extra. The test pack is copied into the
plugin's tests/fixtures, and a CI step compares its SHA256SUMS with the engine checkout's copy, so the plugin tests
need no sibling checkout and the two copies cannot drift. While a train is open the pin is the engine's integration
branch head; it moves to the engine's main commit after the engine PR merges, and only then does the plugin PR
merge. This workflow is where sys.executable stands in for the engine interpreter (section 9).

## 9. Headless parity

The claim: in the same environment and on the same device, these routes write the same file list, cell-identical
tables (run_timestamp aside), the same config_hash and stage hashes, and equal done markers. Two markers are equal
when every key except finished is equal (image_id, config_hash, enabled_stages, stage_hashes, xy_stride, z_range,
failed_stages and any extras): provenance.write_done_marker stamps finished with the time, so two markers are never
byte-identical.
1. The plugin through the subprocess runner.
2. The plugin in process, when the plugin environment is the engine environment.
3. `germquant run IMAGE --pack ID --recipe R` and api.run(image, pack=ID, recipe=R).
4. `germquant project run DIR`.
5. The Snakemake project route.

Why it holds: every route ends in api.run_loaded and process_image. There are no hidden settings: the GUI writes
hashed recipe.yaml and channel_map.yaml, the unhashed post_stages.yaml, and otherwise only --no-<stage> overrides.
post_stages.yaml feeds only the post tables and is recorded in __post.json. The pack record, the runtime block and
the sample row are written but never enter config_hash. The subprocess runner removes the difference between the
plugin environment and the engine environment.

Tests:
- engine tests/test_pack_parity.py and tests/test_project_parity.py on tests/fixtures/packs/test_pack, a classical
  recipe found through $GERMQUANT_PACKS;
- C. elegans parity: `--pack celegans_germline` equals `--config config/config.yaml` (hash 84bdda5a5230,
  identical tables), and each recipe equals its `--profile`;
- the plugin's test_gui_run_writes_the_same_files_as_the_cli, extended to the pack route and the project route with
  both runners, and from comparing six marker keys to comparing every key except finished; in the plugin workflow
  (section 8.10), sys.executable stands in for the engine interpreter.

For the paper (P9): scripts/parity_report.py runs routes 1, 3, 4 and 5 on each pack's real example on the RTX 5090
and writes parity_report.json, where 0 differing cells is expected. scripts/device_report.py runs one project on
CUDA in bfloat16, on the CPU, and on Apple MPS in float32 if a Mac is available. It reports agreement (object F1
between runs, count CCC, spot F1), never identity.

## 10. The training loop on an outside system

1. Run the pack, in a project.
2. Export. export_training_pairs(orientations=None) reads run_geometry. Planar runs export the plane (orientation
   xy, tiled at 1024 px with the offsets recorded). Projected runs export the exact 2D image the model saw (the
   chosen plane or projection, read again at the run's geometry) with labels[0] (orientation proj). 3D runs keep
   xy, xz and yz. The CLI default becomes auto, and the manifest of a 3D export is unchanged (tested on the
   synthetic gonad).
3. Or import labels made elsewhere. `germquant import-labels SRC --format tiff|png|cellacdc-npz|cellpose-npy
   --image IMG --pairs DIR` writes the export layout (image/, labels/, and manifest.csv rows with source
   import:<tool>) from Cellpose GUI, Cell-ACDC, ilastik, Labkit or micro-sam outputs. npz files are read with
   allow_pickle=False. _seg.npy files are read only with --allow-pickle, with a warning. Images listed as test in
   any benchmark split are refused.
4. Correct in the Training dock, as today.
5. Train. `germquant train-nuclei --pairs DIR --out MODELDIR --base cpsam_v2@0f1cc3f7`, split by image or series.
   The card adds dims, pixel_um, object_kind, system, pack {id, version, recipe}, segmented_role and, for public
   data, dataset.
6. Learning curves. `germquant train-curve --pairs DIR --out DIR --sizes 1,2,4,8,16,all --repeats 3 --seed S --base
   ... [--keep-models none|best]` trains nested, seeded subsets by source image and scores each model end to end on
   val right after training, through the benchmark runner with an overlay naming the model. It keeps every card
   with its metrics and sha256, plus curve.csv. --keep-models (default none) decides which weights stay: a
   Cellpose-SAM fine-tune is about 1.2 GB (germline_nuclei_combined is 1,218,646,375 bytes), so keeping all 18
   models of a 6-size, 3-repeat curve would take about 22 GB; best keeps only the model with the highest val score.
7. Use. `germquant project use-model DIR --model training/<round>/models/<name>/<file>` writes the relative path
   into recipe.yaml and the digest into the manifest (a new config_hash, validated: false). run_project refuses
   the model if the file it resolves to has another digest (section 7.5).
8. Promote. `germquant pack derive --from ID@VER --recipe R --model PATH --card CARD --out DIR` writes a derived
   pack with a new id: a recipe that extends the parent by versioned name, a model spec (name@sha8, sha256, card)
   and SHA256SUMS, into the project's packs/ or the user folder. The model then appears in the pack menu. It
   rewrites the card's machine paths (model_file, base_model.path) to relative names plus sha256 and escapes any
   non-ASCII text, so no machine path is published and the ASCII lint passes.
9. Spots are tuned, not trained. `germquant tune-spots --project DIR --points annotations/ [--grid FILE]` runs
   section 6.8 on the user's points and writes best.yaml; applying it writes recipe.yaml.

Yeast cells train the same way; YeaZ fine-tuning is not wired. Paper demonstrations (P9): learning curves on
DeepFoci nuclei and on S-BIAD3750 or BABY yeast cells, spot preset tuning on both, and one timed codeless session
(a lab member who did not write the code, 30 minutes, the weakest outside recipe: correct slices, train, score on
val; corrections.json and the card are deposited, with the time to first result and every blocking issue). Before
the packs reach 1.0.0, an outside pilot lab runs a pack on its own data (P8, O20).

## 11. Tests and acceptance gates

### 11.1 Goldens byte-identical first

Every package that touches the run path (config, model resolution, pipeline, segment, spots, measure, io, batch)
proves first that today's outputs did not move. Only then do its new tests count.
1. Synthetic goldens: tests/test_golden.py synthetic and synthetic_4ch pass unchanged, with sentinels 84bdda5a5230
   and 37d49f43ac38, and so do the P0 pins of section 11.2. test_golden.py alone compares five tables, tolerates
   added columns and reads no TIF or PNG, so the pins cover the column list of every TABLES entry, the label TIF
   arrays and the montage.
2. The 23 real goldens on the RTX 5090 (11 dry-ice and stride variants, 12 RAD-51 gonads):
   `python scripts/make_goldens.py tests/golden/real_manifest.json --out <candidate root>`, then
   `python scripts/regression_diff.py <baseline> <candidate> --strict-cols --strict-ids`, then the file-list check
   of scripts/t5_gate.py. The baseline is the T5 base, made once in P0 with `make_goldens.py --commit <base sha>`.
   The comparison must be exact (atol 0), with no added columns, the same file list and the same TIF arrays.
   On its own, regression_diff walks only the golden's files, passes a renumbered label image unless --strict-ids
   is given, and only notes extra images. So t5_gate also checks that both roots hold the same images and, per
   image, the same file names. The existing comparison against C:/Users/ryane/germquant_golden/8666f07 must also
   stay at 0 failing comparisons. These reruns compare outputs with earlier outputs, never with truth, so they
   score nothing on the sealed split.
3. The fence (section 11.2) is unchanged.

### 11.2 The fence (P0)

tests/test_default_invariants.py with tests/golden/default_invariants.json pins, at the T5 base commit:
- the config_hash of config.yaml, config_ccw77.yaml, config_n2dryice.yaml and every config/profiles/*.yaml,
  through load_config (platform independent);
- the LF-normalised sha256 of every config/channel_maps/*.yaml;
- the literal stage_hashes(cfg, role_to_idx, extras=FIXED, model_sha="absent") of config.yaml with the synthetic
  gonad's roles, where FIXED maps every tool in stages._STAGE_TOOLS to "PINNED". Environment changes then cannot
  move the pin, while payload changes do;
- the key sets of run_manifest.json, <id>__done.json and the top level of <id>__stages.json of the default
  synthetic run, and the qc_flags strings of the synthetic and synthetic_4ch runs;
- the legacy single-plane behaviour: a synthetic (C, 1, Y, X) stack through the patched reader seam, with today's
  bioio rule (spacing (1, 1, 1), spacing_ok False), keeps voxel:UNREADABLE_assumed_isotropic_1um and a volume_um3
  equal to the voxel count; and the fake-BioImage reader test pins the spacing rule. Both run in CI without bioio.
- for both synthetic cases: the exact column list of every TABLES entry (all twelve, not only the five that
  test_golden.py compares), the sha256 of each label TIF's array, and the sha256 of the montage's decoded pixels.
  The montage is drawn by matplotlib with text, so its pixels may differ between ubuntu, windows and the
  workstation. P0 checks this; a pin that differs across platforms runs on the workstation only, inside t5_gate,
  and the PR says so. The same rule applies to the label digests.

Pinned files may not change or disappear; new files may be added. scripts/t5_scope_check.py fails a T5 diff that
modifies or deletes a frozen path: everything under R1 and R8, tests/golden/synthetic/ and
tests/golden/synthetic_4ch/ (named one by one, because a glob such as `synthetic*/` would also freeze the new
golden folders that P2a and P2b add) and tests/golden/real_manifest.json. It runs in CI on T5 branches and inside
scripts/t5_gate.py. t5_gate also runs ruff and `pytest -m "not model"` and, when GERMQUANT_GOLDEN and
GERMQUANT_CANDIDATE are set, the real-golden comparison with --strict-cols --strict-ids, the file-list check and
the workstation-only pins, and prints one PASS/FAIL table for PR descriptions.

### 11.3 Gates

- G1 ruff and `pytest -m "not model"` pass on ubuntu-latest and windows-latest (Python 3.11,
  constraints/ci-cpu-2026-10.txt) and, from P1, on macos-14 (arm64, constraints/ci-mac-arm64-2026-10.txt), within
  the 15-minute job budget; each PR records the CI time it adds. The new Mac file holds the pins of
  ci-cpu-2026-10.txt and takes torch from PyPI, because mac-arm64-2026-09.txt pins neither pytest, ruff nor
  laptrack. From P4 the plugin workflow (section 8.10) also passes on all three.
- G2 synthetic goldens unchanged, with the P0 column, label and montage pins (11.1, item 1).
- G3 the fence unchanged: the P0 pins plus the pins added for new files such as generic_2d.yaml.
- G4 real goldens exact on the RTX 5090 (11.1, item 2: --strict-cols --strict-ids and the file-list check), for P1,
  P2a, P2b, P2c, P5 and P5b, and after any change to the validated .venv. A train's integration rerun may reuse the
  last package rerun if nothing merged after it.
- G5 scope check clean (R1, R8).
- G6 import guards in a fresh interpreter: importing germquant.packs and calling its list_packs and a shallow
  pack_status, and running `python -m germquant.packs list --json`, leave torch, cupy, spotmax and cellpose out of
  sys.modules. germquant.api cannot meet this, because it imports pipeline, which imports torch (a probe on this
  worktree confirmed both: config, provenance and models_registry load none of them, api loads torch). A deep check
  imports germquant.pipeline before any backend discovery, and the AST seal of extbench and tune holds.
- G7 style: no U+2013 or U+2014 in new files; ASCII in packs/, new config files and project writes; SHA256SUMS
  current and every path in it tracked by git; pack lint clean.
- G8 recipe goldens with sentinels in CI for every pack recipe, run with an empty models_dir.
- G9 model goldens exact on the RTX 5090, twice, for every recipe with a model.
- G10 planar invariance (identical planar tables for a file dz of 0.2 and of 5.0) and every refusal test.
- G11 the parity tests of section 9.
- G12 pip freeze of the validated .venv unchanged by the packs-extras work (P8).
- G13 the container smoke step, on the workstation or as a manual workflow run: the image builds, `pack list` shows
  the built-in packs (from P1), and the test pack runs on its OME-TIFF example (from P3). No image is published
  before O7 is settled.

## 12. Packages and build order

### 12.1 Overview

| pkg | name | depends on | size |
|---|---|---|---|
| P0 | fence: pin today's behaviour | none | S, 1.5 days |
| P1 | packs core | P0 | L, 8 to 9 days |
| P2a | projected segmentation and pack knobs | P0 | M, 3.5 days |
| P3 | outside packs 0.1 and example fetching | P1, P2a | M, 5 days |
| P4 | plugin MVP: pack picker, runners and plugin CI | P1, P3 | M, 5 days |
| P2b | single-plane images | P2a | M, 4 days |
| P2c | labels from files | P2a | S, 1.5 days |
| P5 | projects, batch formats, sample sheet, training and tuning engine | P1, P2b, P2c | L, 7 to 8 days |
| P5b | post-stage contract | P5 | M, 4 days |
| P6 | plugin: projects, editing, training, spot tuning, 2D display | P4, P5, P5b | L, 6.5 days |
| P7 | public benchmark harness | P1, P2c (P2b for the 2D sets) | L, 6.5 days plus downloads |
| P8 | pack validation, packs 1.0.0, real-Mac smoke test and pilot | P3, P4, P5, P7 | L, 6 days plus GPU and the pilot |
| P9 | paper artefacts, deposits and release | all | L, 6.5 days plus GPU and one session |

Agent estimates have run optimistic before; the milestone times in section 1.4 include integration and the
golden reruns.

### 12.2 Trains and merge rules

- Train A (M1): P0; then P1 and P2a in parallel worktrees; then P3 and P4. The downloads for P7 start on day 1.
- Train B (M2): P2b, P2c, P5, P5b, P6. P7's code is written in parallel, against stubs.
- Train C (M3 and M4): P7, P8, P9.
- Each train merges its engine packages in dependency order into one integration branch made from main, runs every
  gate that applies so far, and opens one PR to main. G9 and G12 apply from Train C, when P8 adds the model
  goldens and the packs extras. The next train's branch is made from main only after the previous PR has merged,
  checked with `git merge-base --is-ancestor <train head> main`. No stacked PRs. The plugin packages (P4, P6) are
  PRs to napari-germquant; each merges only after its train's engine PR has merged and the plugin CI pin has moved
  to that main commit.
- Shared files are edited in small guarded hunks. P1 and P2a both touch pipeline.process_image, cli.py and
  segment/backends.py, so P2a rebases on P1. P2b and P5 both touch batch.py and cli._batch, so P5 rebases on P2b.
  P5b touches api.run_loaded, cli._batch, timelapse.series and project.run_project (which passes the project's
  post_stages.yaml) after P5, so it rebases on P5. Each rebase reruns t5_gate.

### 12.3 The packages

P0 Fence (S). Tests and scripts only; no production code.
- Scope: section 11.2 (with the column, label and montage pins), plus the T5 baseline of the 23 real goldens made
  at the base commit.
- Acceptance: passes on the base commit unchanged on ubuntu-latest, windows-latest and the workstation; the PR
  records a local, uncommitted experiment showing each pin fire (a key added to the default stage-hash payload, a
  comment edited in a profile, a key added to the manifest, a byte changed in a channel map, a changed
  single-plane flag, a column added to a table outside the five, a relabelled label image, a changed montage
  title) and the scope check fire on an edited frozen file; the PR records whether the label and montage digests
  agree across ubuntu, windows and the workstation; the JSON lists 84bdda5a5230 for
  config.yaml and the sentinels 84bdda5a5230 and 37d49f43ac38; the channel-map digests agree between a CRLF and an
  LF checkout; the baseline exists under C:/Users/ryane/germquant_golden/<base sha>, two of its tags rerun exactly
  (a determinism check), and it compares at 0 failing against 8666f07; t5_gate exits 0.

P1 Packs core (L). Sections 3 and 5.1, including the C. elegans pack.
- Scope: germquant.packs (format, discovery, status, resolve_recipe, pack_record, install, vendor, build-sums,
  verify, lint with its key allow-list, deny-list and outside-recipe rules, pin check) and its torch-free entry
  `python -m germquant.packs` with --version; germquant.packs_golden (the tifffile example reader and the golden
  runner); germquant.datafetch (resumable download with size and sha256 checks); the blob store, its verification
  records and fetch-models with Cellpose-cache adoption; the two extends fallbacks, with recorded anchors and the
  portable stage payload for chains through them; Config.pack_record and Config.anchors; the pack model lookup in
  resolve_model_path (unresolved under the classical method) and training._resolve_base; validate_card on the real
  card keys; the pack and runtime records in run_geometry and the batch-root manifest;
  post_run.segment_device_used; declared_backend_names and resolve_model_name; --pack and --recipe;
  `germquant pack ...`; `python -m germquant`; `germquant --version`; version 0.2.0.dev0; built-in roots limited
  to the checkout and $GERMQUANT_CONFIG_ROOT; .gitattributes (text=auto with binary rules) and the .gitignore
  negations with their re-ignores; a macos-14 CI job and its constraints file; packs/ copied into the Docker
  image; packs/celegans_germline with the base-weight check for its card; the test pack.
- Acceptance:
  1. G1 to G7 pass, G1 including the new macos-14 job (a test that cannot pass there for an environment reason is
     marked with the reason and listed in the PR), and the 23 real goldens are exact.
  2. api.resolve_config(pack="celegans_germline").hash == "84bdda5a5230"; for rad51_foci, segmentation_only,
     n2_sc, ccw77_partition and live_syp2 the pack route's hash equals the --profile hash; on the synthetic gonad
     --pack and --config write cell-identical tables and equal done markers (section 9); the manifests differ only
     by the pack and runtime keys.
  3. The test pack copied under two different temporary roots loads with one config_hash; `extends:
     generic.yaml` from outside the checkout resolves to the checkout file; `extends: pack:test_pack@0.1.0/basic`
     resolves; an unknown pack or recipe raises FileNotFoundError listing the candidates and the packs available;
     a test that loads every config and profile in the checkout shows each resolves to the same file as before;
     the two existing extends tests pass unchanged. A profile that extends pack:celegans_germline@<v>/default with
     one parameter edit gives the same stage hashes when the engine config root moves. The test monkeypatches the
     built-in root to one temporary copy of config/ and packs/celegans_germline, then to another, and passes a
     fixed model digest to stage_hashes. An edited profile on the test pack, with the pack copied under two roots,
     also gives the same stage hashes for every stage. An existing out-of-repo overlay with an absolute
     extends path keeps today's stage hashes.
  4. fetch-models adopts a Cellpose-cache file (a small fake) after size and sha256 checks, adopts --from-file,
     downloads from a local HTTP server with Range support and resumes, discards a mismatch, and does nothing on a
     second call. The blob lands in models_dir()/blobs/<sha256>/, resolve_model_path returns it, and model_sha256
     is its digest. A `pack check` in a new process does not re-hash a blob whose size and mtime match its record,
     and re-hashes after either changes.
  5. A local registry entry with a pack model's name but other weights is refused before segmentation;
     register_local_model refuses `@<8 hex>` names; release and Cellpose built-in names cannot be pack names; a
     declared but missing blob raises the fetch-models message under the cellpose and auto methods, while under
     classical the run passes resolution and records model_sha256 "unused".
  6. pack status reports a missing module (find_spec), a missing model, an engine out of range, a recipe hash
     mismatch and a pixel size outside expects (a warning), each with its fix; `pack check` exits 2 when unusable;
     `run --pack` refuses an unusable pack unless --no-pack-check.
  7. install refuses zip-slip entries, drive letters, symlinks, bad sums, non-ASCII hashed files and a collision
     with a built-in id and version; the SHA256SUMS of celegans_germline verifies on a CRLF Windows checkout and an
     LF Linux checkout; every path in each committed SHA256SUMS is tracked by git
     (`git ls-files --error-unmatch`); a PNG committed under packs/ is byte-identical after a fresh checkout; a
     .DS_Store placed under packs/ is neither tracked nor listed in SHA256SUMS.
  8. In a fresh interpreter, germquant.packs.list_packs, a shallow germquant.packs.pack_status and
     `python -m germquant.packs list --json` leave torch, cupy, spotmax and cellpose out of sys.modules; a
     malformed pack under $GERMQUANT_PACKS is listed as unusable, and a generic-profile run still resolves cpsam.
  9. `batch --pack` writes the pack and runtime records into the batch-root manifest; per-image __stages.json of
     pack runs carries post_run.segment_device_used; default runs carry none of this (the fence holds).
  10. `python -m germquant --help` works and imports torch before anything else.
  11. validate_card accepts every card train_nuclei writes today (model_size; size only as an alias) and marks
      model_file and base_model.path as non-portable; lint refuses a recipe key that is neither set by
      config.yaml nor on the allow-list, and a deny-listed key such as spots.channel; the C. elegans card records
      the base-weight comparison, its method and its result.
  12. With the engine installed non-editable into another project's environment, `pack list` names the missing
      built-in root and both fixes. `extends: generic.yaml` from outside the checkout then raises instead of
      binding to that project's or the working directory's config/ folder. packs.vendor copies a pack version
      that then verifies; the Docker image built from the P1 Dockerfile lists the built-in packs (G13).

P2a Projected segmentation and pack knobs (M). Sections 4.1, 4.2, 4.4 to 4.8 for projected runs.
- Scope: planar.py's projected parts (project_image, the plane selectors including sharpest, segment_planar,
  drop_small_area, 2D-object measurement, check_planar); Stage.planar_ok; cellpose_backend(planar) and
  nuclei._cellpose_2d; planar_ok on cellpose and classical; object_diameter_um through the existing kw; the
  SpotMAX optional keys with SpotMAX's sharpen semantics; render base_plane; PLANAR_COLUMNS with extra_declared;
  run_geometry segment_projection and the flags; a synthetic projected golden; segmentation.nuclei.require_channel
  and qc.min_nuclei (sections 3.10 and 3.12).
- Acceptance:
  1. G1 to G7 pass, the real goldens are exact, and the fence is unchanged.
  2. With Cellpose mocked, the default 3D call is exactly today's (do_3D=True, z_axis=0, channel_axis=None,
     anisotropy=dz/dy, no diameter); object_diameter_um adds only diameter = d / dy; the 2D call is
     model.eval(img2d, channel_axis=None).
  3. The plane picker returns a planted sharpest plane, mid, max and mean; an index out of range raises before
     segmentation.
  4. On a synthetic yeast-like stack (C=2, Z=11, a transmitted-like channel and 3D puncta) with classical
     segmentation and SpotMAX: labels are identical in every plane, the area gate works, volume_um3 is NA,
     area_um2 is present, every spot lies inside a cell, segment:projected_plane is raised, 3D-only stages are
     refused or skip as specified, and the golden is checked in with its sentinel.
  5. check_planar refuses, before any pixel is processed, a segmentation backend without planar_ok in a projected
     run and a 3D-only stage that is on with its roles bound.
  6. A spy test shows that, with the optional SpotMAX keys absent, the preprocess, segmentation, detection and
     feature calls are identical to today's; with sharpen set, detection gets the sharpened image and the feature
     call gets the unsharpened image plus sharp_spots_image, and effect_size then equals the value SpotMAX computes
     on the sharpened image; with spot_radius_z_um set, radii are (rz, r, r) / spacing.
  7. base_plane shows the segmented plane in projected runs, and the default montage is unchanged (the P0 montage
     pin).
  8. With require_channel true, a stack without the segmentation role fails the image before segmentation, with a
     message naming the role; unset, the closure behaves as today (segment:no_<role>_channel_using_first_channel).
     qc.min_nuclei moves the few_nuclei threshold only when set, and the fence's qc_flags strings do not change.

P3 Outside packs 0.1 and example fetching (M). Sections 3.12 and 6.9 (tier R), plus examples.
- Scope: packs/cultured_cells and packs/budding_yeast at 0.1.0 (pack.yaml, SHA256SUMS, card, CITATION.cff,
  CHANGELOG, the frozen pack_base.yaml with the QC settings, recipes, channel maps with no transmitted index
  fallback, the cpsam_v2 card, synthetic examples, recipe goldens, example fetch specs), validated: false; the
  readers extra in the Docker image, pinned as below; datafetch's HTTP-range zip-member reader (ZIP64, with a clear
  error naming any compression it cannot read, such as deflate64); germquant.assemble (stacks per-channel files
  into an OME-TIFF with declared axes, channel names, voxel size and crop, and computes pixels_sha256), with a
  `germquant assemble` command for labs that store one TIFF per channel; `pack fetch-example`; api.fetch_example;
  scripts/make_pack_examples.py.
- Readers pins: bioio 3.5.1, bioio-base 3.5.0, bioio-tifffile 1.3.1 and bioio-ome-tiff 1.4.0 go into a new
  constraints/readers-2026-10.txt. These are the engine .venv versions. No constraints file pins bioio today, and
  the image's desktop-2026-09.txt stays as it is. The engine .venv has no bioio-czi, bioio-lif or bioio-ome-zarr.
  So the image either installs only the two TIFF readers, or pins those three in the same file after one
  workstation check.
- Acceptance:
  1. `germquant pack check` passes for both packs on the workstation after fetch-models, which adopts
     ~/.cellpose/models/cpsam_v2 without a download.
  2. `germquant pack fetch-example <id>` then `germquant run --pack <id>` gives object and spot counts inside each
     example's expect ranges, with the pack record in run_manifest.json; the yeast example finishes on a CPU in
     minutes.
  3. Recipe goldens with sentinels pass in CI on ubuntu, windows and macos-14 with an empty models_dir; the PR
     shows that editing a pack's pack_base.yaml fails them, and that editing generic.yaml or config.yaml on a
     scratch branch moves no outside pack hash.
  4. datafetch refuses size and sha256 mismatches, resumes a partial download, extracts a member from a local
     ZIP64 archive served with Range support, and names a compression method it cannot read; the assemble round
     trip gives the declared axes, names, voxel size, crop and pixels_sha256.
  5. Masks from the bare name cpsam and from cpsam_v2@0f1cc3f7 are identical on CUDA on the yeast example plane,
     and docs/PACKS.md records it.
  6. The DeepFoci member, page order, channel order and voxel size are checked on the first download and written
     into the example spec; if a single member cannot be fetched, the documented fallback is taken (open question
     O9).
  7. pack lint is clean: ASCII, no em or en dashes, no spots.channel key (pack_base.yaml included),
     expected_from_germ_cell false on every instance, the QC settings present; each card says what is provisional
     and which scale effect_size is on; a misspelt key added to pack_base.yaml fails lint.
  8. The container image with the readers extra runs the test pack's OME-TIFF example, mounted through
     $GERMQUANT_PACKS (G13).
  9. A stack with no transmitted channel, run with the yeast pack, fails the image with the require_channel message
     instead of segmenting channel 0.

P4 Plugin MVP: pack picker, runners and plugin CI (M). Sections 8.1 to 8.3, 8.9 and 8.10.
- Scope: EngineRunner with SubprocessRunner and InProcessRunner; the first-run interpreter prompt, the version
  check and per-recipe interpreters; listing and checks through `python -m germquant.packs`; the pack picker and
  card; Get models and Get example through the runner; Run with pack and recipe; Results naming, summary,
  Provenance box and fallback banner; export with --pack; the version pin and feature gate; the plugin's GitHub
  repository and CI workflow (section 8.10).
- Acceptance:
  1. The dock lists the built-in packs and those under $GERMQUANT_PACKS through both runners and shows the card;
     the bench profiles are hidden by default.
  2. An unusable pack is greyed with its reason and fix; Get models through the runner (a fake small model) makes
     it usable; Get example (a local stand-in served by the test) opens the assembled image.
  3. A GUI Run with the test pack, against `germquant run IMAGE --pack test_pack` with the exported argv, gives the
     same file list, cell-identical tables and equal done markers (section 9), through the subprocess runner
     (sys.executable as the engine) and in process.
  4. With object_kind cell, the Labels layer is named "<id> cells" and the summary says cells; the banner shows
     when segmentation_method is classical but the pack expects cellpose.
  5. Export writes --pack commands and an overlay with `extends: pack:test_pack@0.1.0/basic`, in ASCII, with no
     absolute path.
  6. An engine older than 0.2.0.dev0 gets one clear message; the existing plugin tests pass.
  7. Manual check on the workstation: both outside packs run on their real examples through the subprocess runner
     into the engine .venv; screenshots go into the README.
  8. The plugin workflow passes on ubuntu, windows and macos with the engine at a pinned commit, and its test-pack
     copy matches the engine's SHA256SUMS.
  9. With no interpreter configured, the dock asks for one instead of falling back to the in-process runner; a
     different engine version or git sha shows a warning; a recipe mapped to a second interpreter runs there, and
     run_start records it.
  10. A pack card refresh through the subprocess runner loads no torch in the engine process and does not re-hash a
      verified blob.

P2b Single-plane images (M). Sections 4.2, 4.3, 4.7, 4.8, 4.9, 4.10, 4.12 and 4.13 for planar runs.
- Scope: io.planar and planar.decide; Stack.spacing_yx and bioio _spacing_yx; io.voxel_size_um; the planar read
  block (NaN voxel_dz_um, read:planar, done-marker extra); stage_hashes(planar=); is_done(segment_hash_planar=)
  with both hashes in cli._batch and timelapse.series; the SpotMAX planar branch; log and net3d refused; the YX
  writers; `germquant info`; generic_2d.yaml; a synthetic 2D golden; the Z=1 time-lapse test.
- Acceptance:
  1. G1 to G7 pass and the real goldens are exact; the fence is unchanged apart from the added generic_2d pin.
  2. The synthetic_2d golden (generic_2d, classical forced, a calibrated plane): every planted nucleus is found;
     at least 90 percent of the planted puncta with at most 10 percent extras; volume_um3 NA; area_um2 within
     2 percent of the planted area; every z_um 0.0; voxel_dz_um NaN; read:planar raised; the label TIF reads back
     as (Y, X) with PhysicalSizeX and Y and no Z.
  3. Planar tables are identical for a file dz of 0.2 and of 5.0; the SpotMAX planar tables are shown for
     PLANAR_RZ_PX 0.5, 1.0 and 2.0.
  4. Where bioio is installed, a calibrated 256 x 256 OME-TIFF written as CYX and as CZYX with Z=1 runs under
     generic_2d with real dy and dx and no voxel:UNREADABLE flag; with io.planar false the same file reproduces
     today's tables bit for bit (the legacy fence).
  5. Refusals: io.planar true on Z > 1 fails before any stage, naming --z-range and project; log, net3d and
     nuclei3d fail in a planar run before any pixel is processed; under auto a 3D file cut to one plane stays 3D.
  6. Stage hashes carry planar only in planar runs; `batch --resume` on a planar folder skips finished images and
     reprocesses after a recipe edit; the default resume is unchanged.
  7. io.voxel_size_um overrides the spacing with stride applied, and voxel:override is raised only when it is set.
  8. A synthetic Z=1 series runs through process_series and register with z shifts of 0, or planar series are
     refused with a documented message.

P2c Labels from files (S). Section 4.11.
- Scope: src/germquant/io/labels.py (TIFF, OME-TIFF and Cell-ACDC npz readers; a _seg.npy reader used only by
  import-labels behind allow_pickle); segmentation.nuclei.labels_from; the segment branch; the labels sha256 in
  run_geometry and the segment sub-hash; per-image resume hashes in cli._batch and timelapse.series.
- Acceptance:
  1. G1 to G7 pass and the real goldens are exact.
  2. TIFF, OME-TIFF and segm.npz (2D, 3D, with a T axis) are read with stride and z_range; 2D labels are copied
     through z and relabelled consecutively; segmentation_method reads labels_file; the flag appears only then.
  3. The labels sha256 is in run_geometry and the segment sub-hash; editing one labels file makes `batch
     --resume` reprocess that image only.
  4. A missing labels file fails that image with a clear message; a _seg.npy path is refused with the
     import-labels hint; a pickled npz is refused (allow_pickle=False).
  5. A detector-only run of the synthetic gonad with truth labels gives the same spots as a run whose segmentation
     produced the same labels.

P5 Projects, batch formats, sample sheet, training and tuning engine (L). Sections 7, 9 (engine side), 10 and
6.8 (engine side).
- Scope: germquant.project and cli_projects.py; per-run logs under runs/ and the manifest lock; the pack copy at
  project creation (packs.vendor); recipe.yaml only after an edit; write_recipe in ASCII; the sample sheet
  (process_image sample=, --samples, the row sha256 in is_done, the flag); io.input_suffixes, image_stem ids and
  output skipping in batch and the Snakefile; project models by relative path and use-model, with the digest check
  in run_project; export by run geometry; import-labels; the new card fields; train-curve with --keep-models; pack
  derive with portable cards; tune.effect_filter, tune.spots and `germquant tune-spots`; the api additions.
- Acceptance:
  1. G1 to G7 pass and the real goldens are exact; two tags also go through `germquant batch --resume` to show the
     skip path.
  2. project run, run --pack and api.run(pack=) on the test pack write identical file lists, cell-identical
     tables, the same config_hash and stage hashes, and equal done markers (section 9); each run's log under runs/
     holds run_start, one event per image and run_end.
  3. An unmodified project has no recipe.yaml and its hash equals --pack; an edit writes recipe.yaml with a
     versioned pack: line, changes the hash and marks every image stale with the reason.
  4. A project made under root A and copied to root B (different absolute paths, one image by relative path, one
     path with spaces) reports the same config_hash, the same stage hashes for every stage, and status done. The
     same holds for an edited project on celegans_germline (one parameter edit through recipe.yaml), whose stage
     hashes are the same under both roots and with the engine config root moved (section 3.4).
  5. Sample sheet: the per-image tables carry the sheet's design columns; adding a row for another image leaves
     finished images done; editing one row makes only that image stale; a missing image falls back to the regex,
     with the flag only in sheet runs; default runs are unchanged (the fence holds).
  6. Batch on a folder of .nd2 and .ome.tif files without io.input_suffixes processes exactly the .nd2 files that
     the base processes; with [".ome.tif"] the ids drop the double suffix; germquant outputs in the folder are
     never picked up.
  7. add_images refuses a multi-scene file, an nd2 with P > 1 and a duplicate id; --shard 1/3, 2/3 and 3/3 cover
     the image list exactly once, and run as three parallel processes on one project they write three run logs,
     lose no event and never write the manifest, and status merges their logs.
  8. export-training on the synthetic gonad writes a manifest identical to today's; on a projected run it writes
     proj pairs of the segmented plane; on a planar run, xy pairs only, tiled above 1024 px.
  9. import-labels converts TIFF, PNG and Cell-ACDC npz (2D, 3D, with T); it refuses _seg.npy without
     --allow-pickle and refuses test-split images.
  10. train-nuclei (with the seams replaced) and --base cpsam_v2@0f1cc3f7 writes a card with dims, pixel_um,
      object_kind and pack; train-curve subsets are nested and deterministic, --keep-models none keeps cards,
      metrics and digests but no weights, and best keeps one model; use-model writes the relative path; the output
      of pack derive verifies, lists and runs, and its cards hold relative paths and ASCII text only.
  11. tune.effect_filter equals benchmark.tuning.filter_effect_size on synthetic data; tune-spots recovers the
      planted optimum on synthetic 3D and 2D examples; run-once-then-filter equals direct runs at each effect
      size; benchmark.detector is never imported.
  12. Where snakemake is installed, `snakemake -n` lists .ome.tif jobs under image_stem ids in project mode.
  13. A project created on test_pack 0.1.0 from a fake built-in root (no version folder) still opens and runs after
      that root moves test_pack to 0.2.0, because its copy in the project's packs/ resolves. run_project refuses a
      model whose resolved file has another digest than the manifest's, such as a decoy placed under
      $GERMQUANT_MODELS/training/.
  14. collect_project on a folder of dense outside-pack results adds no qc:germline_count_outlier flag, while the
      plain `germquant collect` on the same folder still applies 1.8.

P5b Post-stage contract (M). Section 5.6.
- Scope: germquant.post_stages (discovery from metadata; the post_stages.yaml file in projects and packs, and
  --post-stages on run and batch; the call after process_image in api.run_loaded and cli._batch; the series
  refusal; the new tables; <image_id>__post.json with its sub-hash; resume of post-stages only; collect of post
  tables); api.list_post_stages; needs.post_stages in pack status; docs/POST_STAGES.md; the fixture post-stage.
- Acceptance:
  1. G1 to G7 pass and the real goldens are exact; a run without post_stages.yaml writes no post file and calls
     nothing new (the fence holds).
  2. The fixture post-stage runs after process_image and writes only its declared new table and __post.json; no
     existing table, label image or done marker changes.
  3. Changing its parameters in post_stages.yaml reruns only the post-stage under resume; config_hash, the main
     done marker and is_done are unchanged.
  4. A failing post-stage is recorded in __post.json, and the image stays done.
  5. collect writes batch_post_<name>_<table>.csv; pack status greys a recipe whose needed post-stage is missing,
     from metadata only, and G6 still holds.
  6. The project parity of section 9 holds with the fixture post-stage switched on.
  7. With post-stages set, a time series (`germquant run --t`, or api.run given t) is refused with a message;
     `germquant run`, `germquant batch` and `germquant project run` each call the fixture post-stage exactly once
     per image.

P6 Plugin: projects, editing, training, spot tuning, 2D display (L). Sections 8.4 to 8.8.
- Acceptance:
  1. The New, Save and Open round trip keeps the manifest; statuses match api.project_status.
  2. Run all through the dock, with both runners, against `germquant project run`: the same file list,
     cell-identical tables and equal done markers (section 9).
  3. A codeless test: create a project from the test pack, add the example, fetch models (mocked), run, see
     labels and points, change spots.spot_radius_um (the hash changes and validated: false shows), save and
     reopen; `germquant project run` reproduces the files.
  4. The role editor writes an ASCII channel_map.yaml and recipe.yaml and shows the new hash; a headless run of
     the edited project gives the same tables.
  5. The tunables form writes only declared keys.
  6. Importing a Cell-ACDC npz and a Cellpose _seg.npy (with the allow-pickle choice) creates pairs; Train (seams
     replaced) and Use in project write recipe.yaml; Derive pack lists the new pack.
  7. The spot tuning panel saves points in the truth-point format, runs tune-spots through the runner, shows the
     top five and Apply writes recipe.yaml.
  8. A planar run displays with scale (1.0, dy, dx), the YX spots image gains its axis, and the Points size follows
     the run.
  9. The existing plugin tests pass, and the engine .venv gains no Qt or napari package.
  10. The parameter form lists installed post-stages; switching the fixture post-stage on writes post_stages.yaml
      and leaves config_hash as it is, and a headless run of the project gives the same post table.

P7 Public benchmark harness (L). Sections 6.4 to 6.7, and the comparator scoring.
- Scope: germquant.extbench (datasets, adapters, manifest, split, score, run, report, ledger) and
  cli_extbench.py (`germquant extbench prepare|split|run|score|report|score-external`); benchmarks/ext/<dataset>/
  with dataset.json, files.tsv, split.json and draft protocols; the shared scorer; ceilings; the self-check;
  detector-only mode through labels_from; benchmark runs as projects (through api.run_project once P5 has merged,
  a direct api.run_loaded loop before); the seal test. On the workstation: downloads, conversions, and baseline
  numbers on the tune splits (not gated).
- Acceptance:
  1. The seal test passes and the scope check is clean.
  2. Each adapter converts its stand-in layout in CI; converted files read back with the expected spacing and
     channel names (tifffile in CI, bioio where installed). The spotmax_gt adapter writes MDN1 as smFISH_MDN1,
     which binds to foci by name with no channel:foci:name_unmatched flag, and marks Position_4 detector-only; the
     baby adapter writes each stack as it is and as one plane, and records the plane. The stand-ins cover every
     truth layout read so far; the layouts not yet read are recorded on the first download and get their own
     stand-ins.
  3. The scorer matches hand-computed F1, AP, matching, CCC and bias on toy cases, and equals benchmark.nuclei3d
     and benchmark.spots on shared cases; projected truth is matched in 2D.
  4. The test split refuses without --final, a committed protocol and a clean tree; final_runs.jsonl is
     append-only; a second final run needs --rescore-reason.
  5. On the workstation, the real downloads give the documented counts (S-BIAD3750: 11 MDN1 positions, 10 of them
     end to end, 30 SUT509 and 189 mtDNA positions; the DeepFoci validation subset; BBBC039: 200 images; the BABY
     groups), and every file matches files.tsv.
  6. The DeepFoci self-check reproduces the published inter-expert median Dice (about 0.75 at 1.95 um) within
     the protocol tolerance; the S-BIAD3750 annotator ceilings are computed.
  7. score-external scores the authors' SpotMAX and RS-FISH outputs in S-BIAD3750.
  8. Baseline numbers on the tune splits are written to docs/EXTERNAL_BENCHMARKS.md, marked as not gated.

P8 Pack validation, packs 1.0.0, real-Mac smoke test and pilot (L). Sections 6.6 to 6.10 applied, and 5.4.
- Scope: protocols frozen and committed before any test score, with the PI's gate numbers; contamination checks;
  tuning on the tune splits within budget (plane selector and min_area for yeast; projected against 3D for cells;
  object diameter; spot radius, z radius and effect size, with effect size searched per sharpen value; the z-merge
  ablation); selection on val; recipes frozen at 1.0.0 with sums and goldens regenerated; one final score per
  recipe; validation records; validated only where the gates pass; the smoke tier; the single-plane recipes
  (nuclei_foci_2d scored on BBBC039, cells_2d on the single BABY planes the adapter writes) and cells_spots scored
  on the 5-plane BABY stacks for segmentation; comparators in their own environments; the optional germquant-yeaz
  backend in .venv_packs, with a Mac constraints file for the pack extras; the real-Mac smoke test; the pilot.
- Acceptance:
  1. Git history shows each protocol.md and split.json committed before its first tuning record; tuning stays
     within budget, with every trial logged.
  2. Each final_runs.jsonl holds exactly one final scoring per recipe version.
  3. Validation records exist with metrics, intervals, ceilings, comparators, contamination notes and gate
     verdicts; status validated mirrors the verdict; the packs are at 1.0.0 and their sums current; `pack lint
     --previous` passes.
  4. Model goldens are exact on the RTX 5090, twice; recipe goldens pass in CI; the smoke metrics are pinned.
  5. germquant-yeaz, if built: it is discovered by name in .venv_packs, declares planar_ok, refuses non-planar
     input, records its version and the weights' sha256, passes unit tests with a tiny fake weights file, and
     agrees with Cell-ACDC's YeaZ_v2 on 5 images or documents the difference; pip freeze of the validated .venv is
     unchanged.
  6. G1 to G13 pass.
  7. On a real Apple Silicon Mac (O8): `germquant check-gpu`, then one example per outside pack through the
     plugin's subprocess runner with .venv/bin/python, recorded in docs/MAC_METAL.md. Without a Mac, the README
     keeps its line that the Mac version has not run on a real Mac, and no Mac claim enters the paper.
  8. Before the packs reach 1.0.0, one outside pilot lab (O20) has run a pack on its own data; time to first result
     and every blocking issue are recorded in docs/PILOT.md. If no pilot lab is available, the PI's decision to go
     on without one is recorded there.

P9 Paper artefacts, deposits and release (L). Sections 1.5, 9 and 10 (paper parts).
- Scope: docs/PAPER_REPRODUCE.md with the claim map, a STAR Methods key resources table draft (software with
  versions, models with sha256, datasets with DOIs) and a limitations list that mirrors section 6.10 and risk 17;
  parity_report.json; device_report.json; learning curves, within a written disk and GPU budget; the z-merge
  ablation figure; the timed codeless session by a lab member who did not write the code, with its deposit; figure
  scripts; CITATION.cff for the engine and every pack; pack zips and Zenodo DOIs written back as PATCH releases;
  Zenodo archives of the engine v0.2.0 and plugin 0.2.0 tags; a deposit of benchmarks/ext with the per-image output
  tables and test-item label images, with attribution for the CC BY sets; pixi.toml and environment.yml regenerated
  from the validated pins or removed before the outside-lab tutorials; engine v0.2.0; plugin 0.2.0.
- Acceptance:
  1. parity_report.json shows 0 differing cells across the subprocess-runner GUI route, run --pack, project run
     and Snakemake, on each pack's real example on the RTX 5090.
  2. device_report.json tabulates agreement across devices: the CPU always, MPS if a Mac is available (otherwise
     the Mac claim stays out of the paper).
  3. The learning curves, tuning results, the z-merge ablation and the codeless session can be reproduced from
     their project folders with the commands in PAPER_REPRODUCE.md. The session record names the person's role (not
     the developer), the time to first result and every blocking issue.
  4. Every pack version has a zip with its sha256 and a Zenodo DOI, written back into pack.yaml and CITATION.cff
     as a PATCH release.
  5. Engine v0.2.0 is tagged after G1 to G13 with exact real goldens; plugin 0.2.0 is released with
     germquant>=0.2,<0.3; the tag commit is an ancestor of main.
  6. The engine v0.2.0 and plugin 0.2.0 tags each have a Zenodo DOI; the benchmark deposit has a DOI and carries
     the CC BY attributions; PAPER_REPRODUCE.md holds the key resources table draft and the limitations list.
  7. pixi.toml and environment.yml are regenerated from the validated pins or removed, and the tutorials use the
     README's install route (the engine as an editable checkout, the plugin in its own environment).

### 12.4 Docs

All docs use plain language, Windows and Mac command pairs as in the README, and no em or en dashes.
- New: docs/PACKS.md (with the object_kind naming table), docs/PLANAR.md, docs/PROJECTS.md, docs/TRAINING.md,
  docs/POST_STAGES.md, docs/EXTERNAL_BENCHMARKS.md, docs/PAPER_REPRODUCE.md, docs/PILOT.md,
  docs/tutorials/budding_yeast.md and docs/tutorials/cultured_cells.md (codeless first, then headless),
  packs/<id>/README.md cards, benchmarks/ext/README.md and plugins/germquant-yeaz/README.md.
- Updated: README.md (packs and projects quickstart for outside labs, the install route with the engine as an
  editable checkout and the plugin in its own environment, batch formats, single-plane images),
  docs/GENERIC_PROFILE.md (its "not supported yet" 2D section becomes a pointer to generic_2d.yaml and PLANAR.md,
  plus the cpsam note), docs/BACKENDS.md (planar_ok, declared_backend_names, resolve_model_name, the pack rule,
  germquant-yeaz as the worked example, a pointer to POST_STAGES.md), docs/MAC_METAL.md (the real-Mac smoke result,
  P8), docs/TIMELAPSE.md (a pointer to TRAINING.md), tests/README.md (the fence, the golden tiers and the macOS
  job), docs/PROPOSAL_general_modular_tool.md (a short status section 14), the plugin README (install, interpreter
  choice and CI) and CITATION.cff.

## 13. Risks

1. Frozen texts. Editing config.yaml (84bdda5a5230), config_ccw77.yaml, config_n2dryice.yaml, any existing
   profile or any channel map, even a comment, moves config_hash for every historical run and breaks resume.
   Mitigation: every recipe is a new file; P0 pins every shipped hash and map digest; the scope check refuses the
   edit. Outside packs extend their own frozen pack_base.yaml, so the planned hash bump moves no outside pack hash;
   it moves the C. elegans pack's recipes, which are the frozen files themselves, and user or project files that
   extend generic.yaml.
2. Unconditional provenance. A new key in the default stage-hash payload or in _STAGE_TOOLS, a new manifest or
   done-marker key, or a new qc flag on the default path changes every default sub-hash, makes `batch --resume`
   reprocess everything and moves the qc_flags golden column. Mitigation: conditional keys only, pinned by P0
   with environment-independent tool versions.
3. A placeholder dz leaking into planar numbers. Any planar path that still multiplies by spacing[0] writes a
   value in a fake unit. Mitigation: 3D-only stages are refused, voxel_dz_um is NaN, volume_um3 is NA, and a test
   requires identical planar tables for a file dz of 0.2 and of 5.0.
4. Editing the 3D Cellpose or SpotMAX path by mistake while adding branches. Mitigation: new functions, or an
   if-branch whose else is the old lines verbatim; mocked-call and spy tests; exact real goldens before merge.
5. cuBLAS load order on the RTX 5090. A discovery step or backend that imports CuPy (through SpotMAX) before torch
   silently drops Cellpose to the watershed. Mitigation: packs and projects are static data; declared names come
   from metadata; listing goes through germquant.packs, never germquant.api; a fresh-interpreter test; pack recipes
   use method: cellpose, so a fallback becomes a failure.
6. GUI and headless numbers differing by environment. The plugin venv has the CPU torch wheel and no SpotMAX, and
   it points at a different checkout than the engine .venv. Mitigation: the subprocess runner is the default
   route, and the dock asks for the engine interpreter instead of guessing and checks its version and git sha; the
   runtime block and the run logs record the device, interpreter and engine sha; status warns on a mix; identical
   numbers are claimed only within one environment and device.
7. Environment drift from optional backends (yeaz, a newer Cell-ACDC, bioimageio.core, Spotiflow, StarDist with
   TensorFlow, micro-sam with gurobipy). Mitigation: separate distributions with constraints files, .venv_packs
   and benchmark environments only, the pip freeze check (G12), and exact real goldens before any change to the
   validated .venv.
8. Model resolution relaxed by pack lookups. Mitigation: the lookup runs only where today's code raises; the pin
   check by sha256; reserved `@<8 hex>` names; the refusal of unknown bare names stays.
9. The extends fallback binding a profile to an unintended base: a study file outside the repo with `extends:
   config.yaml` now resolves to the checkout's config.yaml instead of raising. Mitigation: the fallback runs only
   after today's candidates fail, the resolved path is logged, the extends text is hashed, docs/PACKS.md states the
   order, and lint forbids config/config.yaml inside packs.
10. Batch discovery widening: a broad glob on the NAS could pull exported TIFFs into a worm batch. Mitigation:
    io.input_suffixes defaults to [".nd2"], only outside recipes widen it, and discovery skips germquant outputs.
11. Encoding and line endings. Hashed YAML is decoded as cp1252, so bytes such as 0x81 or 0x8D from an accented
    folder name break loading, and core.autocrlf=true changes file bytes on Windows. Mitigation: ASCII writers
    with escapes, install refuses non-ASCII hashed files, LF-normalised digests, and text=auto attributes with
    binary rules under packs/.
12. RAD-51 seal and concurrency. The final scoring work is active, and the public harness shares scorer modules
    with it. Mitigation: the scope check forbids T5 edits there, the AST seal limits imports, nothing names the
    split file, and the smoke tier flags drift if the RAD-51 work changes the shared modules.
13. Benchmark bias and leakage. Mitigation: independent sets (DeepFoci, BABY, Spotiflow), comparators scored by
    the same code, splits and minimums committed before tuning, one final score per pack version, and the harness
    self-check.
14. Licences. Cellpose-SAM weights are BSD-3, but their training data include CC BY-NC sets; YeaZ weights and
    training data carry no licence; SpotMAX's licence reads GPL-3 in its LICENSE and BSD-3 on PyPI; DeepFoci and
    BABY require attribution. Mitigation: licence notes on every card, YeaZ weights fetched from source and never
    redistributed, unlicensed and permission-only data kept out of the paper benchmark, SpotMAX never bundled in a
    pack.
15. Upstream availability. Model and data URLs (Hugging Face, the Helmholtz file box that hosts the YeaZ weights,
    EBI, Zenodo) can move. Mitigation: pinned sha256 and size, commit-pinned Hugging Face URLs, adoption of local
    copies, resumable fetches, and a request to mirror YeaZ weights as release assets.
16. DeepFoci unknowns: the zip member names, page and channel order, whether the foci subsets carry nucleus masks,
    and whether members use a compression Python's zipfile cannot read. Mitigation: P3 checks all of this on the
    first download before any table is frozen, and the fallback is documented.
17. Projected segmentation semantics. A mask copied through z gives column intensities, and spots above or below a
    cell in its column count for it. Mitigation: volume_um3 is NA, area_um2 is added, the flag marks such runs, the
    card says which recipe to use for thick or overlapping samples, and spot assignment is scored against 2D truth.
18. Naming: tables keep nucleus_id, n_nuclei and the nuclei table name even when each row is a whole cell.
    Mitigation: the GUI and cards translate through object_kind, and docs/PACKS.md has a mapping table.
19. Pack and project drift. Mitigation: immutable installed versions, versioned pack: lines, a recomputed hash
    on open, and a copy of the pinned pack version in every project, so a git pull that moves a built-in pack
    cannot strand a project.
20. Untrusted inputs: _seg.npy files are pickles, pack zips can carry path traversal or symlinks, and pack folders
    could hold scripts. Mitigation: allow_pickle=False by default, --allow-pickle only in import-labels, zip-slip
    and symlink refusal, and packs treated strictly as data.
21. Merge conflicts and skipped reruns. Mitigation: small guarded hunks, the fixed merge order, rebase before
    merge, t5_gate after each rebase, and real-golden reruns as merge requirements.
22. Schedule. P6, P7 and P8 are the long poles, a Mac may not be available, the pilot depends on another lab, and
    the plugin depends on napari 0.9 APIs. Mitigation: M1 lands early; the Mac report, the YeaZ backend, the BABY
    and Spotiflow scoring and, if the PI agrees, P5b are the first to move to a follow-up; the macOS CI job and the
    real-Mac smoke test are not. If P5b moves, so do P6 acceptance 10, the post-stage listing in section 8.5, P6's
    dependency on P5b, item 9 of section 1.2 and the words "and post-stages" in claim C1.
23. Git ignore rules dropping pack files. .gitignore ignores `models/` and `*.tif` everywhere, so a plain `git add`
    leaves cards and examples out and every clone fails pack verify. Mitigation: negations at the end of .gitignore
    and a CI test that every path in every SHA256SUMS is tracked.
24. Concurrent project writers. Shards and Snakemake jobs run at once, often on SMB or NFS shares where appends
    from several hosts are not atomic. Mitigation: one log per run, workers never write the manifest, and
    interactive writes take a lock file.
25. Machine paths in hashes and cards. Inherited relative paths are anchored to absolute ones, and train_nuclei
    cards carry absolute model_file and base_model.path values. Mitigation: portable stage payloads for chains
    through the new extends forms, and pack derive rewrites card paths.
26. The wrong channel bound without a QC failure. A name miss falls back to an index with only an advisory flag,
    and an unresolved segmentation role segments channel 0. Mitigation: no index fallback for the yeast transmitted
    role, require_channel in outside packs, and adapters that write channel names the maps bind (section 15f).
27. Post-stage code. A post-stage is third-party code that runs after each image. Mitigation: separate
    distributions outside the validated .venv, discovery from metadata only, new tables only, its own record, and a
    failure that never fails the image.
28. No pilot user. The proposal makes going on depend on a pilot lab actually running the tool. Mitigation: a named
    non-developer for the timed session, an outside pilot lab before 1.0.0, or a recorded PI decision to go on
    without one (O20).

## 14. Open questions

For the project lead:
- O1. Decided 2026-10-11: both are merged (PR #9 at 65b494e, PR #10 at bc62476, merge commits, no rebase), and
  `git merge-base --is-ancestor` confirmed both heads in main. The T5 base is bc62476, whose tree equals the
  spotmax-grid head 3fcb4b5; feat/tranche5 has merged main.
- O2. Point the plugin venv and the engine .venv at the same checkout (main), and set one shared
  $GERMQUANT_MODELS on the workstation?
- O3. Decided 2026-10-11: the project lead left it to convention. The built-in packs (celegans_germline,
  cultured_cells, budding_yeast and the test pack) live in the engine repository under packs/, versioned and
  tested with the engine, as section 3.8 describes; packs made by others install from their own folders or
  distributions through $GERMQUANT_PACKS and the entry-point group. That is the usual split for a tool at this
  stage: shipped defaults with the code, community content outside it.
- O4. Add the readers extra to CI (bioio paths tested in CI, at some install time), or rely on the tifffile seam
  plus workstation tests?
- O5. Fix the cpsam drift in generic.yaml in the planned hash-bump release (number-changing for generic-profile
  runs only)? P1's state-dict comparison (section 3.11) should settle the base weights of germline_nuclei_combined;
  is a card that says "not determined" acceptable if it does not?
- O6. Is the RTX 5090 time for real-golden reruns at P1, P2a, P2b, P2c, P5 and P5b acceptable, or should some be
  batched at train integration? Is the learning-curve budget acceptable: 36 trainings for two curves, and disk only
  for the models --keep-models keeps?
- O7. Who settles the SpotMAX licence question (its LICENSE file reads GPL-3, its package metadata BSD-3)? The
  Dockerfile already installs SpotMAX, so no image, with or without packs, is published until it is settled.
- O8. Is a real Apple Silicon Mac available for the M3 smoke test (P8) and the device report? If not, the Mac claim
  stays out of the paper and the README keeps its warning.
- O9. If a single DeepFoci member cannot be fetched by range requests, is the fallback acceptable: one S-BIAD1494
  stack (CC BY 4.0, unlabelled 3D gamma-H2AX and 53BP1; sizes and channel order not yet checked) as the
  cultured-cell example?
- O10. Which Zenodo account or community holds the pack and weights DOIs, and the engine, plugin and benchmark
  deposit DOIs?
- O21. After 0.2.0, should the engine go to PyPI and the plugin to napari hub? That needs the built-in packs and
  configs shipped as package data with base-dir handling. Version 0.2.0 is installed as an editable checkout
  (section 3.8).
- O22. napari-germquant has no git remote. Should its GitHub repository be public or private? If the engine
  repository is private, the plugin workflow needs a read token to check out the pinned engine commit.

For the PI:
- O11. Confirm the pre-registered minimums for "validated" (section 6.6) before any test score.
- O12. Which outside system leads the paper, yeast or cultured cells? Can a collaborating lab supply 20 to 30
  Rad52 yeast volumes and 20 to 30 RAD51 cultured-cell volumes, each with two annotators? No public set has point
  truth for either, and such a set would also show the training loop.
- O13. Is C.K. Cahoon on idr0063 our PI? It is the only public yeast meiosis 3D plus time spot set.
- O14. May we request permissions: the Cell Tracking Challenge organisers; the YeaZ authors (a licence for weights
  and training data, and permission to mirror the weights); the YeastMate authors?
- O15. Should the 19 C. elegans COSA-1 volumes in S-BIAD3750 (three annotators) become an external worm benchmark
  in this tranche, as a stretch goal, or later?
- O16. Is time-lapse for outside systems needed in the paper, or does the worm live series carry the time-lapse
  claim for now?
- O17. Should per-image channel bindings outside the config hash come later (they need the binding in the done
  marker and in is_done), or is a study channel map with a new hash acceptable to outside labs?
- O18. The post-stage contract (P5b, section 5.6) is the second half of the "plugins in two senses" promise in the
  proposal. Build it in T5 as planned, or move it to the next tranche? It departs from proposal section 3 in two
  ways. There is no local modules folder (R5 forbids running code from a folder). Parameters are declared like
  recipe tunables instead of in a pydantic Params model (no new dependency in the validated environment). Accept
  both?
- O19. Accept the departures from the proposal (its section 11; this design's section 5.5)? Yeast segmentation
  defaults to Cellpose-SAM instead of Cell-ACDC's yeast models; YeaZ comes as an optional direct backend
  (germquant-yeaz); the Cell-ACDC model adapter and bioimage.io backends are left out of T5. Reasons: Cell-ACDC's
  model installers prompt and pip-install on Enter; its StarDist and cellpose_v3 wrappers and Cell-ACDC 1.8.x
  (which needs pandas below 3.0) would change the validated environment; YeaZ weights carry no licence.
  Replacements: the direct yeaz backend now, and bioimage.io as a later package.
- O20. Which outside lab pilots the tool on its own data before the packs reach 1.0.0, and which lab member who did
  not write the code does the timed session (P8, P9)?

Answered 2026-10-23, for Train B (milestone M2: P2b, P2c, P5, P5b and P6). The Train B plan
(C:/Users/ryane/germquant_golden/t5_train_b/PLAN.md, section 7) put its questions Q1 to Q6 to the project lead, and
its PI questions A to C to him as well, under "For you or Cahoon"; each came with the recommendation yes. His reply,
at 2026-10-23 17:36 MDT (23:36 UTC), verbatim: "Yes to all. I will have to get the spo-11 info from her on Monday
she is gone." The second sentence concerns another task (the RAD-51 spo-11 test). So every recommended answer
holds, and Train B proceeds on the plan as written:
- Q1, the time-lapse code: yes. While feat/flow-prior-linking is open, Train B edits nothing in
  src/germquant/timelapse/ or docs/TIMELAPSE.md. In this release single-plane and labels_from time series are
  refused with a clear message, post-stages are refused for time series, a project's sample sheet does not apply to
  time-series images, and the TIMELAPSE.md pointer to TRAINING.md comes in a follow-up (the plan's FU1).
- Q2, the real goldens (the first question of O6): yes. G4 runs once, on the final Train B head. After each package
  merge a 32-minute check on five real tags runs instead, and G4 reruns if anything other than Markdown changes
  after it starts.
- Q3, the engine CI (O4): yes to both parts. The windows-latest tests run in two jobs (the G1 note below), and every
  CI job installs the two bioio TIFF readers from constraints/readers-2026-10.txt, so single-plane OME-TIFF reads
  are tested on all three systems.
- Q4, one write path for projects: yes. P5 adds `germquant project edit` and a torch-free
  `python -m germquant.project` entry (new, add, edit, use-model, rounds, unlock), and the plugin makes every project
  change through the engine (plan D39 and D42).
- Q5, the order with the other workflows (the flow-prior linker, the SC ribbon work and the RAD-51 LoG recipe): yes.
  Whatever is ready before Train B's final G4 merges to main first, and feat/t5-train-b merges main before G4;
  whatever comes later waits for the Train B PR and then updates from main; nothing merges to main between Train
  B's G4 and its PR; and if a PR merged before G4 changes real-golden output, Train B's G4 compares against a new
  baseline made at that main commit.
- Q6: yes to the six decisions below.
- A, the post-stage contract (O18): yes. P5b is built in Train B, with no local modules folder and with parameters
  declared like recipe tunables.
- B, single-plane time series (O16): yes. They are refused in this release; time-lapse stays a worm claim, carried
  by the live GFP::SYP-2 series.
- C, the study channel map (O17): yes. One channel_map.yaml per project, with its own config hash, is acceptable to
  outside labs in T5; per-image channel bindings outside the hash come later.

The six decisions that change this design's text, signed off with Q6 (plan section 8):
- D13 (4.11): the frame of a T axis in a labels file is chosen with a new key, segmentation.nuclei.labels_frame,
  not in the labels_from template. It is an int of at least 0, default 0, read only in labels runs, and it is added
  to recipe_keys and to the lint allow-list.
- D14 (4.11): labels_from works with segmentation.nuclei.project and with io.planar. In a projected run the labels
  file replaces segment_planar's output: a 2D file (or one plane) is copied through z, and the projected measure and
  records stay (volume_um3 NA; area_um2 and equivalent_diameter_um from the 2D labels; the flag
  segment:projected_<mode> and the segment_projection record), with segment:labels_from_file added. The projection
  itself is computed only for the montage's base plane, when the segmentation channel exists, and in a projected run
  a labels file with more than one plane is refused. With io.planar a 2D labels file is used, and check_planar does
  not test the segmentation backend in labels runs. Labels runs have no size gate: min_volume_um3 and min_area_um2
  do not apply.
- D16 (7.8 and 3.13): the project API is new_project, open_project and the rest of 7.8. create_project, which 3.13
  lists, is not added, not even as an alias.
- D22 (section 10, step 6): train-curve scores each model on the held-out pairs with the existing holdout_metrics and
  records the scorer's name. A scorer= seam lets P7 plug in the end-to-end scoring through the benchmark runner.
- D24 (6.8): tune-spots tunes the recipe's spot backend: SpotMAX through run-once-then-filter, as 6.8 says, and any
  other backend, LoG included, with one detection per grid point. The default grids keep every other key at the
  recipe's value: effect_size_min in {0, 0.5, 1, 1.5, 2, 3, 4, 5} for spotmax and min_snr in {1.5, 2, 2.5, 3, 4, 5}
  for log; --max-runs defaults to 50. The objective is per-object Lin's CCC over the annotated objects, rounded to
  3 decimals, then point F1, the order the RAD-51 harness uses.
- D38 (8.7 and 8.8): spot tuning is a plugin dock of its own, and the Time-lapse dock is not changed (8.8's NaN-dz
  rule matters only once planar series are allowed).

G1 note: G1's budget stays 15 minutes per job (11.3). windows-latest took 15.4 minutes on main 8eb73de (its Test step
13.9), so from Train B on its tests run in two CI jobs, shards 0 and 1, each running every other file of the sorted
tests/test_*.py list; ubuntu-latest and macos-14 keep one job (.github/workflows/ci.yml). Each PR records the CI time
it adds, and every Train B hand-off records the Windows test seconds it adds, against about 30 s for a small unit
and 90 s for a P5 unit or P2b-CORE (plan D41). If a shard still passes 15 minutes, the shards are rebalanced before
any test is trimmed.

macOS note (2026-10-30, plan decision D54): from 2026-10-30, G1's macOS job runs on macos-15, the macOS 15 arm64
(Apple Silicon) runner image, instead of macos-14 (.github/workflows/ci.yml). GitHub retires its macOS 14 images on
2026-11-02 and until then fails macos-14 jobs in scheduled brownouts (github.com/actions/runner-images issue 13518:
14:00 to 00:00 UTC on October 29, 30 and 31 and November 1); on 2026-10-30 a brownout cancelled the plugin CI's
macos-14 job in all three attempts of one run. GitHub's ways off macos-14 are macos-latest, macos-15 and macos-26.
macos-15 is the same Apple Silicon runner class one release newer, so constraints/ci-mac-arm64-2026-10.txt, torch
and torchvision from PyPI, golden_rtol 1e-6, one job (shard "all") and the 15-minute budget carry over unchanged;
macos-latest moves to each new macOS on GitHub's schedule (it is macOS 26 today), and macos-26 is the larger step.
Where this design names macos-14 for CI (G1 in 11.3, P1 and P3 in 12.3, the G1 note above), read the macOS job,
macos-15 from this date. The first comment line of constraints/ci-mac-arm64-2026-10.txt still names macos-14,
because Train B edits nothing under constraints/. GitHub keeps the two newest macOS versions, so macos-15 holds until
the images of the next macOS arrive; a later train then moves the label again. The plugin workflow (8.10) takes the
same label in its own commit.

Decisions and contract changes after 2026-10-23 (Train B). The plan's integration agent took decisions D43 to D59
after the answers above and kept the plan's interface contract (its section 9; the record of each change, with its
date and the units it concerns, is its section 9.4); the files are in C:/Users/ryane/germquant_golden/t5_train_b/
(PLAN.md, handoff/T5B-INT_decisions.md). Two of the plan's own decisions are listed first, because docs/LABELS.md
cites them here with T5B-0's D13 and D14:
- D12 (4.11): labels runs (segmentation.nuclei.labels_from set) skip model resolution and the pin check and record
  model_sha256 "unused", as the classical method does; `germquant batch` follows the same rule. run_project skips
  the segmentation model in its digest check with a warning, and project_status reports that model as not needed.
  pack_status does not require the segmentation model of a recipe that sets labels_from itself: a problem with that
  model (missing, the wrong size, another sha256) is a warning, not a blocker.
- D15 (4.11, 4.13): labels_from with a time series is refused wherever D9 refuses a single-plane one: `germquant run
  --t` before any timepoint, run_project before process_series, and process_image when given t. The doc is
  docs/LABELS.md.
- D43 (4.7, 5.1): the planar stage-flag channel is internal to the built-in spotmax backend. A planar call returns
  tables whose attrs["stage_flags"] lists flag names; the pipeline raises each as spots:<name> (for an instance
  spots_<instance>:<name>) and removes the key from every table of every backend before anything is written. 5.1's
  backend contract (planar_ok and the planar keyword) is unchanged, and attrs that a third-party backend leaves are
  dropped. spots:z_merge_off_planar is raised only when the recipe asks for the z-split merge
  (spots.merge_z_columns true). Confirmed by the project lead on 2026-10-30.
- D44 (4.9, 7.8): `germquant batch` passes the single-plane segment hash to its resume check only when the recipe
  sets io.planar, so a default batch makes exactly the base's call.
- D45 (4.13): every time-series refusal of the command line (single-plane, labels_from, a sample sheet,
  post-stages) prints "REFUSED: <message>" on stderr and exits 2 before any timepoint; "ERROR:" stays for usage
  errors.
- D46 (4.11): the per-image labels hash is public: germquant.io.labels.segment_hashes(cfg, image_path, *,
  model_sha="unused") returns (segment_hash, segment_hash_planar, problem), for labels runs only, and
  UNREADABLE_HASH ("labels-unreadable") for a file it cannot hash, so batch, run_project and project_status
  resume labels runs alike.
- D47 (7.8): `germquant project run --device` accepts auto, cuda, mps and cpu.
- D48 (7.7): the engine run log's done_sha256 is the sha256 of the done marker's bytes with CRLF normalised to LF.
- D49 (7.8): `project run`, `status`, `collect` and `export-commands` refuse as the other project commands do (exit
  2, "REFUSED: <message>" on stderr, {"error", "kind"} on stdout with --json); `project run` exits 0, 1 or 2;
  commands built before keep their recorded forms.
- D50: the command line's docstring and the torch-free pack entry's first docstring line name `pack derive`; a
  test ties the latter to the light parser.
- D51: base code that takes plain paths (packs/install.py, segment/training.py's reads) moves to fsutil.long_path
  in a follow-up; until then docs/PROJECTS.md asks for a project folder path under about 170 characters on Windows
  without long paths.
- D52 (4.11): in a labels run the resume reason is "segment stage hash differs (the labels file, or the Cellpose or
  torch version, changed)"; every other run keeps the base's text.
- D53 (8.4, 8.7): the plugin's Project and Spot tuning docks are in napari.yaml, and the one plugin test that pins
  the exact set of docks is widened by the commit that adds a dock (it stays an equality). Confirmed by the project
  lead on 2026-10-30.
- D54: the macOS note above.
- D55 (11.3): the style check of gate G7 runs on tranche 5 branch names only: in CI's T5 step, and in
  tests/test_t5_style_check.py's check of the checkout, which skips on any other branch and on a detached HEAD. The
  project lead's word (2026-10-30): the no en or em dash rule is for prose, not code, so main gets no dash check.
- D56 (10): export-training of a `germquant batch` tree takes each run's config from the batch manifest, so a 3D
  recipe that segments a role other than dna exports that role's channel (the base exported channel 0 as dna); no
  run number changes. Taken in Train B by the project lead's yes of 2026-10-30.
- D57 (8.10): the plugin's run widget test reads config_hash as text, so a hash of digits with one inner e is not
  read as a number. The project lead's yes of 2026-10-30.
- D58: the last unit merges into the P5 package branch went in as one --no-ff octopus merge per round (unit heads
  synced on one package head, each gated and green on its own), with the gate rerun on the merged tree.
- D59: the final review ran in two parts: part A on the P5 package head before its merge into the train (the whole
  train against main), part B on the head G4 ran on, before G4 (every delta since part A).

The interface contract's changes (plan section 9.4, numbered as there):
1. (2026-10-24) the planar stage-flag channel of D43.
2. (2026-10-24) planar.PLANAR_SERIES_REFUSED, planar.one_plane_header and planar.PlanarConfigError (a ValueError),
   and the bioio header key spacing_yx_um ((dy, dx) or None; one plane only).
3. (2026-10-24) the command line's series refusal form of D45.
4. (2026-10-24) germquant.io.labels.segment_hashes and UNREADABLE_HASH (D46).
5. (2026-10-24) the batch resume call of D44.
6. (2026-10-24) the project package as built: `project unlock --json`, Project.config_hash, ProjectError and
   ProjectLockedError, preview_roles' record, `project edit --dry-run`, the refusal form, the word none that clears
   a setting; (2026-10-31) `project add` and `project new --images` refuse an io.input_suffixes that the recipe
   check refuses.
7. (2026-10-24) tune-spots as built: the apply form, the --json records, the exit codes, tuning/<run>/run.json and
   TuningWarning; (2026-10-30) a refusal is one ERROR line and nothing on stdout.
8. (2026-10-24) pack derive as built: germquant.packs' new names, the --json record, derived_from, and the "projects"
   feature entering with run_project.
9. (2026-10-24) the post-stage package's additions: run_post_stages(resume=), list_post_stages(include_unusable=),
   the stage records' inputs, `post_stages list --all`, the write record and its exit codes.
10. (2026-10-24) project runs: D47 to D49.
11. (2026-10-30) import-labels and train-curve as built: channel=, the import record, the pairs manifest's new
    columns, tiles above 1024 px, pixels_sha256 in little-endian byte order, --frame, the exit codes;
    (2026-10-30) train-curve's refusal form; (2026-10-31) D56's export rule, and a channel named by digits is given
    by its index on the command line.
12. (2026-10-30) germquant.api re-exports germquant.tune's TuneRefused, GridError and TuningWarning.
13. (2026-10-30) the sample sheet and batch formats as built: the module's names, the row hash over the row's
    non-empty design cells, the resume check's sample row reasons, the separator, root and duplicate_id rules, the
    refusal forms of --samples and batch.discover_options; (2026-10-31) one image id twice in one format is a
    duplicate_id too, and germquant.api re-exports batch.discover_options and the sheet's SampleSheetError and
    SampleSheetWarning.
14. (2026-10-30) planar.write_labels_yx and planar.write_spots_yx, and the stride reading of Stack.spacing_yx.
15. (2026-10-31) project runs with a sample sheet: collect_project's "sample_sheet", the sheet read once per run, a
    sheet that cannot be used, the series refusal with the project's remedy, reason_code sample_row; (2026-10-31) a
    series that finished before the sheet gave it a row reads stale.
16. (2026-10-31) the post-stage wiring as built: `--post-stages FILE` on run and batch, the list in force (the file
    given, else a project's post_stages.yaml, else the pack's post_stages/<recipe>.yaml), checked before any image;
    api.list_post_stages and post_stages= on api.run and run_loaded; a failed post-stage, and an error that leaves
    the post-stage runner itself, never fail the image; one series refusal function for `run --t` and `project run`,
    which `project status` asks too, with the project's post-stage list; collect's batch_post_<name>_<table>.csv of
    the post-stages in force; the post_stages feature.

## 15. Findings to raise (number-relevant; not fixed in T5)

- a. Since Cellpose 4.2.0 (2026-06-11), the bare name cpsam runs cpsam_v2. nuclei._cellpose passes no
  pretrained_model for cpsam; the installed 4.2.1.1 has MODEL_NAMES [cpsam_v2, cpdino, cpdino-vitb, cpsam] and
  defaults to cpsam_v2; run_manifest records model_sha256 "absent". config/profiles/generic.yaml therefore runs
  cpsam_v2 today. Fixing the branch changes generic-profile masks; C. elegans runs load their model by path and do
  not move. T5 leaves it, and packs pin weights by name and digest. The same doubt applies to any "stock cpsam"
  figure (for example F1 0.17). If that figure's run date is known, the cache dates (cpsam_v2 on 2026-06-12,
  cpsam v1 on 2026-09-24) may tell which weights it used; the state-dict comparison of section 3.11 cannot.
- b. A non-Cellpose segmentation backend still hashes segmentation.nuclei.cellpose_model (the 1.2 GB germline
  model under config.yaml) as model_sha256. Pack recipes avoid it; a general fix is a later, deliberate change.
- c. The nd2 reader reads position 0 of a multi-position file without a flag. T5 refuses such files in projects
  only.
- d. The base weights of germline_nuclei_combined are not recorded: scripts/train_nucleus_model_combined.py builds
  CellposeModel(gpu=True), whose default depends on the Cellpose version installed at training time. The file dates
  point to cpsam_v2 (section 3.11), and P1's state-dict comparison tests it before the card is written.
- e. The workstation's plugin venv and engine .venv run germquant from two different checkouts
  (C:/Users/ryane/germquant_wt_dev and C:/Users/ryane/C_elegans_germline_masking_and_spot_counting), and neither
  is the newest one. Until O2 is settled, a GUI run and a CLI run may not even use the same code. The plugin's
  runner therefore asks for the engine interpreter and checks its version and git sha (section 8.1).
- f. A role whose name does not match binds by its index fallback with only an advisory flag, and an unresolved
  segmentation role makes the segment stage segment channel 0 (pipeline.py, the segment closure). Either way a
  fluorescence channel can be segmented as cells, and in the first case qc_pass stays True. Outside packs opt out
  through segmentation.nuclei.require_channel and a map without an index fallback for the segmentation role; the
  default is left as it is.

## 16. What changed from the three proposals

Kept from the risk-first design: the fence before any feature; packs as static data; LF-normalised SHA256SUMS and
.gitattributes; immutable, versioned installs; versioned pack: lines; an isolated planar module whose gating
refuses instead of skipping; conditional provenance only; the tifffile seam for CI; the subprocess runner; the
never-matching filename regex; planar-aware resume; separate environments for third-party code; the AST seal.

Grafted from the MVP-first design:
- the milestone order, with a codeless outside run (M1) before any tuning;
- the SpotMAX knobs spot_radius_z_um, remove_hot_pixels and sharpen, outside SPOTMAX_DEFAULTS;
- the sharpest plane selector, expects ranges, typed tunables with units;
- pixels_sha256 for examples, and the HTTP-range zip-member fetcher;
- the pack record in the batch-root manifest, and the device actually used;
- digest-qualified model names (`<name>@<8 hex>`), and adoption of the Cellpose cache;
- the sample sheet applied per image, outside the hash, compared by is_done;
- the final-score ledger, gate numbers put to the PI first, and the no-stacked-PR rule.

Grafted from the paper-first design:
- the claim map and paper artefacts (parity and device reports, PAPER_REPRODUCE.md, DOIs, CITATION.cff);
- the content-addressed blob store, and an object diameter in micrometres;
- one shared scorer, the DeepFoci self-check and a strict 3D score next to the lenient one;
- the final-score gate on a committed protocol, and contamination checks;
- labels_from (with the pickle guard and per-image resume hashes), learning curves, the timed session, and the
  spot tuning panel as a core deliverable;
- volume_um3 NA for projected runs, and one single-plane rule for both readers;
- pack lint, including --previous, and module checks with find_spec.

Corrected:
- recipe.yaml exists only after an edit, so an unmodified project hashes exactly like --pack.
- pack status checks modules with find_spec, because the built-in cellpose name is always declared.
- Fence digests use LF-normalised bytes; rad51_split.json and the benchmark harness are guarded by a diff check,
  not pinned, because the RAD-51 work may legitimately change them.
- The yeast map covers the SUT509 and mtDNA names (Bright, phase_contr, mNeon) and leaves out BF, which matches
  BFP. MDN1's gene-named channels are renamed by the adapter (section 17).
- check_planar covers projected runs, and germline and axis are not planar_ok.
- sharpen feeds detection, and the feature step gets both images, as in SpotMAX core. Section 17 corrects the
  reason: SpotMAX computes the effect size on the sharpened image.
- The plugin pins >=0.2.0.dev0,<0.3 during T5.
- The plugin's environment facts (no SpotMAX, CPU torch, a different checkout) make the subprocess runner part of
  M1, not a later option.
- T5 waits for tranche 4 and the SpotMAX-grid work to reach main.

Not grafted: channel bindings outside the hash (is_done cannot see them), slash-qualified pack model references
(they collide with the inherited-path re-anchoring in load_config), moving filter_effect_size out of
benchmark/tuning.py during an active protocol, and byte-level locks that differ between Windows and Linux.

## 17. Revision after the completeness critic

A completeness critic checked the draft against both repositories and reported 24 missing items and 9 wrong ones.
Each was checked again against the code, the deposit file list or a scratch test before it was applied. All 33 held
up. Where a proposed fix did not work or clashed with a rule, another fix was used; those departures are listed
last.

Added (missing items):
1. .gitignore negations and a tracking test, so pack cards and examples reach every clone (3.7, P1, G7, risk 23).
2. Classical runs leave `<name>@<8 hex>` unresolved, so CI recipe goldens need no blob (3.5, 6.9, G8, P1, P3).
3. A minimal post-stage contract as package P5b, and a PI question on its timing (5.6, P5b, P6, O18, risk 27).
4. A macOS CI job, a Mac constraints file for the pack extras and a real-Mac smoke test at M3 (G1, P1, P8, O8).
5. Plugin CI on three systems, starting with a GitHub repository for the plugin (8.10, P4, O22).
6. Packs and the readers extra in the container, a smoke step, and publishing tied to O7 (G13, P1, P3, O7).
7. Exact golden gates: --strict-ids, a file-list check, P0 pins for columns, label TIFs, montage (11.1, 11.2, P0).
8. One log per run under runs/; only interactive commands write the manifest, under a lock (7.2, 7.6, P5, risk 24).
9. Projects keep a copy of their pinned pack version, which survives an engine update (3.8, 7.8, P1, P5, risk 19).
10. Portable stage hashes for the new extends forms, so edited projects stay done on another machine (3.4, P1, P5).
11. A first-run interpreter prompt, a version and git-sha check, and per-recipe interpreters (8.1, P4, 15e).
12. Install route: an editable checkout for 0.2.0, built-in roots only in a checkout, O21 on PyPI (3.8, P9).
13. MDN1 Position_4 detector-only, no transmitted index fallback, require_channel (3.12, 6.1, 6.5, P2a, 15f).
14. Engine, plugin and benchmark deposit DOIs, a key resources table and a limitations list (1.2, 1.5, P9).
15. A non-developer for the timed session and an outside pilot lab before 1.0.0 (1.5, 10, P8, P9, O20, risk 28).
16. A PI question on the departures from the proposal's yeast and backend plan (5.5, O19).
17. train-curve --keep-models, a disk and GPU budget, and portable card paths in pack derive (10, P5, P9, O6).
18. Persisted blob verification and the torch-free entry `python -m germquant.packs` (3.5, 3.9, 8.1, P1, P4).
19. The stacked precondition branches and the order in which they merge (1.6, O1).
20. A state-dict comparison to settle the base weights of germline_nuclei_combined (3.11, 15d, O5).
21. Outside QC settings: expected_from_germ_cell false, the framing check off, qc.min_nuclei (3.12, P2a, P3).
22. The scope check names the two frozen synthetic golden folders instead of a glob (11.2).
23. text=auto attributes with binary rules under packs/ (3.7, P1, risk 11).
24. Outside packs extend a frozen pack_base.yaml, so the config.yaml hash bump never moves them (3.12, P3, risk 1).

Corrected (wrong items):
1. SpotMAX computes the effect size on the sharpened image; tuning and cards follow from that (4.7, 6.8, 16, P2a).
2. The yeast map covers SUT509 and mtDNA, not MDN1; the adapter renames the MDN1 channel (3.12, 6.5, 16, P7).
3. $GERMQUANT_MODELS can shadow a project model; run_project refuses a digest mismatch (7.5, 10, P5).
4. BABY's 5-plane stacks go to the projected recipe; cells_2d gets planes from the adapter (6.1, 6.5, P7, P8).
5. Only germquant.packs can stay free of torch; R5, 3.13 and G6 now say so (P1 acceptance 8).
6. nd2_reader.py gains one Stack field; read_stack and read_nd2_metadata stay untouched (4.10).
7. The yeast recipe drops spots.channel, and lint refuses keys the engine does not read (3.7, 3.12, P1).
8. Done markers are compared on every key except finished (9, P1, P4, P5, P6).
9. Card fields follow the keys train_nuclei really writes, model_size and model_file among them (3.6, P1).

Where the revision departs from the critic's fix:
- .gitignore: a file-level negation fails, since git never re-includes a file under an excluded folder (tested).
- require_channel is a recipe key, not a test of cfg.pack_record, so `--pack` and `--config` behave alike.
- The montage and label pins run in CI only if they agree across platforms; otherwise on the workstation.
- There is no `--version` flag, and `-m germquant` loads torch, so the check runs `-m germquant.packs --version`.
- Post-stages get no local modules folder (R5) and no pydantic (a new dependency in the validated environment).
- The hash bump is handled with a frozen pack base, not a policy of MAJOR bumps after each bump release.
- The state-dict check cannot tell which weights the "stock cpsam F1 0.17" figure used; its run date might (15a).
- spots.channel is removed and linted rather than turned into a new spots.role key.
- Package data for PyPI is left to O21 instead of being built now.

Effect on the schedule: about 13.5 working days, most of it in P1 (2 days), P5b (4 days) and P9 (1.5 days), with
smaller additions elsewhere (section 1.4). No fix moves a C. elegans number: every new behaviour sits behind a new
key, a new extends form, a new command or new tooling.

### 17.1 Verifier pass

A verifier checked the revision against both repositories and the deposit file list. Each issue was checked again
before it was applied:
1. Post-stage settings moved out of hashed recipe text into post_stages.yaml, so a change reruns only the
   post-stage (3.2, 5.6, 7.1, 7.3, 8.5, 9, P5b, P6).
2. Post-stages run once per image, through api.run_loaded and cli._batch; time series refuse them (5.6, 12.2).
3. The cross-root tests move the engine config root and compare every stage hash (P1 and P5 acceptance).
4. The bare-name extends fallback searches only the built-in roots, never the working directory (3.4, 3.8, P1).
5. Config gains anchors next to pack_record (3.4, P1).
6. The key lint reads config.yaml itself, now that outside chains never reach it (3.7, 3.12, P1, P3).
7. collect_project passes the framing factor; the plain `germquant collect` keeps 1.8 (3.12, 7.8, P5).
8. The spotmax_gt adapter maps each truth layout; segmInfo.csv exists for MDN1 and SUT509 only (6.5, P7).
9. The Docker readers pins name the four bioio packages the engine .venv holds (P3).
10. The macOS CI job gets a constraints file with the pytest, ruff and laptrack pins (G1, P1).
11. The pin check runs after the read and before segmentation, and the text now says so (3.5, P1).
12. Every outside recipe sets require_channel; process_image has five hunks; 4.14 names read_stack (3.10, 4.14).
13. Moving P5b names what moves with it (1.4, risk 22); P8 depends on P4; G9 and G12 apply from Train C (12.1,
    12.2, R7).
14. O18 asks the PI to accept the two departures from proposal section 3.
15. .gitignore re-ignores caches and .DS_Store under packs/, checked in a scratch repository (3.7, P1).
16. M1 takes about 4 weeks and the total about 14 to 15 (0, 1.4); three long sentences were split.

Where this pass departs from the verifier's fix:
- config.yaml itself sets spots.channel (line 66), so a lint that accepts every config.yaml key would let it pass.
  packs.lint therefore keeps a deny-list, pack_base.yaml is written without deny-listed keys, and the rule skips
  the C. elegans pack, whose recipes are the frozen files.
- The SUT509 yagya_gt_manual.csv files have not been read either, so 6.5 lists them as unread.
- post_stages.yaml also has a pack form and a --post-stages flag, so every route can switch post-stages on, and
  run_loaded given t refuses them like the series loop.
- The readers pins go into a new constraints/readers-2026-10.txt, because the image's desktop-2026-09.txt pins no
  bioio.
- The non-editable extends test sits in P1 acceptance 12. In a checkout, config_roots() tries the checkout first,
  so only the non-editable case tells the old and new fallback apart.
