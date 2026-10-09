# Proposal: a general modular image-quantification tool built on germquant

Status: proposal for discussion with the PI, written 2026-10-09 from a design pass (codebase map, survey of existing tools, three independent proposals, three judges, one synthesis). Nothing here is implemented. The C. elegans pipeline (germquant) stays as it is; the general tool would be a separate repository that depends on it.

## 1. Recommendation

Keep germquant as the engine, unchanged in name, numbers and launchers, and build the outsider-facing layer as a thin napari plugin in a new repository that depends on germquant pinned at a git tag. Two judges chose this (Stance B); the third chose the new-core design (Stance A) because B's "module" is only a toggle on a fixed stage list and generic features would accrete inside the validated repo. I side with B on sequencing (a GUI in three weeks, no numeric code written) and with the dissenter on the gap: the plugin gets A's entry-point stage contract, so an outside lab adds a measurement without touching germquant, which takes only a few additive, golden-gated patches.

## 2. Repo strategy

Fork? Not for the code. A fork copies the numeric code, the golden gate certifies only one copy, and the two drift at the first bug fix. The same repo is wrong too: napari and Qt would enter the pinned environment the numbers depend on. Mechanics if wanted: the Fork button keeps history but links to upstream and cannot go private; an independent copy is `git clone --bare`, an empty repo, `git push --mirror`. napari-germquant starts from the napari plugin template, without that history.

germquant: finish the RAD-51 goldens (today's NAS run is at 2 of 12; do not start a second), fetch, merge the two branch commits, tag v0.1.0. After that, only additive, golden-gated PRs: `api.py` re-exporting the public functions, the private writers and the synthetic stack builder; `read_image` dispatch, nd2 path untouched, bioio behind it for TIFF, CZI, LIF and OME-Zarr; a generic channel map and profile as new files. New keys read code defaults and the frozen config texts stay unedited, so every historical config_hash reproduces.

napari-germquant holds everything else and pins `germquant @ git+...@v0.1.0`.

Byte-identity: only `process_image` writes result folders, so GUI and CLI runs are indistinguishable on disk; previews write no tables; CI diffs the widget path against a direct `process_image` call on synthetic stacks; the 23 real goldens gate every tag, on the 5090 only.

## 3. What a module is for outsiders

Built-in modules are germquant's 18 stages, documented as cards in CellProfiler's words (image and channel roles in; primary objects, per-object measurements, secondary objects, per-image summary and QC flags out), each with one switch key and one provenance hash. Group them as "objects and puncta, any sample" and "C. elegans germline". Role names stay literal in code; the widget labels them objects, puncta, secondary objects, tissue marker and boundary, and the cards flag the germline-sized defaults an outsider must check.

Added modules are post-stages: one class declaring required and produced kinds (Image with role, Labels with object type, Table keyed by id, Points, Mask, Scalars), a pydantic Params model with unit-named fields, and `run()`. Found through entry points or a local modules folder, they run after `process_image` on what it wrote, append under the append-only schema rule, and get their own sub-hash. Zero germquant edits, so the goldens cannot move.

## 4. What "train models" means

Three model classes, trainers wrapped not written; every model gets a sha256 (already in the segment sub-hash and done marker) and is selected through a profile, never by editing config.yaml.

- Nuclei: `cellpose.train.train_seg` (as in scripts/train_nucleus_model_combined.py), in-process, not finetune.py's CLI subprocess. `export-training` writes slice pairs from a run's own labels so users correct rather than draw. Split by source image, never by slice: the June split leaked adjacent crops and the F1 0.98 has no leak-free provenance. Held-out `segmentation_metrics` against the current model, and a model card (bioimage.io fields, training-set sha256s, the split) before `register-model` accepts it.
- Spots: tune before train, with the cv_detect sweep in the Spots widget. SpotMAX 1.3.1 ships an undocumented 3D U-Net trainer; ask its authors before building. Spotiflow (0.6.5, nothing since mid-2025) is the fallback, off until cross-validated against Imaris counts.
- Masks (granules, envelope, SC): ilastik or Labkit outside the plugin, probability maps imported; no trainer written.

## 5. GUI

napari 0.9 as host, npe2 manifest, magicgui forms, `thread_worker` for the 20-minute gonad. PyQt6 as the one binding (micro-sam requires it, Cell-ACDC's gui extra installs it), in a separate hosts venv; the batch .venv is never touched. Every entry point imports `germquant.pipeline` before spotmax or CuPy (5090 cuBLAS).

First screen: File > Open adds one Image layer per channel with physical scale and flags files with T or P above 1. A Run dock: roles table (red when a required role is unbound), profile and stage checkboxes from `cli_switches()` greyed by `missing_roles`, filename-regex preview, output folder, Run. Results load as Labels and Points with tables as layer features. "Export profile" writes an `extends:` overlay marked `validated: false` plus the exact `germquant batch` command. The pachytene tracer becomes a Shapes path layer writing through `save_trace`; `germquant trace` stays until three re-traces match.

## 6. Reuse, and why not just that tool

napari, npe2, magicgui; bioio readers; Cellpose 4 (GUI, train_seg); Cell-ACDC (3D brush, model folder: replace its stale Cellpose_germlineNuclei folder, whose model path no longer exists, with one defaulting to `get_model()`); the SpotMAX GUI (tuning, INI import); micro-sam (3D annotation); Spotiflow; ilastik and Labkit. None of them composes 3D per-object stages and trains the models those stages use while keeping provenance and a regression gate: Cell-ACDC plus SpotMAX is a fixed segment-then-spots flow with no hashes; CellProfiler is 2D-first; napari alone measures nothing; BiaPy trains but has no per-object stages. Skip napari-assistant and cellpose-napari (old napari), napari-spotiflow (pins torch below the lab's), Allen Segmenter ML (Python below 3.11).

## 7. Phases (focused weeks; two semesters at half time)

0. Freeze and scaffold (1): goldens, tag, api.py PR, plugin repo, constraints, CI, pilot lab named.
1. Open and run (2): read_image PR, reader, Run widget, equivalence test. Milestone, week 3: an outside lab opens its TIFF or CZI stack, maps roles, runs the generic profile, inspects layers.
2. Stage widgets, profile export, tracer, cards (2). Gate: continue only if the pilot lab is running it.
3. Cell-ACDC export and import, model folder, SpotMAX INI (2, parallel with 1 and 2).
4. Extension contract, local modules folder, template, tutorial (2).
5. Nucleus training loop and model cards (3).
6. Spot tuning; Spotiflow only if a dataset fails the threshold path (1 to 3).
7. Pilot run by someone else, first real Mac run, 0.2.0 tags, napari hub (2).

Total 15 to 17.

## 8. Risks

Qt clashes (one binding, separate venv). napari API churn (pin the minor version, thin widgets). Cell-ACDC 1.8.1 and SpotMAX 1.3.3 are out while the lab pins 1.7.1 and 1.3.1; re-verify integration points per bump. A grouped retrain may not beat the shipped model (keep it; the loop is still delivered). Mac masks differ (MPS float32), so parity is 5090-only. Scope creep: time-lapse, tracking and whole-slide stay out. After you graduate, germquant must still work without the plugin.

## 9. Decisions to settle with the PI

- Fork, new repo or same repo? New dependent repo; germquant frozen at v0.1.0.
- First outside user? Name a lab with 3D nuclei plus puncta in a non-nd2 format now; phase 1 is built on their files.
- What does "modular" promise? Toggling validated stages and adding a post-stage without touching germquant; not a graph editor.
- What does "train" promise? Fine-tuning nucleus and spot models from corrections with held-out metrics; nothing else.
- Spot-detector training now or later? Later, when a dataset fails thresholds; ask SpotMAX's authors first.
- PyPI, napari hub, tool paper? Decide at the week-5 gate.
- Priority against SC calibration, COSA-1 and Alpine? Phases 0 to 2 first, then reassess.
- Licence and credit for KoehlerLab labels and others' models? MIT like germquant; cards carry provenance.

## 10. Three sentences for the PI

The lab pipeline stays as it is, name, commands and byte-identical numbers included, and becomes the engine behind a thin napari plugin that other labs can install, point at their own 3D stacks, extend with their own measurement step, and use to fine-tune nucleus and spot models from their corrections with held-out evaluation. We write only glue; the viewer, label editing, spot tuning and training come from napari, Cell-ACDC, SpotMAX, Cellpose and Spotiflow, which other people maintain. About fifteen focused weeks; an outside lab can try it after three; the first thing we need from you is that lab's name.

## 11. Clarified scope (Ryan, 2026-10-09)

The target is a general framework applied to other systems, cultured cells and yeast first: a modular interface
where a user selects models and plugins trained or tuned for their system and runs, and which is trainable on
their own data through napari or another annotation tool. That maps onto the design above with three additions.

System packs. A pack bundles a profile, a channel-role map, a nucleus or cell model (checksum plus a card: what
it was trained on, by whom, base model, date), spot-detection presets, size gates and the default module set.
Packs are what the interface lists; choosing one configures everything and the engine runs unchanged. The
existing model registry (get-model) and profiles (extends) are the seed. Planned packs: C. elegans germline
(what exists today), cultured cells (Cellpose general model plus generic presets), yeast (Cell-ACDC yeast models).

Plugins in two senses. Measurement modules: the existing stages grouped as "any sample: objects and puncta" and
"C. elegans germline", plus the entry-point contract for added steps. Segmentation backends: the segment stage
becomes pluggable (Cellpose, Cell-ACDC's model API including its yeast models, bioimage.io models) instead of
Cellpose-or-classical only.

Training loop. In napari: run the pack, correct the labels layer on a few slices or volumes, save them as
training examples, fine-tune from the pack's model with cellpose.train.train_seg (the call
scripts/train_nucleus_model_combined.py already uses), score on held-out images split by source image,
register the model with a card, and it appears in the pack menu. Labels made in the Cellpose GUI, Cell-ACDC,
ilastik, Labkit or micro-sam import as label images. Spot detection is tuned, not trained: annotate a few
puncta as points, sweep SpotMAX parameters, keep the best.

Additions on the germquant side, all additive and golden-gated: readers beyond nd2 (TIFF, CZI, LIF,
OME-Zarr); a 2D path (cell culture is often 2D or thin stacks; the pipeline is 3D-only today); a generic
channel map and profile; the backend interface for the segment stage. The worm stages do not change.

Time-lapse (Ryan, 2026-10-09, later the same day): the lab needs time-lapse quantification for its own
C. elegans live imaging regardless, so time-lapse is in scope as a module set rather than left to Cell-ACDC.
Per-timepoint quantification reuses every existing stage on each timepoint's stack (a T-aware reader and a loop,
both additive, with provenance per timepoint). Tracking is a pluggable backend, not a tracker of our own:
Cell-ACDC trackers, btrack, ultrack, laptrack, or imported TrackMate results, producing a track id per object per
timepoint appended under the schema rule plus per-track time series. The germline-specific part is drift
correction (registration before linking) and validation of linking against a few hand-followed nuclei in a
packed, moving tissue. First dataset: the lab's own live series (for example the October 2026 GFP::SYP-2 5-minute,
2-hour series on the NAS). Yeast and cell-culture labs get time-lapse through the same backends.

Decisions for the PI: which system first, yeast or mammalian culture, and which lab supplies a few stacks; and whether the time-lapse module is built on the lab's own live data before or after the generic fixed-image path.

## 12. Interface requirement (Ryan, 2026-10-09)

The tool must be usable without code and must also run headless. Both paths drive the same engine and write
the same files, so numbers are identical whichever way a run was started.

Codeless: the napari plugin is the control surface for everything: open a project, choose a system pack,
map channel roles, tick modules, run one image or a batch, review labels and points, correct them, train or
fine-tune a model, register it, and open time-lapse results. No terminal.

Headless: the existing CLI and batch commands (and the Snakemake and container routes) run the same project
with one command on a workstation or a cluster.

Project: a folder with the images (or links to them), the chosen pack and profile, the model versions, and
the run log and outputs. The GUI creates and opens projects; a project exported from the GUI runs headless
unchanged; a headless run opens in the GUI for review. The project file is what makes a codeless workflow
reproducible and what a paper's methods section can cite.

Benchmark data: the lab's own annotated gonads (Imaris spot positions from the .ims files, independent of
export version) plus public datasets for other systems (3D nuclei, spots, yeast), to be listed from the
design passes.

## 13. Rewrite or build on: evidence pass (2026-10-09)

### 1. Rewrite Cellpose and SpotMAX? No.

[Cellpose-SAM](https://www.biorxiv.org/content/10.1101/2025.04.28.651001v1.full) was trained on 22,826 images with 3.34 million ROIs from 18 public datasets, 20 hours on eight H200 GPUs per run, by an HHMI team. At 300 to 600 objects per hour (the [original Cellpose](https://www.biorxiv.org/content/10.1101/2020.02.02.931238v1.full) annotation rate) those labels alone are 5,600 to 11,000 hours, and your fine-tune already recalls 0.957 to 1.0 of Imaris-traced nuclei, so there is no metric left to win. [SpotMAX](https://pmc.ncbi.nlm.nih.gov/articles/PMC13524050/) (Science Advances 2026) is 45k lines plus 130k of Cell-ACDC, validated on 20,692 spots, and its C. elegans COSA-1 case made fewer errors than the expert annotator. germquant calls these engines from two files (segment/nuclei.py, 143 lines; spots/detect.py, 212 lines), about 800 of 6,460 lines. Every benchmark in hand (11 byte-identical goldens, 12-gonad Imaris cross-validation) is bound to pinned cellpose 4.2.1.1 and spotmax 1.3.1; a rewrite zeroes it and invites "why not the published tools?"

### 2. What is yours, and how to frame it

The paper is the layer above the engines: SYP-seeded germline isolation (spot recall 0.995, 14 gonads) and gonad linearization; the RAD-51 recipe (axial PSF-split diagnosis of the 1.6x over-count, z-column merge, triangle threshold plus effect size 3; gonad-mean CCC 0.958 HS, 0.805 noHS); the fine-tuned nucleus model and synthetic generator (stock cpsam collapses on orthogonal reslices, F1 0.17 vs 0.66); lamin-envelope masks; the distance-matched partition coefficient with rotation null and z-shift floor (your most original measurement, robust across 18 configs); the Imaris-matched granule recipe (3187 vs 3186); the stage registry, per-stage provenance hashes and byte-identical golden gate. Frame it as [SMINT](https://pmc.ncbi.nlm.nih.gov/articles/PMC12946760/) did: "the primary contribution is the validated computational and experimental workflow"; the journal judges "the innovation and utility of the method itself" ([launch editorial](https://pmc.ncbi.nlm.nih.gov/articles/PMC9017128/)). Report for one readout, Article for the framework with a second system, Resource for the automated-versus-Imaris comparison.

### 3. Where in-house work could pay off

- A 3D-native nucleus backend. Cellpose-SAM is 2D-trained; your fine-tune gets F1 0.98 in xy but 0.66 orthogonal. Painting gonad volumes (about 600 nuclei each; the only public worm 3D set is [larval somatic tissue](https://mlanthology.org/wacv/2020/weigert2020wacv-starconvex)) for StarDist-3D or a 3D flow model could be worth it. Bar: grouped held-out 3D instance F1 above the current model, as [AIstain](https://pmc.ncbi.nlm.nih.gov/articles/PMC12664894/) had to beat Cellpose 3 and SAM2 head to head.
- The four SpotMAX calls, only if licensing forces it: a clean-room few hundred lines on scikit-image, or [Spotiflow](https://github.com/weigertlab/spotiflow) (BSD-3) as backend. Bar: parity on the cv_detect harness against the 12 Imaris gonads.

### 4. Precedents

In Cell Reports Methods: [SCAMPR](https://pmc.ncbi.nlm.nih.gov/articles/PMC9606134/) (Cellpose), [MCA](https://pmc.ncbi.nlm.nih.gov/articles/PMC12946746/) (Cellpose, StarDist), [FLIM Playground](https://pmc.ncbi.nlm.nih.gov/articles/PMC13494553/), [Convpaint](https://pmc.ncbi.nlm.nih.gov/articles/PMC13030958/), SMINT, and [LAM](https://pmc.ncbi.nlm.nih.gov/articles/PMC9017226/), the midgut analogue of your linearization and staging. Elsewhere: [Cell-ACDC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9356409/) (BMC Biology), [TrackMate 7](https://www.nature.com/articles/s41592-022-01507-1) (Nature Methods), [Celldetective](https://api.elifesciences.org/articles/105302) (eLife, "rather than reinventing existing efficient solutions"), [GoNuclear](https://pmc.ncbi.nlm.nih.gov/articles/PMC11273294/) (Development). Occupied ground to position against: your lab's [Toraason 2021](https://academic.oup.com/genetics/article/217/3/iyab010/6117234) and the Koehler [germlineNuclei model](https://pmc.ncbi.nlm.nih.gov/articles/PMC10750166/).

### 5. Gaps

1. Ground truth unpinned: NAS xlsx copies give about half the June coloc counts (HERM_001 3.45 vs 5.81); the pipeline itself reproduces June exactly. Settle and sha256-pin before quoting any correlation.
2. No 3D instance F1 on your own images (all F1 is 2D, on Koehler data).
3. One date, microscope and genotype; es=3 is a knife edge.
4. No head-to-head against Toraason 2021, stock Cell-ACDC plus SpotMAX, or Imaris on the same images.
5. SC tracer synthetic-only; staging manual with no inter-rater statistic; partition biology confounded by session and exposure.
6. Time-lapse, system packs, napari plugin and training loop are proposal text.
7. Hygiene: no CI, LICENSE file, CITATION.cff, release DOI or demo data; the journal requires code and data deposited at submission, a Limitations section and a STAR key resources table.

### 6. Licences

germquant is MIT; Cellpose is BSD-3 (weights trained on CC-BY-NC data; note that in the model card). SpotMAX is contradictory: LICENSE and README say GPL-3.0, pyproject and PyPI say BSD-3. Treat it as GPL-3 until Padovani answers; fine as an optional extra, it bites only when shipped inside your Docker or Apptainer image. Not a barrier for the paper ([u-track3D](https://github.com/DanuserLab/u-track3D) is GPL MATLAB and published there). Cell-ACDC, Spotiflow, scikit-image and trackpy are BSD-3. Remove the unused AGPL pymupdf from the venv. Update the SpotMAX citation to Science Advances; Spotiflow is not dormant (0.6.5, June 2026).
