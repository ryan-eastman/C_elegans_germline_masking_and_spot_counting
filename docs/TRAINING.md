# The training loop: a segmentation model for your own images

germquant segments nuclei (or, in the yeast pack, whole cells) with Cellpose-SAM. The stock model works on many
samples, and the C. elegans pack ships a model fine-tuned on germline nuclei. When the masks are not good enough
on your images, this loop makes a model from your own corrected labels, and every step records what it used: the
runs the labels came from, the base weights, the training settings and the scores. The loop is the same for
every system (docs/DESIGN_system_packs.md, section 10).

Every command on this page runs in the program's folder, in PowerShell on Windows and in the Terminal on a Mac
(README.md, "Every command in this guide is shown twice"). The commands this page introduces are shown once for
each machine; they differ only in how the program is called and how paths are written. The others are shown once,
as `germquant ...`: type `.venv\Scripts\germquant.exe` for `germquant` on Windows and `.venv/bin/germquant` on a
Mac.

## The nine steps

1. Run the pack on your images, in a project (docs/PROJECTS.md).
2. Export the images the segmenter saw, with its labels: `germquant export-training`.
3. Or import labels made in another program: `germquant import-labels`.
4. Correct the labels, in napari's Training dock or the Cellpose GUI.
5. Train: `germquant train-nuclei`.
6. Measure how much data the model needs: `germquant train-curve`.
7. Use the model in the project: `germquant project use-model` (docs/PROJECTS.md).
8. Promote it into a pack of your own: `germquant pack derive` (docs/PACKS.md).
9. Spots are tuned, not trained: `germquant tune-spots` (docs/TUNING.md).

The time-lapse page (docs/TIMELAPSE.md, "Segmenting without DAPI") walks through the same loop for the live
GFP::SYP-2 series of the lab.

## 1. Run the pack

Run the recipe on a few representative images first, as a project (docs/PROJECTS.md) or with
`germquant run IMAGE --pack ID --recipe R --out results/IMAGE`. The export in step 2 reads every finished run:
a folder with `run_manifest.json` (or `<id>__stages.json` in a batch), `<id>__nuclei_labels.tif` and
`<id>__image_summary.csv`. Pick images that cover your conditions, the dimmest and the most crowded included.

## 2. Export the images the segmenter saw

Windows:
```
.venv\Scripts\germquant.exe export-training "C:\path\to\RESULTS_FOLDER" --out "C:\path\to\PAIRS_FOLDER"
```
Mac:
```
.venv/bin/germquant export-training /path/to/RESULTS_FOLDER --out /path/to/PAIRS_FOLDER
```

The first path is a results folder, a time-lapse series folder or a batch tree. For every finished run under it,
export-training reads the image again from its source file at the xy stride, z range and timepoint the run
recorded, so the exported images sit on exactly the grid the labels were made on, and takes the channel the run
segmented (`segmentation.nuclei.channel`, through the run's channel map). What it writes follows the run's
geometry (`--orientations auto`, the default):

| the run | what export writes | file names |
|---|---|---|
| a 3D run (every C. elegans recipe, `nuclei_foci_3d`) | for each of xy, xz and yz, `--slices-per-image` slices (default 3) drawn at random from those that cut at least `--min-labels` nuclei (default 3); xz and yz slices are stretched along z to the xy pixel size, so every slice is isotropic and a model trained on them runs in 3D with the anisotropy dz / dy | `<id>_xy_0012.tif`, `<id>_xz_0040.tif` |
| a projected run (`segmentation.nuclei.project`, docs/PLANAR.md; the default recipes of the yeast and cultured cells packs) | one pair of orientation `proj`: the 2D image the segmenter saw (the plane it chose, or the maximum or mean projection) with the run's 2D labels. The image is not resampled: it is the one the model saw | `<id>_proj_0003.tif` (the plane), `<id>_proj_0000.tif` (a projection) |

Each pair is `image/<name>.tif` (uint16, the raw intensities; a mean projection is rounded to whole numbers) and
`labels/<name>.tif` (int32, the run's labels). `manifest.csv` lists every pair: the image id, series and
timepoint, orientation, slice index (the plane for a projected run), number of labels, zone, the two file paths,
the source file, the stride and z range, the segmented channel role and index, the pixel size, the z stretch and
the size. `export.json` records the settings and the runs; its `orientations` is the union over the runs, and
when a run is not 3D, `run_orientations` gives each run's geometry and orientations (with the projection record of
a projected run).

- `--orientations xy,yz` names the slices of the 3D runs, in that order (the order decides the seeded draw). A
  projected run gives its `proj` pair whatever the list says: a slice of its labels, which are the 2D labels
  copied through every plane, is no image any model saw.
- The draw is seeded per image (`--seed`, default 0), so adding a run to the folder does not change the slices of
  the others. Exporting again into the same folder replaces the pairs of the images it exports and keeps the rest.
- `--zone-column zone` (the default) tags each pair with the most frequent zone of its nuclei, when the nuclei table
  has that column (the staging stage's); the hold-out metrics are then reported per zone.
- `--max-gb` (default 2) stops the export once that many GB are written.

Before this release the default was `xy,xz,yz` for every run, so a projected run gave xy, xz and yz slices of its
copied labels; it now gives its proj pair. A 3D run gives the same files and the same export.json as before.

Single-plane images (`io.planar`, docs/PLANAR.md) are not exported by geometry in this build yet: their xy pairs,
tiled above 1024 pixels with the offsets recorded, come with the single-plane path.

## 3. Or import labels made elsewhere

Labels made in the Cellpose GUI, Cell-ACDC, ilastik, Labkit or micro-sam can join the same layout:

```
germquant import-labels SRC --format tiff|png|cellacdc-npz|cellpose-npy --image IMG --pairs PAIRS_FOLDER
```

This command is not in this build yet. When it comes, it reads the labels file, refuses a Cellpose `_seg.npy`
file without `--allow-pickle` (such a file can run code when it is opened), refuses an image that a committed
benchmark split holds as test data, and writes `image/`, `labels/` and manifest rows with the source
`import:<tool>`. docs/LABELS.md describes the labels files the engine reads.

## 4. Correct the labels

Open a pair in napari's Training dock, or in napari or the Cellpose GUI by hand (the image, then the labels as a
labels layer), fix the merges and splits, and save the labels back under the same name in `labels/`. Delete a pair
that is not worth correcting (both files); training skips manifest rows whose files are gone. About 30 corrected
pairs across your conditions is a first training set.

## 5. Train

Windows:
```
.venv\Scripts\germquant.exe train-nuclei --pairs "C:\path\to\PAIRS_FOLDER" --out "C:\path\to\MODELS\my_nuclei" --base cpsam_v2@0f1cc3f7
```
Mac:
```
.venv/bin/germquant train-nuclei --pairs /path/to/PAIRS_FOLDER --out /path/to/MODELS/my_nuclei --base cpsam_v2@0f1cc3f7
```

train-nuclei runs Cellpose's own training in process (150 epochs, AdamW, learning rate 1e-5, weight decay 0.1,
batch 8), on the GPU when one is visible (`--cpu` keeps it on the processor). The pairs of whole groups are held
out for scoring, never single slices: `--split-by series` (the default) keeps every timepoint of a series on one
side, `--split-by image` holds out whole images, `--val-frac` (default 0.2) sets the share, `--seed` the draw.

`--base` names the starting weights: a pinned pack model such as `cpsam_v2@0f1cc3f7` (the verified file of the model
store, checked against its pin as a run checks it; `germquant pack fetch-models` puts it there), `cpsam` (Cellpose's
default weights), a registered model name, a model file, or `none` for a randomly initialised network (tests
only).

It writes, in `--out`: the model file (`--name`, default the folder's name; a Cellpose-SAM fine-tune is about
1.2 GB), `training_log.json` (losses per epoch, the split, counts, settings, times) and the model card
`<name>.card.json`. The card holds the base model with its sha256, the sha256 of the weights and of every pair, the
training settings, the hold-out metrics (F1 at IoU 0.5, precision, recall, miss rate and mean matched IoU, overall,
per orientation and per zone, and per pair), and what the data was:

| card field | what it says | where it comes from |
|---|---|---|
| `kind` | the model type, `cellpose` | |
| `dims` | `3d` for a model trained on slices of 3D runs, `2d` on proj pairs, `2d+3d` on both | the pairs' orientations |
| `pixel_um` | the pixel size of the pairs, in micrometres; null when they differ, with `pixel_um_values` listing them | manifest.csv |
| `object_kind` | what one label is: `nucleus` or `cell` | `--kind`, else the pack of the runs, else nucleus |
| `system` | the biological system the pairs show | `--system` (free text), else null |
| `pack` | `{id, version, recipe}` of the pack recipe the runs used; null for runs without a pack, and null with `pack_note` when the runs differ or cannot be read | the runs export.json lists, read again from their results folders |
| `segmented_role` | the channel role the runs segmented, such as `dna` | manifest.csv |
| `dataset` | for public data only: `{id, licence, split_sha256}` | `--dataset ID`: the committed record benchmarks/ext/ID/ (its dataset.json and the sha256 of its split.json) |
| `source_url` | null: a local model has no download address | |

The option is `--kind`, not `--object-kind`, so that `--o` stays short for `--out`. A wrong `--kind` or dataset id
is refused before anything is read or written. The card's `model_file` and `base_model.path` hold this machine's
paths; `pack derive` (step 8) rewrites them before a card is published.

## 6. Measure how much data the model needs

A learning curve trains the same recipe on more and more of your pairs and scores every model on the same
held-out pairs, so you can see whether more corrections would still help.

Windows:
```
.venv\Scripts\germquant.exe train-curve --pairs "C:\path\to\PAIRS_FOLDER" --out "C:\path\to\CURVE_FOLDER" --sizes 1,2,4,8,16,all --repeats 3 --seed 0 --base cpsam_v2@0f1cc3f7
```
Mac:
```
.venv/bin/germquant train-curve --pairs /path/to/PAIRS_FOLDER --out /path/to/CURVE_FOLDER --sizes 1,2,4,8,16,all --repeats 3 --seed 0 --base cpsam_v2@0f1cc3f7
```

- The held-out pairs are the ones train-nuclei holds out with the same `--seed` (whole series, 20 percent). Every
  model of the curve is scored on them, so a curve needs at least two groups.
- `--sizes` counts training groups (source images, or series): `all` is every training group, and a size above
  that is skipped and recorded. For each of the `--repeats` the training groups are put in a seeded order and the
  subset of size k is the first k of it, so within one repeat every subset holds the smaller ones (nested), and
  the same seed gives the same subsets.
- Each model trains from `--base` as train-nuclei trains, with its defaults, into
  `CURVE_FOLDER/models/n<size>_r<repeat>/`, and is scored right after training: the score is the hold-out F1 at
  IoU 0.5 that train-nuclei computes (scorer `holdout_metrics`). Programs can pass their own scorer to
  `germquant.segment.training.train_curve` (`scorer=`); the external benchmarks will plug their end-to-end
  scoring in there, and its name is recorded.
- `--keep-models none` (the default) deletes each model's weights once it is scored; `--keep-models best` keeps
  only the model with the best score (a tie keeps the earlier, smaller one). A 6-size, 3-repeat curve trains 18
  models, about 22 GB if every Cellpose-SAM fine-tune were kept. Every card (with its metrics, its weights' sha256
  and a `curve` record: size, repeat, groups, scorer, score, whether the weights were kept) and every
  `training_log.json` stays.
- `curve.csv` has one row per model: size (and its label, `all` for every group), repeat, the groups, the numbers of
  training and held-out pairs, the score, F1, precision, recall, miss rate, mean matched IoU, the weights' sha256 and
  size, whether they were kept, the model and card paths, the time, and the status (`scored`, or `failed` when every
  slice of a subset holds too few nuclei to train). `curve.json` holds the plan (sizes, subsets, held-out groups,
  scorer, settings), every row and the best model.
- Every model is a full training (150 epochs). The command uses train-nuclei's defaults; from Python,
  `train_curve(..., epochs=..., split_by=..., val_frac=..., cpu=..., min_train_masks=...)` changes them.
- It exits with code 2, and trains nothing, when the pairs cannot make a curve (a single group, no size that fits,
  an unknown base).

## 7. Use the model in the project

Train into the project's round folder (`--out PROJECT/training/r01/models/my_nuclei`), then:

```
germquant project use-model PROJECT --model training/r01/models/my_nuclei/my_nuclei
```

(`--card CARD` names its card.) This writes the model's relative path into the project's recipe.yaml and its
sha256 into the project file, so the project gets a new config hash and is marked not validated; a project run
refuses a model file whose sha256 differs (docs/PROJECTS.md). Outside a project,
`germquant register-model --model MODEL --card CARD` adds the model to the local registry, and a profile names it
in `segmentation.nuclei.cellpose_model` (docs/TIMELAPSE.md, "Register and use").

## 8. Promote it into a pack

```
germquant pack derive --from ID@VER --recipe R --model PATH --card CARD --out PROJECT/packs
```

writes a derived pack: a recipe that extends the parent recipe by its versioned name, the model spec (name, sha256
and card) and checksums. The card's machine paths become relative names with their sha256, so no path of this
computer is published (docs/PACKS.md).

## 9. Spots are tuned, not trained

Spot detection has no model to train here. Instead, mark the spots on a few images (napari's spot tuning dock
writes `annotations/<image_id>__points_<annotator>.csv`) and let `germquant tune-spots --project PROJECT` search the
recipe's spot settings for the best match; applying the result writes recipe.yaml (docs/TUNING.md).

## What is tested

`tests/test_training_export_geometry.py` pins the 3D export of the synthetic gonad (manifest.csv and export.json)
to what the export wrote before this change, and checks the proj pairs of projected runs (the recorded plane, the
projections, the run's stride and z range) and a tree that mixes both. `tests/test_train_curve.py` checks the card
fields (trained from the test pack's stand-in model `tinynet@f5055f8a` through the model store, as CI cannot hold
the 1.2 GB Cellpose-SAM weights; the same command with `cpsam_v2@0f1cc3f7` is an opt-in test,
`GERMQUANT_TRAIN_CARD_BASE`), and the learning curve: nested, seeded subsets, one held-out set, the scorer seam,
and `--keep-models` none and best. `tests/test_training_loop.py` covers the export of time series, the split, the
card's training record and the registry. In all of them Cellpose's training and prediction are replaced by stand-ins,
so the suite never trains the network; a real two-epoch training on the CPU is opt-in, `GERMQUANT_TRAIN_SMOKE=1`.
