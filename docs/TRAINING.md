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
segmented (`segmentation.nuclei.channel`, through the run's channel map). In a batch tree the run folders hold
no `run_manifest.json`: every run of a batch, 3D, projected or single-plane, then takes its role from the batch's
own `run_manifest.json` (the one above the folder that records the run's config hash). A run without such a
manifest is exported from the `dna` role (channel 0, with a warning, when the run's channel map has no `dna` role).
What it writes follows the run's geometry (`--orientations auto`, the default):

| the run | what export writes | file names |
|---|---|---|
| a 3D run (every C. elegans recipe, `nuclei_foci_3d`) | for each of xy, xz and yz, `--slices-per-image` slices (default 3) drawn at random from those that cut at least `--min-labels` nuclei (default 3); xz and yz slices are stretched along z to the xy pixel size, so every slice is isotropic and a model trained on them runs in 3D with the anisotropy dz / dy | `<id>_xy_0012.tif`, `<id>_xz_0040.tif` |
| a projected run (`segmentation.nuclei.project`, docs/PLANAR.md; the default recipes of the yeast and cultured cells packs) | one pair of orientation `proj`: the 2D image the segmenter saw (the plane it chose, or the maximum or mean projection) with the run's 2D labels. The image is not resampled: it is the one the model saw | `<id>_proj_0003.tif` (the plane), `<id>_proj_0000.tif` (a projection) |
| a single-plane run (`io.planar`, docs/PLANAR.md; `config/profiles/generic_2d.yaml`) | the plane the run segmented, as orientation `xy`, with the run's labels. A plane of at most 1024 pixels a side is one pair; a larger plane is cut into the fewest tiles of at most 1024 pixels a side, all of nearly the same size (on each axis they differ by at most one pixel), and a tile without any object is left out. Not resampled, like a proj pair | `<id>_xy_0000.tif` (the plane, or tile 0), `<id>_xy_0003.tif` (tile 3) |

Each pair is `image/<name>.tif` (uint16, the raw intensities; a mean projection is rounded to whole numbers) and
`labels/<name>.tif` (int32, the run's labels). `manifest.csv` lists every pair: the image id, series and
timepoint, orientation, slice index (the plane for a projected run, the tile's number for a single-plane run),
number of labels, zone, the two file paths, the source file, the stride and z range, the segmented channel role
and index, the pixel size, the z stretch and the size. A manifest that holds pairs of a single plane has two more
columns at the end, `tile_y0` and `tile_x0`: where the pair's top left pixel sits in the plane, on the run's grid
(0 and 0 for a whole plane); they are empty for the other pairs, and a manifest of 3D and projected runs only has
exactly the columns above. Tiles are numbered row by row: a plane of 2100 rows and 1500 columns gives tiles 0 to 5,
three rows of two. `export.json` records the settings and the runs; its `orientations` is the union over the runs,
and when a run is not 3D, `run_orientations` gives each run's geometry and orientations (with the projection
record of a projected run, and the tile size and where the pixel size came from for a single-plane run).

The pixel size of a pair is the one the run used. A single-plane run, and any run with `io.voxel_size_um`
(docs/PLANAR.md), records `spacing_source` in its run_manifest.json: its spacing is not the file's, so the export
takes it from the run's image_summary table (voxel_dz_um, voxel_dy_um, voxel_dx_um), and the xz and yz stretch and
`pixel_um` follow the run.

- `--orientations xy,yz` names the slices of the 3D runs, in that order (the order decides the seeded draw). A
  projected run gives its `proj` pair whatever the list says: a slice of its labels, which are the 2D labels
  copied through every plane, is no image any model saw. A single-plane run gives its plane (or tiles) whatever
  the list says.
- The draw is seeded per image (`--seed`, default 0), so adding a run to the folder does not change the slices of
  the others. Exporting again into the same folder replaces the pairs of the images it exports and keeps the rest.
- `--zone-column zone` (the default) tags each pair with the most frequent zone of its nuclei, when the nuclei table
  has that column (the staging stage's); the hold-out metrics are then reported per zone.
- `--max-gb` (default 2) stops the export once that many GB are written.

Before this release the default was `xy,xz,yz` for every run, so a projected run gave xy, xz and yz slices of its
copied labels; it now gives its proj pair. A 3D run gives the same files and the same export.json as before.
Single-plane runs are new in this release.

## 3. Or import labels made elsewhere

Labels made in the Cellpose GUI, Cell-ACDC, ilastik, Labkit or micro-sam can join the same layout, one labels
file and the image it was drawn on at a time:

Windows:
```
.venv\Scripts\germquant.exe import-labels "C:\path\to\LABELS.tif" --format tiff --image "C:\path\to\IMAGE.ome.tif" --pairs "C:\path\to\PAIRS_FOLDER"
```
Mac:
```
.venv/bin/germquant import-labels /path/to/LABELS.tif --format tiff --image /path/to/IMAGE.ome.tif --pairs /path/to/PAIRS_FOLDER
```

| `--format` | the file | read with |
|---|---|---|
| `tiff` | a TIFF or OME-TIFF label image (`.tif`, `.tiff`, `.ome.tif`, `.ome.tiff`): ilastik, Labkit, micro-sam, Fiji | tifffile, with the file's own axes (docs/LABELS.md) |
| `png` | a PNG label image of one plane, grayscale or with indexed colours (the index is the object's id) | Pillow |
| `cellacdc-npz` | a Cell-ACDC `segm.npz`: (Y, X), (Z, Y, X), or with a T axis | numpy, `allow_pickle=False` (docs/LABELS.md) |
| `cellpose-npy` | the masks of a Cellpose GUI `_seg.npy` | only with `--allow-pickle` |

The file name must carry its format (a `.npz` for `cellacdc-npz`, and so on). The labels hold the image's own grid
at full resolution: one plane for an image of one plane, or every plane of the image. They are relabelled 1 to N
in the order of the file's ids, and the pairs are written as export-training writes a run's:

- labels of one plane give the plane as orientation `xy`, as a single-plane run's export does (tiles of at most
  1024 pixels a side above that size, with `tile_y0` and `tile_x0` in the manifest);
- labels of every plane give the xy, xz and yz slices export-training draws for a 3D run: 3 per orientation among
  those that cut at least 3 objects, drawn with seed 0 per image, xz and yz stretched to the xy pixel size with the
  image's voxel size;
- the manifest rows carry the source `import:<tool>` in a `source` column (`import:tiff`, `import:png`,
  `import:cellacdc` or `import:cellpose` by the format; from Python, `tool=` names another program, such as
  `micro-sam`). Pairs export-training wrote leave that column empty. The channel role is left empty, because the
  role the labels were drawn for is not known;
- importing an image of the same id again (the id is the file name without its suffix, plus `_t<NNN>` for a
  timepoint of a time series) replaces its earlier rows (and files of the same names), with a warning when they
  came from another image file, such as an image of the same name in another folder; rows of other images,
  exported or imported, stay. Nothing else in the folder changes: export.json still describes the last export.

Options:
- `--channel C` names the image channel the labels were drawn on, by its index or its name. An image of one channel
  needs none; for an image of several channels the default is channel 0, with a warning that names it.
- `--frame N` picks the frame of a labels file with a T axis (docs/LABELS.md, "Frames"); a file without one holds
  frame 0 only. When the image is a time series, N is also its timepoint (default 0, with a warning that names
  it), and the pairs take that timepoint's image id (`<name>_t<NNN>`), so they stay in the series' group when
  training holds series out. A labels file without a T axis (the labels of one timepoint, as the Cellpose GUI
  or micro-sam saves them) is paired with timepoint N.
- `--allow-pickle` reads a Cellpose `_seg.npy`. Such a file is a Python pickle, and loading a pickle can run code,
  so it is read only with this option, only for a file you made, and with a warning. The other formats never
  unpickle anything; a Cellpose mask can also be saved as a TIFF from the GUI and imported as `tiff`.
- `--json` prints the import record: the image id, the digests of the image file, its pixels and the labels file,
  the dims, the channel, the number of objects, the rows written, how many earlier rows of the image were
  replaced, the benchmark splits compared and the warnings.

Test images are never trained on. The import refuses an image that a committed benchmark split holds as test data
(docs/DESIGN_system_packs.md, section 6.6): every `benchmarks/ext/<dataset>/split.json` under the program's folder,
the working folder or `$GERMQUANT_CONFIG_ROOT` is read with its `converted.json`, and an image whose file sha256
or pixels sha256 is that of a test item is refused, whatever its file name or format. The pixels sha256 is taken
of the image as germquant reads it, every channel and plane at full resolution (the timepoint for a time series),
so the same pixels in another file are still recognised (`germquant.io.pixels.pixels_sha256`). The file sha256
reads the whole image file on every import, so each import from a large time series (one timepoint at a time)
reads its file once in full. With no committed split at all, as in a container, which holds no `benchmarks/`
folder, the import warns `no committed splits found` instead of passing silently; a split whose test items carry
no digests, or a split file that cannot be read, is named in a warning too.

It exits with code 2, and writes nothing, when it refuses: a format the file name does not carry, a `_seg.npy`
without `--allow-pickle`, labels of one plane for an image of several (save the plane the labels were drawn on as
an image of one plane and import against it), labels on another grid than the image, labels without any object,
a channel or frame the image does not have, a test image, a missing file. A refusal prints one `REFUSED:` line on
stderr, and with `--json` also `{"error", "kind"}` on stdout. It exits with code 1 when a write fails part way.
docs/LABELS.md describes the label files the engine reads.

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
| `dims` | `3d` for a model trained on slices of 3D runs (or of imported labels of every plane), `2d` on proj pairs and on planes (single-plane runs, imported labels of one plane), `2d+3d` on both | the pairs' orientations, and the tile columns of the manifest |
| `pixel_um` | the pixel size of the pairs, in micrometres; null when they differ, with `pixel_um_values` listing them | manifest.csv |
| `object_kind` | what one label is: `nucleus` or `cell` | `--kind`, else the object kind the packs of all the runs name (two recipes of one pack included), else nucleus. Imported labels name no object kind: give `--kind cell` for labels of whole cells (Cell-ACDC's yeast cells, for example) |
| `system` | the biological system the pairs show | `--system` (free text), else null |
| `pack` | `{id, version, recipe}` of the pack recipe the runs used, with `"edited": true` when that is a project's edited recipe (docs/PROJECTS.md); null for runs without a pack, and null with `pack_note` when the runs differ (an edited recipe differs from its pack recipe) or cannot be read, or when some labels were imported (they come from no run) | the runs export.json lists, read again from their results folders; the manifest's `source` |
| `segmented_role` | the channel role the runs segmented, such as `dna` | manifest.csv |
| `dataset` | for public data only: `{id, licence, split_sha256}` | `--dataset ID`: the committed record benchmarks/ext/ID/ (its dataset.json and the sha256 of its split.json) |
| `source_url` | null: a local model has no download address | |

`export.json` records the last export into the pairs folder only. Export every run at once (one results tree)
for the card to name its pack; pairs exported into the folder earlier, or copied without their results, give
`pack` null with `pack_note`.

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
- It exits with code 2, and trains nothing, when the pairs cannot make a curve (a single group, no size that fits)
  or the base cannot be used (an unknown base, or a pinned base whose weights fail the pin check).

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
to what the export wrote before this change, from a results folder and from a batch tree, checks that a 3D run of a
batch tree whose recipe segments `central_element` (`live_syp2`) gives that channel's slices, as the same run in a
results folder does, and checks the proj pairs of projected runs (the recorded plane, the projections, the run's
stride and z range, a run of a batch tree), the xy pairs of single-plane runs (the whole plane, the tiles of a plane
of 1030 rows and 1100 columns and their offsets, the stride, the size cap, a batch tree), the spacing of runs with
`io.voxel_size_um`, and trees that mix the geometries.
`tests/test_import_labels.py` imports TIFF, PNG (also with indexed colours) and Cell-ACDC npz files (one plane,
every plane, with a T axis, and a time series, also from labels of one timepoint), refuses a `_seg.npy` without
`--allow-pickle`, refuses a test image of a benchmark split by its file and by its pixels, warns without any split
and for a split it cannot read, pins `pixels_sha256`, checks the merged manifest (an image of the same name in
another folder included) and the card of imported pairs, the exit code of a failed write, and checks, in a fresh
Python, that the import loads torch before it reads an OME-TIFF.
`tests/test_train_curve.py` checks the card fields (trained from the test pack's stand-in model `tinynet@f5055f8a`
through the model store, as CI cannot hold the 1.2 GB Cellpose-SAM weights; the same command with
`cpsam_v2@0f1cc3f7` is an opt-in test, `GERMQUANT_TRAIN_CARD_BASE`), the card's pack and object kind from the run
records, and the learning curve: nested, seeded subsets, one held-out set, the scorer seam, `--keep-models` none
and best, a subset too thin to train, and a base refused with exit 2. `tests/test_training_loop.py` covers the
export of time series, the split, the card's training record and the registry. In all of them Cellpose's training
and prediction are replaced by stand-ins, so the suite never trains the network; a real two-epoch training on the
CPU is opt-in, `GERMQUANT_TRAIN_SMOKE=1`.
