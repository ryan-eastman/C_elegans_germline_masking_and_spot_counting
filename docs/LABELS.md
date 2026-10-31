# Labels from files: count spots in objects someone else segmented

Some labs already have masks: truth masks of a public benchmark, Cell-ACDC or Cellpose GUI output that was
corrected by hand, an ilastik or Labkit export. With `segmentation.nuclei.labels_from` a recipe reads those
masks instead of segmenting. Everything after segmentation runs as usual: measure, spots, QC, the tables. That
is also how detector-only scoring works: the truth masks go in, and only the spot detector is scored.

Nothing changes for a recipe that does not set the key. Its code default is null, no shipped config sets it, and
a run without it reads, records and hashes exactly what it always did.

The design is `docs/DESIGN_system_packs.md`, section 4.11 (and the decisions D12 to D15 recorded in its section
14). The code is `src/germquant/io/labels.py`.

## Switching it on

| key | values | code default | meaning |
|---|---|---|---|
| `segmentation.nuclei.labels_from` | null, or a path template | null | the labels file of each image |
| `segmentation.nuclei.labels_frame` | a whole number, 0 or more | 0 | the frame of a T axis, read only in labels runs |

A template names one file per image with three fields:

| field | value | example for `D:/data/20260101_cells_001.ome.tif` |
|---|---|---|
| `{image_dir}` | the image's folder, absolute | `D:/data` |
| `{image_stem}` | the file name without its format suffix | `20260101_cells_001` |
| `{image_id}` | the image id (for a fixed image, the stem) | `20260101_cells_001` |

```yaml
# my_study.yaml: masks next to the images, as <stem>_masks.tif
extends: generic.yaml
segmentation:
  nuclei:
    labels_from: '{image_dir}/{image_stem}_masks.tif'
```

Quote the template in YAML (a value that starts with `{` is otherwise read as a mapping). A relative template
(one that does not start with `{image_dir}` or an absolute path) is resolved against the folder of the recipe that
is run: the folder above `config/` for a file in a `config/` or `config/profiles/` folder (in this repository, the
repository root), the file's own folder otherwise. Unlike `io.channel_map`, a relative template that a recipe
inherits through `extends` is not anchored to the base recipe's folder, so a recipe that others extend should use
`{image_dir}` or an absolute path. A template without any field names one file for every image. Any other field,
a format spec such as `{image_id:>9}` or an unbalanced brace is refused before any stage runs; write `{{` and `}}`
for a literal brace.

Labels files kept next to the images, as in the example, are image files too. `germquant batch` (and the
Snakemake workflow) with a recipe that lists their suffix in `io.input_suffixes`, as
`config/profiles/generic_2d.yaml` does for `.tif`, takes each `<stem>_masks.tif` for an image of its own, and so
does `germquant project add` given the folder. That image fails before any stage runs, on every run and every
`--resume`, because its own labels file (`<stem>_masks_masks.tif`) does not exist. For a batch, leave the labels
files out with `io.exclude_patterns`, which replaces config.yaml's list, so keep its entries, for example
`exclude_patterns: ["*largeimage*", "*10x*", "*_10x_*", "*_masks.tif"]`, or keep them in a folder the batch does
not search; for a project, add the images by name or by a glob such as `"D:/data/*.ome.tif"` instead of by folder.

Run it like any other recipe.

Windows:
```
.venv\Scripts\germquant.exe run "C:\path\to\IMAGE.nd2" --config "C:\path\to\my_study.yaml" --out "C:\path\to\RESULTS_FOLDER"
```
Mac / Linux:
```
.venv/bin/germquant run /path/to/IMAGE.nd2 --config /path/to/my_study.yaml --out /path/to/RESULTS_FOLDER
```

## Formats

| format | suffixes | read with | axes |
|---|---|---|---|
| TIFF | `.tif`, `.tiff` | tifffile | the file's own axes |
| OME-TIFF | `.ome.tif`, `.ome.tiff` | tifffile | the file's own axes |
| Cell-ACDC `segm.npz` | `.npz` | numpy, `allow_pickle=False` | (Y, X), (Z, Y, X) or (T, Z, Y, X) |

* A label image holds one integer id per pixel, 0 for the background. Integer, boolean and whole-number float
  files are read; fractions, negative ids, more than one channel and colour images are refused.
* A file holds the image's own grid: one plane, or every plane of the image, at full resolution. A file with one
  plane (2D) is copied through z. `--xy-stride` and `--z-range` cut the labels exactly as they cut the image.
* An npz is read with numpy and `allow_pickle=False`: its first array, as Cell-ACDC saves it. An npz whose array
  holds Python objects is refused, and nothing in it is unpickled.
* A 3D array without axis names (a segm.npz, or a plain TIFF that tifffile wrote without metadata) is a z stack
  when its first axis has the image's number of planes. For an image of one plane it is a series of 2D labels
  (T, Y, X), as Cell-ACDC saves a 2D time-lapse. Any other number is refused: a z stack with planes missing or
  added cannot be told apart from a series. A series of 2D labels for an image of several planes needs a named T
  axis: save it as a (T, 1, Y, X) segm.npz or as an OME-TIFF with its axes. A plain TIFF drops an axis of one
  plane, so (T, 1, Y, X) does not survive in it.
* The ids are relabelled consecutively, 1 to N, in the order of the file's ids, after the cut. There is no size
  gate: `segmentation.nuclei.min_volume_um3` and `min_area_um2` do not apply, so every object of the file is an
  object of the run.

## Frames

A label file with a T axis (an OME-TIFF or ImageJ TIFF with T in its axes, a (T, Z, Y, X) segm.npz, or a
(T, Y, X) segm.npz of an image with one plane) holds one label image per frame.
`segmentation.nuclei.labels_frame` picks one; the default is 0. A frame outside the file, or a frame other than 0
for a file without a T axis, is refused. The frame is chosen with this key rather than in the template (decision
D13).

## What a labels run records

Only a labels run records these:

* `segmentation_method` in image_summary reads `labels_file`;
* the flag `segment:labels_from_file`;
* `labels_from` {path, sha256} in the run geometry (run_manifest.json and `<image_id>__stages.json`);
* the labels file's sha256 in the segment stage sub-hash, so it also reaches the completion marker.

No model is resolved, hashed or pin-checked: the run records `model_sha256` as `unused`, as the classical method
does, and a recipe can name a Cellpose model that is not installed on this machine (decision D12). The
segmentation backend records nothing either. `germquant batch` follows the same rule. The pack check that
`--pack` runs before a run (docs/PACKS.md) is a separate step.

`germquant batch --resume` computes the segment hash per image in a labels run, with that image's labels file,
and for a recipe that sets `io.planar` also the single-plane one, which the completion marker of a single-plane
image carries. Editing one labels file therefore reprocesses that image only; every other image is skipped as
usual. Other runners get the same pair from `germquant.io.labels.segment_hashes` (decision D46). The segment
sub-hash also keeps the segment stage's usual tool versions (Cellpose and torch), although neither runs in a
labels run: upgrading either one reprocesses labels-run images under `--resume` too.

## With the 2D settings

* `segmentation.nuclei.project` (a projected run, docs/PLANAR.md): allowed. The labels file replaces the 2D
  segmentation: a file with one plane is copied through z, and the projected measure and records stay
  (volume_um3 NA, area_um2 and equivalent_diameter_um from the 2D labels, the flag `segment:projected_<mode>`, the
  `segment_projection` record). The projection is computed only for the montage's outline panel: of the
  segmentation channel, or of the first channel when that channel is missing and `require_channel` is false (the
  plane `segment_projection` records, as in a segmenting run). A labels file with more than one plane is refused
  in a projected run (decision D14).
* `io.planar` (a single-plane run, docs/PLANAR.md): allowed. On an image of one plane the labels file holds that
  plane (a 2D file, or a file of one plane), and the run is a single-plane run whose objects are the file's: the
  single-plane measure and records stay (`read:planar`, voxel_dz_um NaN, volume_um3 NA, area_um2 and
  equivalent_diameter_um from the labels, the YX label TIFF, `planar` in the run geometry, the stage sub-hashes
  and the completion marker), next to the labels records. A `segmentation.nuclei.project` setting is ignored
  there with the flag `segment:project_ignored_planar`, as in a segmenting single-plane run (decision D7).
* `segmentation.nuclei.require_channel` does not apply: nothing is segmented, so an image without the
  segmentation channel runs, and no `segment:no_<role>_channel_using_first_channel` flag is raised.
* A segmentation backend that cannot segment a 2D image (nuclei3d) is not refused in a labels run, projected or
  single-plane, because it is never called.

## Refused

Each of these fails the image with a message, before any stage runs and before any file is written:

* the labels file of the image does not exist (the message names the file, the template and the image);
* a `.npy` file, such as a Cellpose `_seg.npy` (see below), or any format but the three above;
* a pickled npz;
* a damaged or incomplete TIFF (a truncated copy, or a file still being written), which tifffile reports while
  reading it, and any TIFF that tifffile cannot read;
* labels that do not fit the image: another number of planes, or another size in y and x;
* in a projected run, a labels file with more than one plane;
* a time series: `germquant run IMAGE --t ...` with a labels recipe exits with code 2 before any timepoint, with
  a `REFUSED:` line on stderr, and process_image and `api.run` refuse `t` in a labels run (decision D15). Run the
  image without `--t`.

## Cellpose `_seg.npy` files

The Cellpose GUI saves its masks in `<name>_seg.npy`, a pickled Python dict. A pickle can run code when it is
loaded, so the run path never reads a `.npy` file. For a file you trust, write its masks to a TIFF once and point
`labels_from` at that file:

```python
import numpy as np, tifffile
masks = np.load("IMG_seg.npy", allow_pickle=True).item()["masks"]
tifffile.imwrite("IMG_masks.tif", masks)
```

To train on those masks instead, the training loop's import command (tranche 5 package P5) reads such a file only
with `--allow-pickle` and writes its masks into a training-pairs folder (image and label slices for training, not
a labels file for `labels_from`):

Windows:
```
.venv\Scripts\germquant.exe import-labels "C:\path\to\IMG_seg.npy" --format cellpose-npy --allow-pickle --image "C:\path\to\IMG.nd2" --pairs "C:\path\to\PAIRS_FOLDER"
```
Mac / Linux:
```
.venv/bin/germquant import-labels /path/to/IMG_seg.npy --format cellpose-npy --allow-pickle --image /path/to/IMG.nd2 --pairs /path/to/PAIRS_FOLDER
```
