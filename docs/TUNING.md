# Spot tuning: choosing a recipe's spot settings from your own points

Spots are tuned, not trained. A spot detector has a few settings that decide what counts as a spot (for
SpotMAX the effect-size threshold, for the LoG detector the signal-to-background ratio), and the right
values depend on the stain, the microscope and the sample. Spot tuning searches a grid of such values on
images you have annotated, scores every combination against your points and writes the best one into the
project's recipe, where every later run uses it.

The design is `docs/DESIGN_system_packs.md`, sections 6.8 and 10 (step 9); the code is `src/germquant/tune/`
and the command is `germquant tune-spots`. The napari plugin's Spot tuning dock saves the points and runs the
same command through the engine (design 8.7).

## In short

1. Run the images of a project (`germquant project run DIR`), so each has a finished run in
   `results/<image_id>/`.
2. Click the spots of a few nuclei (or cells) of one or more images, on those runs, and save them as
   `annotations/<image_id>__points_<annotator>.csv` (the plugin's Spot tuning dock does this).
3. Tune:

   ```
   germquant tune-spots --project DIR
   ```

   It prints the five best settings and writes `tuning/<run_id>/`.
4. Apply the best one:

   ```
   germquant tune-spots --project DIR --apply tuning/<run_id>
   ```

   That writes `recipe.yaml`: the project's config hash changes, the project reads `validated: false`,
   and every image becomes stale until it is run again.

## Truth points

A truth-point file is `annotations/<image_id>__points_<annotator>.csv`. The annotator is 1 to 32 letters,
digits, dashes or underscores, for example your initials. The columns, in this order:

| column | meaning |
|---|---|
| `point_id` | 0, 1, 2, ... (each once) |
| `z_um`, `y_um`, `x_um` | the point in micrometres on the image file's own frame: plane 0 of the file is z 0 and pixel 0 is y and x 0, whatever stride and z range the run you annotated used. `z_um` is 0.0 on a single-plane image |
| `zi`, `yi`, `xi` | the voxel of the point on the file's own grid (stride 1, the full z range) |
| `object_id` | the label of the run's nucleus (or cell) at the point: the label of that voxel, or else the most frequent label right around it (one plane and two pixels), as the RAD-51 benchmark maps expert spots; 0 outside every object |
| `projected` | `true` when the point was placed on a projection, so it is matched in 2D |

The file is ASCII with LF line endings. Write it with `germquant.tune.points.write_points`, which takes the
points as napari places them on a run's layers (z counted from the run's first plane, in the run's
micrometres), the run's labels, its voxel size and its geometry, and fills every column:

```python
from germquant.tune.points import points_file_name, write_points

write_points(project_dir / "annotations" / points_file_name("gonad_01", "rje"), zyx_um,
             labels=run_labels, spacing=(dz, dy, dx), xy_stride=2, z_range=(10, 40))
```

`read_points` reads a file back with its types. A tuning run looks the objects up again in the labels of the
run it tunes on, and says so in its warnings when a point now lies in another object than its file records.

What counts as annotated: an object is scored when at least one of your points lies in it. So click every
spot of each nucleus you annotate; a nucleus you leave out is simply not scored (its detections are neither
false positives nor misses). A nucleus with no spot at all cannot be marked in this release.

One annotator per image: a tuning run scores against one set of points. When two people annotated the same
image, keep each person's files in a folder of their own and pass it with `--points`.

## Running a search

```
germquant tune-spots --project DIR [--points annotations/] [--grid FILE] [--max-runs N] [--out DIR] [--json]
```

- `--points`: the folder of truth-point files, relative to the project (default `annotations/`). Files whose
  names are not truth-point names, files of images the project does not hold, and excluded or time-series
  images are skipped with a warning.
- `--grid`: a grid file (below), relative to the project; without it the default grid of the recipe's spot
  backend is searched.
- `--max-runs`: the most detections the search may make over all images (default 50); a larger grid is
  refused before anything runs.
- `--out`: the run folder, relative to the project (default `tuning/<run_id>/`); it must be absent or empty,
  so a tuning run never overwrites another.
- `--json`: print the result (`run_dir`, `best`, `ranked`, `warnings`) as JSON. Standard output then holds
  only that document; the progress, the warnings and whatever a backend prints (SpotMAX's report) go to the
  error stream.

Exit codes: 0 done; 1 the search ran but no setting could be scored, or a detection failed (one `ERROR:`
line names the error); 2 refused before any detection (the project, the grid, `--max-runs`, an image without a
finished run, a truth-point file that cannot be read, such as a `projected` cell that is neither true nor false
or a repeated `point_id`, two annotators for one image, or no truth point in an object of any image). An image
none of whose points lies in an object of its run only gets a warning: it adds nothing to the scores.

Each annotated image needs a finished run in `results/<image_id>/` with its label image (the default
`output.write_label_images: true`). The search runs only the spot detection again, on that run's labels and at
that run's stride and z range, with the call the spots stage itself makes; every key the grid does not name
keeps the value of the config the run was made with. When the recipe has changed since a run (an earlier
Apply, for example), the search says so in a warning and still starts from that run's own settings: run the
images again first when the change matters.

## Grids

A grid file maps dotted `spots.` keys to the values to try:

```yaml
spots.log.min_snr: [1.5, 2, 2.5, 3, 4, 5]
spots.log.diameter_um: [0.4, 0.5]
```

Every combination is one setting (here 12). The keys a built-in backend reads are checked by name and type:

| backend | keys |
|---|---|
| `spotmax` | `spots.spot_radius_um`, `spots.gauss_sigma_um`, `spots.thresholding_method`, `spots.effect_size_metric`, `spots.effect_size_min`, `spots.merge_z_columns`, `spots.z_merge_gap_um`, `spots.z_merge_valley_frac`, `spots.max_spot_candidates`, and the optional pack keys `spots.spot_radius_z_um`, `spots.remove_hot_pixels`, `spots.sharpen` |
| `log` | `spots.log.diameter_um`, `spots.log.min_snr`, `spots.log.min_separation_um`, `spots.log.annulus_inner_um`, `spots.log.annulus_outer_um`, `spots.log.pre_gate_sigmas`, and the shared `spots.merge_z_columns`, `spots.z_merge_gap_um`, `spots.z_merge_valley_frac`, `spots.max_spot_candidates` |
| any other | any `spots.` key, with plain values |

A number takes the type of the key's code default (2 is stored as 2.0 for a decimal key). `spots.backend`,
`spots.enabled` and `spots.instances` are never tuned: they are recipe choices. Only the main spots stage is
tuned; a recipe's extra spot instances (`spots.instances`) are not.

The default grids (used when you pass no `--grid`; no pack ships a grid of its own yet):

| backend | default grid |
|---|---|
| `spotmax` | `spots.effect_size_min` in 0, 0.5, 1, 1.5, 2, 3, 4, 5 |
| `log` | `spots.log.min_snr` in 1.5, 2, 2.5, 3, 4, 5 |

Any other backend needs a grid file.

How many detections a grid means: SpotMAX runs once per group of settings that differ only in
`spots.effect_size_min` (see below), every other backend once per setting, each on every annotated image.
`germquant.tune.grid.n_runs(grid, n_images, backend)` gives the number, which `--max-runs` bounds.

### SpotMAX: run once, then filter

SpotMAX applies `effect_size_min` as its very last step, a plain row filter after the detection, the
features and the z-split merge. So the search runs SpotMAX once per group with the minimum at 0 and applies
each setting's minimum afterwards (`germquant.tune.effect_filter.filter_effect_size`). The tables are
exactly those of a direct run at that minimum (tests/test_tune_spots.py checks it), and a long list of
minima costs one SpotMAX run per image.

`spots.sharpen` is part of the group like every other key. With sharpen on, SpotMAX computes the effect size
on the sharpened image (design 4.7), so a threshold found with one sharpen value never carries over to the
other: search `spots.effect_size_min` together with `spots.sharpen` when you tune both.

## Scoring and the objective

Every setting is scored on every image inside the annotated objects:

- per object, your count of points against the detected count; Lin's concordance correlation coefficient
  (CCC) over the annotated objects of all images pooled (`object_ccc`);
- per point, a one-to-one matching of the detected spots inside those objects with your points, inside a
  cylinder of 0.3 um in xy and 0.4 um in z (the RAD-51 benchmark's tolerances; points marked `projected` are
  matched in 2D); precision, recall and F1 from the counts summed over the images (`point_f1`).

The objective ranks the settings by the pooled object CCC rounded to 3 decimals, then by the pooled point F1
(higher is better; a missing value ranks last, for example the CCC with fewer than two annotated objects or
the F1 of a setting that detects nothing). Settings with the same two values keep the grid's order. This is
the order the RAD-51 harness ranks by.

When every annotated object holds the same number of points (one focus in each, say), the CCC formula is
0 / 0. A setting whose counts equal yours then scores a CCC of 1.0, and every other setting 0.0, so the
setting that finds exactly your counts still ranks first.

## The run folder

A search writes `tuning/<run_id>/` (the run id is the UTC time and 8 hex digits, as run logs use):

| file | content |
|---|---|
| `grid.yaml` | the grid searched, in the `--grid` format, so `--grid tuning/<run_id>/grid.yaml` repeats the search |
| `trials.jsonl` | one line per setting: `trial`, `params`, `backend`, `split`, `shared_detection` (the SpotMAX run it was filtered from), `objective`, `rank`, `pooled` (object_ccc, point_f1, point_precision, point_recall and the counts of objects, points, matches and detections) and `images` (the same per image, with the time of the detection) |
| `ranked.csv` | `rank`, `objective`, `objective_value`, `tie_break`, `tie_break_value`, then one column per grid key, best first; values to 3 decimals |
| `best.yaml` | `objective` (the text), `backend` and `keys` ({dotted key: value}): what `--apply` writes. Not written when no setting could be scored |
| `run.json` | what the search read: each finished run (config hash, done-marker sha256, stride, z range, voxel size), each point file with its sha256 and how many of its points were scored, the tolerances, the grid, the number of detections and how they were made (`shortcut`: run once, then filter, or one detection per set), the engine version and git sha, and the warnings. The sha256 digests are taken with CRLF turned into LF, so a file re-saved with Windows line endings keeps its digest |

Every file is ASCII with LF line endings, written atomically.

## Applying a result

```
germquant tune-spots --project DIR --apply tuning/<run_id> [--plugin-version V] [--json]
```

writes the run's `best.yaml` into the project's `recipe.yaml` through the project's recipe writer, with the
edit origin `tuning` (`--plugin-version`, which the napari plugin passes, goes into the manifest's saved
stamp; the search options `--points`, `--grid`, `--max-runs` and `--out` are refused here): `recipe.yaml` extends the pack's recipe (`extends: pack:<id>@<version>/<recipe>`) and
adds the best keys. A tuning edit may set keys the pack does not offer in its parameter form (such as
`spots.sharpen`), but only the keys of that run's grid, and only when the run tuned the backend the recipe uses
now. The project's config hash changes, `validated` turns false and every image turns stale with the reason
config_hash. In Python: `germquant.tune.apply_tuning(project, "tuning/<run_id>")`, which returns the new
config hash.

## Python

```python
from germquant import tune

res = tune.tune_spots(project)                      # or the project folder; grid=, points=, out=, max_runs=
print(res["best"], res["ranked"][:5], res["warnings"])
tune.apply_tuning(project, res["run_dir"])
```

The tuners themselves work on any finished runs, outside a project too (the benchmark tune splits use them),
also on the runs of `germquant batch`, whose `run_manifest.json` sits at the top of its output folder:

```python
items = [{"image_id": "img_01", "image": "data/img_01.nd2", "result_dir": "results/img_01",
          "points": "annotations/img_01__points_rje.csv"}]
tune.tune_spotmax(items, {"spots.effect_size_min": [1, 2, 3, 4]}, split="tune", out_dir="tuning/a")
tune.tune_grid(items, {"spots.log.min_snr": [2, 3]}, backend="log", out_dir="tuning/b")
```

`split` is recorded with the run; a split named `test` is refused, because tuning never reads a test split
(design 6.6).

Importing `germquant.tune` loads no torch, CuPy, SpotMAX, Cellpose, bioio or dask, so the plugin writes
points in napari's own process. Every tuner imports `germquant.pipeline` (torch) before it reads an image or
loads a backend, so torch claims cuBLAS before SpotMAX or CuPy can (design rule R5). Of the benchmark code,
tune uses only `germquant.benchmark.stats`, `matching`, `nuclei3d` and `spots` and
`germquant.validate.compare` (the seal of design 6.7, checked by tests/test_tune_seal.py).

## Not in this release

- Single-plane (planar) runs are tuned once the single-plane path (Train B package P2b) is in the engine;
  until then such a run is refused with that reason.
- Time series are not tuned.
- The matching tolerances are fixed (0.3 um in xy, 0.4 um in z); per-protocol tolerances for the public
  benchmarks come with their protocols.
- A nucleus without any spot cannot be marked as annotated.
