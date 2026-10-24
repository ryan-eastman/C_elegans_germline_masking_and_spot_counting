# The 2D path: objects segmented on one 2D image

Some samples are not segmented well in 3D. Yeast cells and cultured cells are usually segmented on one
image: a single plane of a transmitted-light channel, or a projection of a nuclear stain. Spots inside
those objects should still be counted in 3D. The 2D path does exactly this, as an opt-in recipe choice.
Nothing changes for a recipe that does not ask for it: every C. elegans run, every shipped config and
every existing output stays byte for byte what it was.

There are two modes, both described here.

| mode | what it is | switched on by | status |
|---|---|---|---|
| projected | a 3D stack whose objects are segmented on one 2D image and copied through every plane; everything after segmentation stays 3D | `segmentation.nuclei.project` | available |
| single plane | an image with one plane, where every quantity is 2D | `io.planar` | available |

The design is `docs/DESIGN_system_packs.md`, section 4. The code is `src/germquant/planar.py`.

## Projected runs

### Switching it on

A projected run is a recipe (a profile) that sets `segmentation.nuclei.project`. Every key below is
read with a code default that reproduces today's behaviour, so the shipped config files never carry
them and their config hashes do not move.

| key | values | code default | meaning |
|---|---|---|---|
| `segmentation.nuclei.project` | `none`, `max`, `mean`, `plane` | `none` | `max` or `mean` segment the maximum or mean projection of the segmentation channel; `plane` segments one plane of it |
| `segmentation.nuclei.project_plane` | `mid`, `sharpest`, or a plane index | `mid` | which plane `plane` uses (see below); ignored by `max` and `mean` |
| `segmentation.nuclei.min_area_um2` | a number | `0.0` (no gate) | the size gate for 2D objects, in square micrometres |
| `segmentation.nuclei.cellpose.object_diameter_um` | a number or null | null | optional: the object diameter Cellpose rescales to (see "Recipe keys for packs") |

The channel that is segmented is still the one `segmentation.nuclei.channel` names (`dna` by default;
a yeast recipe names a transmitted-light role).

An example, the shape of the budding yeast recipe (the numbers are the yeast pack's starting values):

```yaml
# config/profiles/my_yeast_study.yaml (a bare extends name is found next to this file; from outside the
# repository use the absolute path of generic.yaml)
extends: generic.yaml
io:
  channel_map: my_yeast_channels.yaml       # a transmitted role and a foci role
segmentation:
  nuclei:
    method: cellpose
    channel: transmitted
    require_channel: true
    project: plane
    project_plane: sharpest
    min_area_um2: 3.0
spots:
  spot_radius_um: 0.244
  spot_radius_z_um: 1.0
  sharpen: true
  remove_hot_pixels: true
  merge_z_columns: false
```

Run it like any other profile.

Windows:
```
.venv\Scripts\germquant.exe run "C:\path\to\YOUR_STACK.nd2" --profile my_yeast_study --out "C:\path\to\RESULTS_FOLDER"
```
Mac / Linux:
```
.venv/bin/germquant run /path/to/YOUR_STACK.nd2 --profile my_yeast_study --out /path/to/RESULTS_FOLDER
```

### Which plane `plane` uses

- `mid`: the middle plane, index Z // 2 (plane 5 of 11, plane 5 of 10).
- `sharpest`: the plane in best focus. Each plane is smoothed with a Gaussian of sigma 1 pixel, the
  Laplacian of the result is taken (both in float32), and the plane whose Laplacian has the highest
  variance wins. Ties go to the lowest z.
- an index such as `4`: that plane, counting from 0. An index outside the stack stops the run before
  segmentation with a message that gives the number of planes.

Indices count the planes of the stack as read, after `--z-range`, not the planes of the file; `mid` and
`sharpest` also look only at those planes. With `--z-range 10 30` (planes 10 to 29 are read) and
`project_plane: 4`, file plane 14 is segmented. The recorded `segment_projection.z` counts in the same
way, and `z_range` in `run_manifest.json` and `<image_id>__stages.json` gives the offset.

`max` keeps the image's own data type. `mean` is summed in double precision and handed over as float32.
Neither reads `project_plane`, so a value left over from a `plane` recipe is ignored.

### What happens, stage by stage

1. Right after the image is read, before any stage runs, the run is checked (see "Refusals") and the
   plane is resolved.
2. Segment: the segmenter sees one 2D image. `cellpose` runs Cellpose-SAM on the plane (the 2D call,
   `model.eval(image, channel_axis=None)`, with the same model, device choice and retry on the CPU as the
   3D call); `classical` runs the watershed on the plane. Objects smaller than `min_area_um2` are dropped
   (an object with fewer than int(min_area_um2 / (dy * dx)) pixels goes) and the labels are renumbered
   1, 2, 3, and so on. `segmentation.nuclei.min_volume_um3` does not apply to 2D objects. The labels are
   then copied into every plane, so each object is a prism through the whole stack.
3. Measure: intensities are means over each object's column (all planes). `volume_um3` is left empty
   (NA), because a prism is not a cell volume. Two columns are added to the nuclei table, right after its
   declared columns: `area_um2` (the object's pixel count in the 2D labels times dy times dx) and
   `equivalent_diameter_um` (the diameter of a circle with that area). `centroid_z_um` is the middle of
   the stack; it is never empty, so nothing downstream breaks. `n_voxels` counts the column.
4. Spots are counted in 3D, as in every other run. A spot counts for the object whose column it falls in,
   including a spot above or below the cell itself.
5. Write: `<image_id>__nuclei_labels.tif` holds the (Z, Y, X) column labels, so every plane is the same.
   The two new columns appear only in a projected run's nuclei table; no other table changes.
6. Render: the montage draws the outlines on the 2D image the segmenter saw (the chosen plane or the
   projection), not on a maximum projection of the stack, which for a transmitted-light channel shows
   nothing useful.

### What a projected run records

Only a projected run records these; a run without the key records none of them.

- The flag `segment:projected_<mode>` (`segment:projected_plane`, `segment:projected_max`,
  `segment:projected_mean`) in `qc_flags` and in the segment stage's entry of `<image_id>__stages.json`.
- `segment_projection` in `run_manifest.json` and `<image_id>__stages.json`: the mode, the plane
  selector and the plane index, for example `{"mode": "plane", "plane": "sharpest", "z": 3}`. For `max`
  and `mean` both `plane` and `z` are null.
- `segmentation_method` names the segmenter as always (`cellpose`, `classical`, a plugin's name).
- The projection keys sit in the `segmentation` section, so they enter the segment stage's sub-hash and
  `germquant batch --resume` reprocesses an image when they change. The completion marker gains no key.

### Refusals

A projected run stops right after the read, before any pixel is processed and before any file is
written, when something in the recipe cannot work on 2D objects. The error (`PlanarConfigError`, a
`ValueError`) lists every problem with the key that fixes it:

- a stage that only makes sense in 3D is switched on and every channel role it needs is in the image.
  These stages are `germline`, `axis`, `envelope`, `audit`, `staging`, `sc_trace`, `granule`, `coloc`,
  `partition` and `granule_tail` (they are worm biology, and `axis` weights by volume). A stage whose
  channel roles are missing skips itself, as it always has. Switch a stage off in the recipe
  (`axis: {enabled: false}`) or with `--no-<stage>`. `generic.yaml` already switches the worm stages
  off;
- the segmentation backend does not declare that it can segment one plane (`planar_ok`). `cellpose` and
  `classical` can; the experimental `nuclei3d` cannot. Under `method: classical` the backend is never
  called, so it is not checked;
- `segmentation.nuclei.project` holds a value that is not listed above, or, in mode `plane`,
  `project_plane` does or its index lies outside the stack.

A malformed value of one of the recipe keys below stops the run at the same point (see "Recipe keys for
packs").

### Limits

- A prism copied through z gives column intensities, and a spot above or below a cell in its column
  counts for that cell. For thick samples, or cells that lie on top of each other, use a 3D recipe.
- On one plane the classical watershed does not fill holes inside an object (a single plane counts as a
  border in z). It is a baseline, as in 3D.
- A single-plane image (Z = 1) in a recipe without `io.planar` still runs through the 3D path, with the
  1 um voxel guess and the `voxel:UNREADABLE_assumed_isotropic_1um` flag when the file has no z step, exactly
  as before. `io.planar` (next section) makes it a single-plane run.

## Single-plane runs

An image with one plane (a field of cultured cells, a tissue section imaged once) is a 2D image: its objects
have an area, not a volume, and its spots have no depth. A single-plane run treats it so, from the read to the
written images. Every key below is read with a code default that keeps today's behaviour, so the config files
that existed before them never carry them and their config hashes do not move; `generic_2d.yaml` (below) is the
first to set one.

### Switching it on

| key | values | code default | meaning |
|---|---|---|---|
| `io.planar` | `false`, `auto`, `true` | `false` | `auto`: a single-plane run when the image has one plane (below); `true`: a single-plane run, and an image with more than one plane is refused |
| `io.voxel_size_um` | null, a number, `[dy, dx]` or `[dz, dy, dx]` | null | the file's own pixel size in micrometres, before `--xy-stride`, for a file that carries no calibration or a wrong one; it applies to every run kind (below) |

`config/profiles/generic_2d.yaml` is the ready-made recipe. It extends `generic.yaml` (the generic channel map,
the C. elegans stages off; `docs/GENERIC_PROFILE.md`) with `io.planar: auto`, `segmentation.nuclei.min_area_um2:
20.0`, `spots.merge_z_columns: false`, the granule and coloc stages off and a filename regex that matches nothing
(the design columns then read `unknown`; a study profile sets its own regex). Its `io.input_suffixes` lists the
formats a single plane usually comes in; `germquant batch` reads that key from tranche 5 package P5 on, so until
then a batch with this profile finds `.nd2` files only and other formats run one at a time with `germquant run`.

Windows:
```
.venv\Scripts\germquant.exe info "C:\path\to\YOUR_IMAGE.ome.tif"
.venv\Scripts\germquant.exe run "C:\path\to\YOUR_IMAGE.ome.tif" --profile generic_2d --out "C:\path\to\RESULTS_FOLDER"
```
Mac / Linux:
```
.venv/bin/germquant info /path/to/YOUR_IMAGE.ome.tif
.venv/bin/germquant run /path/to/YOUR_IMAGE.ome.tif --profile generic_2d --out /path/to/RESULTS_FOLDER
```

`germquant info` prints a line `planes : single plane, ...` for an image whose header has one plane (no Z axis,
or a Z of 1), with the xy calibration it found in the file, or says that it found none.

### When a run is a single-plane run

- `io.planar: false`, the default: never. A one-plane image runs through the 3D path as it always has (see
  "The legacy golden" below).
- `auto`: when the header has one plane (no Z axis, or a Z of 1) and the stack read has one plane. A 3D file cut
  to one plane with `--z-range` is still a 3D file, so it stays a 3D run. Any other image runs through the 3D
  path, so one recipe serves a folder of planes and stacks.
- `true`: always, and the stack read must have one plane. More than one plane stops the run before any stage
  runs and before any file is written, with a message that names `--z-range` (read one plane) and
  `segmentation.nuclei.project` (segment a projection of the stack instead, with `io.planar` false or auto).

The header test does not use the header's `is_2d`, which means "no Z axis" for an `.nd2` and "Z is 1" for the
other formats.

### The pixel size

dy and dx come from the first of these that exists:

1. `io.voxel_size_um`, multiplied by `--xy-stride` (`spacing_source` "override");
2. the xy calibration a one-plane file carries itself (`spacing_source` "file_yx"). An OME-TIFF without a z step
   reports only PhysicalSizeX and PhysicalSizeY; the 3D reader cannot use that and takes the 1 um guess, a
   single-plane run uses it;
3. the voxel size the reader found, when it found one (an `.nd2` plane; `spacing_source` "file");
4. 1 um per pixel, times the stride, with the flag `voxel:UNREADABLE_xy_assumed_1um` (`spacing_source`
   "assumed").

A TIFF whose resolution unit is "none" (what tifffile writes without a calibration) or an ImageJ TIFF in pixels
counts as uncalibrated. A TIFF from an image editor, at 72 dots per inch, carries a screen resolution rather than
a pixel size: set `io.voxel_size_um` for such files, and check the `germquant info` line first.

The file's z step plays no part in a single-plane run: its tables are identical for a file dz of 0.2 and of 5.0
(tested, gate G10). The tables write `voxel_dz_um` as NaN. Array code receives a z spacing of 1.0
(`planar.PLANAR_DZ`), a placeholder that reaches only code whose result does not depend on dz for one plane.

`io.voxel_size_um` applies to 3D and projected runs as well (decision D6 of the Train B plan). A number or
`[dy, dx]` sets the xy pixel size and keeps the file's dz; `[dz, dy, dx]` sets all three (a dz has no use in a
single-plane run). Every run that sets the key raises `voxel:override`, plus `voxel:override_differs_from_file`
when the file's own calibration of an axis the key sets is more than 1 percent away, and records
`spacing_source` "override". When a 3D run's dz is still the 1 um guess (an override without dz, on a file
without a voxel size), the run also keeps today's flag `voxel:UNREADABLE_assumed_isotropic_1um`.

### What happens, stage by stage

1. Read: the decision and the pixel size above, and the refusals below, before any stage runs.
2. Segment: the plane is segmented by the 2D segmentation that projected runs use (`segment_planar`: Cellpose's
   2D call, or the classical watershed), objects below `segmentation.nuclei.min_area_um2` are dropped and the
   labels are renumbered 1, 2, 3, and so on. `segmentation.nuclei.min_volume_um3` does not apply. A
   `segmentation.nuclei.project` setting is ignored, because one plane needs no projection, and the run raises
   `segment:project_ignored_planar` (decision D7; an unknown mode is still refused).
3. Measure: as in 3D, then `volume_um3` is left empty (NA), `centroid_z_um` is 0.0, and `area_um2` (the pixel
   count times dy times dx) and `equivalent_diameter_um` (the diameter of the circle of that area) follow the
   declared nuclei columns. `n_voxels` is the pixel count.
4. Spots: SpotMAX runs on the plane (`docs/BACKENDS.md`): radii (`PLANAR_RZ_PX`, r / dy, r / dx), the Gaussian
   sigma over the smaller of dy and dx, and no z-split merge. Every `z_um` is 0.0. When the recipe asks for the
   merge (`spots.merge_z_columns: true`, as `config.yaml` sets it), the flag `spots:z_merge_off_planar` (or
   `spots_<name>:z_merge_off_planar` for a spot instance) says that it did not run. `spots.remove_hot_pixels` and
   `spots.sharpen` apply as in 3D; `spots.spot_radius_z_um` has no meaning on one plane and is ignored, with a
   warning in the log. Only `spotmax` can detect spots on one plane.
5. Write: `<image_id>__nuclei_labels.tif` is an int32 OME-TIFF with axes YX and PhysicalSizeX and Y only (no Z),
   `<image_id>__spots.tif` a YX image with a disk at each spot, and masks are YX too. The two 2D columns appear
   in the nuclei table only.
6. Render: the montage as for any run, with the run's own dy for the scale bar.

### What a single-plane run records

Only a single-plane run records these, except `voxel:override`, `voxel:override_differs_from_file` and
`spacing_source` "override", which any run that sets `io.voxel_size_um` records (above); a run without either key
records none of them.

- The flag `read:planar`, and `voxel:UNREADABLE_xy_assumed_1um` or `voxel:override` (with
  `voxel:override_differs_from_file`) as above.
- `voxel_dz_um` NaN in every table.
- In `run_manifest.json` and `<image_id>__stages.json`: `planar: true` and `spacing_source`.
- `"planar": true` in the read, segment, measure and spots stage sub-hashes, and the done-marker extra
  `{"planar": true}`. `germquant batch --resume` computes the planar segment hash once per batch and compares a
  single-plane marker with it, so a finished image is skipped and an edited recipe reprocesses it. A caller that
  gives only the 3D segment hash (the time-lapse loop, an outside caller of `is_done`) reprocesses a single-plane
  image, with the reason "planar flag", rather than skip it unchecked (decision D8).

### Refusals

A single-plane run stops right after the read, before any pixel is processed and before any file is written:

- a malformed `io.planar` or `io.voxel_size_um` (a quoted `"true"` is text, not a boolean, and a list item
  cannot be null) is refused first, by the recipe check: a `ValueError` naming the key, as for the recipe keys
  below;
- `io.planar: true` on an image with more than one plane raises a `PlanarConfigError` of its own (above);
- so does a `segmentation.nuclei.project` value that a projected run would refuse (an unknown mode, or in mode
  `plane` an unknown `project_plane` selector), although a single-plane run then ignores the setting (decision
  D7);
- the checks of stages and backends are collected into one `PlanarConfigError` that lists every problem: those
  a projected run makes (a 3D-only stage switched on with its channel roles bound, a segmentation backend without
  `planar_ok` such as `nuclei3d`), and a spot backend without `planar_ok`, for the main spots instance and for
  every `spots.instances` entry whose role is bound (`log` and `net3d` cannot run on one plane).

`PlanarConfigError` is a `ValueError`, so a caller that catches `ValueError` sees every refusal above.

### Time series

Single-plane time series are refused in this release (decision D9; design section 4.13 allows a refusal). With
`io.planar: true`, or `auto` on an image whose header has one plane, `germquant run IMAGE --t ...` exits with code
2 before any timepoint, `api.run(..., t=...)` raises `PlanarConfigError` before any file is written, and
`process_series` records every timepoint as failed, each with this message:

> single-plane time series are refused in this release: io.planar is true, or auto on an image with one plane,
> and timepoints were asked for (--t, or t=). Run the series with io.planar false, which processes each
> timepoint as a stack of one plane, as before; or save the timepoints as separate single-plane images and run
> those (docs/PLANAR.md, 'Time series')

### The z radius SpotMAX gets on one plane

SpotMAX's spot detection takes a radius in pixels for each axis. On one plane the z radius has no physical
meaning, so a single-plane run gives SpotMAX 1.0 (`spots.detect.PLANAR_RZ_PX`), the value SpotMAX 1.3.1's own
core sets when an image has one plane (design 4.7). `tests/test_planar_spots.py::test_planar_rz_px_table` runs
one synthetic plane (24 planted puncta in 4 nuclei, dy 0.1 um, dx 0.125 um, the spots stage's default SpotMAX
parameters) with 0.5, 1.0 and 2.0:

| PLANAR_RZ_PX | spots | matched to planted puncta | F1 | tables equal to those of 1.0 |
|---|---|---|---|---|
| 0.5 | 24 | 24 | 1.000 | yes |
| 1.0 | 24 | 24 | 1.000 | (the shipped value) |
| 2.0 | 24 | 24 | 1.000 | yes |

The three tables are identical, as expected from SpotMAX's peak finder, which squeezes a (1, Y, X) image and keeps
only the yx minimum distance. 1.0 is the shipped value.

### The legacy golden

`tests/golden/single_plane_legacy/` records what the tranche 5 base commit's code (8eb73de) writes for the
calibrated plane of `tests/planar_helpers.py`, written as an OME-TIFF CYX (no z step) and as an OME-TIFF CZYX
with Z = 1 (with a z step), under `generic_2d.yaml` with `io.planar` false: every table, for each form. It was
made on that commit's code before any single-plane code existed, so `tests/test_planar_bioio.py` compares the new
code with the old one, not with itself: with `io.planar` false the same files give the same tables cell for
cell. The CYX file keeps the 1 um guess and its flag; the CZYX file, which has a z step, runs as a calibrated
stack of one plane. To make it again (only for a deliberate change, with a detached worktree of the base commit
at BASE): from this checkout, `GERMQUANT_UPDATE_GOLDEN=1 PYTHONPATH=BASE/src python -m pytest
tests/test_planar_bioio.py -k legacy`. The update refuses to run on code that has `planar.decide`.

### Limits

- A single-plane run has no time series (above) and no 3D-only stage: germline, axis, envelope, audit,
  staging, sc_trace, granule, coloc, partition and granule_tail need the third dimension or worm biology.
- The classical watershed smooths with a Gaussian of sigma 0.4 um before its threshold, which moves an object's
  edge outward on fine pixels. At the 0.8 um pixel of the synthetic golden it recovers every planted disk
  exactly; a check on the same kind of flat disks (radius 4 to 6 um) gave areas 7 to 22 percent too large at
  pixels of 0.325 to 0.65 um. It is a baseline, as in 3D; Cellpose's 2D call does not share this rule.

## Segmentation backends and 2D objects

A segmentation backend (`docs/BACKENDS.md`) that can segment one plane says so with an attribute,
`planar_ok = True`. In a projected run it then receives the plane as a (1, Y, X) array with the extra
keyword `planar=True`, and returns (1, Y, X) labels. A backend without the attribute never sees that
keyword, and a projected run that selects it is refused. The built-in `cellpose` and `classical`
declare it; `nuclei3d` does not.

## Recipe keys for packs

The same package added four recipe keys that the outside packs need. Each is read with a code default
that keeps today's behaviour, and none of the shipped configs sets one.

- `segmentation.nuclei.require_channel` (default false). When true and the image has no channel for the
  segmentation role, the image fails before segmentation with a message that names the role. Without it
  the segment stage segments channel 0 and raises `segment:no_<role>_channel_using_first_channel`, which
  can silently take a fluorescence channel for the cells. Every outside pack sets it.
- `qc.min_nuclei` (default 10). The threshold of the advisory flag `qc:few_nuclei_<n>`, raised when an
  image has fewer objects. A field of cultured cells holds a handful of nuclei, so the cultured-cell pack
  sets 1; the validated C. elegans path keeps 10. It is passed on only when a recipe sets it.
- `segmentation.nuclei.cellpose.object_diameter_um` (default null). When set, Cellpose is given
  diameter = d / dy pixels (dy after `--xy-stride`) in both the 3D and the 2D call, and Cellpose 4
  rescales the image by 30 / diameter. It is for objects outside Cellpose-SAM's native size range, such
  as nuclei imaged at 100x that span 150 to 300 pixels. Unset, both calls are exactly today's. The
  older `segmentation.nuclei.diameter_um` keeps its meaning: Cellpose-SAM ignores it, the classical
  watershed uses it to space its seeds.
- Three optional SpotMAX keys in the `spots:` section, forwarded only when a recipe sets them:
  `spot_radius_z_um` (the z radius: SpotMAX's radii become (rz, r, r) / spacing instead of r / spacing),
  `remove_hot_pixels` (SpotMAX's hot-pixel removal before the smoothing) and `sharpen` (SpotMAX's
  difference-of-Gaussians sharpening). With `sharpen`, detection runs on the sharpened image and SpotMAX
  computes the spot-versus-background effect size on it, as in SpotMAX's own runs, so `effect_size` and
  `effect_size_min` are then on the sharpened scale: a threshold tuned without `sharpen` (the RAD-51
  value 3.0 among them) does not carry over. The other intensity features stay on the unsharpened image.
  A `spots.instances` entry can set these three keys for its own instance.

These keys are checked strictly in every run that sets one, right after the image is read and before any
stage runs or any file is written (`germquant.recipe_keys`). A malformed value stops the image with a
`ValueError` that names the key:

| key | accepted values |
|---|---|
| `segmentation.nuclei.require_channel`, `spots.remove_hot_pixels`, `spots.sharpen` | `true`, `false` or null |
| `qc.min_nuclei` | a whole number, 0 or more, or null |
| `segmentation.nuclei.cellpose.object_diameter_um`, `spots.spot_radius_z_um` | a positive number of micrometres, or null |
| `segmentation.nuclei.cellpose` | a mapping (holding `object_diameter_um`), or null |

A quoted value is text, so `require_channel: "false"` is refused rather than read as true. The check also
runs under `method: classical` and when spots are off, so a typo shows up in a CPU test run of the recipe
rather than only on the GPU. It matters most for `object_diameter_um`: under `method: auto` an error
inside the Cellpose call makes the segment stage fall back to the classical watershed, so without this
check a typo would silently swap Cellpose masks for watershed masks.

## What does not move

- A recipe without these keys takes exactly today's path: the same Cellpose call (do_3D, z_axis 0, no
  channel axis, anisotropy dz / dy, no diameter), the same SpotMAX calls, the same tables, labels,
  montage, flags, manifest keys, stage hashes and completion markers. The tranche 5 fence
  (`tests/test_default_invariants.py`) and the synthetic goldens check this.
- Tests: `tests/test_planar.py` (the selectors, the area gate, the segmentation control flow and the
  refusals), `tests/test_pack_knobs.py` (Cellpose mocked, SpotMAX spied, `require_channel`,
  `qc.min_nuclei`), `tests/test_recipe_keys.py` (the strict checks of the recipe keys, and that no
  shipped config sets one) and `tests/test_projected_run.py`, which runs a synthetic yeast-like stack
  through a projected recipe (`tests/fixtures/planar/projected_yeast_like.yaml`) and compares it with the
  golden in `tests/golden/synthetic_projected/`.
- A recipe without `io.planar` and `io.voxel_size_um` reads, segments, measures, counts and writes a
  one-plane image exactly as before: `tests/test_planar_bioio.py` compares that run with the legacy golden
  made on the base commit's code, and the fence pins the legacy single-plane rule.
- Single-plane tests: `tests/test_planar_single.py` (the decision, the pixel size, the override, the recipe
  keys, resume's safe default, `germquant info`, the readers), `tests/test_planar_run.py` (the synthetic_2d
  golden in `tests/golden/synthetic_2d/`, whose sentinel is generic_2d.yaml's own hash, dz invariance, the YX
  writers, the refusals, the stage hashes and `batch --resume`), `tests/test_planar_bioio.py` (real one-plane
  OME-TIFFs through bioio, and the legacy golden), `tests/test_planar_series.py` (time series refused) and
  `tests/test_planar_spots.py` (SpotMAX on one plane, the `PLANAR_RZ_PX` table). The test images come from
  `tests/planar_helpers.py`.
