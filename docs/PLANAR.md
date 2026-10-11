# The 2D path: objects segmented on one 2D image

Some samples are not segmented well in 3D. Yeast cells and cultured cells are usually segmented on one
image: a single plane of a transmitted-light channel, or a projection of a nuclear stain. Spots inside
those objects should still be counted in 3D. The 2D path does exactly this, as an opt-in recipe choice.
Nothing changes for a recipe that does not ask for it: every C. elegans run, every shipped config and
every existing output stays byte for byte what it was.

There are two modes. This page describes the first; the second is still to come.

| mode | what it is | switched on by | status |
|---|---|---|---|
| projected | a 3D stack whose objects are segmented on one 2D image and copied through every plane; everything after segmentation stays 3D | `segmentation.nuclei.project` | available |
| single plane | an image with one plane, where every quantity is 2D | `io.planar` | not yet (tranche 5 package P2b) |

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

`max` keeps the image's own data type. `mean` is summed in double precision and handed over as float32.

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
- `segmentation.nuclei.project` or `project_plane` holds a value that is not listed above, or the index
  lies outside the stack.

### Limits

- A prism copied through z gives column intensities, and a spot above or below a cell in its column
  counts for that cell. For thick samples, or cells that lie on top of each other, use a 3D recipe.
- On one plane the classical watershed does not fill holes inside an object (a single plane counts as a
  border in z). It is a baseline, as in 3D.
- Single-plane images (Z = 1) still run through today's 3D path with the 1 um voxel guess and the
  `voxel:UNREADABLE_assumed_isotropic_1um` flag, until package P2b adds `io.planar`.

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

## What does not move

- A recipe without these keys takes exactly today's path: the same Cellpose call (do_3D, z_axis 0, no
  channel axis, anisotropy dz / dy, no diameter), the same SpotMAX calls, the same tables, labels,
  montage, flags, manifest keys, stage hashes and completion markers. The tranche 5 fence
  (`tests/test_default_invariants.py`) and the synthetic goldens check this.
- Tests: `tests/test_planar.py` (the selectors, the area gate, the segmentation control flow and the
  refusals), `tests/test_pack_knobs.py` (Cellpose mocked, SpotMAX spied, `require_channel`,
  `qc.min_nuclei`) and `tests/test_projected_run.py`, which runs a synthetic yeast-like stack through a
  projected recipe (`tests/fixtures/planar/projected_yeast_like.yaml`) and compares it with the golden in
  `tests/golden/synthetic_projected/`.
