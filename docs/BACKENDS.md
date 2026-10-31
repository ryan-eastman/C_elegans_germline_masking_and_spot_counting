# Pluggable backends

The pipeline's model-driven steps (nucleus segmentation today, spot detection next) run behind a
small registry, so that another segmenter can be plugged in without forking the pipeline and a
head-to-head against the current model is a config change, not a rewrite. This document is the
contract a backend has to meet.

A measurement step that runs after the pipeline and writes tables of its own, without replacing a step,
is not a backend but a post-stage: its contract is docs/POST_STAGES.md.

## Nucleus segmentation backends

### Picking one

The config key `segmentation.nuclei.backend` names the segmenter. Its code default is `cellpose`,
and the key is read with that default, so none of the shipped config files changes and the config
hash of every existing run keeps reproducing. Put the key in a profile rather than in
`config/config.yaml`:

```yaml
# config/profiles/my_segmenter.yaml
extends: config.yaml

segmentation:
  nuclei:
    method: auto             # unchanged: auto = the backend, classical watershed on failure
    backend: threshold_cc    # the backend's name (its entry point name, see below)
    threshold_cc:            # keys of the backend's own; the whole section reaches it as cfg_section
      sigma_um: 0.3
```

then `germquant run --profile my_segmenter <file.nd2>` (or `--config config/profiles/my_segmenter.yaml`).
`method` keeps its three values and its meaning; `backend` only says which segmenter `auto` and
`cellpose` run:

| `method` | what runs | on failure |
|---|---|---|
| `auto` (default) | the backend | logged, then the classical watershed; `segmentation_method` says `classical` |
| `cellpose` | the backend, whatever its name | the error propagates and the segment stage fails (it is fatal) |
| `classical` | the classical watershed; the backend is never called | the error propagates |

An unknown backend name raises `UnknownBackendError` under `auto` and `cellpose` alike, with the list
of available names in the message: a misspelt name is a config error, never a silent switch to the
fallback.

### Which channel is segmented

`segmentation.nuclei.channel` names the channel ROLE the segment stage hands to the backend. Every
shipped config says `channel: dna` and the key is read with that code default, so the validated runs and
their goldens segment the DAPI channel exactly as before. A profile for a stain-free live series names
another role (`config/profiles/live_syp2.yaml` says `central_element`, the GFP::SYP-2 channel;
`docs/TIMELAPSE.md`, "Segmenting without DAPI"). When the named role does not resolve in an image the
stage falls back to channel 0 and flags `segment:no_<role>_channel_using_first_channel`. A recipe that
sets `segmentation.nuclei.require_channel: true` (code default false; every outside pack sets it) fails
the image instead, before segmentation, with a message that names the role. The montage's
outline panel sits on the same channel. A backend therefore receives "the channel to segment", which is
the DNA channel unless the profile says otherwise; its first argument keeps the name `dna` for
compatibility.

`segmentation.nuclei.cellpose_model` takes a file path, a Cellpose built-in name (`cpsam`) or the name of
a registered model: a release asset `get-model` has downloaded, or a model added with `germquant
register-model` (`models_registry.register_local_model`), which resolves to the registered file and puts
its sha256 in the run manifest as always.

A file path must name a file. An absolute path is used as it stands. A relative one is looked for under
`$GERMQUANT_MODELS` (with and without its leading `models/`), then next to the config (its base folder),
then from the working directory, where Cellpose itself looks. A relative path that a profile inherits from
a base config in another folder was already made absolute against that folder when the config loaded, so
`$GERMQUANT_MODELS` is not searched for it. A path found nowhere stops the run before segmentation
(`pipeline.resolve_model_path`, `germquant/model_paths.py`): a FileNotFoundError lists the places looked in
and the fixes (fetch the release model with `germquant get-model <name>`, set `GERMQUANT_MODELS` to a folder
that holds `models/<name>`, or give the file's absolute path). Cellpose itself would not stop: it would log
"pretrained model ... not found, using default model" and segment with its stock weights, and the run would
record `model_sha256: absent`. A run that loads no Cellpose model is not checked: method `classical`, another
backend, method `auto` where Cellpose is not installed (the classical watershed runs, as before), or a recipe
that reads its nuclei from labels files (`segmentation.nuclei.labels_from`, docs/LABELS.md: no model is
resolved at all, decision D12). Every
path that resolves behaves as before. To segment with the stock weights on purpose, name them:
`cellpose_model: cpsam`. The check covers config-driven runs (`germquant run`, `run --t`, `batch`,
`api.run`); a direct call of a building block (`api.segment_nuclei`, a backend called with `model=`) hands
its model value to Cellpose unchecked. `germquant finetune --pretrained` takes a model file or a name
Cellpose's command line knows (`cpsam`, `cpsam_v2`, `cpdino`, `cpdino-vitb`) and refuses anything else,
for the same reason.

### What the output says

`segmentation_method` in `<image_id>__image_summary.csv` (the only table that carries it) names the backend that produced the labels
(`cellpose`, `classical`, or the plugin's name), or `classical` after a fallback; it still says
`cellpose` for every validated run. The segment stage's sub-hash in `run_manifest.json` covers the
whole `segmentation` section, so a run with a different backend or different backend keys gets a
different sub-hash. The manifest's `model_sha256` is still the digest of
`segmentation.nuclei.cellpose_model`: a backend that keeps its model under that key gets the
provenance for free. A backend with its own model key declares it through an optional attribute of its
callable, `provenance(cfg_section) -> dict` (added with `nuclei3d`): when the record is not empty the
pipeline adds it to the segment stage sub-hash and to `run_manifest.json` (`segment_backend`), and
`germquant batch --resume` and the time-lapse series compute their segment sub-hash the same way. The
built-in `cellpose` and `classical` declare nothing, so a default run's hashes and manifest are
unchanged; `germquant.segment.backend_provenance(name, cfg_section)` reads the attribute and never
raises.

### The contract

A backend is a plain callable:

```python
def segment(dna, spacing, *, model, diameter_um, device, cfg_section):
    ...
    return labels
```

| argument | type | meaning |
|---|---|---|
| `dna` | `np.ndarray`, shape (Z, Y, X) | the channel to segment, as read from the .nd2 (uint16; float32 for a synthetic stack), not normalised: the DNA (DAPI) channel by default, the role `segmentation.nuclei.channel` names otherwise (GFP::SYP-2 in the live profile) |
| `spacing` | `tuple[float, float, float]` | (dz, dy, dx) voxel size in micrometres from the .nd2 metadata; every physical scale comes from it, never from a constant |
| `model` | `str` | the resolved `segmentation.nuclei.cellpose_model`: a model file path (resolved against the repository or `$GERMQUANT_MODELS`) or a bare name such as `cpsam`. Use it as your model reference or ignore it |
| `diameter_um` | `float` | `segmentation.nuclei.diameter_um`, the expected nucleus diameter |
| `device` | `str` or `None` | `segmentation.nuclei.device`: `cuda`, `mps`, `cpu`, `auto` or `None` (auto). Resolve it with `germquant.device.select_device(device)`, which honours `$GERMQUANT_DEVICE` and falls back when an accelerator is missing |
| `cfg_section` | `dict` | a copy of the whole `segmentation.nuclei` section; read your own keys from it, by convention under a sub-key named after the backend |

Return an integer label image: an `np.ndarray` of shape (Z, Y, X), any integer dtype, no negative
values, 0 for background and one positive label per nucleus. Anything else (another shape, a float
probability map, a bool mask) is rejected with an error that names the backend. `segment_nuclei`
then drops objects below `segmentation.nuclei.min_volume_um3` and relabels consecutively, exactly as
it does for the Cellpose path, so a backend does not filter by size itself. A backend writes nothing
to disk; the label image goes to `<image_id>__nuclei_labels.tif` as always.

The callable is invoked with `dna` and `spacing` positionally and everything else by keyword, so a
backend that declares exactly the signature above works; accepting `**kwargs` as well keeps it
working if a keyword is added later.

Raise on failure rather than return an empty image: under `method: auto` the pipeline then falls
back to the classical watershed and says so in the log and in `segmentation_method`, whereas an
empty label image would run the rest of the pipeline on nothing.

One optional attribute extends the contract for 2D objects (`docs/PLANAR.md`): a backend that can
segment one plane sets `planar_ok = True` on its callable. In a projected run
(`segmentation.nuclei.project`) it then receives the plane as a (1, Y, X) array with the extra keyword
`planar=True` and returns (1, Y, X) labels; the area gate (`segmentation.nuclei.min_area_um2`) and the
relabel follow, and the labels are copied through the stack. A backend without the attribute is never
given that keyword, and a projected run that selects it is refused before any pixel is processed
(`germquant.planar.check_planar`). `cellpose` and `classical` declare it, `nuclei3d` does not.

### The built-ins

| name | implementation |
|---|---|
| `cellpose` | `germquant.segment.backends.cellpose_backend`: `nuclei._cellpose` (Cellpose-SAM in 3D, anisotropy from the spacing, the fine-tuned germline model) with the explicit device from `germquant.device` and one retry on the CPU when the accelerator fails (out of memory, an operator Metal cannot run). In a projected run it makes the 2D call, `nuclei._cellpose_2d` (`model.eval(plane, channel_axis=None)`, the same model and retry). Its one key of its own, `segmentation.nuclei.cellpose.object_diameter_um` (default null), passes Cellpose diameter = d / dy pixels in both calls; unset, both calls are exactly the validated ones. A value that is not a positive number stops the image before any stage, under every method (`germquant.recipe_keys`), because under `method: auto` an error inside the call would only switch the run to the classical fallback |
| `classical` | `germquant.segment.backends.classical_backend`: `nuclei._classical`, Gaussian smoothing, Li threshold and a distance-transform watershed, all spacing-aware. No model and no GPU; a baseline, not a substitute for the fine-tuned model on crowded pachytene nuclei. It runs unchanged on the one plane of a projected run |
| `nuclei3d` | EXPERIMENTAL, never the default: `germquant.segment.backends.nuclei3d_backend`, the in-house 3D U-Net of `germquant.segment.nuclei3d` (below) |

### Built-in `nuclei3d` (experimental in-house 3D model)

Status: EXPERIMENTAL. Selected only by `segmentation.nuclei.backend: nuclei3d` with a model file in
`segmentation.nuclei.nuclei3d.model`; never the default, and not validated beyond the development
scores on the train and val gonads in `docs/BENCHMARK.md` ("nuclei3d"). It needs torch. The model is
trained with `germquant train-nuclei3d`.

The chain:

1. Input: the channel the segment stage hands every backend (DAPI unless the profile says otherwise),
   normalised per image to (raw - p1) / (p99.9 - p1) with the percentiles of the whole volume read on every
   fourth voxel in y and x.
2. A compact residual 3D U-Net (four levels, 16 to 128 channels, an xy-only first pooling for the 1.85
   anisotropy) gives three maps: the foreground probability, the distance to the boundary of the voxel's own
   nucleus in um (the cell distance of Scherr et al. 2020) and the closeness of the nearest neighbouring
   nucleus (1 at a contact, 0 beyond 1 um). Inference is tiled in float32 with an overlap equal to the
   network's receptive radius, so the result does not depend on the tile size beyond float rounding; the
   whole field is inferred, never gated by another segmentation.
3. Instances: the foreground (probability above `fg_threshold`), a seed map S = cell distance x (1 -
   `neighbour_weight` x neighbour closeness), seeds from the h-maxima of S with `h_um`, and a seeded
   watershed of -S inside the foreground, per connected component of the foreground (exact, and the memory
   stays at the size of the largest component).
4. `segment_nuclei` then applies `segmentation.nuclei.min_volume_um3` and the consecutive relabel, as for
   every backend.

Keys under `segmentation.nuclei.nuclei3d`, all read with code defaults:

| key | code default | meaning |
|---|---|---|
| `model` | none (required) | the model file `germquant train-nuclei3d` wrote (`MODELDIR/model.pt`); an absolute path, or one relative to the working directory |
| `device` | the stage's `segmentation.nuclei.device` | `auto` (CUDA when available, honours `GERMQUANT_DEVICE`), `cpu` or `cuda`; Apple MPS runs on the CPU (untested) |
| `tile_zyx` | [96, 384, 384] | inference tile in voxels |
| `fg_threshold` | the model's own, else 0.5 | foreground probability threshold |
| `h_um` | the model's own, else 0.4 | h of the h-maxima seeds, in um of the distance map |
| `neighbour_weight` | the model's own, else 1.0 | how far the neighbour closeness presses the distance map down at contacts (0: the distance map alone) |

`segmentation.nuclei.cellpose_model` and `diameter_um` are not used by this backend (the pipeline still
resolves and hashes the Cellpose model key as for any backend: a missing file under a path records the
digest `absent`, while a bare name must be a registered model or a Cellpose built-in such as `cpsam`).
Under `method: auto` a failure (no model key, a missing file) falls back to the classical
watershed and `segmentation_method` says `classical`; `method: cellpose` makes it an error. A run records
`segmentation_method` `nuclei3d` and the model file's sha256 in the segment stage sub-hash and in
`run_manifest.json` under `segment_backend`. The masks follow the convention of the Imaris surfaces the
model was trained on, which are tighter than the deployed Cellpose masks (about 1.4 times the Imaris volume),
so nucleus volumes, perinuclear shells and everything measured in them change with the backend: compare a
study's two runs stage by stage before using the numbers, and name the backend in the methods text.

```yaml
extends: config.yaml
segmentation:
  nuclei:
    backend: nuclei3d
    nuclei3d:
      model: C:/Users/ryane/germquant_golden/models/nuclei3d_20261010/model.pt
```

### Registering a backend in the same process

```python
from germquant.segment import register_backend, segment_nuclei

register_backend("threshold_cc", segment)
labels, used = segment_nuclei(dna, spacing, method="auto", backend="threshold_cc",
                              cfg_section={"threshold_cc": {"sigma_um": 0.3}})
assert used == "threshold_cc"
```

`register_backend(name, fn, replace=False)` refuses a name that is already taken unless
`replace=True`; `unregister_backend(name)` undoes it (a built-in name goes back to the built-in).
`list_backends()` and `backend_sources()` show what is available and where each backend came from.

### A third-party backend package

An installed package declares its backend as an entry point in the group
`germquant.segment_backends`; the entry point's name is the backend name. Nothing in germquant has
to change: the package is discovered the first time a backend is looked up. A minimal, complete
package follows. The example segments by a threshold plus connected components, enough to show the
plumbing; a real backend would call StarDist, a 3D U-Net or a differently trained Cellpose in the
same place.

```
germquant-threshold-cc/
  pyproject.toml
  src/germquant_threshold_cc/__init__.py
```

`pyproject.toml`:

```toml
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "germquant-threshold-cc"
version = "0.1.0"
requires-python = ">=3.11"
dependencies = ["germquant", "numpy", "scipy", "scikit-image"]

[project.entry-points."germquant.segment_backends"]
threshold_cc = "germquant_threshold_cc:segment"

[tool.setuptools.packages.find]
where = ["src"]
```

`src/germquant_threshold_cc/__init__.py`:

```python
"""A nucleus segmentation backend for germquant: Gaussian smoothing, Otsu threshold, connected
components. Registered as 'threshold_cc' through the germquant.segment_backends entry point."""
import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu


def segment(dna, spacing, *, model, diameter_um, device, cfg_section):
    opts = cfg_section.get("threshold_cc", {})          # this backend's own keys
    sigma_um = float(opts.get("sigma_um", 0.3))
    smooth = ndi.gaussian_filter(dna.astype(np.float32), sigma=[sigma_um / s for s in spacing])
    labels, _n = ndi.label(smooth > threshold_otsu(smooth))
    return labels.astype(np.int32)                      # (Z, Y, X), 0 = background
```

Install it into the same environment as germquant (`pip install -e .` from the package directory)
and check that it is seen:

```
python -c "from germquant.segment import backend_sources; print(backend_sources())"
# {'cellpose': 'built-in', 'classical': 'built-in',
#  'threshold_cc': 'entry point germquant_threshold_cc:segment'}
```

Rules of discovery: a backend whose entry point fails to load (an import error, a missing attribute,
a value that is not callable) is logged as a warning and skipped, never raised, so a broken plugin
cannot stop a run that does not use it. A plugin cannot shadow a built-in or a backend registered
earlier in the process; its entry point is skipped with a warning. Discovery runs once per process;
`discover_entry_points(force=True)` scans again.

### Before trusting a new backend

Compare it with the built-in on the same images before using its numbers: `segmentation_method`
keeps the two runs apart in every table, and the synthetic goldens (`tests/test_golden.py`) certify
only the built-in paths. A different backend is a number-changing choice, so its name belongs in the
profile and in the methods text of anything published from the run.

## Spot detection backends

### The contract

A spot backend is a callable registered under a name in `germquant.spots.backends`:

```python
backend(image, labels, spacing, *, cfg_section, marker) -> (per_spot, per_nucleus)
```

| argument | meaning |
|---|---|
| `image` | the foci channel as a `(Z, Y, X)` numpy array (uint16 or float, as read from the stack) |
| `labels` | the nucleus label image, same shape, 0 = background |
| `spacing` | `(dz, dy, dx)` voxel size in microns; every length a backend uses is converted with it |
| `cfg_section` | the raw `spots:` mapping of the config with the resolved SpotMAX parameters overlaid; a backend reads its own keys from it with code defaults and ignores the rest |
| `marker` | the marker name from the channel map (`RAD-51`, `COSA-1`), written into both tables |

It returns two DataFrames with exactly these columns, in this order:

* per spot: `spot_id, nucleus_id, marker, z_um, y_um, x_um, effect_size, spot_mean_intensity, detector`
  (positions in microns on the image grid, `spot_id` counting from 0, one row per detected spot inside
  a nucleus label);
* per nucleus: `nucleus_id, marker, n_spots, mean_spot_intensity, detector`, one row for EVERY label in
  `labels`, so a nucleus with no spot is a real zero and not a missing row.

`detector` carries the registry name of the backend in every row. `germquant.spots.detect_spots`
dispatches on the name, checks that both tables have the declared columns (missing columns raise, extra
columns are kept after the declared ones, as the append-only schema rule allows) and warns when the
detector column does not carry the name. The constants `PER_SPOT_COLS` and `PER_NUC_COLS` in
`germquant.spots.backends` are the contract in code.

Two optional attributes of the backend callable extend the contract (added with `net3d`; a backend
without them is called exactly as before, so `spotmax` and `log` never see a difference):

* `input_roles`, a tuple of channel roles beyond the foci channel (net3d: `("dna",)`). The pipeline's
  spots stage, its spot instances and the detector bench harness pass the roles the image has as
  `detect_spots(..., channels={role: array})`, and the backend receives them as the keyword `channels`.
  A backend without the attribute is never given that keyword, whatever the caller passes.
  `germquant.spots.backends.backend_input_roles(name)` reads it.
* `provenance(cfg_section) -> dict`, what a run records before any pixel is read (net3d: the model file
  and its sha256). When the dict is not empty the pipeline adds it, with the channel index of every
  extra role, to the spots stage sub-hash and to `run_manifest.json` (`spot_backend`); on the default run
  it is empty and both are unchanged. `backend_provenance(name, cfg_section)` reads it and never raises.

A third optional attribute, `planar_ok = True`, is for single-plane images (`docs/PLANAR.md`): in a single-plane run the backend gets the (1, Y, X) plane and the keyword `planar=True`. `spotmax` declares it; `log` and `net3d` do not, and a single-plane run that selects them is refused before any pixel is processed.

### Selecting a backend

`spots.backend` names the backend for the RAD-51 instance (code default `spotmax`); an entry of
`spots.instances` may carry its own `backend` key, resolved like every other instance key (the
instance value, then `spots.backend`, then the code default). None of the config files is edited for
this; a profile selects it:

```yaml
extends: config.yaml
spots:
  backend: log
  log:
    diameter_um: 0.5
    min_snr: 2.0
```

An instance may also carry its own `log:` mapping; it merges one level deep over the shared one, so a
COSA-1 instance can set a larger `diameter_um` and keep the shared `min_snr`.

What a run records: the QC flag `spots:<backend>_n=<count>` (so the default run still writes
`spots:spotmax_n=...`), the extra flag `spots:unvalidated_backend=<name>` for every backend other than
`spotmax`, the `detector` column of the spots table, and the spots stage sub-hash, which covers the
`spots:` section text and therefore changes when the backend or its keys change. An unknown backend
name raises `ValueError` naming the known backends; inside the pipeline that is recorded as
`spots:FAILED_ValueError` and the run continues with an empty spots table, like any other stage failure.

### Built-in `spotmax` (the validated default)

The SpotMAX chain exactly as the stage always ran it: preprocess, semantic segmentation of where spots
can be, peak detection inside the nucleus masks, per-spot features with the spot-versus-background
effect size, the flood guard, the z-split merge and the effect-size filter. Its code is
`germquant.spots.detect.detect_spots_spotmax`; the registered backend is a thin wrapper that passes the
resolved parameters through, and its outputs are byte-identical to the code before the registry
existed (the synthetic goldens run this path for real). It is the only backend cross-validated against
Imaris (`docs/SPOTMAX_VALIDATION.md`). SpotMAX's minimum-distance peak filter runs as germquant's exact
replica (`germquant.spots.peak_filter`: the same output, without the cost that grew with the square of
the candidate count; `docs/PERFORMANCE.md`).

Parameters, read from the top of the `spots:` section with the pipeline's code defaults:
`spot_radius_um` 0.3, `gauss_sigma_um` 0.08, `thresholding_method` threshold_triangle,
`effect_size_metric` spot_vs_backgr_effect_size_glass, `effect_size_min` 3.0, `merge_z_columns` true,
`z_merge_gap_um` 0.8, `z_merge_valley_frac` 0.8, `max_spot_candidates` 30000.

Three optional keys for pack recipes (`germquant.spots.detect.SPOTMAX_OPTIONAL_KEYS`) are forwarded to
the chain only when the `spots:` section sets them; they have no code default, so a run without them
makes exactly the SpotMAX calls above:

* `spot_radius_z_um`: the z radius, so SpotMAX's radii become (rz, r, r) / spacing instead of r / spacing;
* `remove_hot_pixels`: passed to SpotMAX's `preprocess_image` as `do_remove_hot_pixels`;
* `sharpen`: SpotMAX's own sharpening, as its core runs it. The preprocessed image (hot pixels removed,
  smoothed) is kept; its difference-of-Gaussians copy (`spotmax.filters.DoG_spots` with the radii and the
  nucleus labels) feeds the semantic segmentation and the peak detection, and goes to the feature step as
  `sharp_spots_image` while the feature step's image stays the unsharpened one. SpotMAX computes the
  spot-versus-background effect size on the sharpened image, so with `sharpen` the `effect_size` column
  and `effect_size_min` are on that scale (as in SpotMAX's own runs); a threshold tuned without
  `sharpen`, such as the RAD-51 value 3.0, does not carry over. The other features stay on the
  unsharpened image.

A `spots.instances` entry can set them for its own instance. `remove_hot_pixels` and `sharpen` take
`true` or `false`, and `spot_radius_z_um` a positive number; anything else (a quoted "false" included)
stops the image before any stage, with the key named (`germquant.recipe_keys`). The chain that honours
them is `germquant.spots.spotmax_options`, reached by a branch at the top of `detect_spots_spotmax`, whose
own body is not edited; with neutral values (`spot_radius_z_um` equal to `spot_radius_um`, both flags
false) it gives the default chain's tables.

### Built-in `log` (experimental classical baseline)

Status: EXPERIMENTAL. It is never selected unless `spots.backend: log` is set, it has not been
validated against Imaris counts, and its counts are not comparable with the SpotMAX defaults. It is the
first in-house seed of a detector: a classical Laplacian-of-Gaussian (LoG) blob detector built with
scipy and scikit-image conventions, wired through the same tables, flags and goldens as the validated
path so that a candidate detector can be compared head to head. Its code is
`germquant.spots.log_detect`.

The chain, per nucleus label (each label is processed in its own padded bounding-box crop, so the
memory cost stays small on a full 2600 x 2600 x 90 stack):

1. LoG response `-gaussian_laplace(crop, sigma)` with a per-axis sigma in voxels derived from one
   physical sigma, the standard deviation of a Gaussian whose full width at half maximum is the
   expected focus diameter. Bright blobs give a positive response. Anisotropy is handled by the
   per-axis sigma, not by resampling.
2. Local maxima: a voxel is a candidate when it is the maximum of the response inside an ellipsoid of
   radius `min_separation_um` (a physical radius, so z and xy are treated alike), it lies inside the
   label and its response clears the pre-gate: more than `pre_gate_sigmas` robust standard deviations
   (1.4826 times the median absolute deviation) of the response over the voxels of that nucleus, and
   always more than zero. The pre-gate keeps noise maxima of an empty or noisy nucleus out of the
   flood guard. A greedy pass in physical coordinates settles ties.
3. Flood guard: `_cap_candidates`, the SpotMAX path's guard, keeps the brightest
   `max_spot_candidates` raw voxels with the same loud warning.
4. Signal-to-background gate: the raw peak (the brightest voxel of a core ellipsoid of radius half the
   diameter) divided by the median of the raw image in a local annulus must reach `min_snr`. The
   annulus is restricted to voxels of the same nucleus label when at least 8 are available, so the
   reference is the nucleoplasm and not the dark outside of the nucleus.
5. z-split merge: `_merge_z_columns`, the SpotMAX path's merge, on the integer peak positions, with the
   raw foci channel for the valley test (the shared keys apply). The anisotropic suppression of step 2
   already removes stacked duplicates closer than the separation radius; this is the second guard.
6. Sub-voxel refinement: the intensity-weighted centroid of the background-subtracted core voxels.
7. Nucleus assignment: the label at the refined position; when that is 0 (the centroid drifted just
   outside the mask) the label of the peak voxel, which is never 0 because candidates are only taken
   inside labels. The SpotMAX path has no nearest-label fallback, spots outside every label are
   dropped, so none was added here.

Parameters under `spots.log`, all read with code defaults:

| key | code default | meaning |
|---|---|---|
| `diameter_um` | 0.5 | expected focus diameter; sets the LoG sigma and the core radius |
| `min_snr` | 2.0 | raw peak over local annulus median |
| `min_separation_um` | half the diameter | minimum distance between two maxima |
| `annulus_inner_um` | the diameter | inner radius of the local background annulus |
| `annulus_outer_um` | twice the diameter | outer radius of the local background annulus |
| `pre_gate_sigmas` | 3.0 | candidate pre-gate in robust standard deviations of the in-nucleus LoG response (0: any positive response); a finite number of at least 0, null means the default |

Shared keys read from the top of the `spots:` section: `merge_z_columns`, `z_merge_gap_um`,
`z_merge_valley_frac`, `max_spot_candidates` (same code defaults as above). Unknown keys under
`spots.log` are logged and ignored.

Column meaning for this backend: `spot_mean_intensity` is the mean raw intensity of the core voxels;
`effect_size` is a Glass's delta, the core mean minus the annulus mean in annulus standard deviations
(the same idea as SpotMAX's `spot_vs_backgr_effect_size_glass`, not the same number; infinite when the
annulus has no variance, as on a noise-free synthetic image).

Known limits, to be measured before anyone trusts a count from it: one scale only (no per-spot size, two
touching foci of different sizes resolve by separation alone); bright nucleoplasm texture can pass the
SNR gate; the annulus median is a crude background for a focus at the nuclear edge; the valley test
runs on the raw profile, so noise can keep two stacked peaks apart that SpotMAX would merge. The way to
validate it is the one the SpotMAX path went through: the same gonads with two profiles (`spots.backend`
spotmax and log), `scripts/regression_diff.py` on the two result trees, and the Imaris comparison tools
in `germquant.validate`.

### Built-in `net3d` (experimental learned detector)

Status: EXPERIMENTAL. Selected only by `spots.backend: net3d` with a model file in `spots.net3d.model`;
never the default, not cross-validated against Imaris beyond the development scores on the train and val
gonads in `docs/BENCHMARK.md` ("net3d"). It needs torch. Its code is `germquant.spots.net3d`; the model
is trained with `germquant train-spots` (`docs/BENCHMARK.md` has the recipe and the record).

The chain:

1. Inputs: the foci channel and, for a model trained with DAPI, the DNA channel (through the
   `channels` mapping of the contract above), each normalised per image to (raw - p1) / (p99.9 - p1)
   with the percentiles of the whole volume read on every fourth voxel in y and x (the model file
   records the settings).
2. A compact residual 3D U-Net (four levels, 16 to 128 channels, an xy-only first pooling for the 1.85
   anisotropy of the 0.108 x 0.108 x 0.2 um grid) gives a per-voxel probability of a focus centre.
   Inference is tiled in float32 with an overlap equal to the network's receptive radius, so the result
   does not depend on the tile size beyond float rounding; tiles without a nucleus are skipped.
3. Peaks: voxels above the threshold that are the maximum of the heatmap over an ellipsoid of
   `separation_um` (0.6 um in z, 0.3 um in xy), a greedy ellipsoidal non-maximum suppression, and a
   sub-voxel refinement by the parabola through the log probabilities on each axis.
4. Nucleus assignment: the label at the refined position; when that is background, the nearest
   labelled voxel within `assign_max_um` (0.3 um); otherwise the peak is dropped, as SpotMAX drops
   spots outside every nucleus. `contested` marks a spot with two or more nuclei within that distance.

Keys under `spots.net3d`, all read with code defaults:

| key | code default | meaning |
|---|---|---|
| `model` | none (required) | the model file `germquant train-spots` wrote (`MODELDIR/model.pt`); an absolute path, or one relative to the working directory |
| `threshold` | the model's own | peak probability threshold; by default the value chosen on the val gonads and stored in the model file |
| `device` | `auto` | `auto` (CUDA when available, honours `GERMQUANT_DEVICE`), `cpu` or `cuda` |
| `tile_zyx` | [96, 384, 384] | inference tile in voxels |
| `separation_um` | the model's own, else [0.6, 0.3] | (z, xy) radii of the suppression ellipsoid |
| `assign_max_um` | the model's own, else 0.3 | reach of the nearest-label fallback |

Shared key: `max_spot_candidates` (30000), the flood guard (the most probable peaks are kept, with the
SpotMAX path's warning). The z-split merge keys do not apply.

Columns: `spot_mean_intensity` is the mean raw intensity of a 0.25 um core, `effect_size` the log
backend's Glass's delta against a 0.5 to 1.0 um annulus; after the declared columns come `prob` (the
heatmap peak), `assignment` (`label` or `nearest`) and `contested`. A run records the model file's
sha256 in the spots stage sub-hash and in `run_manifest.json` under `spot_backend`, besides the usual
`spots:net3d_n=` and `spots:unvalidated_backend=net3d` flags. A CPU and a CUDA run agree to float
rounding, not bit for bit, so a peak exactly at the threshold can differ between devices.

```yaml
extends: config.yaml
spots:
  backend: net3d
  net3d:
    model: C:/Users/ryane/germquant_golden/models/net3d_20261010/model.pt
```

### Adding a backend

In process, for a script or a notebook:

```python
from germquant.spots import register_backend

@register_backend("mydetector")
def mydetector(image, labels, spacing, *, cfg_section, marker):
    ...
    return per_spot, per_nucleus   # the columns above, detector == "mydetector"
```

Packaged, so that `spots.backend: mydetector` works in any config once the package is installed, through
the entry-point group `germquant.spot_backends` in the package's `pyproject.toml`:

```toml
[project.entry-points."germquant.spot_backends"]
mydetector = "mypackage.germquant_backend:detect"
```

Entry points are discovered on the first lookup. A failing one (broken metadata, an import error, a
value that is not callable) is logged and skipped, never fatal, and an entry point can not replace a
name that is already registered, so the validated default can not be swapped out from outside. The
registry API is `register_backend(name, fn, replace=False)` (usable as a decorator), `get_backend(name)`,
`list_backends()` and `discover_backends(force=False)`.

A new backend earns its place the way the `log` tests do it (`tests/test_spot_backends.py`): a synthetic
stack with Gaussian foci planted at known positions inside synthetic nuclei, at least 90 percent of them
recovered, at most 10 percent extras, every match assigned to the right nucleus, and the tables checked
against the contract.

### What does not move

The default run never changes: `spots.backend` is read with the code default `spotmax`, the wrapper
calls the same code with the same arguments, the detector column still reads `spotmax` and the flag
text is unchanged, which the synthetic goldens verify cell for cell. The spots table keeps its declared
columns in their order (`germquant.schema.SPOTS`); a backend may not rename, reorder or drop one.

## Backends and system packs

System packs (docs/PACKS.md) name backends and models by name, never by callable, so the engine can tell
whether a pack can run without loading anything.

- `declared_backend_names()`, on both registries (`germquant.segment.backends` and
  `germquant.spots.backends`): the built-ins, the backends registered in this process and the entry-point
  names read from the installed packages' metadata, without loading any of them (so nothing pulls in torch
  or CuPy). `germquant pack check` greys a recipe whose backend is not declared. A declared entry point
  can still fail when it is loaded; `germquant pack check --deep` loads the backends and reports that.
  `germquant.api.list_segment_backends()` and `list_spot_backends()` return the same lists.
- `germquant.models_registry.resolve_model_name(name, expect_sha256=None)`: the helper a third-party
  backend calls for its own model key. It returns the file a registered local or release model points at,
  or, for a pinned pack model `<name>@<8 hex>`, the verified file in the model store, after the pin check
  (the file's sha256 must start with the 8 hex digits and equal the pack's digest; a name of that form that no
  installed pack declares is refused). With `expect_sha256`
  the file's digest must equal it. A pack model whose file is missing raises FileNotFoundError naming
  `germquant pack fetch-models <pack>`.
- The pack rule: a pack recipe sets `segmentation.nuclei.method: cellpose` (with whatever backend it names)
  or `classical`, never `auto`, so a missing Cellpose or a missing plugin fails loudly instead of falling
  back to the classical watershed.
- Post-stages are named the same way: `needs.post_stages` lists the post-stages a recipe needs, and `germquant
  pack check` greys the recipe where one is not offered, reading metadata only (docs/POST_STAGES.md).
