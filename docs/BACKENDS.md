# Pluggable backends

The pipeline's model-driven steps (nucleus segmentation today, spot detection next) run behind a
small registry, so that another segmenter can be plugged in without forking the pipeline and a
head-to-head against the current model is a config change, not a rewrite. This document is the
contract a backend has to meet.

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
stage falls back to channel 0 and flags `segment:no_<role>_channel_using_first_channel`. The montage's
outline panel sits on the same channel. A backend therefore receives "the channel to segment", which is
the DNA channel unless the profile says otherwise; its first argument keeps the name `dna` for
compatibility.

`segmentation.nuclei.cellpose_model` takes a file path, a Cellpose built-in name (`cpsam`) or the name of
a registered model: a release asset `get-model` has downloaded, or a model added with `germquant
register-model` (`models_registry.register_local_model`), which resolves to the registered file and puts
its sha256 in the run manifest as always.

### What the output says

`segmentation_method` in `<image_id>__image_summary.csv` (the only table that carries it) names the backend that produced the labels
(`cellpose`, `classical`, or the plugin's name), or `classical` after a fallback; it still says
`cellpose` for every validated run. The segment stage's sub-hash in `run_manifest.json` covers the
whole `segmentation` section, so a run with a different backend or different backend keys gets a
different sub-hash. The manifest's `model_sha256` is still the digest of
`segmentation.nuclei.cellpose_model`: a backend that keeps its model under that key gets the
provenance for free, a backend with its own model key should log its own digest.

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

### The built-ins

| name | implementation |
|---|---|
| `cellpose` | `germquant.segment.backends.cellpose_backend`: `nuclei._cellpose` (Cellpose-SAM in 3D, anisotropy from the spacing, the fine-tuned germline model) with the explicit device from `germquant.device` and one retry on the CPU when the accelerator fails (out of memory, an operator Metal cannot run) |
| `classical` | `germquant.segment.backends.classical_backend`: `nuclei._classical`, Gaussian smoothing, Li threshold and a distance-transform watershed, all spacing-aware. No model and no GPU; a baseline, not a substitute for the fine-tuned model on crowded pachytene nuclei |

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
Imaris (`docs/SPOTMAX_VALIDATION.md`).

Parameters, read from the top of the `spots:` section with the pipeline's code defaults:
`spot_radius_um` 0.3, `gauss_sigma_um` 0.08, `thresholding_method` threshold_triangle,
`effect_size_metric` spot_vs_backgr_effect_size_glass, `effect_size_min` 3.0, `merge_z_columns` true,
`z_merge_gap_um` 0.8, `z_merge_valley_frac` 0.8, `max_spot_candidates` 30000.

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
