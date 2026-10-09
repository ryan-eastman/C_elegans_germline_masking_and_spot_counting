# Pluggable backends

A backend is one algorithm behind one pipeline stage, selected by a config key that is read with a
code default. The default backend of every stage is the validated one, and its outputs do not move
when another backend is added: the config files are not edited, the schema only grows, and the
synthetic goldens (`tests/test_golden.py`) gate every change. Each stage that can swap its algorithm
has a section below.

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
   radius `min_separation_um` (a physical radius, so z and xy are treated alike), its response is
   positive and it lies inside the label. A greedy pass in physical coordinates settles ties.
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
