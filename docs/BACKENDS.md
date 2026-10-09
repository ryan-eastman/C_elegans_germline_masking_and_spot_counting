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

### What the output says

`segmentation_method` in every table carries the name of the backend that produced the labels
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
| `dna` | `np.ndarray`, shape (Z, Y, X) | the DNA (DAPI) channel as read from the .nd2 (uint16; float32 for a synthetic stack), not normalised |
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

Not yet pluggable: the spot stage runs SpotMAX (`germquant.spots`). This section is reserved for the
spot-detection backend contract, which the spot backend work adds.
