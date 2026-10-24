# Post-stages: adding a measurement without editing germquant

A post-stage is a measurement step that a lab adds as its own installed package. It runs after germquant has
finished an image, reads what that run wrote (the tables, the nucleus label image and the image itself), and
writes new tables of its own. It never edits an existing table, label image or completion marker, it never
changes `config_hash`, and when it fails the failure is recorded and the image stays done.

Post-stages are the second half of the "plugins in two senses" promise: segmentation and spot backends
(docs/BACKENDS.md) replace a step of the pipeline, post-stages add one after it. This page is the contract. The
design is section 5.6 of docs/DESIGN_system_packs.md.

## What a post-stage is

Three pieces, all in the lab's package:

1. An entry point in the group `germquant.post_stages`. Its name is the post-stage's name, its value the function:
   `ring_intensity = "germquant_ring_intensity:measure"`.
2. A declaration file, `germquant_post_stage.yaml`, inside the top-level package that the entry point names. It
   says what the post-stage reads, which tables and columns it writes, and its parameters.
3. The function, `measure(inputs, params) -> {table name: pandas DataFrame}`.

germquant finds a post-stage from the package's metadata and its declaration file alone: listing post-stages,
checking a pack's `needs.post_stages` and validating a `post_stages.yaml` never import the post-stage's code,
torch, CuPy, SpotMAX or Cellpose. The code is imported only when a run calls it, after germquant has loaded its
own pipeline (torch first, so its CUDA libraries load before anything else's).

Post-stages, like third-party backends, are installed into `.venv_packs` or a lab's own environment, never into
the validated `.venv` (rule R6 of the design). No code ever runs from a project or pack folder: a post-stage is
always an installed package.

## The declaration file

`germquant_post_stage.yaml` sits next to the package's `__init__.py`. One file may declare several post-stages
(one entry point each).

```yaml
format: 1
post_stages:
  ring_intensity:                         # the entry point's name
    version: "0.1.0"                      # quoted, a semantic version; part of the sub-hash
    title: "Foci intensity in a shell around each nucleus"     # optional, one line for listings
    reads: {roles: [foci], tables: [], labels: true, image: true}
    writes:
      rings: [nucleus_id, ring_voxels, ring_mean, inside_mean]
    params:
      - {key: radius_um, type: float, unit: um, min: 0.1, max: 5.0, default: 1.0, label: "Shell width"}
```

`reads` (every key optional):

| key | meaning |
|---|---|
| `roles` | the channel roles of the recipe's channel map that the post-stage reads from the image (`dna`, `foci`, ...) |
| `tables` | the main run's tables it reads, by name: `nuclei`, `spots`, `image_summary`, `spots_<instance>`, ... (never another post-stage's tables) |
| `labels` | `true` when it reads the nucleus label image |
| `image` | `true` when it reads the image; then `roles` names the channels, and `roles` without `image: true` is refused |

`writes` maps each table the post-stage writes to its columns, in order. A table is written as
`<image_id>__post_<name>_<table>.csv` with exactly these columns.

`params` declares each parameter like a recipe tunable in a pack (docs/PACKS.md):

| type | value | other keys |
|---|---|---|
| `float` | a number, stored as a float | `min`, `max`, `unit` |
| `int` | an integer (never a float, never true or false) | `min`, `max`, `unit` |
| `bool` | `true` or `false` | |
| `choice` | one of `choices`, with its YAML type (`1` is not `"1"`) | `choices` (required) |
| `str` | text | |

Every parameter needs a `default`, and `label` (text for the parameter form) is optional. Ranges are inclusive.

Names of post-stages, tables and parameters are lower-case snake case of at most 40 characters: letters, digits
and single underscores, starting with a letter (`ring_intensity`, not `RingIntensity`, `ring__intensity` or
`ring_`). That keeps every post table's file name apart from the main run's tables on every file system.
Column and role names start with a letter and hold letters, digits and underscores. A declaration that breaks a
rule makes its post-stage unusable: it is not offered, and `python -m germquant.post_stages list --all` shows
why. So does a package that holds no declaration file, an entry point that names a module instead of a package,
and a name that two installed distributions both declare.

## The function

```python
def measure(inputs, params):
    ...
    return {"rings": dataframe}
```

`params` holds every declared parameter: the value `post_stages.yaml` gives, else the default. `inputs` holds:

| name | what it is |
|---|---|
| `inputs.image_id`, `inputs.out_dir` | the image's id and its results folder (a Path) |
| `inputs.spacing` | (dz, dy, dx) in micrometres as the run recorded it (the `voxel_d*_um` columns), or None |
| `inputs.roles` | the run's channel roles, `{role: channel index or None}` |
| `inputs.tables(name)` | a declared table of the main run, as a new DataFrame each call |
| `inputs.labels()` | the nucleus label image, int32 (Z, Y, X); a single plane is (1, Y, X) |
| `inputs.image(role)` | the channel of a declared role, (Z, Y, X), read through the run's own reader at the run's `xy_stride` and `z_range`, the same grid as the labels |

The arrays are read-only (copy one before changing it), and the image is read once per image however many
post-stages read it. Reading something the declaration does not name (an undeclared table or role, the labels
without `labels: true`) raises an error, which fails that post-stage.

The function returns one DataFrame per declared table, with exactly the declared columns (any order; they are
written in the declared order). germquant writes the files; the function must not write into the results folder
itself. It should import numpy, pandas or anything heavy inside the function, not at the top of its module, so
that importing the package stays cheap.

## Switching post-stages on: post_stages.yaml

A run switches post-stages on through a file of its own, `post_stages.yaml`, a list of `{name, params}`:

```yaml
- name: ring_intensity
  params:
    radius_um: 0.8
```

It sits in a project's root next to `recipe.yaml`, or ships in a pack as `post_stages/<recipe>.yaml`, or is
given to a run as a file; a project's file, or one given to the run, replaces the pack's. An empty list switches
every post-stage off (and so also replaces a pack's file). germquant's config loader never reads the file:
switching a post-stage on or changing its parameters changes neither `config_hash`, nor a stage hash, nor the
completion marker, so it reruns neither the image nor anything but the post-stage.

The file is checked against the post-stages installed in the environment that runs: each entry must name a
usable post-stage, list it once, and give only declared parameters, each of its declared type and inside its
range. Write it with germquant, which checks it first and writes it in one canonical form (ASCII, with YAML
escapes for any other character; LF line endings; the parameter keys sorted; float parameters as floats; an
atomic replace). From Python, the same line works in PowerShell, cmd, bash and zsh:

```
python -c "from germquant.post_stages import write_post_stages_file; write_post_stages_file('post_stages.yaml', [{'name': 'ring_intensity', 'params': {'radius_um': 0.8}}])"
```

or from the command line (see below). A file is also fine written by hand; the runner checks it the same way.

## Running post-stages on a results folder

Until the run commands take post-stages, run them on a finished results folder from Python, with the
environment that has them installed:

```python
from germquant import api
from germquant.post_stages import run_post_stages

cfg = api.resolve_config(profile="generic")       # the config the folder was made with
res = run_post_stages("results/img_001", "img_001", "post_stages.yaml", cfg=cfg, image_path="img_001.nd2")
print(res["ran"], res["skipped"], res["failed"])
```

`entries` is the path of a `post_stages.yaml` or the list itself. `stack=` hands over an image stack already in
memory (then the image is not read again); `resume=False` reruns every post-stage. The call returns the names it
ran, skipped and failed and the record it wrote; it never raises for a post-stage's failure.

It runs only on a folder that a run has finished (one with its completion marker), made with the same
`config_hash` as `cfg`, and never on a timepoint of a time series: post-stages are refused for time series in
this release. Each refusal is recorded for every post-stage.

## What is written: post tables and the record

For each post-stage that ran, one CSV per declared table, `<image_id>__post_<name>_<table>.csv`, and for every
call the record `<image_id>__post.json`:

```json
{
 "format": 1,
 "image_id": "img_001",
 "post_stages_file_sha256": "1fd429c6...",
 "stages": {
  "ring_intensity": {
   "version": "0.1.0",
   "params": {"radius_um": 0.8},
   "sub_hash": "c53e223939a1",
   "done": true,
   "error": null,
   "tables": ["img_001__post_ring_intensity_rings.csv"],
   "finished": "2026-10-24T02:14:36.376158+00:00",
   "inputs": {"read_sub_hash": "0123456789ab", "spacing": [0.5, 0.2, 0.2], "tables": {},
              "roles": {"foci": 1}, "labels": "083de45c...", "image": {"xy_stride": 1, "z_range": null}}
  }
 }
}
```

The `sub_hash` (12 hex digits, like the stage sub-hashes of `run_manifest.json`) covers:
- the post-stage's version, its declared reads and writes, and its parameters (defaults filled in);
- the digest of `post_stages.yaml` (`post_stages_file_sha256`: the sha256 of the file's bytes with CRLF turned
  into LF, or, for a list given from Python, of the text germquant would write from it);
- the digests of its inputs (`inputs`): the main run's read sub-hash, the spacing, the channel index of each
  role it reads, the sha256 of each table it declares (CRLF turned into LF) and of the label image when it
  declares it, and the run's `xy_stride` and `z_range` when it reads the image.

The raw image is not hashed again: the read sub-hash and the tables already change when it does. So a new main
run (whose tables carry a new `run_timestamp`), a changed parameter, a new version of the post-stage or another
`post_stages.yaml` each give a new sub-hash.

On the next call a post-stage is skipped when its record says `done`, its `sub_hash` is the same and its tables
are still there; anything else runs again, the post-stage only. The record lists the post-stages of the current
file: one taken out of `post_stages.yaml` leaves the record on the next call, while its old tables stay on disk.

A failure is written into the record with `done: false` and the reason in `error`, and the other post-stages
still run: a post-stage that is not installed, a parameter outside its range, an exception from the function, a
table or column outside the declaration, or a function that wrote into the results folder itself. Nothing of the
main run changes either way: its tables, label images, stage record and completion marker stay byte for byte as
they were, and the image stays done.

## From the command line

`python -m germquant.post_stages` loads no torch, so it is quick and safe in any environment. It is how
napari-germquant lists post-stages and writes a project's `post_stages.yaml`: through the engine environment,
because only that environment has the post-stages installed.

```
python -m germquant.post_stages list                 # the post-stages this environment offers
python -m germquant.post_stages list --all --json    # also the unusable ones, with the reason; as JSON
```

`write DIR --entries JSON` checks a JSON list of `{name, params}` against this environment's post-stages and
writes `DIR/post_stages.yaml`. It exits 0 when it wrote the file, 2 when it refused (malformed JSON, no such
folder, an entry that cannot run here; every problem is printed and nothing is written) and 1 when the file
could not be written. `--json` prints the outcome (`written`, `path`, `sha256`, `entries`, or `problems`).

macOS and Linux:
```
python -m germquant.post_stages write myproject --entries '[{"name": "ring_intensity", "params": {"radius_um": 0.8}}]'
```

Windows PowerShell 5.1 passes inner double quotes to a program only when they are escaped (PowerShell 7.3
and later take the macOS form):
```
python -m germquant.post_stages write myproject --entries '[{\"name\": \"ring_intensity\", \"params\": {\"radius_um\": 0.8}}]'
```

## Worked example: foci intensity in a shell around each nucleus

A complete package. It measures, for every nucleus, the mean foci intensity in a shell of a given width around
the nucleus and inside it.

```
germquant-ring-intensity/
  pyproject.toml
  src/germquant_ring_intensity/__init__.py
  src/germquant_ring_intensity/germquant_post_stage.yaml
```

`pyproject.toml` (the package data line ships the declaration in the wheel):

```toml
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "germquant-ring-intensity"
version = "0.1.0"
requires-python = ">=3.11"
dependencies = ["germquant", "numpy", "pandas", "scipy"]

[project.entry-points."germquant.post_stages"]
ring_intensity = "germquant_ring_intensity:measure"

[tool.setuptools.packages.find]
where = ["src"]

[tool.setuptools.package-data]
germquant_ring_intensity = ["germquant_post_stage.yaml"]
```

`src/germquant_ring_intensity/germquant_post_stage.yaml`: the declaration shown above, under "The declaration
file".

`src/germquant_ring_intensity/__init__.py`:

```python
"""A germquant post-stage: mean foci intensity in a shell around each nucleus and inside it."""


def measure(inputs, params):
    import numpy as np
    import pandas as pd
    from scipy import ndimage as ndi

    labels = inputs.labels()                          # (Z, Y, X) int32, read-only
    foci = inputs.image("foci").astype(np.float64)    # the foci channel at the run's stride and z range
    spacing = inputs.spacing or (1.0, 1.0, 1.0)       # (dz, dy, dx) in um, as the run recorded it
    radius = params["radius_um"]
    pad = [int(np.ceil(radius / s)) + 1 for s in spacing]
    rows = []
    for nucleus_id, box in enumerate(ndi.find_objects(labels), start=1):
        if box is None:                               # a label number the image does not use
            continue
        box = tuple(slice(max(0, b.start - p), min(n, b.stop + p)) for b, p, n in zip(box, pad, labels.shape))
        inside = labels[box] == nucleus_id
        dist = ndi.distance_transform_edt(~inside, sampling=spacing)
        ring = (dist > 0) & (dist <= radius) & (labels[box] == 0)
        rows.append({"nucleus_id": nucleus_id, "ring_voxels": int(ring.sum()),
                     "ring_mean": float(foci[box][ring].mean()) if ring.any() else float("nan"),
                     "inside_mean": float(foci[box][inside].mean())})
    return {"rings": pd.DataFrame(rows, columns=["nucleus_id", "ring_voxels", "ring_mean", "inside_mean"])}
```

Install it into the environment that will run it, never the validated `.venv`: a lab's own environment, or a
clone of the engine environment such as `.venv_packs` (the name the design uses; below it sits in the
current folder). Install it editable while you develop it, and check that germquant sees it:

Windows:
```
.venv_packs\Scripts\python.exe -m pip install -e germquant-ring-intensity
.venv_packs\Scripts\python.exe -m germquant.post_stages list
```

macOS and Linux:
```
.venv_packs/bin/python -m pip install -e germquant-ring-intensity
.venv_packs/bin/python -m germquant.post_stages list
```

`list` prints `ring_intensity 0.1.0  (germquant-ring-intensity)` and its title. Then switch it on with a
`post_stages.yaml` and run it on a finished results folder, as above. The folder gains
`img_001__post_ring_intensity_rings.csv`:

```
nucleus_id,ring_voxels,ring_mean,inside_mean
1,452,47.89627305223747,51.912286962615326
2,720,49.82049357617895,46.06454626303166
3,688,47.62183161465408,48.86045501652794
```

and the record `img_001__post.json` shown earlier. Change `radius_um` and call it again: only the post-stage
reruns. Like a new backend, a post-stage's numbers are the lab's own: compare it with a known measurement on a
few images before relying on it, and name it and its version in the methods.

## Testing a post-stage

Register the function in the test's process instead of installing it:

```python
from pathlib import Path

import yaml
from germquant.post_stages import register_post_stage, run_post_stages, unregister_post_stage

text = Path("src/germquant_ring_intensity/germquant_post_stage.yaml").read_text(encoding="utf-8")
declaration = yaml.safe_load(text)
register_post_stage("ring_intensity", measure, declaration["post_stages"]["ring_intensity"])
try:
    res = run_post_stages(out_dir, image_id, [{"name": "ring_intensity"}], cfg=cfg, image_path=image)
finally:
    unregister_post_stage("ring_intensity")
```

`register_post_stage` checks the declaration by the same rules as an installed one. A registration lives only
in that process, so a child process (the plugin's engine runner, for example) does not see it.
`germquant.post_stages.testing` holds germquant's own fixture post-stage, `fixture_nuclei` (per nucleus: voxels,
volume, a foci statistic times a factor, the spot count), in both forms: `register_fixture()` registers it in
process, and `write_fixture_distribution(folder)` writes it as an installed distribution (a `.dist-info` folder
with the entry point, plus the package) that any interpreter with `folder` on `PYTHONPATH` discovers. The
engine's and the plugin's tests both use that one helper.

## Reference

| name | what it does |
|---|---|
| `germquant.post_stages.list_post_stages(include_unusable=False)` | the usable post-stages, from metadata: name, version, dist, title, reads, writes, params, source |
| `load_post_stages_file(path)` | the entries of a `post_stages.yaml` (ValueError for a malformed file) |
| `validate_post_stages(entries, declarations=None)` | every problem, against this environment's post-stages or the given declarations |
| `write_post_stages_file(path, entries, declarations=None)` | check, then write the canonical file |
| `post_stages_sha256(path or entries)` | the digest that enters the sub-hashes |
| `run_post_stages(out_dir, image_id, entries, cfg=, image_path=, stack=None, resume=True)` | run them on a finished results folder |
| `register_post_stage(name, fn, declaration, replace=False)`, `unregister_post_stage(name)` | in-process registration |
| `read_post_record(out_dir, image_id)`, `post_table_name(image_id, name, table)` | the record and the table file names |
| `python -m germquant.post_stages list [--json] [--all]`, `write DIR --entries JSON [--json]` | the torch-free command line |
