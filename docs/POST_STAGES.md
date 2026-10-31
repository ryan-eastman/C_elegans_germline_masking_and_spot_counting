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

The `version` is what resume trusts: a post-stage's sub-hash covers its version, never its code. A code change
that can change the numbers needs a new version (or a call with `resume=False`); otherwise, with an editable
install for example, an edit to the function reruns nothing on folders it has already done.

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
| `inputs.spacing` | (dz, dy, dx) in micrometres as the run recorded it (the `voxel_d*_um` columns of its image_summary table, on every route), or None when it recorded none; on a single plane, whose tables record dz as NaN, dz is 1.0, the placeholder germquant's own array code uses there |
| `inputs.roles` | the run's channel roles, `{role: channel index or None}` |
| `inputs.tables(name)` | a declared table of the main run, as a new DataFrame each call |
| `inputs.labels()` | the nucleus label image, int32 (Z, Y, X); a single plane is (1, Y, X) |
| `inputs.image(role)` | the channel of a declared role, (Z, Y, X), read through the run's own reader at the run's `xy_stride` and `z_range`, the same grid as the labels |

The arrays are read-only (copy one before changing it), and the image is read once per image however many
post-stages read it. Reading something the declaration does not name (an undeclared table or role, the labels
without `labels: true`) raises an error, which fails that post-stage.

The function returns one DataFrame per declared table, with exactly the declared columns (any order; they are
written in the declared order). germquant writes the files; the function must not write into the results folder
itself, nor into a subfolder of it. It should import numpy, pandas or anything heavy inside the function, not at
the top of its module, so that importing the package stays cheap.

## Switching post-stages on: post_stages.yaml

A run switches post-stages on through a file of its own, `post_stages.yaml`, a list of `{name, params}`:

```yaml
- name: ring_intensity
  params:
    radius_um: 0.8
```

It sits in a project's root next to `recipe.yaml`, or ships in a pack as `post_stages/<recipe>.yaml`, or is
given to a run with `--post-stages FILE`; a project's file, or one given to the run, replaces the pack's (see
"Running post-stages" below). An empty list switches every post-stage off (and so also replaces a pack's file).
germquant's config loader never reads the file: switching a post-stage on or changing its parameters changes
neither `config_hash`, nor a stage hash, nor the completion marker, so it reruns neither the image nor anything
but the post-stage.

The file is checked against the post-stages installed in the environment that runs: each entry must name a
usable post-stage, list it once, and give only declared parameters, each of its declared type and inside its
range. Two post-stages whose tables would share a file name are refused together: `a` writing table `b_c` and
`a_b` writing table `c` would both write `<image_id>__post_a_b_c.csv`, so switch on only one of them. Write the
file with germquant, which checks it first and writes it in one canonical form (ASCII, with YAML escapes for any
other character; LF line endings; the parameter keys sorted; float parameters as floats; an atomic replace).
From Python, the same line works in PowerShell, cmd, bash and zsh:

```
python -c "from germquant.post_stages import write_post_stages_file; write_post_stages_file('post_stages.yaml', [{'name': 'ring_intensity', 'params': {'radius_um': 0.8}}])"
```

or from the command line (see below). A file is also fine written by hand; the runner checks it the same way.

## Running post-stages

The run commands take post-stages, in the environment that has them installed:

```
germquant run   IMAGE  --config CFG --out OUT --post-stages post_stages.yaml
germquant batch FOLDER --config CFG --out OUT --post-stages post_stages.yaml [--resume]
germquant project run DIR                        # the project's post_stages.yaml
```

and from Python `api.run(IMAGE, ..., post_stages="post_stages.yaml")` or `api.run_loaded(IMAGE, cfg, OUT,
post_stages=[{"name": "ring_intensity", "params": {"radius_um": 0.8}}])`: a path or the list itself.

**Which list a run uses.** `--post-stages FILE` (or `post_stages=`) when given; in a project, its own
`post_stages.yaml`; else, for a pack recipe (`--pack`, `api.run(pack=)` or a project's recipe), the
`post_stages/<recipe>.yaml` the pack ships, when it ships one; else none. So a project's file, or one given to the
run, replaces the pack's, and a file holding `[]` switches every post-stage off. Without any of them nothing
changes: no post file is written, and a run of a recipe file or a profile does not even import
germquant.post_stages.

**Checked before the first image.** The list in force is read and checked against this environment's
post-stages before any image is processed (as `validate_post_stages` checks it). A list that cannot be read, or
holds an entry that cannot run here, refuses the run: `run` and `batch` print one `REFUSED:` line and exit 2
(`REFUSED: post-stages <file>: <problems>` for entries that cannot run here or a file that cannot be opened,
`REFUSED: <file> cannot be read as YAML ...` or `REFUSED: <file>: post-stages must be a list ...` for one that is
not a list of {name, params}), `project run` refuses with its exit 2 before its run log exists. With
`--no-pack-check` a pack's own list that cannot run here is a warning instead, and each image records the failure
in its `__post.json`. api.run and api.run_loaded refuse only a list they cannot read (ValueError, or OSError for a
file that cannot be opened); a post-stage that cannot run there is recorded as failed, image by image.

**When they run.** Right after process_image has written the image's folder, on what it wrote, with the stack the
run read handed over (the image is read once). `germquant run` and `batch` print a line per post-stage (done, up to
date, or FAILED with the reason), and api.run's result gains `post`, what run_post_stages returns. A failure never
fails the image or the batch: it is recorded in `__post.json`, and the image stays done.

**Resume.** `batch --resume` and `project run` (which resumes by default) skip an image whose main run is done
and run only its due post-stages: those whose sub-hash changed (a parameter, the file, a new version of the
post-stage, or an input). That is the rerun "of the post-stage only": the main tables, label images, stage record
and completion marker stay byte for byte as they were; the image is read again only for a post-stage
that reads it. A resumed batch without `--post-stages` runs no post-stage and leaves the post files as they
are. A run that processes the image again (`germquant run`, a batch without `--resume`, or
`project run --no-resume`) gives it a new main run, so every post-stage listed reruns on it.

**Time series.** Post-stages are refused for time series in this release. `germquant run IMAGE --t ...` with a
list that switches a post-stage on (also a pack recipe's own) prints `REFUSED: post-stages are refused for time
series ...` and exits 2 before any timepoint; api.run and api.run_loaded given `t` raise ValueError with the same
message; `project run` fails each series with it (the project's other images run). An empty list switches them
off, so `--post-stages` with a file holding `[]` runs such a series without them.

**Collecting.** `germquant collect` (and `batch --resume`, and `project collect`) stack the post tables into
`batch_post_<name>_<table>.csv` next to the other batch tables, with `image_id` as the first column (a post table
holds only its declared columns) and a sample sheet's own columns joined as in every batch table. Only the tables
in force are stacked: those an image's `__post.json` lists for a post-stage that is done and that ran on the
folder's current main run (its recorded `main_run` and table digests are the folder's now). The tables of a
post-stage switched off or failed since, or made before the main run was redone, stay on disk (nothing of a
user's is deleted) and are left out, each with a warning naming the file and the reason.

**Packs.** A pack recipe that needs a post-stage names it in `needs.post_stages`; `germquant pack check` greys
the recipe where that post-stage is not installed or cannot be used, reading metadata and the declaration file
only, and `run --pack` and `batch --pack` refuse it unless `--no-pack-check` is given. The pack's
`post_stages/<recipe>.yaml` switches post-stages on for that recipe's runs (docs/PACKS.md).

**napari-germquant** lists the post-stages and writes a project's `post_stages.yaml` through the engine
environment (`python -m germquant.post_stages list --json` and `write`, below), because only that environment
has the post-stages installed.

### On a finished results folder, from Python

The same runner works on a folder a run finished earlier:

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
   "sub_hash": "1823ae360f3b",
   "done": true,
   "error": null,
   "tables": ["img_001__post_ring_intensity_rings.csv"],
   "finished": "2026-10-24T03:28:37.963264+00:00",
   "inputs": {"read_sub_hash": "0123456789ab", "main_run": "2026-10-23T00:00:00+00:00",
              "spacing": [0.5, 0.2, 0.2], "tables": {}, "roles": {"foci": 1}, "labels": "083de45c...",
              "image": {"xy_stride": 1, "z_range": null}}
  }
 }
}
```

The `sub_hash` (12 hex digits, like the stage sub-hashes of `run_manifest.json`) covers:
- the post-stage's version, its declared reads and writes, and its parameters (defaults filled in);
- the digest of `post_stages.yaml` (`post_stages_file_sha256`: the sha256 of the file's bytes with CRLF turned
  into LF, or, for a list given from Python, of the text germquant would write from it);
- the digests of its inputs (`inputs`): the main run's read sub-hash and its `run_timestamp` (`main_run`, from
  the image_summary table), the spacing, the channel index of each role it reads, the sha256 of each table it
  declares (CRLF turned into LF) and of the label image when it declares it, and the run's `xy_stride` and
  `z_range` when it reads the image.

The raw image is not hashed again, and the read sub-hash does not cover pixels: a new main run restales every
post-stage through its `run_timestamp`, also one that reads no table. So a new main run, a changed parameter, a
new version of the post-stage or another `post_stages.yaml` each give a new sub-hash.

On the next call a post-stage is skipped when its record says `done`, its `sub_hash` is the same and its tables
are still there; anything else runs again, the post-stage only. The record lists the post-stages of the current
file: one taken out of `post_stages.yaml` leaves the record on the next call, while its old tables stay on disk.
An empty list (every post-stage switched off) leaves a record with no stages, and a folder that never had a
record gets none.

A failure is written into the record with `done: false` and the reason in `error`, and the other post-stages
still run: a post-stage that is not installed, a parameter outside its range, an exception from the function, a
table or column outside the declaration, or a function that wrote into the results folder itself (a subfolder
included). Nothing of the main run changes either way: its tables, label images, stage record and completion
marker stay byte for byte as they were, and the image stays done. Every stage record holds the same keys; `inputs`
is null for a post-stage that was refused before it could run (not installed, a parameter outside its range, a
folder the runner refuses).

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
dependencies = ["numpy", "pandas", "scipy"]

[project.entry-points."germquant.post_stages"]
ring_intensity = "germquant_ring_intensity:measure"

[tool.setuptools.packages.find]
where = ["src"]

[tool.setuptools.package-data]
germquant_ring_intensity = ["germquant_post_stage.yaml"]
```

germquant is not among the dependencies: it is not on PyPI, so pip could only fetch someone else's package of
that name. Install the post-stage into the environment that already runs germquant (below); the function never
imports germquant anyway.

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
    if inputs.spacing is None:                        # a shell in um needs the run's voxel size
        raise ValueError("the run recorded no voxel size")
    spacing = inputs.spacing                          # (dz, dy, dx) in um, as the run recorded it
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
`post_stages.yaml` and run an image with it (`germquant run IMAGE ... --post-stages post_stages.yaml`), or run
it on a finished results folder, as above. The folder gains `img_001__post_ring_intensity_rings.csv`:

```
nucleus_id,ring_voxels,ring_mean,inside_mean
1,452,47.89627305223747,51.912286962615326
2,720,49.82049357617895,46.06454626303166
3,688,47.62183161465408,48.86045501652794
```

and the record `img_001__post.json` shown earlier. Change `radius_um` and run again with `--resume` (a batch)
or call run_post_stages again: only the post-stage reruns.

Like a new backend, a post-stage's numbers are the lab's own: compare it with a known measurement on a few images
before relying on it, and name it and its version in the methods.

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
| `germquant.post_stages.list_post_stages(*, include_unusable=False)` | the usable post-stages, from metadata: name, version, dist, title, reads, writes, params, source |
| `load_post_stages_file(path)` | the entries of a `post_stages.yaml` (ValueError for a malformed file) |
| `validate_post_stages(entries, declarations=None)` | every problem, against this environment's post-stages or the given declarations |
| `write_post_stages_file(path, entries, *, declarations=None)` | check, then write the canonical file |
| `post_stages_sha256(path or entries)` | the digest that enters the sub-hashes |
| `run_post_stages(out_dir, image_id, entries, *, cfg, image_path, stack=None, resume=True)` | run them on a finished results folder |
| `register_post_stage(name, fn, declaration, *, replace=False)`, `unregister_post_stage(name)` | in-process registration |
| `read_post_record(out_dir, image_id)`, `post_table_name(image_id, name, table)` | the record and the table file names |
| `python -m germquant.post_stages list [--json] [--all]`, `write DIR --entries JSON [--json]` | the torch-free command line |
| `germquant run ... --post-stages FILE`, `germquant batch ... --post-stages FILE`, `germquant project run DIR` | the run commands (above) |
| `germquant.api.run(..., post_stages=None)`, `run_loaded(..., post_stages=None)`, `list_post_stages(*, include_unusable=False)` | the same from Python; `list_post_stages` is the listing above |
| `entries_of(path or entries)` | the entries of a list, read (ValueError, or OSError, for one that cannot be read) |
| `post_stages_in_force(given=None, *, project_root=None, pack_folder=None, recipe=None)`, `pack_post_stages_path(pack_folder, recipe)` | which list a run uses, and the list a pack ships for a recipe |
| `post_table_problems(out_dir, image_id, names)` | why collect leaves a post table out (a table it does not name is in force) |
