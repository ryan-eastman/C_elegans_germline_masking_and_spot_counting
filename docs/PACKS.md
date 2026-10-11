# System packs

A system pack is a folder of checked files that you pick by name: the recipes (settings) and the pinned models
for one kind of sample. `germquant run IMAGE --pack celegans_germline` runs the C. elegans recipe; an outside lab
picks a pack made for its own cells. A pack holds no code and no model weights, only text, cards and a small
example, so it is safe to install from anyone.

Packs are new in germquant 0.2. The design is docs/DESIGN_system_packs.md (section 3).

## The packs that ship with the engine

| pack | what for | recipes |
|---|---|---|
| celegans_germline | the C. elegans germline recipes, exactly as the frozen engine config files define them | default, rad51_foci, segmentation_only, n2_sc, ccw77_partition, live_syp2 |

The outside packs (cultured cells, budding yeast) arrive in a later package. The test pack used by the test
suite lives in tests/fixtures/packs/test_pack.

## Using a pack

List what this engine finds, read a pack's card, and check that a recipe can run here:

    germquant pack list
    germquant pack show celegans_germline
    germquant pack check celegans_germline --recipe rad51_foci

`pack check` names everything that stops a recipe and how to fix it: a missing Python package (with the extra
to install), a model that is not in the model store, an engine outside the pack's version range, a pack file
that no longer matches its checksums, a tunable (a setting the parameter form offers) that the recipe does not
set. It exits with code 2 when the recipe cannot run or the pack is not found. With `--image FILE` it also checks the readers extra and warns when
the image's pixel size or z step lies outside what the pack expects; it reads .nd2 and TIFF headers itself, and
for a CZI, LIF or OME-Zarr file it says the size was not checked, because reading those headers loads the
readers. `--deep` also loads the pipeline and the backends, reports the device (CUDA, Apple Metal or CPU) and
reads any image header.

Fetch the models a recipe needs, then run it:

    germquant pack fetch-models celegans_germline
    germquant run GONAD.nd2 --pack celegans_germline --out results/GONAD
    germquant batch NAS_FOLDER --pack celegans_germline --recipe rad51_foci --out NAS_FOLDER_results

On a Mac the same commands work in the Terminal (Applications > Utilities > Terminal) after activating the
environment (`source .venv/bin/activate`); on Windows use the command prompt or PowerShell after
`.venv\Scripts\activate`.

`--pack ID` picks the newest installed version that fits this engine; `--pack ID@VERSION` pins one. `--recipe R`
picks a recipe (the pack's default otherwise). `run --pack` and `batch --pack` refuse a recipe that `pack check`
reports unusable; `--no-pack-check` runs it anyway.

The same numbers: a pack recipe is a profile file that the engine reads like any other, so

    germquant run GONAD.nd2 --pack celegans_germline

and

    germquant run GONAD.nd2 --config config/config.yaml

give the same config_hash (84bdda5a5230), the same stage hashes and the same tables. A pack run only adds two
records to run_manifest.json and to `<image_id>__stages.json`: `pack` (id, version, recipe, where it was found,
the pack's checksum, object kind, dims and status) and `runtime` (the device asked for and the one used, torch,
CUDA, python, platform), plus `post_run.segment_device_used` in the stage record (the device Cellpose finally
ran on, after any retry on the CPU). These never enter a table, the completion marker or a hash, so a `--pack`
run and a `--config` run of the same recipe resume each other.

From Python:

    from germquant import api
    res = api.run("GONAD.nd2", pack="celegans_germline", recipe="rad51_foci", out="results/GONAD")

A program that must stay light (a GUI refreshing a pack card) imports `germquant.packs` instead of
`germquant.api`: listing and checking packs never loads torch, Cellpose, SpotMAX or CuPy. The command
`python -m germquant.packs list --json` (also `check`, `fetch-models`, `install`, `verify`, `lint`) is the
torch-free entry the napari plugin uses; `python -m germquant.packs --version` prints the engine version, git
sha and python.

## Rows of the nuclei table

Table and column names never change, so in a pack whose objects are whole cells, one row of the "nuclei" table
is one cell. The pack card and `pack show` say which:

| object_kind | one row of the nuclei table is | n_nuclei counts | nucleus_id names |
|---|---|---|---|
| nucleus | a nucleus | nuclei | a nucleus |
| cell | a whole cell | cells | a cell |

## Where packs live

`germquant pack list` prints the search order. Packs are found, in this order:

1. the packs of an open project (`<project>/packs/`, from a later package);
2. the built-in packs, `packs/<id>/` in the engine checkout, or under `$GERMQUANT_CONFIG_ROOT`;
3. each folder of `$GERMQUANT_PACKS` (several folders separated by `;` on Windows, `:` on a Mac); the first
   folder is where `germquant pack install` writes;
4. the per-user folder: `%LOCALAPPDATA%\germquant\packs` on Windows, `~/Library/Application Support/germquant/packs`
   on a Mac, `~/.local/share/germquant/packs` on Linux.

Installed packs sit at `<folder>/<id>/<version>/`; a folder of 3 or 4 may also hold a pack you are writing, at
`<folder>/<id>/`. The first copy of an id and version wins, except that a built-in pack cannot be shadowed: a
different copy elsewhere is listed as shadowed and ignored (`pack list --all` shows it). A pack with a broken
manifest is listed as unusable with the reason and never stops a run that does not use it.

The built-in packs are found only in an editable checkout of the engine (`uv pip install -e .` from the checkout,
the supported install for 0.2) or through `$GERMQUANT_CONFIG_ROOT`, as the container sets it. With the engine
installed any other way, `pack list` says so and names both fixes.

Install a pack someone sends you (a folder or a zip):

    germquant pack install budding_yeast-1.0.0.zip

Install checks every file against the pack's SHA256SUMS and refuses archives with paths that leave the pack
folder, symlinks and Windows junctions, non-ASCII recipe and channel-map files, a tunable that the recipe does not
set, and a copy of a built-in pack. Junk (`.DS_Store`, `__pycache__`, compiled Python, a zip's `__MACOSX`
folder) is never copied. An installed version never changes: the same version with other content is refused, and
so is a reinstall over an installed copy whose files no longer match its SHA256SUMS (`--replace` is for pack
developers, and reinstalls such a copy).

## Models

A pack names each model by name and digest, for example `cpsam_v2@0f1cc3f7` (the name, `@`, the first 8 hex
digits of the file's sha256), so the recipe text, and with it the config hash, changes whenever the weights
change. The form `<name>@<8 hex>` is reserved for packs: `germquant register-model` refuses such names. Release
models keep their names (germline_nuclei_combined).

Weights live once per digest in the model store, `models_dir()/blobs/<sha256>/<file>` (models_dir is
`models/models` in the checkout, or `$GERMQUANT_MODELS/models` when that is set, as on the HPC), with a
verification record beside them. `germquant pack fetch-models ID` tries, in order:

1. the file already in the store;
2. Cellpose's own cache (`$CELLPOSE_LOCAL_MODELS_PATH`, else `~/.cellpose/models`), after a size and sha256 check,
   hard-linked when the disk allows and copied otherwise; on the lab workstation this finds cpsam_v2 without a
   download;
3. a file you name with `--from-file NAME=PATH`, after the same check;
4. a download of the pinned address, which resumes after a broken connection and is checked by size and sha256.

A file that fails a check is never registered. Running fetch-models again does nothing. `pack check` trusts a
model whose size and modification time match its verification record and re-reads the whole file only when
either changed, so a GUI can check packs often.

When a run starts, a model named `<name>@<8 hex>` must resolve to weights whose sha256 starts with those 8 hex
digits and equals the pack's digest, wherever the file came from, or the run stops before any pixel is
processed; a name of that form that no installed pack declares is refused too, because 8 hex digits alone do not
pin a file. `germquant train-nuclei --base` applies the same check. A model the pack declares but the store lacks
stops the run with the fetch-models command to run. A recipe run with `segmentation.nuclei.method: classical`
never loads its model, so it needs none.

Third-party backends resolve their own model keys the same way with
`germquant.models_registry.resolve_model_name(name, expect_sha256=None)` (docs/BACKENDS.md).

## Recipes that build on a pack

A profile can extend a pack recipe by name, so a lab's edit stays small and the same on every computer:

    extends: pack:celegans_germline@0.1.0/rad51_foci
    spots:
      effect_size_min: 3.5

Always name the version (`@0.1.0`): an unversioned `pack:<id>/<recipe>` picks the newest version and is meant for
trying things interactively; `pack lint` refuses it inside a pack. `extends` resolves in this order: first the
places it always used (an absolute path; a relative path next to the file, then the repository's config folder;
a bare name in the repository's config folder, then next to the file), and only when none of them exists, a
`pack:` reference, or a bare name such as `generic.yaml` in `config/` and then `config/profiles/` of the engine
checkout (never of the working directory). An existing profile therefore resolves to exactly the file it always
did. The error message lists every place tried and the packs available.

A profile that builds on a pack this way records the paths it inherits (channel map, model, exclusions file)
with their original relative values and hashes those values in its stage hashes, so the same edited recipe gives
the same stage hashes on another computer and finished images stay finished there. Every other profile, the
shipped ones and study files with an absolute `extends` path included, hashes exactly as before.

## Making a pack

The layout (docs/DESIGN_system_packs.md, section 3.2):

    packs/<id>/
      pack.yaml                         the manifest (format 1, ASCII)
      SHA256SUMS                        the digest of every other file
      README.md                         the pack card, in plain language
      CITATION.cff  CHANGELOG.md
      config/profiles/<recipe>.yaml     recipe text (hashed)
      config/profiles/pack_base.yaml    the frozen base every outside recipe extends (hashed with it)
      config/channel_maps/<map>.yaml    channel map (hashed with the recipe)
      models/<model>.card.json          model cards; weights never live in a pack
      example/<name>.ome.tif            a tiny synthetic example, under 2 MB
      example/expected/recipe/<recipe>/ the recipe golden: nuclei.csv, spots.csv, image_summary.csv, config_hash.txt
      validation/<recipe>.json          validation records, written only by the benchmark tool

tests/fixtures/packs/test_pack is a complete small example; scripts/make_test_pack.py shows how its generated
files were made, and `germquant.packs.authoring` writes a frozen base from a profile. Recipes must sit at
`config/profiles/<recipe>.yaml` (the place that anchors their relative paths, such as the channel map, to the
pack folder), and every tunable a recipe lists must be a setting the recipe already has, with a single value:
the parameter form is built from the tunables (`packs.validate_pack` checks both, and install refuses a pack that
fails them).

Steps:

1. Write the recipes and the channel map. Outside recipes extend the pack's own `pack_base.yaml`, written once from
   a resolved profile, so a change to the engine's config.yaml never moves the pack's hashes.
2. Put each recipe's config hash into pack.yaml (`expected_config_hash`); `germquant pack lint` prints any that
   differ.
3. Write the recipe goldens on the synthetic example: `germquant pack golden ID --update`, then compare with
   `germquant pack golden ID`. The recipe tier runs with the classical segmenter and an empty model store.
4. `germquant pack build-sums packs/<id>` and `germquant pack lint packs/<id>`.

`pack lint` checks, for every pack: ASCII text only (no en or em dashes); no config/config.yaml inside the pack;
every expected hash current; SHA256SUMS current; versioned pack: lines; no absolute or machine path in a recipe,
map, card or pack.yaml; every card agrees with its model's digest and size. For outside recipes it also checks
that every worm stage switch is set explicitly, the file-name regex is the neutral `(?!x)x` (which matches
nothing) unless pack.yaml declares the pack's own `filename_regex`, `segmentation.nuclei.method` is cellpose or
classical (never auto, so a missing Cellpose fails loudly), the chain never reaches config/config.yaml, the batch
framing check is off (`qc.germline_outlier_factor` of at least 1e6), `qc.min_nuclei` is set, every spots instance
sets `expected_from_germ_cell: false`, and every key is one the engine reads: a key config/config.yaml sets or one
on the allow-list in `germquant.packs.lint`. Keys config.yaml carries but nothing reads, such as `spots.channel`
(the spots stage always counts the `foci` role), are refused.

Versions follow semver. Before 1.0.0, any change to hashed text or a model digest needs a new minor version;
from 1.0.0 on it needs a new major version. A new recipe or a changed engine range needs at least a minor
version, and a card, docs or metadata change a patch version. `germquant pack lint DIR --previous OLD.zip`
checks the bump. A recipe may say `status: validated` only with a passing validation record.

## Collecting outside results

The plain `germquant collect` keeps the worm's framing check (an image with 1.8 times the median germline count
is flagged). For an outside pack collect with `germquant batch FOLDER --pack ID --out OUT --resume`, which reads
the pack's `qc.germline_outlier_factor`.
