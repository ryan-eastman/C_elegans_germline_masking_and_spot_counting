# System packs

A system pack is a folder of checked files that you pick by name: the recipes (settings) and the pinned models
for one kind of sample. `germquant run IMAGE --pack celegans_germline` runs the C. elegans recipe; an outside lab
picks a pack made for its own cells. A pack holds no code and no model weights, only text, cards and a small
example, so it is safe to install from anyone.

Packs are new in germquant 0.2. The design is docs/DESIGN_system_packs.md (section 3).

## The packs that ship with the engine

| pack | version | what for | recipes |
|---|---|---|---|
| celegans_germline | 0.1.0 | the C. elegans germline recipes, exactly as the frozen engine config files define them | default, rad51_foci, segmentation_only, n2_sc, ccw77_partition, live_syp2 |
| cultured_cells | 0.1.0, provisional | nuclei of cultured cells on the DNA stain, gamma-H2AX and 53BP1 foci per nucleus in 3D | nuclei_foci (default, projected), nuclei_foci_3d |
| budding_yeast | 0.1.0, provisional | whole yeast cells on one transmitted-light plane, spots per cell in 3D | cells_spots (default, projected) |

The two outside packs are provisional: their settings are starting values (docs/DESIGN_system_packs.md, section
3.12), not yet tuned, and no recipe is validated (`status: provisional`, `validated: false`); each card
(packs/<id>/README.md) says what that means and what is provisional. Their recipes extend the pack's own frozen
`pack_base.yaml`, never config/config.yaml, so a change to the engine's config files never moves their config
hashes. The test pack used by the test suite lives in tests/fixtures/packs/test_pack.

The spot threshold is not on the same scale in every recipe. The yeast recipe sets `spots.sharpen`, as the
SpotMAX authors did, so its `effect_size` column and `spots.effect_size_min` (2.5) are on the sharpened
(difference of Gaussians) scale; the cultured cells recipes do not sharpen, so theirs (3.0) are on SpotMAX's
usual, unsharpened scale. A threshold never carries over from one scale to the other.

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
`python -m germquant.packs list --json` (also `show`, `check`, `fetch-models`, `fetch-example`, `install`,
`verify`, `lint` and `derive`) is the torch-free entry the napari plugin uses; `python -m germquant.packs
--version` prints the engine version, git sha and python, and `--version --json` adds the engine's features
(section "Versions and features" below).

`--project DIR` on `pack list`, `show`, `check` and `fetch-models` searches the project's own packs first
(`DIR/packs/`: the copy of the pack version the project was made with, and its derived packs), so a project's
packs list, check and fetch their models like any other:

    germquant pack list --project D:\projects\yeast_2026
    germquant pack check test_pack_derived_1a2b3c4d --project D:\projects\yeast_2026

## Example images

Every pack ships a tiny synthetic example (`example/<name>.ome.tif`, under 2 MB), which its recipe goldens run
on. A pack may also name a real public image, which is never copied into the pack: `germquant pack show ID`
lists the examples, and

    germquant pack fetch-example budding_yeast
    germquant pack fetch-example cultured_cells --out D:\germquant_examples

fetches one (the first real example, or the one `--example NAME` names), checks every downloaded file against the
size and sha256 the pack records, stacks the files into one OME-TIFF with the channel names and voxel size the
pack declares, and checks the result's pixels_sha256 (a digest of the pixels alone, because every OME-TIFF header
carries a fresh UUID). It prints the image path and the command to run it with the pack, and writes
`<example>.json` beside the image: where every byte came from, the data licence and the address of its text, what
was changed (which member or files, the crop, the channel names, the voxel size, the re-encoding as an OME-TIFF)
and the attribution to cite. Share that record with the image: CC BY 4.0 asks a shared copy to name the licence
and say what was changed. A file that fails a check is discarded; running the command again downloads nothing.
It also checks the assembled image's channel names and voxel size against the pack: when a later pack version
corrects them for the same pixels, the image is written again from its pixels (nothing is downloaded) and the
record is rewritten, as it is whenever it no longer matches the pack. The image goes to
`<out>/<pack>/<example>/`, where `--out` defaults to `$GERMQUANT_EXAMPLES`, else
`%LOCALAPPDATA%\germquant\examples` on Windows, `~/Library/Application Support/germquant/examples` on a Mac or
`~/.local/share/germquant/examples` on Linux. A folder inside a pack or a pack root (the engine's packs/ folder,
$GERMQUANT_PACKS, the per-user packs folder) is refused, because git tracks every file under packs/ and public
images are never redistributed. An unusable pack (one whose pack.yaml has problems) fetches nothing. From
Python: `germquant.packs.fetch_example(pack)` (torch-free) or `api.fetch_example(pack)`.

Some public datasets ship as one zip of many gigabytes. An example can name a single member of such a zip: the
fetcher then reads only the end of the archive, its table of contents and that member, with HTTP range requests
(`germquant.datafetch.fetch_zip_member`; ZIP64 archives included). It reads members stored without compression,
deflated or bzip2-compressed; a member compressed another way (deflate64, LZMA, zstd) or encrypted is refused with
a message naming the method, and so is a server that does not answer range requests. The archive is not trusted
to describe itself: a pinned archive size that the server contradicts is named (the deposit changed, so the
example must be pinned again), a table of contents above 64 MB is refused, a member whose recorded sizes cannot be
right is refused before any of its data is read, and a compressed member is inflated only up to the size the
archive records for it. No outside example declares a fallback: when a member cannot be read this way, download
the whole archive from its deposit and stack the member with `germquant assemble` (next section). The DeepFoci
zip stores every member without compression, so the fallback the design kept for it (one S-BIAD1494 stack, open
question O9) was not needed.

| pack | example | data | licence |
|---|---|---|---|
| budding_yeast | sbiad3750_sut509_p1 | BioImage Archive S-BIAD3750, SUT509 smFISH, Position_1: two 14.7 MB TIFFs (Bright, Cy3) | CC0 |
| cultured_cells | deepfoci_valid_258 | DeepFoci, Zenodo 4067741: one member (valid/data_258.tif, 145 MB) of the 10.2 GB validation zip, cropped to its two nuclei | CC BY 4.0 |

Each real example also gives the ranges a run of its recipe should land in (`expect` in pack.yaml: the number of
objects and the spots per object), and `scripts/make_pack_examples.py --run-real` fetches, runs and compares.

## Images with one file per channel

germquant reads one image per run and binds channels to roles by name. When each channel sits in its own TIFF (as
Cell-ACDC and many microscopes keep them), stack them first:

    germquant assemble --channels S_Bright.tif S_Cy3.tif --names Bright Cy3 --spacing 0.24 0.072 0.072 --out S.ome.tif

`--spacing` is the z step and the y and x pixel size in micrometres (a z-stack file often records no z step).
`--axes` declares how each input is laid out, outermost axis first, from the letters C, Z, Y and X: `ZYX` (the
default) is a z-stack of grey pages or a single plane, `ZCYX` a stack whose channels alternate page by page, and
`ZYXC` a stack whose pages hold the channels as samples (an RGB-like page per plane). Nothing is guessed: a file
that does not fit the declared axes is refused with the reason. A file whose own ImageJ or OME metadata records
both its channels and its planes must be declared in that order: `--axes CZYX` for a ZCYX hyperstack (or a
declaration without its channels or planes) is refused, because it would mix the channels up without a sign;
`--force` reads the pages as declared anyway. `--crop-z`, `--crop-y` and `--crop-x START STOP` keep a part of the
image. `--out` must be a new file, never one of the inputs. The command prints the pixels_sha256 of what it wrote.

## cpsam and cpsam_v2@0f1cc3f7

The outside packs name the Cellpose-SAM weights by digest, `cpsam_v2@0f1cc3f7`, while the generic profile names
the bare `cpsam`, which Cellpose 4.2 maps to its default weights, cpsam_v2 (docs/DESIGN_system_packs.md, section
15a). The pack name loads the verified file through `pretrained_model=`, the bare name loads Cellpose's default;
Cellpose builds the network from the weights, not the file name, so the two should give the same masks.
`scripts/compare_cpsam_names.py` checks it on the plane the yeast recipe segments (the sharpest transmitted plane
of the budding_yeast real example). It re-hashes the pack's blob against its pin, records whether the blob is the
same file as Cellpose's own copy (the file the bare name loads), and stops with exit code 2 when the device asked
for cannot run, so a CUDA check never passes on the CPU. On the lab workstation, `germquant pack fetch-models`
hard-linked the blob to ~/.cellpose/models/cpsam_v2, so both names load one file there, through two code paths
(the default load and `pretrained_model=`); the check is about those paths, not about two copies of the weights.
On 2026-10-11, on the CPU, the two gave identical masks: 0 of 490,000 pixels differ, 11 cells each (record
cpsam_names_cpu.json, written by the committed script). The same check on CUDA (`--device cuda`, run through the
GPU lock; its record must say "device": "cuda" and "cuda_available": true) is the acceptance check for the RTX 5090.
It ran on 2026-10-12 at 0521abf: identical masks again, 0 of 490,000 pixels differ, 11 cells each, "device": "cuda",
"cuda_available": true, on the NVIDIA GeForce RTX 5090 (record cpsam_names_cuda.json).

## Rows of the nuclei table

Table and column names never change, so in a pack whose objects are whole cells, one row of the "nuclei" table
is one cell. The pack card and `pack show` say which:

| object_kind | one row of the nuclei table is | n_nuclei counts | nucleus_id names |
|---|---|---|---|
| nucleus | a nucleus | nuclei | a nucleus |
| cell | a whole cell | cells | a cell |

## Where packs live

`germquant pack list` prints the search order. Packs are found, in this order:

1. the packs of a project (`<project>/packs/<id>/<version>/`: the copy of the pack version the project was made
   with, and the packs derived from it), while the project is open or running, and with `--project DIR`;
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

## Derived packs

A model trained in a project (docs/DESIGN_system_packs.md, section 10) can be promoted into a pack of its own,
so it appears in the pack list like any pack model and runs by name on another computer:

    germquant pack derive --from test_pack@0.1.0 --recipe with_model --model PROJECT/training/r01/models/r01_nuclei/r01_nuclei --card PROJECT/training/r01/models/r01_nuclei/r01_nuclei.card.json --out PROJECT/packs

That writes `PROJECT/packs/test_pack_derived_<sha8>/0.1.0/`, where `<sha8>` is the first 8 hex digits of the
model's sha256 (`--id NEW_ID` names it instead):

- `pack.yaml`: version 0.1.0, status provisional; the parent's system, object kind, engine range and expected
  sizes; `derived_from` (the parent, its recipe and its pack_sha256); one recipe; and the model spec
  `<name>@<sha8>` (the name is the card's, else the file's) with its sha256, size and card, and `sources: []`,
  because the weights are not in any download;
- `config/profiles/<recipe>.yaml`: only `extends: pack:<parent>@<version>/<recipe>` and
  `segmentation.nuclei.cellpose_model: <name>@<sha8>`. Both lines are hashed text without a machine path, so the
  recipe's config hash is new and the same on every computer. The recipe keeps the parent's needs (the new
  model in place of the one it replaces), tunables and dims;
- `models/<name>.card.json`: the card made portable. model_file becomes the file's name, base_model.path and
  every other absolute path become bare names (the sha256 values beside them stay), a path quoted in a text is
  shortened the same way, and non-ASCII text is written as JSON escapes. A `derived` record names the parent,
  the recipe, the sha256 of the card it came from and the keys it rewrote;
- `README.md` and `SHA256SUMS`.

A copy of the parent version is written beside it (`PROJECT/packs/test_pack/0.1.0/`), as a new project copies
it, so the derived recipe still resolves after an engine update moves the parent's built-in version. The model
goes into the model store, `models_dir()/blobs/<sha256>/<file>` (hard-linked when the disk allows), with its
verification record, where a run of the derived recipe finds it and checks its pin. Nothing is written before
every check has passed, the pack goes in through the same checks as `pack install`, and deriving again with
the same inputs changes nothing.

`--out` is a folder of packs: a project's `packs/` (searched while the project is open or running, and with
`--project PROJECT`), a `$GERMQUANT_PACKS` folder or the per-user folder (always searched). For any other folder
the command warns that nothing searches it. The engine's own packs/ folder and a folder inside a pack are
refused. So are, each with the reason and the fix: a parent recipe that never runs a Cellpose model
(`segmentation.nuclei.method: classical`, labels from files, or another segmentation backend), a parent recipe
that does not load or no longer hashes to its expected value, a card whose model_sha256 or model_size is not the
file's, an id that is not a pack id or is the parent's or a built-in pack's, and a folder that already holds the
id and version with other content (derive with another `--id`). A refusal exits with code 1.

`--json` prints the record: `ok`, `status` (derived, or already present), `pack`, `id`, `version`, `folder`,
`recipe`, `extends`, `config_hash`, `parent` {pack, folder, copied}, `model` {name, file_name, sha256, size, blob,
placed}, `card`, `rewritten` and `where` (searched, project, or not searched); a refusal prints
`{"ok": false, "error": ...}`. From Python: `germquant.packs.derive_pack(parent, recipe, model, card, out,
pack_id=None)`, which imports nothing heavy.

On another computer the weights are not in the pack. Copy the model file there and run

    germquant pack fetch-models test_pack_derived_<sha8> --project PROJECT --from-file r01_nuclei@<sha8>=PATH

`pack lint` on a derived pack reports what the outside-recipe rules ask of a recipe written by hand, because
those rules count only the keys set inside the pack folder: no recipe golden (the parent's recipe golden runs
the classical segmenter, which never loads this model), and the worm stage switches and the segmentation method,
which the parent's files set. A parent whose recipes are the engine's own files (celegans_germline, `root:
checkout`) adds the outside rules those frozen files were never meant to meet: the neutral file-name regex, the
framing check, `qc.min_nuclei`, `require_channel` and a chain that reaches config/config.yaml. It reports no
ASCII, dash, machine-path or unversioned `pack:` problem.

## Versions and features

`germquant --version` and `python -m germquant.packs --version` print the engine version, git sha and python.
`python -m germquant.packs --version --json` prints them as JSON, with the interpreter (`executable`) and
`features`: what this engine can do that its version string does not say, since the version stays 0.2.0.dev0
through tranche 5. The napari plugin gates its docks on these names.

| feature | what it is |
|---|---|
| packs | system packs (this file) |
| planar | single-plane images |
| labels_from | nucleus labels read from files (segmentation.nuclei.labels_from) |
| projects | project folders (germquant.project) |
| sample_sheet | the sample sheet (`--samples`) |
| import_labels | `germquant import-labels` |
| train_curve | `germquant train-curve` |
| pack_derive | `germquant pack derive`, and `--project DIR` on pack list, show, check and fetch-models |
| tune_spots | `germquant tune-spots` |

A feature is listed when its code is in the engine (`germquant.packs.FEATURE_PROBES`): its module file exists
and, for planar, projects and train_curve, whose module or package exists before the feature does, the module
defines `PLANAR_DZ`, `run_project` (in `germquant/project/run.py`) or `train_curve`. So projects is listed once
`germquant project run` exists, not while the package holds only its first stubs or commands. Nothing is
imported to find out, so the command stays torch-free, and an engine that lacks a feature never lists it.

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
files were made, and scripts/make_pack_examples.py does the same for the two outside packs (their frozen bases,
their synthetic examples, and `--draft-real`, which fetches a real example whose digests are not pinned yet and
prints the sizes, sha256 and pixels_sha256 to write into pack.yaml). `germquant.packs.authoring` writes a frozen
base from a profile, once: `make_pack_examples.py --base` leaves an existing base alone when the regenerated text
is the same and refuses, with the diff, when it differs (after a change to generic.yaml or config.yaml), because
a new base moves every recipe hash of the pack; `--rebase` replaces it, together with a minor version bump and a
CHANGELOG entry. A real example in pack.yaml lists `base_url` (https; plain http to this machine only while
`$GERMQUANT_ALLOW_LOOPBACK_EXAMPLES` is 1, for the stand-in servers of tests, and `pack lint` refuses it in any
case), `sources` (each a `path` under it with its `size` and `sha256`, plus `member`, `zip_size` and `crc32` for one
member of a zip; each saves under its file name, which must be a plain name, unique within the example),
`assemble` (`axes`, `channels`, `spacing_zyx_um` and an optional `crop`), `pixels_sha256`, `licence`,
`licence_spdx` (the SPDX id, such as CC0-1.0 or CC-BY-4.0), `licence_url` (the https address of the licence
text), `attribution` and optionally `recipe` and `expect`; the loader refuses a spec with any of them missing or
malformed. Write every digest and the CRC-32 in quotes: YAML reads an unquoted value of digits only as a number.
Recipes must sit at
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
map, card or pack.yaml; every card agrees with its model's digest and size; the pack's own example exists and
stays under 2 MB, no other file reaches 2 MB, and no TIFF sits in the pack unless it is a declared local example;
a real example's address and its licence_url are https and its `licence_spdx` is one of CC0-1.0, CC-BY-4.0,
CC-BY-3.0, BSD-3-Clause, MIT or Apache-2.0. For outside recipes it also checks
that every worm stage switch is set explicitly, the file-name regex is the neutral `(?!x)x` (which matches
nothing) unless pack.yaml declares the pack's own `filename_regex`, `segmentation.nuclei.method` is cellpose or
classical (never auto, so a missing Cellpose fails loudly), `segmentation.nuclei.require_channel` is true (unless
labels come from files), the chain never reaches config/config.yaml, the batch framing check is off
(`qc.germline_outlier_factor` of at least 1e6), `qc.min_nuclei` is set, every spots instance sets
`expected_from_germ_cell: false`, every recipe has its recipe golden (`example/expected/recipe/<recipe>/` with
nuclei.csv, spots.csv, image_summary.csv and a config_hash.txt equal to the recipe's expected_config_hash), and
every key is one the engine reads: a key config/config.yaml sets or one on the allow-list in
`germquant.packs.lint`. Keys config.yaml carries but nothing reads, such as `spots.channel` (the spots stage
always counts the `foci` role), are refused. `germquant pack golden ID` fails an outside recipe without its golden
and exits 1 when no golden ran at all.

Versions follow semver. Before 1.0.0, any change to hashed text or a model digest needs a new minor version;
from 1.0.0 on it needs a new major version. A new recipe or a changed engine range needs at least a minor
version, and a card, docs or metadata change a patch version. `germquant pack lint DIR --previous OLD.zip`
checks the bump. A recipe may say `status: validated` only with a passing validation record.

## Collecting outside results

The plain `germquant collect` keeps the worm's framing check (an image with 1.8 times the median germline count
is flagged). For an outside pack collect with `germquant batch FOLDER --pack ID --out OUT --resume`, which reads
the pack's `qc.germline_outlier_factor`. `germquant batch` runs the .nd2 files of a folder; folders of TIFF or
OME-TIFF images come with projects in a later germquant package, so run those one at a time with `germquant run`
for now.
