# Projects

A project is a folder that holds one study: the system pack and recipe it uses, the images, any changes to the
recipe, the models trained for it, and its results. The napari plugin creates and edits projects, and the
command line runs the same project unchanged, so a run started from the GUI and one started from a terminal or
a cluster write the same files.

Projects are new in germquant 0.2 (tranche 5, Train B). The design is docs/DESIGN_system_packs.md, section 7;
packs are described in docs/PACKS.md.

    germquant project new STUDY --pack cultured_cells --images "D:/scope/2026-10-20/*.ome.tif"
    germquant project edit STUDY --set spots.effect_size_min=4 --role foci="Cy3 RAD-51" --dry-run
    germquant project edit STUDY --set spots.effect_size_min=4 --role foci="Cy3 RAD-51"

## Layout and manifest

```
STUDY/
  germquant_project.yaml   the manifest (not hashed, ASCII)
  recipe.yaml              hashed; exists only after an edit
  channel_map.yaml         hashed with the recipe; exists only after a role edit
  samples.csv              optional design sheet
  runs/                    one log per run, runs/<run_id>.jsonl
  images/                  optional; images may live anywhere (absolute or UNC paths)
  results/<image_id>/      the results of each image
  annotations/             <image_id>__points_<annotator>.csv, truth points for spot tuning
  training/<round>/        training rounds: pairs, corrections.json, models with their cards
  tuning/<run>/            spot tuning runs: trials.jsonl, ranked.csv, best.yaml
  packs/                   the pinned pack version (copied at creation) and any derived packs
```

`germquant project new DIR --pack ID[@VER] [--recipe R]` copies the pack version into `DIR/packs/<id>/<version>/`
(through the same checks as `germquant pack install`), so the project keeps the exact version it was made with
when the engine moves on. Without `@VER` the newest version that fits this engine is used, and the project pins
it. `--images` adds images at once; `--name` names the project (default: the folder's name). A folder named
`config` cannot hold a project.

The manifest, germquant_project.yaml, is format 1. It records:

| key | what |
|---|---|
| name, id | the project's name and a uuid4 |
| created, saved | `{at, germquant, git_sha, plugin}`: when, and the engine version and git sha that wrote it, plus the plugin version when the plugin asked for the write (`--plugin-version`) |
| pack | `{id, version, recipe, sha256}`: the pinned pack version and the pack_sha256 of the project's copy |
| recipe | `{file, expected_config_hash, overrides, roles, origins}`: recipe.yaml or null; the config hash the project had when it was last written; the stages turned off (`{"spots.enabled": false}`); the channel role choices; where each recipe edit came from (form, tuning, model, or manual for a hand edit) |
| run | `{xy_stride, z_range, device, runner, interpreters}`: the run settings; device and runner null mean the engine's choice; interpreters `{recipe: {python, germquant, git_sha}}` |
| models | the models the project chose: `[{key, ref, sha256, card, card_sha256, source}]` |
| images | `[{id, path, size, mtime, timepoints, excluded, reason}]`; timepoints only for a time series |
| samples | the sample sheet's file name, or null |
| validated | `{value, why, record}`, computed at every write: false after any edit, or while the pack recipe is provisional |

The manifest is written as ASCII (a letter outside ASCII, such as the accented e of a name, is written as a YAML
escape, `Caf\xE9`) with LF line ends,
atomically: a temporary file next to it, flushed to disk, then renamed over the old file, retried for 2 seconds
when Windows reports the file busy (an indexer or a virus scanner). A failure at any point leaves the old file.
Only the commands that edit a project write it: new, add, edit, use-model and rounds --new. A run never writes it.

### Images

`germquant project add DIR PATH|FOLDER|GLOB ...` adds image files, folders and glob patterns:

- A folder is searched, with its subfolders, for the suffixes the recipe lists in io.input_suffixes, or for
  every suffix germquant reads when the recipe lists none (.nd2, .tif, .tiff, .ome.tif, .ome.tiff, .czi, .lif,
  .zarr, .ome.zarr). So a folder of OME-TIFFs can be added to a project on any pack, even one whose base names
  only .nd2 files.
- The search skips germquant's own outputs (`<id>__nuclei_labels.tif`, `<id>__spots.tif` and the other label
  and mask images) and the project's own results/, runs/, training/, tuning/, annotations/ and packs/ folders.
- Each file's header is read. A file is refused, with the reason, when it is missing, has an unsupported suffix,
  is a germquant output or lies in the project's own folders, has a header that cannot be read, holds more than
  one image (a multi-scene OME-TIFF; split it into one file per scene), is an .nd2 with more than one position
  (P > 1; export each position as its own .nd2), or has the id of an image already in the project. The other
  files are still added, and the command exits with code 1 when anything was refused.
- The image id is the file name without its format suffix (`gonad_01.ome.tif` is `gonad_01`) and names the
  results folder, so ids must be unique, and may not differ only in letter case (Windows and macOS would put
  both into one folder).
- An image inside the project folder is recorded by its relative path, so the project can be moved or zipped;
  any other image by its absolute path, a UNC path (`\\server\share\...`) as it is.

## Recipe edits and the hash

An unmodified project has no recipe.yaml. Its recipe is the pack recipe file itself, so its config_hash is the
hash of `germquant run --pack ID@VER --recipe R` and the pack's expected_config_hash, and it can be compared with
every lab that uses that pack version.

Any edit writes recipe.yaml:

```yaml
extends: pack:cultured_cells@0.1.0/nuclei_foci
spots:
  effect_size_min: 4.0
```

The extends line always names the pack version, so installing a newer version of the pack never changes the
project. recipe.yaml holds nothing else (no comment, no project name, no path), so the same edit of the same pack
version gives the same config_hash in every project and on every computer. An edit gives a new hash, and
validated becomes false: the pack's validation does not cover a changed recipe.

`germquant project edit DIR` makes the edits:

| option | what |
|---|---|
| `--set KEY=VALUE` | set a tunable of the pack recipe (the settings the pack declares for the parameter form; `germquant pack show ID` lists them). VALUE is read as YAML: `4` is a number, `true` a boolean, `max` a word |
| `--unset KEY` | drop an edit (or a stage switch, or `io.channel_map`, which drops every role edit) |
| `--role ROLE=CHANNEL` | bind a channel role to an image's exact channel name (see Channel roles); the roles not named keep their binding, and `ROLE=` drops that role's edit |
| `--revert` | drop every recipe edit, role edit and project model: back to the pack recipe |
| `--no-<stage>` | turn a stage off in this project (as `--no-spots` does on `germquant run`; it enters the hash) |
| `--device D`, `--xy-stride N`, `--z-range Z0 Z1`, `--runner R`, `--interpreter RECIPE=PYTHON`, `--samples FILE` | the run settings; the word `none` clears a device, z range, runner, interpreter or sample sheet |
| `--dry-run` | print what would be written (the files, the role preview, the new config_hash) and write nothing |
| `--plugin-version V`, `--json` | the plugin version for the stamps; machine-readable output |

What each kind of edit may set (the engine refuses anything else, naming the key and the rule):

- the parameter form (`--set`): only the tunables the pack recipe declares, with a single value of the declared
  type (a float, an int, true or false, or one of the declared choices) inside the declared minimum and
  maximum;
- spot tuning (`germquant tune-spots --apply`, docs/TUNING.md): spot keys (spots.*) from the tuning run's grid;
- a project model (`project use-model`): segmentation.nuclei.cellpose_model only;
- the role editor (`--role`): the channel map.

A value equal to the pack recipe's own is no edit, and when the last edit is dropped, recipe.yaml is removed and
the project hashes as the pack recipe again. A change of run settings alone (device, stride, interpreter)
never rewrites recipe.yaml. A hand edit of recipe.yaml is kept as an edit (origin manual), and the next edit
rewrites the file in the form above, from the version the manifest pins.

When a project is opened, its hash is computed again. A hash that differs from the one the manifest recorded is
a warning that names the change: for example recipe.yaml extending another pack version than the manifest pins,
or the project's copy of the pack replaced by other content.

## Channel roles

The pipeline finds its channels by role: the DNA stain (`dna`), the foci channel (`foci`) and the others of the
pack's channel map. The pack map names each role by tokens that match the channel names of the images it was made
for (DAPI, Hoechst); a lab whose channels are named after genes or dyes (`Cy3 RAD-51`, `MDN1`) binds the roles
with the role editor:

    germquant project edit STUDY --role dna="Hoechst 405" --role foci="Cy3 RAD-51" --dry-run --json

The edit writes channel_map.yaml: the pack's channel map with each edited role's names replaced by the exact
channel name chosen, its index fallback, required flag and marker kept, and recipe.yaml points at it. A role
matches when its name occurs in a channel's name, ignoring case; a role that matches no channel name falls back
to its index in the pack map, and the run flags it. A channel named with non-ASCII letters (a Greek gamma in
gammaH2AX) is written as a YAML escape and still matches.

`--dry-run` (and the API's preview_roles) shows, for every image of the project, how each role binds with the
edited map: by name, by index, or unbound. The preview is computed by the same code that writes the map and binds
the roles in a run, so it shows exactly what will run. An edit is refused when a required role binds no channel
on an image, or when two roles bind one channel.

The limit in this release: the channel map applies to every image of the project. Per-image bindings outside the
hash are not offered (design 7.3): the stage hashes that decide whether an image is done carry no channel roles,
so a binding kept outside the hash could leave results made with the old binding marked as done.

## Project models

A model trained in the project lives under `training/<round>/models/<name>/`. Rounds are made with

    germquant project rounds STUDY --new            # training/r01/, then r02 and so on
    germquant project rounds STUDY --new fine_tune  # a name of letters, digits, _ and -, up to 40 characters
    germquant project rounds STUDY                  # list them

and `germquant project use-model STUDY --model training/r01/models/m1/m1` puts the model into the recipe:
recipe.yaml names it by its path relative to the project (segmentation.nuclei.cellpose_model), and the manifest
records the file's sha256 and its card (`--card`, by default `<file>.card.json` next to the model, which
`germquant train-nuclei` writes). The model must be a file inside the project; a card that describes another file
is refused. The new hash makes validated false. docs/TRAINING.md describes the training loop.

A run compares the file it resolves with the sha256 in the manifest before it reads any pixel. That matters
because the engine looks a relative model path up under `$GERMQUANT_MODELS` first, so a shared model store holding
a file at the same relative path would otherwise win over the project's own file. To share a model with other
projects, promote it into a pack (`germquant pack derive`, docs/PACKS.md).

## The lock

Every command that edits a project takes the project lock while it writes: germquant_project.lock, a small JSON
file `{host, pid, time, command}` made with an exclusive create, which works the same on Windows, macOS and Linux
and on network shares. A second command that tries to write while the lock is held is refused, and its message
names the holder (computer, process id, time and command): try again when that command has finished.

A lock is removed by the command that made it. A lock left behind by a command that stopped (a crash, a closed
laptop) stays: a lock older than 12 hours is reported as stale, by the refusal and by `germquant project status`,
but it is never removed silently. Remove it, once that process is gone, with the holder's computer name and
process id as the message shows them:

    germquant project unlock STUDY --host LABPC --pid 4242

The command refuses a host or pid that is not the holder's, so it never removes a lock another command has taken
since.

## The light entry

`python -m germquant.project new|add|edit|use-model|rounds|unlock ...` takes the same arguments as `germquant
project ...` and never loads torch, Cellpose, SpotMAX, CuPy or the image readers: it reads .nd2 headers with the
nd2 package and TIFF and OME-TIFF headers with tifffile alone. It starts in a fraction of the time of the full
command line, and it keeps the cuBLAS load order of the RTX 5090 safe (torch must load before CuPy, or Cellpose
silently falls back to the watershed; docs/DESIGN_system_packs.md, rule R5). A CZI, LIF or OME-Zarr header is read
after the pipeline has loaded torch.

The napari plugin makes every project change through this entry, with the engine's own interpreter, so the code
that writes a project is the code that runs it. With `--json` each command prints one JSON object on stdout (an
error as `{"error": ...}`); the exit code is 0 when done, 1 when `add` (or `new --images`) refused some files,
and 2 when the command was refused before anything was written. `--plugin-version V` records the plugin version
in the manifest's stamps.

The Python API is the same code: germquant.project offers new_project, open_project, save_project, add_images,
project_config, preview_roles, write_recipe, revert_recipe, use_model, list_rounds, new_round and break_lock.
Each takes the project (a Project or its folder) explicitly, so it can be called from any thread, and importing
germquant.project loads nothing heavy either.
