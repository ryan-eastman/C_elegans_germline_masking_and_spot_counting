# Projects

A project is a folder that holds one study: the system pack and recipe it uses, the images, any changes to the
recipe, the models trained for it, and its results. The napari plugin creates and edits projects, and the
command line runs the same project unchanged, so a run started from the GUI and one started from a terminal or
a cluster write the same files.

Projects are new in germquant 0.2 (tranche 5, Train B). The design is docs/DESIGN_system_packs.md, section 7;
packs are described in docs/PACKS.md.

Windows (PowerShell, after `.venv\Scripts\activate`):

    germquant project new STUDY --pack cultured_cells --images "D:/scope/2026-10-20/*.ome.tif"
    germquant project edit STUDY --set spots.effect_size_min=4 --role foci="Cy3 RAD-51" --dry-run
    germquant project edit STUDY --set spots.effect_size_min=4 --role foci="Cy3 RAD-51"

Mac (the Terminal, Applications > Utilities > Terminal, after `source .venv/bin/activate`):

    germquant project new STUDY --pack cultured_cells --images "/Volumes/scope/2026-10-20/*.ome.tif"
    germquant project edit STUDY --set spots.effect_size_min=4 --role foci="Cy3 RAD-51" --dry-run
    germquant project edit STUDY --set spots.effect_size_min=4 --role foci="Cy3 RAD-51"

The commands are the same on both; only the paths differ. Keep the quotes around a glob pattern, so germquant
expands it rather than the shell, and around a channel name with a space in it.

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
  tuning/<run>/            spot tuning runs: grid.yaml, trials.jsonl, ranked.csv, best.yaml, run.json
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
Only the commands that edit a project write it: new, add, edit, use-model, rounds --new and tune-spots --apply
(docs/TUNING.md). A run never writes it.

### Images

`germquant project add DIR PATH|FOLDER|GLOB ...` adds image files, folders and glob patterns:

- A folder is searched, with its subfolders, for the suffixes the recipe lists in io.input_suffixes, or for
  every suffix germquant reads when the recipe lists none (.nd2, .tif, .tiff, .ome.tif, .ome.tiff, .czi, .lif,
  .zarr, .ome.zarr). So a folder of OME-TIFFs can be added to a project on any pack, even one whose base names
  only .nd2 files. The recipe's list may not be empty, and each suffix in it must be one of these, written as
  here (lower case, with the dot): every run of the recipe refuses anything else, so a recipe.yaml edited by hand
  to `['.TIF']`, `['tif']` or `['.png']` makes `project add` refuse (exit 2) before anything is added, naming the
  value and the suffixes allowed.
- The search skips germquant's own outputs (`<id>__nuclei_labels.tif`, `<id>__spots.tif` and the other label
  and mask images) and the project's own results/, runs/, training/, tuning/, annotations/ and packs/ folders.
- Each file's header is read. A file is refused, with the reason, when it is missing, has an unsupported suffix,
  is a germquant output or lies in the project's own folders, has a header that cannot be read, holds more than
  one image (a multi-scene OME-TIFF; split it into one file per scene), is an .nd2 with more than one position
  (P > 1; export each position as its own .nd2), or has the id of an image already in the project. The other
  files are still added, and the command exits with code 1 when anything was refused.
- On Windows, a path longer than 260 characters (a deep folder on a NAS) is refused too, naming its length,
  unless Windows long paths are turned on (the LongPathsEnabled setting): the image readers cannot open it.
  Move the folder to a shorter path, or map a drive letter to it. The same holds for a project model.
- The image id is the file name without its format suffix (`gonad_01.ome.tif` is `gonad_01`) and names the
  results folder, so ids must be unique, and may not differ only in letter case (Windows and macOS would put
  both into one folder). An id that no folder can carry on every system is refused too: one ending with a dot
  or a space, and a Windows device name such as CON or COM1.
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
| `--set KEY=VALUE` | set a tunable of the pack recipe (the settings the pack declares for the parameter form; `germquant pack show ID` lists them). VALUE is read as YAML 1.1, as the recipe file is: `4` is a number, `true` a boolean, `max` a word. Write an exponent with a dot and a sign, `1.0e-5` or `2.0e+3` (`1e-5` and `1.0e5` are read as words); `yes` and `no` are booleans, and `010` is the number 8 |
| `--unset KEY` | drop an edit (or a stage switch, or `io.channel_map`, which drops every role edit) |
| `--role ROLE=CHANNEL` | bind a channel role to an image's exact channel name (see Channel roles); the roles not named keep their binding, and `ROLE=` drops that role's edit |
| `--revert` | drop every recipe edit, role edit and project model: back to the pack recipe |
| `--no-<stage>` | turn a stage off in this project (as `--no-spots` does on `germquant run`; it enters the hash) |
| `--device D`, `--xy-stride N`, `--z-range Z0 Z1`, `--runner R`, `--interpreter RECIPE=PYTHON`, `--samples FILE` | the run settings; the word `none` clears a device, z range, runner, interpreter or sample sheet |
| `--dry-run` | print what would be written (the files, the role preview with its problems, the new config_hash) and write nothing; it exits 0 also when the edit itself would be refused, and says why |
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
rewrites the file in the form above, from the version the manifest pins. A typo that leaves recipe.yaml unreadable
as YAML makes every edit refuse with the line and column; fix the file, or go back with `--revert`.

Only a recipe.yaml that the manifest records is used, rewritten or removed. One put into the folder by hand, for
example copied from another project, is not used by runs: opening the project warns about it, no change removes
it, and an edit that would replace it is refused until it is moved out of the folder. Its settings can then be
made with `--set`.

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
the roles in a run, so it shows exactly what will run. An edit is refused (exit 2) while a required role binds
no channel on an image, or two roles bind one channel. The dry run lists those problems (`PROBLEM:` lines, and
`preview.problems` with `--json`) and still exits 0, so the plugin's role editor can show them in red before
anything is applied.

A lab that needs a map the role editor cannot write (several tokens per role, or `match_by: index`) can keep its
own map file in the project and name it in recipe.yaml by hand (`io:` then `channel_map: lab_map.yaml`, a path
inside the project). The project runs with that map, every other edit keeps it, and no change removes it. A role
edit is refused while recipe.yaml names it, since it would replace that map with the role editor's
channel_map.yaml: drop the line first with `--unset io.channel_map` to go back to the pack's map. Likewise a
channel_map.yaml the role editor did not write is never replaced.

The limit in this release: the channel map applies to every image of the project. Per-image bindings outside the
hash are not offered (design 7.3): the stage hashes that decide whether an image is done carry no channel roles,
so a binding kept outside the hash could leave results made with the old binding marked as done.

## Project models

A model trained in the project lives under `training/<round>/models/<name>/`. Rounds are made with

    germquant project rounds STUDY --new            # training/r01/, then r02 and so on
    germquant project rounds STUDY --new fine_tune  # a name of letters, digits, _ and -, up to 40 characters
    germquant project rounds STUDY                  # list them

A round name may not differ from another only in letter case (`R01` next to `r01` is one folder on Windows and
macOS).

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
since. A lock file that cannot be read (a command that stopped while it was making the lock) is reported as
such, whatever its age; once no germquant command is writing the project, delete it by hand.

## The light entry

`python -m germquant.project new|add|edit|use-model|rounds|unlock ...` takes the same arguments as `germquant
project ...` and never loads torch, Cellpose, SpotMAX, CuPy or the image readers: it reads .nd2 headers with the
nd2 package and TIFF and OME-TIFF headers with tifffile alone. It starts in a fraction of the time of the full
command line, and it keeps the cuBLAS load order of the RTX 5090 safe (torch must load before CuPy, or Cellpose
silently falls back to the watershed; docs/DESIGN_system_packs.md, rule R5). A CZI, LIF or OME-Zarr header is read
after the pipeline has loaded torch.

The napari plugin makes every project change through this entry (a tuning result through
`germquant tune-spots --apply`, docs/TUNING.md), with the engine's own interpreter, so the code that writes a
project is the code that runs it. With `--json` each command prints one JSON object on stdout (an
error as `{"error": ..., "kind": ...}`); the exit code is 0 when done, 1 when `add` (or `new --images`) refused
some files, and 2 when the command was refused. A refusal writes nothing. Exit 2 also reports the rare write that
stops part way, when Windows still denies a file after 2 seconds (an indexer or a virus scanner holding it):
recipe.yaml or channel_map.yaml may then be written while the manifest is not, and the message names what was
written. Run the same command again to finish it. `--plugin-version V` records the plugin version in the
manifest's stamps.

The Python API is the same code: germquant.project offers new_project, open_project, save_project, add_images,
project_config, preview_roles, write_recipe, revert_recipe, use_model, list_rounds, new_round and break_lock.
Each takes the project (a Project or its folder) explicitly, so it can be called from any thread, and importing
germquant.project loads nothing heavy either.

## Running a project

    germquant project run STUDY
    germquant project run STUDY --images gonad_01 gonad_02 --no-resume
    germquant project run STUDY --device cpu --json

`project run` runs every image of the project that is not excluded, in the order of the manifest, into
`results/<image_id>/`: the same files `germquant run --pack` writes for one image (the tables, the label images,
the montage, run_manifest.json and the completion marker `<image_id>__done.json`), made by the same code. The
recipe is recipe.yaml once the project was edited, else the pack recipe itself, with the project's stage
switches; an unmodified project therefore runs with the config_hash of `germquant run IMAGE --pack ID@VER
--recipe R`, and its numbers are those of that command.

- `--images ID ...` runs only those images.
- An image already done is skipped: its marker matches the recipe's config_hash, the stages, the xy stride and z
  range, and the segmentation model (a retrained model or a new Cellpose makes it run again). `--no-resume` runs
  every image again.
- `--device auto|cuda|mps|cpu` asks for that device on this computer; without it the project's run.device
  (`project edit --device`) is used, and without that the engine's choice (CUDA when it is available). The run
  sets GERMQUANT_DEVICE for its images and restores it afterwards. The xy stride and z range are the project's
  (`project edit --xy-stride`, `--z-range`).
- An image with timepoints (a time series) runs through the time-lapse loop of `germquant run --t all`
  (docs/TIMELAPSE.md), into `results/<image_id>/t<NNN>/`. A recipe that reads labels from files (labels_from,
  docs/LABELS.md) or runs single-plane images (io.planar true, or auto on an image with one plane, docs/PLANAR.md)
  cannot run a series in this release: such an image fails with that message, before anything of it is read, and
  the other images run.
- With a sample sheet (`project edit STUDY --samples samples.csv`, see Sample sheet below) each image runs with its
  row, as `batch --samples` runs it, and a resume reruns the images whose row changed. A time series that has a row
  fails before anything of it is read (the time-lapse loop gives its timepoints no row: take the row out, or clear
  the sheet with `--samples none`); a series without a row runs as before, its design columns from the file name
  and without the flag `metadata:not_in_sample_sheet`.
- A recipe that reads its nuclei from labels files does not load the segmentation model, so the model is not
  checked; the run says so.
- An image that fails, or cannot be reached (a NAS or a VPN that is not connected: "not reachable"), is named in
  the output and in the run's log, and the run goes on with the next image.

The run is refused, before any image is read and with nothing written (exit 2, the reason after `REFUSED:`, and
`{"error", "kind"}` with `--json`), when:

- the folder is no project, a named image is not in it or is excluded, or `--run-id`, `--shard` or `--device` is
  malformed;
- the recipe does not load, or the pack is not usable here (what `germquant pack check` reports: a missing backend
  or Python module, files that do not match the pack's SHA256SUMS);
- a model the run needs is missing, or cannot be resolved on this computer (a registered model whose file is
  not here, a name that is neither a file nor a Cellpose model), or the file the engine resolves is not the one
  the manifest records (its sha256 differs). This matters because the engine looks a relative model path up
  under `$GERMQUANT_MODELS` first: a shared model store holding a file at the same relative path would
  otherwise run instead of the project's own (see Project models). A pack model that is missing is fetched with
  `germquant pack fetch-models`;
- the sample sheet the project names cannot be opened or read (the message names the file, as `batch --samples`
  does);
- an edit holds the project's lock. A lock older than 12 hours is reported instead, and the run goes on (a run
  never removes a lock).

The exit code is 0 when every image was done or skipped, 1 when an image failed, and 2 for a refusal. A run
never writes the manifest and takes no lock; an edit made while a run goes on applies from the next run (the
run keeps the recipe it started with, and the images it finished then read stale). `--json` prints the run's
record on stdout, `{"run_id", "log", "counts": {"done", "failed", "skipped"}, "exit_code", "interrupted",
"warnings", "message"}`, and the progress lines on stderr. A run that stops part way for a reason that is no
image's own (its log can no longer be written, say) names the image in progress as failed in its log and ends
the log with interrupted true; it exits 1 with one ERROR line, and with `--json` prints `{"error", "kind":
"RunError", "run_id", "log"}`. Ctrl-C stops the run in the image in progress the same way. The plugin's Cancel
stops a run it started as a separate process at once, and a run inside napari before its next image. An image
stopped part way has no completion marker, so the next run does it again, and the images already done are
skipped.

The Python API is the same code: germquant.project.run_project(project, images=None, resume=True, device=None,
shard=None, progress=None, *, via="cli", run_id=None, should_stop=None), also in germquant.api, with
project_status, collect_project and export_commands. `progress(event)` gets every event of the run's log as it
is written, and `should_stop()` is asked before each image; a refusal raises ProjectError, and a run that stops
part way RunError (a RuntimeError naming its run_id and log).

## Shards and SLURM

    germquant project run STUDY --shard 1/3
    germquant project run STUDY --shard 2/3
    germquant project run STUDY --shard 3/3

`--shard K/N` runs the K-th of N slices of the images: every N-th image from the K-th on, in the manifest's
order (after `--images`, when given), so the N shards together run each image exactly once. They can run one
after the other, side by side on one computer, on several computers that share the project folder, or as the
tasks of a SLURM array:

    #SBATCH --array=1-3
    germquant project run /path/to/STUDY --shard ${SLURM_ARRAY_TASK_ID}/3

Each shard writes only its own run log and its own images' results folders, never the manifest, so shards never
write the same file; `germquant project status` merges their logs. Collect the tables (below) once every shard
has finished. `germquant project export-commands STUDY --shards 3` writes these lines for you.

Each shard takes its slice from the manifest when it starts, so add or exclude images only while no shard is
queued or running: a change in between can give an image to two shards, which may then run it at the same
time, or to none. `germquant project status` shows what is left to run.

## Snakemake

    snakemake -s workflow/Snakefile --config project=/path/to/STUDY --cores 1

The workflow's project mode runs one job per image of the project, each `germquant project run STUDY --images
<id> --via snakemake`, then `germquant project collect STUDY`. The jobs are named by the image ids, so
`gonad_01.ome.tif` is `gonad_01`. A job's output is the image's completion marker (the series record
`<id>__series.json` for an image with timepoints): an image whose marker is missing runs again, and after a
recipe edit every job runs again (the project's config_hash is one of the jobs' params). The sample sheet is not
one of them: after a row is edited, added or taken away, `project status` reads that image stale (sample_row) and
`germquant project run STUDY` reruns it; Snakemake does not. Snakemake deletes a
job's output before it runs the job, so every job it runs processes its image again, also an image that `project
status` reads as done, and `--forcerun project_image` runs every image of the project again. An image whose stage
failed keeps its marker, so Snakemake does not try it again; an image that fails ends its job with exit code 1,
and Snakemake then starts no further job unless `--keep-going` is given. To run only what is stale, missing or
failed, use `germquant project run STUDY` (or its shards).

Snakemake comes with the workflow extra (`uv pip install -e ".[workflow]"` from the repository folder), or from an
environment of its own. Its shell steps call the `germquant` command, which must be on PATH (the engine's
environment activated, or its Scripts or bin folder added to PATH). The default mode of the workflow
(`--configfile workflow/config.yaml`, one job per .nd2 of a folder) is unchanged.

## Run logs

Every run writes its own log, `runs/<run_id>.jsonl`. A run id is the run's start time in UTC and eight random hex
digits, `20261024T181503Z_0a1b2c3d`; `--run-id` gives one (the plugin does, to follow the run as it goes), and a
run id that is already used is refused, so no two runs ever write one file. `--via` records who started the run:
cli (the default), gui-subprocess, gui-inprocess or snakemake.

Each line of a log is one JSON object, in ASCII, written and flushed to disk at once, so a crash never leaves a
broken line behind a complete one: `{"event", "run_id", "time", ...}`, with the time in UTC to the millisecond.

| event | fields |
|---|---|
| run_start | via, host, os, interpreter, interpreter_source, engine {germquant, git_sha}, python, constraints, device_requested, device_selected, config_hash, config_hash_source, recipe {kind: project, pack, recipe, edited, file}, xy_stride, disabled_stages, argv, images [{path, image_id, out_dir}], plugin, log_location, project {root, name, id}, shard |
| image_done, image_failed, image_skipped | image_id, out_dir, status (done, failed, skipped), done_sha256, config_hash, stage_hashes, failed_stages, qc_pass, segment_device_used, elapsed_s, error |
| run_end | counts {done, failed, skipped}, interrupted |

These are the events and fields the napari plugin writes for its own runs, plus project and shard, and
interrupted (true when the run was stopped: between images, or part way through one, which is then logged as
failed). done_sha256 is the sha256 of the image's completion marker with Windows line ends (CRLF) read as LF, so
it is the same on every system. A log whose last line is not complete belongs to a run still writing it, or one
that was stopped: readers use its complete lines.

## Status

    germquant project status STUDY
    germquant project status STUDY --json

`project status` computes, for every image, what a run would do with it:

| status | meaning | reason_code |
|---|---|---|
| done | its results match the project's recipe and run settings; a run skips it | |
| stale | its results were made with another recipe or setting; a run makes them again | config_hash (the recipe changed: an edit, or another pack version), stages (a stage was switched on or off), geometry (the xy stride, z range or timepoint), planar_flag (a single-plane result checked without its single-plane segment hash; it runs again), segment_hash (the segmentation model, Cellpose or the labels file changed), sample_row (the image's row of the sample sheet changed, was added, or was taken away) |
| failed | a stage failed last time (failed_stages), or the last run of it failed (its error) | failed_stages, or none |
| missing | no results yet | |
| unreachable | not done, and its file cannot be reached: is the NAS or VPN connected? | |
| excluded | excluded from the project, with the reason | |

A time series that finished before the sample sheet gave it a row reads stale (sample_row), not done: every run
now fails it before anything of it is read (see Running a project), and its reason gives that refusal with the way
out (take the row out of the sheet, or clear the sheet with `--samples none`).

Each image also shows its qc_pass, the run that last processed it (last_run_id), the device its nuclei were
segmented on and the file a run reads (a relative path joined to the project folder). Then:

- the config_hash computed now, the one the manifest recorded, and whether the recipe is edited and validated;
- the pinned pack version, the versions of that pack installed on this computer, and whether the project's copy of
  it is in packs/;
- each model: the project's models (use-model) and the pinned models the pack recipe needs, whether it is present,
  verified (its sha256 is the recorded one) or mismatched (the file the engine resolves, and its digest), and
  whether this recipe needs it;
- the lock: who holds it, since when, and whether it is older than 12 hours;
- every run of the logs, merged by time, with its counts; a run whose log has no run_end (it was stopped without
  one, or is still running elsewhere) is shown as interrupted, with finished false;
- warnings: a config hash that changed since the project was saved (with the pack version change it comes from), the
  project's copy of the pack missing or replaced, a recipe.yaml the manifest does not record, a held or stale lock,
  a model a run would refuse (missing, not the expected file, or a segmentation model the engine cannot resolve
  on this computer, which status reports instead of failing; when the recipe segments with that model, the
  images' segment stage is then not compared), the sample sheet's own warnings, or a sheet a run would refuse
  (the images are then compared as without a sheet), and an environment mix: runs that used more than one device,
  interpreter or engine version (git sha), whose results can differ.

`--json` prints the same as one JSON object, the record germquant.project.project_status returns: format, root,
name, config_hash, expected_config_hash, edited, validated, pack, models, images (one object per image, with the
fields above), warnings, lock and runs.

## Collecting the tables

    germquant project collect STUDY

`project collect` stacks the images' tables into `results/batch_<table>.csv` and rebuilds
`results/batch_summary.csv`, as `germquant collect` does, but for the project's own images only (an image with
timepoints by its timepoint ids, `<id>_t<NNN>`) and with the recipe's QC settings: qc.germline_outlier_factor.
The outside packs set it to 1.0e+9, which switches off the framing check that flags a gonad with far more
germline nuclei than the others; on dense fields of cultured cells or yeast that check would flag nearly every
image. The plain `germquant collect results/` has no recipe and keeps the factor 1.8. With a sample sheet, its
own columns are joined into every table, as `germquant collect --samples` joins them. `--json` prints `{"root",
"results", "outlier_factor", "images", "counts", "summary", "sample_sheet"}` (`sample_sheet` null, or `{"path",
"columns"}` with the columns joined). A project without results is refused.

## Export commands

    germquant project export-commands STUDY
    germquant project export-commands STUDY --shards 4 --out run_study.sh
    germquant project export-commands STUDY --shell ps1 --out run_study.ps1

`project export-commands` prints, or writes to `--out`, the commands that run the project headless: a header
naming the project, its pack, recipe and config_hash, then `project run` (one line per shard with `--shards N`, and
for sh how to run the shards as a SLURM array) and `project collect`. The commands call the engine the project
records for its recipe (`project edit --interpreter RECIPE=PYTHON`), else the interpreter that wrote them, as
`python -m germquant`, so they run the same engine whatever comes first on PATH. `--shell sh` (the default) is for
Linux, macOS and SLURM; `--shell ps1` for Windows PowerShell (its file gets a byte-order mark, so PowerShell 5.1
reads a path with accents right). Everything else a run uses is in the project folder, so the lines name only the
folder.

## Portability

A project folder can be moved, copied or zipped to another computer:

- an image inside the project folder is recorded by its relative path and is found under the new folder; an image
  elsewhere keeps its absolute path, a UNC path (`\\server\share\...`) as it is. An image that cannot be reached
  there (a NAS without its VPN) is reported as unreachable by `project status`, and a run names it, "not reachable",
  and goes on with the others;
- the config_hash and every stage hash are the same in any folder and on any computer: recipe.yaml names the pack
  version, the project's copy of that version travels in packs/, and the hashes take the relative names of the
  pack's files, not where they are. So `project status` on the copy reads the finished images as done, and a run
  skips them;
- a run started on one computer and continued on another (or on a cluster) shows up in status as an environment
  mix: the results are valid, but those made on another device may differ in the last digits.

On Windows without long paths (the LongPathsEnabled setting), keep a project's folder path under about 170
characters: copying the pack into a new project (`project new`, `pack derive`) and the training loop's reads use
ordinary paths, which Windows stops at 260 characters, and the deepest files of a project lie well below its
folder. Images and models past the limit are refused by `project add` and `project use-model`, naming their
length.

## Sample sheet

A sample sheet gives each image its design columns when the file names do not carry them, or carry them wrongly. It
is a CSV file, samples.csv in a project, any file with `--samples`:

    image_id,genotype,sex,treatment,replicate,acquisition_date,slide
    20251105_n2_nohs_HERM_001,syp-2(ok307),male,control,3,2025-11-05,S12
    20251105_n2_nohs_HERM_002,N2,herm,control,4,2025-11-05,S12

- The first line names the columns. Each row names its image by `image_id`, the file name without its format
  suffix (the id that names the results folder: `x` for `x.ome.tif`), or by `file`: a file name, or a path relative
  to the sheet's folder (`/` or `\` between folders, read alike on Windows, Linux and macOS), which then names that
  one file.
- The design columns are `genotype`, `sex`, `germ_cell`, `treatment`, `replicate` and `acquisition_date`, written in
  lower case. A row's values go into every per-image table instead of what the file name gives through
  metadata.filename_regex, exactly as written (spaces around a value are dropped). A cell left empty keeps the file
  name's value. When a row gives `sex` but no `germ_cell`, germ_cell follows the row's sex as it does for a file name
  (herm: oocyte, male: spermatocyte), so the expected synaptonemal complex count follows the corrected sex. A sex
  other than herm or male (H, HERM, M or MALE in any letter case) on a row without germ_cell, and a germ_cell other
  than oocyte or spermatocyte, are kept as written, with a warning: such an image gets no expected synaptonemal
  complex or crossover count (`unknown` and `NA` say the value is not known, and are not warned about).
- Any other column is your own (a slide, a strain note, a scorer). No run reads it: `germquant collect --samples`
  joins it into every `batch_<table>.csv` by image id (empty for an image without a row). A column a table already
  has, such as `n_nuclei` in batch_image_summary.csv, is not joined there, with a warning; batch_summary.csv keeps
  its six columns. A column named like one of the metadata columns every table carries (`file_path`, `config_hash`
  and the others) is refused. Because the join goes by image id, two images in different folders that share an id
  (day1/img_001.nd2 and day2/img_001.nd2) both get the columns of a row that names one of them by its path, though
  the run gives that row's design columns to the named file only: give such images distinct names.
- An image the sheet has no row for keeps the file name's design columns and gets the flag
  `metadata:not_in_sample_sheet` (only in a run that was given a sheet). Rows for images that are not in the run are
  ignored.

Commands, Windows and Mac alike:

    germquant run IMAGE --config CONFIG --out RESULTS --samples samples.csv
    germquant batch FOLDER --config CONFIG --out RESULTS --samples samples.csv [--resume]
    germquant collect RESULTS --samples samples.csv

A plain `germquant batch` writes batch_summary.csv only, so the sheet's own columns reach the batch tables through
`germquant collect --samples` (or `batch --resume --samples`, which collects). In a project the sheet is the file the
manifest names (`germquant project edit STUDY --samples samples.csv`), and a project run applies it image by image as
`batch --samples` does.

**Resuming.** The sheet is not part of the config hash, like the file name it replaces, so adding a row for one
image never reruns another. Each image's row has a digest, row_sha256: the sha256 of the row's design cells written
as canonical JSON. A run with a sheet records it in the image's run geometry (`sample_sheet {source, row_sha256}` in
`<image_id>__stages.json`, and in run_manifest.json where a run writes one per image: `germquant run`, the API and
project runs) and in its completion marker (`sample_row_sha256`; null for an image without a row). `batch --resume` then reruns exactly the images whose row changed: a corrected design cell, a row
added for an image that ran without one, or a row taken away (resuming without the sheet reruns every image that
ran with a row). Editing your own columns reruns nothing. The reason is printed, `reprocessing: sample row differs
(...)`.

**Saving from Excel.** Save the sheet as "CSV UTF-8 (Comma delimited)". The digest depends only on the values, so
saving the same values again, with Windows line ends, a byte-order mark, another column order or quotes around every
cell, reruns nothing. Excel's plain "CSV" on Windows is cp1252: it is read, with a warning, and gives the same
values. A sheet with semicolons between the columns (what Excel saves in some languages) or saved as "Unicode Text"
is refused with a message. Keep dates in ISO form, 2025-11-05 (or 20251105): Excel rewrites a date in the computer's
own format when it saves, so format the column as Text before you type them. A date in another form is kept as it
is, with a warning.

**Time series.** A sheet does not apply to a time series in this release: `germquant run SERIES --t ... --samples
FILE` is refused (exit 2, before any timepoint), because the time-lapse loop gives its timepoints no sample row. Run
the series without the sheet; its design columns then come from the file name.

**From code.** germquant.io.sample_sheet reads and writes sheets without loading anything heavy:
`read_sample_sheet(path)` returns a SampleSheet, whose `row(image_id, path)` is the SampleRow that
`api.run(..., sample=row)` and `process_image(..., sample=row)` take; `write_sample_sheet(path, rows)` writes a list
of dicts as ASCII with LF line ends, atomically, and refuses a value outside ASCII, naming its row and column.
