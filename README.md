# germquant: nucleus segmentation and RAD-51 spot counting in C. elegans germlines

This program takes a confocal image of a worm germline (a `.nd2` file), finds every nucleus in 3D,
decides which nuclei belong to the gonad, and counts the RAD-51 spots inside each one. The answers
come out as spreadsheet files (`.csv`) that open in Excel, plus a picture you can glance at to
confirm it worked and 3D mask files you can load into Imaris if you want to check the segmentation.

If the image also has a PGL-1 channel, the same run measures how much SYP signal sits in the
P granules around each nucleus (the colocalization stage). The full SYP-3 / PGL-1 heat shock
analysis for the ccw77 strain lives in `analysis/coloc/` and has its own README.

You do not need to know how to code to use the basic version. This guide assumes you have never
opened a terminal. Unfamiliar words are explained in "Words you might not know" near the end.

## Which computer, and how the commands are written

The program runs on two kinds of machine:

- The lab workstation: Windows, with an NVIDIA graphics card. Already set up. Fastest.
- An Apple Silicon Mac (M1, M2, M3 or M4 chip): uses the Mac's own graphics chip for the slow part.
  Setting one up takes roughly twenty minutes plus download time (see "Setting up a Mac" near the
  end). It is slower than the workstation, so use it for a few images at a time and send big
  batches to the workstation or the cluster.

The Mac version has been checked with automated tests but has not yet been run on a real Mac. The
Mac times on this page are estimates, not measurements. If you are the first person to run it on a
Mac, run the `check-gpu` command from the set-up section, then one image, and tell Ryan how long it
took and whether anything on this page was wrong.

Every command in this guide is shown twice, once for each machine. They differ only in how the
program is called and how paths are written:

| | Windows (workstation) | Mac |
|---|---|---|
| Where you type | PowerShell (a blue window) | Terminal (a white window) |
| The program | `.venv\Scripts\germquant.exe` | `.venv/bin/germquant` |
| Paths | `C:\Users\you\images\worm1.nd2` | `/Users/you/images/worm1.nd2` |
| Point-and-click file | `quantify.bat` (drag the image onto it) | `quantify.command` (double-click it, then drag the image into the window) |
| The key that runs a command | Enter | Return (the same key) |

Putting a path into a command: delete the placeholder, including its quotes, then drag the file or
folder from File Explorer or Finder into the window. Windows pastes the path with quotes around it
when it contains a space; keep them. A Mac pastes the path with a backslash before each space and
no quotes; leave it exactly as pasted and do not add quotes.

Wherever this page says `you` in a path, it means your own user name on that computer.

## The short version: one image, drag and drop

Copy the image to your own computer first (the Desktop or an images folder). The results folder is
created next to the image, which does not work on the lab NAS share for most accounts and is slow
over the network.

Windows:

1. In the program's folder, find the file `quantify.bat`.
2. Drag your `.nd2` image onto `quantify.bat` and let go.
3. A black window opens and starts working. Leave it open. One image takes about 15 to 25 minutes.
4. When it prints `*** DONE! Your results are in: ***` followed by the results folder and then
   `Press any key to continue`, it is finished. Press any key and the window closes. The results are
   in a new folder next to your image, named `yourimagename_results`.

Mac:

1. In the program's folder, find the file `quantify.command`. The very first time, macOS refuses to
   open it because the file came from the internet. On macOS 15 or newer: double-click it, click
   Done on the warning, open System Settings, Privacy & Security, scroll to the bottom, click Open
   Anyway next to the message about quantify.command, enter your password, then double-click the
   file again and click Open. On macOS 14: right-click the file, choose Open, then click Open. You
   only do this once.
2. Double-click `quantify.command`. A Terminal window opens and asks for the image. Drag your `.nd2`
   file from Finder into that window (its path appears), then press Return. (Dropping the image
   onto the file's icon does nothing on a Mac; the file has to be opened first.)
3. The window starts working. Leave it open. We expect roughly 30 to 60 minutes, but this has not
   been timed yet.
4. When it prints `*** DONE! Your results are in: ***` followed by the results folder and then
   `Press Return to finish`, it is finished. Press Return; Terminal then shows `[Process completed]`
   and leaves the window open, which is normal. Close it with Command-W. The results are in a new
   folder next to your image, named `yourimagename_results`.

That is the whole procedure for a single RAD-51 image. The rest of this page explains the results,
how to run many images, how to run other kinds of images, and what to do when something goes wrong.

## What you need

- A set-up computer: the lab workstation, or a Mac set up as described at the end of this page.
- A copy of the per-gonad `.nd2` file on that computer. The whole-slide overview files (names
  containing `10x` or `largeimage`) will not work; use the individual gonad files.
- Time. About 20 minutes per image on the workstation, longer on a Mac. The computer does the
  work; you wait. Do not start Imaris or another heavy program on the same machine while it runs
  (they compete for the graphics memory).

## Running from a terminal

Use this when the point-and-click file does not work, when you want to choose where the results
go, when your image is not the standard DAPI / SYP / RAD-51 panel, or when you want to process a
whole folder.

1. Open the terminal.
   - Windows: click the Start menu, type `PowerShell`, click "Windows PowerShell". A blue window
     opens.
   - Mac: press Command and Space, type `Terminal`, press Return. A white window opens.
2. Go to the program's folder: type `cd`, then a space, then drag the program's folder from File
   Explorer (or Finder) into the window, then press Enter. The line will look like one of these:

   Windows:
   ```
   cd "C:\Users\you\C_elegans_germline_masking_and_spot_counting"
   ```
   Mac:
   ```
   cd /Users/you/C_elegans_germline_masking_and_spot_counting
   ```
   On the workstation the folder is `C:\Users\ryane\C_elegans_germline_masking_and_spot_counting`.
   If the code was downloaded as a ZIP, the folder name ends in `-main` unless someone renamed it.
3. Run one image. Copy the line for your machine, replace the two placeholder paths with yours
   (delete the placeholder and drag the file or folder in, as described at the top), press Enter.

   Windows:
   ```
   .venv\Scripts\germquant.exe run "C:\path\to\YOUR_IMAGE.nd2" --config config\config.yaml --out "C:\path\to\RESULTS_FOLDER"
   ```
   Mac:
   ```
   .venv/bin/germquant run /path/to/YOUR_IMAGE.nd2 --config config/config.yaml --out /path/to/RESULTS_FOLDER
   ```
   The first path is your image. The path after `--out` is any folder where you want the results
   saved (it is created if it does not exist).
4. Wait. When the blinking prompt comes back, it is finished.

### A whole folder at once

Windows:
```
.venv\Scripts\germquant.exe batch "C:\path\to\FOLDER_OF_IMAGES" --config config\config.yaml --out "C:\path\to\RESULTS_FOLDER"
```
Mac:
```
.venv/bin/germquant batch /path/to/FOLDER_OF_IMAGES --config config/config.yaml --out /path/to/RESULTS_FOLDER
```

This processes every `.nd2` under the folder, including sub-folders, one image at a time, and
skips the `10x` and `largeimage` overviews on its own. It writes one results sub-folder per image
and a `batch_summary.csv` listing them all with their nucleus counts and QC flags. Many images can
take hours; leave the window open. If it is interrupted, add `--resume` to the same line: images
whose results folder already holds a completion marker written with the same config, the same
stages, the same xy stride and the same nucleus model, and with no failed stage, are skipped; the
rest are processed, and the summary is rebuilt from every finished folder of the study. Without
`--resume` the batch reprocesses everything.

To drop particular gonads from a study (a fused carcass, a second germline limb), list their ids in
a JSON file and point `qc.exclusions_file` in the config at it. An entry matches whole name parts
only: `HS_male_07` drops `..._HS_male_07` but not `..._noHS_male_07` or `..._HS_male_070`. The
analysis-style file with `excluded_short_ids` and `excluded_batch` is accepted as is. Both `batch`
and the Snakemake workflow read it, and the excluded files are named in the run manifest.

### Stacking the results of a batch

Windows:
```
.venv\Scripts\germquant.exe collect "C:\path\to\RESULTS_FOLDER"
```
Mac:
```
.venv/bin/germquant collect /path/to/RESULTS_FOLDER
```

writes one `batch_<table>.csv` per table next to `batch_summary.csv`: `batch_nuclei.csv`,
`batch_spots.csv`, `batch_granules.csv`, `batch_coloc.csv`, `batch_image_summary.csv`, plus the
tables of the optional stages (`batch_zones.csv`, `batch_partition.csv`, `batch_sc_tracks.csv`,
`batch_sc_per_nucleus.csv`, `batch_mask_audit.csv`, `batch_granule_tail.csv`). Tables for stages
you did not run are still written, with a header row and no data. Run it whenever you like, for
example after reprocessing a few images.

### Staging pachytene by hand

Automatic pachytene staging was never reliable enough on these gonads, so it is drawn by hand once
per image and reused by every stage that needs it. `trace`, `restage` and the run or batch must all
be given the same `--profile` (or `--config`), because that file decides which traces file is
written and read; the `ccw77_partition` profile turns staging on and names the file.

Windows:
```
.venv\Scripts\germquant.exe trace "C:\path\to\RESULTS_FOLDER" --profile ccw77_partition
```
Mac:
```
.venv/bin/germquant trace /path/to/RESULTS_FOLDER --profile ccw77_partition
```

opens each finished image (DAPI grey, SYP red; scroll for single planes) and you click a line from
the pachytene start to its end, then press Enter. The line is saved in whole-image microns to the
traces file named by `staging.traces_file` in the profile. With staging on, the pipeline then
assigns every germline nucleus a zone (early, mid, late thirds of your line, or pre, post,
off_axis). After you redraw a line, recompute the zones without reprocessing the image:

Windows:
```
.venv\Scripts\germquant.exe restage "C:\path\to\RESULTS_FOLDER" --profile ccw77_partition
```
Mac:
```
.venv/bin/germquant restage /path/to/RESULTS_FOLDER --profile ccw77_partition
```

`trace` also takes optional image ids after the folder, `--traces PATH` (a different traces file,
relative to where you are), `--redo` (revisit traced images) and `--stride N` (display resolution,
default 2).

### Counting a second kind of focus (COSA-1 crossovers)

The spot counter is not tied to RAD-51. A profile can declare further instances of it, each on its
own channel role with its own table and nucleus column, for example COSA-1 crossover foci:

```
spots:
  instances:
    - name: crossover           # -> table spots_crossover, nuclei column n_spots_crossover
      role: crossover_foci      # a role you add to a channel map (do not edit the shipped maps)
      restrict_to_zone: late    # late-pachytene mean when the staging stage has run
      expected_from_germ_cell: true   # 6 per oocyte, 5 per spermatocyte, reported next to the mean
```

Detection parameters default to the `spots` block and can be overridden per instance.

### Choosing what to measure (profiles)

Instead of `--config config/config.yaml` you can name a profile with `--profile`:

| Profile | What it runs |
|---|---|
| `rad51_foci` | nuclei, germline, gonad axis, RAD-51 foci (no PGL-1 or colocalization) |
| `segmentation_only` | nuclei, germline and gonad axis only, no spots or colocalization (for Imaris surfaces) |
| `n2_sc` | adds the SC tracer for fragmentation (see `docs/SC_TRACING.md`) |
| `ccw77_partition` | the SYP-3 / PGL-1 partition analysis: lamin envelopes, hand-traced staging, the Imaris-calibrated granule recipe, on top of `config/config_ccw77.yaml` |
| `generic` | samples that are not a worm germline (cultured cells, yeast, sections): nuclei with the stock Cellpose model and puncta per nucleus, every worm stage off (see `docs/GENERIC_PROFILE.md`) |

Windows:
```
.venv\Scripts\germquant.exe run "C:\path\to\YOUR_IMAGE.nd2" --profile rad51_foci --out "C:\path\to\RESULTS_FOLDER"
```
Mac:
```
.venv/bin/germquant run /path/to/YOUR_IMAGE.nd2 --profile rad51_foci --out /path/to/RESULTS_FOLDER
```

A profile is a short file in `config/profiles/` that switches stages on or off on top of a base
config. Every optional stage can also be switched off on the command line with `--no-<stage>`
(`--no-spots`, `--no-coloc`, `--no-granule`, `--no-axis`, `--no-germline`, `--no-envelope`,
`--no-staging`, `--no-partition`, and so on). The staging stage has to run before the partition
zones or a late-pachytene spot mean can exist. The `alpine` and `slurm` folders in
`config/profiles/` are cluster job settings, not names for `--profile`.

### Samples that are not a worm germline

`--profile generic` runs the nucleus segmentation (stock Cellpose-SAM), the per-nucleus
measurements and the SpotMAX puncta count on any 3D stack with a nuclear stain and an optional
puncta channel, with every *C. elegans* stage (germline isolation, gonad axis, staging, envelope,
SC tracing, partitioning) switched off; `in_germline` then reads True for every nucleus. Its channel
map, `config/channel_maps/generic_nuclei_puncta.yaml`, resolves channels by descriptive name or by
position. `docs/GENERIC_PROFILE.md` says what each output column means outside the worm, how to
adapt the recipe in a study profile (your channel names, spot radius, filename pattern), and why
single-plane (2D) images are not supported yet.

### Segmentation only (no spot counting)

Add `--no-spots` to the run or batch command. It skips the RAD-51 spot detection, so it is much
faster and cannot get stuck on a difficult image. You still get the nucleus table, the 3D nucleus
label image for Imaris, and the montage; the spots table is empty and no spots image is written.

Windows:
```
.venv\Scripts\germquant.exe run "C:\path\to\YOUR_IMAGE.nd2" --config config\config.yaml --out "C:\path\to\RESULTS_FOLDER" --no-spots
```
Mac:
```
.venv/bin/germquant run /path/to/YOUR_IMAGE.nd2 --config config/config.yaml --out /path/to/RESULTS_FOLDER --no-spots
```

`--no-coloc` skips the PGL-1 / SYP colocalization stage in the same way.

### Checking what is in an image

Windows:
```
.venv\Scripts\germquant.exe info "C:\path\to\YOUR_IMAGE.nd2"
```
Mac:
```
.venv/bin/germquant info /path/to/YOUR_IMAGE.nd2
```

prints the channel names and the voxel size without loading the pixels. Use it whenever you are not
sure which config to pick (next section).

### Checking the graphics chip is being used

Windows:
```
.venv\Scripts\germquant.exe check-gpu
```
Mac:
```
.venv/bin/germquant check-gpu
```

On the workstation it names the NVIDIA card. On a Mac the output should contain the words "Apple
Silicon GPU via Metal" and a line starting with "OK". If it says "No accelerator", the nucleus
segmentation will run on the processor and take hours; see "If something breaks".

## Which config file to use

The config file tells the program which channel is which and holds every tunable setting. Pick the
one that matches your image; everything else stays the same. (The files are written with forward
slashes here; on Windows `config/config.yaml` and `config\config.yaml` are the same file.)

| Config | Image type | What runs |
|---|---|---|
| `config/config.yaml` | N2-style DAPI / SYP / RAD-51, with or without a fourth PGL-1 channel | nuclei, germline isolation, RAD-51 spots; colocalization only if a PGL-1 channel is present |
| `config/config_ccw77.yaml` | ccw77 four-colour IF: DAPI / PGL-1::GFP (477) / SYP-3::mCherry (545) / LMN-1 (640) | nuclei, germline isolation, SYP-3 / PGL-1 colocalization with the lamin envelope; there is no RAD-51 channel so spot counting skips itself |
| `config/config_n2dryice.yaml` | the N2 dry-ice test set | as `config.yaml` |

The channel assignments themselves are in `config/channel_maps/`. Channels are matched by name
(the laser line, for example `405`, `477`, `545`, `640`) with a fixed-position fallback. If the
`info` command shows names that are not in the map, or the montage shows the wrong channel in the
wrong place, the map needs editing; ask whoever maintains the tool rather than guessing. Which
wavelength carries which protein changes between experiments (in the N2 panels 477 is SYP; in ccw77
477 is PGL-1 and 545 is SYP-3), so confirm it by eye on a new experiment before trusting any table.

## Your results: where they are and what they mean

Inside the results folder, every file starts with the image name followed by two underscores.

| File | What it is |
|---|---|
| `..._nuclei.csv` | The main result. One row per nucleus. `n_spots` is the RAD-51 count in that nucleus. `in_germline` is True for gonad nuclei and False for gut, debris and other tissue; ignore the False rows. `axis_position_norm` is the position along the gonad (0 = distal tip, 1 = proximal end) when the automatic axis fit succeeded. With a PGL-1 channel there are also `n_granules` and `granule_volume_um3` per nucleus. |
| `..._spots.csv` | One row per RAD-51 spot: 3D position, which nucleus it is in, brightness and effect size. For deeper analysis. |
| `..._granules.csv` | One row per PGL-1 granule (only with a PGL-1 channel). |
| `..._coloc.csv` | The SYP / PGL-1 colocalization metrics, one row per operand (only with a PGL-1 channel). See `docs/COLOCALIZATION.md` for what each column means. |
| `..._image_summary.csv` | One row of totals for the image: nucleus counts, mean spots per nucleus, the headline colocalization numbers, QC pass or fail and the flags behind it, and `segmentation_method`, which must say `cellpose` (see "If something breaks"). |
| `..._montage.png` | The picture to look at every time: each channel, the nucleus outlines and the detected spots side by side. If this looks wrong, the numbers are wrong. |
| `..._nuclei_labels.tif` | The 3D nucleus map, with the voxel size stored in the file. Load it into Imaris as Surfaces to check the segmentation. |
| `..._spots.tif` | The detected spots as a 3D image on the same grid. Load it into Imaris as a channel (Edit, Add Channels) or run Imaris Spots on it (diameter about 0.4 um) to get Spots objects next to your Surfaces. |
| `..._stages.json`, `..._done.json` | Which stages ran, were skipped or failed, with timings; and the marker that says the run finished (batch `--resume` looks for it). |
| `.parquet` copies | The same tables in a compact format for R or Python. |

Optional stages add their own tables next to these (`..._sc_per_nucleus.csv`, `..._zones.csv`,
`..._partition.csv`, `..._granule_tail.csv`, `..._mask_audit.csv`); they are written empty when the
stage is off, so a spreadsheet or script never finds a file missing.

To get the average number of RAD-51 spots per nucleus: open `..._nuclei.csv` in Excel or Numbers,
filter to `in_germline = True`, and average the `n_spots` column.

Every table also carries the image metadata parsed from the file name (date, genotype, treatment,
sex, replicate), the voxel size, the channel map that was used, and the version of the code, so a
result can always be traced back to how it was made.

### Before you trust the numbers

The spot counts have been checked against Imaris counts for N2 worms, with and without heat shock,
and they agree to about one spot per nucleus across 4 to 21 foci per nucleus. They have not been
checked for other genotypes. Before reporting absolute counts on a mutant, count a few gonads in
Imaris and compare; the detection settings were tuned on N2 and may need re-tuning. The details and
the caveats are in `docs/SPOTMAX_VALIDATION.md`.

The `in_germline` decision uses SYP signal and spatial connectivity, never the spot count; on 14
validation gonads the nuclei it kept held about 97 percent of the RAD-51 spots. It can admit sperm
masses and somatic nuclei that touch the gonad; for the ccw77 colocalization analysis that mattered,
and the audit and the fix are described in `analysis/coloc/README.md`.

Numbers from a Mac and from the workstation: the nucleus segmentation runs with slightly different
arithmetic on the two machines (the Mac uses 32-bit numbers on its graphics chip, the workstation
16-bit), so masks can differ a little. Before pooling results from both, run the same gonad on
each and compare the tables (`scripts/regression_diff.py` does this cell for cell). Everything after
segmentation is identical arithmetic on both.

## The colocalization stage

When the channel map resolves a PGL-1 channel, the run also segments the P granules, builds a thin
cytoplasmic shell around each germline nucleus (anchored to the lamin envelope when a lamin channel
is present), and measures how much SYP signal coincides with the granules inside that shell. The
headline columns are `shell_pearson`, `shell_manders_m1` and `shell_manders_m2` in
`..._image_summary.csv`; the per-operand detail is in `..._coloc.csv`. The reasoning behind the
region choice and the meaning of every column are in `docs/COLOCALIZATION.md`.

The stage-resolved analysis of SYP-3 partitioning into P granules under heat shock (partition
coefficient, hand-traced pachytene zones, controls, figures) lives in `analysis/coloc/`, and the
same analysis now runs inside the pipeline as the `ccw77_partition` profile.

## If something breaks

On either machine:

- A red error saying the `.venv` folder is missing: the program is not installed on this
  computer. See the set-up section for your machine.
- "file not found": the image path is wrong. Delete the placeholder and drag the `.nd2` into the
  window to paste it correctly (Windows adds quotes for you; a Mac adds backslashes; leave either
  as it is).
- "no such file or directory" (Mac) or "The term '.venv\Scripts\germquant.exe' is not recognized"
  (Windows) when typing a command: the window is not in the program's folder. Run the `cd` line
  from step 2 of "Running from a terminal" again.
- "Permission denied", "Read-only file system" or "Access is denied" right after the Results line:
  the image is on a folder you cannot write to (the NAS share). Copy the image to your own computer
  and run it from there, or use the terminal route with `--out` pointing at a folder of yours.
- It said DONE but the montage shows blocky, merged, fragmented or missing nuclei: the nucleus
  step did not run the trained model. Open `..._image_summary.csv` and look at the
  `segmentation_method` column; it must say `cellpose`. If it says `classical`, the trained model
  or the graphics set-up is missing and a rough backup method was used; run `check-gpu`, then the
  `get-model` command from the set-up section, and run the image again.
- It found 0 nuclei, or the montage shows the wrong channels: the image's channel names do not
  match the config. Run the `info` command on the file and compare with `config/channel_maps/`.
  Tell whoever maintains the tool.
- The window closed instantly: run the command from the terminal so the error stays on screen.
- The spot step runs for hours on one image: that image is probably flooded with signal or an
  artefact. The run caps the number of candidate spots so it should still finish; if you only need
  the segmentation, rerun with `--no-spots`.
- Anything else: copy the red text and send it to whoever maintains the tool.

Windows:

- A blue box saying "Windows protected your PC" when you double-click `quantify.bat` (this happens
  when the code was downloaded as a ZIP): click "More info", then "Run anyway". If it keeps
  happening, use the terminal route instead of double-clicking.
- "out of memory" mentioning CUDA or the GPU: close other heavy programs, especially Imaris, and
  run it again. Only one image can run at a time on one graphics card.

Mac:

- macOS refuses to open `quantify.command` ("cannot be opened", "not opened", or "Apple could not
  verify"): this is the once-only check described in step 1 of the short version. macOS 15 or
  newer: click Done, then System Settings, Privacy & Security, Open Anyway. macOS 14: right-click
  the file, Open.
- A message that `quantify.command` "could not be executed because you do not have appropriate
  access privileges" (or "Permission denied" from Terminal): in Terminal, in the program's folder,
  run `chmod +x quantify.command` once (step 7 of the Mac set-up).
- `check-gpu` says "No accelerator": torch was installed without Metal support, or the Mac has an
  Intel chip. Redo step 4 of the Mac set-up exactly as written (torch must come from the plain
  install, not the Windows `cu128` address). Intel Macs can still run everything, only much more
  slowly, on the processor.
- The window says "Cellpose failed on mps ... retrying on the CPU" (the message mentions "MPS
  backend out of memory"): the image did not fit in the graphics memory. The run does not stop; it
  carries on using the processor, so the nucleus step can take hours instead of minutes. If you
  would rather not wait, stop it (Control-C), close other applications and run it again, or run
  that image on the workstation. For a quick look only, `--xy-stride 2` runs at half resolution,
  but never use those numbers.
- The first run stops with a message about "command line developer tools" or `xcrun`: click
  Install in the dialog that appears and run the command again when it has finished.
- A run that worked before is now very slow: check that `check-gpu` still says Metal; a macOS
  update can require reinstalling torch (step 4 of the Mac set-up).

## Running many gonads on the university cluster

One gonad is not faster on RMACC Alpine than on the workstation; the point of the cluster is
running a hundred gonads at once, one GPU job per image, while the workstation stays free. The
same container runs on the workstation (RTX 5090), the A100 and the L40. The build, transfer and
SLURM steps are in `docs/HPC_ALPINE.md`; the batch workflow is `workflow/Snakefile`, and the
`alpine` and `slurm` folders in `config/profiles/` are its job settings (passed to Snakemake with
`--workflow-profile`, not to `--profile`).

## Getting the nucleus model (both machines)

The trained nucleus model is a single file, 1.2 GB, named exactly `germline_nuclei_combined` with
no extension. It is not in the download; it is published as a release on the GitHub page of this
project, and the program can fetch and check it for you. In the terminal, in the program's folder,
after the install steps below:

Windows:
```
.venv\Scripts\germquant.exe get-model
```
Mac:
```
.venv/bin/germquant get-model
```

It downloads the file into `models/models/` inside the program's folder, checks it byte for byte,
and prints "verified and installed". Run it again any time; if the file is already there and
correct it says so and downloads nothing. If you were given the file some other way, put it at
`models/models/germline_nuclei_combined` yourself and run the same command to check it.

## Setting up a Windows computer (one time, for a lab tech)

Requirements: 64-bit Windows, an NVIDIA GPU (RTX 30, 40 or 50 series), about 10 GB free disk.

1. Install Python 3.11 from python.org and tick "Add Python to PATH" during the install.
2. Install the `uv` helper: in PowerShell, `pip install uv`.
3. Get this code: download the ZIP from GitHub and unzip it (the folder is named
   `C_elegans_germline_masking_and_spot_counting-main`; move it somewhere sensible and you may drop
   the `-main`), or `git clone` it.
4. In PowerShell, go into the folder (`cd`, a space, drag the folder in, Enter), then run these one
   at a time. The second and third take several minutes and print a lot:
   ```
   uv venv --python 3.11
   uv pip install torch torchvision --index-url https://download.pytorch.org/whl/cu128
   uv pip install -c constraints/desktop-2026-09.txt --extra-index-url https://download.pytorch.org/whl/cu128 -e ".[gpu,sc]"
   uv pip install -c constraints/desktop-2026-09.txt --extra-index-url https://download.pytorch.org/whl/cu128 spotmax cellacdc
   ```
   The constraints file pins every package to the versions the workstation was validated with.
5. Get the nucleus model: `.venv\Scripts\germquant.exe get-model` (see "Getting the nucleus model").
6. In the same window, confirm the GPU is seen: `.venv\Scripts\germquant.exe check-gpu` should
   name your card and print a line starting with "OK".
7. Test the whole thing by dragging a small `.nd2` onto `quantify.bat`.

Linux and cluster installs use the Dockerfile or `apptainer.def`; `environment.yml` and
`pixi.toml` pin the same environment.

## Setting up a Mac (one time)

Requirements: a Mac with an Apple chip (Apple menu, About This Mac, says "Chip: Apple M..."),
macOS 14 or newer, 16 GB of memory or more, about 8 GB free disk. Everything below is typed in
Terminal (Command and Space, type `Terminal`, Return). Remember this has not been run on a real Mac
yet; if a step does not match what you see, tell Ryan what happened.

1. Install Python 3.11: download the macOS installer from python.org and run it, or, if you use
   Homebrew, `brew install python@3.11`. Check with `python3.11 --version`.
2. Install the `uv` helper. If you have Homebrew: `brew install uv`. Otherwise paste this line into
   Terminal, then close Terminal and open it again:
   ```
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```
   Check with `uv --version`.
3. Get this code: download the ZIP from GitHub and unzip it (Safari puts it in Downloads, in a
   folder named `C_elegans_germline_masking_and_spot_counting-main`; move it to your home folder
   and you may drop the `-main`), or `git clone` it. If a dialog offers to install the "command
   line developer tools", click Install and wait for it to finish.
4. Go into the folder (`cd`, a space, drag the folder from Finder into the window, Return), then
   run these one at a time. The second and third take several minutes and print a lot:
   ```
   uv venv --python 3.11
   uv pip install -c constraints/mac-arm64-2026-09.txt -e ".[gpu,sc]"
   uv pip install -c constraints/mac-arm64-2026-09.txt spotmax cellacdc
   ```
   Do not use the `download.pytorch.org/whl/cu128` line from the Windows instructions: that is
   an NVIDIA build and there is no NVIDIA on a Mac. The plain install above brings the Mac
   version of torch, which drives the graphics chip through Metal.
5. Get the nucleus model, in the same window: `.venv/bin/germquant get-model` (see "Getting the
   nucleus model").
6. In the same window (if you closed it, open Terminal and repeat the `cd` from step 4), confirm
   the graphics chip is seen:
   ```
   .venv/bin/germquant check-gpu
   ```
   The output should contain "Apple Silicon GPU via Metal" and a line starting with "OK". If it
   says "No accelerator", see "If something breaks".
7. Make the point-and-click file runnable, once:
   ```
   chmod +x quantify.command
   ```
8. Test the whole thing: double-click `quantify.command` (first time: the once-only macOS check
   described in the short version), drag a small `.nd2` into the Terminal window, press Return.

What to expect on a Mac: the nucleus segmentation is the slow part and should take several times
longer than on the workstation; everything else runs at about the same speed. The RAD-51 spot
counting runs on the processor on every machine. The technical details (device choice, the
`GERMQUANT_DEVICE` setting, what has and has not been tested on a Mac) are in `docs/MAC_METAL.md`.

## Words you might not know

- `.nd2`: the raw image file the Nikon confocal saves. It holds all channels and all z planes.
- Terminal, PowerShell, "the window": where you type commands. Scrolling text is normal.
- `cd`: the command that moves the window into a folder, so that the short names in the commands
  (`.venv`, `config`) are found. Type `cd`, a space, then drag the folder in.
- Path: a file's full address, for example `C:\Users\you\images\worm1.nd2` on Windows or
  `/Users/you/images/worm1.nd2` on a Mac. Drag a file into the window to paste its path.
- Finder: the Mac's file browser (the smiling face in the Dock). File Explorer: the Windows one.
- `.csv`: a spreadsheet file; double-click opens it in Excel (or Numbers on a Mac).
- Segmentation: the computer outlining each nucleus in 3D.
- Label image: a 3D image in which every voxel holds the number of the nucleus it belongs to (0 is
  background). This is what Imaris turns into Surfaces.
- RAD-51 focus (plural foci), spot: the dots being counted; they mark DNA double-strand breaks
  during meiosis.
- Germline isolation: deciding which segmented nuclei are part of the gonad and which are gut,
  debris or other tissue.
- Voxel: one 3D pixel. Here 0.2 um in z and about 0.11 um in x and y; the program reads this from
  each file rather than assuming it.
- GPU, graphics chip, Metal, CUDA: the part of the computer that does the nucleus segmentation
  quickly. NVIDIA cards are driven through CUDA (Windows, cluster), Apple chips through Metal.

## For developers

Pipeline: read `.nd2`, segment nuclei (Cellpose, fine-tuned germline model), isolate the germline,
fit the gonad axis, count spots (SpotMAX), optionally segment P granules and compute SYP / PGL-1
colocalization, and the optional stages (lamin envelope, mask audit, hand-trace staging, SC tracer,
partition coefficient, lit fraction, extra spot instances); write tidy CSV and Parquet tables with
provenance and render a QC montage. The stage list and switches are in `src/germquant/stages.py`;
the code is a Python package under `src/germquant/` (`cli.py`, `pipeline.py`, `device.py`,
`models_registry.py`, and one sub-package per stage: `io`, `segment`, `germline`, `envelope`,
`axis`, `staging`, `spots`, `sc`, `granule`, `coloc`, `partition`, `measure`, `qc`, `render`,
`validate`); `schema.py` declares every output column (append-only). The roadmap that produced this
layout is `docs/ROADMAP_modular_pipeline.md`.

Documents: `docs/ARCHITECTURE.md` (design), `docs/RUNBOOK.md` (full-resolution GPU run and spot
calibration), `docs/SPOTMAX_VALIDATION.md` (how the spot counts were validated against Imaris and
what is still unvalidated), `docs/COLOCALIZATION.md` (the coloc stage), `docs/SC_TRACING.md` (the SC
tracer), `docs/MAC_METAL.md` (Apple Silicon), `docs/ANNOTATION.md` (annotating nuclei and
fine-tuning the Cellpose model), `docs/HPC_ALPINE.md` (cluster runs), `docs/GENERIC_PROFILE.md`
(the generic profile for samples that are not a worm germline), `analysis/coloc/README.md`
(the SYP-3 / PGL-1 analysis layer).

Configuration: `config/config.yaml` holds every tunable value with a comment saying why it is set
the way it is; the `spots` block holds the cross-validated detection parameters. Channel maps are in
`config/channel_maps/`, profiles in `config/profiles/`. Lengths and volumes are always in microns;
the voxel size is read from each file and threaded through every 3D operation. The config text is
hashed into every output, so the shipped config files are not edited casually: new switches are
read with code defaults and profiles overlay them. The trained model is a GitHub release asset
registered in `src/germquant/models_registry.py` (tag, checksum, size); a new model version is a new
release and a new entry there.

Validation and cross-validation tools are in `scripts/` (`cv_*.py`, `validate_*.py`,
`calibrate_spots.py`, `regression_diff.py`, `make_goldens.py`) with the Imaris readers in
`germquant.validate`. Unit tests: `pytest` from the repository root; they run on CPU and need no GPU
or data, and `tests/test_golden.py` guards the numbers against the frozen references. Pinned
environments: `constraints/desktop-2026-09.txt` (workstation), `constraints/mac-arm64-2026-09.txt`
(Apple Silicon).

Other CLI commands: `validate` (compare a pipeline table or label image with hand-scored ground
truth), `prep-training` (export DAPI slices for annotation), `finetune` (fine-tune the Cellpose
model on a labelled folder), `check-gpu`, `get-model`, `collect`, `trace`, `restage`. Run the
program with `--help` for the options.
