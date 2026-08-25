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

## The short version: one image, drag and drop

1. In this folder, find the file `quantify.bat`.
2. Drag your `.nd2` image onto `quantify.bat` and let go.
3. A black window opens and starts working. Leave it open. One image takes about 15 to 25 minutes.
4. When it prints DONE, the results are in a new folder next to your image, named
   `yourimagename_results`.

That is the whole procedure for a single RAD-51 image. The rest of this page explains the results,
how to run many images, how to run other kinds of images, and what to do when something goes wrong.

## What you need

- The lab workstation with the NVIDIA graphics card. It is already set up. On a different
  computer, see "Setting up a new computer" at the end, or ask whoever manages the lab computers.
- A per-gonad `.nd2` file. The whole-slide overview files (names containing `10x` or
  `largeimage`) will not work; use the individual gonad files.
- About 20 minutes per image. The computer does the work; you wait. Do not start Imaris or another
  heavy program on the same machine while it runs (they compete for the graphics card memory).

## Running from PowerShell

Use this when drag and drop does not work, when you want to choose where the results go, when your
image is not the standard DAPI / SYP / RAD-51 panel, or when you want to process a whole folder.

1. Open PowerShell: click the Start menu, type `PowerShell`, click "Windows PowerShell". A black
   window opens.
2. Go to the program's folder. Copy this line, paste it into the window (right-click pastes), press
   Enter:
   ```
   cd "C:\Users\ryane\C_elegans_germline_masking_and_spot_counting"
   ```
   If the code lives somewhere else on your computer, use that folder instead.
3. Run one image. Copy the line below, replace the two paths with yours, press Enter:
   ```
   .venv\Scripts\germquant.exe run "C:\path\to\YOUR_IMAGE.nd2" --config config\config.yaml --out "C:\path\to\RESULTS_FOLDER"
   ```
   The first quoted path is your image. Instead of typing it, drag the `.nd2` file into the
   window and PowerShell pastes the correct path. The path after `--out` is any folder where you
   want the results saved (it is created if it does not exist). Paths with spaces are fine as long
   as they are inside quotes.
4. Wait 15 to 25 minutes. When the blinking prompt comes back, it is finished.

### A whole folder at once

```
.venv\Scripts\germquant.exe batch "C:\path\to\FOLDER_OF_IMAGES" --config config\config.yaml --out "C:\path\to\RESULTS_FOLDER"
```

This processes every `.nd2` under the folder, including sub-folders, one image at a time, and
skips the `10x` and `largeimage` overviews on its own. It writes one results sub-folder per image
and a `batch_summary.csv` listing them all with their nucleus counts and QC flags. Many images can
take hours; leave the window open. If it is interrupted, images that already finished keep their
results folders, but running the same line again reprocesses everything, so point `--out` at a new
folder or move the finished ones aside first.

### Segmentation only (no spot counting)

Add `--no-spots` to either command. It skips the RAD-51 spot detection, so it is much faster and
cannot get stuck on a difficult image. You still get the nucleus table, the 3D nucleus label image
for Imaris, and the montage; the spots table is empty and no spots image is written.

```
.venv\Scripts\germquant.exe run "...\YOUR_IMAGE.nd2" --config config\config.yaml --out "...\RESULTS" --no-spots
```

`--no-coloc` skips the PGL-1 / SYP colocalization stage in the same way.

### Checking what is in an image

```
.venv\Scripts\germquant.exe info "C:\path\to\YOUR_IMAGE.nd2"
```

prints the channel names and the voxel size without loading the pixels. Use it whenever you are not
sure which config to pick (next section).

## Which config file to use

The config file tells the program which channel is which and holds every tunable setting. Pick the
one that matches your image; everything else stays the same.

| Config | Image type | What runs |
|---|---|---|
| `config\config.yaml` | N2-style DAPI / SYP / RAD-51, with or without a fourth PGL-1 channel | nuclei, germline isolation, RAD-51 spots; colocalization only if a PGL-1 channel is present |
| `config\config_ccw77.yaml` | ccw77 four-colour IF: DAPI / PGL-1::GFP (477) / SYP-3::mCherry (545) / LMN-1 (640) | nuclei, germline isolation, SYP-3 / PGL-1 colocalization with the lamin envelope; there is no RAD-51 channel so spot counting skips itself |
| `config\config_n2dryice.yaml` | the N2 dry-ice test set | as `config.yaml` |

The channel assignments themselves are in `config\channel_maps\`. Channels are matched by name
(the laser line, for example `405`, `477`, `545`, `640`) with a fixed-position fallback. If
`germquant info` shows names that are not in the map, or the montage shows the wrong channel in the
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
| `..._coloc.csv` | The SYP / PGL-1 colocalization metrics, one row per operand (only with a PGL-1 channel). See `docs\COLOCALIZATION.md` for what each column means. |
| `..._image_summary.csv` | One row of totals for the image: nucleus counts, mean spots per nucleus, the headline colocalization numbers, QC pass or fail and the flags behind it. |
| `..._montage.png` | The picture to look at every time: each channel, the nucleus outlines and the detected spots side by side. If this looks wrong, the numbers are wrong. |
| `..._nuclei_labels.tif` | The 3D nucleus map. Load it into Imaris as Surfaces to check the segmentation. |
| `..._spots.tif` | The detected spots as a 3D image on the same grid. Load it into Imaris as a channel (Edit, Add Channels) or run Imaris Spots on it (diameter about 0.4 um) to get Spots objects next to your Surfaces. |
| `.parquet` copies | The same tables in a compact format for R or Python. |

To get the average number of RAD-51 spots per nucleus: open `..._nuclei.csv` in Excel, filter to
`in_germline = True`, and average the `n_spots` column.

Every table also carries the image metadata parsed from the file name (date, genotype, treatment,
sex, replicate), the voxel size, the channel map that was used, and the version of the code, so a
result can always be traced back to how it was made.

### Before you trust the numbers

The spot counts have been checked against Imaris counts for N2 worms, with and without heat shock,
and they agree to about one spot per nucleus across 4 to 21 foci per nucleus. They have not been
checked for other genotypes. Before reporting absolute counts on a mutant, count a few gonads in
Imaris and compare; the detection settings were tuned on N2 and may need re-tuning. The details and
the caveats are in `docs\SPOTMAX_VALIDATION.md`.

The `in_germline` decision uses SYP signal and spatial connectivity, never the spot count; on 14
validation gonads the nuclei it kept held about 97 percent of the RAD-51 spots. It can admit sperm
masses and somatic nuclei that touch the gonad; for the ccw77 colocalization analysis that mattered,
and the audit and the fix are described in `analysis\coloc\README.md`.

## The colocalization stage

When the channel map resolves a PGL-1 channel, the run also segments the P granules, builds a thin
cytoplasmic shell around each germline nucleus (anchored to the lamin envelope when a lamin channel
is present), and measures how much SYP signal coincides with the granules inside that shell. The
headline columns are `shell_pearson`, `shell_manders_m1` and `shell_manders_m2` in
`..._image_summary.csv`; the per-operand detail is in `..._coloc.csv`. The reasoning behind the
region choice and the meaning of every column are in `docs\COLOCALIZATION.md`.

The stage-resolved analysis of SYP-3 partitioning into P granules under heat shock (partition
coefficient, hand-traced pachytene zones, controls, figures) is a separate layer built on these
outputs and lives in `analysis\coloc\`.

## If something breaks

- A red error saying the `.venv` folder is missing: the program is not installed on this
  computer. See "Setting up a new computer".
- "running scripts is disabled on this system", or Windows blocked the file: use the PowerShell
  route instead of double-clicking, or ask a lab tech to allow scripts.
- "file not found": the image path is wrong. Drag the `.nd2` into the window to paste it correctly
  and keep the quotes.
- It found 0 nuclei, or the montage shows the wrong channels: the image's channel names do not
  match the config. Run `germquant info` on the file and compare with `config\channel_maps\`. Tell
  whoever maintains the tool.
- The black window closed instantly: run it from PowerShell so the error stays on screen.
- "out of memory" mentioning CUDA or the GPU: close other heavy programs, especially Imaris, and
  run it again. Only one image can run at a time on one graphics card.
- The spot step runs for hours on one image: that image is probably flooded with signal or an
  artefact. The run caps the number of candidate spots so it should still finish; if you only need
  the segmentation, rerun with `--no-spots`.
- Anything else: copy the red text and send it to whoever maintains the tool.

## Running many gonads on the university cluster

One gonad is not faster on RMACC Alpine than on the workstation; the point of the cluster is
running a hundred gonads at once, one GPU job per image, while the workstation stays free. The
same container runs on the workstation (RTX 5090), the A100 and the L40. The build, transfer and
SLURM steps are in `docs\HPC_ALPINE.md`; the batch workflow is `workflow\Snakefile` with the SLURM
profile in `config\profiles\`.

## Setting up a new computer (one time, for a lab tech)

Requirements: 64-bit Windows, an NVIDIA GPU (RTX 30, 40 or 50 series), about 10 GB free disk.

1. Install Python 3.11 from python.org and tick "Add Python to PATH" during the install.
2. Install the `uv` helper: in PowerShell, `pip install uv`.
3. Get this code: download the ZIP from GitHub and unzip it, or `git clone` it.
4. In PowerShell, `cd` into the folder, then run these one at a time:
   ```
   uv venv
   uv pip install -e .
   uv pip install spotmax cellacdc
   uv pip install torch --index-url https://download.pytorch.org/whl/cu128
   ```
5. Put the trained nucleus model at `models\models\germline_nuclei_combined` (about 1.2 GB; it is
   not in the repository, ask Ryan for it).
6. Confirm the GPU is seen: `.venv\Scripts\germquant.exe check-gpu` should print your card.
7. Test the whole thing by dragging a small `.nd2` onto `quantify.bat`.

Linux and cluster installs use the Dockerfile or `apptainer.def`; `environment.yml` and
`pixi.toml` pin the same environment.

## Words you might not know

- `.nd2`: the raw image file the Nikon confocal saves. It holds all channels and all z planes.
- Terminal, PowerShell, "the black window": where you type commands. Scrolling text is normal.
- Path: a file's full address, for example `C:\Users\you\images\worm1.nd2`. Drag a file into the
  window to paste its path.
- `.csv`: a spreadsheet file; double-click opens it in Excel.
- Segmentation: the computer outlining each nucleus in 3D.
- Label image: a 3D image in which every voxel holds the number of the nucleus it belongs to (0 is
  background). This is what Imaris turns into Surfaces.
- RAD-51 focus (plural foci), spot: the dots being counted; they mark DNA double-strand breaks
  during meiosis.
- Germline isolation: deciding which segmented nuclei are part of the gonad and which are gut,
  debris or other tissue.
- Voxel: one 3D pixel. Here 0.2 um in z and about 0.11 um in x and y; the program reads this from
  each file rather than assuming it.

## For developers

Pipeline: read `.nd2`, segment nuclei (Cellpose, fine-tuned germline model), isolate the germline,
fit the gonad axis, count spots (SpotMAX), optionally segment P granules and compute SYP / PGL-1
colocalization, write tidy CSV and Parquet tables with provenance, and render a QC montage. The code
is a Python package under `src\germquant\` (`cli.py`, `pipeline.py`, and one sub-package per stage:
`io`, `segment`, `germline`, `axis`, `spots`, `granule`, `sc`, `coloc`, `measure`, `zones`, `qc`,
`render`, `validate`); `schema.py` declares every output column.

Documents: `docs\ARCHITECTURE.md` (design), `docs\RUNBOOK.md` (full-resolution GPU run and spot
calibration), `docs\SPOTMAX_VALIDATION.md` (how the spot counts were validated against Imaris and
what is still unvalidated), `docs\COLOCALIZATION.md` (the coloc stage), `docs\ANNOTATION.md`
(annotating nuclei and fine-tuning the Cellpose model), `docs\HPC_ALPINE.md` (cluster runs),
`HANDOFF.md` (current state and how to pick the work back up), `analysis\coloc\README.md` (the
SYP-3 / PGL-1 analysis layer).

Configuration: `config\config.yaml` holds every tunable value with a comment saying why it is set
the way it is; the `spots` block holds the cross-validated detection parameters. Channel maps are in
`config\channel_maps\`. Lengths and volumes are always in microns; the voxel size is read from each
file and threaded through every 3D operation.

Validation and cross-validation tools are in `scripts\` (`cv_*.py`, `validate_*.py`,
`calibrate_spots.py`) with the Imaris readers in `germquant.validate`. Unit tests: `pytest` from the
repository root; they run on CPU and need no GPU or data. The trained model and raw data are not in
the repository.

Other CLI commands: `validate` (compare a pipeline table or label image with hand-scored ground
truth), `prep-training` (export DAPI slices for annotation), `finetune` (fine-tune the Cellpose
model on a labelled folder), `check-gpu`. Run `.venv\Scripts\germquant.exe --help` for the options.
