# Time-lapse quantification

Status: the reader ("Reading one timepoint"), the series loop with its resume logic and command line switch
("Running a series"), the live GFP::SYP-2 profile and the nucleus training loop ("Segmenting without DAPI"),
the drift registration ("Registration") and its piecewise variant along the gonad axis ("Piecewise
registration"), the tracking exchange ("Tracking") and the biologist validation ("Validating the linking") are
implemented and tested on synthetic series. The first real series (the October 2026 DLW93 hermaphrodite) has
been run through the loop, the registration, the linking and the piecewise registration (results and their
limits in `C:/Users/ryane/live_results/FIRST_LOOK.md`); the hand-followed validation set has been picked but not
yet followed, so no linking result is validated. The design is [DESIGN_timelapse.md](DESIGN_timelapse.md).
Nothing on this page changes a fixed image: its id, its tables, its stage record and its completion marker
are byte for byte what they were, and tests/test_golden.py still gates that.

The first dataset is the October 2026 GFP::SYP-2 series on the lab NAS: 25 timepoints 300 s apart, 76 (or
66 in the second worm) planes at 0.2 um, two channels (Brightfield, EGFP), 2720 x 2720 pixels at 0.108 um,
56 GB (49 GB) per file. One timepoint is a 2.25 GB stack at full resolution.

## Reading one timepoint

### The keyword

Every reader takes `t`, the index of the timepoint to read, and returns the usual `(C, Z, Y, X)` Stack of
that timepoint with its voxel size, so every stage sees exactly what it sees for a fixed image:

```python
from germquant.io import read_stack, read_image, read_nd2_metadata, read_image_metadata

st = read_stack("series.nd2", t=7)                       # timepoint 7, full resolution
st = read_stack("series.nd2", t=7, xy_stride=4, z_range=(10, 60))   # sub-sampled as always
st = read_image("series.ome.tif", t=7)                   # the same through the format dispatcher
st.t_index                                               # 7 (None when no timepoint was asked for)
```

* `t=None`, the default, is the old behaviour: a T axis takes index 0, as it always did, and a file
  without a T axis is read by exactly the code that read it before. `t=0` on a fixed image is allowed
  and gives the same stack, so a loop can treat a fixed image as a one-timepoint series.
* `t` out of range, or `t` other than 0 on a file without a T axis, raises ValueError before any pixel is
  read, with the number of timepoints in the message. A non-integer `t` is a ValueError too.
* Only that timepoint's planes are read. nd2 chunks its lazy array one frame per (T, Z) with all channels
  in the frame, and the T index is applied before compute, so a 76-plane timepoint touches 76 of the
  1900 frames of the file. `xy_stride` and `z_range` cut the read down further as they do for a fixed image.
* Formats other than `.nd2` (OME-TIFF, CZI, LIF, OME-Zarr through bioio) take the same keyword; the
  timepoint is selected on the lazy array. Without `t`, a file with more than one time point still raises
  UnsupportedImageError, as it did before, and the message now says to pass `t`.

### What the header says about time

`read_nd2_metadata` and `read_image_metadata` return three keys besides the old ones:

| key | meaning | fixed image |
|---|---|---|
| `n_timepoints` | number of timepoints (the size of the T axis) | 1 |
| `period_s` | the programmed interval between timepoints in seconds, from the experiment's TimeLoop (`periodMs`); None for a non-equidistant loop with more than one phase | None |
| `t_seconds` | one value per timepoint: the recorded time stamp of its first plane in seconds since the acquisition started, as NIS-Elements stores it per frame | None |

`t_seconds` is what the camera recorded, not `k * period_s`: in the October series the steps run from
299.5 to 300.5 s (median 300.005 s), and the first stamp is 0.358 s, not 0. It is read from one small
per-frame metadata chunk per timepoint, no pixels. For the bioio formats `period_s` and `t_seconds` are
None for now (not read yet); `n_timepoints` comes from the header. Passing `t` to a metadata call only
checks that the timepoint exists.

`germquant info series.nd2` prints a `timepoints : 25 every 300 s` line for a series and nothing new for
a fixed image.

### Running one timepoint through the pipeline

```python
from germquant import api

res = api.run("series.nd2", profile="n2_sc", out="results/series", t=7)
res["image_id"]          # "series_t007"
```

`process_image(path, cfg, out_dir, t=7)` and `api.run(..., t=7)` run every enabled stage, unchanged, on
that timepoint's stack. Three things differ from a fixed image, and only when `t` is given:

* The image id becomes `<stem>_t<NNN>` (three digits, one underscore): `series.nd2` at `t=7` is
  `series_t007`. The design columns (acquisition date, genotype, treatment, sex, replicate) are parsed
  from the raw file stem first, so they are what they would be for the whole file, and the suffix is added
  afterwards. Every per-image file is named from the suffixed id (`series_t007__nuclei.csv`,
  `series_t007__stages.json`, `series_t007__done.json`, the montage, the label images), so the
  timepoints of one series can share a folder.
* `t_index` (the index) and `t_seconds` (that timepoint's recorded stamp, None when the file carries
  none) are written to the run geometry in `run_manifest.json`, to `<id>__stages.json` and to the
  completion marker `<id>__done.json`. No table gets a new column: `image_id` is the key, as it is for
  every other per-image file.
* `germquant run` has no `--t` switch yet; it arrives with the series loop. Use the API, or
  `process_image`, until then.

tests/test_timelapse_reader.py checks that a timepoint run and a fixed run on the same stack write
identical tables apart from `image_id` and the time stamp, and that the fixed run's manifest, stage
record and marker contain no timepoint fields.

### Cost on the real series

Read-only check on 2026-10-09, workstation reading the 56 GB file over the lab network, nothing copied:

* Header (`read_nd2_metadata`): 0.28 s, no pixel frame read; `n_timepoints` 25, `period_s` 300.0,
  25 time stamps from 0.358 to 7200.367 s. The second worm (49 GB, 66 planes) header: 0.14 s.
* `read_stack(t=7, xy_stride=8)`: 25 s. Exactly the 76 frames of timepoint 7 were read (sequence
  indices 532 to 607, no other timepoint touched), 2.25 GB of pixels at about 90 MB/s over the network;
  the stack in memory is 35 MB at stride 8. Resident memory after the read was 0.36 GB; the peak working
  set was 2.64 GB, which is the memory-mapped file pages of the 76 frames (nd2 copies each frame out of
  the map while it is open, and the map goes away when the read finishes).
* At full resolution one timepoint is 2.25 GB in memory, like a fixed gonad stack; the read time is set
  by the network, not by the stride, since a whole frame is read before it is strided. A 25-timepoint
  loop reads for about 10 minutes over the network; a local copy of the file is the cheaper way to iterate.

## Running a series

```
germquant run SERIES.nd2 --config config/config.yaml --out results/SERIES --t all
germquant run SERIES.nd2 --profile live_syp2 --out results/SERIES --t 0-24 --resume
germquant run SERIES.nd2 --config config/config.yaml --out results/SERIES --t 7 --xy-stride 4
```

`--t` names the timepoints: `all` (every timepoint of the file), one index (`7`), or an inclusive range
(`0-24`). Without `--t` a run is the single-image run it always was: a file with a T axis is read at
its first timepoint, as before, and nothing below applies. `--xy-stride`, `--z-range` and the
`--no-<stage>` switches apply to every timepoint of the series.

What happens, per timepoint k:

1. The reader loads only that timepoint's planes (`read_stack(path, t=k)`; the nd2 dask array is
   indexed before anything is computed, so a 56 GB series never comes into memory as a whole).
2. `process_image` runs every enabled stage on that stack, exactly as it does on a fixed image, into
   the folder `results/SERIES/t<NNN>/` (three digits, `t007`).
3. The image id is `<stem>_t<NNN>`, with a single underscore. The design columns (genotype, sex,
   treatment, replicate, acquisition_date) are parsed from the raw stem first, so they read the same as
   they would for a fixed image with that name. Every per-image file in the folder carries the suffixed
   id: `<stem>_t007__nuclei.csv`, `__spots.csv`, `__image_summary.csv`, `__montage.png`,
   `__nuclei_labels.tif`, `__stages.json`, `__done.json` and the rest of the fixed-image file set, plus
   the folder's own `run_manifest.json`. No per-image table gains a column.
4. `t_index` and `t_seconds` (the acquisition time of the timepoint from the file header, when it has
   one) are written into the run geometry of `run_manifest.json` and `<id>__stages.json` and into the
   completion marker `<id>__done.json`. They are written only when `--t` was given, so the markers of
   a fixed image keep their exact keys.
5. With `timelapse.keep_reference_volume` (code default true) a downsampled copy of the segmentation
   channel, and of the brightfield channel when a role of that name resolves, is kept next to the file set
   as `t<NNN>/<id>__ref_stride<k>.npy` for `germquant register` (see "Registration"). A fixed image never
   writes one.

Next to the timepoint folders the loop keeps two series-level records:

* `results/SERIES/<stem>__series.json`: the timepoints requested, done and failed (each failure with
  its flags), `period_s` and `t_seconds` from the header, the elapsed seconds per timepoint,
  `config_hash`, `pipeline_version`, `git_sha`, the run geometry and the enabled stages. It is
  rewritten after every timepoint, so an interrupted series leaves an honest record of what is on disk
  (`complete` is false until the loop finishes).
* `results/SERIES/<stem>__series_summary.csv` and `.parquet`: one row per requested timepoint, see
  "Series tables".

A timepoint whose run raises (an unreadable frame, a segmentation that blows up) is recorded with an
`EXCEPTION:<type>: <message>` flag and the loop goes on, as `germquant batch` does for an image; its
folder has no completion marker, its summary row has `qc_pass` false and empty counts. A stage that
fails without stopping the image (`<stage>:FAILED_<type>` in `qc_flags`) is a finished timepoint whose
marker lists the failed stage. The exit code of `germquant run --t` is 0 when at least one requested
timepoint finished and 1 when none did.

`--resume` skips every requested timepoint whose completion marker matches the current config hash,
enabled stages, xy stride, z range, timepoint index and segment stage (model file and Cellpose
version), the rule `germquant batch --resume` uses, and reprocesses the rest, including timepoints
whose marker lists a failed stage. The summary rows of skipped timepoints are read back from their
`__image_summary.csv`, so the series summary is complete after a resume; their elapsed times are
carried over from the previous series record.

Self-test, from the design: `germquant run SERIES.nd2 --t 7 --out A` writes the same `A/t007/`
folder, file for file, as `germquant run SERIES.nd2 --t all --out B` writes under `B/t007/` (only the
timestamps and elapsed times differ). The test suite checks this on a synthetic three-timepoint series
with a planted drift (`tests/test_timelapse_series.py`).

`germquant collect results/` stacks the timepoint folders like any images (their ids carry the
`_t<NNN>` suffix, so `batch_nuclei.csv` and `batch_summary.csv` hold every timepoint) and, when any
`<stem>__series_summary.csv` is present under the root, stacks those into `batch_series_summary.csv`.
A results tree without a series gets no such file.

From code, `germquant.api.process_series(path, cfg, out, timepoints="all", xy_stride=1, z_range=None,
resume=False)` is the function the command runs; it returns the series id, the requested, done and
failed timepoints, the record and the summary DataFrame. `api.parse_timepoints("0-24", n_timepoints)`
turns a `--t` value into indices. A single timepoint as a plain image, with the suffixed id but no
`t<NNN>` folder and no series record, is `api.process_image(path, cfg, out, t=7)`.

What the reader provides (the contract shared with the reader package): `read_stack(path, t=None)`
keeps today's behaviour (a T axis is read at index 0); `read_stack(path, t=k)` reads timepoint k and
raises ValueError when k is out of range, or when the file has no T axis and k is not 0 (so `--t 0` on
a fixed image is a one-timepoint series, which is how the loop is checked on a local fixed `.nd2`).
`read_nd2_metadata` reports `n_timepoints` (1 for a fixed image), `period_s` (from the experiment's
TimeLoop) and `t_seconds` (per timepoint, or None). `read_image` / `read_image_metadata` pass `t` on
for the other formats.

Timing: one timepoint of the October series is a full fixed-image run (segmentation of a 76 plane
stack on the GPU, SpotMAX when a `foci` role is mapped), so budget the single-image time per
timepoint; 25 timepoints run overnight. Start with `--t 0 --xy-stride 4` to check the profile and the
channel roles, then `--t all`.

## The series record

`<stem>__series.json` keys: `series_id` (the raw stem), `file_path`, `config_hash`,
`pipeline_version`, `git_sha`, `n_timepoints` (of the file), `period_s`, `t_seconds` (list or null),
`xy_stride`, `z_range`, `enabled_stages`, `resume`, `requested` (indices), `done` (indices with a
completion marker, processed now or skipped by a resume), `failed` (`{"t_index", "flags"}` objects),
`elapsed_s` (seconds per timepoint, keyed by the index as text), `timepoints` (one object per
requested timepoint: `t_index`, `t_seconds`, `image_id`, `out_dir`, `status` done or failed,
`resumed`, `qc_pass`, `qc_flags`, `failed_stages`, `elapsed_s`), `started`, `finished`, `complete`.

## Series tables

`schema.SERIES_TABLES` declares six tables, one file per series, written next to the timepoint folders
as `<stem>__<table>.csv` (and `.parquet`). They live outside `schema.TABLES` on purpose: the stage
tests force every `TABLES` entry through `process_image`, and a timepoint folder holds exactly the
fixed-image file set. The schema rule is the same, append only. Every row carries `SERIES_META`
(`series_id`, `file_path`, the design columns, `n_timepoints`, `period_s`, `pipeline_version`,
`git_sha`, `config_hash`, `run_timestamp`); the join to a timepoint's tables is `image_id`
(`series_id + "_t" + three digits`), and on to a nucleus with `(image_id, nucleus_id)`. The track id
is `time_track_id`; `track_id` already means an SC fragment in `sc_tracks`.

| table | written by | one row per | columns |
|---|---|---|---|
| `series_summary` | this package (`process_series`) | requested timepoint | `t_index`, `t_seconds`, `image_id`, `n_nuclei`, `n_germline_nuclei`, `mean_spots`, `qc_pass`, `qc_flags`, `elapsed_s`, then `sig_<role>_bg` and `sig_<role>_snr` for each channel role the series resolves (the signal_qc stage's background and signal-to-noise), then `SERIES_META` |
| `register_shifts` | registration package (`germquant register`) | timepoint | the step against the reference frame and the accumulated shift to t0 in microns, the brightfield cross-check and its disagreement flag, whether the step was capped, the median nearest-neighbour displacement before and after, stride, upsample factor, method |
| `track_points` | `germquant track-link` / `track-import` | nucleus per timepoint | `time_track_id`, `t_index`, `t_seconds`, `image_id`, `nucleus_id`, the raw centroid (`z_um`, `y_um`, `x_um`), the registered centroid (`reg_z_um`, `reg_y_um`, `reg_x_um` and, appended, the shared-contract names `zr_um`, `yr_um`, `xr_um`), the step from the track's previous point (`step_um` and its components), `is_gap`, `link_source`; then as extras the displacement that was subtracted (`disp_z_um`, `disp_y_um`, `disp_x_um`: raw = registered + displacement), the nucleus measurements joined on `(image_id, nucleus_id)` (`volume_um3`, `in_germline`, `axis_position_um`, `axis_position_aligned_um`, the `dna` / `central_element` mean intensities, `central_element_mean_intensity_norm`, the `sc_*` columns when the sc_trace stage ran) and the fixed-axis positions `axis_position_fixed_um` (registered) and `axis_position_fixed_raw_um` (raw) |
| `tracks` | `germquant track-link` / `track-import` | track | `n_points`, `n_gaps`, `t_start`, `t_end`, `is_full_length`, `path_length_um`, `net_displacement_um`, `speed_um_per_min`, `axis_velocity_um_per_min` (signed, after the t0 sign alignment), `intensity_slope_raw_per_min`, `intensity_slope_norm_per_min`, `sc_length_change_um`, `link_source`; appended: `t_first`, `t_last`, `mean_step_um`, `intensity_slope_raw`, `intensity_slope_norm` (the contract names of the same numbers), `axis_velocity_refit_um_per_min` (the per-frame axis), `axis_velocity_raw_um_per_min` (the fixed axis at the raw centroids, independent of the registration) |
| `track_validation` | `germquant track-validate` | hand-followed point | the biologist's point in raw-frame microns, the label it lands in (`nucleus_id`, NA = segmentation miss) and that label's `time_track_id`, the `outcome` (linked, segmentation_miss, switch, edge_missing, unlinked), the `validation_set` and `zone`; appended, repeated on every point of the hand track: `n_points`, `n_in_label`, `segmentation_recall`, `linking_recall`, `edge_precision`, `switches`, `role` (`drop_duplicates("hand_track_id")` is the per-track table) |
| `register_piecewise` | `germquant register --mode piecewise` | timepoint per axial segment | the segment (`segment`, `n_segments`, its extent `axis_start_um` to `axis_end_um` on the origin's fixed axis, its centre `axis_centre_um`, `n_nuclei_t0`), the cumulative displacement `dz_um`, `dy_um`, `dx_um` and the step `step_d*_um` with `capped` and `filled`, the brightfield cross-check of the same crop, the per-segment quality numbers, the crop box; see "Piecewise registration" |

## Tracking

Linking the nuclei of a series across its timepoints is a file exchange, never a change to the per-timepoint
results (design section 5): the tracking commands read a finished series folder and write only series-level
files next to the `t<NNN>/` folders. A tracker change never reruns segmentation, and nothing is ever
resampled: every per-timepoint table and label image stays in its own raw frame, and the registration enters
as coordinates only. The shared contract between the time-lapse packages:

* the registration package writes `<stem>__register_shifts.csv` with the cumulative shift of every timepoint
  to t0 (`t_index`, `dz_um`, `dy_um`, `dx_um`; a positive shift means the sample moved by that much). A
  **registered coordinate is the raw centroid minus the cumulative shift of its timepoint**. Without the
  table the series is treated as already aligned (zero shifts) and every record says so (`shift_source`
  `none`). Label images are placed on the registered grid by the whole-voxel rounding of the shift;
* `germquant register --mode piecewise` writes `<stem>__register_piecewise.csv` instead: a cumulative
  displacement per axial segment of the germline per timepoint. With `track-link --registration piecewise` a
  **registered coordinate is the raw centroid minus the displacement interpolated along the gonad axis at the
  nucleus's own axis position** ("Piecewise registration" below), and each label is placed on the registered
  grid by the whole-voxel rounding of its own nucleus's displacement. The export, the import and the strips use
  the rigid table only;
* the loop keeps a downsampled copy of the segmentation channel per timepoint,
  `t<NNN>/<id>__ref_stride<k>.npy` (uint16, Z x Y/k x X/k, `timelapse.reference_stride`, code default 8,
  written when `timelapse.keep_reference_volume` is true, which is the default for a series run and never
  for a fixed image). The export and the strips use it as their image; without it they say so (the export
  writes no image, the strips need `--image file`).

What a command reads from the series folder: `<stem>__series.json` (the timepoints listed as done, their
`t_seconds`, `period_s`, the file path), `t<NNN>/<id>__nuclei.csv` (centroids, `in_germline`,
`axis_position_um`, the mean intensities, the `sc_*` columns, the voxel size and the channel map),
`t<NNN>/<id>__nuclei_labels.tif`, `t<NNN>/<id>__image_summary.csv` (`total_germline_length_um`, used when a
timepoint's axis is reversed) and the shift table. A timepoint that failed or was not run is a gap in the
series, not a frame of zeros. Every command takes `--stem` when a folder holds several series, and
`--config` / `--profile` for the `timelapse.*` keys; without them the keys are read from the series'
`run_manifest.json`, and every key has a code default (table below), so no config file needs editing.

### Export: `germquant track-export`

```
germquant track-export results/SERIES --out exchange/SERIES [--crop germline|all] [--stride 2] [--image ref|file|none]
```

writes a Cell-ACDC position folder, `exchange/SERIES/Position_1/Images/`, with

* `<stem>_segm.npz`: a `(T, Z, Y, X)` integer array (uint16 while the largest nucleus id fits) of the
  germline nucleus labels (`in_germline` True; `--crop all` keeps every label) of every timepoint on disk,
  placed on the registered grid, cropped to the union bounding box of the germline across the timepoints
  (`--crop all`: the whole, shifted field) and sub-sampled by `--stride` in y and x (z is kept whole). The
  label values are the `nucleus_id`s of each timepoint, so the array is an untracked segmentation in
  Cell-ACDC's sense: frame-local ids, the tracker's job is to make them consistent;
* `<stem>_<channel>.tif`: the `central_element` channel on the same grid (`--channel-name`, default
  `central_element`), uint16, for display. `--image ref` (default) builds it from the reference volumes with
  a nearest-neighbour resample onto the export grid (coarse: the reference stride is 8, so a nucleus is a few
  pixels), `--image file` reads each timepoint of the series file at the export stride (a full frame per
  timepoint, 25 s each over the lab network for the October series), `--image none` writes no image;
* `<stem>_metadata.csv`: SizeT, SizeZ, SizeY, SizeX, TimeIncrement (`period_s`), the physical voxel size of
  the export grid, `basename` and `channel_0_name`, in Cell-ACDC's `Description,values` form;
* `exchange/SERIES/export.json`: the crop box in registered voxels of the t0 grid, the stride, the exported
  shape and dtype, the shifts applied per timepoint (microns and the rounded voxels), the shift source, the
  frame index of every `t_index`, the label id map of every timepoint (exported value to `nucleus_id`,
  identity but explicit) and the files written.

The germline box is found in a first pass over the label images (each read once, nothing kept), the labels
are placed in a second pass; the exported array is the only large object in memory. An estimate above
`timelapse.track.export_max_gb` (2 GB) is warned about; raise the stride. On the October series the
germline crop at stride 2 is about 45 MB per timepoint before compression.

### Linking: `germquant track-link`

```
germquant track-link results/SERIES --backend overlap|overlap_inhouse|laptrack|none [--gap 2] [--max-dist-um 3.0] [--min-overlap 0.4] [--stride 2] [--all-nuclei] [--registration rigid|piecewise]
```

links the germline nuclei (`--all-nuclei`: every label) across the timepoints on disk and writes
`<stem>__track_points.csv` and `<stem>__tracks.csv` (and `.parquet`) next to the timepoint folders, plus
`<stem>__track_link.json` with the backend, its parameters, the registration and shift source, the counts and
the axis sign per timepoint.

`--registration` (key `timelapse.track.registration`, code default `rigid`) picks the registered frame the
backends link in. `rigid` subtracts each timepoint's cumulative shift from `<stem>__register_shifts.csv` (zero
shifts without the table, `shift_source` `none`). `piecewise` subtracts each nucleus's own displacement,
interpolated along the gonad axis from `<stem>__register_piecewise.csv` (`germquant register --mode
piecewise`; the command stops with exit code 1 and says so when the table is missing); the overlap backends then
move each label by the whole-voxel rounding of its own displacement (`backend_info.placement` in the record),
and the record carries `registration_info` (the segments, the origin, the timepoints with an estimate). The
backends:

* `overlap` (default): Cell-ACDC's own overlap tracker, run headless. `cellacdc.trackers.CellACDC.
  CellACDC_tracker.track_frame`, the function Cell-ACDC's tracker loop calls, is called frame by frame on the
  registered label grid (stride `--stride`, default 2), so only two frames are in memory; it imports without
  a Qt binding in the numeric venv. The rule: a nucleus takes the previous nucleus it covers most, measured
  as intersection over the previous nucleus' volume, when that is at least `--min-overlap` (0.4, Cell-ACDC's
  `IoA_thresh`) and no other nucleus covers that previous nucleus better; every other nucleus starts a new
  track. The first timepoint's track ids are its nucleus ids. When Cell-ACDC cannot be imported the command
  falls back to `overlap_inhouse` and says so in the record (`note`) and on the console;
* `overlap_inhouse`: the same rule in NumPy (`track.overlap_assign`), kept for environments without the
  `spots` extra; the test suite checks both give the same tracks on the synthetic series;
* `laptrack`: LAP linking (laptrack, BSD; `pip install laptrack`) of the registered centroids with the
  squared Euclidean cost, cutoff `--max-dist-um` squared (3.0 um), gap closing over `--gap` frames and
  splitting and merging off. `--gap N` is laptrack's `gap_closing_max_frame_count`: a track ending at t may
  join one starting at t + d when 1 <= d <= N, so the default 2 bridges one missing timepoint (a nucleus the
  segmentation lost for a frame, or a whole timepoint that failed); 0 turns gap closing off. This is what
  closes the single-frame gaps the overlap tracker cannot;
* `none`: no linking; the points table is written with `time_track_id` empty and `link_source` `none`, so
  every measurement can be plotted per timepoint without tracks.

The nearest-neighbour spacing of the packed rows is 3.3 to 4.7 um median and 1.2 um at the 10th percentile
(design section 10): keep `--max-dist-um` below the spacing, and look at the strips before trusting a
larger value.

**track_points** has one row per nucleus per timepoint: `time_track_id`, `t_index`, `t_seconds`, `image_id`,
`nucleus_id`, the raw centroid, the registered centroid (`reg_*_um` and the contract names `zr_um`, `yr_um`,
`xr_um`), `step_um` and its components (the registered displacement from the track's previous point, empty
on a track's first point), `is_gap` (True on a point that follows one or more missing timepoints in its
track), `link_source` (the backend, `imported:<tool>`, or `none` for a nucleus without a track), then the
nucleus measurements joined on `(image_id, nucleus_id)`: `volume_um3`, `n_voxels`, `in_germline`,
`axis_position_norm`, `axis_position_um`, `zone`, `n_spots`, the `sc_*` columns, `dna_mean_intensity` and
`central_element_mean_intensity` when the nuclei table has them, plus the derived columns
`axis_position_aligned_um` (the axis position after the sign alignment below),
`central_element_mean_intensity_norm` (the SYP mean divided by the median over the germline nuclei of the
same timepoint, the bleaching-free quantity), `axis_position_fixed_um` and `axis_position_fixed_raw_um` (the
first timepoint's axis at the registered and at the raw centroid, "One reference axis" below) and `disp_z_um`,
`disp_y_um`, `disp_x_um` (the displacement subtracted from the raw centroid: the timepoint's cumulative shift
under the rigid registration, the nucleus's own interpolated displacement under the piecewise one).

**tracks** has one row per track: `n_points`, `n_gaps` (timepoints missing inside the span), `t_start` /
`t_first`, `t_end` / `t_last`, `is_full_length` (a point at every timepoint on disk), `path_length_um`,
`net_displacement_um`, `mean_step_um`, `speed_um_per_min` (path length over the recorded time span),
`axis_velocity_um_per_min` (the least-squares slope of `axis_position_fixed_um` against minutes, signed:
positive is towards the proximal end in t0's orientation; see "One reference axis" below),
`axis_velocity_refit_um_per_min` (the same slope on `axis_position_aligned_um`, the per-frame axis; kept
for comparison), `axis_velocity_raw_um_per_min` (the slope of `axis_position_fixed_raw_um`: the motion along
the gonad in the microscope frame, whatever the registration; it includes any drift of the whole sample, and a
z jump of the whole field moves it a little through the axis fit), `intensity_slope_raw_per_min` /
`intensity_slope_raw` (the slope of `central_element_mean_intensity`) and `intensity_slope_norm_per_min` /
`intensity_slope_norm` (the slope of the median-normalised intensity; the control baseline for the
heat-shock series), `sc_length_change_um` (last minus first `sc_total_length_um` within the track, NA
without the sc_trace stage), `link_source`. Minutes come from the recorded `t_seconds`, else `t * period_s`,
else the timepoint index.

**Which velocity is the flow.** Every geometric column of the tracks (steps, path length, net displacement,
speed, `axis_velocity_um_per_min`) is measured in the registered frame the linking used. A rigid registration
on GFP follows the mean motion of the germline, so these columns are the motion relative to that mean; a
piecewise registration follows the germline segment by segment, so they are the motion relative to the local
flow (small by construction: that is why it links better). The flow itself is then
`axis_velocity_raw_um_per_min`, or the cumulative displacement per segment in `<stem>__register_piecewise.csv`.
On a series where the worm does not drift (the brightfield stays put), the raw frame is the worm's frame.

**Axis sign alignment.** The axis stage orients each timepoint on its own, with a heuristic (the small,
dense nuclei are distal), and that orientation can flip between frames. For every timepoint the nuclei
linked to t0 are taken and the Spearman rank correlation of their `axis_position_um` at t0 and at that
timepoint is computed: a negative correlation means the orientation is reversed, and the frame's positions
are flipped as `L - p`, with `L` the timepoint's `total_germline_length_um` from its image_summary (its
largest axis position when that is missing). A timepoint with fewer than three nuclei linked to t0 inherits
the previous timepoint's sign. The signs are in `<stem>__track_link.json` (`axis_sign_by_t`) and the command
prints the flipped timepoints. The sign is exact; the per-frame origin is not, which is why the velocity
uses one reference axis instead.

**One reference axis.** The axis stage refits the axis at every timepoint and puts its origin at that
frame's most distal labelled nucleus, so `axis_position_um` of a fixed point in the worm moves from frame
to frame. On the October series at stride 4 a fixed point at x = 100 um sat anywhere from 54 to 83 um
along the per-frame axis, drifting by -0.05 to -0.09 um per minute, and the refitted-axis velocity read
about half the true flow (0.087 against 0.174 um per minute on 135 full-length tracks). The noise does not
average out. `axis_position_fixed_um` therefore applies the FIRST timepoint's axis to every frame: the
position at a registered point is a local linear fit of the first timepoint's `axis_position_um` over
its 16 nearest germline nuclei, evaluated at the point (exact where the axis is locally straight,
extrapolated with the local gradient past the last first-timepoint nucleus). Without axis positions at
the first timepoint, a centerline is fitted to its nuclei with the axis stage's curve fit instead.
tests/test_timelapse_tracking.py plants a flipped frame and checks the velocity; tests/test_tranche3_fixes.py
plants a curved gonad with a jittering per-frame origin and checks that the fixed-axis velocity recovers
the true flow while the refitted one does not.

### Importing a tracker's result: `germquant track-import`

```
germquant track-import results/SERIES --from trackmate spots.csv            (or the TrackMate model .xml)
germquant track-import results/SERIES --from imaris points.csv --time-base 1
germquant track-import results/SERIES --from cellacdc exchange/SERIES/Position_1/Images/<stem>_segm.npz
```

Every importer produces the same two tables with `link_source` `imported:<tool>`; a nucleus the import does
not cover keeps its row with an empty `time_track_id` and `link_source` `none`.

* TrackMate: the spots CSV (`ID`, `TRACK_ID`, `POSITION_X/Y/Z`, `FRAME`; the description and unit rows
  TrackMate writes under the header are skipped; spots without a track are dropped and counted), or the
  model XML (spots from `AllSpots`, tracks from the `AllTracks` edges, `FilteredTracks` honoured when
  present). Positions are physical units (microns).
* Imaris: a points table with `id, t, z, y, x` (aliases accepted: `track_id` / `TrackID`, `t_index` /
  `frame` / `time`, `z_um`...), in raw-frame microns. `--time-base 1` for Imaris's 1-based time index.
* Cell-ACDC: the tracked `segm.npz` (Cell-ACDC overwrites the exported file in place). The export's grid is
  rebuilt from `export.json` (found next to the file or given with `--export-json`) and every exported label
  takes the tracked value covering most of its voxels. A label whose majority is 0 was deleted in the GUI
  (no track); a tracked value claimed by two of our labels in one frame is a merge (the larger overlap keeps
  it, the other gets no track); a label split between tracked values is logged as a split. Objects the
  biologist added have no counterpart in our labels and are counted, not measured. The counts are in the
  record (`backend_info`) and on the console.

Points are mapped to our nuclei by looking the voxel up in that timepoint's label image; `--snap-um D` lets
a point just outside every label take the nearest nucleus centroid within D. Two imported tracks landing
in the same nucleus at one timepoint are a conflict: the first keeps it, the count is logged. `--space
export` says the coordinates were measured on the exported stack (TrackMate or Imaris opened the
exported tif): they are converted back through `export.json` (registered = box origin + coordinate, raw =
registered + the timepoint's shift).

### Strips: `germquant track-strips`

```
germquant track-strips results/SERIES --out strips/SERIES [--size-um 8] [--planes 3] [--image ref|file] [--max-tracks N]
```

writes one PNG per track, `<stem>__track<NNNNN>_strip.png`, with one tile per timepoint of the series: the
`central_element` channel around the nucleus the track names at that timepoint (its raw centroid in that
frame, so the tracked nucleus is always centred), `--size-um` wide, max-projected over `--planes` planes
around the centroid, intensity scaled per track (bleaching shows as dimming along the strip), the timepoint
index under each tile and a dash where the track has no point. An identity swap shows as a jump in texture
and neighbours between two tiles. `--image ref` (default) uses the reference volumes (stride 8: coarse but
honest about what was tracked), `--image file` reads the series file at `--stride` (default 2) for full
texture, one frame per timepoint.

### Keys and code defaults

| key | code default | read by |
|---|---|---|
| `timelapse.track.export_stride` | 2 | track-export (`--stride`), track-strips with `--image file` |
| `timelapse.track.export_crop` | `germline` | track-export (`--crop`) |
| `timelapse.track.export_max_gb` | 2.0 | track-export (a warning) |
| `timelapse.track.max_dist_um` | 3.0 | track-link laptrack (`--max-dist-um`) |
| `timelapse.track.gap_frames` | 2 | track-link laptrack (`--gap`) |
| `timelapse.track.min_overlap` | 0.4 | track-link overlap backends (`--min-overlap`) |
| `timelapse.track.link_stride` | 2 | track-link overlap backends (`--stride`) |
| `timelapse.track.strip_size_um` | 8.0 | track-strips (`--size-um`) |
| `timelapse.track.strip_planes` | 3 | track-strips (`--planes`) |
| `timelapse.track.strip_tile_px` | 64 | track-strips (tile side in output pixels) |
| `timelapse.track.registration` | `rigid` | track-link (`--registration`) |
| `timelapse.register.piecewise_segments` | 4 | register `--mode piecewise` (`--segments`) |
| `timelapse.register.piecewise_margin_um` | 6.0 | register `--mode piecewise` (`--margin-um`) |
| `timelapse.register.piecewise_taper` | 0.25 | register `--mode piecewise` (`--taper`) |

A command line flag wins over the key; the key wins over the code default. No shipped config file names
any of them, so every config hash is unchanged.

From code: `germquant.timelapse.export_cellacdc`, `link_series` (with `registration=`), `import_tracks`,
`track_strips`, `register_series_piecewise` and `load_piecewise`, and `find_series` / `load_shifts` /
`object_table` in `germquant.timelapse.layout` for reading a series folder and its registered coordinates.

## Validating the linking

Build the hand set before tuning anything, or it stops being blind (design section 10).

### 1. Pick the nuclei: `scripts/pick_validation_nuclei.py`

```
python scripts/pick_validation_nuclei.py results/SERIES [--seed 0] [--n 16 --tz 6 --pach 7 --edge 3] [--out picks.csv]
```

chooses, at the first timepoint of the series, 16 germline nuclei (`timelapse.validate.n_nuclei`): 6 in
the transition zone (`n_transition`), 7 in pachytene (`n_pachytene`) and 3 at the edge of the field
(`n_edge`), with a seeded generator (`seed`, 0), so the same call picks the same nuclei. Zones come from the
staging stage's `zone` column when it ran (`pre` is the transition zone, `early` / `mid` / `late` pachytene);
without a trace they come from `axis_position_um` quantiles (transition zone up to the 0.35 quantile,
`tz_quantile`; pachytene up to 0.90, `pach_quantile`; `centroid_x_um` orders the nuclei when there is no
axis). Edge nuclei lie within `edge_um` (5 um) of the xy border of the field; the most proximal nuclei
stand in when there are too few. The picks alternate `holdout` and `constraint`: tune the linking on the
constraint nuclei, judge it on the holdout ones, which nobody looks at until the verdict.

The worksheet `<stem>__validation_picks.csv` is written into the series folder (`--out` elsewhere): a
comment line, then one row per pick with `hand_track_id`, `role`, `zone`, `t_index` (0), `z_um`, `y_um`,
`x_um` (raw-frame microns), `nucleus_id`, `image_id`, `axis_position_um`. The biologist opens the series in
Imaris, finds each nucleus at t0 from its coordinates (check one known nucleus first: Imaris positions are
microns from the dataset origin, which Bio-Formats puts at 0 for an .nd2), follows it through the series
and appends one row per later timepoint with `hand_track_id`, `t_index`, `z_um`, `y_um`, `x_um`; the
other columns stay empty. A nucleus that leaves the field or cannot be followed simply stops. The file may
also be an Imaris points export renamed to those five columns.

### 2. Score the linking: `germquant track-validate`

```
germquant track-validate results/SERIES --points <stem>__validation_picks.csv [--snap-um 0.3] [--set all|constraint|holdout]
```

Every hand point is looked up in that timepoint's label image. A point outside every label is a
**segmentation miss**, not a linking error (the segmentation lost the nucleus, or the biologist's point
missed it; `--snap-um` lets a point just outside take the nearest nucleus within that distance). The
in-label points of a hand track are mapped to `time_track_id` through `<stem>__track_points.csv` (so the
verdict is about the last `track-link` or `track-import`), and consecutive in-label points form the hand
edges. Per hand track:

* `n_points`, `n_in_label`, `segmentation_recall` = n_in_label / n_points;
* `linking_recall` = hand edges the pipeline kept in one track / hand edges. Because the edges join
  consecutive in-label points, a dropout the pipeline bridged counts as kept, and a track that fragmented at
  a dropout counts against it;
* `edge_precision` = pipeline edges leaving a hand point that land on the hand nucleus at the next
  timepoint the track has / pipeline edges whose target timepoint has a hand point in a label;
* `switches` = hand edges where the pipeline's identity jumped: the track at the earlier point continues on
  a different nucleus, or the track at the later point came from a different nucleus (an identity swap);
  `edges_missing` = hand edges the pipeline broke without another nucleus involved (a fragmented or
  unlinked track). Both lower the linking recall; only switches are the identity errors the acceptance is
  about;
* `role` (holdout | constraint) and `zone`, from the worksheet, else the picks file in the series folder,
  else holdout.

Per point, `outcome` is `linked`, `segmentation_miss`, `switch`, `edge_missing` or `unlinked` (in a label
that no track covers). The verdict, printed and written to `<stem>__track_validation.json`:

| criterion | threshold (key, code default) | judged on |
|---|---|---|
| switches | `timelapse.validate.max_switches_holdout`, 0 | the holdout tracks (the tracks scored when none is holdout) |
| linking recall | `timelapse.validate.min_linking_recall`, 0.95 | pooled over the holdout tracks (all tracks when none is marked holdout) |
| segmentation recall | `timelapse.validate.min_segmentation_recall`, 0.85 | the pachytene zone when the hand set has zones, else all points |

`--set constraint` scores only the constraint tracks (`--set holdout` only the holdout ones; `all`, the default,
every track): that is how linking settings (a backend, the rigid or the piecewise registration, a cutoff) are
compared while the holdout tracks stay unseen; the chosen setting is then judged once with `--set all` (or
`holdout`). A subset run judges the switches and the linking recall on the tracks it scores, prints which set it
scored, and writes only `<stem>__track_validation_<set>.json`; the `track_validation` table is written by a run
over every track.

The exit code is 0 on PASS and 2 on FAIL (an error is 1). `<stem>__track_validation.csv` holds one row per
hand point (schema `track_validation`: the point, the label it landed in and its track, the outcome, the
validation set and zone) with the per-track scores repeated on every row of the track;
`drop_duplicates("hand_track_id")` gives the per-track table, which the JSON also carries. The design's
acceptance for tuning: zero switches on the holdout nuclei, linking recall at or above 0.95, pachytene
segmentation recall at or above 0.85; when the overlap backend fails only on dropouts, `laptrack` with gap
closing is the next thing to try, and the strips show what happened at any flagged point.

tests/test_timelapse_tracking.py scores a perfect hand set (every metric 1, PASS), a set where two hand
tracks exchange nuclei half way (one switch each, recall 3/4, FAIL) and a set with a point in the background
(a segmentation miss that leaves the linking recall at 1), checks that a fragmented track is a missing
edge, not a switch, and that `--set constraint` scores the constraint tracks alone (their switch and recall,
the subset record, the refusal of a set the hand file does not hold).

## Segmenting without DAPI: the training loop

A live worm carries no DNA stain. The October series has two channels, transmitted light and GFP::SYP-2,
so the nuclei have to be segmented on the synaptonemal complex signal itself, and the shipped nucleus model
(trained on DAPI) fragments SC bundles into several objects. The loop below gets from that first run to a
model trained on the series' own corrected labels, with provenance at every step. Everything in it is CPU
tested on a synthetic series (`tests/test_live_profile.py`, `tests/test_training_loop.py`); the real training
run on the GFP data has not been done yet.

### The live channel map and profile

`config/channel_maps/live_gfp_syp2.yaml` has two roles and no `dna`:

| role | names (substring, case-insensitive) | fallback | required | consumed by |
|---|---|---|---|---|
| `central_element` | EGFP, GFP, 488, 477 | channel 1 (flagged) | yes | the segment stage (see the channel key), germline isolation, the signal QC, the series reference volume |
| `brightfield` | Brightfield, BF, TD, Trans | none: by name only | no | nothing; recorded by the signal QC (`sig_brightfield_*`), given per-nucleus intensity columns and a montage panel, and kept as the registration cross-check volume |

`config/profiles/live_syp2.yaml` extends `config/config.yaml` with that map and sets
`segmentation.nuclei.channel: central_element`, switches spots, granule, coloc, envelope, partition and
staging off, keeps the signal QC on with its flags off, and sets `timelapse.keep_reference_volume: true`.
Germline isolation and the axis stay on: the GFP channel is the SYP marker the germline stage seeds on.
The nuclei table of a run made with it carries `central_element_*` and `brightfield_*` intensity columns
and the three `dna_*_intensity` columns as NaN (a map without a `dna` role writes them empty rather than
leaving them out, so a reader that expects them finds them); `dapi_present` is False. The stages that need
RAD-51, PGL-1 or lamin are recorded as `skipped` / `disabled` in `<id>__stages.json`.

```
germquant run SERIES.nd2 --profile live_syp2 --out results/SERIES --t 0 --xy-stride 4    # check the roles
germquant run SERIES.nd2 --profile live_syp2 --out results/SERIES --t all
```

### The channel key

`segmentation.nuclei.channel` names the channel ROLE the segment stage segments. Every shipped config
carries `channel: dna` and the key is read with the code default `dna`, so every frozen recipe resolves
exactly the pixels it always did and the goldens do not move. When the named role does not resolve in an
image the stage falls back to channel 0 and writes the flag `segment:no_<role>_channel_using_first_channel`
(for `dna` that is the text the fixed path always wrote). The montage draws the nucleus outlines over the
segmented channel. The key is part of the `segmentation` section, so a profile that changes it gets a new
segment sub-hash and a new config hash. `docs/BACKENDS.md` has the backend side of it: the first argument a
backend receives is this channel, whatever its role.

### 1. Export slices to correct

```
germquant export-training results/SERIES --out training/live_syp2 --orientations xy,xz,yz --slices-per-image 3 --seed 0
```

`export-training` takes one results folder, a series folder (every `t<NNN>/` folder is a run) or a batch
tree, and for every finished run (a folder with `run_manifest.json`, `<id>__nuclei_labels.tif` and
`<id>__image_summary.csv`):

* re-reads the stack from the source file named in the run's tables, at the `xy_stride`, `z_range` and
  `t_index` recorded in `run_manifest.json`, so the slices sit on exactly the grid the labels were made on;
* takes the channel the run segmented (`segmentation.nuclei.channel` of the recorded config, resolved
  through the run's `channel_map_json`);
* for each orientation draws `--slices-per-image` slices at random (seeded per image, so adding a run to
  the folder does not change the slices of the others) from the slices that cut at least `--min-labels`
  nuclei (default 3; when fewer qualify, from every slice that cuts one);
* writes `image/<id>_<orient>_<index>.tif` (uint16, the raw intensities) and `labels/<id>_<orient>_<index>.tif`
  (int32, the run's labels). `xz` and `yz` slices are resampled along z to the xy pixel size (linear for
  the image, nearest neighbour for the labels), so every exported slice is isotropic. A 2D model trained on
  isotropic slices applies in 3D with `anisotropy = dz / dy`, which is how the segment stage runs Cellpose,
  and the xz and yz slices teach it the z profile of a nucleus that xy slices alone cannot;
* records every pair in `manifest.csv`: `image_id`, `series_id` (the id without its `_t<NNN>` suffix),
  `t_index`, `orientation`, `slice_index`, `n_labels`, `zone`, the two file paths, the source file, the
  geometry, the channel role and index, `pixel_um`, `z_zoom`, `height`, `width`; and the settings in
  `export.json`. `--zone-column zone` (the default) tags each slice with the most frequent zone among its
  nuclei when the nuclei table has that column (the staging stage's); without it the tag is empty and the
  biologist adds zones by hand in the manifest. Running the export again into the same folder replaces the
  rows (and files) of the images exported, and keeps the rest.
* stops at `--max-gb` (default 2) of written slices: a 56 GB series must not become a 56 GB training folder.

Include the last timepoint (t24 in the October series): bleaching makes it the hardest frame.

### 2. Correct the labels

Open a pair in napari (`napari training/live_syp2/image/X.tif training/live_syp2/labels/X.tif`, the second
layer as labels) or in the Cellpose GUI, fix the merges and splits, and save the labels back under the same
name in `labels/`. A pair that is not worth correcting is deleted (both files); `train-nuclei` skips the
manifest rows whose files are gone. About 30 corrected slices across zones and timepoints is the design's
first training set; the zone tag in `manifest.csv` is what lets the hold-out metrics be reported per zone.

### 3. Train

```
germquant train-nuclei --pairs training/live_syp2 --out models/live_syp2_nuclei --epochs 150 --split-by series --val-frac 0.2 --seed 0
```

`train-nuclei` runs `cellpose.train.train_seg` in process, as `scripts/train_nucleus_model_combined.py`
did (AdamW, learning rate 1e-5, weight decay 0.1, batch 8, 150 epochs), on the GPU when one is visible
(`--cpu` keeps it off every accelerator). The hold-out set is whole groups, never slices: `--split-by series`
(the default) holds out whole series, so two neighbouring timepoints of one worm never straddle the split;
`--split-by image` holds out whole images (for fixed images the two are the same). At least one group is
held out whenever there are two or more; a single group trains alone and the card says there are no
hold-out metrics. `--base` is the starting weights: `cpsam` (the default, Cellpose's stock model), a
registered model name (`germline_nuclei_combined` once `get-model` has fetched it), a model file, or
`none` for a randomly initialised network (smoke tests only). Cellpose drops training slices with fewer
than `--min-train-masks` nuclei (default 5); the command warns about them beforehand.

It writes, in `--out`:

* the model file `<name>` (`--name`, default the folder's name), float32 Cellpose-SAM weights, 1.2 GB;
* `training_log.json`: the losses per epoch, the split (groups and files on each side), counts per
  orientation, every setting, the device, start and end times;
* `<name>.card.json`, the model card: the base model (name, path, sha256), the sha256 of the model file,
  the training set with the sha256 of every image and label file, its image id, series, orientation, slice
  index, zone and split side, the orientations and zones, the epochs and the other settings, the date, the
  germquant, Cellpose and torch versions, and the hold-out metrics: `germquant.validate.segmentation_metrics`
  on every validation pair (F1 at IoU 0.5 from the pooled true and false positives and misses, precision,
  recall, the miss rate fn / n_gt, the mean matched IoU) overall, per orientation and per zone, plus the
  per-pair rows.

The design's gate, read off the card before the model is used for numbers: F1 at IoU 0.5 of 0.85 or better
on pachytene slices and 0.75 on transition-zone slices, and a miss rate under 10 percent.

### 4. Register and use

```
germquant register-model --model models/live_syp2_nuclei/live_syp2_nuclei --card models/live_syp2_nuclei/live_syp2_nuclei.card.json
germquant get-model live_syp2_nuclei          # verifies the registered file (sha256), downloads nothing
```

`register-model` adds the model to `local_models.json` next to the release assets (`models/models/` under
the checkout, or `$GERMQUANT_MODELS/models/`): name, absolute path, sha256, size, the card's path and
summary. The sha256 is computed from the file and must agree with the card's `model_sha256`, so a card
cannot be attached to other weights; a release name or a Cellpose built-in name is refused; registering a
name again replaces its entry. Then a profile names it:

```yaml
extends: live_syp2.yaml
segmentation:
  nuclei:
    cellpose_model: live_syp2_nuclei
```

`segmentation.nuclei.cellpose_model` accepts a registered name: a bare name that is neither a file next to
the config nor under `$GERMQUANT_MODELS` is looked up in the registry (local models, then downloaded
release assets) and resolves to the registered file; `cpsam` and the other Cellpose names still pass
through. The run's `run_manifest.json` records the file's digest as `model_sha256` and the segment stage's
sub-hash changes, so `--resume` reprocesses every timepoint made with the old model. Rerun the series with
the new profile, then `export-training` again if a second round of corrections is wanted: the loop is
cheap to repeat because the labels to correct are the model's own.

### What is tested, and what is not

`tests/test_live_profile.py` runs the synthetic two-channel series through the profile end to end on the
CPU (the channel key, the NaN dna columns, the signal QC of both roles, the planted drift in the nucleus
centroids). `tests/test_training_loop.py` covers the export (file names, dtypes, the isotropic resampling,
the manifest, the zone tag, the seeded draw, the size cap), the split rule, the training wrapper with the
Cellpose calls replaced by fakes (the losses, the sha256 of every pair, the per-orientation F1, the card),
and the registry (register, verify, resolve into the segment stage). A real training run of a randomly
initialised network for two epochs on the CPU is opt-in, `GERMQUANT_TRAIN_SMOKE=1`, because it takes about
90 seconds and writes 1.2 GB. The GFP model itself, its gate numbers and the rerun of the October series
are the work of the first real loop.

## Registration

Drift between timepoints is measured, never corrected in the images: `germquant register` writes a shift
table, and the tracking package subtracts the shifts from the nucleus centroids. Nothing the stages see is
resampled, so every per-timepoint table stays what the fixed-image pipeline would have written.

### The reference volumes

`germquant run --t` keeps, next to each timepoint's fixed-image file set, a downsampled copy of the
segmentation channel: `t<NNN>/<image_id>__ref_stride<k>.npy`, a uint16 array of shape Z x Y/k x X/k (every
plane, every k-th row and column of the stack as read, so a run with `--xy-stride 2` and k = 8 keeps every
16th pixel of the file; a float stack is rounded and clipped). When the channel map resolves a role named
`brightfield`, that channel is kept the same way as `<image_id>__ref_brightfield_stride<k>.npy`. Two config
keys, both read with code defaults, so no config file changes:

| key | code default | meaning |
|---|---|---|
| `timelapse.reference_stride` | 8 | k, the xy stride of the reference volumes (4 for a finer estimate at four times the disk space) |
| `timelapse.keep_reference_volume` | true | whether a series run keeps them; a fixed image (`run` without `--t`) never writes one |

The volumes come out of the stack `process_image` has already read (its `on_stack` callback), so the 56 GB
file is not read twice; at stride 8 the October series costs 18 MB per timepoint per channel. The
segmentation channel is the role `segmentation.nuclei.channel` names (code default `dna`), channel 0 when
that role does not resolve (the segment stage's own fallback). Each timepoint's entry in
`<stem>__series.json` records them under `reference_volumes`: `stride`, `role`, `channel_index`,
`role_resolved`, `spacing_um` (dz, dy x k, dx x k), `shape`, `dtype` and `files` (role to path relative to
the series folder); the record's top level carries `keep_reference_volume` and `reference_stride`. A
failure to write a volume is logged and recorded in that entry's `error`; the timepoint itself goes on.

A series processed before the volumes existed (or with the key off) is not lost: `germquant run SERIES.nd2
--t all --out ... --resume` with the key on skips every finished timepoint as before and, when a resumed
timepoint lacks its volumes, reads that timepoint back from the file at the run's geometry and writes them
(no stage runs; `regenerated: true` in the entry). On the real series that is a read of about 25 s per
timepoint over the network instead of a night of segmentation.

### The command

```
germquant register results/SERIES
germquant register results/SERIES --reference t0 --upsample 20 --max-step-um 2 --channel brightfield
germquant register results/SERIES --rigid
```

`SERIES_DIR` is the `--out` folder of the series run. The command needs no config: it reads the series
record (or, without one, the `t<NNN>/` folders and their file sets), the reference volumes, the per-timepoint
`__nuclei.csv` tables and the voxel size of the volumes (from the record entry, else the `__image_summary.csv`
voxel columns times the stride).

For each timepoint, in `t_index` order, the translation between its reference volume and the reference
frame's is estimated with `skimage.registration.phase_cross_correlation` at `--upsample` 10 (the code
default): the estimate is quantised to a tenth of a reference-volume voxel, 0.086 um at stride 8 on the
October series, below the full-resolution pixel. The reference frame is the previous registrable timepoint
(`--reference previous`, the default; the steps are accumulated to the first timepoint with a volume), or the
first timepoint itself (`--reference t0`: every estimate is independent, the step is the difference of
consecutive estimates, and an estimate that went wrong does not propagate). A step longer than
`--max-step-um` (code default 5.0 um; 0 switches the cap off) is shortened to that length, its direction kept,
and flagged `capped`; in `previous` mode a capped step stays in the cumulative shift of every later timepoint,
which is the cost of chaining and is why the quality columns and the brightfield cross-check exist.

`--channel ROLE` registers another kept volume (the segmentation channel's by default; `brightfield` when one
was kept; the other kept volume becomes the cross-check). `--normalization none|phase`: by default the
cross-power spectrum is not whitened (plain cross-correlation with skimage's upsampled peak refinement; the
`method` column says `cross_correlation`). Nuclei at stride 4 to 8 are smooth blobs with little texture, and
whitening (`phase`, the textbook phase correlation, method `phase_cross_correlation`) then amplifies the
empty high frequencies: on the synthetic series of the tests it biased a step by up to 0.6 reference voxel
where the plain correlation was exact to the quantisation, with and without camera noise. `phase` stays
available for textured volumes. `--taper FRACTION` (0, off, by default) subtracts each volume's mean and
multiplies it by a Tukey window that tapers that fraction of every axis (0.25: the outer eighth on each
side) before the correlation; see the caveats below for when to use it. `--rigid` is described below.

Sign convention, shared by every time-lapse package: a positive shift means the sample moved by that much in
the raw frame (a feature at position p at t0 sits at p + shift at t), so a registered coordinate is the raw
centroid minus the cumulative shift of its timepoint. `phase_cross_correlation` reports the opposite (the
shift that brings the moving volume back onto the reference) and the command negates it. Microns, in the raw
frame, (z, y, x) order.

### The outputs

`<stem>__register_shifts.csv` and `.parquet` (`schema.SERIES_TABLES["register_shifts"]`), one row per
timepoint, the declared columns first, then `SERIES_META` (the design columns from the series summary, the
run's `config_hash`, this registration's version, git sha and timestamp):

| column | meaning |
|---|---|
| `t_index`, `t_seconds`, `image_id` | the timepoint (`t_seconds` from the series record, NaN without one) |
| `reference_t_index` | the frame this timepoint was registered against (empty for the origin) |
| `channel_role` | the role whose volume was registered (`segmentation` when no record names it) |
| `dz_um`, `dy_um`, `dx_um` | the cumulative shift to the first timepoint; the tracking package subtracts it. `cum_shift_z_um`, `cum_shift_y_um`, `cum_shift_x_um` carry the same numbers (the names declared first) |
| `step_dz_um`, `step_dy_um`, `step_dx_um` | the step against the reference frame, after the cap; `shift_z_um`, `shift_y_um`, `shift_x_um` are the same numbers |
| `capped` | the step hit `--max-step-um` |
| `bf_shift_z_um`, `bf_shift_y_um`, `bf_shift_x_um` | the brightfield (cross-check) step, never capped |
| `brightfield_dz_um`, `brightfield_dy_um`, `brightfield_dx_um` | the brightfield cumulative shift |
| `brightfield_disagreement_um`, `bf_disagreement_um` | the distance between the brightfield step and the uncapped primary step (the same number twice); NaN without a brightfield volume |
| `disagreement_flag` | that distance is above 1 um (`register_series(disagreement_um=...)` from code) |
| `quality_before_um`, `quality_after_um` | the median nearest-neighbour displacement of the nucleus centroids (germline nuclei when `in_germline` marks any) between this timepoint and its reference frame, raw and after the step is subtracted; `nn_disp_before_um`, `nn_disp_after_um` are the same numbers. NaN without nuclei tables |
| `error` | why a timepoint has no estimate (`no reference volume`, a shape mismatch after a stride change, ...); empty otherwise. The next registrable timepoint is registered against the last one that was, so one bad frame costs one row, not the chain |
| `xy_stride`, `upsample_factor`, `method` | the reference volumes' stride, the upsampling, `cross_correlation` or `phase_cross_correlation` |

`<stem>__register_report.json`: the parameters (reference mode, upsample, normalization, cap, disagreement
threshold, channel and cross-check roles, stride, voxel size, origin timepoint), per-step summaries (step, raw
step when capped, cumulative, brightfield step and cumulative, disagreement, quality before and after), the
lists of capped and disagreeing timepoints, the errors, the median quality before and after, the largest
step, the cumulative shift at the last timepoint, and the rigid fits when asked for.

The quality numbers are a sanity check, not a loss: after a good registration what remains is the real flow
of the nuclei plus segmentation jitter, so `quality_after_um` should drop well below `quality_before_um` and
stay well under the nearest-neighbour spacing (3.3 to 4.7 um in the goldens). On the lab series the flow
itself is a few voxels per 5 minutes, about the size of the drift; the brightfield, which drifts but does not
flow, is what tells them apart when it disagrees.

Two caveats. A feature the field edge cuts (tissue leaving the frame, a bright object at the border) does
not move between frames and pulls a correlation towards zero shift; the test gonad of
tests/test_timelapse_series.py has such a blob on purpose and registers a quarter short in x against t0 and
half short when two such steps are chained, while the sub-voxel series without one registers to the
quantisation. On a field as small as that toy (80 x 20 reference voxels) no remedy works, because a window
skews the content it keeps; on the 340 x 340 reference volumes of the lab series `--taper 0.25` is the
standard remedy and costs nothing when no edge is cut (on the sub-voxel series it changes the recovery by
at most a tenth of a reference voxel over four steps). When `quality_after_um` does not drop well below
`quality_before_um`, or the brightfield disagrees, try `--taper 0.25` and `--reference t0` before anything
else. And a uniform flow of every nucleus is indistinguishable from drift in the GFP channel; the
brightfield cross-check is the only handle on it.

### Level 2: a rigid fit (experimental)

`--rigid` adds, per step, a Kabsch rigid fit on nucleus centroids: the centroids of the timepoint are
translated by the step, linked to the reference frame's by mutual nearest neighbour within 0.7 times the
reference frame's median nearest-neighbour spacing, and the proper rotation plus residual translation that
maps the linked pairs onto each other is stored in the report (`rigid`: `n_linked`, `link_radius_um`,
`rotation_deg`, `rotation_axis` in (z, y, x) order, `rotation_matrix`, `translation_um`, `rmsd_before_um`,
`rmsd_after_um`; the fit fields are null with fewer than six links). The shift table is unchanged: nothing
consumes the rotation yet, it is the diagnostic the design asked for before deciding whether a rotation stage
is needed (DESIGN_timelapse.md section 4). A pure translation fits to about 0 degrees; a planted 3 degree
turn per step is recovered within 0.5 degree on the synthetic series.

### Piecewise registration (`--mode piecewise`)

```
germquant register results/SERIES --mode piecewise [--segments 4] [--margin-um 6] [--taper 0.25] [--channel ROLE] [--max-step-um 5] [--upsample 10]
germquant track-link results/SERIES --backend laptrack --registration piecewise
```

Why: on the first real series the worm does not drift (the brightfield steps are 0.17 um or less), but the
germline slides inside it non-rigidly, the proximal half up to 3.6 um per 5 minutes while the distal end
barely moves. One translation per frame cannot describe that; registering on GFP breaks the distal links, on
brightfield the proximal ones (`C:/Users/ryane/live_results/FIRST_LOOK.md`, section 7).

What it does (`germquant.timelapse.piecewise`):

1. At the origin (the first timepoint with a reference volume) the germline nuclei are placed on the fixed axis
   of the tracks ("One reference axis": the origin's `axis_position_um` smoothed by a local linear fit over the
   16 nearest nuclei; a fitted centerline when the axis stage did not run) and split into K segments of equal
   nucleus count (K = `--segments`, key `timelapse.register.piecewise_segments`, code default 4; at least 5
   nuclei per segment). A segment's centre is the median axis position of its nuclei.
2. For every later timepoint and every segment, the translation between the previous registrable frame and this
   one is estimated exactly as a rigid step is (the same `displacement`: plain cross-correlation by default,
   upsampled peak refinement to 1/`--upsample` of a reference voxel, the same `--max-step-um` cap with the
   direction kept and the `capped` flag), on a crop of the two reference volumes: the y and x extent of the
   segment's origin nuclei moved to where the segment sat in the previous frame, widened by `--margin-um` (key
   `timelapse.register.piecewise_margin_um`, 6.0 um), every plane in z. The crop is mean-subtracted and tapered
   with a Tukey window (`--taper`, key `timelapse.register.piecewise_taper`, 0.25), because a crop always cuts
   the gonad where the segment ends and an untapered cut edge pulls the estimate towards zero. The reference
   frame is always the previous one (`--reference t0` is refused). The steps are accumulated per segment.
3. A segment whose own estimate fails (a crop that left the field, a flat volume) takes the step interpolated
   from the other segments and is flagged `filled` with the reason in `error`; a timepoint where every segment
   fails is a set of rows with an error, and the chain continues from the last registrable frame.
4. The other kept volume (brightfield when GFP is registered) is registered on the same crops as the
   cross-check, never capped: `brightfield_step_d*_um`, `brightfield_d*_um`, `brightfield_disagreement_um` and
   `disagreement_flag` (above 1 um) per segment.
5. A nucleus at axis position a gets the displacement D(t, a) interpolated linearly between the segment centres
   (constant beyond the outermost centres), and its registered coordinate is its raw centroid minus D(t, a).
   Because a is the nucleus's position on the origin's axis after registration, it is found by a fixed-point
   iteration (it converges in a few steps; the displacement varies slowly along the axis). Nothing is resampled.

`<stem>__register_piecewise.csv` (and `.parquet`), one row per timepoint per segment, the declared columns of
`schema.REGISTER_PIECEWISE` then `SERIES_META`: `t_index`, `t_seconds`, `image_id`, `reference_t_index`,
`origin_t_index`, `channel_role`, `segment`, `n_segments`, `axis_start_um`, `axis_end_um`, `axis_centre_um`,
`n_nuclei_t0`, the cumulative displacement `dz_um`, `dy_um`, `dx_um` (positive = the segment moved by that
much in the raw frame), the step `step_dz_um`, `step_dy_um`, `step_dx_um`, `capped`, `filled`, the brightfield
columns, `n_nuclei` and `quality_before_um` / `quality_after_um` (the median nearest-neighbour displacement of
the segment's nuclei at this timepoint against the reference frame's, raw and after the registration of both
frames), the crop box `crop_y0`, `crop_y1`, `crop_x0`, `crop_x1` in reference voxels, `error`, `xy_stride`,
`upsample_factor`, `method`, `taper`, `margin_um`. `<stem>__register_piecewise_report.json` holds the
parameters, the segments, the per-step and per-segment summaries, the capped, disagreeing and filled segment
steps, the errors, the quality medians per segment and the cumulative displacement at the last timepoint. The
rigid `<stem>__register_shifts.csv` is not touched; `track-link --registration piecewise` reads the piecewise
table (`load_piecewise`), rebuilding the fixed axis from the origin's nuclei table exactly as the registration
built it.

What the registered frame means here: it follows the germline itself, segment by segment. Linking in it sees
nearly stationary nuclei, which is the point, but the tracks' speeds and `axis_velocity_um_per_min` are then the
motion relative to the local flow. The flow is the registration (the cumulative displacement per segment) or
`axis_velocity_raw_um_per_min` ("Which velocity is the flow" above).

Two limits, both seen on the October series (FIRST_LOOK.md section 7): the reference volumes at the default
`timelapse.reference_stride` 8 are coarse for a segment crop (on a stride-4 run a reference voxel is 3.47 um,
about one nucleus, and the segment estimates jump by whole voxels), so keep `timelapse.reference_stride: 2` for
a series meant for the piecewise registration; and the brightfield cross-check is not a test of the segment
estimate, because the gonad has texture in transmitted light too, so the brightfield crop over a sliding
segment partly follows the germline and partly the still body wall.

### Tests

tests/test_timelapse_register.py builds a synthetic series of 16 gaussian nuclei and a brightfield channel
with a planted sub-voxel drift per step (0.4 planes, 1.3 and 2.7 voxels), a zero-mean shear flow of the
nuclei, and writes what a series run leaves behind. It checks that the cumulative shift is recovered within
0.2 reference voxel at stride 4 and 0.5 at stride 8 (the tolerances are in reference-volume voxels, the unit
phase correlation works in), against the previous frame and against t0; that a planted 1.6 um jump is capped
to a 1 um step and flagged while the brightfield, which jumped too, raises no disagreement; that a jump in
the nuclei alone raises the disagreement flag at that step and only there; that the quality drops after
every step; that a series without reference volumes fails with a message naming
`timelapse.keep_reference_volume`; that a missing volume costs one row and the chain continues; that
`--rigid` finds no rotation on a translation and recovers a planted 3 degree turn; and, on the
three-timepoint gonad of test_timelapse_series.py, that `germquant run --t all` writes the volumes next to
the fixed-image file set and records them, that a fixed image writes none, that `--resume` reads missing
volumes back without running a stage, and that `germquant register` runs on the loop's output.

tests/test_timelapse_piecewise.py plants the non-rigid flow of the October series on a synthetic gonad: 60
nuclei in two jittered rows along x, a brightfield body that does not move, and a shear flow along the axis
(each nucleus moves u = U s / L per timepoint, 0 at the distal end and U = 4.6 um at the proximal end, more than
half the 4 um spacing). It checks the segment layout (equal counts, centres), the interpolation (linear between
centres, constant beyond) and the fixed-point registration on a straight axis; that each segment's step is the
shear at its centre within 0.3 um and the cumulative displacement within 0.3 um per step; that the brightfield
cross-check reads zero and flags exactly the sliding segments; that the per-segment residual after
registration is under 0.35 um and under half the raw one; that the interpolated displacement of every nucleus
follows the planted flow (within 0.3 um per step between the outermost centres, within the clamp error beyond
them); that laptrack after the rigid registration joins different nuclei (at least 8 wrong links; 75 on this
series) while after the piecewise registration all 60 tracks are full length with no wrong link, for laptrack
and for the overlap rule on per-nucleus shifted labels; that the raw-frame axis velocity of every track equals
its planted flow while the registered-frame one stays small; the command line (`--mode piecewise`, the
options it refuses, the segment count from the series' run_manifest.json and from `--segments`,
`--registration piecewise` without the table exits 1, the default stays rigid); the missing-volume and
too-few-nuclei errors; and that `register_piecewise` is appended last to `schema.SERIES_TABLES`.

## What stays unchanged

The fixed-image path is byte-identical: `process_image` without `t` makes the same reader call, the
same image id, the same tables and the same marker keys as before; `tests/test_golden.py` runs the
synthetic goldens unchanged, and `config/config.yaml`, the study configs and the channel maps are not
edited. Every new key is read with a code default; the new schema entries are appended, not inserted.
`segmentation.nuclei.channel` is read with the code default `dna`, the value every shipped config carries.
