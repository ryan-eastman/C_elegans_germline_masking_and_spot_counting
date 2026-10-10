# Time-lapse quantification

Status: the reader ("Reading one timepoint"), the loop over a series with its resume logic and command line
switch ("Running a series") and the drift registration ("Registration") are implemented; tracking and the
biologist validation are the next packages and are described in [DESIGN_timelapse.md](DESIGN_timelapse.md).
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

`schema.SERIES_TABLES` declares five tables, one file per series, written next to the timepoint folders
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
| `register_shifts` | `germquant register` (timelapse.register; see "Registration") | timepoint | the cumulative shift to the first timepoint (`dz_um`, `dy_um`, `dx_um`) and the step against the reference frame in microns, the brightfield cross-check and its disagreement flag, whether the step was capped, the median nearest-neighbour displacement before and after the step, an error text, stride, upsample factor, method |
| `track_points` | tracking package | object per timepoint | `time_track_id`, `t_index`, `t_seconds`, `image_id`, `nucleus_id`, raw and registered centroids, the step from the previous point, `is_gap`, `link_source`; nucleus measurements are joined on `(image_id, nucleus_id)` |
| `tracks` | tracking package | track | `n_points`, `n_gaps`, `t_start`, `t_end`, `is_full_length`, path length, net displacement, speed, signed axis velocity after the t0 sign alignment, SYP intensity slopes raw and median-normalised, SC length change, `link_source` |
| `track_validation` | `germquant track-validate` | hand-followed point | the biologist's point in raw-frame microns, the label it lands in and its track, the outcome (linked, segmentation_miss, switch, edge_missing), the validation set (holdout or constraint) and zone |

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

## What stays unchanged

The fixed-image path is byte-identical: `process_image` without `t` makes the same reader call, the
same image id, the same tables and the same marker keys as before; `tests/test_golden.py` runs the
synthetic goldens unchanged, and `config/config.yaml`, the study configs and the channel maps are not
edited. Every new key is read with a code default; the new schema entries are appended, not inserted.
