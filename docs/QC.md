# Quality control

What the pipeline records about each image beyond the biology, and the checks it runs on those
records. Every check writes a flag string into the `qc_flags` column of `image_summary` (also in
`batch_summary.csv` and in the per-stage record `<image_id>__stages.json`). A flag is advisory
unless it names a missing required channel, a role collision or an image with no nuclei: only those
set `qc_pass` to false. The flags described here never fire on a run that does not switch them on,
so the `qc_flags` of every existing run and of the synthetic goldens (`tests/test_golden.py`) are
unchanged.

## The few-nuclei flag

`qc:few_nuclei_<n>` is raised, as an advisory flag, when an image has between 1 and 9 segmented
objects; `qc:NO_NUCLEI` (a QC fail) when it has none. A recipe can move the threshold with
`qc.min_nuclei` (code default 10, read only when a recipe sets it): a field of cultured cells holds a
handful of nuclei, so the cultured-cell pack sets 1, and the flag would otherwise mark almost every
image. No shipped config sets the key, so the C. elegans runs keep the threshold of 10
(`docs/PLANAR.md`, "Recipe keys for packs"). The value must be a whole number, 0 or more; anything else
(2.7, -1, `true`, a quoted "10") stops the image before any stage runs, with the key named.

## Acquisition record

### Where the record comes from

A Nikon `.nd2` carries a text description (`nd2.ND2File(path).text_info["description"]`) with one
`Plane #n:` block per channel. Each block names the channel (`Name: 545`) and lists the camera
settings, the X-Light filter wheels and the Celesta laser line that was on. The `acquisition`
stage (config key `acquisition.read_exposures`, off by default, on in the `ccw77_partition`
profile) parses every block into one record per channel with these keys:

| key | from the description | type |
|---|---|---|
| `exposure_ms` | `Exposure: 200 ms` (an exposure given in seconds is converted) | number |
| `power` | `Line:5; ExW:545; Power: 96.0; On` (the line named like the channel, or the only line listed) | number |
| `em_filter` | `X-Light, FilterChanger(EM Wheel): 3 (595/31)` | string, the label `595/31` |
| `ex_filter` | `X-Light, FilterChanger(EX Wheel): 1 (MXR00543 - ... (Ex))` | string, the label |
| `dichroic` | `X-Light, FilterChanger(DM Wheel): 1 (MXR00543 - ... (Di))` | string, the label |
| `camera_mode` | `Camera Mode: 12-bit (CMS)` | string |
| `bit_depth` | `Bit Depth: 12-bit` | integer |
| `binning` | `Binning: 1x1` | string |
| `conversion_gain` | `Conversion Gain: Sensitivity` | string |
| `readout_speed` | `Readout Speed: 100.0MHz` | string |
| `modality` | `Modality: Widefield Fluorescence, Spinning Disk Confocal` | string |

A line the description does not carry gives a null value. The filter keys hold the filter label
inside the parentheses, not the wheel position in front of it: `3 (595/31)` is recorded as
`595/31`, because the label is what a person recognises and what an expected-settings entry names.
Only an `.nd2` has this description; for TIFF, CZI, LIF and OME-Zarr input the stage is skipped
with the reason `not an .nd2`.

### What is written

The record is mapped onto the channel roles of the channel map by channel name (an exact match,
then a substring match either way, so a channel called `Widefield 545` still finds the `545`
plane). `image_summary` gets, per resolved role:

| column | meaning |
|---|---|
| `exp_ms_<role>` | exposure in ms (the original column of the stage, unchanged) |
| `laser_pow_<role>` | laser power as the software reports it, percent (unchanged) |
| `em_filter_<role>` | emission filter label |
| `camera_mode_<role>` | camera mode |
| `bit_depth_<role>` | bit depth |

and one fixed column, `acquisition_json`, holding the complete record of every channel in the file
(channel names as keys, every key of the table above) as compact JSON. The per-role columns take
their names from the channel map, so they are written after the declared columns like the other
per-role columns of the pipeline; `acquisition_json` is declared in `schema.py` and is therefore
present, empty, on every run that does not switch the stage on. A value the description lacks is an
empty cell (`NaN` in pandas, `NA` in R once the empty string is read as missing); a channel the
description has no plane for leaves every column of its role empty.

Reading the JSON back:

```python
import json, pandas as pd
isum = pd.read_csv("results/<image_id>__image_summary.csv")
record = json.loads(isum["acquisition_json"].iloc[0])      # {"405": {...}, "477": {...}, "545": {...}, "640": {...}}
record["545"]["em_filter"]                                  # "595/31"
```

```r
library(jsonlite)
isum <- read.csv("results/<image_id>__image_summary.csv", na.strings = c("", "NA"))
record <- fromJSON(isum$acquisition_json[1])
record[["545"]]$em_filter
```

### Looking before running

`germquant info FILE.nd2` prints the table for one file, no pixels loaded:

```
file       : C:\Users\ryane\ccw80_local\20260716_ccw80_I01F_488gtbp1_555syp3_647lmn1_HS_male018.nd2
sizes      : {'Z': 44, 'C': 4, 'Y': 2720, 'X': 2720}
dtype      : uint16   2D=False
voxel µm   : dz=0.2000 dy=0.1083 dx=0.1083  (anisotropy z/xy=1.85)
channels   : ['405', '477', '545', '640']
acquisition (from the .nd2 description; the full record goes to image_summary.acquisition_json):
  channel    exposure_ms  power  em_filter    camera_mode
  405                  7     49  438/24       12-bit (CMS)
  477                 50     91  511/20       12-bit (CMS)
  545                200     96  595/31       12-bit (CMS)
  640                 30     67  685/40       12-bit (CMS)
```

This is the fastest way to find the emission filters and exposures a session actually used before
writing them into a profile.

### Expected settings: the check

Config section `acquisition.expected`, read with a code default of "absent", meaning no check. It
is a mapping from channel role to the settings the study expects for that role:

| key | meaning | flag raised |
|---|---|---|
| `em_filter` | one accepted emission filter label, or a list of them | `acq:<role>_em_filter_mismatch` |
| `exposure_ms_max` | the longest exposure the study allows, ms | `acq:<role>_exposure_above_max` |
| `power_max` | the highest laser power the study allows | `acq:<role>_power_above_max` |
| `camera_mode` | one accepted camera mode, or a list | `acq:<role>_camera_mode_mismatch` |

String values are compared after removing whitespace and ignoring case (`12-bit (CMS)` equals
`12-BIT(cms)`). A value the description does not carry counts as a mismatch for the two string
checks, because the expectation could not be confirmed; the two maxima stay silent when the value
is absent, since nothing is above a maximum. Equal to the maximum is not above it. A role of the
section that did not resolve in an image (no lamin channel in a 3-colour file) is not checked. The
flags are advisory: `qc_pass` does not change, the image is processed as usual, and the flags sit
in `qc_flags` next to the per-role columns that show the offending values.

The section is validated before any pixel is read: a role the channel map does not define, a key
other than the four above or a value of the wrong type raises a `ValueError` naming the entry, so a
typo cannot silently check nothing. The check needs the stage: with `acquisition.read_exposures`
off (or a non-nd2 input) a warning says the check was skipped.

### Worked example: the ccw77 pack

The lab found that one session of ccw77 images had a dead mCherry::SYP-3 fluorophore imaged at 97
percent laser power and 300 ms, and that the 545 channel of that session was taken through the
wrong emission filter, so the 647 lamin stain could reach it: the partition coefficient floor came
out above 1 as a result. The expected settings of the pack (the ccw80 and good ccw77 sessions, read
with `germquant info`) are 438/24 for 405 (DAPI), 511/20 for 477 (PGL-1), 595/31 for 545 (SYP-3)
and 685/40 for 640 (lamin), all at `12-bit (CMS)`. A study profile that records and checks them:

```yaml
# config/profiles/<your study>.yaml: the August recipe plus the acquisition check. A new profile is a
# new recipe with its own config_hash; the shipped profiles are not edited.
extends: ccw77_partition.yaml

acquisition:
  read_exposures: true              # already on in ccw77_partition; stated here because the check needs it
  expected:
    dna:
      em_filter: 438/24
    granule:
      em_filter: 511/20
      exposure_ms_max: 100
    central_element:                # SYP-3 on the 545 line: the channel the bad session got wrong
      em_filter: 595/31
      exposure_ms_max: 250
      power_max: 96.5
      camera_mode: 12-bit (CMS)
    lamin:
      em_filter: [685/40, 700/75]   # a list accepts either filter
```

Run with `germquant batch E:/<folder> --profile <your study> --out results/`. An image from the bad
session then carries

```
acq:central_element_em_filter_mismatch;acq:central_element_exposure_above_max;acq:central_element_power_above_max
```

in `qc_flags`, with `em_filter_central_element = 685/40`, `exp_ms_central_element = 300` and
`laser_pow_central_element = 97` beside it, while an image from a good session carries none of the
`acq:` flags. The flags are the handle for the downstream analysis: filter `batch_image_summary.csv`
on `qc_flags` not containing `acq:central_element_` before pooling partition coefficients, or keep
the flagged images as a labelled group. Whether the lamin signal did reach the 545 channel in those
images is the question the signal QC below is for; this record only says the settings were wrong.

## Signal QC

Stage `signal_qc` (`src/germquant/signal_qc.py`), on by default (`signal_qc.enabled`), runs right
after `measure`, before any stage trusts a channel. `--no-signal_qc` switches it off for one run. It
needs no particular channel role: every role the channel map resolved is measured, and roles the
image does not carry leave their columns NaN. It adds only new columns and its own table, and its
flags are off unless a profile enables them, so the existing columns and the `qc_flags` text of
every run made before it existed are unchanged.

### What it answers

Two questions per image. Is each channel measurable, that is, does the tissue stand above the camera
noise and below saturation? And does one channel leak into another, that is, do two channels that
stain different structures nevertheless rise and fall together voxel by voxel in the region where
the leak would matter? The ccw77 case is the model for both: a 545 channel holding nothing but noise
and whatever the filter admitted from the 647 lamin stain fails the first test (its SNR is that of
pure noise) and the second (it correlates with the lamin channel in the perinuclear ring, exactly
where the P granules sit).

### Sampling and regions

Everything is computed on a deterministic xy sub-sample of the stack (every `signal_qc.stride`-th
row and column, every z plane), so a 2048 x 2048 x 40 stack costs seconds and the result is the same
on every run. The nucleus labels, sub-sampled the same way, define the regions through one
anisotropy-aware Euclidean distance transform (the voxel spacing of the sub-sampled grid is used, so
the dilations are in microns in every direction):

| region | definition | used for |
|---|---|---|
| tissue | labels dilated by `signal_qc.tissue_dilation_um` (code default 2.0 um) | the voxels OUTSIDE it are the background sample; the 99th percentile inside it is the signal level |
| nuclei | the labels | cross-talk inside the nuclei |
| perinuclear_ring | labels dilated by `signal_qc.ring_um` (code default 1.0 um), minus the labels | cross-talk where the nuclear envelope and the P granules sit |
| cytoplasm | tissue minus ring minus nuclei | cross-talk in the cytoplasm proper |

### Per-channel columns (image_summary)

One block per channel role, named `sig_<role>_<metric>` for the roles `dna`, `central_element`,
`foci`, `granule` and `lamin` (a role a custom channel map adds beyond these rides along as an extra
column after the declared ones). All intensities are raw camera units of the channel.

| column | meaning |
|---|---|
| `sig_<role>_offset` | 0.5th percentile of the whole sub-sample: the camera offset |
| `sig_<role>_bg` | background level: the median of the voxels outside the tissue; when fewer than `signal_qc.min_bg_voxels` (code default 1000) lie outside, the 5th percentile of the channel |
| `sig_<role>_bg_source` | `outside_tissue`, or `p5` when the fallback was used (then `noise` is a lower bound, because the darkest 5 percent is a sample cut from above) |
| `sig_<role>_noise` | 1.4826 x the median absolute deviation of the background sample: the robust sigma of the camera plus shot noise |
| `sig_<role>_p99` | 99th percentile inside the tissue: the signal level |
| `sig_<role>_snr` | (p99 - bg) / noise; NaN when the noise is 0 |
| `sig_<role>_sat_frac` | fraction of ALL sampled voxels at or above the saturation level (rule below) |
| `sig_<role>_bg_per_ms` | (bg - offset) / exposure_ms; NaN unless the acquisition stage supplied `exp_ms_<role>` |
| `sig_<role>_p99_per_ms` | (p99 - offset) / exposure_ms; NaN likewise |
| `sig_<role>_bg_per_ms_pct` | bg_per_ms / laser power; NaN unless `laser_pow_<role>` is known as well |

The normalized levels are what make sessions comparable. A background of 400 counts means nothing on
its own; a background of 1.3 counts per ms per percent of laser power is a property of the sample
and the filter set, and a session that sits three times higher was imaged through the wrong filter,
with a dirty objective, or with a stain that did not take.

Reading the SNR: the 99th percentile of pure Gaussian noise sits 2.33 sigma above its mean, so a
channel with no live fluorophore reports an SNR near 2.3 regardless of how hard it was driven. The
code default `signal_qc.min_snr` of 3.0 sits just above that floor. A stained channel in a usable
image reports tens to hundreds.

Saturation level rule. The level is `2 ** bit_depth - 1` when a bit depth is known: the acquisition
record's `bit_depth_<role>` field when the acquisition stage supplies one (it reads the nd2
description, where the camera reports `Bit Depth: 12-bit`), else the `signal_qc.bit_depth` config
key (code default none; a profile for the 12-bit camera can set `signal_qc.bit_depth: 12`). Without
a bit depth the level is the maximum of the voxel dtype: the in-memory integer dtype (65535 for the
uint16 nd2 pixels), else the integer dtype the file was read as, else the float maximum, in which
case `sat_frac` is 0 by construction. A 12-bit camera writing into uint16 files saturates at 4095,
not 65535, so on such data `sat_frac` under-reports until a bit depth is supplied.

### The cross-talk table (`<image_id>__channel_crosstalk.csv` and `.parquet`)

One row per unordered pair of resolved channel roles and region, on the same sub-sample, stacked by
`germquant collect` into `batch_channel_crosstalk.csv` like every other table. The shared metadata
block (image_id, genotype, sex, voxel size, config hash, ...) is stamped on every row.

| column | meaning |
|---|---|
| `role_a`, `role_b` | the two roles, in channel-map order |
| `region` | `nuclei`, `perinuclear_ring` or `cytoplasm` |
| `n_voxels` | sampled voxels in the region |
| `pearson_r` | Pearson correlation of the two channels over those voxels; NaN when either is constant or the region is empty |
| `slope` | least-squares slope of the offset-subtracted `role_b` on the offset-subtracted `role_a`: the leak fraction when `role_b` is a scaled copy of `role_a` |

Two stains of different structures correlate weakly (|r| below about 0.3 in the ring, where both
have structure). A channel that is a scaled copy of another, through bleed-through or because the
same fluorophore is seen through two filters, gives r near 1 with a slope equal to the leak
fraction. Note that the DNA channel correlates with every nuclear stain inside the nuclei by
biology, so the `nuclei` row of a `dna` pair is the least informative one; the `perinuclear_ring`
row is where a lamin leak into a granule or SYP channel shows.

### Flags

Off by default (`signal_qc.flags.enabled`, code default false), because a new flag on an existing
run would change its `qc_flags` column and break the golden comparison. A profile that enables them
gets, with the code-default thresholds:

| flag | fires when |
|---|---|
| `signal:<role>_low_snr` | `sig_<role>_snr` < `signal_qc.min_snr` (3.0) |
| `signal:<role>_saturated` | `sig_<role>_sat_frac` > `signal_qc.max_sat_frac` (0.001) |
| `crosstalk:<role_a>_<role_b>_<region>` | `pearson_r` > `signal_qc.crosstalk.max_r` (0.8) |

A NaN metric never fires. The flags are advisory: they appear in `qc_flags` and in the stage's
record in `<image_id>__stages.json`, and do not set `qc_pass` to false.

### Configuration keys and code defaults

Every key is read with a code default, so the shipped configs are not edited and their hashes keep
reproducing; a profile sets only what it changes. The stage's sub-hash in `run_manifest.json` and
`<image_id>__stages.json` covers the `signal_qc` section.

```yaml
signal_qc:
  enabled: true              # the stage switch (--no-signal_qc on the command line)
  stride: 4                  # xy sub-sampling step; every z plane is kept
  tissue_dilation_um: 2.0    # labels dilated by this = tissue; outside it = background sample
  ring_um: 1.0               # labels dilated by this, minus labels = perinuclear ring
  min_bg_voxels: 1000        # fewer background voxels -> bg falls back to the 5th percentile (bg_source p5)
  bit_depth: null            # saturation level 2**bit_depth - 1 when the acquisition record gives none
  min_snr: 3.0               # flags only
  max_sat_frac: 0.001        # flags only
  crosstalk:
    max_r: 0.8               # flags only
  flags:
    enabled: false           # the flags above are written only when true
```

### What it does not do

It does not correct anything: a channel that fails is reported, not fixed, and the stages that
consume it still run. The partition stage's "not measurable" verdict when the SYP signal sits at the
floor is a separate item that reads these numbers. The background sample assumes part of the field
is empty; a stack cropped tightly to the tissue triggers the `p5` fallback and its lower-bound
noise. Cross-talk by correlation cannot tell bleed-through from true colocalization of two stains
on the same structure; the `slope` and the comparison across sessions and filter sets do that, and
the acquisition section's filter names say which sessions to compare.
