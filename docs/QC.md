# Quality control

What the pipeline records about each image beyond the biology, and the checks it runs on those
records. Every check writes a flag string into the `qc_flags` column of `image_summary` (also in
`batch_summary.csv` and in the per-stage record `<image_id>__stages.json`). A flag is advisory
unless it names a missing required channel, a role collision or an image with no nuclei: only those
set `qc_pass` to false. The flags described here never fire on a run that does not switch them on,
so the `qc_flags` of every existing run and of the synthetic goldens (`tests/test_golden.py`) are
unchanged.

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

To be written by the signal QC package: per-channel background and signal-to-noise normalised by
exposure and laser power, the cross-talk check between channels in the envelope shell, and the
partition stage's "not measurable" verdict when the SYP signal sits at the floor.
