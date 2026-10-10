# Quality control of the acquisition and the signal

Checks that protect the biology from a bad acquisition. They exist because a session of ccw77 images
carried a dead mCherry::SYP-3 fluorophore imaged at 97 percent laser power and 300 ms, through the
wrong emission filter, and nothing in the pipeline noticed: every stage assumed its channel carried
signal, and the partition coefficient floor came out above 1. The acquisition section (exposure,
power and emission filter per channel from the nd2 description, and the mismatch flags against a
pack's expected settings) is written by the acquisition package; this file holds the signal section.

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
