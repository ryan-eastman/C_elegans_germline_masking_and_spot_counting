# Partition measurability: was there a SYP signal to partition?

The partition stage (`partition.enabled`, switched on by the `ccw77_partition` profile) reports how
much more SYP-3 sits in the P granules than in the cytoplasm at the same distance from the nuclear
envelope. That ratio means something only when the SYP channel carried a live fluorophore. Two
`image_summary` columns and one opt-in QC flag record whether it did.

## The columns and the flag

| Column or flag | Written when | Value |
|---|---|---|
| `partition_syp_snr` | the partition stage ran | the SYP signal-to-noise the signal_qc stage measured (its `image_summary` field `sig_central_element_snr`); NaN when that stage did not run or the SYP role did not resolve |
| `partition_measurable` | the partition stage ran | `partition_syp_snr >= partition.min_syp_snr` when the gate is on and the ratio is known; otherwise True |
| `partition:low_syp_snr` (in `qc_flags`, and on the partition entry of `<image_id>__stages.json`) | the gate is on and the image sits below it | this image's partition numbers are a floor artefact, not a measurement |

Both columns are NA when the partition stage did not run (switched off, skipped because the granule
stage did not run, or failed), like every other `partition_*` column.

The gate is `partition.min_syp_snr`. Its code default is 0.0, which means off: on the default path
the two columns are written, `partition_measurable` reads True, and no flag fires, so the `qc_flags`
text of every existing run is unchanged. The flag is a profile choice.

Every partition number is still written for a flagged image. The flag, not a blank, is the signal:
a reader can see what a dead channel produces (the floor and its controls) and filter on
`partition_measurable` when pooling.

## Switching the gate on

Add the key to a study profile. It is deliberately not in `config/profiles/ccw77_partition.yaml`:
the text of that profile is part of the `config_hash` of every run made with it, a commented line
changes the hash as much as a live one, and its real goldens must keep reproducing. A profile that
extends it, inside the repository or anywhere else, is a new recipe with its own hash:

```yaml
# config/profiles/ccw77_partition_gated.yaml
extends: ccw77_partition.yaml
partition:
  min_syp_snr: 3.0    # flag partition:low_syp_snr when the SYP signal-to-noise is below 3
```

The key is live in this new file (a commented-only block would load as an empty section). 3.0 is a
starting point, not a validated threshold: set it
from the distribution of `partition_syp_snr` across sessions whose fluorophore is known to be live
and the session known to be dead, and keep the value in the profile so it enters the hash.

The ratio itself comes from the signal_qc stage (per-channel background and signal-to-noise,
normalised by exposure and laser power, built alongside this check as part of the acquisition and
channel QC work), which writes `sig_central_element_snr` after the measure stage. With that stage
off the ratio is NaN and the gate reports the image as measurable, because nothing says otherwise: a
threshold without the signal_qc stage is not a check. Turn both on together.

## Why the flag exists: the dead mCherry::SYP-3 session

One session of ccw77 images (mCherry::SYP-3 in the 545 channel, with a far-red 647 lamin stain in a
separate channel) was taken with a fluorophore that was no longer live, imaged at 97 percent laser
power and 300 ms exposure, and all of its images were recorded through the wrong emission filter, so
the lamin stain may have bled into the 545 channel. Nothing in the tables said so. The partition
coefficient for those images came out on a floor above 1, consistently, and read like a weak
enrichment.

A channel with no live fluorophore at full power and a long exposure holds camera noise and whatever
the mismatched filter admits. A lamin stain admitted this way lands on the nuclear envelope, exactly
where the perinuclear granules sit. The partition coefficient compares granule voxels with cytoplasm
voxels at the same distance from the envelope, so both sides see the same bleed, but the ratio of two
noise-dominated means with the background subtracted is not centred on 1, and on these images it
settled slightly above it. The controls do not rescue it:

- the 180 degree rotation null rotates the granule mask in-plane and measures the same noise field,
  so it sits on the same floor;
- the z-shift floor shifts the mask 6 planes, samples the same blur, and sits there too;
- `partition_coef_specific` (the coefficient divided by the z-shift floor) is a ratio of two floors,
  near 1.

What separates a live channel from a dead one is the signal itself: the SYP level inside the nuclei
(the synaptonemal complex) against the background of the same channel, normalised by exposure and
power so that a channel pushed to 97 percent and 300 ms is not mistaken for a bright one. That is
`sig_central_element_snr`, and that is what `partition_measurable` gates on.

The August tables offer one indirect check: the partition coefficient computed with lamin-defined
shells and with DAPI-defined shells differs by a median of under 0.01, which neither confirms nor
refutes the bleed. The direct tests are the 545-versus-640 voxel correlation in the envelope band of
the raw images (the signal_qc cross-talk check) and the emission filter names in the nd2 descriptions
(the acquisition stage). The measurability gate is what keeps a dead channel out of a pooled
partition estimate in the meantime.

## Related

- `docs/COLOCALIZATION.md`: the coloc stage and the perinuclear shell.
- `docs/ROADMAP_modular_pipeline.md`: the append-only schema and golden rules these columns follow.
- `src/germquant/partition/measurable.py` is the check; `tests/test_partition_qc.py` its tests.
