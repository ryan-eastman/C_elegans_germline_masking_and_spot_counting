# The generic profile: nuclei plus puncta in samples that are not a worm germline

`config/profiles/generic.yaml` runs the pipeline on any 3D stack that has a nuclear stain and,
optionally, a puncta channel: cultured cells, yeast, tissue sections, organoids. It overlays
`config/config.yaml` with the channel map `config/channel_maps/generic_nuclei_puncta.yaml` and
switches the C. elegans stages off. The worm path is untouched: `config/config.yaml`, the worm
channel maps and the config hash of every run made with them are exactly what they were; the
generic profile is a new recipe with its own `config_hash`.

Windows:
```
.venv\Scripts\germquant.exe run "C:\path\to\YOUR_IMAGE.nd2" --profile generic --out "C:\path\to\RESULTS_FOLDER"
.venv\Scripts\germquant.exe batch "C:\path\to\FOLDER" --profile generic --out "C:\path\to\RESULTS_FOLDER"
```
Mac / Linux:
```
.venv/bin/germquant run /path/to/YOUR_IMAGE.nd2 --profile generic --out /path/to/RESULTS_FOLDER
.venv/bin/germquant batch /path/to/FOLDER --profile generic --out /path/to/RESULTS_FOLDER
```

Run `germquant info YOUR_IMAGE.nd2` first: it prints the channel names and the voxel size, which is
what the channel map matches against.

## What runs

| Stage | In the generic profile | Why |
|---|---|---|
| read | runs | resolves `dna` / `foci` / `granule` from the channel names, then by position |
| segment | runs with the stock Cellpose-SAM model (`cpsam`) | the shipped fine-tuned model is a worm germline nucleus model |
| measure | runs | volume, centroid and per-channel intensities for every nucleus |
| germline | off | needs the SYP marker; without it every nucleus is the analysis set |
| axis | off | a gonad is a tube with a distal and a proximal end; a dish is not |
| spots | runs when a `foci` channel resolves | SpotMAX puncta per nucleus, the same detector as RAD-51 |
| granule, coloc | skip themselves | they need the `granule` and `central_element` roles; the generic map has no central_element |
| staging, envelope, audit, sc_trace, partition, granule_tail | off | worm biology: pachytene zones, the lamin envelope, the synaptonemal complex, SYP-3 / PGL-1 partitioning |
| qc, write, render | run | QC flags, tables, label images, montage |

`<image_id>__stages.json` records every stage as `ran`, `skipped` (with the reason: `disabled`, or
`missing channel role(s) [...]` for granule) or `failed`. The stock Cellpose model is downloaded
by Cellpose on first use. Without a GPU Cellpose runs on the processor, slowly; where Cellpose or
torch is not installed at all, `segmentation.nuclei.method: auto` falls back to the classical
watershed, and `image_summary.segmentation_method` records which one ran.

## Channel roles

| Role | Required | Name tokens (case-insensitive substring) | Position fallback | Output label (`marker`) |
|---|---|---|---|---|
| `dna` | yes | DAPI, Hoechst, 405, DNA, SYTOX, DRAQ5, TO-PRO, nuclei, nuclear, H2B | channel 0 | DNA |
| `foci` | no | foci, puncta, spot | channel 1 | foci |
| `granule` | no | granule, condensate | channel 2 | granule |

The tokens are generic words rather than laser lines because which wavelength carries the puncta
is experiment-specific. Channels named by laser line (`405`, `488`, `561`, `640`) therefore resolve
by position, and every role that fell back to its position raises a `channel:<role>:name_unmatched_
used_index_N` QC flag (advisory; `qc_pass` is unaffected). A single-channel image resolves `dna`
and leaves the optional roles unset without a flag; a two-channel image never resolves `granule`.

To silence the flag, or to label the spots table with the protein name, copy
`generic_nuclei_puncta.yaml`, add your channel names to `names` and your protein to `marker`, and
point `io.channel_map` of a study profile at the copy (next section). Do not edit the shipped maps:
their text is part of the config hash of every run made with them. The `central_element` role is
deliberately absent; a study whose second marker plays the SYP part (a marker that defines the
tissue and whose cytoplasmic aggregates are compared with the granules) can add it to the copy and
the granule and coloc stages run as they do for the worm, with the worm's meaning.

## Outputs, and what the worm columns mean here

Every table of `src/germquant/schema.py` is written on every run, empty for the stages that did
not run, so consumers never hit a missing file. The columns that carry a worm meaning read as
follows under the generic profile (the comments in `schema.py` say the same):

- `nuclei`: `in_germline` is True on every row (the germline stage is off, so every segmented
  nucleus is the analysis set; it is never NA). `axis_position_norm` and `axis_position_um` are NA
  on every row (no axis). `n_spots` is the puncta count per nucleus when a `foci` channel resolved.
  The worm-only columns (`n_granules`, `sc_*`, `envelope_*`, `zone`, ...) are NA. After the declared
  columns come the per-role intensities (`dna_mean_intensity`, `foci_mean_intensity`,
  `granule_integrated_intensity`, ...) for every resolved role.
- `spots`: one row per punctum with its nucleus, position in microns, SpotMAX effect size and the
  `marker` from the channel map.
- `image_summary`: `n_germline_nuclei` equals `n_nuclei`; `mean_spots` is the mean over all nuclei;
  `total_germline_length_um` and `axis_confidence` are NaN. With no `foci` channel the flag reads
  `spots:disabled_segmentation_only` (historical wording that also covers "no puncta channel"; it
  changes only in the planned hash-bump release).
- `batch_summary.csv` (from `batch` and `collect`): `n_germline` equals `n_nuclei`. The framing
  flag `qc:germline_count_outlier_<n>_vs_median<m>_review_framing` (a nucleus count above 1.8 times
  the batch median) was designed to catch two gonad arms in one frame; on a plate it marks the
  densest fields and is advisory only.
- Images: `<image_id>__nuclei_labels.tif` (int32 label image with the voxel size in the OME tags),
  `<image_id>__spots.tif` (a blob per punctum on the same grid) and `<image_id>__montage.png`
  (one panel per resolved channel plus the nucleus outlines with the puncta circled).
- Filename metadata: `genotype`, `sex`, `treatment`, `replicate` and `acquisition_date` come from
  `metadata.filename_regex`, which encodes the worm naming convention. A filename it does not match
  gets the `metadata.defaults` (`unknown` / `NA`) and `image_id` is the file stem; supply your own
  regex with the same named groups in a study profile to fill them.

## Adapting the recipe: a study profile

A study profile is a small YAML overlay. In `config/profiles/` it is selectable with
`--profile <name>`; anywhere else, pass it with `--config`. Chained `extends` is supported, so the
overlay sits on top of `generic.yaml`, which sits on top of `config.yaml`:

```yaml
# config/profiles/hela_gh2ax.yaml (a bare name is found next to this file; from outside the
# repository use the absolute path of generic.yaml)
extends: generic.yaml
io:
  channel_map: config/channel_maps/hela_dapi_gh2ax.yaml   # your copy of generic_nuclei_puncta.yaml
spots:
  spot_radius_um: 0.25          # the expected radius of your puncta (0.3 um is a RAD-51 focus)
  effect_size_min: 3.0          # tune on a few images against a hand count: docs/SPOTMAX_VALIDATION.md
segmentation:
  nuclei:
    min_volume_um3: 20.0        # debris gate; 4 um3 is sized for germline nuclei
metadata:
  filename_regex: '(?P<date>\d{8})_(?P<genotype>[A-Za-z0-9]+)_(?P<treatment>[a-z0-9]+)_(?P<replicate>\d+)$'
```

Settings that were calibrated on the worm and are inherited unchanged: `spots.spot_radius_um`,
`spots.effect_size_min`, `spots.thresholding_method` (`threshold_triangle`), the z-column merge
(`spots.merge_z_columns`, `z_merge_gap_um`, sized for a confocal axial PSF of about 0.8 um) and
`segmentation.nuclei.min_volume_um3`. They are starting points, not validated values, for any
other sample; `docs/SPOTMAX_VALIDATION.md` describes how the worm values were cross-validated.

To count puncta in the third channel as well, declare a second spot instance on the `granule`
role; it gets its own table and nucleus column:

```yaml
spots:
  instances:
    - name: pbodies                 # -> table spots_pbodies, nuclei column n_spots_pbodies
      role: granule
      expected_from_germ_cell: false
```

## Single-plane (2D) images: not supported yet

There is no `generic_2d.yaml` because the pipeline does not support single-plane images today.
What was checked (2026-10-09; a synthetic two-channel plane, classical segmenter, CPU):
`process_image` finishes without an error on a Z=1 stack, but every quantity keeps its 3D meaning.
`volume_um3` is the pixel area times a z step that does not exist, the classical segmenter seeds
and gates objects in 3D, the Cellpose path would be called with `do_3D=True` on one plane, and the
spot detector receives a z radius. A 2D path needs, in order:

1. `src/germquant/io/nd2_reader.py`: `read_nd2_metadata` already reports `is_2d` (no Z axis in the
   file), but `read_stack` pads a Z axis of length 1 and keeps whatever z calibration the file
   carries (the `nd2` reader returns the channel `axesCalibration` triple, a placeholder for a
   plane). `Stack` needs a `planar` flag, and the spacing must keep the real dy / dx while marking
   dz as not a length, instead of `pipeline.process_image` raising `voxel:UNREADABLE`.
2. `src/germquant/segment/nuclei.py`: `_cellpose` runs `model.eval(..., do_3D=True, z_axis=0,
   anisotropy=...)`; a plane needs `do_3D=False` on the (Y, X) image. `_classical` and
   `_drop_small` need an area gate (`min_area_um2`) in place of `min_volume_um3`, which multiplies
   by dz.
3. `src/germquant/measure/regionprops3d.py`: `volume_um3 = n_voxels * dz * dy * dx`; a plane needs
   `area_um2` and `n_pixels`, appended to `schema.NUCLEI` (append-only: `volume_um3` stays), and
   centroids without a z.
4. `src/germquant/spots/detect.py`: `radii_px = spot_radius_um / spacing` hands SpotMAX a z radius
   and `_merge_z_columns` assumes stacked peaks; SpotMAX's 2D mode wants a 2D image and labels, and
   the z merge must be skipped.
5. `src/germquant/pipeline.py` (`_save_labels`, `_save_spots_image`, `_save_mask_image`) and
   `src/germquant/spots/export.py`: the TIFs are tagged ZYX with `PhysicalSizeZ`; a plane should be
   written as YX.
6. `src/germquant/pipeline.py`: read `meta["is_2d"]` once and thread the planar flag to the stages
   above. `src/germquant/stages.py` needs no change: the 3D-only worm stages are off in the generic
   profile already. `src/germquant/render.py` already shows a plane as itself (`_mip`).
7. `src/germquant/schema.py`: the appended 2D columns; `tests/test_end_to_end.py`: a Z=1 fixture;
   `tests/test_golden.py`: a golden case for the 2D path.
8. Then `config/profiles/generic_2d.yaml` is a one-key overlay on `generic.yaml`.

A thin stack (a few planes) is a 3D stack as far as the pipeline is concerned and runs today,
with the caveat that Cellpose's 3D flows and the z-column merge are meant for a well-sampled z.
