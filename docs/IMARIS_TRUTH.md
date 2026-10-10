# The Imaris truth reader: `germquant imaris-truth`

The lab scores RAD-51 in Imaris: automatic Spots and Surfaces, curated by threshold inside a
hand-drawn box, saved as an `.ims` project (HDF5). This command reads such a project directly and
writes the expert truth as plain files in OUR coordinate frame, so the benchmark can compare the
pipeline's spots and nuclei with the expert's on the same image, voxel for voxel. It is the first of
three packages of the RAD-51 cross-validation (B1 writes the truth files, B3 consumes them).

Nothing in the pipeline changes. The reader lives in `germquant.validate.imaris_ims` next to the
older `read_imaris_ims` (unchanged), and the command is wired into `germquant.cli` like `get-model`.

    .venv\Scripts\germquant.exe imaris-truth "<project>.ims" --out <folder> --nd2 "<the same image>.nd2"
    .venv\Scripts\germquant.exe imaris-truth "<project>.ims" --out <folder> --no-masks     # tables only, fast

`--nd2` names the raw stack the project was made from; its voxel size defines our frame. Without
it a sibling `<stem>.nd2` next to the `.ims` is used, and without that Imaris's own voxel size (the
params JSON then carries a warning, because the two frames differ, see below). Reading an `.ims`
needs `h5py` (`uv pip install -e ".[imaris]"`; the `spots` extra already brings it in). The file is
opened read only.

From Python: `export_imaris_truth(ims, out_dir, nd2=..., masks=True)` returns the params dict.

## The truth file contract

For a project named `<stem>.ims` the output folder holds four files. The column lists below are the
contract; the columns after them are extras that a consumer may ignore (select columns by name, the
contract columns come first and in this order).

`<stem>__imaris_spots.csv`, one row per Imaris spot:

| column | meaning |
| --- | --- |
| `spot_id` | Imaris's spot ID (the `ID` field of the Spot table) |
| `x_um`, `y_um`, `z_um` | position in OUR frame: voxel index times the nd2 voxel size |
| `xi`, `yi`, `zi` | integer voxel indices, `floor((imaris - ExtMin) / imaris_voxel)`; not clipped |
| `radius_um` | the Imaris radius (one value per export, or per spot when region growing was on) |
| `quality` | the Quality statistic; NaN when the export does not store it (none of the 12 files does) |
| `spot_class` | the Imaris class label(s), e.g. `Class A`; several joined with `;`; empty when none |
| `imaris_x_um`, `imaris_y_um`, `imaris_z_um` | the raw Imaris positions (stage micrometres, float32) |
| extra `imaris_mask_label` | Imaris's own "Intensity Max" of the mask channel for the spot: the lab's spot-to-nucleus assignment (0 outside every nucleus; empty when the statistic is absent) |
| extra `spots_object` | the Imaris object name (`Spots 1`) |

`<stem>__imaris_nuclei.csv`, one row per Imaris surface:

| column | meaning |
| --- | --- |
| `surface_id` | Imaris's surface ID (SurfaceModel `ID`) |
| `label` | the value of this nucleus in `<stem>__imaris_labels.tif`; 0 when it has no mask |
| `x_um`, `y_um`, `z_um` | centre of mass in OUR frame (voxel index times the nd2 voxel size) |
| `xi`, `yi`, `zi` | voxel indices of the centre of mass |
| `volume_um3` | Imaris's surface volume |
| `n_voxels_mask` | voxels carrying `label` in the label image; empty when there is no mask |
| `imaris_x_um`, `imaris_y_um`, `imaris_z_um` | the raw centre of mass (stage micrometres) |
| extra `mask_centroid_dist_um` | distance between the surface's centre of mass and the centroid of its mask blob (the match quality) |
| extra `surface_class` | the Imaris class label(s) of the surface |
| extra `surfaces_object` | the Imaris object name (`Surfaces 1`) |

`<stem>__imaris_labels.tif`: an int32 ZYX label image of the IMAGE size (not Imaris's padded
array), 0 where no traced surface, written by the pipeline's own label writer (OME-TIFF, zlib,
`PhysicalSize` tags in our voxel size), so `tifffile.imread` gives the same array shape and dtype as
the pipeline's `<id>__nuclei_labels.tif`. Absent when no mask channel could be identified or when
`--no-masks` was given; `mask.reason_no_labels` in the params JSON says why.

`<stem>__imaris_params.json`: everything about the export. Voxel sizes of both frames
(`voxel_imaris_xyz_um`, `voxel_ours_xyz_um`, `our_voxel_source`), `ext_min_xyz_um`,
`ext_max_xyz_um`, `shape_xyz`, `data_shape_zyx_padded`, `channel_names`, `imaris_version`; the full
`CreationParameters` text of every Spots and Surfaces object (`creation_parameters`, keyed by object
name) and per object its ROI boxes, class label groups and class counts, the names of the stored
statistics, and whether region growing was on (`objects`, plus the file-level `region_growing` and
`quality_present` flags); `roi_boxes` and `roi_check`; `n_spots`, `n_surfaces`,
`spot_classes_present`, `surface_classes_present`; and `mask`: which channel served as the mask,
how it was identified (`identification`, with the score of every channel), the matching strategy
and tolerance, `n_surfaces_with_mask`, the match distances, how many values were shared and how
they were split, the orphan values dropped, the centroid sanity check
(`n_centroid_inside_own_label`, `n_centroid_in_cavity`, `n_centroid_in_other_label`), the
agreement with Imaris's Intensity Max statistic, `labels_bbox_xyz_voxel` (the traced region) and
the label image's shape. `warnings` lists everything worth a look.

## The two frames

Imaris keeps objects in stage micrometres. The image starts at `DataSetInfo/Image ExtMin{0,1,2}`
and Imaris's voxel size is `(ExtMax - ExtMin) / size`, which is ROUNDED: 0.108 um here where the nd2
says 0.108333 um. Voxel `i` of an axis covers `ExtMin + i * v .. ExtMin + (i + 1) * v`, so an object
position `p` lies in voxel `floor((p - ExtMin) / v)` with Imaris's `v`. Our tables place voxel `i`
at `i * v_nd2`. Converting an Imaris position to our frame therefore goes through the voxel index
(`x_um = xi * v_nd2`), never through "Imaris micrometres minus ExtMin": that shortcut drifts by
`i * (0.108333 - 0.108)`, almost a micrometre at the far edge of a 2600-pixel field (measured on
the real spots: mean 0.37 to 0.50 um, up to 0.86 um). Spot positions are sub-voxel (their fractional
voxel coordinates are uniformly distributed), so `floor` gives the containing voxel and `x_um` is
that voxel's position in the pipeline's convention; the sub-voxel detail stays in `imaris_*`.

## What an .ims holds (the 12 N2 RAD-51 projects, Imaris 10.2 and 11.0)

Checked with h5py on the files under `C:/Users/ryane/n2_rad51_local/imaris/`:

* `DataSet/ResolutionLevel 0/TimePoint 0/Channel k/Data`: uint16, gzip, chunks (8, 256, 256), and
  PADDED to chunk multiples: (72, 2816, 2816) for an image of 2600 x 2600 x 71. The image is
  `[:Z, :Y, :X]`, the padding is zero. Channels 0 to 2 are 405 (DAPI), 477 (SYP), 545 (RAD-51);
  channel 3 is `Masked 405`, the lab's mask channel (see below). The reader reads it in z-chunks
  of the file's chunk size straight into an int32 array and never loads a whole array as float.
* `Scene8/Content/Points0` (`Spots 1`): the `Spot` table (ID, PositionX/Y/Z as float32, Radius),
  statistics (`StatisticsType` names a statistic and, through `Factor`, its channel; `StatisticsValue`
  holds one row per object and statistic, keyed by the object ID), class labels, and the creation
  parameters. The `CreationParameters` DATASET holds the full XML; the group attribute of the same
  name is a shorter copy, so the dataset is read first.
* Quality is not stored: the spots were thresholded on Quality in the wizard (the threshold is in
  the creation parameters, `mLowerThresholdManualValue`), but Quality is in none of the 32 stored
  statistics of any of the 12 files, so `quality` is NaN throughout. Stored per spot are the
  intensity statistics of every channel, position, diameter, volume, distances to neighbours and
  to the image border, and `Shortest Distance to Surfaces`.
* Classes: `LabelGroupNames` (`Set 1`), `LabelValues` (`Class A`, `Class B`), `LabelSets`
  (objects grouped by their label combination, cumulative end indices into `LabelSetLabelIDs` and
  `LabelSetObjectIDs`). The decoding reproduces Imaris's own `Number of labelled Spots per Class`
  counts exactly. These classes are an AUTOMATIC classification the wizard applied, by `Average
  Distance To 3 Nearest Neighbours` with a 7.85 um threshold for the spots and `Distance from
  Origin` for the surfaces, not a manual curation. A surface can carry both labels.
* `Scene8/Content/MegaSurfaces0` (`Surfaces 1`): `SurfaceModelInfo` (centre of mass, volume, area,
  ellipsoid axes), `SurfaceModel` (ID and the axis-aligned bounding box `ExtendMin/Max`), the same
  statistics machinery, and the creation parameters: smoothing 0.4 um, a manual threshold per
  gonad, seeded region growing with a 2.5 um seed diameter, and a `Number of Voxels` filter (for
  example 4514.5 to 22647). The surface voxel masks themselves live in the proprietary block tree
  and are not read.
* The ROI box of the Surfaces wizard is recorded (`<bpRegionOfInterest mMinXYZ=... mMaxXYZ=...>`,
  voxel indices, inclusive) and the reader reports it, but on both real files it does NOT contain
  the surfaces: 65 of 146 and 31 of 88 surface centroids lie inside it. It records a wizard state,
  not the traced region. The traced region is `mask.labels_bbox_xyz_voxel`, the bounding box of the
  kept labels; a warning says so whenever the box misses surfaces.
* Region growing: two HS male exports have `mEnableRegionGrowing="true"` on the spots
  (`20251105_n2_hs_MALE _.ims` with 772,485 spots and only one class label, and
  `20251105_n2_hs_MALE _002.ims`, 13,951 spots, saved by Imaris 11.0). Their radii vary per spot.
  The flag is `region_growing` in the params JSON.
* The mask channel: Imaris filled every surface with its own value in `Masked 405` (the lab's
  coloc-channel trick). The values run 1..N in roughly surface order but are neither the surface IDs
  nor the row indices, and the channel holds a few values no surface claims: nuclei deleted after
  masking (whole blobs of 13,000 to 22,000 voxels) and small fragments (5 to 364 voxels). Some
  nuclei have a cavity in the fill where the centre of mass lands.

## How the nucleus masks are recovered

1. Identify the mask channel with a pixel probe: for every channel, read the plane of each surface
   centroid and look at the 7 x 7 patch around the centroid voxel. In a mask channel the patch holds
   exactly one distinct nonzero value for most surfaces (a cavity at the centroid still counts) and
   those values differ between surfaces; a fluorescence channel varies from voxel to voxel. The
   channel with the highest fraction of such patches (at least half) wins; the score of every
   channel is in `mask.identification`. On the real files channel 3 scores 145 of 146 and 87 of 88
   with as many distinct values; the fluorescence channels score 0.
   Imaris's per-surface intensity statistics cannot do this job: the Intensity Min of the mask
   channel is 0 for every surface (Imaris counts boundary voxels outside the fill), so "Min equals
   Max" never holds. The Intensity Max statistic is still read, as a cross-check (below) and as the
   label source for `--no-masks`.
2. Read the channel as int32 labels, in z-chunks. When the values are not unique per surface (a
   binary mask), connected components of the nonzero channel become the blobs instead.
3. Match each surface inside Imaris's own bounding box of that surface (widened by two voxels): of
   the values present there, take the one whose voxels have their centroid nearest to the surface's
   centre of mass, within 1 um. A neighbour reaching into the box has its centroid near the box
   edge and is not chosen. On the real files the chosen blob is 0.015 um (median) to 0.036 um (max)
   from the centre of mass and the next blob is more than 2 um away.
4. A value claimed by two surfaces is split: first into connected components (one value on two
   separate nuclei, which happens with random mask values), then, for a component holding several
   centroids, by nearest surface centroid. The second surface gets a fresh label above the largest
   value. The real files needed no split.
5. Values no surface claims are set to 0 and listed as orphans (`mask.orphan_values`,
   `n_orphan_voxels`). Then `n_voxels_mask` is counted per label and the label under each surface's
   centroid voxel is checked: its own label, 0 (a cavity in the fill) or another surface's label.

## Validation on the real files (2026-10-09, workstation, GPU hidden)

| | `20251105_N2_nohs_HERM _001` (noHS) | `20251105_n2_hs_HERM _` (HS) | `20251105_n2_hs_MALE _` (HS, region growing) |
| --- | --- | --- | --- |
| spots, surfaces | 1222, 146 | 5152, 88 | 772,485, 64 |
| surfaces with a mask | 146 of 146 | 88 of 88 | 64 of 64 |
| label image | (71, 2600, 2600) int32, 2.6 MB | (67, 2600, 2600) int32, 2.3 MB | (67, 2600, 2600) |
| match distance, median / max | 0.016 / 0.036 um | 0.015 / 0.028 um | |
| centroid in own label / cavity / other label | 142 / 4 / 0 | 75 / 13 / 0 | 63 / 1 / 0 |
| orphan values dropped (voxels) | 5 (14,747; one whole nucleus, four fragments) | 4 (79,599; four whole nuclei) | 7 |
| mask voxels times voxel volume over Imaris volume | 1.06 (1.035 to 1.11) | 1.05 (1.03 to 1.08) | |
| label equals Imaris's Intensity Max statistic | 143 of 146 | 85 of 88 | |
| spot's voxel label equals `imaris_mask_label` | 1199 of 1222 (98.1 %) | 5118 of 5152 (99.3 %) | |
| wall time, peak working set | 35 s, 4.4 GB | 31 s, 4.2 GB | 40 s |

Every label in the nuclei table is unique and present in the label image, and the image holds no
label without a nuclei row. The six surfaces whose recovered label differs from Imaris's Intensity
Max statistic all have a neighbour with a larger value touching them: one a traced nucleus 4.9 um
away, five an orphan blob or fragment of a nucleus deleted after masking. The statistic takes the
maximum over the surface's voxels and so reports the neighbour; the centroid match reports the
surface's own fill. This is the collision the lab's coloc-max join is known to suffer from (about
3 % of nuclei, see `imaris_xlsx.py`), so the label image, not the statistic, is the truth for
per-nucleus work. The spot disagreements are all boundary spots whose centre voxel is outside every
nucleus while Imaris's 0.2 um sphere reaches one (ours 0, Imaris nonzero; never a different
nucleus). `--no-masks` on the HS file took 23 s and 0.14 GB and gave the statistic's label for
every surface, equal to the full run for 85 of 88.

The peak memory is the int32 label volume (1.9 GB for 71 x 2600 x 2600) plus its copy inside the
TIFF writer; the mask channel itself is never held as float.

## Notes for the benchmark (B3)

* Per-nucleus truth: join spots to nuclei through the label image (`labels[zi, yi, xi]`, with the
  spot's `xi, yi, zi`) and the nuclei table's `label`. `imaris_mask_label` is the lab's own rule and
  differs for boundary spots (about 1 to 2 %), never by assigning another nucleus.
* Only a subset of the nuclei of a gonad is traced; spots outside every traced nucleus (58 % of the
  spots on the noHS file, 93 % on the HS file) are real detections in untraced nuclei and
  cytoplasm, not errors. Compare our nuclei to the truth only where a traced surface exists.
* Orphan blobs are not truth (they were curated away). The recorded ROI box is not a spatial
  filter; use `mask.labels_bbox_xyz_voxel`.
* `quality` is NaN for this data set; `spot_class` is an automatic neighbour-distance class.
* The exports with region growing have per-spot radii and, for `20251105_n2_hs_MALE _.ims`, 772,485
  spots (the detector flooded); its spots table is large and its class labels are empty.

## Tests

`tests/test_imaris_truth.py` writes a tiny synthetic HDF5 in the .ims layout (padded arrays, Spot
and SurfaceModelInfo compound tables, statistics with Quality and the mask-channel Intensity Max,
class label sets, creation parameters with an ROI, extent offsets) and checks the frame conversion,
the mask recovery (a shared value on two nuclei split into two labels, an orphan blob dropped, a
surface without a mask), the four written files and the params JSON, the `--no-masks` path, the
binary-mask fallback, an export without surface statistics, the Imaris-voxel fallback, the nd2
header path, and the CLI wiring.
