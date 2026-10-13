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
  The flag is `region_growing` in the params JSON. (Corrected 2026-10-12, protocol changelog entry 11: in
  `20251105_n2_hs_MALE _.ims` the radii do not vary. Its 772,485 spots are unthresholded, unclassified
  candidates, every one at radius 0.2 um, saved under creation parameters that say region growing on and
  Quality 13.358; its curated 7,628 spots are in the project's 2025-12-31 statistics export, see "The
  recorded earlier export" below.)
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

## The frame check and the xlsx spot source (protocol changelog entry 9)

    germquant imaris-truth "<project>.ims" --out <folder> --nd2 "<the same image>.nd2" --frame-check
    germquant imaris-truth "<project>.ims" --out <folder> --nd2 "<image>.nd2" --spots-xlsx "<export>.xlsx" \
        --xlsx-expect-sha256 <hex> --xlsx-expect-spots <n> --xlsx-expect-surfaces <n>

**Why.** Imaris stores the extents as text with six significant digits. On the 12 registered projects that is
exact (2600 x 0.108 = 280.8), but a field whose stage position exceeds 10000 um in magnitude loses the second
decimal: the 2026-01-05 males read ExtMin 16249.9 and ExtMax 16543.7 on x, so the stored-extent voxel is 0.1080147
um on that axis and 0.108 on the other, and ExtMin itself is known only to +-0.05 um, half a voxel. Spot positions
are sub-voxel, so the frame decides which voxel a spot is in.

**The test** (`frame_check`, `frame_agreement`). Imaris stores per spot the Intensity Center of every channel: the
value of the voxel that holds the spot's centre, an integer. For a candidate frame (ExtMin, voxel), a spot agrees
when the image of its source channel (from the .nd2 given with `--nd2`, else the .ims channel) at
floor((p - ExtMin) / v) equals its Intensity Center, and when the mask channel there equals its Intensity Center of
the mask channel. Positions are float32, so a spot within two float32 steps of a voxel boundary also agrees with
the voxel across it (`agreement`; the plain comparison is `agreement_strict`). A frame passes with at least 99.5
percent of the spots overall and, for the spot channel, in each third of the field in x and in y. The frames are
tried in order: `stored` (the reader's rule), `rounded` (the voxel rounded to 4 decimals), then a search that moves
ExtMin in 0.001 um steps inside the rounding interval of its text, for either voxel (per axis to the highest
agreement; ties go to the higher strict agreement, then to the smaller shift, so the search does not stop at the
edge of the agreeing range where a spot sits on a voxel boundary). The first that passes
converts every Imaris position of the export (spots and surfaces): `ext_min_xyz_um` and `voxel_imaris_xyz_um` in
the params are then that frame (the stored ones are kept under `*_stored`), and `frame_check` lists every
candidate with its agreement overall and per third. When nothing passes, the stored frame is kept,
`frame_check.chosen` is None, and the benchmark manifest marks the spot truth unusable. An export whose Spots do not
store the Intensity Center cannot be tested (`frame_check.untestable`).

On the eight registered train and val gonads the stored frame agrees for every spot on both channels; the strict
comparison misses 0 to 330 spots per gonad (0 to 0.27 percent), every one within 0.003 voxel of a boundary
(docs/BENCHMARK.md, "Entry 9").

**The read guard** (protocol changelog entry 10). With `--frame-check` or `--spots-xlsx`, the .ims, the .nd2 and
the xlsx are checked before anything is opened (`benchmark_split.guard_truth_read`): a file whose name is the stem
of a registered gonad (`benchmark_split.REGISTERED_STEMS`) passes the role guard, so a train or val project is read
and a test gonad's project is refused; any other file must not be sealed (`SEALED_NAME_MARKERS`: a 2025-11-05
export, the lab's aggregate tables). `--pre-flight` allows a test gonad, for the frame check of each test gonad in
the pre-flight of docs/FINAL_SCORING.md only (the Intensity Center test, no metric).

**What the benchmark does with it.** The manifest marks an addition's spot truth usable only when the export ran
the frame check and it chose a frame (an export without `--frame-check` is not usable), and then by the decision of
a changelog entry or a section 4.4 agreement of at least 99.35 percent (`benchmark_split.addition_spot_truth_usable`;
docs/BENCHMARK.md, "Entry 9").

**The xlsx spot source** (`--spots-xlsx`; entry 9 (e)). For a project whose saved Spots object is not the curated
set (20260105_n2_hs_1_MALE_1 holds 1,226,263 saved spots; its 2026-01-13 export the curated 6,512), the spots come
from the Imaris statistics export (`imaris_xlsx.read_export_objects`: the Position, Diameter, Intensity Center,
Intensity Max and Volume sheets and the Overall totals). The export is verified first, and every criterion is in
`spots_source.verification` of the params: V1 the sha256 is the expected one; V2 the Overall totals and the
Position sheet hold the expected numbers of spots and surfaces; V3 the export's surfaces equal the .ims Surfaces
one to one (by ID, or by position when the ID sets differ): position within 0.005 um, volume within 1e-5 relative,
the same Intensity Max of the mask channel; V4 every spot lies inside the .ims extent; V5 the frame test above
passes on the export's own Intensity Center values; V6 every radius (Diameter X and Y over two) is 0.2 um; V7 the
section 4.4 agreement (`truth_assignment_agreement`, the rule of `scripts/check_truth_assignment.py`) is at least
99.35 percent. `passed` is true only when all seven pass. The spots table then holds the export's spots, placed
with the frame V5 chose; `imaris_mask_label` is the export's Intensity Max of the mask channel (the lab's coloc
value); Quality and the class are empty (the export holds no per-spot class); `spots_object` is "Spots 1 (xlsx)".
Reported, not a criterion: the fraction of the export's spots within 0.005 um of a spot of the saved object. The
nuclei and the label image come from the .ims as always. A failed verification still writes the files, with
`passed` false and a warning; the manifest then marks the spot truth unusable and the gonad stays nucleus truth
(a failure of V7 alone, the section 4.4 check, can be decided otherwise only by a changelog entry, as for a gonad
whose spots come from the .ims; protocol changelog entry 10 (d)).

The params also carry `creation_parameter_summary` (`benchmark_split.creation_parameter_summary`): the Spots and
Surfaces settings entry 9 (d)(4) compares with sections 4.1 and 4.2 of the protocol.

## The recorded earlier export (protocol changelog entry 11)

    germquant imaris-truth "<project>.ims" --out <folder> --nd2 "<image>.nd2" --spot-source tests/golden/rad51_spot_sources.json

For a registered TRAIN gonad whose saved Spots object is not the thresholded, classified set but whose thresholded,
classified spots survive in an Imaris statistics export of the same project (n2rad51_HS_male_000: 772,485
unthresholded candidates in the .ims, a subset of 7,628 with classes in its 2025-12-31 export
`20251105_n2_hs_MALE_3.xlsx`; the .ims was saved again on 2026-01-21, after that export), `--spot-source` takes the
gonad's entry from the record of entry 11 (`benchmark_split.load_spot_sources`, checked against the split record,
`--split-file`, by default `rad51_split.json` next to the record) and builds the spot truth from the export without
replacing the .ims as the file of record. "Curated" in the record and in entry 11 is an inference from the
thresholded subset and from the class counts the project's own classifier reproduces; the Quality threshold of
the 7,628 is stored in neither file, and the settings the export was made with are not recorded (protocol
changelog entry 12):

* the guard: the project must be the gonad's registered one (by its stem) and have the size the source records
  (before it is opened); the export is admitted only as the file pinned in code
  (`benchmark_split.SPOT_SOURCE_EXPORTS`: name and size without reading a byte, then the role guard, a train
  gonad, then its sha256, all before anything parses it; entry 12); the 2025-12-31 exports of the other heat-shock
  males stay sealed, and `--pre-flight` is refused;
* the spots: the .ims spots the export names, matched one to one by position within 0.005 um (the export
  renumbers the IDs and writes positions to 0.001 um), with their .ims IDs, positions and radius, in the .ims's
  order, converted with the frame the check of entry 9 (d)(3) chooses on the .ims's own spots (the frame check is
  implied); the class from the project's own one-dimensional classification (`classification_rule`: Class A below
  the threshold stored in the Spots creation parameters, Class B at or above) applied to the Average Distance To
  3 Nearest Neighbours recomputed among the matched .ims positions (`average_distance_to_neighbours`). Those
  positions are float32 and lie within 0.0013 um of the positions the export was made from, so the recomputed
  distances are not the export's own values unrounded; on n2rad51_HS_male_000 they differ from the written values by
  -0.00097 to +0.00059 um. They are therefore checked spot by spot against the export (E8 below); `imaris_mask_label`
  the export's per-spot Intensity Max of the mask channel; `spots_object` "Spots 1 (export subset)";
* the checks, in `spots_source.verification`: E1 the export's sha256 is the pinned one and the .ims's the one the
  source records (entry 12; hashed before the .ims is opened); E2 the Overall totals and class counts and the
  Position sheet hold the recorded numbers, and the image's voxel count is the .ims's; E3 the export's surfaces
  equal the .ims Surfaces one to one (position within 0.005 um, volume within 0.0005 um^3, half the export's last
  digit, plus two float32 steps; the mask Intensity Max equal); E4 every truth spot inside the image; E5 the
  Intensity Center test on the export's own values chooses the frame of the .ims spots; E6 every spot at the
  recorded diameter within 1e-6 um (0.4 x 0.4 x 0.9 um: no region-grown sizes. Protocol section 4.1 has no rule
  that reads the export settings from the diameter, and the settings are not recorded: the registered .ims's own
  flag says region growing on while its saved spots are also constant); E7 the position match is one to one, every
  export spot within 0.005 um of exactly one .ims spot; E8 the classification reproduces the export's Overall class
  counts and, spot by spot, the classes the export alone forces (`classes_from_written_values`, entry 12: a spot
  whose written distance lies more than 0.0005 um plus two float32 steps from the threshold takes the class of
  that value, and the Overall counts must give the remaining spots one class; on n2rad51_HS_male_000 one spot is
  left, export ID 386332, written 4.115 against 4.11504, and the counts make it Class B, which the recomputed
  distance 4.115152 also gives); E9 the export's mask Intensity Max equals the .ims's for every matched spot (the
  independent check of the mask label: the truth's `imaris_mask_label` is the export's, so the per-surface coloc
  join of `scripts/check_truth_assignment.py --spot-sources` compares the export with itself); E10 the section 4.4
  agreement at least 99.35 percent. `passed` is true only when all ten pass. The params of entry 11 (written at
  e187f8b) carry the earlier E1 and E8 (no .ims digest, totals only) and the earlier wording of E6's rule;
* the files: the four truth files as always (the nuclei and the label image from the .ims, which must equal the
  registered export of the project: `scripts/check_truth_assignment.py --spot-sources` compares them) and
  `spot_source_id_map.csv` (export ID, .ims ID, match distance, the recomputed and the exported 3-neighbour
  distance, the class). `region_growing` in the params is that of the spots used (false when E6 passes); the flag
  of the saved Spots object is kept as `region_growing_saved_spots_object`.

The manifest (`scripts/benchmark_manifest.py --spot-sources --spot-sources-truth-root`) reads this export for the
gonad and marks its spot truth usable by `benchmark_split.spot_source_truth_usable`: E1 to E9 and the frame check,
which no decision overrides, then a changelog decision, else the section 4.4 agreement of the `--truth-check` table
of at least 99.35 percent.

## Tests

`tests/test_imaris_truth.py` writes a tiny synthetic HDF5 in the .ims layout (padded arrays, Spot
and SurfaceModelInfo compound tables, statistics with Quality and the mask-channel Intensity Max,
class label sets, creation parameters with an ROI, extent offsets) and checks the frame conversion,
the mask recovery (a shared value on two nuclei split into two labels, an orphan blob dropped, a
surface without a mask), the four written files and the params JSON, the `--no-masks` path, the
binary-mask fallback, an export without surface statistics, the Imaris-voxel fallback, the nd2
header path, and the CLI wiring.

The frame check and the xlsx spot source (entry 9): the frame rule on synthetic spots (the stored frame; a voxel
that drifts over the field, caught in the far third and fixed by the rounded voxel; an ExtMin that lost 0.03 um,
found by the search; values that belong to no frame; the float32 boundary allowance), the export with
`frame_check_on` on a synthetic project that stores the Intensity Center (the same truth as without it), an
export without the statistic (untestable), and the xlsx spot source on a synthetic export written with openpyxl:
all seven criteria pass, and a wrong sha256, a moved surface, a wrong radius and a spot outside the extent each
fail their criterion. Added with protocol changelog entry 10: V2 (one spot more expected), the volume part of V3
(volumes 5e-5 off), V5 (Intensity Centers off by one) and V7 (the 99.35 percent minimum) each fail alone; the pass
rule's thirds and mask channel; the voxel rounded to four decimals (0.10801 on a 0.108 um grid of 2720 voxels); an
export whose stored x extent lost 0.03 um, converted with the frame the search chose (every spot in its true
voxel, the nuclei in the chosen frame), with and without the float32 allowance (the search's tie-break); the coloc
value taken from the export's Intensity Max of the mask channel, not its Intensity Center; and the refusals of the
read guard (a test gonad's project or .nd2, a sealed xlsx or .nd2) before any file is opened, with `allow_test`
letting a test-named project through.

The recorded earlier export (entry 11): the classification parsed from a creation-parameter text and the class
boundary (4.115 rounded falls below 4.11504, 4.115243 above), the 3-neighbour distance, the one-to-one position
match (a spot taken twice, an ambiguous second neighbour, a spot beyond the tolerance), the volume rule of E3; on
the synthetic project under n2rad51_HS_male_000's registered stem with an export that renumbers the IDs, rounds the
positions and writes one distance within half a digit of the threshold on the wrong side: every criterion passes,
the truth spots carry the .ims IDs and the recomputed classes, and the nuclei equal the plain export's; E1, E2, E3,
E5, E6, E7, E8, E9 and E10 each fail alone on a defect of their own. Added with protocol changelog entry 12: E1 on
an .ims digest that is not the recorded one; E2 alone when the record expects other class counts than the
export's Overall; E8 alone on a compensating swap of two written distances and on a sheet that says 99 for every
spot; E4 alone in `export_subset_spots` on a 1 x 1 x 1 image; `classes_from_written_values` (one spot within half
a digit resolved by the counts, two left open, counts that do not fit, a missing value, values just beyond the
tolerance). Refused before anything is opened: another gonad's project, an export of another size than the pin, an
export with the pinned name and size but other bytes (hashed, never parsed), a project of another size than the
source records, and two spot sources at once.
