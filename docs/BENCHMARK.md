# Benchmark harness: our output against the Imaris expert truth

`germquant.benchmark` scores the pipeline against the expert's Imaris annotations gonad by gonad,
with two scripts on top of it:

* `scripts/benchmark_spots.py`: the RAD-51 spot detector. Per-nucleus counts (Pearson, Lin's CCC,
  bias) and one-to-one spot matching (precision, recall, F1).
* `scripts/benchmark_nuclei3d.py`: the 3D nucleus segmentation against the expert's partial label
  image (F1 at IoU 0.5, AP over IoU 0.5 to 0.95, splits, merges, z-extent error).

Both take the same manifest, both write a per-gonad CSV, a per-nucleus CSV and a pooled JSON with
cluster-bootstrap intervals by gonad. Nothing here changes a pipeline stage: the package reads the
frozen outputs under the golden tree and the truth files, and the numeric contract of
`germquant run` is untouched.

The truth is the lab's Imaris work on the twelve N2 RAD-51 cross-validation gonads of 2025-11-05
(no heat shock and heat shock, hermaphrodites and males): threshold-curated automatic Spots and
Surfaces inside a hand-drawn ROI box, 90 to 150 traced nuclei per gonad (1198 in all) and 1000 to
5000 spots per gonad. Only a SUBSET of the nuclei was traced, so every nucleus-level score below is
either recall-like or restricted to the traced subset, and the partial-annotation rule of the nucleus
benchmark exists for that reason.

## The manifest

A JSON list of entries, or an object whose `entries` list holds them (`scripts/benchmark_manifest.py`
builds it; `germquant.benchmark.manifest.load_manifest` reads both forms). One entry per gonad:

| key | meaning |
|---|---|
| `tag` | the golden tag, e.g. `n2rad51_noHS_herm_001`; the row key of every output |
| `split` | a grouping label (`noHS`, `HS`, or a cross-validation fold); pooled scores are also given per split when there is more than one |
| `results_dir` | the pipeline output folder holding `<image>__nuclei.csv`, `<image>__spots.csv` and `<image>__nuclei_labels.tif` |
| `truth_dir` | the folder holding the truth files for this gonad (below) |
| `nd2` | the source stack; read only for the voxel size when the tables lack their `voxel_d{z,y,x}_um` columns |
| `ims` | the Imaris project the truth was exported from; its file stem names the truth files, and the spot script can fall back to reading it directly |

Optional: `stem` overrides the truth file stem. Unknown keys are kept in `Entry.extra` and ignored.

```json
{"entries": [
  {"tag": "n2rad51_noHS_herm_001", "split": "noHS",
   "results_dir": "C:/Users/ryane/germquant_golden/8666f07/n2rad51_noHS_herm_001",
   "truth_dir": "C:/Users/ryane/germquant_golden/imaris_truth/n2rad51_noHS_herm_001",
   "nd2": "C:/Users/ryane/n2_rad51_local/20251105_N2_noHS/20251105_N2_nohs_HERM _001.nd2",
   "ims": "C:/Users/ryane/n2_rad51_local/imaris/noHS_ims2/20251105_N2_nohs_HERM _001.ims"}
]}
```

## The truth file contract (consumed here, produced by the Imaris reader package)

For an .ims named `<stem>.ims` the truth folder holds:

* `<stem>__imaris_spots.csv`: `spot_id`, `x_um`, `y_um`, `z_um` (OUR frame, see below), `xi`, `yi`,
  `zi` (integer voxel indices, `floor(imaris_coord / imaris_voxel)`), `radius_um`, `quality` (NaN when
  absent), `spot_class` (string or empty; Imaris's Class A / Class B), `imaris_x_um`, `imaris_y_um`,
  `imaris_z_um` (the Imaris values).
* `<stem>__imaris_nuclei.csv`: `surface_id`, `label` (the mask value), `x_um`, `y_um`, `z_um` (our
  frame), `xi`, `yi`, `zi`, `volume_um3`, `n_voxels_mask` (NaN without a mask), `imaris_x_um`,
  `imaris_y_um`, `imaris_z_um`.
* `<stem>__imaris_labels.tif`: int32 ZYX label image the size of the image, 0 where no traced
  surface, written like the pipeline's own label images (OME PhysicalSize tags). Absent when no mask
  channel could be identified; the reason is then in the params JSON.
* `<stem>__imaris_params.json`: the voxel sizes of both frames, `ext_min`, `shape`, the
  CreationParameters text of every Spots and Surfaces object, the ROI box if found, which channel
  served as the mask and how it was identified, `n_spots`, `n_surfaces`, the spot classes present and
  the region-growing flag.

The harness reads the ROI box from the params JSON with `roi_box_from_params`: a `roi_box` (or `roi`,
`region_of_interest`) entry written as `{"min_xyz": [x, y, z], "max_xyz": [x, y, z]}` in voxel
indices, or a pair of triplets; failing that, the `mMinXYZ="x y z" mMaxXYZ="x y z"` text of the
Surfaces CreationParameters anywhere in the JSON. Voxel sizes are read from `voxel_zyx_um_ours` /
`nd2_voxel_zyx_um` and `voxel_xyz_um_imaris` / `imaris_voxel_xyz_um` (zyx or xyz spellings).

**The frame.** Imaris reports voxel CENTRES scaled by its rounded voxel size (0.108 um where the .nd2
says 0.108333), so matching in raw micrometres drifts by almost a micrometre across the 2600 pixel
field. The truth files therefore carry voxel indices computed with Imaris's own voxel size
(`floor(coord / voxel)`) and coordinates equal to index times the .nd2 voxel size, which is exactly how
the pipeline forms its own `x_um`/`y_um`/`z_um`. Our voxel size is taken from the pipeline tables'
`voxel_d{z,y,x}_um` columns, else the .nd2 header, else the params JSON (`per_gonad.csv` records the
source in `voxel_source`).

**Fallback path.** When the truth files are absent and `--truth` is `auto` (the default) or `ims`,
`benchmark_spots.py` builds the truth tables itself from the .ims (`truth_from_ims`, using
`germquant.validate.imaris_ims.read_imaris_ims` and the conversion above). That path has no label
image, no `quality` and no `spot_class`, so the detector-only and no-Class-B scores stay empty and the
traced subset uses the centroid rule; `truth_source` says `ims_fallback` for those gonads. The
nucleus benchmark has no fallback: it needs the truth label image.

## Running

```powershell
$env:CUDA_VISIBLE_DEVICES = "-1"
python scripts/benchmark_spots.py --manifest benchmark_manifest.json --out benchmark_spots_out
python scripts/benchmark_nuclei3d.py --manifest benchmark_manifest.json --out benchmark_nuclei3d_out
```

Options of `benchmark_spots.py`: `--split` (only entries with that split), `--tol-xy 0.3` and
`--tol-z 0.4` (the matching cylinder, um), `--surface-tol 2.5` (centroid distance for the traced
subset without a truth label image), `--overlap-frac 0.5` (volume fraction for the traced subset with
one), `--traced-rule auto|overlap|centroid|none`, `--max-candidates` and `--max-dense-cells` (the
Hungarian budgets, below), `--n-boot 1000`, `--seed 0`, `--truth auto|files|ims`, `--write-matches`
(adds `per_spot_matches.csv`).

Options of `benchmark_nuclei3d.py`: `--split`, `--tags`, `--group-by split|condition` (what the
per-gonad `split` column and the pooling by split hold), `--coverage auto|roi|bbox|none|protocol` (the
coverage region of the partial-annotation rule; `protocol` is the pre-registered rule, below),
`--dilate-um 2.0`, `--majority-frac 0.5` (splits and merges of the box rules), the protocol's parameters
`--min-truth-voxels 1000`, `--coverage-xy-um 1.0`, `--coverage-z-planes 1`, `--scope-frac 0.5` and
`--split-frac 0.2`, `--max-dense-cells`, `--n-boot`, `--seed`. Which prediction is scored: by default the
pipeline's own `<image>__nuclei_labels.tif` (for the goldens, the deployed Cellpose-SAM labels);
`--pred-dir DIR` scores `DIR/<tag>__labels.tif`; `--model MODEL.pt` runs the nuclei3d model (below) through
`segment_nuclei` on the cached DAPI channel (`--cache-root`, `--device`, `--tile`, `--min-volume-um3 4.0`,
and `--fg-threshold`, `--h-um`, `--neighbour-weight` to override the model's instance parameters) and
writes the labels to `<out>/labels/<tag>__labels.tif`, so `--pred-dir` can score them again. A test gonad
(the pre-registered split's, or one the manifest calls test) is skipped and listed under `skipped` unless
`--final` is given with a `--reason`.

A gonad whose outputs or truth are missing is skipped, not fatal: it is listed under `skipped` in
`pooled.json` with the reason. The exit code is 0 when at least one gonad was scored.

## The spot benchmark

**Expert spots into our nuclei.** Each expert spot is looked up in OUR nucleus label image at its
voxel (`zi`, `yi`, `xi`). Where that voxel is background, the most frequent non-zero label in a
neighbourhood of 1 voxel in z and 2 voxels in xy is taken (the prototype's rule), so a spot on the edge
of a nucleus still counts. The per-nucleus expert count is set against our own `n_spots`.

**Per-nucleus agreement.** Pearson's r, Lin's CCC and the bias (mean of ours minus expert, positive
when we over-count) over all germline nuclei (`in_germline`), and over the traced subset: with a truth
label image, our nuclei whose label has more than half of its voxels inside one truth label
(`--overlap-frac`); without one, our nuclei whose centroid lies within 2.5 um of a truth surface
centroid (`--surface-tol`). `traced_rule` in `per_gonad.csv` says which rule was used.

**Spot matching.** Our spots and the expert's are matched one to one inside an anisotropic cylinder:
dxy <= 0.3 um and |dz| <= 0.4 um by default (about three xy pixels and two z planes on the
0.108 x 0.108 x 0.2 um grid). Candidate pairs come from a KD-tree on the coordinates scaled by the
tolerances; the cost of a pair is `sqrt((dxy / tol_xy)^2 + (dz / tol_z)^2)`, so a pair at the edge of
the cylinder costs about 1. The matching is the Hungarian algorithm
(`scipy.optimize.linear_sum_assignment`) run on each connected component of the candidate graph
(`germquant.benchmark.matching.assign_sparse`): no candidate pair joins two components, so this is
exactly the global optimum, and the dense matrices stay tiny (a component is usually one or two spots
on each side) instead of 30000 by 30000. A component whose matrix would exceed `--max-dense-cells`
(4 million cells by default) is matched greedily (nearest pair first) and `n_greedy_components` counts
it; more candidate pairs than `--max-candidates` (5 million) switches the gonad to greedy matching.
`matching_method` records `hungarian`, `hungarian+greedy` or `greedy` per gonad, and
`matching_notes` says why.

**Three scores.**

* `precision`, `recall`, `f1`: our spots against all expert spots.
* `no_b_*`: without the expert's Class B spots, when the truth has a `spot_class` column. In these
  projects the class is Imaris's own filter on "Average Distance To 3 Nearest Neighbours" (Class B are
  the isolated spots, more than about 8 um from their neighbours). Class B spots and the spots we
  matched to them leave both counts, so finding them is neither rewarded nor punished.
* `det_*` (detector-only, with a truth label image): expert spots inside the traced nuclei against
  our spots inside the same traced nuclei, matched afresh. This removes our nucleus segmentation and
  the untraced part of the gonad from the comparison and scores the detector alone.

### Outputs of `benchmark_spots.py`

`per_gonad.csv`, one row per gonad:

| column | meaning |
|---|---|
| `tag`, `split` | from the manifest |
| `truth_source` | `truth files` or `ims_fallback` |
| `voxel_dz_um`, `voxel_dy_um`, `voxel_dx_um`, `voxel_source` | our frame's voxel size and where it came from |
| `n_nuclei_total`, `n_germline_nuclei` | our nuclei, and those with `in_germline` |
| `traced_rule`, `n_traced` | the traced-subset rule used and how many germline nuclei it marks |
| `n_expert_spots`, `n_expert_spots_class_b` | expert spots, and how many are Class B |
| `n_expert_spots_in_our_labels`, `n_expert_spots_in_germline_labels` | expert spots that land in any of our nuclei, and in a germline nucleus |
| `n_our_spots`, `n_our_spots_in_germline` | our spots, and those assigned to a germline nucleus |
| `mean_expert_all`, `mean_ours_all`, `pearson_r_all`, `ccc_all`, `bias_all` | per-nucleus agreement over all germline nuclei |
| `mean_expert_traced`, `mean_ours_traced`, `pearson_r_traced`, `ccc_traced`, `bias_traced` | the same over the traced subset |
| `n_matched`, `n_ours`, `n_expert`, `precision`, `recall`, `f1` | spot matching against all expert spots |
| `matching_method`, `n_candidate_pairs`, `n_greedy_components`, `matching_notes` | how the matching was solved |
| `has_spot_class`, `no_b_n_matched`, `no_b_n_ours`, `no_b_n_expert`, `no_b_precision`, `no_b_recall`, `no_b_f1` | the no-Class-B variant (empty without a class column) |
| `det_n_matched`, `det_n_ours`, `det_n_expert`, `det_precision`, `det_recall`, `det_f1` | the detector-only variant (empty without a truth label image) |
| `tol_xy_um`, `tol_z_um`, `surface_tol_um`, `overlap_frac` | the settings |
| `elapsed_s` | seconds for the gonad |

`per_nucleus.csv`, one row per germline nucleus: `tag`, `split`, `nucleus_id`, `centroid_x_um`,
`centroid_y_um`, `centroid_z_um`, `expert_n`, `ours_n`, `traced`.

`per_spot_matches.csv` (with `--write-matches`): `tag`, `our_spot_id`, `expert_spot_id`, `cost`,
`dxy_um`, `dz_um` for every matched pair of the all-spots matching.

`pooled.json`: `n_gonads`, `tags`, `n_expert_spots`, `n_our_spots`; `per_nucleus_all` and
`per_nucleus_traced` (`n_nuclei`, `n_gonads`, `mean_expert`, `mean_ours`, and `ccc`, `pearson_r`,
`bias` each as `{point, lo, hi, n_nan}`); `spots_all`, `spots_without_class_b` and
`spots_detector_only` (`n_gonads`, summed `n_matched`, `n_ours`, `n_expert`, and `precision`, `recall`,
`f1` each as `{point, lo, hi, n_nan}` from the summed counts, plus `mean_per_gonad_f1`);
`matching_methods` and `truth_sources` per tag; `settings`; `bootstrap` (`n_boot`, `seed`, unit, the
percentiles); `skipped`; `manifest`; and `by_split` with the same blocks per split when there is more
than one. NaN is written as `null`.

## The nucleus benchmark

**Overlaps.** A joint histogram of (our label, truth label) pairs is built plane by plane
(`nuclei3d.overlap_table`), together with every object's voxel count, z extent and centroid, so no
dense n_pred by n_truth matrix is formed and the two 2600 x 2600 x 71 images are never duplicated.
IoU is `intersection / (size_pred + size_truth - intersection)` for the pairs that share a voxel.

**Matching.** At each IoU threshold the pairs at or above it are matched one to one by the same
component-wise Hungarian algorithm (cost 1 - IoU), which maximises the number of matches first and the
summed IoU second.

**Partial-annotation rule.** A prediction that overlaps no truth object AND whose centroid lies outside
the truth coverage region is ignored: not a false positive, not a true positive. The region is the ROI
box from the params JSON when present (`coverage_rule` = `roi_box`), else the bounding box of all truth
labels dilated by 2 um (`truth_bbox_dilated_2um`); `--coverage bbox` forces the latter, `--coverage
none` switches the rule off. `n_pred_ignored_outside_coverage` and `n_pred_no_overlap_in_coverage` are
reported so the rule can be judged. Keep in mind that the ROI box holds nuclei the expert did not keep
(the truth is threshold-curated inside the box), so precision inside the box is a lower bound; recall
of the traced nuclei is the number to trust.

**The pre-registered rule (`--coverage protocol`).** `docs/BENCHMARK_PROTOCOL.md` (section 5.3), written
after the rule above, defines the scope differently, and a nucleus score reported under the protocol uses
this mode (`coverage_rule` starts with `protocol_`). Truth labels of fewer than 1000 voxels are removed
first (`n_truth_fragments_removed`, `truth_fragments_removed`; none of the eight train and val gonads has
one but n2rad51_noHS_male_000, whose 44-voxel label is removed); region C is the union of the remaining
truth masks dilated by 1.0 um in xy and one plane in z (an exact cylindrical dilation); a prediction is
in scope when at least half of its voxels lie in C, and out-of-scope predictions are ignored
(`n_pred_ignored_outside_coverage`), so the ROI box is never used. Splits and merges count in-scope
predictions that cover at least 20 percent of a truth nucleus's volume (`split_rule`
`in_scope_cover_0.2_of_truth`): a truth nucleus with two or more such predictions is split, a prediction
covering two or more truth nuclei is a merge. Matching, AP, recall and the z-extent error are the same as
above. On the deployed labels of the two val gonads the two rules disagree a lot: the ROI box of
n2rad51_HS_herm_001 holds 142 unmatched Cellpose nuclei that the protocol leaves out of scope (table
under "nuclei3d").

**Scores.** `f1_iou50` with `precision_iou50`, `recall_iou50`, `tp_iou50`, `fp_iou50`, `fn_iou50`
(TP are matched pairs, FP the counted predictions left unmatched, FN the truth objects left unmatched);
`ap_iou50_95`, the mean over the ten thresholds 0.50, 0.55, ..., 0.95 of TP / (TP + FP + FN) (the
DSB-2018 / StarDist definition; there are no confidence scores to rank by), with `tp_iouNN`,
`fp_iouNN`, `fn_iouNN` and `f1_iouNN` at every threshold; `mean_matched_iou_iou50`; `n_splits`: truth
objects that hold the majority (more than `--majority-frac` of the volume) of two or more predictions;
`n_merges`: predictions that hold the majority of two or more truth objects; `z_extent_error_mean_um`,
`z_extent_error_median_um`, `z_extent_error_mean_abs_um` over the pairs matched at IoU 0.5, where the
error is the predicted z span minus the truth z span and a span is `(zmax - zmin + 1) * dz`
(`n_z_pairs` pairs).

### Outputs of `benchmark_nuclei3d.py`

`per_gonad.csv`: `tag`, `split`, `n_pred`, `n_truth`, `n_pred_ignored_outside_coverage`,
`n_pred_counted`, `n_pred_no_overlap_in_coverage`, `coverage_rule`, `coverage_min_zyx`,
`coverage_max_zyx`, `n_splits`, `n_merges`, the per-threshold `tp_iouNN` / `fp_iouNN` / `fn_iouNN` /
`f1_iouNN`, `precision_iou50`, `recall_iou50`, `mean_matched_iou_iou50`, `ap_iou50_95`,
`matching_method`, `n_z_pairs`, the three `z_extent_error_*_um`, `voxel_dz_um`, `voxel_dy_um`,
`voxel_dx_um`, then (appended with the protocol mode) `z_extent_error_mean_slices`,
`z_extent_error_mean_abs_slices`, `split_rule`, `n_truth_fragments_removed`, `truth_fragments_removed`,
and from the script `voxel_source`, `roi_box`, `elapsed_s`, `pred_source`, `role` (plus `segment_s`,
`labels_file` and `segmentation_method` with `--model`).

`per_nucleus.csv`, one row per TRUTH nucleus: `tag`, `split`, `truth_label`, `truth_voxels`,
`truth_volume_um3`, `truth_z_span_um`, `matched_iou50`, `pred_label` (0 when unmatched), `iou`,
`pred_voxels`, `pred_z_span_um`, `z_extent_error_um`, `n_pred_majority_inside`, `is_split` (under the
active split rule), `n_pred_split_rule` (the count that rule uses), `z_extent_error_slices`.

`pooled.json`: `n_gonads`, `tags`, summed `n_pred`, `n_pred_counted`,
`n_pred_ignored_outside_coverage`, `n_truth`, `tp_iou50`, `fp_iou50`, `fn_iou50`; `precision_iou50`,
`recall_iou50`, `f1_iou50` and `ap_iou50_95` as `{point, lo, hi, n_nan}` from the summed counts;
`mean_per_gonad_f1_iou50`, `mean_per_gonad_ap`; `n_splits`, `n_merges`; `z_extent_error_um` (pooled
mean with interval) and `n_z_pairs`; `coverage_rules` and `matching_methods` per tag; `bootstrap`;
`settings`; `skipped`; `manifest`; `by_split`.

## Statistics (`germquant.benchmark.stats`)

`pearson`, `lin_ccc` (population moments, the convention of `germquant.validate.compare`) and `bias`
(ours minus expert) with `agreement_metrics` bundling them; `precision_recall_f1` from counts.

`cluster_bootstrap(groups, statistic, n_boot=1000, seed=0)` resamples whole gonads with replacement
and recomputes the pooled statistic each time; the 2.5 and 97.5 percentiles of the resampled values
are the interval. The gonad is the unit of replication because nuclei and spots inside one gonad share
its staining, focus and expert, so an interval over nuclei would be far too narrow. With few gonads
(seven tonight) the intervals are wide and a single gonad can move them; with one gonad the interval
collapses onto the point estimate. `pooled_agreement`, `pooled_detection` and `pooled_mean` are the
statistics the scripts use. Everything is deterministic for a given seed.

`paired_comparison(a, b)` compares two detectors scored on the same gonads: per-gonad differences of a
metric, the mean difference with a bootstrap interval, how many gonads favour each side, and the
Wilcoxon signed-rank p-value when at least five differences are non-zero. Use it on two `per_gonad.csv`
files from two runs (same manifest, same order):

```python
import pandas as pd
from germquant.benchmark.stats import paired_comparison
a = pd.read_csv("run_new/per_gonad.csv").set_index("tag")
b = pd.read_csv("run_old/per_gonad.csv").set_index("tag").loc[a.index]
print(paired_comparison(a["f1"], b["f1"], label_a="new", label_b="old"))
```

## Caveats

* The truth is partial and threshold-curated: inside the ROI box Imaris segmented automatically and
  the expert kept the surfaces that passed a Number-of-Voxels filter, so an untraced nucleus is not a
  false positive and a spot the expert did not keep is not necessarily noise. Recall-like numbers and
  the traced subset are the ones to read first; the pooled F1 is an audit, not a verdict.
* The expert's spots are threshold-curated automatic Imaris Spots, not hand-placed points, so a part
  of every disagreement is Imaris's detector against ours rather than an expert against ours.
* Two HS male projects were exported with spot region growing on, which moves the reported spot
  centre towards the intensity-weighted centre of the grown region; `region_growing` in the params
  JSON marks them.
* The Imaris mask channel gives every surface a distinct value, but value collisions have been seen
  (about 3 percent of nuclei on some gonads in the xlsx exports); the truth label image inherits the
  reader's handling of them.
* Everything is scored in the pipeline frame. A gonad whose tables lack voxel columns, whose .nd2 is
  unreadable and whose params JSON has no voxel size cannot be scored.

## First run on the goldens (2026-10-09)

`scripts/benchmark_spots.py` on the frozen outputs under `C:/Users/ryane/germquant_golden/8666f07/`
(pipeline commit 8666f07), truth by the .ims FALLBACK path (the truth files of the Imaris reader
package did not exist yet, so no Class B, no detector-only score, traced subset by the centroid rule),
cylinder 0.3 / 0.4 um, 1000 bootstrap resamples, seed 0. Every gonad was solved exactly by the
Hungarian algorithm (1057 to 1499 candidate pairs; no greedy fallback), in about one second per gonad
once the label image was loaded. Eight of the twelve gonads had outputs: `n2rad51_noHS_male_002` has
none (it floods the detector until a performance fix) and the three HS male gonads had not been run
yet; all four are listed under `skipped` in the pooled JSON. The outputs are not committed; the
command above reproduces them.

| gonad | germline nuclei | traced | expert spots | our spots | r (all) | CCC (all) | bias (all) | r (traced) | CCC (traced) | bias (traced) | P | R | F1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| noHS_herm_000 | 1243 | 107 | 1845 | 1934 | 0.945 | 0.743 | +0.71 | 0.912 | 0.609 | +2.37 | 0.501 | 0.525 | 0.512 |
| noHS_herm_001 | 595 | 146 | 1222 | 1653 | 0.942 | 0.783 | +1.07 | 0.912 | 0.675 | +2.37 | 0.553 | 0.748 | 0.636 |
| noHS_herm_002 | 814 | 93 | 1300 | 1521 | 0.916 | 0.786 | +0.61 | 0.851 | 0.657 | +1.28 | 0.596 | 0.697 | 0.642 |
| noHS_male_000 | 513 | 143 | 930 | 1440 | 0.944 | 0.767 | +1.22 | 0.872 | 0.563 | +3.11 | 0.488 | 0.755 | 0.592 |
| noHS_male_001 | 442 | 101 | 1256 | 1519 | 0.940 | 0.728 | +1.66 | 0.857 | 0.538 | +3.21 | 0.429 | 0.518 | 0.469 |
| HS_herm_000 | 607 | 88 | 5152 | 1552 | 0.889 | 0.864 | +0.32 | 0.887 | 0.838 | +0.90 | 0.562 | 0.169 | 0.260 |
| HS_herm_001 | 1167 | 95 | 2753 | 2600 | 0.836 | 0.643 | +1.34 | 0.822 | 0.423 | +3.93 | 0.326 | 0.308 | 0.317 |
| HS_herm_002 | 635 | 84 | 2364 | 1773 | 0.938 | 0.784 | +1.32 | 0.903 | 0.574 | +3.27 | 0.464 | 0.348 | 0.397 |

Bias is ours minus expert, spots per nucleus. Pooled over the eight gonads (6016 germline nuclei, 857
traced), intervals are the 2.5 and 97.5 bootstrap percentiles by gonad:

* per nucleus, all germline nuclei: CCC 0.758 [0.711, 0.805], r 0.898 [0.864, 0.939], bias +0.99
  [+0.74, +1.25] (expert 1.31, ours 2.30 spots per nucleus);
* per nucleus, traced subset: CCC 0.609 [0.534, 0.691], r 0.848 [0.808, 0.889], bias +2.58 [+1.96,
  +3.12] (expert 3.06, ours 5.65);
* spots: 6683 matched of our 13992 and the expert's 16822: precision 0.478 [0.416, 0.539], recall
  0.397 [0.278, 0.605], F1 0.434 [0.346, 0.549]; mean per-gonad F1 0.478;
* by split: noHS F1 0.566 [0.506, 0.631] (P 0.513, R 0.632), CCC all 0.764; HS F1 0.314 [0.260,
  0.397] (P 0.429, R 0.248), CCC all 0.751.

The matched pairs sit at a mean dxy of 0.014 um and dz of 0.037 um, and 82.5 percent of them are in
the same z plane: the index-wise frame conversion puts a matched spot on the expert's voxel. The
prototype (`rad51_ims_benchmark.py`, a 0.5 um sphere, greedy) gave the same per-nucleus r and CCC to
three decimals on the four gonads it covered and F1 within 0.002, so the tighter cylinder and the exact
matching change nothing there.

What the numbers say so far: the per-nucleus correlation is high on every gonad (r 0.84 to 0.95) but
we count about one spot per germline nucleus more than the expert, and 2.6 more per traced nucleus;
the HS gonads carry far more expert spots (2364 to 5152) than we detect (1552 to 2600), so recall drops
to 0.17 to 0.35 there. Whether that is the expert's threshold or our detector is what the detector-only
score is for, and it needs the truth label image. These are the first numbers of the harness, not a
calibration verdict.

## Detector bench harness

The scripts above score outputs that already exist. The detector bench harness produces new ones: it
runs a spot detection backend (SpotMAX with other settings, the experimental `log` backend, an
in-house detector later) on exactly the inputs the deployed spots stage saw, scores the result the same
way as the baseline, and keeps a record of every trial, which the tuning budget of
`docs/BENCHMARK_PROTOCOL.md` (section 6.4) counts. Three pieces:

* `scripts/bench_cache.py` writes the cache: the inputs of the golden runs' spots stage, one folder per
  gonad.
* `germquant.benchmark.detector`: `load_cache(tag)`, `run_backend(tag, backend, params)`,
  `score(tag, per_spot, per_nucleus)`, `pool_scores`, the test guard and the trial log helpers.
* `scripts/bench_detector.py` runs and scores every usable gonad of the train or the val split and
  appends the trial record.

`scripts/bench_reproduce_golden.py` is the check that the harness reproduces the deployed pipeline (below).

### The cache

`python scripts/bench_cache.py [--verify-sha256] [--tags ...] [--overwrite]` reads the benchmark
manifest (`--manifest`, default `C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/manifest.json`)
and, for every entry whose results folder holds the golden outputs (nuclei table, label image, spots
table, image summary, run manifest), writes `<cache root>/<tag>/` (`--cache-root`, default
`C:/Users/ryane/germquant_golden/bench_cache`). Entries without outputs are skipped and listed; a
complete cache is kept unless `--overwrite`.

| file | content |
|---|---|
| `foci.npy`, `dna.npy` | the foci and DNA channels, (Z, Y, X), full resolution, in the reader's dtype (uint16 for the lab's .nd2 stacks), read with `germquant.io.read.read_image`, the reader the pipeline uses. The channel indices come from the golden run's `channel_map_json` (its image summary), never from a constant. |
| `labels.npy` | int32 (Z, Y, X): exactly the label image the golden run handed to `detect_spots` |
| `nuclei.csv` | `nucleus_id`, `in_germline`, the three centroid columns and the golden `n_spots` of every row of the golden nuclei table |
| `meta.json` | written last (a folder without it is an unfinished cache and `load_cache` refuses it): `image_id`, `stem`, `spacing_zyx_um` (the reader's voxel size, exact), `shape_zyx`, `role_indices`, `channel_names`, `marker` and where it came from, `spots_config` (the golden run's `spots:` section, from its run_manifest.json), `nd2` (path, size, the sha256 of the split record and whether it was verified), `labels` (source, content, label count), `n_nuclei`, `n_germline`, `germline_label_ids`, `arrays` (file, dtype, shape and a sha256 of the data in C order for each array), the golden run's git sha, config hash, timestamp, tool versions, spot count and files, the role, condition and sex, and when and with which code the cache was made |

**Which label image.** In `src/germquant/pipeline.py`, at the pinned commit 8666f07 and on the current
branch alike, the spots stage calls `detect_spots(stack.data[spots_idx], labels, spacing, ...)` with
`labels` as the nucleus segmentation returned it. The germline stage only splits the rows of the nuclei
table (`in_germline`) and never edits the label image, and the run saves that same array as
`<image>__nuclei_labels.tif` (int32). So the spots stage saw every segmented nucleus, germline and
off-gonad alike, and `labels.npy` is that TIFF unchanged. The germline ids go into meta.json for
scoring; the off-gonad nuclei keep their spots, as they do in the pipeline.

Arrays are written through `numpy.lib.format.open_memmap` eight z planes at a time (`--z-chunk`)
under a temporary name that is renamed when complete. `--verify-sha256` hashes every .nd2 and refuses
one that is not the registered file. The test gonads are cached like the others: the cache holds
images and labels and no metric, and scoring a test gonad stays behind the guard below.

Built on 2026-10-10 with `--verify-sha256` (all eleven stacks matched their registered sha256): eleven
gonads, 40.94 GB on disk (40,939,186,156 bytes; 3.84 GB per no-heat-shock gonad of 71 planes, 3.62 GB
per heat-shock gonad of 67 planes), 11 to 15 s per gonad. `n2rad51_noHS_male_002` was skipped: its
golden outputs were still being frozen. `cache_index.json` at the cache root lists what each run did.

### Running and scoring (library)

```python
from germquant.benchmark import detector as bd

per_spot, per_nucleus = bd.run_backend("n2rad51_noHS_herm_000", "spotmax", {"effect_size_min": 2.5})
row, per_nuc, matches = bd.score("n2rad51_noHS_herm_000", per_spot, per_nucleus)
```

`run_backend` calls `germquant.spots.detect_spots` exactly as the spots stage does: the cached foci
channel and label image (read into memory), the voxel size, the marker, the nine SpotMAX keyword
arguments built from the `spots:` section with the stage's own expressions and code defaults
(`pipeline_spot_kwargs`), `backend` and `cfg_section`. `params` override keys of the golden section
first (`effective_section`); a nested mapping such as `log:` merges one level deep, the rule the
pipeline applies to a spot instance, so `{"log": {"min_snr": 3}}` changes one key of the log backend
and `{}` is the deployed configuration. `params` may not carry `backend`.

`score` reuses `germquant.benchmark.spots.benchmark_gonad` with the detector's spots and per-nucleus
counts in place of the golden ones, the cached label image as our nuclei and the truth files the
manifest lists. Every column of the spot benchmark's `per_gonad.csv` is there (matching against all
expert spots, `no_b_*` without Class B, `det_*` detector-only in the traced nuclei, the per-nucleus
agreement over all germline nuclei and over the traced subset), plus:

| column | meaning |
|---|---|
| `role` | train, val or test (the per-gonad `split` column holds the condition, noHS or HS, so `pooled.json` pools by condition under `by_split`) |
| `n_expert_spots_class_a` | expert spots that are not Class B |
| `class_a_mean_expert_all`, `class_a_pearson_r_all`, `class_a_ccc_all`, `class_a_bias_all`, and the same with `_traced` | the per-nucleus agreement with the expert count taken from Class A spots only; our count is unchanged (protocol section 5.1: the truth count changes with the class, the method's does not) |
| `backend`, `detect_s`, `score_s`, `n_spots_reconcile`, `n_unknown_nucleus_ids` | the backend, seconds spent detecting and scoring, whether the per-nucleus counts sum to the spot count, and nucleus ids the detector reported that the cached nuclei table does not hold (0 expected) |

`per_nucleus.csv` gains `role` and `expert_n_class_a`. `pool_scores` returns the blocks of
`spots.pool` with the cluster bootstrap by gonad plus `per_nucleus_all_class_a` and
`per_nucleus_traced_class_a`, overall and per condition.

### `scripts/bench_detector.py`

```powershell
$env:CUDA_VISIBLE_DEVICES = "-1"
python scripts/bench_detector.py --backend spotmax --params '{\"effect_size_min\": 2.5}' --split val --out <dir>
python scripts/bench_detector.py --backend log --params log_params.json --split train --jobs 2 --out <dir>
```

`--params` is a JSON object or the path of a JSON or YAML file (`{}` by default). `--split` is train or
val; `--tags` restricts the run to gonads of that split (a tag of another split is refused). Gonads
whose spot truth is unusable (`spot_truth_usable` false in the manifest: the two heat-shock males
exported with Imaris region growing) and gonads without a cache are skipped and listed. `--jobs N` runs
N gonads at once in separate processes (a SpotMAX worker peaked at 10.4 to 12.1 GB of working set on
these gonads, measured on 2026-10-10; the cached foci channel and labels alone are 2.9 GB), so pick N
by the free memory. Other options: `--manifest`, `--cache-root`, `--trials-log`,
`--no-trial-log`, `--n-boot` (1000) and `--seed` (20251105, the protocol's), `--tol-xy` and `--tol-z`
(0.3 and 0.4 um), `--verbose` (keeps SpotMAX's progress bars), `--detector-only` (a second run of the
backend with the truth masks as its nuclei, for the detector-only count of protocol section 5.1; see
"Leave-one-gonad-out") and `--note JSON` (a value stored as `note` in the trial record and pooled.json).

Outputs in `--out`: `per_gonad.csv`, `per_nucleus.csv`, `per_truth_nucleus.csv` (the section 5.1 partition,
"Leave-one-gonad-out"), `per_spot_matches.csv`, `pooled.json` (the pooled blocks above with `backend`, `params`,
`params_sha256`, `split`, `skipped`, `git_sha`, `manifest`, `cache_root` and the bootstrap settings) and
`detections/<tag>__spots.csv`, the detector's own spot table per gonad (with `--detector-only` also
`detections_detector_only/<tag>__spots.csv`).

**The trial log.** Every run appends one JSON line to `trials.jsonl` next to the manifest
(`C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/trials.jsonl`): `utc`, `git_sha` (with
`-dirty` when the tree had uncommitted changes), `backend`, `params`, `params_sha256` (sha256 of the
key-sorted compact JSON of `params`), `split`, `tags`, `skipped`, `final`, `pooled` (the seven pooled
blocks with their intervals), `pooled_by_condition`, `per_gonad` (the headline numbers of each gonad,
for paired comparisons), `bootstrap`, `out_dir`, `manifest`, `cache_root`, the versions of germquant,
SpotMAX, numpy, scipy, scikit-image and pandas, and `elapsed_s`. It is the record the protocol's
tuning budget (at most 60 configurations per competitor, each listed with its val score) is counted
from.

**The test guard.** `--split test`, or a test gonad among the selected ones, is refused unless
`--final` is given. The test gonads are those the pre-registration names
(`germquant.validate.benchmark_split.PREREGISTERED_ROLES`) and those the manifest calls test; a
manifest that disagrees with the pre-registration is refused outright. `--final` also needs `--reason`,
refuses a code tree with uncommitted changes unless `--allow-dirty`, refuses a third test scoring of the
same backend (protocol section 10), and appends its record to `final_runs.jsonl` instead of
`trials.jsonl`. The library refuses too: `run_backend` and `score` raise `TestSetGuardError` for a test
gonad unless called with `allow_test=True`, before anything is loaded, and `score` refuses a gonad with
unusable spot truth unless `include_unusable=True`.

### Harness validation

`scripts/bench_reproduce_golden.py --tags n2rad51_noHS_herm_000 n2rad51_HS_herm_000` ran SpotMAX with
each golden run's own `spots:` section from the cache and compared the result with the golden
`<image>__spots.parquet` (exact binary floats) and `<image>__spots.csv` (parsed with round-trip
precision) column by column in row order (spot_id, nucleus_id, marker, z_um, y_um, x_um, effect_size,
spot_mean_intensity, detector) and the per-nucleus counts with the golden nuclei table. Both gonads are
identical bit for bit, against both files: 1934 of 1934 spots (1274 nuclei) and 1552 of 1552 spots (662
nuclei), every column, every count, 85 to 97 s per gonad (on the uncommitted tree and again at commit
957e307). The two ran one after the other in one process, so nothing carries over between gonads.

Two wider checks came with the first runs below. The SpotMAX detections that `bench_detector.py` wrote
for all six train and val gonads (`detections/<tag>__spots.csv`) equal the golden spot tables column for
column. And the 50 columns that the bench `per_gonad.csv` shares with the earlier
`benchmark_spots.py` output (`C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/spots/per_gonad.csv`,
read for these six gonads only) are equal with zero tolerance, so a score from the harness is the
baseline's score. A second SpotMAX run of both splits gave the same numbers to the last digit.

### First runs (2026-10-10, commit 957e307)

Four runs, logged as lines 3 to 6 of `trials.jsonl` (lines 1 and 2 are the same SpotMAX runs made
before the code was committed, with the same numbers, logged with git sha `3923789...-dirty`). Outputs
under `C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/detector_bench/<backend>_<config>_<split>_957e307/`.
SpotMAX with `params {}` is the deployed configuration; `log` with `params {}` is the experimental
Laplacian-of-Gaussian backend with its code defaults (diameter 0.5 um, min SNR 2.0), never tuned on
these gonads. Train holds the four gonads with usable spot truth (noHS_herm_000, noHS_male_000,
noHS_male_001, HS_herm_000), val its two hermaphrodites (noHS_herm_001, HS_herm_001); the two
region-growing heat-shock males were skipped. Intervals: cluster bootstrap by gonad, 1000 resamples,
seed 20251105; with four and two gonads they are wide, and a per-condition number below rests on a
single heat-shock gonad.

| | SpotMAX, train | log, train | SpotMAX, val | log, val |
|---|---|---|---|---|
| germline nuclei / traced | 2805 / 425 | 2805 / 425 | 1762 / 208 | 1762 / 208 |
| per nucleus, germline, A+B: CCC | 0.783 [0.740, 0.831] | 0.861 [0.785, 0.926] | 0.695 [0.643, 0.783] | 0.926 [0.880, 0.957] |
| per nucleus, germline, A+B: bias | +0.87 [+0.58, +1.42] | +0.60 [+0.30, +1.11] | +1.25 [+1.07, +1.34] | +0.16 [-0.14, +0.73] |
| per nucleus, traced, A+B: CCC | 0.636 [0.548, 0.751] | 0.732 [0.547, 0.924] | 0.590 [0.436, 0.703] | 0.872 [0.821, 0.932] |
| per nucleus, traced, A+B: bias | +2.48 [+1.64, +3.12] | +1.73 [+0.80, +2.72] | +2.88 [+2.30, +3.76] | +0.72 [-0.17, +1.30] |
| per nucleus, germline, Class A: CCC | 0.783 [0.740, 0.831] | 0.858 [0.784, 0.920] | 0.676 [0.613, 0.783] | 0.922 [0.879, 0.949] |
| per nucleus, traced, Class A: CCC | 0.642 [0.549, 0.758] | 0.733 [0.546, 0.917] | 0.574 [0.411, 0.703] | 0.869 [0.821, 0.905] |
| spots A+B: F1 | 0.409 [0.295, 0.558] | 0.493 [0.388, 0.625] | 0.428 [0.317, 0.636] | 0.558 [0.427, 0.734] |
| spots A+B: precision / recall | 0.495 / 0.348 | 0.644 / 0.400 | 0.414 / 0.443 | 0.766 / 0.439 |
| matched / ours / expert | 3193 / 6445 / 9183 | 3672 / 5701 / 9183 | 1762 / 4253 / 3975 | 1747 / 2282 / 3975 |
| spots without Class B: F1 | 0.436 | 0.526 | 0.494 | 0.711 |
| spots detector-only (traced nuclei): F1 | 0.621 | 0.767 | 0.597 | 0.824 |
| detection time per gonad | 68 to 106 s | 14 to 29 s | 77 to 127 s | 18 to 33 s |

Per condition (`by_split` in pooled.json), germline CCC and bias, SpotMAX against log: train noHS
(three gonads) 0.752, +1.02 against 0.828, +0.77; train HS (HS_herm_000) 0.864, +0.32 against 0.938,
-0.03; val noHS 0.783, +1.07 against 0.880, +0.73; val HS (HS_herm_001) 0.643, +1.34 against 0.957,
-0.14. Spot F1 in the heat-shock gonads stays low for both (0.260 and 0.357 on HS_herm_000, 0.317 and
0.427 on HS_herm_001) because most heat-shock expert spots lie outside our germline nuclei (1334 of
5152 and 969 of 2753 fall in them); inside the traced nuclei the detector-only F1 is 0.636 against 0.916
and 0.432 against 0.830.

What this does and does not say:

* The untuned `log` defaults agree with the Imaris truth better than the deployed SpotMAX recipe on
  every pooled train and val number above, mostly by counting fewer spots per nucleus (smaller
  positive bias) at a higher precision. These are development numbers on six gonads; no test gonad was
  scored, and nothing here is a win condition of the protocol.
* There is no usable heat-shock male in train or val (both training ones were exported with region
  growing). On the two heat-shock hermaphrodites the deployed SpotMAX recipe OVER-counts per germline
  nucleus (+0.32 and +1.34) and `log` is close to zero (-0.03 and -0.14). The heat-shock undercount of
  the deployed recipe, the in-house detector's main target, is therefore not visible in train or val as
  they stand; seeing it during development needs the re-exported heat-shock male spot sets.
* Val holds hermaphrodites only, so anything selected on val is selected without a male.
* Imaris Spots detects blobs of a set diameter after a Gaussian background subtraction, which is close in
  spirit to a Laplacian of Gaussian, so a LoG detector may agree with it partly by construction;
  agreement with this truth is not yet agreement with a person (protocol section 9).

## net3d

`germquant.spots.net3d` is an EXPERIMENTAL learned RAD-51 foci detector: a compact residual 3D U-Net
that turns the RAD-51 and DAPI channels into a per-voxel heatmap of focus centres, trained on the lab's
Imaris spots. It is the spot backend `net3d` (never the default; `docs/BACKENDS.md` has its keys and the
two optional attributes it added to the backend contract), it is trained with `germquant train-spots`,
and it is scored with the detector bench harness above. Everything here is development on the train and
val gonads of `docs/BENCHMARK_PROTOCOL.md`. The training code never opens a test gonad (every gonad it
loads passes the test guard first) and no test gonad has been scored.

### What it learns from

* Train: the four train gonads with usable spot truth, n2rad51_HS_herm_000, n2rad51_noHS_herm_000,
  n2rad51_noHS_male_000 and n2rad51_noHS_male_001. The other two train gonads, the heat-shock males 000
  and 002, were exported with Imaris region growing and are left out (their spot truth is unusable).
  Training therefore sees one heat-shocked gonad, a hermaphrodite, and no heat-shocked male.
* Val (early stopping, model selection and the peak threshold): n2rad51_HS_herm_001 and
  n2rad51_noHS_herm_001, both hermaphrodites.
* No heat-shock male is in train or val. The heat-shock undercount of the deployed SpotMAX recipe, the
  in-house detector's main target, cannot be seen during development (see "First runs" above: on the
  two heat-shock hermaphrodites SpotMAX over-counts); it needs the re-exported heat-shock male spot sets.
* The roles come from the split file (`--split-file tests/golden/rad51_split.json`); a split file that
  moves a pre-registered gonad, or that disagrees with the manifest, is refused.

Inputs come from the bench cache (`foci.npy`, `dna.npy`, the golden run's `labels.npy` and germline ids)
and the truth files the manifest names. Each channel is normalised per gonad to
`(raw - p1) / (p99.9 - p1)`, with the percentiles of the whole volume on every fourth voxel in y and x,
so the input never depends on a segmentation. On these stacks p1 is 87 to 89 in both channels and p99.9
of the RAD-51 channel is 208 to 343 on the train gonads, 149 (HS_herm_001) and 223 (noHS_herm_001) on val.

### Targets, loss region and weights (`net3d/data.py`)

* Target: a Gaussian of peak 1 at every Class A spot (sigma 0.25 um in z and 0.12 um in xy, about 1.25
  and 1.1 voxels, the size of the confocal image of a focus; overlapping Gaussians take their maximum),
  centred on the Imaris position converted to the image grid (`(coordinate - ExtMin) / imaris_voxel -
  0.5`, whose rounding is the truth files' `zi, yi, xi`).
* Loss region: the golden run's germline nucleus labels dilated by 0.3 um (exact ellipsoidal dilation).
  The truth spots were detected over the whole field with a per-gonad Quality threshold set by eye
  (protocol section 4.1), and the counts that matter are per nucleus; outside this region the weight
  is 0 and the network is unconstrained there (it does put low probabilities on background, which the
  nucleus assignment drops; see the peak finding below).
* Inside it: weight 1 within the traced truth coverage (the protocol's region C: the truth masks dilated
  by 1.0 um in xy and one plane in z) and 0.5 elsewhere in the germline, where the same automatic
  detection holds but nobody looked nucleus by nucleus. The traced coverage is about 0.5 to 0.65 percent of a
  stack's voxels, the rest of the germline region 1.2 to 4.5 percent.
* Class B spots (Imaris's "isolated" class, not false positives) get weight 0 inside an ellipsoid of
  0.5 um in z and 0.3 um in xy, so the network is neither taught to find nor to reject them.
* Every weight is multiplied by `1 + 9 * target`, which spends the gradient near the foci without moving
  any voxel's optimum (each voxel's weighted cross-entropy is still minimised by its own target). The
  loss is the weighted binary cross-entropy of the logits, divided by the batch's summed weight.
* Patches of 32 x 128 x 128 voxels (6.4 x 13.9 x 13.9 um), eight per step: a gonad drawn uniformly,
  then half the time a Class A spot inside the loss region plus up to a third of a patch of offset,
  otherwise a germline nucleus centroid plus up to half a patch; a patch may hang over the stack's edge
  by a quarter (zero input, zero weight there). Augmentation: y and x flips, a y/x transpose and a
  log-uniform gain per channel (RAD-51 0.75 to 1.33, DAPI 0.8 to 1.25). Each sample draws from its own
  generator seeded with (seed, epoch, index).

### The network (`net3d/model.py`)

Four levels of one residual block each (two 3 x 3 x 3 convolutions, batch normalisation, SiLU; 16, 32,
64 and 128 channels), an xy-only first pooling (1 x 2 x 2, which makes the next level roughly
isotropic at 0.2 x 0.217 x 0.217 um) and 2 x 2 x 2 after it, transposed-convolution upsampling, a 1 x 1 x 1
head whose bias starts at -4; 1,421,889 parameters. Pooling is a reshape and a max and upsampling a
transposed convolution with kernel equal to stride, so the CUDA backward pass is deterministic. Batch
normalisation (an affine map in eval mode) keeps a voxel's output independent of how the volume is
tiled. The receptive radius is at most 30 voxels in z and 58 in y and x (`receptive_radius`), and a test
checks that a perturbed voxel never reaches farther.

### Training and model selection (`net3d/train.py`)

AdamW (learning rate 1e-3, weight decay 1e-4), 200 warm-up steps then a cosine decay to 5 percent over
the maximum of 40 epochs of 200 steps, bf16 autocast on CUDA, seeds for python, numpy and torch,
`torch.use_deterministic_algorithms(True)` with cuDNN benchmarking off and `CUBLAS_WORKSPACE_CONFIG`
set. Four launches with the same seed gave bit-identical losses and validation numbers for every epoch
they shared.

After every epoch the val gonads go through the backend's own path (tiled float32 inference, peak
finding, nucleus assignment in the cached labels) and are scored with the bench harness's functions:
expert spots mapped into our labels, Lin's CCC of the per-nucleus counts over the traced subset
(Class A+B, the protocol's primary partition as the harness computes it: our nuclei more than half
inside a truth mask) and over all germline nuclei, and the Hungarian spot matching in the 0.3 um by
0.4 um cylinder for F1. The expensive parts (the traced-subset overlap, the expert spots' labels) are
computed once per gonad; at the first validation the numbers are compared with
`germquant.benchmark.detector.score` itself and training stops if they differ (the difference was 0 in
every run). The two val gonads are pooled (nuclei concatenated, spot counts summed) at each of 18
thresholds, 0.05 to 0.90.

Selection (`select_config`): every scored (epoch, threshold) whose pooled traced CCC is within 0.01 of
the best one so far is eligible, and the eligible configuration with the best pooled F1 wins (then the
higher CCC, the earlier epoch, the lower threshold). This rule is a declared deviation from protocol
section 6.4 (protocol changelog entry 3), adopted on val after the first launches; the 0.01 band
borrows the scale of the section 10 stop rule and is not a tie band the protocol defines. The first launches used a strict CCC maximum with F1 only for
exact ties; it jumped between thresholds 0.05 and 0.60 on differences of 0.001 (it picked epoch 5 at
threshold 0.05, germline CCC 0.959, F1 0.608, over epoch 4 at 0.45 with traced CCC 0.979 and F1 0.627),
so the rule was changed before the run reported below; this was decided on val, never on test. Early
stopping: 8 validations in a row without a change of the selected configuration. Every epoch's weights
are kept (`checkpoints/epoch_NNN.pt`, sha256 in the log), and `model.pt` holds the selected epoch's
weights with the selected threshold, which the backend uses unless `spots.net3d.threshold` says otherwise.

### Peak finding and assignment (`net3d/peaks.py`, `net3d/infer.py`)

Inference runs in float32 on tiles whose origins are multiples of the network's size divisor and whose
overlap is the receptive radius, each tile contributing only its core: tiled and untiled heatmaps agree
to float rounding (tested), tiles without a nucleus are skipped. Peaks are the voxels above the
threshold that are the maximum over an ellipsoid of 0.6 um in z and 0.3 um in xy, settled by a greedy
ellipsoidal suppression (for these local maxima only exact ties are left to settle), refined to
sub-voxel positions by the parabola through the log probabilities on each axis, and assigned to the label
at the refined position or, failing that, to the nearest label within 0.3 um (`contested` when two labels
are that close). Peaks with no label within reach are dropped.

The network's output outside the loss region is unconstrained, and on the heat-shock val gonad at epoch 3
of a training run, 47.6 million voxels (10 percent of the stack) were above 0.05, with 1.93 million
local maxima, of which 1269 were assigned at that threshold. The peak finder handles this exactly: a
3 x 3 x 3 maximum filter prefilters dense blocks (28 s instead of 49 s for the identical candidate
set), suppression only loops over candidates that have a partner, and a coarse grid of blocks within
reach of a label drops the far peaks before refinement and assignment (3210 instead of 240015 peaks to
assign, the identical 1269 assigned). Three launches of the same seeded run, with the original peak
finder, with the pair-only suppression and a coarse grid, and with the final version, gave the same
validation numbers at every threshold, pooled and per gonad: epochs 1 to 4 between the first two, and
epochs 1 to 10 between the last two (8280 values, none different), with bit-identical losses. The
interrupted launches' training logs are kept in
`C:/Users/ryane/germquant_golden/models/net3d_20261010/earlier_launches/` (launch 1 stopped after one
epoch to change the GPU cache handling, launches 2 and 3 for the peak finder and the selection rule). At a
working threshold (0.3) the whole chain takes about 3 s per gonad.

Timing on the RTX 5090 (shared with other jobs): 39 to 44 s per training epoch, about 4.3 GB of
allocated GPU memory in training and about 12.5 GB reserved during the validation inference; inference
of one val gonad takes 10 to 25 s. On the CPU (4 threads, competing with other jobs) the inference of
one heat-shock stack took 566 s.

### Running it

```powershell
$env:PYTORCH_CUDA_ALLOC_CONF = "expandable_segments:True"
germquant train-spots --cache C:/Users/ryane/germquant_golden/bench_cache `
    --split-file tests/golden/rad51_split.json --out C:/Users/ryane/germquant_golden/models/net3d_20261010 `
    --epochs 40 --seed 0 --device cuda --min-free-gb 12 --prep-dir C:/Users/ryane/germquant_golden/models/net3d_prep
python scripts/bench_detector.py --backend net3d --params net3d_val_params.json --split val --out <dir>
```

where `net3d_val_params.json` holds `{"net3d": {"model": ".../model.pt", "device": "cuda"}}`. Other
options: `--patch Z,Y,X` (multiples of 4, 8, 8), `--batch-size`, `--steps-per-epoch`, `--lr`,
`--patience`, `--val-every`, `--no-dna` (the RAD-51-only ablation), `--thresholds`, `--tile`,
`--select-tolerance`, `--base-channels`, `--levels`, `--manifest`, `--wait-minutes`, and for retraining a
selected model on other gonads (protocol sections 6.1 and 6.2; "Leave-one-gonad-out") `--recipe-from MODELDIR`
(its settings, seed and schedule, its selected epoch count and threshold, no validation), `--train-tags`
(train on exactly these train or val gonads, in this order), `--val-tags` (validate on these; the option alone
for none), `--stop-after-epochs` and `--threshold` (a run's epoch count within the schedule, and the threshold
stored for a run without val gonads). `--min-free-gb`
waits (polling every 60 s) until nvidia-smi reports that much free memory; on Windows,
`torch.cuda.mem_get_info` reports the whole card as free while other processes hold most of it, so it is
not used for this. The loss-region maps (one uint8 volume per train gonad, 0.45 GB each) are computed
once into `--prep-dir` and reused by later runs. Outputs in MODELDIR: `model.pt`, `training_log.json`
(settings, gonads, every epoch's loss, timing, memory, checkpoint sha256 and the full validation sweep,
the selection), `model_card.json` (training, val and excluded gonads, truth classes, epochs, the
selection, the val metrics, sha256, git sha, limitations) and `checkpoints/`.

### First training run and its val scores (2026-10-10)

The model: `C:/Users/ryane/germquant_golden/models/net3d_20261010/model.pt`, sha256
`8cdb710bd2f0fe7acef49c6117d7113ad61e224405960d2847b34c8327a4c4f7`, trained by commit 87d77f1 with seed
0 (17:19 to 17:40 UTC, 1268 s on the RTX 5090). It ran 12 epochs and stopped early (the selection did
not change in 8 validations); the selected configuration is epoch 4 at threshold 0.40, out of 216 scored
(12 epochs times 18 thresholds; 32 were within 0.01 of the best traced CCC). At epoch 4 the learning rate
was still 0.99 of its peak (the cosine runs over the maximum of 40 epochs).

Scored with `scripts/bench_detector.py --backend net3d` on the val split (trial log line 79, commit
5181bed, CUDA; outputs in
`C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/detector_bench/net3d_20261010_val_5181bed/`).
The harness's numbers are the training-time validation numbers exactly (1909 matched of our 2089 and
the expert's 3975 spots in both). SpotMAX (golden config) and `log` (defaults) are lines 3 and 6, the runs
of the "First runs" section above, on the same two gonads. Intervals: cluster bootstrap by gonad, 1000
resamples, seed 20251105; with two gonads they only show how far apart the two gonads are.

| val, pooled over 2 gonads | SpotMAX (golden config) | log (defaults) | net3d |
|---|---|---|---|
| per nucleus, germline, A+B: CCC | 0.695 [0.643, 0.783] | 0.926 [0.880, 0.957] | 0.982 [0.981, 0.982] |
| per nucleus, germline, A+B: bias | +1.25 [+1.07, +1.34] | +0.16 [-0.14, +0.73] | +0.04 [+0.03, +0.08] |
| per nucleus, traced, A+B: CCC | 0.590 [0.436, 0.703] | 0.872 [0.821, 0.932] | 0.977 [0.974, 0.977] |
| per nucleus, traced, A+B: bias | +2.88 [+2.30, +3.76] | +0.72 [-0.17, +1.30] | -0.00 [-0.11, +0.06] |
| per nucleus, germline, Class A: CCC | 0.676 | 0.922 | 0.973 |
| per nucleus, traced, Class A: CCC | 0.574 | 0.869 | 0.969 |
| spots A+B: F1 | 0.428 [0.317, 0.636] | 0.558 [0.427, 0.734] | 0.630 [0.494, 0.854] |
| spots A+B: precision / recall | 0.414 / 0.443 | 0.766 / 0.439 | 0.914 / 0.480 |
| matched / ours / expert | 1762 / 4253 / 3975 | 1747 / 2282 / 3975 | 1909 / 2089 / 3975 |
| spots without Class B: F1 | 0.494 | 0.711 | 0.805 |
| spots detector-only (traced nuclei): F1 | 0.597 | 0.824 | 0.945 |
| detection time per gonad | 77 to 127 s (CPU) | 18 to 33 s (CPU) | 20 to 28 s (CUDA) |

| val gonad: germline CCC / traced CCC / F1 (our spots, expert spots) | SpotMAX | log | net3d |
|---|---|---|---|
| n2rad51_noHS_herm_001 | 0.783 / 0.703 / 0.636 (1653, 1222) | 0.880 / 0.821 / 0.734 (1453, 1222) | 0.981 / 0.975 / 0.854 (1065, 1222) |
| n2rad51_HS_herm_001 | 0.643 / 0.436 / 0.317 (2600, 2753) | 0.957 / 0.932 / 0.427 (829, 2753) | 0.982 / 0.974 / 0.494 (1024, 2753) |

The net3d per-nucleus bias is +0.08 (germline) and +0.06 (traced) on the no-heat-shock gonad and +0.03
and -0.11 on the heat-shock one. Of its spots, 1035 and 985 lie inside a label and 30 and 39 were
assigned by the nearest-label fallback; 13 and 65 are contested.

How much the selection on val could flatter it: at the selected epoch every threshold from 0.25 to 0.55
gives a traced CCC of at least 0.963 and an F1 of at least 0.615; at threshold 0.40 every epoch from 3
to 11 gives a traced CCC of 0.94 to 0.98 and an F1 of 0.58 to 0.63. Epoch 12 fell to 0.887 at 0.40 (0.932
at its own best threshold, 0.10), and the best threshold of an epoch moved between 0.10 and 0.60 from
epoch to epoch.

**Device check.** The same model on the CPU (4 threads, `--no-trial-log`, outputs in
`.../detector_bench/net3d_20261010_val_cpu_check_5181bed/`) gives the identical spots on both gonads
(1065 and 1024, the same nuclei, positions within 6e-7 um and probabilities within 1e-6 of the CUDA
run) and so the identical scores; detection took 288 s and 537 s per gonad against 20 s and 28 s.

**Ablation without DAPI.** The same recipe and seed with the RAD-51 channel only (`--no-dna`;
`C:/Users/ryane/germquant_golden/models/net3d_nodna_20261010/model.pt`, sha256
`7eb7563decaf3d3f3178c716eda2a9118ff7714dba1547aa1c9e908c7a5d6e7f`, commit 5181bed, whose training code
is that of 87d77f1; 16 epochs, 1159 s) selected epoch 8 at threshold 0.10 out of 288 configurations. On
val (trial log line 96, commit 9ec3f3b): germline CCC 0.983 [0.977, 0.986] with bias +0.04, traced CCC
0.980 [0.974, 0.983] with bias +0.05, Class A CCC 0.975 (germline) and 0.973 (traced), F1 0.635 [0.498,
0.853] (precision 0.925, recall 0.483; 1921 matched of 2076 and 3975), without Class B 0.813,
detector-only 0.951; per gonad 0.977 / 0.974 / 0.853 (1105 spots) without heat shock and 0.986 / 0.983 /
0.498 (971 spots) with it. Every CCC and F1 is within 0.01 of the two-channel model's (the precision, 0.925
against 0.914, differs by 0.011).

What this does and does not say:

* On the two val hermaphrodites net3d agrees with the Imaris truth far better than the deployed SpotMAX
  recipe and the untuned `log` baseline on every pooled number: a per-nucleus CCC of 0.98 over the
  germline nuclei and over the traced subset, with biases of +0.04 and -0.005 spots per nucleus (SpotMAX
  +1.25 and +2.88, `log` +0.16 and +0.72), a spot precision of 0.91 and a detector-only F1 of 0.945.
* These are development numbers. net3d's epoch and threshold were chosen on these two gonads among 216
  configurations; SpotMAX is the deployed configuration and `log` its never-tuned defaults (other `log`
  settings are being tried on train, trial log lines 7 onward). The plateau above suggests the optimism
  of the selection is a few hundredths of CCC, well below the gap to the baselines, but only the test
  set can show it. No test gonad was opened.
* The heat-shock undercount, the in-house detector's main target, is not tested: there is no heat-shock
  male in train or val. On the heat-shock hermaphrodite net3d counts 1024 spots where the expert has 2753
  over the whole field; per germline nucleus its bias is +0.03, and most of the other expert spots lie
  outside our nuclei (969 of the 2753 fall in germline labels, "First runs" above), which caps the
  whole-field recall (0.339) and F1 (0.494) of any detector that counts inside nuclei.
* DAPI brings no measurable gain on val: the RAD-51-only network matches the two-channel one within
  0.01 on every CCC and F1. With two val gonads that means no gain could be seen, not that there is none.
* The selected epoch is early (4 of 12) with the learning rate at its peak, and each epoch's best
  threshold moved between 0.10 and 0.60, so the network's calibration was not settled from epoch to
  epoch; the ablation picked threshold 0.10. A schedule that anneals sooner is the obvious next run; it was
  not tried here.
* A learned detector imitates its truth, including the expert's per-gonad Imaris Quality thresholds
  (6.8 to 102.3), which the per-gonad normalisation lets it partly mimic. Agreement with Imaris is not
  agreement with a person (protocol section 9).
* Before any test scoring, the protocol (section 6.1) asks for the selected configuration to be retrained
  on train plus val; that, the leave-one-gonad-out analysis of section 6.2 and the changelog entry of
  section 10 have not been done. (Package E1 did the first two, "Leave-one-gonad-out" below; the
  pre-scoring entry of section 10 is step 4.2 of `docs/FINAL_SCORING.md`.)

## Tuning on train and val

Package D3, 2026-10-10. The two classical detectors that the bench harness runs were tuned on the train
gonads and chosen between on the val gonads: the experimental `log` backend, and SpotMAX, the deployed
detector, re-tuned as the competitor that section 8 of `docs/BENCHMARK_PROTOCOL.md` requires. No test
gonad was run, scored or looked at. Every parameter set below is a record in `trials.jsonl` (split train
or val, with a `tuning` block naming its search and stage), and every output folder is under
`C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/tuning_d3/`.

### Rules, fixed before the first trial

The plan was written to `tuning_d3/search_plan.json` before any D3 trial ran.

* Gonads. Train: the four with usable spot truth, n2rad51_noHS_herm_000, noHS_male_000, noHS_male_001
  and HS_herm_000 (2805 germline nuclei, 425 traced, 9183 expert spots). Val: n2rad51_noHS_herm_001 and
  HS_herm_001 (1762 germline nuclei, 208 traced, 3975 expert spots). The two heat-shock males of train,
  HS_male_000 and HS_male_002, were skipped because their spots were exported with Imaris region growing.
  **No heat-shock male took part in any choice below, and val holds hermaphrodites only**, so the final
  choices were made without seeing any male on val and without seeing a heat-shock male at all.
* Objective. The pooled per-nucleus Lin's CCC over the germline nuclei of the split, the truth counted
  from Class A+B spots and the method counted end-to-end (the detector's own spots and per-nucleus counts
  on the deployed nucleus masks): `per_nucleus_all` in pooled.json, the harness's headline per-nucleus
  number. It is rounded to 3 decimals and the pooled spot F1 (Class A+B, 0.3 um in xy by 0.4 um in z)
  breaks ties. The traced-subset CCC, the Class A numbers, the F1 without Class B and the detector-only
  F1 are reported and steered nothing. The harness's traced subset (our nuclei with more than half of
  their volume inside one Imaris surface, counted in our labels) is not the primary partition of
  protocol section 5.1, which takes the Imaris masks themselves as the nuclei; neither was used to select.
  (Since E1 the harness computes the section 5.1 partition as well; see "Leave-one-gonad-out".)
* Budgets. At most 60 parameter sets for `log` (search `d3_log`) and 24 for SpotMAX (search
  `d3_spotmax`); `bench_tune.py --budget` refuses a run that would exceed them.
* Selection. The best 5 `log` sets on train went to val and the best of the 5 on val by the same
  objective is the final `log` configuration; the best 3 SpotMAX sets likewise. The deployed SpotMAX
  recipe is the reference on val.

### The runner: `scripts/bench_tune.py`

`bench_detector.py` runs one parameter set per call and reloads each gonad's cache (about 3 GB) every
time. `bench_tune.py` (library `germquant.benchmark.tuning`) takes a list of parameter sets (a JSON or
YAML file of `{name, params}`, `params` exactly as `bench_detector.py --params`), keeps each gonad in one
worker process that loads its cache once, runs every set of its share with `detector.run_backend`, then
releases the foci channel, loads the truth once and scores every set with `detector.score`. Each set is
pooled with `detector.pool_scores`, written to `<out-root>/<name>/` with the files of
`bench_detector.py` (per_gonad.csv, per_nucleus.csv, per_spot_matches.csv, pooled.json,
detections/), and logged with `detector.trial_record` plus a `tuning` block (search, stage, name,
objective, whether the detection was shared, wall time). `<out-root>/summary.csv` and `summary.json`
rank the sets by the objective. `--top-from SUMMARY ... --top K` takes the K best sets of earlier
summaries (how the finalists went to val); `--budget N` with `--search` counts the distinct parameter
sets already logged under that search; `--jobs` sets the worker processes and `--chunks` the tasks per
gonad. There is no `--final`: a test gonad is refused by the script and by `run_configs`.

```powershell
$env:CUDA_VISIBLE_DEVICES = "-1"
python scripts/bench_tune.py --backend log --configs logA.json --split train --search d3_log --stage A `
    --budget 60 --jobs 4 --out-root <dir>
python scripts/bench_tune.py --backend log --top-from <A>/summary.json <B>/summary.json --top 5 `
    --split val --search d3_log --stage val --out-root <dir>
```

**Shared SpotMAX detections.** `detect_spots_spotmax` applies `effect_size_min` as its last step, a row
filter after detection, the features and the z-split merge (a missing effect size passes; nothing is
filtered when no spot has one). Sets that differ only in `effect_size_min` therefore share everything
before it: the runner makes one SpotMAX run per such group with `effect_size_min` 0 and applies each
set's minimum afterwards (`tuning.filter_effect_size`, which rebuilds both tables as
`detect._finish_tables` does). The equality is tested on synthetic data, and on the real gonads the
deployed recipe obtained this way (triangle, 0.3 um, 0.08 um, filtered at 3.0) is bit for bit the
golden spots table, every column and every per-nucleus count, on all six train and val gonads (1934,
1440, 1519 and 1552 spots on the train gonads, 1653 and 2600 on the val gonads; `detector.compare_with_golden`).
Those runs also give the same 66 per-gonad columns as D1's `bench_detector.py` runs of the deployed recipe
(trials 3 and 4), and the `log` run with the code defaults gives D1's detections byte for byte on the four
train gonads, so the `pre_gate_sigmas` key left the default path unchanged on real data. `--no-shortcut`
runs every set on its own.

### Trial counts

88 trial records, all with split train or val: `log` 55 parameter sets on train (budget 60) and its 5
finalists on val; SpotMAX 24 parameter sets on train (budget 24) and its 3 finalists plus the deployed
recipe on val. Wall time about 56 minutes for `log` and 66 minutes for SpotMAX, run with at most four
worker processes in all (other jobs shared the workstation; the SpotMAX stages used one or two workers
to stay within memory). No configuration was cut. All records carry git sha `1a95863...-dirty`: the
uncommitted changes at the time were this documentation, the tests and the profiles; the detector,
harness and runner code that ran is the committed 1a95863. The parameter lists of every stage are in
`tuning_d3/configs/` (logA, logB, logC, smS1, smS2, smVal); the outputs are in `tuning_d3/log_A_train`,
`log_B_train`, `log_C_train`, `log_val`, `spotmax_S1_train`, `spotmax_S2_train` and `spotmax_val`, each
with its `summary.csv`, and the console logs sit next to them.

### `log` on train

| stage | what varied (everything else at the code defaults) | sets | wall time |
|---|---|---|---|
| A | `diameter_um` 0.35, 0.42, 0.5, 0.6, 0.7 by `min_snr` 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 | 30 | 18.5 min, 4 workers |
| B | at the stage-A best (0.5 um, SNR 2.5): `min_separation_um` 0.25 (the default, half the diameter), 0.375 and 0.5 by `pre_gate_sigmas` 2, 3, 4 by `merge_z_columns` on and off, less the stage-A point | 17 | 23.2 min, 2 workers |
| C | at the best of A and B (pre-gate 2, merge on): `diameter_um` 0.46, 0.5, 0.54 by `min_snr` 2.25, 2.5, 2.75, less the centre | 8 | 9.2 min, 2 workers |

Stage A, pooled germline CCC / bias on train (bias is ours minus expert, foci per germline nucleus; 0.5
um with SNR 2.0 is the code default):

| `diameter_um` \ `min_snr` | 1.5 | 2.0 | 2.5 | 3.0 | 3.5 | 4.0 |
|---|---|---|---|---|---|---|
| 0.35 | 0.113 / +14.09 | 0.462 / +3.09 | 0.640 / +1.30 | 0.769 / +0.35 | 0.775 / -0.26 | 0.651 / -0.66 |
| 0.42 | 0.131 / +8.20 | 0.799 / +0.90 | 0.903 / +0.12 | 0.864 / -0.28 | 0.745 / -0.60 | 0.592 / -0.82 |
| 0.5 | 0.161 / +6.04 | 0.861 / +0.60 | **0.934 / -0.04** | 0.879 / -0.36 | 0.760 / -0.62 | 0.622 / -0.80 |
| 0.6 | 0.213 / +4.39 | 0.854 / +0.46 | 0.923 / -0.12 | 0.869 / -0.39 | 0.754 / -0.62 | 0.626 / -0.80 |
| 0.7 | 0.272 / +3.13 | 0.836 / +0.31 | 0.895 / -0.20 | 0.839 / -0.44 | 0.730 / -0.65 | 0.622 / -0.80 |

* The SNR gate is the parameter that matters. At SNR 1.5 every diameter floods (13,503 to 44,842 spots
  against the expert's 9183); the per-nucleus bias falls through zero near SNR 2.5 and is negative from
  3.0 up at every diameter but 0.35 um. 0.5 um is the best diameter at every SNR from 2.0 to 3.0.
* The flood guard (`max_spot_candidates`, 30000 candidates) fired in 30 stage-A runs: every gonad at 0.35
  um and noHS_herm_000 at 0.42 um. Those counts are floors in the pipeline's sense.
* Stage B: at SNR 2.5 the pre-gate and the z merge change the spot counts by less than 1 percent (the
  SNR gate removes what they let through; the six sets at the default separation are within 0.0005 in
  CCC); a separation of 0.75 or 1.0 diameters lowers the CCC to 0.932 and 0.920.
* Stage C: no neighbour beats the centre (the best is 0.931, at 0.54 um).
* By condition, at the best setting the three gonads without heat shock are nearly unbiased (+0.12 per
  germline nucleus) and the heat-shock hermaphrodite is undercounted (-0.61). Its bias reaches zero near
  SNR 2.0, theirs between 2.5 and 3.0: one global threshold cannot fit both. The Imaris Quality
  thresholds the annotator set were lower for the heat-shock projects (24.2 and 17.4 against 42 to 69).

The five best of the 55 on train (all one point of the search, 0.5 um and SNR 2.5, with different
pre-gate and merge settings), with the untuned defaults and the deployed SpotMAX recipe for reference:

| set | germline CCC [95% CI] | bias | traced CCC | traced bias | germline CCC, Class A | spot F1 [95% CI] | P / R | F1 without B | detector-only F1 | matched / ours / expert |
|---|---|---|---|---|---|---|---|---|---|---|
| logB_d0.5_s2.5_k0.5_g2_m1 | 0.934 [0.892, 0.971] | -0.04 | 0.936 | +0.15 | 0.934 | 0.516 [0.347, 0.753] | 0.869 / 0.367 | 0.561 | 0.882 | 3369 / 3877 / 9183 |
| logA_d0.5_s2.5 | 0.934 [0.891, 0.971] | -0.04 | 0.935 | +0.14 | 0.934 | 0.516 [0.347, 0.754] | 0.870 / 0.367 | 0.561 | 0.883 | 3367 / 3870 / 9183 |
| logB_d0.5_s2.5_k0.5_g4_m1 | 0.934 [0.890, 0.971] | -0.04 | 0.935 | +0.14 | 0.934 | 0.516 [0.346, 0.755] | 0.871 / 0.366 | 0.561 | 0.883 | 3363 / 3860 / 9183 |
| logB_d0.5_s2.5_k0.5_g2_m0 | 0.934 [0.892, 0.971] | -0.03 | 0.934 | +0.17 | 0.934 | 0.515 [0.347, 0.753] | 0.866 / 0.367 | 0.561 | 0.880 | 3369 / 3891 / 9183 |
| logB_d0.5_s2.5_k0.5_g3_m0 | 0.934 [0.891, 0.971] | -0.03 | 0.934 | +0.16 | 0.934 | 0.515 [0.346, 0.754] | 0.867 / 0.367 | 0.561 | 0.880 | 3367 / 3884 / 9183 |
| `log` defaults (rank 33 of 55) | 0.861 [0.785, 0.926] | +0.60 | 0.732 | +1.73 | 0.858 | 0.493 [0.388, 0.625] | 0.644 / 0.400 | 0.526 | 0.767 | 3672 / 5701 / 9183 |
| SpotMAX deployed | 0.783 [0.740, 0.831] | +0.87 | 0.636 | +2.48 | 0.783 | 0.409 [0.295, 0.558] | 0.495 / 0.348 | 0.436 | 0.621 | 3193 / 6445 / 9183 |

Names: `k` is the separation in diameters, `g` the pre-gate in robust standard deviations, `m` the z merge
(1 on, 0 off).

### `log` on val

| set | germline CCC [95% CI] | bias | traced CCC | traced bias | germline CCC, Class A | spot F1 [95% CI] | P / R | F1 without B | detector-only F1 | matched / ours / expert |
|---|---|---|---|---|---|---|---|---|---|---|
| **logB_d0.5_s2.5_k0.5_g2_m0** (selected) | 0.874 [0.792, 0.963] | -0.37 | 0.866 | -0.74 | 0.880 | 0.482 [0.249, 0.825] | 0.950 / 0.323 | 0.660 | 0.835 | 1282 / 1349 / 3975 |
| logB_d0.5_s2.5_k0.5_g3_m0 | 0.874 [0.792, 0.962] | -0.37 | 0.866 | -0.74 | 0.880 | 0.481 [0.249, 0.825] | 0.951 / 0.322 | 0.659 | 0.835 | 1280 / 1346 / 3975 |
| logB_d0.5_s2.5_k0.5_g2_m1 | 0.872 [0.787, 0.963] | -0.37 | 0.866 | -0.74 | 0.878 | 0.482 [0.250, 0.826] | 0.956 / 0.323 | 0.661 | 0.835 | 1282 / 1341 / 3975 |
| logA_d0.5_s2.5 | 0.872 [0.787, 0.962] | -0.37 | 0.866 | -0.74 | 0.878 | 0.482 [0.250, 0.826] | 0.957 / 0.322 | 0.661 | 0.835 | 1280 / 1338 / 3975 |
| logB_d0.5_s2.5_k0.5_g4_m1 | 0.871 [0.786, 0.961] | -0.37 | 0.864 | -0.75 | 0.877 | 0.481 [0.249, 0.825] | 0.957 / 0.321 | 0.660 | 0.832 | 1277 / 1334 / 3975 |
| `log` defaults (D1, trial 6) | 0.926 [0.880, 0.957] | +0.16 | 0.872 | +0.72 | 0.922 | 0.558 [0.427, 0.734] | 0.766 / 0.439 | 0.711 | 0.824 | 1747 / 2282 / 3975 |
| SpotMAX deployed | 0.695 [0.643, 0.783] | +1.25 | 0.590 | +2.88 | 0.676 | 0.428 [0.317, 0.636] | 0.414 / 0.443 | 0.494 | 0.597 | 1762 / 4253 / 3975 |

The selected `log` configuration is `diameter_um` 0.5, `min_snr` 2.5, `pre_gate_sigmas` 2.0 and
`merge_z_columns` off: tied with the pre-gate 3 variant at 0.874 (0.8744 against 0.8740) and ahead on
F1 (0.4816 against 0.4811).

* The five finalists are within 0.0033 of each other in val CCC and 0.0012 in val F1: the val step could
  not tell them apart, and the choice among them means no more than the rule that made it.
* Against the deployed SpotMAX recipe on val the selected set has the higher CCC (0.874 against 0.695),
  the smaller bias (-0.37 against +1.25), the higher F1 (0.482 against 0.428, from precision 0.950
  against 0.414 at recall 0.323 against 0.443) and the higher detector-only F1 (0.835 against 0.597).
* The untuned `log` defaults beat every finalist on val (CCC 0.926, F1 0.558). The difference is the
  heat-shock gonad (table under "What the heat-shock gonads show"): at SNR 2.5 it is undercounted
  (-0.48 per germline nucleus), at the default SNR 2.0 it is not (-0.14). Train chose SNR 2.5 because
  three of its four gonads are without heat shock. Under the rule the defaults (rank 33 on train) were
  not a finalist; this is reported, not acted on.
* The pre-gate of 2 lets more than 30000 candidates through on noHS_herm_000 (35,024) and HS_herm_001
  (49,611), so those runs log the flood warning. With the cap lifted (no metric computed) noHS_herm_000
  gave the same 1130 spots and HS_herm_001 425 instead of 424 (two spots replaced by three), so the
  counts are not floors in practice. The pre-gate 3 sets stay under the cap on all six gonads (at most
  20,827 candidates).

### SpotMAX on train

| stage | what varied | sets | wall time |
|---|---|---|---|
| S1 | `thresholding_method` triangle, otsu, li by `effect_size_min` 2, 3, 4, 5, at the deployed `spot_radius_um` 0.3 and `gauss_sigma_um` 0.08 | 12 | 31.6 min, 2 workers |
| S2 | `threshold_triangle` (the S1 best, at effect size 4): (`spot_radius_um`, `gauss_sigma_um`) (0.25, 0.08), (0.35, 0.08), (0.3, 0.05), (0.3, 0.12) by `effect_size_min` 3.5, 4.0, 4.5 | 12 | 27.9 min, 1 worker |

With the shared detections the 24 sets took 7 SpotMAX runs per gonad (3 in S1, 4 in S2). The
`threshold_li` run on noHS_herm_000 alone took 29 minutes (24,599 spots before the effect-size filter;
the flood guard fired once in S1, most likely in this run).

S1, pooled germline CCC / bias on train:

| `thresholding_method` \ `effect_size_min` | 2 | 3 | 4 | 5 |
|---|---|---|---|---|
| threshold_triangle | 0.694 / +1.31 | 0.783 / +0.87 (deployed) | **0.839 / +0.39** | 0.810 / -0.01 |
| threshold_otsu | 0.691 / +0.61 | 0.787 / +0.26 | 0.809 / -0.06 | 0.738 / -0.35 |
| threshold_li | 0.120 / +7.49 | 0.338 / +3.53 | 0.635 / +1.46 | 0.751 / +0.46 |

Otsu starves the heat-shock gonad (bias -0.88 to -1.41 on HS_herm_000 at every effect size) and li floods
the gonads without heat shock at low effect sizes (bias +9.25 on them at 2; 25,196 spots on the four
train gonads against the expert's 9183), the two failure modes the deployed config's comments describe.
The five best of the 24, and the deployed recipe (rank 13):

| set | germline CCC [95% CI] | bias | traced CCC | traced bias | germline CCC, Class A | spot F1 [95% CI] | P / R | F1 without B | detector-only F1 | matched / ours / expert |
|---|---|---|---|---|---|---|---|---|---|---|
| smS2_tri_r0.3_g0.12_e3.5 | 0.855 [0.830, 0.881] | +0.26 | 0.764 | +0.89 | 0.856 | 0.394 [0.266, 0.563] | 0.581 / 0.298 | 0.426 | 0.654 | 2741 / 4715 / 9183 |
| smS2_tri_r0.35_g0.08_e4 | 0.841 [0.803, 0.881] | +0.19 | 0.767 | +0.58 | 0.839 | 0.377 [0.258, 0.533] | 0.571 / 0.281 | 0.405 | 0.630 | 2583 / 4522 / 9183 |
| smS1_tri_r0.3_g0.08_e4 | 0.839 [0.813, 0.874] | +0.39 | 0.745 | +1.13 | 0.838 | 0.396 [0.274, 0.560] | 0.555 / 0.308 | 0.424 | 0.638 | 2825 / 5090 / 9183 |
| smS2_tri_r0.3_g0.12_e4 | 0.835 [0.758, 0.891] | +0.03 | 0.734 | +0.29 | 0.836 | 0.376 [0.251, 0.547] | 0.613 / 0.272 | 0.407 | 0.646 | 2495 / 4071 / 9183 |
| smS2_tri_r0.25_g0.08_e4 | 0.833 [0.797, 0.880] | +0.27 | 0.725 | +0.87 | 0.831 | 0.377 [0.257, 0.540] | 0.552 / 0.286 | 0.404 | 0.622 | 2628 / 4764 / 9183 |
| deployed (smS1_tri_r0.3_g0.08_e3) | 0.783 [0.740, 0.831] | +0.87 | 0.636 | +2.48 | 0.783 | 0.409 [0.295, 0.558] | 0.495 / 0.348 | 0.436 | 0.621 | 3193 / 6445 / 9183 |

The best set has the largest smoothing searched (0.12 um) and, at that smoothing, the smallest effect
size (3.5), so the SpotMAX optimum may lie outside the grid; the 24-set budget did not allow following
it. In S2 the lightest smoothing (0.05 um) was the worst choice at every effect size.

### SpotMAX on val

| set | germline CCC [95% CI] | bias | traced CCC | traced bias | germline CCC, Class A | spot F1 [95% CI] | P / R | F1 without B | detector-only F1 | matched / ours / expert |
|---|---|---|---|---|---|---|---|---|---|---|
| **smS2_tri_r0.3_g0.12_e3.5** (selected) | 0.784 [0.721, 0.891] | +0.60 | 0.717 | +1.20 | 0.760 | 0.444 [0.347, 0.630] | 0.506 / 0.395 | 0.528 | 0.635 | 1571 / 3104 / 3975 |
| smS2_tri_r0.35_g0.08_e4 | 0.775 [0.732, 0.846] | +0.55 | 0.756 | +0.87 | 0.752 | 0.420 [0.338, 0.581] | 0.488 / 0.369 | 0.500 | 0.605 | 1466 / 3006 / 3975 |
| smS1_tri_r0.3_g0.08_e4 | 0.724 [0.656, 0.850] | +0.86 | 0.668 | +1.62 | 0.699 | 0.424 [0.325, 0.620] | 0.449 / 0.402 | 0.497 | 0.603 | 1598 / 3558 / 3975 |
| deployed (smS1_tri_r0.3_g0.08_e3) | 0.695 [0.643, 0.783] | +1.25 | 0.590 | +2.88 | 0.676 | 0.428 [0.317, 0.636] | 0.414 / 0.443 | 0.494 | 0.597 | 1762 / 4253 / 3975 |

The selected SpotMAX configuration is `threshold_triangle`, `spot_radius_um` 0.3, `gauss_sigma_um`
0.12, `effect_size_min` 3.5: the best of the three finalists on val, and better than the deployed recipe
on every pooled val number except recall (CCC 0.784 against 0.695, bias +0.60 against +1.25, F1 0.444
against 0.428). It still overcounts on both val gonads (+0.30 and +0.75 per germline nucleus). By section 7
of the protocol the SpotMAX reference is the better of the deployed and the re-tuned recipe chosen on val,
which on these numbers is the re-tuned one (changelog entry 2 leaves open whether this search or the
section 8 grid supplies the re-tuned competitor). On traced-subset CCC, the closest proxy the harness
computes for the primary metric, the val rule would pick smS2_tri_r0.35_g0.08_e4 instead (0.756
against 0.717), so the reference must be chosen under the primary metric before test. On val the selected `log` set leads the selected
SpotMAX set in CCC (0.874 against 0.784) and F1 (0.482 against 0.444); these are development numbers on
two gonads, not a win condition.

### What the heat-shock gonads show

There is no usable heat-shock male in train or val, so the only heat-shock gonads seen are two
hermaphrodites, HS_herm_000 (train) and HS_herm_001 (val). Per gonad, our spots in germline nuclei
against the expert's spots that fall in our germline nuclei, the means per germline nucleus, and the
agreement:

| detector | gonad | ours | ours / expert in germline nuclei | mean ours / expert | CCC | bias | traced CCC | traced bias | spot F1 | detector-only F1 |
|---|---|---|---|---|---|---|---|---|---|---|
| SpotMAX deployed | HS_herm_001 (val) | 2600 | 2529 / 969 | 2.17 / 0.83 | 0.643 | +1.34 | 0.436 | +3.76 | 0.317 | 0.432 |
| SpotMAX re-tuned | HS_herm_001 (val) | 1909 | 1847 / 969 | 1.58 / 0.83 | 0.721 | +0.75 | 0.536 | +2.23 | 0.347 | 0.539 |
| `log` defaults | HS_herm_001 (val) | 829 | 811 / 969 | 0.69 / 0.83 | 0.957 | -0.14 | 0.932 | -0.17 | 0.427 | 0.830 |
| `log` tuned | HS_herm_001 (val) | 424 | 414 / 969 | 0.35 / 0.83 | 0.792 | -0.48 | 0.571 | -1.18 | 0.249 | 0.504 |
| SpotMAX deployed | HS_herm_000 (train) | 1552 | 1527 / 1334 | 2.52 / 2.20 | 0.864 | +0.32 | 0.843 | +0.85 | 0.260 | 0.636 |
| SpotMAX re-tuned | HS_herm_000 (train) | 1138 | 1118 / 1334 | 1.84 / 2.20 | 0.816 | -0.36 | 0.772 | -0.39 | 0.231 | 0.618 |
| `log` defaults | HS_herm_000 (train) | 1351 | 1319 / 1334 | 2.17 / 2.20 | 0.938 | -0.02 | 0.978 | +0.29 | 0.357 | 0.916 |
| `log` tuned | HS_herm_000 (train) | 986 | 975 / 1334 | 1.61 / 2.20 | 0.875 | -0.59 | 0.946 | -0.66 | 0.299 | 0.844 |

* The deployed recipe does not undercount these two heat-shock gonads; it overcounts them, by +0.32 and
  +1.34 foci per germline nucleus (+0.85 and +3.76 in the traced nuclei). The heat-shock undercount
  that the nine-gonad baseline attributes to the deployed recipe is not visible in train or val; it
  comes from gonads this package did not look at, and seeing it during development needs the heat-shock
  male spot sets re-exported from Imaris without region growing.
* Most expert spots under heat shock lie outside our germline nuclei: 1334 of 5152 and 969 of 2753 fall
  inside them. Whole-field spot recall and F1 are therefore low for every detector on these gonads, and
  the per-nucleus numbers are the informative ones.
* Tuning towards the gonads without heat shock moves both detectors towards undercounting heat shock:
  the tuned `log` set undercounts both heat-shock gonads (-0.59 and -0.48); the re-tuned SpotMAX
  undercounts the training one (-0.36) and overcounts the val one less than the deployed recipe (+0.75
  against +1.34); the untuned `log` defaults are nearly unbiased on both (-0.02 and -0.14). The
  annotator's thresholds were lower under heat shock, and one global threshold chosen on a training set
  that is three quarters no-heat-shock tends to undercount heat shock against that truth. A detector
  meant to remove the heat-shock undercount needs heat-shock males with usable truth in train or val;
  today train has none usable and val holds no male at all.

### The two profiles

* `config/profiles/bench_log_tuned.yaml`: `spots.backend: log`, `spots.merge_z_columns: false`,
  `spots.log.diameter_um: 0.5`, `spots.log.min_snr: 2.5`, `spots.log.pre_gate_sigmas: 2.0`.
* `config/profiles/bench_spotmax_tuned.yaml`: `spots.thresholding_method: threshold_triangle`,
  `spots.spot_radius_um: 0.3`, `spots.gauss_sigma_um: 0.12`, `spots.effect_size_min: 3.5`.

Both extend `config/config.yaml` and change only `spots:` keys (`tests/test_bench_tuning.py` checks it,
and that the overlay used as bench parameters over the base section gives exactly the profile's section;
the golden runs' `spots:` section of all six cached train and val gonads equals that of
`config/config.yaml`). Neither is the default, neither has been run on a test gonad, and
`config/config.yaml`, its hash and the deployed recipe are unchanged. Select one with
`germquant run X.nd2 --profile bench_log_tuned` (or `bench_spotmax_tuned`); a `log` run carries the QC
flag `spots:unvalidated_backend=log`.

### What this does not say

* These are development numbers on four train and two val gonads. The intervals are cluster bootstraps by
  gonad (1000 resamples, seed 20251105); with two val gonads they rest on one gonad per condition.
* Nothing here is a win condition of the protocol, and the test split was not touched. The searches
  differ from the pre-registered design (the SpotMAX grid, val scores for finalists only, the germline
  partition as the objective); protocol changelog entry 2 records the differences, and each of them is
  post hoc (protocol rule 3).
* The objective is the germline partition of the harness, not the traced-subset partition of protocol
  section 5.1, and the harness's traced subset is not that partition either.
* No heat-shock male informed any choice, and val holds hermaphrodites only.
* Imaris Spots is close in spirit to a Laplacian of Gaussian, so part of the `log` backend's agreement
  with this truth may come by construction; agreement with this truth is not agreement with a person
  (protocol section 9).

## Leave-one-gonad-out

Package E1, 2026-10-10: the protocol's primary partition (section 5.1) in the harness, leave-one-gonad-out
estimates for the five detectors (section 6.2), and net3d retrained on train plus val (section 6.1). No test
gonad was run, scored or looked at: every run below is on the six train and val gonads with usable spot
truth, the harness's test guard was never passed (no `--final` run; `final_runs.jsonl` does not exist), and
the plan for the one-time test scoring is `docs/FINAL_SCORING.md`, which nothing here executed. Every E1 output
folder is under `C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/logo_e1/` (the SpotMAX grid at the
end of this section writes under `prereg_grid/`), every model under `C:/Users/ryane/germquant_golden/models/`,
and every scoring run is a record in `trials.jsonl`.

### The protocol's primary partition (section 5.1)

The harness's "traced" subset (our germline nuclei with more than half of their volume inside one truth mask,
counted in OUR labels) is not the partition section 5.1 makes primary, and reviewers said so. Section 5.1
takes the truth masks themselves as the nuclei, the same for every method and independent of any method's
segmentation. `germquant.benchmark.spots` now computes exactly that, next to the existing partitions:

* Nuclei: the recovered Imaris masks (`<stem>__imaris_labels.tif`) with fragments below 1000 voxels removed
  (section 4.2; `truth_mask_lut`). The truth export had already split collided values and dropped the
  unclaimed ones (its params JSON: `n_shared_values` 0 and 2 to 8 orphan values per gonad on the six gonads),
  so the fragment rule is the only step left; it removes one label among the train and val gonads, a
  44-voxel label of n2rad51_noHS_male_000.
* Assignment (`assign_to_masks`), for truth spots and the method's spots alike: the mask label at the centre
  voxel; if that voxel is background, the most frequent non-zero label of the 3 x 3 x 3 neighbourhood (ties
  to the smaller label); otherwise unassigned. Centre voxels: the truth files' `zi, yi, xi`; for the method,
  round(coordinate / voxel), the convention of section 4.3. The spot region R of section 5 (the masks
  dilated by one voxel in every axis) is read as the 3 x 3 x 3 cube, and then a centre voxel is in R exactly
  when the rule assigns it (a test checks this voxel by voxel), so the R restriction needs no extra step.
* Counts: truth count = truth spots assigned to the nucleus, from Class A+B (primary), Class A and Class B
  spots; method count = the method's end-to-end spots assigned to it (every spot of its final list, whatever
  nucleus of its own it was given). With `--detector-only`, `bench_detector.py` runs the spot stage a second
  time with the truth masks in place of the deployed nuclei (`run_backend(..., labels=...)`) and reports that
  detector-only count next to the end-to-end one, as section 5.1 asks.
* Statistics: Pearson r, Lin's CCC, bias (method minus truth) and mean absolute error (`stats.mean_abs_error`,
  now also in the other agreement blocks), per gonad (`tm_` columns of `per_gonad.csv`: `tm_n_nuclei`,
  `tm_n_fragments_removed`, assigned and unassigned truth and method spots, `tm_ccc`, `tm_pearson_r`,
  `tm_bias`, `tm_mae`, `tm_mean_expert`, `tm_mean_ours`, the same under `tm_class_a_` and `tm_class_b_`, and
  `tm_det_` for the detector-only count) and pooled with the cluster bootstrap by gonad
  (`per_nucleus_truth_masks`, `_class_a`, `_class_b`, `_detector_only` in pooled.json and the trial records).
  `per_truth_nucleus.csv` holds one row per truth nucleus (`truth_label`, `n_voxels`, `expert_n`,
  `expert_n_class_a`, `expert_n_class_b`, `ours_n`, `ours_n_detector_only`).
* Where: `bench_detector.py`, `bench_tune.py` (summary column `tm_ccc`) and `benchmark_spots.py` all report
  it; `stats.paired_cluster_bootstrap` is the paired comparison of section 6.3 (both methods' pooled
  statistic on the same resample of gonads).
* Spot F1 inside R (section 5.2): the method's spots and the truth spots whose centre voxel lies in R (the
  rule above) are matched one to one in the anisotropic cylinder of the whole-field score. Per gonad the
  `tm_spots_` columns hold the end-to-end score, `tm_spots_no_b_` Class A alone (Class B spots in R and the
  method spots matched to them leave both counts), `tm_spots_b_recall` the Class B recall and
  `tm_det_spots_` the detector-only score; pooled.json has them as `spots_truth_masks`,
  `spots_truth_masks_class_a` and `spots_truth_masks_detector_only`. This is the W2 metric of
  docs/FINAL_SCORING.md; the whole-field F1 (`spots_all`) stays for comparison. It was added in commit
  e19a029, after the leave-one-gonad-out runs below, so their held-out tables have no `tm_spots_` columns
  and every spot F1 before "Results of the pre-registered SpotMAX grid" is the whole-field one; that subsection
  adds the F1 inside R from the same runs scored again.

On the six gonads it covers 676 truth nuclei (train: 88, 107, 143 and 97 for HS_herm_000, noHS_herm_000,
noHS_male_000 and noHS_male_001; val: 95 and 146 for HS_herm_001 and noHS_herm_001) and assigns 2054 of the
13158 truth spots (357 of 5152 in HS_herm_000, 162 of 2753 in HS_herm_001: most heat-shock truth spots lie
outside the traced nuclei). Truth check: the section 4.4 check, run with `scripts/check_truth_assignment.py` against the lab's mask-channel
assignment (`imaris_mask_label`, Imaris's per-spot maximum of the Masked channel, the value the xlsx coloc counts
come from; no xlsx export of these gonads is registered), agrees on 99.35 to 99.99 percent of the spots of all
eight train and val gonads. The disagreements are boundary spots (0 to 14 geometric-only and 1 to 22
coloc-only per gonad) and spots whose coloc value is a removed mask fragment (0 to 58); only HS_male_002, a
region-growing gonad, has spots assigned to two different nuclei (3). No truth-export bug was found. Table:
`benchmarks/n2rad51_8666f07/truth_assignment_check_trainval.csv`. The test gonads are checked in the pre-flight
of docs/FINAL_SCORING.md.

### The folds

Section 6.2 over the train and val gonads; for spot detectors the two region-growing heat-shock males drop
out, which leaves six gonads and six folds (`germquant.benchmark.logo.logo_folds`, `bench_logo.py folds`). A
fold tunes or trains on the other five; its own selection uses the designated val gonads that are among them
and, in place of a held-out val gonad, the first train gonad of the same condition by tag order.

| held out | role | fold train gonads | fold val gonads |
|---|---|---|---|
| n2rad51_HS_herm_000 | train | noHS_herm_000, noHS_male_000, noHS_male_001 | HS_herm_001, noHS_herm_001 |
| n2rad51_HS_herm_001 | val | noHS_herm_000, noHS_male_000, noHS_male_001 | noHS_herm_001, HS_herm_000 (substitute) |
| n2rad51_noHS_herm_000 | train | HS_herm_000, noHS_male_000, noHS_male_001 | HS_herm_001, noHS_herm_001 |
| n2rad51_noHS_herm_001 | val | HS_herm_000, noHS_male_000, noHS_male_001 | HS_herm_001, noHS_herm_000 (substitute) |
| n2rad51_noHS_male_000 | train | HS_herm_000, noHS_herm_000, noHS_male_001 | HS_herm_001, noHS_herm_001 |
| n2rad51_noHS_male_001 | train | HS_herm_000, noHS_herm_000, noHS_male_000 | HS_herm_001, noHS_herm_001 |

The two heat-shock folds tune on the three no-heat-shock train gonads only: with one heat-shock gonad left in
train after the region-growing exclusion, holding out either heat-shock hermaphrodite leaves none to tune on.

### net3d: the frozen recipe on every fold

E1 fixes the recipe rather than the selection: each fold retrains net3d on its five gonads with the
recipe of the selected model `net3d_20261010` (same architecture, data and loss settings, normalisation,
seed 0, AdamW with 200 warm-up steps and the cosine decay over the maximum of 40 epochs of 200 steps), stops
after the selected 4 epochs and keeps the selected threshold 0.40, with no validation and no selection inside
the fold. `germquant.spots.net3d.train.recipe_from` reads all of that from the model's `training_log.json`, and
`germquant train-spots --recipe-from MODELDIR --train-tags ...` (or `bench_logo.py net3d-fold`) runs it; the
five gonads are taken in the recipe's order (its train gonads, then its val gonads, the held-out one left out).

The check that the recipe is frozen: retrained by `--recipe-from` on its own four train gonads, without the
validation the original run made after every epoch, the model's weights after epoch 4 equal those of the
selected checkpoint `checkpoints/epoch_004.pt` bit for bit (all 134 tensors), every epoch's training loss and
learning rate equal the original run's (epoch losses 0.048756, 0.031216, 0.029856, 0.030222), and the stored
threshold is 0.40 (`models/net3d_recipe_check_20261010/`, 163 s on the RTX 5090, commit eca6f12; the file
sha256 differs from the original `model.pt` only through its metadata). So the folds below differ from the
selected model only in their training gonads.

What this estimate does and does not cover. The recipe's epoch and threshold were chosen on the two
designated val gonads (protocol changelog entry 3). In the four folds that hold out a train gonad, the held-out
gonad never informed the recipe; in the two folds that hold out a val gonad (HS_herm_001, noHS_herm_001), it
did, so those two held-out scores are optimistic for net3d, and the estimate as a whole leaves out the
variance of net3d's own selection step. Section 6.2 asks for selection inside each fold; repeating net3d's
selection per fold (twelve epochs times eighteen thresholds on each fold's val gonads) was not done.

### The tuned classical detectors: the D3 search replayed on every fold

"Rerun the selection on the five non-held-out gonads within the same search space and budget": the D3
searches were staged and data dependent (tuning_d3/search_plan.json): for `log` a 30-set grid (stage A),
17 sets around its best (stage B), up to 8 around the best of A and B (stage C), the 5 best of all on train
scored on val and the best of those on val selected; for SpotMAX a 12-set grid (S1), 12 sets around its best
(S2), the 3 best on train to val. The objective is D3's: pooled germline per-nucleus CCC rounded to three
decimals, then pooled spot F1. `germquant.benchmark.logo.replay_log_search` and `replay_spotmax_search`
replay the whole procedure on a fold's own gonads: the grid ranked on the fold's train gonads, the
neighbourhoods generated around the fold's own bests by the plan's rules, the finalists ranked on the fold's
val gonads. Scores come from the run folders (`index_runs`), so a set already scored on a gonad is never run
again, and the replay names any (set, gonad) pair it still needs. Replayed on the D3 split itself (train the
four train gonads, val the two val gonads) it gives D3's stage bests, finalists and selections exactly, and the
generated stage lists equal the D3 configuration files name for name and parameter for parameter.

What had to be run:

* Every D3 set on the two val gonads (D3 had scored only its finalists there): 55 `log` sets (trial lines
  97 to 151, search `d3_log`, stage `logo_val`, `logo_e1/val_scoring/log/`) and 24 SpotMAX sets (lines 174 to
  197, search `d3_spotmax`, `logo_e1/val_scoring/spotmax/`). The searches' distinct-set counts stay 55 and 24.
  The val scores of the D3 finalists came out identical to D3's (5 `log` and 4 SpotMAX sets: pooled points and
  intervals, every per-gonad number).
* In the four folds whose train gonads include the heat-shock hermaphrodite HS_herm_000, the replay follows the
  D3 path (stage-A best 0.5 um at SNR 2.5, stage-B best pre-gate 2 with z merge on, so stage C is D3's; S1 best
  triangle at effect size 4, so S2 is D3's). In the two heat-shock folds, whose train gonads are the three
  no-heat-shock train gonads, it does not. For `log` the stage-B best is separation 0.75 diameters with pre-gate
  4 (`logB_d0.5_s2.5_k0.75_g4_m1`), so stage C is centred there: 8 new sets (`logo_e1/fold_search/log_C_train3noHS/`,
  search `logo_log_train3noHS`, trial lines 152 to 159), one of which reached the finalists and was scored on
  those folds' val gonads (lines 160 and 161, the second on the substitute HS_herm_000). For SpotMAX the S1 best is
  `threshold_otsu` at effect size 4, so S2 is centred there: 12 new sets (`logo_e1/fold_search/spotmax_S2_train3noHS/`,
  search `logo_spotmax_train3noHS`, lines 198 to 209), two of which reached the finalists and were scored on those
  folds' val gonads (lines 210 to 213). Each fold's own search stays within D3's budget (55 `log` sets, 24
  SpotMAX sets).

The per-fold selections (`logo_e1/selection_log.json`, `logo_e1/selection_spotmax.json`, with every
finalist's train and val objective):

| held out | `log` selected | SpotMAX selected |
|---|---|---|
| n2rad51_HS_herm_000 | logB_d0.5_s2.5_k0.75_g2_m1 | smS1_otsu_r0.3_g0.08_e4 |
| n2rad51_HS_herm_001 | logB_d0.5_s2.5_k0.75_g2_m1 | smS1_otsu_r0.3_g0.08_e4 |
| n2rad51_noHS_herm_000 | logC_d0.5_s2.25_k0.5_g2_m1 | smS2_tri_r0.3_g0.12_e3.5 (the D3 selection) |
| n2rad51_noHS_herm_001 | logC_d0.5_s2.25_k0.5_g2_m1 | smS2_tri_r0.3_g0.12_e3.5 (the D3 selection) |
| n2rad51_noHS_male_000 | logC_d0.46_s2.25_k0.5_g2_m1 | smS2_tri_r0.3_g0.12_e3.5 (the D3 selection) |
| n2rad51_noHS_male_001 | logB_d0.5_s2.5_k0.5_g2_m0 (the D3 selection) | smS2_tri_r0.3_g0.12_e3.5 (the D3 selection) |

`log` picks another set than D3's in five of six folds: the same diameter and SNR with a wider separation (0.75
diameters, z merge on) in the two heat-shock folds, SNR 2.25 instead of 2.5 with z merge on (and 0.46 um in one)
in three others. SpotMAX picks the D3 set wherever the fold's tuning gonads include a heat-shock gonad, and an
otsu set where they do not.

### Deployed SpotMAX and the log defaults

Neither is tuned, so a fold changes nothing and the held-out score of a gonad is its per-gonad score. Both
were re-run per gonad with the same command as the tuned detectors (`bench_logo.py heldout --params {}`), so
that every held-out score carries the section 5.1 partition and the detector-only count. The end-to-end
numbers equal the earlier runs of the same parameters (D1, trial lines 3 to 6: thirteen headline numbers per
gonad, none different), and the deployed recipe's detections equal the golden spot tables of all six gonads
(`detector.compare_with_golden`, every column, every row). Likewise every re-tuned SpotMAX and tuned `log`
held-out score equals the earlier run of the same parameters on the same gonad wherever one exists.

### Held-out scoring

Every held-out gonad of every detector is one `bench_detector.py --detector-only` run on that gonad alone
(`bench_logo.py heldout` and `net3d-fold` call it), with the fold in the trial record's `note`
(`{"logo": {"detector", "held_out", "fold_train", "fold_val", "substitute", ...}}`; net3d's also names the fold
model and its sha256, the tuned detectors' the selected set and the finalists). Outputs:
`logo_e1/heldout/<detector>/<held-out tag>/`; trial lines 162 to 167 (tuned `log`), 168 to 173 (`log` defaults),
214 to 219 (SpotMAX re-tuned), 220 to 225 (SpotMAX deployed) and 226 to 231 (net3d). The pooled estimates and
the paired comparisons are `bench_logo.py report` (`logo_e1/report/`: `logo_report.json`, `logo_per_gonad.csv`,
`logo_pooled.csv`, `logo_paired.csv`): the six held-out gonads pooled with the cluster bootstrap by gonad (1000
resamples, seed 20251105), and every pair of detectors compared on the same resamples (protocol section 6.3).

Provenance. E1 added trial lines 97 to 231 (135 records, all split train or val). Lines 97 to 99 carry git sha
eca6f12 with `-dirty` (the sha is read when a record is appended, and at that moment the tree held an
uncommitted edit of `scripts/bench_logo.py` and its test); every other line carries 240cc72, clean. The
scoring code that ran (`germquant.benchmark.spots`, `detector.score`, `tuning`) is that of 9b206ba in every
record, unchanged in both later commits. Every net3d model was trained at commit 240cc72.

### The commands

With `$L = "C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/logo_e1"`, `$T` the `tuning_d3` folder
next to it and `$M = "C:/Users/ryane/germquant_golden/models"`; `CUDA_VISIBLE_DEVICES=-1` except for net3d, whose
commands ran under the workstation's GPU lock (`gpu_lock.py run --name ... -- COMMAND`).

```powershell
# every D3 set on the val gonads ($L/configs/*_all.json are the D3 configuration files concatenated)
python scripts/bench_tune.py --backend log --configs $L/configs/log_all.json --split val --search d3_log `
    --stage logo_val --budget 60 --jobs 4 --chunks 2 --out-root $L/val_scoring/log
python scripts/bench_tune.py --backend spotmax --configs $L/configs/spotmax_all.json --split val `
    --search d3_spotmax --stage logo_val --budget 24 --jobs 1 --chunks 4 --out-root $L/val_scoring/spotmax
# the sets the replay called for in the heat-shock folds (written with logo.log_stage_c / spotmax_stage_s2)
python scripts/bench_tune.py --backend log --configs $L/configs/log_C_train3noHS.json --split train `
    --tags n2rad51_noHS_herm_000 n2rad51_noHS_male_000 n2rad51_noHS_male_001 --search logo_log_train3noHS `
    --stage C --budget 13 --jobs 3 --out-root $L/fold_search/log_C_train3noHS
python scripts/bench_tune.py --backend log --configs $L/configs/log_valstep_train3noHS.json --split val `
    --search logo_log_train3noHS --stage val --budget 13 --out-root $L/fold_search/log_val_train3noHS
python scripts/bench_tune.py --backend log --configs $L/configs/log_valstep_train3noHS.json --split train `
    --tags n2rad51_HS_herm_000 --search logo_log_train3noHS --stage val_substitute --budget 13 `
    --out-root $L/fold_search/log_valsub_train3noHS
python scripts/bench_tune.py --backend spotmax --configs $L/configs/spotmax_S2_train3noHS.json --split train `
    --tags n2rad51_noHS_herm_000 n2rad51_noHS_male_000 n2rad51_noHS_male_001 --search logo_spotmax_train3noHS `
    --stage S2 --budget 12 --jobs 1 --chunks 4 --out-root $L/fold_search/spotmax_S2_train3noHS
python scripts/bench_tune.py --backend spotmax --configs $L/configs/spotmax_valstep_train3noHS.json --split val `
    --search logo_spotmax_train3noHS --stage val --budget 12 --out-root $L/fold_search/spotmax_val_train3noHS
python scripts/bench_tune.py --backend spotmax --configs $L/configs/spotmax_valstep_train3noHS.json --split train `
    --tags n2rad51_HS_herm_000 --search logo_spotmax_train3noHS --stage val_substitute --budget 12 `
    --out-root $L/fold_search/spotmax_valsub_train3noHS
# the replayed selections (each was first run on the D3 and val scorings alone: it named the sets and gonads
# still missing and exited with code 3; the fold_search runs above are those sets)
python scripts/bench_logo.py select --backend log --runs $T/log_A_train $T/log_B_train $T/log_C_train $T/log_val `
    $L/val_scoring/log $L/fold_search/log_C_train3noHS $L/fold_search/log_val_train3noHS `
    $L/fold_search/log_valsub_train3noHS --out $L/selection_log.json
python scripts/bench_logo.py select --backend spotmax --runs $T/spotmax_S1_train $T/spotmax_S2_train $T/spotmax_val `
    $L/val_scoring/spotmax $L/fold_search/spotmax_S2_train3noHS $L/fold_search/spotmax_val_train3noHS `
    $L/fold_search/spotmax_valsub_train3noHS --out $L/selection_spotmax.json
# held-out scoring
python scripts/bench_logo.py heldout --detector log_tuned --backend log --selection $L/selection_log.json `
    --out-root $L/heldout/log_tuned
python scripts/bench_logo.py heldout --detector log_defaults --backend log --params "{}" --out-root $L/heldout/log_defaults
python scripts/bench_logo.py heldout --detector spotmax_retuned --backend spotmax --selection $L/selection_spotmax.json `
    --out-root $L/heldout/spotmax_retuned
python scripts/bench_logo.py heldout --detector spotmax_deployed --backend spotmax --params "{}" `
    --out-root $L/heldout/spotmax_deployed
# net3d: the recipe check, the six folds (each under the GPU lock) and the train-plus-val model
germquant train-spots --cache C:/Users/ryane/germquant_golden/bench_cache --split-file tests/golden/rad51_split.json `
    --out $M/net3d_recipe_check_20261010 --recipe-from $M/net3d_20261010 `
    --train-tags n2rad51_HS_herm_000 n2rad51_noHS_herm_000 n2rad51_noHS_male_000 n2rad51_noHS_male_001 `
    --device cuda --prep-dir $M/net3d_prep --min-free-gb 8
python scripts/bench_logo.py net3d-folds --recipe-from $M/net3d_20261010 --models-root $M/net3d_logo_20261010 `
    --out-root $L/heldout/net3d --log-dir $L/logs --gpu-lock C:/Users/ryane/germquant_golden/gpu_lock.py `
    --prep-dir $M/net3d_prep --min-free-gb 8
germquant train-spots --cache C:/Users/ryane/germquant_golden/bench_cache --split-file tests/golden/rad51_split.json `
    --out $M/net3d_trainval_20261010 --recipe-from $M/net3d_20261010 `
    --train-tags n2rad51_HS_herm_000 n2rad51_noHS_herm_000 n2rad51_noHS_male_000 n2rad51_noHS_male_001 `
                 n2rad51_HS_herm_001 n2rad51_noHS_herm_001 --device cuda --prep-dir $M/net3d_prep --min-free-gb 8
# the pooled estimates and the paired comparisons
python scripts/bench_logo.py report --heldout net3d=$L/heldout/net3d log_tuned=$L/heldout/log_tuned `
    log_defaults=$L/heldout/log_defaults spotmax_retuned=$L/heldout/spotmax_retuned `
    spotmax_deployed=$L/heldout/spotmax_deployed --out $L/report
```

### The net3d fold models

`C:/Users/ryane/germquant_golden/models/net3d_logo_20261010/fold_<held-out tag>/` (model.pt, model_card.json with
the recipe, training_log.json, checkpoints), all trained at commit 240cc72 on the RTX 5090 under the GPU lock,
158 to 191 s each (the first fold also built the loss-region maps of the two val gonads, 35 s):

| fold (held out) | trained on | model.pt sha256 |
|---|---|---|
| n2rad51_HS_herm_000 | noHS_herm_000, noHS_male_000, noHS_male_001, HS_herm_001, noHS_herm_001 | 969fec3edb8d980b8e69e2c98f7837de617fcbfee4f9876e5c0d098c073fa548 |
| n2rad51_HS_herm_001 | HS_herm_000, noHS_herm_000, noHS_male_000, noHS_male_001, noHS_herm_001 | 4fa6d1edb9cf6bad95469b47f93f9e7f1a4113ec505aae0e7cb1a70c246a3844 |
| n2rad51_noHS_herm_000 | HS_herm_000, noHS_male_000, noHS_male_001, HS_herm_001, noHS_herm_001 | 932a5dac3bb879301910bd09d381d2dbf7740f018abeff2785bb7b7df0fe05a3 |
| n2rad51_noHS_herm_001 | HS_herm_000, noHS_herm_000, noHS_male_000, noHS_male_001, HS_herm_001 | 1aa2176f3b547de91dfb0aa38fc8fbf46af103e29dc7d0727189d7dd9f78a37c |
| n2rad51_noHS_male_000 | HS_herm_000, noHS_herm_000, noHS_male_001, HS_herm_001, noHS_herm_001 | aff880e6e90853059c2f7bb4f8e0f063b23271f43c155b63cd08b470a8cf8ecb |
| n2rad51_noHS_male_001 | HS_herm_000, noHS_herm_000, noHS_male_000, HS_herm_001, noHS_herm_001 | 963da323640fa764945af2d50b89a3a53074368f16dbbc0507285971bb0e71b8 |

### Results

Development numbers on the six train and val gonads; each held-out gonad is scored by a configuration or model
that its own fold chose or trained without it (with net3d's caveat above). Cluster bootstrap by gonad, 1000
resamples, seed 20251105; truth from Class A+B unless stated; end-to-end counts on the deployed nucleus masks.
`logo_e1/report/` holds every number (`logo_pooled.csv`, `logo_per_gonad.csv`, `logo_paired.csv`,
`logo_report.json`).

| pooled over the 6 held-out gonads | SpotMAX deployed | SpotMAX re-tuned (D3 replayed) | log defaults | log tuned (D3 replayed) | net3d (frozen recipe) |
|---|---|---|---|---|---|
| per nucleus, truth masks (5.1), A+B: CCC | 0.601 [0.499, 0.718] | 0.742 [0.655, 0.812] | 0.777 [0.623, 0.920] | 0.893 [0.818, 0.955] | 0.878 [0.834, 0.935] |
| per nucleus, truth masks, A+B: bias | +2.66 [+1.95, +3.28] | +0.37 [-0.72, +1.17] | +1.33 [+0.50, +2.20] | +0.24 [-0.50, +0.95] | +0.23 [-0.58, +0.88] |
| per nucleus, truth masks, A+B: MAE | 2.79 [2.28, 3.31] | 1.48 [1.20, 1.85] | 1.43 [0.76, 2.22] | 0.84 [0.50, 1.26] | 0.85 [0.49, 1.27] |
| per nucleus, truth masks, Class A: CCC | 0.599 [0.489, 0.723] | 0.752 [0.674, 0.817] | 0.776 [0.625, 0.911] | 0.895 [0.826, 0.954] | 0.881 [0.839, 0.929] |
| per nucleus, truth masks, detector-only count: CCC | 0.727 [0.674, 0.783] | 0.746 [0.588, 0.842] | 0.824 [0.703, 0.931] | 0.893 [0.831, 0.938] | 0.879 [0.834, 0.935] |
| per nucleus, germline, A+B: CCC | 0.750 [0.688, 0.814] | 0.821 [0.734, 0.881] | 0.882 [0.820, 0.934] | 0.898 [0.850, 0.947] | 0.892 [0.826, 0.954] |
| per nucleus, germline, A+B: bias | +1.02 [+0.69, +1.35] | -0.02 [-0.51, +0.41] | +0.43 [+0.08, +0.85] | -0.05 [-0.40, +0.31] | +0.16 [-0.35, +0.56] |
| per nucleus, traced (harness), A+B: CCC | 0.620 [0.524, 0.728] | 0.751 [0.652, 0.830] | 0.776 [0.624, 0.923] | 0.889 [0.815, 0.950] | 0.867 [0.820, 0.926] |
| spots A+B: F1 | 0.415 [0.320, 0.562] | 0.368 [0.219, 0.567] | 0.513 [0.414, 0.647] | 0.501 [0.359, 0.727] | 0.512 [0.364, 0.705] |
| spots without Class B: F1 | 0.454 [0.363, 0.572] | 0.425 [0.269, 0.583] | 0.573 [0.474, 0.685] | 0.576 [0.432, 0.747] | 0.580 [0.423, 0.767] |
| matched / ours / expert spots | 4955 / 10698 / 13158 | 3505 / 5883 / 13158 | 5419 / 7983 / 13158 | 4737 / 5764 / 13158 | 5131 / 6881 / 13158 |

Nuclei in every column: 676 truth masks, 4567 germline nuclei, 633 nuclei in the harness's traced subset.

Per held-out gonad, truth-mask CCC / germline CCC / spot F1:

| held-out gonad | SpotMAX deployed | SpotMAX re-tuned | log defaults | log tuned | net3d |
|---|---|---|---|---|---|
| HS_herm_000 | 0.843 / 0.864 / 0.260 | 0.587 / 0.622 / 0.133 | 0.979 / 0.938 / 0.357 | 0.951 / 0.863 / 0.299 | 0.825 / 0.789 / 0.241 |
| HS_herm_001 | 0.391 / 0.643 / 0.317 | 0.734 / 0.796 / 0.213 | 0.911 / 0.957 / 0.427 | 0.550 / 0.769 / 0.248 | 0.959 / 0.984 / 0.490 |
| noHS_herm_000 | 0.611 / 0.743 / 0.512 | 0.743 / 0.844 / 0.532 | 0.693 / 0.823 / 0.564 | 0.901 / 0.932 / 0.629 | 0.848 / 0.801 / 0.553 |
| noHS_herm_001 | 0.679 / 0.783 / 0.636 | 0.835 / 0.891 / 0.630 | 0.818 / 0.880 / 0.734 | 0.963 / 0.973 / 0.828 | 0.952 / 0.969 / 0.834 |
| noHS_male_000 | 0.545 / 0.767 / 0.592 | 0.746 / 0.893 / 0.604 | 0.480 / 0.755 / 0.626 | 0.772 / 0.879 / 0.729 | 0.816 / 0.892 / 0.732 |
| noHS_male_001 | 0.533 / 0.728 / 0.469 | 0.742 / 0.871 / 0.462 | 0.826 / 0.909 / 0.624 | 0.964 / 0.975 / 0.686 | 0.921 / 0.938 / 0.651 |

Per held-out gonad, truth-mask bias / traced CCC / F1 without Class B / our spots:

| held-out gonad | SpotMAX deployed | SpotMAX re-tuned | log defaults | log tuned | net3d |
|---|---|---|---|---|---|
| HS_herm_000 | +0.91 / 0.843 / 0.292 / 1552 | -2.33 / 0.574 / 0.158 / 574 | +0.34 / 0.978 / 0.402 / 1351 | -0.68 / 0.942 / 0.345 / 950 | -1.70 / 0.812 / 0.287 / 721 |
| HS_herm_001 | +3.79 / 0.436 / 0.377 / 2600 | -0.77 / 0.792 / 0.353 / 537 | -0.15 / 0.932 / 0.656 / 829 | -1.08 / 0.569 / 0.423 / 407 | -0.15 / 0.959 / 0.733 / 958 |
| noHS_herm_000 | +2.35 / 0.611 / 0.515 / 1934 | +1.44 / 0.738 / 0.535 / 1556 | +1.47 / 0.678 / 0.566 / 1679 | +0.62 / 0.890 / 0.632 / 1306 | +0.91 / 0.803 / 0.556 / 1861 |
| noHS_herm_001 | +2.34 / 0.703 / 0.652 / 1653 | +0.67 / 0.855 / 0.650 / 1195 | +1.24 / 0.821 / 0.754 / 1453 | +0.21 / 0.964 / 0.856 / 1086 | -0.17 / 0.950 / 0.863 / 1000 |
| noHS_male_000 | +3.13 / 0.564 / 0.603 / 1440 | +1.03 / 0.775 / 0.617 / 995 | +2.98 / 0.486 / 0.635 / 1554 | +1.59 / 0.784 / 0.743 / 1219 | +1.28 / 0.818 / 0.747 / 1219 |
| noHS_male_001 | +3.30 / 0.530 / 0.473 / 1519 | +1.34 / 0.755 / 0.466 / 1026 | +1.24 / 0.813 / 0.630 / 1117 | -0.01 / 0.959 / 0.694 / 796 | +0.65 / 0.910 / 0.657 / 1122 |

Paired differences (the first detector minus the second, both pooled on the same 1000 resamples of the six
gonads; an asterisk marks an interval that excludes zero):

| pair | truth masks CCC | germline CCC | traced CCC | spot F1 | F1 without Class B | truth masks MAE |
|---|---|---|---|---|---|---|
| net3d - log tuned | -0.015 [-0.084, +0.073] | -0.006 [-0.088, +0.082] | -0.021 [-0.097, +0.078] | +0.011 [-0.049, +0.112] | +0.004 [-0.052, +0.104] | +0.012 [-0.370, +0.399] |
| net3d - log defaults | +0.101 [-0.051, +0.239] | +0.010 [-0.069, +0.087] | +0.092 [-0.064, +0.234] | -0.001 [-0.065, +0.072] | +0.007 [-0.056, +0.087] | -0.583 [-1.211, +0.162] |
| net3d - SpotMAX re-tuned | +0.136 [+0.090, +0.204] * | +0.071 [+0.004, +0.137] * | +0.117 [+0.065, +0.193] * | +0.144 [+0.081, +0.212] * | +0.155 [+0.085, +0.246] * | -0.632 [-0.792, -0.446] * |
| net3d - SpotMAX deployed | +0.277 [+0.139, +0.409] * | +0.142 [+0.018, +0.250] * | +0.247 [+0.114, +0.371] * | +0.097 [+0.027, +0.185] * | +0.126 [+0.043, +0.259] * | -1.939 [-2.631, -1.160] * |
| log tuned - log defaults | +0.116 [-0.014, +0.217] | +0.016 [-0.072, +0.101] | +0.113 [-0.022, +0.222] | -0.012 [-0.080, +0.084] | +0.003 [-0.075, +0.087] | -0.595 [-1.068, +0.043] |
| log tuned - SpotMAX re-tuned | +0.151 [+0.040, +0.278] * | +0.077 [+0.012, +0.151] * | +0.138 [+0.014, +0.274] * | +0.133 [+0.087, +0.175] * | +0.151 [+0.113, +0.188] * | -0.643 [-1.164, -0.199] * |
| log tuned - SpotMAX deployed | +0.292 [+0.215, +0.360] * | +0.148 [+0.078, +0.200] * | +0.268 [+0.201, +0.337] * | +0.085 [+0.021, +0.169] * | +0.122 [+0.070, +0.192] * | -1.951 [-2.429, -1.484] * |
| log defaults - SpotMAX re-tuned | +0.035 [-0.148, +0.254] | +0.061 [-0.051, +0.189] | +0.025 [-0.167, +0.261] | +0.145 [+0.062, +0.200] * | +0.148 [+0.061, +0.224] * | -0.049 [-0.951, +0.757] |
| log defaults - SpotMAX deployed | +0.176 [+0.042, +0.337] * | +0.132 [+0.047, +0.223] * | +0.155 [+0.024, +0.312] * | +0.097 [+0.068, +0.126] * | +0.119 [+0.070, +0.185] * | -1.357 [-2.219, -0.680] * |
| SpotMAX re-tuned - SpotMAX deployed | +0.141 [-0.043, +0.292] | +0.071 [-0.074, +0.155] | +0.130 [-0.063, +0.277] | -0.047 [-0.102, +0.008] | -0.029 [-0.092, +0.026] | -1.308 [-2.024, -0.511] * |

### Which differences survive

Plainly, with six gonads and the paired cluster bootstrap:

* net3d and the tuned `log` each beat both SpotMAX configurations on every number in the table: the truth-mask,
  germline and traced CCC, both spot F1 variants and the truth-mask mean absolute error; all twenty-four of those
  paired intervals exclude zero. Against the SpotMAX re-tuning procedure the truth-mask CCC gains are +0.136
  (net3d) and +0.151 (tuned `log`); against the deployed recipe +0.277 and +0.292.
* net3d against the tuned `log`: no difference survives on any number (truth-mask CCC -0.015 [-0.084, +0.073]).
* net3d and the tuned `log` against the `log` defaults: no difference survives (truth-mask CCC +0.101 and +0.116,
  both intervals reaching below zero). The defaults are far off on one male gonad (noHS_male_000, truth-mask
  CCC 0.480) and the best or close to it on the two heat-shock hermaphrodites (0.979 and 0.911), so six gonads
  do not settle it.
* The `log` defaults beat the deployed recipe on every number, and the re-tuned SpotMAX on the two spot F1
  variants only.
* SpotMAX re-tuned against deployed: only the smaller mean absolute error survives (1.48 against 2.79 foci per
  truth nucleus); the CCC gain (+0.141 [-0.043, +0.292]) does not, and its F1 is lower (-0.047 [-0.102, +0.008]).

Sensitivity to net3d's two optimistic folds (`logo_e1/report_train_heldout/`, `bench_logo.py report --tags` with
the four train gonads: the folds where the recipe never saw the held-out gonad): net3d's truth-mask CCC is 0.845
[0.816, 0.890] and the tuned `log`'s 0.889 [0.808, 0.957]; net3d minus tuned `log` is -0.044 [-0.115, +0.029] in
truth-mask CCC, and on those four gonads the tuned `log` is ahead in germline CCC (-0.056 [-0.108, -0.005]) and spot
F1 (-0.037 [-0.065, -0.018]). net3d still beats SpotMAX re-tuned (+0.154 [+0.084, +0.231]) and deployed (+0.219
[+0.076, +0.348]) in truth-mask CCC there.

### The train-plus-val model (protocol section 6.1)

`C:/Users/ryane/germquant_golden/models/net3d_trainval_20261010/model.pt`, sha256
`d091b9ef58850783387a4f5aa0fea5ecc9310abc027145771936fecd0de001c0`: the frozen recipe of net3d_20261010 (sha256
8cdb710b..., commit 87d77f1; epoch 4 of the 40-epoch schedule, threshold 0.40, seed 0) retrained on the six train
and val gonads with usable spot truth (HS_herm_000, noHS_herm_000, noHS_male_000, noHS_male_001, HS_herm_001,
noHS_herm_001, in that order), no validation; commit 240cc72, 159 s on the RTX 5090 under the GPU lock (epoch
losses 0.047437, 0.025822, 0.024363, 0.024501). Its `model_card.json` records the sha256, the training gonads, the
epochs (4 run of the 40-epoch schedule), the threshold, the git sha, the recipe (source model, its sha256 and git
sha, its train and val gonads, the selected epoch and threshold, the selection rule) and limitations written from
the gonads it saw. It is the net3d of `docs/FINAL_SCORING.md`. It was not scored on anything: every gonad with
usable spot truth outside the test set is one of its training gonads.

### What this does and does not say

* These are development estimates on six gonads, two of them heat-shocked hermaphrodites and none a heat-shock
  male. The intervals are wide, and one gonad moves a pooled number by several hundredths.
* The tuned classical detectors' estimates include their whole selection (the replayed search), so they are the
  out-of-sample performance of the D3 procedure, judged by its own germline objective. net3d's estimate does not
  include its selection, and two of its six folds saw their held-out gonad through the recipe; on the other four
  it trails the tuned `log` by a few hundredths. A selection by the primary metric was not replayed for any of
  these five detectors (the SpotMAX grid procedure below is the first that was).
* The `log` backend and Imaris Spots are both blob detectors, so part of `log`'s agreement may come by
  construction; agreement with this truth is not agreement with a person (protocol section 9).
* The two heat-shock folds tune on no-heat-shock gonads only, which is where the replayed searches leave the D3
  path and where the tuned detectors undercount heat shock (truth-mask bias -0.68 and -1.08 for the tuned `log`,
  -2.33 and -0.77 for SpotMAX re-tuned): one global threshold chosen without heat-shock data undercounts heat shock
  against this truth.
* No test gonad was scored, and nothing here is a win condition of the protocol. The section 4.4 cross-check of
  the geometric assignment against the lab's mask-channel assignment agrees on 99.35 to 99.99 percent of the
  spots of the eight train and val gonads (see "Leave-one-gonad-out", truth check).

### The pre-registered SpotMAX grid (the code and the planned run)

Decision 2 of `docs/FINAL_SCORING.md` offers the grid of protocol section 8 item 2 in place of the D3 search as
the re-tuned SpotMAX competitor. The code for it is in place. This subsection was written before the run; the
results are in the next one.

`germquant.benchmark.logo.prereg_spotmax_grid` lists the 60 sets in the protocol's order: `thresholding_method`
(otsu, li, triangle, yen), then `effect_size_min` (2.0 to 4.0 in steps of 0.5), then `z_merge_gap_um` (0.6, 0.8,
1.0). A set changes only those three keys. Every other `spots:` key keeps the golden run's value, so
`pg_tri_e3.0_z0.8` is the deployed recipe.

The rule is that of sections 6.1 and 6.4. Every set is scored on the selection gonads, and the best on the
primary metric is selected. There is no train stage, so every set is a finalist. There is no early stop either:
all 60 sets are ranked. `primary_objective_on` is the primary metric, pooled as `bench_detector.py` pools it. It
is the truth-mask CCC of section 5.1 (Class A+B truth, end-to-end counts; the `per_nucleus_truth_masks` block)
over the pooled nuclei of the selection gonads, rounded to three decimals. Ties go to the spot F1 inside R of
section 5.2 (`spots_truth_masks`), then to the grid's order. A NaN value ranks last.

`select_grid` ranks the sets and lists the (set, gonad) pairs still to run. A pair is still to run when it has
no score, or when its run folder lacks the truth-mask counts (per_truth_nucleus.csv or the `tm_spots_`
columns). So a set is never passed over only because a folder lacks them. `replay_grid` applies the rule to one
fold's selection gonads (its val gonads from `logo_folds`). The run index reads the truth-mask counts only when
it is built with `index_runs(..., truth_masks=True)`. The D3 replay therefore reads and returns exactly what it
did before.

`bench_logo.py grid-configs` writes the 60 sets as a `bench_tune.py --configs` file. `bench_logo.py select
--backend spotmax --search prereg-grid` applies the rule on every fold and on the designated split, from the run
folders given with `--runs`. The designated split is the split's own train and val gonads with usable spot
truth (`designated_split`). The selection file has the layout of the D3 selection files, plus `search`, the grid
(`configs`), the designated split's selection (`designated`) and the missing pairs merged by set (`to_run`).
`heldout --selection` reads it as it reads a D3 selection. While any pair is missing, `select` exits with
code 3. Run folders that disagree, or a fold without selection gonads, are refused with code 2.

On each gonad SpotMAX runs once per method and z gap: 12 runs, since the five effect sizes share one
detection. The train scores do not enter the designated split's choice. The folds use those of two train
gonads: in the two folds that hold out a val gonad, n2rad51_HS_herm_000 or n2rad51_noHS_herm_000 takes that
gonad's place. `select` run on the val folder alone lists exactly these two under `to_run`.
n2rad51_noHS_male_000 and n2rad51_noHS_male_001 are never a selection gonad, so the planned train run scores
them for the record only. Adding `--tags n2rad51_HS_herm_000 n2rad51_noHS_herm_000` to that run would save 24
of the plan's 72 SpotMAX runs and change no selection. A SpotMAX worker peaks at 10.4 to 12.1 GB of working set
(see `scripts/bench_detector.py` above), and the workstation has 63 GB. The commands therefore use two workers,
as the D3 SpotMAX stages did (E1 used one). Raise `--jobs` only when the free memory allows about 12 GB per
worker; it changes no score.

**What is known before the run.** Nine grid sets are the same SpotMAX configuration as D3 S1 sets: methods otsu,
li and triangle at effect sizes 2.0, 3.0 and 4.0 with z gap 0.8 (for example `pg_otsu_e4.0_z0.8` and
`smS1_otsu_r0.3_g0.08_e4`). Their effective `spots:` sections are identical (`detector.effective_section`), but
they are recorded under another `params_sha256` because the keys differ. E1 scored every D3 set on the two val
gonads under the primary partition. The detection and scoring code have not changed since then, apart from the
new inside-R columns. So the nine sets check the run: their val truth-mask CCC should equal E1's
(`logo_e1/val_scoring/spotmax/summary.csv`). One of them, otsu at effect size 4.0, reached 0.822 there, and the
deployed recipe, also one of the nine, reached 0.556 (decision 2 of `docs/FINAL_SCORING.md`). If the run
reproduces E1, the designated grid choice therefore has a val truth-mask CCC of at least 0.822.

**How the code reads section 6.1.** Section 6.1 says that parameter-only competitors are tuned on train plus
val by the grid rule of 6.4. Here every set gets a val score (6.4), and the set with the best val score is
frozen. Train scores are used only where a train gonad takes a held-out val gonad's place in a fold.

The planned run, with `$G = "C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/prereg_grid"` and `$L`
the `logo_e1` folder:

```powershell
python scripts/bench_logo.py grid-configs --out $G/configs.json
python scripts/bench_tune.py --backend spotmax --configs $G/configs.json --split val --search prereg_spotmax_grid `
    --stage val --budget 60 --jobs 2 --chunks 4 --out-root $G/val
python scripts/bench_tune.py --backend spotmax --configs $G/configs.json --split train --search prereg_spotmax_grid `
    --stage train --budget 60 --jobs 2 --chunks 4 --out-root $G/train
python scripts/bench_logo.py select --backend spotmax --search prereg-grid --runs $G/val $G/train `
    --out $G/selection.json
python scripts/bench_logo.py heldout --detector spotmax_grid --backend spotmax --selection $G/selection.json `
    --out-root $G/logo/heldout/spotmax_grid
python scripts/bench_logo.py report --heldout net3d=$L/heldout/net3d log_tuned=$L/heldout/log_tuned `
    log_defaults=$L/heldout/log_defaults spotmax_retuned=$L/heldout/spotmax_retuned `
    spotmax_deployed=$L/heldout/spotmax_deployed spotmax_grid=$G/logo/heldout/spotmax_grid --out $G/logo/report
```

### Results of the pre-registered SpotMAX grid (2026-10-11)

Protocol changelog entry 6 declared this run before it started, and entry 7 records it. Every run is on the six
train and val gonads with usable spot truth. No test gonad was run, scored or opened: no `--final` run was made,
`final_runs.jsonl` does not exist, and no record in `trials.jsonl` names a test gonad. The numbers come from the
run folders under `$G` and `$L` (as above) and under `tuning_d3/` next to them, read with the repository's own
functions: `logo.index_runs` with `truth_masks=True`, `select_grid`, `replay_grid`, `pool_heldout`, `paired` and
`detector.compare_with_golden`. The analysis that produced them ran no detector. The detector runs reported here
are the grid run and the re-runs of "The F1 inside R, and the E1 runs scored again". CCC is the pooled truth-mask
CCC of section 5.1 (Class A+B truth, end-to-end counts). F1 inside R is the spot F1 of section 5.2, and the
whole-field F1 is the older `spots_all` score. Bias and MAE are in foci per truth nucleus.

**What ran.** The planned run above, as written (two workers, four chunks), from the worktree at commit abc49ca
(entry 6, committed at 00:59:15 MDT on 2026-10-11) on a clean tree. A driver script ran the commands in that order
with `CUDA_VISIBLE_DEVICES=-1` and logged every step to `$G/run.log`. The run went from 00:59:51 to 02:18:22 MDT,
and every step exited with code 0: the 60 sets (`grid-configs`), `bench_tune.py` on val (24 SpotMAX runs, 1293 s)
and on train (48 runs, 2740 s; the two region-growing males skipped), `select` (4 s), the six held-out runs (60 to
98 s each) and `report` (20 s). Trial lines 232 to 351 hold the grid: 120 records (60 val, 60 train), 60 distinct
sets under the search label `prereg_spotmax_grid`, the effect-size shortcut on and complete in every record,
SpotMAX 1.3.1. Lines 352 to 357 hold the held-out runs. All carry abc49ca without `-dirty`. No run failed, so the
failure rule of entry 6 (c) was not needed. With the 24 D3 and 12 E1 SpotMAX sets, 96 named SpotMAX sets have now
been scored, 87 of them distinct in effect and 77 of those on val, as entry 6 (d) counted. Counted per competitor,
77 sets scored on val is over the 60 of section 6.4; counting the grid alone is a reading of entry 6 (d), for the
project lead to confirm. The grid sets equal in effect to earlier sets are exactly the nine of entry 6 (b).
`bench_tune.py` ranks the sets by its own built-in rule (the germline CCC, then the whole-field F1: the D3
objective) in its console summary and in the `rank` and `objective_ccc` columns of its `summary.csv`, and names
that objective in its trial records. By that rule `pg_otsu_e3.0_z0.8` is first on val (germline CCC 0.863; the
selection is 10th) and `pg_tri_e4.0_z0.8` is first on train. That order selects nothing. The selection is
`select_grid`'s, in `$G/selection.json`.

**The checks of entry 6 (c).** All passed. They were made from the files after the driver had run every step,
`select` and the held-out runs included, so they confirmed the selection afterwards instead of gating it.

* The deployed recipe. The golden `spots:` section of all six gonads holds the deployed values that entry 6 (b)
  lists, and `pg_tri_e3.0_z0.8` leaves it unchanged. Its detections equal the golden spot tables on all six
  gonads (`detector.compare_with_golden` against the parquet and the CSV tables: every compared column, row order
  included, and no per-nucleus count different). That is 1552, 2600, 1934, 1653, 1440 and 1519 spots on
  HS_herm_000, HS_herm_001, noHS_herm_000, noHS_herm_001, noHS_male_000 and noHS_male_001.
* The nine sets of entry 6 (b). Their val tables equal those of E1 (trial lines 174 to 184), and their train
  tables those of D3 stage S1 (lines 54 to 64), in every common column of per_gonad, per_nucleus and the spot
  tables; on val per_truth_nucleus is equal too. Only `detect_s` and `score_s` differ. Val CCC, the same in E1
  and in the grid, at effect sizes 2, 3 and 4: triangle 0.515, 0.556, 0.635; otsu 0.766, 0.794, 0.822; li 0.327,
  0.410, 0.544. Train germline CCC, the same in D3 and in the grid: triangle 0.694, 0.783, 0.839; otsu 0.691,
  0.787, 0.809; li 0.120, 0.338, 0.635.
* The held-out runs. Each of the six, made without the shortcut, equals the grid run of the same set on the same
  gonad in per_gonad, per_nucleus, per_truth_nucleus, the spot tables and the parameters. Only the run times and
  the detector-only columns differ, and the grid makes no detector-only run. The HS_herm_001 run
  (`pg_tri_e4.0_z1.0`) is the first check of a 1.0 um gap without the shortcut. li, yen and the 0.6 um gap still
  rest on the order of steps alone.
* Also equal to the grid's tables: the E1 held-out runs made without the shortcut, that is the deployed recipe on
  all six gonads (trial lines 220 to 225) and the otsu set on both heat-shock hermaphrodites (lines 214 and 215).
* Every set has a score on both val gonads, and no value is undefined.

**The selection.** `pg_otsu_e4.0_z0.8` (`threshold_otsu`, `effect_size_min` 4.0, `z_merge_gap_um` 0.8), over
241 truth nuclei: val CCC 0.822 (0.822104 unrounded), F1 inside R 0.612 (402 matched, 642 of ours and 671 truth
spots in R), bias -0.12, MAE 1.17, whole-field F1 0.365. `pg_otsu_e4.0_z1.0` ties it on both keys. The two give
byte-identical detections on all six gonads, so listing order decides between two runs that are the same in
effect. Recomputed from the run folders, the selection and all 60 val scores equal what `select` wrote. The
unrounded CCC (summary column `tm_ccc`) picks the same set. It reorders ranks 3 to 5 (`pg_otsu_e4.0_z0.6`,
0.82040, moves above `pg_otsu_e3.5_z0.8`, 0.81992). It also reverses three triangle pairs that tie at three
decimals and that the rule orders by F1 inside R: at effect sizes 4.0, 3.5 and 2.5 the 1.0 um gap moves above the
0.8 um gap (ranks 16 and 17, 19 and 20, 27 and 28).

The ranking on val follows. All 15 otsu sets take ranks 1 to 15. In `$G/selection.json` the `val_objective` of the
`designated` entry lists all 60 sets in rank order with both keys, and `$G/analysis/tables.md` lists them with
their train scores. The `bench_tune.py` summaries hold every set's val and train CCC, bias, MAE and whole-field F1.
The F1 inside R, the tie-break, is not in them or in the grid's trial records, only in each set folder's
`pooled.json` and `per_gonad.csv` and in `selection.json`. Train columns pool the four train gonads (435 truth
nuclei):

| rank | set | val CCC | unrounded | val F1 inside R | val bias | train CCC | train bias |
|---|---|---|---|---|---|---|---|
| 1 | pg_otsu_e4.0_z0.8 (selected) | 0.822 | 0.82210 | 0.612 | -0.12 | 0.711 | +0.40 |
| 2 | pg_otsu_e4.0_z1.0 | 0.822 | 0.82210 | 0.612 | -0.12 | 0.711 | +0.40 |
| 3 | pg_otsu_e3.5_z0.8 | 0.820 | 0.81992 | 0.620 | +0.06 | 0.699 | +0.80 |
| 4 | pg_otsu_e3.5_z1.0 | 0.820 | 0.81992 | 0.620 | +0.06 | 0.699 | +0.80 |
| 5 | pg_otsu_e4.0_z0.6 | 0.820 | 0.82040 | 0.610 | -0.10 | 0.703 | +0.44 |
| 6 | pg_otsu_e3.5_z0.6 | 0.815 | 0.81525 | 0.618 | +0.09 | 0.690 | +0.85 |
| 7 | pg_otsu_e3.0_z0.8 | 0.794 | 0.79371 | 0.624 | +0.23 | 0.637 | +1.39 |
| 8 | pg_otsu_e3.0_z1.0 | 0.794 | 0.79371 | 0.624 | +0.23 | 0.637 | +1.39 |
| 9 | pg_otsu_e3.0_z0.6 | 0.789 | 0.78852 | 0.621 | +0.26 | 0.629 | +1.44 |
| 10 | pg_otsu_e2.5_z0.8 | 0.777 | 0.77733 | 0.630 | +0.32 | 0.548 | +2.01 |
| 16 | pg_tri_e4.0_z0.8 (best triangle) | 0.635 | 0.63459 | 0.606 | +1.76 | 0.736 | +1.20 |
| 23 | pg_tri_e3.0_z0.8 (the deployed recipe) | 0.556 | 0.55647 | 0.599 | +2.91 | 0.626 | +2.53 |
| 24 | pg_li_e4.0_z0.8 (best li) | 0.544 | 0.54447 | 0.574 | +2.52 | 0.477 | +2.91 |
| 46 to 60 | the 15 yen sets | 0.071 | 0.07066 | 0.049 | -2.68 | 0.219 to 0.313 | -2.39 to -2.25 |

**Early stop and the six-gonad pick.** The early stop of section 10, read as in entry 6 (d), would have stopped
at the 20th set in listing order, `pg_li_e2.5_z0.8`. The best rounded val CCC so far was 0.822, against 0.815
before the last ten sets: a gain of 0.007. That is the fifth li set, before any triangle set (from position 31)
or yen set (from 46). The sets scored to that point give the same selection, `pg_otsu_e4.0_z0.8`. Integer
thousandths and floating point give the same stop. On the six train and val gonads pooled (676 truth nuclei;
reported only, entry 6 (c)), `select_grid` picks the same set again: CCC 0.747, F1 inside R 0.610, bias +0.21,
MAE 1.50, again tied with `pg_otsu_e4.0_z1.0`. Third is `pg_otsu_e4.0_z0.6` (0.741), and the deployed recipe is
18th (0.601).

**The selection per gonad, train and val.** Entry 6 (b) asks for the selection's train scores next to it. Per
gonad, CCC / bias / F1 inside R, with the deployed recipe for comparison:

| gonad | split | truth nuclei | pg_otsu_e4.0_z0.8 | deployed recipe (pg_tri_e3.0_z0.8) |
|---|---|---|---|---|
| HS_herm_000 | train | 88 | 0.587 / -2.33 / 0.401 | 0.843 / +0.91 / 0.635 |
| noHS_herm_000 | train | 107 | 0.759 / +1.03 / 0.676 | 0.611 / +2.35 / 0.659 |
| noHS_male_000 | train | 143 | 0.725 / +1.38 / 0.648 | 0.545 / +3.13 / 0.634 |
| noHS_male_001 | train | 97 | 0.816 / +0.72 / 0.637 | 0.533 / +3.30 / 0.549 |
| pooled train | | 435 | 0.711 / +0.40 / 0.609 | 0.626 / +2.53 / 0.621 |
| HS_herm_001 | val | 95 | 0.734 / -0.77 / 0.454 | 0.391 / +3.79 / 0.430 |
| noHS_herm_001 | val | 146 | 0.809 / +0.30 / 0.650 | 0.679 / +2.34 / 0.684 |
| pooled val | | 241 | 0.822 / -0.12 / 0.612 | 0.556 / +2.91 / 0.599 |

The heat-shock bias is what this table is for. The selection undercounts both heat-shock hermaphrodites. On
HS_herm_000 it counts 1.73 foci per truth nucleus where the truth has 4.06 (the deployed recipe 4.97); on
HS_herm_001 it counts 0.94 against 1.71. The deployed recipe overcounts both. On the four gonads without heat shock
(three train, one val) the selection overcounts by 0.30 to 1.38 foci per nucleus, the deployed recipe by 2.34 to
3.30. Pooled over the four train gonads the selection reaches CCC 0.711 [0.621, 0.791], bias +0.40 [-1.03, +1.26]
and MAE 1.69, against 0.626 [0.545, 0.748], +2.53 [+1.69, +3.20] and 2.71 for the deployed recipe. On those four
gonads pooled the selection ranks 4th of 60 by the rule; `pg_tri_e4.0_z0.8` and `pg_tri_e4.0_z1.0` lead there with
0.736 (reported only; train scores enter no designated selection).

**The section 7 reference.** Section 7 takes the better on val of SpotMAX deployed and SpotMAX re-tuned. If the
project lead confirms `pg_otsu_e4.0_z0.8` as SpotMAX re-tuned (decision 2 of `docs/FINAL_SCORING.md`), it is the
reference: val CCC 0.822 against 0.556 for the deployed recipe (+0.266). That difference is in sample, since the
set was chosen on these two gonads; out of sample the grid procedure minus the deployed recipe is +0.083 [-0.079,
+0.197] (below). Both are reported on test either way.

**The D3 candidates (post hoc, entry 2).** Entry 6 (a) asks for them either way. On the same val score (E1's
scoring of all 24 D3 sets, `logo_e1/val_scoring/spotmax/`): the D3 selection `smS2_tri_r0.3_g0.12_e3.5` reaches
0.695 (bias +1.32, MAE 1.63) and entry 2's alternative `smS2_tri_r0.35_g0.08_e4` 0.723 (+1.00, 1.68). The best
D3 set, `smS1_otsu_r0.3_g0.08_e4`, is `pg_otsu_e4.0_z0.8` in effect (0.822). The two triangle candidates lie
outside the grid (smoothing 0.12 um, radius 0.35 um). By val CCC they would sit between grid ranks 15 (0.760)
and 16 (0.635). Their F1 inside R is not known: their E1 val tables predate the `tm_spots_` columns.

**The candidate cap.** Three SpotMAX runs reached the cap (`max_spot_candidates` 30000), all in the train run, at
run.log lines 380, 426 and 461: 32720 candidates each, the brightest 30000 kept and 2720 dropped. The warning
names no set or gonad. From the log order, each falls inside one of the three `threshold_li` runs on
noHS_herm_000 (z gaps 0.6, 0.8 and 1.0 um), one warning per run. Those runs gave 25045, 24599 and 24520 spots
before the effect-size filter, where every other li run gave at most 13619, and D3's li run on that gonad gave
24599 with one such warning. The three runs served the 15 li sets on noHS_herm_000, so those 15 train counts are
floors. They were scored as they are (entry 6 (c)). They enter no designated selection, which uses val only. In
the one fold that selects on noHS_herm_000 (noHS_herm_001 held out) the li sets rank 31st or lower (best 0.347,
against 0.772 for the selected set). No val run and no held-out run reached the cap.

**Leave-one-gonad-out (entry 6 (e)).** `replay_grid` applied the rule to each fold's selection gonads, and every
fold selection equals `$G/selection.json`. Each held-out gonad was then scored afresh, without the shortcut:

| held out | selection gonads | selected | selection CCC / F1 inside R | held out: CCC / bias / F1 inside R / whole-field F1 |
|---|---|---|---|---|
| HS_herm_000 | HS_herm_001, noHS_herm_001 | pg_otsu_e4.0_z0.8 | 0.822 / 0.612 | 0.587 / -2.33 / 0.401 / 0.133 |
| HS_herm_001 | noHS_herm_001, HS_herm_000 (substitute) | pg_tri_e4.0_z1.0 | 0.766 / 0.648 | 0.472 / +2.84 / 0.488 / 0.326 |
| noHS_herm_000 | HS_herm_001, noHS_herm_001 | pg_otsu_e4.0_z0.8 | 0.822 / 0.612 | 0.759 / +1.03 / 0.676 / 0.495 |
| noHS_herm_001 | HS_herm_001, noHS_herm_000 (substitute) | pg_otsu_e4.0_z0.8 | 0.772 / 0.614 | 0.809 / +0.30 / 0.650 / 0.584 |
| noHS_male_000 | HS_herm_001, noHS_herm_001 | pg_otsu_e4.0_z0.8 | 0.822 / 0.612 | 0.725 / +1.38 / 0.648 / 0.594 |
| noHS_male_001 | HS_herm_001, noHS_herm_001 | pg_otsu_e4.0_z0.8 | 0.822 / 0.612 | 0.816 / +0.72 / 0.637 / 0.453 |

Five folds select `pg_otsu_e4.0_z0.8`. The fold that holds out HS_herm_001 selects on noHS_herm_001 and
HS_herm_000, and there a triangle set wins narrowly (0.766, against 0.758 for the best otsu set). It fits
HS_herm_000 (bias +0.17) and overcounts the held-out HS_herm_001 (+2.84).

Pooled over the six held-out gonads (cluster bootstrap by gonad, 1000 resamples, seed 20251105), with the five
detectors of entry 4 (b) taken from their runs scored again (see "The F1 inside R" below; `$G/logo/report_f1R/`):

| pooled over the 6 held-out gonads | truth masks CCC | bias | MAE | germline CCC | whole-field F1 | F1 inside R |
|---|---|---|---|---|---|---|
| SpotMAX grid (replayed) | 0.684 [0.551, 0.788] | +0.72 [-0.46, +1.70] | 1.78 [1.38, 2.30] | 0.772 [0.690, 0.876] | 0.366 [0.234, 0.543] | 0.603 [0.512, 0.655] |
| SpotMAX deployed | 0.601 [0.499, 0.718] | +2.66 [+1.95, +3.28] | 2.79 [2.28, 3.31] | 0.750 [0.688, 0.814] | 0.415 [0.320, 0.562] | 0.614 [0.543, 0.661] |
| SpotMAX re-tuned (D3 replayed) | 0.742 [0.655, 0.812] | +0.37 [-0.72, +1.17] | 1.48 [1.20, 1.85] | 0.821 [0.734, 0.881] | 0.368 [0.219, 0.567] | 0.627 [0.521, 0.678] |
| log defaults | 0.777 [0.623, 0.920] | +1.33 [+0.50, +2.20] | 1.43 [0.76, 2.22] | 0.882 [0.820, 0.934] | 0.513 [0.414, 0.647] | 0.786 [0.726, 0.859] |
| log tuned (D3 replayed) | 0.893 [0.818, 0.955] | +0.24 [-0.50, +0.95] | 0.84 [0.50, 1.26] | 0.898 [0.850, 0.947] | 0.501 [0.359, 0.727] | 0.848 [0.772, 0.900] |
| net3d (frozen recipe) | 0.878 [0.834, 0.935] | +0.23 [-0.58, +0.88] | 0.85 [0.49, 1.27] | 0.892 [0.826, 0.954] | 0.512 [0.364, 0.705] | 0.858 [0.793, 0.909] |

Inside R the grid procedure's precision is 0.545 and its recall 0.674 (the deployed recipe 0.470 and 0.883,
net3d 0.828 and 0.890).

Paired differences with the grid procedure (the other detector minus the grid procedure, as
`report_f1R/logo_paired.csv` lists them, both pooled on the same 1000 resamples; an asterisk marks an interval
that excludes zero):

| other detector minus SpotMAX grid | truth masks CCC | germline CCC | truth masks MAE | whole-field F1 | F1 inside R |
|---|---|---|---|---|---|
| net3d | +0.194 [+0.101, +0.334] * | +0.120 [+0.001, +0.227] * | -0.923 [-1.535, -0.515] * | +0.146 [+0.108, +0.204] * | +0.255 [+0.207, +0.317] * |
| log tuned | +0.209 [+0.100, +0.342] * | +0.126 [+0.053, +0.182] * | -0.935 [-1.392, -0.535] * | +0.134 [+0.050, +0.199] * | +0.245 [+0.184, +0.326] * |
| log defaults | +0.093 [-0.118, +0.342] | +0.110 [-0.044, +0.244] | -0.340 [-1.423, +0.499] | +0.146 [+0.089, +0.191] * | +0.183 [+0.086, +0.335] * |
| SpotMAX re-tuned (D3 replayed) | +0.058 [-0.028, +0.177] | +0.049 [-0.008, +0.109] | -0.291 [-0.898, +0.115] | +0.002 [-0.045, +0.032] | +0.024 [+0.002, +0.041] * |
| SpotMAX deployed | -0.083 [-0.197, +0.079] | -0.022 [-0.122, +0.092] | +1.016 [+0.344, +1.518] * | +0.049 [+0.005, +0.098] * | +0.011 [-0.046, +0.085] |

The F1 inside R between the five detectors of entry 4 (b), which their first held-out tables could not give:

| pair | F1 inside R | F1 inside R, Class A | F1 inside R, detector-only |
|---|---|---|---|
| net3d - log tuned | +0.010 [-0.058, +0.094] | +0.009 [-0.053, +0.083] | +0.026 [-0.038, +0.109] |
| net3d - log defaults | +0.072 [-0.044, +0.139] | +0.079 [-0.029, +0.141] | +0.072 [-0.033, +0.130] |
| net3d - SpotMAX re-tuned | +0.231 [+0.188, +0.299] * | +0.230 [+0.188, +0.295] * | +0.230 [+0.183, +0.310] * |
| net3d - SpotMAX deployed | +0.244 [+0.170, +0.331] * | +0.252 [+0.182, +0.344] * | +0.194 [+0.136, +0.251] * |
| log tuned - log defaults | +0.063 [-0.033, +0.114] | +0.070 [-0.016, +0.114] | +0.046 [-0.050, +0.094] |
| log tuned - SpotMAX re-tuned | +0.221 [+0.154, +0.314] * | +0.221 [+0.155, +0.311] * | +0.204 [+0.133, +0.311] * |
| log tuned - SpotMAX deployed | +0.235 [+0.193, +0.286] * | +0.243 [+0.197, +0.300] * | +0.168 [+0.116, +0.213] * |
| log defaults - SpotMAX re-tuned | +0.158 [+0.065, +0.321] * | +0.151 [+0.063, +0.305] * | +0.158 [+0.065, +0.342] * |
| log defaults - SpotMAX deployed | +0.172 [+0.097, +0.281] * | +0.173 [+0.096, +0.288] * | +0.122 [+0.061, +0.196] * |
| SpotMAX re-tuned - SpotMAX deployed | +0.014 [-0.073, +0.069] | +0.022 [-0.059, +0.086] | -0.036 [-0.148, +0.016] |

Per held-out gonad, F1 inside R:

| held-out gonad | SpotMAX deployed | SpotMAX re-tuned | SpotMAX grid | log defaults | log tuned | net3d |
|---|---|---|---|---|---|---|
| HS_herm_000 | 0.635 | 0.401 | 0.401 | 0.914 | 0.850 | 0.699 |
| HS_herm_001 | 0.430 | 0.454 | 0.488 | 0.832 | 0.516 | 0.923 |
| noHS_herm_000 | 0.659 | 0.696 | 0.676 | 0.767 | 0.882 | 0.846 |
| noHS_herm_001 | 0.684 | 0.686 | 0.650 | 0.827 | 0.925 | 0.939 |
| noHS_male_000 | 0.634 | 0.671 | 0.648 | 0.684 | 0.801 | 0.846 |
| noHS_male_001 | 0.549 | 0.617 | 0.637 | 0.784 | 0.908 | 0.879 |

**Which differences survive.** With six gonads and the paired cluster bootstrap:

* net3d and the tuned `log` beat the grid procedure on every number of the paired table with it (and on the traced
  CCC), as they beat both other SpotMAX configurations.
* The `log` defaults beat the grid procedure on every spot F1 (whole field and inside R, with or without Class B,
  end-to-end and detector-only) and on no CCC, bias or MAE.
* The D3 procedure (SpotMAX re-tuned) against the grid procedure: of the end-to-end numbers only its higher F1
  inside R survives, +0.024 [+0.002, +0.041] (Class A +0.031 [+0.006, +0.057]). Its CCC gain, +0.058 [-0.028,
  +0.177], does not. On the detector-only count it has the higher CCC (+0.041 [+0.011, +0.074]) and the smaller MAE
  (-0.210 [-0.348, -0.088]), and its detector-only spot F1 in the traced nuclei is higher too (+0.024 [+0.004,
  +0.039]).
* The deployed recipe against the grid procedure: no difference in any CCC (truth masks -0.083 [-0.197, +0.079]) or
  in the end-to-end F1 inside R (+0.011 [-0.046, +0.085]). The grid procedure has the smaller bias and the smaller
  mean absolute error (deployed minus grid: bias +1.945 [+1.443, +2.515], MAE +1.016 [+0.344, +1.518], on the
  germline nuclei +0.316 [+0.071, +0.569]). The deployed recipe has the higher whole-field F1 (+0.049 [+0.005,
  +0.098]) and the higher detector-only F1 inside R (+0.065 [+0.014, +0.159]).
* On the F1 inside R, net3d, the tuned `log` and the `log` defaults each beat both SpotMAX configurations of entry
  4 (b). No difference survives among the three, or between the two SpotMAX configurations (+0.014 [-0.073,
  +0.069]; on the whole-field F1 that pair was -0.047 [-0.102, +0.008]).

**The F1 inside R, and the E1 runs scored again.** Entry 6 (e) noted that pairs on the F1 inside R need that block
in the report and the held-out runs of entry 4 scored again. Commit ca8537d, whose parent is abc49ca, adds
`spots_truth_masks`, `spots_truth_masks_class_a` and `spots_truth_masks_detector_only` to `logo.DETECTION_BLOCKS`
(the `tm_spots_`, `tm_spots_no_b_` and `tm_det_spots_` counts), so that `bench_logo.py report` pools and pairs
them. A held-out folder without those columns (one written before commit e19a029) is left out of those blocks, and
its other numbers do not change. The commit changes only `logo.py`, `bench_logo.py` and `tests/test_bench_logo.py`.
It was prepared during the train run, after the val run had written its scores (01:23 MDT), and entered the tree
after the grid's last step: every grid record and the grid's report carry abc49ca without `-dirty`. It was
committed at 02:19:08 MDT.

At ca8537d the five detectors' held-out gonads were then run again with the same parameter sets and fold models
(commands below), from 02:19:10 to 02:30:19 MDT, every step with exit code 0 (`$G/rescore.log`; the last held-out
run ended at 02:29:59, then the report ran; trial lines 358 to 387, 30 records, six per detector, clean). All 30
runs reproduce their E1 originals (lines 162 to 173 and 214 to 231). In per_gonad.csv, per_nucleus.csv,
per_truth_nucleus.csv and per_spot_matches.csv every E1 column but the run times is identical, and so are the spot
tables (end-to-end and detector-only) and the parameter hashes. net3d, on the GPU, gave the same detections bit for
bit. Its fold model files, hashed after the run, still have the sha256 of E1's records (lines 226 to 231); the
re-run records do not carry that sha256. The run times differ (`detect_s` and `score_s` in all 30 runs,
`detect_detector_only_s` in 28). The new runs add 20 per_gonad columns (`tm_spots_*`, `tm_spots_no_b_*`,
`tm_spots_n_expert_class_b`, `tm_spots_b_recall`, `tm_det_spots_*`). The report over them, `$G/logo/report_f1R/`,
equals the E1 report (`logo_e1/report/`) in all 145 pooled and 290 paired rows and adds the F1 inside R. It also
equals all 174 pooled and 435 paired rows of the abc49ca report (`$G/logo/report/`). As in E1, the tuned `log`
reached its backend's candidate cap on HS_herm_001 (43967 candidates) and noHS_herm_000 (35024).

**What this does and does not say.**

* If the project lead confirms the recommendation of entry 6 and the readings of its (c) and (d), the grid's
  selection becomes SpotMAX re-tuned as section 8 item 2 defines it (the grid, chosen on val) and, by section 7,
  the reference. The readings that can be recomputed do not change the set: the six-gonad pick, the unrounded CCC
  and the 20 sets before the early stop all give `pg_otsu_e4.0_z0.8`. The rule picked the set that entry 6
  expected, and its scores per gonad show the heat-shock undercount that decision 2 named.
* No single grid set fits both heat-shock hermaphrodites against this truth. Every grid set's bias has the same
  sign on the two, and every set misses one of them by at least 1.7 foci per truth nucleus (the closest,
  `pg_otsu_e2.0_z0.6`: -1.70 on HS_herm_000, -0.64 on HS_herm_001). This was computed after the run; entry 6 did
  not declare it. In the folds, otsu won whenever HS_herm_001 was a selection gonad, and a triangle set won, by
  0.008, in the one fold where HS_herm_000 took its place.
* The leave-one-gonad-out estimate measures the grid procedure, not one set: five folds chose `pg_otsu_e4.0_z0.8`
  and one a triangle set. Out of sample the procedure's CCC (0.684) lies between the deployed recipe's (0.601) and
  the D3 procedure's (0.742), and neither difference survives. The D3 procedure's leads on the F1 inside R and on
  the detector-only CCC and MAE do survive (see "Which differences survive"). The largest per-gonad difference
  between the two procedures is the fold that holds out HS_herm_001, where the grid procedure picked a triangle set
  (held-out CCC 0.472 against 0.734). The selection's val CCC of 0.822 is in sample: the set was chosen on those
  two gonads.
* The grid procedure's estimate is not blind to these gonads in one respect: entry 6 wrote its readings (val only,
  rounding, tie-breaks) knowing the scores of nine of its sets on all six gonads, the selection's among them (entry
  6 (b) and (f)).
* net3d and the tuned `log` beat all three SpotMAX configurations on the F1 inside R as well, the W2 metric, by
  +0.221 to +0.255.
* These are development estimates on six gonads, two of them heat-shock hermaphrodites and none a heat-shock male.
  The intervals are wide, and nothing here is a win condition. A re-export of the two heat-shock males (decision 1)
  would add two folds to this estimate and leave the designated selection as it is (entry 6 (c)).
* No test gonad was run, scored or opened. Which configuration is the reference stays the project lead's choice
  (decision 2; protocol changelog entry 7).

**The commands.** The grid ran the planned commands above as written. The E1 runs were scored again at ca8537d
with `$R = "$G/logo/heldout_rescored"`. The five `heldout` runs ran side by side, with `CUDA_VISIBLE_DEVICES=-1`
except for net3d, which ran under the GPU lock. `$G/logo/net3d_rescore_selection.json` names each fold's E1
model (`models/net3d_logo_20261010/fold_<held-out tag>/model.pt`) in the layout of a selection file, with no
selection.

```powershell
python scripts/bench_logo.py heldout --detector log_defaults --backend log --params "{}" `
    --out-root $R/log_defaults
python scripts/bench_logo.py heldout --detector log_tuned --backend log --selection $L/selection_log.json `
    --out-root $R/log_tuned
python scripts/bench_logo.py heldout --detector spotmax_deployed --backend spotmax --params "{}" `
    --out-root $R/spotmax_deployed
python scripts/bench_logo.py heldout --detector spotmax_retuned --backend spotmax `
    --selection $L/selection_spotmax.json --out-root $R/spotmax_retuned
python C:/Users/ryane/germquant_golden/gpu_lock.py run --name "net3d held-out re-scoring (entry 6 e)" -- `
    python scripts/bench_logo.py heldout --detector net3d --backend net3d `
    --selection $G/logo/net3d_rescore_selection.json --out-root $R/net3d
python scripts/bench_logo.py report --heldout net3d=$R/net3d log_tuned=$R/log_tuned `
    log_defaults=$R/log_defaults spotmax_retuned=$R/spotmax_retuned spotmax_deployed=$R/spotmax_deployed `
    spotmax_grid=$G/logo/heldout/spotmax_grid --out $G/logo/report_f1R
```

## nuclei3d

`germquant.segment.nuclei3d` is an EXPERIMENTAL in-house 3D nucleus segmentation model, phase 2 of
`docs/DESIGN_inhouse_components.md` (Component 2): a compact residual 3D U-Net with three heads, trained on
the lab's recovered Imaris nucleus masks with a partial-annotation loss, whose instances come from the
h-maxima of its distance map and a seeded watershed. It is the segmentation backend `nuclei3d` (never the
default; `docs/BACKENDS.md` has its keys), it is trained with `germquant train-nuclei3d`, and it is scored
with `scripts/benchmark_nuclei3d.py` under the protocol's scope rule (`--coverage protocol`, "The nucleus
benchmark" above). Everything here is development on the train and val gonads of
`docs/BENCHMARK_PROTOCOL.md`: the training code never opens a test gonad (every gonad it loads passes the
test guard first), the benchmark script skips test gonads unless `--final`, and no metric of any method
was computed on a test gonad in this package.

### What it learns from

* Train: all six train gonads, n2rad51_HS_herm_000 (88 traced nuclei), HS_male_000 (64), HS_male_002 (82),
  noHS_herm_000 (107), noHS_male_000 (144) and noHS_male_001 (97): 582 traced nuclei. The two heat-shock
  males were exported with Imaris spot region growing, which ruins their spot truth (protocol section 4.1)
  but leaves their nucleus surfaces intact, so unlike net3d this model trains on them.
* Val (early stopping, model selection and the instance-step grid): n2rad51_HS_herm_001 (95 traced nuclei)
  and noHS_herm_001 (146), both hermaphrodites.
* The truth: `<stem>__imaris_labels.tif` of each gonad (`docs/IMARIS_TRUTH.md`), every surface-matched
  label kept whatever its size. The traced masks hold 5,780 to 31,560 voxels per nucleus on the train and
  val gonads (median 11,300 to 19,700 per gonad) and span 14 to 32 planes in z (2.8 to 6.4 um; the axial
  blur makes them about as tall as they are wide), except one 44-voxel label in noHS_male_000.
* The deployed Cellpose-SAM labels of the golden runs (the bench cache's `labels.npy`) decide only where the
  loss is silent; no target value comes from them.
* Input: the DAPI channel of the bench cache, normalised per gonad to (raw - p1) / (p99.9 - p1) with the
  percentiles of the whole stack on every fourth voxel in y and x (p1 is 87 to 89 on every train and val
  gonad, p99.9 471 to 756).

### Targets, loss region and weights (`nuclei3d/data.py`)

* A deployed Cellpose label with at least half of its voxels inside traced masks is the traced nucleus's
  own label; every other label is an untraced nucleus. On the six train gonads, of the 743 Cellpose labels
  that touch a traced mask, 152 have less than 5 percent of their voxels inside one, 7 between 5 and 20
  percent, none between 20 and 40, 16 between 40 and 50 (treated as untraced, the safe side for a label
  that merges a traced nucleus with an untraced one) and 568 at least half; a matched pair puts about 0.7
  there, since the Cellpose masks are about 1.4 times the Imaris volume.
* Classes and loss weights: traced (inside a traced mask) weight 1, target foreground; shell (inside a
  traced nucleus's own Cellpose label, outside its mask) weight 1, target background; ignore (inside an
  untraced Cellpose label or within 0.5 um of one, an exact ellipsoidal dilation, and outside every traced
  mask) weight 0; background (everything else, where neither a traced mask nor a Cellpose label is
  present) weight 1. The shell's weight is a choice the design's rule did not settle (it names traced
  masks, background outside every Cellpose label and the untraced labels): the annotator traced that
  nucleus and its surface ends at the mask, so the truth calls the shell background, and teaching it
  teaches the Imaris boundary. `--w-traced-shell 0` makes it unknown instead; that variant was not
  trained. Sampled on every fourth voxel in y and x, the train gonads are 93.7 to 97.6 percent background,
  1.9 to 5.9 percent ignore, 0.19 to 0.37 percent traced and 0.09 to 0.13 percent shell.
* Targets inside a traced mask (0 everywhere else): foreground 1; the cell distance, the Euclidean
  distance (with the voxel size) from the voxel to the nearest voxel outside its own nucleus, the edge of
  the stack counting as outside (at most 1.74 um on noHS_male_000); the neighbour closeness
  (1 - d / 1.0 um) ** 3 clipped to [0, 1], d the distance to the nearest voxel of another nucleus (another
  traced mask, or a Cellpose label that holds less than half of this mask, standing in for untraced
  neighbours), which tops out at 0.71 at a contact (d is at least one voxel) and reaches 0 at 1 um. On
  noHS_male_000, 101 of the 144 traced nuclei have a neighbour within 1 um.
* Loss: per voxel, the binary cross-entropy of the foreground logit plus the squared error of the cell
  distance (um) plus the squared error of the neighbour closeness, times the weight, summed over the batch
  and divided by its summed weight.
* Patches of 48 x 160 x 160 voxels (9.6 x 17.3 x 17.3 um, a whole nucleus deep), four per step; a gonad
  drawn uniformly, then the centre: a traced nucleus (probability 0.6, plus up to a third of the patch), a
  germline nucleus of the golden run (0.3, plus up to half the patch) or anywhere in the stack; a patch may
  hang over the stack's edge by a quarter (zero input, zero weight there). Augmentation: y and x flips, a
  y/x transpose and a log-uniform gain of 0.8 to 1.25. Each sample draws from its own generator seeded with
  (seed, epoch, index). The target maps are computed once per gonad (18 to 30 s) into the prep folder.

### The network (`nuclei3d/model.py`)

The encoder and decoder of net3d, reused: four levels of one residual block each (two 3 x 3 x 3
convolutions, batch normalisation, SiLU; 16, 32, 64 and 128 channels), an xy-only first pooling (1 x 2 x
2) and 2 x 2 x 2 after it, transposed-convolution upsampling, and a 1 x 1 x 1 head with three outputs
(foreground logit with its bias starting at -2, cell distance, neighbour closeness); 1,421,475 parameters,
a receptive radius of at most 30 voxels in z and 58 in y and x.

### Training, validation and selection (`nuclei3d/train.py`)

AdamW (learning rate 1e-3, weight decay 1e-4), 250 warm-up steps then a cosine decay to 5 percent over the
maximum of 30 epochs of 250 steps, bf16 autocast on CUDA, seeds for python, numpy and torch and the
deterministic settings of net3d (two CPU runs with the same seed give identical weights, tested).

After every epoch both val gonads go through the backend's own path on the whole stack: tiled float32
inference, the instance step with fixed parameters (foreground 0.5, h 0.4 um, neighbour weight 1), the
pipeline's size filter (4 um3), then `germquant.benchmark.nuclei3d.benchmark_gonad` with the protocol's
scope rule. Only the foreground components that touch region C go through the instance step, which
changes no metric (a component that never touches C can only produce out-of-scope predictions) and saves
time; the counts of all predictions in the training log cover those components only. The two gonads are
pooled by summing TP, FP and FN at every IoU threshold. Selection: the best pooled F1 at IoU 0.5 (the
protocol's primary nucleus metric), then the best pooled AP, then the earlier epoch; early stopping after
8 validations without a new best. Every epoch's weights are kept with their sha256. For the selected epoch,
the instance step is then scored on the same two gonads over a grid of h (0.2, 0.4, 0.6, 0.8 um), the
foreground threshold (0.4, 0.5, 0.6) and the neighbour weight (0, 1), 24 configurations, chosen by the same
rule and written into `model.pt`. Together with the validated epochs these are the configurations the
protocol's tuning budget counts (section 6.4, at most 60).

### The instance step and inference (`nuclei3d/infer.py`)

The whole stack is inferred in float32 tiles of 96 x 384 x 384 voxels overlapping by the receptive radius,
each tile contributing its core, so the result does not depend on the tiling beyond float rounding
(tested); no tile is skipped and no other segmentation gates the inference. Foreground: probability above
the threshold. Seed map: S = cell distance x (1 - neighbour weight x closeness), 0 outside the foreground.
Seeds: the h-maxima of S (skimage, full connectivity), each connected group one seed. A seeded watershed
of -S inside the foreground (face connectivity), per connected component of the foreground in its bounding
box, which is exact and keeps the memory at the largest component. The pipeline's `_drop_small` follows in
`segment_nuclei`.

### Running it

```powershell
$env:PYTORCH_CUDA_ALLOC_CONF = "expandable_segments:True"
germquant train-nuclei3d --cache C:/Users/ryane/germquant_golden/bench_cache `
    --split-file tests/golden/rad51_split.json --out C:/Users/ryane/germquant_golden/models/nuclei3d_20261010 `
    --epochs 30 --steps-per-epoch 250 --batch-size 4 --patch 48,160,160 --seed 0 --patience 8 `
    --post-grid "h_um=0.2,0.4,0.6,0.8;fg_threshold=0.4,0.5,0.6;neighbour_weight=0,1" `
    --prep-dir C:/Users/ryane/germquant_golden/models/nuclei3d_prep --device cuda
python scripts/benchmark_nuclei3d.py --manifest <manifest> --split val --coverage protocol --group-by condition `
    --seed 20251105 --model C:/Users/ryane/germquant_golden/models/nuclei3d_20261010/model.pt --device cuda `
    --out <dir>
```

Other options of `train-nuclei3d`: `--lr`, `--warmup-steps`, `--val-every`, `--base-channels`, `--levels`,
`--tile`, `--w-traced-shell`, `--fg-threshold`, `--h-um`, `--neighbour-weight` (the instance step of the
validation), `--min-volume-um3`, `--manifest`, `--min-free-gb`, `--wait-minutes`. Outputs in MODELDIR:
`model.pt` (the selected epoch's weights with the instance-step parameters), `training_log.json` (settings,
gonads, every epoch's losses, timing, memory, checkpoint sha256 and validation, the selection, the grid),
`model_card.json` and `checkpoints/`.

### First training run and its val scores (2026-10-10)

The model: `C:/Users/ryane/germquant_golden/models/nuclei3d_20261010/model.pt`, sha256
`11920205258dceece3e2bd6053621cbce390177f2baa70b9f40eefc0980bbb83` (5,735,871 bytes), trained by commit
02d8389 with seed 0 (21:31 to 22:45 UTC, 4458 s on the RTX 5090, of which 1553 s were the instance-step
grid). It ran 16 epochs and stopped early (no new best in 8 validations); the selected epoch is 8, with the
learning rate still at 0.87 of its peak. The grid kept the training's own instance-step parameters
(foreground 0.5, h 0.4 um, neighbour weight 1). 40 configurations were scored on val in all (16 epochs, 24
grid points), within the protocol's budget of 60. Timing: 56 to 57 s per training epoch (5.0 GB of GPU
memory allocated, 7.2 GB reserved), 114 to 162 s per validation of both gonads (13.4 to 13.7 GB reserved during
the inference): 29 to 33 s of inference per gonad, 14 to 45 s for the instance step and 6 to 8 s for the scoring. The
console log is `train_console.log` next to the model; a one-epoch smoke run of the same code on the real
data, made before it, is kept in `C:/Users/ryane/germquant_golden/models/nuclei3d_work/dryrun/`.

Pooled val scores by epoch (protocol rule, F1 at IoU 0.5 / AP over 0.5 to 0.95; the per-gonad F1 of the
heat-shock and the no-heat-shock hermaphrodite after the bar):

| epoch | 1 | 2 | 3 | 4 | 5 | 6 | 7 | **8** | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| F1 | 0.895 | 0.940 | 0.895 | 0.915 | 0.940 | 0.872 | 0.904 | **0.968** | 0.956 | 0.908 | 0.960 | 0.947 | 0.891 | 0.940 | 0.956 | 0.923 |
| AP | 0.437 | 0.502 | 0.531 | 0.571 | 0.570 | 0.529 | 0.472 | **0.582** | 0.583 | 0.500 | 0.602 | 0.622 | 0.548 | 0.597 | 0.596 | 0.586 |
| F1, HS_herm_001 | 0.879 | 0.844 | 0.747 | 0.790 | 0.857 | 0.720 | 0.821 | **0.923** | 0.906 | 0.802 | 0.900 | 0.888 | 0.774 | 0.872 | 0.917 | 0.828 |
| F1, noHS_herm_001 | 0.905 | 0.997 | 0.979 | 0.986 | 0.990 | 0.956 | 0.952 | **0.997** | 0.986 | 0.969 | 0.997 | 0.983 | 0.959 | 0.983 | 0.979 | 0.979 |

The no-heat-shock gonad is at 0.95 to 0.997 from epoch 2 on; the heat-shock gonad swings between 0.72 and
0.92 from one epoch to the next, almost all of it in recall (0.61 to 0.88), and that swing decides the
pooled F1. The selected epoch is the top of a noisy curve: the median pooled F1 of epochs 8 to 16 is
0.947, about 0.02 below the selected 0.968. AP kept rising after the selected epoch (0.622 at epoch 12).

The instance-step grid at the selected epoch (pooled over both val gonads; the 24 rows are in the training
log):

| h (um) | foreground 0.4 / 0.5 / 0.6, neighbour weight 1: F1 | the same with neighbour weight 0: F1 | splits / merges (0.5, weight 1) |
|---|---|---|---|
| 0.2 | 0.920 / 0.912 / 0.912 | 0.938 / 0.934 / 0.922 | 18 / 0 |
| 0.4 | 0.964 / **0.968** / 0.968 | 0.958 / 0.962 / 0.962 | 0 / 1 |
| 0.6 | 0.881 / 0.885 / 0.885 | 0.872 / 0.877 / 0.877 | 0 / 1 |
| 0.8 | 0.785 / 0.794 / 0.790 | 0.785 / 0.794 / 0.790 | 0 / 5 |

h is the parameter that matters: 0.2 um splits nuclei (16 to 18 splits, precision 0.88 to 0.91), 0.6 um and
above loses nuclei (recall 0.82 to 0.83 at 0.6, 0.71 at 0.8). The foreground threshold moves F1 by at most 0.016.
The neighbour head helps a little at h 0.4 (229 against 227 matched nuclei) and hurts at h 0.2, where its
damping of the distance map at contacts adds seeds.

### Against the deployed Cellpose-SAM on val

`scripts/benchmark_nuclei3d.py --coverage protocol` on the two val gonads: the deployed labels are the
8666f07 goldens' `<image>__nuclei_labels.tif` (Cellpose-SAM fine-tuned as `germline_nuclei_combined`, the
pipeline's own path); nuclei3d is the frozen model through `segment_nuclei` on the whole stack (`--model`,
CUDA, the size filter of 4 um3). Outputs under
`C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/nuclei3d_val/` (`cellpose_golden_protocol`,
`nuclei3d_20261010_protocol`, and `cellpose_golden_auto`, `nuclei3d_20261010_auto` for the earlier rule).
The whole-stack metrics of nuclei3d equal the training-time validation exactly (TP, FP, FN, AP and the
z-extent error on both gonads). Intervals: cluster bootstrap by gonad, 1000 resamples, seed 20251105; with
two gonads they are the range of the two.

| val, protocol rule, pooled over 2 gonads | deployed Cellpose-SAM | nuclei3d |
|---|---|---|
| F1 at IoU 0.5 | 0.863 [0.863, 0.863] | 0.968 [0.923, 0.997] |
| precision / recall at IoU 0.5 | 0.863 / 0.863 | 0.987 / 0.950 |
| matched / in-scope predictions / truth | 208 / 241 / 241 | 229 / 232 / 241 |
| AP over IoU 0.5 to 0.95 | 0.338 [0.311, 0.357] | 0.582 [0.493, 0.645] |
| splits / merges (20 percent rule) | 0 / 0 | 0 / 1 |
| z-extent error of matched pairs, mean | +0.51 um [+0.44, +0.61] | +0.09 um [-0.14, +0.48] |
| predictions in the whole stack / ignored as out of scope | 1887 / 1646 | 1003 / 771 |
| segmentation time per gonad on the RTX 5090 | about 19 min of Cellpose-SAM network time (golden run logs) | 50 to 72 s (inference, instances, size filter) |

| val gonad | method | in scope / truth | TP / FP / FN | F1 | recall | AP | mean matched IoU | z-extent error, mean (mean absolute) |
|---|---|---|---|---|---|---|---|---|
| n2rad51_noHS_herm_001 | Cellpose-SAM | 146 / 146 | 126 / 20 / 20 | 0.863 | 0.863 | 0.357 | 0.737 | +0.44 um (0.48), +2.2 planes |
| n2rad51_noHS_herm_001 | nuclei3d | 145 / 146 | 145 / 0 / 1 | 0.997 | 0.993 | 0.645 | 0.828 | -0.14 um (0.30), -0.7 planes |
| n2rad51_HS_herm_001 | Cellpose-SAM | 95 / 95 | 82 / 13 / 13 | 0.863 | 0.863 | 0.311 | 0.715 | +0.61 um (0.62), +3.0 planes |
| n2rad51_HS_herm_001 | nuclei3d | 87 / 95 | 84 / 3 / 11 | 0.923 | 0.884 | 0.493 | 0.801 | +0.48 um (0.57), +2.4 planes |

Per gonad, nuclei3d minus Cellpose: F1 +0.134 and +0.060, AP +0.287 and +0.182, recall +0.130 and +0.021,
mean matched IoU +0.092 and +0.087 (no-heat-shock, heat-shock). Under the earlier rule (the ROI box) the
same labels give F1 0.694 and 0.514 (Cellpose, pooled 0.610, AP 0.205) against 0.824 and 0.680 (nuclei3d,
pooled 0.765, AP 0.396): the box holds 91 and 142 unmatched Cellpose nuclei and 61 and 68 unmatched nuclei3d
ones, untraced nuclei the protocol's rule leaves out.

**Device check.** The same model on the CPU (4 threads; `nuclei3d_20261010_cpu_check`) gives the
identical label image for n2rad51_noHS_herm_001, voxel for voxel (476 nuclei, 7,748,740 labelled voxels),
in 898 s of segmentation instead of 50 s.

Where the two differ, from the val label images:

* Cellpose finds every traced nucleus and draws it larger. Its in-scope predictions equal the truth count on
  both gonads (146 and 95), and each of the 33 traced nuclei it misses at IoU 0.5 lies 95 to 100 percent
  inside one Cellpose nucleus with an IoU of 0.39 to 0.50: near misses of the boundary convention (its
  masks are about 1.4 times the Imaris volume, protocol section 4.2, and run 0.44 to 0.61 um longer in z),
  not missed nuclei. nuclei3d learned the Imaris convention it is scored against, which is most of its gain
  in F1, AP and matched IoU. The decisive comparison, Cellpose fine-tuned on the same masks, is below.
* nuclei3d merges crowded nuclei, which the protocol's merge count cannot see. On the heat-shock gonad each
  of the 11 traced nuclei it misses lies 94 to 100 percent inside a nuclei3d object 2.5 to 10.2 times its
  volume, 27 to 54 planes tall where the traced nucleus spans 18 to 28, that covers at least 20 percent of
  2 to 9 deployed Cellpose nuclei: merges, mostly along z, with untraced neighbours. 7 of these 11 objects
  have less than half of their volume in region C and are ignored by the scope rule, so they cost recall
  but no precision, and the merge count (20 percent of two or more traced nuclei) sees one. Over the whole
  stacks nuclei3d makes 476 and 527 objects where Cellpose makes 617 and 1270, and 95 and 211 of its objects
  cover at least 20 percent of two or more Cellpose nuclei (none and 5 of them in scope). Outside the
  traced region its training had no supervision on untraced nuclei (zero weight there), and these numbers
  say that it under-segments dense untraced regions, the heat-shock gonad most of all. Its masks are not
  ready to replace the deployed ones in the pipeline, whatever the traced-subset scores say.

### The fine-tuned Cellpose competitor (on the same masks)

The design's decisive competitor: Cellpose fine-tuned on the masks nuclei3d learns from. Run as the
stretch of this package, once, with no tuning:

* Slices: `scripts/export_imaris_training_pairs.py` wrote 144 slice pairs from the six train gonads (8 per
  gonad and orientation, xy, xz and yz, each cutting at least 3 traced nuclei, seed 0; xz and yz stretched
  along z to the xy pixel size; each cropped to its traced labels plus 48 pixels; 2152 nucleus sections in
  all) to `C:/Users/ryane/germquant_golden/models/cellpose_ft_imaris_pairs_20261010/`. The labels are the
  traced masks; the voxels nuclei3d leaves out of its loss (untraced Cellpose labels plus 0.5 um) were filled
  with the gonad's 1st DAPI percentile, since Cellpose's loss has no ignore label (9.1 to 28.4 million voxels
  per gonad). That fill changes the image the network sees (an untraced nucleus becomes a dark hole), a limit
  of this competitor.
* Training: `germquant train-nuclei --pairs <pairs> --out C:/Users/ryane/germquant_golden/models/cellpose_ft_imaris_20261010
  --base <the deployed germline_nuclei_combined file> --epochs 100 --split-by series --min-train-masks 3
  --batch-size 8 --learning-rate 1e-5 --seed 0` on the RTX 5090 (992 s; one series group, so all 144 pairs
  trained and nothing was held out). Model `cellpose_ft_imaris_20261010`, sha256
  `a4334d0016f3ce548c92ddb722011be908c922d576e5d740c1405c714561c38c`; base model sha256
  `e4e32320cfc321a145ccbc0c52138eded50e37b87130f34bae2296c88b148f57`; Cellpose 4.2.1.1.
* Scoring: `scripts/benchmark_nuclei3d.py --cellpose-model <model> --coverage protocol --split val`, the
  pipeline's own Cellpose path in 3D on the cached DAPI with the size filter of 4 um3 (1140 to 1180 s of
  segmentation per gonad); outputs `cellpose_ft_imaris_protocol` and `cellpose_ft_imaris_auto` next to the
  others.

| val, protocol rule, pooled over 2 gonads | deployed Cellpose-SAM | nuclei3d | Cellpose fine-tuned on the masks |
|---|---|---|---|
| F1 at IoU 0.5 | 0.863 [0.863, 0.863] | 0.968 [0.923, 0.997] | 0.974 [0.967, 0.984] |
| precision / recall | 0.863 / 0.863 | 0.987 / 0.950 | 0.952 / 0.996 |
| TP / FP / FN | 208 / 33 / 33 | 229 / 3 / 12 | 240 / 12 / 1 |
| AP over IoU 0.5 to 0.95 | 0.338 [0.311, 0.357] | 0.582 [0.493, 0.645] | 0.590 [0.554, 0.656] |
| splits / merges | 0 / 0 | 0 / 1 | 1 / 0 |
| mean matched IoU, noHS / HS | 0.737 / 0.715 | 0.828 / 0.801 | 0.798 / 0.840 |
| z-extent error, mean | +0.51 um | +0.09 um | -0.25 um [-0.32, -0.15] |
| F1, noHS_herm_001 / HS_herm_001 | 0.863 / 0.863 | 0.997 / 0.923 | 0.967 / 0.984 |
| AP, noHS_herm_001 / HS_herm_001 | 0.357 / 0.311 | 0.645 / 0.493 | 0.554 / 0.656 |
| objects in the whole stack, noHS / HS (median volume) | 617 / 1270 (35.5 / 39.4 um3) | 476 / 527 (34.2 / 50.3 um3) | 647 / 1234 (23.2 / 27.2 um3) |
| earlier rule (ROI box): pooled F1 / AP | 0.610 / 0.205 | 0.765 / 0.396 | 0.688 / 0.338 |
| segmentation time per gonad (RTX 5090) | about 19 min (network) | 50 to 72 s | 1140 to 1180 s |

On the traced nuclei the fine-tuned Cellpose ties with nuclei3d: pooled F1 0.974 against 0.968 and AP 0.590
against 0.582, differences far below the protocol's minimum effect of 0.03, and the two split the gonads
(nuclei3d ahead on the no-heat-shock gonad by 0.030 in F1, the fine-tuned Cellpose ahead on the heat-shock
gonad by 0.061, where it finds all 95 traced nuclei). The one traced nucleus it misses has a best IoU of
0.26. It keeps the deployed model's object count over the whole stacks (647 and 1234 against 617 and 1270)
where nuclei3d makes 476 and 527, while the total foreground volume of the two new models is the same on
the heat-shock gonad (37,332 and 37,351 um3 against 55,136 for the deployed labels): both draw the Imaris
boundary, and nuclei3d packs it into fewer, larger objects. Against the deployed labels as a reference
(not a truth), 193 of the fine-tuned Cellpose's 1234 heat-shock objects also cover at least 20 percent of
two or more deployed nuclei (7 of 647 without heat shock), against 211 of nuclei3d's 527 (95 of 476), so
neither is free of merges there, but nuclei3d merges far more of its objects. Under the design's stop rule
(a fine-tuned competitor at least the minimum gain ahead on val ends the iteration) this is a tie on the
protocol's metric, not a stop; on the whole-stack object counts the fine-tuned Cellpose is the better
segmenter, and on runtime nuclei3d is about 20 times faster.

### What this does and does not say

* These are development numbers on two val gonads, both hermaphrodites. The epoch and the instance-step
  parameters were chosen on them among 40 configurations, so the scores are optimistic: the selected epoch
  is the top of a curve whose later epochs have a median pooled F1 of 0.947, and the heat-shock gonad's F1
  swings by 0.2 from epoch to epoch. No test gonad was opened or scored, by nuclei3d or by any other method.
* On the traced nuclei, scored by the protocol's rule, nuclei3d agrees with the Imaris masks better than the
  deployed Cellpose-SAM on both val gonads (pooled F1 0.968 against 0.863, AP 0.582 against 0.338), but
  Cellpose fine-tuned on the same masks does as well (0.974 and 0.590): the gain over the deployed model
  comes from the training masks, the Imaris boundary convention, not from the architecture.
* The protocol's partial-annotation scope rule is lenient on objects that reach far outside the traced
  region, and the whole-stack diagnostic above shows nuclei3d under-segmenting untraced crowded nuclei. The
  traced-subset metrics cannot show this; the RAD-51 counts on its masks (not computed), a complete
  annotation of a crop (the design's napari crops) or the blind human subset can.
* W3 of the protocol is a test-set criterion after retraining on train plus val (section 6.1); neither was
  done. Nor was leave-one-gonad-out (section 6.2), the `--w-traced-shell 0` variant, or a second seed.
* No heat-shock male is in val; the train set holds both heat-shock males (their surfaces are usable), so
  unlike the spot detectors this model saw heat-shocked males in training.
* Runtime is a clear difference: about a minute per gonad against about twenty for Cellpose-SAM in 3D on the
  same GPU.

## Tests

`tests/test_benchmark.py`: synthetic label volumes and spot sets with known answers for every metric
(a split, a merge, an out-of-coverage prediction under the bounding-box rule and the ROI rule, the
z-extent error, the AP arithmetic), the chunked Hungarian against the dense reference (and the classic
case where greedy matching loses a pair), the greedy fallbacks and their flags, bootstrap intervals
that contain the point estimate and repeat for the same seed, the manifest and truth loaders, the
frame conversion, both scripts run on a temporary manifest in a subprocess, and the .ims fallback
adapter with a stubbed reader.

`tests/test_bench_detector.py`: the detector bench harness. A real pipeline run on a synthetic uint16
gonad is cached and re-run: the cached label image equals the one the spots stage received,
`run_backend` hands `detect_spots` exactly the pipeline's arguments (image, labels, spacing, the nine
keyword arguments, marker, backend, `cfg_section`) and reproduces its spots table (log backend always,
SpotMAX when installed); the keyword arguments follow the stage for non-default values too; the
params merge equals `pipeline.spot_backend_section`. Tiny planted-foci gonads with truth files run the
runner end to end with the log backend (scores, Class A columns, conditions under `by_split`), the
trial log (one line per run, the params sha256), the test guard in the script and in the library, the
final-run log and its limit of two (on synthetic test-role gonads and temporary logs only), and the
process pool against the sequential run.

`tests/test_net3d.py`: net3d on the CPU with tiny networks (skips itself without torch). The network's
shapes and receptive-radius bound (a perturbed voxel never reaches farther), tiled inference equal to
untiled inference in xy and z and with the tile skip, the normalisation, peak finding on a synthetic
heatmap (log-parabola refinement, ellipsoidal suppression, the threshold property the validation sweep
relies on), the suppression against the plain greedy and the candidate scan against its definition
(dense and sparse blocks, exact ties), the reach grid against assigning every peak, nucleus assignment
with the nearest-label fallback and the contested flag, the selection rule, the exact dilations, the
targets, weights and masked loss, reproducible sampling, the Imaris-to-grid conversion of the truth
positions, the backend contract with a randomly initialised model (the DNA channel through `channels`,
a foci-only model, the provenance hook, thresholds from zero to above every probability, deterministic output), channels
reaching only backends that ask for them, the default run's stage hashes and manifest unchanged, a
pipeline run with `spots.backend: net3d`, and `germquant train-spots` for one epoch on tiny planted gonads
(the test gonad and the region-growing gonad are never opened, the same seed gives the same weights, the
validation scorer agrees with the harness, a split file that moves a gonad is refused, and
`scripts/bench_detector.py --backend net3d` scores the trained model from the cache).

`tests/test_bench_tuning.py`: the parameter search runner. The SpotMAX effect-size filter applied to a
run made with `effect_size_min` 0 reproduces runs made with minima from 0 to above the largest effect
size, table for table, and `run_configs` gives the same tables and scores with the shared detection as
without it (SpotMAX when installed); the ranking rule (rounding, the F1 tie-break, missing values last)
and the chunking that keeps a shared SpotMAX group in one task; on the planted-foci gonads,
`bench_tune.py` writes one trial record per parameter set with its tuning block and the same per-gonad
rows, pooled blocks and detections as `bench_detector.py` for the same parameters, the process pool with
two tasks per gonad equals the sequential run, the budget guard counts distinct parameter sets per search,
`--top-from` sends the best sets to val, and test gonads are refused by the script and the library. The
two tuned profiles extend `config/config.yaml`, change only `spots:` keys and, used as bench parameters
over the base section, give exactly the profile's section. `tests/test_spot_backends.py` covers the
`spots.log.pre_gate_sigmas` key (default, validation, that it reaches the candidate step and gates).

`tests/test_bench_partition.py`: the section 5.1 partition. The assignment rule on tiny label images (centre
voxel, the 3 x 3 x 3 mode with ties to the smaller label, unassigned beyond, the image edge), fragments removed
everywhere through the lookup table, the assignment equal to the spot region R voxel by voxel, the
per-truth-nucleus counts for Class A+B, A and B with assigned and unassigned spots and the detector-only column,
the `tm_` statistics against `agreement_metrics` and the mean absolute error, the paired cluster bootstrap (the
same draws as `cluster_bootstrap`, a method against itself, a better method), the independence of the partition
from our own segmentation on a planted gonad, and `bench_detector.py --detector-only --note` end to end (the
detector-only run on truth masks equal to the cached nuclei gives the same spots, the pooled blocks, the trial
record's blocks and note).

`tests/test_bench_logo.py`: leave-one-gonad-out. The fold rule (designated val gonads, the same-condition
substitute by tag order, no test gonad and no gonad without usable spot truth, the real split's six folds), the
D3 stage generators, the staged replay on a synthetic score index (stage by stage, naming what is missing, the
selection following val among the finalists), the index refusing two disagreeing scores of one gonad, and
`bench_logo.py folds`, `heldout` and `report` end to end on the planted gonads (one trial record per held-out
gonad with its fold, the pooled truth-mask CCC equal to the concatenated tables, the paired difference, `report
--tags` from one folder or per-gonad folders, and the pooled F1 inside R equal to the F1 of the summed
`tm_spots_` counts, with held-out tables that lack those columns left out of that block). The pre-registered
SpotMAX grid: its 60 sets, their order and names, the deployed recipe among them and only `spots:` keys
changed; the primary objective equal to the point
estimates `detector.pool_scores` gives over the same gonads (a gonad without truth nuclei or with NaN counts
inside R included) and to `bench_tune.py`'s own pooling of the folders it writes (log backend); the truth-mask
counts in the run index (read only on request, refused when two runs disagree); the grid rule (the rounding, the
F1 and order tie-breaks, NaN values last, no selection while a pair is missing or lacks the truth-mask counts, no
test gonad); each fold's own selection gonads on the real split; and `grid-configs` and `select --search
prereg-grid` end to end (a file `bench_tune.py` accepts, exit code 3 while pairs are missing, a selection that
`heldout --selection` reads, the D3 replay unchanged). On the real split's folds, `to_run` merges the two
substitutes of the val folds, the substitutes alone complete the selection, and disagreeing folders or a split
without val gonads are refused with exit code 2.

`tests/test_net3d_recipe.py`: the frozen recipe. Retraining by `--recipe-from` on the selected model's own
gonads gives the selected checkpoint bit for bit with the same learning rates, the card and the log record the
recipe, the train-plus-val run trains on the val gonad and never opens the test or region-growing gonad,
`bench_logo.py net3d-folds` trains each fold on the recipe's other gonads and scores its held-out gonad, and the
refusals (a test gonad, an unusable gonad, a gonad both trained and validated on, no val gonad without a fixed
epoch count and threshold, a threshold with val gonads, `--recipe-from` with `--threshold`).

`tests/test_nuclei3d.py`: nuclei3d and the protocol scope, on the CPU with tiny networks (skips itself
without torch). The network's shapes, heads and receptive-radius bound; the loss (the weights, a
zero-weight voxel changes nothing whatever the network says there, the normalisation, an all-zero weight);
the targets on synthetic spheres (the cell distance equals the Euclidean distance transform inside each
traced mask, the neighbour closeness is high at a contact and zero beyond the radius, a nucleus's own
larger Cellpose label is not its neighbour, the class map: traced, shell, ignore inside an untraced
Cellpose label and its 0.5 um margin, background) and their cache round trip; reproducible, weighted
patches with zero weight off the stack; the instance step splitting two touching spheres (with and
without the neighbour damping), one sphere staying one, an h above every dome or an empty foreground giving
nothing; tiled inference equal to untiled inference in xy and z, and the seed-map path equal to the maps
path; the backend contract with a randomly initialised model (segment_nuclei, the config keys over the
model file's parameters, the classical fallback under `auto`, a missing model file, a file that is not a
nuclei3d model, the provenance and its sha256, nothing for cellpose and classical); the default run's stage
hashes and segment provenance unchanged and the segment sub-hash changed by nuclei3d; a pipeline run with
`segmentation.nuclei.backend: nuclei3d` (`segmentation_method` and `segment_backend` in the manifest);
`germquant train-nuclei3d` for one epoch on tiny planted gonads with a post grid (the test gonad is never
opened, the same seed gives the same weights, the log, card and model file agree, a split file that moves
a gonad is refused) and `scripts/benchmark_nuclei3d.py --model` and `--pred-dir` on the result (the test
gonad skipped, the training-time grid scores equal to the script's); with an oracle in place of the
network, the validation scorer equal to the script and perfect; `--cellpose-model` with a stubbed Cellpose;
the protocol mode (fragments removed, the half-volume scope rule at its edge, the 20 percent splits and
merges, the z-extent error in slices); and `scripts/export_imaris_training_pairs.py` (train gonads only,
the truth masks as labels, one series so `train-nuclei` trains on every pair).
