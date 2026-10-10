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

Options of `benchmark_nuclei3d.py`: `--split`, `--coverage auto|roi|bbox|none` (the coverage region
of the partial-annotation rule), `--dilate-um 2.0`, `--majority-frac 0.5` (splits and merges),
`--max-dense-cells`, `--n-boot`, `--seed`.

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
`voxel_dx_um`, `voxel_source`, `roi_box`, `elapsed_s`.

`per_nucleus.csv`, one row per TRUTH nucleus: `tag`, `split`, `truth_label`, `truth_voxels`,
`truth_volume_um3`, `truth_z_span_um`, `matched_iou50`, `pred_label` (0 when unmatched), `iou`,
`pred_voxels`, `pred_z_span_um`, `z_extent_error_um`, `n_pred_majority_inside`, `is_split`.

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

## Tests

`tests/test_benchmark.py`: synthetic label volumes and spot sets with known answers for every metric
(a split, a merge, an out-of-coverage prediction under the bounding-box rule and the ROI rule, the
z-extent error, the AP arithmetic), the chunked Hungarian against the dense reference (and the classic
case where greedy matching loses a pair), the greedy fallbacks and their flags, bootstrap intervals
that contain the point estimate and repeat for the same seed, the manifest and truth loaders, the
frame conversion, both scripts run on a temporary manifest in a subprocess, and the .ims fallback
adapter with a stubbed reader.
