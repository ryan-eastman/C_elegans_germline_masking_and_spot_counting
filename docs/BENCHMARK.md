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
(0.3 and 0.4 um), `--verbose` (keeps SpotMAX's progress bars).

Outputs in `--out`: `per_gonad.csv`, `per_nucleus.csv`, `per_spot_matches.csv`, `pooled.json` (the
pooled blocks above with `backend`, `params`, `params_sha256`, `split`, `skipped`, `git_sha`,
`manifest`, `cache_root` and the bootstrap settings) and `detections/<tag>__spots.csv`, the detector's
own spot table per gonad.

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
  The truth spots were detected over the whole field, but the per-gonad Quality threshold was set while
  looking at the nuclei, and the counts that matter are per nucleus; outside this region the weight is 0
  and the network is unconstrained there (it does put low probabilities on background, which the
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
set. Two launches with the same seed gave bit-identical losses and validation numbers.

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
higher CCC, the earlier epoch, the lower threshold). 0.01 is the protocol's scale for a change of the
primary metric that counts (section 10). The first launches used a strict CCC maximum with F1 only for
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
assign, the identical 1269 assigned). On the same run the validation numbers of four epochs were
recomputed with the old and the new peak finder: identical at every threshold, pooled and per gonad.
At a working threshold (0.3) the whole chain takes about 3 s per gonad.

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
`--select-tolerance`, `--base-channels`, `--levels`, `--manifest`, `--wait-minutes`. `--min-free-gb`
waits (polling every 60 s) until nvidia-smi reports that much free memory; on Windows,
`torch.cuda.mem_get_info` reports the whole card as free while other processes hold most of it, so it is
not used for this. The loss-region maps (one uint8 volume per train gonad, 0.45 GB each) are computed
once into `--prep-dir` and reused by later runs. Outputs in MODELDIR: `model.pt`, `training_log.json`
(settings, gonads, every epoch's loss, timing, memory, checkpoint sha256 and the full validation sweep,
the selection), `model_card.json` (training, val and excluded gonads, truth classes, epochs, the
selection, the val metrics, sha256, git sha, limitations) and `checkpoints/`.

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
