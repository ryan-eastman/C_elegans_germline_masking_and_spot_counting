# RAD-51 benchmark protocol (pre-registration)

Registered 2026-10-09. No benchmark result has been computed under this protocol; the changelog at the
end says what had been looked at before it was written. The split record is
`tests/golden/rad51_split.json`, its code is `germquant.validate.benchmark_split`,
`scripts/make_rad51_split.py` rebuilds the record and `scripts/benchmark_manifest.py` writes
`benchmarks/<name>/manifest.json`, the input list for the scoring step.

## 1. What this document is for

The pipeline counts RAD-51 foci per nucleus in the C. elegans germline. The lab's reference for those
counts is Imaris: automatic Spots and Surfaces, curated with hand-set thresholds, on 12 N2 gonads
imaged on 2025-11-05. This protocol fixes, before any result exists, which gonads may be used for
training and tuning and which are held out, what the truth is and what it is not, how agreement is
measured, what would count as an improvement over the deployed pipeline, which competitors are run,
and when to stop. The point is that the numbers reported later cannot have been shaped by a choice of
metric, split or threshold made after seeing them.

Rules of the document:

1. A result is any metric of section 5 computed on any test gonad for any method, the deployed
   pipeline included.
2. Until the first result, the sections above the changelog may be edited; each edit gets a line in
   the changelog.
3. From the first result on, nothing above the changelog is edited. A change of method, metric or
   data is appended to the changelog with the date, the reason, and whether a test result had been
   seen at that point. A change made after a test result is labelled post hoc there and in every
   write-up that uses it.
4. The split never changes. A gonad leaves the benchmark only for a defect on the truth side (a
   corrupted file, an annotation that does not belong to the image), recorded in the changelog,
   never because of a score.

## 2. The data

Twelve gonads from one imaging session (2025-11-05), N2, three channels (405 DAPI, 477 SYP, 545
RAD-51), 2600 x 2600 x 71 voxels without heat shock and 2600 x 2600 x 67 with heat shock, voxel
0.108333 x 0.108333 x 0.2 um, 16-bit Nikon `.nd2` stacks (2.88 GB and 2.72 GB each). Each stack has
an Imaris project (`.ims`, 1.5 to 2.8 GB) with the lab's Spots and Surfaces and a fourth channel,
`Masked 405`, that holds the surface masks (section 4.2). The size and sha256 of all 24 files are in
the split record; `check_split_files` (and the test
`test_split_paths_exist_with_the_registered_sizes`) confirm that the data on disk is the data that
was registered.

| tag | condition | sex | image id (`.nd2` stem) | Imaris spots | surfaces | spot Quality threshold | region growing | role |
|---|---|---|---|---|---|---|---|---|
| n2rad51_noHS_herm_000 | noHS | herm | `20251105_N2_nohs_HERM _` | 1845 | 107 | 42.5 | no | train |
| n2rad51_noHS_herm_001 | noHS | herm | `20251105_N2_nohs_HERM _001` | 1222 | 146 | 42.1 | no | val |
| n2rad51_noHS_herm_002 | noHS | herm | `20251105_N2_nohs_HERM _002` | 1300 | 93 | 60.6 | no | test |
| n2rad51_noHS_male_000 | noHS | male | `20251105_N2_nohs_MALE _` | 930 | 144 | 68.7 | no | train |
| n2rad51_noHS_male_001 | noHS | male | `20251105_N2_nohs_MALE _001` | 1256 | 97 | 43.3 | no | train |
| n2rad51_noHS_male_002 | noHS | male | `20251105_N2_nohs_MALE _002` | 817 | 130 | 102.3 | no | test |
| n2rad51_HS_herm_000 | HS | herm | `20251105_n2_hs_HERM _` | 5152 | 88 | 24.2 | no | train |
| n2rad51_HS_herm_001 | HS | herm | `20251105_n2_hs_HERM _001` | 2753 | 95 | 17.4 | no | val |
| n2rad51_HS_herm_002 | HS | herm | `20251105_n2_hs_HERM _002` | 2364 | 84 | 16.9 | no | test |
| n2rad51_HS_male_000 | HS | male | `20251105_n2_hs_MALE _` | 772485 | 64 | 13.4 | yes | train |
| n2rad51_HS_male_001 | HS | male | `20251105_n2_hs_MALE _001` | 10411 | 68 | 8.1 | no | test |
| n2rad51_HS_male_002 | HS | male | `20251105_n2_hs_MALE _002` | 13951 | 82 | 6.8 | yes | train |

The spot counts, surface counts, Quality thresholds and region-growing flags were read from the
creation parameters stored in each `.ims` (`Scene8/Content/Points0/CreationParameters`,
`mEnableRegionGrowing`) with h5py; the split record carries them per gonad. 1198 surfaces in all.

The deployed pipeline, for every comparison below, is the frozen output tree at commit 8666f07
(`C:/Users/ryane/germquant_golden/8666f07/<tag>/`, `config/config.yaml` with config hash
84bdda5a5230, nucleus model `germline_nuclei_combined`, SpotMAX with `threshold_triangle` and
`effect_size_min` 3.0).

## 3. The split and why

| role | gonads |
|---|---|
| test (4) | n2rad51_noHS_herm_002, n2rad51_noHS_male_002, n2rad51_HS_herm_002, n2rad51_HS_male_001 |
| val (2) | n2rad51_noHS_herm_001, n2rad51_HS_herm_001 |
| train (6) | n2rad51_noHS_herm_000, n2rad51_noHS_male_000, n2rad51_noHS_male_001, n2rad51_HS_herm_000, n2rad51_HS_male_000, n2rad51_HS_male_002 |

Why this split:

- Every condition-by-sex cell appears in the test set once, so the headline result covers both
  conditions and both sexes and no cell can be left out by tuning.
- Replicate 002 is the test gonad in three cells. In the heat-shock male cell the two other gonads
  (000 and 002) were exported with Imaris region growing on, which ruins their spot truth (772,485
  and 13,951 spots); 001 is the only heat-shock male with a clean spot export, so it is the test
  gonad there, and the two region-growing gonads go to train, where the defect can only cost training
  signal, never a test score.
- n2rad51_noHS_male_002 is in the test set although the deployed pipeline cannot process it today
  (the spot detector floods; section 10 says how such a gonad is reported). Its truth is sound, and a
  gonad that is hard for the current pipeline is exactly the one that must not be quietly left out of
  the test.
- Val holds one no-heat-shock and one heat-shock hermaphrodite, so thresholds and early stopping see
  both conditions. It holds no male: anything tuned on val is tuned on hermaphrodites, and the males
  in train are the only male signal a learned component sees before test. That is a limitation of a
  12-gonad set, stated here rather than discovered later.
- Train holds the remaining six, including both region-growing gonads (usable as nucleus truth, not
  as spot truth).

With 4 test gonads every confidence interval in section 6 is wide. The leave-one-gonad-out analysis
over the 8 train and val gonads is the supporting evidence; it does not replace the test set.

## 4. The truth and its caveats

### 4.1 Spots

Imaris automatic Spots on channel 545 (RAD-51): estimated diameter 0.4 x 0.4 x 0.9 um, background
subtraction on, detected over the whole field (the Spots creation parameters carry no region of
interest, and the spot positions span the full 2600 x 2600 field), then kept by a per-gonad manual
Quality threshold set by eye (6.8 to 102.3 across the twelve files, table above). Every spot is stored
with radius 0.2 um. The Quality values themselves are not stored in the `.ims` statistics tables, only
the threshold that was applied, so the existing projects cannot be re-thresholded; that needs Imaris.

Each spot carries a Class, A or B, from Imaris's one-dimensional classification on the average
distance to the 3 nearest neighbouring spots: A below a per-gonad threshold (3.2 to 12.6 um), B
above it. Class B therefore means isolated, not false: a nucleus with a single focus in a sparse
region produces a Class B spot. In n2rad51_noHS_herm_001 that is 1142 A and 80 B of 1222 spots; in
n2rad51_HS_herm_000, 3981 A and 1171 B of 5152.

Consequences:

- The spot truth is a thresholded automatic detection, not a hand-verified count. It is the lab's
  reference, so it is what the benchmark measures agreement with, and the human blind subset (section
  9) says how far two people agree with each other and with it.
- Counts and matching are reported for Class A+B (primary), Class A alone and Class B alone
  (secondary). A method that agrees with A+B but not with A alone has a different sensitivity to
  isolated foci, which is worth knowing.
- The two region-growing exports (n2rad51_HS_male_000 and n2rad51_HS_male_002) are excluded from
  every spot metric until they are re-exported from Imaris without region growing (the re-export gets
  a new size and hash in the split record and a changelog line). Their surfaces are unaffected and
  stay in the nucleus truth.

### 4.2 Surfaces (nuclei) and their masks

Imaris automatic Surfaces on channel 405 (DAPI) inside a hand-drawn region of interest, with surface
smoothing 0.4 um, background subtraction (largest sphere 2.5 um), a per-gonad manual absolute
intensity threshold, seeded region growing (estimated diameter 2.5 um) and a Number of Voxels
filter. The result is 64 to 146 surfaces per gonad, 1198 in all: a subset of the nuclei in the field
(the deployed pipeline finds 380 to 1240 germline nuclei per gonad). The annotator traced enough
nuclei for the lab's statistics, not every nucleus. The region-of-interest box recorded in the
Surfaces creation parameters is not a reliable record of the final extent (in n2rad51_noHS_herm_001
the box ends at y = 1442 while the surfaces reach y = 2234), so the coverage of the truth is defined
by the masks themselves, never by that box.

Masks: Imaris wrote a fourth channel, `Masked 405`, into every project, filled with one constant value
per surface (1 to 152 in n2rad51_noHS_herm_001; the lab's coloc-channel trick, by which the last
channel's per-object maximum identifies the surface a spot sits in). The truth export recovers each
nucleus mask from that channel and pairs it with its `SurfaceModelInfo` centroid and volume. Caveats
the export records per gonad: fill values can collide (two surfaces with the same value, about 3
percent of surfaces in the lab's xlsx exports), and a few small fragments carry values of their own
(n2rad51_noHS_herm_001 has 151 distinct values for 146 surfaces, four of them 5 to 364 voxels). The
scoring code ignores fragments below 1000 voxels and treats a collided value as one object only when
its connected components cannot be told apart. The voxel fill is a few percent larger than the mesh
volume in `SurfaceModelInfo` (voxelised fill against triangulated surface).

The deployed pipeline's masks are about 1.4 times the Imaris surface volume
(`docs/SPOTMAX_VALIDATION.md`), so IoU against this truth is depressed for every method by the
truth's tight threshold. The nucleus win condition is a relative comparison between methods against
the same truth, which that systematic difference leaves intact.

### 4.3 The coordinate frame

Imaris coordinates are voxel centres scaled by Imaris's rounded voxel size (0.108 um in xy, 0.2 um in
z); the `.nd2` metadata says 0.108333 x 0.108333 x 0.2 um. Over a 2600-pixel field the difference
grows to 0.87 um, three times the xy matching radius, so matching in raw micrometres fails
systematically on one side of the field. The rule: convert an Imaris coordinate to a voxel index with
Imaris's own voxel size, index = floor((coordinate - ExtMin) / imaris_voxel), and place it in the
pipeline frame as index times the `.nd2` voxel size. Pipeline spot coordinates are index times the
`.nd2` voxel size already (index = round(coordinate / voxel)). Everything in section 5 is computed in
voxel indices or in the pipeline frame derived from them. The lesson comes from the prototype named
in the changelog.

### 4.4 What the truth is not

Not a complete annotation of the field (surfaces), not a hand-verified count (spots), and not
independent of the per-gonad thresholds an annotator chose. The lab's canonical per-nucleus number
is the xlsx coloc method (`germquant.validate.imaris_xlsx`), which assigns spots to surfaces by the
mask-channel value. The traced-subset counts of section 5.1 reproduce that assignment geometrically
(spot voxel inside the recovered mask); the two are checked against each other once per gonad before
any scoring, and a disagreement is a truth-export bug to fix, not a result.

## 5. Metrics

Common definitions:

- Truth masks: the recovered Imaris nucleus masks of a gonad, fragments removed.
- Truth coverage region C: the union of the truth masks dilated by 1.0 um in xy and one slice in z.
  A predicted nucleus with less than half of its volume inside C is ignored by the instance metrics
  (it is probably an untraced nucleus, which the truth cannot judge).
- Spot region R: the union of the truth masks dilated by one voxel in every axis. Only truth spots
  and predicted spots whose centre voxel lies in R take part in spot matching and in the
  traced-subset counts.
- Spot to nucleus assignment: a spot belongs to the mask label at its centre voxel; if that voxel is
  background, to the most frequent non-zero label in the 3 x 3 x 3 voxel neighbourhood; otherwise it
  is unassigned. Unassigned truth spots are counted and reported per gonad.
- Two ways of producing a method's spots. End-to-end: the method's final spot list as its pipeline
  produces it, with its own nucleus masks if it uses any. Detector-only: the method's spot stage run
  on the raw image with the truth masks in place of its own nuclei. For mask-free detectors
  (Spotiflow, RS-FISH, Imaris) the two coincide.

### 5.1 Per-nucleus RAD-51 count agreement

Unit: nucleus. For each nucleus, truth count = truth spots assigned to it; method count = the
method's spots assigned to it. Statistics over the nuclei of a gonad set: Pearson r, Lin's
concordance correlation coefficient (CCC), bias (mean of method minus truth) and mean absolute
error, pooled over all nuclei of the set and also per gonad.

Two partitions:

- Traced subset (primary): the truth masks are the nuclei. Method count = the method's end-to-end
  spots assigned to the truth mask, with the detector-only count reported next to it. This partition
  is the same for every method and independent of any method's segmentation, so it is the one the
  win condition uses.
- All germline nuclei (secondary): the deployed pipeline's germline masks (`nuclei.csv` rows with
  `in_germline` true and `nuclei_labels.tif` of the 8666f07 goldens) are the nuclei for every spot
  method, 380 to 1240 per gonad; truth count = truth spots assigned to the mask, method count = the
  method's spots assigned to the mask. This shows detector agreement over the whole gonad, including
  the nuclei the annotator did not trace, on one fixed partition. A method with its own masks is
  additionally reported on its own masks, labelled as such; those numbers are not comparable across
  methods.

Class: truth counts under A+B (primary), A alone and B alone.

### 5.2 Spot matching

One-to-one matching between the method's spots and the truth spots inside R by the Hungarian
algorithm. A pair is admissible when the xy distance is at most 0.3 um and the z distance at most
0.4 um (an anisotropic cylinder, both distances in the pipeline frame); the cost of an admissible
pair is sqrt((dxy / 0.3)^2 + (dz / 0.4)^2) and inadmissible pairs cannot be matched. Precision =
matched / method spots in R; recall = matched / truth spots in R; F1 is their harmonic mean.
Reported end-to-end and detector-only, for Class A+B (primary), A and B, per gonad and pooled.

### 5.3 Nucleus instances (3D)

Predicted nuclei in scope: those with at least half of their volume inside C. Matching: one-to-one
by IoU (Hungarian, maximising total IoU) at threshold t; TP = matched pairs with IoU at least t, FN =
truth nuclei without a match, FP = in-scope predictions without a match.

- F1 at IoU 0.5 (primary), and recall at IoU 0.5 on its own (precision depends on the ignore rule,
  recall does not).
- AP = the mean over t in 0.50, 0.55, ..., 0.95 of TP / (TP + FP + FN).
- Splits: truth nuclei overlapped by two or more in-scope predictions each covering at least 20
  percent of the truth volume. Merges: in-scope predictions covering at least 20 percent of the
  volume of two or more truth nuclei.
- Z-extent error: for the pairs matched at IoU 0.5, the predicted nucleus's z span in slices minus
  the truth's, reported as mean signed and mean absolute error, in slices and in um (0.2 um per
  slice).
- Mean IoU of the matched pairs.

### 5.4 Reporting

Every table has one row per gonad and a pooled row, and every test gonad appears in every table,
with "no output" where a method produced nothing (section 10). Software versions, config hashes,
model weight sha256 and the split file sha256 (all in the manifest or the run manifests) are recorded
with each result, under `benchmarks/<name>/` next to the manifest.

## 6. Training and evaluation design

### 6.1 Roles

Train is for fitting and tuning; val is for model selection, thresholds and early stopping; test is
scored once, by the frozen configuration. After selection on val, a learned component is retrained on
train plus val with the selected settings, and that model is the one scored on test. Parameter-only
competitors are tuned on train plus val by the grid rule of 6.4.

### 6.2 Leave-one-gonad-out

Over the 8 train and val gonads: 8 folds, each trains or tunes on 7 and scores the held-out gonad
with the metrics of section 5. Selection inside a fold uses only its 7 gonads (the designated val
gonads when they are among the 7, otherwise one gonad of the same condition chosen by tag order). The
8 out-of-sample gonads are pooled with their own bootstrap interval. Leave-one-gonad-out is the
variance estimate and the guard against a lucky test split: a test result whose direction it
contradicts is reported as fragile.

For spot detectors the two region-growing gonads have no usable spot truth, so spot-side training
and leave-one-gonad-out run on 6 gonads (4 train, 2 val); the nucleus side uses all 8.

### 6.3 Confidence intervals

Cluster bootstrap by gonad: resample the gonads of the set with replacement (a gonad drawn twice
contributes its nuclei or spots twice), recompute the pooled statistic, 1000 resamples, seed
20251105, percentiles 2.5 and 97.5 as the 95 percent interval. A paired comparison (method A against
method B) computes both statistics on the same resample and takes the difference; the interval of
the difference decides the win conditions. With 4 test gonads an interval of a difference that
excludes zero requires the gain to hold on essentially every gonad; that is deliberately
conservative.

### 6.4 Tuning budget

At most 60 configurations per competitor, each listed in the result record with its val score; the
configuration with the best val score on the primary metric is frozen before test is scored.

## 7. Win conditions (pre-stated)

Evaluated on the test set with the primary settings (traced subset, Class A+B, end-to-end) against
the SpotMAX reference, which is the better of SpotMAX deployed and SpotMAX re-tuned on val (chosen on
val, not on test; both are reported on test).

- W1 (counts): CCC(new) minus CCC(SpotMAX reference) is at least 0.03 and the paired bootstrap
  interval of the difference excludes zero.
- W2 (spots): F1(new) is not lower than F1(SpotMAX reference) by more than 0.02 (point estimates; the
  interval is reported).
- W3 (nuclei, only when the new pipeline changes the nucleus segmentation): F1 at IoU 0.5 of the new
  segmentation minus that of the deployed Cellpose model `germline_nuclei_combined` is at least 0.03
  and the paired interval excludes zero.

Decision: W1 and W2 together decide whether the new detector replaces SpotMAX; W3 decides the
nucleus model on its own. W1 without W2 (better counts, worse localisation) is a mixed result and
does not replace the detector; the write-up says so.

## 8. Competitors

1. SpotMAX deployed: the pipeline at commit 8666f07 (section 2). Its outputs are the goldens;
   nothing is re-run.
2. SpotMAX re-tuned: the same detector and the deployed masks, grid over thresholding method
   (`threshold_otsu`, `threshold_li`, `threshold_triangle`, `threshold_yen`), `effect_size_min`
   (2.0, 2.5, 3.0, 3.5, 4.0) and `z_merge_gap_um` (0.6, 0.8, 1.0): 60 configurations, chosen by the
   val rule of 6.1 and 6.4.
3. Spotiflow zero-shot: the published general model with its default probability threshold on the
   RAD-51 channel (3D, or slice-wise with z merging if the 3D model does not apply; which one is
   decided and written in the changelog before test is scored). Spotiflow fine-tuned: from the same
   weights, on the 4 train gonads with usable spot truth (Class A+B spots over the whole field as
   point labels), val for early stopping.
4. RS-FISH: difference of Gaussians with RANSAC localisation; sigma and intensity threshold over a
   grid of at most 60 on train plus val; anisotropy 0.2 / 0.108333.
5. Imaris with one global Quality threshold: the same Spots parameters re-run in Imaris on all 12
   gonads with a single Quality threshold for all, chosen on train plus val by the primary metric.
   It needs a re-export from Imaris, because the existing projects store only the applied threshold.
   This competitor measures what Imaris gives without per-gonad hand tuning.
6. The human blind two-scorer subset of section 9.

## 9. Human blind subset and its size

Two scorers, blind to each other and to every method's output, count RAD-51 foci in a fixed sample
of truth nuclei, each shown as a cropped z-stack of the RAD-51 channel with the nucleus mask
outline. Sample: 160 nuclei, 40 per test gonad, drawn with seed 20251105 from the traced subset,
stratified by truth-count tertile within the gonad so that sparse and dense nuclei are both covered.

Size: the quantity estimated is a concordance of about 0.8 between two counts. With the Fisher z
approximation the standard error of a correlation is 1 / sqrt(n - 3); a 95 percent interval
half-width of 0.05 around r = 0.8 corresponds to a z-interval of 0.158, so n = (1.96 / 0.158)^2 + 3
= 157 nuclei, rounded up to 160.

The scorers' counts give (a) the inter-scorer CCC, the ceiling any method is compared with, (b) each
scorer's CCC against the Imaris truth, the truth's own agreement with a person, and (c) each method's
CCC against the mean of the two scorers on these nuclei. Disagreements between the scorers are data
and are not resolved.

## 10. Stop rule

- Development uses train and val only. Before test is scored, the method's code commit, config hash
  and weight sha256 are written into the changelog; a method still being changed is not scored on
  test.
- The test set is scored at most twice per method family: once with the pre-registered method, and
  at most once more after a revision declared in the changelog, with its reason, before the second
  scoring. Both results are reported wherever either is. A third attempt is not allowed on this test
  set; it needs a new test set from a new imaging session, registered the same way.
- A method that produces no output for a test gonad (as the deployed pipeline does for
  n2rad51_noHS_male_002 today) gets a row saying so in every table; the gonad is not dropped and
  nothing is substituted. The win conditions are evaluated over the test gonads where both methods
  produced output, labelled with that number (3 of 4, say). A win declared on fewer than all 4 test
  gonads is provisional until the missing gonad is processed or its failure is explained in the
  changelog.
- Negative result: if W1 fails after the allowed attempts, the SpotMAX reference stays deployed and
  the negative result is recorded in the changelog and in the write-up.
- Tuning stops at the budget of 6.4, or earlier when the val primary metric has not improved by
  more than 0.01 over the last 10 configurations.

## 11. Changelog (append only)

Entries are appended, never edited. Each one says whether a test result had been computed at that
point.

1. 2026-10-09. Protocol written and split registered (`tests/golden/rad51_split.json`, 12 gonads,
   sha256 of all 24 files). No result has been computed under this protocol. Pipeline outputs
   available at registration (8666f07 goldens, as the manifest built that night records):
   n2rad51_noHS_herm_000, _001, _002, n2rad51_noHS_male_000, _001, n2rad51_HS_herm_000, _001, _002
   complete (8 of 12); n2rad51_noHS_male_002 failed (detector flood, no outputs until a performance
   fix); n2rad51_HS_male_000 running; n2rad51_HS_male_001 and _002 not yet run. The Imaris truth
   export (the `<stem>__imaris_*` files) is being written in parallel and no truth file existed yet;
   the file names in `benchmark_split.TRUTH_FILES` are the expected contract and will be confirmed
   here when it lands. Prior exposure, disclosed: before this
   protocol, a prototype outside the repository (`C:/Users/ryane/germquant_golden/rad51_ims_benchmark.py`)
   computed exploratory agreement numbers for the deployed pipeline against the Imaris spots
   (per-nucleus r and CCC on the deployed masks, greedy spot matching at 0.5 um) on the June
   `results_cv` outputs of all 12 gonads and on the 8666f07 goldens of the four no-heat-shock gonads
   finished at the time (noHS herm 000, 001, 002 and noHS male 000). That includes test gonads. Those
   numbers established the coordinate-frame rule of 4.3 and showed the deployed pipeline in the
   per-gonad CCC range 0.5 to 0.8; the split was fixed by the cell structure and the region-growing
   defect, and the margins of section 7 were set by the project lead before this document, not
   derived from those numbers. No competitor and no new method has been scored.
2. 2026-10-10. Development tuning of the classical detectors on train and val (package D3;
   `docs/BENCHMARK.md`, "Tuning on train and val"). No test gonad was run or scored for it and no
   per-gonad test number was looked at. Disclosed: the pooled agreement of the deployed recipe over the
   nine gonads that had outputs and usable spot truth, three of them test gonads (per-nucleus CCC 0.84,
   bias +0.77 foci per nucleus, spot F1 0.39), had been computed before this entry with the benchmark
   harness and was known when the searches were designed; by rule 1 that computation was the first
   result, and this entry is the first to record it. What the searches did differently from sections 6
   and 8 (each of (a) to (d) is post hoc under rule 3: it was made after the first result): (a) the SpotMAX re-tuned competitor was searched over `thresholding_method`
   (`threshold_triangle`, `threshold_otsu`, `threshold_li`; `threshold_yen` not searched),
   `effect_size_min`, `spot_radius_um` and `gauss_sigma_um` (`z_merge_gap_um` not searched), 24
   configurations in two stages, instead of the 60-configuration grid of section 8 item 2, to fit
   SpotMAX's compute cost; (b), for the same compute cost and because the package was designed as a
   search on train with finalists sent to val, both searches (SpotMAX, and the experimental `log` backend with 55
   configurations in three stages) ran on the four train gonads with usable spot truth, and only their
   best configurations (3 for SpotMAX, 5 for `log`) were scored on val, where section 6.4 lists a val
   score for every configuration and section 6.1 tunes parameter-only competitors on train plus val;
   the early stop of section 10 was not used; (c) the selection objective was the pooled per-nucleus
   CCC over the germline nuclei of the deployed masks (end-to-end, Class A+B; the harness's
   `per_nucleus_all`), with the pooled spot F1 as tie-break, not the traced-subset partition of section
   5.1, because the germline partition has about six times as many nuclei as the traced subset and the
   harness does not compute the section 5.1 primary partition; (d), to make the pre-gate searchable,
   the `log` backend's candidate pre-gate became a config key, `spots.log.pre_gate_sigmas`
   (code default 3.0, the former constant). Which SpotMAX configuration is the re-tuned competitor of
   section 7, this search's or one from the section 8 grid, is to be decided and written here before
   test is scored; on the harness's traced-subset CCC, the closest proxy it computes for the primary
   metric, the val rule would pick smS2_tri_r0.35_g0.08_e4 (0.756) over the selected
   smS2_tri_r0.3_g0.12_e3.5 (0.717), so the reference has to be chosen under the primary metric.
3. 2026-10-10. net3d, the learned 3D U-Net detector (package D2; `docs/BENCHMARK.md`, "net3d"). No
   test gonad was run or scored for it. It was trained on the four train gonads with usable spot truth
   and selected on the two val gonads. Its selection rule is declared here as a deviation from 6.4,
   post hoc under rule 3 and adopted on val after launches 1 to 3: every (epoch, threshold) scored on
   val (216 for the frozen model, 12 epochs times 18 thresholds; 288 for the no-DAPI ablation; more in
   the interrupted launches) whose pooled traced-subset CCC is within 0.01 of the best is eligible, and
   the eligible one with the best pooled spot F1 is frozen. Each (epoch, threshold) pair counts as one
   configuration against the 6.4 budget, so net3d exceeded it. The literal 6.4 rule (best val score
   on the primary metric) would have frozen epoch 5 at threshold 0.05 (val germline CCC 0.959, F1
   0.608) instead of epoch 4 at 0.40 (germline CCC 0.982, F1 0.630). Frozen model: sha256
   8cdb710bd2f0fe7acef49c6117d7113ad61e224405960d2847b34c8327a4c4f7, trained at commit 87d77f1 with
   seed 0 (`C:/Users/ryane/germquant_golden/models/net3d_20261010/`). Section 6.1 asks for the selected
   configuration to be retrained on train plus val before test is scored; that has not been done.
4. 2026-10-10. Package E1: the section 5.1 partition in the harness, leave-one-gonad-out estimates and net3d
   retrained on train plus val (`docs/BENCHMARK.md`, "Leave-one-gonad-out"). No test gonad was run, scored or
   looked at, and no test result was computed: no `--final` run was made, `final_runs.jsonl` does not exist,
   and `docs/FINAL_SCORING.md`, the plan for the one-time test scoring, was not executed. (a) The primary
   partition of section 5.1 is implemented as written there (`germquant.benchmark.spots.truth_mask_table`; `tm_`
   columns, `per_truth_nucleus.csv` and the `per_nucleus_truth_masks` blocks of `bench_detector.py`,
   `bench_tune.py` and `benchmark_spots.py`): the truth masks with fragments below 1000 voxels removed are the
   nuclei (one 44-voxel label of n2rad51_noHS_male_000 is removed among the train and val gonads; the truth export
   had already split collided values), every truth spot and every method spot is assigned by the section 5 rule
   (the spot region R read as the masks dilated by a 3 x 3 x 3 cube, which equals the set of assigned centre
   voxels), and Pearson r, Lin's CCC, bias and mean absolute error are reported for Class A+B, A and B truth, per
   gonad and pooled; `bench_detector.py --detector-only` re-runs the spot stage with the truth masks as its nuclei
   for the detector-only count. The harness's earlier "traced" subset is kept and is not this partition. The
   spot F1 of section 5.2 (matching restricted to R) is the `spots_truth_masks` block and the `tm_spots_` columns;
   the whole-field F1 is kept for comparison. Truth check: the section 4.4 check, run with `scripts/check_truth_assignment.py` against the lab's mask-channel
   assignment (`imaris_mask_label`, Imaris's per-spot maximum of the Masked channel, the value the xlsx coloc counts
   come from; no xlsx export of these gonads is registered), agrees on 99.35 to 99.99 percent of the spots of all
   eight train and val gonads. The disagreements are boundary spots (0 to 14 geometric-only and 1 to 22
   coloc-only per gonad) and spots whose coloc value is a removed mask fragment (0 to 58); only HS_male_002, a
   region-growing gonad, has spots assigned to two different nuclei (3). No truth-export bug was found. Table:
   `benchmarks/n2rad51_8666f07/truth_assignment_check_trainval.csv`. The test gonads are checked in the pre-flight
   of docs/FINAL_SCORING.md.
   (b) Leave-one-gonad-out over the six train and val gonads with usable spot truth: six folds (for a held-out
   val gonad the first same-condition train gonad by tag order joins the remaining val gonad), each held-out
   gonad scored by `bench_detector.py --detector-only` with its fold in the trial record's `note` (trial lines 162
   to 173 and 214 to 231). net3d: each fold retrained on its five gonads by the frozen recipe of net3d_20261010
   (its settings, seed and 40-epoch schedule, stopped after the selected 4 epochs, threshold 0.40, no
   validation); retraining the recipe on its own four train gonads reproduced the selected weights bit for bit.
   This departs from section 6.2, which selects inside each fold: the recipe's epoch and threshold were chosen on
   the two designated val gonads, so the two folds that hold one of them out are optimistic for net3d, and its
   estimate leaves out the variance of its selection. Tuned log and re-tuned SpotMAX: the whole D3 search replayed
   on each fold's own gonads (its stages, budgets of 55 and 24 sets, finalists and objective, which is the
   germline partition, not the primary metric); first all 55 log and 24 SpotMAX D3 sets were scored on the two val
   gonads (lines 97 to 151 and 174 to 197); in the two folds whose train gonads are the three no-heat-shock train
   gonads the replay left the D3 path, and the sets it called for were run (lines 152 to 161 and 198 to 213). The
   deployed SpotMAX recipe and the log defaults: their per-gonad scores, re-run with the same harness. Pooled over
   the six held-out gonads, primary partition (Class A+B, end-to-end), CCC [95 percent cluster-bootstrap
   interval]: SpotMAX deployed 0.601 [0.499, 0.718], SpotMAX re-tuned 0.742 [0.655, 0.812], log defaults 0.777
   [0.623, 0.920], log tuned 0.893 [0.818, 0.955], net3d 0.878 [0.834, 0.935]. Paired differences whose interval
   excludes zero: net3d and the tuned log better than both SpotMAX configurations on every partition, on spot F1
   with and without Class B and in mean absolute error (net3d minus SpotMAX re-tuned +0.136 [+0.090, +0.204] in
   primary CCC); the log defaults better than the deployed recipe on every number and than the re-tuned SpotMAX
   on F1; the re-tuned SpotMAX better than the deployed recipe in mean absolute error only. No difference survives
   between net3d and the tuned log (-0.015 [-0.084, +0.073]); on the four folds that hold out a train gonad the
   tuned log leads net3d in germline CCC and F1. These are development estimates; nothing here is a win
   condition. (c) net3d retrained on the six train and val gonads by the same frozen recipe (protocol section
   6.1):
   `C:/Users/ryane/germquant_golden/models/net3d_trainval_20261010/model.pt`, sha256
   d091b9ef58850783387a4f5aa0fea5ecc9310abc027145771936fecd0de001c0, trained at commit 240cc72, with its
   `model_card.json`; it is the net3d of `docs/FINAL_SCORING.md` and was not scored. Every choice in (b) and (c) was
   made after the first result of entry 2 and is post hoc under rule 3. No test result was computed.

5. 2026-10-10. nuclei3d, the in-house 3D nucleus model (package E2; `docs/BENCHMARK.md`, "nuclei3d").
   No test gonad was run or scored for it, and no metric of any method was computed on a test gonad in
   this package. (a) The scope rule of section 5.3 is now implemented in the nucleus harness as
   `coverage="protocol"` (`germquant.benchmark.nuclei3d`; `scripts/benchmark_nuclei3d.py --coverage
   protocol`): truth fragments under 1000 voxels removed, region C, a prediction in scope with at least
   half of its volume in C, splits and merges at 20 percent of the truth volume. Until this entry the
   harness implemented only its earlier rule (the ROI box when present, else the truth bounding box
   dilated by 2 um; majority splits and merges), which section 4.2 rules out as a record of the traced
   extent; that rule stays the script's default so earlier outputs reproduce, and every nucleus number of
   this entry uses the protocol rule. The script now skips test gonads unless `--final` is given with a
   `--reason`. (b) What was trained: a residual 3D U-Net with foreground, cell-distance and neighbour
   heads, on the six train gonads (582 traced nuclei; the region growing of the two heat-shock males
   affects their spots, not their surfaces) against the recovered Imaris masks, with zero loss weight
   inside untraced deployed Cellpose labels plus a 0.5 um margin; selected on the two val gonads by the
   pooled section 5.3 F1 at IoU 0.5, then AP, then the earlier epoch, with fixed instance-step parameters
   during training and then a grid of 24 instance-step configurations for the selected epoch, which kept
   the training's parameters. 40 configurations were scored on val (16 epochs and the 24 grid points),
   within the 60 of section 6.4, each listed with its val score in the training log. Frozen model: sha256
   11920205258dceece3e2bd6053621cbce390177f2baa70b9f40eefc0980bbb83 (epoch 8), trained at commit 02d8389
   with seed 0 (`C:/Users/ryane/germquant_golden/models/nuclei3d_20261010/`). For comparison the
   deployed Cellpose labels (the 8666f07 goldens) were scored on the same two val gonads under the same
   rule. (c) Choices that neither this protocol nor the design fixed, made before training: the region of
   a traced nucleus's own Cellpose label outside its traced mask has weight 1 and target background (the
   design's weight rule leaves that region unnamed); a Cellpose label with at least half of its voxels
   inside traced masks is a traced nucleus's own label; the neighbour closeness reaches 0 at 1 um (power
   3); the validation runs the instance step only on foreground components that touch region C, which
   changes no metric. (d) Not done: retraining on train plus val (6.1), leave-one-gonad-out (6.2), the
   RAD-51 counts on the new masks; W3 is a test-set criterion and has not been evaluated. (e) The
   design's decisive competitor was run once, without tuning: the deployed Cellpose model fine-tuned on 144
   slices of the same six train gonads' traced masks (`scripts/export_imaris_training_pairs.py`, untraced
   nuclei filled with background; `germquant train-nuclei`, 100 epochs; sha256
   a4334d0016f3ce548c92ddb722011be908c922d576e5d740c1405c714561c38c,
   `C:/Users/ryane/germquant_golden/models/cellpose_ft_imaris_20261010/`), one configuration scored on the
   two val gonads by the same rule. On val it ties with nuclei3d (pooled F1 at IoU 0.5 0.974 against
   0.968, AP 0.590 against 0.582; deployed Cellpose 0.863 and 0.338), a difference far below the 0.03 of
   section 7, so the design's stop rule is not triggered; `docs/BENCHMARK.md` reports the whole-stack
   diagnostic in which nuclei3d makes far fewer, larger objects than both Cellpose models (527 against
   1270 deployed and 1234 fine-tuned on the heat-shock gonad), merging crowded untraced nuclei that the
   traced-subset metrics cannot see.
   (f) Disclosed: two output folders of the nucleus harness, `benchmarks/n2rad51_8666f07/nuclei3d` and
   `nuclei3d_bbox`, were written on 2026-10-10 at 02:47 and 02:49 UTC (20:47 and 20:49 local time on
   2026-10-09), three minutes after the manifest, under the earlier rule. They hold test gonads: rows for
   n2rad51_HS_herm_002, n2rad51_HS_male_001 and n2rad51_noHS_herm_002 (deployed Cellpose-SAM nuclei against
   the recovered masks). They were written by the coordinating session, which ran benchmark_nuclei3d.py and
   benchmark_spots.py over every gonad with pipeline outputs that evening as the first look at the deployed
   pipeline; their per-gonad console lines, test gonads included, were displayed in that session and some
   were reported to the project lead (for example the heat-shock male test gonad's 15.7 against 12.9 foci per
   nucleus, and per-gonad nucleus F1 lines). By rule 1 these are the first nucleus results on test, and the
   spot lines extend what entry 2 discloses. They are results of the deployed pipeline only; no method was
   tuned or selected on them, and nobody in packages D1 to E3 opened their test rows. Each of the folders
   benchmarks/n2rad51_8666f07/spots, nuclei3d and nuclei3d_bbox now carries a README saying it holds test rows.

6. 2026-10-11. The project lead's reply to the five decisions of `docs/FINAL_SCORING.md` (section 1), and the
   SpotMAX grid of section 8 item 2, declared before the grid was run (nine of its sets, in effect, had been run
   earlier under D3 names; see (b)). No test gonad was run, scored or looked at for this entry, no `--final` run
   was made and `final_runs.jsonl` does not exist. The test numbers seen before it are the deployed pipeline's
   (entries 1, 2 and 5 (f)) and those of the pre-protocol sweep named in (b).
   (a) A message to the project lead (Ryan Eastman) asked him to merge a pull request and to settle the five
   decisions, and listed his open tasks, the re-export of the two heat-shock males among them. He replied on
   2026-10-10 at 23:18 MDT (05:18 UTC on 2026-10-11), verbatim: "Yes to everything and merged". Decisions 3 (both
   SpotMAX versions use the family's two test scorings) and 4 (the missing competitors of section 8 are scored
   later) were put to him as "accept that" and are accepted as written. Decision 5 (the post hoc choices stand) is
   accepted with a reading: in `docs/FINAL_SCORING.md` its list names "the SpotMAX grid not run" and it confirms
   the methods of section 2 as they stand, the D3 row of SpotMAX re-tuned among them, and both are left to
   decision 2. Decisions 1 and 2 are choices, which the reply does not settle by itself; what follows for them is a
   reading, not his words. Decision 1 (wait for the re-export, or score now) is read as waiting, which keeps both
   options open: nothing is scored on test until he confirms one explicitly. Decision 2 (the SpotMAX reference; the
   best D3 set on val starves the heat-shock gonads): the protocol's own grid is run now, on train and val only,
   as (b) to (e) declare it, and the recommendation is that its section 6.4 selection becomes SpotMAX re-tuned.
   That settles which rule supplies the set, not the heat-shock concern: the selection is most likely the otsu set
   of decision 2 (see (b)). He confirms the recommendation and the readings of (c) and (d) before the pre-scoring
   entry; if he picks another set, that entry records it and its reason as post hoc. The grid's ranking and the D3
   candidates are reported either way, the D3 candidates as post hoc (entry 2).
   (b) The grid: `thresholding_method` (`threshold_otsu`, `threshold_li`, `threshold_triangle`, `threshold_yen`) by
   `effect_size_min` (2.0, 2.5, 3.0, 3.5, 4.0) by `z_merge_gap_um` (0.6, 0.8, 1.0), 60 sets, each run on the
   deployed nucleus masks. A set's parameters are exactly these three keys over the golden `spots:` section, so
   every other key keeps its deployed value (`spot_radius_um` 0.3, `gauss_sigma_um` 0.08, `merge_z_columns` true,
   `z_merge_valley_frac` 0.8, `effect_size_metric` `spot_vs_backgr_effect_size_glass`, `max_spot_candidates`
   30000). Listing order: method outermost as listed, then `effect_size_min`, then `z_merge_gap_um`, both
   ascending. Names: `pg_<otsu|li|tri|yen>_e<value>_z<value>`, one decimal each; the 38th, `pg_tri_e3.0_z0.8`, is
   in effect the deployed recipe. Known before this entry: (1) nine sets (tri, otsu and li at e2.0, e3.0 and e4.0,
   z0.8) are in effect the D3 stage-S1 sets `smS1_<tri|otsu|li>_r0.3_g0.08_e<2|3|4>`, with train scores from D3
   (germline partition, among trial lines 54 to 64) and val scores from E1 (section 5.1 partition, among lines 174
   to 184); entry 4 (b) adds section 5.1 scores for the deployed recipe on all six train and val gonads and for the
   otsu set on both heat-shock hermaphrodites (held-out runs, lines 214, 215 and 220 to 225; truth-mask bias -2.33
   on HS_herm_000). Among the nine are the deployed recipe (val truth-mask CCC 0.556) and `pg_otsu_e4.0_z0.8`
   (0.822, the best of the nine), the otsu set of decision 2. (2) Before the protocol, a sweep recorded in
   `docs/SPOTMAX_VALIDATION.md` (section 4, June 2026) ran otsu, li and triangle by `effect_size_min` on the six
   gonads without heat shock, two of them now test gonads (noHS_herm_002 and noHS_male_002), and gives per-gonad
   results on them; the deployed recipe was chosen from it. Entry 1 does not name it. Decision 2 is not taken on
   that record (`docs/FINAL_SCORING.md`, section 5). If the nine sets reproduce, the selection is
   `pg_otsu_e4.0_z0.8` unless one of the 51 sets not yet scored ranks above it. The rule is applied as written
   whatever it does to the heat-shock gonads, and the selection's train scores are reported with it.
   (c) Scored: all 60 sets on the two val gonads and on the four train gonads with usable spot truth;
   `bench_tune.py` skips the two region-growing heat-shock males (section 4.1) and refuses a test gonad. Section
   6.1 tunes parameter-only competitors "on train plus val", while 6.4 and section 8 item 2 select by the val
   score. This entry selects on the two val gonads, as 6.4 is written, and reports what `select_grid` picks on the
   six train and val gonads pooled; this reading was made knowing (b). Train scores serve that report and the folds
   of (e) and enter no selection on the designated split, so a re-export under decision 1 would add two folds
   to (e) and leave that selection as it is. The effect-size shortcut of D3 and E1 is used (`effect_size_min` is
   the last step of `detect_spots_spotmax`, a row filter, so one SpotMAX run per method and gap, with the minimum
   at 0, serves all five minima; `tuning.filter_effect_size`, 12 runs per gonad). On real data, runs without it
   gave the same results for the deployed recipe (the golden spot tables, all six train and val gonads), the otsu
   set at 4.0 and two D3 triangle sets (`docs/BENCHMARK.md`); li, yen and the gaps 0.6 and 1.0 rest on the order of
   steps alone. Checks, where a difference is a bug that stops the selection until its cause is found and recorded
   in this changelog: `pg_tri_e3.0_z0.8` gives the golden spot tables (`detector.compare_with_golden`), the nine
   sets give their earlier counts, and each held-out run of (e), made without the shortcut, gives its grid score. A
   run that fails is repeated with fewer workers and nothing else changed; a set that fails again is recorded in
   this changelog and ranks last. A run that reaches the candidate cap is scored as it is, and the cap warning is
   reported.
   (d) Selection (6.1 and 6.4), by `select_grid` on the val gonads: the pooled val truth-mask CCC (section 5.1
   partition of entry 4 (a), Class A+B truth, end-to-end counts, pooled over the nuclei of both gonads as
   `bench_detector.py` pools it; summary column `tm_ccc`) rounded to 3 decimals, higher first; ties by the pooled
   val spot F1 inside R (section 5.2, block `spots_truth_masks`), higher first; then listing order. Rounding and
   the tie-breaks go beyond the letter of 6.4 (the best val score): if the unrounded CCC (that column) picks
   another set, both are reported, as in entry 3. An undefined value ranks last, as does a set that failed
   twice (c); otherwise nothing is selected until every set has a score on both val gonads. The first set is the
   grid's selection. If the project lead confirms the recommendation of (a), it is frozen as SpotMAX re-tuned, and
   it is then also the reference of section 7, since the deployed recipe is one of the 60. The ranking and the
   `bench_tune.py` summaries (every set's val and train scores) are the record of 6.4. Budget (6.4): 60 distinct
   sets under the search label `prereg_spotmax_grid` (`--budget 60`; the train run repeats them). The 24 D3 and 12
   E1 SpotMAX sets (searches `d3_spotmax`, `logo_spotmax_train3noHS`) are exploratory and post hoc: reported, and
   not candidates if the grid is confirmed. With them 96 named SpotMAX sets will have been scored, 87 of them
   distinct in effect and 77 of those on val; counted per competitor that is over the 60 of 6.4, as it was for
   net3d (entry 3), so counting the grid alone is a reading. Early stop (section 10): read as governing sequential
   searches, where a set is chosen or run after earlier ones are scored; a fixed grid of exactly the budget is run
   in full. Where the rule would have stopped in listing order (the first set, from the 11th on, at which the best
   rounded val CCC so far is at most 0.01 above the best before the last ten) is reported with what it would have
   selected, computed from the val scores that `select` writes. With the known otsu scores that point most likely
   falls among the li sets, before any triangle or yen set; this reading was made knowing that.
   (e) Leave-one-gonad-out (6.2, the six folds of entry 4 (b)): `replay_grid` applies the rule of (d) to each
   fold's selection gonads, the designated val gonads among its five and, for a held-out val gonad, the first
   same-condition train gonad by tag order (HS_herm_000 for HS_herm_001, noHS_herm_000 for noHS_herm_001), so the
   folds that hold out a train gonad select as the designated split does. Unlike the D3 replay, every fold selects
   by the primary metric. Each held-out gonad is scored afresh (`bench_detector.py --detector-only`) and pooled
   with the five detectors of entry 4 (b), every pair on the same resamples (6.3). The section 5.2 F1 inside R is
   not among the blocks `bench_logo.py report` pools at e19a029 (`logo.DETECTION_BLOCKS`), and the held-out runs of
   entry 4 have no columns for it, so pairs on it need that block added and those runs scored again.
   (f) Code, added with this entry. `germquant.benchmark.logo`: `PREREG_SPOTMAX_METHODS`,
   `PREREG_SPOTMAX_EFFECT_SIZES`, `PREREG_SPOTMAX_Z_GAPS` and `prereg_spotmax_grid()` for (b);
   `primary_objective_on(idx, sha, tags)`, on a run index built by `index_runs(..., truth_masks=True)`, and
   `select_grid(idx, configs, tags)` for (d); `replay_grid(idx, configs, train, val)` for (e), with the returns of
   `replay_spotmax_search`. `scripts/bench_logo.py`: `grid-configs --out FILE`, and `select --search prereg-grid`,
   which writes the six fold selections and the designated one in the D3 replay's layout, exit code 3 included (the
   D3 replay stays the default). Outputs, under G =
   `C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/prereg_grid`: `G/configs.json`, the `bench_tune.py`
   runs on val and then train (`G/val`, `G/train`), `G/selection.json`, the held-out runs
   (`G/logo/heldout/spotmax_grid`) and the report with the five detectors of entry 4 (`G/logo/report`). The runs
   start after this entry is committed, so every grid trial record carries the sha of a commit that holds it. The
   grid and its val rule are the pre-registered ones (sections 8, 6.1 and 6.4), so running them is not a change,
   and the folds of (e) are those of entry 4 (b). What this entry adds (the val-only reading, the rounding and
   tie-breaks, the failure rule, the budget count, the early-stop reading) comes after the first result of entry 2
   and is post hoc under rule 3, as in entries 2 and 4; it was written knowing the scores of the nine sets of (b).
   No test result was computed for this entry.
   (g) Disclosed: the real-golden regression checks made before each pull request reran the deployed pipeline
   on test gonads. For tranche 1 (candidate 75448d8) n2rad51_HS_herm_002, n2rad51_HS_male_001 and
   n2rad51_noHS_herm_002 finished (n2rad51_noHS_male_002 did not); for tranche 2, n2rad51_HS_male_001; for
   tranches 3 and 4, n2rad51_HS_male_001 and n2rad51_noHS_male_002. Each rerun was compared only with its
   pinned 8666f07 outputs (every table cell for cell, every label image, the config hash). No truth file was
   opened and no metric against truth was computed. The check logs print the pipeline's own nucleus and spot
   counts for those gonads, numbers the 8666f07 goldens of entry 1 already hold. These reruns test the code of
   the deployed pipeline, not a method under development, and nothing was tuned or selected on them; no
   earlier entry recorded them.

7. 2026-10-11. The SpotMAX grid of entry 6, run as that entry declares it, and the section 5.2 F1 inside R added to
   the leave-one-gonad-out pairs (`docs/BENCHMARK.md`, "Results of the pre-registered SpotMAX grid", which gives
   every number). No test gonad was run, scored or opened for this entry: no `--final` run was made,
   `final_runs.jsonl` does not exist, and no record of `trials.jsonl` names a test gonad. No test result was
   computed.
   (a) The run. Commit abc49ca, which holds entry 6 (committed at 00:59:15 MDT on 2026-10-11), on a clean tree,
   with the planned commands of `docs/BENCHMARK.md` as written (a driver ran them in order with
   `CUDA_VISIBLE_DEVICES=-1`), from 00:59:51 to 02:18:22 MDT; every step exited with code 0 (`G/run.log`, with G as
   in entry 6 (f)). Trial lines 232 to 351 hold the grid (60 val and 60 train records, 60 distinct sets under
   `prereg_spotmax_grid`, the effect-size shortcut on and complete in every record, SpotMAX 1.3.1) and lines 352 to
   357 the six held-out runs, all at abc49ca without `-dirty`. No run failed, so the failure rule of entry 6 (c)
   was not used. The counts of entry 6 (d) are as it gave them: 60 distinct grid sets, and with the D3 and E1 sets
   96 named SpotMAX sets, 87 distinct in effect, 77 of those on val. Counted per competitor, that is over the 60 of
   section 6.4; counting the grid alone is the reading of entry 6 (d) that the project lead is to confirm. The grid
   sets equal in effect to earlier sets are exactly the nine of entry 6 (b). `bench_tune.py` ranks the sets by its
   own germline rule in its console summary and in the rank columns of its summary.csv (val: `pg_otsu_e3.0_z0.8`
   first, the selection 10th), and names that objective in its trial records; that order selects nothing, and the
   selection is `select_grid`'s (`G/selection.json`).
   (b) The checks of entry 6 (c) all passed. `pg_tri_e3.0_z0.8` leaves the golden `spots:` section unchanged and
   gives the golden spot tables on all six gonads (`detector.compare_with_golden`, parquet and CSV: every compared
   column, and no per-nucleus count different). The nine sets of entry 6 (b) give their earlier tables, on val
   those of E1 (lines 174 to 184) and on train those of D3 stage S1 (lines 54 to 64), in every common column but
   `detect_s` and `score_s`. Each held-out run of entry 6 (e), made without the shortcut, equals the grid run of
   its set on its gonad in every table, apart from the run times and the detector-only columns that the grid does
   not compute. One of them (`pg_tri_e4.0_z1.0` on n2rad51_HS_herm_001) is the first such check of a 1.0 um gap;
   li, yen and the 0.6 um gap still rest on the order of steps alone. The E1 held-out runs of the deployed recipe
   and of the otsu set (lines 214, 215 and 220 to 225), also made without the shortcut, equal the grid's tables
   too. Three train runs reached the candidate cap (`threshold_li` on n2rad51_noHS_herm_000, one per z gap; 32720
   candidates, 2720 dropped). The warning names no set; the order of the run log places it. They were scored as
   they are, so the 15 li train counts on that gonad are floors. No val run and no held-out run reached the cap.
   (c) The selection of entry 6 (d) is `pg_otsu_e4.0_z0.8` (`threshold_otsu`, `effect_size_min` 4.0,
   `z_merge_gap_um` 0.8): val truth-mask CCC 0.822 (0.822104 unrounded), val F1 inside R 0.612. `pg_otsu_e4.0_z1.0`
   ties it on both keys with byte-identical detections on all six gonads, and listing order decides. Every set has
   a score on both val gonads, and no value is undefined. The unrounded CCC picks the same set; it reorders ranks 3
   to 5 and three triangle pairs tied at three decimals. The deployed recipe ranks 23rd of 60 (0.556). Reported
   only: `select_grid` on the six train and val gonads pooled picks the same set (CCC 0.747), and the early stop of
   section 10, read as in entry 6 (d), would have stopped at the 20th set, `pg_li_e2.5_z0.8` (best so far 0.822
   against 0.815 before the last ten), with the same selection. The selection's train scores, as entry 6 (b) asks:
   pooled over the four train gonads CCC 0.711 [0.621, 0.791] (4th of 60 there; `pg_tri_e4.0_z0.8` leads with
   0.736) and bias +0.40 foci per truth nucleus. That pooled bias combines -2.33 on n2rad51_HS_herm_000 (CCC 0.587;
   the deployed recipe +0.91, CCC 0.843) with +0.72 to +1.38 on the three train gonads without heat shock. On the
   val gonads, in sample: -0.77 on n2rad51_HS_herm_001 (the deployed recipe +3.79) and +0.30 on
   n2rad51_noHS_herm_001. Put plainly, the selection undercounts both heat-shock hermaphrodites (1.73 foci per
   truth nucleus where the truth has 4.06, and 0.94 against 1.71) and overcounts the four gonads without heat
   shock. The deployed recipe overcounts all six; on n2rad51_HS_herm_000 it is better in CCC (0.843 against 0.587)
   and in F1 inside R (0.635 against 0.401). If the selection is confirmed as SpotMAX re-tuned, it is the reference
   of section 7 (val CCC 0.822 against 0.556 for the deployed recipe, in sample). The D3 candidates, post hoc
   (entry 2), on the same val score: the D3 selection 0.695 and entry 2's alternative 0.723; the best D3 set is the
   selection in effect.
   (d) Leave-one-gonad-out, as entry 6 (e) declares it: every fold selection equals `G/selection.json`. Five folds
   select `pg_otsu_e4.0_z0.8`. The fold that holds out n2rad51_HS_herm_001 selects on n2rad51_noHS_herm_001 and
   n2rad51_HS_herm_000 and picks `pg_tri_e4.0_z1.0` (0.766), which overcounts the held-out gonad (bias +2.84).
   Pooled over the six held-out gonads, the grid procedure reaches a truth-mask CCC of 0.684 [0.551, 0.788] and an
   F1 inside R of 0.603 [0.512, 0.655]. In the pairs with it, net3d and the tuned log are better in truth-mask,
   germline and traced CCC, truth-mask mean absolute error and both spot F1 (truth-mask CCC +0.194 and +0.209, F1
   inside R +0.255 and +0.245). The log defaults are better on every spot F1 and on no CCC, bias or mean absolute
   error. The D3 procedure is ahead on the point estimates of both primary numbers. Its F1 inside R (+0.024
   [+0.002, +0.041]), its detector-only CCC (+0.041 [+0.011, +0.074]) and its detector-only mean absolute error
   (-0.210 [-0.348, -0.088]) have intervals that exclude zero; its truth-mask CCC (+0.058 [-0.028, +0.177]) does
   not. The deployed recipe differs in no CCC (truth masks -0.083 [-0.197, +0.079]) and not in the F1 inside R
   (+0.011 [-0.046, +0.085]); it has the larger bias and mean absolute error and the higher whole-field F1. These
   are development estimates; nothing here is a win condition.
   (e) The F1 inside R in the pairs, which entry 6 (e) said needs the block added and the held-out runs of entry 4
   scored again. Commit ca8537d adds `spots_truth_masks`, `spots_truth_masks_class_a` and
   `spots_truth_masks_detector_only` to `logo.DETECTION_BLOCKS` and changes only `logo.py`, `bench_logo.py` and
   `tests/test_bench_logo.py`; a held-out folder without the `tm_spots_` columns is left out of those blocks. It
   was prepared during the train run, after the val run had written every val score (01:23:14 MDT;
   `G/logo/net3d_rescore_selection.json` was written at 01:25:12). It entered the tree only after the grid's report
   (02:18:22): every grid record and that report carry abc49ca without `-dirty`. It was committed at 02:19:08 MDT.
   The five detectors of entry 4 (b) were then run again on their held-out gonads at ca8537d, with the same
   parameter sets and fold models (trial lines 358 to 387; 02:19:10 to 02:29:59 MDT, then the report until
   02:30:19; `G/rescore.log`). All 30 runs equal their E1 originals (lines 162 to 173 and 214 to 231) in every
   result column of every table and in their spot tables and parameter hashes; only the run times differ. The six
   net3d re-run records do not carry the fold model sha256 that E1's records carry; the fold model files, hashed
   after the run, match lines 226 to 231. The report over them (`G/logo/report_f1R`) equals the E1 report in all
   145 pooled and 290 paired rows and the abc49ca report (`G/logo/report`) in all 174 and 435, and adds the F1
   inside R. Nothing in (e) changes a method, a metric or a selection.
   (f) Deviations from entry 6: none in the grid, the rule, the data or the scoring. Two things that entry 6 left
   open are recorded here. The driver ran every step without stopping, so the checks of (c) were made from the
   files afterwards: they confirmed the selection instead of gating each step (the held-out check needs the
   held-out runs in any case), and none failed. The code, the 30 re-runs and the report of (e) are not steps of
   entry 6, which named them as needed; they change no method, metric definition, data or selection.
   (g) Decision 2 stays the project lead's. He confirms or overrides the recommendation of entry 6 and its readings
   of (c) and (d): selection on the val gonads only, the rounding and tie-breaks, the failure rule, the budget
   counted for the grid alone, and the early stop read as not applying to a fixed grid. Entry 6 made these readings
   before the run, and he has not confirmed them yet. The readings that can be recomputed do not change the set:
   the six-gonad pick, the unrounded CCC and the 20 sets before the early stop all give `pg_otsu_e4.0_z0.8`, the
   tie-breaks only put it before a byte-identical set, and no run failed. The recommendation of entry 6 is kept,
   for him to confirm or override. Section 8 item 2 defines SpotMAX re-tuned as this grid chosen by the val rule.
   The grid ran and selected as declared, and every check passed. The selection is the set entry 6 expected, with
   the heat-shock undercount that decision 2 already named. Another set chosen now would be chosen after these
   results. What the results add for his choice: the undercount on both heat-shock hermaphrodites, known from the
   held-out runs of entry 4 (b), reproduces. Computed after the run, and not declared in entry 6: every grid set
   misses one of the two by at least 1.7 foci per truth nucleus. Out of sample the grid procedure's truth-mask CCC
   (0.684) lies between the deployed recipe's (0.601) and the D3 procedure's (0.742), and neither difference
   excludes zero. The D3 procedure is ahead of it on the F1 inside R and on the detector-only CCC and mean absolute
   error, with intervals that exclude zero (see (d)). The largest per-gonad difference is the fold that holds out
   n2rad51_HS_herm_001, where the grid procedure picked a triangle set (held-out CCC 0.472 against 0.734). If he
   picks another set, the pre-scoring entry records it and its reason as post hoc (entry 6 (a)). Decision 1 is
   unchanged: nothing is scored on test until he confirms one of its options.
