# Final scoring: the one-time test scoring of the RAD-51 detectors

Written 2026-10-10 (package E1). **No step of this plan has been run.** It is the plan for the one test
scoring that `docs/BENCHMARK_PROTOCOL.md` allows per method family (sections 7, 8 and 10): the decisions to
take first, the frozen methods, the exact commands in order, how the win conditions are read, and the
changelog text to append before and after. No test gonad was scored, run or looked at while it was written.

The test gonads (protocol section 3): n2rad51_HS_herm_002, n2rad51_HS_male_001, n2rad51_noHS_herm_002,
n2rad51_noHS_male_002.

## 1. Decisions the project lead takes before anything runs

### Decision 1. Wait for the re-exported heat-shock male spots, or score now

* The two heat-shock males of train, n2rad51_HS_male_000 and n2rad51_HS_male_002, were exported from Imaris
  with region growing on; their spot truth is unusable (protocol section 4.1) and no heat-shock male has
  informed any choice: not the D3 searches, not net3d's training or selection, not the leave-one-gonad-out
  estimates. Val holds hermaphrodites only.
* The deployed recipe's heat-shock undercount, the in-house detector's main target, is not visible in train
  or val: on the two heat-shock hermaphrodites the deployed SpotMAX recipe over-counts (docs/BENCHMARK.md,
  "What the heat-shock gonads show"). The only heat-shock male with clean spot truth is n2rad51_HS_male_001,
  a test gonad. Scoring now therefore tests the main claim on one gonad that no method has been shaped on.
* Waiting means: re-export both projects without region growing, record their new size and sha256 in the
  split record with a changelog line (section 4.1), rebuild their truth files, and then decide whether they
  enter training. If they do, every learned or tuned method changes (net3d retrained on eight gonads, the
  classical searches rerun or at least rescored, the leave-one-gonad-out estimates redone on eight gonads):
  a declared revision before any test scoring, which costs a day of compute and the Imaris work. If they
  only enter the leave-one-gonad-out estimate (as evidence, not as training data), nothing frozen changes.
* Scoring now means the write-up states that no heat-shock male informed any choice and that the
  heat-shock-male result rests on one test gonad.

Status (2026-10-11, protocol changelog entry 6): the project lead replied "Yes to everything and merged" to the
message that listed the five decisions. That does not settle this choice by itself. It is read as waiting
for the re-export, the reading that keeps both options open: no test gonad is scored until he confirms one
option explicitly.

### Decision 2. Which SpotMAX configuration is the reference of section 7

Section 7: the reference is "the better of SpotMAX deployed and SpotMAX re-tuned on val (chosen on val, not
on test; both are reported on test)", under the primary metric (section 5.1 partition, Class A+B,
end-to-end). Changelog entry 2 left open which search supplies the re-tuned competitor. Package E1 scored
all 24 D3 SpotMAX sets on the two val gonads with the primary partition (trial log lines 174 to 197,
`logo_e1/val_scoring/spotmax/summary.csv`). The candidates, pooled over n2rad51_HS_herm_001 and
n2rad51_noHS_herm_001 (241 truth nuclei; with two gonads an interval only shows how far apart they are):

| candidate | parameters | val truth-mask CCC [interval] | bias | MAE | val spot F1 | how it would be chosen |
|---|---|---|---|---|---|---|
| SpotMAX deployed | triangle, radius 0.3 um, smoothing 0.08 um, effect size 3.0 | 0.556 [0.391, 0.679] | +2.91 | 2.93 | 0.428 | the 8666f07 recipe |
| D3 selection | triangle, 0.3, 0.12, 3.5 | 0.695 [0.512, 0.835] | +1.32 | 1.63 | 0.444 | the D3 rule (germline objective, three finalists on val) |
| entry 2's alternative | triangle, 0.35, 0.08, 4.0 | 0.723 [0.652, 0.762] | +1.00 | 1.68 | 0.420 | the best D3 finalist on val by the harness's older traced subset |
| best D3 set on val under the primary metric | otsu, 0.3, 0.08, 4.0 | 0.822 [0.734, 0.822] | -0.12 | 1.17 | 0.365 | section 6.4 read literally over the D3 search space (every set has a val score now) |

Every re-tuned candidate beats the deployed recipe on val under the primary metric, so by section 7 the
reference is a re-tuned configuration; which one is the decision. Points to weigh: the otsu set is the
literal 6.4 choice within the D3 space but has the lowest F1 of the four, in D3 every otsu set starved the
heat-shock training gonad (bias -0.88 to -1.41 foci per germline nucleus on n2rad51_HS_herm_000), and held out
in the leave-one-gonad-out estimate that gonad got a truth-mask bias of -2.33 from it; the leave-one-gonad-out
estimate of "SpotMAX re-tuned" measures the D3 procedure, which picked the D3 selection in four of six folds and
the otsu set in the two heat-shock folds; none of the re-tuned candidates was chosen under the primary metric by
a rule written before the search. Parameter files: `final_params/spotmax_retuned_primary.json` (sha256 84e27f45...) for the otsu set,
`spotmax_retuned_d3.json` for the D3 selection. In the leave-one-gonad-out estimate (docs/BENCHMARK.md,
"Leave-one-gonad-out") the re-tuning procedure (the D3 search replayed per fold) reached a truth-mask CCC of
0.742 [0.655, 0.812] against 0.601 [0.499, 0.718] for the deployed recipe; that difference, +0.141 [-0.043,
+0.292], does not exclude zero.

The section 8 item 2 grid (thresholding method in otsu, li, triangle, yen; `effect_size_min` 2.0 to 4.0;
`z_merge_gap_um` 0.6, 0.8, 1.0; 60 configurations) has never been run. Running it now is allowed (development
on train and val, val score for every configuration, section 6.4), takes several hours of SpotMAX, and would
replace the D3 search as the re-tuned competitor if the project lead prefers the pre-registered grid.

Status (2026-10-11, protocol changelog entry 6): a choice, which the reply "Yes to everything and merged" does not
settle by itself. Entry 6 declares the section 8 item 2 grid before it is run: 60 sets
(`pg_<method>_e<effect size>_z<gap>`, search `prereg_spotmax_grid`) scored on the four usable train gonads and the
two val gonads, selected on val by the pooled truth-mask CCC rounded to 3 decimals, then the spot F1 inside R,
then listing order. Recommendation: that selection becomes SpotMAX re-tuned, and by section 7 the reference, once
the project lead confirms it and the readings of entry 6 (c) and (d); the three D3 rows above stay reported as post
hoc. The grid contains the otsu row above (`pg_otsu_e4.0_z0.8`, in effect the same set), so if its val score
reproduces, the selection is that set unless a grid set not scored yet ranks above it on val. The choice is made on
train and val only (section 5); entry 6 (b) lists what was known before, including a pre-protocol sweep on two
test gonads.

Status (2026-10-11, protocol changelog entry 7): still the project lead's choice. The grid ran as entry 6 declared
it, on train and val only, and every check of entry 6 (c) passed. Its selection is `pg_otsu_e4.0_z0.8`
(`threshold_otsu`, effect size 4.0, z gap 0.8), in effect the otsu row of the table above: val truth-mask CCC 0.822
against 0.556 for the deployed recipe (rank 23 of 60). By section 7 it is the reference once he confirms it and the
readings of entry 6 (c) and (d); the readings that can be recomputed (val only, rounding, early stop) leave the set
unchanged. Its scores per gonad show the heat-shock concern above: truth-mask bias -2.33 on the train gonad
n2rad51_HS_herm_000 (CCC 0.587; the deployed recipe +0.91, CCC 0.843) and -0.77 on the val gonad
n2rad51_HS_herm_001, against +0.72 to +1.38 on the three train gonads without heat shock and +0.30 on the val gonad
n2rad51_noHS_herm_001. Computed after the run: every grid set misses one of the two heat-shock hermaphrodites by at
least 1.7 foci per truth nucleus. Two of the four test gonads are heat-shocked, and no heat-shock male informed
this choice. If the selection undercounts them as it undercounts both heat-shock hermaphrodites here, it is a weak
reference on half the test set, which makes W1 and W2 easier to meet. In the leave-one-gonad-out estimate the grid
procedure reaches a truth-mask CCC of 0.684 [0.551, 0.788], against 0.601 for the deployed recipe and 0.742 for the
D3 procedure, and neither difference excludes zero. The D3 procedure is ahead on the F1 inside R (+0.024 [+0.002,
+0.041]) and on the detector-only CCC (+0.041 [+0.011, +0.074]) and mean absolute error (-0.210 [-0.348, -0.088]),
intervals that exclude zero (docs/BENCHMARK.md, "Results of the pre-registered SpotMAX grid"). The recommendation
of entry 6 is kept: section 8 item 2 defines SpotMAX re-tuned as this grid chosen by the val rule, the grid ran and
selected as declared, and the selection is the set entry 6 expected, with the heat-shock undercount this decision
already named. Another set chosen now would be chosen after these results. If he picks another set, the pre-scoring
entry records it and its reason as post hoc.

### Decision 3. The limit of two test scorings per method family

`scripts/bench_detector.py --final` refuses a third test scoring of the same backend (protocol section 10:
"at most twice per method family"). SpotMAX deployed and SpotMAX re-tuned are both backend `spotmax`, so
scoring both (section 7 asks for both) uses both scorings of that family and leaves no declared revision for
either. Either accept that (the SpotMAX family is the reference, not a method under development), or decide
before scoring that they are two families and give the guard a family label (a code change, declared in the
pre-scoring changelog entry). This document assumes the first.

Status (2026-10-11, protocol changelog entry 6): accepted as written by the reply "Yes to everything and merged"
(the question put it as "accept that"); SpotMAX deployed and SpotMAX re-tuned use the two test scorings of the
`spotmax` family.

### Decision 4. The competitors of section 8 that do not exist yet

Spotiflow (zero-shot and fine-tuned), RS-FISH, Imaris with one global Quality threshold and the human blind
subset (section 9) are not implemented. The win conditions compare a new detector with the SpotMAX reference
only, so they do not block W1 and W2; but each later competitor is then scored after the first results are
known, which the write-up must say.

Status (2026-10-11, protocol changelog entry 6): accepted as written; these competitors are scored after the first
test results, and the write-up says so.

### Decision 5. The post hoc choices that stand

Changelog entry 2 (the D3 searches: germline objective, finalists only on val, the SpotMAX grid not run) and
entry 3 (net3d's selection rule) are post hoc under protocol rule 3, and so is the leave-one-gonad-out design
of entry 4. The write-up labels them. Confirm that the methods below are scored as they stand.

Status (2026-10-11, protocol changelog entry 6): accepted with a reading. The choices of entries 2, 3 and 4
stand and are labelled post hoc, and the methods below are scored as they stand except SpotMAX re-tuned:
decision 2 leaves that row open, and "the SpotMAX grid not run" no longer holds once entry 6 runs the grid.
The SpotMAX re-tuned row of section 2 stands only if decision 2 keeps the D3 selection.

## 2. The frozen methods

| competitor | backend | parameters (`--params` file) | weights | what it is |
|---|---|---|---|---|
| SpotMAX deployed | spotmax | `final_params/spotmax_deployed.json`: `{}` (sha256 e3566b3a...) | none | the 8666f07 recipe: the golden runs' `spots:` section, config hash 84bdda5a5230 |
| SpotMAX re-tuned (D3 selection; depends on decision 2) | spotmax | `final_params/spotmax_retuned_d3.json` (sha256 67db0cf7...) | none | `threshold_triangle`, radius 0.3 um, smoothing 0.12 um, effect size 3.5 (profile `bench_spotmax_tuned`) |
| tuned log | log | `final_params/log_tuned.json` (sha256 3a65e844...) | none | diameter 0.5 um, SNR 2.5, pre-gate 2, z merge off (profile `bench_log_tuned`) |
| net3d (train plus val) | net3d | `final_params/net3d_trainval.json` (sha256 2f9097e6...) | `models/net3d_trainval_20261010/model.pt`, sha256 d091b9ef58850783387a4f5aa0fea5ecc9310abc027145771936fecd0de001c0 (commit 240cc72) | the frozen recipe of net3d_20261010 (epoch 4 of the 40-epoch schedule, threshold 0.40) retrained on the six train and val gonads with usable spot truth (protocol section 6.1) |

`final_params/` is `C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/final_params/`; its five files
(the four above and the decision 2 alternative `spotmax_retuned_primary.json`) were written on 2026-10-10 with the
sha256 given here and under decision 2. SpotMAX re-tuned is not frozen until decision 2 is settled: the row above
is the D3 selection, and the grid selection of protocol changelog entry 6 replaces it if the project lead confirms
the recommendation. The grid selected `pg_otsu_e4.0_z0.8` (entry 7): `threshold_otsu`, `effect_size_min` 4.0 and
`z_merge_gap_um` 0.8, the `selected_params` under `designated` in `prereg_grid/selection.json`. Laid over the
golden `spots:` section (`detector.effective_section`), these three keys give the same section as
`spotmax_retuned_primary.json`: the two parameter sets run the same SpotMAX configuration under different parameter
hashes. A confirmed grid selection still gets its own three-key file, as below: `spotmax_retuned_primary.json`
belongs to a D3 set, and entry 6 (d) makes the D3 sets non-candidates once the grid is confirmed. If decision 2
picks another re-tuned configuration, write its file there first
(`{"thresholding_method": ..., "spot_radius_um": ..., "gauss_sigma_um": ..., "effect_size_min": ...}`; for a grid
set, its three keys `{"thresholding_method": ..., "effect_size_min": ..., "z_merge_gap_um": ...}`, the
`selected_params` under `designated` in `prereg_grid/selection.json`), use it in step 4.4, and name it in the
pre-scoring entry.

## 3. The order

1. Pre-flight (no metric is computed).
2. The pre-scoring changelog entry (protocol section 10: code commit, config hash and weight sha256 are
   written before test is scored), committed.
3. SpotMAX deployed, then the golden check of its detections.
4. SpotMAX re-tuned: the configuration that decision 2 settles (under the recommendation of protocol changelog
   entry 6, the grid selection `pg_otsu_e4.0_z0.8`, run and selected as entry 7 records, which must still be
   confirmed and written to `final_params/` before step 1).
5. Tuned log.
6. net3d (train plus val), on the GPU under the workstation's lock.
7. The comparison (pooling only) and the win conditions.
8. The post-scoring changelog entry, committed.

Every method is frozen before the first is scored, so no test number can shape another method; the
reference is scored first only because everything else is compared with it.

## 4. The commands

Run from the repository root of the commit named in the pre-scoring entry, in PowerShell with the
repository's environment active (`python` is its interpreter and imports this checkout's `germquant`).

### 4.1 Pre-flight

```powershell
git status --porcelain                 # must print nothing (--final refuses uncommitted code)
git rev-parse HEAD                     # goes into the pre-scoring entry
$env:CUDA_VISIBLE_DEVICES = "-1"
# the one test gonad without a bench cache: images and labels only, no metric
python scripts/bench_cache.py --tags n2rad51_noHS_male_002 --verify-sha256
# the golden run of that gonad must be the deployed pipeline (git sha 8666f07..., config hash 84bdda5a5230)
python -c "import json; m = json.load(open(r'C:/Users/ryane/germquant_golden/8666f07/n2rad51_noHS_male_002/run_manifest.json')); print(m['git_sha'], m['config_hash'])"
# the net3d weights are the frozen ones (d091b9ef...)
python -c "from germquant.provenance import file_sha256; print(file_sha256(r'C:/Users/ryane/germquant_golden/models/net3d_trainval_20261010/model.pt'))"
# the parameter files are the frozen ones (sha256 as in section 2 and under decision 2; a new file for the
# grid selection of decision 2 has its sha256 printed here, and the pre-scoring entry names it)
Get-FileHash -Algorithm SHA256 C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/final_params/*.json
# protocol section 4.4 on the test gonads: the geometric truth assignment against the lab's mask channel;
# agreement should match the train and val gonads (99.35 to 99.99 percent); there the differences were
# boundary spots, spots whose mask-channel value is a removed mask fragment (the larger group), and 3 spots
# that one region-growing gonad assigns to different nuclei
python scripts/check_truth_assignment.py --tags n2rad51_HS_herm_002 n2rad51_HS_male_001 n2rad51_noHS_herm_002 n2rad51_noHS_male_002
# no test scoring has happened
Test-Path C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/final_runs.jsonl   # must be False
```

The golden run of n2rad51_noHS_male_002 was made by the pinned commit 8666f07 (finished 2026-10-10 after 205.5
minutes in the flooded spot stage), so the deployed competitor on that gonad is the 8666f07 pipeline; the
check above confirms it. The same gonad rerun by the code with the KD-tree peak filter is byte-identical to it.

### 4.2 The pre-scoring changelog entry

Append to section 11 of docs/BENCHMARK_PROTOCOL.md (fill the brackets), then commit:

> N. [date]. Test scoring declared before any test result (protocol section 10). Methods, frozen at commit
> [git sha]: SpotMAX deployed (backend spotmax, params `{}`, the golden `spots:` section, config hash
> 84bdda5a5230); SpotMAX re-tuned [configuration and params file sha256]; the SpotMAX reference of section 7
> is [deployed or re-tuned], chosen on val under the primary metric ([val numbers]); tuned log (backend log,
> params file sha256 3a65e844...); net3d (backend net3d, model.pt sha256 [sha256], the frozen recipe of
> net3d_20261010 retrained on the six train and val gonads with usable spot truth). Decisions 1 to 5 of
> docs/FINAL_SCORING.md: [one line each]. No test result has been computed at this point.

### 4.3 SpotMAX deployed

```powershell
$env:CUDA_VISIBLE_DEVICES = "-1"
$F = "C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07"
python scripts/bench_detector.py --backend spotmax --params $F/final_params/spotmax_deployed.json `
    --split test --final --detector-only --jobs 2 `
    --reason "one-time test scoring (protocol section 10): SpotMAX deployed, the 8666f07 recipe" `
    --out $F/final/spotmax_deployed
```

Two SpotMAX workers need about 25 GB of memory; n2rad51_noHS_male_002 may take much longer than the others
(its detector flooded before the peak-filter fix). Then check, with no metric, that the detections are the
golden spot tables (protocol section 8 item 1: the deployed competitor's outputs are the goldens):

```powershell
@'
import pandas as pd
from germquant.benchmark import detector as bd
from germquant.benchmark.manifest import load_manifest
E = {e.tag: e for e in load_manifest(bd.DEFAULT_MANIFEST)}
for t in ["n2rad51_HS_herm_002", "n2rad51_HS_male_001", "n2rad51_noHS_herm_002", "n2rad51_noHS_male_002"]:
    run = pd.read_csv(f"C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/final/spotmax_deployed/detections/{t}__spots.csv", float_precision="round_trip")
    gold = bd.read_golden_spots(E[t].results_dir, bd._entry_stem(E[t]))
    print(t, bd.compare_with_golden(run, pd.DataFrame(), gold)["identical"])
'@ | python -
```

Every line must say True. A False means the scored detections are not the goldens; stop and record it.

### 4.4 SpotMAX re-tuned

```powershell
python scripts/bench_detector.py --backend spotmax --params $F/final_params/spotmax_retuned_d3.json `
    --split test --final --detector-only --jobs 2 `
    --reason "one-time test scoring (protocol section 10): SpotMAX re-tuned, the decision 2 choice" `
    --out $F/final/spotmax_retuned
```

This is the second and last `--final` run the guard allows for backend `spotmax` (decision 3). The `--params` file
depends on decision 2: `spotmax_retuned_d3.json` is the D3 selection. If the project lead confirms the grid
selection of protocol changelog entry 6, its own file (written to `final_params/` first, section 2) takes its place
here, and the pre-scoring entry names it.

### 4.5 Tuned log

```powershell
python scripts/bench_detector.py --backend log --params $F/final_params/log_tuned.json `
    --split test --final --detector-only --jobs 2 `
    --reason "one-time test scoring (protocol section 10): tuned log (D3 selection)" `
    --out $F/final/log_tuned
```

### 4.6 net3d (train plus val)

```powershell
Remove-Item Env:CUDA_VISIBLE_DEVICES -ErrorAction SilentlyContinue
$env:PYTORCH_CUDA_ALLOC_CONF = "expandable_segments:True"     # a no-op on Windows, kept for the workstation rule
python C:/Users/ryane/germquant_golden/gpu_lock.py run --name "final test scoring: net3d_trainval" -- `
    python scripts/bench_detector.py --backend net3d --params $F/final_params/net3d_trainval.json `
    --split test --final --detector-only `
    --reason "one-time test scoring (protocol section 10): net3d retrained on train plus val by the frozen recipe" `
    --out $F/final/net3d_trainval
```

### 4.7 The comparison

Pooling only; every number comes from the four runs above. New methods first, so that every paired
difference reads new minus reference:

```powershell
$env:CUDA_VISIBLE_DEVICES = "-1"
python scripts/bench_logo.py report `
    --tags n2rad51_HS_herm_002 n2rad51_HS_male_001 n2rad51_noHS_herm_002 n2rad51_noHS_male_002 `
    --heldout net3d_trainval=$F/final/net3d_trainval log_tuned=$F/final/log_tuned `
              spotmax_retuned=$F/final/spotmax_retuned spotmax_deployed=$F/final/spotmax_deployed `
    --out $F/final/report
```

`report` refuses a folder that lacks a gonad. A method that produced no output for a gonad (its run lists the
gonad under `skipped`) gets a row saying so in every table (protocol section 10); the win conditions are then
evaluated over the gonads where both methods produced output: rerun `report` with `--tags` holding those
gonads for that pair and label the result with their number (3 of 4, say).

Reading the win conditions (protocol section 7) from `final/report/logo_paired.csv`, with R the reference of decision 2:

* W1 (counts): the row with `pair` = "net3d_trainval - R", `block` = truth_masks, `stat` = ccc: `diff` at least
  0.03 and `lo` above 0 (the paired cluster bootstrap by gonad, 1000 resamples, seed 20251105). The same row
  for "log_tuned - R" says whether the tuned log wins.
* W2 (spots): `block` = spots_truth_masks (protocol section 5.2: matching restricted to R, the truth masks;
  spots_all is the whole-field score, kept for comparison only), `stat` = f1: `diff` not below -0.02 (point estimates; the interval is
  reported). `logo.DETECTION_BLOCKS` holds this block since commit ca8537d (`"spots_truth_masks": "tm_spots_"`,
  with `spots_truth_masks_class_a` and `spots_truth_masks_detector_only`), so `bench_logo.py report` pools and
  pairs it; a held-out folder without the `tm_spots_` columns (written before commit e19a029) is left out of these
  blocks (protocol changelog entry 7). If `logo_paired.csv` has no spots_truth_masks row for a pair, W2 cannot be
  read for it: stop and record it.
* W3 does not apply: no method changes the nucleus segmentation.
* Fragility (section 6.2): for every pair, compare the direction of the test difference with the
  leave-one-gonad-out difference of docs/BENCHMARK.md ("Leave-one-gonad-out"); a test result whose direction
  the leave-one-gonad-out estimate contradicts is reported as fragile. The leave-one-gonad-out directions, for
  reference: net3d and the tuned log are above both SpotMAX configurations on the truth-mask CCC, the germline
  and traced CCC and spot F1 (the paired intervals exclude zero); net3d against the tuned log shows no
  difference (truth-mask CCC -0.015 [-0.084, +0.073]); SpotMAX re-tuned against deployed shows none except a
  smaller mean absolute error. These SpotMAX re-tuned directions are those of the D3 replay; if decision 2 makes
  the grid selection the re-tuned competitor, its pairs are read from the grid replay of protocol changelog entries
  6 and 7 instead. Both replays, with the section 5.2 F1 inside R for every pair, are in
  `prereg_grid/logo/report_f1R/logo_paired.csv` (docs/BENCHMARK.md, "Results of the pre-registered SpotMAX grid").
  There net3d and the tuned log are above all three SpotMAX configurations on the F1 inside R (+0.221 to +0.255,
  the intervals exclude zero) and better than the grid procedure on every CCC, the truth-mask mean absolute error
  (lower) and both spot F1. On the F1 inside R no difference survives between net3d and the tuned log (+0.010
  [-0.058, +0.094]) or between SpotMAX re-tuned and deployed (+0.014 [-0.073, +0.069]). The grid procedure minus
  the deployed recipe (the report lists "spotmax_deployed - spotmax_grid", with the opposite sign): truth-mask CCC
  +0.083 [-0.079, +0.197], F1 inside R -0.011 [-0.085, +0.046], mean absolute error -1.016 [-1.518, -0.344],
  whole-field F1 -0.049 [-0.098, -0.005].
* Report next to the primary numbers: Class A and Class B truth (truth_masks_class_a; per_gonad.csv
  `tm_class_b_*`), the detector-only count (truth_masks_detector_only), the germline and traced partitions,
  and spot F1 without Class B, every one per gonad and pooled (section 5.4), with the run records of
  final_runs.jsonl (versions, git sha, params sha256) and the split file sha256.

### 4.8 The post-scoring changelog entry

Append to section 11 (fill the brackets), then commit:

> N+1. [date]. Test scored once per method (final_runs.jsonl lines [n to m], commit [git sha]): SpotMAX
> deployed, SpotMAX re-tuned, tuned log, net3d (train plus val), on [4 of 4] test gonads [list any gonad
> without output]. Primary metric (section 5.1 partition, Class A+B, end-to-end), pooled CCC [95 percent
> interval]: SpotMAX deployed [..], SpotMAX re-tuned [..], tuned log [..], net3d [..]. W1 net3d against the
> reference [R]: difference [..] [interval]: [met / not met]; W2: F1 difference [..]: [met / not met]. Tuned
> log against R: W1 [..], W2 [..]. Decision (section 7): [the new detector replaces SpotMAX / the SpotMAX
> reference stays deployed / mixed result, not replaced]. Fragile results (section 6.2): [none / which].
> These are the first test results of the new methods; any change from here on is post hoc (rule 3) and a
> second scoring needs a revision declared here first (section 10).

## 5. What must not happen

* No method is changed, tuned or re-selected after any test number is seen; a revision is declared in the
  changelog before its (second and last) scoring.
* No test gonad is scored with `--no-trial-log`, with `scripts/benchmark_spots.py`, or outside `--final`.
* The SpotMAX reference is not chosen, or changed, on test numbers.
