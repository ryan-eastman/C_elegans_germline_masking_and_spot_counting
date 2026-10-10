# Performance notes

Findings about where the pipeline spends its time, and the fixes, each with the measurement behind it.
A fix listed here changes no output: the default path stays byte-identical, and the evidence for that
is part of each entry.

## The SpotMAX spot stage on a flooded gonad (2026-10-10)

### Symptom

On `n2rad51_noHS_male_002` the spot stage (`germquant.spots.detect_spots`, backend `spotmax`) ran for
hours on one core, so the first freeze of the 8666f07 goldens produced no output for this gonad
(docs/BENCHMARK_PROTOCOL.md, changelog entry 1). The 30000 candidate cap of the flood guard
(`spots.max_spot_candidates`) did not help, because the time was spent before the cap is applied.

### Profile

`cProfile` around `detect_spots` with the `.nd2` foci channel (545 nm, index 2, read with
`germquant.io.read_stack`), the June label image
(`results_cv/20251105_N2_nohs_MALE _002/..._nuclei_labels.tif`, 570 nuclei; the pinned-commit golden
labels did not exist yet) and the spots section of the 8666f07 golden run manifest
(`threshold_triangle`, `effect_size_min` 3.0, `max_spot_candidates` 30000), capped at 20 minutes:

| step | time (s) |
|---|---|
| Gaussian preprocessing (`spotmax.pipe.preprocess_image`) | 5.5 |
| semantic segmentation (`spots_semantic_segmentation`) | 32.6 |
| the rest of `spotmax.features.find_local_peaks` before the filter (`peak_local_max` 3.2 of it) | 3.9 |
| `spotmax.filters.filter_valid_points_min_distance` on 479,901 peak candidates | 1157.4 at the cap, unfinished |

The minimum-distance filter took 96 percent of the capped run. It sorts the candidates by intensity
and walks them in that order; for every candidate it keeps, it builds masked arrays over ALL
candidates (a subtraction, a division, a square, a sum and a comparison) to find the ones inside the
kept one's ellipsoid. The cost is therefore the number of kept candidates times the number of
candidates. In the 1157.4 s it ran before the cap it completed 25,070 iterations, 46.2 ms each
(`numpy.ma` division 491.6 s, sum 379.7 s, subtraction 174.8 s, the loop body itself 104.6 s). The
candidate cap of the stage comes after this step, which is why it never helped.

### The old runtime, extrapolated

The loop runs once per kept candidate. The exact replica below keeps 207,340 of the 479,901
candidates (the same rows SpotMAX keeps, see the evidence), so SpotMAX's loop needs 207,340
iterations: 207,340 x 46.2 ms = 9,572 s for the filter alone. Timing SpotMAX's function without the
profiler on a random subset of 20,000 of the same candidates (18,929 kept) took 30.2 s and 29.1 s in
two runs, 77 to 80 ns per (kept candidate x candidate) pair, which predicts 7,650 to 7,930 s at full
size; the profiled figure includes the profiler and the larger arrays. With the 42.6 s before the
filter and the 164 s after it (measured with the new code, which does exactly the same work there),
the old stage needs about 7,900 to 9,800 s: 2.2 to 2.7 hours on one core, for this label image.

### The fix

`germquant.spots.peak_filter.filter_valid_points_min_distance` is an exact replica of SpotMAX's
function (SpotMAX 1.3.1, the pinned version): the same `argsort` call for the intensity order (so
tied intensities come out in the same order), the same greedy walk, and the same floating point
ellipsoid test (`np.square(diff / min_distance).sum(axis=1) < 1`) on the same integer difference
vectors. Only the search for the candidates near a kept one changes: a KD-tree on the coordinates
divided by `min_distance` proposes every pair closer than 1 + 1e-6 in those units, a superset of the
pairs the test can accept, and the SpotMAX test decides each pair. Inputs outside the domain this was
verified on (points that are not a 2D numpy array of signed integers, an intensity array of another
length, a radius that is not a positive finite number, coordinates above 1e8 radii) go to SpotMAX's
own function unchanged.

`detect_spots_spotmax` runs SpotMAX's `spot_detection` inside `fast_min_distance_filter()`, which
installs the replica in `spotmax.filters` for that call only and puts SpotMAX's function back
afterwards (nested and concurrent calls are reference counted). Nothing else in the chain changes.

### Evidence that the output is identical

1. Unit tests (`tests/test_peak_filter.py`): the replica against the installed SpotMAX function on
   120 random 2D and 3D inputs with tied intensities, duplicate points and nine radius settings
   (both return forms), on pairs that sit exactly on the ellipsoid boundary, the fallback domain, the
   context manager, and the whole SpotMAX stage with and without the replica on a flooded synthetic
   volume (hundreds of candidates, the cap hit), compared with `assert_frame_equal(check_exact=True)`.
2. Real candidates: random subsets of 5,000 and 20,000 of the 479,901 `n2rad51_noHS_male_002`
   candidates, plus 400 further random inputs, all on the replica's fast path: identical rows in
   identical order (on the committed code, SpotMAX 2.4 s and 29.1 s, replica 0.003 s and 0.016 s).
   On all 479,901 candidates, an independent exact reference that shares no neighbour search with
   the replica (a dictionary of integer coordinates and the 47 integer offsets that pass the SpotMAX
   ellipsoid test, with the same sort and walk) keeps the same 207,340 rows in the same order.
3. Three gonads end to end, old code (commit 3923789) against new code, each on the golden label
   image that the pipeline handed to the spots stage (`<stem>__nuclei_labels.tif`; its label set
   equals the `nucleus_id` set of `<stem>__nuclei.csv`) and the `.nd2` foci channel, with the golden
   spots section. Both tables were equal cell for cell, floats compared bit for bit, and both runs
   reproduced the golden `<stem>__spots.parquet` and `<stem>__spots.csv` (per-spot columns) and the
   `n_spots` column of `<stem>__nuclei.parquet` with no mismatch:

   | gonad | candidates into the filter | kept | spot rows | nucleus rows | old stage (s) | new stage (s) |
   |---|---|---|---|---|---|---|
   | n2rad51_noHS_herm_000 | 7,806 | 3,729 | 1,934 | 1,274 | 80.4 | 82.1 |
   | n2rad51_HS_herm_000 | 28,437 | 13,048 | 1,552 | 662 | 99.4 | 66.8 |
   | n2rad51_HS_male_001 | 53,883 | 25,303 | 4,858 | 442 | 326.6 | 164.6 |

   Only identity was checked; no agreement metric was computed on any gonad (n2rad51_HS_male_001 is
   a test gonad of the benchmark protocol). The stage times were measured on a shared workstation
   with other jobs running, one gonad at a time.
4. The synthetic goldens (`tests/test_golden.py`) pass unchanged.

### The new runtime

The new code on `n2rad51_noHS_male_002` with the same inputs as the profile, without the profiler:

| step | time (s) |
|---|---|
| Gaussian preprocessing | 5.0 |
| semantic segmentation | 23.7 |
| peak detection (`spot_detection`; the replica filter 0.7 s of it) | 6.6 |
| flood guard (54,487 candidates inside nuclei, capped at 30,000) | 0.0 |
| per-spot features on 30,000 candidates (`spots_calc_features_and_filter`) | 157.5 |
| z-split merge | 3.1 |
| tables | 3.3 |
| the whole stage | 199.8 |

About 3.3 minutes instead of 2.2 to 2.7 hours. The per-spot feature step (SpotMAX, about 5 ms per
candidate) is now the largest part; the cap bounds it at about 2.6 minutes.

### What did not change

- The gonad still floods. With the June labels, 54,487 candidates lie inside nuclei, more than the
  30,000 cap, so the stage logs `FLOODED: 54487 spot candidates inside nuclei exceed the cap
  max_spot_candidates=30000; keeping the brightest 30000 and dropping 24487. ...` and the count of
  this gonad is a floor, not a measurement. The warning is a log line only: the default path adds no
  QC flag and no column for it (flags enter `qc_flags` and would change the golden text).
- The candidate counts themselves come from SpotMAX's semantic segmentation (the triangle threshold
  over the aggregated nuclei); the flood is a property of the image and the recipe, which this fix
  does not touch.
- Memory: a full 2600 x 2600 x 71 stack peaks at about 11 GB in the spot stage (the profiled process
  peaked at 10.9 GB working set).
