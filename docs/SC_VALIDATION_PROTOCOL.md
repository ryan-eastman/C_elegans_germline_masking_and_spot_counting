# SC ribbon validation protocol (pre-registration)

Registered 2026-10-23 (changelog entry 1). When this was written, no SC output of any method had been compared
with any SC truth; entry 1 lists everything that had been looked at. The method under test is `sc_straighten`
v1 (`DESIGN.md`). The split record is `C:/Users/ryane/sc_ribbon_work/sc_split.json` (sha256
db328c1cf4bc2b08dda32027a0e2567a128700a26facac4f7dc47738f67ae9eb, LF line endings), built by
`scripts/make_sc_split.py` from the inventory metadata and the bytes of the truth copies. The build copies it to
`tests/golden/sc_split.json` and this file to `docs/SC_VALIDATION_PROTOCOL.md`, both byte-identical. The file is
ASCII (rule R9) and follows the layout of `docs/BENCHMARK_PROTOCOL.md` (the RAD-51 protocol).

## 1. What this document is for

The stage traces SC ribbons (SYP) in 3D, straightens them and measures them. This protocol fixes, before any
result exists: which images may be used for development and which are sealed; what the truth is and what it is
not; how truth is matched to our nuclei and ribbons; the metrics; the success criteria and the decisions they
drive; and who scores the test set, and when. The point is that the numbers reported later cannot have been shaped
by a choice of split, metric or threshold made after seeing them.

Rules of the document:

1. A result is any metric of section 6 computed for any method (sc_straighten, a baseline or a comparator) on any
   test image, on the live filament file, or on any sealed synthetic seed.
2. Until the first result, the sections above the changelog may be edited; each edit gets a changelog line that
   says whether development results motivated it.
3. From the first result on, nothing above the changelog is edited. A change of method, metric or data is
   appended to the changelog with the date, the reason, and whether a test result had been seen; a change made
   after a test result is labelled post hoc there and in every write-up that uses it.
4. The split never changes. A gonad leaves the protocol only for a defect on the truth side (a corrupted file, an
   annotation that does not belong to the image, zone letters missing or out of order), recorded in the
   changelog before any metric is computed, never because of a score.

## 2. The truth that exists, and what it is not

### 2.1 T1: SC zone truth (the primary real-data truth)

- **What it is.** For each gonad, one hand-placed Imaris Measurement Points polyline along the gonad axis (7 to 21
  points, named A, B, C, ... in row order), and a master table that names, per gonad, the points where the SC zones
  change: `SC_fragstart`, `SC_fullstart`, `SC_fraglate`, `Pachy_end` (and `gonad_start`, the first point). The lab's
  R code (`20260121_SC_zone_analysis.R`, the Libuda-lab gonad linearization zone method) labels the polyline's
  segments: None before SC_fragstart, SCfrag up to SC_fullstart, SCfull up to SC_fraglate, LPfrag up to Pachy_end,
  with lengths in x and y only.
- **Sets.** N2, January 2026 tracing (master `20260120_measurement_points_SC.xlsx`): 39 gonads from 2025-10-21 (6),
  2025-11-05 (12) and 2026-01-05 (21). spo-11 IR (DLW251, 12) and CCW69 (9 listed, 8 with exports): reserved.
  DLW190 (42): excluded (38 images missing, staining notes, SC marker not checked).
- **Resolution.** Consecutive points are about 20 um apart (median 19.8 um, IQR 13.9 to 28.1, maximum 93.5); a zone
  boundary is known only to the nearest drawn point, several nuclear diameters. Section 4.2 scores accordingly.
- **What it is not.** No nucleus and no ribbon was measured. It says nothing about SC length per nucleus or per
  chromosome. It is one annotator's judgement of SYP morphology along the gonad; the only repeatability data are
  the DLW190 January and April tracings of the same gonads.
- **The lab's derived tables are not truth here.** Read as code: the R summary takes each zone's start as the
  cumulative length at the END of the zone's first segment, one point after the named point, while the segment
  labels start the zone at the named point. This protocol uses the named points.

### 2.2 T2: live Imaris Filament truth (descriptive check only)

- One file has saved geometry: `2025_12_05_dlw193_mch_syp3_male1_nohs image 1_twonucleiSCtraced.ims` (E: copy;
  identical object tables to the NAS copy `..._midpachytene.ims`), DLW193 mCherry::SYP-3, live, male, no heat shock,
  time index 60 of 61: two mid-pachytene male nuclei, 5 unbranched Filaments (autopath) each, plus 3 single-vertex
  stubs, which are ignored. The traced time point exists as
  `E:/2025_12_05_dlw193_mch_syp3_male1_nohs image 1_twonucleiSCtraced_time61_allz.tif` (65 planes stored as ImageJ
  frames; read with z 0.223077 um and xy 0.108 um set explicitly).
- **What it is not.** Fixed-image truth: the strain, marker, modality and SNR differ from the N2 IF images, and it
  is 10 ribbons in 2 nuclei. Autopath vertices sit on voxel centres. The three other Filament files are empty
  (re-save requested, section 11).

### 2.3 T3: synthetic truth

Exact-render phantoms from `germquant.sc_straighten.synth`, with these fixed distributions (a change after
registration is a protocol edit under rule 2):

- Grid: voxel 0.2 x 0.108333 x 0.108333 um (z, y, x).
- PSF: Gaussian, sigma_xy 0.149 um and sigma_z 0.365 um (FWHM 0.35 and 0.86 um, Design A's measurement on 232 flat
  strands of herm_000), summed over truth points 0.005 um apart; the truth is never rounded to voxels.
- Intensities: camera offset 100 counts; nucleoplasm U(30, 90) counts above offset inside the nucleus (Design B's
  measured 130 to 190 counts); each ribbon's peak above the nucleoplasm after the PSF is the scene peak x
  U(0.7, 1.3); scene peak 430 counts unless the family says otherwise.
- Noise: Gaussian with variance 4.04 x I - 394 counts^2 at the noise-free absolute intensity I (offset included),
  floored at 25 counts^2 (the herm_000 fit; the floor matches a read noise of about 5 to 6 counts).
- Geometry: smooth ribbons with curvature radius at least 1.0 um. Oocyte nuclei: radius U(2.2, 2.6) um, 6 ribbons
  of length U(4.5, 7.0) um, one of them U(4.5, 5.7) um (the X). Spermatocyte nuclei: 5 ribbons. Ribbons lie in a
  shell between 0.55 and 0.95 of the radius (peripheral, as at pachytene), in random orientations. Labels are
  spheres of radius R + U(-0.2, +0.4) um.
- Families: iso, dense, gaps, foci, pairs, noise, snr (section 6.1).
- Development may use any seed below 7,000,000. Sealed seeds (scored once, section 10): iso 7100000 to 7100199,
  dense 7200000 to 7200099, gaps 7300000 to 7300099, foci 7400000 to 7400099, pairs 7500000 to 7500099, noise
  7600000 to 7600049, snr 7700000 to 7700099 (also in sc_split.json). The generator is the code at the scoring
  commit.
- What it is not: real geometry. Decomposition rates on real nuclei can differ; that is what T4 is for.

### 2.4 T4: independent hand traces on fixed images (requested, not yet existing)

Defined in section 11. Until T4 exists, no per-ribbon or per-nucleus length claim is made for fixed images beyond
"validated on synthetic data".

## 3. The split

### 3.1 The rule (applied before any SC output was compared with any SC truth)

1. 2025-11-05: the RAD-51 train and val gonads are development; the four RAD-51 test gonads are SC test (the
   RAD-51 seal: test, never development).
2. 20251021_n2_nohs_HERM_2 is development: example rows of its export were printed during the inventory.
3. The four 2026-01-05 heat-shock males with local copies are development.
4. Every other N2 gonad is free. Within its (session, condition, sex) cell the free gonads are ordered by
   sha256("sc-split-v1:" + gonad id) and take, in that order, the roles of the cell's registered sequence:

   | cell | free gonads | role sequence |
   |---|---|---|
   | 2025-10-21 noHS herm | 2 | test, reserved |
   | 2025-10-21 noHS male | 3 | dev, test, reserved |
   | 2026-01-05 HS herm | 7 | test, dev, reserved, test, dev, reserved, reserved |
   | 2026-01-05 HS male | 5 | test, reserved, test, reserved, reserved |
   | 2026-01-05 noHS herm | 2 | dev, test |
   | 2026-01-05 noHS male | 3 | dev, test, reserved |

5. spo-11 IR and CCW69 gonads are reserved; DLW190 gonads are excluded.
6. The live Filament file is test (check S3).
7. Synthetic: development seeds and sealed seeds as in 2.3.

### 3.2 The registered roles (N2)

Test, 12 gonads (3 per condition x sex cell, from two or three sessions each):

| gonad (Madeleine's id) | image stem | session | cond. | sex | rule | notes |
|---|---|---|---|---|---|---|
| 20251105_n2_hs_HERM_2 | `20251105_n2_hs_HERM _002` | 20251105 | HS | herm | 1 | n2rad51_HS_herm_002 (RAD-51 test) |
| 20260105_n2_hs_2_HERM_2 | `20260105_n2_hs_2_HERM_2` | 20260105 | HS | herm | 4 | |
| 20260105_n2_hs_HERM_1 | `20260105_n2_hs_HERM_1` | 20260105 | HS | herm | 4 | |
| 20251105_n2_hs_MALE_1 | `20251105_n2_hs_MALE _001` | 20251105 | HS | male | 1 | n2rad51_HS_male_001 (RAD-51 test) |
| 20260105_n2_hs_2_MALE_1 | `20260105_n2_hs_2_MALE_1` | 20260105 | HS | male | 4 | |
| 20260105_n2_hs_2_MALE_2 | `20260105_n2_hs_2_MALE_2` | 20260105 | HS | male | 4 | |
| 20251021_n2_nohs_HERM_3 | `20251021_n2_nohs_HERM_3` | 20251021 | noHS | herm | 4 | SC on 545 |
| 20251105_N2_nohs_HERM_2 | `20251105_N2_nohs_HERM _002` | 20251105 | noHS | herm | 1 | n2rad51_noHS_herm_002 (RAD-51 test) |
| 20260105_n2_nohs_HERM_4 | `20260105_n2_nohs_HERM_4` | 20260105 | noHS | herm | 4 | |
| 20251021_n2_nohs_MALE_3 | `20251021_n2_nohs_MALE_3` | 20251021 | noHS | male | 4 | SC on 545 |
| 20251105_N2_nohs_MALE_2 | `20251105_N2_nohs_MALE _002` | 20251105 | noHS | male | 1 | n2rad51_noHS_male_002 (RAD-51 test) |
| 20260105_n2_nohs_MALE_2 | `20260105_n2_nohs_MALE_2` | 20260105 | noHS | male | 4 | |

Development, 18 gonads (4 in each cell, 6 in the heat-shock male cell):

| gonad (Madeleine's id) | image stem | session | cond. | sex | rule | notes |
|---|---|---|---|---|---|---|
| 20251105_n2_hs_HERM_1 | `20251105_n2_hs_HERM _001` | 20251105 | HS | herm | 1 | n2rad51_HS_herm_001 (RAD-51 val) |
| 20251105_n2_hs_HERM_3 | `20251105_n2_hs_HERM _` | 20251105 | HS | herm | 1 | n2rad51_HS_herm_000 (RAD-51 train) |
| 20260105_n2_hs_1_HERM_1 | `20260105_n2_hs_1_HERM_1` | 20260105 | HS | herm | 4 | |
| 20260105_n2_hs_HERM_2 | `20260105_n2_hs_HERM_2` | 20260105 | HS | herm | 4 | |
| 20251105_n2_hs_MALE_2 | `20251105_n2_hs_MALE _002` | 20251105 | HS | male | 1 | n2rad51_HS_male_002 (RAD-51 train) |
| 20251105_n2_hs_MALE_3 | `20251105_n2_hs_MALE _` | 20251105 | HS | male | 1 | n2rad51_HS_male_000 (RAD-51 train) |
| 20260105_n2_hs_1_MALE_1 | `20260105_n2_hs_1_MALE_1` | 20260105 | HS | male | 3 | local copy |
| 20260105_n2_hs_MALE_1 | `20260105_n2_hs_MALE_1` | 20260105 | HS | male | 3 | local copy |
| 20260105_n2_hs_MALE_2 | `20260105_n2_hs_MALE_2` | 20260105 | HS | male | 3 | local copy |
| 20260105_n2_hs_MALE_3 | `20260105_n2_hs_MALE_3` | 20260105 | HS | male | 3 | local copy |
| 20251021_n2_nohs_HERM_2 | `20251021_n2_nohs_HERM_2` | 20251021 | noHS | herm | 2 | SC on 545; seen |
| 20251105_N2_nohs_HERM_1 | `20251105_N2_nohs_HERM _001` | 20251105 | noHS | herm | 1 | n2rad51_noHS_herm_001 (RAD-51 val) |
| 20251105_N2_nohs_HERM_3 | `20251105_N2_nohs_HERM _` | 20251105 | noHS | herm | 1 | n2rad51_noHS_herm_000 (RAD-51 train); seen by three agents |
| 20260105_n2_nohs_HERM_1 | `20260105_n2_nohs_HERM_1` | 20260105 | noHS | herm | 4 | |
| 20251021_n2_nohs_MALE_2 | `20251021_n2_nohs_MALE_2` | 20251021 | noHS | male | 4 | SC on 545 |
| 20251105_N2_nohs_MALE_1 | `20251105_N2_nohs_MALE _001` | 20251105 | noHS | male | 1 | n2rad51_noHS_male_001 (RAD-51 train) |
| 20251105_N2_nohs_MALE_3 | `20251105_N2_nohs_MALE _` | 20251105 | noHS | male | 1 | n2rad51_noHS_male_000 (RAD-51 train) |
| 20260105_n2_nohs_MALE_3 | `20260105_n2_nohs_MALE_3` | 20260105 | noHS | male | 4 | |

Reserved, 9 N2 gonads (sealed, not scored under this version): 20260105_n2_hs_2_HERM_1, 20260105_n2_hs_2_HERM_4,
20260105_n2_hs_HERM_3, 20260105_n2_hs_1_MALE_2, 20260105_n2_hs_1_MALE_3, 20260105_n2_hs_2_MALE_3,
20251021_n2_nohs_HERM_4, 20251021_n2_nohs_MALE_4, 20260105_n2_nohs_MALE_1. Also reserved: the 12 spo-11 IR gonads
and the 9 CCW69 gonads. Excluded: the 42 DLW190 gonads (their January and April tracings serve only as the human
repeatability reference of 6.2).

### 3.3 Why this split

- Every condition x sex cell has three test gonads from at least two imaging sessions, so the headline covers both
  conditions, both sexes and the session-to-session variation, and no cell can be dropped by tuning.
- The four RAD-51 test gonads stay test. Their `.ims` files are never opened; their polylines enter the pipeline
  frame through the nd2 stage rule (4.1), which reproduced the eight sibling `.ims` headers of that session to
  0.02 um.
- Development holds the eight RAD-51 train and val gonads (their germquant goldens exist, herm_000 has been seen
  by three agents) and gets six heat-shock males, because heat-induced late SC fragmentation in spermatocytes is the
  lab's main question and the RAD-51 split had none in val.
- 2025-10-21 has the SC on 545 (channels swapped): one gonad per sex of that session is development, so the session
  channel map is built and checked there before any test gonad of the session is run.
- The reserved gonads keep a clean second test set from the same sessions (and the IR and CCW69 sets an external
  one) for a later version or for a third attempt, which this test set does not allow (section 10).
- With 12 test gonads the intervals are wide; development results are reported next to the test results as
  supporting evidence, never instead of them.

### 3.4 What development may and may not use

- May: the 18 development images (germquant runs for the ten without one, sc_straighten runs, viewing, any
  statistic), their zone truth through the guarded loader from registration on, development synthetic seeds, the
  inventory metadata.
- May not, before the scoring of section 10: run sc_straighten on a test or reserved image; read a test or reserved
  zone export or their rows of a master table; read the vertices or statistics of the live Filament file; generate
  a sealed seed; look at the germquant outputs of a test image beyond a completion and channel-map check. The code
  enforces this (DESIGN 2.4): the command refuses sealed images without `--allow-sealed`, and the truth loaders
  raise `SealedError` without `final=True`.

## 4. Frames and matching rules

### 4.1 Zone polyline into the pipeline frame

- Read the gonad's export, sheet `Position` (one header row skipped, as the R code does), columns `Position X`,
  `Position Y`, `Name`; points in row order. Row k must be named by the k-th letter (A, B, ...); otherwise it is a
  truth defect (rule 4).
- Continuous index `i = (coord - ExtMin) / v_imaris - 0.5` (Imaris positions are voxel centres) and pipeline
  um `= i x v_nd2`, separately in x and y; `v_imaris = (ExtMax - ExtMin) / size` from the `.ims` header attributes
  (`germquant.validate.imaris_ims.read_geometry`; no dataset is read), `v_nd2` from the run (0.108333 um). For the
  four RAD-51 test gonads the `.ims` is never opened: ExtMin xy = the nd2 stage xy of the first frame and
  `v_imaris` xy = the nd2 voxel rounded to 0.001 um (0.108), the rule of `sc_ribbon_work/scripts/match_sc_to_images.py`
  with the stage read by `scripts/image_extents.py`. On every development gonad both rules are computed and must
  agree within 0.05 um in xy before scoring.
- Arc length s along the converted polyline, in xy only (the R code's convention). It is 0.31 % longer than the R
  code's lengths (0.108333 / 0.108); positions are compared in this frame and zone lengths normalised.

### 4.2 Zone label of each nucleus, analysis population, guard bands

- **Nuclei.** Rows of the germquant run's `nuclei.csv` with `in_germline` true whose sc_straighten status is
  `traced` or `no_sc`. Failed nuclei are excluded and counted per gonad; a gonad with more than 5 % failed nuclei is
  flagged in every table.
- **Projection.** `(s_i, r_i) = staging.zones.project_to_polyline((centroid_x_um, centroid_y_um), polyline)`.
- **Off-axis cut.** c = 20 um; when at least 20 nuclei have r <= 20 um and 0 <= s <= L, c = min(20, 2.5 x their
  median r) (the staging rule). Keep r_i <= c and 0 <= s_i <= s(Pachy_end).
- **Segment and zone.** With points p_1 ... p_n and the point indices f, u, l, e of SC_fragstart, SC_fullstart,
  SC_fraglate and Pachy_end: the nucleus lies in segment b (s(p_b) <= s_i < s(p_(b+1)); s_i = s(p_n) counts as
  segment n - 1), and its zone is None if b < f, SCfrag if f <= b < u, SCfull if u <= b < l, LPfrag if
  l <= b <= e. Letters must satisfy f <= u <= l <= e, else a truth defect.
- **Boundaries.** A point k in {f, u, l} is a boundary when 1 < k < n and the zones on both sides of it exist.
- **Guard band** of boundary k: s(p_k) - 0.5 (s(p_k) - s(p_(k-1))) < s < s(p_k) + 0.5 (s(p_(k+1)) - s(p_k)), half
  the neighbouring point spacing on each side. Nuclei in a guard band are left out of S2a and S2b and kept for S2c.

### 4.3 Live Filament truth

- Truth ribbons: the 10 Filaments of time index 60 with at least 2 vertices; nucleus A = filament ids 100000060 to
  100000064, nucleus B = 100000068 to 100000072 (the inventory's grouping, confirmed at scoring by the two vertex
  clusters).
- Frame: `index = (p - ExtMin) / v - 0.5` with v = (0.108, 0.108, 0.223077) um (x, y, z) and ExtMin from the
  `.ims` header attributes; um = index x v, which indexes the TIFF directly.
- Territories: for each nucleus, the ball around the centroid of its truth vertices with radius = the largest
  vertex distance from that centroid + 0.5 um, split between A and B by the plane bisecting their centroids; labels
  1 and 2. Only the centres and extents of the truth enter our run, never their lengths.
- Our run: `straighten_image` on the TIFF (the one channel as central_element) with these labels, the frozen
  parameters, no spot table, germ_cell spermatocyte (expected 5). The truth ribbons are also measured by the same
  code, once as drawn (no refinement, ends as drawn) and once refined (re-centring and the end rule).
- Tolerance between a truth sample and our samples: sqrt((d_xy / 0.3)^2 + (d_z / 0.6)^2) <= 1, distances in um.
  Length recall = truth length within tolerance of our strands / truth length; length precision = our strand length
  within tolerance of truth / our strand length in the territory.
- One-to-one matching within each nucleus: Hungarian assignment between truth ribbons and our strands of at least
  1.0 um, cost = mean symmetric centreline distance (both resampled at 0.05 um), admissible at 0.5 um or less. A
  matched pair is **matchable** when our strand covers at least 0.8 of the ribbon and at least 0.8 of the strand
  lies on the ribbon (within tolerance).

### 4.4 Synthetic matching

As 4.3, with the exact truth polylines. Gaps: a detected gap matches a truth gap when their s intervals, projected
onto the truth ribbon, overlap; one-to-one by largest overlap. Foci: the synthetic spot table is the truth foci
moved by Gaussian localisation noise (sigma 0.03 um in xy, 0.06 um in z); a focus is scored on the strand matched to
its truth ribbon.

### 4.5 Checks before scoring (truth side and run side)

1. Every test zone export, the N2 master table and the R script match the sha256 in sc_split.json.
2. Zone letters exist among the point names, satisfy f <= u <= l <= e, and the names follow A, B, ... in row order.
3. Every converted polyline lies inside its image's xy extent within 1 um.
4. The live file shows 10 ribbons of 2 or more vertices and 3 single-vertex stubs at time index 60.
5. Every test image has a finished germquant run at xy stride 1 whose central_element channel is the session's SC
   channel (477; 545 for 2025-10-21, through a session channel map kept outside the repo, built and checked on the
   two 2025-10-21 development gonads first).

A failure of 1 to 4 is a truth defect (rule 4, changelog first). A failure of 5 is fixed by rerunning germquant
before scoring; a gonad that cannot be run gets a "no output" row (section 10).

## 5. Methods and comparators

### 5.1 The method under test

sc_straighten v1 at the frozen commit with the frozen resolved `sc_straighten` section (params_hash) and the
frozen estimator values of 5.2, all written into the pre-scoring changelog entry (section 9).

### 5.2 The zone boundary estimator (pre-specified)

- Per-nucleus calls: full = `piece_lw_mean_um` >= theta_full; present = `total_length_um` >= theta_present (columns
  of `<id>__sc_strand_nuclei.csv`, DESIGN 5.7).
- Over the analysis nuclei of a gonad (guard bands included), sorted by s: f(s_j) = the mean call over nuclei with
  |s_i - s_j| <= 5 um, defined at each nucleus position s_j with at least 3 nuclei in the window.
- SC_fullstart estimate: the smallest s_j with f_full >= 0.5 there and at every defined position in
  [s_j, s_j + 20 um], with at least 3 defined positions in that range; none: not found.
- SC_fragstart estimate: the same with f_present.
- SC_fraglate estimate (only when SC_fullstart was found): the smallest s_j > (SC_fullstart estimate + 20 um) with
  f_full < 0.5 there and at every defined position in [s_j, min(s_j + 20 um, s_end)], with at least 3 defined
  positions; none: no late fragmentation (s_end).
- theta_full in {1.0, 1.5, ..., 6.0} um and theta_present in {0.5, 1.0, ..., 10.0} um are chosen on development
  gonads only, each to maximise the mean per-gonad balanced accuracy of its call (outside guard bands; full against
  SCfull versus SCfrag and LPfrag; present against any SC zone versus None); ties go to the smaller value.

### 5.3 Comparators (reported, never gating)

- Intensity baseline: `central_element_mean_intensity` of the run's nuclei.csv (raw mean SYP over the label), in
  place of our score in S2a and S2b, with the paired difference.
- Each other nucleus-tier column of DESIGN 5.7 (n_pieces, gaps_per_10um, longest_piece_um, on_strand_frac,
  profile_cv_corr) in S2a, labelled exploratory.
- Human repeatability: on the DLW190 gonads traced in January and again in April (35 with both), the difference of
  each boundary's normalised position s(p_k) / s(p_e) between the two tracings. Computed by the orchestrator at
  scoring time; a reference for S2c, not a criterion.

## 6. Metrics

### 6.1 S1: synthetic (sealed seeds at scoring; development seeds during development)

Families: **iso** 200 isolated ribbons, length U(2, 8) um, peak U(300, 800) counts, tilt to the xy plane U(0, 60)
deg for 160 of them and U(60, 75) deg for 40 (reported apart); **dense** 100 nuclei (50 oocyte, 50 spermatocyte);
**gaps** 100 isolated ribbons with 1 or 2 gaps of U(0.3, 1.0) um (no signal along the gap); **foci** 100 isolated
ribbons with 1 to 3 foci at known s (Gaussian spots of peak 300 counts in a foci channel); **pairs** 25 each of an X
crossing (angle U(40, 90) deg), a T contact, a parallel pair 0.5 um apart in xy and a pair stacked 1.0 um apart in
z; **noise** 50 nuclei with nucleoplasm only; **snr** 100 dense nuclei, 25 at each peak of 200, 300, 430 and 700
counts.

- S1a (iso, tilt up to 60 deg): length error (ours minus truth) mean and 95 % limits of agreement (mean +/- 1.96 SD);
  median centreline distance; share of ribbons with exactly one matched strand of 1 um or more.
- S1b (dense): per-nucleus total length relative error (mean, SD) and CCC of totals; coverage (truth length within
  tolerance of our strands / truth length); off-truth length (our length not within tolerance / truth length);
  strands per nucleus and the matchable fraction (descriptive).
- S1c (gaps of 0.3 um or more): precision and recall of gap detection; mean absolute gap length error of matched
  gaps.
- S1d (foci): |s error| of assigned foci (median and 90th percentile); share assigned to the correct ribbon.
- S1e (all families): `complete` precision, the share of `complete` strands that are matchable with |length error|
  <= 0.4 um; `complete` recall (descriptive).
- S1f (noise): share of nuclei with a total under 1.0 um.
- S1g (pairs, snr; descriptive): correct continuation at crossings, `contact` ends at T contacts, merging of close
  pairs; the recovery curve; **the SNR floor**: the lowest scene peak whose dense mean total error is within +/-5 %
  (used by S4).

### 6.2 S2: zones (the 12 test gonads)

- **S2a (primary).** Per gonad, AUC of `piece_lw_mean_um` for SCfull nuclei against fragmented nuclei (SCfrag and
  LPfrag), outside guard bands; a gonad is evaluable with at least 10 nuclei in each class. Statistic: the mean AUC
  over evaluable test gonads. Ties count one half.
- **S2b.** Per gonad, AUC of `total_length_um` for nuclei in any SC zone against None nuclei (at least 10 each); the
  mean over evaluable gonads.
- **S2c.** SC_fullstart: a hit when the estimate lies in [s(p_(u-1)), s(p_(u+1))]; the hit rate over all test gonads (not
  found counts as a miss); the median |estimate - s(p_u)| over found estimates, in um and as a fraction of s(p_e).
- Secondary (reported): S2a split into SCfull against SCfrag and SCfull against LPfrag; the comparators of 5.3 with
  paired differences; SC_fragstart and SC_fraglate hit rates with the same interval rule; the SCfull zone length,
  (s(p_l) - s(p_u)) / s(p_e), against ours, as bias and CCC over gonads; breakdowns per condition x sex cell and per
  session; nuclei left out per reason; the human repeatability of 5.3.

### 6.3 S3: live Filament check (descriptive)

- S3a: each truth ribbon measured as drawn against its stored Imaris `Segment Length` (relative difference).
- S3b: per nucleus, length recall and precision of our strands, and our nucleus total against the refined truth
  total (relative error).
- S3c: the number of matchable ribbons; their length error (ours minus refined truth) and, beside it, against the
  as-drawn truth.
- S3d: refined truth against as-drawn truth per ribbon (the end rule and the voxel-centre polyline).

### 6.4 S4: fixed-image hand traces (registered now, scored when T4 exists)

Per-nucleus total length (nuclei whose peak is above the SNR floor of S1g): bias in percent, CCC, Bland-Altman.
`complete` strands: length bias, 95 % limits of agreement, CCC. Decomposition: `complete` precision and recall, the
matchable fraction. Centreline distance and end error of matchable pairs. Gaps, if the traces mark breaks as
separate filaments. Matching as in 4.3 with the hand traces in place of the live truth, inside our nucleus labels.

### 6.5 Reporting

Every table has one row per test gonad and a pooled row; every test gonad appears in every table, with "no output"
where a method produced nothing. Each result carries the commit, params_hash, split sha256 and generator commit,
under `C:/Users/ryane/germquant_golden/benchmarks/sc_ribbon_v1/final/`.

## 7. Statistics

- Zone metrics: cluster bootstrap by gonad (`germquant.benchmark.stats.cluster_bootstrap`), 1000 resamples, seed
  20261023, percentiles 2.5 and 97.5; for a mean of per-gonad AUCs the evaluable gonads are resampled. Paired
  comparisons use the same resamples and the interval of the difference.
- Synthetic metrics: bootstrap over scenes, 1000 resamples, seed 20261023; limits of agreement by the normal
  approximation with their bootstrap interval. The synthetic criteria of 8.1 use the point estimates; the
  intervals are reported.
- S3: values only (2 nuclei, 10 ribbons); no interval is claimed.

## 8. Success criteria (pre-stated) and decisions

### 8.1 Criteria

| id | criterion |
|---|---|
| S1a | mean length error within +/-0.10 um; 95 % limits within +/-0.40 um; median centreline distance <= 0.05 um; one strand per ribbon in >= 95 % |
| S1b | mean total relative error within +/-5 %, SD <= 7 %; CCC >= 0.80; coverage >= 0.95; off-truth <= 0.10 |
| S1c | gap precision >= 0.80 and recall >= 0.80; gap length error <= 0.15 um |
| S1d | 90th percentile of the absolute s error <= 0.15 um; correct ribbon >= 0.95 |
| S1e | `complete` precision >= 0.90 |
| S1f | >= 95 % of noise-only nuclei with a total under 1.0 um |
| S2a (primary) | mean per-gonad AUC >= 0.80 and the lower 95 % bound >= 0.70, over at least 8 evaluable test gonads (fewer: inconclusive) |
| S2b | mean per-gonad AUC >= 0.85 and the lower bound >= 0.75 |
| S2c | SC_fullstart hit rate >= 0.70 (9 of 12) and median error <= 20 um |
| S3a | every ribbon within 1 % (a failure is an importer defect: S3 stops and is reported, and a fix is a revision) |
| S3b | in both nuclei, recall >= 0.80 and total within +/-15 % |
| S3c | >= 5 of 10 ribbons matchable, median absolute length error <= 0.5 um |
| S4 (later) | total bias within +/-5 % and CCC >= 0.80 above the SNR floor; `complete` strands bias within +/-0.10 um, 95 % limits within +/-0.40 um; `complete` precision >= 0.90 |

### 8.2 Decisions

- **D1, nucleus-tier continuity validated:** S2a and S2b pass, and S1b, S1c and S1f pass. Then `piece_lw_mean_um`,
  `n_pieces`, `gaps_per_10um`, `total_length_um` (as SC presence) and `on_strand_frac` may be used to compare zones,
  conditions and sexes on fixed N2 images from this microscope; the docs say so, and a post-scoring commit that
  changes no number replaces the image flag `sc_straighten:uncalibrated` by `sc_straighten:zones_validated`.
- **D2, zone estimator:** D1 and S2c pass. Then the per-gonad estimator may become a production output (automated
  zone lengths) in a later version, labelled with the measured error.
- **D3, lengths:** absolute per-ribbon and per-nucleus lengths are reported as "validated on synthetic data" when
  S1a and S1b pass, and as validated on fixed images only after S4 passes. Per-chromosome statistics use only
  `complete` or `verified` strands in every case.
- **D4, in-pipeline stage:** the phase-2 registry stage (DESIGN 1.3) is built only after D1.
- **D5, failure:** if S2a fails, the continuity metrics are not recommended for zone analysis. At most one revision
  may be declared and scored (section 10); after that, a new test set is needed (the reserved gonads, by a
  changelog entry made before any sc_straighten output of them exists).
- S3 never changes D1 to D4; it only says whether the live modality is worth a later version.

## 9. Development, tuning budget and freeze

- Tuning on real data: at most 40 configurations evaluated on development zone metrics, each a line of
  `benchmarks/sc_ribbon_v1/dev/dev_trials.jsonl` (params_hash, the changed keys, development S2a, S2b, S2c and the
  synthetic development gates). Stop when the development S2a statistic has not improved by more than 0.01 over
  the last 10 configurations. Synthetic sweeps on development seeds are not limited.
- Readiness before the freeze: CI green; the development synthetic gates of 8.1 reported; development S2 computed
  on every development gonad that has a germquant run (at least the eight of 2025-11-05); theta_full and
  theta_present chosen by 5.2; a robustness report: the change of the development S2a statistic when each of
  follow.seed_min, follow.lo, link.max_gap_um and gaps.q moves by +/-20 % (reported, not a criterion).
- Freeze: a pre-scoring changelog entry records the commit, params_hash and the resolved section, theta_full,
  theta_present, the split record sha256, and the germquant run folder and config hash of every test image. A
  method still being changed is not scored.

## 10. Scoring (the orchestrator, once, later)

- Who and when: the orchestrator, after the freeze entry is committed and pushed, outside the build workflow, and
  only on the project lead's go. Not in the workflow that builds and tunes the method.
- Steps:
  1. Pre-flight: the checks of 4.5. Defects go to the changelog before anything else runs.
  2. `scripts/benchmark_sc.py zones --final --reason TEXT`: sc_straighten on the 12 test images with the frozen
     parameters (`--allow-sealed`), then S2 and the comparators. `filaments --final`: S3. `synthetic --final`: S1 on
     the sealed seeds. Each run appends a line to `benchmarks/sc_ribbon_v1/final_runs.jsonl`.
  3. The report under `benchmarks/sc_ribbon_v1/final/`.
  4. A post-scoring changelog entry: the results, the decisions of 8.2, and anything unusual (labelled post hoc
     where it is).
- The test set is scored at most twice per method family: once with the pre-registered method and at most once
  more after a revision declared in the changelog, with its reason, before the second scoring. Both results are
  reported wherever either is. A third attempt needs a new test set.
- A method that produces no output for a test image gets a "no output" row in every table; criteria are evaluated
  over the images with output and labelled with their number (11 of 12, say); a pass on fewer than all images is
  provisional until the missing one is processed or its failure explained in the changelog.
- After the first scoring a changelog entry may set `"sealed": false` in the split record for the test gonads, so
  the lab can process them; the reserved gonads stay sealed.

## 11. Requests for new truth, and how it enters

- **R1: the three empty Filament files** (live DLW193; 1, 5 and 4 traced nuclei by name): ask Madeleine to re-save
  them with the Filaments committed (check that the Filaments object is not empty before closing Imaris). When they
  arrive they are development ribbon truth for the live modality, registered by a changelog entry with their sha256
  before any comparison.
- **R2: independent 3D hand traces on fixed images (T4, for S4).** Imaris FilamentTracer as in Cahoon et al. 2023,
  each nucleus traced before its gonad's sc_straighten output has been seen by the tracer and never by correcting
  our output. Per condition x sex cell: at least 10 nuclei on development gonads and 10 on test gonads (15 each
  preferred), at most 5 per gonad, mid-pachytene (SCfull zone), the whole nucleus inside the stack, every SC of the
  nucleus traceable (6 in oocytes, 5 in spermatocytes; otherwise the nucleus is skipped), spread along the zone. A
  second person traces 20 of these nuclei for the inter-annotator ceiling. The files and their sha256 enter by a
  changelog entry before any comparison; development traces then serve calibration, and test traces are scored at
  the next scoring.

## 12. Changelog (append only)

Entries are appended, never edited. Each one says whether a test result had been computed at that point.

1. 2026-10-23. Protocol written and split registered (`sc_split.json`, sha256
   db328c1cf4bc2b08dda32027a0e2567a128700a26facac4f7dc47738f67ae9eb; 102 gonads: N2 18 development, 12 test,
   9 reserved; IR 12 and CCW69 9 reserved; DLW190 42 excluded; the live Filament file test). No SC output of any
   method had been compared with any SC truth, and no test result existed. Looked at before writing, all of it
   recorded in the inventory, the code brief and the two designs: (a) n2rad51_noHS_herm_000, development: its SYP
   channel and untuned outputs of the sc_trace tracer (code brief), Design A's and Design B's probes (overlays of
   about 120 nuclei, untuned statistics); (b) 20251021_n2_nohs_HERM_2, development: example rows of its export
   (Overall, first rows of Position, Distance, Angle and Intensity Center), no zone letter; (c) every truth file:
   sheet names, headers, row counts, point counts and names, point bounding boxes for image matching, the master
   tables' metadata columns (zone letters read in code only, never printed); (d) the live Filament file: object
   tables, vertex counts, the voxel bounding boxes of its two nuclei and the statistic names (no length or
   intensity value); (e) the R scripts as code; (f) the per-gonad metadata of the inventory. Seals: none of the four
   RAD-51 test `.ims` files was opened, and no RAD-51 xlsx export, Imaris spot truth or lab aggregate was opened.
2. 2026-10-23, the build workflow (branch feat/sc-ribbon-quant), before any result: no test result had been
   computed, no test or reserved image was run, no test or reserved truth was read and no sealed seed was
   generated. Nothing above the changelog was edited. Implementation notes, recorded so the scoring commit can be
   read against this protocol:
   (a) The method is the package `germquant.sc.straighten` (command `germquant sc-straighten`); the generator of
   section 2.3 is `germquant.sc.straighten.synth`.
   (b) Generator details that section 2.3 leaves open, as coded: the iso, gaps and foci ribbons are planar
   circular arcs of curvature radius U(1.2, 4.0) um (raised so that an arc spans at most 0.9 pi), the chord tilted
   by the stated tilt toward z, uniform azimuth, the bow in the image plane; a gap leaves at least 0.8 um of ribbon
   on either side and two gaps occur only on ribbons of 4 um or more; foci lie at least 0.4 um from the ends and
   0.6 um apart; the X crossing and the T contact lie in one image plane; a pairs scene has one enclosing nucleus;
   nucleus radii and labels as in 2.3; the image is rounded to uint16.
   (c) Motivated by development results (synthetic development seeds): in the dense and snr families a ribbon
   that comes closer than 0.25 um in 3D to a ribbon already placed in its nucleus is drawn again (at most 60
   draws; then the draw with the largest clearance): two chromosomes' SCs do not pass through each other. Without
   this rule 13.5 percent of the dense truth length lay within half a PSF FWHM of another ribbon.
   (d) Motivated by development results (synthetic development seeds, then development image
   n2rad51_noHS_herm_000 for the multiplicity reference only; no zone truth had been compared): defaults changed
   from DESIGN section 8, tunable keys follow.seed_min 4.0 to 2.5, link.max_skew_um 0.3 to null, refine.knot_um
   0.3 to 0.5, refine.lambda 1.0 to 0.1, gaps.q 0.35 to 0.6; and, among keys the design fixes, follow.momentum 0.7
   to 0.5, link.dir_window_um measured from the tip [0.3, 1.2], new keys follow.max_turn_deg 30 (pure prediction
   when the local line turns away by more than 30 degrees in one step, or inside another strand's claim tube),
   link.tip_um 0.3 and link.max_gap_contact_um 1.5 (end matching from anchors inside the drifting tips),
   dedupe.ghost_max_um 1.5, dedupe.ghost_tz_min 0.7 and dedupe.ghost_reach_um [1.0, 0.35] (axial shadows of ribbon
   ends and crossings removed), ends.max_extension_um 0.5, refine.lambda_imported 0.01, gaps.window_frac 0.9 and
   gaps.ring_um 0.45 (a gap's length is the missing signal above a local background over the dip, divided by the
   median), and the multiplicity reference is the nucleus's own median strand signal. The QC code low_snr has a
   separate fail-level code low_snr_fail. docs/SC_STRAIGHTEN.md gives the reasons. These are development values;
   the freeze entry of section 9 records the values that are scored.
3. 2026-10-23, the build workflow, before any result (no test result, no test or reserved image or truth read,
   no sealed seed). Development trial 1 (commit 8d23105, the defaults of entry 2) was scored on the eight
   2025-11-05 development gonads (dev_trials.jsonl). Then, motivated by synthetic development results only (gap
   recall 0.68 and ribbons broken at gaps): link.unlinked_cost 40 to 60 (tunable), tips that meet or overlap are
   matched by their collinearity, and accounting drops strands that lie for at least dedupe.shadow_frac (0.8, a new
   fixed-group key) of their length inside the contamination reach of a longer strand. Development trial 2 uses
   these values.
4. 2026-10-23, the build workflow, before any result. Development trial 2 (commit f7dc9a3) scored. Then, motivated
   by synthetic development seeds 2000 to 2049 (isolated ribbons over-measured by 1 to 3 um where the trace ran
   along the axial PSF tail into the territory edge, or joined an axial shadow): axial shadows are also removed
   before the end matching, and a trace that reaches the territory or stack edge after its ridge has faded ends
   free (trimmed to its last strong point). Development trial 3 uses these values.
5. 2026-10-24, the review round of the build workflow, before any result: no test result had been computed, no
   test or reserved image was run, no test or reserved truth was read and no sealed seed was generated.
   Development trial 3 (commit 3e88403) had been scored on the eight 2025-11-05 development gonads. A correctness
   review and a protocol review of the branch followed; their findings were checked on development seeds and
   development images, and these changes were made before development trial 4:
   (a) Method, motivated by synthetic development seeds and by the trial-3 runs of the eight development images
   (strand geometry and end types only; no zone truth motivated any of them): T contacts are resolved after the
   axial shadows are removed, and a `contact` end with no surviving strand within `follow.contact_free_um` (a new
   fixed-group key, [0.6, 0.4] um axial and lateral) is free again (2161 of the 6357 contact ends of the trial-3
   runs had no partner); linked pieces whose tips overlap are joined without a backward step (835 of 23,419
   trial-3 strands had a folded frame, 788 of them linked); an end fit whose edge sits at the inner bound of its
   window is fitted again from an anchor `ends.inside_um` deeper, and the first fit stands when that one fails;
   gaps: a candidate is a run below `gaps.window_frac` of the strand median, its window is widened by
   `gaps.tail_sigmas` (new, 1.0) PSF sigmas along the strand, the dip is judged against a straight baseline
   between the profile medians within `gaps.plateau_um` (new, 0.4 um) outside the window, and the length is the
   missing fraction of that baseline (tilted gaps measured up to 30 percent short with the strand median as the
   reference); foci: the closest strand and the ambiguity test use the distance to the strand proper, so a strand
   that ends short of a focus and only points at it neither claims it nor makes it ambiguous; QC:
   `label_merged_suspect` compares a label with the median volume of the traced nuclei within
   `qc.merged_label_window_um` (new, 15 um) along the gonad axis, `noise_fit_poor` judges the fit over the noise
   bins above the darkest decile of the sampled voxels, and the input fallback to a parent `golden_manifest.json`
   for the run geometry has the image flag `geometry_from_golden_manifest`. Keys introduced in earlier rounds
   without a changelog line are recorded here: `follow.contact_um` [0.3, 0.15], `noise.var_floor` 1.0 and
   `refine.se_curve_um` 0.10 (a constant in place of the calibration table of DESIGN 4.10).
   (b) Checked and not adopted: the outward end direction from a line fitted 0.6 to 1.5 um inside the end (the
   length bias of isolated ribbons tilted 30 to 60 degrees changed by 0.01 um on development seeds); a straight or
   Hermite bridge between the flanks of a gap as the path of the gap measurement (tilted gaps came out shorter,
   and along the true centreline they were short as well: the cause was the reference level, fixed in (a)).
   (c) Scoring path (code only; no rule of this protocol changed): `zones --final` uses the frozen theta_full and
   theta_present of the freeze entry from the command line (before, it would have chosen them on the scored test
   gonads), each gonad's registered frame rule (the eight test gonads registered with the .ims header rule would
   have used the nd2 stage rule), the frame agreement on the image corners with the voxel included, and runs
   sc-straighten itself on the test images; `--allow-sealed` needs a reason and logs a start and an end line per
   sealed image; every --final run logs a start line before it touches sealed data and an end line after it; the
   pre-flight checks 1, 3 and 5 of section 4.5 are coded; every gonad in scope gets a row, with "no output" and
   its reason where it could not be scored, plus a pooled row; sealed synthetic scoring has a path
   (`calibrate_sc_straighten.py --final --reason`), refused without a reason; image stems are compared with the
   sealed stems without regard to case, and imported geometry from the registered live Filament file is refused
   without `--allow-sealed`.
   (d) Reporting added under rule 2 (secondary, never gating): S1a is also reported in the tilt bins 0 to 30, 30
   to 60 and 60 to 75 degrees and S1e per family; at the sealed scoring S1b and S1e are also reported on the
   generator as registered in section 2.3 (dense and snr ribbons without the 0.25 um clearance of entry 2(c));
   the development S1 gates are also reported on the fresh development seeds 3000 to 3049.
   (e) Development data added: the four 2026-01-05 heat-shock male development gonads have germquant runs
   (germquant_golden/additions/8666f07), so development trial 4 scores them with the eight 2025-11-05 gonads
   (section 9: every development gonad with a germquant run).
   (f) Open, for the project lead, before the freeze: the guard bands of section 4.2 (half the neighbouring point
   spacing) remove every nucleus of a zone one segment long, so on the eight 2025-11-05 development gonads SCfrag
   kept no nucleus outside the guard bands in 6 of 8 and S2a was evaluable on 3; a guard of a quarter spacing
   would make all 8 evaluable (nucleus counts only; no AUC has been computed with another guard). The scoring
   code takes `--guard-frac` (default the registered 0.5; `--final` refuses another value unless a changelog
   entry sets it). Also open: the population of S4 when no SNR floor exists (no scene peak met the 5 percent
   rule of S1g on development seeds).
   (g) The protocol file `C:/Users/ryane/sc_ribbon_work/VALIDATION_PROTOCOL.md` stays the registration snapshot
   (sha256 d56b68859b8a500f6aad964fbb3c8cb9d96437c3108fe957fb8760dd125a2a01); this repository copy is the living
   protocol, and its changelog is appended here.
6. 2026-10-24, the review round, before any result (no test result, no test or reserved image or truth read, no
   sealed seed). Development trial 4 (the runs of commit 6b7be22 on the twelve development gonads with germquant
   runs) was scored; its first scoring counted the six development gonads without a run as SC_fullstart misses, so
   it was scored again at commit cf43042 with the rates over the twelve gonads with output (both lines are in
   dev_trials.jsonl, the second marked as the rescore). Commit cf43042 also reads the .nd2 channels plane by plane,
   which changes memory only (identical outputs checked on two development gonads). No parameter was changed after
   trial 4. The development S1 gates were reported on development seeds 2000 to 2049, on the fresh seeds 3000 to
   3049 and, for S1b and S1e, on the registered generator (entry 5 (d)).
7. 2026-10-24, the v2 development round, written before any number below was computed (no test result had been
   computed, no test or reserved image was run, no test or reserved truth was read and no sealed seed was
   generated). An analysis of the continuity score on the zone truth of the twelve development gonads with
   trial-4 runs (the outputs of commit 6b7be22; nothing is run again), motivated by development results: on S2a
   the raw SYP intensity (mean AUC 0.97) separated SCfull from fragmented nuclei better than piece_lw_mean_um
   (0.83, per gonad 0.45 to 0.98, 7 evaluable gonads), and 39 to 46 percent of the traced length runs steeper than
   60 degrees. Looked at before this entry: the trial-4 zone tables (the per-gonad AUCs of the nucleus columns and
   of the intensity in zones_per_gonad.csv) and, without zone truth, the trial-4 strand tables: frac_steep60 is
   bimodal (28 percent of the 35,840 strands under 0.1, 35 percent above 0.9), steep strands are short and dim
   (median 1.5 um, median integrated signal 107 against 186), 83 percent of the strands carry no detected gap, and
   in the traced nuclei with 30 um of SC or more the median piece_lw_mean_um is 2.6 um (a full ribbon is 4.5 to
   7.0 um in the phantoms of 2.3). The rules below were fixed after that look and before any zone label was joined
   to them. Scripts and logs: C:/Users/ryane/sc_ribbon_work/v2/ (scripts/cont_*.py, logs/cont_*).
   (a) Candidates, computed from the trial-4 strands, samples and gaps tables (pieces: each strand split at the s
   windows of its gaps, as in v1; this reconstruction reproduces piece_lw_mean_um within 2e-5 um in every
   nucleus). A strand is steep when frac_steep60 > 0.5 (the steep_z rule, qc.steep_frac_review), in-plane
   otherwise. C1: piece_lw_mean_um over the pieces of in-plane strands only (0 without one). C2: on_strand_frac
   times the share of the strands' integrated signal (in-strand samples outside gaps, clipped at 0) that lies on
   pieces of 3.0 um or more of in-plane strands (0 without strands): the share of the nucleus's SYP on long
   in-plane SC; 3.0 um is half a full-length ribbon and is not tuned. C3: the break density, pieces of 0.25 um or
   more of in-plane strands per um of their summed length, scored as its negative (a nucleus with under 1 um of
   such pieces counts as the most fragmented); pieces are counted rather than detected gaps because most breaks
   end a strand. R0 is piece_lw_mean_um.
   (b) Reported for R0, C1, C2, C3 and the intensity comparator (central_element_mean_intensity): the per-gonad
   S2a AUC with the registered rules (4.2, 6.2), the mean, SD, minimum and maximum over the evaluable gonads, the
   cluster-bootstrap interval of the mean and of the paired differences against R0 and against the comparator
   (section 7: 1000 resamples, seed 20261023); SCfull against SCfrag and against LPfrag; all of it again with
   guard_frac 0.25 and 0 (for decision 1 of the project lead; the registered 0.5 stays unless a later entry
   changes it).
   (c) Diagnostics (they explain and never choose): per gonad and zone, the median intensity above the camera
   offset of 100 counts, the steep share of the traced length, the in-plane and steep strand counts against the
   expected count, and the pieces per um; within SCfull and within the fragmented zones, the Spearman correlation
   of each score with the intensity; the S2a AUC of each score's within-gonad rank after a least-squares fit on
   the rank of the intensity over the S2a nuclei (classes pooled, the class not used), and of the intensity after
   the same fit on each score (what each carries beyond the other); the S2a AUC of each score averaged over the
   nuclei within 5 um along the gonad axis (the window of the 5.2 estimator; a measure of per-nucleus noise, never
   an S2a statistic).
   (d) The recommendation rule, fixed now: a candidate may replace piece_lw_mean_um as the S2a statistic only when
   its mean AUC under the registered guard exceeds R0's by at least 0.03 with a paired bootstrap interval above 0,
   it is below R0 in at most two evaluable gonads, and its paired difference keeps its sign with guard_frac 0.25;
   between qualifying candidates the one closest to v1 (C1, then C3, then C2) wins unless another is ahead by more
   than 0.03 (C1 keeps the unit of theta_full and the grid of 5.2). The choice is in-sample (7 evaluable gonads
   under the registered guard): it is checked on the six development gonads without runs once they have runs,
   before the freeze, and no candidate is added after that check. Each candidate counts as one configuration
   against the budget of section 9.
8. 2026-10-24, the v2 development round (steep strands), written before any number of the candidate rules below
   was computed (no test result had been computed, no test or reserved image was run, no test or reserved truth was
   read and no sealed seed was generated). Motivated by development results: on the twelve development gonads with
   trial-4 runs, 46 percent of the traced length runs steeper than 60 degrees, against 13.4 percent for randomly
   oriented ribbons. Looked at before this entry (scripts and logs: C:/Users/ryane/sc_ribbon_work/v2/,
   scripts/steep_*.py and logs/steep_*; outputs under v2/steep/): (i) from the trial-4 tables and, for the zone of
   each nucleus, the zone truth of the development gonads through the guarded loader: the steep share of the traced
   length per gonad (no heat shock 26 to 33 percent, heat shock 44 to 62 percent) and per zone (None 74 percent,
   SCfrag 58, SCfull 35, LPfrag 38), its |t_z| histogram and the steep strand lengths; (ii) 3D Gaussian fits on the
   raw counts of the development images: the RAD-51 foci of the runs' spot tables (foci channel), the brightest SYP
   maxima of each nucleus and the object under each strand's brightest sample, and a Gaussian fitted to each
   strand's core profile along s; (iii) development phantoms that are not registered families (seeds 5000 to 5459,
   the v1 code of commit 5cb3cb6): nuclei with PSF-blurred point puncta and with spheres of 0.25 to 0.8 um, dense
   nuclei with 6 point puncta at least 0.6 um from the ribbons, two flat arcs crossing at height differences of 0.2
   to 1.0 um, and isolated arcs of 1.5 to 6 um tilted 60 to 90 degrees. v1 makes a strand of 85 percent of the
   isolated point puncta, 99 percent of them steep (half-maximum length median 1.0 um, in-plane extent 0.07 um),
   and the Gaussian fitted to their core profile has sigma 0.385 um (95th percentile 0.55); every tilted arc gave
   0.83 um or more; the crossings gave no steep strand.
   (a) Candidate rules, fixed now and evaluated as development diagnostics only (no package code changes in this
   entry; adopting a rule needs a later entry). A strand is steep when frac_steep60 > 0.5. B1 (blob rule, after
   measurement): a steep strand with length_xy_um <= 0.5 um whose core profile over s from -0.5 to length_um +
   0.5 um is fitted by least squares with a Gaussian plus a constant of sigma times tz_abs_mean <= 1.5 times the
   image's axial PSF sigma (psf.sigma_z of the record) is a blob: it leaves the strand set and the nucleus totals
   (it would be reported as a punctum). B1 is also reported with 1.3 and 1.8 times sigma_z and with 0.3 and 0.75 um.
   B2 (minimum in-plane extent): a steep strand with length_xy_um under 0.3 um (and, separately, 0.5 um) is
   removed. B3 (the PSF extent law): a steep strand with length_xy_um <= 0.5 um whose traced z extent is at most
   2 c sqrt(2 ln(S / follow.lo)) + 0.3 um, with c = sqrt(sigma_z^2 + fields.sigma_um[0]^2) and S its largest ridge
   strength, is removed. B4 (puncta mask): a strand is removed when the 3D Gaussian fitted to the raw counts around
   its brightest sample (window +/-1.2 um axial, +/-0.54 um lateral) is compact (lateral sigmas at most 0.25 um,
   their ratio at most 1.5, axial sigma at most 1.5 times sigma_z, relative residual RMS at most 0.2) and 80 percent
   of the strand lies inside that object's ellipsoid (3 axial sigmas, 0.3 um lateral). B5 (anisotropic ridge scale,
   a tracer change; phantoms and two development images only): fields.sigma_um [0.2, 0.17, 0.17], which makes a
   PSF-blurred point isotropic in voxel units. B6 (deconvolution before tracing; phantoms only): 20
   Richardson-Lucy iterations with the Gaussian PSF of section 2.3 on the image before v1.
   (b) Reported for each rule, choosing nothing that is scored: on the phantom families of (iii), at seeds 5000 to
   5459 and on the fresh seeds 5500 to 5959, the share of punctum strands removed and the share of ribbon strands
   and ribbon length removed (tilted arcs by tilt bin, dense nuclei, crossings); on the twelve development gonads,
   per gonad and zone, the share of strands and of length removed, the steep share of the remaining traced length
   against 13.4 percent, and the strands per traced nucleus against the expected count. No S1 gate and no S2
   metric is computed in this entry; a rule adopted later counts as a configuration under section 9 and gets its
   own entry before its zone metrics are computed.
9. 2026-10-24, the v2 development round (steep strands), written before the number below was computed (no test
   result had been computed, no test or reserved image was run, no test or reserved truth was read and no sealed
   seed was generated). The rules of entry 8 were evaluated as registered; B5 was run on phantoms only (its phantom
   result already ruled it out, so no development image was run with it). Development diagnostics added after
   entry 8, none of which changes the method: a family of short tilted arcs (0.6 to 1.5 um, 60 to 90 degrees,
   seeds 6000 to 6119), dense phantoms at scene peaks of 60 to 430 counts (seeds 6100 to 6139), averaged axial
   profiles of compact RAD-51 foci and SYP maxima, the depth dependence of their axial sigma, and the orientation of
   elongated SYP objects fitted in the raw images (full-covariance Gaussians in an isotropic window, calibrated on
   dense phantoms, seeds 6300 to 6359). Motivated by those development results: the image's PSF record fell back to
   sigma_z 0.365 um on 10 of the 12 images, while the compact RAD-51 foci of each session, scaled from 595 to
   511 nm, give 0.380 um (2025-11-05 no heat shock), 0.428 um (2025-11-05 heat shock) and 0.464 um (2026-01-05).
   B1s, fixed now and reported like B1 (development diagnostics only, choosing nothing scored): B1 with sigma_z the
   session value above in place of the record's, at 1.3 and 1.5 times, with length_xy_um <= 0.5 um.
10. 2026-10-24, the v2 development round, the results of entry 7 (no test result had been computed, no test or
   reserved image was run, no test or reserved truth was read and no sealed seed was generated). Computed as entry 7
   fixed it (C:/Users/ryane/sc_ribbon_work/v2/logs/cont_zone_analysis.log and .json): S2a under the registered
   guard, 7 evaluable gonads, R0 0.828 (0.671 to 0.933), C1 0.802, C2 0.745, C3 0.776, intensity 0.965; paired
   differences against R0 C1 -0.026 (-0.058 to -0.000), C2 -0.082 (-0.136 to -0.035), C3 -0.052 (-0.091 to -0.012);
   with guard_frac 0.25, 12 evaluable, R0 0.802 (0.728 to 0.859), C1 0.789, C2 0.779, C3 0.752, intensity 0.910;
   with 0, R0 0.776 and intensity 0.874. No candidate met rule 7 (d): piece_lw_mean_um stays the S2a statistic and
   no method or metric changed; the candidates are three lines of dev_trials.jsonl. The AUCs with guard_frac 0.25
   and 0 were seen on development gonads before decision 1. Explanatory checks made after entry 7, choosing nothing:
   axis profiles and strand end classes of the development gonads; SYP projections of random nuclei of
   20251105_N2_nohs_MALE_1 and 20251105_N2_nohs_MALE_3; development phantoms that are not registered families (seeds
   2000 to 2049, scene peaks 200 to 700 counts, the dense nucleus of each seed with continuous ribbons and the same
   nucleus with every ribbon cut into pieces of 1 to 2 um by gaps of 0.4 to 0.8 um), where R0 separated continuous
   from cut ribbons with AUC 0.99 to 1.00 at every peak and the mean intensity over the label with 0.55 to 0.61.
11. 2026-10-24, the v2 development round (its rules), written before any number of the v2 changes was computed (no
   test result had been computed, no test or reserved image was run, no test or reserved truth was read and no
   sealed seed was generated). Motivated by development results: entries 7 to 10, and a failure taxonomy of the
   dense phantoms made without a method change and without zone truth (the code of commit 5cb3cb6 with hooks on the
   phase-A steps, development seeds 2000 to 2049 and 3000 to 3049, both generators;
   C:/Users/ryane/sc_ribbon_work/v2/dense_taxonomy/report_v1.json). The round follows
   C:/Users/ryane/sc_ribbon_work/v2/PLAN_V2.md (sha256 in v2/logs/prereg_entry11.txt). Fixed now:
   (a) Seeds. Development seeds used before this round: 1000 to 1023, 2000 to 2049, 3000 to 3049, 5000 to 5459, 5500
   to 5959, 6000 to 6139 and 6300 to 6359. Changes of this round are designed on 2000 to 2049 and 3000 to 3049,
   tuned on 10000 to 10099 (any number of runs) and confirmed on 11000 to 11099 (once per change and per revision);
   the development snapshot of the configuration proposed for the freeze is computed once, on 12000 to 12099, which
   no run touches before. Development families that are not registered use the same blocks; they live outside the
   generator of 2.3, whose scenes stay identical for every seed (a test pins them).
   (b) A change is adopted on synthetic evidence: its targeted metric, named in its entry, improves on the tuning
   seeds by the amount the entry states and moves the same way on the confirmation seeds by at least half of that
   amount, and on the confirmation seeds no criterion metric of 8.1 gets worse than before the change, on the same
   scenes, by more than: S1a mean length error 0.03 um, 95 percent limits 0.05 um, median centreline distance 0.01
   um, one strand per ribbon 0.02; S1b mean total error 1.5 percentage points, SD 0.5 points, CCC 0.03, coverage
   0.01, off-truth 0.01; S1c precision and recall 0.03, gap length error 0.03 um; S1d 90th percentile 0.02 um,
   correct ribbon 0.02; S1e precision 0.02; S1f 0.02; the median time per nucleus on n2rad51_noHS_herm_000 stays
   within 1.5 times that of v1 and no process exceeds 8 GB. A parameter that a change introduces is set from physics
   or the literature where its entry says so, otherwise chosen on the tuning seeds by the rule its entry states;
   none is chosen on development images. Every change has a switch whose other value gives the configuration before
   it; a change that fails is withdrawn by its switch, and a change gets at most two revisions, each with its own
   entry, before it is withdrawn. The criteria of 8.1 themselves are judged on the snapshot seeds.
   (c) Zone truth enters this round at most three times: Z1 after the per-image PSF and the punctum class; Z2 after
   the decomposition changes (contact stop, end matching, junction end fit, complete guard and, if made, junction
   re-pairing); Z3 after a realism-driven change, if one is made. A checkpoint runs sc-straighten with the
   cumulative configuration on the twelve development gonads with trial-4 runs and computes S2 by the registered
   rules (theta_full and theta_present chosen again by 5.2) under guard_frac 0.5 and 0.25, with the paired
   differences against development trial 4 (cluster bootstrap, section 7); it is one line of dev_trials.jsonl. A
   checkpoint can only veto: when the paired S2a difference against trial 4 is below -0.03 under both guards, or the
   S2b mean is more than 0.01 below that of trial 4, the changes it covers are withdrawn and that configuration is
   not scored on zone truth again; at Z1 the per-image PSF and the punctum class and flag stay while punctum strands
   return to every total, count and continuity column as in v1; at Z2 the complete guard stays (no S2 input depends
   on the strand class). No two variants are compared on zone truth and no parameter is chosen on it, the registered
   grid of 5.2 excepted. Statistics by zone are computed only at the checkpoints.
   (d) Truth-free statistics of the development images (flag rates, strands per nucleus, steep shares, end classes,
   profile statistics, PSF estimates) may be computed at any time; they may calibrate development phantoms (the
   realism entry of this round) and set no method parameter directly.
   (e) Hold-out. The six development gonads without a run (20251021_n2_nohs_HERM_2, 20251021_n2_nohs_MALE_2,
   20260105_n2_hs_1_HERM_1, 20260105_n2_hs_HERM_2, 20260105_n2_nohs_HERM_1, 20260105_n2_nohs_MALE_3) get germquant
   runs now, the 2025-10-21 session channel map being built and checked on the first two, and no sc-straighten run
   until the configuration proposed for the freeze is fixed. Then the trial-4 configuration and that configuration
   run on them once; S2 is computed with theta_full and theta_present of the twelve (so S2c is out of sample there)
   and with the paired differences; the result is one line of dev_trials.jsonl. When the paired S2a difference over
   the gonads evaluable under guard_frac 0.25 is below -0.05, the Z1 configuration is proposed for the freeze
   instead; with fewer than three evaluable gonads the check is reported and decides nothing. Nothing is changed
   after this check.
   (f) Not changed in this round: the QC flag thresholds, the S1 criteria of 8.1 and the S2 rules; the definition of
   piece_lw_mean_um, except that punctum strands are not pieces; the registered synthetic families and
   sc_split.json.
   (g) Reporting added under rule 2 (secondary, never gating), at every checkpoint and at the scoring: the S2a AUC
   of piece_lw_mean_um's within-gonad rank after a least-squares fit on the rank of the intensity (entry 7 (c));
   n_puncta, puncta_signal_frac and steep_length_frac as exploratory comparators (5.3); and docs/SC_STRAIGHTEN.md
   states that the intensity comparator measures the amount of SYP, not its continuity (entry 10).
   (h) The package VERSION becomes 2.0. The method under test of 5.1 is the configuration named by the freeze entry;
   no test result exists for any version, so its scoring is the first of section 10.
12. 2026-10-24, the v2 development round (per-image PSF), written before any number of this change was computed (no
   test result had been computed, no test or reserved image was run, no test or reserved truth was read and no
   sealed seed was generated). Motivated by development results (entry 9 and the trial-4 records): the PSF of 10 of
   the 12 development images fell back to sigma_z 0.365 um (0 to 15 flat isolated bright stretches against
   psf.min_stretches 20; 83 and 97 stretches on the two images with an estimate, 0.3675 and 0.3825 um), while the
   compact RAD-51 foci of the runs' spot tables, scaled from 595 to 511 nm, gave 0.380, 0.428 and 0.464 um by
   session.
   Method. (a) The flat-stretch estimate of v1 (DESIGN 4.9) with psf.min_stretches N_s. (b) A foci estimate: every
   focus of the run's `spots` table (RAD-51) in a nucleus that phase A processed is fitted, in the foci channel, by
   weighted least squares with a 3D Gaussian of diagonal covariance (axes z, y, x) plus a constant (amplitude A,
   constant B, centre, sigma_z, sigma_y, sigma_x; the centre within 0.4 um axially and 0.3 um laterally of the spot,
   sigma_z in [0.08, 2.0] um, the lateral sigmas in [0.05, 1.0] um), the weights 1 / var(I) of the image's noise
   model at the observed counts (floored at 25 counts^2, as in the diagnostics of entry 9), on the raw counts of the
   window of +/-psf.foci_window_um (1.2 um axial, 0.54 um lateral, rounded to whole voxels) around the spot's voxel.
   A fit is kept when the window lies inside the nucleus crop (crop_pad_um around the label, inside the stack), no
   other spot of the table lies within psf.foci_isolation_um (1.0 um, 3D), both lateral sigmas are at most
   psf.foci_max_sigma_xy_um (0.25 um) and their ratio at most psf.foci_max_ratio (1.5), sigma_z lies more than 1
   percent inside its bounds, the relative residual RMS (the RMS of the unweighted residuals over A) is at most
   psf.foci_max_rel_rms (0.2) and A is at least psf.foci_min_snr (8) times its standard error (the fit covariance
   scaled by the reduced chi-square). With at least psf.foci_min kept fits (N_f) and both emission wavelengths (the
   central_element and foci channels of the source metadata; without both there is no foci estimate) the foci
   estimate is the median of the kept sigma_z and the median of the kept sqrt(sigma_y sigma_x), each times the ratio
   of the central_element to the foci emission wavelength. (c) The variance of each estimate is that of a bootstrap
   over its items (stretches, or kept fits; 200 resamples, numpy seed 0); with both estimates the PSF is, for each
   axis, their inverse-variance weighted mean (source both), with one it is that one (source stretches or foci),
   with none it is the fallback with the image flag psf_fallback as in v1 (source fallback). One PSF per image, no
   depth dependence. With psf.foci "off" the record and every output are those of v1 (sources estimated and
   fallback). The record gives both estimates with their item counts and standard errors, the kept and rejected foci
   by reason, the wavelengths and the source; new image flag psf_from_foci (info) when the foci enter the PSF. New
   keys (fixed group): psf.foci ("on"), psf.foci_min (N_f), psf.foci_window_um [1.2, 0.54], psf.foci_isolation_um
   1.0, psf.foci_max_sigma_xy_um 0.25, psf.foci_max_ratio 1.5, psf.foci_max_rel_rms 0.2, psf.foci_min_snr 8;
   psf.min_stretches becomes N_s. The package VERSION becomes 2.0 with this change (entry 11 (h)).
   Development family psf_check (not registered; its own family code; seeds as entry 11 (a)): one image of six dense
   nuclei of 2.3 in a row (oocyte and spermatocyte alternating, the 0.25 um clearance of entry 2 (c)) and three
   isolated arcs as in the iso family but tilted U(30, 75) degrees, each in its own enclosing nucleus; the SYP
   channel is rendered with the PSF of the scene, sigma_z 0.33, 0.38, 0.43 or 0.48 um (seed modulo 4) and sigma_xy
   U(0.135, 0.165) um, and the nucleoplasm is blurred with it; a foci channel holds RAD-51-like point spots in the
   dense nuclei, a Poisson number of mean 0.15 per um3 of the nucleus, uniform in its sphere, peak U(150, 450)
   counts above the nucleoplasm, rendered with the scene PSF times 595/511 in every axis; the spot table is the
   truth with the localisation noise of 4.4; the wavelengths are 511 nm (central_element) and 595 nm (foci);
   everything else as in 2.3. N_s and N_f are chosen on psf_check at the tuning seeds: the items of each source (the
   flat stretches, the kept foci fits) are pooled over the scenes of one sigma_z level, and for n in {5, 8, 12, 20}
   (stretches) and {15, 25, 40, 60} (foci) 500 random draws of n items without replacement (numpy seed 0) give the
   estimate of that source alone and its relative error of sigma_z against the level; N_s and N_f are the smallest n
   at which at every level the median of the relative error lies within +/-3 percent and the 90th percentile of its
   absolute value is at most 8 percent, else the largest n of the set. Targeted metrics, at the tuning and then once
   at the confirmation seeds: (i) the same error criterion for the combined estimate (N_s stretches and N_f foci
   drawn together, 500 draws, combined as in (c)) at every level; (ii) on the arcs of the psf_check scenes at
   sigma_z 0.43 and 0.48 um the mean length error lies closer to 0 with this change than with the configuration
   before it (psf.foci "off", psf.min_stretches 20). The change is adopted when (i) and (ii) hold at the
   confirmation seeds and the rules of entry 11 (b) hold. Reported, choosing nothing: the estimates, sources and
   standard errors of the twelve development images at Z1, and the difference of the two sources where both exist.
13. 2026-10-24, the v2 development round (per-image PSF, departures from entry 12), written before any tuning or
   confirmation number of entry 12 was computed (no test result had been computed, no test or reserved image was
   run, no test or reserved truth was read and no sealed seed was generated). Motivated by the design seeds of
   psf_check (2000 to 2049, entry 11 (a); C:/Users/ryane/sc_ribbon_work/v2/design/ and scripts/design_psf*.py), with
   the code of entry 12 before its commit:
   (a) With the constant background of entry 12 (b), 55 to 73 percent of the isolated foci whose window lay inside
   the crop failed the lateral-sigma rule, and the pooled median sigma_z of the kept fits was 2.6 to 5.9 percent
   above the scene PSF at the four levels. Most foci lie within 1 um of the blurred edge of the nucleoplasm, which a
   constant cannot absorb: fits of isolated foci on a flat background (60 per case, two PSFs, peaks 150 to 450) were
   within 1.5 percent and all compact, and joint fits of the neighbouring spots inside the window (sharing the
   sigmas) changed nothing. The fit of entry 12 (b) therefore has a planar background, B plus a gradient in z, y and
   x over the window (three more parameters), and is otherwise unchanged; on the design seeds 87 to 96 percent of
   those foci then pass and the pooled median sigma_z of the kept fits lies within -1.8 to +3.0 percent of the scene
   PSF.
   (b) Both sources over-estimate the lateral sigma on the design seeds (the stretches by 7 to 10 percent, the foci
   with the planar background by 6 to 11 percent), while the fallback lateral sigma (0.15 um) is the mean of the
   psf_check range and the value rendered in 2.3. The image PSF therefore takes only sigma_z from the estimates of
   entry 12; its sigma_xy stays psf.fallback_sigma_um[1], the combination of entry 12 (c) is that of sigma_z alone,
   and the lateral estimates of both sources are recorded and choose nothing. The flag psf_from_foci, the sources,
   N_s and N_f, the selection on the tuning seeds and the targeted metrics of entry 12 are unchanged.
14. 2026-10-24, the v2 development round (per-image PSF, revision 1 of entries 12 and 13), written before any number
   of this revision was computed (no test result had been computed, no test or reserved image was run, no test or
   reserved truth was read and no sealed seed was generated). N_s and N_f were chosen by the rule of entry 12 on the
   psf_check tuning seeds (10000 to 10099, commit b0cf9ac;
   C:/Users/ryane/sc_ribbon_work/v2/logs/step1_psf_select_b0cf9ac_tuning.json): the stretch estimate alone had a
   median error of sigma_z of +2.1 to +5.3 percent at every n (+4.5 to +5.3 percent at sigma_z 0.33 um, +3.3 percent
   at 0.38 um), so the rule was met at no n and N_s is 20 (the largest of the set); the foci estimate met the rule
   at every level from n = 15 (medians -1.8 to +2.8 percent, 90th percentiles of the absolute error at most 5.6
   percent), so N_f is 15; the combined estimate of entry 12 (c) with N_s 20 and N_f 15 failed targeted metric (i)
   at sigma_z 0.33 um (median +4.2 percent) and met it at the other levels. Targeted metric (ii) was not computed.
   Revision 1, motivated by these tuning-seed results: the image PSF takes sigma_z from the foci estimate when it
   exists (source foci), else from the stretch estimate (source stretches), else the fallback (source fallback); the
   inverse-variance combination of entry 12 (c) is dropped, and both estimates are still recorded with their
   difference where both exist. psf.min_stretches is 20 and psf.foci_min 15; everything else of entries 12 and 13 is
   unchanged. Targeted metric (i) is now the error criterion of entry 12 for the image PSF this rule makes from N_s
   stretches and N_f foci drawn together (500 draws) at every level, and (ii) is unchanged; the adoption rule of
   entry 12 holds for this revision, on the tuning and then once on the confirmation seeds.
15. 2026-10-24, the v2 development round (per-image PSF, revision 2 of entries 12 to 14), written before any
   confirmation number of the change was computed (no test result had been computed, no test or reserved image was
   run, no test or reserved truth was read and no sealed seed was generated). Revision 1 (commit 8b744af) on the
   tuning seeds (C:/Users/ryane/sc_ribbon_work/v2/logs/step1_8b744af_tuning.json and
   step1_psf_select_8b744af_tuning.json): every psf_check scene took the foci estimate (25 to 29 kept fits per
   scene), with per-image relative errors of sigma_z of median +2.9, +1.1, -0.5 and -1.5 percent at sigma_z 0.33,
   0.38, 0.43 and 0.48 um (90th percentiles of the absolute error 4.9, 2.9, 1.4 and 3.1 percent), against +10.6,
   -4.0, -15.1 and -24.0 percent with the fallback; targeted metric (i) as entry 14 defined it, on draws of N_s 20
   stretches and N_f 15 foci, failed at 0.33 um only (median +3.1 percent); targeted metric (ii) cannot measure the
   change: the mean length errors of the arcs were 0.1231 um with and without it at 0.43 um and 0.1723 and 0.1721 um
   at 0.48 um, because the lengths of the end rule do not depend on the image PSF in these scenes (the fitted edge
   sigma stays inside its PSF-scaled bounds). The criteria of 8.1 on the registered families were identical (their
   scenes have no wavelengths and fewer than 20 stretches, so they keep the fallback). Revision 2 changes nothing in
   the method of revision 1. Its targeted metrics, stated now: (i) the image PSF that the pipeline makes for each
   psf_check scene (sigma_z of its record) has, at every level, a median relative error within +/-3 percent and a
   90th percentile of its absolute error of at most 8 percent; (ii) is dropped; the tolerances and the cost rule of
   entry 11 (b) stay. The tuning values of this (i) are those above, seen before this entry, so the confirmation
   seeds decide alone: the change is adopted when (i) and the tolerances hold there (once, seeds 11000 to 11099) and
   the cost rule holds on n2rad51_noHS_herm_000; otherwise it is withdrawn (psf.foci off) and the round continues
   with the PSF of v1. No further revision of the per-image PSF is made in this round.
16. 2026-10-24, the v2 development round (puncta), written before any number of this change was computed (no test
   result had been computed, no test or reserved image was run, no test or reserved truth was read and no sealed
   seed was generated). The per-image PSF of entries 12 to 15 was withdrawn by the rule of entry 15: on the
   confirmation seeds (11000 to 11099, commit b58f6d3;
   C:/Users/ryane/sc_ribbon_work/v2/logs/step1_b58f6d3_confirm.json and step1_b58f6d3_confirm_psf_errors.json) the
   image PSF of the psf_check scenes had median relative errors of sigma_z of +3.4, +0.9, -0.4 and -1.8 percent at
   sigma_z 0.33, 0.38, 0.43 and 0.48 um (90th percentiles of the absolute error 5.8, 2.9, 1.7 and 3.3 percent; the
   fallback of v1 +10.6, -4.0, -15.1 and -24.0 percent), so the restated targeted metric (i) failed at 0.33 um; the
   criteria of 8.1 on the registered families were identical, and the cost rule held on n2rad51_noHS_herm_000
   (median time per nucleus 0.172 s against 0.173 s, 3.0 GB; there the foci gave sigma_z 0.363 um and the flat
   stretches of v1 0.3675 um). psf.foci is off from commit 4cd1452, which keeps the code behind the switch: the
   image PSF is that of v1. The hold-out of entry 11 (e): the germquant run of 20251021_n2_nohs_HERM_2 (the session
   channel map of v2/channel_maps, SC on 545) finished and passed the completion and channel checks (central_element
   is the nd2 channel 545; Cellpose segmentation, 935 nuclei), but its process peaked at 42.6 GB against the 8 GB
   per-process budget of the shared workstation, so the germquant runs of the other five hold-out gonads were not
   made; H1 will run on the hold-out gonads that have a germquant run and, with fewer than three evaluable, is
   reported and decides nothing.
   Motivated by development results (entries 8 and 9): on development phantoms v1 makes a strand of 85 percent of
   isolated PSF-blurred point puncta (99 percent of them steep) and of 81 percent of aggregates of 0.25 to 0.8 um,
   and dense nuclei with 6 puncta gain 14 percent of total length; rule B1 at 1.3 times sigma_z with length_xy_um <=
   0.5 um removed 83 to 86 percent of the point puncta, 72 to 73 percent of the aggregates and 69 percent of the
   dense-nucleus puncta, none of 60 tilted arcs of 1.5 to 6 um, 1 of 45 steep arcs of 0.9 to 1.2 um and none of 33
   of 1.2 to 1.5 um (29 of 42 of 0.6 to 0.9 um); on the twelve development images B1s removed 6.6 to 29.7 percent of
   the strands of a gonad. No S1 gate and no S2 metric has been computed with these rules.
   Method (B1s as registered in entries 8 and 9, with the image's sigma_z of the PSF of v1 (the per-image PSF being
   withdrawn)): after measurement, a traced strand with frac_steep60 > 0.5 and length_xy_um <= punctum.max_xy_um
   whose core profile over s from -0.5 um to length_um + 0.5 um (every 0.05 um), fitted by least squares with a
   Gaussian plus a constant, has sigma times tz_abs_mean at most punctum.sigma_ratio times the image's sigma_z gets
   strand_class punctum and the strand flag punctum (info); the class order is verified, punctum, fragment,
   merged_suspect, complete, partial. Punctum strands stay in the strand, sample and strip outputs and are left out
   of total_length_um, total_length_xy_um, n_strands, n_main, n_fragments, n_complete, n_merged_suspect,
   count_status, complete_nucleus, the pieces and their columns (n_pieces, mean_piece_um, longest_piece_um,
   piece_lw_mean_um), n_gaps, gaps_per_10um, gap_fraction, multiplicity_frac, the multiplicity reference of the
   nucleus and the image-level references of the multiplicity and of double_width, the strand tubes of
   on_strand_frac and coverage_proxy (their voxels count as territory) and the foci assignment; the PSF stretches
   never hold one (a stretch is flat over 1 um). They still count as neighbours in the contamination of other
   strands and in the nucleus's snr_core_median, as in v1. New nucleus columns: n_puncta; puncta_signal_frac (the
   SYP above background inside the punctum tubes and outside the strand tubes over that in the territory, so that
   on_strand_frac + puncta_signal_frac <= 1); steep_length_frac (the share of total_length_um steeper than 60
   degrees to the image plane). New strand column punctum_sigma_z_um (the fitted sigma times tz_abs_mean; NaN where
   the rule does not apply). New keys (fixed group): punctum.enabled (true; false gives v1), punctum.exclude (true;
   false keeps the class and the flag and returns punctum strands to every column above as in v1: the fallback of
   entry 11 (c) at Z1), punctum.max_xy_um 0.5, punctum.sigma_ratio 1.3. The S1 scoring
   (germquant.benchmark.sc_ribbons.score_scene) scores the strands the nucleus tier counts: under punctum.exclude
   the punctum strands are left out of the totals, the coverage, the off-truth length and the matching (v1 has no
   such class, so no v1 score changes).
   New development families, not registered, with the family codes of the entry-8 generator
   (v2/scripts/steep_phantoms.py, so that a seed gives the same scene there): puncta_pt and puncta_agg (a nucleus of
   2.3 holding 8 PSF-blurred points, or spheres of diameter U(0.25, 0.8) um, at least 1.2 um apart, and no ribbon)
   and dense_puncta (the dense nucleus of 2.3 plus 6 points at least 0.6 um from every ribbon), puncta peaks
   log-uniform in [100, 1500] counts above the nucleoplasm; dense_puncta_clean is each dense_puncta scene with the
   same draws and the same noise but without the puncta. A strand belongs to a punctum when half or more of its
   in-strand samples lie within 1.5 um axially and 0.4 um laterally of the punctum's centre and under half within
   the protocol tolerance of a ribbon (the attribution of entry 8). Targeted metrics, on the tuning and then once on
   the confirmation seeds, against the configuration before this change (punctum.enabled false) on the same scenes:
   (a) the excess of the dense_puncta nucleus total (total_length_um) over that of dense_puncta_clean, the mean over
   the scene pairs, falls by at least 60 percent (confirmation: by at least 30 percent); (b) at least 80 percent of
   the point puncta and 65 percent of the aggregates leave no strand outside the punctum class (confirmation: the
   same levels); (c) on the registered families no strand lying half or more within the tolerance of an iso, gaps or
   foci ribbon becomes punctum, and punctum strands lying on dense ribbons hold at most 0.5 percent of the dense
   truth length (confirmation: the same). Adopted when (a) to (c) and the rules of entry 11 (b) hold. Reported,
   choosing nothing: the share of puncta inside a strand of two or more members (0.3 um of its samples within the
   punctum's ellipsoid). Z1 (entry 11 (c)) follows this entry.
17. 2026-10-24, the v2 development round (puncta, revision 1 of entry 16), written before any number of this
   revision was computed (no test result had been computed, no test or reserved image was run, no test or reserved
   truth was read and no sealed seed was generated). Entry 16 (commit 7c3a25e) on the tuning seeds against
   punctum.enabled false (C:/Users/ryane/sc_ribbon_work/v2/logs/step2_7c3a25e_tuning.json): (a) failed, the
   dense_puncta excess fell by 36 percent (6.47 to 4.13 um per nucleus, of which 2.16 um on strands attributed to
   puncta but not of class punctum and 1.97 um on ribbon strands that differ from the clean scenes); (b) held, 83
   percent of the point puncta and 78 percent of the aggregates left no strand outside the class; (c) held, no
   strand on an iso, gaps or foci ribbon and none on a dense ribbon became punctum; the criteria of 8.1 on the
   registered families were identical; 38 percent of the dense-nucleus puncta and 22 to 24 percent of the isolated
   ones lay inside strands of two or more members. That share is above the trigger of PLAN_V2.md section 5 for the
   deferred item "B1 before end matching" (10 percent), and revision 1 is that item: in phase A, after the tracer,
   the axial shadows and the T contacts and before the end matching, a traced strand whose polyline has frac_steep60
   > 0.5 (tangents over +/-0.1 um) and length_xy_um <= punctum.max_xy_um, and whose core profile (the core ellipse
   of the measurement, background-subtracted raw counts) along the polyline extended straight by 0.5 um at both
   ends, fitted by the Gaussian of entry 16, has sigma times its mean |t_z| at most punctum.sigma_ratio times the
   fallback sigma_z (the image PSF is estimated after phase A), takes no part in the end matching and rejoins the
   strands after it; the punctum class of entry 16 is still decided on the measured strand. New key (fixed group):
   punctum.before_linking (true; false gives entry 16). Everything else of entry 16, its targeted metrics (a) to (c)
   and their confirmation levels, is unchanged; they are computed again on the tuning seeds against the
   configuration before entry 16, then once on the confirmation seeds.
18. 2026-10-24, the v2 development round (puncta, revision 2 of entries 16 and 17), written before any confirmation
   number of the puncta change was computed (no test result had been computed, no test or reserved image was run, no
   test or reserved truth was read and no sealed seed was generated). Revision 1 (commit 55e7d05) on the tuning
   seeds against the configuration before entry 16
   (C:/Users/ryane/sc_ribbon_work/v2/logs/step2r1_55e7d05_tuning.json): (a) failed again, the dense_puncta excess
   fell by 54 percent (6.47 to 2.96 um per nucleus: 1.69 um on strands attributed to puncta but not of class
   punctum, mostly steep strands longer than punctum.max_xy_um in the image plane or whose profile fit ran to its
   bound, and 1.27 um on ribbon strands that differ from the clean scenes); (b) held (90 percent of the point puncta
   and 80 percent of the aggregates left no strand outside the class); (c) held; the criteria of 8.1 on the
   registered families moved by at most 0.06 percentage points (dense mean total error -6.92 to -6.98 percent, SD
   8.59 to 8.56 percent, CCC 0.535 to 0.536). The removal of puncta from the nucleus tier is therefore withdrawn:
   the remaining excess sits on strands that the rule cannot tell from steep ribbons with the constants that the
   plan holds fixed (PLAN_V2.md 2.5) and on ribbon strands, which no punctum rule removes. Revision 2 keeps the
   class and the flag as a label only: punctum.exclude false (punctum strands count in every column as in v1, the
   fallback of entry 11 (c)) and punctum.before_linking false (the end matching as in v1), with the rule of entry 16
   on the measured strand, the columns n_puncta, puncta_signal_frac (here the share of the territory's SYP inside
   the punctum tubes, which overlap the strand tubes) and steep_length_frac, and the S1 scoring of all strands. No
   length, continuity or intensity column of v1 changes; the strand class changes for punctum strands only (they
   leave partial, fragment or merged_suspect, so n_merged_suspect can fall; a steep strand is never complete in v1).
   Targeted metrics of revision 2, once on the confirmation seeds against the configuration before entry 16: (b) at
   least 80 percent of the point puncta and 65 percent of the aggregates are labelled with no strand outside the
   class (their tuning values for this configuration are those of entry 16 on the tuning seeds, 83 and 78 percent,
   seen before this entry); (c) no strand on an iso, gaps or foci ribbon and at most 0.5 percent of the dense ribbon
   length is labelled punctum; and the tolerances of entry 11 (b). The change is adopted when they hold, else
   withdrawn (punctum.enabled false). Z1 is not run: with the per-image PSF withdrawn and the puncta a label only,
   no input of S2 differs from trial 4, so zone truth is not used here; the realism probe measures the trial-4 runs
   (runs/6b7be22) instead of Z1 runs.
19. 2026-10-24, the v2 development round (realism probe; diagnostics only, no method change), written before any
   number below was computed (no test result had been computed, no test or reserved image was run, no test or
   reserved truth was read and no sealed seed was generated). The puncta change of entries 16 to 18 was withdrawn by
   the rule of entry 18: revision 2 (commit 215e24d) on the confirmation seeds (11000 to 11099;
   C:/Users/ryane/sc_ribbon_work/v2/logs/step2r2_215e24d_confirm.json) labelled 84 percent of the point puncta and
   79 percent of the aggregates with no strand outside the class and 0.05 percent of the dense ribbon length, but
   one strand lying on a gaps ribbon (seed 11049, a 1.26 um vertical strand on a ribbon tilted 28 degrees,
   length_xy_um 0.48) was labelled punctum, so targeted metric (c) failed; the criteria of 8.1 were identical;
   punctum.enabled is false from commit bb07ede, which keeps the code behind the switch. The configuration is
   therefore that of v1 (trial 4) for every output, and Z1 was not run (entry 18).
   Motivated by development results (entries 7 to 10): in the traced development nuclei with 30 um of SC or more the
   median piece_lw_mean_um is 2.6 um (entry 7) and the in-plane strands of SCfull nuclei have a median length of 2.4
   to 2.5 um with 77 to 88 percent of their ends free (entry 10), while continuous dense phantoms give a
   piece_lw_mean_um of 5.2 um for 5.8 um of truth; 39 to 46 percent of the real traced length runs steeper than 60
   degrees, 91 percent in the dimmest octile of nuclei, against 16 to 23 percent on dense phantoms with white noise.
   The phantoms of 2.3 have one PSF, a uniform nucleoplasm, a constant brightness along each ribbon and no puncta.
   New development family real_dense, not registered (its own family code; seeds as entry 11 (a)): the dense nucleus
   of 2.3 (the 0.25 um clearance of entry 2 (c)) at a scene peak of 200, 300, 430 or 700 counts (seed modulo 4, as
   the snr family), with four features drawn from separate random streams so that a feature switched off leaves the
   others' draws unchanged, each also switched on alone (variants base, psf, brightness, texture, puncta and all):
   (a) the axial PSF sigma drawn uniformly from the twelve per-image values of the compact RAD-51 foci of entry 9
   (v2/steep/psf, the median sigma_z of each image times 511/595; the lateral sigma 0.149 um of 2.3); (b) each
   ribbon's brightness times exp(G(s) - var/2), G a stationary Gaussian process along s with an exponential
   covariance whose variance and correlation length are the medians over the twelve trial-4 runs of the estimates on
   in-plane (frac_steep60 <= 0.1), isolated (frac_contaminated <= 0.1) strands of 3 um or more without gaps in
   traced nuclei with total_length_um >= 30 um (the log of the core profile over its median, every 0.1 um, the mean
   noise variance (core_se / core)^2 removed; the correlation length from the slope of the log autocovariance at 0.2
   to 1.0 um); (c) a nucleoplasmic texture: unit white noise blurred with the scene PSF inside the nucleus, scaled
   so that after a second Gaussian filter of the PSF size its standard deviation over the nucleus is sqrt(r) times
   the nucleus's nucleoplasm level, r the median over the twelve trial-4 images of (var_s - var_n) / m^2 in at most
   120 traced no_sc nuclei per image (var_s the variance of the PSF-smoothed counts over the territory eroded by two
   voxels, var_n the noise model's variance after the same filter, m the mean count above the camera offset of 100);
   (d) puncta, half points and half spheres of diameter U(0.25, 0.8) um, peaks log-uniform in [100, 1500] counts, at
   least 0.6 um from every ribbon and 1.2 um apart, their number Poisson with the mean count per traced nucleus of
   rule B1s (entry 9, with the sigma_z of each image's trial-4 PSF record and the core-profile fits of v2/steep/psf)
   in the tercile of the nuclei's SYP brightness (central_element_mean_intensity) that the scene peak maps to (200:
   lowest, 300 and 430: middle, 700: highest). None of these measurements reads a zone label. Reported, choosing
   nothing, for the configuration of v1 on the tuning seeds (10000 to 10099) per variant: S1b (mean, SD and CCC of
   the total error), strands and main strands per nucleus against the expected count, the end classes, the steep
   share of the traced length (pooled and by scene peak), the median length of in-plane strands (frac_steep60 <=
   0.5) and the median piece_lw_mean_um; beside the same statistics of the trial-4 runs in traced nuclei with
   total_length_um >= 30 um (the steep share also by octile of SYP brightness). A feature reproduces a real symptom
   when switching it on alone moves the phantom statistic at least half of the way from its base value to the real
   one. From here on the synthetic panels of this round also report real_dense with all four features (descriptive).
20. 2026-10-24, the v2 development round (contact stop), written before any number of this change was computed (no
   test result had been computed, no test or reserved image was run, no test or reserved truth was read and no
   sealed seed was generated). The realism probe of entry 19 (commit ee4d056, tuning seeds, the configuration of v1;
   C:/Users/ryane/sc_ribbon_work/v2/logs/realism_probe_ee4d056_tuning.json) found that no feature switched on alone
   reproduces a real symptom: the texture moved the phantom statistics 17 to 29 percent of the way to the real ones
   (steep share of the traced length 0.17 to 0.21 against 0.40 real; median in-plane strand 5.36 to 4.52 um against
   2.49; median piece_lw_mean_um 5.03 to 4.33 um against 2.63; main strands over the expected count 0.93 to 1.23
   against 2.19), the brightness modulation 30 percent for piece_lw_mean_um and at most 10 percent otherwise, the
   PSF and the puncta at most 12 percent; all four together moved piece_lw_mean_um 61 percent of the way (3.58 um),
   the steep share 43 percent (0.27) and the in-plane strand length 37 percent (4.29 um), and the dense S1b of that
   variant was much worse (SD 18.3 percent, CCC 0.13, against 10.4 percent and 0.32 without the features). So no
   realism-driven change is made in this round (step 9 of the plan is skipped) and Z3 is not run; real_dense with
   all four features is reported beside the panels. Motivated by development results (the taxonomy of entry 11): 91
   percent of the 662 contact stops of the tracer on the dense phantoms hit a strand lying at more than 32 degrees
   to the trace (median 64), because the stop compares the trace with the local Hessian direction, which inside an
   unresolved junction is the trace's own; such stops are involved in 378 of the 470 tracer breaks and in 132 of the
   193 final breaks. Method: every claimed voxel keeps the unit tangent of the strand that claimed it, taken at the
   claiming point (central differences of the traced polyline); a trace that enters another strand's claim tube
   stops as contact only when |t . t_claim| > follow.contact_cos (0.85, unchanged), whatever the local lineness, and
   otherwise continues by pure prediction, as inside any claim tube in v1. New key (fixed group):
   follow.contact_test ("claim_tangent"; "local" gives v1). Targeted metrics (dense phantoms with the taxonomy of
   v2/scripts/dense_taxonomy.py, and the pairs family), against the configuration before this change on the same
   scenes: breaks per nucleus after the tracer and in the final strands fall by at least 30 percent (confirmation:
   15 percent); switches (a strand leaving one ribbon for another, the final strands) per nucleus rise by at most 10
   percent (confirmation: the same bound); the dense S1b SD falls (confirmation: does not rise); the share of X
   crossings of the pairs family with both ribbons matched does not fall by more than 0.02. The tuning seeds run the
   dense and pairs families; the confirmation seeds run every registered family for the tolerances of entry 11 (b),
   once; the cost rule is checked on n2rad51_noHS_herm_000.
21. 2026-10-24, the v2 development round (end matching), written before any number of this change was computed (no
   test result had been computed, no test or reserved image was run, no test or reserved truth was read and no
   sealed seed was generated). The contact stop of entry 20 (commit 40fe390) failed two of its targeted metrics on
   the tuning seeds against the configuration before it
   (C:/Users/ryane/sc_ribbon_work/v2/logs/step4_40fe390_tuning.json): breaks after the tracer fell by 18 percent
   only (4.88 to 4.01 per nucleus; at least 30 required) and switches in the final strands rose by 13.5 percent
   (3.34 to 3.79 per nucleus; at most 10 allowed), while the final breaks fell by 39 percent (1.99 to 1.21), the
   dense S1b SD fell (8.59 to 8.47 percent) and the X crossings of the pairs family kept both ribbons as often
   (0.96); the dense S1b mean and CCC, which no metric of entry 20 targets, moved from -6.92 to -1.37 percent and
   from 0.535 to 0.699. Tracer events on the same scenes (an instrumented copy of the tracer,
   v2/scripts/trace_events.py; logs/trace_events_s4pre_tuning.json and trace_events_s4_tuning.json) show why:
   contact stops fell from 675 to 163, but of the 414 passes through the claim tube of a strand on another ribbon by
   a trace whose own ribbon goes on for more than 1 um beyond it, 44 went on along their own ribbon, 246 stopped
   within 0.9 um after the tube and 104 inside it; at the exits of the 246 the trace direction lay a median 39
   degrees off the ribbon (20 degrees for the passes that went on), 0.21 um laterally and 0.36 um axially from it,
   so pure prediction does not regain the ribbon behind an unresolved junction; the added switches are spurs at
   strand ends (tracer 1.51 to 2.31 per nucleus, final 1.22 to 1.39; on a single T contact the stem now runs through
   the bar and stops beyond it, three seeds) and switches through junctions after the end matching (2.12 to 2.40).
   One revision was explored on the design seeds (2000 to 2049 and 3000 to 3049) with code that was never committed
   (v2/logs/dev/): no weak run inside another strand's tube and, for a trace that stops before regaining its own
   ridge, a cut back to the tube entry as a contact end; against v1, tracer breaks fell by 26 and 23 percent, final
   breaks by 8 and 18 percent (entry 20 there: 32 and 57 percent), switches by 17 and 22 percent, and the dense mean
   total error returned to that of v1 (-5.3 percent on both blocks against -4.9 and -4.8). Neither variant meets the
   tracer target, and regaining a ribbon 39 degrees off the predicted direction needs a search inside junctions,
   which is the question of the junction re-pairing step; entry 20 is therefore withdrawn without a revision:
   follow.contact_test is "local" by default again (the switch and its tests stay), and the decomposition continues
   from that of v1. Motivated by development results (the taxonomy of entry 11, made with the decomposition of v1):
   of the 411 tracer breaks with facing ends, 37 percent were linked correctly, 30 percent were rejected by the 50
   degree rule although the two anchors lay nearly in line (median offsets 0.11 um axial and 0.07 um lateral over
   0.8 um: inside the junction the end direction hooks toward the neighbour), 22 percent had turned onto another
   ribbon there and 8 percent lost the right partner to a wrong one; strands that join two ribbons end to end have a
   median length of 10.3 um and strands that turn at a T contact 8.3 um, while no correctly decomposed strand
   exceeds 7.5 um (truth maximum 7.0 um). Method: (a) when points of an end's last 1.0 um lie within one PSF FWHM
   (axial and lateral) of another strand's points, the end's anchor and outward direction (link.tip_um and
   link.dir_window_um, as in v1) are measured on the strand from the first point beyond the last of them, unless
   less than 0.3 um of the strand lies beyond it (then as in v1); the gap is still counted from the tip (the
   anchor's arc offset from the tip includes the shift); the FWHM is that of psf.fallback_sigma_um (0.86 um axial,
   0.35 um lateral), because the image PSF is estimated after phase A; (b) after the matching, a chain longer than
   link.max_chain_um is cut at its most expensive link, again until no chain is longer or no link is left; cut ends
   are free and get the end rule. link.max_chain_um is 9.0 um, set without phantoms: the longest pachytene SC of
   wild-type C. elegans is about 7 um (the field review of the build, section 0) and the registered iso, gaps and
   foci ribbons reach 8 um, so no registered ribbon can reach the bound; it is held fixed in this round although
   link.* is tunable. New keys: link.junction_dirs (true; false gives v1; fixed group) and link.max_chain_um (9.0;
   null gives v1). Targeted metrics (the dense family, the taxonomy script, which takes each end's anchor and
   direction as the end matching computed them and counts only the links that the chain bound keeps), against the
   configuration before this change on the same scenes: angle rejections at breaks with facing ends fall by at least
   a third (confirmation: a sixth); wrong-ribbon links fall by at least 30 percent (confirmation: 15 percent);
   correct links do not fall (confirmation: the same); strands longer than 9 um that hold more than one ribbon fall
   by at least half (confirmation: a quarter). The tuning seeds run the dense, gaps and iso families; the
   confirmation seeds run every registered family for the tolerances of entry 11 (b), once; the cost rule is checked
   on n2rad51_noHS_herm_000.
22. 2026-10-24, the v2 development round (end matching, revision 1 of entry 21), written before any tuning or
   confirmation number of this revision was computed (no test result had been computed, no test or reserved image
   was run, no test or reserved truth was read and no sealed seed was generated). Entry 21 (commit f39cb85) on the
   tuning seeds against the configuration before it
   (C:/Users/ryane/sc_ribbon_work/v2/logs/step5_f39cb85_tuning.json) failed three of its four targeted metrics and
   every S1b tolerance: angle rejections fell by 18 percent (127 to 104; a third required), wrong-ribbon links rose
   from 107 to 181 and correct links fell from 179 to 100, while strands longer than 9 um holding more than one
   ribbon fell from 0.64 to 0.07 per nucleus; the dense S1b mean went from -6.92 to -17.19 percent, the SD from 8.59
   to 12.29 percent, the CCC from 0.535 to 0.210 and the coverage from 0.958 to 0.864. Two causes: with an anchor
   measured beyond a junction, the anchor estimate of the gap (the distance between the anchors less their offsets
   from the tips) is 0 for tips that lie microns apart, and such pairs fall under the overlap rule that judges only
   the collinearity of the two directions (on the first 20 tuning scenes the wrong links had a median distance of
   1.4 um between their tips, 90th percentile 3.0 um, with estimated gaps of 0); and the chain bound and the cycle
   cuts then removed 56 of the links that the matching had made correctly at tracer breaks (1 in v1), because
   chimeric chains grew past 9 um. Explored on the design seeds (2000 to 2049 and 3000 to 3049) with code that was
   never committed (v2/logs/dev/), against the configuration before entry 21, with the taxonomy script applying the
   rule of the end matching to the pairs it diagnoses: with the gap of such pairs taken between the tips and their
   angles to the line between the anchors, wrong-ribbon links fell from 45 to 18 and from 48 to 23 with the chain
   bound and to 25 and 35 without it; correct links went from 91 and 101 to 86 and 82 with the bound (85 and 86 when
   it cuts at the largest turn instead) and to 93 and 101 without it; angle rejections went from 63 and 68 to 65 and
   72 either way, so measuring the directions beyond the junction does not rescue the pairs that the 50 degree rule
   rejects; without the bound the dense S1b mean went from -4.9 to -3.9 and from -4.8 to -3.2 percent, the SD from
   8.1 to 7.2 and from 9.0 to 8.3 percent and the CCC from 0.579 to 0.669 and from 0.594 to 0.643. Revision 1,
   motivated by these results: (a) of entry 21 stays, and for a pair in which at least one end was measured beyond a
   junction the gap is the distance between the two tips and the angles are those of the two outward directions to
   the line between the two anchors (the overlap rule of v1 is not applied to such pairs; anchors closer than 0.05
   um are judged by the collinearity of the directions as in v1); every other pair is judged as in v1; (b) of entry
   21 is withdrawn: link.max_chain_um is null by default (the switch and its tests stay). New key (fixed group):
   link.junction_gap ("tips"; "anchors" gives entry 21). The angle-rejection metric of entry 21 is dropped, because
   the design seeds show that (a) does not act on those rejections, and the long-chimera metric goes with (b); both
   are reported. Targeted metrics of revision 1 (the dense family, the taxonomy script as above), against the
   configuration before entry 21 on the same scenes: wrong-ribbon links fall by at least 30 percent (confirmation:
   15 percent); correct links do not fall (confirmation: the same); the dense S1b SD falls (confirmation: does not
   rise). The tuning seeds run the dense, gaps and iso families; the confirmation seeds run every registered family
   for the tolerances of entry 11 (b), once; the cost rule is checked on n2rad51_noHS_herm_000.
23. 2026-10-24, the v2 development round (complete class), written before any tuning or confirmation number of this
   change was computed (no test result had been computed, no test or reserved image was run, no test or reserved
   truth was read and no sealed seed was generated). Revision 1 of the end matching (entry 22, commit ecf52a3) on
   the tuning seeds against the configuration before entry 21
   (C:/Users/ryane/sc_ribbon_work/v2/logs/step5r1_ecf52a3_tuning.json): wrong-ribbon links fell by 31 percent (107
   to 74; met), correct links fell from 179 to 172 (failed) and the dense S1b SD rose from 8.59 to 8.81 percent
   (failed); the mean went from -6.92 to -6.15 percent and the CCC from 0.535 to 0.564, angle rejections from 127 to
   130 and strands longer than 9 um holding more than one ribbon from 0.64 to 0.40 per nucleus. The end matching
   change is withdrawn without a second revision (link.junction_dirs is false by default; the switches and their
   tests stay): directions measured beyond a junction rescue none of the pairs that the 50 degree rule rejects, and
   what remains, fewer wrong links at the price of correct ones, gives no rule that the design seeds pointed to. The
   junction end fit of PLAN_V2.md step 6 was implemented as the plan states it and run on the design seeds (2000 to
   2049 and 3000 to 3049) with code that was never committed
   (C:/Users/ryane/sc_ribbon_work/v2/logs/dev/step6_junction_endfit_explored.diff): against v1, the mean shortfall
   at T and end-to-end ends went from 0.39 and 0.38 um to 0.42 and 0.44 um, the dense S1b mean from -4.9 and -4.8 to
   -7.4 and -7.4 percent, the CCC from 0.579 and 0.594 to 0.527 and 0.533, and the pooled S1e precision from 0.910
   and 0.888 to 0.862 and 0.857; with the bound on e at the path's closest approach to the neighbour, at that point
   plus half the lateral FWHM, or with no bound, the shortfall stayed between 0.38 and 0.42 um and the mean between
   -5.6 and -7.2 percent: the neighbour term takes up the ribbon's own end signal. The junction end fit is therefore
   not proposed, and no entry of the round carries it. With the decomposition of v1 the dense S1b on the tuning
   seeds has SD 8.59 percent and CCC 0.535, so the junction re-pairing of PLAN_V2.md step 8 follows this entry.
   Motivated by development results (the taxonomy of entry 11): 27 of the 77 complete strands of the dense phantoms
   were wrong (10 spurs at T or end-to-end contacts, 8 chimeras, 5 pieces of broken ribbons with both ends fitted, 4
   end errors); another strand's end lay within 2 PSF FWHM of the switch point for 92 percent of the crossing turns,
   55 percent of the T turns and 39 percent of the end-to-end joins. Method: besides the rule of v1, a complete
   strand needs (i) length_um <= classes.max_complete_um (9.0 um: the longest pachytene SC of wild-type C. elegans
   is about 7 um, the field review of the build, section 0, and the registered iso, gaps and foci ribbons reach 8
   um; entry 21 set the same bound on chains, which entry 22 withdrew; here it bounds the class only), (ii) no end
   within contamination_um of another strand, and (iii) no end of another strand whose outward direction lies within
   45 degrees of the strand's local tangent (either sense) within classes.complete_guard_fwhm PSF FWHM (axial and
   lateral, the image PSF) of a point of the strand more than 0.5 um from its own ends; the other strands are those
   of the whole image. classes.complete_guard_fwhm is chosen on the tuning seeds from {1.0, 1.5, 2.0}: the smallest
   value with pooled S1e precision >= 0.92 and dense S1e precision >= 0.85 there, else 2.0; the choice is repeated
   by the same rule after any later change of this round that alters the strands, and confirmed once on the
   confirmation seeds with the final configuration. A check of the code on the design seeds 2000 to 2049 (every
   registered family, code not committed) gave pooled S1e precision 0.896 in v1 and 0.933, 0.940 and 0.947 at 1.0,
   1.5 and 2.0 FWHM, dense precision 0.707 and 0.792, 0.826 and 0.882, and an unchanged iso recall (0.92). New keys:
   classes.complete_guard (true; false gives v1), classes.complete_guard_fwhm, classes.max_complete_um (9.0; null
   gives v1). Targeted metrics on the confirmation seeds: pooled S1e precision >= 0.90 and dense S1e precision >=
   0.80; complete recall is reported by family and stays on iso at least 0.9 times that of v1; and the tolerances of
   entry 11 (b). No length, continuity or intensity column changes: the class of a strand that fails the guard
   becomes partial instead of complete, and with it n_complete and n_partial, complete_nucleus and
   total_length_complete_um; no input of S2 depends on the class.
24. 2026-10-24, the v2 development round (complete guard choice; junction re-pairing), written before any
   confirmation number of the complete guard was computed (no test result had been computed, no test or reserved
   image was run, no test or reserved truth was read and no sealed seed was generated). The rule of entry 23 on the
   tuning seeds (commit 5a3481b, every registered family;
   C:/Users/ryane/sc_ribbon_work/v2/logs/step7_5a3481b_tuning_s7_f10.json, _f15.json and _f20.json): pooled S1e
   precision 0.893 in v1 and 0.925, 0.939 and 0.943 with classes.complete_guard_fwhm 1.0, 1.5 and 2.0; dense
   precision 0.649 and 0.696, 0.674 and 0.643; no value reaches 0.85 on dense, so the rule sets 2.0. Complete recall
   at 2.0 against v1: iso 0.95 and 0.95, gaps 0.28 and 0.28, foci 1.00 and 1.00, pairs 0.585 and 0.33, dense 0.091
   and 0.033, snr 0.082 and 0.040; the tolerances of entry 11 (b) held at every value. The junction re-pairing of
   PLAN_V2.md step 8 (the trigger holds: dense S1b SD 8.59 percent and CCC 0.535 on the tuning seeds) was
   implemented as the plan states it, with three additions that keep the traced strands where the pairing agrees
   with them (a pair that restores the traced continuation keeps the traced points of the zone; an unpaired branch
   where the strand itself ended in the zone keeps its traced end; any other unpaired branch runs straight on along
   its direction to the point closest to the other strands of the zone and ends contact), and run on the design
   seeds (2000 to 2049 and 3000 to 3049, the dense family) with code that was never committed
   (C:/Users/ryane/sc_ribbon_work/v2/logs/dev/): with link.turn_cost_per_deg 1.0 and link.offset_cost 40, turns and
   end-to-end joins went from 2.34 and 2.44 to 2.36 and 2.64 per nucleus (at least 30 percent fewer required), the
   dense S1b SD from 8.09 and 9.02 to 7.98 and 8.62 percent and the CCC from 0.579 and 0.594 to 0.555 and 0.603; at
   the corners of the grid (1.5 and 80, 0.75 and 20) and with the turn of a pair limited to 50 or 30 degrees, turns
   and joins changed by -1 to +13 percent. The re-pairing joined pieces across zones (breaks before the end matching
   4.06 and 4.16 to 3.68 and 3.60 per nucleus) but made switches there (through-switches before the end matching
   0.58 and 0.66 to 0.96 and 1.10), and the turns and joins of the final strands are mostly made after it, by the
   end matching. The junction re-pairing is therefore not proposed and no entry of the round carries it. The
   configuration proposed for this round is therefore v1 with the complete guard at 2.0 FWHM, every other v2 switch
   at its v1 value; the guard is confirmed once now, on the confirmation seeds (11000 to 11099, every registered
   family) against v1, by the targeted metrics and tolerances of entry 23, with the cost rule on
   n2rad51_noHS_herm_000. Z2 is not run: no decomposition change was adopted and no input of S2 depends on the
   strand class, so the S2 inputs of the proposed configuration are those of development trial 4. Step 9 does not
   apply (entry 20: no probe feature reproduced a real symptom), so Z3 is not run either.
25. 2026-10-24, the v2 development round (the configuration proposed for the freeze), written before any snapshot,
   flag-rate or hold-out number of this configuration was computed (no test result had been computed, no test or
   reserved image was run, no test or reserved truth was read and no sealed seed was generated). The complete guard
   of entries 23 and 24 (commit a540680, classes.complete_guard_fwhm 2.0) on the confirmation seeds (11000 to 11099,
   every registered family) against v1 (C:/Users/ryane/sc_ribbon_work/v2/logs/step7_a540680_confirm.json): pooled
   S1e precision 0.883 to 0.938 (at least 0.90 required; met), dense S1e precision 0.672 to 0.741 (at least 0.80
   required; failed), iso complete recall 0.91 and 0.91 (met), and the tolerances of entry 11 (b) held; complete
   recall fell on pairs (0.595 to 0.33), dense (0.078 to 0.036) and snr (0.101 to 0.063). The cost rule held on
   n2rad51_noHS_herm_000 (median time per nucleus 0.154 s against 0.158 s, peak 3.0 GB; complete strands 228 against
   471). The guard is withdrawn by the rule of entry 23: classes.complete_guard is false by default (the switch and
   its tests stay). No change of the round is adopted. The configuration proposed for the freeze is therefore that
   of v1: every v2 switch at the value that gives v1 (V1_VALUES of germquant/sc/straighten/params.py: psf.foci off,
   punctum.enabled false, follow.contact_test local, link.junction_dirs false, link.max_chain_um null,
   classes.complete_guard false), with the package VERSION 2.0 of entry 11 (h); the nucleus columns n_puncta,
   puncta_signal_frac (both 0 with the punctum class off) and steep_length_frac are added, and every v1 column
   equals that of v1 (tests/test_sc_straighten_v2_compat.py). Step 10 of PLAN_V2.md follows with this configuration:
   the snapshot on seeds 12000 to 12099 (this configuration and v1, by their switches, on the same scenes), the flag
   rates of the twelve development gonads (truth-free), the cost and peak memory on the largest 2026-01-05 stack,
   and H1 on the hold-out gonad that has a germquant run (20251021_n2_nohs_HERM_2; theta_full 2.0 um and
   theta_present 10.0 um, those of development trial 4 on the twelve, since the S2 inputs of this configuration are
   those of trial 4), which with one gonad is reported and decides nothing; then the results entry.
26. 2026-10-24, the v2 development round (results), written after the snapshot, the flag rates and H1 of entry 25
   were computed (no test result had been computed, no test or reserved image was run, no test or reserved truth was
   read and no sealed seed was generated). Snapshot (once, seeds 12000 to 12099, commit a711fe5;
   C:/Users/ryane/sc_ribbon_work/dev_eval/summary_a711fe5.json and v2/logs/snapshot_a711fe5.json): the configuration
   proposed for the freeze and v1 gave the same results in every family. Criteria of 8.1: S1a met (isolated ribbons
   up to 60 degrees: mean length error +0.076 um, 95 percent limits -0.069 to +0.221 um, median centreline distance
   0.017 um, one strand per ribbon in all); S1b not met (dense nuclei: mean total error -6.68 percent, SD 8.69
   percent, CCC 0.493; coverage 0.962 and off-truth 0.010 met); S1c met (gap precision 1.00, recall 0.91, gap length
   error 0.034 um); S1d met (90th percentile of the s error 0.064 um, all on the right ribbon); S1e not met
   (complete precision 0.890 pooled; 0.672 in dense nuclei); S1f met (all noise-only nuclei under 1 um). With the
   registered generator the dense totals were -15.3 percent (SD 9.9 percent, CCC 0.221) and the dense complete
   precision 0.775. Development families, descriptive: psf_check complete precision 0.79; of the isolated point
   puncta 83 percent and of the aggregates 78 percent left a strand (the punctum class being off), and the puncta of
   dense nuclei 49.5 percent; real_dense complete precision 0.36. Zone checkpoints: none was run (entries 18 and
   24), because no change that reaches an input of S2 was adopted; the development S2 of the proposed configuration
   is therefore that of trial 4 (S2b 0.987 over twelve gonads, S2a 0.83 over seven, S2c within one point spacing in
   11 of 12), in-sample as PLAN_V2.md 2.5 states (S2c with thresholds chosen on the same gonads). Real data,
   truth-free: the twelve development gonads run at a711fe5 (C:/Users/ryane/sc_ribbon_work/runs/a711fe5) have the
   flag, status and class rates of trial 4 to the fourth decimal (dev_eval/flag_rates_a711fe5.json: 35,840 strands,
   4,740 traced nuclei) and nucleus columns identical to those of runs/6b7be22 (status, n_strands, n_main,
   total_length_um, piece_lw_mean_um, n_pieces, longest_piece_um, gaps_per_10um and on_strand_frac of every
   nucleus); four workers took 0.8 to 2.4 min per gonad, and the main process peaked at 2.5 to 3.0 GB on the
   2025-11-05 stacks and 3.8 to 3.9 GB on the 2026-01-05 stacks (106 x 2720 x 2720). H1 (one line of
   dev_trials.jsonl; runs/h1_a711fe5 and runs/h1_a711fe5_v1; v2/zones/H1_h1_a711fe5/checkpoint.json): only
   20251021_n2_nohs_HERM_2 has a germquant run (entry 16), so H1 decides nothing; with theta_full 2.0 um and
   theta_present 10.0 um of trial 4 on the twelve, S2b was 0.989 under the registered guard and 0.985 under
   guard_frac 0.25 (intensity comparator 0.991 and 0.989); S2a was not evaluable under the registered guard and was
   0.739 under guard_frac 0.25 (intensity comparator 0.670); SC_fullstart fell within one point spacing, 0.98 um
   from the boundary (out of sample); the two configurations gave identical S2 (paired differences 0); the
   sc-straighten process peaked at 3.5 GB. Nothing is changed after H1. The robustness report of section 9, the
   thresholds on all eighteen development gonads and the freeze entry wait for Ryan's decisions (PLAN_V2.md section
   6).
27. 2026-10-24, the v2 development round (an independent review of entries 10 to 26: disclosures, corrections and
   code-only changes; no method change), written before any number of the changes below was computed (no test result
   had been computed, no test or reserved image was run, no test or reserved truth was read and no sealed seed was
   generated). After entry 26 an independent review checked the round (its scripts and outputs are kept in
   C:/Users/ryane/sc_ribbon_work/v2/review/). With its own code it recomputed the dense S1b of the snapshot from the
   v1 package of commit 5cb3cb6 (mean -6.68 percent, SD 8.69 percent, CCC 0.493), S1a of the isolated ribbons
   (+0.076 um, 95 percent limits -0.069 to +0.221 um) and H1 through the guarded loader, all as in entry 26, and
   found every shared column of the twelve runs/6b7be22 and runs/a711fe5 tables identical. Nothing was chosen.
   Disclosed and corrected now, from numbers already computed:
   (a) At 2026-10-24T17:13:02Z, between entries 13 and 14, v2/scripts/v2_zones.py was checked by scoring the trial-4
   runs (runs/6b7be22) against themselves on the twelve development gonads under guard_frac 0.5 and 0.25
   (v2/zones/T4check_check, v2/logs/zones_T4check.log). It reproduced dev_trials.jsonl line 5 and the numbers of
   entries 7 and 10 (S2a 0.828 and 0.802, the intensity comparator 0.965 and 0.910, the residual AUC of entry 7 (c)
   0.616 and 0.609) and the comparators of trial 4 under the registered guard, and it gave for the first time, under
   guard_frac 0.25: S2b 0.984 (intensity comparator 0.986), the S2a AUCs of on_strand_frac 0.857, longest_piece_um
   0.799, n_pieces 0.773 and gaps_per_10um 0.596, and theta_full 2.0 um and theta_present 10.0 um chosen by 5.2,
   with SC_fullstart within one point spacing in 11 of 12 gonads. It compared no configuration and chose nothing,
   and it is neither a checkpoint nor a line of dev_trials.jsonl; but entry 11 (c) says that statistics by zone are
   computed only at the checkpoints, and no entry recorded it.
   (b) The S1b mean total error of v1 met its +/-5 percent band only on the design seeds (-4.9 and -4.8 percent on
   2000 to 2049 and 3000 to 3049, 50 scenes each); on every fresh block of 100 scenes it failed (tuning -6.92,
   confirmation -8.33 and snapshot -6.68 percent), and the pooled S1e precision there was 0.893, 0.883 and 0.890.
   The v1 of the snapshot and the reference configuration of H1 were the code of version 2.0 with V1_VALUES, not the
   v1 commits; tests/test_sc_straighten_v2_compat.py and the review's recompute tie them to commit 5cb3cb6. The
   zones_summary.json files of H1 name the threshold source "chosen on these development gonads", while v2_zones.py
   fixed theta_full 2.0 um and theta_present 10.0 um, those of trial 4, as entry 25 states.
   (c) Order of the round: design-seed runs preceded entries 13 and 21 to 24, and the junction end fit and the
   junction re-pairing of PLAN_V2.md steps 6 and 8 were explored on the design seeds before any entry named them,
   with variants that the plan does not state (entries 23 and 24 record them); the contact stop (entry 20) and the
   end matching (entries 21 and 22) failed on the tuning seeds and were never confirmed. Entries 15, 18 and 22
   restated or dropped targeted metrics after the tuning numbers of their change were seen; the confirmation seeds
   or the tuning seeds decided each, and none was adopted. A later round counts such a revision as a new change.
   Entry 11 and decision 1 of the PR were not approved by the project lead before the round (PLAN_V2.md section 6),
   so no freeze entry is written before the project lead states whether entries 10 to 26 and the later entries of
   this round stand. PLAN_V2.md, part of the registration of entry 11, is outside git; its sha256
   c19b155a03f2fa63d5376b2c181f92932a978a123508b070ca64ba75c7e0be32 has not changed since entry 11
   (v2/logs/prereg_entry11.txt to prereg_entry26.txt).
   (d) B6 of entry 8 ran 20 unregularised Richardson-Lucy iterations on counts that still carried the camera offset,
   with zero-padded borders, and the unadapted pipeline after it; its result (isolated totals 36 times too long) is
   inconclusive about deconvolution as such. The realism probe of entry 19 rendered the brightness modulation and
   the puncta weaker than the real statistics they were calibrated to (the review measured the brightness statistic
   of entry 19 (b) on the phantoms; entry 28), so the skip of PLAN_V2.md step 9 in entry 20 rests on an
   under-calibrated probe.
   (e) The first full test suite of the hand-off of entry 26 ran with the GPU visible, so its small torch tests ran
   on the GPU outside C:/Users/ryane/germquant_golden/gpu_lock.py (one failure, tests/test_net3d_recipe.py, which
   also fails on the code of main with the shared environment and passes with CUDA hidden; the counts of entry 26
   come from a rerun with CUDA hidden); from here on every suite runs with CUDA_VISIBLE_DEVICES=-1. The germquant
   run of 20251021_n2_nohs_HERM_2 peaked at 42.6 GB against the 8 GB per-process budget (entry 16), so the
   2025-10-21 session channel map is checked on HERM_2 only, while protocol 4.5 (5) and entry 11 (e) want both
   development gonads of that session checked before the test gonads 20251021_n2_nohs_HERM_3 and
   20251021_n2_nohs_MALE_3 are run; 20251021_n2_nohs_MALE_2 and the four 2026-01-05 hold-out gonads wait for a
   memory window that the project lead approves. Line 9 of dev_trials.jsonl (H1) holds bare NaN values, which JSON
   does not allow; it is rewritten with null in their place and nothing else changes (the file before is kept as
   v2/logs/dev_trials_before_entry27.jsonl).
   (f) Code-only changes, no method change and no rule of this protocol changed (every output of sc-straighten at
   the defaults stays that of entry 25 apart from the code digest; tests/test_sc_straighten_v2_compat.py and the
   fence tests check it): tests/test_sc_straighten_v2_compat.py fails off the workstation as well when the number of
   strands or the strip index of a synthetic case differs from v1 (before, every difference off the workstation was
   a skip, so CI could not fail on a v1 regression); scripts/benchmark_sc.py zones reports what entry 11 (g)
   registered for the scoring, the per-gonad and mean S2a AUCs of n_puncta, puncta_signal_frac and steep_length_frac
   and the residual AUC of piece_lw_mean_um on the intensity (germquant.benchmark.sc_ribbons.residual_auc, the
   function of v2/scripts/v2_zones.py), with their bootstrap intervals, keeps the threshold tables of 5.2 in
   zones_summary.json when it chooses the thresholds, fixes them outside --final when --theta-full and
   --theta-present are given (source "fixed (command line)"), and writes null for every non-finite number in its
   JSON files and lines; v2/scripts/v2_zones.py passes the thresholds that way and uses the package's residual AUC;
   link.junction_dirs is read as the other switches are (a string "off" is off), and params.resolve refuses a value
   of follow.contact_test other than local or claim_tangent, of link.junction_gap other than anchors or tips, and of
   the on/off switches psf.foci, punctum.enabled, punctum.exclude, punctum.before_linking, link.junction_dirs and
   classes.complete_guard other than true, false, on, off, yes or no; germquant.sc.straighten.synth_dev refuses
   final=True (its families have no sealed seeds before decision 9 of the PR); inputs.discover reads no nuclei.csv
   whose file name or folder name is a sealed image stem or RAD-51 tag (compared without regard to case), and the
   runner refuses such an image as sealed before any of its files is read. No zone statistic is computed with these
   changes before the next checkpoint or the scoring.
28. 2026-10-24, the v2 development round (realism probe recalibrated; diagnostics only, no method change), written
   before any number below was computed (no test result had been computed, no test or reserved image was run, no
   test or reserved truth was read and no sealed seed was generated). Entry 27 made no method change. Motivated by
   the review of entry 27: with the brightness estimator of entry 19 (b) (v2/scripts/realism_measure.py) on the
   tuning seeds 10000 to 10059 and the configuration of v1, real_dense_brightness gave 0.042 (7 strands) in phantom
   nuclei with total_length_um >= 30 um and 0.031 (26 strands) with >= 20 um, against 0.026 and 0.019 for
   real_dense_base and 0.1224 on the trial-4 runs (v2/review/realism_brightness_check.json). Entry 19 (b) used the
   statistic of blurred, measured real profiles as the variance and correlation length of the modulation of the
   unblurred ribbon, so the rendered feature is far weaker than the real statistic; the puncta means of entry 19 (d)
   count rule B1s with the sigma_z of the trial-4 PSF records (mostly the fallback 0.365 um), and B1s labels only
   part of the puncta of the phantoms; and the texture scale of entry 19 (c) was never checked with its estimator on
   phantoms. Whatever this finds, no method change is made in this round (entry 11 (e): nothing is changed after
   H1); a feature that reproduces a real symptom is carried to the next round with its numbers. Nothing here reads
   zone truth.
   Method, fixed now. The development family real_dense of entry 19 takes a calibration argument and can be drawn
   without ribbons and puncta (its families and scenes of entry 19 stay identical); scripts
   v2/scripts/realism_recal.py; the configuration of v1 (the defaults); the tuning seeds 10000 to 10099; four worker
   processes.
   (a) Brightness. For each latent variance v in {0.12, 0.25, 0.5, 1.0} and latent correlation length in {0.1, 0.2,
   0.3} um, the real_dense_brightness scenes with that modulation, and the statistic of entry 19 (b)
   (realism_measure.brightness: in-plane, isolated strands of 3 um or more without gaps; the log of the core profile
   over its median, its variance less the mean noise variance, the median over the strands; the correlation length
   from the slope of the log autocovariance at 0.2 to 1.0 um) over the traced phantom nuclei with total_length_um >=
   20 um (the cut at 30 um of the real data keeps 6 or 7 strands per 60 scenes, and every phantom nucleus holds
   full-length SC, which is what that cut selects among real nuclei); the same at >= 30 um is reported. The pair
   kept is the one whose measured variance and correlation length lie nearest to the real 0.1224 and 0.3126 um by
   sqrt(ln(v_m / 0.1224)^2 + ln(l_m / 0.3126)^2); a pair with fewer than 10 strands, a variance at or below 0 or no
   correlation length is not eligible, and with no eligible pair the calibration of entry 19 stays.
   (b) Texture. Scenes of the real_dense geometry, nucleoplasm, texture and noise without ribbons and puncta, with
   the latent texture variance of entry 19 (c) at 0, 0.1082 and 0.2164 (relative to the squared nucleoplasm level);
   per scene the statistic of entry 19 (c) (realism_measure.texture: the variance of the counts smoothed with the
   PSF of the image record over the territory of the label, reach 0.4 um, eroded by two voxels, less the variance
   that the image's fitted noise model gives after the same filter, over the squared mean count above the camera
   offset of 100), the median over the scenes of each level. A least-squares line through the three medians gives
   the latent value whose statistic equals the real 0.1082, clipped to [0, 0.5]; with a slope at or below 0 the
   value of entry 19 stays.
   (c) Puncta. The B1s counts per traced nucleus of the twelve trial-4 runs (realism_measure.b1s_counts on
   v2/steep/psf) recomputed with each image's compact-foci sigma_z (the values of entry 19 (a)) in place of that of
   its PSF record, their means by the terciles of SYP brightness of entry 19 (d), each divided by the B1s yield of
   the phantoms: the punctum-class strands per rendered punctum, pooled over the real_dense_puncta scenes run with
   punctum.enabled true and punctum.exclude false (the label of entry 18; there the image record has the scene's
   PSF).
   (d) The probe. Families real_dense_r2 (all four features) and real_dense_r2_base, _psf, _brightness, _texture and
   _puncta: the scenes of entry 19 with the calibration of (a) to (c) and the PSF feature of entry 19 (a), run once
   on the tuning seeds (v2/scripts/v2_s1_panel.py, name probe_r2) and reported with the statistics of entry 19 and
   its rule: a feature reproduces a real symptom when switching it on alone moves the phantom statistic at least
   half of the way from its base value to the real one; the same for all four features together, and their S1b, are
   reported, descriptive.
29. 2026-10-24, the v2 development round (results of entry 28; diagnostics only, no method change), written after
   the numbers below were computed (no test result had been computed, no test or reserved image was run, no test or
   reserved truth was read and no sealed seed was generated). Computed as entry 28 fixed it on the tuning seeds
   10000 to 10099 with the configuration of v1 (commit 9726612 for (a) to (c), e53c814 for (d);
   C:/Users/ryane/sc_ribbon_work/v2/logs/realism_recal_a_9726612_tuning.json, realism_recal_b_9726612_tuning.json,
   realism_recal_c_9726612_tuning.json, realism_probe_r2_e53c814_tuning.json and realism_probe_r2_mechanism.json).
   (a) Brightness. The statistic on the phantoms was noisy and not monotonic in the latent values (8 to 43 strands
   per pair at total_length_um >= 20 um, 0 to 10 at >= 30 um; measured variances 0.034 to 0.131); at the latent
   variance of entry 19 (0.12) the phantoms measured 0.034 to 0.042, as the review found. The nearest pair was
   latent variance 0.5 and latent correlation length 0.1 um (measured 0.102 and 0.315 um from 18 strands; log
   distance 0.18; the next was 1.0 and 0.1 um, 0.31).
   (b) Texture. Without texture the phantom statistic was already 0.141 (median over 100 scenes), above the real
   0.108, and 0.266 and 0.377 at latent 0.108 and 0.216 (slope 1.09), so the latent texture is 0: by this estimator
   the texture of entry 19 made the phantoms about 2.5 times as textured as the real nuclei, and the real excess
   variance does not exceed what the phantoms' geometry alone gives (the blurred edge of the nucleoplasm inside the
   territory).
   (c) Puncta. With each image's compact-foci sigma_z, the B1s counts per traced nucleus rose from 1.069, 0.869 and
   0.634 (entry 19 (d), reproduced exactly with the records' sigma_z) to 1.478, 1.223 and 0.896 by tercile of SYP
   brightness; on real_dense_puncta B1s labelled 38 strands for 89 rendered puncta (yield 0.43, far below the 72 to
   86 percent of isolated puncta of entry 16), so the recalibrated means are 3.46, 2.86 and 2.10 puncta per nucleus.
   (d) The probe (real_dense_r2 families, 100 scenes each; real values 0.403 for the steep share, 2.49 um for the
   median in-plane strand, 2.63 um for the median piece_lw_mean_um and 2.19 for main strands over the expected
   count, all in traced nuclei with 30 um of SC or more). The recalibrated brightness modulation alone moved the
   median piece_lw_mean_um from 5.03 to 3.65 um, 58 percent of the way to the real value, so by the rule of entry 19
   it reproduces that real symptom; it moved the median in-plane strand 17 percent of the way (5.36 to 4.86 um), the
   steep share 7 percent and main strands over the expected count 4 percent. Its mechanism is detected gaps in
   continuous ribbons, not broken strands: gaps per traced nucleus rose from 1.76 to 4.29 (median gaps_per_10um 0.47
   to 1.46) and pieces from 7.0 to 9.9, while strands per nucleus went from 5.26 to 5.60 and the median total from
   28.8 to 27.8 um (the dense ribbons hold no true gap). The recalibrated puncta moved the steep share 28 percent of
   the way (0.174 to 0.239), main strands over the expected count 16 percent and piece_lw_mean_um 16 percent; the
   texture (latent 0) left the scenes those of the base; the PSF feature was that of entry 19 (at most 12 percent).
   All four features together moved the steep share 68 percent of the way (0.330) and piece_lw_mean_um 73 percent
   (3.27 um), the in-plane strand 35 percent (4.37 um) and main strands over the expected count 31 percent (1.32);
   their dense S1b was mean -5.4 percent, SD 16.0 percent and CCC 0.184 (entry 19: -3.2 percent, 18.3 percent and
   0.133). So the skip of PLAN_V2.md step 9 in entry 20 was not supported: with the probe calibrated as entry 28
   fixed, the brightness feature reproduces the low continuity of real full SC, through false gaps at dim stretches,
   and the likely form of step 9 named by the plan (a weak-run and gap rule that keeps a ribbon through dim
   stretches) would have been the change to make. As entry 28 fixed, no method change is made in this round (entry
   11 (e)); the finding is carried to the next round, where such a rule is designed on these phantoms and judged by
   the rules of a new entry, and it qualifies S2a in the meantime: on development data part of the low
   piece_lw_mean_um of full SC may be dim stretches read as gaps. The robustness report of section 9, the thresholds
   on all development gonads and the freeze entry wait for the project lead's decisions.
30. 2026-10-24, the v2 development round (the CI check of entry 27 (f), departure; no method change), written before
   any number of the change below was computed (no test result had been computed, no test or reserved image was run,
   no test or reserved truth was read and no sealed seed was generated). CI at commit b1610c7: the Linux and Windows
   jobs passed; both macOS arm64 jobs failed one case of tests/test_sc_straighten_v2_compat.py, gaps/entry2/2004,
   whose strip index differed from v1 while its number of strands (2) did not. On the Linux and macOS runners the
   table digests of most cases differ from those of the workstation (the floating point of another platform; on the
   Windows runner they all match), and no case differed in its number of strands on any runner. A strip runs every
   0.05 um along its strand, so a last-digit change of an end can add a column, and an end refitted from deeper
   inside can move by more than that. The structural check of entry 27 (f) is therefore the number of strands of
   every synthetic case alone, on every platform; the strip index stays with the table digests among the workstation
   pins (GERMQUANT_WORKSTATION_PINS=1), which the hand-off runs on the workstation.
31. 2026-10-30, the v3 development round (its rules, and the disclosures of its interrupted start), written before
   any number of the v3 changes was computed on the tuning, confirmation or snapshot seeds of this round, on a
   development image or on zone truth (no test result had been computed, no test or reserved image was run, no test
   or reserved truth was read and no sealed seed was generated). Motivated by development results (entries 8, 9, 16
   to 19, 28 and 29): rule B1 labels most PSF-blurred SYP puncta and small aggregates of the development phantoms
   and almost no tilted ribbon longer than about 1 um, but on real images its reference, the image PSF, is mostly
   the fallback sigma_z (0.365 um, against 0.380 to 0.464 um from the compact RAD-51 foci of entry 9), and its
   confirmation failed on one strand of a ribbon (entry 19); brightness variation along the ribbons at the level
   measured on the trial-4 strands reproduces, on the recalibrated phantoms, 58 percent of the fall of the median
   piece_lw_mean_um from continuous phantoms to real full SC, through gaps read at dim stretches (entry 29). The
   round makes two changes, each behind switches whose defaults give the outputs of v2 (which equal those of v1 in
   every column of v1): the punctum class B1 with the image's measured axial sigma (entry 32) and dim stretches that
   are not gaps (entry 33). The long vertical structures that hold most of the steep length in the heat-shock and
   no-SC zones are not addressed: without hand traces (protocol R2) and a bead PSF (decision 9 of the PR) they
   cannot be told from real tilted SC, so they are reported (the steep share by zone), not changed.
   (a) The interrupted start (resume note C:/Users/ryane/sc_ribbon_work/v3/RESUME_NOTE.md). A first build of this
   round began on 2026-10-25 (UTC) after entry 30 and stopped at about 01:55Z the same day, before any entry of the
   round was committed; its draft of this entry (v3/entries/entry31.txt) was never committed and is replaced by this
   one. Before any entry, with code that was never committed (sc_straighten code digest d143d3b93f68), it computed
   on development seeds: (i) at 00:45Z to 00:48Z, v3/scripts/design_gaps.py on the design seeds 2000 to 2049 and
   3000 to 3049 at the defaults (the features of every detected gap of five families, against the truth gaps of the
   gaps family; v3/design/gaps_design_v1.json); (ii) at 01:12Z to 01:51Z, five synthetic panels on the design seeds
   2000 to 2049 (17 families each): the defaults, the punctum class of entry 16 with punctum.exclude, the punctum
   rule of entry 32 as then coded, the dim-stretch rule of entry 33 as then coded, and the punctum rule without its
   shape test (stopped at 679 of 850 scenes), and two reports of them (v3/design/report_d_v2.json and
   report_design4.json). All of it is discarded: no target, parameter, decision or table of this round uses it, and
   the resumed round did not open it. The values that the uncommitted code gave the dim-stretch rule after (i)
   (gaps.dim_fill 0.15, gaps.dim_z 2.0) are not kept as such: entry 33 chooses them on the tuning seeds by a rule
   fixed there; the other new constants are set in entries 32 and 33 from earlier entries. The golden
   tests/golden/sc_straighten/v2_compat.json made then (with the package code of acbee94, as its builder required)
   was recomputed on 2026-10-30 from a clean copy of acbee94 with identical digests in all 19 cases and the code
   digest deedb8a7e9fd (v3/logs/resume_verify_v2_compat_golden.json); it is kept. Before any number, the resumed
   round added to the uncommitted code the guards and the measured axial sigma of the dim-stretch rule and the
   development family gaps_psf (entry 33), and tests.
   (b) Seeds. Development seeds used before this round: 1000 to 1023, 2000 to 2049, 3000 to 3049, 5000 to 5459, 5500
   to 5959, 6000 to 6139, 6300 to 6359, 10000 to 10099, 11000 to 11099 and 12000 to 12099 (2000 to 2049 and 3000 to
   3049 again by (a)). This round tunes on 13000 to 13099, confirms on 14000 to 14099 (once per change and per
   revision) and computes the snapshot of the configuration it proposes once, on 15000 to 15099, with v2 on the same
   scenes; no run has touched 13000 to 15099, and no design-seed run is made in this round after this entry. Scripts
   and outputs: C:/Users/ryane/sc_ribbon_work/v3/ (scripts/v3_panel.py, v3_targets.py, v3_zones.py, v3_cost.py;
   panel/, logs/).
   (c) Adoption. A change is adopted when its targeted metrics hold on the tuning seeds and at their confirmation
   levels on the confirmation seeds, the tolerances of entry 11 (b) hold on the confirmation seeds against v2 on the
   same scenes (scored over the strands the nucleus tier counts, so punctum strands are out of the S1 scoring when
   they are out of the tier), the cost rule of entry 11 (b) holds on n2rad51_noHS_herm_000, and the checkpoint of
   (d) does not veto it. A parameter that a change introduces is set from physics or earlier entries where its entry
   says so, or chosen on the tuning seeds by the rule its entry states, never on development images or zone truth. A
   change gets at most one revision, with its own entry before any number of the revision (a new change in the sense
   of entry 27 (c)); otherwise it is withdrawn by its switch. Whatever is adopted, the package defaults stay those
   of v2 in this round: an adopted change enters the configuration proposed for the freeze by its switch values
   (params.V3_PUNCTUM, params.V3_GAPS), and a default changes only by a later entry after the project lead's
   decisions 2 and 3 of the PR.
   (d) Zone truth enters once, at checkpoint Z-v3, after both changes are decided on synthetic evidence:
   sc-straighten runs on the twelve development gonads with trial-4 runs with v2 (the defaults), with each change
   alone and with both, whatever the synthetic stage decided (so that the real-data effect of each rule is
   reported), and S2 is computed by the registered rules (theta_full and theta_present chosen again by 5.2 for each
   configuration) under guard_frac 0.5 and 0.25, with the paired differences against v2 (cluster bootstrap, section
   7) and the residual AUC of piece_lw_mean_um after the intensity (entry 7 (c)) with its paired difference. The
   checkpoint can only veto (entry 11 (c)): a change that passed the synthetic stage is withdrawn when, in the
   configuration with that change alone, the paired S2a difference against v2 is below -0.03 under both guards or
   the S2b mean is more than 0.01 below that of v2; when both changes pass and neither is vetoed alone, the
   configuration with both is checked by the same rule, and if it is vetoed, the change with the lower paired S2a
   difference alone under guard_frac 0.25 is withdrawn. No parameter is chosen and nothing is revised after the
   checkpoint. Each configuration is one line of dev_trials.jsonl (four lines; with H1-v3 the budget of section 9
   then holds 15 of 40).
   (e) Hold-out H1-v3, after the checkpoint: v2 and the configuration this round proposes run once on the
   development hold-out gonads that have a germquant run (20251021_n2_nohs_HERM_2 and those of (f)); S2 with the
   thresholds of the twelve of each configuration (S2c out of sample), the paired differences and the residual AUC;
   one line of dev_trials.jsonl per configuration. With three or more gonads evaluable under guard_frac 0.25, a
   paired S2a difference below -0.05 withdraws the changes of the proposal (v2 is proposed instead); with fewer it
   is reported and decides nothing. If no change is adopted, H1-v3 runs v2 and both changes and decides nothing.
   Nothing is changed after H1-v3.
   (f) The germquant runs (xy stride 1) of the five development hold-out gonads that entry 27 (e) left for a memory
   window above the 8 GB budget (20251021_n2_nohs_MALE_2, 20260105_n2_hs_1_HERM_1, 20260105_n2_hs_HERM_2,
   20260105_n2_nohs_HERM_1 and 20260105_n2_nohs_MALE_3; entry 11 (e), protocol 4.5 (5)) are made under a gate that
   the orchestrator set for this round: one at a time; each started only when Windows reports at least 40 GiB of
   free physical memory right before its start (PowerShell Get-CimInstance Win32_OperatingSystem FreePhysicalMemory,
   read inside the GPU lock; otherwise the lock is released and the gate asked again 10 minutes later); the GPU only
   through C:/Users/ryane/germquant_golden/gpu_lock.py; the peak memory of each run logged
   (v3/scripts/heavy_germquant.py and mem_gate_run.py, v3/logs/heavy_germquant.jsonl). No approval of decision 4 of
   the PR (a memory window above 8 GB) by the project lead is on record; the gate is the orchestrator's. State when
   this entry was written: 20251021_n2_nohs_MALE_2 (2026-10-25, 25 min, peak working set 45.0 GiB, peak commit 54.5
   GiB) and 20260105_n2_hs_1_HERM_1 (27 min, 50.3 and 65.6 GiB) finished and passed the completion, channel and
   segmentation checks (central_element is the nd2 channel 545 and 477, the sessions' SC channels; Cellpose; 585 and
   532 nuclei), which completes the 2025-10-21 channel-map check of protocol 4.5 (5) on both development gonads of
   that session; 20260105_n2_hs_HERM_2 was interrupted by the stop and left no output, and was started again on
   2026-10-30 at 12:08Z (43.2 GiB free), before this entry; the other two follow. The gate guards the start only:
   during the two finished runs the free physical memory fell to 0.33 and 0.01 GiB. No sc-straighten run is made on
   these gonads before H1-v3, and their outputs are looked at only by those checks.
   (g) Not changed in this round: the S1 criteria of 8.1, the S2 rules, the registered families and sc_split.json,
   the QC thresholds and every switch of v2; the package VERSION stays 2.0 (the record's parameters name the rules).
   With the switches of entries 32 and 33 off, every column of every table, the strips and the classes equal those
   of the v2 head (commit acbee94; tests/test_sc_straighten_v3_compat.py against
   tests/golden/sc_straighten/v2_compat.json), and every v1 column equals v1
   (tests/test_sc_straighten_v2_compat.py); the params_hash changes with the new keys. No output column is added;
   the record gains the keys punctum_reference, dim_reference and dim_stretches only when their switches are on.
   (h) Reporting added under rule 2 (secondary, never gating), at the checkpoint and at the scoring: the steep share
   of the counted traced length by zone, and the puncta, gaps and dim stretches per traced nucleus by zone.
32. 2026-10-30, the v3 development round (the punctum class with the image's measured axial sigma), written before
   any number of this change was computed on the seeds of entry 31 (b), on a development image or on zone truth (no
   test result had been computed, no test or reserved image was run, no test or reserved truth was read and no
   sealed seed was generated). Motivated by development results (entries 8, 9 and 16 to 19): rule B1 at 1.3 times
   sigma_z labelled 83 to 86 percent of PSF-blurred point puncta, 72 to 73 percent of aggregates of 0.25 to 0.8 um
   and none of 60 tilted arcs of 1.5 to 6 um (entry 16), but the confirmation of entry 18 labelled one strand of a
   gaps ribbon (a 1.26 um vertical strand on a ribbon tilted 28 degrees), and its reference on real images, the
   image PSF, fell back to 0.365 um on 10 of the 12 development images while their compact RAD-51 foci give 0.380 to
   0.464 um (entry 9), so that on wide-PSF images the rule misses puncta.
   Method (switches; params.V3_PUNCTUM: punctum.enabled true, punctum.exclude true, punctum.sigma_source foci_depth,
   punctum.one_gaussian true; the defaults are those of v2). After measurement, a traced strand with frac_steep60
   above 0.5 and length_xy_um at most punctum.max_xy_um (0.5 um) whose core profile over s from -0.5 um to length_um
   + 0.5 um (every 0.05 um), fitted by least squares with a Gaussian plus a constant (entry 16), (i) is one
   Gaussian: the fit exists, its sigma lies more than 1 percent inside its bounds (0.05 and 5.0 um), its centre lies
   on the strand (0 <= mu <= length_um), its amplitude is positive and the RMS of its residuals is at most
   punctum.max_rel_rms (0.2, the threshold of the foci fits of entry 12) of the amplitude, and (ii) has sigma times
   tz_abs_mean at most punctum.sigma_ratio (1.3, entries 8, 9 and 16) times the reference axial sigma at the
   strand's mean depth, is of class punctum with the strand flag punctum (info). The reference
   (psf.axial_reference): the foci estimate of entries 12 and 13 on the image's compact isolated RAD-51 foci (the
   fit with the planar background and its rules; at least psf.foci_min, 15, kept fits; both emission wavelengths),
   the median kept sigma_z times the ratio of the central_element to the foci emission wavelength; with at least
   psf.ref_depth_min_fits (30, twice psf.foci_min) kept fits whose depths (5th to 95th percentile) span at least
   psf.ref_depth_min_span_um (3.0 um, over which the depth dependence of entry 9, about 20 percent per 10 um,
   changes sigma_z by about 6 percent, twice the per-image error of the foci estimate in entry 16), a depth term:
   the Theil-Sen slope of the kept fits' sigma_z on their depth (times the same ratio), evaluated at the strand's
   depth clipped to the 5th to 95th percentile of the foci's depths, the value clipped to the 10th to 90th
   percentile of their scaled sigma_z; without a foci estimate the reference is the image PSF's sigma_z (that of v1,
   as in entry 16). The image PSF itself stays that of v1 (psf.foci off); the record names the reference
   (punctum_reference). Punctum strands stay in the strand, sample and strip outputs and leave the nucleus tier as
   under punctum.exclude of entry 16 (every total, count, piece and continuity column, the multiplicity references,
   the strand tubes of on_strand_frac and coverage_proxy, the foci assignment); n_puncta and puncta_signal_frac
   report them; the S1 scoring leaves them out (entry 16). The end matching is that of v1 (punctum.before_linking
   false).
   New development family puncta_psf (family code 109; not registered; the seeds of entry 31 (b)): one image of two
   layers 6.5 um apart in z, each with three dense nuclei of 2.3 holding 3 puncta each, one nucleus with 6 puncta
   and no ribbon, and two steep isolated arcs (tilt U(60, 90) degrees, length U(0.9, 4.0) um, peak U(300, 800)
   counts) in their own nuclei; puncta half points and half spheres of diameter U(0.25, 0.8) um, peaks log-uniform
   in [100, 1500] counts, at least 0.6 um from every ribbon and 1.2 um apart; the axial PSF the psf_check level
   (0.33, 0.38, 0.43 or 0.48 um by seed modulo 4) at the middle depth with a relative depth gradient U(-0.02, 0.02)
   per um, the lateral sigma U(0.135, 0.165) um; a foci channel of RAD-51-like spots in every nucleus but the arcs'
   (as psf_check, rendered with the local PSF times 595/511) with its spot table and the wavelengths 511 and 595 nm.
   Targeted metrics, on the tuning seeds and then once on the confirmation seeds, against v2 on the same scenes
   (v3/scripts/v3_targets.py): (a) yield: at least 80 percent of the point puncta of puncta_pt and 65 percent of the
   aggregates of puncta_agg leave no strand outside the punctum class (entry 16 (b); there the reference is the
   image PSF), and on puncta_psf at least 65 percent of the puncta at every PSF level; (b) reference: on puncta_psf
   the median relative error of the reference at the middle depth against the rendered sigma_z there lies within
   +/-5 percent at every level; (c) ribbons: no strand lying half or more within the tolerance of an iso, gaps, foci
   or pairs ribbon is of class punctum and no matched strand of such a ribbon is (entries 16 (c) and 19); the
   punctum-labelled strand length lying within the tolerance of the ribbons of dense, dense_reg, snr, snr_reg,
   dense_puncta and real_dense_r2, and of the dense ribbons of puncta_psf, is at most 0.5 percent of that family's
   truth length; on puncta_psf at most 1 percent of the steep arcs of 2.0 um or more and at most 5 percent of those
   of 1.5 to 2.0 um carry a punctum strand lying half or more on them or as their matched strand (shorter arcs are
   reported); (d) excess: the excess of the mean nucleus total of dense_puncta over that of dense_puncta_clean falls
   by at least 30 percent (confirmation: 15 percent; entry 16 asked 60 percent and reached 36 percent without the
   shape test); (e) the tolerances of entry 11 (b) on the confirmation seeds and the cost rule. The confirmation
   levels of (a) to (c) are the same. Reported, choosing nothing: the reference's error at the puncta's depths with
   the depth term and with the constant estimate, the sources of the reference by level, and, at the checkpoint, the
   reference of the twelve development images against their image PSF.
