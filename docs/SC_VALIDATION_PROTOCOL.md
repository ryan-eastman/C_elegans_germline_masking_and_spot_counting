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
