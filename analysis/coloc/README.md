# SYP-3 / P-granule colocalization analysis (ccw77)

Analysis layer for the SYP-3 (mCherry knock-in) partitioning into PGL-1::GFP P granules under heat
shock, built on top of the germquant pipeline outputs (Cellpose nuclei + germline isolation).

## Layout

| path | what |
|---|---|
| `chload.py` | canonical nd2 loader; asserts channel order (405 DAPI / 477 PGL-1 / 545 SYP-3 / 640 lamin) and voxel size on every read; `parse_iid`, `is_excluded` |
| `exclusions.json` | single source of truth for gonads excluded after the per-gonad mask-contamination audit |
| `scripts/pc_lamin_worker.py` | whole-gonad partition coefficient, lamin-envelope masking (+ DAPI masking for comparison), rotation and z-shift controls |
| `scripts/trace_pachytene.py` | pop-up tool to draw the pachytene region per gonad (scroll = z-planes, m = max projection); writes `staging/zones` |
| `scripts/rezone_all.py` | re-derive zones from saved traces after any change to the zoning rule |
| `scripts/pc_zone_worker.py` | early / mid / late pachytene partition coefficients from the traced zones |
| `scripts/granule_tail.py`, `scripts/tail_by_distance.py` | per-granule SYP-3 excess: fraction of granules that "light up", by distance from the envelope, with GFP bleed-through calibration (v4 drops no-envelope labels first) |
| `scripts/amount_metrics.py`, `denominator_check.py`, `exposure_sim.py` | method robustness: absolute-amount metrics, cytoplasm-denominator and background variants, synthetic re-exposure |
| `scripts/pubstyle.py`, `scripts/fig_pub_*.py` | publication figures (Arial, 180 mm, Okabe-Ito, PDF + 600 dpi PNG) |
| `scripts/qc_mask_audit.py`, `qc_mask_zoom.py` | per-gonad nucleus-mask audit: all-depth overlay, lamin-only nucleus detection (missed rings), no-envelope test for every label; zoom with flagged ids |
| `scripts/nucleus_filter.py`, `pc_filter_compare.py`, `build_filter_compare.py` | per-nucleus filters (no lamin envelope; outside the traced territory) and the three-way whole-gonad sensitivity table `results/pc_filter_compare.csv` (staged rows untouched) |
| `scripts/fig_pub_pachytene_pooled.py`, `fig_pub_partition.py` | figure 1 (primary): partition coefficient within the pooled hand-traced pachytene region (`pach_*` columns of `pc_zone_all.csv`); figure 1S: the whole-gonad version (all `in_germline` labels, supplementary) |
| `scripts/fig_grant_assets.py`, `fig_grant_pipeline.py` | grant pipeline figure: real crops + editable SVG (`figures/grant/pipeline_figure.svg`, 180 x 48 mm artboard, spare parts outside it) |
| `scripts/make_staged_review.py` | builds the pre-commit review page from the staged diff |
| `scripts/crescent_axis.py` | shared data-location and crop-cache helpers (`find_run`, `find_nd2`, `load_crops`, `lamin_labels`) used by every worker; the rest of the file is the automated staging attempt, validated not accurate enough (~5 rows error) and superseded by hand tracing |
| `results/` | per-gonad tables (acquisition metadata, PC tables, zone rows, robustness checks) |
| `staging/` | hand traces (`pachytene_traces.json`), per-nucleus zones, human-verified landmarks |
| `figures/` | publication figures (`figpub1` pooled pachytene, `figpub1S` whole gonad, `figpub2..6`), the grant pipeline figure (`grant/`), and the segmentation-validation and robustness figures (`fig71`, `fig72`) |

## Key results (clean 13 gonads, lamin masking)

* Heat shock raises the granule-specific SYP-3 partition coefficient in males (pooled hand-traced
  pachytene region 1.085 v 1.311, P = 0.016, figure 1; early pachytene P = 0.016, figure 2; whole gonad
  P = 0.032, figure 1S); hermaphrodites are underpowered (n = 2 v 2, floor P = 0.333), not null.
* The metric that matches the images is the fraction of P granules holding SYP-3 above half the nuclear
  level (no-envelope objects removed, `granule_tail_v4.csv`): 22% in HS males, 6% in unshocked males
  (P = 0.032), 14% v 1% in hermaphrodites (n = 2 v 2).
* Imaging session is the dominant covariate (the 8-Jul session amplified the HS response in both sexes);
  exposure, background choice, cytoplasm-denominator choice and GFP bleed-through were each tested and
  ruled out as drivers.

## Mask audit (2026-08-24)

The pipeline's `in_germline` flag is DAPI-only. A per-label audit of the clean 13 gonads
(`results/mask_audit/`, 13 independent visual reviews in `visual_review.json`) found two kinds of
non-germline object inside the whole-gonad nucleus set:

* objects with no LMN-1 envelope (sperm / spermatid masses, condensed debris): 26% of nuclear volume in
  noHS_male_003, 23% in HS_male_010, 17% in HS_male_009, 8% in HS_herm_011, 5% in HS_male_008, under
  2.5% elsewhere. The no-envelope test (`nucleus_filter.py`) is reliable in the proximal/sperm regions
  (0 to 2 borderline false removals per gonad) but flags 6 and 16 dim distal mitotic nuclei in
  HS_herm_08 and noHS_herm_004.
* ring-bearing somatic nuclei (seminal vesicle / vas deferens, intestinal, spermathecal) outlined as
  germline in 9 of 13 gonads; a lamin test cannot separate these, and the territory filter only removes
  blocks physically disconnected from the tube (3 males).

The stage-resolved results (hand-traced zones) are unaffected: 4 no-envelope objects received a zone
across all 13 gonads. The whole-gonad tables are the exposed ones. Sensitivity (male, noHS n = 4 v HS
n = 5, granule-specific PC): all labels 1.158 v 1.333 (P = 0.032); no-envelope removed 1.152 v 1.336
(P = 0.032); plus off-trace territories removed 1.141 v 1.337 (P = 0.032); pooled hand-traced pachytene
region 1.085 v 1.311 (P = 0.016). Hermaphrodites (n = 2 v 2) move by +0.18 to +0.22 in every version
(P floor 0.333). The pooled pachytene region (`pach_*` columns) is therefore figure 1, and the whole-gonad version is
figure 1S. Figures 3 to 5 (imaging-session covariate) stay on the whole-gonad values because they
also show the untraced, excluded gonads for context. The per-granule lit-fraction table and figure 6 now use
the no-envelope-filtered envelope set (`granule_tail_v4.csv`; unfiltered it was male noHS 8.2% v HS
24.7%, filtered 5.8% v 21.9%, P = 0.032 both).

## Caveats

* Scripts use absolute paths under `C:\Users\ryane\` (data on the local SSD / NAS); adjust
  `RUN_DIRS`, `E_DIRS` and `CA` when running elsewhere.
* `staging/cache` (per-gonad crops, ~60 GB) is not committed; workers rebuild it from the nd2s.
* HS and no-HS gonads of both sexes should be imaged in one session at fixed exposures before the sex
  comparison is made (see `PIPELINE_AND_ACQUISITION_SPEC.md`).
