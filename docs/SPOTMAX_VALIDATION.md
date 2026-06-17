# SpotMAX spot-counting: integration + Imaris validation

**Date:** 2026-06-17 · **Pipeline:** `segment → isolate germline → SpotMAX spots → Imaris`
**Test image:** `20251105_N2_nohs_HERM_001` (N2 control, 2600×2600×71, voxel 0.108×0.108×0.2 µm)

This documents the switch from the v1 `blob_log` foci detector to **SpotMAX**, the cross-gonad
calibration of the detector against your Imaris ground truth, and the honest limits of that
calibration. Every headline was independently re-derived and **adversarially verified** by two
multi-agent reviews that caught and corrected several of my own conclusions (see §5).

---

## 1. What changed

- New module **`germquant.spots.detect_spots`** wraps SpotMAX's validated chain (preprocess →
  semantic segmentation → peak detection → 160 per-spot features incl. effect size) into tidy
  per-spot + per-nucleus tables. Runs **headless** (no GUI/prompts).
- `pipeline.py` gained a gated **`spots`** stage that replaces `foci`; schema gained a `spots`
  table + `n_spots`/`mean_spots`; the montage overlays spots.
- `config.yaml` defaults to the PI's minimal path: **`spots` on; `foci`/`sc`/`zones` off** (via
  flags — nothing deleted, all recoverable by flipping `enabled: true`).
- Segmentation still uses **your trained model** `germline_nuclei_combined` (not stock cpsam).

## 2. Validation vs Imaris (same image)

The `.ims` for this exact gonad (`E:\Madeleine\N2\20251105_N2_noHS\…HERM _001.ims`) is the ground
truth: **146 nuclei (Surfaces) + 1222 RAD-51 spots (Spots)**. It is the same gonad as
`20251105_N2_nohs_HERM_1.xlsx` (146/146 surfaces agree to ~1e-3 µm in position and volume — not a
coincidental count tie). The lab's canonical per-nucleus number for it (xlsx coloc method) = **5.81**.

> **Important caveat:** Imaris masked only a *subset* of nuclei (146), "enough to get usable data".
> So nucleus **count, precision, and F1 are meaningless** — our extra nuclei are real, not errors.
> Only **recall of the masked subset** and **per-nucleus RAD-51** are valid comparisons.

| Metric | Result |
|---|---|
| **Segmentation recall** of the 146 Imaris nuclei | **1.000** (every one found within 1.5 µm) |
| **Per-nucleus RAD-51 correlation** (ours vs Imaris) | **Pearson r = 0.864**, Spearman ~0.88, p<0.0005 |
| **Total spots** after the z-fix (§3), confirmed end-to-end | **1262 ≈ Imaris 1222** (+3.3 %) |
| **Per-nucleus RAD-51** (146 matched, de-duped) | ours 4.66 — between Imaris's own inside-sphere 2.92 and coloc 5.81 |

The whole-pipeline run with the z-fix produced **617 nuclei → 595 germline, 1262 spots, qc_pass**
(down from 2033 pre-fix; max foci/nucleus 33 → 21 as z-splits collapsed). The per-nucleus number
sits between Imaris's two assignment methods — their ~2× disagreement (inside-sphere 2.92 vs coloc
5.81 *for the same gonad*) is now the dominant uncertainty, not our detector. Hence §5.

## 3. The over-count was Z-axis spot-splitting (not dim noise) — FIXED

Uncalibrated, SpotMAX detected 2033 spots (~1.6× Imaris). The cause is **not** spurious dim peaks —
it is **z-axis spot-splitting**: the confocal axial PSF (~0.6–0.8 µm) is far wider than the z
spot-radius (1.5 voxels = 0.3 µm), so one RAD-51 focus is detected as **2–3 stacked peaks at the
same (x,y) voxel**. 60.7 % of in-nucleus spots share an exact (x,y) column.

**Fix:** `detect_spots(..., merge_z_columns=True, z_merge_gap_um=0.8)` collapses peaks sharing an
(x,y) voxel within `z_merge_gap_um` (keeps the brightest), while leaving genuinely distinct foci
that are farther apart in z. Result: **2033 → 1262 spots ≈ Imaris's 1222** (confirmed end-to-end). On by default.

(`effect_size_min` does **not** subsume this: z-split duplicate peaks are *bright* — they survive an
effect-size cut — so the merge is required. Verified on HERM_001: merge on → 5.94/nuc; merge off →
9.35/nuc, ~1.57× over.)

## 4. Cross-gonad validation — the detector needed a threshold change (the central result)

Tuning on one gonad is exactly the trap, so we ran the full pipeline on **all 6 N2 no-HS gonads**
(20251105; 3 HERM, 3 MALE) that have Imaris GT, matched each gonad's traced surfaces to our nuclei,
and compared to the lab's coloc per-nucleus counts. This **overturned the single-gonad picture**:

- **Segmentation is robust** — recall 0.96–1.0 across all 6 gonads.
- **Spot detection with SpotMAX's default per-image otsu threshold is NOT** — it *starves* some
  images (HERM_002 → 1.04/nuc vs coloc 4.66) and *floods* others (MALE_002 → 109.7/nuc vs 8.15).
  Pooled cross-gonad CCC ≈ **0.04**. otsu adapts to each image's intensity histogram and swings the
  wrong way per image.

Sweeping thresholding method × `effect_size_min` over all 6 gonads, the robust choice is
**`threshold_triangle` + `effect_size_min = 3.0`** (`threshold_li` over-floods):

| gonad | ours (triangle/es3) | Imaris coloc |
|---|---|---|
| HERM_ | 4.98 | 3.72 |
| HERM_001 | 5.94 | 5.81 |
| HERM_002 | 3.74 | 4.66 |
| MALE_ | 6.80 | 7.57 |
| MALE_001 | 5.86 | 4.35 |
| MALE_002 | 9.46 | 8.15 |

| metric | triangle/es3 | otsu/es3 | otsu/es0 (old default) |
|---|---|---|---|
| pooled per-nucleus CCC (n=712) | **0.483** | 0.441 | 0.040 |
| **gonad-mean CCC** | **0.805** | 0.711 | 0.040 |
| gonad-mean MAE | **0.98** | 1.35 | — |
| pooled bias | +0.39 | −0.45 | +18.7 |

### Honest caveats (verified)
- **One date / condition / session.** All 6 are 20251105 / N2 / no-HS (only sex varies). **No HS, no
  other genotype, no second scope/date.** Re-CV before trusting the absolute numbers there.
- **It's a modest, not dominant, win.** Triangle beats otsu in **5 of 6** leave-one-gonad-out folds
  (dropping HERM_002, otsu/es3 edges ahead). Call it "best available," not "robustly dominant."
- **`effect_size_min = 3` is a knife-edge, not a plateau** — CCC drops sharply at es=2 (over-counts)
  and es≥4 (under-counts); es=3 wins partly by zeroing the mean bias on *these* gonads. Expect to
  re-tune per condition.
- **Per-nucleus CCC ≈ 0.48 is modest**; the cleaner story is **gonad-mean agreement (CCC 0.80, MAE
  ~1 spot/nucleus)**. Residuals are mixed-sign scatter (4 over / 2 under), *not* a correctable
  constant bias — do **not** apply a global correction.
- Our cellpose masks are ~1.4× the Imaris surface volume; only **recall + per-nucleus** are valid
  comparisons (Imaris masked a subset).
- **Confirmed end-to-end:** the committed triangle/es3 config, run through the full pipeline on
  HERM_001, reproduces the CV exactly (1653 spots, matching `cv_detect` grid) — so it's wired
  correctly. The *generalization* caveats above (single date/condition) still stand.

## 5. How this was verified

Two independent background workflows (4 investigators → 4 adversarial verifiers → synthesis each).
The first **refuted** my initial "no spot-splitting" conclusion (my <0.3 µm test was blind to
z-splitting since the z-voxel is 0.2 µm) and fixed a bug list. The second **confirmed** triangle/es3
as the best config but **corrected** my overclaims (the 5-of-6-fold dependency; "~1.5" not "within
1.5" since MALE_001 = 1.505; otsu's flood is 109.7 not 130) and found two more bugs now fixed: the
GT coloc reader double-counted spots when two Imaris surfaces share a mask label (~3% of nuclei), and
a `fillna` mismatch between the CV script and production. Net: conclusion stands, numbers tightened.

## 6. Open calibration work (priority order)

1. **Add HS gonads** — highest priority; RAD-51 density + image statistics differ most under heat, so
   the threshold/effect-size interaction is most likely to need re-tuning there.
2. **≥1 non-N2 genotype + a second date/scope**; re-validate per-nucleus **recall and CCC** there.
3. **Pin one ground truth per metric.** The `.ims` (1222 spots) and xlsx-raw (3086) disagree ~2.5×
   for the same gonad; the lab's canonical reference is the **xlsx coloc** method.
4. Tooling to do all of the above already exists: `scripts/cv_run.py` (batch),
   `scripts/build_cv_manifest.py` (auto-pair GT), `scripts/cv_analyze.py` / `cv_sweep.py` /
   `cv_detect.py` (per-gonad + pooled agreement, threshold × effect-size sweeps).

## 7. How to run

```powershell
$py = ".venv\Scripts\python.exe"
# full pipeline on one .nd2 (writes tidy tables + montage + label mask):
& $py -m germquant.cli run "<image.nd2>" --config config/config.yaml --out results_dir
# or the bundled real-data runner:
& $py scripts/run_real.py
# validate the spots/segmentation vs the Imaris .ims for the same image:
& $py scripts/validate_same_image.py
# Imaris GT distribution for the N2 controls (xlsx coloc method):
& $py scripts/imaris_gt_counts.py
```

Readers: `germquant.validate.imaris_xlsx` (lab coloc method, canonical per-nucleus counts) and
`germquant.validate.imaris_ims` (Spots + Surfaces in image-relative µm, same frame as our output).

## 8. Cellpose retraining feasibility

There are **no ready-to-use nucleus label masks** on the E: drive — but 188 `.ims` projects hold
Imaris nucleus Surfaces, and 179 `.nd2` raw stacks. The catch: Imaris stores surfaces in a
**proprietary block-tree format** (`/Scene8/Content/MegaSurfaces0`), *not* voxel masks or simple
meshes, so exact extraction is hard/unreliable.

**Recommended path (cleanest):** export Surface label masks directly from Imaris (you have it):
1. In Imaris, open a gonad with traced nucleus Surfaces.
2. Surfaces → **Mask All** / *Convert to Labels* → export the label channel as a 3D TIFF.
3. Save the matching raw DAPI channel as a 3D TIFF (or keep the `.nd2`).
4. Do this for ~20–50 gonads → `(raw.tif, nuclei_labels.tif)` pairs.
5. Retrain with `germquant finetune` (or `scripts/train_nucleus_model_combined.py`), **splitting by
   animal/gonad, not by slice**, and **keep the existing model** (`germline_nuclei_combined`) — the
   new one must beat it on held-out gonads before replacing it.

This is **optional**: the current model already scored recall 1.000 here and F1 0.98 historically.
More real labels would help generalization but are not blocking the SpotMAX work.
