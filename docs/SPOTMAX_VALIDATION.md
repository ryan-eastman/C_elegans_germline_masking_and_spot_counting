# SpotMAX spot-counting: integration + Imaris validation

**Date:** 2026-06-17 · **Pipeline:** `segment → isolate germline → SpotMAX spots → Imaris`
**Test image:** `20251105_N2_nohs_HERM_001` (N2 control, 2600×2600×71, voxel 0.108×0.108×0.2 µm)

This documents the switch from the v1 `blob_log` foci detector to **SpotMAX**, and what we learned
validating it against your Imaris ground truth. Findings here were independently re-derived and
**adversarially verified** (a 9-agent review caught and corrected one wrong conclusion — see §4).

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
| **Per-nucleus RAD-51 correlation** (ours vs Imaris) | **Pearson r = 0.856**, Spearman 0.879, p<0.0005 |
| **Total spots** after the z-fix (§3) | **~1269 ≈ Imaris 1222** |

## 3. The over-count was Z-axis spot-splitting (not dim noise) — FIXED

Uncalibrated, SpotMAX detected 2033 spots (~1.6× Imaris). The cause is **not** spurious dim peaks —
it is **z-axis spot-splitting**: the confocal axial PSF (~0.6–0.8 µm) is far wider than the z
spot-radius (1.5 voxels = 0.3 µm), so one RAD-51 focus is detected as **2–3 stacked peaks at the
same (x,y) voxel**. 60.7 % of in-nucleus spots share an exact (x,y) column.

**Fix:** `detect_spots(..., merge_z_columns=True, z_merge_gap_um=0.8)` collapses peaks sharing an
(x,y) voxel within `z_merge_gap_um` (keeps the brightest), while leaving genuinely distinct foci
that are farther apart in z. Result: **2033 → ~1269 spots ≈ Imaris's 1222.** On by default.

**Do NOT use `effect_size_min ≈ 4.0`.** It can hit the same per-nucleus number arithmetically, but
it is overfit to one gonad/stage, collapses the whole-germline mean, and *masks* the z-splitting
defect instead of fixing it. `effect_size_min` stays **0.0**.

## 4. How this was verified

A background workflow ran 4 independent investigators → 4 adversarial verifiers → 1 synthesizer.
It **confirmed** gonad identity, canonical 5.81, recall 1.000, and r=0.856; and **refuted** my
initial "no spot-splitting" conclusion (my <0.3 µm test was structurally blind to z-splitting,
since the z-voxel is 0.2 µm). It also produced the bug list now fixed in this commit (silent
`except`, n_spots off-gonad invariant, radius-aware `spots_per_surface`, hard-coded 5.84→5.81).

## 5. Open calibration work (before trusting absolute RAD-51 numbers)

1. **Run the pipeline on all N2 gonads that have xlsx GT** (HERM_1/2/3, MALE_1/2/3, 20251021 set).
   Today our per-spot data exists only for HERM_001, so cross-validation is impossible.
2. **Leave-one-gonad-out**: confirm `merge_z_columns` (and any threshold) generalizes — the GT
   per-nucleus count itself ranges 3.72–8.15 across gonads, so single-gonad tuning is fragile.
3. **Pin one ground truth per metric.** The `.ims` (1222) and xlsx-raw (3086) disagree ~2.5× for
   the same gonad; the lab's canonical reference is the **xlsx coloc** method (848 assigned → 5.81).

## 6. How to run

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

## 7. Cellpose retraining feasibility

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
