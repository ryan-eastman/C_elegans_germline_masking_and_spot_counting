# Architecture — *C. elegans* germline RAD-51 spot quantification

> Design of record for **this** pipeline: `segment nuclei → isolate germline → count spots (SpotMAX)
> → validate vs Imaris`. (The original fork also had SC-tracing / zone-calling / blob_log-foci
> stages; those were removed — recover them from the upstream `C_elegans_ml` repo if ever needed.)

**Scope:** confocal Nikon `.nd2` Z-stacks → headless, reproducible per-image (or batch) run → tidy
CSV/Parquet (one row per nucleus / per spot) + QC montage. Headline readout: **RAD-51 foci per
nucleus**, cross-validated against the lab's Imaris counts.

**This dataset (confirmed from the files):** 3-channel uint16 stacks, Z ≈ 67–91, 2600×2600, voxel
**0.108 × 0.108 × 0.20 µm** (z anisotropic ≈ 1.85× xy). Channels named by laser line → `405` =
**DAPI** (DNA), `477` = **SYP** (synaptonemal complex; used to find the gonad), `545` = **RAD-51**
(the foci we count). Filenames encode genotype (N2), sex (HERM/MALE), treatment (HS/noHS), replicate.

---

## 1. Why open-source Python, with Imaris only as ground truth

Build the load-bearing, unattended pipeline in open Python on the RTX 5090; keep **Imaris off the
critical path**, used only to **validate** the automated counts. Imaris is Windows/macOS-only (no
Linux), its "Batch Process" is GUI point-and-click, and its scripting bridge needs the GUI + legacy
Python — none of which fits "point at a folder → run unattended". The lab's Imaris exports (per-spot
`Spots` + per-nucleus `Surfaces`, counted by the coloc-channel-max method) are the **ground truth** we
calibrate and cross-validate the open detector against (see [SPOTMAX_VALIDATION.md](SPOTMAX_VALIDATION.md)).

## 2. End-to-end stages

Voxel spacing `(dz,dy,dx)` is read once at stage 0 and is a **required argument** to Cellpose
`anisotropy`, the SpotMAX spot radii, and `regionprops` voxel-volume. A single forgotten `spacing`
silently corrupts every length/volume — so it is threaded through every 3D op, never hardcoded.

| # | Stage | Tool / method |
|---|-------|---------------|
| 0 | Read `.nd2` + metadata | **`nd2`** (Talley Lambert) — voxel size, channel names. Pull `(dz,dy,dx)` once. |
| 0b | Channel-role resolution | config name→role map (`config/channel_maps/`), validated at load; missing expected channel = QC flag. |
| 1 | Nucleus segmentation (3D) | **Cellpose** `do_3D, anisotropy=dz/dxy` with a **fine-tuned germline model** (`germline_nuclei_combined`); recall 0.96–1.0 vs Imaris-traced nuclei. Stock cpsam fragments crowded pachytene nuclei → fine-tuning is essential. |
| 2 | Measure nuclei | `regionprops` (volume, centroid, per-channel intensity), spacing-aware. |
| 3 | Germline isolation | drop nuclei segmented OUTSIDE the gonad (gut, debris) via **SYP intensity + spatial connectivity** (`syp_seeded_cc`) — never the spot count, so spot recall stays an independent check. Adds `in_germline`; excluded nuclei kept (flagged) for a complete record. |
| 4 | Axis linearization | principal-curve centerline through germline nucleus centroids → per-nucleus distal→proximal position (`axis_position_norm/um`). |
| 5 | Spot counting | **SpotMAX** (`germquant.spots`): preprocess → semantic segmentation (`threshold_triangle`) → peak detection → 160 per-spot features incl. spot-vs-background **effect size** → effect-size filter + **z-axis spot-split merge**. Cross-validated config (§ SPOTMAX_VALIDATION). |
| 6 | Tidy output | long CSV + Parquet, one row per nucleus / per spot / per image; every row carries voxel size + git SHA + config hash. |
| 7 | QC / render | overlay montage (nucleus outlines + spots) + QC flags (0 objects, missing channel). |

## 3. Method grounding

- **Segmentation:** crowded pachytene nuclei (6 chromosomes balled together) defeat zero-shot
  Cellpose, which shatters one nucleus into fragments. Fine-tuning on real hand-labels fixes this
  ([microPub 2023](https://micropublication.org/journals/biology/micropub.biology.001062/)); our
  model holds whole nuclei across `do_3D` orthogonal reslices (the regime stock collapses on).
- **Germline isolation:** SYP-seeded connected component keeps the whole gonad tube (including the
  SYP-negative distal tip, because it's spatially contiguous) while dropping disconnected gut/debris.
- **Spot counting (SpotMAX):** a generalist 3D spot detector that takes **external** nucleus masks
  and over-proposes candidates, then filters by per-spot effect size (spot vs local background) —
  far more sensitive than the old `blob_log`. The detection threshold and effect-size cut are
  **cross-validated against Imaris coloc counts across 6 gonads** (the otsu default is unstable
  per-image; `threshold_triangle` + `effect_size_min=3` is robust). Details + caveats:
  [SPOTMAX_VALIDATION.md](SPOTMAX_VALIDATION.md).
- **Ground truth (Imaris):** the lab fills each nucleus Surface with a random mask value, takes each
  spot's max in that channel, and groups spots by value → spots-per-nucleus. Readers in
  `germquant.validate` reproduce this from the `.xlsx` (coloc) and read the `.ims` (Spots/Surfaces)
  directly. Imaris masks only a *subset* of nuclei, so only **recall** + **per-nucleus** counts are
  valid comparisons (not nucleus count / precision / F1).

## 4. Reproducibility & GPU

- **GPU (verified):** RTX 5090 = **sm_120**; needs **PyTorch ≥2.7 (cu128)** wheels (first stable with
  sm_120). `germquant check-gpu` asserts `torch.cuda.get_device_capability()==(12,0)`.
- **Env:** `uv` venv on Windows (this box) / `pixi` for Linux+HPC; the trained model lives at
  `models/models/germline_nuclei_combined` (gitignored, ~1.2 GB).
- **Provenance:** every run writes `run_manifest.json` (git SHA, tool versions, resolved config,
  per-file voxel size); every output row carries `git_sha` + `config_hash`.
- **Tests:** `pytest` (CPU) covers the spacing-critical maths, segmentation/germline/spots units, and
  a synthetic full-pipeline run — no GPU or data needed.
