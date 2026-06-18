# Runbook — full-resolution GPU run + spot calibration

The real quantification happens on the RTX 5090 (or HPC) at full resolution with the fine-tuned
Cellpose model. (For a no-coding, one-image run, use `quantify.bat` — see the README. This runbook is
for full-res batches and calibration.)

## 0. One-time setup (native Windows on the 5090, with `uv` — verified)
```powershell
uv venv --python 3.11
uv pip install -e .
uv pip install --index-url https://download.pytorch.org/whl/cu128 "torch>=2.7"   # Blackwell sm_120
uv pip install "cellpose>=3.1" spotmax cellacdc                                   # segmenter + spot detector
.venv\Scripts\germquant.exe check-gpu     # MUST print CUDA cap (12, 0) -> 5090/sm_120 + cu128 torch
```
If `check-gpu` fails, the torch wheel didn't match Blackwell — reinstall torch from the `cu128` index.
The fine-tuned model must exist at `models\models\germline_nuclei_combined` (~1.2 GB, gitignored).
Each gonad stack is ~3 GB in RAM; ≥32 GB system RAM recommended.

## 1. Sanity-check one stack end-to-end (full res, GPU)
```powershell
.venv\Scripts\germquant.exe info "E:\...\20251105_n2_nohs_HERM _001.nd2"   # confirm voxel + channels
.venv\Scripts\germquant.exe run  "E:\...\20251105_n2_nohs_HERM _001.nd2" --config config\config.yaml --out results_test
```
Then open `results_test\...__montage.png` and eyeball that nuclei outlines + spots look right, and
`...__image_summary.csv` for `n_germline_nuclei` + `mean_spots`.

## 2. Spot-detection config (already cross-validated; re-tune only for new conditions)
The `spots:` block in `config/config.yaml` is **cross-validated against Imaris on N2 no-HS**:
`thresholding_method: threshold_triangle`, `effect_size_min: 3.0`, `merge_z_columns: true`.

| Symptom | Knob | Direction |
|---|---|---|
| Spots flooded with noise on some images | `effect_size_min` | raise (3 → 4–5) — but re-check against Imaris |
| Real foci missed (counts too low) | `thresholding_method` | `threshold_triangle` (sensitive) vs `threshold_otsu` (stricter) |
| One focus counted as 2–3 (z-stacking) | `merge_z_columns` / `z_merge_gap_um` | keep on; widen the gap if needed |
| Nuclei merged / fragmented | retrain Cellpose (docs/ANNOTATION.md) | — |

**Do not trust new conditions (HS, other genotypes) without re-running the cross-validation below —**
the otsu/effect-size optimum shifts with image statistics. See [SPOTMAX_VALIDATION.md](SPOTMAX_VALIDATION.md).

## 3. Batch a whole folder
```powershell
.venv\Scripts\germquant.exe batch "E:\...\folder_of_gonads" --config config\config.yaml --out results_batch
```
Skips the `10x`/`largeimage` overviews, mirrors the input tree, and writes `batch_summary.csv`
(one row per image + QC flags). Snakemake (`workflow/Snakefile`) does the same on HPC/SLURM.

## 4. Validate / re-calibrate against Imaris
```powershell
# per-gonad Imaris RAD-51/nucleus from the xlsx coloc method:
.venv\Scripts\python.exe scripts\imaris_gt_counts.py
# our output vs the .ims for ONE image (recall + paired per-nucleus):
.venv\Scripts\python.exe scripts\validate_same_image.py
# cross-gonad CV (segment a folder, auto-pair the GT, sweep threshold x effect-size):
.venv\Scripts\python.exe scripts\cv_run.py            --ndir <nd2_dir> --out results_cv
.venv\Scripts\python.exe scripts\build_cv_manifest.py --ndir <nd2_dir> --xdir <xlsx_dir> --results results_cv --xlsx-glob "<glob>" --out cv_manifest.json
.venv\Scripts\python.exe scripts\cv_detect.py         --manifest cv_manifest.json --ndir <nd2_dir> --out results_cv_detect
```
`cv_detect` reports per-gonad + pooled agreement (CCC, bias) for each `(threshold, effect_size_min)`,
so you can confirm — or re-tune — the config for a new condition. To retrain the segmentation model
(more real labels), see [ANNOTATION.md](ANNOTATION.md) + `scripts/train_nucleus_model_combined.py`.

## Gotchas
- **Anisotropy:** never disable it; voxel spacing `(dz,dy,dx)` is read from each `.nd2` and threaded
  into every 3D op. The tests assert non-isotropic spacing.
- **GPU memory:** one big stack at a time; close Imaris/other GPU apps first.
- **Blackwell + TensorFlow:** stick to the PyTorch path (Cellpose); TF tools lag on sm_120.
