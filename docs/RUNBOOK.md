# Runbook — full-resolution GPU run + tuning

The Mac smoke tests proved the plumbing at low resolution with the classical fallback. The
real quantification happens here, on the RTX 5090 (or HPC), at full resolution with Cellpose.

## 0. One-time setup (on the 5090 workstation)
```bash
git clone https://github.com/ryan-eastman/C_elegans_ml.git
cd C_elegans_ml
pixi install                       # or: conda env create -f environment.yml && pip install -e ".[gpu,sc,viz]"
pixi run -e gpu check-gpu          # MUST print "CUDA cap (12, 0)" — confirms 5090/sm_120 + cu128 torch
```
If `check-gpu` fails, the torch wheel didn't match Blackwell — reinstall torch from the
`cu128` index (see ARCHITECTURE.md §4). Each gonad stack is ~3 GB in RAM; ≥32 GB system RAM
recommended.

## 1. Sanity-check one stack end-to-end (full res, GPU)
```bash
germquant info "/nas/.../20251105_n2_nohs_HERM _001.nd2"     # confirm voxel/channels
germquant run "/nas/.../20251105_n2_nohs_HERM _001.nd2" --config config/config.yaml --out results_test/
open results_test/**/**__montage.png                          # eyeball seg + foci + zones
```

## 2. Tune three thresholds against full-res reality (config/config.yaml — no code edits)
Inspect the montage + the `*__image_summary.csv` and adjust:

| Symptom | Knob | Direction |
|---|---|---|
| RAD-51 foci over-counted (control should be sparse) | `foci.threshold_rel` | raise (e.g. 0.1 → 0.2–0.3) |
| Foci missed | `foci.min_sigma_um` / `max_sigma_um` | widen around the real spot size |
| No transition zone found (all "pachytene") | `zones.polarization_threshold` (in code) / check whole gonad is in frame | lower threshold; TZ needs the distal tip imaged |
| Nuclei merged/over-segmented | fine-tune Cellpose (docs/ANNOTATION.md) + `segmentation.nuclei.diameter_um` | set diameter to measured mean |
| SC length implausible / fragments wrong | `sc.intensity_percentile`, `sc.ridge_sigmas_um`, `sc.min_fragment_length_um` | tune ridge to the SYP filament width |

## 3. Batch the whole NAS folder
Edit `workflow/config.yaml` → `input_root` / `output_root`, then:
```bash
# workstation:
pixi run -e gpu snakemake -s workflow/Snakefile --configfile workflow/config.yaml --cores 8
# HPC (SLURM, 1 GPU/job) — fill in config/profiles/slurm first:
snakemake -s workflow/Snakefile --configfile workflow/config.yaml --workflow-profile config/profiles/slurm
```
Output mirrors the input tree; `batch_summary.csv` lists every image + QC flags.

## 4. Stats & figures (R / Positron)
```bash
pixi run Rscript analysis_R/analyze.R /nas/madeleine_results
# -> figures/{sc_fragmentation,sc_length,rad51_foci,zone_length}.png + group_summary.csv
```
Open `analysis_R/germquant_analysis.Rproj` in Positron to iterate interactively.

## 5. Validate against ground truth (makes it publishable)
- Segmentation: `germquant validate --seg-pred ... --seg-truth ... ` → F1 / IoU.
- Counts/lengths: put your hand scores in a CSV keyed by `image_id,nucleus_id`, then
  `germquant validate --pred X__sc_per_nucleus.csv --truth manual.csv \
       --key image_id nucleus_id --pred-col n_fragments --truth-col n_fragments_manual --out validation/`
  → CCC, Bland-Altman plot.
- SC tracing: trace ~15–25 nuclei in Imaris/Fiji-SNT as the cross-check (ARCHITECTURE.md Risk 1).

## Gotchas
- **Anisotropy**: never disable it; the pipeline asserts non-isotropic spacing in tests.
- **TZ/pachytene length** needs the *whole gonad arm* (distal tip → late pachytene) in frame.
- **No-DAPI zoning** (synapsis-state) is wild-type only and must be calibrated vs DAPI gonads.
- **Blackwell + TensorFlow**: StarDist/CARE lag on sm_120 — stick to the PyTorch path (Cellpose).
