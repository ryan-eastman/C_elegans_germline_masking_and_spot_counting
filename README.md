# C. elegans germline — nucleus segmentation + RAD-51 spot counting

Reproducible, **headless** 3D quantification of RAD-51 (or other) foci in the *C. elegans* meiotic
germline from confocal `.nd2` z-stacks, **cross-validated against Imaris**.

**Pipeline:** `read .nd2 → segment nuclei (fine-tuned Cellpose) → isolate the gonad → count spots
inside each nucleus (SpotMAX) → tidy per-nucleus + per-spot tables (CSV/Parquet) + QC montage`.
Counting is headless and reproducible; the GUI / Imaris is for *viewing*, not counting.

## Status
- **Segmentation — production-ready.** Fine-tuned Cellpose model (`germline_nuclei_combined`);
  instance **recall 0.96–1.0** vs Imaris-traced nuclei across 6 gonads.
- **RAD-51 spot counting — cross-validated on N2 no-HS.** `threshold_triangle` + effect-size
  filtering + a z-axis spot-split fix track the lab's Imaris coloc counts (**gonad-mean CCC 0.80,
  ~1 spot/nucleus MAE**). Full story + honest caveats: [docs/SPOTMAX_VALIDATION.md](docs/SPOTMAX_VALIDATION.md).
  Heat-shock / other-genotype CV is the next step (tooling is in `scripts/cv_*.py`).
- **Voxel-size aware throughout** — voxel size is read from each `.nd2` and threaded into every 3D op.

## Install
```bash
pip install -e .                        # core (CPU): read, measure, render, batch
pip install -e ".[gpu,viz,workflow]"    # RTX 5090 / HPC: Cellpose-SAM (cu128 torch), napari, snakemake
pip install spotmax cellacdc            # the SpotMAX spot detector (+ its framework)
```
The spot pipeline also needs the trained model at `models/models/germline_nuclei_combined`
(gitignored, ~1.2 GB; rebuild via `scripts/train_nucleus_model_combined.py`). For a pinned
environment use [pixi.toml](pixi.toml) / [environment.yml](environment.yml); the GPU container is
[Dockerfile](Dockerfile) → [apptainer.def](apptainer.def).

## Quickstart
```bash
# 1. inspect a file's channels + voxel size (no pixels loaded)
germquant info "STACK.nd2"

# 2. run one stack: segment -> germline -> SpotMAX spots -> tidy tables + montage
germquant run "STACK.nd2" --config config/config.yaml --out results/

# 3. batch a whole folder (mirrors the tree, writes batch_summary.csv)
germquant batch /path/to/nd2s --config config/config.yaml --out results/
```

## Outputs (per image — tidy, long-format, one row per object)
| file | one row per | key columns |
|---|---|---|
| `*__nuclei.csv` | segmented nucleus | `nucleus_id`, `volume_um3`, centroid (µm), `in_germline`, **`n_spots`**, axis position |
| `*__spots.csv` | detected spot | `nucleus_id`, position (µm), **`effect_size`**, intensity |
| `*__image_summary.csv` | image | nuclei, germline nuclei, **`mean_spots`**, QC flags |
| `*__nuclei_labels.tif` | — | 3D integer label mask (QC / Imaris import) |
| `*__montage.png` | — | QC montage (nucleus outlines + spots) |

Every row carries a shared metadata block (`image_id`, genotype, sex, treatment, voxel size,
`git_sha`, `config_hash`) so analysis in R/Positron joins on `image_id` (+ `nucleus_id`).

## Configure
The **`spots:`** block in [config/config.yaml](config/config.yaml) holds the cross-validated detection
parameters — `thresholding_method: threshold_triangle`, `effect_size_min: 3.0`,
`merge_z_columns: true` (the z-axis spot-split fix). Channel→role mapping lives in
`config/channel_maps/` (the N2 set: `405`→DAPI, `477`→SYP, `545`→RAD-51). The legacy
**SC / zones / blob_log-foci** stages are present but **disabled** (`enabled: false`); this fork's
pipeline is segment → spots → Imaris.

## Validate against Imaris
```bash
python scripts/imaris_gt_counts.py        # per-gonad RAD-51/nucleus from the xlsx coloc method
python scripts/validate_same_image.py     # our pipeline output vs the .ims for one image
# cross-gonad CV (segment a folder of gonads, then sweep threshold x effect-size vs Imaris GT):
python scripts/cv_run.py            --ndir <nd2_dir> --out results_cv
python scripts/build_cv_manifest.py --ndir <nd2_dir> --xdir <xlsx_dir> --results results_cv \
                                    --xlsx-glob "<glob>" --out cv_manifest.json
python scripts/cv_detect.py         --manifest cv_manifest.json --ndir <nd2_dir> --out results_cv_detect
```
Readers: `germquant.validate.imaris_xlsx` (lab coloc method) and `germquant.validate.imaris_ims`
(Spots + Surfaces in image-relative µm). Imaris masks only a *subset* of nuclei, so only **recall**
and **per-nucleus** counts are valid comparisons — not nucleus count / precision / F1.

## Layout
- `src/germquant/` — the package: `io/` (nd2 reader), `segment/` (Cellpose), `germline/` (gonad
  isolation), `spots/` (SpotMAX), `measure/`, `validate/` (Imaris readers + agreement stats),
  `pipeline.py`, `cli.py`.
- `scripts/` — `run_real.py`, `cv_*.py` (cross-validation), `validate_*` / `*imaris*` (validation),
  `train_nucleus_model_combined.py` (retrain the segmentation model).
- `docs/` — **SPOTMAX_VALIDATION.md** (calibration + caveats), ARCHITECTURE.md, RUNBOOK.md, ANNOTATION.md.
- `config/` — `config.yaml` + `channel_maps/`. · `tests/` — `pytest` (CPU, no GPU/data needed).
