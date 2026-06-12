# celegans-germline-quant

Reproducible, headless 3D quantification of the *C. elegans* meiotic germline from confocal
`.nd2` z-stacks.

**Point it at a folder (e.g. on a NAS) → it reconstructs 3D objects → writes a mirrored output
folder of tidy results (for R / Positron) + publication-quality renders.**

## What it measures
- **Synaptonemal complex (SC) per nucleus** — traces the SYP/axis filament in 3D and reports
  total length (µm), per-fragment length distribution, and **fragment count** (the heat
  phenotype: heat *fragments* the SC in spermatocytes but not oocytes).
- **Transition-zone & pachytene zone lengths** along the distal→proximal germline axis
  (DAPI-crescent morphology; no-DAPI synapsis-state fallback for other experiments).
- **RAD-51 foci per nucleus** — 3D blob detection, assigned to segmented nuclei.
- **Granules / generic 3D objects** — optional module (count, volume, surface area, intensity).

Every readout is **spacing-aware**: voxel size (`0.108 × 0.108 × 0.20 µm`, z anisotropic ≈1.85×)
is read from each `.nd2` and threaded into every 3D operation. See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

## Status
Core pipeline **runs end-to-end on real `.nd2` data** (read → segment → axis → zones → SC →
foci → measure → tidy CSV/Parquet → QC montage → provenance manifest). The science modules
(SC tracing, zone-calling, foci thresholds, Cellpose segmentation) are **first-pass and need
full-resolution GPU runs + ground-truth validation** before publication — tracked in
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) §8.

## Install
```bash
# core (CPU) — enough to read, segment (classical), measure, render, batch
pip install -e .

# on the RTX 5090 / HPC, add GPU + SC tracing:
pip install -e ".[gpu,sc,viz]"   # Cellpose-SAM (cu128 torch), skan, napari
```
(A pinned `pixi`/conda env + Apptainer recipe for the 5090 and SLURM lands next — see §6 todo.)

## Quickstart
```bash
# 1. peek at a file's channels + voxel size (no pixels loaded)
germquant info "data/raw_examples/.../20251105_n2_nohs_HERM _001.nd2"

# 2. run one stack (full res on GPU; --xy-stride downsamples for a laptop smoke test)
germquant run STACK.nd2 --config config/config.yaml --out results/

# 3. batch a whole NAS folder -> mirrored output tree + batch_summary.csv
germquant batch /nas/madeleine_images --config config/config.yaml --out /nas/madeleine_results
```

## Configure
- [config/config.yaml](config/config.yaml) — all parameters (segmentation, SC, foci, zones, axis, render, output).
- [config/channel_maps/](config/channel_maps/) — map fluor → biological role per experiment.
  The N2 set uses [n2_dapi_syp3_rad51.yaml](config/channel_maps/n2_dapi_syp3_rad51.yaml)
  (`405`→DAPI, `477`→SYP-3, `545`→RAD-51).

## Outputs (per image, tidy long-format)
`*__nuclei`, `*__sc_tracks`, `*__sc_per_nucleus`, `*__foci`, `*__zones`, `*__image_summary`
(CSV + Parquet) · `*__nuclei_labels.tif` · `*__montage.png` · `run_manifest.json`.
Every row carries `image_id, genotype, sex, germ_cell, treatment, replicate, voxel_d{z,y,x}_um,
git_sha, config_hash, run_timestamp` — join on `image_id` (+ `nucleus_id`) in R, facet by
`sex`/`treatment`/`zone`.

## Layout
```
config/        # YAML: global params + per-experiment channel maps
src/germquant/ # io, segment, axis, zones, sc, foci, measure, render, qc, provenance, pipeline, cli
workflow/      # Snakemake batch orchestration (NAS → mirrored outputs)   [next]
analysis_R/    # R / Positron project consuming the tidy outputs           [next]
tests/         # spacing/anisotropy + filename-metadata unit tests
data/          # raw .nd2 (gitignored)
docs/ARCHITECTURE.md   # the verified, cited design of record
```
