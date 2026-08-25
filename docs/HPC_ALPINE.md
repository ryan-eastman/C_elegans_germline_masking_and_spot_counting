# Running germquant on RMACC Alpine (SLURM, one GPU job per gonad)

**Why:** a single gonad is *not* faster on Alpine — the workstation RTX 5090 (Blackwell) beats
Alpine's A100 per-GPU for this Cellpose/SpotMAX workload. The win is **parallelism + freeing the
desktop**: queue one 1-GPU SLURM job per `.nd2` and they run simultaneously, so a 100-gonad batch
finishes in ~one gonad's wall-time instead of days, and the 5090 is yours again. Analysis
(coloc/EDA/figures) stays on the desktop — only the **segmentation** needs the GPU fan-out.

The `cu128` torch wheel already includes the A100 (sm_80) and L40 (sm_89) kernels, so **the same
container image is portable** across the 5090, A100, and L40 — no separate build.

## One-time setup

### 1. Build the image (on any Linux box with Docker, or on Alpine with Apptainer)
```bash
# option A — build the Docker image, convert to .sif:
docker build -t germquant:gpu .
apptainer build germquant.sif docker-daemon://germquant:gpu
# option B — build the .sif straight from the published base on Alpine (no Docker):
#   apptainer build germquant.sif apptainer.def
```
The image does NOT contain the 1.2 GB trained model (gitignored) or the raw data — those are
bind-mounted at run time (below).

### 2. Stage to Alpine scratch (Globus or rsync)
```bash
DST=/scratch/alpine/$USER
rsync -a germquant.sif                         $DST/
rsync -a models/models/germline_nuclei_combined $DST/models/models/   # the trained Cellpose model
rsync -a /path/to/raw_nd2/                      $DST/data/             # or Globus for big transfers
```

### 3. Put the container wrapper on PATH as `germquant`
```bash
mkdir -p ~/bin && ln -sf "$PWD/scripts/germquant-apptainer" ~/bin/germquant
export PATH=~/bin:$PATH
export GERMQUANT_SIF=$DST/germquant.sif
export GERMQUANT_MODELS=$DST/models          # -> bound to /opt/germquant/models in the container
module load apptainer
germquant check-gpu                           # sanity: should print "OK — A100 ..." on a GPU node
```

### 4. Fill in your allocation
Edit `config/profiles/alpine/config.yaml`: set `slurm_account` to your allocation, and pick the
partition (`al40` = L40, usually least contended → best for "queue and walk away"; `aa100` = A100).

## Run a batch (queue one GPU job per .nd2)
Point `workflow/config.yaml` at the staged data + output dirs, then:
```bash
snakemake -s workflow/Snakefile --configfile workflow/config.yaml \
          --workflow-profile config/profiles/alpine
```
Snakemake submits one `process_one` SLURM job per gonad (up to `jobs: 20` at once), mirrors the
input tree to the output tree, and writes `batch_summary.csv` at the end. **Resume is automatic** —
re-running skips gonads with an existing `_done/<idx>.done`, so a walltime kill just picks up where
it left off.

## Notes / gotchas
- **Config paths are relative to the working dir** (the channel map is loaded relative to the config
  file). Run `snakemake` from the repo root inside `$DST`, or bind the repo in (`GERMQUANT_BIND`).
- **Model path:** `config.yaml` points at `models/models/germline_nuclei_combined` (relative). The
  wrapper binds `$GERMQUANT_MODELS` to `/opt/germquant/models`, so keep the model under
  `$GERMQUANT_MODELS/models/germline_nuclei_combined` and run from a dir where that relative path
  resolves — or set an absolute `cellpose_model` in a copy of the config.
- **Validate first:** run 1–2 gonads on Alpine and confirm the outputs match the desktop
  (`shell_pearson`, nucleus count) before firing the full batch.
- **Skip `ami100`** (AMD MI100) — it needs a ROCm torch build, not this cu128 image.
- **NAS vs scratch:** Alpine can't see the lab NAS; stage the `.nd2` to `/scratch/alpine` first
  (Globus for the big multi-GB transfers).
