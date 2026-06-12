#!/usr/bin/env bash
# Full-resolution GPU run + ground-truth validation for germquant (RTX 5090 / HPC).
#
# Usage:
#   pixi run -e gpu scripts/run_and_validate.sh <nd2_folder> <results_dir> [truth_dir]
# e.g.
#   pixi run -e gpu scripts/run_and_validate.sh /nas/madeleine_images results/ data/ground_truth
#
# Prereqs (once, on each machine):
#   pixi install                 # resolves + writes pixi.lock (commit it)
#   pixi run -e gpu check-gpu    # must print: CUDA cap (12, 0)  -> 5090 / sm_120 + cu128 torch
#
# Ground-truth files expected under <truth_dir> (see data/README.md):
#   per_nucleus.csv : image_id,nucleus_id,sc_total_length_um,n_fragments,n_foci   (Imaris)
#   zones.csv       : image_id,transition_zone_um,pachytene_um                    (optional)
# If Imaris nucleus numbering can't be matched to germquant's, drop nucleus_id and the
# validate step will simply report n=0 matched — fall back to comparing per-image means.
set -euo pipefail

NDIR="${1:?usage: run_and_validate.sh <nd2_folder> <results_dir> [truth_dir]}"
RES="${2:?results dir}"
TRUTH="${3:-}"
CFG="${GQ_CONFIG:-config/config.yaml}"

echo "== 0. GPU sanity =="
germquant check-gpu 2>/dev/null || python -c "import torch; cc=torch.cuda.get_device_capability(); print('CUDA cap', cc); assert cc==(12,0), cc"

echo "== 1. full-resolution batch (xy-stride 1) =="
germquant batch "$NDIR" --config "$CFG" --out "$RES"
echo "   -> $RES (see batch_summary.csv + per-image __montage.png; eyeball QC flags first)"

if [[ -z "$TRUTH" ]]; then
  echo "No truth dir given — run complete. Add ground truth and re-run to validate."
  exit 0
fi

echo "== 2. concat per-image tidy tables across the results tree =="
python - "$RES" <<'PY'
import sys, os, glob, pandas as pd
res = sys.argv[1]
def collect(tbl):
    fs = [f for f in glob.glob(os.path.join(res, "**", f"*__{tbl}.csv"), recursive=True)
          if os.sep + "_done" + os.sep not in f]
    return (pd.concat([pd.read_csv(f) for f in fs], ignore_index=True), len(fs)) if fs else (None, 0)
for tbl in ("sc_per_nucleus", "nuclei"):
    df, n = collect(tbl)
    if df is not None:
        df.to_csv(os.path.join(res, f"combined_{tbl}.csv"), index=False)
        print(f"   combined_{tbl}.csv  <- {n} images")
# zones: long (one row per zone) -> wide per image to match truth's TZ/pachytene columns
zd, n = collect("zones")
if zd is not None and {"image_id", "zone", "length_um"} <= set(zd.columns):
    w = zd.pivot_table(index="image_id", columns="zone", values="length_um", aggfunc="first").reset_index()
    w = w.rename(columns={"transition_zone": "transition_zone_um", "pachytene": "pachytene_um"})
    w.to_csv(os.path.join(res, "combined_zones_wide.csv"), index=False)
    print(f"   combined_zones_wide.csv  <- {n} images")
PY

V="$RES/validation"
echo "== 3. validate vs ground truth  (CCC + Bland-Altman -> $V/*) =="
val() {  # pred_csv truth_csv pred_col truth_col out_subdir
  germquant validate --pred "$1" --truth "$2" --key image_id nucleus_id \
    --pred-col "$3" --truth-col "$4" --out "$V/$5" || echo "   [skip] $5 (check columns / overlap)"
}

if [[ -f "$TRUTH/per_nucleus.csv" ]]; then
  val "$RES/combined_sc_per_nucleus.csv" "$TRUTH/per_nucleus.csv" n_fragments        n_fragments        sc_fragments
  val "$RES/combined_sc_per_nucleus.csv" "$TRUTH/per_nucleus.csv" sc_total_length_um sc_total_length_um sc_length
  val "$RES/combined_nuclei.csv"         "$TRUTH/per_nucleus.csv" n_foci             n_foci             rad51_foci
else
  echo "   (no $TRUTH/per_nucleus.csv — skipping per-nucleus validation)"
fi

if [[ -f "$TRUTH/zones.csv" && -f "$RES/combined_zones_wide.csv" ]]; then
  # zones are per-image: validate with image_id as the only key
  germquant validate --pred "$RES/combined_zones_wide.csv" --truth "$TRUTH/zones.csv" \
    --key image_id --pred-col transition_zone_um --truth-col transition_zone_um --out "$V/tz_length" || echo "   [skip] tz_length"
  germquant validate --pred "$RES/combined_zones_wide.csv" --truth "$TRUTH/zones.csv" \
    --key image_id --pred-col pachytene_um --truth-col pachytene_um --out "$V/pachytene_length" || echo "   [skip] pachytene_length"
fi

# Segmentation F1 / IoU (optional) — needs a ground-truth label image (e.g. an Imaris-exported
# or hand-painted mask) aligned to one stack's __nuclei_labels.tif:
#   germquant validate --seg-pred "$RES/.../X__nuclei_labels.tif" \
#       --seg-truth "$TRUTH/X_truth_labels.tif" --iou 0.5 --out "$V/segmentation"

echo "Done. Agreement stats (CCC, Pearson, bias, MAE) + Bland-Altman PNGs in $V/"
echo "Then build figures:  pixi run Rscript analysis_R/analyze.R $RES"
