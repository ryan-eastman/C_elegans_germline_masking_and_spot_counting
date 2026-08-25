#!/bin/bash
# Lane driver for the lamin-mask batch: runs pc_lamin_worker.py on each image_id passed as an
# argument, skipping any gonad whose JSON row already exists (resumable after VPN drops).
PY="C:/Users/ryane/C_elegans_germline_masking_and_spot_counting/.venv/Scripts/python.exe"
CA="C:/Users/ryane/coloc_analysis"
cd "$CA"
mkdir -p "$CA/pc_lamin_rows"
for iid in "$@"; do
  if [ -f "$CA/pc_lamin_rows/$iid.json" ] && ! grep -q '"lamin": null' "$CA/pc_lamin_rows/$iid.json"; then
    echo "SKIP (done): $iid"
    continue
  fi
  echo "START $iid  $(date +%H:%M:%S)"
  PYTHONIOENCODING=utf-8 "$PY" scripts/pc_lamin_worker.py "$iid" >> "$CA/pc_lamin_rows/lane.log" 2>&1 \
    && echo "OK  $iid  $(date +%H:%M:%S)" || echo "FAIL $iid  $(date +%H:%M:%S)"
done
echo "LANE DONE $(date +%H:%M:%S)"
