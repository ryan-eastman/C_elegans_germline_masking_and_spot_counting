#!/bin/bash
PY="C:/Users/ryane/C_elegans_germline_masking_and_spot_counting/.venv/Scripts/python.exe"
CA="C:/Users/ryane/coloc_analysis"; cd "$CA"
mkdir -p "$CA/staging"
for iid in "$@"; do
  if [ -f "$CA/staging/${iid}_landmarks.json" ] && [ -f "$CA/staging/cache/${iid}.npz" ]; then echo "SKIP $iid"; continue; fi
  echo "START $iid $(date +%H:%M:%S)"
  PYTHONIOENCODING=utf-8 "$PY" scripts/crescent_axis.py "$iid" >> "$CA/staging/batch.log" 2>&1 && echo "OK $iid" || echo "FAIL $iid"
done
echo "DONE $(date +%H:%M:%S)"
