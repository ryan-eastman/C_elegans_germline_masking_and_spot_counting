#!/bin/bash
# After Ryan re-traces the pachytene region: rezone from the saved polylines, recompute the zone PCs for
# every traced, non-excluded ccw77 gonad, rebuild the tables and the trace-dependent figures.
# Usage: bash scripts/rerun_after_retrace.sh          (log: staging/rerun_after_retrace.log)
PY="C:/Users/ryane/C_elegans_germline_masking_and_spot_counting/.venv/Scripts/python.exe"
CA="C:/Users/ryane/coloc_analysis"
cd "$CA"
LOG="$CA/staging/rerun_after_retrace.log"
export PYTHONIOENCODING=utf-8 MPLBACKEND=Agg
echo "START $(date +%H:%M:%S)" > "$LOG"
mkdir -p pc_zone_rows_v3_pre_retrace && cp pc_zone_rows/*.json pc_zone_rows_v3_pre_retrace/ 2>/dev/null
"$PY" scripts/rezone_all.py >> "$LOG" 2>&1 || { echo "rezone FAILED" >> "$LOG"; exit 1; }
IDS=$("$PY" -c "
import json, sys
sys.path.insert(0, r'C:\Users\ryane\coloc_analysis'); import chload
t = json.load(open('staging/pachytene_traces.json', encoding='utf-8'))
print(' '.join(k for k, v in t.items() if 'ccw77' in k and v.get('status') == 'traced' and not chload.is_excluded(k)))")
echo "zone worker on: $IDS" >> "$LOG"
for iid in $IDS; do
  echo "zone $iid $(date +%H:%M:%S)" >> "$LOG"
  "$PY" scripts/pc_zone_worker.py "$iid" 2>&1 | grep -v -i warning | tail -1 >> "$LOG"
done
"$PY" scripts/collect_tables.py >> "$LOG" 2>&1
"$PY" scripts/build_clean_table.py >> "$LOG" 2>&1
"$PY" scripts/build_pachytene_pooled.py >> "$LOG" 2>&1
"$PY" scripts/build_filter_compare.py >> "$LOG" 2>&1
"$PY" scripts/fig_pub_pachytene_pooled.py >> "$LOG" 2>&1
"$PY" scripts/fig_pub_stage.py >> "$LOG" 2>&1
"$PY" scripts/fig_grant_assets.py 2>&1 | grep -v -i warning >> "$LOG"
"$PY" scripts/fig_grant_pipeline.py >> "$LOG" 2>&1
echo "RERUN DONE $(date +%H:%M:%S)" >> "$LOG"
