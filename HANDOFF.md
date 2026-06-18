# HANDOFF — current state + how to pick back up

**Snapshot date:** 2026-06-17 · **Repo:** this fork (origin = github.com/ryan-eastman/C_elegans_germline_masking_and_spot_counting)
**Everything below is committed + pushed** (latest commit `5e9bc0f`). A restart cannot lose committed code.

---

## 1. If the computer restarts / moves to the lab

1. **Reconnect the E: drive** (the 2TB external). Every validation run reads raw images + Imaris
   ground truth from `E:\Madeleine\...`. Nothing validation-related works without it.
2. **Resume the Claude session:** open a terminal **in `C:\Users\ryane`** and run `claude --resume`
   (or `claude --continue`), pick this session, and say "pick up the handoff".
3. **Background jobs do NOT survive a restart.** Any run that was in flight is dead; partial per-gonad
   results already on disk are kept, but re-launch the run (commands in §4). No data is lost.

All commands below use the project's Python:

    cd "C:\Users\ryane\C_elegans_germline_masking_and_spot_counting"
    $py = ".venv\Scripts\python.exe"

---

## 2. What's DONE (validated + committed)

- **Pipeline:** `read .nd2 -> segment nuclei (trained Cellpose) -> isolate germline -> linearize axis
  -> count RAD-51 spots (SpotMAX) -> tidy CSV + montage`. Config: `config/config.yaml` `spots:` block.
- **Cross-validated config:** `thresholding_method: threshold_triangle`, `effect_size_min: 3.0`,
  `merge_z_columns: true`.
- **Validated vs Imaris (N2):** no-HS gonad-mean CCC **0.80**; HS gonad-mean CCC **0.96** (generalizes,
  no re-tuning). Full story: `docs/SPOTMAX_VALIDATION.md`.
- **Repo is lab-usable:** beginner README + `quantify.bat` (drag a .nd2 on it). Legacy SC/zones/foci
  stages stripped; all old-pipeline references audited out; **25 tests pass** (`$py -m pytest`).
- **QC:** `germquant batch` flags germline-count outliers (two-arm / framing) in `batch_summary.csv`.
- **Flood guard v1 (partial, committed):** the spots stage caps candidate peaks at
  `spots.max_spot_candidates` (30000) before the feature step. Helps, BUT the 2026-06-18 re-test showed
  the wedge is actually UPSTREAM — inside SpotMAX's own detection/segmentation, which is single-threaded
  (~1 of 16 cores) — so the cap alone is NOT sufficient. The real fix is in §3.

## 3. Mutant validation — BLOCKED on one robustness fix (DO THIS FIRST)

- **DLW188 (syp-2):** `HERM _1` runs fine (849 nuclei, 4016 spots, ~21 min). `HERM _2` is a **bad image**
  (almost certainly the bleed-through this dataset is known for): it floods SpotMAX's detector, which
  then grinds **single-threaded for hours and hangs the whole batch** — there's no per-image escape
  hatch. Confirmed twice (overnight 8 h wedge; a 2026-06-18 re-test sat 14 min in detection on 1 core,
  then killed). The cap (§2) doesn't help because the stall is upstream of it.
- **So the batch is NOT safe to re-run as-is** — one bad gonad stalls everything. The validated **N2
  result is unaffected**; this only blocks the *mutant breadth* runs.

### PRIORITY FIX (before re-running any mutant batch) — two layers
1. **Cheap up-front flood check** in/around `germquant.spots.detect_spots`: one NumPy pass over the
   RAD-51 channel inside the nucleus masks — if the bright-voxel fraction is abnormally high (calibrate
   the threshold against N2 + DLW188 `HERM_1`, which are normal), set flag `spots:flooded_skipped` and
   **skip** the spots stage (return empty) instead of entering SpotMAX's slow chain. O(voxels), seconds.
2. **Hard per-gonad timeout** as a backstop — run `detect_spots` in a child process and terminate after
   N minutes (flag `spots:timeout`) so NOTHING (detection or features) can ever stall the batch.
   (Windows = spawn, so pass the arrays in; for one gonad that's fine.)
Then re-run §4: bad gonads get flagged + skipped, good ones processed. `_cap_candidates` stays as a
third line of defence. Also worth: lower `max_spot_candidates` (~8–10k) so even kept floods stay fast
single-threaded.
- The rest of the mutant validation (DLW190, CCW68) is in §6, unchanged.

## 4. Re-launch the validation runs (if interrupted)

    # DLW188 syp-2 (the one that was running):
    & $py scripts\cv_run.py --ndir "E:\Madeleine\SYP2_het_mutants\20251021_DLW188_noHS" --out results_cv_dlw188_nohs
    & $py scripts\cv_run.py --ndir "E:\Madeleine\SYP2_het_mutants\20251022_DLW188_HS"   --out results_cv_dlw188_hs

## 5. Analyze a finished validation run (does the config match Imaris for that genotype?)

    # 1. auto-pair each gonad's .ims to its xlsx ground truth (by surface count):
    & $py scripts\build_cv_manifest.py --ndir "E:\Madeleine\SYP2_het_mutants\20251021_DLW188_noHS" `
        --xdir "E:\Madeleine\202511_imaris_export_het_mutants_Crest" `
        --results results_cv_dlw188_nohs --xlsx-glob "20251021_dlw188_nohs_*.xlsx" --out cv_manifest_dlw188_nohs.json
    # 2. report per-gonad + pooled agreement (CCC, bias) vs Imaris coloc:
    & $py scripts\cv_analyze.py --manifest cv_manifest_dlw188_nohs.json --tol 2.5 --out cv_results_dlw188_nohs
    # (repeat with the _hs paths + "20251022_dlw188_hs_*.xlsx" for the heat-shock batch)

## 6. Still TODO (the rest of "validate all genotypes")

- **DLW188:** analyze when its run finishes (§5). [genotype #1 of the mutants]
- **DLW190 (syp-3 het):** raw images are NESTED in `imsfiles/` subfolders with different naming, e.g.
  `E:\Madeleine\SYP3_het_mutants\20251022_DLW190_noHS\20251022_dlw190_nohs_imsfiles\`. `cv_run.py` globs
  one folder non-recursively, so point `--ndir` at the actual subfolder that holds the per-gonad
  `.nd2` (confirm the .nd2 are there + named with HERM/MALE; the dated folder itself only had the
  subfolder). xlsx GT: `20251022_dlw190_*.xlsx` (12 files) in the same export dir.
- **CCW68 (syp-3 + msh-5 het):** has raw `.nd2`/`.ims` but its GT is only the *processed* R-output
  xlsx (`...processed.xlsx` in `msh5_data/`), which has different sheets than the raw Imaris export the
  readers parse. Needs a small reader adapter (or the raw exports) before it can be validated.
- **CCW63 (syp-5 het):** xlsx GT exists but NO raw images on the drive — can't validate.

## 7. Data locations on E: (must be connected)

- Raw images + `.ims`: `E:\Madeleine\N2\...` (N2), `E:\Madeleine\SYP2_het_mutants\...` (dlw188),
  `E:\Madeleine\SYP3_het_mutants\...` (dlw190), `E:\Madeleine\SYP3_het_MSH5_mutants\...` (ccw68).
- Imaris xlsx ground truth: `E:\Madeleine\202511_imaris_export_het_mutants_Crest\` (raw exports for N2,
  dlw188, dlw190); processed exports under `E:\Madeleine\20251227_..._processed_..._SYP2_555_RAD51\`.

## 8. Run the pipeline on a new image (no validation needed)

- Easiest (lab members): drag a `.nd2` onto `quantify.bat`. See `README.md`.
- CLI: `& $py -m germquant.cli run "IMAGE.nd2" --config config\config.yaml --out results\`
- Whole folder: `& $py -m germquant.cli batch "FOLDER" --config config\config.yaml --out results\`
