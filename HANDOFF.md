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
  stages stripped; all old-pipeline references audited out; **21 tests pass** (`$py -m pytest`).
- **QC:** `germquant batch` flags germline-count outliers (two-arm / framing) in `batch_summary.csv`.

## 3. What's RUNNING / IN PROGRESS

- **DLW188 (syp-2 het) mutant validation** — background, ~6 h. Outputs to `results_cv_dlw188_nohs/`
  (9 gonads) then `results_cv_dlw188_hs/` (8 gonads). Check progress:
  `Get-Content results_cv_dlw188_nohs\run.log | Select-String "=== done|COMPLETE"`.
  If it finished, analyze it (§5). If it died on restart, re-run it (§4).

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
