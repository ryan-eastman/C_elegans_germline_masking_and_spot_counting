# PGL-1 (p-granule) surfacing + SYP↔PGL-1 colocalization

Measures **direct voxel/object overlap** between the SYP signal ("the SC channel") and PGL-1
p-granules, on a 4-channel image (DAPI / SYP / RAD-51 / PGL-1). Built to test whether **cytoplasmic
SYP aggregates coincide with perinuclear P-granules**. RAD-51 spot counting is unchanged and runs
alongside as before.

## Why the region matters (read this first)
SYP's **SC ribbon is inside the nucleus**; **P-granules dock just OUTSIDE the nuclear envelope**
(perinuclear cytoplasm). So overlap measured inside the nucleus mask is ~0 *by construction*, and a
whole-image Pearson is dominated by empty background. This pipeline therefore:

- restricts every metric to a **region R** = the germline-nucleus union **dilated by a perinuclear
  shell** (`coloc.region_dilation_um`, default 1.5 µm), which includes the cytoplasm where granules and
  cytoplasmic SYP aggregates live;
- compares PGL-1 against **two SYP operands**:
  - **`syp_aggregate`** — cytoplasmic SYP blobs segmented in the shell (the **headline**: this is where
    aggregate↔granule coincidence would show up);
  - **`sc_ribbon`** — the intranuclear SC ribbon (a **control**: it should *not* coincide with
    perinuclear granules; if it does, suspect segmentation bleed).

## What it computes (per image, per operand — the `coloc` table)
Headline (object/mask overlap, the defensible "how much do they coincide" numbers):
- `dice`, `jaccard` of the SYP mask ∩ PGL-1 mask
- `frac_granules_overlapping_sc` — fraction of PGL-1 granules touching the SYP mask
  (`coloc.object_overlap_min_frac` sets how much of a granule must be inside to count)
- `overlap_volume_um3`, `mean/median_granule_to_sc_um` (nearest-SYP distance, informative even at 0
  overlap)
- `overlap_pvalue` / `overlap_zscore` — a **translation null**: is the observed overlap more than you'd
  get by randomly repositioning the granule mask in R (`coloc.n_random` shuffles)? This is what makes
  "genuine overlap" defensible.

Intensity coloc:
- `manders_m1` (fraction of SYP intensity inside the PGL-1 mask), `manders_m2` (symmetric) —
  mask-restricted.
- `pearson_r` — **diagnostic only** (inflated by the shared sparse background; never the headline).
- `costes_threshold_*` — only if `coloc.costes: true`.

Per object (`granules` table): one row per PGL-1 granule — centroid, volume, intensity, nearest germline
nucleus, and per-operand `overlaps_*` / `overlap_frac_*` / `nearest_*_um`. Per nucleus, `nuclei` gains
`n_granules` and `granule_volume_um3`.

## Outputs
- Tables: `*__granules.csv`, `*__coloc.csv`, plus `n_granules` on `*__nuclei.csv` and headline coloc
  fields on `*__image_summary.csv`.
- Imaris-loadable TIFs (same voxel grid as the nucleus labels): `*__granules_labels.tif` (→ Surfaces),
  `*__sc_ribbon.tif`, `*__syp_aggregate.tif`.
- Montage: a `coloc: SYP-agg ∩ PGL-1` overlay panel (SYP red / PGL-1 green / overlap yellow).

## Config knobs (`config/config.yaml`)
```yaml
granule: {thresholding_method: threshold_triangle, min_volume_um3: 0.03, max_volume_um3: 8.0}
sc:      {enabled: true, ridge_hyst_low_pct: 45, ridge_hyst_high_pct: 80}
coloc:   {enabled: true, sc_operand: both, region: perinuclear_shell,
          region_dilation_um: 1.5, object_overlap_min_frac: 0.0, n_random: 100, costes: false}
```
The channel map is `config/channel_maps/n2_dapi_syp_rad51_pgl1.yaml` — set the real PGL-1 channel name
(laser line) there so it name-matches instead of the index-3 fallback. `--no-coloc` skips the whole
stage. A 3-channel (no PGL-1) image skips it automatically and RAD-51 output is identical.

## Validating against Imaris (no hand-painting of voxels)
1. **Coloc numbers** → Imaris **Colocalization** module on SYP × PGL-1 within a germline ROI
   (draw a Surface as the mask so the denominator matches R). Export Manders M1/M2 + Pearson.
2. **Granule detection** → Imaris **Surfaces** (or Spots) on the PGL-1 channel; correct a few; export the
   Statistics (`.ims` or `.xlsx`).
3. Score it:
   ```
   python scripts/validate_coloc.py --granules <img>__granules.csv --coloc <img>__coloc.csv \
       --imaris-ims <img>.ims --imaris-coloc-json imaris_coloc.json --out validation_coloc
   ```
   (`imaris_coloc.json` = `{"manders_m1":.., "manders_m2":.., "pearson_r":..}`.) Reports granule
   recall/precision + volume CCC and puts our Manders/Pearson next to Imaris'.
4. **Tune first** without GT if useful: `scripts/coloc_param_sweep.py` shows how the numbers move with the
   granule threshold, volume floor, and shell thickness.

## Known limits
- SC ribbon **length** (optional, needs `skan`, `pip install -e .[sc]`) is usable; per-nucleus SC
  **fragment count** is NOT recoverable in 3D at pachytene density — don't build on it.
- Manders/overlap depend on the granule threshold and `region_dilation_um`; calibrate them against
  Imaris before quoting absolute numbers (build-now, tune-later).
- The two operands are *different SYP pools*; report which one you mean.
