# SYP-3 in p-granules: pipeline + acquisition spec

_Working doc for the ccw77 SYP-3 (mCherry) x PGL-1 (GFP) p-granule partition study. Written 2026-07-14 during the pipeline rebuild. Everything here is aimed at making the measurement defensible on properly-imaged slides; Hallie's current images are exploratory and used only to build/validate the pipeline._

## 1. The bug we found (and fixed)

The channel map used to process every gonad on disk had **SYP and PGL swapped**. It assigned:

- `central_element` (SYP) -> channel index 1
- `granule` (PGL) -> channel index 2

But by three independent checks the truth is the opposite:

| index | emission | fluorophore | morphology | localization | identity |
|---|---|---|---|---|---|
| 0 | 438 nm | DAPI | chromatin | nuclei | DNA |
| 1 | **511 nm (GFP)** | PGL-1::GFP | **puncta** | **on nuclear rim** | **PGL** |
| 2 | **595 nm (mCherry)** | mCherry::syp-3 | **threads** | **inside nuclei** | **SYP** |
| 3 | 681 nm | far-red | rings | envelope | lamin |

So `index1 = PGL`, `index2 = SYP`. The config channel-map **file** has since been corrected (it now resolves PGL->477/index1 and SYP->545/index2 by name), but **all processed outputs on disk predate that fix** and are stale. The pipeline *code* is correct (`syp = central_element`, `pgl = granule`); only the map was wrong, so a re-run fixes every output with no code change.

Guard added: `chload.py` is a single canonical loader that asserts the .nd2 channel order (by emission wavelength) on every read and throws if it ever changes, so this cannot silently recur.

Note: the 0622 and 0708 sets share the identical channel order. The "syp1" vs "syp3" in filenames is a typo; the genotype (mCherry::syp-3) is the source of truth, so the 545 channel is SYP-3 in both.

## 2. Why this is a hard measurement (the honest part)

The biological question is whether SYP-3 (a synaptonemal-complex protein, normally intranuclear) concentrates in cytoplasmic/perinuclear p-granules, and whether heat shock changes that. The partition coefficient (PC = intensity inside condensate / intensity in surrounding cytoplasm) is the standard metric and p-granules are a textbook condensate, so the framing is sound. The threats to defensibility are the data, not the method:

1. **SYP is intranuclear and bright; granules are perinuclear.** Even on confocal there is residual PSF spread (~250 nm) from the bright nuclear SC into immediately adjacent voxels, exactly where docked granules sit. Distance-matching reduces this bias; a thin (~1 PSF) exclusion rim at the envelope removes most of the rest.
2. **Exposure confound (current data only):** SYP was shot at 90 vs 200 ms covarying with noHS/HS. The PC is a within-gonad ratio so it survives better than raw intensity, but cross-condition intensity claims are not trustworthy on this data. Fixed by acquisition discipline.
3. **No calibration:** without a known-in-granule protein and a known-cytoplasmic protein imaged the same way, we do not know what PC value counts as "enriched" in this setup.

Bottom line: a **null** result (PC ~ 1, SYP not in granules) is defensible and matches the Imaris eyeball. A **positive** result is only defensible with the acquisition controls in section 3.

## 3. Acquisition spec for the new slides (this is what makes it publishable)

1. **Confocal, fixed exposure / laser / gain across ALL conditions and both channels.** Kills the 90-vs-200 ms confound. Same settings noHS and HS, male and herm.
2. **One SYP marker** (SYP-3), consistent across the whole set.
3. **Calibration controls, imaged identically (the single highest-value add):**
   - a **positive**: a protein known to be in p-granules (GLH-1 or PGL-3) -> defines "enriched" PC.
   - a **negative**: a diffuse cytoplasmic protein -> defines PC ~ 1 baseline.
4. **Nyquist z-sampling** so granules are resolved in z (current z-step 0.2 um is okay; keep or improve).
5. **Autofluorescence / bleed-through check:** a no-fluorophore or single-label control to confirm the 545 channel has no PGL bleed and vice versa. The RAD-51-in-red slide (pgl-1::gfp, no SYP) is the strongest version of this.
6. **TetraSpeck bead stack on the same scope, same channels/settings (HIGH VALUE, two-for-one).** It gives (a) the **chromatic-aberration registration** (x/y/z shift + warp between GFP and mCherry, which we currently cannot measure, our data-driven estimate is confounded by biology; lateral looks small centrally, axial is unknown and is the axis we care about), AND (b) an **empirical PSF** for deconvolution (strictly better than the Gaussian PSF currently used). Image beads across the field so the edge warp is captured.
7. Note laser lines / dyes in the filename or metadata consistently so the channel map never has to be guessed.

Imaging is spinning-disk confocal, Nikon Plan Apo lambda D 60x/1.42 oil (n=1.515), 50 um pinholes, voxel 0.108x0.108x0.2 um. The 50 um pinhole (~0.83 um in sample space) is wider than the Airy disk, which is why axial out-of-focus bleed is significant and deconvolution is warranted.

## 4. Analysis pipeline (state as of tonight)

- **Load:** `chload.py`, channel order asserted. PGL=idx1, SYP=idx2.
- **Nuclei:** Cellpose 3D on DAPI (already good; reused from disk, not re-run). Nuclei are excluded from everything downstream so intranuclear SYP cannot contaminate.
- **Germline:** existing germline mask; cytoplasm region = germline tissue minus nuclei.
- **Granule detection:** moving to **SpotMAX** (the pipeline's validated spot detector, already used for RAD-51), pointed at the PGL channel inside the germline-cytoplasm mask. The hand-rolled LoG/threshold detector worked on clean gonads but missed dim puncta and could not be tuned to catch all of them; SpotMAX is the right tool. [being validated tonight]
- **Measurement (partition coefficient):** for each gonad, PC = (SYP inside granules - bg) / (SYP in cytoplasm - bg), background-subtracted.
  - **Distance-matched** variant: compare SYP in granule vs non-granule cytoplasm at the same distance from the nearest nucleus (removes the perinuclear-blur bias). This is a custom control, not a field standard, but defensible.
  - **Rotation null:** rotate the granule mask 180 deg and recompute; PC should exceed the null if enrichment is granule-specific. Custom specificity check.
- **Per-gonad QC:** a PGL granule-vs-cytoplasm contrast flag auto-drops diffuse/garbage gonads (e.g. HS herm 0622) so they never enter the stats.

## 5. Open decisions / to confirm with Hallie

- Confirm SYP-3 (not SYP-1) is the 545 knock-in in both sets.
- Which calibration control proteins are available in the strain collection (GLH-1 / PGL-3 / a cytoplasmic marker).
- Whether to re-run the existing pipeline on all gonads with the corrected map (fixes the stale object/coloc outputs; CPU-only, no Cellpose needed) or just carry the corrected scratch analysis forward.

## 6. RESULT of the corrected run (fig51, n=22, 2026-07-15)

All 22 gonads processed, all passed QC. SYP-3 distance-matched partition coefficient by condition:

| condition | n | PC (dist-matched) | rotation null |
|---|---|---|---|
| male noHS | 5 | 1.50 | 0.87 |
| male HS | 9 | 1.69 | 0.98 |
| herm noHS | 4 | 1.45 | 0.93 |
| herm HS | 4 | 1.45 | 0.95 |

- SYP-3 reads **~1.5x enriched** in p-granules in **every** condition, above the rotation null (~0.9).
- **No heat-shock effect:** male HS vs noHS p=0.11, herm p=0.89.
- The enrichment being **uniform and present in noHS** is the tell. If this were HS-induced biology, HS would exceed noHS. It doesn't. The most likely cause is a **systematic artifact: GFP(PGL) -> mCherry(SYP) spectral bleed-through** (a PGL-bright voxel reads as SYP-bright regardless of real SYP). The rotation null sitting just below 1 is consistent with this.
- **Conclusion: no defensible evidence of heat-shock-induced SYP-3 partitioning into p-granules from this data.** This matches the Imaris eyeball (SYP intranuclear, PGL perinuclear, separate compartments) and the long-standing skepticism about "coloc in noHS."
- The **single-label bleed-through control** (section 3.5) is the one experiment that resolves this. If, on fixed-exposure confocal with that control, the enrichment persists above the bleed-through baseline AND tracks heat shock, then it's real. Until then it is not.

## 7. What is NOT yet done

- Test the bleed-through hypothesis directly (needs the single-label control image; not in current data).
- Port the SpotMAX granule detector + PC into the pipeline proper (`granule/segment.py`, `coloc/metrics.py`).
- Optional: re-run the full germquant pipeline on all gonads with the corrected channel map to refresh the stale object/coloc CSVs (CPU-only, no Cellpose).
- Add the thin-rim (~1 PSF) envelope exclusion as a sensitivity check.
