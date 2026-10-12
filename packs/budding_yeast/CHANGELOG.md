# Changelog of the budding yeast pack

## 0.1.0 (2026-10-11)

- First version, provisional (not validated): the recipe cells_spots (dims projected) with the starting values
  of docs/DESIGN_system_packs.md, section 3.12: Cellpose-SAM (cpsam_v2@0f1cc3f7) on the sharpest plane of the
  transmitted channel, whole-cell masks copied through z, spots counted per cell in 3D by SpotMAX with the
  authors' SUT509 settings (sharpen on, so effect_size_min 2.5 is on the sharpened scale).
- The frozen base config/profiles/pack_base.yaml, written once from the generic profile (every worm stage
  switch off, the neutral file-name regex, qc.germline_outlier_factor 1e9, qc.min_nuclei 10).
- The channel map yeast_transmitted_foci.yaml, whose transmitted role has no index fallback.
- The model card of cpsam_v2; the synthetic example and its recipe golden; the real example
  sbiad3750_sut509_p1 (BioImage Archive S-BIAD3750, CC0), fetched with `germquant pack fetch-example`.
