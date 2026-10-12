# Changelog of the cultured cells pack

## 0.1.0 (2026-10-11)

- First version, provisional (not validated): the recipes nuclei_foci (the default, dims projected) and
  nuclei_foci_3d (dims 3d) with the starting values of docs/DESIGN_system_packs.md, section 3.12: Cellpose-SAM
  (cpsam_v2@0f1cc3f7) on the DNA maximum projection or in 3D, foci counted per nucleus in 3D by SpotMAX, with a
  second focus channel (53BP1) as the spots instance foci2 (sharpen off: effect_size_min 3.0 is on SpotMAX's
  usual, unsharpened scale).
- The frozen base config/profiles/pack_base.yaml, written once from the generic profile (every worm stage
  switch off, the neutral file-name regex, qc.germline_outlier_factor 1e9, qc.min_nuclei 1).
- The channel map cultured_dna_foci.yaml.
- The model card of cpsam_v2; the synthetic example and its two recipe goldens; the real example
  deepfoci_valid_258 (DeepFoci, Zenodo 4067741, CC BY 4.0), one member of a 10.2 GB zip fetched with HTTP range
  requests by `germquant pack fetch-example`.
