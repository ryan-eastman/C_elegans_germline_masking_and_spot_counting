# Changelog of the test pack

## 0.1.0 (2026-10-11)

- First version: recipes basic (classical, LoG) and with_model (Cellpose with the stand-in tinynet@f5055f8a),
  the frozen base written from config/profiles/generic.yaml, a synthetic example and recipe goldens.
- Before any release (tranche 5, package P3): both recipes set segmentation.nuclei.require_channel: true, which
  pack lint now asks of every outside recipe (new hashes 02bb83cf74f9 and 81e38c5cf371; the golden tables did
  not change).
