# Validation of the C. elegans germline pack

The recipes of this pack are the engine's frozen config files, validated where the engine documents it:

- docs/SPOTMAX_VALIDATION.md: the SpotMAX spot counts against Imaris truth, and the cross-gonad calibration of
  the spot settings (threshold_triangle, effect_size_min 3.0, the z-column merge).
- docs/BENCHMARK.md: the RAD-51 benchmark harness and its protocol (docs/BENCHMARK_PROTOCOL.md).
- The 23 real goldens (tests/golden/real_manifest.json; 11 dry-ice and stride variants, 12 RAD-51 gonads),
  which every tranche-5 package must reproduce exactly on the RTX 5090.

This folder holds no validation record in the pack format (validation/<recipe>.json, written only by the
public benchmark tool), so every recipe reads "provisional". Nothing here computes or reads a score on the
sealed RAD-51 test split.
