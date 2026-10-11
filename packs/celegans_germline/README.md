# C. elegans germline pack (celegans_germline 0.1.0)

The germquant recipes for the C. elegans meiotic germline, offered as a system pack. Each recipe IS a frozen
config file of the engine checkout (`root: checkout`), never a copy, so

    germquant run GONAD.nd2 --pack celegans_germline --out results/GONAD

writes exactly what `germquant run GONAD.nd2 --config config/config.yaml --out results/GONAD` writes, with the
same config_hash (84bdda5a5230) and stage hashes. Only the run manifest gains two records: the pack and the
runtime (device, torch, python).

## What it is for, and not for

For: fixed, dissected C. elegans gonads imaged as confocal z-stacks (Nikon .nd2) with DAPI and SYP, plus
RAD-51 foci, PGL-1 granules or lamin where a recipe uses them; and live GFP::SYP-2 time-lapse series.

Not for: other organisms or samples. Use an outside pack (cultured cells, budding yeast) or the generic profile.
The worm recipes isolate the germline by SYP, fit a gonad axis and parse the worm file-name convention; on other
samples those stages give wrong answers.

## Recipes

| recipe | dims | profile | stages beyond nuclei | channel naming |
|---|---|---|---|---|
| default | 3d | config/config.yaml | germline, axis, RAD-51 spots, P granules and coloc | n2_dapi_syp_rad51_pgl1 (405, 477, 545, 647) |
| rad51_foci | 3d | config/profiles/rad51_foci.yaml | germline, axis, RAD-51 spots | as default |
| segmentation_only | 3d | config/profiles/segmentation_only.yaml | germline (for Imaris surfaces) | as default |
| n2_sc | 3d | config/profiles/n2_sc.yaml | germline, axis, spots, SC tracing (needs skan) | as default |
| ccw77_partition | 3d | config/profiles/ccw77_partition.yaml | envelope, audit, staging, granules, partition, lit fraction | ccw77_dapi_syp_pgl_lamin |
| live_syp2 | 3d+t | config/profiles/live_syp2.yaml | per-timepoint nuclei on GFP::SYP-2, signal QC | live_gfp_syp2 (no DNA stain) |

bench_log_tuned, bench_spotmax_tuned and config_n2dryice.yaml (a fixation test, regression use only) are not
recipes of this pack.

Status: every recipe reads "provisional" in the pack format, because the pack format calls a recipe
validated only with a validation record written by the public benchmark tool. The engine's own validation of
these recipes is documented in the files listed under Validation.

## Object kind

nucleus: one row of the nuclei table is one nucleus.

## Model

`germline_nuclei_combined` (release asset, sha256 e4e32320...8f57, 1,218,646,375 bytes; card:
`models/germline_nuclei_combined.card.json`): a June 2026 Cellpose-SAM fine-tune on 48 real hand-labelled
germline nucleus slices plus about 49 synthetic slices. Fetch it with `germquant get-model` or
`germquant pack fetch-models celegans_germline`. Its base weights are not recorded; the file dates point to
cpsam_v2, and the card records the tensor-by-tensor comparison with both cached bases (method, and the result
once it has been run on the workstation; until then the base reads "not determined").

Licence: Cellpose-SAM weights are BSD-3-Clause, but Cellpose's own training data include CC BY-NC sets, so a
model fine-tuned from them inherits that question for commercial use.

## How the presets were chosen

The spot settings (threshold_triangle, effect_size_min 3.0, the z-column merge) were cross-validated against
Imaris counts on six N2 gonads (docs/SPOTMAX_VALIDATION.md). The germline isolation keeps the gonad's spots
(mean spot recall about 0.97 on 14 real gonads, scripts/validate_germline_select.py), and the nucleus model
scored F1 0.98 on held-out real slices against 0.87 for stock cpsam (both from the comments of
config/config.yaml).

## Validation

See `validation/README.md`: it points at docs/SPOTMAX_VALIDATION.md, docs/BENCHMARK.md and the 23 real goldens.
This pack never computes or reads a test-split score of the RAD-51 benchmark.

## Limitations

The spot threshold was tuned on N2 without heat shock; re-check it on other genotypes and treatments. The axis
and germline stages assume one gonad arm in the field. SC lengths are uncalibrated (the sc:uncalibrated flag).

## Citation

See CITATION.cff. Pack text licence: MIT.
