# Test pack (test_pack 0.1.0)

A tiny outside pack for the germquant engine and plugin tests. It is not for real images.

## What it is for

- Checking the pack machinery end to end in CI: discovery through `$GERMQUANT_PACKS`, `germquant pack check`,
  `fetch-models`, `run --pack`, `batch --pack`, the recipe golden, lint and install.
- A stand-in for the outside packs: it extends its own frozen base (`config/profiles/pack_base.yaml`), switches
  every worm stage off, uses the neutral filename regex `(?!x)x` and the outside QC settings, and its recipes
  set `segmentation.nuclei.require_channel: true`, as every outside recipe does.

## Recipes

| recipe | dims | segmentation | spots | model |
|---|---|---|---|---|
| basic (default) | 3d | classical watershed | LoG (classical) | none |
| with_model | 3d | Cellpose | LoG (classical) | tinynet@f5055f8a, a stand-in |

`basic` needs only the core dependencies, so it runs on any machine. `with_model` names a pinned pack model
that is 64 bytes of text (see `models/tinynet.card.json`, field `stand_in_content`, plus one line feed): it
exercises `pack status`, `fetch-models` and the pin check, but Cellpose cannot load it. Its recipe golden
forces the classical segmenter, which never resolves the model, so the golden runs with an empty model store.

## Object kind

One row of the nuclei table is one nucleus.

## Example

`example/test_pack_synthetic.ome.tif`: a deterministic two-channel field (DAPI and foci), 12 x 64 x 96 voxels
at 0.5 x 0.2 x 0.2 um, five nuclei with two puncta each, written by `scripts/make_test_pack.py`. The recipe
goldens live under `example/expected/recipe/<recipe>/`.

## Limitations

Synthetic data only; no validation record; the spot detector is the experimental LoG baseline (every run
carries the flag `spots:unvalidated_backend=log`).

## Citation

Not citable: test data. Licence MIT.
