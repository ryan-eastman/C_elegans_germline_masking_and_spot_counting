# Budding yeast pack (budding_yeast 0.1.0)

Whole budding yeast cells are found on one transmitted-light plane, and spots (smFISH probes, fluorescent
protein foci) are counted per cell in 3D. One command runs it:

    germquant run STACK.ome.tif --pack budding_yeast --out results/STACK

## Status: provisional, not validated

Version 0.1.0 is provisional. Every setting of the recipe is a starting value (docs/DESIGN_system_packs.md,
section 3.12), seeded from the SpotMAX authors' own settings and not yet tuned. No recipe is validated: the
pack has no validation record, the recipe reads `status: provisional` and `validated: false`, and no number it
gives has been compared with ground truth yet. What is provisional:

- every value in config/profiles/cells_spots.yaml: the plane picker (sharpest), the smallest cell (3 um2), the
  spot radii, the Gaussian smoothing, the spot strength threshold and the sharpening;
- the expect ranges of the real example below (what a run is expected to give, not a score);
- the use of the stock Cellpose-SAM weights for yeast cells (no yeast fine-tune).

Tuning on a tune split of S-BIAD3750, validation with frozen protocols and confidence intervals, and model
goldens on the RTX 5090 come with version 1.0.0.

## What it is for, and not for

For: z-stacks of budding yeast with a transmitted-light channel (phase contrast, bright field or DIC) and one
spot channel, with pixels of about 0.05 to 0.15 um and a z step of about 0.15 to 0.5 um (`pack check --image`
warns outside these ranges).

Not for: stacks without a transmitted-light channel (the image fails, see "Channel naming"); single-plane
images (a single-plane recipe comes in a later version); tracking cells over time; thick samples or cells that
lie on top of each other (the masks are columns through the whole stack); nuclear segmentation.

## Recipes

| recipe | dims | what it does | channels it needs |
|---|---|---|---|
| cells_spots (default) | projected | Cellpose-SAM segments the sharpest plane of the transmitted channel; the whole-cell masks are copied through every plane; SpotMAX counts spots per cell in 3D | transmitted and foci (dna optional) |

## Object kind: one row of the nuclei table is one cell

Table and column names never change, so in this pack:

| object_kind | one row of the nuclei table is | n_nuclei counts | nucleus_id names |
|---|---|---|---|
| cell | a whole cell | cells | a cell |

Because the masks are 2D cells copied through z, `volume_um3` is empty (a column is not a cell volume) and the
nuclei table gains `area_um2` and `equivalent_diameter_um` from the 2D mask (docs/PLANAR.md).

## Channel naming

The channel map, config/channel_maps/yeast_transmitted_foci.yaml, binds roles by a part of the channel name
(case does not matter):

| role | names | if no name matches |
|---|---|---|
| transmitted | phase, contr, PhC, Bright, DIC, trans | nothing: the channel is missing |
| foci | smFISH, FISH, Quasar, Cy3, Cy5, mNeon, GFP, YFP, foci, spots | channel 1, with an advisory flag |
| dna (optional) | DAPI, Hoechst | nothing |

The transmitted role never binds by position. A stack without such a channel gets the QC failure
`channel:transmitted:MISSING_REQUIRED`, and because the recipe sets `segmentation.nuclei.require_channel`, the
image fails with a message naming the role instead of segmenting channel 0, which could be a fluorescence
channel. The name BF is left out on purpose, because it would also match BFP.

Labs that name channels after genes (MDN1, ACT1) will find that no name matches: bind the roles with the
role editor of the napari plugin (its project dock, a later plugin release), or until then point a study profile
at a copy of this map with their names added. Never edit the map in the pack: its text is part of the config
hash.

## The spot threshold and its scale

`sharpen` is on, as in the SpotMAX authors' runs, so SpotMAX computes the spot-versus-background effect size on
the sharpened (difference of Gaussians) image: the `effect_size` column and `spots.effect_size_min` (2.5) are on
the sharpened scale. A threshold chosen without sharpening, such as the C. elegans RAD-51 value 3.0, does not
carry over; change `effect_size_min` and `sharpen` together.

## Model

`cpsam_v2@0f1cc3f7` (card: models/cpsam_v2.card.json): the stock Cellpose-SAM weights cpsam_v2, sha256
0f1cc3f7...c667, 1,233,586,851 bytes, the weights Cellpose 4.2 loads by default. They were not trained for this
pack. Weights licence: BSD-3-Clause; the Cellpose training data include CC BY-NC sets, so commercial use is an
open question. YeaZ yeast images are part of Cellpose-SAM's training data (sampled at 1 to 2 percent), so yeast
segmentation scores of this model are not independent of YeaZ data.

    germquant pack fetch-models budding_yeast

adopts Cellpose's own copy (~/.cellpose/models/cpsam_v2) after a size and sha256 check, or downloads the pinned
Hugging Face commit. On the CPU, the bare name `cpsam` gave masks identical to `cpsam_v2@0f1cc3f7` on the plane
this recipe segments in the real example (docs/PACKS.md records the check); the same check on CUDA is pending.

## How the starting values were chosen

From the authors' SpotMAX settings for the SUT509 smFISH positions of S-BIAD3750: spot radii (1.0, 0.244,
0.244) um, a Gaussian sigma of 0.75 pixels at 0.072 um (0.054 um), the triangle threshold, an effect size
threshold of 2.5 on the sharpened image, hot-pixel removal, and no merging of z columns. The deposit's own
cells were segmented by Cell-ACDC on one plane of the transmitted channel, which is why the recipe segments one
plane; it takes the plane in best focus and drops objects under 3 um2.

## Examples

- `synthetic`: example/yeast_synthetic.ome.tif, 13 x 160 x 288 voxels at 0.24 x 0.072 x 0.072 um: six cells
  with 20 spots, one of them weak (its effect size on the sharpened scale, about 3.2, lies just above the
  threshold of 2.5), two small bright specks between the cells that the segmenter does not find, and a small
  object of about 1.9 um2 that it finds and the 3 um2 area gate drops. The recipe golden runs it with the
  classical segmenter (no model needed) and compares every table cell with example/expected/recipe/cells_spots/,
  so a change of the spot threshold or of the smallest cell shows in the tables.
- `sbiad3750_sut509_p1`: SUT509 smFISH, Position_1 of the BioImage Archive deposit S-BIAD3750 (CC0; Padovani et
  al. 2026, Science Advances). It is never redistributed:

      germquant pack fetch-example budding_yeast
      germquant run IMAGE --pack budding_yeast --out results/sut509_p1

  The first command downloads two 14.7 MB TIFFs (Bright and Cy3), checks their size and sha256, stacks them into
  one OME-TIFF (2 channels, 15 planes of 700 x 700 pixels, 0.24 x 0.0720635 x 0.0720635 um) and checks its
  pixels_sha256; it prints the image path to run. On the lab workstation's CPU (2026-10-11) the run took 40
  seconds and found 11 cells and 8 spots; the deposit's two annotators marked 6 and 9 spots in this position.
  The expect ranges are 8 to 15 cells and 0.3 to 1.5 spots per cell.

## Validation

None yet (see "Status"). The plan for 1.0.0: S-BIAD3750 for spots (against both annotators) and BABY for cells,
with splits and protocols committed before any score, 95 percent cluster-bootstrap intervals, and the scores of
SpotMAX's own outputs as comparators.

## Limitations

- A cell's mask is a column through the stack: a spot above or below the cell in its column counts for it.
- Training-data overlap: YeaZ images are in Cellpose-SAM's training data, the S-BIAD3750 masks are
  hand-corrected YeaZ output, and SpotMAX's AI was trained on annotations from the same groups.
- germquant does not run SpotMAX's spotFIT step, so the authors' deposited SpotMAX outputs are a comparator, not
  a target to match exactly.
- QC for fields that are not gonads: the advisory flag `qc:few_nuclei_<n>` fires below 10 cells
  (`qc.min_nuclei`), and the batch framing check is off (`qc.germline_outlier_factor` 1e9). `germquant batch
  FOLDER --pack budding_yeast --out OUT` runs the .nd2 files of a folder and collects them with the pack's QC
  settings (add `--resume` to collect again); folders of TIFF or OME-TIFF images come with projects in a later
  germquant package, so run those one at a time with `germquant run` for now. The plain `germquant collect` keeps
  the worm's framing check.

## Citation

See CITATION.cff: cite germquant and this pack, Cellpose-SAM and SpotMAX, and for the example S-BIAD3750. Pack
text licence: MIT.
