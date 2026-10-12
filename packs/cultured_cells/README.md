# Cultured cells pack (cultured_cells 0.1.0)

Nuclei of cultured cells are found on the DNA stain, and DNA repair foci (gamma-H2AX, and 53BP1 as a second
focus channel) are counted per nucleus in 3D. One command runs it:

    germquant run STACK.ome.tif --pack cultured_cells --out results/STACK

## Status: provisional, not validated

Version 0.1.0 is provisional. Every setting of the two recipes is a starting value
(docs/DESIGN_system_packs.md, section 3.12), not yet tuned. No recipe is validated: the pack has no validation
record, both recipes read `status: provisional` and `validated: false`, and no number they give has been
compared with ground truth yet. What is provisional:

- every value in config/profiles/nuclei_foci.yaml and nuclei_foci_3d.yaml: the projection, the smallest nucleus
  (20 um2, or 60 um3 in 3D), the focus radii, the Gaussian smoothing and the focus strength threshold;
- running Cellpose-SAM at the image's own scale: at 100x the nuclei are larger than Cellpose-SAM's usual object
  size (see "Limitations"); the object diameter is one of the settings tuned for 1.0.0;
- the expect ranges of the real example below (what a run is expected to give, not a score).

Tuning on a tune split of DeepFoci, validation with frozen protocols and confidence intervals, and model goldens
on the RTX 5090 come with version 1.0.0.

## What it is for, and not for

For: fixed cultured cells (cell lines, primary cultures) stained for DNA and one or two focus markers, imaged as
z-stacks with pixels of about 0.04 to 0.35 um and a z step of about 0.1 to 1 um (`pack check --image` warns
outside these ranges).

Not for: tissue sections; live-cell tracking; whole-cell (cytoplasm) segmentation; single-plane images (a
single-plane recipe comes in a later version); heavily overlapping nuclei in the projected recipe (use
nuclei_foci_3d).

## Recipes

| recipe | dims | what it does | channels it needs |
|---|---|---|---|
| nuclei_foci (default) | projected | Cellpose-SAM segments the maximum projection of the DNA channel; objects under 20 um2 are dropped; the masks are copied through every plane; SpotMAX counts foci per nucleus in 3D, the second focus channel as the spots instance foci2. Runs on a CPU. | dna and foci (foci2 optional) |
| nuclei_foci_3d | 3d | Cellpose-SAM segments the DNA channel in 3D (with the anisotropy dz / dy); objects under 60 um3 are dropped; foci as in nuclei_foci. A GPU is recommended. | dna and foci (foci2 optional) |

## Object kind: one row of the nuclei table is one nucleus

| object_kind | one row of the nuclei table is | n_nuclei counts | nucleus_id names |
|---|---|---|---|
| nucleus | a nucleus | nuclei | a nucleus |

The main foci are in `<image_id>__spots.csv` (the nuclei column `n_spots`), the second focus channel in
`<image_id>__spots_foci2.csv` (the nuclei column `n_spots_foci2`). In nuclei_foci the masks are 2D nuclei copied
through z, so `volume_um3` is empty and the nuclei table gains `area_um2` and `equivalent_diameter_um`
(docs/PLANAR.md).

## Channel naming

The channel map, config/channel_maps/cultured_dna_foci.yaml, binds roles by a part of the channel name (case does
not matter), then by position:

| role | names | if no name matches |
|---|---|---|
| dna | DAPI, Hoechst, DNA, SiR-DNA, DRAQ5, H2B | channel 0 |
| foci | gH2AX, gammaH2AX, H2AX, RAD51, foci, puncta | channel 1 |
| foci2 (optional) | 53BP1, foci2 | channel 2 |

A channel bound by position raises the advisory flag `channel:<role>:name_unmatched_used_index_<n>`. Laser-line
names (405, 488, 561, 640) are left out on purpose, so such channels bind by position. The recipes set
`segmentation.nuclei.require_channel`, so an image whose DNA channel cannot be bound fails with a message naming
the role. Labs whose channels carry other names (a protein name, a dye): bind the roles with the role editor of
the napari plugin (its project dock, a later plugin release), or until then point a study profile at a copy of
this map with their names added. Never edit the map in the pack: its text is part of the config hash.

## The focus threshold and its scale

`sharpen` is off in both recipes, so `effect_size` and `spots.effect_size_min` (3.0) are on SpotMAX's usual,
unsharpened scale (the spot-versus-background effect size on the smoothed image). If you switch `sharpen` on,
the effect size is computed on the sharpened image and the threshold must be chosen again.

## Model

`cpsam_v2@0f1cc3f7` (card: models/cpsam_v2.card.json): the stock Cellpose-SAM weights cpsam_v2, sha256
0f1cc3f7...c667, 1,233,586,851 bytes, the weights Cellpose 4.2 loads by default. They were not trained for this
pack. Weights licence: BSD-3-Clause; the Cellpose training data include CC BY-NC sets (and LIVECell, TissueNet and
the Cellpose data sets), so commercial use is an open question.

    germquant pack fetch-models cultured_cells

adopts Cellpose's own copy (~/.cellpose/models/cpsam_v2) after a size and sha256 check, or downloads the pinned
Hugging Face commit.

## How the starting values were chosen

Focus settings: a radius of 0.3 um in xy and 0.6 um in z, a Gaussian sigma of 0.08 um, the triangle threshold and
an effect size threshold of 3.0, the values validated for C. elegans RAD-51 foci (docs/SPOTMAX_VALIDATION.md),
with z columns merged. Nuclei: the maximum projection of the DNA stain, which suits flat cultured nuclei and runs
on a CPU, and a size gate of 20 um2 (a disk about 5 um across), or 60 um3 in 3D (a sphere about 4.9 um across).

## Examples

- `synthetic`: example/cultured_synthetic.ome.tif, four nuclei with 18 gamma-H2AX and 14 53BP1 foci, 16 x 144 x
  176 voxels at 0.3 x 0.2 x 0.2 um. The recipe goldens run it with the classical segmenter (no model needed) and
  compare every table cell with example/expected/recipe/<recipe>/.
- `deepfoci_valid_258`: one field of the DeepFoci foci validation subset (Zenodo 10.5281/zenodo.4067741, CC BY 4.0;
  Vicar et al. 2021, Computational and Structural Biotechnology Journal). It is never redistributed:

      germquant pack fetch-example cultured_cells
      germquant run IMAGE --pack cultured_cells --out results/deepfoci_258

  The first command reads one member, valid/data_258.tif (145 MB), out of the 10.2 GB zip with HTTP range
  requests, checks its size, CRC-32 and sha256, crops it to the two nuclei of the field and writes one OME-TIFF
  (3 channels 53BP1, gH2AX and DAPI, 48 planes of 480 x 896 pixels, 0.3 x 0.0647 x 0.0647 um), then checks its
  pixels_sha256. pack.yaml records how the member, the page order (ZYXC), the channel order and the voxel size
  were checked on the first download. No fallback example is declared: every member of the DeepFoci zip is
  stored without compression, so the member is always readable by range requests, and the design's fallback for
  a member that cannot be fetched (one S-BIAD1494 stack, open question O9) was not needed. If the deposit ever
  changes so that the member cannot be read, the fetch stops with a message naming why; download the whole
  archive from Zenodo and assemble the member with `germquant assemble` (docs/PACKS.md). On the lab workstation's CPU (2026-10-11) the run took 43 seconds and
  found 2 nuclei with 110 gamma-H2AX and 67 53BP1 foci; one of the two nuclei was only partly segmented (about
  97 of 157 um2), see "Limitations". The expect ranges are 1 to 3 nuclei and 35 to 80 foci per nucleus.

## Validation

None yet (see "Status"). The plan for 1.0.0: DeepFoci for nuclei and foci (its published splits, its two
validation experts as the ceiling, its 1.95 um matching radius plus a strict 3D score), with protocols committed
before any score and 95 percent cluster-bootstrap intervals.

## Limitations

- At 100x, nuclei span 150 to 300 pixels, larger than Cellpose-SAM's usual object size, and at the image's own
  scale Cellpose-SAM can clip or split them: on the DeepFoci example one nucleus came out at about 97 of 157 um2.
  Setting `segmentation.nuclei.cellpose.object_diameter_um: 12` in a study profile segmented both nuclei whole on
  that example (and ran five times faster); it changes the config hash and is the setting to tune for 1.0.0.
- In nuclei_foci a nucleus mask is a column through the stack: a focus above or below the nucleus in its column
  counts for it. Use nuclei_foci_3d for thick or overlapping nuclei.
- Anisotropy: typical 100x stacks are far coarser in z than in xy. nuclei_foci_3d runs Cellpose-SAM in 3D with
  the anisotropy dz / dy, which is about 4.6 on the DeepFoci voxel (0.3 / 0.0647 um); nothing has been checked
  at that ratio yet, so treat nuclei_foci_3d masks as provisional, and prefer nuclei_foci for flat cultured
  nuclei.
- Training-data overlap: the Cellpose-SAM training collection may include images from DeepFoci or BBBC039, the
  public sets these recipes are to be scored on. That overlap is checked (Cellpose-SAM's training lists against
  both sets) before any score is reported, and a score on overlapping images will say so.
- The DeepFoci test labels were clicked on 2D projections and matched with a lenient radius; DeepFoci scores
  will say so next to the numbers.
- QC for fields that are not gonads: with `qc.min_nuclei` 1 the advisory flag `qc:few_nuclei_<n>` never fires for
  a field with at least one nucleus (a DeepFoci field holds about 5), and the batch framing check is off
  (`qc.germline_outlier_factor` 1e9). `germquant batch FOLDER --pack cultured_cells --out OUT` runs the .nd2 files
  of a folder and collects them with the pack's QC settings (add `--resume` to collect again); folders of TIFF or
  OME-TIFF images come with projects in a later germquant package, so run those one at a time with `germquant run`
  for now. The plain `germquant collect` keeps the worm's framing check.

## Citation

See CITATION.cff: cite germquant and this pack, Cellpose-SAM and SpotMAX, and for the example the DeepFoci paper
and deposit (CC BY 4.0 asks for attribution). Pack text licence: MIT.
