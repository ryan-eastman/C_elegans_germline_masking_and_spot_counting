# Design: in-house foci detector and 3D nucleus model, with their benchmarks

Status: design from a 7-agent pass on 2026-10-09 (data and hooks inventory, detector and 3D-model survey, two designs, two judges, one synthesis). Not implemented. Both judges expect ties against the deployed tools; the plan is written so that a tie is still a publishable Resource result and a win is a pre-registered one.

Both judges expect ties, and both "first benchmark" claims are pre-empted (SpotMAX COSA-1, Toraason 2021, NIS3D, WORM, BlastoSPIM), so frame the paper as a Resource: nuclei and foci benchmarked on the same gonads, partial-annotation training at both stages, drop-in backends under a byte-identical default, and whether the swap changes the biology.

## Component 1: foci (log3d, net3d, assignment, harness)

**Design.** Package spots/backends/ with one detect() contract and a shared peaks.py: ellipsoidal NMS, SNR or probability gate, quadratic sub-voxel refinement, assign_nucleus with a 0.3 um fallback and a contested flag. log3d: anisotropic scale-normalised LoG, separable torch convolutions in 2 GB z chunks. net3d: RAD-51 plus DAPI residual 3D U-Net, Gaussian heatmap target, masked MSE, deterministic tiled inference. Corrections: effect_size must use the SpotMAX definition or a new column name; Spotiflow is bench-only, not a backend, unless it passes a compatibility test in the pinned venv.

**Training data.** Truth from the .ims only. Corrected totals: 7,370 noHS and 34,631 HS spots, 24,362 in two HS males. HS MALE_ and HS MALE_002 were run with region growing on: out of spot training and test until re-exported, leaving about 28k spots in 10 standard gonads. Class A/B is blocking (Class B may be over half the HS HERM spots): report A-only and A+B. Known region: germline Cellpose labels dilated 0.3 um (weight 1), traced masks (2), elsewhere 0.

**Benchmark and win.** Hungarian matching in an anisotropic cylinder (0.3 um xy, 0.4 um z); per-nucleus counts detector-only (Imaris mask) and end-to-end (our masks); CCC with cluster-bootstrap CIs. Leave-one-gonad-out is primary; test gonads confirm. Competitors: SpotMAX deployed and re-tuned on train plus val (its parameters touched test gonads), Spotiflow synth_3d and smfish_3d zero-shot and fine-tuned, RS-FISH, Imaris with one global Quality, and a blind two-scorer human subset sized by a power calculation. Win: end-to-end CCC gain over SpotMAX with the paired CI excluding zero and a pre-stated minimum gain, F1 not lower by more than 0.02, human agreement at least Imaris's own; a tie ships as an option.

**Integration.** spots.method defaulting to spotmax (config_hash 84bdda5a5230 untouched), one branch in _spots(); sub-hash payload and _STAGE_TOOLS additions conditional on method != spotmax, default sub-hash pinned in a test. Golden gate: 0 failing on all 12 n2rad51 tags plus dry-ice, and a synthetic golden with PSF-shaped foci.

**Effort: 16 weeks.** W1-2 pre-registration. W3-4 harness, both SpotMAX columns. W5-7 peaks.py, log3d, dispatch, tests, goldens. W8-9 competitors in .venv_bench, runtime and GPU memory. W10-13 net3d: loader, leave-one-gonad-out, learning curve, ablations, tolerance. W14-15 human anchor, rule applied once. W16 figures, docs, PR.

**Risks.** Thresholded automatic truth (15x Quality spread) and learning Imaris merges; Class B shrinking the HS HERM truth; fine-tuned Spotiflow winning.

**Novelty claim.** A RAD-51-specific, whole-gonad, nucleus-resolved foci benchmark scoring per-nucleus count agreement, a boundary-aware assignment step, and a nucleus-conditioned detector with partial-annotation masking as a drop-in backend.

## Component 2: nuclei3d

**Design.** Residual 3D U-Net (plain torch), xy-only first pooling, three heads: foreground, cell distance in um, neighbour distance (Scherr 2020); per-voxel supervision weights; h-maxima seeds, watershed, then the untouched _drop_small; deterministic inference. Corrections: 32 x 256 x 256 patches are 1.5 to 2 nuclei deep in z; nnU-Net v2's ignore label is prior art for the masked loss.

**Training data.** read_imaris_masks: crop padding, split collisions by connected components, match to SurfaceModelInfo (centroid within 0.1 um), inherit Class. Keep every surface-matched component regardless of size: six gonads had a 10-voxel floor, so a blanket 5 um3 cut deletes real surfaces. Weights: traced 1, background 1 where neither traced nor Cellpose-labelled, 0 inside untraced Cellpose labels plus margin. Complete napari crops (three held-out regions, one ccw77, one blind re-annotation by Madeleine), evaluation only.

**Benchmark and win.** instances3d.py: Hungarian IoU matching, F1 and AP over IoU 0.5 to 0.95, split and merge counts, z-extent error, CTC partial-support precision. Competitors: deployed Cellpose-SAM, Cellpose-SAM fine-tuned on identical masks (decisive), stock cpsam, classical watershed, StarDist-3D and EmbedSeg trained on the lab masks on Alpine, the 2.5D twin and ablations; public set if licensed. Win: paired F1 at IoU 0.5 and matched IoU gain of at least 0.03 with CI excluding zero, downstream RAD-51 CCC not worse; a tie with a runtime win is a tie.

**Integration.** backends.py registry, one additive branch after the legacy names, conditional stage_hashes payload, profile, models_registry entry, default sub-hash pinned. Never the default: it changes volumes and perinuclear shells downstream, documented per stage.

**Effort: 16 weeks.** W1-2 reader, zarrs, harness, goldens scored, pre-registration. W3-5 model, training, Cellpose fine-tune, full-stack RAM and watershed timing. W6-8 val-only iteration, learning curve, ablations, crops. W9-10 wiring, tests, golden diff, docs. W11-13 StarDist-3D, EmbedSeg, public set, ccw77 and dry-ice readouts. W14-15 final evaluation, downstream rerun, figures. W16 slack.

**Risks.** Tie at 1.85x anisotropy; Imaris convention is not the envelope; per-gonad subset bias; pseudo-label gate inherits Cellpose misses; 582 training nuclei is marginal.

**Novelty claim.** A compact anisotropy-aware 3D nucleus model trained from partial Imaris surfaces with a pseudo-label-gated ignore mask, a partial-ground-truth instance harness, and the first benchmark pairing nuclei and foci on the same gonads.

## Order and combined timeline

Unified split, frozen first: test noHS HERM_002, noHS MALE_002, HS HERM_002, HS MALE_001; val noHS HERM_001, HS HERM_001; train the rest. Phase 0, weeks 1 to 3, shared: one Imaris reader for spots and masks, manifest, per-gonad CreationParameters table, pre-registration docs, Madeleine's answers, NAS .ims check, per-stage runtime with current code (the goldens have no __stages.json), Cellpose tolerance. Phase 1, weeks 4 to 19: foci (truth readier). Phase 2, weeks 18 to 33: nuclei, overlapping the human-anchor weeks. Phase 3, weeks 34 to 38: paper and release. About nine months. Stop rule: if fine-tuned Spotiflow or Cellpose beats the in-house model on val by the minimum gain, stop iterating and report.

## Where the judges disagreed and my call

- C1 novelty 4 versus 6, risk 7 versus 5: the COSA-1 pre-emption stands, log3d likely ties, net3d uncertain; claim narrowed, a tie is publishable.
- Spotiflow flag: the README says --3d True; judge 1 is right.
- LoG memory and stage runtimes: budget judge 2's 2 GB; measure the stage split before quoting it.
- Human anchor: power calculation over a fixed count, two scorers either way.
- EmbedSeg: run it under CC BY-NC for benchmarking, no code into the backend.
- Conditional hash payloads: judge 2 is right, applied to both stages.

## Reviewer demands

1. One pre-registered split and frozen manifest for both components.
2. Per-gonad CreationParameters table; region-growing gonads out until re-exported.
3. Truth called threshold-curated Imaris detections; power-sized human anchor.
4. Narrowed novelty claims with prior art cited.
5. SpotMAX re-tuned on train plus val, both columns reported.
6. Missing competitors run head to head on the lab data.
7. Leave-one-gonad-out with cluster-bootstrap CIs and pre-stated minimum effects.
8. Convention separated from architecture; downstream conclusions tested with CIs.
9. Complete crops with inter-annotator agreement; per-gonad subset bias.
10. Determinism and runtime with current code, two runs per backend.
11. Release with DOI, consent, licence, CI, STAR methods, Limitations.
12. Factual slips corrected (do_3D, KIT-GE rank, 2002.02857, get-model, totals, filters, collisions).
