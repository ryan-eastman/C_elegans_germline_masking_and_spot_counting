# SC tracing for fragmentation (the `sc_trace` stage)

Restored 2026-09-10 from the June 2026 upstream tracer (`C_elegans_ml`, branch `imaris-validation`,
`src/germquant/sc/skeleton.py`, byte-identical) as an optional stage of the modular pipeline
(docs/ROADMAP_modular_pipeline.md, step 6). Switch: `sc.trace.enabled` (off by default; profile
`n2_sc` turns it on; `--no-sc_trace` switches it off on the command line).

## What it measures, per germline nucleus

The SYP channel is cropped per nucleus, resampled to isotropic voxels, ridge-filtered (Sato, physical
scales `sc.ridge_sigmas_um`), turned into a continuous strand mask by hysteresis thresholding inside
the nucleus (`sc.ridge_hyst_low_pct` / `sc.ridge_hyst_high_pct`), smoothed (in-plane closing, holes
under 27 voxels filled, components under 30 voxels dropped), skeletonized in 3D and measured with skan.

| column (nuclei table) | meaning | trust |
|---|---|---|
| `sc_total_length_um` | summed skeleton length of all SC fragments | recoverable: on synthetic truth 34.7 v 35.2 um before mask smoothing, about 7 percent shorter after, per-nucleus correlation about 0.6 |
| `sc_n_fragments_lb` | number of disconnected skeleton components | a LOWER BOUND, not the biological count: the six SCs overlap in 3D at pachytene density, components cap near 2 to 3, four counting methods all correlate below 0.15 with truth |
| `sc_fragmentation_index` | SYP intensity coefficient of variation within the nucleus | noisy per nucleus (corr 0.21), separates control from heat at the population level (Cohen's d about 0.57); report per gonad or per condition, never as a per-nucleus classifier |
| `sc_expected_n_tracks` | 6 for oocytes, 5 for spermatocytes (X univalent) | from the sample name |

Tables `sc_tracks` (one row per fragment: length, branches, junctions, tortuosity) and
`sc_per_nucleus` (the full per-nucleus aggregate) are written on every run, empty when the stage is off.
`image_summary` gains `mean_sc_total_length_um`, `mean_sc_fragmentation_index`,
`mean_sc_n_fragments_lb` over germline nuclei. Pachytene-restricted means come from the staging stage
(hand traces) once that lands.

## Calibration status

Synthetic ground truth only (`scripts/validate_sc_tracer.py`, `scripts/synth_nuclei.py`). No Imaris
SC ground truth exists yet; the E: drive exports are RAD-51 spots, nucleus surfaces and centerlines.
Until Filament traces of 15 to 25 pachytene nuclei per condition are compared (roadmap step 13:
`germquant validate --mode sc`), every run carries the `sc:uncalibrated` flag and absolute lengths
should not be quoted. Calibration will be recorded as a bias with a concordance coefficient, never
baked in as a length multiplier.

## Two mask cores, deliberately separate

`sc.enabled` still means the July ribbon CONTROL operand inside the coloc stage
(`germquant.sc.surface.surface_sc_ribbon`: size gate 10 voxels, no hole filling); it has nothing to do
with tracing and is known to leak into the perinuclear shell. `sc.trace.enabled` is the tracer. They
are not unified until the tracer is calibrated; the rename of `sc.enabled` to
`coloc.sc_ribbon_operand` waits for the deliberate config-hash bump release.

## Ribbon straightening: a separate, opt-in command

`germquant sc-straighten` (docs/SC_STRAIGHTEN.md) traces each SC ribbon in 3D on finished results,
straightens it along its centreline and measures it (length, profile, gaps, width, foci positions). It
is a post-run command that writes only new files; this stage, its tables and its switch are unchanged.
On real data this stage's `sc_total_length_um` follows nucleus volume and background rather than SC
content (golden n2rad51_noHS_herm_000, untuned, 2026-10-23): the percentile thresholds keep a fixed
volume fraction of every nucleus (about 43 percent of the eroded interior), the dimmest nuclei get the
longest skeletons and a nucleus filled with background-level noise yields about 190 um; read it as
uncalibrated.
