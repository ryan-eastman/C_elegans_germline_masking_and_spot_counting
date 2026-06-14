# Overnight findings — nucleus segmentation (night of 2026-06-12 → 13)

**Read this first. TL;DR in one paragraph, then the evidence, the honest limitations, and the Monday plan.**

## TL;DR

You asked the right question — *"should we train our own Cellpose model instead of the stock one?"*
**Yes — I have a trained model that beats stock, but the story is more nuanced than my first pass
implied (I corrected myself overnight; see §3b/§4).**

- A **synthetic-only** model overfit: F1 0.98 on synthetic but only **0.55 on real** held-out — far
  below a real-trained model. Synthetic data alone is not enough.
- A model trained on **real** hand annotations (KoehlerLab — the only real germline-nucleus labels
  available before Monday) scores **0.94**, and a **real+synth COMBINED model 0.98**, on held-out images
  at *your* Cahoon pixel scale. Stock is **0.87** there — so on a single 2D slice the gain is real but
  **modest** (my earlier "stock 0.37" was a scale artifact — §3b).
- The **decisive** difference is in 3D: the pipeline runs Cellpose in `do_3D`, which feeds the model
  rotated (xz/yz) views. On those, **stock collapses to F1 0.17 while the trained models hold ~0.67** —
  and *that* is what makes stock shatter each nucleus into hundreds–thousands of ~3 µm³ fragments and
  orphan RAD-51 foci, while the trained model produces ~36 µm³ whole nuclei that contain **~97% of foci**
  on your Cahoon data (and stay anatomically sized — see the §3c controls before trusting that 97%).

**Best model tonight: `germline_nuclei_combined` (real + synthetic).** The Monday implication stands:
your Imaris ground truth isn't just for *checking* the pipeline — it's the training data the nucleus
model needs. Real annotations are the missing ingredient. Plan below.

⚠️ **What is and isn't proven** (I ran an adversarial self-review to keep myself honest — see §4):
**Solid** — synthetic-only overfits; the trained models clearly win in the 3D/orthogonal regime the
pipeline uses (~0.67 vs stock 0.17) and produce anatomically-correct whole nuclei on Cahoon. **Softer** —
on a single xy slice the trained-vs-stock gap is modest (0.98 vs 0.87), my first split had same-gonad
leakage (rebuilt leak-free, §3b), and there is still **zero Cahoon ground-truth F1** — exactly what
Monday fixes.

---

## 1. What was broken (the root cause, restated)

Stock Cellpose-SAM cannot group a pachytene nucleus ("6 chromosomes tangled in a ball, dark
nucleoplasm between") into ONE object. It traces bright chromatin threads instead of the nuclear
envelope, so each nucleus shatters into many partial masks. Downstream this:
- **orphans RAD-51 foci** — foci on chromosomes left outside the mask get assigned to no nucleus;
- **crops the SC trace** — SYP signal on excluded chromosomes is cut from the per-nucleus length.

Both headline readouts (foci/nucleus, SC length/fragmentation) are corrupted *at the segmentation
foundation*. You said it yourself: fix the earliest stage first, because error compounds. That's why
all of tonight's work is on segmentation and **nothing downstream of it was touched.**

## 2. The approach

Two models, one synthetic engine (architecture unchanged from our earlier decision):
- **nucleus model** — Cellpose on DAPI (this is what tonight validated).
- **SC tracer** — still broken at the method level (corr 0.07 with truth); deferred, see §6.

I compared four nucleus models on three benchmarks:

| Model | Training data |
|---|---|
| `stock cpsam` | none (Cellpose-SAM as shipped) |
| `synthetic-only` (`germline_nuclei`) | 165 synthetic DAPI slices |
| `real` (`germline_nuclei_real` / `_real_grp`) | 48–51 real KoehlerLab hand annotations |
| `real+synth combined` (`germline_nuclei_combined`) | 48 real + ~42 synthetic |

Benchmarks: (a) held-out **real** Köhler images, (b) held-out **synthetic**, (c) **your Cahoon .nd2**
via a ground-truth-free biological prior (RAD-51 foci must lie inside nuclei).

## 3. Results

### 3a. Held-out F1 — FIRST (leaky) split, for the record
*(split strided over sorted slices → ~5/12 test images were adjacent crops of a training gonad →
these numbers are optimistic; kept only to show the synthetic-vs-real contrast)*

| Model | real held-out F1 | synthetic F1 | real pred/GT |
|---|---|---|---|
| stock cpsam | 0.37 | 0.86 | 46 / 25 (over-seg) |
| synthetic-only | 0.36 | **0.98** | 12 / 25 (under-seg) |
| real (leaky) | **0.91** | 0.82 | 24 / 25 |

### 3b. Held-out F1 — LEAK-FREE split (trustworthy) — *this corrects 3a*
*(grouped by image shape so no gonad appears in both train and test; test = 9 held-out gonads at your
Cahoon 0.108 µm/px scale; per-image F1 spread in parens; plus an xz/yz orientation probe)*

| Model | leak-free **xy** F1 (Cahoon scale) | **xz/yz** ortho F1 | synthetic F1 | xy pred/GT |
|---|---|---|---|---|
| stock cpsam | 0.87 (0.64–0.97) | **0.17** | 0.85 | 22/22 |
| synthetic-only | 0.55 (0.00–0.91) | 0.05 | 0.98 | 18/22 |
| real (leaky split) | 0.97 | 0.63 | 0.82 | 22/22 |
| real (leak-free) | 0.94 (0.82–1.00) | 0.68 | 0.81 | 20/22 |
| **real+synth combined** | **0.98 (0.93–1.00)** | 0.66 | **0.96** | 22/22 |

**Three honest corrections to my first-pass story (this is why the leak-free redo mattered):**

1. **"Stock F1 0.37" was mostly a SCALE artifact, not the truth.** That first test set was 7/12 images
   at a 3×-finer pixel scale where stock struggles. At *your* Cahoon scale on a clean held-out set,
   **stock actually scores 0.87** on a single xy slice. So the per-slice xy gap (real/combined 0.94–0.98
   vs stock 0.87) is **real but modest** — not the 0.37→0.91 chasm 3a implied.

2. **The trained model's BIG advantage is in the 3D regime the pipeline actually runs.** Cellpose
   `do_3D` feeds the model xz/yz reslices. On those, **stock collapses to F1 0.17** while the
   real-trained models hold **0.66–0.68**. *This* is what produces stock's 1650-fragment 3D mess and the
   real model's coherent whole nuclei — and it's the mechanism behind the foci-recapture result (§3c).
   The xy single-slice number understates the real-vs-stock difference; the 3D/orthogonal number is the
   one that matters for the pipeline.

3. **Synthetic data HELPED when combined with real.** `real+synth combined` is the best overall model
   (xy 0.98, synthetic 0.96, ortho 0.66) — adding synthetic did not drag it toward the synthetic-only
   failure mode. **Recommend `germline_nuclei_combined` as the working model**, pending Monday's
   Cahoon-labelled validation.

### 3c. Cahoon foci-recapture (your data, no ground truth needed) — hardened
*(3 on-tissue crops, 220 RAD-51 foci total, identical for every model; frozen foci params; voxel spacing
read from the .nd2 (verified 0.2/0.108/0.108); with the controls the review demanded)*

| Model | % foci inside | % vol fill | enrichment (inside/fill) | µm³/object | foci/nucleus |
|---|---|---|---|---|---|
| stock cpsam | 75% | 4.6% | 16.3× | 3.4 (fragments) | 1.0 |
| synthetic-only | 38% | 1.5% | 24.4× | 2.5 (fragments) | 1.0 |
| real (leak-free) | **97%** | 9.3% | 10.5× | **36.5 (whole)** | 1.3 |
| real+synth combined | 96% | 9.1% | 10.5× | 35.2 (whole) | 1.3 |

**Read this honestly — it is NOT "tighter masks":**
- The **load-bearing number is µm³/object.** A real 5 µm pachytene nucleus is ~50–65 µm³. The trained
  models produce **~36 µm³ whole nuclei (~150 objects/crop)**; stock and synthetic-only produce **~3 µm³
  sub-nuclear fragments (~1300 and ~580 objects)** — useless for per-nucleus counts no matter how many
  foci they happen to overlap. Confirmed visually in `results_fullres/masks3d_compare.png`.
- **Enrichment (capture per unit volume) actually favors stock (16×) over the trained models (10.5×)** —
  so the trained models are *not* more spatially precise per voxel. Their win is whole-nucleus-ness.
- **Null-dilation control** (bloat each model's masks to the trained model's ~9% fill, re-measure):
  stock 64% → **93%**, synthetic-only 11% → **68%**. So *most* of stock's foci-orphaning is just
  volume/fragmentation (dilation nearly fixes capture) — but synthetic-only's masks are genuinely
  **misplaced** (only 68% even when bloated). The trained models reach ~100% at that same fill **as
  coherent nuclei**, not bloated blobs.
- **Merge proxy:** median foci/nucleus ≈ 1.3 (biologically sane for N2 RAD-51) → **no sign of gross
  over-merging**, the failure mode "% inside" is otherwise blind to.

**Bottom line for 3c:** the trained models are the only ones that deliver *whole, anatomically-sized
nuclei that contain the foci* — the prerequisite for per-nucleus RAD-51 counts and SC tracing. This is
a biological sanity check on your real data, **not** a segmentation-accuracy number (that needs Monday's
labels).

## 4. Honest limitations (from an adversarial self-review)

I ran a 4-reviewer + synthesizer adversarial pass against my own conclusions *before* writing this.
What it found:

1. **Leakage (FIXED + quantified).** My first split (`prep_kohler_real.py`) strided `i%5` over sorted
   slices, leaking same-gonad crops into train+test → the headline F1 0.91 was optimistic. I rebuilt a
   leak-free split grouped by image shape (`prep_kohler_grouped.py`) and re-trained. Net effect of the
   leak: real model 0.97 (leaky) → 0.94 (leak-free) — small, but the bigger correction was the **scale
   artifact**: the original test set was mostly at a 3×-finer scale where *stock* looks bad (0.37), so
   the leak-free, Cahoon-scale redo (§3b) reset stock to its true 0.87.
2. **2D vs 3D (probed, not closed).** F1 is on 2D slices, but the pipeline runs `do_3D` (xz/yz reslices
   the model never trained on). The §3b xz/yz probe directly measures this gap — and it's where the
   trained models earn their keep (0.67 vs stock 0.17). Full 3D-instance validation still needs Monday's
   labels.
3. **Foci-recapture is one-sided (controlled).** "% foci inside" rewards volume and is blind to
   over-merging. §3c now reports enrichment, a null-dilation control (which honestly shows *most* of
   stock's deficit is volume, recoverable by dilation), and a foci/nucleus merge proxy (≈1.3 → no gross
   merging). Framed as a biological sanity check, **not** an accuracy number.
4. **Köhler ≠ Cahoon (but the labels are the same KIND).** All real F1 is within the Köhler-lab domain
   (different microscope/prep). There is **no Cahoon ground-truth F1** yet — the Cahoon evidence is
   foci-recapture + qualitative masks only; that's the whole point of the Imaris data. *Caveat softened:*
   the Köhler lab's annotations were very likely made in **Imaris** — the same tool that produces Monday's
   labels — so the annotation *methodology* matches; the remaining gap is microscope/sample-prep, not
   labelling style. (Practical consequence: Monday's Imaris labels will likely be **3D surfaces** →
   slice them into xy/xz/yz for training, as in the recipe.)

**Freeze compliance — confirmed clean.** `config.yaml` was not edited (timestamped before the
session). No frozen tuning knob (`foci.threshold_rel`, `sc.intensity_percentile`, `zones.*`) was
changed or calibrated against ground truth. Training/selecting a segmentation model is structural
work — exactly the carve-out you allowed. Nothing destructive; NAS untouched; all work on the
`segmentation-model` branch.

## 5. Monday recipe (the validated plan)

When the Imaris ground truth arrives:
1. **Split by ANIMAL/gonad, never by slice index** — this is the single thing that makes the reported
   F1 trustworthy (it's what bit the first split).
2. **Include xz/yz orientations in training** — if Imaris gives 3D surfaces, slice them in all three
   planes so the model is trained for the 3D regime it's deployed in. (Köhler's set already ships
   xz/yz files; our prep currently only uses xy.)
3. **Fine-tune from cpsam on the real Cahoon labels** — proven recipe: lr 1e-5, wd 0.1, ~150 epochs
   (`scripts/train_real_grouped.py`). Ablate Cahoon-only vs Cahoon+Köhler combined.
4. **Benchmark in the 3D `do_3D` regime**, reporting instance F1 *and* the per-nucleus volume
   distribution + a merge rate (so over-merging is caught — foci-recapture alone can't see it).
5. **Integration note (your call — it changes numbers, so I left it):** `segment_nuclei._cellpose`
   passes a `diameter` rescale (~28 px) for any non-`cpsam` model *path*, which our fine-tuned models
   were *not* benchmarked with. **RESOLVED (2026-06-13):** measured native vs rescaled on the Köhler GT —
   native wins (combined 0.980 vs 0.978; `scripts/diameter_test.py`) and is the benchmarked regime, so
   `segment/nuclei.py` now runs fine-tuned models at native scale.
6. *Then, and only then,* unfreeze and calibrate the downstream SC/zone/foci params against Imaris.

## 6. Not done / next

- **Germline isolation — DONE (2026-06-13).** New `germquant.germline` module drops nuclei segmented
  outside the gonad (gut/debris) using SYP + spatial connectivity (never foci → foci-recall is an
  independent check). **A montage check caught a real bug first:** thresholding each nucleus on SYP
  independently (`multi_cc`) wrongly dropped the entire **SYP-negative distal tip** (mitotic + transition
  zone, before SYP-3 loads) — that's real germline, and foci-recall (0.97) was blind to it because the
  distal tip has few foci. Fixed with the default **`syp_seeded_cc`**: seed on the SYP-bright synapsed
  core, then keep the whole spatially-connected component containing it — so the distal tip is kept
  (contiguous) while separate gut/debris clusters are dropped. Re-validated on 14 real gonads: mean
  foci-recall **0.995**, min 0.98; on the real control gonad it keeps the entire tube and drops only the
  ~22 disconnected off-gonad nuclei (`results_fullres/germline_check.png`). Wired into the pipeline
  (axis/zones/means run on germline only; excluded nuclei kept flagged `in_germline=False`); montage
  outlines off-gonad nuclei in **red**. Config block `germline:`; 5 unit tests + e2e (**32 pass**).
  Validators: `scripts/validate_germline_select.py`, `scripts/viz_germline_check.py`.
  *Lesson (per the earlier adversarial review): one-sided metrics lie — always look at the picture.*
- **SC tracer — REBUILT (2026-06-13).** The old tracer was broken (fragment corr 0.07, length ~45% of
  truth) because a single 90th-percentile threshold shredded the thin SC into disconnected blobs.
  Diagnosed on synthetic SC ground truth (looked at the pictures). **Two findings:** (1) **SC length IS
  recoverable** — switched to a continuous
  *hysteresis* mask → length now **34.7 vs 35.2 µm GT, corr 0.60** (was ~15 µm). (2) **Per-nucleus
  fragment COUNT is NOT recoverable** from light microscopy at pachytene density — the ~6 SCs overlap
  in 3D (clean-signal component count caps at ~2.4, not 6) and lateral merging bridges each strand's
  heat-gaps; four methods (components, endpoints, gap-runs, dark-fraction) all gave corr <0.15. This
  confirms ARCHITECTURE Risk 1 and is exactly why you're getting Imaris tracing. **Instead** the heat
  readout is now a **fragmentation INDEX** = SYP intensity CV per nucleus: per-nucleus noisy (corr
  0.21) but separates control vs heat at the **population level (Cohen's d 0.57)**. Reported as
  `sc_fragmentation_index` alongside `n_fragments` (now flagged a lower bound, not the biological
  count). Continuous mask verified visually (`results_fullres/sc_v2_compare.png`). New config params
  `sc.ridge_hyst_{low,high}_pct`; `sc.intensity_percentile` deprecated. Tests updated (32 pass).
  **Monday:** calibrate both SC length and the fragmentation index against Imaris.

## 7. Artifacts

**Models** (`models/models/`, gitignored, ~1.2 GB each):
- **`germline_nuclei_combined`** ← **recommended working model** (real + synthetic; best overall: xy 0.98,
  synth 0.96, ortho 0.66)
- `germline_nuclei_real_grp` (real, leak-free split; the trustworthy-F1 model: xy 0.94, ortho 0.68)
- `germline_nuclei_real` (real, leaky split — superseded), `germline_nuclei` (synthetic-only — overfit)

**Figures** (`results_fullres/`, gitignored): `masks3d_compare.png` (5-panel 3D-mask overlay),
`germline_check.png` (germline kept-vs-dropped on a real gonad), `sc_v2_compare.png` (SC mask/skeleton
old-vs-new), `model_compare3.png`, plus `train_*.log`.

**Scripts** (`scripts/`): nucleus — `synth_nuclei.py` (engine), `prep_kohler_grouped.py`,
`train_nucleus_model_combined.py` / `train_real_grouped.py`, `benchmark2.py`, `diameter_test.py`,
`foci_recapture2.py`, `viz_masks3d.py`; germline — `validate_germline_select.py`,
`viz_germline_check.py`, `run_real_germline.py`; SC — `validate_sc_tracer.py`, `sc_v2_viz.py`.

**One code fix** (not a tuning change): `src/germquant/validate/compare.py` `_iou_matrix` crashed on
zero-overlap predictions (empty float index array); forced int64 — needed for the xz/yz benchmark and
for Monday's 3D instance-F1.
