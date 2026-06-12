# Architecture — *C. elegans* germline 3D quantification

> Design of record. Produced from a web-grounded, adversarially-verified research pass
> (15 agents, ~540k tokens). Citations are inline. Decisive throughout; inference flagged.

**Scope:** confocal Nikon `.nd2` Z-stacks → unattended NAS-folder batch → tidy CSV/Parquet
for R/Positron + publication renders. Core readouts: **per-nucleus SC length & fragmentation**,
**transition-zone vs pachytene zone length**, **RAD-51 foci/nucleus**; granule 3D metrics and
crossover foci (COSA-1/MSH-5/ZHP-3) as optional modules.

**This dataset (confirmed from the files):** 3-channel uint16 stacks, Z≈67–91, 2600×2600,
voxel **0.108 × 0.108 × 0.20 µm** (z anisotropic ≈1.85× xy). Channels named by laser line →
`405`/em438 = **DAPI** (DNA), `477`/em511 = **SYP-3** (central element / SC), `545`/em595 =
**RAD-51** (DNA-damage foci). Filenames encode design: genotype (N2), sex (HERM=oocytes /
MALE=spermatocytes), treatment (HS / noHS), replicate. Biology (Cahoon lab): heat fragments
the SC in spermatocytes but not oocytes — sexually dimorphic.

---

## 1. Verdict: Imaris vs open Python

**Build the load-bearing, unattended, provenance-tracked pipeline entirely in open-source
Python on the RTX 5090 + SLURM. Keep Imaris off the critical path** — used only for (a) one
validation cross-check on SC tracing and (b) optional hero figures.

Driven by an architectural fact, not preference: **Imaris is Windows/macOS-only (no Linux)**,
so it cannot run headless on a SLURM node, and its "Batch Process" is GUI point-and-click, not
a CLI ([System Requirements](https://imaris.oxinst.com/support/system-requirements);
[Imaris on NIH HPC](https://hpc.nih.gov/apps/imaris.html)). Its ImarisXT/ICE bridge needs the
GUI running plus legacy Python 2.7 ([cvbi python-XTensions](https://cvbi.github.io/python-XTensions/)).
That defeats "point at NAS → run unattended."

Open tooling matches Imaris on every quantitative readout. The adversarial check found
"Imaris is *the* standard for SC tracing" **overstated**: peer-reviewed *C. elegans* SC work
also traces with open **Fiji Simple Neurite Tracer (SNT)**
([Sci Adv 2024, PMC11753403](https://pmc.ncbi.nlm.nih.gov/articles/PMC11753403/)), and `skan`
skeleton lengths validate against Fiji AnalyzeSkeleton
([Nunez-Iglesias et al., PeerJ, PMC5816961](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5816961/)).
Imaris keeps a genuine edge only in **semi-automated 3D tracing of tangled SC in intact
whole-mount pachytene nuclei** ([Hurlock et al., JCB 2020, PMC7199856](https://pmc.ncbi.nlm.nih.gov/articles/PMC7199856/);
[Cahoon et al., eLife 2023](https://elifesciences.org/articles/84538)); the leading open
dedicated SC tool is explicitly 2D-spread-only and *not* for intact *C. elegans* nuclei
([PMC9894712](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9894712/)). → validate the open
tracer against an Imaris/SNT subset (Risk 1).

---

## 2. End-to-end stages

Voxel spacing `(dz,dy,dx)` is read once at stage 0 and is a **required argument** to Cellpose
`anisotropy`, `skan` `spacing`, `marching_cubes` `spacing`, `regionprops` voxel-volume, and
foci-radius conversion. A single forgotten `spacing` silently corrupts every length/area/volume.

| # | Stage | Tool / method | Why |
|---|-------|---------------|-----|
| 0 | Read `.nd2` + metadata | **`nd2`** (Talley Lambert); bioio-nd2 optional | exposes `voxel_size()`, channel names, `to_dask()`; aicsimageio archived → bioio. Pull `(dz,dy,dx)` once, thread everywhere. |
| 0b | Channel-role resolution | config name→role map, validated at load | channels vary; detect DAPI presence; missing expected channel = QC flag |
| 1 | (opt) Deconvolution | **RedLionfish** GPU Richardson-Lucy ([PMID 39309225](https://pubmed.ncbi.nlm.nih.gov/39309225/)) | sharpens SC for skeletons. Keep raw; never before intensity readouts |
| 1b | (opt) Denoise | **Cellpose3 restoration** / Noise2Void ([PMC11903308](https://pmc.ncbi.nlm.nih.gov/articles/PMC11903308/)) | better segmentation; never before raw-intensity readouts |
| 2 | Nucleus segmentation 3D | **Cellpose-SAM** `do_3D, anisotropy=dz/dxy`; StarDist-3D alt; ilastik fallback | generalist 3D; crowded germline needs fine-tuning (Jaccard 0.47→0.78, [microPub 2023](https://micropublication.org/journals/biology/micropub.biology.001062/)) |
| 3 | Germline axis + distal→proximal order | **Gonad Linearization** (Libuda lab, [PMC8045727](https://pmc.ncbi.nlm.nih.gov/articles/PMC8045727/), GPL-3) reimpl | draw axis, project nuclei ⟂, normalize 0→1. Axis is hand-drawn upstream → see Risk 3 |
| 4 | Zone calling (TZ vs pachytene) | DAPI crescent per row **(boundary = ≥2 crescent nuclei / ≥60%, NOT 80%)**; no-DAPI synapsis-state classifier | [MacQueen & Villeneuve PMC312723](https://pmc.ncbi.nlm.nih.gov/articles/PMC312723/); no-DAPI route is novel → §3.4 |
| 5 | Per-nucleus SC tracing | crop nucleus → `skimage.filters.sato` → threshold → `skeletonize` (3D) → `skan.summarize(spacing=…)`; SNT/Imaris on hard nuclei | µm length on anisotropic stacks; report per-track **and** per-nucleus sum |
| 6 | Crossover/DSB foci | `skimage.feature.blob_log` (simple) / **big-fish**, **RS-FISH** (sub-pixel 3D) | COSA-1/MSH-5/ZHP-3 late-pachytene only; RAD-51 across zones; radius in µm |
| 7 | Granule 3D segmentation | threshold (Li) + watershed, or Cellpose; **pyclesperanto** GPU watershed | field-standard granule metrics |
| 8 | Measurement | `regionprops` (volume, intensity, centroid) + `marching_cubes(spacing=…)` → `mesh_surface_area` | reproduces Imaris Surfaces stats in citable code |
| 9 | Tidy output | long CSV + Parquet (`pyarrow`), one row per object/nucleus/zone/image | §5; every row carries voxel size + git SHA |
| 10 | Renders / QC | **napari** + napari-animation; **Blender Microscopy Nodes** hero figs; overlay montages | scriptable GPU renders as pipeline artifacts; montages double as supp figures |

---

## 3. Method grounding (publication credibility)

- **SC length:** trace central-element (SYP) or axis (HTP-3/HIM-3) filament in 3D per cropped
  pachytene nucleus, µm. WT **oocyte = 6 SCs** (5 autosomal + XX), **spermatocyte = 5** (X
  univalent) — scope by sex. Literature convention is **per-chromosome** length; emit per-track
  distribution **and** per-nucleus sum ([Cahoon eLife 2023, PMC10611432](https://pmc.ncbi.nlm.nih.gov/articles/PMC10611432/)).
  *For this project the headline is SC **fragmentation**: number of discrete tracks/fragments per
  nucleus + fragment-length distribution + total length — the heat phenotype.*
- **Transition zone:** crescent/polarized DAPI chromatin. Boundary = distal-most row with **≥2
  crescent nuclei** (count) or **≥60% crescent** (percentage) — **not 80%**. Report rows **and**
  µm (via spacing) **and** fraction-of-germline ([Phillips/Dernburg MMB 2009]; [WormBook Meiosis, NBK430708](https://www.ncbi.nlm.nih.gov/books/NBK430708/)).
- **Pachytene:** completed synapsis — thick parallel DAPI tracks flanking continuous SYP; subdivide
  into equal early/mid/late thirds of normalized span. SYP-1 asymmetric disassembly at
  pachytene→diplotene brackets but doesn't cleanly mark the proximal end.
- **No-DAPI fallback (validated in principle, 2 mandatory refinements):** loading order HTP-3 →
  HIM-3 → HTP-1/2 (axis), then SYP-1/2/3/4 (central) ([Severson Genes Dev 2009, PMC2720254](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2720254/)).
  TZ = axis loaded + SYP incomplete/polarized; pachytene = full co-tracking. **(1)** threshold on
  *axial/filamentous* HTP-3, not total signal. **(2)** wild-type only; in synapsis mutants
  "incomplete SYP" is not uniquely diagnostic. Calibrate against DAPI gonads of the same strain.
- **Crossover foci (bonus):** COSA-1/MSH-5/ZHP-3 **late pachytene only**; WT oocyte ≈6, sperm ≈5.
  Report %-nuclei-with-N-foci across thirds.
- **RAD-51/DSB foci:** canonical 6–7 zone foci/nucleus; zoning is a prerequisite for valid foci
  staging — the modules are coupled.

---

## 4. Reproducibility & orchestration

- **Env:** `pixi` (`pixi.toml` + `pixi.lock`), multi-platform (`linux-64` 5090/HPC, `osx-arm64` M4).
- **GPU (verified):** RTX 5090 = **sm_120**; needs CUDA 12.8+ wheels. **PyTorch ≥2.7 (cu128)**
  first stable with sm_120 ([issue #164342](https://github.com/pytorch/pytorch/issues/164342)).
  Pin identical torch on workstation + HPC; assert `torch.cuda.get_device_capability()==(12,0)` at
  build. Prefer PyTorch tools (Cellpose) over TensorFlow (StarDist/CARE) on Blackwell; isolate any
  TF tool in its own container.
- **Containers:** one Docker recipe (`nvidia/cuda:12.8` base + explicit cu128 torch wheel) → run on
  SLURM via **Apptainer** `.sif` (`--nv`).
- **Orchestration:** **Snakemake ≥8**, one job per `.nd2`, mirror NAS input tree → output tree,
  `--executor slurm`, resources in a profile.
- **Provenance:** Snakemake HTML report + a per-run `run_manifest.json` (git SHA, `pixi.lock` hash,
  tool versions, resolved config, per-file voxel size).
- **QC + validation:** overlay/montage PNGs (raw + mask + trace/zone + scale bar) + QC CSV that
  flags failures (0 objects, implausible length, missing channel). Validate segmentation with
  **Metrics Reloaded** (detection F1 + per-instance DSC/IoU).
- **Archival:** OME-Zarr intermediates (optional); `.nd2` stays read-only source.

## 5–8

See [README](../README.md) §schema for the tidy data model, and §8 of the design pass for the
top risks (SC tracing in dense 3D; novel no-DAPI zoning; non-automated axis; crowded-germline
segmentation fine-tuning; Blackwell CUDA; anisotropy). Repo layout is realized under `src/germquant/`.
