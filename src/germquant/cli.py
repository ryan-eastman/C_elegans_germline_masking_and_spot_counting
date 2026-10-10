"""germquant command line.

  germquant info  STACK.nd2
  germquant run   STACK.nd2 --config config/config.yaml --out results/ [--xy-stride 4 --z-range 20 40]
  germquant run   SERIES.nd2 --config config/config.yaml --out results/ --t all|7|0-24 [--resume]
  germquant register results/SERIES [--reference previous|t0] [--upsample 10] [--max-step-um 5] [--channel ROLE] [--rigid]
  germquant register results/SERIES --mode piecewise [--segments 4] [--margin-um 6] [--taper 0.25]
  germquant batch /nas/folder --config config/config.yaml --out /nas/folder_results

`info` and `run` also take TIFF / OME-TIFF, CZI, LIF and OME-Zarr files when the optional `readers`
extra is installed (uv pip install -e ".[readers]" from the repository folder); `batch` works on .nd2 folders.
`run --t` processes the named timepoints of a time-lapse one by one into results/t<NNN>/ and writes the
series record and summary next to them (docs/TIMELAPSE.md); without --t a run is the single-image run.

The nucleus training loop (docs/TIMELAPSE.md, "Segmenting without DAPI: the training loop"):

  germquant export-training RESULTS_OR_SERIES_DIR --out PAIRS [--orientations xy,xz,yz] [--slices-per-image 3]
  germquant train-nuclei --pairs PAIRS --out MODELDIR [--base cpsam] [--epochs 150] [--split-by series] [--cpu]
  germquant register-model --model MODELDIR/NAME --card MODELDIR/NAME.card.json [--name NAME]

The EXPERIMENTAL net3d RAD-51 foci detector (docs/BENCHMARK.md, "net3d"), trained on the bench cache against
the Imaris spots of the train gonads, selected on the val gonads, never on test:

  germquant train-spots --cache DIR --split-file tests/golden/rad51_split.json --out MODELDIR [--epochs N] [--patch Z,Y,X] [--seed S]

`register` estimates the drift of a finished series from the reference volumes that loop kept and writes the
register_shifts table and a report next to the series record; nothing is resampled. `--mode piecewise` estimates
one translation per axial segment of the germline instead (the register_piecewise table), for a germline that
slides inside the worm non-rigidly.

Tracking a finished series (docs/TIMELAPSE.md, "Tracking" and "Validating the linking"):
  germquant track-export   SERIES_DIR --out DIR [--crop germline|all] [--stride 2]   (a Cell-ACDC position folder)
  germquant track-link     SERIES_DIR --backend overlap|overlap_inhouse|laptrack|none [--gap 2] [--max-dist-um 3.0]
                           [--registration rigid|piecewise]
  germquant track-import   SERIES_DIR --from trackmate|imaris|cellacdc FILE
  germquant track-strips   SERIES_DIR --out DIR [--size-um 8]
  germquant track-validate SERIES_DIR --points hand.csv
"""
from __future__ import annotations

# --- cuBLAS load-order guard (must be the FIRST import) ---------------------------------------
# torch (cuBLAS 12.8) and CuPy (cuBLAS 12.9, via the nvidia-*-cu12 wheels SpotMAX/CuPy pull in)
# each ship their own cublas64_12.dll. On Windows the first one imported claims that DLL name for
# the whole process; if CuPy wins, torch's batched GEMMs fail on the RTX 5090 (Blackwell/sm_120)
# with CUBLAS_STATUS_INVALID_VALUE and Cellpose-SAM silently falls back to classical watershed.
# germquant.exe enters here, so importing torch first (before any SpotMAX/CuPy import) is the one
# place that guarantees torch's cuBLAS loads first. See pipeline.py for the same guard (defence in
# depth for `from germquant.pipeline import ...` used by scripts/notebooks).
try:
    import torch  # noqa: F401  (side effect: claim cublas64_12.dll before CuPy can)
except Exception:
    pass

import argparse
import fnmatch
import logging
import sys
from pathlib import Path

from . import provenance
from .config import find_profile, load_config
from .io import image_format, read_image_metadata


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="germquant", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    pi = sub.add_parser("info", help="print the channels + voxel size of an image (no pixels loaded)")
    pi.add_argument("nd2", help="an .nd2 (or .tif / .ome.tif / .czi / .lif / .zarr with the readers extra)")

    pr = sub.add_parser("run", help="process a single image (.nd2; other formats with the readers extra)")
    pr.add_argument("nd2", help="an .nd2 (or .tif / .ome.tif / .czi / .lif / .zarr with the readers extra)")
    _add_config_args(pr)
    pr.add_argument("--out", required=True)
    pr.add_argument("--xy-stride", type=int, default=1, help="downsample xy for a quick test")
    pr.add_argument("--z-range", type=int, nargs=2, default=None, metavar=("Z0", "Z1"))
    pr.add_argument("--t", default=None, metavar="all|N|A-B",
                    help="time-lapse: process these timepoints ('all', one index such as 7, or an inclusive "
                         "range such as 0-24) one by one into OUT/t<NNN>/, each as a fixed image under the "
                         "id <stem>_t<NNN>, and write <stem>__series.json and <stem>__series_summary.csv to OUT "
                         "(docs/TIMELAPSE.md). Without --t the file's first timepoint is the image, as before.")
    pr.add_argument("--resume", action="store_true",
                    help="with --t: skip timepoints whose completion marker matches this config, geometry and "
                         "segment stage (an interrupted series picks up where it stopped)")
    _add_stage_switches(pr)

    pb = sub.add_parser("batch", help="process every .nd2 under a folder, mirroring the tree")
    pb.add_argument("folder")
    _add_config_args(pb)
    pb.add_argument("--out", required=True)
    pb.add_argument("--xy-stride", type=int, default=1)
    pb.add_argument("--resume", action="store_true",
                    help="skip images whose results folder holds a completion marker written with the "
                         "same config and the same stages (an interrupted batch picks up where it stopped)")
    _add_stage_switches(pb)

    pc = sub.add_parser("collect", help="stack the per-image tables under a results folder into batch_<table>.csv "
                                        "and rebuild batch_summary.csv (safe to rerun any time)")
    pc.add_argument("results_root")

    prg = sub.add_parser("register", help="time-lapse drift as a shift table: phase correlation of the reference "
                                          "volumes `run --t` kept per timepoint, accumulated to the first timepoint, "
                                          "written as <stem>__register_shifts.csv plus a report (docs/TIMELAPSE.md); "
                                          "no image is resampled")
    prg.add_argument("series_dir", help="the --out folder of a `germquant run SERIES.nd2 --t ...` run")
    prg.add_argument("--reference", choices=["previous", "t0"], default="previous",
                     help="register each timepoint against the previous one and accumulate (default), or "
                          "directly against the first timepoint")
    prg.add_argument("--upsample", type=int, default=10, metavar="N",
                     help="phase correlation upsampling: estimates are quantised to 1/N of a reference-volume voxel "
                          "(default 10)")
    prg.add_argument("--max-step-um", type=float, default=5.0, metavar="UM",
                     help="cap on one step in microns; a longer step is shortened to it and flagged capped "
                          "(default 5.0; 0 = no cap)")
    prg.add_argument("--channel", default=None, metavar="ROLE",
                     help="channel role whose reference volume is registered (default: the segmentation channel's; "
                          "'brightfield' when a brightfield volume was kept)")
    prg.add_argument("--normalization", choices=["none", "phase"], default="none",
                     help="cross-power spectrum normalization: 'none' (plain cross-correlation, the default: exact "
                          "on smooth nuclei at stride 4 to 8) or 'phase' (whitened; for textured volumes)")
    prg.add_argument("--taper", type=float, default=None, metavar="FRACTION",
                     help="Tukey window fraction applied to both volumes (0.25 tapers the outer eighth of each axis): "
                          "quiets tissue the field edge cuts, a stationary feature that pulls the estimate towards zero. "
                          "Default 0 (off) for --mode rigid; for --mode piecewise the timelapse.register.piecewise_taper "
                          "key, code default 0.25 (a segment crop always cuts the gonad)")
    prg.add_argument("--rigid", action="store_true",
                     help="EXPERIMENTAL: also fit a rigid rotation (Kabsch) on nucleus centroids linked by nearest "
                          "neighbour after the translation; stored in the report only, the shift table is unchanged "
                          "(--mode rigid only)")
    prg.add_argument("--mode", choices=["rigid", "piecewise"], default="rigid",
                     help="rigid: one translation per timepoint, <stem>__register_shifts.csv (default); piecewise: one "
                          "translation per axial segment of the germline per timepoint, interpolated along the gonad "
                          "axis for every nucleus, <stem>__register_piecewise.csv (the rigid table is not touched)")
    prg.add_argument("--segments", type=int, default=None, metavar="K",
                     help="piecewise: number of segments of equal nucleus count along the axis (default: the "
                          "timelapse.register.piecewise_segments key, code default 4)")
    prg.add_argument("--margin-um", type=float, default=None, metavar="UM",
                     help="piecewise: how far a segment's crop reaches beyond its nuclei in y and x (default: the "
                          "timelapse.register.piecewise_margin_um key, code default 6.0)")
    gcfg = prg.add_mutually_exclusive_group(required=False)
    gcfg.add_argument("--config", help="config yaml for the timelapse.register.* keys (default: the series' run_manifest.json)")
    gcfg.add_argument("--profile", help="profile name instead of --config")

    pt = sub.add_parser("trace", help="draw the pachytene region on each finished image (pop-up); saves the "
                                      "polylines in whole-image microns to the traces file the staging stage reads")
    pt.add_argument("results_root")
    pt.add_argument("image_ids", nargs="*", help="only these images (default: every finished image without a trace)")
    _add_config_args(pt)
    pt.add_argument("--traces", help="traces JSON, relative to the working directory (default: staging.traces_file "
                                     "of the config, relative to the repo)")
    pt.add_argument("--redo", action="store_true", help="also show images that already have a trace")
    pt.add_argument("--stride", type=int, default=2, help="xy stride of the display (2 = half resolution)")

    ps = sub.add_parser("restage", help="recompute <image_id>__zones.csv for finished images from the traces file "
                                        "(cheap geometry; nothing else is touched)")
    ps.add_argument("results_root")
    ps.add_argument("image_ids", nargs="*")
    _add_config_args(ps)
    ps.add_argument("--traces", help="traces JSON, relative to the working directory (default: staging.traces_file "
                                     "of the config, relative to the repo)")

    pv = sub.add_parser("validate", help="compare pipeline output to hand-scored ground truth")
    pv.add_argument("--pred", help="pipeline CSV (counts/lengths mode)")
    pv.add_argument("--truth", help="ground-truth CSV")
    pv.add_argument("--key", nargs="+", default=["image_id", "nucleus_id"], help="join key columns")
    pv.add_argument("--pred-col")
    pv.add_argument("--truth-col")
    pv.add_argument("--seg-pred", help="predicted label image .tif (segmentation mode)")
    pv.add_argument("--seg-truth", help="ground-truth label image .tif")
    pv.add_argument("--iou", type=float, default=0.5)
    pv.add_argument("--out", default="validation")

    pp = sub.add_parser("prep-training", help="export DAPI z-slices from .nd2 for annotation")
    pp.add_argument("folder")
    pp.add_argument("--config", required=True)
    pp.add_argument("--out", required=True)
    pp.add_argument("--n-slices", type=int, default=3)
    pp.add_argument("--xy-stride", type=int, default=1)

    pf = sub.add_parser("finetune", help="fine-tune a germline Cellpose model (GPU)")
    pf.add_argument("labeled_dir")
    pf.add_argument("--out-model", required=True)
    pf.add_argument("--pretrained", default="cpsam")
    pf.add_argument("--epochs", type=int, default=100)
    pf.add_argument("--print-only", action="store_true", help="print the command, don't run")

    sub.add_parser("check-gpu", help="pre-flight: confirm an accelerator is visible (a CUDA GPU with compute "
                                     "capability >= 7.0, or an Apple Silicon GPU through Metal)")

    pm = sub.add_parser("get-model", help="download the trained Cellpose nucleus model (1.2 GB) from the GitHub "
                                          "release into models/models/ and verify its checksum")
    pm.add_argument("name", nargs="?", default="germline_nuclei_combined")
    pm.add_argument("--force", action="store_true", help="download even if a verified copy exists")

    pit = sub.add_parser("imaris-truth", help="read the expert truth out of an Imaris .ims project (Spots, Surfaces, "
                                              "the traced nucleus masks) into the truth files the benchmark compares "
                                              "against (docs/IMARIS_TRUTH.md)")
    pit.add_argument("ims", help="an Imaris .ims project (HDF5); read only")
    pit.add_argument("--out", required=True, help="folder for <stem>__imaris_{spots,nuclei}.csv, "
                                                  "<stem>__imaris_labels.tif and <stem>__imaris_params.json")
    pit.add_argument("--nd2", help="the .nd2 the project was made from; its voxel size defines our frame "
                                   "(default: a sibling <stem>.nd2, else Imaris's own voxel size)")
    pit.add_argument("--no-masks", action="store_true",
                     help="skip the mask channel (no label image; much faster, tables only)")

    pe = sub.add_parser("export-training",
                        help="write paired 2D training slices (image/ and labels/ TIFs plus manifest.csv) from the "
                             "nucleus labels of every finished run under a results or series folder, for correction "
                             "in napari or the Cellpose GUI (docs/TIMELAPSE.md)")
    pe.add_argument("results_dir", help="a results folder, a time-lapse series folder (its t<NNN>/ subfolders are "
                                        "the runs) or a batch tree")
    pe.add_argument("--out", required=True, help="the training folder: image/, labels/, manifest.csv, export.json")
    pe.add_argument("--orientations", default="xy,xz,yz",
                    help="comma-separated subset of xy,xz,yz; xz and yz slices are resampled along z to the xy "
                         "pixel size so every slice is isotropic (default: all three)")
    pe.add_argument("--slices-per-image", type=int, default=3, metavar="N",
                    help="slices per run and orientation, drawn at random (seeded) from slices that cut at least "
                         "--min-labels nuclei (default 3)")
    pe.add_argument("--seed", type=int, default=0, help="seed of the slice draw (default 0)")
    pe.add_argument("--zone-column", default="zone",
                    help="nuclei-table column whose most frequent value among a slice's nuclei becomes its zone "
                         "tag in the manifest (default zone, the staging stage's; empty when absent)")
    pe.add_argument("--min-labels", type=int, default=3, help="a slice must cut at least this many nuclei (default 3)")
    pe.add_argument("--max-gb", type=float, default=2.0,
                    help="stop once this many GB of slices are written (default 2; a 56 GB series must not "
                         "become a 56 GB training folder)")

    ptn = sub.add_parser("train-nuclei",
                         help="train a Cellpose nucleus model on corrected slice pairs (cellpose.train.train_seg in "
                              "process), split by image or series, never by slice; writes the model file, "
                              "training_log.json and a model card with hold-out F1 (docs/TIMELAPSE.md)")
    ptn.add_argument("--pairs", required=True, help="the folder export-training wrote (image/, labels/, manifest.csv)")
    ptn.add_argument("--base", default="cpsam",
                     help="starting weights: a registered model name, a model file, a Cellpose built-in name "
                          "(default cpsam) or 'none' for a randomly initialised network (smoke tests only)")
    ptn.add_argument("--out", required=True, metavar="MODELDIR",
                     help="output folder: MODELDIR/<name> (the model), training_log.json, <name>.card.json")
    ptn.add_argument("--name", help="model name (default: the name of MODELDIR)")
    ptn.add_argument("--epochs", type=int, default=150, help="training epochs (default 150, as the combined model)")
    ptn.add_argument("--split-by", choices=["image", "series"], default="series",
                     help="hold out whole source images, or whole series so neighbouring timepoints never straddle "
                          "the split (default series; for fixed images both are the same)")
    ptn.add_argument("--val-frac", type=float, default=0.2, help="fraction of groups held out (default 0.2, at least one)")
    ptn.add_argument("--seed", type=int, default=0, help="seed of the split (default 0)")
    ptn.add_argument("--cpu", action="store_true", help="train on the CPU even when an accelerator is available")
    ptn.add_argument("--batch-size", type=int, default=8, help="images per step (default 8)")
    ptn.add_argument("--learning-rate", type=float, default=1e-5, help="AdamW learning rate (default 1e-5)")
    ptn.add_argument("--min-train-masks", type=int, default=5,
                     help="Cellpose drops training slices with fewer nuclei than this (default 5)")

    pts = sub.add_parser("train-spots",
                         help="train the EXPERIMENTAL net3d RAD-51 foci detector (a residual 3D U-Net heatmap) on the "
                              "bench cache against the Imaris Class A spots of the train gonads with usable spot "
                              "truth; early stopping and the peak threshold on the val gonads, scored with the bench "
                              "harness; never touches a test gonad. Writes MODELDIR/model.pt, training_log.json and "
                              "model_card.json (docs/BENCHMARK.md, net3d)")
    pts.add_argument("--cache", required=True, help="the bench cache root (scripts/bench_cache.py writes it)")
    pts.add_argument("--split-file", required=True,
                     help="the pre-registered split record with the train / val / test roles "
                          "(tests/golden/rad51_split.json)")
    pts.add_argument("--out", required=True, metavar="MODELDIR",
                     help="output folder (keep it outside the repository)")
    pts.add_argument("--epochs", type=int, default=40, help="maximum epochs (default 40; early stopping may end sooner)")
    pts.add_argument("--patch", default="32,128,128",
                     help="training patch Z,Y,X in voxels, multiples of 4,8,8 (default 32,128,128)")
    pts.add_argument("--seed", type=int, default=0, help="seed of the patch draw, the initialisation and torch (default 0)")
    pts.add_argument("--manifest", default=None,
                     help="benchmark manifest naming the truth folders and the usable spot truth (default: "
                          "germquant.benchmark.detector.DEFAULT_MANIFEST)")
    pts.add_argument("--batch-size", type=int, default=8, help="patches per step (default 8)")
    pts.add_argument("--steps-per-epoch", type=int, default=200, help="steps per epoch (default 200)")
    pts.add_argument("--lr", type=float, default=1e-3, help="AdamW peak learning rate (default 1e-3)")
    pts.add_argument("--patience", type=int, default=8,
                     help="stop after this many validations without improvement (default 8)")
    pts.add_argument("--val-every", type=int, default=1, help="validate every N epochs (default 1)")
    pts.add_argument("--no-dna", action="store_true", help="RAD-51 channel only (the ablation without DAPI)")
    pts.add_argument("--device", choices=["auto", "cuda", "cpu"], default="auto",
                     help="auto = CUDA when available (default)")
    pts.add_argument("--min-free-gb", type=float, default=0.0,
                     help="on CUDA, wait until the GPU has this many GB free, polling every 60 s (default 0: no wait)")
    pts.add_argument("--wait-minutes", type=float, default=120.0, help="give up waiting for the GPU after this long")
    pts.add_argument("--prep-dir", default=None,
                     help="folder of the per-gonad loss-region maps (default MODELDIR/prep); share it between runs")
    pts.add_argument("--base-channels", type=int, default=16, help="width of the first level (default 16)")
    pts.add_argument("--levels", type=int, default=4, help="U-Net levels (default 4)")
    pts.add_argument("--thresholds", default=None,
                     help="comma-separated validation threshold grid (default 0.05 to 0.90 in steps of 0.05)")
    pts.add_argument("--tile", default="96,384,384", help="inference tile Z,Y,X for the validation (default 96,384,384)")
    pts.add_argument("--select-tolerance", type=float, default=0.01,
                     help="configurations whose val traced CCC is within this of the best count as tied and the "
                          "best F1 among them is selected (default 0.01, the protocol's meaningful change)")

    prm = sub.add_parser("register-model",
                         help="add a locally trained model to the model registry (local_models.json next to the "
                              "release assets), so a profile can name it in segmentation.nuclei.cellpose_model and "
                              "`get-model NAME` verifies it")
    prm.add_argument("--model", required=True, help="the model file `train-nuclei` wrote")
    prm.add_argument("--card", help="its model card JSON (<name>.card.json next to the model file)")
    prm.add_argument("--name", help="registry name (default: the card's name, else the file name)")

    _add_track_parsers(sub)

    args = p.parse_args(argv)
    _force_utf8_stdio()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    return {
        "info": lambda: _info(args.nd2),
        "run": lambda: _run(args),
        "batch": lambda: _batch(args),
        "collect": lambda: _collect(args),
        "register": lambda: _register(args),
        "trace": lambda: _trace(args),
        "restage": lambda: _restage(args),
        "validate": lambda: _validate(args),
        "prep-training": lambda: _prep_training(args),
        "finetune": lambda: _finetune(args),
        "check-gpu": lambda: _check_gpu(),
        "get-model": lambda: _get_model(args),
        "imaris-truth": lambda: _imaris_truth(args),
        "export-training": lambda: _export_training(args),
        "train-nuclei": lambda: _train_nuclei(args),
        "train-spots": lambda: _train_spots(args),
        "register-model": lambda: _register_model(args),

        "track-export": lambda: _track_export(args),
        "track-link": lambda: _track_link(args),
        "track-import": lambda: _track_import(args),
        "track-strips": lambda: _track_strips(args),
        "track-validate": lambda: _track_validate(args),
    }[args.cmd]()


def _add_track_parsers(sub) -> None:
    """The tracking exchange of a time-lapse series (germquant.timelapse.track / track_validate). Every command
    takes the series folder `germquant run --t` wrote; --config / --profile are optional here: without them the
    timelapse.* keys are read from the series' run_manifest.json, and every key has a code default."""
    def common(parser, stem_help="the series stem when the folder holds several (default: the one <stem>__series.json)"):
        parser.add_argument("series_dir", help="the results folder of `germquant run SERIES.nd2 --t ...` (t<NNN>/ folders inside)")
        parser.add_argument("--stem", default=None, help=stem_help)
        g = parser.add_mutually_exclusive_group(required=False)
        g.add_argument("--config", help="config yaml for the timelapse.* keys (default: the series' run_manifest.json)")
        g.add_argument("--profile", help="profile name instead of --config")

    pe = sub.add_parser("track-export", help="write the germline nucleus labels (registered, cropped, strided) and the "
                                             "central_element channel as a Cell-ACDC Position_1/Images folder + export.json")
    common(pe)
    pe.add_argument("--out", required=True, help="folder to create Position_1/Images in")
    pe.add_argument("--crop", choices=["germline", "all"], default=None,
                    help="germline: in_germline labels cropped to their union bounding box (default); all: every label, whole field")
    pe.add_argument("--stride", type=int, default=None, help="xy sub-sampling of the export (default 2)")
    pe.add_argument("--image", choices=["ref", "file", "none"], default="ref",
                    help="image stack source: the reference volumes the loop kept (ref, default), the series file at the "
                         "export stride (file; one full frame read per timepoint), or no image (none)")
    pe.add_argument("--channel-name", default="central_element", help="Cell-ACDC channel name of the image file")

    pl = sub.add_parser("track-link", help="link the nuclei across timepoints and write <stem>__track_points.csv and "
                                           "<stem>__tracks.csv next to the timepoint folders")
    common(pl)
    pl.add_argument("--backend", choices=["overlap", "overlap_inhouse", "laptrack", "none"], default="overlap",
                    help="overlap: Cell-ACDC's overlap tracker run headless (falls back to overlap_inhouse, the same rule in "
                         "NumPy, when it cannot be imported); laptrack: LAP linking of the registered centroids with gap "
                         "closing; none: the points table only (measurements per timepoint, no tracks)")
    pl.add_argument("--gap", type=int, default=None, help="laptrack: bridge a track end at t to a start at t + d for 1 <= d <= GAP "
                                                          "(default 2 = one missing timepoint; 0 = off)")
    pl.add_argument("--max-dist-um", type=float, default=None, help="laptrack: largest registered displacement linked (default 3.0)")
    pl.add_argument("--min-overlap", type=float, default=None,
                    help="overlap backends: least intersection over the previous nucleus' volume (default 0.4, Cell-ACDC's)")
    pl.add_argument("--stride", type=int, default=None, help="overlap backends: xy stride of the label grid (default 2)")
    pl.add_argument("--all-nuclei", action="store_true", help="link every label, not only in_germline nuclei")
    pl.add_argument("--registration", choices=["rigid", "piecewise"], default=None,
                    help="rigid: subtract each timepoint's shift from <stem>__register_shifts.csv; piecewise: subtract "
                         "each nucleus's displacement interpolated along the gonad axis from <stem>__register_piecewise.csv "
                         "(germquant register --mode piecewise). Default: the timelapse.track.registration key, code "
                         "default rigid")

    pi = sub.add_parser("track-import", help="turn a tracker's result (TrackMate CSV / XML, an Imaris points table, a "
                                             "Cell-ACDC tracked segm.npz) into the same two tables, link_source imported:<tool>")
    common(pi)
    pi.add_argument("--from", dest="source", required=True, choices=["trackmate", "imaris", "cellacdc"])
    pi.add_argument("file", help="the tracker's file")
    pi.add_argument("--space", choices=["raw", "export"], default="raw",
                    help="raw: coordinates in raw-frame microns (default); export: measured on the exported stack")
    pi.add_argument("--export-json", default=None, help="the export.json of the exported stack (default: found next to FILE)")
    pi.add_argument("--time-base", type=int, default=0, help="index of the first timepoint in the file (Imaris: 1)")
    pi.add_argument("--snap-um", type=float, default=0.0,
                    help="a point just outside every label takes the nearest nucleus within this distance (default 0 = off)")
    pi.add_argument("--all-nuclei", action="store_true", help="keep every label in the points table, not only in_germline nuclei")

    ps = sub.add_parser("track-strips", help="one PNG per track with a tile per timepoint around the tracked nucleus: the "
                                             "non-programmer's check that an identity did not jump")
    common(ps)
    ps.add_argument("--out", required=True)
    ps.add_argument("--size-um", type=float, default=None, help="tile width in microns (default 8)")
    ps.add_argument("--planes", type=int, default=None, help="planes max-projected around the centroid (default 3)")
    ps.add_argument("--image", choices=["ref", "file"], default="ref", help="reference volumes (default) or the series file")
    ps.add_argument("--stride", type=int, default=None, help="with --image file: xy stride of the read (default 2)")
    ps.add_argument("--max-tracks", type=int, default=None, help="only the first N tracks")

    pv = sub.add_parser("track-validate", help="score the linking against hand-followed nuclei (scripts/pick_validation_nuclei.py "
                                               "picks them; the biologist completes the worksheet in Imaris)")
    common(pv)
    pv.add_argument("--points", required=True, help="the completed worksheet: hand_track_id, t_index, z_um, y_um, x_um (raw-frame microns)")
    pv.add_argument("--snap-um", type=float, default=0.0,
                    help="a hand point just outside every label takes the nearest nucleus within this distance (default 0 = off)")
    pv.add_argument("--set", dest="subset", choices=["all", "constraint", "holdout"], default="all",
                    help="score only the hand tracks of this role: constraint to compare linking settings while the "
                         "holdout tracks stay unseen, holdout (or all, the default) to judge the chosen setting once; "
                         "a subset run writes only <stem>__track_validation_<set>.json")


def _track_cfg(args):
    """The config the timelapse.* keys are read from: --config / --profile when given, else the series' own
    run_manifest.json config (a plain dict), else None (every key falls back to its code default)."""
    import json

    if getattr(args, "config", None) or getattr(args, "profile", None):
        return load_config(_resolve_config_path(args))
    man = Path(args.series_dir) / "run_manifest.json"
    try:
        with open(man, encoding="utf-8") as fh:
            return json.load(fh).get("config") or None
    except (OSError, ValueError):
        return None


def _track_export(args) -> int:
    from .timelapse.track import export_cellacdc

    rec = export_cellacdc(args.series_dir, args.out, crop=args.crop, stride=args.stride, image=args.image,
                          channel_name=args.channel_name, stem=args.stem, cfg=_track_cfg(args))
    shape = "x".join(str(v) for v in rec["shape"])
    print(f"{rec['series_id']}: {len(rec['timepoints'])} timepoint(s) exported as (T, Z, Y, X) = {shape}, stride {rec['stride']}, "
          f"crop {rec['crop']}, shifts from {rec['shift_source']}")
    print(f"  labels -> {rec['files']['segm']}" + (f"; image -> {rec['files']['image']}" if rec["files"]["image"] else "; no image written"))
    print(f"  -> {Path(args.out) / 'export.json'}")
    return 0


def _track_link(args) -> int:
    from .timelapse.track import link_series

    from .timelapse.piecewise import PiecewiseRegistrationMissing

    try:
        res = link_series(args.series_dir, backend=args.backend, gap=args.gap, max_dist_um=args.max_dist_um,
                          min_overlap=args.min_overlap, stride=args.stride, germline_only=not args.all_nuclei,
                          stem=args.stem, cfg=_track_cfg(args), registration=args.registration)
    except PiecewiseRegistrationMissing as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    rec = res["record"]
    print(f"{rec['series_id']}: backend {rec['backend']} ({rec['link_source']}): {rec['n_objects']} nucleus-timepoints over "
          f"{len(rec['timepoints'])} timepoint(s) -> {rec['n_tracks']} track(s), {rec['n_full_length']} full length; "
          f"{rec['registration']} registration, shifts from {rec['shift_source']}")
    if rec.get("note"):
        print(f"  note: {rec['note']}")
    flipped = [t for t, s in rec.get("axis_sign_by_t", {}).items() if int(s) < 0]
    if flipped:
        print(f"  axis orientation flipped relative to t0 at timepoint(s) {', '.join(flipped)} (aligned)")
    print(f"  -> {res['paths'].get('track_points')}\n  -> {res['paths'].get('tracks')}")
    return 0


def _track_import(args) -> int:
    from .timelapse.track import import_tracks

    res = import_tracks(args.series_dir, args.source, args.file, space=args.space, time_base=args.time_base,
                        export_json=args.export_json, snap_um=args.snap_um, germline_only=not args.all_nuclei, stem=args.stem)
    rec = res["record"]
    info = rec.get("backend_info", {})
    print(f"{rec['series_id']}: imported {args.source} from {Path(args.file).name}: {rec['n_tracks']} track(s), "
          f"{rec['n_objects_unlinked']} of {rec['n_objects']} nucleus-timepoints without a track")
    print("  " + ", ".join(f"{k} {v}" for k, v in info.items()))
    print(f"  -> {res['paths'].get('track_points')}\n  -> {res['paths'].get('tracks')}")
    return 0


def _track_strips(args) -> int:
    from .timelapse.track import track_strips

    written = track_strips(args.series_dir, args.out, size_um=args.size_um, planes=args.planes, image=args.image,
                           stride=args.stride, max_tracks=args.max_tracks, stem=args.stem, cfg=_track_cfg(args))
    print(f"{len(written)} strip(s) -> {args.out}")
    return 0


def _track_validate(args) -> int:
    from .timelapse.track_validate import format_verdict, validate_tracks

    res = validate_tracks(args.series_dir, args.points, stem=args.stem, cfg=_track_cfg(args), snap_um=args.snap_um,
                          subset=args.subset)
    print(format_verdict(res))
    print(f"  -> {res['paths'].get('track_validation') or res['paths'].get('record')}")
    return 0 if res["verdict"]["pass"] else 2


def _add_config_args(parser) -> None:
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("--config", help="config yaml (config/config.yaml, or any file; may use `extends:`)")
    g.add_argument("--profile", help="name of a profile in config/profiles/<NAME>.yaml (a small overlay "
                                     "on config/config.yaml that switches the stages a study needs)")


def _resolve_config_path(args) -> Path:
    if getattr(args, "profile", None):
        try:
            return find_profile(args.profile)       # the same lookup germquant.api.run(profile=...) uses
        except FileNotFoundError as e:
            raise SystemExit(str(e)) from None
    return Path(args.config)


def _add_stage_switches(parser) -> None:
    """One generated ``--no-<stage>`` per switchable stage in `germquant.stages.STAGES` (so
    ``--no-spots`` and ``--no-coloc`` keep working and new optional stages get a switch for free)."""
    from .stages import cli_switches

    for st in cli_switches():
        parser.add_argument(f"--no-{st.name}", action="store_true",
                            help=st.legacy_cli_help or f"skip the {st.name} stage ({st.help})")


def _apply_switches(cfg, args) -> None:
    from .stages import apply_cli_switches

    for msg in apply_cli_switches(cfg, args):
        print(msg)


def _force_utf8_stdio() -> None:
    """Windows consoles default to cp1252; our status glyphs (✓ ⚠ ≥ µ) would raise
    UnicodeEncodeError *after* the work is done, reporting a success as a crash. Reconfigure
    stdout/stderr to UTF-8 (replacing anything truly unencodable) so output never aborts a run.
    """
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")  # py3.7+ TextIOWrapper
        except (AttributeError, ValueError):  # already wrapped / not reconfigurable
            pass


def _info(nd2: str) -> int:
    m = read_image_metadata(nd2)       # read_nd2_metadata for an .nd2, the bioio header for other formats
    print(f"file       : {nd2}")
    print(f"sizes      : {m['sizes']}")
    print(f"dtype      : {m['dtype']}   2D={m['is_2d']}")
    dz, dy, dx = m["spacing"]
    print(f"voxel µm   : dz={dz:.4f} dy={dy:.4f} dx={dx:.4f}  (anisotropy z/xy={dz/dy:.2f})")
    print(f"channels   : {m['channel_names']}")
    if int(m.get("n_timepoints", 1) or 1) > 1:       # a time series: say so (fixed images print nothing here)
        period = m.get("period_s")
        every = f" every {period:g} s" if period else ""
        print(f"timepoints : {m['n_timepoints']}{every}  (run one with t=<index>; docs/TIMELAPSE.md)")
    if image_format(nd2) == "nd2":
        _print_acquisition(nd2)
    return 0


def _print_acquisition(nd2: str) -> None:
    """The per-channel acquisition table from the .nd2 description (exposure, laser power, emission filter,
    camera mode): the record the acquisition stage writes and the one acquisition.expected is checked
    against (docs/QC.md). Other formats carry no such description, so nothing is printed for them."""
    from .io.acquisition import format_table, parse_description, read_description

    try:
        per_name = parse_description(read_description(nd2))
    except Exception as e:  # noqa: BLE001 - the header fields above are already printed; say why this part is not
        print(f"acquisition: description unreadable ({e})")
        return
    if not per_name:
        print("acquisition: the description lists no per-channel planes")
        return
    print("acquisition (from the .nd2 description; the full record goes to image_summary.acquisition_json):")
    print(format_table(per_name))


def _run(args) -> int:
    from .api import run

    cfg = load_config(_resolve_config_path(args))
    _apply_switches(cfg, args)
    if getattr(args, "t", None) is not None:
        return _run_series(args, cfg)
    if getattr(args, "resume", False):
        raise SystemExit("--resume applies to a time series: give --t as well (for example --t all)")
    # germquant.api.run is the one path to a results folder: a GUI or a notebook calling it with the same
    # settings writes the same files (process_image writes run_manifest.json itself)
    res = run(args.nd2, config=cfg, out=args.out, xy_stride=args.xy_stride, z_range=args.z_range)
    mark = "✓" if res["qc_pass"] else "⚠"
    print(f"\n{mark} {res['image_id']}: {res['n_nuclei']} nuclei, qc_pass={res['qc_pass']}")
    if res["qc_flags"]:
        print("  flags:", "; ".join(res["qc_flags"]))
    print(f"  -> {res['out_dir']}")
    # QC status is recorded in the output tables; a flagged image is NOT a process
    # failure (so unattended Snakemake batches complete). Only exceptions are errors.
    return 0


def _run_series(args, cfg) -> int:
    """``germquant run SERIES.nd2 --t ...``: the per-timepoint loop (timelapse.series.process_series). A
    timepoint that raises is recorded in the series record and the loop goes on; the exit code is 0 when at
    least one requested timepoint finished and 1 when none did."""
    from .timelapse import process_series

    res = process_series(args.nd2, cfg, args.out, timepoints=args.t, xy_stride=args.xy_stride,
                         z_range=args.z_range, resume=bool(getattr(args, "resume", False)))
    rec = res["record"]
    n_req, n_done, n_failed = len(res["requested"]), len(res["done"]), len(res["failed"])
    n_pass = sum(1 for e in rec["timepoints"] if e.get("status") == "done" and e.get("qc_pass"))
    n_resumed = sum(1 for e in rec["timepoints"] if e.get("resumed"))
    print(f"\n{res['series_id']}: {n_req} timepoint(s) requested of {res['n_timepoints']}, {n_done} done"
          + (f" ({n_resumed} already on disk)" if n_resumed else "") + f", {n_failed} failed, {n_pass} passed QC")
    for f in res["failed"]:
        print(f"  t{f['t_index']:03d} FAILED: {'; '.join(f['flags'])}")
    print(f"  -> {res['out_dir']}  ({Path(res['record_path']).name}, {Path(res['summary_path']).name})")
    return 0 if n_done else 1


def _batch(args) -> int:
    import pandas as pd

    from . import batch as B
    from .pipeline import process_image

    cfg = load_config(_resolve_config_path(args))
    _apply_switches(cfg, args)
    root = Path(args.folder)
    out_root = Path(args.out)
    glob = cfg.get("io.input_glob", "**/*.nd2")
    excludes = cfg.get("io.exclude_patterns", [])
    exclusions = B.load_exclusions(cfg.get("qc.exclusions_file"), cfg.base_dir)

    files, dropped = B.discover_files(root, glob, excludes, exclusions)
    for f, why in dropped:
        if why.startswith("exclusions_file"):
            print(f"excluded by {why}: {f.name}")
    if not files:
        print(f"No .nd2 files matched {glob} under {root}", file=sys.stderr)
        return 1

    # image-independent provenance, computed once per batch: the segment sub-hash lets --resume notice a
    # retrained model or a Cellpose upgrade (per-image geometry lives in <id>__stages.json / __done.json)
    from .pipeline import resolve_model_path
    from .stages import stage_hashes as stage_hashes_fn

    model_sha = (provenance.file_sha256(resolve_model_path(cfg))
                 if str(cfg.get("segmentation.nuclei.method", "auto")) != "classical" else "unused")
    segment_hash = stage_hashes_fn(cfg, {}, provenance.tool_versions(), model_sha)["segment"]
    prov = provenance.write_manifest(out_root, config_hash=cfg.hash, config=cfg.as_dict(),
                                     extra={"n_files": len(files), "input_root": str(root),
                                            "excluded": [str(f) for f, _ in dropped],
                                            "xy_stride": int(args.xy_stride), "model_sha256": model_sha,
                                            "enabled_stages": B.enabled_stage_names(cfg)})
    print(f"Processing {len(files)} files -> {out_root}")
    summaries = []
    resume = bool(getattr(args, "resume", False))
    for i, f in enumerate(files, 1):
        rel = f.relative_to(root).parent
        out_dir = out_root / rel / f.stem
        print(f"[{i}/{len(files)}] {f.name}")
        if resume:
            done, why = B.is_done(out_dir, f.stem, cfg, xy_stride=args.xy_stride, segment_hash=segment_hash)
            if done:
                print("    already done (completion marker matches); skipping")
                continue
            if why != "no completion marker":
                print(f"    reprocessing: {why}")
        try:
            res = process_image(f, cfg, out_dir, xy_stride=args.xy_stride, prov=prov)
            isum = res["tables"]["image_summary"]
            n_germ = int(isum["n_germline_nuclei"].iloc[0]) if "n_germline_nuclei" in isum else 0
            summaries.append({"image_id": res["image_id"], "n_nuclei": res["n_nuclei"],
                              "n_germline": n_germ, "qc_pass": res["qc_pass"],
                              "qc_flags": ";".join(res["qc_flags"]), "out_dir": res["out_dir"]})
        except Exception as e:  # noqa: BLE001
            logging.exception("FAILED %s", f.name)
            summaries.append({"image_id": f.stem, "n_nuclei": 0, "n_germline": 0, "qc_pass": False,
                              "qc_flags": f"EXCEPTION:{e}", "out_dir": str(out_dir)})

    # Framing QC: flag images whose germline count is a strong outlier vs the batch median (a robust,
    # threshold-free proxy for "two gonad arms / extra tissue / fuller distal capture in frame", worth
    # eyeballing the montage; the axis/position readout for such gonads is unreliable). Spot COUNTS are
    # unaffected, so this is advisory, not a failure. The same rule serves `collect`.
    factor = float(cfg.get("qc.germline_outlier_factor", 1.8))
    from .fsutil import long_path
    if resume:
        # a resumed batch: every completed folder of THIS study on disk, plus this pass's failures
        failed = [s for s in summaries if not s["qc_pass"] and str(s["qc_flags"]).startswith("EXCEPTION")]
        counts = B.collect(out_root, image_ids={f.stem for f in files}, extra_rows=failed, outlier_factor=factor)
        print(f"\nDone. Resumed batch collected: {counts}. Summary -> {out_root / 'batch_summary.csv'}")
        return 0
    summaries = B.apply_framing_qc(summaries, factor)
    pd.DataFrame(summaries, columns=B.SUMMARY_COLS).to_csv(long_path(out_root / "batch_summary.csv"), index=False)
    n_pass = sum(s["qc_pass"] for s in summaries)
    print(f"\nDone. {n_pass}/{len(files)} passed QC. Summary -> {out_root / 'batch_summary.csv'}")
    return 0


def _traces_path(cfg, args) -> Path:
    tf = getattr(args, "traces", None)
    if tf:
        return Path(tf).resolve()                  # a command-line path is relative to the working directory
    tf = cfg.get("staging.traces_file", "staging/pachytene_traces.json")
    return Path(tf) if Path(tf).is_absolute() else cfg.base_dir / tf   # the config value: relative to the repo


def _trace(args) -> int:
    from .staging import load_traces
    from .staging.tracer import results_images, run_gui

    cfg = load_config(_resolve_config_path(args))
    tf = _traces_path(cfg, args)
    items = results_images(args.results_root)
    have = load_traces(tf)
    if args.image_ids:
        items = [it for it in items if it[0] in set(args.image_ids)]
    elif not args.redo:
        items = [it for it in items if have.get(it[0], {}).get("status") not in ("traced", "skipped")]
    if not items:
        print("nothing to trace (use --redo to revisit traced images)")
        return 0
    print(f"{len(items)} image(s) to trace -> {tf}")
    run_gui(items, cfg, tf, stride=int(args.stride),
            off_axis_um=float(cfg.get("staging.off_axis_um") or 20.0))
    return 0


def _restage(args) -> int:
    from .staging.tracer import restage

    cfg = load_config(_resolve_config_path(args))
    done = restage(args.results_root, _traces_path(cfg, args), cfg, args.image_ids or None)
    print(f"restaged {len(done)} image(s)")
    return 0


def _register(args) -> int:
    """``germquant register SERIES_DIR``: the drift of a time-lapse as a shift table (timelapse.register,
    docs/TIMELAPSE.md "Registration"). Exit code 1, with the reason, when the folder is not a series or holds no
    reference volumes; a timepoint that cannot be registered is a row with an error, not a failure. ``--mode
    piecewise`` registers each axial segment of the germline (timelapse.piecewise) instead."""
    from .timelapse.register import DEFAULTS as RIGID_DEFAULTS
    from .timelapse.register import register_series

    if getattr(args, "mode", "rigid") == "piecewise":
        return _register_piecewise(args)
    for opt, given in (("--segments", args.segments), ("--margin-um", args.margin_um)):
        if given is not None:
            print(f"ERROR: {opt} applies to --mode piecewise only", file=sys.stderr)
            return 1
    try:
        res = register_series(args.series_dir, reference=args.reference, upsample=args.upsample,
                              max_step_um=args.max_step_um, channel=args.channel, rigid=args.rigid,
                              normalization=args.normalization,
                              taper=RIGID_DEFAULTS["taper"] if args.taper is None else args.taper)
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    rep, par = res["report"], res["report"]["parameters"]
    print(f"\n{res['series_id']}: {rep['n_timepoints']} timepoint(s), {rep['n_registered']} step(s) registered on "
          f"{par['channel_role']} ({rep['method']}, stride {par['stride']}, upsample {par['upsample_factor']}, "
          f"reference {par['reference']})")
    if rep.get("cumulative_last_um"):
        dz, dy, dx = rep["cumulative_last_um"]
        cap = f"{par['max_step_um']:g} um" if par.get("max_step_um") else "none"
        print(f"  cumulative shift at t{rep['last_t_index']:03d}: dz={dz:+.2f} dy={dy:+.2f} dx={dx:+.2f} um; "
              f"largest step {rep['max_step_um']:.2f} um; {rep['n_capped']} capped (cap {cap})")
    q = rep["quality"]
    if q.get("median_before_um") is not None and q.get("median_after_um") is not None:
        print(f"  quality: median nearest-neighbour displacement {q['median_before_um']:.2f} -> "
              f"{q['median_after_um']:.2f} um (before -> after the step, over {q['n_steps']} step(s))")
    if par.get("cross_check_role"):
        med = rep.get("median_disagreement_um")
        med_txt = f"{med:.2f} um" if med is not None else "n/a"
        print(f"  {par['cross_check_role']} cross-check: median disagreement {med_txt}, {rep['n_disagreement']} step(s) "
              f"above {par['disagreement_um']:g} um")
    else:
        print("  no brightfield reference volume: cross-check skipped")
    if rep.get("rigid") is not None:
        angles = [r["rotation_deg"] for r in rep["rigid"] if r.get("rotation_deg") is not None]
        print(f"  rigid fit (experimental, report only): {len(angles)} step(s) fitted, largest rotation "
              f"{max(angles) if angles else 0.0:.2f} deg")
    for k, why in rep.get("errors", {}).items():
        print(f"  t{int(k):03d}: {why}")
    print(f"  -> {res['table_path']}, {res['report_path']}")
    return 0


def _register_piecewise(args) -> int:
    """``germquant register SERIES_DIR --mode piecewise``: one translation per axial segment per timepoint
    (timelapse.piecewise, docs/TIMELAPSE.md "Piecewise registration")."""
    from .timelapse.piecewise import register_series_piecewise

    if args.reference != "previous":
        print("ERROR: --mode piecewise registers every timepoint against the previous one (--reference previous)",
              file=sys.stderr)
        return 1
    if args.rigid:
        print("ERROR: --rigid (the experimental Kabsch fit) applies to --mode rigid only", file=sys.stderr)
        return 1
    try:
        cfg = _track_cfg(args)
        res = register_series_piecewise(args.series_dir, segments=args.segments, margin_um=args.margin_um,
                                        taper=args.taper, upsample=args.upsample, max_step_um=args.max_step_um,
                                        channel=args.channel, normalization=args.normalization, cfg=cfg)
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    rep, par = res["report"], res["report"]["parameters"]
    print(f"\n{res['series_id']}: piecewise registration, {par['segments']} segment(s) along the axis "
          f"({par['axis_source']}), {rep['n_timepoints']} timepoint(s), {rep['n_registered']} step(s) registered on "
          f"{par['channel_role']} ({rep['method']}, stride {par['stride']}, upsample {par['upsample_factor']}, "
          f"margin {par['margin_um']:g} um, taper {par['taper']:g})")
    for s, q in zip(rep["segments"], rep["quality"]["by_segment"]):
        cum = rep["cumulative_last_um"][s["segment"]]
        before, after = q.get("median_before_um"), q.get("median_after_um")
        qtxt = (f"median nearest-neighbour displacement {before:.2f} -> {after:.2f} um" if before is not None and
                after is not None else "no quality numbers")
        dtxt = (f"; {par['cross_check_role']} disagreement {q['median_disagreement_um']:.2f} um"
                if q.get("median_disagreement_um") is not None else "")
        print(f"  segment {s['segment']} (axis {s['axis_start_um']:.1f} to {s['axis_end_um']:.1f} um, {s['n_nuclei_t0']} "
              f"nuclei): cumulative at t{rep['last_t_index']:03d} dz={cum[0]:+.2f} dy={cum[1]:+.2f} dx={cum[2]:+.2f} um; "
              f"{qtxt}{dtxt}")
    print(f"  {len(rep['capped'])} segment step(s) capped, {len(rep['brightfield_disagreement'])} above the "
          f"{par['disagreement_um']:g} um disagreement threshold, {len(rep['filled'])} filled from neighbours")
    for k, why in rep.get("errors", {}).items():
        print(f"  t{int(k):03d}: {why}")
    print(f"  -> {res['table_path']}, {res['report_path']}")
    return 0


def _collect(args) -> int:
    from . import batch as B

    counts = B.collect(args.results_root)
    for k, v in counts.items():
        print(f"  {k}: {v} rows")
    print(f"  -> {Path(args.results_root) / 'batch_<table>.csv'} and batch_summary.csv")
    return 0


def _validate(args) -> int:
    import json

    import pandas as pd

    from .validate import bland_altman_plot, compare_table, segmentation_metrics

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if args.seg_pred and args.seg_truth:
        import tifffile

        pred = tifffile.imread(args.seg_pred)
        gt = tifffile.imread(args.seg_truth)
        m = segmentation_metrics(pred, gt, iou_threshold=args.iou)
        (out / "segmentation_metrics.json").write_text(json.dumps(m, indent=2))
        print(f"Segmentation @ IoU≥{args.iou}: F1={m['f1']:.3f}  precision={m['precision']:.3f}  "
              f"recall={m['recall']:.3f}  mean_IoU={m['mean_iou']:.3f}  (TP={m['tp']} FP={m['fp']} FN={m['fn']})")
        return 0

    if not (args.pred and args.truth and args.pred_col and args.truth_col):
        print("counts mode needs --pred --truth --pred-col --truth-col (and --key)", file=sys.stderr)
        return 1
    paired, stats = compare_table(
        pd.read_csv(args.pred), pd.read_csv(args.truth), args.key, args.pred_col, args.truth_col
    )
    paired.to_csv(out / "paired.csv", index=False)
    (out / "agreement.json").write_text(json.dumps(stats, indent=2))
    if stats.get("n", 0) >= 2:
        bland_altman_plot(paired["pred"], paired["truth"], out / "bland_altman.png",
                          title=f"{args.pred_col} vs {args.truth_col}")
    print(f"n={stats.get('n')}  CCC={stats.get('ccc', float('nan')):.3f}  "
          f"Pearson={stats.get('pearson_r', float('nan')):.3f}  bias={stats.get('bias_mean_diff', float('nan')):.2f}  "
          f"MAE={stats.get('mae', float('nan')):.2f}")
    print(f"  -> {out}")
    return 0


def _prep_training(args) -> int:
    from .segment.finetune import prep_training_data

    cfg = load_config(args.config)
    root = Path(args.folder)
    excludes = cfg.get("io.exclude_patterns", [])
    files = [f for f in sorted(root.rglob("*.nd2"))
             if not any(fnmatch.fnmatch(f.name.lower(), pat.lower()) for pat in excludes)]
    if not files:
        print(f"No .nd2 under {root}", file=sys.stderr)
        return 1
    written = prep_training_data(files, cfg.channel_map, args.out,
                                 n_slices=args.n_slices, xy_stride=args.xy_stride)
    print(f"Wrote {len(written)} DAPI slices to {args.out}. Annotate them in the Cellpose GUI "
          f"(see docs/ANNOTATION.md), then `germquant finetune`.")
    return 0


def _finetune(args) -> int:
    from .segment.finetune import finetune_cellpose

    finetune_cellpose(args.labeled_dir, args.out_model, pretrained=args.pretrained,
                      n_epochs=args.epochs, run=not args.print_only)
    return 0


def _get_model(args) -> int:
    from .models_registry import main as get_model_main

    return get_model_main([args.name] + (["--force"] if args.force else []))


def _export_training(args) -> int:
    from .segment.training import export_training_pairs

    orientations = tuple(o.strip() for o in str(args.orientations).split(",") if o.strip())
    manifest = export_training_pairs(args.results_dir, args.out, orientations=orientations,
                                     slices_per_image=args.slices_per_image, seed=args.seed,
                                     zone_column=args.zone_column or None, min_labels=args.min_labels,
                                     max_bytes=float(args.max_gb) * 1e9)
    n_img = manifest["image_id"].nunique() if len(manifest) else 0
    per = manifest.groupby("orientation").size().to_dict() if len(manifest) else {}
    print(f"{len(manifest)} slice pair(s) from {n_img} run(s) -> {args.out}  "
          f"({', '.join(f'{k}: {v}' for k, v in sorted(per.items()))})")
    print(f"  correct labels/<name>.tif in napari (napari {args.out}/image/X.tif {args.out}/labels/X.tif) or the "
          f"Cellpose GUI, keep the file names, then `germquant train-nuclei --pairs {args.out} --out MODELDIR`")
    return 0


def _train_nuclei(args) -> int:
    from .segment.training import train_nuclei

    res = train_nuclei(args.pairs, args.out, base=args.base, name=args.name, epochs=args.epochs,
                       split_by=args.split_by, val_frac=args.val_frac, seed=args.seed, cpu=args.cpu,
                       learning_rate=args.learning_rate, batch_size=args.batch_size,
                       min_train_masks=args.min_train_masks)
    sp = res["split"]
    print(f"\n{res['name']}: trained on {sp['n_train']} slice(s) from {len(sp['train_groups'])} {sp['by']} group(s), "
          f"{sp['n_val']} held out from {len(sp['val_groups'])}, device {res['device']}")
    overall = res["metrics"].get("overall")
    if overall:
        print(f"  hold-out at IoU 0.5: F1 {overall['f1']:.3f}  precision {overall['precision']:.3f}  "
              f"recall {overall['recall']:.3f}  miss rate {overall['miss_rate']:.3f}  "
              f"(TP {overall['tp']} FP {overall['fp']} FN {overall['fn']})")
        for o, m in res["metrics"].get("per_orientation", {}).items():
            print(f"    {o}: F1 {m['f1']:.3f} on {m['n_pairs']} slice(s)")
    elif res["metrics"].get("error"):
        print(f"  hold-out metrics missing: {res['metrics']['error']}")
    print(f"  model -> {res['model_path']} (sha256 {res['model_sha256'][:12]}...)")
    print(f"  card  -> {res['card_path']}")
    print(f"  register it:  germquant register-model --model \"{res['model_path']}\" --card \"{res['card_path']}\"")
    return 0


def _train_spots(args) -> int:
    from .spots.net3d.train import DEFAULT_THRESHOLDS, parse_patch, train_spots

    thresholds = (DEFAULT_THRESHOLDS if not args.thresholds else
                  tuple(float(t) for t in str(args.thresholds).split(",") if t.strip()))
    res = train_spots(args.cache, args.split_file, args.out, manifest=args.manifest, epochs=args.epochs,
                      patch=parse_patch(args.patch), seed=args.seed, batch_size=args.batch_size,
                      steps_per_epoch=args.steps_per_epoch, lr=args.lr, patience=args.patience,
                      val_every=args.val_every, use_dna=not args.no_dna, device=args.device,
                      base_channels=args.base_channels, levels=args.levels, thresholds=thresholds,
                      tile_zyx=parse_patch(args.tile), prep_dir=args.prep_dir, min_free_gb=args.min_free_gb,
                      wait_minutes=args.wait_minutes, select_tolerance=args.select_tolerance)
    g = res["gonads"]
    print(f"\nnet3d: trained on {[r['tag'] for r in g['train']]} ({res['device']}), validated on "
          f"{[r['tag'] for r in g['val']]}; left out {[r['tag'] for r in g['excluded']]}; test gonads never opened")
    vb = res.get("val_best") or {}
    print(f"  selected epoch {res['best_epoch']} of {res['epochs_run']} ({res['stop_reason']}), threshold "
          f"{res['best_threshold']}: val traced CCC {vb.get('ccc_traced')}, germline CCC {vb.get('ccc_all')}, "
          f"F1 {vb.get('f1')}")
    print(f"  model -> {res['model_path']} (sha256 {res['model_sha256'][:12]}...)")
    print(f"  card  -> {res['card_path']}")
    print(f"  score it: python scripts/bench_detector.py --backend net3d --params "
          f"'{{\"net3d\": {{\"model\": \"{Path(res['model_path']).as_posix()}\"}}}}' --split val --out DIR")
    return 0


def _register_model(args) -> int:
    from .models_registry import register_main

    argv = ["--model", args.model]
    if args.card:
        argv += ["--card", args.card]
    if args.name:
        argv += ["--name", args.name]
    return register_main(argv)


def _imaris_truth(args) -> int:
    from .validate.imaris_ims import export_imaris_truth

    p = export_imaris_truth(args.ims, args.out, nd2=args.nd2, masks=not args.no_masks)
    m = p["mask"]
    line = f"{p['stem']}: {p['n_spots']} spots, {p['n_surfaces']} surfaces"
    if m.get("labels_file"):
        line += (f", {m['n_surfaces_with_mask']}/{p['n_surfaces']} surfaces with a mask from channel "
                 f"{m['channel']} ('{m['channel_name']}'), labels {tuple(m['labels_shape_zyx'])}")
    else:
        line += f", no label image ({m.get('reason_no_labels')})"
    print(line)
    if m.get("labels_file"):
        print(f"  centroid check: {m['n_centroid_inside_own_label']} inside own label, {m['n_centroid_in_cavity']} on a "
              f"cavity (0), {m['n_centroid_in_other_label']} inside another label; "
              f"{m['n_orphan_values']} unclaimed mask value(s) dropped")
    for w in p["warnings"]:
        print(f"  warning: {w}")
    print(f"  -> {args.out}")
    return 0


def _check_gpu() -> int:
    """Pre-flight accelerator check for a real run. Passes on any CUDA GPU adequate for Cellpose-SAM
    (compute capability >= 7.0: the workstation RTX 5090, the Alpine A100 / L40) and on an Apple
    Silicon GPU through Metal (MPS). Returns nonzero when torch is missing, only the CPU is available,
    or the CUDA GPU is too old. Honours $GERMQUANT_DEVICE like the pipeline does."""
    from . import device as D

    try:
        import torch
    except Exception as e:  # noqa: BLE001
        print(f"torch not importable ({e}); install germquant[gpu].", file=sys.stderr)
        return 1
    chosen = D.select_device()
    print(D.describe())
    if chosen == "cuda":
        cap = torch.cuda.get_device_capability()
        if cap < (7, 0):
            print(f"GPU compute capability {tuple(cap)} < 7.0 is too old for Cellpose-SAM.", file=sys.stderr)
            return 1
        known = {(12, 0): "RTX 5090 (Blackwell)", (9, 0): "H100 (Hopper)",
                 (8, 9): "L40/L40S (Ada)", (8, 0): "A100 (Ampere)"}
        print(f"OK, {known.get(tuple(cap), 'CUDA GPU')}: usable for Cellpose-SAM + SpotMAX.")
        return 0
    if chosen == "mps":
        print("OK, Apple Silicon GPU via Metal: Cellpose-SAM runs on it in float32 (operators Metal lacks "
              "fall back to the CPU); SpotMAX runs on the CPU. Expect a run to take longer than on the RTX 5090.")
        return 0
    print("No accelerator: neither CUDA nor Metal (MPS) is available to torch. Cellpose will run on the CPU, "
          "which is very slow; check the torch install (see docs/MAC_METAL.md or the README).", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
