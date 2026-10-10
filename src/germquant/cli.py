"""germquant command line.

  germquant info  STACK.nd2
  germquant run   STACK.nd2 --config config/config.yaml --out results/ [--xy-stride 4 --z-range 20 40]
  germquant run   SERIES.nd2 --config config/config.yaml --out results/ --t all|7|0-24 [--resume]
  germquant batch /nas/folder --config config/config.yaml --out /nas/folder_results

`info` and `run` also take TIFF / OME-TIFF, CZI, LIF and OME-Zarr files when the optional `readers`
extra is installed (uv pip install -e ".[readers]" from the repository folder); `batch` works on .nd2 folders.
`run --t` processes the named timepoints of a time-lapse one by one into results/t<NNN>/ and writes the
series record and summary next to them (docs/TIMELAPSE.md); without --t a run is the single-image run.
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

    args = p.parse_args(argv)
    _force_utf8_stdio()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    return {
        "info": lambda: _info(args.nd2),
        "run": lambda: _run(args),
        "batch": lambda: _batch(args),
        "collect": lambda: _collect(args),
        "trace": lambda: _trace(args),
        "restage": lambda: _restage(args),
        "validate": lambda: _validate(args),
        "prep-training": lambda: _prep_training(args),
        "finetune": lambda: _finetune(args),
        "check-gpu": lambda: _check_gpu(),
        "get-model": lambda: _get_model(args),
        "imaris-truth": lambda: _imaris_truth(args),
    }[args.cmd]()


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
