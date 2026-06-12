"""germquant command line.

  germquant info  STACK.nd2
  germquant run   STACK.nd2 --config config/config.yaml --out results/ [--xy-stride 4 --z-range 20 40]
  germquant batch /nas/folder --config config/config.yaml --out /nas/folder_results
"""
from __future__ import annotations

import argparse
import fnmatch
import logging
import sys
from pathlib import Path

from . import provenance
from .config import load_config
from .io import read_nd2_metadata


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="germquant", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    pi = sub.add_parser("info", help="print .nd2 channels + voxel size (no pixels loaded)")
    pi.add_argument("nd2")

    pr = sub.add_parser("run", help="process a single .nd2")
    pr.add_argument("nd2")
    pr.add_argument("--config", required=True)
    pr.add_argument("--out", required=True)
    pr.add_argument("--xy-stride", type=int, default=1, help="downsample xy for a quick test")
    pr.add_argument("--z-range", type=int, nargs=2, default=None, metavar=("Z0", "Z1"))

    pb = sub.add_parser("batch", help="process every .nd2 under a folder, mirroring the tree")
    pb.add_argument("folder")
    pb.add_argument("--config", required=True)
    pb.add_argument("--out", required=True)
    pb.add_argument("--xy-stride", type=int, default=1)

    args = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")

    if args.cmd == "info":
        return _info(args.nd2)
    if args.cmd == "run":
        return _run(args)
    if args.cmd == "batch":
        return _batch(args)
    return 1


def _info(nd2: str) -> int:
    m = read_nd2_metadata(nd2)
    print(f"file       : {nd2}")
    print(f"sizes      : {m['sizes']}")
    print(f"dtype      : {m['dtype']}   2D={m['is_2d']}")
    dz, dy, dx = m["spacing"]
    print(f"voxel µm   : dz={dz:.4f} dy={dy:.4f} dx={dx:.4f}  (anisotropy z/xy={dz/dy:.2f})")
    print(f"channels   : {m['channel_names']}")
    return 0


def _run(args) -> int:
    from .pipeline import process_image

    cfg = load_config(args.config)
    out = Path(args.out)
    prov = provenance.write_manifest(out, config_hash=cfg.hash, config=cfg.as_dict())
    z_range = tuple(args.z_range) if args.z_range else None
    res = process_image(args.nd2, cfg, out, xy_stride=args.xy_stride, z_range=z_range, prov=prov)
    mark = "✓" if res["qc_pass"] else "⚠"
    print(f"\n{mark} {res['image_id']}: {res['n_nuclei']} nuclei, qc_pass={res['qc_pass']}")
    if res["qc_flags"]:
        print("  flags:", "; ".join(res["qc_flags"]))
    print(f"  -> {res['out_dir']}")
    # QC status is recorded in the output tables; a flagged image is NOT a process
    # failure (so unattended Snakemake batches complete). Only exceptions are errors.
    return 0


def _batch(args) -> int:
    import pandas as pd

    from .pipeline import process_image

    cfg = load_config(args.config)
    root = Path(args.folder)
    out_root = Path(args.out)
    glob = cfg.get("io.input_glob", "**/*.nd2")
    excludes = cfg.get("io.exclude_patterns", [])

    files = [
        f for f in sorted(root.glob(glob))
        if not any(fnmatch.fnmatch(f.name.lower(), pat.lower()) for pat in excludes)
    ]
    if not files:
        print(f"No .nd2 files matched {glob} under {root}", file=sys.stderr)
        return 1

    prov = provenance.write_manifest(out_root, config_hash=cfg.hash, config=cfg.as_dict(),
                                     extra={"n_files": len(files), "input_root": str(root)})
    print(f"Processing {len(files)} files -> {out_root}")
    summaries = []
    for i, f in enumerate(files, 1):
        rel = f.relative_to(root).parent
        out_dir = out_root / rel / f.stem
        print(f"[{i}/{len(files)}] {f.name}")
        try:
            res = process_image(f, cfg, out_dir, xy_stride=args.xy_stride, prov=prov)
            summaries.append({"image_id": res["image_id"], "n_nuclei": res["n_nuclei"],
                              "qc_pass": res["qc_pass"], "qc_flags": ";".join(res["qc_flags"]),
                              "out_dir": res["out_dir"]})
        except Exception as e:  # noqa: BLE001
            logging.exception("FAILED %s", f.name)
            summaries.append({"image_id": f.stem, "n_nuclei": 0, "qc_pass": False,
                              "qc_flags": f"EXCEPTION:{e}", "out_dir": str(out_dir)})

    pd.DataFrame(summaries).to_csv(out_root / "batch_summary.csv", index=False)
    n_pass = sum(s["qc_pass"] for s in summaries)
    print(f"\nDone. {n_pass}/{len(files)} passed QC. Summary -> {out_root / 'batch_summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
