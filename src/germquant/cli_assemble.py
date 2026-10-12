"""`germquant assemble`: stack one TIFF per channel into one OME-TIFF (germquant.assemble; docs/PACKS.md).

  germquant assemble --channels A.tif B.tif --names Bright Cy3 --spacing DZ DY DX --out X.ome.tif
                     [--axes ZYX] [--crop-z Z0 Z1] [--crop-y Y0 Y1] [--crop-x X0 X1] [--force] [--json]

For labs (and tools such as Cell-ACDC) that keep each channel in its own file: germquant reads one image per
run and binds channels to roles by name, so the assembled file carries the channel names and the voxel size.
Nothing heavy is imported (numpy and tifffile).
"""
from __future__ import annotations

import json
import sys


def add_parsers(sub) -> None:
    """The `germquant assemble` command (cli.main calls this once)."""
    pa = sub.add_parser("assemble", help="stack one TIFF per channel (or a multi-channel TIFF with a declared page "
                                         "order) into one OME-TIFF with channel names and voxel size")
    pa.add_argument("--channels", nargs="+", required=True, metavar="TIFF",
                    help="the input files, one per channel (or files holding several channels, see --axes)")
    pa.add_argument("--names", nargs="+", required=True, metavar="NAME",
                    help="one name per output channel, in input order (the names the pack's channel map binds)")
    pa.add_argument("--spacing", nargs=3, type=float, required=True, metavar=("DZ", "DY", "DX"),
                    help="the voxel size in micrometres: z step, then the y and x pixel size")
    pa.add_argument("--out", required=True, help="the OME-TIFF to write (for example X.ome.tif)")
    pa.add_argument("--axes", default="ZYX",
                    help="the declared axes of each input, outermost first, from C, Z, Y and X (default ZYX: a "
                         "z-stack of grey pages, or one plane). ZCYX: channels alternate page by page; ZYXC: each "
                         "page holds the channels as samples")
    for ax in ("z", "y", "x"):
        pa.add_argument(f"--crop-{ax}", nargs=2, type=int, default=None, metavar=("START", "STOP"),
                        help=f"keep {ax} pixels START to STOP-1 of the input (default: all)")
    pa.add_argument("--force", action="store_true",
                    help="read the pages as --axes declares even where the file's own ImageJ or OME metadata records "
                         "its channels and planes in another order")
    pa.add_argument("--json", action="store_true", help="print the result record as JSON")


def dispatch(args) -> int:
    from .assemble import AssembleError, assemble

    crop = {ax: list(getattr(args, f"crop_{ax}")) for ax in ("z", "y", "x") if getattr(args, f"crop_{ax}")}
    try:
        rec = assemble(args.channels, args.names, args.spacing, args.out, axes=args.axes, crop=crop or None,
                       force=args.force)
    except (AssembleError, FileNotFoundError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    if args.json:
        print(json.dumps(rec, indent=1))
        return 0
    c, z, y, x = rec["shape"]
    print(f"wrote {rec['path']}: {c} channel(s) {rec['channel_names']}, {z} plane(s) of {y} x {x} pixels "
          f"({rec['dtype']}), voxel {rec['spacing_zyx_um']} um" + (f", crop {rec['crop']}" if rec["crop"] else ""))
    print(f"  pixels_sha256 {rec['pixels_sha256']}")
    return 0
