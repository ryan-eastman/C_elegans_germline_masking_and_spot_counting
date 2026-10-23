"""`germquant sc-straighten` and `germquant sc-view` (germquant.sc.straighten; docs/SC_STRAIGHTEN.md).

  germquant sc-straighten RESULTS_ROOT [IMAGE_IDS ...] [--config C | --profile P] [--set KEY=VALUE ...]
                          [--out DIR] [--geometry auto|import|both] [--paths-from ims:F|swc:F|csv:F ...]
                          [--refine-imported on|off] [--ends fit|as_drawn] [--nuclei germline|pachytene|all]
                          [--nucleus-ids N ...] [--workers N] [--force] [--no-resume] [--no-strips] [--no-overview]
                          [--source-map OLD=NEW ...] [--allow-sealed] [--dry-run]
  germquant sc-view RESULTS_DIR IMAGE_ID [--nucleus N] [--strips] [--stride K]

A post-run command: it reads finished results (nuclei.csv, the label TIFF, spot tables, the run's geometry and
channel map) and the source image, and writes only new files. Exit code 0 when at least one image finished (or
was up to date), 1 when none did, 2 for usage errors.
"""
from __future__ import annotations

import sys


def add_parsers(sub) -> None:
    ps = sub.add_parser("sc-straighten", help="trace, straighten and measure every SC ribbon of finished results "
                                              "(new files only; docs/SC_STRAIGHTEN.md)")
    ps.add_argument("results_root")
    ps.add_argument("image_ids", nargs="*", help="image ids or result-folder names (default: every finished image)")
    g = ps.add_mutually_exclusive_group(required=False)
    g.add_argument("--config", help="a config yaml; only its sc_straighten section is read")
    g.add_argument("--profile", help="a profile name; only its sc_straighten section is read")
    ps.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE",
                    help="override sc_straighten keys, e.g. follow.seed_min=3.5 (enters params_hash)")
    ps.add_argument("--out", help="write under this folder, mirroring the tree (required inside a golden tree)")
    ps.add_argument("--geometry", choices=["auto", "import", "both"], default="auto",
                    help="auto: trace the ribbons; import: measure imported paths only; both")
    ps.add_argument("--paths-from", nargs="*", default=[], metavar="KIND:FILE",
                    help="imported geometry: ims:FILE (Imaris Filaments), swc:FILE or csv:FILE")
    ps.add_argument("--refine-imported", choices=["on", "off"], default="on")
    ps.add_argument("--ends", choices=["fit", "as_drawn"], default="fit", help="end rule for imported paths")
    ps.add_argument("--nuclei", choices=["germline", "pachytene", "all"], default=None)
    ps.add_argument("--nucleus-ids", nargs="*", type=int, default=None)
    ps.add_argument("--workers", type=int, default=None)
    ps.add_argument("--force", action="store_true", help="rerun even when the record is up to date")
    ps.add_argument("--no-resume", action="store_true", help="the same as --force")
    ps.add_argument("--no-strips", action="store_true")
    ps.add_argument("--no-overview", action="store_true")
    ps.add_argument("--source-map", nargs="*", default=[], metavar="OLD=NEW",
                    help="rewrite the recorded source path prefix OLD to NEW (moved or NAS images)")
    ps.add_argument("--allow-sealed", action="store_true",
                    help="run images of the sealed SC validation split (only scripts/benchmark_sc.py --final)")
    ps.add_argument("--dry-run", action="store_true", help="list what would run")

    pv = sub.add_parser("sc-view", help="open sc-straighten outputs in napari (needs the napari-germquant environment)")
    pv.add_argument("results_dir")
    pv.add_argument("image_id")
    pv.add_argument("--nucleus", type=int, default=None)
    pv.add_argument("--strips", action="store_true", help="also show the straightened strips of --nucleus")
    pv.add_argument("--stride", type=int, default=2, help="xy stride of the whole-image SYP layer")


def _params(args):
    from .sc.straighten import params as P

    cfg = None
    if getattr(args, "config", None) or getattr(args, "profile", None):
        from .config import find_profile, load_config

        path = find_profile(args.profile) if args.profile else args.config
        cfg = load_config(path)
    return P.from_config(cfg, list(args.set or []))


def dispatch(args) -> int:
    if args.cmd == "sc-view":
        return _view(args)
    try:
        params = _params(args)
    except (KeyError, ValueError, FileNotFoundError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
    smap = []
    for s in args.source_map or []:
        if "=" not in s:
            print(f"ERROR: --source-map wants OLD=NEW, got {s!r}", file=sys.stderr)
            return 2
        old, new = s.split("=", 1)
        smap.append((old, new))
    if args.geometry != "auto" and not args.paths_from:
        print("ERROR: --geometry import|both needs --paths-from", file=sys.stderr)
        return 2
    from .sc.straighten.runner import run_results

    nuclei = args.nuclei or params["nuclei"]
    res = run_results(args.results_root, args.image_ids or None, params=params, out=args.out,
                      geometry=args.geometry, paths_from=args.paths_from, refine_imported=args.refine_imported == "on",
                      ends=args.ends, nuclei=nuclei, nucleus_ids=args.nucleus_ids, workers=args.workers,
                      force=bool(args.force or args.no_resume), strips=False if args.no_strips else None,
                      overview=False if args.no_overview else None, source_map=smap,
                      allow_sealed=args.allow_sealed, dry_run=args.dry_run)
    ok = [r for r in res if r.status in ("done", "skipped", "dry_run")]
    return 0 if ok else 1


def _view(args) -> int:
    try:
        import napari  # noqa: F401
    except Exception:  # noqa: BLE001
        print("sc-view needs napari, which is not in this environment: run it from the napari-germquant "
              "environment (pip install napari there, or the viz extra).", file=sys.stderr)
        return 1
    import json
    from pathlib import Path

    import napari

    from .sc.straighten import layers as LY

    rd = Path(args.results_dir)
    rec = json.loads((rd / f"{args.image_id}__sc_straighten.json").read_text(encoding="utf-8"))
    v = napari.Viewer(title=f"sc-straighten: {args.image_id}")
    try:
        from .sc.straighten.inputs import read_roles

        src = rec["inputs"]["source"]["path"]
        ce_idx = rec["inputs"]["channel_map"]["central_element"]
        arr = read_roles(src, {"central_element": ce_idx}, z_range=rec["inputs"].get("z_range"))["central_element"]
        k = max(1, int(args.stride))
        import pandas as pd

        nuc = pd.read_csv(Path(rec["inputs"]["nuclei_csv"]["path"]), nrows=1)
        sp = tuple(float(nuc[c].iloc[0]) for c in ("voxel_dz_um", "voxel_dy_um", "voxel_dx_um"))
        v.add_image(arr[:, ::k, ::k], name="SYP", scale=(sp[0], sp[1] * k, sp[2] * k), colormap="gray")
    except Exception as e:  # noqa: BLE001 - the strands still show without the image
        print(f"(no SYP layer: {e})")
    for data, kw, kind in LY.strand_layers(rd, args.image_id, [args.nucleus] if args.nucleus else None):
        getattr(v, f"add_{kind}")(data, **kw)
    if args.strips and args.nucleus:
        v2 = napari.Viewer(title=f"strips: {args.image_id} nucleus {args.nucleus}")
        for data, kw, kind in LY.strip_layers(rd, args.image_id, args.nucleus):
            getattr(v2, f"add_{kind}")(data, **kw)
    napari.run()
    return 0
