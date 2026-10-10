"""Benchmark the RAD-51 spot detector against the Imaris expert spots over a manifest of gonads.

    python scripts/benchmark_spots.py --manifest benchmark_manifest.json --out benchmark_spots_out
        [--split noHS] [--tol-xy 0.3] [--tol-z 0.4] [--surface-tol 2.5] [--overlap-frac 0.5]
        [--traced-rule auto|overlap|centroid] [--max-candidates N] [--max-dense-cells N]
        [--n-boot 1000] [--seed 0] [--truth auto|files|ims] [--write-matches]

The manifest is the one scripts/benchmark_manifest.py builds: entries with tag, split, results_dir,
truth_dir, nd2, ims (germquant.benchmark.manifest documents the shape). For every gonad the pipeline
outputs are read from results_dir and the truth from truth_dir (<stem>__imaris_spots.csv and friends,
the truth file contract of the Imaris reader package). When those files are absent and --truth is
"auto" or "ims", the truth tables are built straight from the .ims here (``truth_from_ims``, the
FALLBACK path: no label image, no quality, no spot class, so the detector-only and no-Class-B scores
stay empty and the traced subset uses the centroid rule); per_gonad.csv says which path each gonad
took in ``truth_source``.

Outputs in --out: per_gonad.csv (one row per gonad), per_nucleus.csv (one row per germline nucleus
with expert_n, ours_n and traced), pooled.json (pooled scores with cluster-bootstrap intervals by
gonad, overall and by split, plus the skipped gonads and the settings). --write-matches adds
per_spot_matches.csv with every matched (our spot, expert spot) pair. docs/BENCHMARK.md explains
every column and key. Exit code 0 when at least one gonad was scored, 1 otherwise.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from germquant.benchmark import spots as bs  # noqa: E402
from germquant.benchmark.manifest import (  # noqa: E402
    Truth,
    find_pipeline_outputs,
    imaris_um_to_voxel_index,
    load_labels,
    load_manifest,
    load_truth,
    our_voxel_zyx,
    voxel_index_to_um,
)
from germquant.benchmark.matching import DEFAULT_MAX_CANDIDATES, DEFAULT_MAX_DENSE_CELLS  # noqa: E402
from germquant.benchmark.report import write_outputs  # noqa: E402


def truth_from_ims(ims_path: str | Path, voxel_zyx) -> Truth:
    """FALLBACK truth builder, used only when the truth files of the Imaris reader package are absent.

    Builds the truth tables straight from the .ims with ``read_imaris_ims`` (image-relative Imaris
    micrometres) and the documented frame conversion: voxel index = floor(imaris_um / imaris_voxel)
    with Imaris's own rounded voxel size, our micrometres = index * the .nd2 voxel size. Produces the
    contract columns with what the reader gives: no label image, quality NaN, spot_class empty,
    label and n_voxels_mask NaN. ``imaris_*`` hold the image-relative Imaris values (Imaris minus
    ExtMin), which is what ``read_imaris_ims`` returns.
    """
    from germquant.validate.imaris_ims import read_imaris_ims

    spots, surfaces, geom = read_imaris_ims(ims_path)
    vox_im_xyz = np.asarray(geom["voxel_xyz_um"], float)
    vox_our_xyz = np.asarray(voxel_zyx, float)[::-1]

    def convert(df: pd.DataFrame):
        um_im = df[["x_um", "y_um", "z_um"]].to_numpy(float).reshape(-1, 3)
        idx = imaris_um_to_voxel_index(um_im, vox_im_xyz)
        return idx, voxel_index_to_um(idx, vox_our_xyz), um_im

    sidx, sum_, sraw = convert(spots)
    tspots = pd.DataFrame({
        "spot_id": np.arange(len(spots)), "x_um": sum_[:, 0], "y_um": sum_[:, 1], "z_um": sum_[:, 2],
        "xi": sidx[:, 0], "yi": sidx[:, 1], "zi": sidx[:, 2],
        "radius_um": spots["radius_um"].to_numpy(float) if len(spots) else np.zeros(0),
        "quality": np.full(len(spots), np.nan), "spot_class": [""] * len(spots),
        "imaris_x_um": sraw[:, 0], "imaris_y_um": sraw[:, 1], "imaris_z_um": sraw[:, 2]})
    nidx, num, nraw = convert(surfaces)
    tnuclei = pd.DataFrame({
        "surface_id": np.arange(len(surfaces)), "label": np.full(len(surfaces), np.nan),
        "x_um": num[:, 0], "y_um": num[:, 1], "z_um": num[:, 2],
        "xi": nidx[:, 0], "yi": nidx[:, 1], "zi": nidx[:, 2],
        "volume_um3": surfaces["volume_um3"].to_numpy(float) if len(surfaces) else np.zeros(0),
        "n_voxels_mask": np.full(len(surfaces), np.nan),
        "imaris_x_um": nraw[:, 0], "imaris_y_um": nraw[:, 1], "imaris_z_um": nraw[:, 2]})
    params = {"truth_source": "ims_fallback", "ims": str(ims_path),
              "voxel_xyz_um_imaris": [float(v) for v in vox_im_xyz],
              "voxel_zyx_um_ours": [float(v) for v in voxel_zyx],
              "shape_xyz": [int(v) for v in geom["shape_xyz"]],
              "ext_min_xyz": [float(v) for v in geom["ext_min_xyz"]],
              "n_spots": int(len(spots)), "n_surfaces": int(len(surfaces))}
    return Truth(Path(ims_path).stem, tspots, tnuclei, None, params, source="ims_fallback")


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, help="benchmark manifest JSON")
    ap.add_argument("--out", required=True, help="output folder")
    ap.add_argument("--split", default=None, help="only entries with this split")
    ap.add_argument("--tol-xy", type=float, default=bs.DEFAULT_TOL_XY_UM, help="cylinder radius in xy (um)")
    ap.add_argument("--tol-z", type=float, default=bs.DEFAULT_TOL_Z_UM, help="cylinder half-height in z (um)")
    ap.add_argument("--surface-tol", type=float, default=bs.DEFAULT_SURFACE_TOL_UM,
                    help="centroid distance (um) for the traced subset without a truth label image")
    ap.add_argument("--overlap-frac", type=float, default=bs.DEFAULT_OVERLAP_FRAC,
                    help="fraction of our nucleus inside one truth label for the traced subset")
    ap.add_argument("--traced-rule", choices=["auto", "overlap", "centroid", "none"], default="auto")
    ap.add_argument("--max-candidates", type=int, default=DEFAULT_MAX_CANDIDATES,
                    help="candidate pairs per gonad before greedy matching is used instead")
    ap.add_argument("--max-dense-cells", type=int, default=DEFAULT_MAX_DENSE_CELLS,
                    help="largest dense matrix (cells) solved exactly per connected component")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--truth", choices=["auto", "files", "ims"], default="auto",
                    help="truth files from truth_dir, the .ims fallback, or files first then the fallback")
    ap.add_argument("--write-matches", action="store_true", help="also write per_spot_matches.csv")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    entries = load_manifest(args.manifest)
    if args.split is not None:
        entries = [e for e in entries if e.split == args.split]
    rows, per_nuc, matches, skipped = [], [], [], []
    for e in entries:
        t0 = time.time()
        try:
            outputs = find_pipeline_outputs(e.results_dir)
        except (FileNotFoundError, ValueError) as err:
            skipped.append({"tag": e.tag, "reason": f"pipeline outputs: {err}"})
            print(f"  skip {e.tag}: {err}")
            continue
        if outputs.labels_path is None:
            skipped.append({"tag": e.tag, "reason": "no __nuclei_labels.tif in results_dir"})
            print(f"  skip {e.tag}: no nucleus label image")
            continue
        truth = None
        if args.truth in ("auto", "files"):
            truth = load_truth(e.truth_dir, e.truth_stem)
        if truth is None and args.truth in ("auto", "ims"):
            if e.ims is None or not Path(e.ims).exists():
                skipped.append({"tag": e.tag, "reason": "no truth files and no .ims to fall back on"})
                print(f"  skip {e.tag}: no truth files and no .ims")
                continue
            voxel, _ = our_voxel_zyx(outputs, e, None)
            truth = truth_from_ims(e.ims, voxel)
        if truth is None:
            skipped.append({"tag": e.tag, "reason": f"no truth files for stem {e.truth_stem!r} in {e.truth_dir}"})
            print(f"  skip {e.tag}: no truth files")
            continue
        voxel, vsrc = our_voxel_zyx(outputs, e, truth)
        labels = load_labels(outputs.labels_path)
        truth_labels = load_labels(truth.labels_path) if truth.has_labels else None
        row, pn, mt = bs.benchmark_gonad(
            outputs, truth, labels, voxel, tag=e.tag, split=e.split, truth_labels=truth_labels,
            tol_xy_um=args.tol_xy, tol_z_um=args.tol_z, surface_tol_um=args.surface_tol,
            overlap_frac=args.overlap_frac, traced_rule=args.traced_rule,
            max_candidates=args.max_candidates, max_dense_cells=args.max_dense_cells)
        row["voxel_source"] = vsrc
        row["elapsed_s"] = round(time.time() - t0, 1)
        rows.append(row)
        per_nuc.append(pn)
        matches.append(mt)
        print(f"  {e.tag} [{truth.source}]: nuclei all r={row['pearson_r_all']:.3f} ccc={row['ccc_all']:.3f} "
              f"bias={row['bias_all']:+.2f} (expert {row['mean_expert_all']:.2f} / ours {row['mean_ours_all']:.2f}) | "
              f"traced n={row['n_traced']} r={row['pearson_r_traced']:.3f} ccc={row['ccc_traced']:.3f} | "
              f"spots P/R/F1 {row['precision']:.3f}/{row['recall']:.3f}/{row['f1']:.3f} "
              f"({row['n_matched']} of ours {row['n_our_spots']}, expert {row['n_expert_spots']}; "
              f"{row['matching_method']}) | {row['elapsed_s']} s")
        del labels, truth_labels
    if not rows:
        print("no gonad could be scored")
        write_outputs(args.out, pd.DataFrame(), pd.DataFrame(), {"n_gonads": 0, "skipped": skipped})
        return 1
    per_gonad = pd.DataFrame(rows)
    per_nucleus = pd.concat(per_nuc, ignore_index=True)
    pooled = bs.pool(per_gonad, per_nucleus, n_boot=args.n_boot, seed=args.seed)
    pooled["skipped"] = skipped
    pooled["manifest"] = str(args.manifest)
    pooled["settings"].update({"max_candidates": args.max_candidates, "max_dense_cells": args.max_dense_cells,
                               "traced_rule": args.traced_rule, "truth": args.truth})
    extra = {"per_spot_matches": pd.concat(matches, ignore_index=True)} if args.write_matches else None
    paths = write_outputs(args.out, per_gonad, per_nucleus, pooled, extra)
    pa, pt, sa = pooled["per_nucleus_all"], pooled["per_nucleus_traced"], pooled["spots_all"]
    print(f"POOLED {len(rows)} gonads: per-nucleus all n={pa['n_nuclei']} ccc={pa['ccc']['point']:.3f} "
          f"[{pa['ccc']['lo']:.3f}, {pa['ccc']['hi']:.3f}] r={pa['pearson_r']['point']:.3f} "
          f"bias={pa['bias']['point']:+.2f}; traced n={pt['n_nuclei']} ccc={pt['ccc']['point']:.3f} "
          f"[{pt['ccc']['lo']:.3f}, {pt['ccc']['hi']:.3f}]; spots F1={sa['f1']['point']:.3f} "
          f"[{sa['f1']['lo']:.3f}, {sa['f1']['hi']:.3f}] P={sa['precision']['point']:.3f} R={sa['recall']['point']:.3f}")
    print("outputs ->", paths["pooled"].parent)
    return 0


if __name__ == "__main__":
    sys.exit(main())
