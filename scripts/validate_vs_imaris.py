#!/usr/bin/env python
"""Validate the pipeline against REAL Imaris ground truth — segmentation F1, and per-nucleus
agreement (CCC / Bland-Altman / correlation) for SC length, RAD-51 foci, and the fragmentation
index. Reuses germquant.validate (agreement_stats, bland_altman_plot, segmentation_metrics).

INPUTS (the export spec):
  * pipeline output dir with <image>__nuclei.csv and <image>__nuclei_labels.tif
  * Imaris GT 3D label TIFF (one integer per nucleus, same dims as DAPI)
  * Imaris GT CSV: columns centroid_x_um,centroid_y_um,centroid_z_um + any of
    sc_length_um, n_sc_fragments, n_rad51_foci  (nuclei matched pipeline<->GT by centroid)

    python scripts/validate_vs_imaris.py --selftest          # validate the harness on synthetic
    python scripts/validate_vs_imaris.py --results results_germline --image <img> \
        --gt-labels imaris_labels.tif --gt-csv imaris_stats.csv
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

from germquant.validate import agreement_stats, bland_altman_plot, segmentation_metrics

# pipeline column -> GT column for each per-nucleus readout
READOUTS = [
    ("SC length (um)", "sc_total_length_um", "sc_length_um"),
    ("RAD-51 foci", "n_foci", "n_rad51_foci"),
    ("fragmentation index vs GT frag count", "sc_fragmentation_index", "n_sc_fragments"),
]
_XYZ = ["centroid_x_um", "centroid_y_um", "centroid_z_um"]


def match_by_centroid(pipe, gt, tol_um=3.0):
    """One-to-one nearest-centroid match (<= tol_um). Returns (pipe_idx, gt_idx) aligned arrays."""
    from scipy.spatial import cKDTree
    pc = pipe[_XYZ].to_numpy(float)
    gc = gt[_XYZ].to_numpy(float)
    dist, idx = cKDTree(pc).query(gc, distance_upper_bound=tol_um)
    order = np.argsort(dist)                       # resolve collisions nearest-first
    used, pi, gi = set(), [], []
    for g in order:
        if not np.isfinite(dist[g]) or idx[g] in used:
            continue
        used.add(idx[g])
        pi.append(idx[g])
        gi.append(g)
    return np.array(pi, int), np.array(gi, int)


def report(pipe, gt, pipe_labels, gt_labels, out_dir, tag):
    print(f"\n===== {tag} =====")
    if pipe_labels is not None and gt_labels is not None:
        m = segmentation_metrics(pipe_labels, gt_labels, iou_threshold=0.5)
        print(f"SEGMENTATION  F1@0.5={m['f1']:.3f}  matched-IoU={m['mean_iou']:.3f}  "
              f"pred={m['n_pred']} gt={m['n_gt']} (tp{m['tp']} fp{m['fp']} fn{m['fn']})")

    pi, gi = match_by_centroid(pipe, gt)
    print(f"matched {len(pi)}/{len(gt)} GT nuclei to pipeline nuclei by centroid (<=3um)")
    for label, pcol, gcol in READOUTS:
        if pcol not in pipe.columns or gcol not in gt.columns:
            continue
        p = pipe.iloc[pi][pcol].to_numpy(float)
        t = gt.iloc[gi][gcol].to_numpy(float)
        ok = np.isfinite(p) & np.isfinite(t)
        if ok.sum() < 2:
            continue
        s = agreement_stats(p[ok], t[ok])
        print(f"  {label:<38} n={s['n']:4d}  CCC={s.get('ccc', float('nan')):.3f}  "
              f"r={s.get('pearson_r', float('nan')):.3f}  bias={s['bias_mean_diff']:+.2f}  MAE={s['mae']:.2f}")
        if out_dir:
            bland_altman_plot(p[ok], t[ok], os.path.join(out_dir, f"ba_{pcol}.png"),
                              title=f"{tag}: {label}")
    return len(pi)


def _selftest():
    """Route the synthetic SC ground truth through the harness (no real data needed)."""
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
    from synth_nuclei import DX, DY, DZ, degrade, synthesize_scene
    from germquant.measure import measure_objects
    from germquant.sc.skeleton import trace_sc
    sp = (DZ, DY, DX)
    dapi, syp, label, sc_gt = synthesize_scene(seed=3, frag_rate=0.5)
    meas = measure_objects(label, {"dna": dapi}, sp).rename(columns={"label": "nucleus_id"})
    _, per = trace_sc(degrade(syp, seed=30), label, sp)
    pipe = meas.merge(per[["nucleus_id", "sc_total_length_um"]], on="nucleus_id", how="left")
    gtpn = sc_gt.groupby("nucleus_id").agg(sc_length_um=("sc_length_um", "sum")).reset_index()
    # build a GT table with centroids (from the perfect labels) — stands in for the Imaris export
    gt = meas.merge(gtpn, on="nucleus_id")[_XYZ + ["sc_length_um"]]
    print("SELF-TEST: synthetic SC ground truth through the Imaris harness (labels identical -> F1 1.0;"
          " SC length CCC should match validate_sc_tracer ~0.6)")
    report(pipe, gt, label, label, None, "selftest")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--results")
    ap.add_argument("--image")
    ap.add_argument("--gt-labels")
    ap.add_argument("--gt-csv")
    a = ap.parse_args()
    if a.selftest:
        _selftest()
        return

    import tifffile
    ncsv = glob.glob(f"{a.results}/{a.image}*__nuclei.csv")[0]
    pipe = pd.read_csv(ncsv)
    if "in_germline" in pipe.columns:                       # compare only the gonad nuclei
        pipe = pipe[pipe["in_germline"] != False]           # noqa: E712
    lbl = glob.glob(f"{a.results}/{a.image}*__nuclei_labels.tif")
    pipe_labels = tifffile.imread(lbl[0]) if lbl else None
    gt_labels = tifffile.imread(a.gt_labels) if a.gt_labels else None
    gt = pd.read_csv(a.gt_csv)
    report(pipe, gt, pipe_labels, gt_labels, a.results, a.image)


if __name__ == "__main__":
    main()
