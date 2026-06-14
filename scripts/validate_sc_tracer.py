#!/usr/bin/env python
"""Validate the v2 SC tracer against synthetic ground truth with PERFECT nucleus labels (isolating
the tracer). Reports what is recoverable: SC LENGTH (reliable) and the FRAGMENTATION INDEX
(population-level control-vs-heat). Per-nucleus fragment COUNT is a lower bound (strands merge in 3D).

    python scripts/validate_sc_tracer.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from synth_nuclei import DX, DY, DZ, degrade, synthesize_scene  # noqa: E402

from germquant.sc.skeleton import trace_sc  # noqa: E402

SP = (DZ, DY, DX)


def cohens_d(a, b):
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((b.mean() - a.mean()) / (sp + 1e-12))


def main():
    print("LENGTH (the reliable readout): detected vs GT, per-nucleus correlation")
    idx = {}
    for cond, frag in (("CONTROL", 0.0), ("HEAT", 1.0)):
        dapi, syp, label, gt = synthesize_scene(seed=3, frag_rate=frag)
        sypimg = degrade(syp, seed=30)
        _, det = trace_sc(sypimg, label, SP)
        gtpn = gt.groupby("nucleus_id").agg(gt_len=("sc_length_um", "sum")).reset_index()
        m = gtpn.merge(det[["nucleus_id", "sc_total_length_um", "sc_fragmentation_index"]],
                       on="nucleus_id", how="left").fillna(0.0)
        cl = m[["gt_len", "sc_total_length_um"]].corr().iloc[0, 1]
        print(f"  {cond:8} GT_len {m.gt_len.mean():5.1f}um  detected {m.sc_total_length_um.mean():5.1f}um  "
              f"corr={cl:.2f}   frag_index mean={m.sc_fragmentation_index.mean():.3f}")
        idx[cond] = m.sc_fragmentation_index.values

    d = cohens_d(idx["CONTROL"], idx["HEAT"])
    print(f"\nFRAGMENTATION INDEX (population control-vs-heat): "
          f"control {idx['CONTROL'].mean():.3f} -> heat {idx['HEAT'].mean():.3f}  Cohen's d = {d:+.2f}")
    print("  (per-nucleus noisy; population-level proxy for the heat phenotype — calibrate vs Imaris.)")


if __name__ == "__main__":
    main()
