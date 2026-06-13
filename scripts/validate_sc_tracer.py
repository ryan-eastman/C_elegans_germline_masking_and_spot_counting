#!/usr/bin/env python
"""Validate the pipeline's SC tracer against synthetic ground truth, with PERFECT nucleus labels
(isolating the tracer from the nucleus model). Does ridge+skeletonize recover the known
fragment counts and SC lengths? Tests control (intact) and heat (fragmented).

    python scripts/validate_sc_tracer.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from synth_nuclei import DX, DY, DZ, degrade, synthesize_scene  # noqa: E402

from germquant.sc.skeleton import trace_sc  # noqa: E402

SP = (DZ, DY, DX)


def main():
    for cond, frag in (("CONTROL", 0.0), ("HEAT", 0.7)):
        dapi, syp, label, gt = synthesize_scene(seed=3, frag_rate=frag)
        sypimg = degrade(syp, seed=30)
        _, det = trace_sc(sypimg, label, SP, intensity_percentile=90.0, min_fragment_length_um=0.5)
        gtpn = gt.groupby("nucleus_id").agg(gt_frags=("n_fragments", "sum"),
                                            gt_len=("sc_length_um", "sum")).reset_index()
        m = gtpn.merge(det[["nucleus_id", "n_fragments", "sc_total_length_um"]],
                       on="nucleus_id", how="left").fillna(0.0)
        cf = m[["gt_frags", "n_fragments"]].corr().iloc[0, 1]
        cl = m[["gt_len", "sc_total_length_um"]].corr().iloc[0, 1]
        print(f"=== {cond} (n={len(m)} nuclei) ===")
        print(f"  fragments/nucleus:  GT {m.gt_frags.mean():5.1f}   detected {m.n_fragments.mean():5.1f}   corr={cf:.2f}")
        print(f"  SC length/nucleus:  GT {m.gt_len.mean():5.1f}um detected {m.sc_total_length_um.mean():5.1f}um corr={cl:.2f}")
        print(f"  %% nuclei traced: {100 * (m.n_fragments > 0).mean():.0f}%%")


if __name__ == "__main__":
    main()
