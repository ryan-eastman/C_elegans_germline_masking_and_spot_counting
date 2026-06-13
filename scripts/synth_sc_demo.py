#!/usr/bin/env python
"""Validate the SYP/SC channel + ground truth: same nuclei, control (intact SC) vs heat
(fragmented SC), with the KNOWN per-nucleus fragment counts. This is the ground truth the SC
tracer will be validated against -- the heat phenotype with answers we control.

    python scripts/synth_sc_demo.py
"""
import os
import sys

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from synth_nuclei import degrade, synthesize_scene  # noqa: E402


def norm(a):
    lo, hi = np.percentile(a, 1), np.percentile(a, 99.5)
    return np.clip((a - lo) / (hi - lo + 1e-9), 0, 1)


def main():
    dapi, sypC, label, gtC = synthesize_scene(seed=1, frag_rate=0.0)   # control: intact SCs
    _, sypH, _, gtH = synthesize_scene(seed=1, frag_rate=0.7)          # heat: same nuclei, broken
    dimg = degrade(dapi, seed=10)
    scimg = degrade(sypC, seed=11)
    shimg = degrade(sypH, seed=12)
    z = dimg.shape[0] // 2

    def mip(v):
        return v[max(0, z - 3):z + 4].max(0)

    pcC = gtC.groupby("nucleus_id").n_fragments.sum()
    pcH = gtH.groupby("nucleus_id").n_fragments.sum()
    print("SC ground truth (per nucleus, frag pieces):")
    print("  CONTROL: mean=%.1f median=%.0f  (= n chromosomes, each 1 intact SC)" % (pcC.mean(), pcC.median()))
    print("  HEAT:    mean=%.1f median=%.0f  (SCs broken into pieces)" % (pcH.mean(), pcH.median()))
    print("  per-chromosome SC length: mean=%.1f um (the SC tracer target)" % gtC.sc_length_um.mean())

    ys, xs = slice(70, 280), slice(70, 280)
    fig, ax = plt.subplots(1, 3, figsize=(18, 6))
    panels = [(mip(dimg), "DAPI (chromatin) — nucleus model"),
              (mip(scimg), "SYP / SC — CONTROL intact (%.0f frag/nuc)" % pcC.median()),
              (mip(shimg), "SYP / SC — HEAT fragmented (%.0f frag/nuc)" % pcH.median())]
    for a, v, t in zip(ax, *zip(*panels)):
        a.imshow(norm(v[ys, xs]), cmap="gray"); a.set_title(t); a.axis("off")
    fig.tight_layout(); fig.savefig("results_fullres/synth_sc.png", dpi=130, bbox_inches="tight")
    print("wrote results_fullres/synth_sc.png")


if __name__ == "__main__":
    main()
