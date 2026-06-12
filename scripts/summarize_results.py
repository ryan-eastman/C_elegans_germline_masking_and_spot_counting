#!/usr/bin/env python
"""Quick Python-side summary of a germquant results tree: per-image readouts + a condition
(sex x treatment) comparison for the heat phenotype. Complements analysis_R/ (which makes the
publication figures). SC numbers are only meaningful once sc.intensity_percentile is calibrated
against ground truth — see docs/RUNBOOK.md.

    python scripts/summarize_results.py <results_dir>
"""
from __future__ import annotations

import glob
import os
import sys

import pandas as pd


def _per_image(results_dir: str) -> pd.DataFrame:
    rows = []
    for nf in glob.glob(os.path.join(results_dir, "**", "*__nuclei.csv"), recursive=True):
        if os.sep + "_done" + os.sep in nf:
            continue
        n = pd.read_csv(nf)
        if n.empty:
            continue
        meta = n.iloc[0]
        scf = nf.replace("__nuclei.csv", "__sc_per_nucleus.csv")
        sc = pd.read_csv(scf) if os.path.exists(scf) else pd.DataFrame(columns=["n_fragments", "sc_total_length_um"])
        traced = sc[sc["n_fragments"] > 0] if "n_fragments" in sc else sc
        rows.append({
            "image_id": meta.get("image_id"),
            "sex": meta.get("sex"), "treatment": meta.get("treatment"),
            "germ_cell": meta.get("germ_cell"), "n_nuclei": len(n),
            "mean_foci": round(n["n_foci"].mean(), 2) if "n_foci" in n else float("nan"),
            "sc_traced_pct": round(100 * (sc["n_fragments"] > 0).mean()) if len(sc) else 0,
            "mean_frags_traced": round(traced["n_fragments"].mean(), 2) if len(traced) else 0.0,
            "mean_sclen_traced_um": round(traced["sc_total_length_um"].mean(), 2) if len(traced) else 0.0,
        })
    return pd.DataFrame(rows).sort_values(["treatment", "sex", "image_id"], ignore_index=True)


def _by_condition(per_image: pd.DataFrame) -> pd.DataFrame:
    g = per_image.groupby(["germ_cell", "treatment"], dropna=False)
    out = g.agg(
        n_gonads=("image_id", "size"),
        nuclei=("n_nuclei", "mean"),
        foci=("mean_foci", "mean"),
        sc_frags=("mean_frags_traced", "mean"),
        sc_len_um=("mean_sclen_traced_um", "mean"),
    ).round(2).reset_index()
    return out


def main(results_dir: str) -> int:
    pi = _per_image(results_dir)
    if pi.empty:
        print(f"No __nuclei.csv found under {results_dir}")
        return 1
    print("=== per gonad ===")
    print(pi.to_string(index=False))
    print("\n=== by condition (germ_cell x treatment) — the heat-phenotype comparison ===")
    print(_by_condition(pi).to_string(index=False))
    print("\nNote: SC fragment counts/lengths are PRELIMINARY until sc.intensity_percentile is")
    print("calibrated vs Imaris ground truth. The HS-vs-noHS contrast within each sex is the test:")
    print("heat should raise SC fragmentation in spermatocytes (male) but not oocytes (herm).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else "results"))
