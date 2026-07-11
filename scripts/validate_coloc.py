"""Validate PGL-1 granule + SYP<->PGL-1 coloc output against Imaris ground truth.

Run AFTER you have (a) a pipeline result for a 4-channel image (its `*__granules.csv` and
`*__coloc.csv`) and (b) Imaris ground truth for the SAME image:
  * PGL-1 Surfaces  -> an .ims project, OR a CSV with columns z_um,y_um,x_um[,volume_um3]
  * SYP x PGL-1 Colocalization module stats -> a small JSON: {"manders_m1":.., "manders_m2":..,
    "pearson_r":..}  (any subset)

Example
-------
  python scripts/validate_coloc.py \
      --granules results/<img>__granules.csv --coloc results/<img>__coloc.csv \
      --imaris-ims E:/Ryan/<img>.ims  --imaris-coloc-json imaris_coloc.json \
      --out validation_coloc

It prints a summary and writes report.json (+ a granule-volume Bland-Altman when there are matches).
No ground truth is needed to RUN the pipeline; this only scores it once you have Imaris exports.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from germquant.validate import compare_coloc_metrics, compare_granules


def _load_imaris_granules(args) -> pd.DataFrame:
    if args.imaris_granules_csv:
        df = pd.read_csv(args.imaris_granules_csv)
        need = {"z_um", "y_um", "x_um"}
        if not need.issubset(df.columns):
            raise SystemExit(f"--imaris-granules-csv needs columns {need}; got {list(df.columns)}")
        return df
    if args.imaris_ims:
        from germquant.validate.imaris_ims import read_imaris_ims

        _spots, surfaces, _geom = read_imaris_ims(args.imaris_ims)
        # PGL-1 granules are exported as Imaris Surfaces (read as x/y/z_um + volume_um3)
        return surfaces
    return pd.DataFrame(columns=["z_um", "y_um", "x_um", "volume_um3"])


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--granules", required=True, help="our *__granules.csv")
    p.add_argument("--coloc", required=True, help="our *__coloc.csv")
    p.add_argument("--imaris-ims", help="Imaris .ims with PGL-1 Surfaces")
    p.add_argument("--imaris-granules-csv", help="alt to --imaris-ims: CSV z_um,y_um,x_um[,volume_um3]")
    p.add_argument("--imaris-coloc-json", help="JSON of Imaris Coloc-module metrics")
    p.add_argument("--operand", default="syp_aggregate", choices=["syp_aggregate", "sc_ribbon"])
    p.add_argument("--max-match-um", type=float, default=1.5)
    p.add_argument("--out", default="validation_coloc")
    args = p.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    our_g = pd.read_csv(args.granules)
    our_c = pd.read_csv(args.coloc)
    ref_g = _load_imaris_granules(args)
    imaris_metrics = json.loads(Path(args.imaris_coloc_json).read_text()) if args.imaris_coloc_json else None

    gstats = compare_granules(our_g, ref_g, max_match_um=args.max_match_um)
    metrics = compare_coloc_metrics(our_c, imaris_metrics, operand=args.operand)

    report = {"granules": gstats, "coloc_metrics": metrics.to_dict(orient="records")}
    (out / "report.json").write_text(json.dumps(report, indent=2, default=float))
    metrics.to_csv(out / "coloc_metrics_compare.csv", index=False)

    print("== PGL-1 granule detection vs Imaris ==")
    print(f"  ours={gstats['n_our']}  imaris={gstats['n_ref']}  matched={gstats['n_matched']}  "
          f"recall={gstats['recall_of_ref']:.3f}  precision={gstats['precision']:.3f}  "
          f"mean_dist={gstats['mean_match_dist_um']:.3f} um")
    vs = gstats.get("volume_stats") or {}
    if vs.get("n", 0) >= 2:
        print(f"  granule volume: CCC={vs.get('ccc', float('nan')):.3f}  "
              f"bias={vs.get('bias_mean_diff', float('nan')):.3f} um^3  MAE={vs.get('mae', float('nan')):.3f}")
    print("\n== SYP x PGL-1 coloc metrics (ours vs Imaris) ==")
    print(metrics.to_string(index=False))
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
