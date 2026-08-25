"""Sweep the coloc-sensitive parameters on ONE real 4-channel .nd2 and report how the headline
SYP<->PGL-1 numbers move. Use this Monday (before/after you have Imaris GT) to see which knobs the
result is sensitive to and to pick a defensible operating point.

Swept knobs (the ones that actually move the numbers):
  * granule.thresholding_method   (triangle | otsu | li)
  * granule.min_volume_um3        (speckle floor)
  * coloc.region_dilation_um      (perinuclear shell thickness — P-granules sit just outside nuclei)

Spot detection is skipped by default (--with-spots to keep it) since coloc doesn't need RAD-51, and
the null is shrunk (--n-random) for speed. Output: a tidy sweep.csv + a printed table.

Example
-------
  python scripts/coloc_param_sweep.py DATA/<img>.nd2 --config config/config.yaml \
      --out sweep_out --xy-stride 2
"""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import pandas as pd

from germquant.config import load_config
from germquant.pipeline import process_image


def _headline(res, dil, thr, minv):
    coloc = res["tables"]["coloc"]
    granules = res["tables"]["granules"]
    agg = coloc[coloc["sc_operand"] == "syp_aggregate"] if not coloc.empty else coloc
    a = agg.iloc[0] if len(agg) else {}
    g = lambda k: float(a[k]) if k in getattr(a, "index", []) else float("nan")  # noqa: E731
    return {
        "threshold": thr, "min_volume_um3": minv, "region_dilation_um": dil,
        "n_granules": int(len(granules)),
        "frac_granules_overlapping_syp_agg": g("frac_granules_overlapping_sc"),
        "manders_m1": g("manders_m1"), "manders_m2": g("manders_m2"),
        "overlap_volume_um3": g("overlap_volume_um3"), "overlap_pvalue": g("overlap_pvalue"),
        "qc_flags": ";".join(res.get("qc_flags", [])),
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("nd2")
    p.add_argument("--config", required=True)
    p.add_argument("--out", default="coloc_sweep")
    p.add_argument("--xy-stride", type=int, default=1)
    p.add_argument("--thresholds", nargs="+", default=["threshold_triangle", "threshold_otsu"])
    p.add_argument("--min-vol", nargs="+", type=float, default=[0.03, 0.1])
    p.add_argument("--dilation", nargs="+", type=float, default=[1.0, 1.5, 2.0])
    p.add_argument("--n-random", type=int, default=50)
    p.add_argument("--with-spots", action="store_true", help="also run SpotMAX (slower; not needed)")
    args = p.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    combos = list(itertools.product(args.thresholds, args.min_vol, args.dilation))
    print(f"Sweeping {len(combos)} combinations on {Path(args.nd2).name} ...")

    rows = []
    for i, (thr, minv, dil) in enumerate(combos, 1):
        cfg = load_config(args.config)
        cfg._data["segmentation"]["nuclei"]["method"] = cfg.get("segmentation.nuclei.method", "auto")
        cfg.set("spots.enabled", bool(args.with_spots))
        cfg.set("granule.thresholding_method", thr)
        cfg.set("granule.min_volume_um3", minv)
        cfg.set("coloc.region_dilation_um", dil)
        cfg.set("coloc.n_random", args.n_random)
        combo_out = out / f"combo_{i:02d}"
        print(f"[{i}/{len(combos)}] threshold={thr} min_vol={minv} dilation={dil}")
        try:
            res = process_image(args.nd2, cfg, combo_out, xy_stride=args.xy_stride, prov=None)
            rows.append(_headline(res, dil, thr, minv))
        except Exception as e:  # noqa: BLE001
            print(f"    FAILED: {type(e).__name__}: {e}")
            rows.append({"threshold": thr, "min_volume_um3": minv, "region_dilation_um": dil,
                         "n_granules": -1, "qc_flags": f"EXCEPTION:{e}"})

    df = pd.DataFrame(rows)
    df.to_csv(out / "sweep.csv", index=False)
    print("\n" + df.to_string(index=False))
    print(f"\n-> {out / 'sweep.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
