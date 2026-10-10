#!/usr/bin/env python
"""Pick the nuclei a biologist follows by hand to validate the linking of a time-lapse series
(docs/TIMELAPSE.md, "Validating the linking").

    python scripts/pick_validation_nuclei.py results/SERIES [--seed 0] [--n 16 --tz 6 --pach 7 --edge 3] [--out picks.csv]

A stratified, seeded choice of nuclei at the first timepoint of the series: by default 16, of which 6 in the
transition zone, 7 in pachytene and 3 at the edge of the field (zones from the staging stage when it ran, else
from axis_position_um quantiles), alternately marked holdout (never looked at while tuning) and constraint
(used to tune). The worksheet <stem>__validation_picks.csv is written into the series folder (or --out): one
row per pick with hand_track_id, role, zone, t_index, z_um, y_um, x_um (raw-frame microns), nucleus_id. The
biologist opens the series in Imaris, finds each nucleus at t0 from its coordinates, follows it and appends
one row per later timepoint (hand_track_id, t_index, z_um, y_um, x_um). `germquant track-validate SERIES_DIR
--points <the completed file>` then scores the linking. The same seed gives the same picks, so pick BEFORE
tuning, or the set stops being blind.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from germquant.timelapse.track_validate import DEFAULTS, pick_validation_nuclei  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("series_dir", help="the results folder of `germquant run SERIES.nd2 --t ...`")
    ap.add_argument("--stem", default=None, help="the series stem when the folder holds several")
    ap.add_argument("--n", type=int, default=None, help=f"nuclei to pick (default {DEFAULTS['n_nuclei']})")
    ap.add_argument("--tz", type=int, default=None, help=f"of which transition zone (default {DEFAULTS['n_transition']})")
    ap.add_argument("--pach", type=int, default=None, help=f"of which pachytene (default {DEFAULTS['n_pachytene']})")
    ap.add_argument("--edge", type=int, default=None, help=f"of which at the edge of the field (default {DEFAULTS['n_edge']})")
    ap.add_argument("--edge-um", type=float, default=None,
                    help=f"a nucleus within this distance of the xy border is an edge candidate (default {DEFAULTS['edge_um']})")
    ap.add_argument("--seed", type=int, default=None, help=f"random seed (default {DEFAULTS['seed']})")
    ap.add_argument("--out", default=None, help="write the worksheet here (default: <series_dir>/<stem>__validation_picks.csv)")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    df = pick_validation_nuclei(a.series_dir, n=a.n, n_tz=a.tz, n_pach=a.pach, n_edge=a.edge, seed=a.seed,
                                edge_um=a.edge_um, stem=a.stem, out=a.out)
    counts = df["zone"].value_counts().to_dict()
    roles = df["role"].value_counts().to_dict()
    print(f"{len(df)} nuclei picked at t{int(df['t_index'].iloc[0]):03d} (seed {df.attrs['seed']}): "
          + ", ".join(f"{k} {v}" for k, v in counts.items()) + "; " + ", ".join(f"{k} {v}" for k, v in roles.items()))
    print(f"zones from {df.attrs['zone_source']}")
    print(df[["hand_track_id", "role", "zone", "nucleus_id", "z_um", "y_um", "x_um"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
