"""Development calibration of sc_straighten on synthetic DEVELOPMENT seeds (docs/SC_VALIDATION_PROTOCOL.md, 2.3:
seeds below 7,000,000 only; the sealed ranges are refused by germquant.sc.straighten.synth).

    python scripts/calibrate_sc_straighten.py --out DIR [--grid follow.seed_min=2.5,3,3.5,4 follow.lo=1,1.5 ...]
           [--families iso dense gaps foci pairs noise snr] [--n 30] [--seed0 1000] [--procs 6]

Every configuration of the grid (the cartesian product of the given values; tunable keys only, DESIGN section 8)
runs the same development scenes of every family; one JSON line per configuration goes to DIR/calibration.jsonl
with the S1 metrics of each family (germquant.benchmark.sc_ribbons.summarise_scores). Synthetic sweeps on
development seeds are not limited by the protocol's tuning budget (section 9).
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
from multiprocessing import get_context
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from germquant.sc.straighten import params as PRM  # noqa: E402

TUNABLE = PRM.TUNABLE_PREFIXES


def _one(task):
    fam, seed, sets = task
    import warnings

    warnings.simplefilter("ignore")
    from germquant.benchmark import sc_ribbons as B
    from germquant.sc.straighten import synth
    from germquant.sc.straighten.pipeline import straighten_image

    P = PRM.resolve(None, sets)
    sc = synth.make_scene(fam, seed)
    roles = {k: v for k, v in sc.channels.items() if k != "central_element"}
    spots = {"spots": sc.spots} if len(sc.spots) else None
    t0 = time.perf_counter()
    res = straighten_image(sc.channels["central_element"], sc.labels, sc.spacing, sc.nuclei, params=P, roles=roles,
                           spot_tables=spots, image_id=f"{fam}_{seed}")
    s = B.score_scene(sc, res)
    s.extra["elapsed_s"] = time.perf_counter() - t0
    s.extra["kind"] = sc.info.get("pair", "")
    return fam, seed, tuple(sets), s


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--grid", nargs="*", default=[], metavar="KEY=V1,V2,...")
    ap.add_argument("--fixed", nargs="*", default=[], metavar="KEY=VALUE", help="applied to every configuration")
    ap.add_argument("--families", nargs="*", default=["iso", "dense", "gaps", "foci", "pairs", "noise"])
    ap.add_argument("--n", type=int, default=30, help="scenes per family")
    ap.add_argument("--seed0", type=int, default=1000)
    ap.add_argument("--procs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--allow-fixed", action="store_true",
                    help="also sweep keys the design fixes (an exploratory sweep; recorded as such in the output)")
    a = ap.parse_args(argv)
    keys, vals = [], []
    for g in a.grid:
        k, v = g.split("=", 1)
        if not a.allow_fixed and not any(k == t or (t.endswith(".") and k.startswith(t)) for t in TUNABLE):
            raise SystemExit(f"{k} is not a tunable key ({TUNABLE})")
        keys.append(k)
        vals.append(v.split(","))
    configs = [list(a.fixed) + [f"{k}={x}" for k, x in zip(keys, combo)] for combo in itertools.product(*vals)] \
        if keys else [list(a.fixed)]
    seeds = list(range(a.seed0, a.seed0 + a.n))
    tasks = [(fam, s, tuple(c)) for c in configs for fam in a.families for s in seeds]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    from germquant.benchmark import sc_ribbons as B

    by = {}
    t0 = time.time()
    with get_context("spawn").Pool(a.procs) as pool:
        for k, (fam, seed, sets, s) in enumerate(pool.imap_unordered(_one, tasks, chunksize=2)):
            by.setdefault(sets, {}).setdefault(fam, []).append(s)
            if (k + 1) % 50 == 0:
                print(f"{k + 1}/{len(tasks)} scenes, {time.time() - t0:.0f} s", flush=True)
    with open(out / "calibration.jsonl", "a", encoding="utf-8") as fh:
        for sets, fams in by.items():
            rec = {"sets": list(sets), "params_hash": PRM.resolve(None, list(sets)).hash, "seeds": [seeds[0], seeds[-1]],
                   "exploratory_fixed_keys": bool(a.allow_fixed), "families": {}}
            for fam, scores in fams.items():
                scores.sort(key=lambda s: s.seed)
                rec["families"][fam] = B.summarise_scores(scores)
                if fam == "pairs":
                    for kind in ("x_crossing", "t_contact", "parallel_xy", "stacked_z"):
                        sub = [s for s in scores if s.extra.get("kind") == kind]
                        rec["families"][f"pairs_{kind}"] = {
                            "n": len(sub),
                            "both_matched": float(sum(all(m["matched"] for m in s.matches) for s in sub) / max(len(sub), 1)),
                            "total_rel_err": float(sum((s.our_length - s.truth_length) / s.truth_length for s in sub)
                                                   / max(len(sub), 1))}
                rec["families"][fam]["elapsed_s_mean"] = float(sum(s.extra["elapsed_s"] for s in scores) / len(scores))
            fh.write(json.dumps(rec, default=lambda o: None if o != o else o) + "\n")
            d = rec["families"].get("dense", {})
            i = rec["families"].get("iso", {}).get("tilt_le60", {})
            g = rec["families"].get("gaps", {})
            print(f"{' '.join(sets) or 'defaults'}: dense {d.get('total_rel_err_mean', float('nan')):+.3f} "
                  f"sd {d.get('total_rel_err_sd', float('nan')):.3f} cov {d.get('coverage_mean', float('nan')):.3f} "
                  f"off {d.get('off_truth_frac_mean', float('nan')):.3f} | iso {i.get('len_err_mean', float('nan')):+.3f} "
                  f"sd {i.get('len_err_sd', float('nan')):.3f} one {i.get('one_strand_frac', float('nan')):.2f} | "
                  f"gaps rec {g.get('gap_recall', float('nan')):.2f} prec {g.get('gap_precision', float('nan')):.2f}",
                  flush=True)
    print(f"done in {time.time() - t0:.0f} s -> {out / 'calibration.jsonl'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
