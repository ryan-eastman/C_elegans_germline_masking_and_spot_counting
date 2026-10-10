"""Check the detector bench harness: SpotMAX with the golden run's own spots config, run from the bench
cache, must reproduce the golden spots table exactly.

    python scripts/bench_reproduce_golden.py --tags n2rad51_noHS_herm_000 n2rad51_HS_herm_000
        [--manifest PATH] [--cache-root DIR] [--out report.json] [--backend spotmax]

For every tag: germquant.benchmark.detector.run_backend(tag, "spotmax", {}) (the cached foci channel,
label image and voxel size, the golden spots: section) and a column by column comparison with the golden
<image>__spots.parquet (exact binary floats; the CSV, parsed with round-trip precision, when there is no
parquet): spot_id, nucleus_id, marker, z_um, y_um, x_um, effect_size, spot_mean_intensity, detector, in
row order, plus the per-nucleus counts against the golden nuclei table's n_spots; the same comparison
against the CSV is reported next to it. Prints the result per tag and writes a JSON report with --out. No truth file is read: this is a reproducibility check, not a
score, and it refuses test gonads all the same. Exit code 0 when every tag is identical.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import pandas as pd  # noqa: E402

from germquant.benchmark import detector as bd  # noqa: E402
from germquant.benchmark.manifest import load_manifest  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tags", nargs="+", required=True)
    ap.add_argument("--manifest", default=str(bd.DEFAULT_MANIFEST))
    ap.add_argument("--cache-root", default=str(bd.DEFAULT_CACHE_ROOT))
    ap.add_argument("--backend", default="spotmax", help="the backend the golden run used")
    ap.add_argument("--out", default=None, help="write the report as JSON here")
    args = ap.parse_args(argv)
    os.environ.setdefault("TQDM_DISABLE", "1")
    entries = {e.tag: e for e in load_manifest(args.manifest)}
    roles = {t: bd.entry_role(e) for t, e in entries.items()}
    try:
        bd.guard_test(args.tags, roles=roles)
    except bd.TestSetGuardError as e:
        print(f"REFUSED: {e}")
        return 2
    report, all_ok = {}, True
    for tag in args.tags:
        t0 = time.time()
        cache = bd.load_cache(tag, args.cache_root, arrays=("foci", "labels"), mmap=False)
        per_spot, per_nuc = bd.run_backend(tag, args.backend, {}, cache=cache)
        golden_dir = cache.meta["golden"]["results_dir"]
        gp = bd.golden_paths(golden_dir, cache.meta["stem"])
        gnuc = pd.read_csv(gp["nuclei"])
        res = bd.compare_with_golden(per_spot, per_nuc, bd.read_golden_spots(golden_dir, cache.meta["stem"]), gnuc)
        res["golden_source"] = "parquet" if gp["spots_parquet"].is_file() else "csv"
        # the CSV as well (parsed with round-trip precision, so it holds the exact floats pandas wrote)
        res_csv = bd.compare_with_golden(per_spot, per_nuc, pd.read_csv(gp["spots"], float_precision="round_trip"),
                                         gnuc)
        res["csv"] = {"identical": res_csv["identical"], "columns": res_csv["columns"]}
        res["identical"] = bool(res["identical"] and res_csv["identical"])
        res["elapsed_s"] = round(time.time() - t0, 1)
        report[tag] = res
        all_ok = all_ok and res["identical"]
        diffs = {c: r["n_different"] for c, r in res["columns"].items() if r.get("compared") and r["n_different"]}
        print(f"{tag}: {'IDENTICAL' if res['identical'] else 'DIFFERENT'} ({res['n_ours']} spots, golden "
              f"{res['n_golden']}; per-nucleus counts differing: {res.get('nuclei', {}).get('n_count_different')}; "
              f"columns differing ({res['golden_source']}): {diffs or 'none'}; CSV identical: "
              f"{res_csv['identical']}) [{res['elapsed_s']} s]", flush=True)
        del cache, per_spot, per_nuc
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
