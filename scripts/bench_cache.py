"""Build the detector bench cache: the exact inputs of the golden runs' spots stage, one folder per gonad.

    python scripts/bench_cache.py [--manifest C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/manifest.json]
        [--split-file tests/golden/rad51_split.json] [--cache-root C:/Users/ryane/germquant_golden/bench_cache]
        [--tags TAG ...] [--z-chunk 8] [--verify-sha256] [--overwrite]

For every manifest entry whose results folder holds the golden pipeline outputs (nuclei table, label
image, spots table, image summary, run manifest) it writes <cache-root>/<tag>/ with foci.npy and dna.npy
(the two channels at full resolution, read with germquant's own reader; channel indices from the golden
run's channel_map_json), labels.npy (int32, exactly the label image the golden run handed to
detect_spots: the full nucleus segmentation, which the run also saved as <image>__nuclei_labels.tif),
nuclei.csv (nucleus_id, in_germline, centroids, golden n_spots) and meta.json (written last; see
germquant.benchmark.detector for every field). Arrays are written through numpy's open_memmap a few z
planes at a time. Entries without outputs are skipped and listed; an existing complete cache is kept
unless --overwrite. The test gonads are cached like the others: the cache holds images and labels,
no metric, and scoring a test gonad stays behind the guard of bench_detector.py.

At the end it prints the size on disk per gonad and in total, and writes <cache-root>/cache_index.json.
Exit code 0 unless a gonad failed with an error (a skip is not an error).
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from germquant.benchmark import detector as bd  # noqa: E402
from germquant.benchmark.manifest import load_manifest  # noqa: E402


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", default=str(bd.DEFAULT_MANIFEST), help="benchmark manifest JSON")
    ap.add_argument("--split-file", default=str(REPO / "tests" / "golden" / "rad51_split.json"),
                    help="the pre-registered split record (supplies the .nd2 path, size and sha256)")
    ap.add_argument("--cache-root", default=str(bd.DEFAULT_CACHE_ROOT), help="cache root folder")
    ap.add_argument("--tags", nargs="*", default=None, help="only these tags")
    ap.add_argument("--z-chunk", type=int, default=bd.DEFAULT_Z_CHUNK, help="z planes per write")
    ap.add_argument("--verify-sha256", action="store_true",
                    help="hash every .nd2 and refuse one that is not the registered file")
    ap.add_argument("--overwrite", action="store_true", help="rebuild caches that already exist")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    entries = load_manifest(args.manifest)
    if args.tags:
        known = {e.tag for e in entries}
        unknown = [t for t in args.tags if t not in known]
        if unknown:
            print(f"unknown tags (not in the manifest): {unknown}")
            return 2
        entries = [e for e in entries if e.tag in set(args.tags)]
    split = None
    if args.split_file and Path(args.split_file).is_file():
        split = json.loads(Path(args.split_file).read_text(encoding="utf-8"))
    else:
        print(f"note: no split record at {args.split_file}; .nd2 paths come from the manifest, no sha256")
    root = Path(args.cache_root)
    root.mkdir(parents=True, exist_ok=True)
    results, failed = [], []
    for e in entries:
        t0 = time.time()
        try:
            r = bd.cache_gonad(e, root, split_record=split, z_chunk=args.z_chunk,
                               verify_sha256=args.verify_sha256, overwrite=args.overwrite)
        except Exception as err:  # noqa: BLE001 - one bad gonad must not stop the others; reported below
            traceback.print_exc()
            r = {"tag": e.tag, "status": "failed", "reason": f"{type(err).__name__}: {err}"}
            failed.append(e.tag)
        r["elapsed_s"] = round(time.time() - t0, 1)
        results.append(r)
        extra = f" ({r['reason']})" if r.get("reason") else ""
        size = f" {r['bytes'] / 1e9:.2f} GB" if r.get("bytes") else ""
        print(f"  {e.tag}: {r['status']}{size}{extra} [{r['elapsed_s']} s]", flush=True)

    sizes = bd.cache_size_bytes(root)
    index = {"created_utc": bd._utc_now(), "manifest": args.manifest, "split_file": args.split_file,
             "cache_format": bd.CACHE_FORMAT, "this_run": results, "complete": sizes["per_tag"],
             "total_bytes": sizes["total"]}
    (root / "cache_index.json").write_text(json.dumps(index, indent=1) + "\n", encoding="utf-8")
    skipped = [r for r in results if r["status"] == "skipped"]
    print(f"cache root {root}: {len(sizes['per_tag'])} complete gonad caches, "
          f"{sizes['total'] / 1e9:.2f} GB on disk ({sizes['total']} bytes)")
    for tag, b in sizes["per_tag"].items():
        print(f"  {tag}: {b / 1e9:.3f} GB")
    if skipped:
        print("skipped: " + "; ".join(f"{r['tag']} ({r['reason']})" for r in skipped))
    if failed:
        print(f"FAILED: {failed}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
