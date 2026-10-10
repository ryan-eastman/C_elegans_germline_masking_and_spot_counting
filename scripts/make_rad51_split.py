"""Build tests/golden/rad51_split.json, the pre-registered split record (docs/BENCHMARK_PROTOCOL.md).

    python scripts/make_rad51_split.py [--ims-root C:/Users/ryane/n2_rad51_local/imaris]
        [--hash-cache tests/golden/rad51_split.json] [--out tests/golden/rad51_split.json]
        [--no-metadata]

Reads the n2rad51 entries of tests/golden/real_manifest.json for the .nd2 paths, pairs each with its
Imaris project (<ims-root>/<noHS_ims2|HS_ims2>/<stem>.ims), hashes both files (sha256; about a minute
for the 24 files when they are in the disk cache, several minutes cold), reads the .nd2 sizes and voxel
size and the Imaris spot and surface counts, Quality threshold and region-growing flag from the .ims,
and writes the split with the roles fixed in germquant.validate.benchmark_split.PREREGISTERED_ROLES.

The roles never change. Re-running this script only refreshes the file facts; a changed size or hash
means the data changed, which the protocol's changelog must record before any result is computed on
the changed file.

--hash-cache reuses the digests of an earlier record (a split file, or a plain JSON mapping
{path: {"size_bytes": n, "sha256": hex}}) for files whose size has not changed, so a re-run after a
note or metadata change does not re-read 50 GB. --no-metadata skips the .nd2 and .ims reads (hashes
only; region_growing is then null).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from germquant.validate import benchmark_split as bs  # noqa: E402

DEFAULT_IMS_ROOT = "C:/Users/ryane/n2_rad51_local/imaris"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--real-manifest", default=str(REPO / "tests" / "golden" / "real_manifest.json"))
    ap.add_argument("--ims-root", default=DEFAULT_IMS_ROOT, help="folder holding noHS_ims2/ and HS_ims2/")
    ap.add_argument("--hash-cache", help="earlier split file or {path: {size_bytes, sha256}} JSON to reuse digests from")
    ap.add_argument("--out", default=str(REPO / "tests" / "golden" / "rad51_split.json"))
    ap.add_argument("--registered-on", default=bs.REGISTERED_ON, help="registration date written into the file")
    ap.add_argument("--no-metadata", action="store_true", help="hashes only; skip the .nd2 and .ims reads")
    args = ap.parse_args(argv)

    cache = bs.load_hash_cache(args.hash_cache) if args.hash_cache else None
    split = bs.build_split(args.real_manifest, args.ims_root, hash_cache=cache,
                           registered_on=args.registered_on, read_metadata=not args.no_metadata)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(split, indent=1) + "\n", encoding="utf-8")

    print(f"{'tag':28s} {'role':5s} {'nd2 GB':>7s} {'ims GB':>7s} {'spots':>7s} {'surf':>5s} {'RG':>5s}")
    for g in split["gonads"]:
        ims = g["ims"]
        rg = ims.get("region_growing")
        print(f"{g['tag']:28s} {g['role']:5s} {g['nd2']['size_bytes'] / 1e9:7.3f} {ims['size_bytes'] / 1e9:7.3f} "
              f"{str(ims.get('n_spots', '')):>7s} {str(ims.get('n_surfaces', '')):>5s} "
              f"{'true' if rg is True else 'false' if rg is False else 'null':>5s}")
    print(f"wrote {out} ({len(split['gonads'])} gonads; roles "
          + ", ".join(f"{r}={len(v)}" for r, v in split["roles"].items()) + ")")
    return 0


if __name__ == "__main__":
    sys.exit(main())
