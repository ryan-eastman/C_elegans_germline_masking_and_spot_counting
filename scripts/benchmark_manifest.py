"""Build benchmarks/<name>/manifest.json: the pre-registered split joined with what exists on disk.

    python scripts/benchmark_manifest.py --truth-root <dir> [--name n2rad51]
        [--split tests/golden/rad51_split.json] [--results-root C:/Users/ryane/germquant_golden/8666f07]
        [--out benchmarks] [--truth-file masks={stem}__imaris_labels.tif ...]

One entry per gonad: the split role, condition and sex, the .nd2 and .ims (exists, size matches the
registered size), the pipeline output files under <results-root>/<tag>/ (nuclei table, nuclei_labels.tif,
spots table, image summary, run manifest) and the Imaris truth files under <truth-root>
(<stem>__imaris_<kind>), each with its path and an exists flag, plus outputs_complete, truth_complete and
spot_truth_usable (false for the gonads whose spots were exported with Imaris region growing). Nothing is
computed from the data: the manifest is the input list for the scoring step and a record of what was
available when it was built (created_utc, and the sha256 of the split file it was built from).

--truth-file overrides one expected truth file name (kind=pattern, {stem} is the image id); files the
export wrote under any other name still appear in each entry's "found" list.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from germquant.validate import benchmark_split as bs  # noqa: E402

DEFAULT_RESULTS_ROOT = "C:/Users/ryane/germquant_golden/8666f07"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--split", default=str(REPO / "tests" / "golden" / "rad51_split.json"))
    ap.add_argument("--results-root", default=DEFAULT_RESULTS_ROOT, help="pipeline outputs, one folder per tag")
    ap.add_argument("--truth-root", required=True, help="folder with the <stem>__imaris_* truth files")
    ap.add_argument("--name", help="benchmark name (default: the split's name)")
    ap.add_argument("--out", default=str(REPO / "benchmarks"), help="root for <name>/manifest.json")
    ap.add_argument("--truth-file", action="append", default=[], metavar="KIND=PATTERN",
                    help="override an expected truth file name, e.g. masks={stem}__imaris_labels.tif")
    args = ap.parse_args(argv)

    overrides = {}
    for item in args.truth_file:
        kind, sep, pattern = item.partition("=")
        if not sep or not kind or not pattern:
            ap.error(f"--truth-file expects KIND=PATTERN, got {item!r}")
        overrides[kind] = pattern

    split = bs.load_split(args.split)
    manifest = bs.build_manifest(split, args.results_root, args.truth_root, name=args.name,
                                 truth_files=overrides or None, split_path=args.split)
    out = bs.write_manifest(manifest, args.out)

    print(f"{'tag':28s} {'role':5s} {'outputs':8s} {'truth':8s} {'spot truth':10s}")
    for e in manifest["gonads"]:
        outputs = "complete" if e["outputs_complete"] else "missing"
        truth = "complete" if e["truth_complete"] else ("partial" if e["truth"]["found"] else "missing")
        usable = "usable" if e["spot_truth_usable"] else "excluded"
        print(f"{e['tag']:28s} {e['role']:5s} {outputs:8s} {truth:8s} {usable:10s}")
    s = manifest["summary"]
    print(f"wrote {out}: {s['n_gonads']} gonads, outputs complete for {len(s['outputs_complete'])}, "
          f"truth complete for {len(s['truth_complete'])}, spot truth usable for {len(s['spot_truth_usable'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
