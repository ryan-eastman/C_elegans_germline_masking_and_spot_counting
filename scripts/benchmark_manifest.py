"""Build benchmarks/<name>/manifest.json: the pre-registered split joined with what exists on disk.

    python scripts/benchmark_manifest.py --truth-root <dir> [--name n2rad51]
        [--split tests/golden/rad51_split.json] [--results-root C:/Users/ryane/germquant_golden/8666f07]
        [--out benchmarks] [--truth-file masks={stem}__imaris_labels.tif ...]
        [--additions tests/golden/rad51_additions.json --additions-results-root DIR
         --truth-check CSV --spot-truth-decisions tests/golden/rad51_additions_decisions.json]
        [--manifest-file NAME.json]

One entry per gonad: the split role, condition and sex, the .nd2 and .ims (exists, size matches the
registered size), the pipeline output files under <results-root>/<tag>/ (nuclei table, nuclei_labels.tif,
spots table, image summary, run manifest) and the Imaris truth files under <truth-root>
(<stem>__imaris_<kind>), each with its path and an exists flag, plus outputs_complete, truth_complete and
spot_truth_usable (false for the gonads whose spots were exported with Imaris region growing). Nothing is
computed from the data: the manifest is the input list for the scoring step and a record of what was
available when it was built (created_utc, and the sha256 of the split file it was built from).

--truth-file overrides one expected truth file name (kind=pattern, {stem} is the image id); files the
export wrote under any other name still appear in each entry's "found" list.

--additions appends the gonads of the additions record of protocol changelog entry 9 after the registered
ones, which stay exactly as they are (their pipeline outputs under --additions-results-root/<tag>/, default
the results root; their truth files under the truth root). An addition's spot_truth_usable follows
benchmark_split.addition_spot_truth_usable (protocol changelog entries 9 (d) and 10 (c)): the export's own
checks (the frame check of the truth params, the xlsx verification), then the decision of a changelog entry
(--spot-truth-decisions, the committed decisions file), else the section 4.4 agreement of --truth-check (the
--out table of scripts/check_truth_assignment.py) of at least 99.35 percent. Without either option no
addition sourced from its .ims is usable. Both files are recorded with their sha256 (spot_truth_inputs).
--manifest-file names the file written in <out>/<name>/ (default manifest.json): a revised
manifest goes next to the registered one under another name, so trials.jsonl and final_runs.jsonl, which the
scoring scripts find next to the manifest, stay the same files.
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
    ap.add_argument("--additions", default=None, help="the additions record of protocol changelog entry 9")
    ap.add_argument("--additions-results-root", default=None,
                    help="pipeline outputs of the additions, one folder per tag (default: --results-root)")
    ap.add_argument("--truth-check", default=None, metavar="CSV",
                    help="with --additions: the section 4.4 table of scripts/check_truth_assignment.py --out")
    ap.add_argument("--spot-truth-decisions", default=None, metavar="JSON",
                    help="with --additions: the decisions file of a changelog entry on the additions' spot "
                         "truth (tests/golden/rad51_additions_decisions.json, entry 10)")
    ap.add_argument("--manifest-file", default="manifest.json",
                    help="file name of the manifest in <out>/<name>/ (default manifest.json)")
    args = ap.parse_args(argv)

    overrides = {}
    for item in args.truth_file:
        kind, sep, pattern = item.partition("=")
        if not sep or not kind or not pattern:
            ap.error(f"--truth-file expects KIND=PATTERN, got {item!r}")
        overrides[kind] = pattern

    if (args.truth_check or args.spot_truth_decisions) and not args.additions:
        ap.error("--truth-check and --spot-truth-decisions need --additions")
    split = bs.load_split(args.split)
    additions = bs.load_additions(args.additions) if args.additions else None
    manifest = bs.build_manifest(split, args.results_root, args.truth_root, name=args.name,
                                 truth_files=overrides or None, split_path=args.split, additions=additions,
                                 additions_results_root=args.additions_results_root, additions_path=args.additions,
                                 truth_check=args.truth_check, spot_truth_decisions=args.spot_truth_decisions)
    out = bs.write_manifest(manifest, args.out, args.manifest_file)

    print(f"{'tag':28s} {'role':5s} {'outputs':8s} {'truth':8s} {'spot truth':10s}")
    for e in manifest["gonads"]:
        outputs = "complete" if e["outputs_complete"] else "missing"
        truth = "complete" if e["truth_complete"] else ("partial" if e["truth"]["found"] else "missing")
        usable = "usable" if e["spot_truth_usable"] else "excluded"
        reason = f"  {e['spot_truth_reason']}" if e.get("addition") else ""
        print(f"{e['tag']:28s} {e['role']:5s} {outputs:8s} {truth:8s} {usable:10s}{reason}")
    s = manifest["summary"]
    print(f"wrote {out}: {s['n_gonads']} gonads, outputs complete for {len(s['outputs_complete'])}, "
          f"truth complete for {len(s['truth_complete'])}, spot truth usable for {len(s['spot_truth_usable'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
