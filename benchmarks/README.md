# Benchmarks

One folder per benchmark, `benchmarks/<name>/`, holding the manifest the scoring step reads and,
later, the results written against it. The protocol that governs the RAD-51 benchmark is
`docs/BENCHMARK_PROTOCOL.md`; the split it scores is `tests/golden/rad51_split.json`.

`manifest.json` is rebuilt from the split and from what exists on disk, so it is not committed (it
carries workstation paths, exists flags and a timestamp):

    python scripts/benchmark_manifest.py --truth-root <folder with the <stem>__imaris_* files>
        [--results-root C:/Users/ryane/germquant_golden/8666f07] [--name n2rad51] [--out benchmarks]

Each entry pairs one gonad's role in the split with its `.nd2`, its `.ims`, the pipeline output files
under `<results-root>/<tag>/` and the Imaris truth files under `--truth-root`, each with an exists
flag, plus `outputs_complete`, `truth_complete` and `spot_truth_usable` (false for the two heat-shock
male gonads whose spots were exported with Imaris region growing). The split file's sha256 is recorded
in the manifest header, so a result can always be traced to the exact split it was scored against.
