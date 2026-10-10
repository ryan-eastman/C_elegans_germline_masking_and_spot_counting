# The test suite

How the suite is laid out, what the markers mean, and how to run it on the workstation, on a Mac
and on GitHub Actions. Everything here runs on the CPU and needs no microscopy data: the synthetic
gonads are built inside the fixtures.

## What is in here

- One module per stage or concern (`test_axis.py`, `test_coloc.py`, `test_granule.py`,
  `test_partition.py`, `test_sc.py`, `test_staging.py`, `test_spots.py`, ...), plus the config and
  profile machinery (`test_profiles.py`), the stage registry (`test_stages.py`), batch bookkeeping
  (`test_batch.py`) and the model download command (`test_models_registry.py`, which fakes the
  network).
- `test_imaris_truth.py`: the Imaris truth reader (`germquant imaris-truth`) on a tiny synthetic
  HDF5 written in the `.ims` layout; skips itself when `h5py` is not installed.
- `test_end_to_end.py`: the whole single-image pipeline on synthetic 3-channel and 4-channel gonads
  with the `.nd2` reader monkeypatched. Its fixtures (`patched_reader`, `patched_reader_4ch`) are
  imported by several other modules.
- `test_timelapse_reader.py`: the T-aware reader (`docs/TIMELAPSE.md`). A fake `ND2File` with nd2's
  one-chunk-per-frame layout records which planes a read touches, so the tests show that `t=<index>`
  reads one timepoint's planes and nothing else, that `t=None` is the old behaviour, that out-of-range
  values are refused before any pixel is read, and that the header reports `n_timepoints`, `period_s`
  and `t_seconds`. The bioio branch runs on a fake `BioImage` and, when `bioio-ome-tiff` is installed,
  on an OME-TIFF with a T axis written by tifffile. `process_image(t=...)` is checked to suffix the
  image id, record the timepoint in the markers only, and write the same numbers as a fixed run on the
  same stack.
- `test_golden.py`: the numeric gate. `test_synthetic_golden` runs the same two synthetic gonads
  with the classical nucleus segmenter and SpotMAX and compares every output table cell for cell
  with the references under `golden/synthetic/` and `golden/synthetic_4ch/`. A `config_hash.txt`
  sentinel next to each reference pins the text of `config/config.yaml` and the channel map, so an
  edit to either fails the test even when no number moved. The output schema is append-only: a new
  column is allowed, a renamed, removed or reordered column is not. `test_real_golden` is the opt-in
  real-data comparison (see the environment variables below). The goldens were produced on the
  workstation (Windows, the pinned `constraints/desktop-2026-09.txt` environment) and are checked
  with `check_exact=True`.
- `test_benchmark_split.py`: the pre-registered RAD-51 benchmark split (`golden/rad51_split.json`,
  `docs/BENCHMARK_PROTOCOL.md`): the 12 tags and their train/val/test roles, the sha256 helper, the
  Imaris creation-parameter reader on a synthetic HDF5 project, and the manifest builder
  (`scripts/benchmark_manifest.py`) on a temporary tree.
- `test_gpu.py`: the real Cellpose-SAM 3D path. See the markers.
- `test_benchmark.py`: the benchmark harness (`germquant.benchmark`, see `docs/BENCHMARK.md`) on
  synthetic label volumes and spot sets with known answers for every metric (a split, a merge, an
  out-of-coverage prediction, the z-extent error, the chunked Hungarian against the dense result, the
  bootstrap intervals), plus both benchmark scripts run on a temporary manifest in a subprocess.

## Markers

Registered in `pyproject.toml` (`[tool.pytest.ini_options] markers`); CI passes `--strict-markers`,
so a misspelt marker is an error rather than a silently empty selection.

- `gpu`: needs a CUDA or Apple Metal accelerator. `tests/test_gpu.py` skips itself when none is
  present, and `pytest.importorskip("torch")` / `pytest.importorskip("cellpose")` at the top of the
  module turn a missing package into a skip instead of a collection error. On a machine WITH a GPU
  the test runs and occupies the card for a few seconds. To keep a run on the CPU (for example while
  a batch job owns the GPU) set `CUDA_VISIBLE_DEVICES=-1`; an empty value does not work on Windows,
  which drops empty environment variables.
- `model`: needs large model weights, either the 1.2 GB nucleus model that `germquant get-model`
  downloads or the Cellpose-SAM weights that Cellpose fetches on first use. CI deselects it with
  `-m "not model"`. Today only `test_gpu.py` carries it; any future test that loads the fine-tuned
  model must carry it too, so CI never downloads a gigabyte.

Opt-ins driven by environment variables rather than markers:

- `GERMQUANT_GOLDEN=<golden tree>` and `GERMQUANT_CANDIDATE=<run tree>`: `test_real_golden` runs
  `scripts/regression_diff.py` on two `germquant run` output trees and fails on any difference.
- `GERMQUANT_UPDATE_GOLDEN=1`: `test_synthetic_golden` rewrites the synthetic references instead of
  checking them. Only for a deliberately number-changing step, run from the code that is being
  certified; see the docstring of `test_golden.py` and `scripts/make_goldens.py`.

Other automatic skips: the extended-path check in `test_fsutil.py` runs only on Windows, its
POSIX twin only elsewhere, `test_staging.py` skips its comparison with the August staging tables
when they are not on the machine, and `test_benchmark_split.py` skips the check that the 24
registered `.nd2` and `.ims` files exist with their registered sizes when the N2 RAD-51 data is not
on the machine.

## Running the suite locally

Workstation (Windows, the pinned `.venv`), from the repository root:

    .venv\Scripts\python.exe -m pytest -q

PowerShell, keeping the GPU out of it:

    $env:CUDA_VISIBLE_DEVICES = "-1"
    .venv\Scripts\python.exe -m pytest -q

Mac or Linux, in the activated environment: `pytest -q`. The whole suite is 136 tests and takes
about 75 s on the workstation (134 pass, 2 skip, or 3 skip with the GPU hidden). Useful selections:

    pytest tests/test_golden.py -q        # just the numeric gate
    pytest -q -m "not model"              # what CI runs
    pytest -q -m gpu                      # the accelerator path, on a machine that has one

From a git worktree whose editable install still points at the main checkout, put the worktree's
`src` first so the tests import the code under test:

    $env:PYTHONPATH = "<worktree>\src"
    python -c "import germquant; print(germquant.__file__)"    # must print the worktree path
    python -m pytest tests -q -p no:cacheprovider

Lint, the same command CI runs: `ruff check src tests scripts` (`uvx ruff check src tests scripts`
works without installing ruff). The rule set is written out in `pyproject.toml` under
`[tool.ruff.lint]` so a ruff release cannot move the gate.

## Continuous integration

`.github/workflows/ci.yml` runs on every push to `main` and to `feat/**` branches, on pull requests
and on demand (the "Run workflow" button). One job per OS, `ubuntu-latest` and `windows-latest`,
with Python 3.11:

1. `uv` installs torch and torchvision as CPU wheels from the PyTorch CPU index, then germquant with
   the `gpu`, `spots`, `sc` and `dev` extras (Cellpose, SpotMAX and Cell-ACDC, skan, pytest and
   ruff) against `constraints/ci-cpu-2026-10.txt`: the validated workstation versions minus the
   CUDA-only packages. The uv cache directory is saved between runs.
2. `ruff check src tests scripts`.
3. `pytest -q -rs -m "not model" --strict-markers`: the full suite including the synthetic golden
   tests; `test_gpu.py` reports itself as skipped.

`QT_QPA_PLATFORM=offscreen` and `CUDA_VISIBLE_DEVICES=-1` are set for the whole job. So is `GERMQUANT_GOLDEN_RTOL=1e-12`: the synthetic goldens are compared exactly on the
workstation that made them, but the GitHub runners have other CPUs and sum floats in a different vector order,
which moved `coloc.pearson_r` by one part in 1e14 on the first CI run; the tolerance allows that and nothing
more (leave the variable unset locally so any real change fails). The README
badge shows the result of the latest run on the default branch.

Rehearsal on the workstation (2026-10-09, the same commands in a fresh `uv` environment with a cold
cache): torch and torchvision 39 s, the rest of the packages 31 s, ruff clean, pytest about 120 s.
Three and a half minutes in all, well inside the budget; a GitHub runner is slower, but not by
enough to matter. That rehearsal also found that `load_config` used to read the config and
channel-map files with the machine's locale codec, so the `config_hash` (and the golden sentinel)
came out differently in UTF-8 mode, on Linux and on macOS than on the workstation (cp1252). The
codec is now pinned to the workstation's (`fsutil.CONFIG_TEXT_ENCODING`), and the golden tests pass
with and without `PYTHONUTF8=1`. Outputs made before the pin on a Mac, on Linux or in UTF-8 mode carry the
other hash for the same file (config.yaml: d20d4321faad instead of 84bdda5a5230); treat them as a different config version.

What CI does and does not prove: it shows that a clean machine can install the pinned environment
and that the numbers are reproduced exactly on it. The references were made on the workstation, so
a last-digit difference on another CPU or operating system would surface as a golden failure. That
is a finding, not noise: first check whether a pin moved, then decide whether the difference is a
platform effect that deserves a documented tolerance. Until the first green run on `ubuntu-latest`,
Linux is unrehearsed: the install was resolved for Linux with `uv pip compile`, but the suite itself
has only been run on Windows.

## The real golden image sets

Two sets of real gonads are frozen under `C:/Users/ryane/germquant_golden/8666f07/` (tags in
`tests/golden/real_manifest.json`):

* `n2dryice_*` (7 gonads plus 4 stride-4 variants): the 4-channel N2 set of 2026-06-25 (DAPI, SYP-3, PGL-1,
  lamin). That imaging session was a fixation test (dry-ice freeze-crack), not the lab's standard preparation.
  It is a regression input only: it shows whether a code change moves the numbers on a fixed input. Do not use
  it as biological reference data, as training or benchmark data, or as representative images in a paper.
* `n2rad51_*` (12 gonads): the 3-channel N2 RAD-51 cross-validation set of 2025-11-05 (DAPI, SYP, RAD-51),
  no heat shock and heat shock, hermaphrodites and males, imaged with the standard preparation and paired
  with Imaris annotations. These are both regression inputs and the benchmark set for spot detection.
  Their train/val/test split, with the size and sha256 of every `.nd2` and `.ims`, is
  `tests/golden/rad51_split.json` (rebuilt by `scripts/make_rad51_split.py`); the benchmark protocol is
  `docs/BENCHMARK_PROTOCOL.md`.

