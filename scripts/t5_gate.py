"""Tranche 5 gate: the checks a T5 package must pass, ending in one PASS/FAIL table for the PR description.

    python scripts/t5_gate.py [--base REF] [--no-workstation] [--require-golden] [--skip STEP ...]
                              [--markdown OUT.md] [--json OUT.json]

docs/DESIGN_system_packs.md, sections 11.1 to 11.3 (gates G1 to G5). The steps, in order:

  ruff         ruff check src tests scripts with the ruff version CI installs (G1)
  pytest       python -m pytest -m "not model" on this checkout's src: the unit suite, the synthetic goldens
               and the fence, tests/test_default_invariants.py (G1, G2, G3)
  scope        scripts/t5_scope_check.py: no frozen path modified or deleted, no existing pin changed (G3, G5)
  golden       scripts/regression_diff.py $GERMQUANT_GOLDEN $GERMQUANT_CANDIDATE --strict-cols --strict-ids:
               the real goldens exact, no added column, no renumbered label image (G4)
  filelist     both roots hold the same images and, per image, the same file names (regression_diff walks
               only the golden's files and only notes extra images; design 11.1 item 2)
  workstation  tests/test_default_invariants.py with GERMQUANT_WORKSTATION_PINS=1: the pins that may differ
               between platforms (scope "workstation" in the pins file, today the montage pixels), enforced
               on the lab workstation (G2)
  legacy       scripts/regression_diff.py $GERMQUANT_GOLDEN_LEGACY $GERMQUANT_CANDIDATE: the comparison with
               the older golden tree C:/Users/ryane/germquant_golden/8666f07 stays at 0 failing
               comparisons (design 11.1 item 2, G4)

golden and filelist run when GERMQUANT_GOLDEN and GERMQUANT_CANDIDATE are set, legacy when
GERMQUANT_GOLDEN_LEGACY and GERMQUANT_CANDIDATE are set. Setting only one of GERMQUANT_GOLDEN and
GERMQUANT_CANDIDATE, or GERMQUANT_GOLDEN_LEGACY without them, is an error (exit 2), so a misspelt variable
cannot turn the real-golden comparison into a quiet SKIP. workstation always runs (about 10 seconds), because
G2 counts the montage pin for every package; --no-workstation skips it on a machine other than the lab
workstation. A step that does not run is SKIP; with --require-golden (a G4 run) a skipped golden,
filelist, workstation or legacy step is a FAIL. The candidate tree is made beforehand with
scripts/make_goldens.py (hours on the RTX 5090): this gate only compares trees and never runs the real
goldens.

The test steps run with PYTHONPATH=<this checkout>/src (so a worktree tests its own code), with
CUDA_VISIBLE_DEVICES=-1 unless it is set, and without GERMQUANT_GOLDEN_RTOL or GERMQUANT_UPDATE_GOLDEN, so
the synthetic goldens are compared exactly and never rewritten.

Verdict and exit code: PASS (0) when ruff, pytest and scope ran and no step failed; FAIL (1) when a step
failed; INCOMPLETE (3) when nothing failed but ruff, pytest or scope was skipped with --skip; 2 when the
golden variables are set only in part.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CI_CONSTRAINTS = REPO / "constraints" / "ci-cpu-2026-10.txt"     # the file .github/workflows/ci.yml installs
CHECKS = {
    "ruff": "ruff check src tests scripts",
    "pytest": 'pytest -m "not model" (suite, synthetic goldens, fence)',
    "scope": "t5_scope_check (R1, R8, 11.2, G3)",
    "golden": "real goldens, --strict-cols --strict-ids",
    "filelist": "same images and file names",
    "workstation": "workstation-only pins (GERMQUANT_WORKSTATION_PINS=1)",
    "legacy": "older golden tree, 0 failing comparisons",
}
STEPS = tuple(CHECKS)
CORE_STEPS = ("ruff", "pytest", "scope")                 # skipping one of these makes the run INCOMPLETE
GOLDEN_STEPS = ("golden", "filelist", "workstation", "legacy")   # a SKIP here is a FAIL with --require-golden
_CLEARED = ("GERMQUANT_GOLDEN", "GERMQUANT_CANDIDATE", "GERMQUANT_GOLDEN_LEGACY", "GERMQUANT_GOLDEN_RTOL",
            "GERMQUANT_UPDATE_GOLDEN", "GERMQUANT_WORKSTATION_PINS")


@dataclass
class Result:
    step: str
    status: str            # PASS | FAIL | SKIP
    detail: str
    seconds: float = 0.0

    @property
    def check(self) -> str:
        return CHECKS[self.step]


def _env(**extra: str) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if k not in _CLEARED}
    src = str(REPO / "src")
    env["PYTHONPATH"] = src + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env.setdefault("CUDA_VISIBLE_DEVICES", "-1")
    env["PYTHONIOENCODING"] = "utf-8"
    env.update(extra)
    return env


def _run(cmd: list[str], env: dict[str, str] | None = None) -> tuple[int, list[str]]:
    """Run ``cmd`` in the checkout, echoing its output; returns (exit code, output lines)."""
    print(f"\n$ {' '.join(cmd)}", flush=True)
    proc = subprocess.Popen(cmd, cwd=REPO, env=env or _env(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, encoding="utf-8", errors="replace")
    lines: list[str] = []
    assert proc.stdout is not None
    for line in proc.stdout:
        sys.stdout.write(line)
        lines.append(line.rstrip("\n"))
    return proc.wait(), lines


def _last(lines: list[str], pattern: str, default: str) -> str:
    rx = re.compile(pattern)
    hits = [ln.strip().strip("= ") for ln in lines if rx.search(ln)]
    return hits[-1] if hits else default


def _status(rc: int) -> str:
    return "PASS" if rc == 0 else "FAIL"


def ruff_pin() -> str | None:
    """The ruff version CI lints with (the ruff== line of the CI constraints file), or None."""
    try:
        text = CI_CONSTRAINTS.read_text(encoding="utf-8")
    except OSError:
        return None
    m = re.search(r"^ruff==(\S+)\s*$", text, re.M)
    return m.group(1) if m else None


def ruff_command() -> tuple[list[str], str]:
    """(command, what it is) for the CI version of ruff: a ruff on PATH only when it is that version,
    otherwise uvx running the pinned version, otherwise this interpreter's ruff module."""
    pin = ruff_pin()
    on_path = shutil.which("ruff")
    if on_path:
        r = subprocess.run([on_path, "--version"], capture_output=True, text=True)
        version = r.stdout.split()[-1] if r.returncode == 0 and r.stdout.split() else "?"
        if pin is None or version == pin:
            return [on_path], f"ruff {version} on PATH"
    if shutil.which("uvx"):
        return (["uvx", f"ruff@{pin}"], f"uvx ruff@{pin}") if pin else (["uvx", "ruff"], "uvx ruff, no pin found")
    return [sys.executable, "-m", "ruff"], "python -m ruff"


def step_ruff() -> Result:
    cmd, what = ruff_command()
    rc, lines = _run([*cmd, "check", "src", "tests", "scripts"])
    found = _last(lines, r"All checks passed|Found \d+ error", f"exit code {rc}")
    return Result("ruff", _status(rc), f"{found} ({what})")


def step_pytest() -> Result:
    rc, lines = _run([sys.executable, "-m", "pytest", "-m", "not model", "-q", "-rs", "-p", "no:cacheprovider",
                      "tests"])
    return Result("pytest", _status(rc), _last(lines, r"\d+ (passed|failed|errors?)\b", f"exit code {rc}"))


def step_scope(base: str) -> Result:
    with tempfile.TemporaryDirectory() as tmp:
        report = Path(tmp) / "scope.json"
        rc, lines = _run([sys.executable, str(REPO / "scripts" / "t5_scope_check.py"), "--base", base,
                          "--json", str(report)])
        rep = json.loads(report.read_text(encoding="utf-8")) if report.exists() else None
    if rep is None:
        return Result("scope", "FAIL", _last(lines, r"ERROR", f"exit code {rc}, no report"))
    detail = (f"base {rep['base'][:7]}, {rep['n_changed']} changed paths, {rep['fence']['pins_checked']} pins "
              f"checked, {len(rep['violations'])} frozen paths or pins touched")
    if rep["violations"]:
        detail += ": " + ", ".join(f"{v['path']} {v['change']}" for v in rep["violations"][:5])
    return Result("scope", _status(rc), detail)


def step_regression_diff(step: str, golden: str, candidate: str, strict: bool) -> Result:
    cmd = [sys.executable, str(REPO / "scripts" / "regression_diff.py"), golden, candidate]
    if strict:
        cmd += ["--strict-cols", "--strict-ids"]
    rc, lines = _run(cmd)
    return Result(step, _status(rc), _last(lines, r"golden images checked|no golden images|does not exist",
                                           f"exit code {rc}"))


def image_file_lists(root: Path) -> dict[str, set[str]]:
    """{<folder>/<image id>: its file names} for every folder that holds <image id>__*.csv files (the image
    keys of regression_diff.py, whose image id is the text before the first "__"), plus
    {<folder>/(other files): ...} for the folder's files without "__" (run_manifest.json) and subfolders."""
    out: dict[str, set[str]] = {}
    for folder in sorted({f.parent for f in root.rglob("*__*.csv")}):
        rel = folder.relative_to(root).as_posix()
        for f in folder.iterdir():
            name = f.name + ("/" if f.is_dir() else "")
            group = f.name.split("__")[0] if ("__" in f.name and f.is_file()) else "(other files)"
            out.setdefault(group if rel == "." else f"{rel}/{group}", set()).add(name)
    return out


def compare_file_lists(golden: Path, candidate: Path) -> list[str]:
    """Problems: images present on one side only, and per image the file names that differ."""
    g, c = image_file_lists(golden), image_file_lists(candidate)
    if not g:
        return [f"no images (*__*.csv) under {golden}"]
    problems = [f"missing from the candidate: {k}" for k in sorted(set(g) - set(c))]
    problems += [f"only in the candidate: {k}" for k in sorted(set(c) - set(g))]
    for key in sorted(set(g) & set(c)):
        lost, added = sorted(g[key] - c[key]), sorted(c[key] - g[key])
        if lost or added:
            problems.append(f"{key}: missing {lost}, added {added}")
    return problems


def step_filelist(golden: str, candidate: str) -> Result:
    gp, cp = Path(golden), Path(candidate)
    if not gp.is_dir() or not cp.is_dir():
        return Result("filelist", "FAIL", f"not a folder: {golden if not gp.is_dir() else candidate}")
    problems = compare_file_lists(gp, cp)
    for p in problems:
        print(f"  filelist: {p}")
    if problems:
        return Result("filelist", "FAIL", f"{len(problems)} differences, first: {problems[0]}")
    n = len([k for k in image_file_lists(gp) if not k.endswith("(other files)")])
    return Result("filelist", "PASS", f"{n} images, the same file names on both sides")


def step_workstation() -> Result:
    rc, lines = _run([sys.executable, "-m", "pytest", "tests/test_default_invariants.py", "-m", "not model",
                      "-q", "-rs", "-p", "no:cacheprovider"], env=_env(GERMQUANT_WORKSTATION_PINS="1"))
    return Result("workstation", _status(rc), _last(lines, r"\d+ (passed|failed|errors?)\b", f"exit code {rc}"))


def _header() -> str:
    def git(*args: str) -> str:
        r = subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else "?"

    dirty = "-dirty" if git("status", "--porcelain", "--untracked-files=no") else ""
    when = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    return (f"t5_gate: {git('rev-parse', '--short=7', 'HEAD')}{dirty} on {git('rev-parse', '--abbrev-ref', 'HEAD')}, "
            f"{when}, Python {platform.python_version()}, {platform.platform(terse=False)}")


def table(results: list[Result], header: str) -> str:
    rows = [header, "", "| step | check | result | detail | seconds |", "|---|---|---|---|---|"]
    rows += [f"| {r.step} | {r.check} | {r.status} | {r.detail.replace('|', '/')} | "
             f"{'' if r.status == 'SKIP' else r.seconds} |" for r in results]
    return "\n".join(rows)


def verdict(results: list[Result]) -> tuple[str, int]:
    """(verdict line, exit code): FAIL beats INCOMPLETE beats PASS."""
    failed = [r.step for r in results if r.status == "FAIL"]
    if failed:
        return f"FAIL ({', '.join(failed)})", 1
    skipped = [r.step for r in results if r.status == "SKIP" and r.step in CORE_STEPS]
    if skipped:
        return f"INCOMPLETE (skipped: {', '.join(skipped)})", 3
    return "PASS", 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default="auto", help="base for the scope check (default: auto)")
    ap.add_argument("--no-workstation", action="store_true",
                    help="skip the workstation-only pins (on a machine other than the lab workstation)")
    ap.add_argument("--require-golden", action="store_true",
                    help="count a skipped golden, filelist, workstation or legacy step as a FAIL")
    ap.add_argument("--skip", nargs="*", default=[], choices=STEPS, metavar="STEP",
                    help=f"steps to skip: {', '.join(STEPS)}")
    ap.add_argument("--markdown", help="also write the table and the verdict here")
    ap.add_argument("--json", help="also write the results here")
    args = ap.parse_args(argv)
    try:
        sys.stdout.reconfigure(errors="replace")
    except (AttributeError, ValueError):
        pass

    golden, candidate = os.environ.get("GERMQUANT_GOLDEN"), os.environ.get("GERMQUANT_CANDIDATE")
    legacy = os.environ.get("GERMQUANT_GOLDEN_LEGACY")
    if bool(golden) != bool(candidate) or (legacy and not candidate):
        missing = "GERMQUANT_CANDIDATE" if not candidate else "GERMQUANT_GOLDEN"
        print(f"t5_gate: ERROR: {missing} is not set but its partner is. Set GERMQUANT_GOLDEN and "
              f"GERMQUANT_CANDIDATE together (GERMQUANT_GOLDEN_LEGACY only with them), or none of them.")
        return 2
    no_real = None if (golden and candidate) else "GERMQUANT_GOLDEN and GERMQUANT_CANDIDATE not set"
    plan = {
        "ruff": (step_ruff, None),
        "pytest": (step_pytest, None),
        "scope": (lambda: step_scope(args.base), None),
        "golden": (lambda: step_regression_diff("golden", golden, candidate, strict=True), no_real),
        "filelist": (lambda: step_filelist(golden, candidate), no_real),
        "workstation": (step_workstation, "skipped with --no-workstation" if args.no_workstation else None),
        "legacy": (lambda: step_regression_diff("legacy", legacy, candidate, strict=False),
                   None if (legacy and candidate) else "GERMQUANT_GOLDEN_LEGACY and GERMQUANT_CANDIDATE not set"),
    }
    results: list[Result] = []
    for step in STEPS:
        fn, skip_reason = plan[step]
        if step in args.skip:
            skip_reason = "skipped with --skip"
        if skip_reason:
            status = "FAIL" if (args.require_golden and step in GOLDEN_STEPS) else "SKIP"
            results.append(Result(step, status, skip_reason))
            continue
        t0 = time.perf_counter()
        try:
            res = fn()
        except Exception as e:  # noqa: BLE001 - one broken step must not hide the others
            res = Result(step, "FAIL", f"{type(e).__name__}: {e}")
        res.seconds = round(time.perf_counter() - t0, 1)
        results.append(res)

    line, rc = verdict(results)
    text = table(results, _header()) + f"\n\nt5_gate: {line}"
    print("\n" + text)
    if args.markdown:
        Path(args.markdown).write_text(text + "\n", encoding="utf-8")
    if args.json:
        Path(args.json).write_text(json.dumps({"verdict": line, "exit_code": rc,
                                               "steps": [{**asdict(r), "check": r.check} for r in results]},
                                              indent=1), encoding="utf-8")
    return rc


if __name__ == "__main__":
    sys.exit(main())
