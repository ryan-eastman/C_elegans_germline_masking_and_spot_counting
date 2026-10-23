"""scripts/t5_style_check.py (gate G7; the Train B plan, decision D3 and task 10.1).

On throwaway git repositories: an en or em dash in a line the diff adds fails, in code, docs and workflows,
while a dash that was there at the base passes; a non-ASCII character fails in a text file added under
config/, packs/, tests/fixtures/ or tests/golden/, in any YAML file, and in the lines added to such a file,
and passes elsewhere; binary files are skipped (git's detection, a NUL byte, the suffix, a .gitattributes
rule); untracked files count without --head and ignored ones never; the base is the scope check's, so
main's own lines are never charged to the branch; and the frozen R8 paths are never opened. The last test
runs the check on this checkout when its full history and origin/main are present, so t5_gate's pytest step
enforces G7 on the workstation; CI's shallow checkout skips it, and the T5 step of the workflow runs the
script on the pushed commit instead.

Every non-ASCII character below is built with chr(), so this file itself stays ASCII.
"""
from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")

EN, EM, MU, E_ACUTE = chr(0x2013), chr(0x2014), chr(0x00B5), chr(0x00E9)
MICRO = "non-ASCII U+00B5 MICRO SIGN"
ACUTE = "non-ASCII U+00E9 LATIN SMALL LETTER E WITH ACUTE"


def _load(name: str):
    path = ROOT / "scripts" / f"{name}.py"
    mod = sys.modules.get(name)
    if mod is not None and Path(getattr(mod, "__file__", "")).resolve() == path.resolve():
        return mod
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


scope = _load("t5_scope_check")
style = _load("t5_style_check")

BASE = {
    "src/germquant/pipeline.py": f"X = 1  # an old {EN} dash, there at the base\n",
    "docs/README.md": "text\n",
    "tests/golden/pins.json": '{"a": "' + MU + 'm"}\n',          # an old non-ASCII line, there at the base
    ".gitignore": "*.log\n",
    ".gitattributes": "*.dat binary\n",
    "docs/BENCHMARK_PROTOCOL.md": "protocol\n",
    "docs/FINAL_SCORING.md": "scoring\n",
    "tests/golden/rad51_split.json": "{}\n",
    "src/germquant/validate/benchmark_split.py": "S = 1\n",
    "src/germquant/benchmark/detector.py": "D = 1\n",
}


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@example.com", "-c", "core.autocrlf=false",
                           *args], cwd=repo, capture_output=True, text=True, check=True).stdout.strip()


def _write(repo: Path, rel: str, data) -> None:
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(data.encode("utf-8") if isinstance(data, str) else data)


def _commit(repo: Path, message: str) -> str:
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", message)
    return _git(repo, "rev-parse", "HEAD")


@pytest.fixture(scope="module")
def template(tmp_path_factory):
    """The base layout, committed once per module on main."""
    r = tmp_path_factory.mktemp("style_template") / "repo"
    r.mkdir()
    _git(r, "init", "-q")
    _git(r, "checkout", "-q", "-b", "main")
    for rel, text in BASE.items():
        _write(r, rel, text)
    return r, _commit(r, "base")


@pytest.fixture
def repo(template, tmp_path):
    """(a private copy of the template repository, its base commit)."""
    src, base = template
    shutil.copytree(src, tmp_path / "repo")
    return tmp_path / "repo", base


def _found(rep: dict) -> list[tuple]:
    return sorted((f["path"], f["line"], f["rule"], f["what"]) for f in rep["findings"])


def test_both_rules_on_the_working_tree(repo):
    r, base = repo
    # staged changes reach the check through git's patch
    _write(r, "src/germquant/pipeline.py", f"X = 1  # an old {EN} dash, there at the base\nY = 2  # a new {EN} one\n")
    _write(r, "docs/new.md", f"title\nan {EM} em dash\n")
    _write(r, ".github/workflows/ci.yml", f"name: CI\n# 5 {MU}s\n")                       # (b): YAML anywhere
    _write(r, "config/profiles/new.yaml", f"io:\n  unit: {MU}m\n")                          # (b)
    _write(r, "tests/golden/pins.json", '{"a": "' + MU + 'm"}\n{"b": 1}\n{"c": "caf' + E_ACUTE + '"}\n')
    _write(r, "docs/with space.md", f"++ b/not_a_header\n+++ b/nor_this\nx {EN}\n")         # header look-alikes
    _write(r, "docs/crlf.md", f"a\r\nb {EM}\r\n")
    _write(r, "docs/no_eol.md", f"x {EN}")                                                   # no final newline
    _write(r, "tests/golden/synthetic_2d/labels.tif", b"II*\x00" + EN.encode() + b"\x00" * 16)   # binary (NUL)
    _write(r, "tests/fixtures/x/flat.png", EN.encode() + b" no NUL byte, but a .png\n")     # binary (suffix)
    _git(r, "add", "-A")
    # untracked files are read from the working tree
    _write(r, "packs/p/notes.txt", f"{MU}m\n")                                               # (b)
    _write(r, "tests/fixtures/x/data.csv", f"a\n{MU}\n")                                     # (b)
    _write(r, "docs/any.yml", f"k: {E_ACUTE}\n")                                             # (b): YAML anywhere
    _write(r, "docs/notes.md", f"{MU}m is fine in docs\n")                                   # passes
    _write(r, "src/germquant/io/new.py", f"UNIT = '{MU}m'\n")                                # passes
    _write(r, "tests/fixtures/x/helper.py", f"UNIT = '{MU}m'\n")                             # passes: not a text suffix
    _write(r, "docs/scratch.md", f"draft {EM}\n")                                            # (a), untracked
    _write(r, "run.log", f"{EN} ignored\n")                                                  # .gitignore
    _write(r, "tests/golden/x/raw.dat", f"{EN}\n")                                           # .gitattributes binary
    _write(r, "tests/golden/x/blob.bin", b"\x00" + EN.encode())                              # binary (NUL)
    rep = style.check(r, base)
    assert rep["status"] == "FAIL"
    assert _found(rep) == sorted([
        ("src/germquant/pipeline.py", 2, "a", "U+2013 EN DASH"),          # the old dash on line 1 passes
        ("docs/new.md", 2, "a", "U+2014 EM DASH"),
        (".github/workflows/ci.yml", 2, "b", MICRO),
        ("config/profiles/new.yaml", 2, "b", MICRO),
        ("tests/golden/pins.json", 3, "b", ACUTE),                        # the old line 1 passes
        ("docs/with space.md", 3, "a", "U+2013 EN DASH"),
        ("docs/crlf.md", 2, "a", "U+2014 EM DASH"),
        ("docs/no_eol.md", 1, "a", "U+2013 EN DASH"),
        ("packs/p/notes.txt", 1, "b", MICRO),
        ("tests/fixtures/x/data.csv", 2, "b", MICRO),
        ("docs/any.yml", 1, "b", ACUTE),
        ("docs/scratch.md", 1, "a", "U+2014 EM DASH"),
    ])
    assert [s["path"] for s in rep["skipped"]] == sorted([
        "tests/golden/synthetic_2d/labels.tif", "tests/fixtures/x/flat.png", "tests/golden/x/raw.dat",
        "tests/golden/x/blob.bin"])
    assert "run.log" not in rep["files"] and "docs/notes.md" in rep["files"]
    # with --head the commit is checked, and nothing was committed after the base
    rep = style.check(r, base, head="HEAD")
    assert rep["status"] == "PASS" and rep["n_files"] == 0


def test_a_commit_the_report_and_the_exit_codes(repo, tmp_path):
    r, base = repo
    assert style.main(["--repo", str(r), "--base", base, "--head", "HEAD"]) == 0
    _write(r, "src/germquant/new_module.py", f'"""A docstring with an {EM} em dash."""\n')
    head = _commit(r, "a dash in code")
    out = tmp_path / "style.json"
    assert style.main(["--repo", str(r), "--base", base, "--head", "HEAD", "--json", str(out)]) == 1
    rep = json.loads(out.read_text(encoding="utf-8"))
    assert rep["status"] == "FAIL" and rep["head"] == head and rep["base"] == base
    assert _found(rep) == [("src/germquant/new_module.py", 1, "a", "U+2014 EM DASH")]
    assert style.main(["--repo", str(r), "--base", "no-such-ref"]) == 2


def test_the_auto_base_charges_only_the_branch(repo):
    """The scope check's base: main's own lines are not the branch's, also when the recorded T5 base has
    diverged from main (a line then counts only when it is added against both bases)."""
    r, _ = repo
    _git(r, "checkout", "-q", "-b", "feat/tranche5")
    _write(r, "docs/DESIGN.md", "design\n")
    side = _commit(r, "a side commit that main never takes")
    _write(r, "tests/golden/default_invariants.json", json.dumps({"base_commit": side}) + "\n")
    _commit(r, "P0")
    _git(r, "checkout", "-q", "main")
    _write(r, "docs/README.md", f"text\nmain's own {EN} line\n")
    _commit(r, "main adds a dash of its own")
    _git(r, "checkout", "-q", "feat/tranche5")
    _git(r, "merge", "-q", "--no-ff", "--no-edit", "main")
    rep = style.check(r, "auto", head="HEAD")
    assert rep["base"] == side and rep["second_base"] == _git(r, "rev-parse", "main")
    assert rep["status"] == "PASS", rep["findings"]
    _write(r, "docs/README.md", f"text\nmain's own {EN} line\nthe branch's own {EM} line\n")
    assert _found(style.check(r, "auto")) == [("docs/README.md", 3, "a", "U+2014 EM DASH")]


def test_the_frozen_r8_paths_are_never_opened(repo, monkeypatch):
    r, base = repo
    for p in style.R8_FILES:
        _write(r, p, f"edited {EM}\n")              # the scope check fails these edits; this check never reads them
    _write(r, "src/germquant/benchmark/detector.py", f"D = 2  # {EN}\n")
    _write(r, "src/germquant/benchmark/new_module.py", f"N = 1  # {EN}\n")      # untracked, in the sealed folder
    _git(r, "add", "docs/BENCHMARK_PROTOCOL.md")
    calls, reads = [], []
    real_git, real_read = style._git_bytes, style.read_untracked
    monkeypatch.setattr(style, "_git_bytes", lambda repo, *a: calls.append(a) or real_git(repo, *a))
    monkeypatch.setattr(style, "read_untracked", lambda repo, p: reads.append(p) or real_read(repo, p))
    rep = style.check(r, base)
    assert rep["status"] == "PASS" and rep["files"] == [] and rep["skipped"] == [] and reads == []
    on_paths = [a for a in calls if a[0] in ("diff", "ls-files")]
    assert len(on_paths) == 3 and all(set(style.R8_EXCLUDES) <= set(a) for a in on_paths)
    assert set(style.R8_FILES) == {p for p, rule in scope.FROZEN_FILES.items() if rule == "R8"}


def test_text_suffixes_are_those_of_the_packs_format():
    from germquant.packs.format import TEXT_SUFFIXES

    assert style.TEXT_SUFFIXES == TEXT_SUFFIXES


def _checkout_problem() -> str | None:
    """None when this checkout can run the check on itself, otherwise why it cannot."""
    def git(*args: str) -> subprocess.CompletedProcess:
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True)

    if git("rev-parse", "--is-inside-work-tree").stdout.strip() != "true":
        return "not a git checkout"
    if git("rev-parse", "--is-shallow-repository").stdout.strip() != "false":
        return "a shallow clone, as in CI: the history the base needs is missing (CI's T5 step runs the check)"
    if git("rev-parse", "--verify", "--quiet", "origin/main^{commit}").returncode != 0:
        return "no origin/main in this checkout"
    return None


def test_this_checkout_passes_the_style_check():
    """G7 on the workstation: the working tree of this checkout against the scope check's auto base."""
    problem = _checkout_problem()
    if problem:
        pytest.skip(problem)
    rep = style.check(ROOT)
    assert rep["status"] == "PASS", "\n".join(f"{f['path']}:{f['line']}: ({f['rule']}) {f['what']}: {f['text']}"
                                              for f in rep["findings"])
