"""scripts/t5_scope_check.py and scripts/t5_gate.py (docs/DESIGN_system_packs.md 11.2, gates G3 and G5).

The scope check runs on a throwaway git repository with the frozen layout: an edit, a deletion or a rename
of a frozen file fails, a new file in a frozen folder passes, any change in the sealed folder fails,
uncommitted and untracked changes count, the automatic base prefers the recorded T5 base while main has
not reached it, and a recorded base that has diverged from main charges neither main's own changes nor
the side commits to the branch. The fence files may not be deleted, and the pins file is append-only,
both while the base predates it and once main holds it. The gate's verdicts run with stub steps, and its
file-list check runs on small output trees.
"""
from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod            # dataclasses look their module up while the class is built
    spec.loader.exec_module(mod)
    return mod


scope = _load("t5_scope_check")
gate = _load("t5_gate")

LAYOUT = {
    "config/config.yaml": "a: 1\n",
    "config/profiles/generic.yaml": "extends: config.yaml\n",
    "config/profiles/slurm/config.yaml": "jobs: 1\n",
    "config/channel_maps/map.yaml": "roles: {}\n",
    "tests/golden/synthetic/nuclei.csv": "x\n1\n",
    "tests/golden/synthetic_4ch/nuclei.csv": "x\n2\n",
    "tests/golden/real_manifest.json": "[]\n",
    "src/germquant/benchmark/detector.py": "X = 1\n",
    "src/germquant/pipeline.py": "Y = 1\n",
    "docs/README.md": "text\n",
}
PINS_PATH = "tests/golden/default_invariants.json"
FENCE_PATH = "tests/test_default_invariants.py"
PINS = {
    "format": 1,
    "about": "the fence",
    "base_commit": "",
    "config_hashes": {"config/config.yaml": "84bdda5a5230"},
    "runs": {"synthetic": {"manifest_keys": ["config", "config_hash"],
                           "tif_arrays": {"nuclei_labels.tif": {"sha256": "aa", "scope": "all"}},
                           "montage_pixels": {"sha256": "bb", "scope": "workstation"}}},
}


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@example.com", "-c", "core.autocrlf=false",
                           *args], cwd=repo, capture_output=True, text=True, check=True).stdout.strip()


def _write(repo: Path, rel: str, text: str) -> None:
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8", newline="\n")


def _commit_all(repo: Path, message: str) -> str:
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", message)
    return _git(repo, "rev-parse", "HEAD")


def _make_repo(r: Path) -> Path:
    r.mkdir(parents=True)
    _git(r, "init", "-q")
    _git(r, "checkout", "-q", "-b", "main")
    for rel, text in LAYOUT.items():
        _write(r, rel, text)
    _commit_all(r, "base")
    return r


@pytest.fixture
def repo(tmp_path):
    return _make_repo(tmp_path / "repo")


def _paths(rep: dict, key: str) -> list[str]:
    return sorted(v["path"] for v in rep[key])


def test_clean_and_allowed_changes_pass(repo):
    base = _git(repo, "rev-parse", "HEAD")
    _write(repo, "src/germquant/pipeline.py", "Y = 2\n")                 # not frozen
    _write(repo, "config/profiles/generic_2d.yaml", "extends: generic.yaml\n")   # a new recipe file
    _write(repo, "tests/golden/synthetic_2d/nuclei.csv", "x\n3\n")       # a new golden folder
    rep = scope.check(repo, base)
    assert rep["status"] == "PASS" and rep["violations"] == []
    assert _paths(rep, "allowed") == ["config/profiles/generic_2d.yaml"]


@pytest.mark.parametrize("change, path, rule", [
    ("edit", "config/config.yaml", "R1"),
    ("edit", "config/profiles/generic.yaml", "R1"),
    ("edit", "config/profiles/slurm/config.yaml", "R1"),
    ("edit", "config/channel_maps/map.yaml", "R1"),
    ("delete", "tests/golden/synthetic/nuclei.csv", "11.2"),
    ("edit", "tests/golden/synthetic_4ch/nuclei.csv", "11.2"),
    ("edit", "tests/golden/real_manifest.json", "11.2"),
    ("edit", "src/germquant/benchmark/detector.py", "R8"),
    ("add", "src/germquant/benchmark/new_module.py", "R8"),       # the sealed folder takes no new file
])
def test_frozen_changes_fail_in_the_working_tree(repo, change, path, rule):
    base = _git(repo, "rev-parse", "HEAD")
    if change == "delete":
        (repo / path).unlink()
    else:
        _write(repo, path, "changed\n")
    rep = scope.check(repo, base)                                     # uncommitted (and untracked) changes count
    assert rep["status"] == "FAIL", rep
    assert [(v["path"], v["rule"]) for v in rep["violations"]] == [(path, rule)]


def test_committed_rename_counts_as_a_deletion(repo):
    base = _git(repo, "rev-parse", "HEAD")
    _git(repo, "mv", "config/channel_maps/map.yaml", "config/channel_maps/map_v2.yaml")
    _git(repo, "commit", "-q", "-m", "rename")
    rep = scope.check(repo, base, head="HEAD")
    assert rep["status"] == "FAIL"
    assert [(v["path"], v["change"]) for v in rep["violations"]] == [("config/channel_maps/map.yaml", "deleted")]
    assert _paths(rep, "allowed") == ["config/channel_maps/map_v2.yaml"]
    assert scope.main(["--repo", str(repo), "--base", base, "--head", "HEAD"]) == 1


def test_auto_base_prefers_the_recorded_t5_base_until_main_reaches_it(repo):
    _write(repo, "docs/BENCHMARK_PROTOCOL.md", "protocol\n")
    main_tip = _commit_all(repo, "base with protocol")
    _git(repo, "checkout", "-q", "-b", "feat/tranche5")
    _write(repo, "docs/BENCHMARK_PROTOCOL.md", "protocol, entry 7\n")    # legitimate RAD-51 work, not T5
    t5_base = _commit_all(repo, "protocol entry")
    _write(repo, PINS_PATH, json.dumps({"base_commit": t5_base}))
    _write(repo, "tests/test_fence.py", "X = 1\n")
    _commit_all(repo, "P0")
    rep = scope.check(repo, "auto", head="HEAD")
    assert rep["base"] == t5_base and "recorded" in rep["base_reason"] and rep["second_base"] is None
    assert rep["status"] == "PASS"
    # an explicit base that includes the protocol commit sees it
    assert scope.check(repo, main_tip, head="HEAD")["status"] == "FAIL"
    # once main holds the T5 base and moves on, the merge base with main is the later one
    _git(repo, "checkout", "-q", "main")
    _git(repo, "merge", "-q", "--ff-only", "feat/tranche5")
    _write(repo, "docs/README.md", "text 2\n")
    _commit_all(repo, "main moves on")
    _git(repo, "checkout", "-q", "-b", "t5/p1")
    _write(repo, "src/germquant/packs.py", "P = 1\n")
    _commit_all(repo, "P1")
    rep = scope.check(repo, "auto", head="HEAD")
    assert rep["base"] == _git(repo, "rev-parse", "main") and "merge base" in rep["base_reason"]
    assert rep["status"] == "PASS" and rep["n_changed"] == 1


def test_a_diverged_recorded_base_charges_only_the_branch(repo):
    """The recorded T5 base is a side commit that main never took. Main then updates an R8 file on its own
    and the T5 branch merges main: that update is not the branch's change. A frozen file the branch edits
    itself still fails, even the same R8 file."""
    _write(repo, "docs/BENCHMARK_PROTOCOL.md", "protocol\n")
    _commit_all(repo, "protocol")
    _git(repo, "checkout", "-q", "-b", "feat/tranche5")
    _write(repo, "docs/DESIGN.md", "design\n")
    side = _commit_all(repo, "design, docs only (never merged into main)")
    _write(repo, PINS_PATH, json.dumps({"base_commit": side}))
    _commit_all(repo, "P0")
    _git(repo, "checkout", "-q", "main")
    _write(repo, "docs/BENCHMARK_PROTOCOL.md", "protocol, entry 8\n")    # legitimate RAD-51 work on main
    _commit_all(repo, "protocol entry 8")
    _git(repo, "checkout", "-q", "feat/tranche5")
    _git(repo, "merge", "-q", "--no-ff", "--no-edit", "main")
    rep = scope.check(repo, "auto", head="HEAD")
    assert rep["base"] == side and rep["second_base"] == _git(repo, "rev-parse", "main")
    assert "diverged" in rep["base_reason"]
    assert rep["status"] == "PASS", rep["violations"]
    _write(repo, "config/config.yaml", "a: 2\n")
    _write(repo, "docs/BENCHMARK_PROTOCOL.md", "protocol, entry 8, then a T5 edit\n")
    _commit_all(repo, "T5 edits two frozen files")
    rep = scope.check(repo, "auto", head="HEAD")
    assert rep["status"] == "FAIL"
    assert _paths(rep, "violations") == ["config/config.yaml", "docs/BENCHMARK_PROTOCOL.md"]


def test_unresolvable_base_exits_2(repo):
    assert scope.main(["--repo", str(repo), "--base", "no-such-ref"]) == 2


# ---------------------------------------------------------------------------------------------------
# The fence itself (gate G3): the fence files stay, and existing pins never change
# ---------------------------------------------------------------------------------------------------

SETUPS = ("base predates the pins", "base holds the pins")


@pytest.fixture(scope="module")
def t5_templates(tmp_path_factory):
    """P0 (the pins file and the fence module) on feat/tranche5, then a package branch feat/t5-p1, built
    once per setup. With "base predates the pins" main has not taken P0, so the pins are compared with the
    commit that added them; with "base holds the pins" main took P0 first, so they are compared with the
    base."""
    out = {}
    for setup in SETUPS:
        r = _make_repo(tmp_path_factory.mktemp("t5_template") / "repo")
        base = _git(r, "rev-parse", "HEAD")
        _git(r, "checkout", "-q", "-b", "feat/tranche5")
        _write(r, PINS_PATH, json.dumps({**PINS, "base_commit": base}, indent=1) + "\n")
        _write(r, FENCE_PATH, "X = 1\n")
        p0 = _commit_all(r, "P0")
        if setup == "base holds the pins":
            _git(r, "checkout", "-q", "main")
            _git(r, "merge", "-q", "--ff-only", "feat/tranche5")
        _git(r, "checkout", "-q", "-b", "feat/t5-p1")
        out[setup] = (r, p0)
    return out


@pytest.fixture(params=SETUPS)
def t5_branch(request, t5_templates, tmp_path):
    """A private copy of one template repository."""
    template, p0 = t5_templates[request.param]
    repo = tmp_path / "repo"
    shutil.copytree(template, repo)
    return SimpleNamespace(repo=repo, p0=p0, setup=request.param)


def _edit_pins(repo: Path, edit) -> None:
    pins = json.loads((repo / PINS_PATH).read_text(encoding="utf-8"))
    edit(pins)
    _write(repo, PINS_PATH, json.dumps(pins, indent=1) + "\n")


# name -> (edit of the pins, given the P0 commit; expected status; the pin that fails, or the allowed change)
PIN_EDITS = {
    "a pin value changed": (lambda p, p0: p["config_hashes"].update({"config/config.yaml": "000000000000"}),
                            "FAIL", "config_hashes/config/config.yaml"),
    "a pin removed": (lambda p, p0: p["runs"]["synthetic"].pop("manifest_keys"),
                      "FAIL", "runs/synthetic/manifest_keys"),
    "a key added to a pinned list": (lambda p, p0: p["runs"]["synthetic"]["manifest_keys"].append("runtime"),
                                     "FAIL", "runs/synthetic/manifest_keys"),
    "base_commit moved": (lambda p, p0: p.update(base_commit=p0), "FAIL", "base_commit"),
    "a pin added": (lambda p, p0: p["config_hashes"].update({"config/profiles/generic_2d.yaml": "123456789abc"}),
                    "PASS", None),
    "the free text edited": (lambda p, p0: p.update(about="the fence, reworded"), "PASS", None),
    "a scope relaxed": (lambda p, p0: p["runs"]["synthetic"]["montage_pixels"].update(scope="all"),
                        "PASS", "runs/synthetic/montage_pixels/scope"),
    "a scope made workstation-only": (
        lambda p, p0: p["runs"]["synthetic"]["tif_arrays"]["nuclei_labels.tif"].update(scope="workstation"),
        "PASS", "runs/synthetic/tif_arrays/nuclei_labels.tif/scope"),
}


@pytest.mark.parametrize("name", list(PIN_EDITS))
def test_existing_pins_are_append_only(t5_branch, name):
    edit, status, where = PIN_EDITS[name]
    repo = t5_branch.repo
    _edit_pins(repo, lambda pins: edit(pins, t5_branch.p0))
    for head in (None, "HEAD"):                                       # uncommitted, then committed
        if head:
            _commit_all(repo, name)
        rep = scope.check(repo, "auto", head=head)
        assert rep["status"] == status, rep["violations"]
        failing = [f"{PINS_PATH}:{where}"] if status == "FAIL" else []
        assert _paths(rep, "violations") == failing
        scoped = [f"{PINS_PATH}:{where}"] if "scope" in name else []
        assert _paths(rep, "allowed") == scoped
        assert rep["fence"]["pins_reference"] is not None and rep["fence"]["pins_checked"] == 8   # "about" is free text
        if name != "base_commit moved":
            reason = "the base" if t5_branch.setup == "base holds the pins" else "first added"
            assert reason in rep["fence"]["pins_reference_reason"]


@pytest.mark.parametrize("path", [FENCE_PATH, PINS_PATH])
def test_the_fence_files_may_not_be_deleted_or_renamed(t5_branch, path):
    repo = t5_branch.repo
    (repo / path).unlink()
    rep = scope.check(repo, "auto")                                   # an uncommitted deletion
    assert rep["status"] == "FAIL"
    assert [(v["path"], v["change"]) for v in rep["violations"]] == [(path, "deleted")]
    _git(repo, "checkout", "--", path)
    moved = path.replace("tests/", "tests/old_")
    (repo / moved).parent.mkdir(parents=True, exist_ok=True)
    _git(repo, "mv", path, moved)                                     # a committed rename
    _commit_all(repo, "move the fence")
    rep = scope.check(repo, "auto", head="HEAD")
    assert rep["status"] == "FAIL"
    assert [(v["path"], v["change"]) for v in rep["violations"]] == [(path, "deleted")]


# ---------------------------------------------------------------------------------------------------
# scripts/t5_gate.py
# ---------------------------------------------------------------------------------------------------

@pytest.fixture
def stub_gate(monkeypatch):
    """t5_gate with stub steps (no ruff, pytest or git runs) and no golden variables; returns the list of
    steps that ran."""
    ran: list[str] = []

    def stub(step):
        def run(*args):
            ran.append(step)
            return gate.Result(step, "PASS", "stub")
        return run

    for step in ("ruff", "pytest", "scope", "workstation"):
        monkeypatch.setattr(gate, f"step_{step}", stub(step))
    monkeypatch.setattr(gate, "_header", lambda: "t5_gate: test")
    for var in ("GERMQUANT_GOLDEN", "GERMQUANT_CANDIDATE", "GERMQUANT_GOLDEN_LEGACY"):
        monkeypatch.delenv(var, raising=False)
    return ran


def test_gate_runs_the_workstation_pins_by_default(stub_gate):
    assert gate.main([]) == 0
    assert stub_gate == ["ruff", "pytest", "scope", "workstation"]
    stub_gate.clear()
    assert gate.main(["--no-workstation"]) == 0
    assert stub_gate == ["ruff", "pytest", "scope"]


def test_gate_verdicts(stub_gate, tmp_path):
    md, js = tmp_path / "gate.md", tmp_path / "gate.json"
    # skipping ruff, pytest or scope never ends in PASS
    assert gate.main(["--skip", "pytest", "--markdown", str(md)]) == 3
    assert md.read_text(encoding="utf-8").rstrip().endswith("t5_gate: INCOMPLETE (skipped: pytest)")
    # a G4 run: a skipped golden step is a FAIL, the comparison with the older golden tree (legacy) included
    assert gate.main(["--require-golden", "--json", str(js)]) == 1
    rep = json.loads(js.read_text(encoding="utf-8"))
    assert {s["step"]: s["status"] for s in rep["steps"]} == {
        "ruff": "PASS", "pytest": "PASS", "scope": "PASS", "golden": "FAIL", "filelist": "FAIL",
        "workstation": "PASS", "legacy": "FAIL"}
    assert rep["verdict"] == "FAIL (golden, filelist, legacy)" and rep["exit_code"] == 1


@pytest.mark.parametrize("names", [("GERMQUANT_GOLDEN",), ("GERMQUANT_CANDIDATE",), ("GERMQUANT_GOLDEN_LEGACY",),
                                   ("GERMQUANT_GOLDEN", "GERMQUANT_GOLDEN_LEGACY")])
def test_gate_refuses_half_set_golden_variables(stub_gate, monkeypatch, tmp_path, names):
    for name in names:
        monkeypatch.setenv(name, str(tmp_path))
    assert gate.main([]) == 2
    assert stub_gate == []                                            # nothing ran


@pytest.mark.parametrize("path_version, expected", [(None, "uvx"), ("0.0.1", "uvx"), ("pin", "path")])
def test_gate_lints_with_the_ci_ruff_version(monkeypatch, path_version, expected):
    pin = gate.ruff_pin()
    assert pin and pin[0].isdigit(), "constraints/ci-cpu-2026-10.txt has no ruff== line"
    version = pin if path_version == "pin" else path_version
    monkeypatch.setattr(gate.shutil, "which", lambda name: {"uvx": "uvx", "ruff": "/fake/ruff" if version else None}
                        .get(name))
    monkeypatch.setattr(gate.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=0, stdout=f"ruff {version}\n"))
    cmd, _ = gate.ruff_command()
    assert cmd == (["/fake/ruff"] if expected == "path" else ["uvx", f"ruff@{pin}"])


def _tree(root: Path, files: list[str]) -> Path:
    for rel in files:
        _write(root, rel, "x\n")
    return root


def test_file_list_check(tmp_path):
    files = ["tag/img__nuclei.csv", "tag/img__nuclei_labels.tif", "tag/run_manifest.json",
             "tag2/other__spots.csv", "golden_manifest.json"]
    g = _tree(tmp_path / "g", files)
    assert gate.compare_file_lists(g, _tree(tmp_path / "same", files)) == []
    extra = _tree(tmp_path / "extra", [*files, "tag/img__new_table.csv"])
    assert gate.compare_file_lists(g, extra) == ["tag/img: missing [], added ['img__new_table.csv']"]
    fewer = _tree(tmp_path / "fewer", [f for f in files if not f.startswith("tag2/")])
    assert gate.compare_file_lists(g, fewer) == ["missing from the candidate: tag2/other"]
    no_manifest = _tree(tmp_path / "no_manifest", [f for f in files if f != "tag/run_manifest.json"])
    assert gate.compare_file_lists(g, no_manifest) == ["missing from the candidate: tag/(other files)"]
    (tmp_path / "empty").mkdir()
    assert gate.compare_file_lists(tmp_path / "empty", g) == [f"no images (*__*.csv) under {tmp_path / 'empty'}"]
