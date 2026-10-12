"""Tranche 5 scope check: fail a diff that modifies or deletes a frozen path or changes an existing fence pin.

    python scripts/t5_scope_check.py [--base auto|REF] [--head REF] [--json OUT] [--list]

docs/DESIGN_system_packs.md, section 2 (rules R1 and R8), section 11.2 (the fence) and gate G3. Frozen paths:

  * R1, frozen texts: config/config.yaml, config/config_ccw77.yaml, config/config_n2dryice.yaml, and every
    file that exists at the base under config/profiles/ (generic.yaml, the bench profiles and the Snakemake
    profiles included) and under config/channel_maps/. New recipes are new files, so a file ADDED to those
    folders passes.
  * R8, the RAD-51 seal: tests/golden/rad51_split.json, src/germquant/validate/benchmark_split.py,
    docs/BENCHMARK_PROTOCOL.md, docs/FINAL_SCORING.md and the folder src/germquant/benchmark/, which is
    sealed: nothing in it may be modified, deleted or added.
  * 11.2: tests/golden/synthetic/ and tests/golden/synthetic_4ch/, named one by one (a glob such as
    synthetic*/ would also freeze the golden folders that later packages add; a new file inside these two
    passes), and tests/golden/real_manifest.json.

A rename counts as a deletion. For these paths the check reads only file NAMES from git and never opens a
frozen file.

The fence itself (gate G3: "the P0 pins plus the pins added for new files"):

  * tests/golden/default_invariants.json and tests/test_default_invariants.py may not be deleted or renamed;
  * the pins file is append-only. Every pin (every leaf of the JSON; a list is one leaf) of the old copy
    must be in the new copy with the same value, and new pins pass. Two exceptions: the free text under
    "about", and a "scope" leaf, which says where a pin is enforced, not what it pins: it may switch
    between "all" and "workstation" and is then listed as allowed, so the PR shows it (design 11.2).
    base_commit is a pin too, so a diff cannot move its own base;
  * the old copy is the pins file at the base. While the base predates the file (tranche 5 work before
    the train that holds P0 has reached main), it is the file as written by the commit that first added
    it to the head's history. A pin that a later package adds is protected once it is in the base.

What is compared with the base:
  * without --head: the working tree (committed, staged, unstaged and untracked files), so an uncommitted
    edit fails too;
  * with --head REF: the commit REF (CI passes the pushed commit or the pull request head).

The base: --base auto (the default) takes the later of two commits, the T5 base recorded in
tests/golden/default_invariants.json (when it is an ancestor of the head) and the merge base of the head
with origin/main (or main). The first keeps the diff to tranche 5 while the base has not reached main; the
second keeps changes that main made on its own out of a branch's diff (the RAD-51 work may legitimately
update the split record on main). If the two have diverged (neither is an ancestor of the other), a path
counts only when it differs from both, so neither main's own changes nor the commits beside the recorded
base are charged to the branch. --base REF names the base explicitly.

T5 branches (where CI runs this check, .github/workflows/ci.yml): a branch, or a pull request head or base,
whose name contains "tranche5", contains a path segment that starts with "t5" (feat/t5-p1,
feat/t5/train-a) or starts with "t5". Every T5 package and integration branch must be named that way.

Exit codes: 0 nothing frozen touched, 1 a frozen path modified, deleted or added to a sealed folder, a
fence file deleted or an existing pin changed or removed, 2 the base, the head or the old copy of the pins
file could not be read.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PINS_FILE = "tests/golden/default_invariants.json"
FENCE_MODULE = "tests/test_default_invariants.py"
MAIN_REFS = ("origin/main", "main")

# path -> the rule that freezes it
FROZEN_FILES = {
    "config/config.yaml": "R1",
    "config/config_ccw77.yaml": "R1",
    "config/config_n2dryice.yaml": "R1",
    "tests/golden/rad51_split.json": "R8",
    "src/germquant/validate/benchmark_split.py": "R8",
    "docs/BENCHMARK_PROTOCOL.md": "R8",
    "docs/FINAL_SCORING.md": "R8",
    "tests/golden/real_manifest.json": "11.2",
}
# every file in these folders at the base is frozen; a new file may be added
FROZEN_FOLDERS = {
    "config/profiles/": "R1",
    "config/channel_maps/": "R1",
    "tests/golden/synthetic/": "11.2",
    "tests/golden/synthetic_4ch/": "11.2",
}
# nothing in these folders may be modified, deleted or added
SEALED_FOLDERS = {
    "src/germquant/benchmark/": "R8",
}
# the fence (gate G3): never deleted or renamed; the pins file is also append-only
FENCE_FILES = (PINS_FILE, FENCE_MODULE)
PIN_SCOPES = ("all", "workstation")     # the values a "scope" leaf may switch between
FREE_TEXT = ("about",)                  # top-level keys of the pins file that are prose, not pins

_KIND = {"M": "modified", "D": "deleted", "A": "added", "T": "type changed", "U": "unmerged",
         "C": "copied", "R": "renamed"}


class GitError(RuntimeError):
    pass


def _git(repo: Path, *args: str, check: bool = True) -> subprocess.CompletedProcess:
    r = subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, encoding="utf-8",
                       errors="replace")
    if check and r.returncode != 0:
        raise GitError(f"git {' '.join(args)}: {(r.stderr or r.stdout).strip()}")
    return r


def _commit(repo: Path, ref: str) -> str | None:
    r = _git(repo, "rev-parse", "--verify", "--quiet", f"{ref}^{{commit}}", check=False)
    return r.stdout.strip() if r.returncode == 0 and r.stdout.strip() else None


def _is_ancestor(repo: Path, a: str, b: str) -> bool:
    return _git(repo, "merge-base", "--is-ancestor", a, b, check=False).returncode == 0


def _exists(repo: Path, commit: str, path: str) -> bool:
    return _git(repo, "cat-file", "-e", f"{commit}:{path}", check=False).returncode == 0


def _show(repo: Path, commit: str, path: str) -> str | None:
    r = _git(repo, "show", f"{commit}:{path}", check=False)
    return r.stdout if r.returncode == 0 else None


def _head_text(repo: Path, head: str | None, path: str) -> str | None:
    """The file at ``head`` (a commit), or in the working tree when ``head`` is None; None if absent."""
    if head is not None:
        return _show(repo, head, path)
    p = repo / path
    return p.read_text(encoding="utf-8") if p.is_file() else None


def recorded_base(repo: Path, head: str | None) -> str | None:
    """base_commit from the pins file: the committed one at ``head``, or the working-tree copy."""
    text = _head_text(repo, head, PINS_FILE)
    if not text:
        return None
    try:
        return json.loads(text).get("base_commit") or None
    except (ValueError, AttributeError):
        return None


def resolve_base(repo: Path, base: str, head_commit: str, head: str | None) -> tuple[str, str, str | None]:
    """(base commit, why it was chosen, second base). The second base is None except when the recorded T5
    base and the merge base with main have diverged; a path then counts only when it differs from both."""
    if base != "auto":
        sha = _commit(repo, base)
        if sha is None:
            raise GitError(f"cannot resolve the base {base!r} to a commit")
        return sha, f"--base {base}", None
    candidates: list[tuple[str, str]] = []
    pinned = recorded_base(repo, head)
    if pinned:
        sha = _commit(repo, pinned)
        if sha is not None and _is_ancestor(repo, sha, head_commit):
            candidates.append((sha, f"the T5 base recorded in {PINS_FILE}"))
    for ref in MAIN_REFS:
        main = _commit(repo, ref)
        if main is None:
            continue
        mb = _git(repo, "merge-base", head_commit, main, check=False).stdout.strip()
        if mb:
            candidates.append((mb, f"the merge base with {ref}"))
        break
    if not candidates:
        raise GitError("no base found (no recorded T5 base that is an ancestor of the head, no main branch); "
                       "pass --base REF")
    if len(candidates) == 1 or candidates[0][0] == candidates[1][0]:
        return candidates[0][0], candidates[0][1], None
    (p, p_why), (m, m_why) = candidates
    if _is_ancestor(repo, m, p):
        return p, f"{p_why}, later than {m_why}", None
    if _is_ancestor(repo, p, m):
        return m, f"{m_why}, later than the recorded T5 base", None
    return p, (f"{p_why}, which has diverged from {m_why} ({m[:12]}): a path counts only when it differs "
               f"from both"), m


def _split_z(text: str) -> list[str]:
    return [t for t in text.split("\0") if t]


def changed_paths(repo: Path, base: str, head: str | None) -> list[tuple[str, str]]:
    """[(status letter, path)] of every path that differs between the base and the head (a commit), or
    the working tree when ``head`` is None (untracked files count as added)."""
    args = ["diff", "--name-status", "--no-renames", "--ignore-submodules", "-z", base]
    if head is not None:
        args.append(head)
    tokens = _split_z(_git(repo, *args, "--").stdout)
    out = [(tokens[i][0], tokens[i + 1]) for i in range(0, len(tokens) - 1, 2)]
    if head is None:
        untracked = _split_z(_git(repo, "ls-files", "--others", "--exclude-standard", "-z").stdout)
        out += [("A", p) for p in untracked]
    return out


def files_at(repo: Path, commit: str, paths) -> set[str]:
    return set(_split_z(_git(repo, "ls-tree", "-r", "--name-only", "-z", commit, "--", *paths).stdout))


def _folder(path: str, folders: dict[str, str]) -> str | None:
    return next((f for f in folders if path.startswith(f)), None)


def classify(changes: list[tuple[str, str]], base_files: set[str]) -> tuple[list[dict], list[dict]]:
    """(violations, allowed additions). ``base_files`` lists the frozen and sealed folders at the base."""
    bad: list[dict] = []
    allowed: list[dict] = []
    for status, path in changes:
        kind = _KIND.get(status, status)
        if path in FROZEN_FILES:
            bad.append({"path": path, "change": kind, "rule": FROZEN_FILES[path], "why": "frozen file"})
            continue
        sealed = _folder(path, SEALED_FOLDERS)
        if sealed:
            bad.append({"path": path, "change": kind, "rule": SEALED_FOLDERS[sealed],
                        "why": f"sealed folder {sealed}"})
            continue
        folder = _folder(path, FROZEN_FOLDERS)
        if folder:
            if status == "A" and path not in base_files:
                allowed.append({"path": path, "change": kind, "rule": FROZEN_FOLDERS[folder],
                                "why": f"new file in {folder}"})
            else:
                bad.append({"path": path, "change": kind, "rule": FROZEN_FOLDERS[folder],
                            "why": f"file that existed at the base in {folder}"})
    return bad, allowed


def _leaves(obj, prefix: tuple[str, ...] = ()):
    """(key path, value) of every leaf of a JSON value: a list is one leaf, an empty object has none."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from _leaves(v, (*prefix, str(k)))
    else:
        yield prefix, obj


def _short(value) -> str:
    s = json.dumps(value)
    return s if len(s) <= 60 else s[:57] + "..."


def compare_pins(old, new) -> tuple[list[dict], list[dict], int]:
    """(violations, allowed scope switches, number of old pins checked): every leaf of ``old`` must be in
    ``new`` with the same value, except the free text and a "scope" leaf switching within PIN_SCOPES."""
    new_leaves = dict(_leaves(new))
    bad: list[dict] = []
    allowed: list[dict] = []
    n = 0
    for path, value in _leaves(old):
        if path and path[0] in FREE_TEXT:
            continue
        n += 1
        where = f"{PINS_FILE}:{'/'.join(path)}"
        if path not in new_leaves:
            bad.append({"path": where, "change": "removed", "rule": "G3",
                        "why": "existing pin, pins are only ever added"})
            continue
        got = new_leaves[path]
        if got == value:
            continue
        if path[-1:] == ("scope",) and value in PIN_SCOPES and got in PIN_SCOPES:
            allowed.append({"path": where, "change": f"scope {value} -> {got}", "rule": "11.2",
                            "why": "where the pin is enforced, not its value"})
            continue
        bad.append({"path": where, "change": f"changed from {_short(value)} to {_short(got)}", "rule": "G3",
                    "why": "existing pin, pins are only ever added"})
    return bad, allowed, n


def _first_added(repo: Path, head_commit: str, path: str) -> str | None:
    """The first commit in the head's history that added ``path``, or None."""
    r = _git(repo, "log", "--diff-filter=A", "--format=%H", "--reverse", head_commit, "--", path, check=False)
    shas = r.stdout.split()
    return shas[0] if shas else None


def check_fence(repo: Path, bases: list[str], head_commit: str, head: str | None) -> tuple[list[dict], list[dict],
                                                                                          dict]:
    """Gate G3 on the fence files: (violations, allowed scope switches, what the pins were compared with).
    ``bases`` are the base and, for diverged bases, the second base: the first that holds a fence file
    gives its old copy, otherwise the commit that first added it to the head's history does."""
    bad: list[dict] = []
    allowed: list[dict] = []
    info = {"files": list(FENCE_FILES), "pins_reference": None, "pins_reference_reason": None,
            "pins_checked": 0}
    for path in FENCE_FILES:
        ref = next((b for b in bases if _exists(repo, b, path)), None)
        why = "the base"
        if ref is None:
            ref = _first_added(repo, head_commit, path)
            why = "the commit that first added it (the base predates the file)"
        if ref is None:
            continue                       # never part of the head's history: nothing to protect yet
        new_text = _head_text(repo, head, path)
        if new_text is None:
            bad.append({"path": path, "change": "deleted", "rule": "G3",
                        "why": "the fence; it may not be deleted or renamed"})
            continue
        if path != PINS_FILE:
            continue
        try:
            old = json.loads(_show(repo, ref, path) or "")
        except ValueError as e:
            raise GitError(f"cannot read the old copy of {path} at {ref[:12]}: {e}") from e
        try:
            new = json.loads(new_text)
        except ValueError as e:
            bad.append({"path": path, "change": "unreadable", "rule": "G3", "why": f"not valid JSON ({e})"})
            continue
        b, a, n = compare_pins(old, new)
        bad += b
        allowed += a
        info.update(pins_reference=ref, pins_reference_reason=why, pins_checked=n)
    return bad, allowed, info


def check(repo: Path = REPO, base: str = "auto", head: str | None = None) -> dict:
    """Run the check; returns the report (``report["status"]`` is PASS or FAIL). Raises GitError."""
    head_commit = _commit(repo, head or "HEAD")
    if head_commit is None:
        raise GitError(f"cannot resolve the head {head or 'HEAD'!r} to a commit")
    base_sha, why, second = resolve_base(repo, base, head_commit, head)
    folders = [*FROZEN_FOLDERS, *SEALED_FOLDERS]
    base_files = files_at(repo, base_sha, folders)
    missing = sorted(set(FROZEN_FILES) - files_at(repo, base_sha, list(FROZEN_FILES)))
    changes = changed_paths(repo, base_sha, head)
    if second:
        # diverged bases: a path counts only when it differs from both, and a file that existed at
        # either base is an existing file
        also = {p for _, p in changed_paths(repo, second, head)}
        changes = [(s, p) for s, p in changes if p in also]
        base_files |= files_at(repo, second, folders)
    bad, allowed = classify(changes, base_files)
    fence_bad, fence_allowed, fence = check_fence(repo, [base_sha, *([second] if second else [])], head_commit,
                                                  head)
    bad += fence_bad
    allowed += fence_allowed
    return {"status": "FAIL" if bad else "PASS", "base": base_sha, "base_reason": why, "second_base": second,
            "head": head_commit if head else "working tree", "n_changed": len(changes),
            "violations": bad, "allowed": allowed, "base_lacks": missing, "fence": fence,
            "frozen": {"files": FROZEN_FILES, "folders": FROZEN_FOLDERS, "sealed": SEALED_FOLDERS}}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default="auto", help="base commit or ref (default: auto, see above)")
    ap.add_argument("--head", help="compare this commit instead of the working tree")
    ap.add_argument("--repo", default=str(REPO), help=argparse.SUPPRESS)
    ap.add_argument("--json", help="write the report here")
    ap.add_argument("--list", action="store_true", help="also print every frozen path at the base")
    args = ap.parse_args(argv)
    repo = Path(args.repo)
    try:
        rep = check(repo, args.base, args.head)
    except GitError as e:
        print(f"t5_scope_check: ERROR: {e}")
        return 2
    pins_bad = [v for v in rep["violations"] if v["path"].startswith(PINS_FILE + ":")]
    print(f"t5_scope_check: base {rep['base'][:12]}: {rep['base_reason']}")
    print(f"  compared with {('commit ' + rep['head'][:12]) if args.head else 'the working tree'}")
    print(f"  {rep['n_changed']} changed paths, {len(rep['violations']) - len(pins_bad)} touching a frozen path")
    fence = rep["fence"]
    if fence["pins_reference"]:
        print(f"  pins file: {fence['pins_checked']} existing pins compared with {fence['pins_reference'][:12]}, "
              f"{fence['pins_reference_reason']}; {len(pins_bad)} changed or removed")
    for m in rep["base_lacks"]:
        print(f"  warning: the base has no {m}; is it a germquant commit at or after the T5 base?")
    for a in rep["allowed"]:
        print(f"  allowed: {a['path']} ({a['change']}: {a['why']}, {a['rule']})")
    for b in rep["violations"]:
        print(f"  FROZEN: {b['path']} {b['change']} ({b['why']}, rule {b['rule']})")
    if args.list:
        for f in sorted(files_at(repo, rep["base"], [*FROZEN_FOLDERS, *SEALED_FOLDERS]) | set(FROZEN_FILES)):
            print(f"  frozen at the base: {f}")
        print(f"  fence (G3): {PINS_FILE} (append-only), {FENCE_MODULE} (never deleted)")
    if args.json:
        Path(args.json).write_text(json.dumps(rep, indent=1), encoding="utf-8")
    print(rep["status"] + ("" if rep["status"] == "PASS" else
                           ": a tranche 5 diff may not modify or delete a frozen path or change an existing "
                           "pin (design 11.2, gate G3)"))
    return 0 if rep["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
