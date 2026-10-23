"""Tranche 5 style check (gate G7): no en or em dash in the lines a diff adds, and pure ASCII where it is required.

    python scripts/t5_style_check.py [--base auto|REF] [--head REF] [--json OUT]

docs/DESIGN_system_packs.md, gate G7 (section 11.3: "no U+2013 or U+2014 in new files; ASCII in packs/, new
config files and project writes") and rule R9 (hashed text is pure ASCII); the Train B plan, decision D3. Two
rules:

  (a) no U+2013 (EN DASH) or U+2014 (EM DASH) in a line the diff adds, in any text file: code, docs,
      workflows and data files alike;
  (b) pure ASCII in the text files the diff adds under config/, packs/, tests/fixtures/ and tests/golden/,
      and in every *.yaml or *.yml file it adds anywhere. A text file here is one whose suffix is on
      TEXT_SUFFIXES of germquant.packs.format (.yaml, .yml, .json, .md, .csv, .txt, .cff), read from that
      module's source, so this check needs only git and the standard library. Where the diff modifies such
      a file, the lines it adds must be ASCII.

Skipped: deleted files; binary files (git's own binary detection, which follows the binary rules of
.gitattributes, a NUL byte in the first 8000 bytes of an untracked file, and the suffixes of BINARY_SUFFIXES,
such as .tif, .ome.tif, .npz and .png); and the frozen R8 paths of scripts/t5_scope_check.py (the RAD-51
seal), which this check never opens: every git command it runs on paths excludes them by pathspec, and it
never reads them from the working tree.

What is compared, and with which base, is as in scripts/t5_scope_check.py, whose base resolution this check
imports: without --head the working tree (committed, staged, unstaged and untracked files), with --head REF
the commit REF; --base auto (the default) takes the later of the T5 base recorded in
tests/golden/default_invariants.json and the merge base with origin/main. When those two bases have diverged,
a line counts only when it is added against both, so neither main's own changes nor the side commits are
charged to the branch. CI runs this check after the scope check in the T5 step of .github/workflows/ci.yml;
tests/test_t5_style_check.py also runs it on the checkout, so t5_gate's pytest step enforces G7 on the
workstation.

Exit codes, as the scope check's: 0 nothing found, 1 a finding, 2 the base or the head could not be read.
"""
from __future__ import annotations

import argparse
import ast
import importlib.util
import json
import subprocess
import sys
import unicodedata
from collections import Counter
from pathlib import Path, PurePosixPath

try:
    import t5_scope_check as scope          # run as a script: scripts/ is on sys.path
except ModuleNotFoundError:                 # loaded by file path (the tests): the scope check sits next to it
    _spec = importlib.util.spec_from_file_location(
        "t5_scope_check", Path(__file__).resolve().with_name("t5_scope_check.py"))
    scope = importlib.util.module_from_spec(_spec)
    sys.modules["t5_scope_check"] = scope
    _spec.loader.exec_module(scope)

GitError = scope.GitError
REPO = Path(__file__).resolve().parents[1]
FORMAT_SOURCE = REPO / "src" / "germquant" / "packs" / "format.py"
DASHES = {chr(0x2013): "U+2013 EN DASH", chr(0x2014): "U+2014 EM DASH"}   # by code point: no dash in this file
ASCII_FOLDERS = ("config/", "packs/", "tests/fixtures/", "tests/golden/")   # rule (b), for TEXT_SUFFIXES files
ASCII_EVERYWHERE = (".yaml", ".yml")                                        # rule (b), wherever they are
# skipped as binary whatever git says: the binary rules of .gitattributes (.tif, .tiff, .png, .npz, .zip, .h5;
# .ome.tif ends in .tif) and the other image, array, model, archive and office formats of this repository
BINARY_SUFFIXES = (".tif", ".tiff", ".png", ".npz", ".zip", ".h5", ".npy", ".nd2", ".ims", ".czi", ".lif",
                   ".jpg", ".jpeg", ".gif", ".bmp", ".pdf", ".pt", ".pth", ".pkl", ".gz", ".parquet",
                   ".xlsx", ".xls", ".docx", ".pptx")
SNIFF = 8000                        # the bytes git looks at for a NUL when it decides a file is binary
R8_FILES = tuple(sorted(p for p, rule in scope.FROZEN_FILES.items() if rule == "R8"))
R8_FOLDERS = tuple(sorted(f for f, rule in scope.SEALED_FOLDERS.items() if rule == "R8"))
R8_EXCLUDES = tuple([f":(top,exclude,literal){p}" for p in R8_FILES] +
                    [f":(top,exclude,literal){f}" for f in R8_FOLDERS])
RULES = {"a": "no U+2013 or U+2014 in a line the diff adds (any text file)",
         "b": "pure ASCII in text files added under " + ", ".join(ASCII_FOLDERS) +
              " and in every added *.yaml or *.yml"}


def _text_suffixes(source: Path = FORMAT_SOURCE) -> tuple[str, ...]:
    """TEXT_SUFFIXES of germquant.packs.format, read from its source (no germquant import)."""
    tree = ast.parse(source.read_text(encoding="utf-8"), filename=str(source))
    for node in tree.body:
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, ast.AnnAssign):
            targets = [node.target]
        else:
            continue
        if node.value is not None and any(isinstance(t, ast.Name) and t.id == "TEXT_SUFFIXES" for t in targets):
            return tuple(ast.literal_eval(node.value))
    raise RuntimeError(f"no TEXT_SUFFIXES assignment in {source}")


TEXT_SUFFIXES = _text_suffixes()


def is_r8(path: str) -> bool:
    return path in R8_FILES or path.startswith(R8_FOLDERS)


def _suffix(path: str) -> str:
    return PurePosixPath(path).suffix.lower()


def binary_by_name(path: str) -> bool:
    return _suffix(path) in BINARY_SUFFIXES


def ascii_rule_applies(path: str) -> bool:
    """Rule (b): a text file under one of ASCII_FOLDERS, or a YAML file anywhere."""
    suffix = _suffix(path)
    return suffix in ASCII_EVERYWHERE or (suffix in TEXT_SUFFIXES and path.startswith(ASCII_FOLDERS))


def _git_bytes(repo: Path, *args: str) -> bytes:
    """Run git in ``repo`` and return its raw output; a failure raises GitError."""
    r = subprocess.run(["git", "-c", "core.quotePath=false", *args], cwd=repo, capture_output=True)
    if r.returncode != 0:
        msg = (r.stderr or r.stdout).decode("utf-8", "replace").strip()
        raise GitError(f"git {' '.join(args)}: {msg}")
    return r.stdout


def _split_z(raw: bytes) -> list[str]:
    return [t.decode("utf-8", "surrogateescape") for t in raw.split(b"\0") if t]


_ESCAPES = {ord("a"): 7, ord("b"): 8, ord("t"): 9, ord("n"): 10, ord("v"): 11, ord("f"): 12, ord("r"): 13,
            ord('"'): 34, ord("\\"): 92}


def _unquote(raw: bytes) -> bytes:
    """A path as git prints it in a patch header: C-style quoted when it holds a special character."""
    if not (raw.startswith(b'"') and raw.endswith(b'"') and len(raw) >= 2):
        return raw
    body, out, i = raw[1:-1], bytearray(), 0
    while i < len(body):
        c = body[i]
        if c == 92 and i + 1 < len(body):                       # a backslash escape
            nxt = body[i + 1]
            octal = body[i + 1:i + 4]
            if len(octal) == 3 and all(48 <= d <= 55 for d in octal):      # \ooo, one byte
                out.append(int(octal, 8))
                i += 4
                continue
            out.append(_ESCAPES.get(nxt, nxt))
            i += 2
            continue
        out.append(c)
        i += 1
    return bytes(out)


def _header_path(raw: bytes) -> str | None:
    """The path of a "+++ " header line (without the "+++ "), or None for /dev/null."""
    if raw.endswith(b"\t"):                     # git ends a name that holds a space with a tab
        raw = raw[:-1]
    raw = _unquote(raw)
    if raw == b"/dev/null":
        return None
    if not raw.startswith(b"b/"):
        raise GitError(f"unexpected patch header path {raw!r}")
    return raw[2:].decode("utf-8", "surrogateescape")


def _hunk_counts(line: bytes) -> tuple[int, int, int]:
    """(old count, new start, new count) of a "@@ -a,b +c,d @@" hunk header."""
    try:
        old, new = line.split(b" ")[1:3]
        old_n = int(old.split(b",")[1]) if b"," in old else 1
        start, _, n = new[1:].partition(b",")
        return old_n, int(start), int(n) if n else 1
    except (ValueError, IndexError) as e:
        raise GitError(f"unreadable hunk header {line!r}") from e


def parse_patch(raw: bytes) -> dict[str, list[tuple[int, bytes]]]:
    """{path: [(line number in the new file, added line)]} from `git diff -U0` output. The hunk counts say
    how many lines each hunk holds, so an added line that starts with "++" is never read as a header."""
    lines = raw.split(b"\n")
    out: dict[str, list[tuple[int, bytes]]] = {}
    path: str | None = None
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(b"diff --git "):
            path = None
            i += 1
            while i < len(lines) and not lines[i].startswith((b"@@ ", b"diff --git ")):
                if lines[i].startswith(b"+++ "):
                    path = _header_path(lines[i][4:])
                i += 1
            continue
        if line.startswith(b"@@ ") and path is not None:
            old_n, lineno, new_n = _hunk_counts(line)
            i += 1
            seen_old = seen_new = 0
            while i < len(lines) and (seen_old < old_n or seen_new < new_n or lines[i].startswith(b"\\")):
                body = lines[i]
                if body.startswith(b"+"):
                    out.setdefault(path, []).append((lineno, body[1:]))
                    lineno += 1
                    seen_new += 1
                elif body.startswith(b"-"):
                    seen_old += 1
                elif body.startswith(b" "):                   # context lines (none with -U0)
                    lineno += 1
                    seen_old += 1
                    seen_new += 1
                elif not body.startswith(b"\\"):              # "\ No newline at end of file" is skipped
                    break
                i += 1
            continue
        i += 1
    return out


def _diff_args(base: str, head: str | None) -> list[str]:
    return [base] + ([head] if head is not None else [])


def tracked_changes(repo: Path, base: str, head: str | None) -> tuple[dict[str, list[tuple[int, bytes]]],
                                                                      set[str], set[str]]:
    """(added lines per text file, binary files, every changed path) between the base and the head, or the
    working tree when ``head`` is None. Deleted files and the R8 paths are left out by git itself."""
    common = ["--no-renames", "--ignore-submodules", "--diff-filter=d"]
    spec = ["--", ".", *R8_EXCLUDES]
    numstat = _split_z(_git_bytes(repo, "diff", "--numstat", "-z", *common, *_diff_args(base, head), *spec))
    changed: set[str] = set()
    binary: set[str] = set()
    for entry in numstat:
        added, _, rest = entry.partition("\t")
        _, _, path = rest.partition("\t")
        changed.add(path)
        if added == "-" or binary_by_name(path):
            binary.add(path)
    patch = _git_bytes(repo, "diff", "-U0", "--no-color", "--no-ext-diff", "--no-textconv", "--src-prefix=a/",
                       "--dst-prefix=b/", *common, *_diff_args(base, head), *spec)
    added_lines = {p: v for p, v in parse_patch(patch).items() if p not in binary}
    stray = sorted(set(added_lines) - changed)
    if stray:
        raise GitError(f"the patch names paths that the numstat listing does not: {stray[:3]}")
    return added_lines, binary, changed


def read_untracked(repo: Path, path: str) -> bytes:
    """The bytes of an untracked file (never called for an R8 path)."""
    return (repo / path).read_bytes()


def untracked_changes(repo: Path) -> tuple[dict[str, list[tuple[int, bytes]]], set[str], set[str]]:
    """(lines per text file, binary files, every path) of the untracked files that git does not ignore."""
    paths = [p for p in _split_z(_git_bytes(repo, "ls-files", "--others", "--exclude-standard", "-z", "--", ".",
                                            *R8_EXCLUDES)) if not is_r8(p)]
    binary = {p for p in paths if binary_by_name(p)}
    binary |= attr_binary(repo, [p for p in paths if p not in binary])
    lines: dict[str, list[tuple[int, bytes]]] = {}
    for p in paths:
        if p in binary:
            continue
        data = read_untracked(repo, p)
        if b"\0" in data[:SNIFF]:
            binary.add(p)
            continue
        rows = data.split(b"\n")
        if rows and rows[-1] == b"":
            rows.pop()
        lines[p] = [(n, row) for n, row in enumerate(rows, 1)]
    return lines, binary, set(paths)


def attr_binary(repo: Path, paths: list[str]) -> set[str]:
    """The paths that a binary rule of .gitattributes covers (the binary macro unsets the diff attribute)."""
    if not paths:
        return set()
    r = subprocess.run(["git", "check-attr", "-z", "--stdin", "diff"], cwd=repo, capture_output=True,
                       input=b"".join(p.encode("utf-8", "surrogateescape") + b"\0" for p in paths))
    if r.returncode != 0:
        raise GitError(f"git check-attr: {(r.stderr or r.stdout).decode('utf-8', 'replace').strip()}")
    out = _split_z(r.stdout)                    # path, attribute, value, path, attribute, value, ...
    return {out[k] for k in range(0, len(out) - 2, 3) if out[k + 2] == "unset"}


def _intersect(first: dict[str, list[tuple[int, bytes]]],
               second: dict[str, list[tuple[int, bytes]]]) -> dict[str, list[tuple[int, bytes]]]:
    """The lines of ``first`` that ``second`` adds too (per path, counted as a multiset of line texts)."""
    out: dict[str, list[tuple[int, bytes]]] = {}
    for path, rows in first.items():
        pool = Counter(text for _, text in second.get(path, []))
        kept = []
        for n, text in rows:
            if pool[text] > 0:
                pool[text] -= 1
                kept.append((n, text))
        if kept:
            out[path] = kept
    return out


def _excerpt(text: str, at: int, width: int = 70) -> str:
    """ASCII excerpt of a line around position ``at`` (non-ASCII shown as escapes)."""
    lo = max(0, at - width // 2)
    piece = text[lo:lo + width].strip()
    return piece.encode("ascii", "backslashreplace").decode("ascii")


def _char_name(ch: str) -> str:
    return f"U+{ord(ch):04X} {unicodedata.name(ch, 'UNNAMED')}"


def line_findings(path: str, lineno: int, raw: bytes) -> list[dict]:
    """The findings of rules (a) and (b) on one added line of ``path``."""
    raw = raw[:-1] if raw.endswith(b"\r") else raw
    text = raw.decode("utf-8", "replace")
    out = []
    for ch, name in DASHES.items():
        at = text.find(ch)
        if at >= 0:
            out.append({"path": path, "line": lineno, "rule": "a", "what": name, "text": _excerpt(text, at)})
    if ascii_rule_applies(path) and any(b > 127 for b in raw):
        at = next(k for k, ch in enumerate(text) if ord(ch) > 127)
        try:
            raw.decode("utf-8")
            what = f"non-ASCII {_char_name(text[at])}"
        except UnicodeDecodeError:
            what = "non-ASCII bytes (not UTF-8)"
        out.append({"path": path, "line": lineno, "rule": "b", "what": what, "text": _excerpt(text, at)})
    return out


def check(repo: Path = REPO, base: str = "auto", head: str | None = None) -> dict:
    """Run the check; returns the report (``report["status"]`` is PASS or FAIL). Raises GitError."""
    head_commit = scope._commit(repo, head or "HEAD")
    if head_commit is None:
        raise GitError(f"cannot resolve the head {head or 'HEAD'!r} to a commit")
    base_sha, why, second = scope.resolve_base(repo, base, head_commit, head)
    lines, binary, changed = tracked_changes(repo, base_sha, head)
    if head is None:
        u_lines, u_binary, u_paths = untracked_changes(repo)
        lines.update(u_lines)
        binary |= u_binary
        changed |= u_paths
    if second:
        # diverged bases: a line counts only when it is added against both (as a path does in the scope check)
        s_lines, _, s_changed = tracked_changes(repo, second, head)
        if head is None:
            s_lines.update(u_lines)
            s_changed |= u_paths
        lines = _intersect(lines, s_lines)
        binary = {p for p in binary if p in s_changed}
    findings = [f for path in sorted(lines) for n, raw in lines[path] for f in line_findings(path, n, raw)]
    return {"status": "FAIL" if findings else "PASS", "base": base_sha, "base_reason": why, "second_base": second,
            "head": head_commit if head else "working tree", "rules": RULES,
            "n_files": len(lines), "n_lines": sum(len(v) for v in lines.values()),
            "files": sorted(lines), "findings": findings,
            "skipped": [{"path": p, "why": "binary"} for p in sorted(binary)],
            "never_opened": {"files": list(R8_FILES), "folders": list(R8_FOLDERS), "rule": "R8"},
            "text_suffixes": list(TEXT_SUFFIXES), "ascii_folders": list(ASCII_FOLDERS)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default="auto", help="base commit or ref (default: auto, as the scope check)")
    ap.add_argument("--head", help="check this commit instead of the working tree")
    ap.add_argument("--repo", default=str(REPO), help=argparse.SUPPRESS)
    ap.add_argument("--json", help="write the report here")
    args = ap.parse_args(argv)
    try:
        rep = check(Path(args.repo), args.base, args.head)
    except GitError as e:
        print(f"t5_style_check: ERROR: {e}")
        return 2
    print(f"t5_style_check: base {rep['base'][:12]}: {rep['base_reason']}")
    print(f"  compared with {('commit ' + rep['head'][:12]) if args.head else 'the working tree'}")
    print(f"  {rep['n_lines']} added lines in {rep['n_files']} text files checked, {len(rep['skipped'])} binary "
          f"files skipped, the R8 paths never opened")
    for f in rep["findings"]:
        print(f"  FOUND: {f['path']}:{f['line']}: ({f['rule']}) {f['what']}: {f['text']}")
    if args.json:
        Path(args.json).write_text(json.dumps(rep, indent=1), encoding="utf-8")
    print(rep["status"] + ("" if rep["status"] == "PASS" else
                           f": {len(rep['findings'])} findings; rule (a): {RULES['a']}; rule (b): {RULES['b']} "
                           f"(gate G7)"))
    return 0 if rep["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
