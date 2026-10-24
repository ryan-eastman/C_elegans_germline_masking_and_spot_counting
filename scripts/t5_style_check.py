"""Tranche 5 style check (gate G7): no en or em dash in the lines a diff adds, and pure ASCII where it is required.

    python scripts/t5_style_check.py [--base auto|REF] [--head REF] [--json OUT]

docs/DESIGN_system_packs.md, gate G7 (section 11.3: "no U+2013 or U+2014 in new files; ASCII in packs/, new
config files and project writes") and rule R9 (hashed text is pure ASCII); the Train B plan, decision D3. Two
rules:

  (a) no U+2013 (EN DASH) or U+2014 (EM DASH) in a line the diff adds, in any text file: code, docs,
      workflows and data files alike. A line that is not UTF-8 is also read as cp1252, the ANSI code page
      that Windows PowerShell 5.1's Set-Content and Add-Content write by default, where the two dashes are
      the bytes 0x96 and 0x97. A UTF-16 or UTF-32 text file (one that starts with its byte-order mark) holds
      NUL bytes, so git takes it for binary and its lines cannot be checked: it is reported instead, to be
      written as UTF-8;
  (b) pure ASCII in the text files the diff adds under config/, packs/, tests/fixtures/ and tests/golden/,
      and in every *.yaml or *.yml file it adds anywhere. A text file here is one whose suffix is on
      TEXT_SUFFIXES of germquant.packs.format (.yaml, .yml, .json, .md, .csv, .txt, .cff), read from that
      module's source, so this check needs only git and the standard library. Where the diff modifies such
      a file, the lines it adds must be ASCII.

Skipped: deleted files; binary files (git's own binary detection, which follows the binary rules of
.gitattributes, a NUL byte in the first 8000 bytes of an untracked file, and the suffixes of BINARY_SUFFIXES,
such as .tif, .ome.tif, .npz and .png), except the UTF-16 and UTF-32 text files of rule (a); an untracked
folder that holds a git repository of its own (git lists only its name); and the frozen R8 paths of
scripts/t5_scope_check.py (the RAD-51 seal), which this check never opens: every git command it runs on
paths excludes them by pathspec, and it never reads them from the working tree.

What is compared, and with which base, is as in scripts/t5_scope_check.py, whose base resolution this check
imports: without --head the working tree (committed, staged, unstaged and untracked files), with --head REF
the commit REF; --base auto (the default) takes the later of the T5 base recorded in
tests/golden/default_invariants.json and the merge base with origin/main. When those two bases have diverged,
a line counts only when it is added against both, so neither main's own changes nor the side commits are
charged to the branch. CI runs this check after the scope check in the T5 step of .github/workflows/ci.yml;
tests/test_t5_style_check.py also runs it on the checkout, so t5_gate's pytest step enforces G7 on the
workstation.

Exit codes, as the scope check's: 0 nothing found, 1 a finding, 2 the base or the head could not be read
(a git error, or a working-tree file that cannot be read).
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
# a line that is not UTF-8 is read in this code page too: Windows PowerShell 5.1's Set-Content and Add-Content
# write it by default, and it stores the two dashes as the single bytes 0x96 and 0x97
LEGACY_CODEC = "cp1252"
# the byte-order marks of UTF-16 and UTF-32 text, longest first: git takes such a file for binary (its ASCII
# characters hold NUL bytes), so rule (a) reports the file instead of its lines
WIDE_BOMS = ((b"\xff\xfe\x00\x00", "UTF-32 LE"), (b"\x00\x00\xfe\xff", "UTF-32 BE"),
             (b"\xff\xfe", "UTF-16 LE"), (b"\xfe\xff", "UTF-16 BE"))
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
RULES = {"a": "no U+2013 or U+2014 in a line the diff adds (any text file, also as cp1252 bytes; a UTF-16 or "
              "UTF-32 text file is reported, since git cannot diff it)",
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


def wide_encoding(start: bytes) -> str | None:
    """The UTF-16 or UTF-32 encoding whose byte-order mark ``start`` begins with, or None."""
    return next((name for bom, name in WIDE_BOMS if start.startswith(bom)), None)


def read_file(repo: Path, path: str, limit: int | None = None) -> bytes:
    """The bytes of a working-tree file, or its first ``limit`` bytes. Never called for an R8 path: the
    listings that name the paths exclude them."""
    with open(repo / path, "rb") as fh:
        return fh.read() if limit is None else fh.read(limit)


def _wide_texts(repo: Path, head: str | None, binary: set[str]) -> dict[str, str]:
    """{path: encoding} of the changed files that git takes for binary only by its NUL test (not by their
    suffix or a .gitattributes rule) and that start with a UTF-16 or UTF-32 byte-order mark: text that git
    cannot diff. The first bytes come from ``head``, or from the working tree when ``head`` is None."""
    sniff = [p for p in sorted(binary) if not binary_by_name(p)]
    by_rule = attr_binary(repo, sniff)
    out: dict[str, str] = {}
    for p in sniff:
        if p in by_rule:
            continue
        start = (_git_bytes(repo, "cat-file", "blob", f"{head}:{p}")[:4] if head is not None
                 else read_file(repo, p, 4))
        enc = wide_encoding(start)
        if enc:
            out[p] = enc
    return out


def tracked_changes(repo: Path, base: str, head: str | None, sniff: bool = True) -> tuple[
        dict[str, list[tuple[int, bytes]]], dict[str, str], dict[str, str], set[str]]:
    """(added lines per text file, {skipped path: why}, {UTF-16 or UTF-32 text file: encoding}, every changed
    path) between the base and the head, or the working tree when ``head`` is None. Deleted files and the R8
    paths are left out by git itself. ``sniff=False`` leaves the wide-text test out (every binary is skipped)."""
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
    wide = _wide_texts(repo, head, binary) if sniff else {}
    patch = _git_bytes(repo, "diff", "-U0", "--no-color", "--no-ext-diff", "--no-textconv", "--src-prefix=a/",
                       "--dst-prefix=b/", *common, *_diff_args(base, head), *spec)
    added_lines = {p: v for p, v in parse_patch(patch).items() if p not in binary}
    stray = sorted(set(added_lines) - changed)
    if stray:
        raise GitError(f"the patch names paths that the numstat listing does not: {stray[:3]}")
    return added_lines, {p: "binary" for p in binary if p not in wide}, wide, changed


def untracked_changes(repo: Path) -> tuple[dict[str, list[tuple[int, bytes]]], dict[str, str], dict[str, str],
                                           set[str]]:
    """(lines per text file, {skipped path: why}, {UTF-16 or UTF-32 text file: encoding}, every path) of the
    untracked files that git does not ignore. git lists an untracked folder that holds a repository of its
    own by its name and a final "/": it is skipped, never read."""
    paths = [p for p in _split_z(_git_bytes(repo, "ls-files", "--others", "--exclude-standard", "-z", "--", ".",
                                            *R8_EXCLUDES)) if not is_r8(p)]
    skipped = {p: "nested repository" for p in paths if p.endswith("/")}
    files = [p for p in paths if p not in skipped]
    binary = {p for p in files if binary_by_name(p)}
    binary |= attr_binary(repo, [p for p in files if p not in binary])
    lines: dict[str, list[tuple[int, bytes]]] = {}
    wide: dict[str, str] = {}
    for p in files:
        if p in binary:
            continue
        data = read_file(repo, p)
        if b"\0" in data[:SNIFF]:
            enc = wide_encoding(data)
            if enc:
                wide[p] = enc
            else:
                binary.add(p)
            continue
        rows = data.split(b"\n")
        if rows and rows[-1] == b"":
            rows.pop()
        lines[p] = [(n, row) for n, row in enumerate(rows, 1)]
    skipped.update((p, "binary") for p in binary)
    return lines, skipped, wide, set(paths)


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
    """The findings of rules (a) and (b) on one added line of ``path``. A line that is not UTF-8 is searched
    for the dashes as cp1252 too."""
    raw = raw[:-1] if raw.endswith(b"\r") else raw
    text = raw.decode("utf-8", "replace")
    try:
        raw.decode("utf-8")
        utf8 = True
        texts = ((text, ""),)
    except UnicodeDecodeError:
        utf8 = False
        texts = ((text, ""), (raw.decode(LEGACY_CODEC, "replace"),
                              f", as the {LEGACY_CODEC} byte of a non-UTF-8 line"))
    out = []
    for ch, name in DASHES.items():
        hit = next(((t, t.find(ch), note) for t, note in texts if ch in t), None)
        if hit:
            t, at, note = hit
            out.append({"path": path, "line": lineno, "rule": "a", "what": name + note, "text": _excerpt(t, at)})
    if ascii_rule_applies(path) and any(b > 127 for b in raw):
        at = next(k for k, ch in enumerate(text) if ord(ch) > 127)
        what = f"non-ASCII {_char_name(text[at])}" if utf8 else "non-ASCII bytes (not UTF-8)"
        out.append({"path": path, "line": lineno, "rule": "b", "what": what, "text": _excerpt(text, at)})
    return out


def wide_finding(path: str, encoding: str) -> dict:
    """Rule (a) on a UTF-16 or UTF-32 text file: git takes it for binary, so its lines cannot be checked."""
    bom = next(b for b, name in WIDE_BOMS if name == encoding)
    return {"path": path, "line": 1, "rule": "a",
            "what": f"{encoding} text, which git takes for binary, so its lines cannot be checked; write it as UTF-8",
            "text": "byte-order mark " + bom.hex(" ").upper()}


def check(repo: Path = REPO, base: str = "auto", head: str | None = None) -> dict:
    """Run the check; returns the report (``report["status"]`` is PASS or FAIL). Raises GitError, or OSError
    when a working-tree file cannot be read."""
    head_commit = scope._commit(repo, head or "HEAD")
    if head_commit is None:
        raise GitError(f"cannot resolve the head {head or 'HEAD'!r} to a commit")
    base_sha, why, second = scope.resolve_base(repo, base, head_commit, head)
    lines, skipped, wide, changed = tracked_changes(repo, base_sha, head)
    if head is None:
        u_lines, u_skipped, u_wide, u_paths = untracked_changes(repo)
        lines.update(u_lines)
        skipped.update(u_skipped)
        wide.update(u_wide)
        changed |= u_paths
    if second:
        # diverged bases: a line counts only when it is added against both (as a path does in the scope check)
        s_lines, _, _, s_changed = tracked_changes(repo, second, head, sniff=False)
        if head is None:
            s_lines.update(u_lines)
            s_changed |= u_paths
        lines = _intersect(lines, s_lines)
        skipped = {p: w for p, w in skipped.items() if p in s_changed}
        wide = {p: enc for p, enc in wide.items() if p in s_changed}
    findings = [f for path in sorted(lines) for n, raw in lines[path] for f in line_findings(path, n, raw)]
    findings += [wide_finding(p, enc) for p, enc in wide.items()]
    findings.sort(key=lambda f: (f["path"], f["line"]))      # stable: a line's own findings keep their order
    return {"status": "FAIL" if findings else "PASS", "base": base_sha, "base_reason": why, "second_base": second,
            "head": head_commit if head else "working tree", "rules": RULES,
            "n_files": len(lines), "n_lines": sum(len(v) for v in lines.values()),
            "files": sorted(lines), "findings": findings,
            "skipped": [{"path": p, "why": w} for p, w in sorted(skipped.items())],
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
    except (GitError, OSError) as e:
        print(f"t5_style_check: ERROR: {e}")
        return 2
    print(f"t5_style_check: base {rep['base'][:12]}: {rep['base_reason']}")
    print(f"  compared with {('commit ' + rep['head'][:12]) if args.head else 'the working tree'}")
    n_skipped = Counter(s["why"] for s in rep["skipped"])
    print(f"  {rep['n_lines']} added lines in {rep['n_files']} text files checked, {n_skipped['binary']} binary "
          f"files skipped" + (f", {n_skipped['nested repository']} nested repositories skipped"
                              if n_skipped["nested repository"] else "") + ", the R8 paths never opened")
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
