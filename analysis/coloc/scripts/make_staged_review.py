"""Build the pre-commit review page (.staged_review.html, untracked, in .git/info/exclude) for the
C_elegans_germline_masking_and_spot_counting repo from whatever is currently staged. Nothing is
committed here. The proposed commit message lives in coloc_analysis/COMMIT_MSG.txt (git commit -F)."""
import html
import os
import re
import subprocess

REPO = r"C:/Users/ryane/C_elegans_germline_masking_and_spot_counting"
OUT = os.path.join(REPO, ".staged_review.html")
MSG = r"C:/Users/ryane/coloc_analysis/COMMIT_MSG.txt"
NOTE = ("Nothing has been committed. Branch <b>{branch}</b>. {n} files staged. Since the last review: (1) the "
        "whole-gonad worker rerun with the rotation-null fix finished; only the <code>rot</code> control column moved. "
        "(2) A per-label mask audit (README, section 'Mask audit') found sperm / debris and somatic nuclei inside the "
        "whole-gonad nucleus set; the hand-traced stage results are unaffected. The audit tools, a three-way "
        "sensitivity table (<code>results/pc_filter_compare.csv</code>) and the pooled-pachytene columns / figure 1b are "
        "staged as ADDITIONS; existing figures and tables are unchanged. OPEN DECISION: which whole-gonad version is "
        "primary (recommended: pooled hand-traced pachytene, figure 1b). (3) Grant pipeline figure under "
        "<code>figures/grant</code>.")


def git(*a):
    return subprocess.run(["git", *a], cwd=REPO, capture_output=True, text=True, encoding="utf-8",
                          errors="replace").stdout

branch = git("rev-parse", "--abbrev-ref", "HEAD").strip()
files = [ln.split(chr(9)) for ln in git("diff", "--cached", "--name-status").splitlines() if ln.strip()]
numstat = {ln.split(chr(9))[2]: ln.split(chr(9))[:2] for ln in git("diff", "--cached", "--numstat").splitlines() if ln.strip()}
binaries = {p for p, (a, d) in numstat.items() if a == "-" or p.lower().endswith(".svg")}   # svg: show, do not diff

if os.path.exists(MSG):
    msg = open(MSG, encoding="utf-8").read()
else:                                   # first run: lift the message out of the previous page
    old = open(OUT, encoding="utf-8").read()
    msg = html.unescape(re.search(r"<h2>Proposed commit message</h2><pre>(.*?)</pre>", old, re.S).group(1))
    open(MSG, "w", encoding="utf-8").write(msg)

STAT = {"A": "add", "M": "mod", "D": "del", "R": "mod"}
rows = []
for f in files:
    st, path = f[0][0], f[-1]
    kind = "binary" if path in binaries else "code"
    rows.append(f'<tr class="{STAT.get(st, "mod")}"><td>{st}</td><td><a href="#f{len(rows)}">{html.escape(path)}</a></td><td>{kind}</td></tr>')

def colour(line):
    e = html.escape(line)
    if line.startswith("+++") or line.startswith("---"):
        return f'<span class="f">{e}</span>'
    if line.startswith("@@"):
        return f'<span class="k">{e}</span>'
    if line.startswith("+"):
        return f'<span class="a">{e}</span>'
    if line.startswith("-"):
        return f'<span class="r">{e}</span>'
    return f'<span class="h">{e}</span>' if line.startswith("diff ") or line.startswith("index ") else e

parts = []
for i, f in enumerate(files):
    st, path = f[0][0], f[-1]
    parts.append(f'<h3 id="f{i}">{st} {html.escape(path)}</h3>')
    if path in binaries:
        if path.lower().endswith((".png", ".svg")):
            parts.append(f'<img src="{html.escape(path)}" style="max-width:100%;border:1px solid #ddd">')
        else:
            parts.append('<p class="h">binary file</p>')
        continue
    d = git("diff", "--cached", "--no-color", "--", path)
    lines = d.splitlines()
    cap = 200 if path.lower().endswith((".csv", ".json")) else 1500     # data files: head only
    if len(lines) > cap:
        lines = lines[:cap] + [f"... {len(d.splitlines()) - cap} more lines not shown"]
    parts.append("<pre>" + chr(10).join(colour(ln) for ln in lines) + "</pre>")

page = ["<!doctype html><html><head><meta charset=\"utf-8\"><title>Staged commit review</title>",
        "<style>body{font-family:Segoe UI,Arial,sans-serif;margin:24px;max-width:1400px} pre{background:#f6f8fa;padding:12px;overflow-x:auto;font-size:12px;line-height:1.35}",
        "table{border-collapse:collapse;font-size:13px} td{padding:2px 10px;border-bottom:1px solid #eee} tr.add td:first-child{color:#1a7f37;font-weight:bold} tr.mod td:first-child{color:#9a6700;font-weight:bold} tr.del td:first-child{color:#cf222e;font-weight:bold}",
        ".a{color:#1a7f37} .r{color:#cf222e} .h{color:#6e7781} .f{color:#0969da;font-weight:bold} .k{color:#8250df} h2{margin-top:32px} h3{margin-top:28px;font-family:Consolas,monospace;font-size:13px} .note{background:#fff8c5;padding:10px 14px;border-left:4px solid #d4a72c}</style></head><body>",
        "<h1>Staged commit: review before committing</h1>",
        f'<div class="note">{NOTE.format(branch=branch, n=len(files))}</div>',
        f"<h2>Proposed commit message</h2><pre>{html.escape(msg)}</pre>",
        f"<h2>Files ({len(files)})</h2><table><tr><th>status</th><th>path</th><th>type</th></tr>{''.join(rows)}</table>",
        "<h2>Diffs</h2>", *parts, "</body></html>"]
open(OUT, "w", encoding="utf-8").write(chr(10).join(page))
print(f"wrote {OUT}: {len(files)} files, {os.path.getsize(OUT) // 1024} KB, branch {branch}")
