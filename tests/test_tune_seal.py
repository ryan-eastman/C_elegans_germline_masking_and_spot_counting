"""The tune seal and import order (design 6.7 and 2, rules R5 and R8; plan D28; gate G6; P5 acceptance 11).

* AST: germquant.tune (and the tune-spots command module) imports, of the benchmark and validation code, only
  germquant.benchmark.stats, matching, nuclei3d and spots and germquant.validate.compare, never through
  importlib or __import__, and its sources never name the split file of the RAD-51 benchmark;
* every tuner imports germquant.pipeline (torch) as its first statement (D28);
* in a fresh interpreter: importing germquant.tune loads no torch, CuPy, SpotMAX, Cellpose, bioio or dask.array;
  a tiny LoG tune and a tiny SpotMAX tune leave germquant.benchmark.detector (and tuning and logo) out of
  sys.modules, load germquant.pipeline and torch before any spot backend, and import torch before CuPy,
  SpotMAX, dask and bioio. In this suite's own process the detector is already loaded (tests/test_bench_tuning.py
  imports it at collection), so these checks run only in child interpreters.
"""
from __future__ import annotations

import ast
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest
import tune_helpers as H
from packs_helpers import HEAVY, PYTHON, ROOT, subprocess_env

TUNE = ROOT / "src" / "germquant" / "tune"
SOURCES = sorted(TUNE.glob("*.py")) + [ROOT / "src" / "germquant" / "cli_tune.py"]
ALLOWED = {"germquant.benchmark.stats", "germquant.benchmark.matching", "germquant.benchmark.nuclei3d",
           "germquant.benchmark.spots", "germquant.validate.compare"}
SEALED = ("germquant.benchmark", "germquant.validate")
READERS = ("bioio", "dask.array")


def _split_file_name() -> str:
    """The stem of the RAD-51 split record, taken from the scope check's frozen list, so this file does not
    hold the name itself (rule R8)."""
    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        import t5_scope_check as scope
    finally:
        sys.path.pop(0)
    names = [Path(p).stem for p, rule in scope.FROZEN_FILES.items() if rule == "R8" and p.endswith(".json")]
    assert len(names) == 1, names
    return names[0]


def _targets(tree: ast.AST, package: str):
    """(line, module) for every import in ``tree``: an ``import a.b`` gives a.b; ``from a import b`` gives a.b
    for each name (b may be a submodule) and a itself; relative imports are resolved against ``package``."""
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield node.lineno, alias.name
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = package.split(".")[: len(package.split(".")) - node.level + 1]
                mod = ".".join(base + ([node.module] if node.module else []))
            else:
                mod = node.module or ""
            yield node.lineno, mod
            for alias in node.names:
                yield node.lineno, f"{mod}.{alias.name}"


def _sealed_problems(path: Path, package: str | None = None) -> list[str]:
    package = package or ("germquant.tune" if path.parent == TUNE else "germquant")
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    problems = []
    seen_from: dict[int, set] = {}
    for line, target in _targets(tree, package):
        seen_from.setdefault(line, set()).add(target)
    for line, targets in seen_from.items():
        for target in targets:
            if not any(target == s or target.startswith(s + ".") for s in SEALED):
                continue
            if target in SEALED:
                # "from germquant.benchmark import spots" names the package and the submodule; the package alone
                # (import germquant.benchmark, from germquant import benchmark) is not allowed
                if not any(t.startswith(target + ".") for t in targets):
                    problems.append(f"{path.name}:{line} imports the package {target}")
                continue
            module = next((a for a in ALLOWED if target == a or target.startswith(a + ".")), None)
            if module is None:
                problems.append(f"{path.name}:{line} imports {target}")
    for node in ast.walk(tree):
        name = node.id if isinstance(node, ast.Name) else node.attr if isinstance(node, ast.Attribute) else None
        if name in ("import_module", "__import__"):
            problems.append(f"{path.name}:{node.lineno} imports by name ({name})")
    return problems


def test_ast_allow_list(tmp_path):
    assert len(SOURCES) >= 7 and all(p.is_file() for p in SOURCES), SOURCES
    problems = [msg for p in SOURCES for msg in _sealed_problems(p)]
    assert not problems, problems
    # the checker itself catches every forbidden form (probes in tmp_path, read as modules of germquant.tune)
    bad = ("from germquant.benchmark import detector\n", "import germquant.benchmark.tuning\n",
           "from ..benchmark import logo\n", "from ..validate import benchmark_split\n",
           "from germquant import benchmark\n", "import importlib\nimportlib.import_module('x')\n",
           "from germquant.benchmark.detector import score\n", "import germquant.validate\n")
    good = ("from germquant.benchmark import spots, stats\n", "from ..benchmark.spots import match_spots\n",
            "from germquant.validate.compare import agreement_stats\n", "from .grid import check_grid\n")
    for i, text in enumerate(bad + good):
        p = tmp_path / f"probe_{i}.py"
        p.write_text(text, encoding="ascii")
        assert bool(_sealed_problems(p, "germquant.tune")) == (text in bad), text


def test_no_split_file_name_in_tune_sources():
    needle = _split_file_name()
    for p in SOURCES:
        assert needle not in p.read_text(encoding="utf-8"), p


def test_every_tuner_imports_the_pipeline_first():
    """D28: the functions that read images or load backends begin with ``import germquant.pipeline``."""
    tree = ast.parse((TUNE / "spots.py").read_text(encoding="utf-8"))
    funcs = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
    for name in ("tune_spots", "_tune"):
        body = [s for s in funcs[name].body if not (isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant))]
        first = body[0]
        assert isinstance(first, ast.Import) and [a.name for a in first.names] == ["germquant.pipeline"], name
    for name in ("tune_spotmax", "tune_grid"):            # thin wrappers: everything happens in _tune
        body = [s for s in funcs[name].body if not (isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant))]
        assert len(body) == 1 and isinstance(body[0], ast.Return), name
        call = body[0].value
        assert isinstance(call, ast.Call) and isinstance(call.func, ast.Name) and call.func.id == "_tune", name


_PROBE = ("import sys, atexit\n"
          "atexit.register(lambda: sys.stderr.write('HEAVY=' + ','.join(m for m in {heavy!r} if m in sys.modules)"
          " + '\\n'))\n")


def _child(code: str, heavy=HEAVY, args=()) -> tuple[str, subprocess.CompletedProcess]:
    r = subprocess.run([PYTHON, "-c", _PROBE.format(heavy=tuple(heavy)) + code, *args], capture_output=True,
                       text=True, timeout=600, env=subprocess_env(), cwd=ROOT)
    m = re.search(r"HEAVY=(.*)", r.stderr)
    assert m, f"no probe line; rc={r.returncode} stderr={r.stderr[-3000:]}"
    return m.group(1).strip(), r


def test_import_tune_stays_light():
    code = ("import germquant.tune, germquant.tune.points, germquant.tune.grid, germquant.tune.apply\n"
            "import germquant.cli_tune\n"
            "assert 'germquant.pipeline' not in sys.modules and 'germquant.benchmark' not in sys.modules\n")
    heavy, r = _child(code, heavy=HEAVY + READERS)
    assert r.returncode == 0, r.stderr[-3000:]
    assert heavy == "", f"importing germquant.tune loaded {heavy}"


def test_fresh_interpreter_tune_keeps_the_detector_out_and_imports_torch_first(tmp_path):
    pytest.importorskip("bioio")
    pytest.importorskip("torch")
    # a small planted field as an OME-TIFF for the child, which reads it through bioio; the finished run of the
    # same pixels is made here through the .nd2 reader seam, and the truth points are saved on it
    field, seed, shape = H.TINY
    pix, truth, _ = H.planted_pixels(field, seed, shape)
    image = H.write_ome_tiff(tmp_path / "tiny.ome.tif", pix)
    stack = H.planted_stack("tiny")[0]
    assert (stack.data == pix).all()
    rd = tmp_path / "results" / "tiny"
    with pytest.MonkeyPatch.context() as mp:
        H.patch_planted_reader(mp, {"tiny.nd2": stack})
        H.finished_run(tmp_path / "tiny.nd2", H.basic_recipe(), rd)
    points = H.write_truth(tmp_path / "tiny__points_ryan.csv", truth, rd, "tiny")
    item = {"image_id": "tiny", "image": str(image), "result_dir": str(rd), "points": str(points)}
    have_spotmax = __import__("importlib").util.find_spec("spotmax") is not None
    code = ("import json\n"
            "import germquant.spots.detect as D\n"
            "import germquant.tune as T\n"
            "assert 'torch' not in sys.modules and 'germquant.pipeline' not in sys.modules\n"
            "loads = []\n"
            "def wrap(name, _orig=D.get_backend):\n"
            "    loads.append([name, 'germquant.pipeline' in sys.modules, 'torch' in sys.modules])\n"
            "    return _orig(name)\n"
            "D.get_backend = wrap\n"
            "item, out, sm = json.loads(sys.argv[1]), sys.argv[2], sys.argv[3] == '1'\n"
            "best = [T.tune_grid([item], {'spots.log.min_snr': [2.0, 3.0]}, backend='log', out_dir=out + '/log')['best']]\n"
            "if sm:\n"
            "    best.append(T.tune_spotmax([item], {'spots.effect_size_min': [0.0, 29.0]}, out_dir=out + '/sm')['best'])\n"
            "order = list(sys.modules)\n"
            "print(json.dumps({'loads': loads, 'best': best,\n"
            "                  'sealed': [m for m in ('germquant.benchmark.detector', 'germquant.benchmark.tuning',\n"
            "                                         'germquant.benchmark.logo') if m in sys.modules],\n"
            "                  'pos': {m: order.index(m) for m in ('torch', 'germquant.pipeline', 'cupy', 'spotmax',\n"
            "                                                     'cellpose', 'dask', 'bioio') if m in sys.modules}}))\n")
    _, r = _child(code, args=(json.dumps(item), str(tmp_path / "tuning"), "1" if have_spotmax else "0"))
    assert r.returncode == 0, r.stderr[-3000:]
    out = json.loads(r.stdout.strip().splitlines()[-1])
    assert out["sealed"] == [], f"a tune loaded {out['sealed']}"
    assert out["best"][0]["keys"] == {"spots.log.min_snr": 3.0}
    if have_spotmax:
        assert out["best"][1]["keys"] == {"spots.effect_size_min": 29.0}
    names = [name for name, _, _ in out["loads"]]
    assert names[:2] == ["log", "log"] and (not have_spotmax or names[-1] == "spotmax"), names
    assert all(pipe and torch for _, pipe, torch in out["loads"]), out["loads"]     # the pipeline came first
    pos = out["pos"]
    assert all(pos["torch"] < p for m, p in pos.items() if m not in ("torch", "germquant.pipeline")), pos
    assert all(pos["germquant.pipeline"] < pos[m] for m in ("bioio", "dask", "cupy", "spotmax") if m in pos), pos
