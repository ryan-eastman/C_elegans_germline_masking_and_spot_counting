"""Import guards and entry points (design 3.8, 3.9, 3.13 and gate G6; P1 acceptance 8, 10 and 12).

* In a fresh interpreter, germquant.packs.list_packs, a shallow germquant.packs.pack_status and
  ``python -m germquant.packs list --json`` leave torch, cupy, spotmax and cellpose out of sys.modules; a
  malformed pack under $GERMQUANT_PACKS is listed as unusable, and a generic-profile run still resolves cpsam.
  So do ``python -m germquant.packs derive``, ``list --json --project DIR`` and ``--version --json`` (Train B,
  P5-DERIVE), which also leave bioio and dask.array out, and a pack status of a recipe that needs a post-stage
  (Train B, P5b), which reads the post-stage from metadata and never imports it.
* ``python -m germquant --help`` works and imports torch before anything else of germquant.
* ``germquant --version`` and ``python -m germquant.packs --version`` print the engine version, git sha and
  python.
* With the engine installed non-editable inside another project (a copy of the package in a fake
  site-packages under a folder with its own pyproject.toml and config/), ``pack list`` names the missing
  built-in root and both fixes, and ``extends: generic.yaml`` from outside the checkout raises instead of
  binding to that project's or the working directory's config/ folder.
"""
from __future__ import annotations

import importlib.util
import json
import os
import re
import shutil
import subprocess

import pytest
from packs_helpers import FIXTURE_PACKS, HEAVY, PYTHON, ROOT, pack_env, subprocess_env  # noqa: F401

from germquant import cli, packs, pipeline
from germquant.config import load_config

_PROBE = ("import sys, atexit\n"
          "atexit.register(lambda: sys.stderr.write('HEAVY=' + ','.join(m for m in {heavy!r} if m in sys.modules)"
          " + '\\n'))\n")
# reading an image header through bioio imports dask.array, and dask imports CuPy whenever CuPy is installed:
# CI has no CuPy, so the image probes also refuse bioio and dask.array themselves
READERS = ("bioio", "dask.array")
EXAMPLE = FIXTURE_PACKS / "test_pack" / "example" / "test_pack_synthetic.ome.tif"


def _heavy_after(code: str, heavy=HEAVY, **env) -> tuple[str, subprocess.CompletedProcess]:
    full = _PROBE.format(heavy=tuple(heavy)) + code
    r = subprocess.run([PYTHON, "-c", full], capture_output=True, text=True, timeout=300,
                       env=subprocess_env(**env), cwd=ROOT)
    m = re.search(r"HEAVY=(.*)", r.stderr)
    assert m, f"no probe line; rc={r.returncode} stderr={r.stderr[-2000:]}"
    return m.group(1).strip(), r


def test_packs_api_stays_free_of_torch_in_a_fresh_interpreter(tmp_path):
    broken = tmp_path / "env_packs" / "broken" / "0.1.0"
    broken.mkdir(parents=True)
    (broken / "pack.yaml").write_text("format: 1\nid: [\n", encoding="ascii")
    paths = f"{FIXTURE_PACKS}{os.pathsep}{tmp_path / 'env_packs'}"
    code = ("import germquant.packs as P\n"
            "rows = P.list_packs(all=True)\n"
            "assert any(r['state'] == 'unusable' and r['id'] == 'broken' for r in rows), rows\n"
            "st = P.pack_status('test_pack')\n"
            "assert st['usable'], st\n"
            "st2 = P.pack_status('test_pack', recipe='with_model')\n"
            "assert not st2['usable']\n")
    heavy, r = _heavy_after(code, GERMQUANT_PACKS=paths, GERMQUANT_MODELS=tmp_path / "models")
    assert r.returncode == 0, r.stderr[-2000:]
    assert heavy == "", f"germquant.packs loaded {heavy}"


def test_python_m_germquant_packs_list_json_stays_free_of_torch(tmp_path):
    code = ("import runpy, sys\n"
            "sys.argv = ['germquant.packs', 'list', '--json', '--all']\n"
            "try:\n"
            "    runpy.run_module('germquant.packs', run_name='__main__', alter_sys=True)\n"
            "except SystemExit as e:\n"
            "    assert not e.code, e.code\n")
    heavy, r = _heavy_after(code, GERMQUANT_PACKS=FIXTURE_PACKS, GERMQUANT_MODELS=tmp_path / "models")
    assert r.returncode == 0, r.stderr[-2000:]
    assert heavy == "", f"python -m germquant.packs loaded {heavy}"
    listing = json.loads(r.stdout)
    assert listing["builtin_root_found"] is True
    assert {"celegans_germline", "test_pack"} <= {p["id"] for p in listing["packs"]}


def test_derive_and_the_project_listing_stay_free_of_torch(tmp_path):
    """Gate G6 for `pack derive` and `--project` (Train B P5-DERIVE): in one fresh interpreter, `python -m
    germquant.packs derive` into a project's packs/ folder, `list --json --all --project DIR` and `--version
    --json` (the features) load neither torch, CuPy, SpotMAX or Cellpose nor bioio or dask.array, and the listing
    shows the derived pack and the parent's copy as the project's."""
    proj = tmp_path / "proj"
    mdir = proj / "training" / "r01" / "models" / "m1"
    mdir.mkdir(parents=True)
    model = mdir / "m1"
    model.write_bytes(b"a stand-in model, not weights\n")
    card = mdir / "m1.card.json"
    card.write_text(json.dumps({"name": "m1", "model_file": str(model), "base_model": {"name": "none"}}),
                    encoding="ascii")
    code = ("import json, runpy, sys\n"
            "def run(argv):\n"
            "    sys.argv = ['germquant.packs', *argv]\n"
            "    try:\n"
            "        runpy.run_module('germquant.packs', run_name='__main__', alter_sys=True)\n"
            "    except SystemExit as e:\n"
            "        assert not e.code, (argv, e.code)\n"
            f"run(['derive', '--from', 'test_pack@0.1.0', '--recipe', 'with_model', '--model', {str(model)!r},\n"
            f"     '--card', {str(card)!r}, '--out', {str(proj / 'packs')!r}, '--json'])\n"
            f"run(['list', '--json', '--all', '--project', {str(proj)!r}])\n"
            "run(['--version', '--json'])\n")
    heavy, r = _heavy_after(code, heavy=HEAVY + READERS, GERMQUANT_PACKS=FIXTURE_PACKS,
                            GERMQUANT_MODELS=tmp_path / "models")
    assert r.returncode == 0, r.stderr[-2000:]
    assert heavy == "", f"pack derive, the --project listing or the features loaded {heavy}"
    dec, text, docs = json.JSONDecoder(), r.stdout.strip(), []
    while text:
        doc, end = dec.raw_decode(text)
        docs.append(doc)
        text = text[end:].strip()
    derived, listing, version = docs
    assert derived["ok"] and derived["pack"].startswith("test_pack_derived_")
    sources = {(p["id"], p["source"]) for p in listing["packs"]}
    assert {(derived["id"], "project"), ("test_pack", "project")} <= sources
    assert version["features"][0] == "packs" and "pack_derive" in version["features"]


def test_shallow_check_with_an_image_stays_free_of_torch_and_cupy(tmp_path):
    """Gate G6 with an image: the expects check reads the TIFF header with tifffile, never through bioio (whose
    dask import loads CuPy before torch where CuPy is installed, the cuBLAS hazard of rule R5)."""
    code = ("import germquant.packs as P, importlib.util\n"
            f"st = P.pack_status('test_pack', image={str(EXAMPLE)!r})\n"
            "if importlib.util.find_spec('bioio') is not None:\n"
            "    assert st['usable'] and st['warnings'] == [], st\n")
    heavy, r = _heavy_after(code, heavy=HEAVY + READERS, GERMQUANT_PACKS=FIXTURE_PACKS,
                            GERMQUANT_MODELS=tmp_path / "models")
    assert r.returncode == 0, r.stderr[-2000:]
    assert heavy == "", f"a shallow pack_status with an image loaded {heavy}"


def test_python_m_germquant_packs_check_with_an_image_stays_free_of_torch_and_cupy(tmp_path):
    code = ("import runpy, sys\n"
            f"sys.argv = ['germquant.packs', 'check', 'test_pack', '--image', {str(EXAMPLE)!r}, '--json']\n"
            "try:\n"
            "    runpy.run_module('germquant.packs', run_name='__main__', alter_sys=True)\n"
            "except SystemExit as e:\n"
            "    print('EXIT', e.code, file=sys.stderr)\n")
    heavy, r = _heavy_after(code, heavy=HEAVY + READERS, GERMQUANT_PACKS=FIXTURE_PACKS,
                            GERMQUANT_MODELS=tmp_path / "models")
    assert r.returncode == 0 and re.search(r"EXIT (0|2)\b", r.stderr), r.stderr[-2000:]
    assert json.loads(r.stdout)["pack"] == "test_pack"
    assert heavy == "", f"python -m germquant.packs check --image loaded {heavy}"


def test_needs_post_stages_check_stays_free_of_torch(tmp_path):
    """Gate G6 for needs.post_stages (Train B P5b, acceptance 5): in one fresh interpreter, the shallow pack_status
    of the post-stage pack (plan D35: its recipe basic needs the fixture post-stage) greys the recipe while the
    fixture is not installed and passes it once the on-disk fixture distribution is on sys.path. It reads metadata
    and the declaration file only: the fixture's package is never imported, and neither torch, CuPy, SpotMAX,
    Cellpose, bioio nor dask.array is loaded."""
    from post_helpers import POST_PACK_REF, make_post_pack

    from germquant.post_stages.testing import FIXTURE_NAME, FIXTURE_PACKAGE, write_fixture_distribution

    make_post_pack(tmp_path / "packs")
    dist = write_fixture_distribution(tmp_path / "dist")
    code = ("import importlib, sys\n"
            "import germquant.packs as P\n"
            f"st = P.pack_status({POST_PACK_REF!r}, recipe='basic')\n"
            f"assert not st['usable'] and any({FIXTURE_NAME!r} in b['what'] for b in st['blockers']), st\n"
            f"sys.path.append({str(dist)!r})\n"
            "importlib.invalidate_caches()\n"
            f"st = P.pack_status({POST_PACK_REF!r}, recipe='basic')\n"
            "assert st['usable'], st['blockers']\n"
            f"assert {FIXTURE_PACKAGE!r} not in sys.modules\n")
    heavy, r = _heavy_after(code, heavy=HEAVY + READERS, GERMQUANT_PACKS=tmp_path / "packs",
                            GERMQUANT_MODELS=tmp_path / "models")
    assert r.returncode == 0, r.stderr[-2000:]
    assert heavy == "", f"pack status of needs.post_stages loaded {heavy}"


def test_deep_check_imports_torch_before_any_backend(tmp_path):
    """Gate G6: a deep check imports germquant.pipeline (torch) before it loads any backend, so torch claims
    cuBLAS before CuPy, SpotMAX or Cellpose can. Each backend load is wrapped to record what was imported by
    then, and sys.modules keeps the order in which imports started."""
    if importlib.util.find_spec("torch") is None:
        pytest.skip("torch is not installed here")
    code = ("import json, sys\n"
            "import germquant.packs as P\n"
            "from germquant.segment import backends as SB\n"
            "from germquant.spots import backends as SPB\n"
            "assert 'torch' not in sys.modules and 'germquant.pipeline' not in sys.modules\n"
            "loads = []\n"
            "for mod in (SB, SPB):\n"
            "    def wrap(name, _orig=mod.get_backend):\n"
            "        loads.append([name, 'germquant.pipeline' in sys.modules, 'torch' in sys.modules])\n"
            "        return _orig(name)\n"
            "    mod.get_backend = wrap\n"
            "st = P.pack_status('test_pack', deep=True)\n"
            "order = list(sys.modules)\n"
            "print(json.dumps({'device': st['device'], 'usable': st['usable'], 'loads': loads,\n"
            "                  'pos': {m: order.index(m) for m in ('torch', 'cupy', 'spotmax', 'cellpose', 'dask')\n"
            "                          if m in sys.modules}}))\n")
    _, r = _heavy_after(code, GERMQUANT_PACKS=FIXTURE_PACKS, GERMQUANT_MODELS=tmp_path / "models")
    assert r.returncode == 0, r.stderr[-2000:]
    out = json.loads(r.stdout.strip().splitlines()[-1])
    assert out["usable"] and out["device"] in ("cpu", "cuda", "mps")
    assert [name for name, _, _ in out["loads"]] == ["classical", "log"]          # the basic recipe's backends
    assert all(pipe and torch for _, pipe, torch in out["loads"]), out["loads"]    # the pipeline came first
    pos = out["pos"]
    assert "torch" in pos and all(pos["torch"] < p for m, p in pos.items() if m != "torch"), pos


def test_malformed_pack_does_not_break_a_generic_run(pack_env, monkeypatch, tmp_path):  # noqa: F811
    root = tmp_path / "env_packs" / "broken" / "0.1.0"
    root.mkdir(parents=True)
    (root / "pack.yaml").write_text("not: [valid\n", encoding="ascii")
    monkeypatch.setenv("GERMQUANT_PACKS", str(tmp_path / "env_packs"))
    packs.refresh()
    assert any(r["state"] == "unusable" for r in packs.list_packs(all=True))
    cfg = load_config(ROOT / "config" / "profiles" / "generic.yaml")
    assert pipeline.resolve_model_path(cfg) == "cpsam"


def test_python_m_germquant_help_imports_torch_first():
    pytest.importorskip("torch")
    r = subprocess.run([PYTHON, "-X", "importtime", "-m", "germquant", "--help"], capture_output=True, text=True,
                       timeout=300, env=subprocess_env(), cwd=ROOT)
    assert r.returncode == 0 and "usage: germquant" in r.stdout, r.stderr[-2000:]
    names = [ln.rsplit("|", 1)[-1].strip() for ln in r.stderr.splitlines() if ln.startswith("import time:")]
    first_torch = names.index("torch")
    ours = [i for i, n in enumerate(names) if n.startswith("germquant.") and n not in ("germquant.__main__",)]
    assert ours, "no germquant submodule was imported"
    # importtime lists a module when its import finishes: torch (cli.py's first import) completes before any
    # other germquant module, cli itself included
    assert first_torch < min(ours), [names[i] for i in ours[:5]]


def test_version_entries(capsys):
    from germquant import __version__

    with pytest.raises(SystemExit) as e:
        cli.main(["--version"])
    assert e.value.code == 0
    line = capsys.readouterr().out.strip()
    assert line.startswith(f"germquant {__version__} git ") and " python " in line
    assert __version__ == "0.2.0.dev0"
    r = subprocess.run([PYTHON, "-m", "germquant.packs", "--version", "--json"], capture_output=True, text=True,
                       timeout=300, env=subprocess_env(), cwd=ROOT)
    info = json.loads(r.stdout)
    assert info["germquant"] == __version__ and info["git_sha"] and info["python"] and info["executable"]


def _non_editable_layout(tmp_path):
    """A copy of the package in a fake site-packages inside another project that has its own pyproject.toml
    and decoy config/ folders (one in the project, one in the working directory)."""
    other = tmp_path / "otherproject"
    site = other / ".venv" / "Lib" / "site-packages"
    shutil.copytree(ROOT / "src" / "germquant", site / "germquant",
                    ignore=shutil.ignore_patterns("__pycache__"))
    (other / "pyproject.toml").write_text("[project]\nname = 'napari-like'\n", encoding="ascii")
    for decoy in (other / "config", tmp_path / "cwd" / "config"):
        (decoy / "profiles").mkdir(parents=True)
        (decoy / "generic.yaml").write_text("io: {channel_map: nowhere.yaml}\n", encoding="ascii")
        (decoy / "profiles" / "generic.yaml").write_text("io: {channel_map: nowhere.yaml}\n", encoding="ascii")
    study = tmp_path / "study"
    study.mkdir()
    (study / "x.yaml").write_text("extends: generic.yaml\n", encoding="ascii")
    return site, study


def test_non_editable_install_names_the_missing_root_and_never_binds_elsewhere(tmp_path):
    site, study = _non_editable_layout(tmp_path)
    env = subprocess_env()
    env["PYTHONPATH"] = str(site)                     # the copied package, not this checkout's src
    env.pop("GERMQUANT_CONFIG_ROOT", None)
    env["GERMQUANT_PACKS"] = str(tmp_path / "no_packs_here")
    cwd = tmp_path / "cwd"
    r = subprocess.run([PYTHON, "-m", "germquant.packs", "list"], capture_output=True, text=True, timeout=300,
                       env=env, cwd=cwd)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "No built-in packs" in r.stdout and "editable checkout" in r.stdout and "GERMQUANT_CONFIG_ROOT" in r.stdout
    code = ("import germquant, sys\n"
            f"assert germquant.__file__.startswith({str(site)!r}), germquant.__file__\n"
            "from germquant.config import load_config, builtin_roots\n"
            "assert builtin_roots() == [], builtin_roots()\n"
            "try:\n"
            f"    load_config({str(study / 'x.yaml')!r})\n"
            "except FileNotFoundError as e:\n"
            "    print('RAISED', e)\n"
            "else:\n"
            "    print('BOUND')\n")
    r = subprocess.run([PYTHON, "-c", code], capture_output=True, text=True, timeout=300, env=env, cwd=cwd)
    assert r.returncode == 0, r.stderr[-2000:]
    assert r.stdout.startswith("RAISED") and "generic.yaml" in r.stdout
    # with GERMQUANT_CONFIG_ROOT pointing at a checkout (the container's fix), the same file resolves there
    env["GERMQUANT_CONFIG_ROOT"] = str(ROOT)
    r = subprocess.run([PYTHON, "-c", code.replace("assert builtin_roots() == [], builtin_roots()\n", "")],
                       capture_output=True, text=True, timeout=300, env=env, cwd=cwd)
    assert r.stdout.startswith("BOUND"), r.stdout + r.stderr[-1500:]
