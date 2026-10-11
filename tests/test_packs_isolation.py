"""Import guards and entry points (design 3.8, 3.9, 3.13 and gate G6; P1 acceptance 8, 10 and 12).

* In a fresh interpreter, germquant.packs.list_packs, a shallow germquant.packs.pack_status and
  ``python -m germquant.packs list --json`` leave torch, cupy, spotmax and cellpose out of sys.modules; a
  malformed pack under $GERMQUANT_PACKS is listed as unusable, and a generic-profile run still resolves cpsam.
* ``python -m germquant --help`` works and imports torch before anything else of germquant.
* ``germquant --version`` and ``python -m germquant.packs --version`` print the engine version, git sha and
  python.
* With the engine installed non-editable inside another project (a copy of the package in a fake
  site-packages under a folder with its own pyproject.toml and config/), ``pack list`` names the missing
  built-in root and both fixes, and ``extends: generic.yaml`` from outside the checkout raises instead of
  binding to that project's or the working directory's config/ folder.
"""
from __future__ import annotations

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


def _heavy_after(code: str, **env) -> tuple[str, subprocess.CompletedProcess]:
    full = _PROBE.format(heavy=HEAVY) + code
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
