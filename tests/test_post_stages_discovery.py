"""Post-stage discovery, declarations, post_stages.yaml and the torch-free commands (design 5.6; plan 4.5, phase
CORE; decisions D31, D34 and D39; plugin needs I12 to I14).

* Discovery reads the entry-point metadata and the declaration file germquant_post_stage.yaml only: the fixture's
  module is never imported, and listing loads no torch (fresh interpreters, gate G6).
* A declaration is validated; a broken one lists its post-stage as unusable, with the reason.
* In-process registration (tests, notebooks) behaves like an installed post-stage.
* The on-disk fixture distribution (germquant.post_stages.testing.write_fixture_distribution) is found by a child
  interpreter through PYTHONPATH.
* ``python -m germquant.post_stages write DIR --entries JSON`` validates against the post-stages installed in its
  own environment and writes DIR/post_stages.yaml as ASCII with LF; a refusal exits 2 and writes nothing.
* post_stages.yaml round-trips, its digest is taken on LF-normalised bytes, and a list of entries digests like the
  file write_post_stages_file makes of it.
"""
from __future__ import annotations

import copy
import hashlib
import importlib
import json

import pytest
from post_helpers import HEAVY, READERS, fixture_post_stage, heavy_after, run_child  # noqa: F401

from germquant import post_stages as PS
from germquant.post_stages import declaration as D
from germquant.post_stages import testing as T

# the parts of a row that come from the declaration file (dist and source depend on how it was installed)
_DECL_KEYS = ("name", "version", "title", "reads", "writes", "params")


def _row(rows, name):
    hits = [r for r in rows if r["name"] == name]
    assert len(hits) == 1, [r["name"] for r in rows]
    return hits[0]


def test_discovery_reads_only_metadata_and_the_declaration(tmp_path):
    """A child interpreter with the fixture distribution on PYTHONPATH lists it with its declaration, and neither
    the fixture's module nor torch, CuPy, SpotMAX, Cellpose, bioio, dask.array, numpy or pandas is loaded."""
    site = T.write_fixture_distribution(tmp_path / "site")
    code = ("import json, sys\n"
            "import germquant.post_stages as PS\n"
            "rows = PS.list_post_stages()\n"
            "print(json.dumps(rows))\n"
            f"assert {T.FIXTURE_PACKAGE!r} not in sys.modules, 'discovery imported the post-stage module'\n")
    heavy, r = heavy_after(code, heavy=HEAVY + READERS + ("numpy", "pandas"), pythonpath=[site])
    assert r.returncode == 0, r.stderr[-3000:]
    assert heavy == "", f"post-stage discovery loaded {heavy}"
    row = _row(json.loads(r.stdout), T.FIXTURE_NAME)
    want = D.normalise_declaration(T.FIXTURE_NAME, T.fixture_declaration())
    assert {k: row[k] for k in _DECL_KEYS} == {k: want[k] for k in _DECL_KEYS}
    assert row["dist"] == T.FIXTURE_DIST
    assert row["source"] == f"entry point {T.FIXTURE_PACKAGE}:measure"
    assert row["version"] == T.FIXTURE_VERSION
    assert row["reads"] == {"roles": ["foci"], "tables": ["nuclei", "spots"], "labels": True, "image": True}
    assert row["writes"] == {T.FIXTURE_TABLE: list(T.FIXTURE_COLUMNS)}
    params = {p["key"]: p for p in row["params"]}
    assert params["factor"]["type"] == "float" and params["factor"]["default"] == 1.0
    assert params["statistic"]["choices"] == ["mean", "max"]
    assert params["fail"] == {"key": "fail", "type": "bool", "default": False, "label": params["fail"]["label"]}


def _bad(edit):
    d = copy.deepcopy(T.fixture_declaration())
    edit(d)
    return D.declaration_problems(T.FIXTURE_NAME, d)


_P = "params"


@pytest.mark.parametrize("edit, fragment", [
    (lambda d: d.pop("version"), "version"),
    (lambda d: d.update(version=1.0), "version"),
    (lambda d: d.update(version="one"), "version"),
    (lambda d: d.update(colour="red"), "unknown key"),
    (lambda d: d.update(title=3), "title"),
    (lambda d: d.update(reads=[]), "reads"),
    (lambda d: d["reads"].update(tables=["post_other_counts"]), "post_"),
    (lambda d: d["reads"].update(tables=["Nuclei"]), "table name"),
    (lambda d: d["reads"].update(tables=["nuclei", "nuclei"]), "twice"),
    (lambda d: d["reads"].update(image=False), "reads.roles"),
    (lambda d: d["reads"].update(roles=[]), "reads.image"),
    (lambda d: d["reads"].update(roles=["foci", "bad role"]), "role"),
    (lambda d: d["reads"].update(labels="yes"), "reads.labels"),
    (lambda d: d["reads"].update(extra=True), "unknown key"),
    (lambda d: d.update(writes={}), "writes"),
    (lambda d: d.update(writes={"counts": []}), "column"),
    (lambda d: d.update(writes={"counts": ["a", "a"]}), "twice"),
    (lambda d: d.update(writes={"counts": ["bad column"]}), "column"),
    (lambda d: d.update(writes={"bad__name": ["a"]}), "table name"),
    (lambda d: d.update(writes={"trailing_": ["a"]}), "table name"),
    (lambda d: d.update(params={}), "params"),
    (lambda d: d[_P].append({"key": "factor", "type": "float", "default": 1.0}), "twice"),
    (lambda d: d[_P].append({"key": "Bad-Key", "type": "float", "default": 1.0}), "key"),
    (lambda d: d[_P].append({"key": "k", "type": "colour", "default": 1}), "type"),
    (lambda d: d[_P].append({"key": "k", "type": "float"}), "default"),
    (lambda d: d[_P].append({"key": "k", "type": "float", "min": 2, "max": 1, "default": 1.5}), "min"),
    (lambda d: d[_P].append({"key": "k", "type": "float", "min": 0, "max": 1, "default": 5}), "default"),
    (lambda d: d[_P].append({"key": "k", "type": "float", "min": True, "default": 1.0}), "min"),
    (lambda d: d[_P].append({"key": "k", "type": "float", "default": float("nan")}), "default"),
    (lambda d: d[_P].append({"key": "k", "type": "int", "default": 1.5}), "default"),
    (lambda d: d[_P].append({"key": "k", "type": "int", "min": 0.5, "default": 1}), "min"),
    (lambda d: d[_P].append({"key": "k", "type": "choice", "default": "a"}), "choices"),
    (lambda d: d[_P].append({"key": "k", "type": "choice", "choices": ["a", "a"], "default": "a"}), "choices"),
    (lambda d: d[_P].append({"key": "k", "type": "choice", "choices": ["a"], "default": "b"}), "default"),
    (lambda d: d[_P].append({"key": "k", "type": "choice", "choices": [1, 2], "default": 1.0}), "default"),
    (lambda d: d[_P].append({"key": "k", "type": "bool", "default": 0}), "default"),
    (lambda d: d[_P].append({"key": "k", "type": "bool", "min": 0, "default": True}), "min"),
    (lambda d: d[_P].append({"key": "k", "type": "str", "default": 3}), "default"),
    (lambda d: d[_P].append({"key": "k", "type": "float", "unit": 5, "default": 1.0}), "unit"),
    (lambda d: d[_P].append({"key": "k", "type": "float", "default": 1.0, "defualt": 2}), "unknown key"),
])
def test_the_declaration_is_validated(edit, fragment):
    probs = _bad(edit)
    assert probs, "a broken declaration passed"
    assert any(fragment in p for p in probs), probs
    for p in probs:                                   # problems are ASCII text that names the post-stage
        assert p.isascii() and T.FIXTURE_NAME in p, p


def test_a_sound_declaration_and_the_declaration_file():
    good = T.fixture_declaration()
    assert D.declaration_problems(T.FIXTURE_NAME, good) == []
    assert any("name" in p for p in D.declaration_problems("Bad-Name", good))
    assert any("name" in p for p in D.declaration_problems("double__underscore", good))
    assert any("name" in p for p in D.declaration_problems("x" * 41, good))
    doc = {"format": 1, "post_stages": {T.FIXTURE_NAME: good}}
    decl, probs = D.declaration_from_document(doc, T.FIXTURE_NAME)
    assert probs == [] and decl == good
    for bad, fragment in (([], "mapping"), ({"format": 2, "post_stages": {}}, "format"),
                          ({"format": 1}, "post_stages"), ({"format": 1, "post_stages": []}, "post_stages"),
                          ({"format": 1, "post_stages": {}, "extra": 1}, "unknown key"),
                          ({"format": 1, "post_stages": {"other": good}}, "declares no post-stage named")):
        decl, probs = D.declaration_from_document(bad, T.FIXTURE_NAME)
        assert decl is None and any(fragment in p for p in probs), (bad, probs)
    # a registration is checked by the same rules
    broken = copy.deepcopy(good)
    broken["writes"] = {}
    with pytest.raises(ValueError, match="writes"):
        PS.register_post_stage("broken_stage", lambda inputs, params: {}, broken)
    assert "broken_stage" not in {r["name"] for r in PS.list_post_stages()}


def _dist(site, name, entry_points, files):
    return T.write_distribution(site, name=name, version="0.1.0", entry_points=entry_points, files=files)


def test_broken_declarations_are_listed_as_unusable(tmp_path, monkeypatch):
    """Installed post-stages whose declaration cannot be used are listed (with include_unusable) with the reason,
    never offered as usable, and refused by validate_post_stages. Nothing is imported to find that out."""
    site = tmp_path / "site"
    good_text = T.declaration_text(T.FIXTURE_NAME)
    _dist(site, "lab-no-file", {"no_file": "lab_no_file:run"}, {"lab_no_file/__init__.py": "def run(i, p):\n    return {}\n"})
    _dist(site, "lab-bad-yaml", {"bad_yaml": "lab_bad_yaml:run"},
          {"lab_bad_yaml/__init__.py": "", f"lab_bad_yaml/{D.DECLARATION_FILE}": "format: 1\npost_stages: [\n"})
    _dist(site, "lab-other-name", {"other_name": "lab_other:run"},
          {"lab_other/__init__.py": "", f"lab_other/{D.DECLARATION_FILE}": good_text})
    _dist(site, "lab-module", {"module_only": "lab_module_only:run"}, {"lab_module_only.py": "def run(i, p):\n    return {}\n"})
    _dist(site, "lab-twice-a", {"twice": "lab_twice_a:run"},
          {"lab_twice_a/__init__.py": "", f"lab_twice_a/{D.DECLARATION_FILE}": T.declaration_text("twice")})
    _dist(site, "lab-twice-b", {"twice": "lab_twice_b:run"},
          {"lab_twice_b/__init__.py": "", f"lab_twice_b/{D.DECLARATION_FILE}": T.declaration_text("twice")})
    _dist(site, "lab-bad-name", {"Bad-Name": "lab_bad_name:run"},
          {"lab_bad_name/__init__.py": "", f"lab_bad_name/{D.DECLARATION_FILE}": good_text})
    _dist(site, "lab-not-utf8", {"not_utf8": "lab_not_utf8:run"},
          {"lab_not_utf8/__init__.py": "", f"lab_not_utf8/{D.DECLARATION_FILE}": b"format: 1\n# caf\xe9\n"})
    _dist(site, "lab-good", {"good_stage": "lab_good.sub:run"},
          {"lab_good/__init__.py": "", "lab_good/sub.py": "def run(i, p):\n    return {}\n",
           f"lab_good/{D.DECLARATION_FILE}": T.declaration_text("good_stage")})
    monkeypatch.syspath_prepend(str(site))
    importlib.invalidate_caches()
    rows = {r["name"]: r for r in PS.list_post_stages(include_unusable=True)}
    expect = {"no_file": "no germquant_post_stage.yaml", "bad_yaml": "cannot be read",
              "other_name": "declares no post-stage named", "module_only": "a module, not a package",
              "twice": "more than one distribution", "Bad-Name": "name", "not_utf8": "UTF-8"}
    for name, fragment in expect.items():
        assert rows[name]["usable"] is False, rows[name]
        assert any(fragment in p for p in rows[name]["problems"]), (name, rows[name]["problems"])
    assert rows["good_stage"]["usable"] is True and rows["good_stage"]["problems"] == []
    assert rows["good_stage"]["dist"] == "lab-good" and rows["good_stage"]["source"] == "entry point lab_good.sub:run"
    usable = {r["name"] for r in PS.list_post_stages()}
    assert "good_stage" in usable and not usable & set(expect)
    probs = PS.validate_post_stages([{"name": "no_file"}])
    assert probs and "not usable" in probs[0] and "no germquant_post_stage.yaml" in probs[0], probs
    assert PS.validate_post_stages([{"name": "good_stage"}]) == []


def test_in_process_registration():
    name = T.register_fixture()
    try:
        row = _row(PS.list_post_stages(), name)
        assert row["source"] == "registered" and row["dist"] is None
        assert row["version"] == T.FIXTURE_VERSION and row["writes"] == {T.FIXTURE_TABLE: list(T.FIXTURE_COLUMNS)}
        with pytest.raises(ValueError, match="already registered"):
            T.register_fixture()
        T.register_fixture(replace=True)
        with pytest.raises(TypeError):
            PS.register_post_stage("not_callable", 3, T.fixture_declaration())
        assert PS.validate_post_stages([{"name": name, "params": {"factor": 2.5, "statistic": "max"}}]) == []
        fn, decl = PS.load_post_stage(name)
        assert callable(fn) and decl["version"] == T.FIXTURE_VERSION
        for entries, fragment in (
                ([{"name": name, "params": {"factor": 99.0}}], "factor"),
                ([{"name": name, "params": {"factor": "2"}}], "factor"),
                ([{"name": name, "params": {"statistic": "median"}}], "statistic"),
                ([{"name": name, "params": {"min_voxels": 1.5}}], "min_voxels"),
                ([{"name": name, "params": {"radius": 1.0}}], "radius"),
                ([{"name": name}, {"name": name}], "more than once"),
                ([{"name": "nope"}], "not installed"),
                ([{"params": {}}], "name"),
                ([{"name": name, "params": [1]}], "params"),
                ([{"name": name, "colour": 1}], "unknown key"),
                ({"name": name}, "list")):
            probs = PS.validate_post_stages(entries)
            assert probs and any(fragment in p for p in probs), (entries, probs)
    finally:
        PS.unregister_post_stage(name)
    assert name not in {r["name"] for r in PS.list_post_stages()}
    with pytest.raises(LookupError):
        PS.unregister_post_stage(name)
    with pytest.raises(LookupError, match="not installed"):
        PS.load_post_stage(name)


def test_the_on_disk_fixture_is_found_in_a_subprocess(tmp_path):
    """I14: the fixture written as a dist-info plus a package is discovered from metadata by any interpreter that
    has the folder on PYTHONPATH, here through the command line the plugin's runner calls."""
    site = T.write_fixture_distribution(tmp_path / "site")
    assert (site / f"{T.FIXTURE_PACKAGE}-{T.FIXTURE_VERSION}.dist-info" / "entry_points.txt").is_file()
    assert (site / T.FIXTURE_PACKAGE / D.DECLARATION_FILE).is_file()
    r = run_child(["-m", "germquant.post_stages", "list", "--json"], pythonpath=[site])
    assert r.returncode == 0, r.stderr[-3000:]
    row = _row(json.loads(r.stdout), T.FIXTURE_NAME)
    assert row["dist"] == T.FIXTURE_DIST and row["version"] == T.FIXTURE_VERSION
    assert row["params"] == D.normalise_declaration(T.FIXTURE_NAME, T.fixture_declaration())["params"]
    # the suite's own interpreter does not have it: the child found it through PYTHONPATH alone
    assert T.FIXTURE_NAME not in {r["name"] for r in PS.list_post_stages()}


def test_list_json_loads_no_torch(tmp_path):
    """I12 and gate G6: ``python -m germquant.post_stages list --json`` stays free of torch, CuPy, SpotMAX,
    Cellpose, bioio and dask.array."""
    site = T.write_fixture_distribution(tmp_path / "site")
    code = ("import runpy, sys\n"
            "sys.argv = ['germquant.post_stages', 'list', '--json']\n"
            "try:\n"
            "    runpy.run_module('germquant.post_stages', run_name='__main__', alter_sys=True)\n"
            "except SystemExit as e:\n"
            "    assert not e.code, e.code\n")
    heavy, r = heavy_after(code, heavy=HEAVY + READERS, pythonpath=[site])
    assert r.returncode == 0, r.stderr[-3000:]
    assert heavy == "", f"python -m germquant.post_stages list loaded {heavy}"
    assert T.FIXTURE_NAME in {row["name"] for row in json.loads(r.stdout)}


def test_write_validates_in_its_own_environment(tmp_path):
    """D39: the plugin writes post_stages.yaml through the engine environment, because only that environment has
    the post-stages installed. The fixture sits on the child's PYTHONPATH only: the suite's own validation refuses
    the entries, the child accepts them (and loads no torch), and an unknown name or an out-of-range parameter
    exits 2 without touching the file."""
    site = T.write_fixture_distribution(tmp_path / "site")
    proj = tmp_path / "proj"
    proj.mkdir()
    entries = [{"name": T.FIXTURE_NAME, "params": {"statistic": "max", "factor": 2}}]
    assert any("not installed" in p for p in PS.validate_post_stages(entries))
    code = ("import runpy, sys\n"
            f"sys.argv = ['germquant.post_stages', 'write', {str(proj)!r}, '--entries', {json.dumps(entries)!r},"
            " '--json']\n"
            "try:\n"
            "    runpy.run_module('germquant.post_stages', run_name='__main__', alter_sys=True)\n"
            "except SystemExit as e:\n"
            "    assert not e.code, e.code\n")
    heavy, r = heavy_after(code, heavy=HEAVY + READERS, pythonpath=[site])
    assert r.returncode == 0, r.stderr[-3000:]
    assert heavy == "", f"python -m germquant.post_stages write loaded {heavy}"
    out = json.loads(r.stdout)
    target = proj / PS.POST_STAGES_FILE
    data = target.read_bytes()
    assert out["written"] is True and out["path"] == str(target)
    assert data.isascii() and b"\r" not in data
    assert out["sha256"] == hashlib.sha256(data).hexdigest()
    # the file holds the entries as given, with the float parameter as a float and the keys sorted
    assert PS.load_post_stages_file(target) == [{"name": T.FIXTURE_NAME, "params": {"factor": 2.0,
                                                                                    "statistic": "max"}}]
    assert out["entries"] == PS.load_post_stages_file(target)
    for bad, fragment in (([{"name": "nope"}], "not installed"),
                          ([{"name": T.FIXTURE_NAME, "params": {"factor": 99}}], "factor"),
                          ("not json", "JSON")):
        text = bad if isinstance(bad, str) else json.dumps(bad)
        r = run_child(["-m", "germquant.post_stages", "write", proj, "--entries", text, "--json"],
                      pythonpath=[site])
        assert r.returncode == 2, (bad, r.returncode, r.stderr[-2000:])
        assert fragment in r.stderr, r.stderr
        refusal = json.loads(r.stdout)
        assert refusal["written"] is False and any(fragment in p for p in refusal["problems"])
        assert target.read_bytes() == data, "a refused write touched the file"
    r = run_child(["-m", "germquant.post_stages", "write", tmp_path / "missing", "--entries", "[]"],
                  pythonpath=[site])
    assert r.returncode == 2 and "not a folder" in r.stderr, r.stderr


@pytest.mark.usefixtures("fixture_post_stage")
def test_post_stages_file_round_trip(tmp_path):
    """post_stages.yaml: written ASCII with LF (non-ASCII text escaped), read back unchanged, refused when
    malformed; its digest is taken on LF-normalised bytes, and a list digests like the file written from it."""
    name = T.FIXTURE_NAME
    path = tmp_path / PS.POST_STAGES_FILE
    entries = [{"name": name, "params": {"statistic": "max", "factor": 3}}]
    PS.write_post_stages_file(path, entries)
    data = path.read_bytes()
    assert data.isascii() and b"\r" not in data
    loaded = PS.load_post_stages_file(path)
    assert loaded == [{"name": name, "params": {"factor": 3.0, "statistic": "max"}}]
    assert PS.post_stages_sha256(path) == hashlib.sha256(data).hexdigest() == PS.post_stages_sha256(entries)
    path.write_bytes(data.replace(b"\n", b"\r\n"))
    assert PS.post_stages_sha256(path) == hashlib.sha256(data).hexdigest(), "the digest is not LF-normalised"
    assert PS.load_post_stages_file(path) == loaded
    # a text parameter with non-ASCII characters is written with YAML escapes and read back as it was
    word = "caf" + chr(0xE9) + " " + chr(0x3B3) + "H2AX"
    decls = {"text_stage": {"version": "1.0.0", "reads": {"tables": ["nuclei"]}, "writes": {"t": ["a"]},
                            "params": [{"key": "note", "type": "str", "default": ""}]}}
    p2 = tmp_path / "text" / PS.POST_STAGES_FILE
    p2.parent.mkdir()
    PS.write_post_stages_file(p2, [{"name": "text_stage", "params": {"note": word}}], declarations=decls)
    assert p2.read_bytes().isascii()
    assert PS.load_post_stages_file(p2) == [{"name": "text_stage", "params": {"note": word}}]
    with pytest.raises(ValueError, match="not installed"):
        PS.write_post_stages_file(p2, [{"name": "nope"}])
    assert PS.load_post_stages_file(p2)[0]["params"]["note"] == word, "a refused write changed the file"
    # an empty file, and an explicit empty list, switch every post-stage off
    (tmp_path / "empty.yaml").write_text("", encoding="ascii")
    assert PS.load_post_stages_file(tmp_path / "empty.yaml") == []
    PS.write_post_stages_file(tmp_path / "none.yaml", [])
    assert PS.load_post_stages_file(tmp_path / "none.yaml") == []
    for text, fragment in (("name: x\n", "list"), ("- {params: {}}\n", "name"), ("- {name: x, params: [1]}\n",
                                                                                 "params"),
                           ("- {name: x, colour: 1}\n", "unknown key"), ("- [\n", "cannot be read"),
                           ("- name: 3\n", "name")):
        bad = tmp_path / "bad.yaml"
        bad.write_text(text, encoding="ascii")
        with pytest.raises(ValueError, match=fragment):
            PS.load_post_stages_file(bad)
