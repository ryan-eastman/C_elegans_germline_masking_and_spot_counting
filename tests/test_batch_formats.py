"""Batch over outside formats (docs/DESIGN_system_packs.md 7.8; PLAN.md 4.4 "P5-SHEET" and decision D20; P5
acceptance 6).

* A default batch discovers exactly what 8eb73de discovers: the .nd2 files of a folder that also holds .ome.tif
  files, a folder named like an image and germquant outputs, under the ids Path.stem gives (io.formats.image_stem
  equals Path.stem for an .nd2). The function of 8eb73de is kept below, verbatim, as the reference.
* With io.input_suffixes the batch discovers those suffixes, matched on the whole format suffix (x.ome.tif is
  .ome.tif, never .tif), and keys out_dir, resume and the collect ids on image_stem, so x.ome.tif is the image x.
* germquant's own output images are never discovered, and with io.input_suffixes set nothing under an output root
  inside the input folder is either (D20: a default batch keeps the base's set, an .nd2 under the output root
  included); an output root at or above the input folder keeps every input (review S3).
* io.input_suffixes is checked strictly: a list of the suffixes germquant reads (recipe_keys), in batch before any
  file is discovered and in process_image before any stage.
* Where a snakemake executable is available (D25: a throwaway environment, never the .venv), a dry run of the
  Snakefile prints the default mode's .nd2 commands, the .ome.tif jobs of a recipe with io.input_suffixes under
  their image_stem ids, and the pack mode's `germquant run --pack` commands.
"""
from __future__ import annotations

import fnmatch
import json
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest
from packs_helpers import TEST_PACK, subprocess_env
from test_end_to_end import patched_reader  # noqa: F401  (fixture: the 3-channel synthetic gonad)

from germquant import batch, cli, pipeline, provenance, recipe_keys
from germquant.batch import Exclusions
from germquant.config import load_config
from germquant.io.formats import SUPPORTED_SUFFIXES, image_stem
from germquant.project.images import OUTPUT_SUFFIXES

ROOT = Path(__file__).resolve().parents[1]


def base_discover_files(root, glob="**/*.nd2", exclude_patterns=None, exclusions=None):
    """batch.discover_files as it is at 8eb73de (src/germquant/batch.py lines 83 to 103), verbatim but for the name:
    the reference a default batch must keep."""
    root = Path(root)
    excl = exclusions if isinstance(exclusions, Exclusions) else Exclusions([str(x) for x in (exclusions or [])])
    keep: list[Path] = []
    dropped: list[tuple[Path, str]] = []
    for f in sorted(root.glob(glob)):
        if f.suffix.lower() != ".nd2":
            continue
        name = f.name.lower()
        pat = next((p for p in (exclude_patterns or []) if fnmatch.fnmatch(name, p.lower())), None)
        if pat is not None:
            dropped.append((f, f"exclude_pattern {pat}")); continue
        hit = excl.matches(f.stem)
        if hit is not None:
            dropped.append((f, f"exclusions_file {hit}")); continue
        keep.append(f)
    return keep, dropped


# a study folder: .nd2 stand-ins (upper case, a space, a sub folder, an excluded gonad, a 10x overview, an .nd2
# inside the output root), outside formats, germquant's own outputs, a folder named like an image, an OME-Zarr
# folder and other files
FILES = ["a_HS_male_07.nd2", "B_noHS_herm_01.ND2", "c_10x_overview.nd2", "d.nd2", "e f.nd2",
         "sub/g.nd2", "sub dir/h.nd2", "i.ome.tif", "j.tif", "k.TIFF", "l.ome.tiff", "sub/m.ome.tif",
         "a_HS_male_07__nuclei_labels.tif", "d__spots.tif", "d__granules_labels.tif", "d__montage.png",
         "n.nd2.bak", "notes.txt", "o.czi", "results/r.nd2", "results/d/d__nuclei_labels.tif", "results/q.ome.tif",
         "p.ome.zarr/.zattrs", "p.ome.zarr/0/0.0.0"]
FOLDERS = ["folder.nd2", "w.tif"]                      # folders that only look like images
# names whose suffix rules are subtle, for discovery alone (a batch would name a results folder "x." after one)
TRICKY = ["x..nd2", ".nd2", "z.Nd2", "y.nd2.nd2"]
# what 8eb73de discovers with the glob **/*.nd2 and no exclusion; pathlib's glob ignores letter case on Windows only
BASE_ND2 = ["a_HS_male_07.nd2", "c_10x_overview.nd2", "d.nd2", "e f.nd2", "folder.nd2", "results/r.nd2",
            "sub dir/h.nd2", "sub/g.nd2"] + (["B_noHS_herm_01.ND2"] if os.name == "nt" else [])


def _tree(tmp_path: Path, tricky: bool = False) -> Path:
    root = tmp_path / "in"
    for name in FILES + (TRICKY if tricky else []):
        p = root / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"")
    for name in FOLDERS:
        (root / name).mkdir(parents=True)
    return root


def _default_set(root: Path) -> list[Path]:
    """What a `germquant batch` of 8eb73de processes under config.yaml's io keys: the base function's files."""
    cfg = load_config(ROOT / "config" / "config.yaml")
    return base_discover_files(root, cfg.get("io.input_glob"), cfg.get("io.exclude_patterns"))[0]


def _rel(root: Path, files) -> list[str]:
    return [f.relative_to(root).as_posix() for f in files]


def _study(tmp_path: Path, **io_keys) -> Path:
    """A study profile outside the repository: config.yaml, the classical segmenter, and the given io keys."""
    lines = [f"extends: {(ROOT / 'config' / 'config.yaml').as_posix()}", "segmentation:", "  nuclei:",
             "    method: classical"]
    if io_keys:
        lines.append("io:")
        lines += [f"  {k}: {json.dumps(v)}" for k, v in io_keys.items()]
    study = tmp_path / "study.yaml"
    study.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return study


def _fake_process(calls: list, tables: bool = True):
    """process_image with today's exact keywords (D40): a call with any other keyword fails. It records the file and
    the out_dir, and writes the image summary and a completion marker under the image id process_image gives (the
    file name without its format suffix)."""
    def fake(nd2, cfg, out_dir, *, xy_stride=1, z_range=None, prov=None):
        iid = image_stem(nd2)
        calls.append((Path(nd2).name, Path(out_dir)))
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        if tables:
            pd.DataFrame([{"image_id": iid, "n_nuclei": 3, "n_germline_nuclei": 2, "qc_pass": True, "qc_flags": ""}]) \
                .to_csv(out_dir / f"{iid}__image_summary.csv", index=False)
        provenance.write_done_marker(out_dir, iid, config_hash=cfg.hash, enabled_stages=batch.enabled_stage_names(cfg),
                                     stage_hashes={"segment": _seg_hash(cfg)}, xy_stride=xy_stride)
        return {"image_id": iid, "n_nuclei": 3, "qc_pass": True, "qc_flags": [], "out_dir": str(out_dir),
                "tables": {"image_summary": pd.DataFrame([{"n_germline_nuclei": 2}])}, "stages": {}}
    return fake


def _seg_hash(cfg):
    from germquant.stages import stage_hashes

    return stage_hashes(cfg, {}, provenance.tool_versions(), "unused")["segment"]


def _spy_discover(monkeypatch) -> list[dict]:
    seen: list[dict] = []
    real = batch.discover_files

    def spy(*a, **k):
        seen.append(k)
        return real(*a, **k)

    monkeypatch.setattr(batch, "discover_files", spy)
    return seen


# ---------------------------------------------------------------------------------------------------
# acceptance 6: the default set, the .ome.tif ids, the outputs
# ---------------------------------------------------------------------------------------------------

def test_default_suffixes_match_base_nd2_set(tmp_path, monkeypatch, capsys):
    tricky = _tree(tmp_path / "tricky", tricky=True)
    for glob in ("**/*.nd2", "**/*", "*.nd2", "**/*.ND2", "*"):
        for excludes, excl in ((None, None), (["*10x*", "*largeimage*"], ["HS_male_07"])):
            assert batch.discover_files(tricky, glob, excludes, excl) == base_discover_files(tricky, glob, excludes,
                                                                                             excl), (glob, excludes)
    root = _tree(tmp_path)
    assert sorted(_rel(root, batch.discover_files(root)[0])) == sorted(BASE_ND2)
    files, dropped = batch.discover_files(root, "**/*.nd2", ["*10x*"], ["HS_male_07"])
    assert files == base_discover_files(root, "**/*.nd2", ["*10x*"], ["HS_male_07"])[0]
    assert {f.name: why.split()[0] for f, why in dropped} == {"a_HS_male_07.nd2": "exclusions_file",
                                                               "c_10x_overview.nd2": "exclude_pattern"}
    # through `germquant batch` with a recipe that sets no io.input_suffixes: the same files, in the same order, under
    # the same ids and results folders (Path.stem), and discover_files is called with today's arguments only (D20)
    calls: list = []
    monkeypatch.setattr(pipeline, "process_image", _fake_process(calls))
    seen = _spy_discover(monkeypatch)
    out = tmp_path / "out"
    assert cli.main(["batch", str(root), "--config", str(_study(tmp_path)), "--out", str(out)]) == 0
    assert seen == [{}]
    want = _default_set(root)
    assert len(want) == len(BASE_ND2) - 1                   # the 10x overview is excluded by config.yaml
    assert [c[0] for c in calls] == [f.name for f in want]
    assert [c[1] for c in calls] == [out / f.relative_to(root).parent / f.stem for f in want]
    assert all(image_stem(f) == f.stem for f in want)
    summ = pd.read_csv(out / "batch_summary.csv", dtype=str)
    assert list(summ["image_id"]) == [f.stem for f in want]
    # --resume finds every marker under the same ids and collects the same study
    calls.clear()
    capsys.readouterr()
    assert cli.main(["batch", str(root), "--config", str(_study(tmp_path)), "--out", str(out), "--resume"]) == 0
    assert calls == [] and capsys.readouterr().out.count("already done") == len(want)
    assert sorted(pd.read_csv(out / "batch_summary.csv", dtype=str)["image_id"]) == sorted(f.stem for f in want)


def test_ome_tif_ids_drop_double_suffix(tmp_path, monkeypatch, capsys):
    root = _tree(tmp_path)
    assert image_stem("x.ome.tif") == "x" and image_stem("x.ome.tiff") == "x" and image_stem("x.nd2") == "x"
    files, _ = batch.discover_files(root, "**/*", suffixes=(".ome.tif",), skip_root=root / "results")
    assert _rel(root, files) == ["i.ome.tif", "sub/m.ome.tif"]
    # a suffix is the whole format suffix: .tif never takes an OME-TIFF, and a list may name several
    files, _ = batch.discover_files(root, "**/*", suffixes=(".tif", ".tiff"))
    assert sorted(_rel(root, files)) == ["j.tif", "k.TIFF", "w.tif"]
    files, _ = batch.discover_files(root, "**/*", suffixes=(".ome.zarr", ".czi", ".nd2"))
    assert sorted(_rel(root, files)) == sorted(["B_noHS_herm_01.ND2", "a_HS_male_07.nd2", "c_10x_overview.nd2",
                                                "d.nd2", "e f.nd2", "folder.nd2", "o.czi", "p.ome.zarr",
                                                "results/r.nd2", "sub dir/h.nd2", "sub/g.nd2"])
    # the exclusions match the image id, so an OME-TIFF is excluded by its id like an .nd2
    files, dropped = batch.discover_files(root, "**/*", None, ["m"], suffixes=(".ome.tif",))
    assert _rel(root, files) == ["i.ome.tif", "results/q.ome.tif"] and [d[1] for d in dropped] == ["exclusions_file m"]
    # `germquant batch` with io.input_suffixes [".ome.tif"]: out_dir, the marker, resume and collect use the id x
    calls: list = []
    monkeypatch.setattr(pipeline, "process_image", _fake_process(calls))
    seen = _spy_discover(monkeypatch)
    study = _study(tmp_path, input_glob="**/*", input_suffixes=[".ome.tif"])
    out = tmp_path / "out"
    assert cli.main(["batch", str(root), "--config", str(study), "--out", str(out)]) == 0
    assert seen == [{"suffixes": (".ome.tif",), "skip_root": out}]
    assert calls == [("i.ome.tif", out / "i"), ("q.ome.tif", out / "results" / "q"), ("m.ome.tif", out / "sub" / "m")]
    assert (out / "i" / "i__done.json").is_file() and not (out / "i.ome").exists()
    assert list(pd.read_csv(out / "batch_summary.csv")["image_id"]) == ["i", "q", "m"]
    decisions: list = []
    real_done = batch.is_done
    monkeypatch.setattr(batch, "is_done", lambda *a, **k: decisions.append(a[:2]) or real_done(*a, **k))
    calls.clear()
    capsys.readouterr()
    assert cli.main(["batch", str(root), "--config", str(study), "--out", str(out), "--resume"]) == 0
    assert calls == [] and capsys.readouterr().out.count("already done") == 3
    assert decisions == [(out / "i", "i"), (out / "results" / "q", "q"), (out / "sub" / "m", "m")]
    assert sorted(pd.read_csv(out / "batch_summary.csv")["image_id"]) == ["i", "m", "q"]
    # an image that fails is named by its id too
    monkeypatch.setattr(pipeline, "process_image", lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")))
    assert cli.main(["batch", str(root), "--config", str(study), "--out", str(tmp_path / "out2")]) == 0
    summ = pd.read_csv(tmp_path / "out2" / "batch_summary.csv")
    assert list(summ["image_id"]) == ["i", "q", "m"] and summ["qc_flags"].str.startswith("EXCEPTION").all()


def test_outputs_never_discovered(tmp_path, monkeypatch):
    root = _tree(tmp_path)
    every = tuple(SUPPORTED_SUFFIXES)
    files, _ = batch.discover_files(root, "**/*", suffixes=every, skip_root=root / "results")
    names = _rel(root, files)
    assert not [n for n in names if "__" in n], names                       # no label, spot or mask image
    assert not [n for n in names if n.startswith("results/")], names        # nothing under the output root
    assert "j.tif" in names and "i.ome.tif" in names and "d.nd2" in names
    # every output image name germquant writes, .tif and .tiff, at any depth, in any letter case
    out_names = [f"z{s}" for s in OUTPUT_SUFFIXES] + [f"Z{s.upper()}F" for s in OUTPUT_SUFFIXES]
    for n in out_names:
        (root / "sub" / n).write_bytes(b"")
    files, _ = batch.discover_files(root, "**/*", suffixes=every)
    assert not [f.name for f in files if f.name in out_names]
    # the default set holds no output (no germquant output ends in .nd2), and keeps the base's .nd2 under the output
    # root (D20: the output root is skipped only for a recipe that sets io.input_suffixes)
    files, _ = batch.discover_files(root)
    assert "results/r.nd2" in _rel(root, files) and not [f for f in files if "__" in f.name]
    # `germquant batch` with the output root inside the input folder: its label images and its own OME-TIFF are
    # never taken for inputs, and a second batch into the same root finds the same images
    calls: list = []
    monkeypatch.setattr(pipeline, "process_image", _fake_process(calls))
    study = _study(tmp_path, input_glob="**/*", input_suffixes=[".ome.tif", ".tif"])
    out = root / "results"
    assert cli.main(["batch", str(root), "--config", str(study), "--out", str(out)]) == 0
    first = [c[0] for c in calls]
    assert first == ["i.ome.tif", "j.tif", "m.ome.tif", "w.tif"]
    calls.clear()
    assert cli.main(["batch", str(root), "--config", str(study), "--out", str(out)]) == 0
    assert [c[0] for c in calls] == first


def test_an_output_root_at_or_above_the_input_folder_keeps_the_inputs(tmp_path, monkeypatch):
    """Review S3: with io.input_suffixes, results next to the images (--out the input folder) or around them (--out
    a parent folder) still find the images; only an output root inside the input folder is skipped (D20), and
    germquant's own images are never taken."""
    study = tmp_path / "study"
    images = study / "images"
    images.mkdir(parents=True)
    for n in ("a.ome.tif", "b.ome.tif"):
        (images / n).write_bytes(b"")
    (images / "results" / "a").mkdir(parents=True)
    (images / "results" / "a" / "c.ome.tif").write_bytes(b"")          # under an output root inside the input
    (images / "results_old").mkdir()
    (images / "results_old" / "d.ome.tif").write_bytes(b"")             # beside it, its name begins as the root's
    every, outside = ["a.ome.tif", "b.ome.tif", "c.ome.tif", "d.ome.tif"], ["a.ome.tif", "b.ome.tif", "d.ome.tif"]
    monkeypatch.chdir(tmp_path)                                         # a relative --out, as cli._batch passes it
    for out, want in ((images, every), (study, every), (tmp_path, every), (tmp_path / "elsewhere", every),
                      (images / "results", outside), (images / "results" / "a", outside)):
        for skip in (out, str(out) + os.sep, os.path.relpath(out, tmp_path)):
            files, _ = batch.discover_files(images, "**/*", suffixes=(".ome.tif",), skip_root=skip)
            assert [f.name for f in files] == want, (out, skip)
    shutil.rmtree(images / "results")
    shutil.rmtree(images / "results_old")
    calls: list = []
    monkeypatch.setattr(pipeline, "process_image", _fake_process(calls))
    prof = _study(tmp_path, input_glob="**/*", input_suffixes=[".ome.tif", ".tif"])
    for out in (images, study):
        assert cli.main(["batch", str(images), "--config", str(prof), "--out", str(out)]) == 0
        first = [c[0] for c in calls]
        assert first == ["a.ome.tif", "b.ome.tif"], out
        for name, out_dir in calls:
            (out_dir / f"{image_stem(name)}__nuclei_labels.tif").write_bytes(b"")
        calls.clear()
        assert cli.main(["batch", str(images), "--config", str(prof), "--out", str(out)]) == 0
        assert [c[0] for c in calls] == first, out
        calls.clear()


def test_two_formats_of_one_id_in_one_folder_keep_the_first(tmp_path, monkeypatch, capsys):
    """Review pass 2, S1 (reviews/P5_pass2_sheet.md): with io.input_suffixes listing several formats, x.nd2 and its
    export x.tif in one folder have one image id, which names one results folder. The first by name is kept and
    every other one is dropped with the reason duplicate_id and printed, so no two runs write one results folder;
    the same id in another folder is another image, and a default batch (.nd2 only) never meets the rule."""
    root = tmp_path / "in"
    for name in ("img.nd2", "img.ome.tif", "img.tif", "sub/img.tif", "two.ome.tif", "TWO.tif", "x.nd2", "y.tif"):
        p = root / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b"")
    files, dropped = batch.discover_files(root, "**/*", suffixes=(".nd2", ".ome.tif", ".tif"))
    names = _rel(root, files)
    twins = {"two.ome.tif", "TWO.tif"}            # which of the two sorts first is the system's letter case rule
    assert [n for n in names if n not in twins] == ["img.nd2", "sub/img.tif", "x.nd2", "y.tif"]
    assert len(twins & set(names)) == 1
    why = {f.name: w for f, w in dropped}
    assert set(why) == {"img.ome.tif", "img.tif"} | (twins - set(names))
    assert why["img.tif"].startswith("duplicate_id img: img.nd2 in this folder has the same image id")
    assert all(w.startswith("duplicate_id ") for w in why.values())
    # the default set: .nd2 files only, one suffix per id, so nothing is dropped for it
    assert _rel(root, batch.discover_files(root)[0]) == ["img.nd2", "x.nd2"]
    assert not [w for _, w in batch.discover_files(root, "**/*")[1] if w.startswith("duplicate_id")]
    # `germquant batch`: one run per results folder, the others printed as left out; a resume finds every marker
    calls: list = []
    monkeypatch.setattr(pipeline, "process_image", _fake_process(calls))
    study = _study(tmp_path, input_glob="**/*", input_suffixes=[".nd2", ".ome.tif", ".tif"])
    out = tmp_path / "out"
    assert cli.main(["batch", str(root), "--config", str(study), "--out", str(out)]) == 0
    printed = capsys.readouterr().out
    assert [c[0] for c in calls] == [Path(n).name for n in names]
    assert len({c[1] for c in calls}) == len(calls)
    assert "left out img.ome.tif: duplicate_id img" in printed and "left out img.tif: duplicate_id img" in printed
    calls.clear()
    assert cli.main(["batch", str(root), "--config", str(study), "--out", str(out), "--resume"]) == 0
    assert calls == [] and capsys.readouterr().out.count("already done") == len(names)
    # one format, two names that differ in letter case (two files only where the file system tells case apart): the
    # base keeps both, and so does the rule, which compares formats
    case = tmp_path / "case"
    case.mkdir()
    (case / "Q.nd2").write_bytes(b"")
    (case / "q.nd2").write_bytes(b"")
    want = sorted(p.name for p in case.iterdir())
    assert sorted(f.name for f in batch.discover_files(case, "**/*", suffixes=(".nd2", ".tif"))[0]) == want
    # a file dropped by a pattern holds no id (review pass 2 part 2a, N2): with the .nd2 files left out by
    # io.exclude_patterns, the OME-TIFF of the pair runs
    files, dropped = batch.discover_files(root, "**/*", ["*.nd2"], suffixes=(".nd2", ".ome.tif", ".tif"))
    why = {f.name: w for f, w in dropped}
    assert "img.ome.tif" in [f.name for f in files] and why["img.nd2"] == "exclude_pattern *.nd2"
    assert why["img.tif"].startswith("duplicate_id img: img.ome.tif in this folder")
    # one id twice in one format, where the file system tells letter case apart (Linux: x.tif and x.TIF, listed here
    # as such a system lists them): one results folder, so the later is left out (review pass 2 part 2a, N1); the
    # default suffixes keep 8eb73de's files, both z.ND2 and z.nd2 (D20)
    listing = [case / "y.TIF", case / "y.tif", case / "z.ND2", case / "z.nd2"]
    monkeypatch.setattr(type(case), "glob", lambda self, pattern: iter(listing))
    kept, gone = batch.discover_files(case, "**/*", suffixes=(".nd2", ".tif"))
    assert [f.name for f in kept] == ["y.TIF", "z.ND2"] and [f.name for f, _ in gone] == ["y.tif", "z.nd2"]
    assert all(w.startswith("duplicate_id ") for _, w in gone)
    assert [f.name for f in batch.discover_files(case, "**/*")[0]] == ["z.ND2", "z.nd2"]


# ---------------------------------------------------------------------------------------------------
# io.input_suffixes: the strict check
# ---------------------------------------------------------------------------------------------------

def test_the_suffix_list_is_the_readers_list():
    assert recipe_keys.INPUT_SUFFIXES == tuple(SUPPORTED_SUFFIXES)
    assert "io.input_suffixes" in recipe_keys.IO_KEYS and "io.input_suffixes" not in recipe_keys.CHECKED_KEYS
    # the one shipped profile that sets it is valid
    assert recipe_keys.problems(load_config(ROOT / "config" / "profiles" / "generic_2d.yaml")) == []


@pytest.mark.parametrize("value", [None, [".nd2"], (".ome.tif", ".tif"), list(SUPPORTED_SUFFIXES), [".tif", ".tif"]])
def test_input_suffixes_accepted(value):
    assert recipe_keys.input_suffixes_problem(value) is None


@pytest.mark.parametrize("value, hint", [
    (".ome.tif", "a single suffix is still a list"), ([], "got []"), ([".TIF"], "lower case"), (["tif"], "with the dot"),
    ([".png"], "'.png'"), ([1], "1"), ([None], "None"), ({"a": ".nd2"}, "{'a': '.nd2'}"), ("", "got ''"),
    ([".ome.tif", ".jpg"], "'.jpg'"),
])
def test_input_suffixes_refused(value, hint):
    problem = recipe_keys.input_suffixes_problem(value)
    assert problem and problem.startswith("io.input_suffixes must be"), problem
    assert hint in problem, problem


def test_a_malformed_input_suffixes_is_refused_before_any_file(tmp_path, monkeypatch, capsys):
    root = _tree(tmp_path)
    seen = _spy_discover(monkeypatch)
    calls: list = []
    monkeypatch.setattr(pipeline, "process_image", _fake_process(calls))
    out = tmp_path / "out"
    rc = cli.main(["batch", str(root), "--config", str(_study(tmp_path, input_suffixes=".ome.tif")), "--out", str(out)])
    err = capsys.readouterr().err
    assert rc == 2 and seen == [] and calls == [] and not out.exists()
    assert "REFUSED: " in err and "io.input_suffixes must be" in err
    # the helper both cli._batch and the Snakefile use: {} for a recipe that does not set the key
    assert batch.discover_options(load_config(ROOT / "config" / "config.yaml"), out) == {}
    assert batch.discover_options(load_config(ROOT / "config" / "profiles" / "generic_2d.yaml"), out) == {
        "suffixes": (".ome.tif", ".ome.tiff", ".tif", ".tiff", ".czi", ".nd2"), "skip_root": out}
    with pytest.raises(ValueError, match="io.input_suffixes must be"):
        batch.discover_options(load_config(_study(tmp_path, input_suffixes=[".png"])), out)


def test_process_image_refuses_a_malformed_input_suffixes(tmp_path, patched_reader):  # noqa: F811
    cfg = load_config(_study(tmp_path, input_suffixes=["tif"]))
    out = tmp_path / "out"
    with pytest.raises(ValueError, match=r"recipe refused before any stage ran: io\.input_suffixes must be"):
        pipeline.process_image("20251105_n2_nohs_HERM_001.nd2", cfg, out, prov=None)
    assert not list(out.glob("*.csv")) and not (out / "run_manifest.json").exists()


# ---------------------------------------------------------------------------------------------------
# the Snakefile (D25: snakemake in a throwaway environment only)
# ---------------------------------------------------------------------------------------------------

def _snakemake() -> str:
    exe = os.environ.get("GERMQUANT_SNAKEMAKE") or shutil.which("snakemake")
    if not exe:
        pytest.skip("snakemake is not installed here: it is checked in a throwaway environment, never in the .venv "
                    "(PLAN.md D25); set GERMQUANT_SNAKEMAKE to its snakemake executable to run this test")
    return exe


def _dry_run(exe: str, home: Path, **config) -> str:
    """`snakemake -n -p` from the repository folder (the Snakefile's configfile is relative to it), leaving nothing
    behind there; the printed commands."""
    (home / "AppData" / "Local").mkdir(parents=True, exist_ok=True)
    left = {name: (ROOT / name).exists() for name in (".snakemake", "snakemake")}
    try:
        r = subprocess.run([exe, "-n", "-p", "-s", str(ROOT / "workflow" / "Snakefile"), "--cores", "1", "--config",
                            *[f"{k}={v}" for k, v in config.items()]], capture_output=True, text=True, cwd=ROOT,
                           env=subprocess_env(GERMQUANT_PACKS=TEST_PACK.parent, USERPROFILE=home, HOME=home,
                                              LOCALAPPDATA=home / "AppData" / "Local"), timeout=600)
    finally:
        for name, was in left.items():
            if not was:
                shutil.rmtree(ROOT / name, ignore_errors=True)
    assert r.returncode == 0, (r.stdout + r.stderr)[-3000:]
    return r.stdout + r.stderr


def test_snakemake_dry_run_formats_and_pack_mode(tmp_path):
    exe = _snakemake()
    root = _tree(tmp_path)
    home = tmp_path / "home"
    out = tmp_path / "out"
    # default mode: the .nd2 jobs only, under Path.stem ids, as at 8eb73de
    text = _dry_run(exe, home, input_root=root, output_root=out, germquant_config=_study(tmp_path).as_posix())
    jobs = [line for line in text.splitlines() if line.strip().startswith("germquant run")]
    want = _default_set(root)
    assert len(jobs) == len(want) and all("--config" in j and ".ome" not in j for j in jobs), jobs
    for f in want:
        assert any(f'"{f}"' in j and f'--out "{out / f.relative_to(root).parent / f.stem}"' in j for j in jobs), f
    # io.input_suffixes: the OME-TIFF jobs under their image_stem ids, nothing from the output root
    study = _study(tmp_path, input_glob="**/*", input_suffixes=[".ome.tif"])
    text = _dry_run(exe, home, input_root=root, output_root=root / "results", germquant_config=study.as_posix())
    jobs = [line for line in text.splitlines() if line.strip().startswith("germquant run")]
    assert len(jobs) == 2 and all(".ome.tif" in j.split("--out")[0] for j in jobs), jobs
    assert {j.split("--out")[1].split('"')[1] for j in jobs} == {str(root / "results" / "i"),
                                                                  str(root / "results" / "sub" / "m")}
    # pack mode: germquant_pack ID@VER with recipe R runs `germquant run --pack`, discovering with that recipe's keys;
    # germquant_config, which workflow/config.yaml always sets, is not read for it (here it would match no file;
    # review N4)
    nothing = _study(tmp_path, input_glob="no_such_folder/*.nd2")
    text = _dry_run(exe, home, input_root=root, output_root=out, germquant_pack="test_pack@0.1.0", recipe="basic",
                    germquant_config=nothing.as_posix())
    jobs = [line for line in text.splitlines() if line.strip().startswith("germquant run")]
    assert jobs and all('--pack "test_pack@0.1.0" --recipe "basic"' in j and "--config" not in j for j in jobs), jobs
