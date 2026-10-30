"""The sample sheet (docs/DESIGN_system_packs.md 7.4; PLAN.md 4.4 "P5-SHEET" and decision D21; P5 acceptance 5;
docs/PROJECTS.md, "Sample sheet").

* The synthetic_sheet golden (tests/golden/synthetic_sheet/, its sentinel pinned in tests/golden/default_invariants.json):
  a two-row sheet (samples.csv, written by write_sample_sheet) on two synthetic images, one of them left out of it.
  The image in the sheet takes the sheet's design columns in every per-image table (germ_cell follows the sheet's
  sex); the one left out keeps the file name's (metadata.filename_regex) and raises metadata:not_in_sample_sheet.
  Everything else equals the synthetic golden, a run without a sheet. records.json pins the row sha256 and the
  sheet records of run_geometry and the done marker; batch_image_summary.csv the sheet's extra column joined at collect.
* `germquant batch --samples --resume`: a new row for another image keeps every image done; editing one row's design
  cell reruns that image only; a row added for the image left out reruns it; an Excel re-save (CRLF, a byte-order
  mark, another column order, padded and quoted cells, cp1252) reruns nothing.
* The row hash (D21) is the sha256 of the canonical JSON of the row's design cells; the file is read as utf-8-sig,
  else cp1252 with a warning; a non-ISO acquisition_date is kept as given, with a warning; write_sample_sheet writes
  ASCII with LF line ends atomically and refuses a non-ASCII value, naming its cell.
* Default runs are unchanged: no sample keyword reaches process_image or is_done unless a sheet is given (D40), and no
  sheet record, flag or column is written without one. `run --t ... --samples` is refused before any timepoint.
* batch.is_done compares the row sha256 (reasons starting "sample row"), and in a labels run its segment reason names
  the labels file (D52).

Regenerate the golden only for a deliberate change of this path: ``GERMQUANT_UPDATE_GOLDEN=1 pytest
tests/test_sample_sheet.py``, as for tests/test_golden.py.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import logging
import os
import shutil
import warnings
from pathlib import Path

import pandas as pd
import pytest
from labels_helpers import TEMPLATE, assert_tables_equal, read_table
from project_helpers import heavy_after
from test_end_to_end import CHANNELS, SPACING, _synthetic_stack

from germquant import api, batch, cli, pipeline, planar, provenance, schema
from germquant.config import load_config
from germquant.io import labels as L
from germquant.io import sample_sheet as S
from germquant.io.sample_metadata import parse_sample

ROOT = Path(__file__).resolve().parents[1]
GOLD = ROOT / "tests" / "golden" / "synthetic_sheet"
SYNTHETIC = ROOT / "tests" / "golden" / "synthetic"
IN = "20251105_n2_nohs_HERM_001.nd2"            # has a row in the golden sheet
OUT_OF = "20251021_n2_hs_MALE_2.nd2"             # left out of it: the file name's regex gives its design columns
IN_ID, OUT_ID = Path(IN).stem, Path(OUT_OF).stem
FLAG = "metadata:not_in_sample_sheet"
DESIGN = ("genotype", "sex", "germ_cell", "treatment", "replicate", "acquisition_date")
GOLDEN_ROWS = [
    {"image_id": IN_ID, "genotype": "syp-2(ok307)", "sex": "male", "treatment": "control", "replicate": "3",
     "acquisition_date": "2025-11-05", "slide": "S12"},
    {"image_id": "20251105_n2_nohs_HERM_009", "genotype": "N2", "sex": "herm", "treatment": "heat", "replicate": "1",
     "acquisition_date": "2025-11-05", "slide": "S13"},
]
UPDATE = bool(os.environ.get("GERMQUANT_UPDATE_GOLDEN"))
E_ACUTE = chr(0xE9)                              # built with chr, so this source stays ASCII


def _patch_reader(mp) -> None:
    """process_image reads the 3-channel synthetic gonad for any path (tests/test_end_to_end.py patched_reader)."""
    stack = _synthetic_stack()
    meta = {"sizes": {"C": 3, "Z": 20, "Y": 64, "X": 160}, "spacing": SPACING, "spacing_ok": True,
            "channel_names": CHANNELS, "dtype": "float32", "is_2d": False}
    mp.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: meta)
    mp.setattr(pipeline, "read_stack", lambda *a, **k: stack)


def _golden_cfg():
    """config.yaml with the classical segmenter, set as tests/test_golden.py sets it (its config_hash, the synthetic
    golden's 84bdda5a5230: the sheet is not hashed)."""
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    return cfg


def _fast_cfg():
    """config.yaml, classical, without the spots stage and the montage: the quickest run that writes every record."""
    cfg = _golden_cfg()
    cfg.set("spots.enabled", False)
    cfg.set("render.montage", False)
    return cfg


def _json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _write_csv(path: Path, rows: list[dict], *, columns=None, encoding="utf-8", newline="\n", quote_all=False,
               pad="") -> Path:
    """A sheet as a spreadsheet program might save it: any encoding, line end, column order, quoting and padding."""
    columns = list(columns or rows[0])
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator=newline, quoting=csv.QUOTE_ALL if quote_all else csv.QUOTE_MINIMAL)
    w.writerow([f"{pad}{c}{pad}" for c in columns])
    for r in rows:
        w.writerow([f"{pad}{r.get(c, '')}{pad}" if r.get(c, "") != "" else "" for c in columns])
    Path(path).write_bytes(buf.getvalue().encode(encoding))
    return Path(path)


def _canonical(fields: dict) -> str:
    return hashlib.sha256(json.dumps(fields, sort_keys=True, ensure_ascii=True, separators=(",", ":"))
                          .encode("ascii")).hexdigest()


# ---------------------------------------------------------------------------------------------------
# the synthetic_sheet golden and acceptance 5 on its runs
# ---------------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def golden_runs(tmp_path_factory):
    """The golden: the two synthetic images through process_image with their rows of the golden sheet, then collect
    with the sheet. The sheet is written by write_sample_sheet when the golden is regenerated."""
    work = tmp_path_factory.mktemp("sheet_golden")
    sheet_path = GOLD / "samples.csv"
    if UPDATE:
        GOLD.mkdir(parents=True, exist_ok=True)
        S.write_sample_sheet(sheet_path, GOLDEN_ROWS)
    sheet = S.read_sample_sheet(sheet_path)
    cfg = _golden_cfg()
    out = work / "results"
    runs = {}
    with pytest.MonkeyPatch.context() as mp:
        _patch_reader(mp)
        for key, name in (("in_sheet", IN), ("not_in_sheet", OUT_OF)):
            iid = Path(name).stem
            row = sheet.row(iid, work / name)
            res = pipeline.process_image(work / name, cfg, out / iid, prov=None, sample=row)
            runs[key] = {"id": iid, "row": row, "out": out / iid, "res": res}
    counts = batch.collect(out, sample_sheet=sheet)
    return {"sheet": sheet, "cfg": cfg, "out": out, "runs": runs, "counts": counts}


def _records(golden_runs) -> dict:
    rec = {}
    for key, run in golden_runs["runs"].items():
        d, iid = run["out"], run["id"]
        stages, manifest = _json(d / f"{iid}__stages.json"), _json(d / "run_manifest.json")
        assert manifest["sample_sheet"] == stages["sample_sheet"]
        geo = stages["sample_sheet"]
        rec[key] = {"image_id": iid, "row_sha256": run["row"].row_sha256,
                    "sample_sheet": {"source": Path(geo["source"]).name, "row_sha256": geo["row_sha256"]},
                    "marker_sample_row_sha256": provenance.read_done_marker(d, iid)["sample_row_sha256"],
                    "read_flags": [f for f in stages["stages"]["read"]["flags"] if f.startswith("metadata:")]}
    return rec


def test_synthetic_sheet_golden(golden_runs):
    out, cfg = golden_runs["out"], golden_runs["cfg"]
    files = {f"{key}/image_summary.csv": run["out"] / f"{run['id']}__image_summary.csv"
             for key, run in golden_runs["runs"].items()}
    files["batch_image_summary.csv"] = out / "batch_image_summary.csv"
    records = _records(golden_runs)
    if UPDATE:
        for rel, f in files.items():
            (GOLD / rel).parent.mkdir(parents=True, exist_ok=True)
            read_table(f).to_csv(GOLD / rel, index=False)
        (GOLD / "records.json").write_bytes((json.dumps(records, indent=1, sort_keys=True) + "\n").encode("ascii"))
        (GOLD / "config_hash.txt").write_text(cfg.hash + "\n", encoding="utf-8")
        pytest.skip("synthetic_sheet golden regenerated")
    sentinel = (GOLD / "config_hash.txt").read_text(encoding="utf-8").strip()
    assert cfg.hash == sentinel, "config/config.yaml or its channel map changed: the golden's config_hash no longer holds"
    for rel, f in files.items():
        ref = pd.read_csv(GOLD / rel, low_memory=False, float_precision="round_trip")
        assert_tables_equal(read_table(f), ref, rel)
    assert records == _json(GOLD / "records.json")
    # what the golden shows: the row's design columns (germ_cell from its sex) and the regex fallback with the flag
    got_in = pd.read_csv(GOLD / "in_sheet" / "image_summary.csv", dtype=str).iloc[0]
    assert {c: got_in[c] for c in DESIGN} == {"genotype": "syp-2(ok307)", "sex": "male", "germ_cell": "spermatocyte",
                                             "treatment": "control", "replicate": "3",
                                             "acquisition_date": "2025-11-05"}
    got_out = pd.read_csv(GOLD / "not_in_sheet" / "image_summary.csv", dtype=str).iloc[0]
    assert {c: got_out[c] for c in DESIGN} == {"genotype": "N2", "sex": "male", "germ_cell": "spermatocyte",
                                              "treatment": "heat", "replicate": "2", "acquisition_date": "20251021"}
    assert FLAG in got_out["qc_flags"].split(";") and FLAG not in got_in["qc_flags"].split(";")
    stacked = pd.read_csv(GOLD / "batch_image_summary.csv", dtype=str, keep_default_na=False).set_index("image_id")
    assert list(stacked.columns)[-1] == "slide" and stacked.loc[IN_ID, "slide"] == "S12"
    assert stacked.loc[OUT_ID, "slide"] == ""
    assert records["in_sheet"]["row_sha256"] == _canonical({k: v for k, v in GOLDEN_ROWS[0].items() if k in DESIGN})
    assert records["not_in_sheet"]["row_sha256"] is None and records["not_in_sheet"]["read_flags"] == [FLAG]


def test_sheet_changes_only_the_design_columns(golden_runs):
    """Both runs equal the synthetic golden (the image without a sheet) in every other column."""
    for key, run in golden_runs["runs"].items():
        for table in ("image_summary", "nuclei", "spots"):
            got = read_table(run["out"] / f"{run['id']}__{table}.csv")
            ref = pd.read_csv(SYNTHETIC / f"{table}.csv", low_memory=False, float_precision="round_trip")
            keep = [c for c in ref.columns if c not in DESIGN + ("image_id", "qc_flags")]
            assert_tables_equal(got[keep], ref[keep], f"{key} {table}", append_only=False)
        flags = set(read_table(run["out"] / f"{run['id']}__image_summary.csv")["qc_flags"].iloc[0].split(";"))
        base = set(pd.read_csv(SYNTHETIC / "image_summary.csv")["qc_flags"].iloc[0].split(";"))
        assert flags == (base | {FLAG} if key == "not_in_sheet" else base)


def test_design_columns_from_sheet(golden_runs):
    run = golden_runs["runs"]["in_sheet"]
    want = {"genotype": "syp-2(ok307)", "sex": "male", "germ_cell": "spermatocyte", "treatment": "control",
            "replicate": "3", "acquisition_date": "2025-11-05"}
    written = 0
    for table in schema.TABLES:
        df = pd.read_csv(run["out"] / f"{run['id']}__{table}.csv", dtype=str, keep_default_na=False)
        assert set(want) <= set(df.columns), table
        if len(df):
            written += 1
            assert (df["image_id"] == IN_ID).all(), table
            for col, value in want.items():
                assert (df[col] == value).all(), (table, col)
    assert written >= 3
    # germ_cell came from the sheet's sex: the file name says HERM (oocyte), the row says male
    assert parse_sample(IN, golden_runs["cfg"].get("metadata.filename_regex"))["germ_cell"] == "oocyte"
    # the records of a sheet run: run_geometry (manifest and stage record) and the done marker
    row = run["row"]
    assert row.row_sha256 == _canonical({k: GOLDEN_ROWS[0][k] for k in DESIGN if k in GOLDEN_ROWS[0]})
    assert row.fields == {k: GOLDEN_ROWS[0][k] for k in DESIGN if k in GOLDEN_ROWS[0]}
    stages = _json(run["out"] / f"{IN_ID}__stages.json")
    assert stages["sample_sheet"] == {"source": os.path.abspath(GOLD / "samples.csv"), "row_sha256": row.row_sha256}
    assert _json(run["out"] / "run_manifest.json")["sample_sheet"] == stages["sample_sheet"]
    assert provenance.read_done_marker(run["out"], IN_ID)["sample_row_sha256"] == row.row_sha256
    assert FLAG not in run["res"]["qc_flags"] and FLAG not in stages["stages"]["read"]["flags"]
    # the sheet is not hashed: the config hash is the synthetic golden's, and the image with a row and the image
    # without one (the same pixels and recipe) have the same stage hashes
    pins = _json(ROOT / "tests" / "golden" / "default_invariants.json")
    assert stages["config_hash"] == pins["golden_sentinels"]["synthetic"]
    left = golden_runs["runs"]["not_in_sheet"]
    assert provenance.read_done_marker(run["out"], IN_ID)["stage_hashes"] == provenance.read_done_marker(
        left["out"], left["id"])["stage_hashes"]


def test_missing_image_falls_back_with_flag_only_in_sheet_runs(golden_runs, tmp_path, monkeypatch):
    run, cfg = golden_runs["runs"]["not_in_sheet"], golden_runs["cfg"]
    regex = parse_sample(OUT_OF, cfg.get("metadata.filename_regex"), cfg.get("metadata.defaults"))
    assert regex == {"image_id": OUT_ID, "genotype": "N2", "sex": "male", "germ_cell": "spermatocyte",
                     "treatment": "heat", "replicate": "2", "acquisition_date": "20251021"}
    summary = pd.read_csv(run["out"] / f"{OUT_ID}__image_summary.csv", dtype=str).iloc[0]
    assert {c: summary[c] for c in regex} == regex
    assert run["row"] == S.SampleRow({}, os.path.abspath(GOLD / "samples.csv"), None)
    assert FLAG in run["res"]["qc_flags"] and FLAG in summary["qc_flags"].split(";")
    stages = _json(run["out"] / f"{OUT_ID}__stages.json")
    assert FLAG in stages["stages"]["read"]["flags"]
    assert stages["sample_sheet"] == {"source": os.path.abspath(GOLD / "samples.csv"), "row_sha256": None}
    marker = provenance.read_done_marker(run["out"], OUT_ID)
    assert "sample_row_sha256" in marker and marker["sample_row_sha256"] is None
    assert run["res"]["qc_pass"] is True                     # the flag is advisory
    # the same image without a sheet: no flag, no record
    _patch_reader(monkeypatch)
    res = pipeline.process_image(OUT_OF, _fast_cfg(), tmp_path / "plain", prov=None)
    assert FLAG not in res["qc_flags"]
    assert "sample_sheet" not in _json(tmp_path / "plain" / f"{OUT_ID}__stages.json")
    assert "sample_row_sha256" not in provenance.read_done_marker(tmp_path / "plain", OUT_ID)


def test_default_runs_unchanged(tmp_path, monkeypatch, capsys):
    """D40 and the fence: without a sheet nothing new is passed on, written or flagged."""
    _patch_reader(monkeypatch)
    out = tmp_path / "plain"
    res = pipeline.process_image(IN, _fast_cfg(), out, prov=None)
    stages, manifest = _json(out / f"{IN_ID}__stages.json"), _json(out / "run_manifest.json")
    assert "sample_sheet" not in stages and "sample_sheet" not in manifest and FLAG not in res["qc_flags"]
    assert sorted(provenance.read_done_marker(out, IN_ID)) == ["config_hash", "enabled_stages", "failed_stages",
                                                                "finished", "image_id", "stage_hashes", "xy_stride",
                                                                "z_range"]
    # api.run and run_loaded pass sample on only when it is given
    seen = []

    def fake_process_image(path, cfg, out_dir, *, xy_stride=1, z_range=None, prov=None, **kw):
        seen.append(kw)
        return {"image_id": "x", "n_nuclei": 0, "qc_pass": True, "qc_flags": [], "out_dir": str(out_dir),
                "tables": {}, "stages": {}}

    monkeypatch.setattr(api, "process_image", fake_process_image)
    cfg = load_config(ROOT / "config" / "config.yaml")
    row = S.SampleRow({"genotype": "N2"}, "s.csv", _canonical({"genotype": "N2"}))
    api.run("g.nd2", config=cfg, out=tmp_path / "a")
    api.run("g.nd2", config=cfg, out=tmp_path / "b", sample=row)
    api.run_loaded("g.nd2", cfg, tmp_path / "c")
    api.run_loaded("g.nd2", cfg, tmp_path / "d", sample=row)
    assert seen == [{}, {"sample": row}, {}, {"sample": row}]
    # review N3: without a row, api.run hands run_loaded today's keywords (a fake of today's exact signature would
    # raise on any other)
    loaded = []

    def exact_loaded(image, cfg, out=None, *, xy_stride=None, z_range=None, t=None):
        loaded.append(str(image))
        return {}

    with monkeypatch.context() as mp:
        mp.setattr(api, "run_loaded", exact_loaded)
        api.run("g.nd2", config=cfg, out=tmp_path / "e")
    assert loaded == ["g.nd2"]
    # `germquant batch` without --samples: process_image and is_done get today's keywords (a fake with today's exact
    # signature would raise on any other), and the batch writes what it always has
    src = tmp_path / "in"
    src.mkdir()
    for n in (IN, OUT_OF):
        (src / n).write_bytes(b"")
    calls, decisions = [], []

    def exact(nd2, cfg, out_dir, *, xy_stride=1, z_range=None, prov=None):
        calls.append(Path(nd2).name)
        return pipeline_real(nd2, cfg, out_dir, xy_stride=xy_stride, z_range=z_range, prov=prov)

    pipeline_real = pipeline.process_image
    monkeypatch.setattr(pipeline, "process_image", exact)
    real_done = batch.is_done
    monkeypatch.setattr(batch, "is_done", lambda *a, **k: decisions.append(sorted(k)) or real_done(*a, **k))
    study = _study(tmp_path)
    bout = tmp_path / "out"
    assert cli.main(["batch", str(src), "--config", str(study), "--out", str(bout)]) == 0
    assert calls == [OUT_OF, IN]
    assert cli.main(["batch", str(src), "--config", str(study), "--out", str(bout), "--resume"]) == 0
    assert decisions == [["segment_hash", "xy_stride"]] * 2
    assert "sample sheet" not in capsys.readouterr().out
    # a default collect joins nothing
    cols = pd.read_csv(bout / "batch_image_summary.csv").columns
    assert "slide" not in cols and list(pd.read_csv(bout / "batch_summary.csv").columns) == batch.SUMMARY_COLS


def test_samples_with_t_refused(tmp_path, monkeypatch, capsys):
    """D21: `run --t ... --samples` is refused before any timepoint (the series loop passes no sample row), with
    "REFUSED:" and exit 2 (D45), after the labels and single-plane refusals."""
    import germquant.timelapse as T

    series: list = []
    monkeypatch.setattr(T, "process_series", lambda *a, **k: series.append(a) or {})
    out = tmp_path / "out"
    missing = tmp_path / "not_read.csv"                       # never opened: the refusal comes first
    rc = cli.main(["run", "series.nd2", "--config", str(ROOT / "config" / "config.yaml"), "--out", str(out),
                   "--t", "0-1", "--samples", str(missing)])
    err = capsys.readouterr().err
    assert rc == 2 and series == [] and f"REFUSED: {S.SERIES_REFUSAL}" in err
    assert not out.exists() or not any(out.iterdir())
    base = load_config(ROOT / "config" / "config.yaml")
    assert cli._series_refusal(base, "s.nd2") is None
    assert cli._series_refusal(base, "s.nd2", samples=str(missing)) == S.SERIES_REFUSAL
    flat = load_config(_study(tmp_path, planar="true"))
    assert cli._series_refusal(flat, "s.nd2", samples=str(missing)) == planar.PLANAR_SERIES_REFUSED
    labels = load_config(_study(tmp_path, labels=True))
    assert cli._series_refusal(labels, "s.nd2", samples=str(missing)) == L.SERIES_REFUSAL
    # process_image refuses a sheet row with t before any pixel is read
    reads: list = []
    monkeypatch.setattr(pipeline, "read_nd2_metadata", lambda *a, **k: reads.append(a) or {})
    monkeypatch.setattr(pipeline, "read_stack", lambda *a, **k: reads.append(a))
    row = S.SampleRow({"sex": "male"}, "s.csv", _canonical({"sex": "male"}))
    with pytest.raises(S.SampleSheetError, match="time series"):
        pipeline.process_image(IN, base, tmp_path / "t", prov=None, t=0, sample=row)
    with pytest.raises(S.SampleSheetError, match="time series"):
        api.run_loaded(IN, base, tmp_path / "t", t=0, sample=row)
    assert reads == []
    # and anything but a SampleRow is refused there too, before the read
    with pytest.raises(TypeError, match="SampleRow"):
        pipeline.process_image(IN, base, tmp_path / "t", prov=None, sample={"sex": "male"})
    with pytest.raises(S.SampleSheetError, match="row_sha256"):
        pipeline.process_image(IN, base, tmp_path / "t", prov=None, sample=S.SampleRow({"sex": "male"}, "s", "0" * 64))
    with pytest.raises(S.SampleSheetError, match="design columns"):
        pipeline.process_image(IN, base, tmp_path / "t", prov=None, sample=S.SampleRow({"strain": "x"}, "s", None))
    assert reads == []


def test_run_takes_the_row_of_the_image(tmp_path, monkeypatch, capsys):
    seen = []

    def fake_run(image, **kw):
        seen.append((Path(image).name, kw))
        return {"image_id": Path(image).stem, "n_nuclei": 1, "qc_pass": True, "qc_flags": [], "out_dir": "x"}

    monkeypatch.setattr(api, "run", fake_run)
    sheet = _write_csv(tmp_path / "samples.csv", [{"image_id": IN_ID, "genotype": "N2", "sex": "herm"}])
    cfg = str(ROOT / "config" / "config.yaml")
    img = tmp_path / IN
    assert cli.main(["run", str(img), "--config", cfg, "--out", str(tmp_path / "o"), "--samples", str(sheet)]) == 0
    assert seen[-1][1]["sample"] == S.read_sample_sheet(sheet).row(IN_ID, img)
    assert seen[-1][1]["sample"].fields == {"genotype": "N2", "sex": "herm"}
    capsys.readouterr()
    assert cli.main(["run", str(tmp_path / OUT_OF), "--config", cfg, "--out", str(tmp_path / "o"),
                     "--samples", str(sheet)]) == 0
    assert seen[-1][1]["sample"].row_sha256 is None
    assert "not in the sample sheet" in capsys.readouterr().out
    assert cli.main(["run", str(img), "--config", cfg, "--out", str(tmp_path / "o")]) == 0
    assert "sample" not in seen[-1][1]
    # an OME-TIFF takes the row of its image id, the file name without .ome.tif (review pass 2, N1)
    ome = _write_csv(tmp_path / "ome.csv", [{"image_id": "gonad_7", "genotype": "him-8(e1489)"}])
    assert cli.main(["run", str(tmp_path / "gonad_7.ome.tif"), "--config", cfg, "--out", str(tmp_path / "o"),
                     "--samples", str(ome)]) == 0
    assert seen[-1][1]["sample"].fields == {"genotype": "him-8(e1489)"}
    # an unreadable sheet is refused before the run starts (exit 2, REFUSED, D49)
    bad = tmp_path / "bad.csv"
    bad.write_bytes(b"genotype,sex\nN2,herm\n")
    n = len(seen)
    with pytest.raises(SystemExit) as e:
        cli.main(["run", str(img), "--config", cfg, "--out", str(tmp_path / "o"), "--samples", str(bad)])
    assert e.value.code == 2 and len(seen) == n
    assert "REFUSED: " in capsys.readouterr().err


# ---------------------------------------------------------------------------------------------------
# `germquant batch --samples --resume` (acceptance 5)
# ---------------------------------------------------------------------------------------------------

BATCH = ["20251105_n2_nohs_HERM_001.nd2", "20251105_n2_nohs_HERM_002.nd2", "20251021_n2_hs_MALE_2.nd2"]
BATCH_ROWS = [
    {"image_id": "20251105_n2_nohs_HERM_001", "genotype": "unc-119(ed3)", "sex": "herm",
     "treatment": f"t{E_ACUTE}moin", "replicate": "1", "notes": "first"},
    {"image_id": "20251105_n2_nohs_HERM_002", "genotype": "unc-119(ed3)", "sex": "herm", "treatment": "control",
     "replicate": "2", "notes": "second"},
    {"image_id": "20251105_n2_nohs_HERM_077", "genotype": "N2", "sex": "male", "treatment": "heat", "replicate": "7",
     "notes": "not in the folder"},
]


def _study(tmp_path: Path, planar: str | None = None, labels: bool = False) -> Path:
    """A study profile outside the repository: config.yaml, the classical segmenter, no spots, no montage."""
    lines = [f"extends: {(ROOT / 'config' / 'config.yaml').as_posix()}", "segmentation:", "  nuclei:",
             "    method: classical"]
    if labels:
        lines.append(f"    labels_from: '{TEMPLATE}'")
    lines += ["spots:", "  enabled: false", "render:", "  montage: false"]
    if planar is not None:
        lines += ["io:", f"  planar: {planar}"]
    study = tmp_path / ("study_planar.yaml" if planar else "study_labels.yaml" if labels else "study.yaml")
    study.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return study


@pytest.fixture(scope="module")
def batch_runs(tmp_path_factory):
    """One `germquant batch --samples` of three stand-ins: two with rows, one left out of the sheet."""
    work = tmp_path_factory.mktemp("sheet_batch")
    src = work / "in"
    src.mkdir()
    for n in BATCH:
        (src / n).write_bytes(b"")
    study = _study(work)
    sheet = _write_csv(work / "samples.csv", BATCH_ROWS)
    out = work / "out"
    with pytest.MonkeyPatch.context() as mp:
        _patch_reader(mp)
        assert cli.main(["batch", str(src), "--config", str(study), "--out", str(out), "--samples", str(sheet)]) == 0
    return {"src": src, "study": study, "sheet": sheet, "out": out}


def _resume(batch_runs, tmp_path, monkeypatch, capsys, sheet: Path) -> tuple[list[str], str]:
    """`batch --resume --samples SHEET` on a copy of the fixture's results: (the images processed, the output)."""
    out = tmp_path / "out"
    shutil.copytree(batch_runs["out"], out)
    _patch_reader(monkeypatch)
    calls: list[str] = []
    real = pipeline.process_image

    def spy(nd2, *a, **k):
        calls.append(Path(nd2).name)
        return real(nd2, *a, **k)

    monkeypatch.setattr(pipeline, "process_image", spy)
    capsys.readouterr()
    assert cli.main(["batch", str(batch_runs["src"]), "--config", str(batch_runs["study"]), "--out", str(out),
                     "--resume", "--samples", str(sheet)]) == 0
    return calls, capsys.readouterr().out


def test_batch_with_samples_applies_rows(batch_runs):
    out = batch_runs["out"]
    for r in BATCH_ROWS[:2]:
        iid = r["image_id"]
        s = pd.read_csv(out / iid / f"{iid}__image_summary.csv", dtype=str).iloc[0]
        assert (s["genotype"], s["treatment"], s["replicate"]) == (r["genotype"], r["treatment"], r["replicate"])
        assert FLAG not in str(s["qc_flags"])
        assert provenance.read_done_marker(out / iid, iid)["sample_row_sha256"] == _canonical(
            {k: v for k, v in r.items() if k in DESIGN})
    left = Path(BATCH[2]).stem
    s = pd.read_csv(out / left / f"{left}__image_summary.csv", dtype=str).iloc[0]
    assert s["genotype"] == "N2" and FLAG in s["qc_flags"].split(";")
    assert provenance.read_done_marker(out / left, left)["sample_row_sha256"] is None


def test_new_row_for_other_image_keeps_done(batch_runs, tmp_path, monkeypatch, capsys):
    rows = BATCH_ROWS + [{"image_id": "20251105_n2_nohs_HERM_078", "genotype": "N2", "sex": "herm",
                          "treatment": "heat", "replicate": "8", "notes": "another"}]
    calls, printed = _resume(batch_runs, tmp_path, monkeypatch, capsys, _write_csv(tmp_path / "s.csv", rows))
    assert calls == [] and printed.count("already done") == 3


def test_row_added_for_the_image_left_out_reruns_it(batch_runs, tmp_path, monkeypatch, capsys):
    left = Path(BATCH[2]).stem
    rows = BATCH_ROWS + [{"image_id": left, "genotype": "him-8(e1489)", "sex": "male", "treatment": "heat",
                          "replicate": "2", "notes": "added"}]
    calls, printed = _resume(batch_runs, tmp_path, monkeypatch, capsys, _write_csv(tmp_path / "s.csv", rows))
    assert calls == [BATCH[2]] and printed.count("already done") == 2
    assert "reprocessing: sample row differs" in printed
    s = pd.read_csv(tmp_path / "out" / left / f"{left}__image_summary.csv", dtype=str).iloc[0]
    assert s["genotype"] == "him-8(e1489)" and FLAG not in str(s["qc_flags"])


def test_edited_row_restales_only_that_image(batch_runs, tmp_path, monkeypatch, capsys):
    rows = [dict(r) for r in BATCH_ROWS]
    rows[1]["genotype"] = "unc-119(ed3); him-5(e1490)"        # a design cell of the second image
    rows[0]["notes"] = "edited note"                          # an extra column: joined at collect, not run on
    calls, printed = _resume(batch_runs, tmp_path, monkeypatch, capsys, _write_csv(tmp_path / "s.csv", rows))
    assert calls == [BATCH[1]] and printed.count("already done") == 2
    assert "reprocessing: sample row differs (this image's row in the sample sheet changed)" in printed
    iid = rows[1]["image_id"]
    s = pd.read_csv(tmp_path / "out" / iid / f"{iid}__image_summary.csv", dtype=str).iloc[0]
    assert s["genotype"] == "unc-119(ed3); him-5(e1490)"
    stacked = pd.read_csv(tmp_path / "out" / "batch_image_summary.csv", dtype=str).set_index("image_id")
    assert stacked.loc[rows[0]["image_id"], "notes"] == "edited note"     # the resumed batch collects with the sheet


@pytest.mark.parametrize("form", ["crlf_bom_reordered_padded_quoted", "cp1252"])
def test_excel_resave_keeps_every_image_done(batch_runs, tmp_path, monkeypatch, capsys, form):
    columns = list(reversed(list(BATCH_ROWS[0])))
    if form == "cp1252":
        sheet = _write_csv(tmp_path / "s.csv", BATCH_ROWS, encoding="cp1252", newline="\r\n")
    else:
        sheet = _write_csv(tmp_path / "s.csv", BATCH_ROWS, columns=columns, encoding="utf-8-sig", newline="\r\n",
                           quote_all=True, pad=" ")
    calls, printed = _resume(batch_runs, tmp_path, monkeypatch, capsys, sheet)
    assert calls == [] and printed.count("already done") == 3
    assert ("cp1252" in printed) == (form == "cp1252")


# ---------------------------------------------------------------------------------------------------
# batch.is_done: the row sha256, and the labels reason (D52)
# ---------------------------------------------------------------------------------------------------

def _marker(out, iid, cfg, extra=None, segment="abc"):
    provenance.write_done_marker(out, iid, config_hash=cfg.hash, enabled_stages=batch.enabled_stage_names(cfg),
                                 stage_hashes={"segment": segment}, extra=extra)


def test_is_done_compares_the_sample_row(tmp_path):
    cfg = load_config(ROOT / "config" / "config.yaml")
    x, y = "1" * 64, "2" * 64
    _marker(tmp_path, "row", cfg, {"sample_row_sha256": x})
    assert batch.is_done(tmp_path, "row", cfg, sample_row_sha256=x) == (True, "done")
    assert batch.is_done(tmp_path, "row", cfg, segment_hash="abc", sample_row_sha256=x) == (True, "done")
    assert batch.is_done(tmp_path, "row", cfg, sample_row_sha256=y) == (
        False, "sample row differs (this image's row in the sample sheet changed)")
    done, why = batch.is_done(tmp_path, "row", cfg)                  # a sheet removed, or the image's row deleted
    assert done is False and why.startswith("sample row: the image ran with a sample sheet row")
    _marker(tmp_path, "left", cfg, {"sample_row_sha256": None})      # a sheet run, the image not in the sheet
    assert batch.is_done(tmp_path, "left", cfg) == (True, "done")
    assert batch.is_done(tmp_path, "left", cfg, sample_row_sha256=x) == (
        False, "sample row differs (the image had no row in the sample sheet when it ran)")
    _marker(tmp_path, "plain", cfg)                                  # a run without a sheet: today's marker
    assert batch.is_done(tmp_path, "plain", cfg) == (True, "done")
    assert batch.is_done(tmp_path, "plain", cfg, sample_row_sha256=x) == (
        False, "sample row differs (the image ran without a sample sheet)")
    # the row is checked before the segment hash, also for a single-plane marker
    _marker(tmp_path, "flat", cfg, {"planar": True, "sample_row_sha256": x})
    assert batch.is_done(tmp_path, "flat", cfg, segment_hash="zzz", segment_hash_planar="abc",
                         sample_row_sha256=x) == (True, "done")
    assert batch.is_done(tmp_path, "flat", cfg, segment_hash="zzz", segment_hash_planar="abc",
                         sample_row_sha256=y)[1].startswith("sample row differs")
    # a stale marker of another kind reports that first
    assert batch.is_done(tmp_path, "row", cfg, xy_stride=4, sample_row_sha256=y)[1].startswith("xy_stride")


def test_labels_resume_reason_names_the_labels_file(tmp_path):
    base = load_config(ROOT / "config" / "config.yaml")
    old = "segment stage hash differs (model file or Cellpose version changed)"
    new = "segment stage hash differs (the labels file, or the Cellpose or torch version, changed)"
    _marker(tmp_path, "a", base)
    assert batch.is_done(tmp_path, "a", base, segment_hash="zzz") == (False, old)
    labels = load_config(_study(tmp_path, labels=True))
    _marker(tmp_path, "b", labels)
    assert batch.is_done(tmp_path, "b", labels, segment_hash="zzz") == (False, new)
    assert batch.is_done(tmp_path, "b", labels, segment_hash="abc") == (True, "done")
    flat = load_config(_study(tmp_path, planar="true"))
    _marker(tmp_path, "c", flat, {"planar": True})
    assert batch.is_done(tmp_path, "c", flat, segment_hash="x", segment_hash_planar="zzz") == (False, old)
    flat.set("segmentation.nuclei.labels_from", TEMPLATE)
    _marker(tmp_path, "d", flat, {"planar": True})
    assert batch.is_done(tmp_path, "d", flat, segment_hash="x", segment_hash_planar="zzz") == (False, new)
    assert batch.is_done(tmp_path, "d", flat, segment_hash="x")[1].startswith("planar flag")


# ---------------------------------------------------------------------------------------------------
# reading and writing the sheet (D21)
# ---------------------------------------------------------------------------------------------------

def test_row_hash_is_the_canonical_json_of_the_design_cells(tmp_path):
    sheet = S.read_sample_sheet(_write_csv(tmp_path / "s.csv", [
        {"image_id": "a", "genotype": " N2 ", "sex": "male", "treatment": "", "notes": "x"},
        {"image_id": "b", "genotype": "", "sex": "", "treatment": "", "notes": "only a note"}]))
    a = sheet.row("a")
    assert a.fields == {"genotype": "N2", "sex": "male"}                       # stripped; an empty cell is left out
    assert a.row_sha256 == _canonical({"genotype": "N2", "sex": "male"}) == S.row_sha256(a.fields)
    assert a.row_sha256 == hashlib.sha256(b'{"genotype":"N2","sex":"male"}').hexdigest()
    assert sheet.row("b").fields == {} and sheet.row("b").row_sha256 == hashlib.sha256(b"{}").hexdigest()
    assert sheet.row("c") == S.SampleRow({}, os.path.abspath(tmp_path / "s.csv"), None)
    assert sheet.source == os.path.abspath(tmp_path / "s.csv") and sheet.key == "image_id"
    assert sheet.image_ids == ["a", "b"] and sheet.extra_columns == ["notes"]
    assert sheet.columns == ["image_id", "genotype", "sex", "treatment", "notes"]
    assert sheet.rows == [{"image_id": "a", "genotype": "N2", "sex": "male", "treatment": "", "notes": "x"},
                          {"image_id": "b", "genotype": "", "sex": "", "treatment": "", "notes": "only a note"}]
    # a non-ASCII value is escaped in the canonical JSON, so the hash does not depend on the file's encoding
    row = S.read_sample_sheet(_write_csv(tmp_path / "u.csv", [{"image_id": "a", "treatment": f"t{E_ACUTE}moin"}]))
    assert row.row("a").row_sha256 == hashlib.sha256(b'{"treatment":"t\\u00e9moin"}').hexdigest()


def test_excel_resaves_keep_every_row_hash(tmp_path):
    rows = [{"image_id": "a", "genotype": "N2", "sex": "male", "treatment": f"t{E_ACUTE}moin", "notes": "x"},
            {"image_id": "b", "genotype": "unc-119(ed3)", "sex": "herm", "replicate": "2", "notes": "y"}]
    cols = ["image_id", "genotype", "sex", "treatment", "replicate", "notes"]
    base = S.read_sample_sheet(_write_csv(tmp_path / "base.csv", rows, columns=cols))
    assert base.row("b").fields == {"genotype": "unc-119(ed3)", "sex": "herm", "replicate": "2"}
    forms = {
        "crlf": _write_csv(tmp_path / "crlf.csv", rows, columns=cols, newline="\r\n"),
        "bom": _write_csv(tmp_path / "bom.csv", rows, columns=cols, encoding="utf-8-sig"),
        "reordered": _write_csv(tmp_path / "re.csv", rows, columns=list(reversed(cols))),
        "padded_quoted": _write_csv(tmp_path / "pq.csv", rows, columns=cols, quote_all=True, pad="  "),
        "cp1252": _write_csv(tmp_path / "cp.csv", rows, columns=cols, encoding="cp1252", newline="\r\n"),
    }
    for name, path in forms.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", S.SampleSheetWarning)
            sheet = S.read_sample_sheet(path)
        for iid in ("a", "b"):
            assert sheet.row(iid).row_sha256 == base.row(iid).row_sha256, (name, iid)
            assert sheet.row(iid).fields == base.row(iid).fields, (name, iid)
    # one edited design cell changes that row's hash only; an edited extra cell changes none
    edited = [dict(r) for r in rows]
    edited[1]["replicate"] = "3"
    edited[0]["notes"] = "z"
    sheet = S.read_sample_sheet(_write_csv(tmp_path / "ed.csv", edited, columns=cols))
    assert sheet.row("a").row_sha256 == base.row("a").row_sha256
    assert sheet.row("b").row_sha256 != base.row("b").row_sha256


def test_cp1252_reads_with_a_warning(tmp_path):
    path = _write_csv(tmp_path / "s.csv", [{"image_id": "a", "genotype": f"caf{E_ACUTE}"}], encoding="cp1252")
    with pytest.warns(S.SampleSheetWarning, match="cp1252"):
        sheet = S.read_sample_sheet(path)
    assert sheet.encoding == "cp1252" and sheet.row("a").fields == {"genotype": f"caf{E_ACUTE}"}
    assert len(sheet.warnings) == 1 and "CSV UTF-8" in sheet.warnings[0]
    utf8 = _write_csv(tmp_path / "u.csv", [{"image_id": "a", "genotype": f"caf{E_ACUTE}"}])
    with warnings.catch_warnings():
        warnings.simplefilter("error")                       # UTF-8 reads without a warning
        assert S.read_sample_sheet(utf8).encoding == "utf-8-sig"


@pytest.mark.parametrize("value, iso", [("2025-11-05", True), ("20251105", True), ("11/5/2025", False),
                                        ("2025-13-05", False), ("5 Nov 2025", False), ("2025-11-5", False)])
def test_non_iso_date_is_kept_with_a_warning(tmp_path, value, iso):
    path = _write_csv(tmp_path / "s.csv", [{"image_id": "a", "acquisition_date": value}])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        sheet = S.read_sample_sheet(path)
    assert sheet.row("a").fields == {"acquisition_date": value}               # kept as given
    found = [str(w.message) for w in caught if issubclass(w.category, S.SampleSheetWarning)]
    assert found == sheet.warnings
    if iso:
        assert found == []
    else:
        assert len(found) == 1 and value in found[0] and "ISO" in found[0] and "row 2" in found[0]


def test_write_then_read_round_trips(tmp_path, monkeypatch):
    rows = [{"image_id": "a", "genotype": "N2", "sex": "male", "notes": "has, a comma"},
            {"image_id": "b", "genotype": "unc-119(ed3)", "sex": "", "notes": 'says "hi"'}]
    path = S.write_sample_sheet(tmp_path / "s.csv", rows)
    assert path == tmp_path / "s.csv"
    data = path.read_bytes()
    assert b"\r" not in data and data.decode("ascii").startswith("image_id,genotype,sex,notes\n")
    sheet = S.read_sample_sheet(path)
    assert sheet.rows == rows and sheet.columns == list(rows[0]) and sheet.warnings == []
    # numbers and None are written as text and an empty cell
    S.write_sample_sheet(tmp_path / "n.csv", [{"image_id": "c", "replicate": 3, "notes": None}])
    assert S.read_sample_sheet(tmp_path / "n.csv").rows == [{"image_id": "c", "replicate": "3", "notes": ""}]
    # the write is atomic: a failure leaves the old file and no temporary file
    def boom(*a, **k):
        raise OSError("disk full")

    monkeypatch.setattr(S, "_replace", boom)
    with pytest.raises(OSError, match="disk full"):
        S.write_sample_sheet(path, [{"image_id": "z"}])
    assert path.read_bytes() == data and sorted(p.name for p in tmp_path.iterdir()) == ["n.csv", "s.csv"]


@pytest.mark.parametrize("rows, match", [
    ([{"image_id": "a", "genotype": f"caf{E_ACUTE}"}], "row 2 (a), column genotype"),     # row 1 is the header
    ([{"image_id": "a", f"not{E_ACUTE}s": "x"}], "column"),
    ([{"genotype": "N2"}], "image_id or file"),
    ([{"image_id": "a"}, {"image_id": "a"}], "twice"),
    ([{"image_id": ""}], "no image_id"),
    ([{"image_id": "a", "replicate": object()}], "row 2 (a), column replicate"),
    ([], "no rows"),
])
def test_write_refuses(tmp_path, rows, match):
    with pytest.raises(S.SampleSheetError, match=match.replace("(", r"\(").replace(")", r"\)")):
        S.write_sample_sheet(tmp_path / "s.csv", rows)
    assert not (tmp_path / "s.csv").exists()


@pytest.mark.parametrize("text, match", [
    (b"", "empty"),
    (b"genotype,sex\nN2,herm\n", "image_id or file"),
    (b"image_id,genotype\na,N2\na,N3\n", "rows 2 and 3"),
    (b"image_id,genotype\n,N2\n", "row 2 has no image_id"),
    (b"image_id,Genotype\na,N2\n", "lower case"),
    (b"image_id,genotype,genotype\na,N2,N3\n", "twice"),
    (b"image_id,genotype\na,N2,extra\n", "more cells"),
    (b"image_id,,genotype\na,x,N2\n", "no name"),
    (b"image_id;genotype\na;N2\n", "semicolon"),
    ("image_id,genotype\na,N2\n".encode("utf-16"), "UTF-16"),
    (b"image_id,genotype\na,\x81\x8d\n", "cannot be read"),
    (b"image_id,file_path\na,x\n", "file_path"),
])
def test_read_refusals(tmp_path, text, match):
    path = tmp_path / "s.csv"
    path.write_bytes(text)
    with pytest.raises(S.SampleSheetError, match=match) as e:
        S.read_sample_sheet(path)
    assert str(path) in str(e.value) or os.path.abspath(path) in str(e.value)


def test_read_tolerates_what_spreadsheets_add(tmp_path):
    """Empty trailing rows and empty unnamed trailing columns, as Excel saves them, are ignored."""
    path = tmp_path / "s.csv"
    path.write_bytes(b"image_id,genotype,,\r\na,N2,,\r\n,,,\r\n,,,\r\n")
    sheet = S.read_sample_sheet(path)
    assert sheet.columns == ["image_id", "genotype"] and sheet.rows == [{"image_id": "a", "genotype": "N2"}]


def test_file_column_keys_rows(tmp_path):
    (tmp_path / "imgs").mkdir()
    path = _write_csv(tmp_path / "s.csv", [{"file": "x.ome.tif", "genotype": "N2"},
                                          {"file": "imgs/y.nd2", "genotype": "N3"}])
    sheet = S.read_sample_sheet(path)
    assert sheet.key == "file" and sheet.image_ids == ["x", "y"] and sheet.extra_columns == []
    assert sheet.row("x").fields == {"genotype": "N2"}                       # a bare name: its image id
    assert sheet.row("x", tmp_path / "anywhere" / "x.ome.tif").fields == {"genotype": "N2"}
    assert sheet.row("y", tmp_path / "imgs" / "y.nd2").fields == {"genotype": "N3"}   # a path: that file only
    assert sheet.row("y", tmp_path / "other" / "y.nd2").row_sha256 is None
    assert sheet.row("y").fields == {"genotype": "N3"}
    # with both columns the image_id is the key and file an extra column
    both = S.read_sample_sheet(_write_csv(tmp_path / "b.csv", [{"image_id": "x", "file": "x.ome.tif"}]))
    assert both.key == "image_id" and both.extra_columns == ["file"]
    with pytest.raises(S.SampleSheetError, match="rows 2 and 3"):
        S.read_sample_sheet(_write_csv(tmp_path / "d.csv", [{"file": "a/x.nd2"}, {"file": "b/x.ome.tif"}]))


@pytest.mark.parametrize("cells, warned", [
    ({"sex": "hermaphrodite"}, "sex 'hermaphrodite'"), ({"sex": "Hermaphrodite"}, "sex 'Hermaphrodite'"),
    ({"sex": "female"}, "sex 'female'"), ({"germ_cell": "Oocyte"}, "germ_cell 'Oocyte'"),
    ({"sex": "herm", "germ_cell": "oocytes"}, "germ_cell 'oocytes'"),
    ({"sex": "Herm"}, None), ({"sex": "M"}, None), ({"sex": "unknown"}, None), ({"sex": "NA"}, None),
    ({"germ_cell": "spermatocyte"}, None), ({"germ_cell": "unknown"}, None),
    ({"sex": "hermaphrodite", "germ_cell": "oocyte"}, None), ({"genotype": "N2"}, None),
])
def test_a_germ_cell_no_count_can_be_made_from_is_warned(tmp_path, monkeypatch, capsys, cells, warned):
    """Review S1: a sex that parse_sample's map does not know (on a row without germ_cell), or a germ_cell other than
    oocyte or spermatocyte, would give the image no expected SC or crossover count without a word. It is kept as
    given (the same fields and row hash), with a warning that names the row and the image; the CLI prints it."""
    path = _write_csv(tmp_path / "s.csv", [{"image_id": "a", **cells}])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        sheet = S.read_sample_sheet(path)
    assert sheet.row("a").fields == cells and sheet.row("a").row_sha256 == _canonical(cells)
    found = [str(w.message) for w in caught if issubclass(w.category, S.SampleSheetWarning)]
    assert found == sheet.warnings
    if warned is None:
        assert found == []
        return
    assert len(found) == 1 and warned in found[0], found
    assert "row 2 (image a)" in found[0] and "no expected count" in found[0]
    # the value stays as written: germ_cell is what parse_sample's rule makes of the sex, or the cell itself
    sample = {"sex": "herm", "germ_cell": "oocyte"}
    S.apply_row(sample, sheet.row("a"))
    assert sample["germ_cell"] == cells.get("germ_cell", "unknown")
    # `germquant run --samples` prints it, as it prints the date warning
    monkeypatch.setattr(api, "run", lambda image, **kw: {"image_id": "a", "n_nuclei": 0, "qc_pass": True,
                                                         "qc_flags": [], "out_dir": "x"})
    assert cli.main(["run", str(tmp_path / "a.nd2"), "--config", str(ROOT / "config" / "config.yaml"), "--out",
                     str(tmp_path / "o"), "--samples", str(path)]) == 0
    assert f"warning: {found[0]}" in capsys.readouterr().out


def test_a_file_path_with_backslashes_names_the_same_file_everywhere(tmp_path, monkeypatch):
    """Review S2: a `file` cell written with Windows separators names the same image and file on Linux and macOS,
    where pathlib.Path is a POSIX path (stood in for by PurePosixPath, so the test bites on every system). The cell
    is kept as written in SampleSheet.rows, and the key is not hashed."""
    from pathlib import PurePosixPath

    import germquant.io.formats as F

    cells = [{"file": "day1\\img_001.nd2", "genotype": "N2"}]
    sheet = S.read_sample_sheet(_write_csv(tmp_path / "s.csv", cells))
    assert sheet.image_ids == ["img_001"] and sheet.rows == cells
    assert sheet.row("img_001", tmp_path / "day1" / "img_001.nd2").fields == {"genotype": "N2"}
    assert sheet.row("img_001", tmp_path / "day2" / "img_001.nd2").row_sha256 is None
    assert sheet.row("img_001").row_sha256 == _canonical({"genotype": "N2"})
    monkeypatch.setattr(F, "Path", PurePosixPath)
    monkeypatch.setattr(S, "Path", PurePosixPath)
    text = "file,genotype\nday1\\img_001.nd2,N2\nday2/img_002.nd2,N3\n"
    posix = S._parse(text, "/study/samples.csv", PurePosixPath("/study"), "utf-8-sig")
    assert posix.image_ids == ["img_001", "img_002"]
    assert posix.row("img_001", PurePosixPath("/study/day1/img_001.nd2")).fields == {"genotype": "N2"}
    assert posix.row("img_002", PurePosixPath("/study/day2/img_002.nd2")).fields == {"genotype": "N3"}
    assert posix.rows[0]["file"] == "day1\\img_001.nd2"


def test_the_collect_hint_quotes_its_paths(tmp_path, monkeypatch, capsys):
    """Review N7: after a plain `batch --samples` whose sheet has columns of its own, the printed `germquant collect`
    line quotes both paths, so it can be pasted as one command when a folder name holds a space."""
    src = tmp_path / "in put"
    src.mkdir()
    (src / IN).write_bytes(b"")
    rows = []

    def fake(nd2, cfg, out_dir, *, xy_stride=1, z_range=None, prov=None, sample=None):
        rows.append(sample)
        return {"image_id": IN_ID, "n_nuclei": 1, "qc_pass": True, "qc_flags": [], "out_dir": str(out_dir),
                "tables": {"image_summary": pd.DataFrame([{"n_germline_nuclei": 1}])}, "stages": {}}

    monkeypatch.setattr(pipeline, "process_image", fake)
    sheet = _write_csv(tmp_path / "my sheet.csv", [{"image_id": IN_ID, "genotype": "N2", "slide": "S1"}])
    out = tmp_path / "out dir"
    assert cli.main(["batch", str(src), "--config", str(_study(tmp_path)), "--out", str(out), "--samples",
                     str(sheet)]) == 0
    assert [r.fields for r in rows] == [{"genotype": "N2"}]
    assert f'`germquant collect "{out}" --samples "{sheet}"`' in capsys.readouterr().out


def test_importing_the_module_is_light():
    heavy, r = heavy_after("import germquant.io.sample_sheet\n")
    assert r.returncode == 0, r.stderr[-2000:]
    assert heavy == "", heavy


# ---------------------------------------------------------------------------------------------------
# collect --samples
# ---------------------------------------------------------------------------------------------------

def test_collect_joins_extra_columns(golden_runs, tmp_path, caplog, capsys):
    out = tmp_path / "results"
    shutil.copytree(golden_runs["out"], out)
    plain = batch.collect(out)                                  # without the sheet: no extra column anywhere
    no_sheet = {t: list(pd.read_csv(out / f"batch_{t}.csv", nrows=0).columns) for t in ("nuclei", "image_summary")}
    assert all("slide" not in cols for cols in no_sheet.values())
    # the sheet's extra columns are appended to every stacked table; the design columns stay the run's own
    rows = [dict(GOLDEN_ROWS[0], genotype="edited after the run", n_nuclei="99"), dict(GOLDEN_ROWS[1])]
    sheet = _write_csv(tmp_path / "s.csv", rows)
    with caplog.at_level(logging.WARNING, logger="germquant.batch"):
        counts = batch.collect(out, sample_sheet=sheet)
    assert counts == plain
    nuclei = pd.read_csv(out / "batch_nuclei.csv", dtype=str, keep_default_na=False)
    assert list(nuclei.columns) == no_sheet["nuclei"] + ["slide", "n_nuclei"]
    assert set(nuclei.loc[nuclei["image_id"] == IN_ID, "slide"]) == {"S12"}
    assert set(nuclei.loc[nuclei["image_id"] == OUT_ID, "slide"]) == {""}
    assert set(nuclei["genotype"]) == {"syp-2(ok307)", "N2"}
    summ = pd.read_csv(out / "batch_image_summary.csv", dtype=str, keep_default_na=False)
    assert list(summ.columns) == no_sheet["image_summary"] + ["slide"]        # n_nuclei is a column of its own there
    assert summ.set_index("image_id").loc[IN_ID, "n_nuclei"] != "99"
    assert "the sample sheet's column n_nuclei is already a column of batch_image_summary.csv" in caplog.text
    assert list(pd.read_csv(out / "batch_summary.csv").columns) == batch.SUMMARY_COLS
    for t in schema.TABLES:
        cols = list(pd.read_csv(out / f"batch_{t}.csv", nrows=0).columns)
        assert "slide" in cols, t
    # `germquant collect ROOT --samples FILE` does the same; an unreadable sheet is refused and nothing is written
    before = (out / "batch_nuclei.csv").read_bytes()
    assert cli.main(["collect", str(out), "--samples", str(sheet)]) == 0
    assert (out / "batch_nuclei.csv").read_bytes() == before
    bad = tmp_path / "bad.csv"
    bad.write_bytes(b"genotype\nN2\n")
    (out / "batch_nuclei.csv").unlink()
    capsys.readouterr()
    with pytest.raises(SystemExit) as e:
        cli.main(["collect", str(out), "--samples", str(bad)])
    assert e.value.code == 2 and "REFUSED: " in capsys.readouterr().err and not (out / "batch_nuclei.csv").exists()
