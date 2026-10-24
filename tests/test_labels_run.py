"""Labels from files through process_image and `germquant batch` (docs/DESIGN_system_packs.md 4.11 and 3.10 item 4;
decisions D12 to D15 of the Train B plan; docs/LABELS.md).

* A detector-only run of the synthetic gonad of tests/test_end_to_end.py with its planted truth labels is the
  synthetic_labels golden (tests/golden/synthetic_labels/, its config_hash sentinel pinned in
  tests/golden/default_invariants.json), and it equals, table for table, a run whose segmentation produced the
  same labels (acceptance 5).
* segmentation_method reads labels_file and the flag segment:labels_from_file is raised only in labels runs;
  run_geometry holds labels_from {path, sha256} and the sha256 enters the segment sub-hash (acceptance 3).
* `batch --resume` computes the segment hash per image, so editing one labels file reprocesses that image only; a
  missing labels file fails that image with a clear message (acceptance 3, 4).
* Labels runs resolve, hash and pin-check no model (model_sha256 "unused"), and cli._batch follows the same rule
  (D12); require_channel does not apply; in a projected run the projected measure and records stay and a labels
  file with more than one plane is refused (D14); labels with t are refused (D15).

Regenerate the golden only for a deliberate change of this path: ``GERMQUANT_UPDATE_GOLDEN=1 pytest
tests/test_labels_run.py``, as for tests/test_golden.py.
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import planar_helpers as H
import pytest
from labels_helpers import (
    TEMPLATE,
    assert_tables_equal,
    disks,
    patch_reader,
    read_table,
    truth_labels,
    write_ome,
    write_tiff,
)
from test_end_to_end import _synthetic_stack
from test_projected_run import CELL_R, CELLS, DEBRIS, DEBRIS_R, IMAGE, IMAGE_ID, NX, NY, NZ, SP, Z_FOCUS, recipe
from test_projected_run import yeast_like_stack

from germquant import api, batch, cli, models_registry, pipeline, provenance, schema, stages
from germquant.config import load_config
from germquant.io import labels as L
from germquant.io.nd2_reader import Stack
from germquant.segment import backends as SB

ROOT = Path(__file__).resolve().parents[1]
GOLD = ROOT / "tests" / "golden" / "synthetic_labels"
GOLD_TABLES = ("nuclei", "spots", "image_summary")
ND2 = "20251105_n2_nohs_HERM_001.nd2"
ND2_ID = "20251105_n2_nohs_HERM_001"
FLAG = "segment:labels_from_file"
PLANE = f"{H.STEM}.nd2"            # the calibrated single plane of tests/planar_helpers.py, through its reader seam


def _cfg(base: Path = ROOT / "config" / "config.yaml", **overrides):
    cfg = load_config(base)
    for key, value in overrides.items():
        cfg.set(key, value)
    return cfg


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def gonad_runs(tmp_path_factory):
    """The synthetic gonad run twice with every default stage (spots included): detector-only through
    labels_from (the golden run), and segmented by a registered backend that returns the same truth labels."""
    work = tmp_path_factory.mktemp("labels")
    truth = truth_labels()
    labels_file = write_tiff(work / f"{ND2_ID}_truth_labels.tif", truth)
    image = work / ND2                     # never read from disk: the synthetic reader stands in
    with pytest.MonkeyPatch.context() as mp:
        patch_reader(mp, _synthetic_stack())
        cfg_l = _cfg(**{"segmentation.nuclei.labels_from": TEMPLATE})
        res_l = pipeline.process_image(image, cfg_l, work / "labels", prov=None)
        mp.setattr(SB, "_REGISTRY", dict(SB._REGISTRY))
        mp.setattr(SB, "_SOURCES", dict(SB._SOURCES))
        mp.setitem(SB._STATE, "discovered", True)

        def truth_backend(dna, spacing, *, model, diameter_um, device, cfg_section):
            assert dna.shape == truth.shape
            return truth.copy()

        SB.register_backend("truth", truth_backend)
        cfg_s = _cfg(**{"segmentation.nuclei.method": "cellpose", "segmentation.nuclei.backend": "truth"})
        res_s = pipeline.process_image(image, cfg_s, work / "segmented", prov=None)
    return {"truth": truth, "labels_file": labels_file, "cfg_l": cfg_l, "res_l": res_l, "out_l": work / "labels",
            "cfg_s": cfg_s, "res_s": res_s, "out_s": work / "segmented"}


def test_synthetic_labels_golden(gonad_runs):
    out, cfg = gonad_runs["out_l"], gonad_runs["cfg_l"]
    files = {t: out / f"{ND2_ID}__{t}.csv" for t in GOLD_TABLES}
    if os.environ.get("GERMQUANT_UPDATE_GOLDEN"):
        GOLD.mkdir(parents=True, exist_ok=True)
        for t, f in files.items():
            read_table(f).to_csv(GOLD / f"{t}.csv", index=False)
        (GOLD / "config_hash.txt").write_text(cfg.hash + "\n", encoding="utf-8")
        pytest.skip("synthetic_labels golden regenerated")
    sentinel = (GOLD / "config_hash.txt").read_text(encoding="utf-8").strip()
    assert cfg.hash == sentinel, ("config/config.yaml, its channel map or the golden's labels_from template "
                                  "changed: the golden's config_hash no longer reproduces")
    for t, f in files.items():
        ref = pd.read_csv(GOLD / f"{t}.csv", low_memory=False, float_precision="round_trip")
        assert_tables_equal(read_table(f), ref, t)
    # the golden is a real detector-only run: the planted truth, one RAD-51 focus in each germline nucleus
    nuclei = pd.read_csv(GOLD / "nuclei.csv")
    assert len(nuclei) == int(gonad_runs["truth"].max()) == 15
    assert int(nuclei["in_germline"].sum()) == 12 and (nuclei.loc[nuclei["in_germline"], "n_spots"] >= 1).all()


def test_method_and_flag_only_in_labels_runs(gonad_runs):
    res_l, res_s = gonad_runs["res_l"], gonad_runs["res_s"]
    assert res_l["tables"]["image_summary"]["segmentation_method"].iloc[0] == "labels_file"
    assert FLAG in res_l["qc_flags"] and FLAG in res_l["stages"]["segment"]["flags"]
    assert res_s["tables"]["image_summary"]["segmentation_method"].iloc[0] == "truth"
    assert FLAG not in res_s["qc_flags"] and not any(f.startswith("segment:") for f in res_s["qc_flags"])
    assert res_l["stages"]["segment"]["status"] == "ran"
    # require_channel's fallback flag is never raised in a labels run: nothing is segmented
    assert not any(f.startswith("segment:no_") for f in res_l["qc_flags"])


def test_labels_sha_in_geometry_and_segment_hash(gonad_runs):
    out, cfg, path = gonad_runs["out_l"], gonad_runs["cfg_l"], gonad_runs["labels_file"]
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest, record = _json(out / "run_manifest.json"), _json(out / f"{ND2_ID}__stages.json")
    done = _json(out / f"{ND2_ID}__done.json")
    for rec in (manifest, record):
        assert rec["labels_from"] == {"path": str(path), "sha256": sha}
        assert rec["model_sha256"] == "unused" and "segment_backend" not in rec
    tools = provenance.tool_versions()
    with_sha = stages.stage_hashes(cfg, {}, tools, "unused", labels_sha256=sha)
    without = stages.stage_hashes(cfg, {}, tools, "unused")
    assert manifest["stage_hashes"]["segment"] == record["stage_hashes"]["segment"] == done["stage_hashes"]["segment"]
    assert manifest["stage_hashes"]["segment"] == with_sha["segment"] != without["segment"]
    assert {k: v for k, v in with_sha.items() if k != "segment"} == {k: v for k, v in without.items() if k != "segment"}
    # no new marker key: the labels sha256 reaches resume through the segment hash
    assert sorted(done) == ["config_hash", "enabled_stages", "failed_stages", "finished", "image_id",
                            "stage_hashes", "xy_stride", "z_range"]
    # the segmented run records none of it
    seg = _json(gonad_runs["out_s"] / "run_manifest.json")
    assert "labels_from" not in seg and seg["model_sha256"] == "absent"


def test_detector_only_equals_segmented_run(gonad_runs):
    """Acceptance 5: the truth labels read from a file give the tables of a run whose segmentation produced them."""
    import tifffile

    out_l, out_s = gonad_runs["out_l"], gonad_runs["out_s"]
    lab_l = tifffile.imread(out_l / f"{ND2_ID}__nuclei_labels.tif")
    np.testing.assert_array_equal(lab_l, tifffile.imread(out_s / f"{ND2_ID}__nuclei_labels.tif"))
    np.testing.assert_array_equal(lab_l, gonad_runs["truth"])
    for t in ("nuclei", "spots", "granules", "coloc"):
        assert_tables_equal(read_table(out_l / f"{ND2_ID}__{t}.csv"), read_table(out_s / f"{ND2_ID}__{t}.csv"), t,
                            append_only=False)
    assert len(pd.read_csv(out_l / f"{ND2_ID}__spots.csv")) == 12
    isum_l, isum_s = (read_table(o / f"{ND2_ID}__image_summary.csv") for o in (out_l, out_s))
    assert isum_l["segmentation_method"].iloc[0] == "labels_file" and isum_s["segmentation_method"].iloc[0] == "truth"
    flags_l, flags_s = (set(str(df["qc_flags"].iloc[0]).split(";")) for df in (isum_l, isum_s))
    assert flags_l == flags_s | {FLAG}
    keep = [c for c in isum_s.columns if c not in ("segmentation_method", "qc_flags")]
    assert_tables_equal(isum_l[keep], isum_s[keep], "image_summary", append_only=False)


def _study(tmp_path, names, **extra) -> tuple[Path, Path]:
    """A folder of .nd2 stand-ins (the synthetic reader reads them) and a study profile outside the repository with
    labels_from set, a pinned pack model that is not installed, and the slow stages off."""
    src = tmp_path / "in"
    src.mkdir()
    for n in names:
        (src / n).write_bytes(b"")
    lines = [f"extends: {(ROOT / 'config' / 'profiles' / 'generic.yaml').as_posix()}",
             "segmentation:", "  nuclei:", "    method: cellpose", "    cellpose_model: cpsam_v2@0f1cc3f7",
             f"    labels_from: '{TEMPLATE}'", "spots:", "  enabled: false", "render:", "  montage: false"]
    lines += [f"{k}: {v}" for k, v in extra.items()]
    study = tmp_path / "study.yaml"
    study.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return src, study


def _no_model(monkeypatch):
    """D12: a labels run neither resolves nor pin-checks a model (both would raise for this study's model)."""
    def resolve(cfg):
        raise AssertionError("resolve_model_path called in a labels run")

    def pin(*a, **k):
        raise AssertionError("check_model_pin called in a labels run")

    monkeypatch.setattr(pipeline, "resolve_model_path", resolve)
    monkeypatch.setattr(models_registry, "check_model_pin", pin)


def test_batch_resume_reprocesses_only_the_edited_labels_image(tmp_path, monkeypatch, capsys):
    names = ["20251105_n2_nohs_HERM_001.nd2", "20251105_n2_nohs_HERM_002.nd2", "20251105_n2_nohs_HERM_003.nd2"]
    src, study = _study(tmp_path, names)
    truth = truth_labels()
    for n in names[:2]:                                    # the third image has no labels file
        write_tiff(src / f"{Path(n).stem}_truth_labels.tif", truth)
    patch_reader(monkeypatch, _synthetic_stack())
    _no_model(monkeypatch)
    calls: list[str] = []
    real = pipeline.process_image

    def spy(nd2, *a, **k):
        calls.append(Path(nd2).name)
        return real(nd2, *a, **k)

    monkeypatch.setattr(pipeline, "process_image", spy)
    decisions: list[dict] = []
    real_done = batch.is_done
    monkeypatch.setattr(batch, "is_done", lambda *a, **k: decisions.append(k) or real_done(*a, **k))
    out = tmp_path / "out"
    assert cli.main(["batch", str(src), "--config", str(study), "--out", str(out)]) == 0
    assert calls == names
    summ = pd.read_csv(out / "batch_summary.csv").set_index("image_id")
    missing = src / "20251105_n2_nohs_HERM_003_truth_labels.tif"
    assert summ.loc["20251105_n2_nohs_HERM_003", "qc_flags"].startswith("EXCEPTION:labels file not found")
    assert str(missing) in summ.loc["20251105_n2_nohs_HERM_003", "qc_flags"]
    assert summ.loc["20251105_n2_nohs_HERM_001", "n_nuclei"] == 15
    # cli._batch's model rule is process_image's (D12): "unused" in the batch manifest, and its per-image segment
    # hash (germquant.io.labels.segment_hashes, D46) is the one process_image wrote into the marker; a recipe
    # without io.planar has no planar hash (D44)
    assert _json(out / "run_manifest.json")["model_sha256"] == "unused"
    first = Path(names[0]).stem
    marker = _json(out / first / f"{first}__done.json")
    assert L.segment_hashes(load_config(study), src / names[0]) == (marker["stage_hashes"]["segment"], None, None)
    # --resume: the per-image segment hashes cli._batch computes equal process_image's, so the two are skipped
    calls.clear()
    capsys.readouterr()
    assert cli.main(["batch", str(src), "--config", str(study), "--out", str(out), "--resume"]) == 0
    printed = capsys.readouterr().out
    assert calls == [names[2]] and printed.count("already done") == 2
    assert "labels file not found" in printed                # the missing file is named before the retry
    # today's is_done keywords (D40, D44): no segment_hash_planar for a recipe that does not set io.planar; each
    # image its own hash, and the unreadable one a value no marker holds
    assert [sorted(d) for d in decisions] == [["segment_hash", "xy_stride"]] * 3
    assert decisions[0]["segment_hash"] == marker["stage_hashes"]["segment"] != decisions[2]["segment_hash"]
    assert decisions[2]["segment_hash"] == L.UNREADABLE_HASH
    # editing one labels file reprocesses that image only
    edited = truth.copy()
    edited[edited == 3] = 0
    write_tiff(src / "20251105_n2_nohs_HERM_002_truth_labels.tif", edited)
    calls.clear()
    capsys.readouterr()
    assert cli.main(["batch", str(src), "--config", str(study), "--out", str(out), "--resume"]) == 0
    printed = capsys.readouterr().out
    assert calls == names[1:] and printed.count("already done") == 1
    assert "segment stage hash differs" in printed
    summ = pd.read_csv(out / "batch_summary.csv").set_index("image_id")
    assert summ.loc["20251105_n2_nohs_HERM_002", "n_nuclei"] == 14 and summ.loc["20251105_n2_nohs_HERM_001", "n_nuclei"] == 15


def test_segment_hashes(tmp_path):
    """germquant.io.labels.segment_hashes (D46): one image's segment sub-hash in a labels run, its labels file's
    sha256 in the payload; the single-plane one too only when the recipe sets io.planar (D44); "labels-unreadable",
    a value no marker holds, with the reason, for a file that cannot be found or hashed; refused for a recipe that is
    not a labels run, whose segment hash is one per run."""
    image = tmp_path / ND2
    labels_file = write_tiff(tmp_path / f"{ND2_ID}_truth_labels.tif", np.ones((4, 4), np.int32))
    sha = hashlib.sha256(labels_file.read_bytes()).hexdigest()
    tools = provenance.tool_versions()
    cfg = _cfg(**{"segmentation.nuclei.labels_from": TEMPLATE})
    assert L.segment_hashes(cfg, image) == (stages.stage_hashes(cfg, {}, tools, "unused", labels_sha256=sha)["segment"],
                                            None, None)
    assert L.segment_hashes(cfg, image, model_sha="m")[0] == \
        stages.stage_hashes(cfg, {}, tools, "m", labels_sha256=sha)["segment"]
    for value in ("auto", True, False):
        cfg_p = _cfg(**{"segmentation.nuclei.labels_from": TEMPLATE, "io.planar": value})
        flat = stages.stage_hashes(cfg_p, {}, tools, "unused", labels_sha256=sha)["segment"]
        plane = stages.stage_hashes(cfg_p, {}, tools, "unused", planar=True, labels_sha256=sha)["segment"]
        assert flat != plane
        assert L.segment_hashes(cfg_p, image) == (flat, None if value is False else plane, None), value
    # a file that cannot be found or hashed
    other = tmp_path / "20251105_n2_nohs_HERM_009.nd2"
    seg, seg_planar, problem = L.segment_hashes(cfg, other)
    assert (seg, seg_planar) == (L.UNREADABLE_HASH, None) and problem.startswith("labels file not found")
    assert str(tmp_path / "20251105_n2_nohs_HERM_009_truth_labels.tif") in problem
    cfg_p = _cfg(**{"segmentation.nuclei.labels_from": TEMPLATE, "io.planar": "auto"})
    assert L.segment_hashes(cfg_p, other)[:2] == (L.UNREADABLE_HASH, L.UNREADABLE_HASH)
    seg, _, problem = L.segment_hashes(_cfg(**{"segmentation.nuclei.labels_from": "{image_dir}/{image_stem}.npy"}),
                                       image)
    assert seg == L.UNREADABLE_HASH and "import-labels" in problem
    with pytest.raises(ValueError, match="labels_from"):
        L.segment_hashes(_cfg(), image)


def test_missing_labels_file_fails_that_image(tmp_path, monkeypatch):
    patch_reader(monkeypatch, _synthetic_stack())
    out = tmp_path / "r"
    out.mkdir()
    (out / f"{ND2_ID}__done.json").write_text("{}", encoding="utf-8")         # a stale marker of an older run
    cfg = _cfg(ROOT / "config" / "profiles" / "generic.yaml", **{"segmentation.nuclei.labels_from": TEMPLATE})
    with pytest.raises(FileNotFoundError) as err:
        pipeline.process_image(tmp_path / ND2, cfg, out, prov=None)
    msg = str(err.value)
    assert "labels file not found" in msg and str(tmp_path / f"{ND2_ID}_truth_labels.tif") in msg
    assert TEMPLATE in msg and ND2_ID in msg
    assert not list(out.iterdir())                       # refused before any stage: nothing written, marker cleared
    # a _seg.npy, and a file that does not fit the image, are refused the same way, before any stage
    cfg.set("segmentation.nuclei.labels_from", "{image_dir}/{image_stem}_seg.npy")
    with pytest.raises(L.LabelsError, match="import-labels"):
        pipeline.process_image(tmp_path / ND2, cfg, out, prov=None)
    write_tiff(tmp_path / f"{ND2_ID}_small.tif", np.ones((80, 160), np.int32))
    cfg.set("segmentation.nuclei.labels_from", "{image_dir}/{image_stem}_small.tif")
    with pytest.raises(L.LabelsError, match=r"\(80, 160\).*\(20, 80, 320\)"):
        pipeline.process_image(tmp_path / ND2, cfg, out, prov=None)
    (tmp_path / f"{ND2_ID}_folder.tif").mkdir()                              # a folder under a labels name
    cfg.set("segmentation.nuclei.labels_from", "{image_dir}/{image_stem}_folder.tif")
    with pytest.raises(L.LabelsError, match="labels file cannot be read"):
        pipeline.process_image(tmp_path / ND2, cfg, out, prov=None)
    assert not list(out.iterdir())


def test_no_model_is_resolved_in_labels_runs(tmp_path, monkeypatch):
    """D12: model_sha256 "unused"; no resolution, no pin check, no segmentation backend record."""
    patch_reader(monkeypatch, _synthetic_stack())
    _no_model(monkeypatch)
    write_tiff(tmp_path / f"{ND2_ID}_truth_labels.tif", truth_labels())
    cfg = _cfg(ROOT / "config" / "profiles" / "generic.yaml",
               **{"segmentation.nuclei.labels_from": TEMPLATE, "segmentation.nuclei.method": "cellpose",
                  "segmentation.nuclei.cellpose_model": "cpsam_v2@0f1cc3f7",
                  "segmentation.nuclei.backend": "nuclei3d", "spots.enabled": False, "render.montage": False})
    res = pipeline.process_image(tmp_path / ND2, cfg, tmp_path / "r", prov=None)
    manifest = _json(tmp_path / "r" / "run_manifest.json")
    assert manifest["model_sha256"] == "unused" and "segment_backend" not in manifest
    assert res["n_nuclei"] == 15 and res["tables"]["image_summary"]["segmentation_method"].iloc[0] == "labels_file"
    assert pipeline.segment_backend_provenance(cfg) == {}


def test_labels_frame_stride_and_z_range_reach_the_run(tmp_path, monkeypatch):
    """The run passes labels_frame (D13), xy_stride and z_range to the labels (acceptance 2): frame 1 of a
    (T, Z, Y, X) segm.npz, read at xy_stride 2 and z_range (6, 15), is the truth cut as the stack is cut."""
    import tifffile

    full = _synthetic_stack()

    def read_cut(path, *, xy_stride=1, z_range=None, **kw):          # read_stack's slicing (nd2_reader)
        data = full.data if z_range is None else full.data[:, z_range[0]:z_range[1]]
        dz, dy, dx = full.spacing
        return dataclasses.replace(full, data=np.ascontiguousarray(data[:, :, ::xy_stride, ::xy_stride]),
                                   spacing=(dz, dy * xy_stride, dx * xy_stride), downsample_xy=xy_stride)

    patch_reader(monkeypatch, full)                                  # the header holds the full 20 planes
    monkeypatch.setattr(pipeline, "read_stack", read_cut)
    truth = truth_labels()
    np.savez_compressed(str(tmp_path / f"{ND2_ID}_segm.npz"), np.stack([np.zeros_like(truth), truth]))
    cfg = _cfg(ROOT / "config" / "profiles" / "generic.yaml",
               **{"segmentation.nuclei.labels_from": "{image_dir}/{image_stem}_segm.npz",
                  "segmentation.nuclei.labels_frame": 1, "spots.enabled": False, "render.montage": False})
    out = tmp_path / "cut"
    res = pipeline.process_image(tmp_path / ND2, cfg, out, prov=None, xy_stride=2, z_range=(6, 15))
    assert res["n_nuclei"] == 15                                     # frame 0 is empty: frame 1 was read
    np.testing.assert_array_equal(tifffile.imread(out / f"{ND2_ID}__nuclei_labels.tif"), truth[6:15, ::2, ::2])


def test_labels_in_a_projected_run(tmp_path, monkeypatch):
    """D14: the labels file replaces segment_planar's output in a projected run; the projected measure and records
    stay, the projection is computed only for the montage, and a file with more than one plane is refused."""
    import tifffile

    seen = []
    real = pipeline.make_montage

    def spy(*args, **kw):
        seen.append(kw)
        return real(*args, **kw)

    monkeypatch.setattr(pipeline, "make_montage", spy)
    stack, _ = yeast_like_stack()
    patch_reader(monkeypatch, stack)
    # the eight cells and the two specks of debris, which a segmenting projected run's area gate drops
    truth2d = disks((NY, NX), CELLS, CELL_R)
    debris = disks((NY, NX), DEBRIS, DEBRIS_R, first_id=len(CELLS) + 1)
    truth2d = np.where(debris > 0, debris, truth2d)
    template = "{image_dir}/{image_stem}_masks.tif"
    write_tiff(tmp_path / f"{IMAGE_ID}_masks.tif", truth2d)
    cfg = recipe(**{"segmentation.nuclei.labels_from": template, "spots.enabled": False})
    out = tmp_path / "p"
    res = pipeline.process_image(tmp_path / IMAGE, cfg, out, prov=None)
    labels = tifffile.imread(out / f"{IMAGE_ID}__nuclei_labels.tif")
    assert labels.shape == (NZ, NY, NX) and all(np.array_equal(labels[z], truth2d) for z in range(NZ))
    assert res["n_nuclei"] == len(CELLS) + len(DEBRIS)            # no size gate: the debris below min_area_um2 stays
    nuclei = pd.read_csv(out / f"{IMAGE_ID}__nuclei.csv")
    cols = list(nuclei.columns)
    assert cols[len(schema.NUCLEI):len(schema.NUCLEI) + 2] == ["area_um2", "equivalent_diameter_um"]
    assert nuclei["volume_um3"].isna().all()
    counts = np.bincount(truth2d.ravel())
    np.testing.assert_allclose(nuclei.set_index("nucleus_id")["area_um2"].sort_index(), counts[1:] * SP[1] * SP[2])
    np.testing.assert_allclose(nuclei["equivalent_diameter_um"], 2 * np.sqrt(nuclei["area_um2"] / np.pi))
    assert "segment:projected_plane" in res["qc_flags"] and FLAG in res["qc_flags"]
    assert {"segment:projected_plane", FLAG} <= set(res["stages"]["segment"]["flags"])
    manifest = _json(out / "run_manifest.json")
    assert manifest["segment_projection"] == {"mode": "plane", "plane": "sharpest", "z": Z_FOCUS}
    assert manifest["labels_from"]["path"] == str(tmp_path / f"{IMAGE_ID}_masks.tif")
    assert res["tables"]["image_summary"]["segmentation_method"].iloc[0] == "labels_file"
    np.testing.assert_array_equal(seen[-1]["base_plane"], stack.data[0][Z_FOCUS])   # the montage's base plane
    # MDN1 Position_4's case: no transmitted channel. require_channel does not apply, nothing is projected
    no_trans = Stack(data=stack.data[::-1].copy(), spacing=SP, channel_names=["Cy3", "Cy5"], path="x.nd2",
                     spacing_ok=True)
    patch_reader(monkeypatch, no_trans)
    res2 = pipeline.process_image(tmp_path / IMAGE, cfg, tmp_path / "nt", prov=None)
    assert res2["n_nuclei"] == len(CELLS) + len(DEBRIS) and not any(f.startswith("segment:no_") for f in res2["qc_flags"])
    assert seen[-1]["base_plane"] is None
    assert _json(tmp_path / "nt" / "run_manifest.json")["segment_projection"]["z"] is None
    # with require_channel false the record resolves its plane on channel 0, as a segmenting run does, and the
    # montage shows that very plane
    cfg_c0 = recipe(**{"segmentation.nuclei.labels_from": template, "spots.enabled": False,
                       "segmentation.nuclei.require_channel": False})
    pipeline.process_image(tmp_path / IMAGE, cfg_c0, tmp_path / "c0", prov=None)
    z_c0 = _json(tmp_path / "c0" / "run_manifest.json")["segment_projection"]["z"]
    assert z_c0 is not None
    np.testing.assert_array_equal(seen[-1]["base_plane"], no_trans.data[0][z_c0])
    # a segmentation backend without planar_ok is not refused: nothing is segmented
    patch_reader(monkeypatch, stack)
    cfg3 = recipe(**{"segmentation.nuclei.labels_from": template, "spots.enabled": False,
                     "segmentation.nuclei.method": "cellpose", "segmentation.nuclei.backend": "nuclei3d"})
    assert pipeline.process_image(tmp_path / IMAGE, cfg3, tmp_path / "n3", prov=None)["n_nuclei"] == 10
    # a labels file with more than one plane is refused in a projected run, before any stage
    write_ome(tmp_path / f"{IMAGE_ID}_masks3d.ome.tif", np.repeat(truth2d[None], NZ, axis=0), "ZYX")
    cfg4 = recipe(**{"segmentation.nuclei.labels_from": "{image_dir}/{image_stem}_masks3d.ome.tif",
                     "spots.enabled": False})
    out4 = tmp_path / "refused"
    with pytest.raises(L.LabelsError, match="more than one plane"):
        pipeline.process_image(tmp_path / IMAGE, cfg4, out4, prov=None)
    assert out4.is_dir() and not list(out4.iterdir())


def test_labels_plus_t_is_refused(tmp_path, monkeypatch, capsys):
    """D15: a labels_from time series is refused in cli._run_series before the loop, in api.run and in
    process_image."""
    patch_reader(monkeypatch, _synthetic_stack())
    write_tiff(tmp_path / f"{ND2_ID}_truth_labels.tif", truth_labels())
    cfg = _cfg(ROOT / "config" / "profiles" / "generic.yaml", **{"segmentation.nuclei.labels_from": TEMPLATE})
    with pytest.raises(L.LabelsError, match="time series"):
        pipeline.process_image(tmp_path / ND2, cfg, tmp_path / "t", prov=None, t=0)
    with pytest.raises(L.LabelsError, match="time series"):
        api.run(tmp_path / ND2, config=cfg, out=tmp_path / "api", t=0)
    from germquant import timelapse

    def no_series(*a, **k):
        raise AssertionError("process_series called")

    monkeypatch.setattr(timelapse, "process_series", no_series)
    study = tmp_path / "s.yaml"
    study.write_text(f"extends: {(ROOT / 'config' / 'profiles' / 'generic.yaml').as_posix()}\n"
                     f"segmentation:\n  nuclei:\n    labels_from: '{TEMPLATE}'\n", encoding="utf-8")
    capsys.readouterr()
    rc = cli.main(["run", str(tmp_path / ND2), "--config", str(study), "--out", str(tmp_path / "cli"), "--t", "0-1"])
    err = capsys.readouterr().err
    assert rc == 2 and f"REFUSED: {L.SERIES_REFUSAL}" in err
    assert not (tmp_path / "cli").exists()
    # the same image without --t runs
    res = api.run(tmp_path / ND2, config=_cfg(ROOT / "config" / "profiles" / "generic.yaml",
                                              **{"segmentation.nuclei.labels_from": TEMPLATE, "spots.enabled": False,
                                                 "render.montage": False}), out=tmp_path / "fixed")
    assert res["n_nuclei"] == 15


def test_one_series_refusal_for_labels_and_single_planes(tmp_path, monkeypatch, capsys):
    """D45: cli._series_refusal holds every refusal of `run --t`, the labels one first because it needs no header
    read, and _run_series prints each as "REFUSED: <message>" with exit code 2 before any timepoint (the planar one
    of P2b too, which printed "ERROR:" before the sync)."""
    from germquant import planar, timelapse

    def no_header(*a, **k):
        raise AssertionError("the header was read")

    with monkeypatch.context() as mp:
        mp.setattr(pipeline, "read_input_metadata", no_header)
        cfg = H.generic_2d(**{"segmentation.nuclei.labels_from": TEMPLATE})       # io.planar auto, and labels
        assert cli._series_refusal(cfg, tmp_path / PLANE) == L.SERIES_REFUSAL
    H.patch_reader(monkeypatch, series=2)

    def no_series(*a, **k):
        raise AssertionError("process_series called")

    monkeypatch.setattr(timelapse, "process_series", no_series)
    study = tmp_path / "series_2d.yaml"
    study.write_text(f"extends: {(ROOT / 'config' / 'profiles' / 'generic_2d.yaml').as_posix()}\n"
                     "segmentation:\n  nuclei:\n    method: classical\nspots:\n  enabled: false\n", encoding="utf-8")
    assert cli._series_refusal(load_config(study), PLANE) == planar.PLANAR_SERIES_REFUSED
    capsys.readouterr()
    rc = cli.main(["run", PLANE, "--config", str(study), "--out", str(tmp_path / "s"), "--t", "0-1"])
    err = capsys.readouterr().err
    assert rc == 2 and f"REFUSED: {planar.PLANAR_SERIES_REFUSED}" in err and "ERROR" not in err
    assert not (tmp_path / "s").exists()
    # a labels recipe that is also single-plane gets the labels message (the first refusal), in the same form
    study.write_text(study.read_text(encoding="utf-8").replace(
        "method: classical\n", f"method: classical\n    labels_from: '{TEMPLATE}'\n"), encoding="utf-8")
    rc = cli.main(["run", PLANE, "--config", str(study), "--out", str(tmp_path / "s"), "--t", "0-1"])
    err = capsys.readouterr().err
    assert rc == 2 and f"REFUSED: {L.SERIES_REFUSAL}" in err and planar.PLANAR_SERIES_REFUSED not in err


# ---- labels on single-plane images (io.planar, decision D14; after the sync on merge #1) ----------------------

def _plane_labels() -> np.ndarray:
    """int32 (Y, X) truth labels of the calibrated plane of tests/planar_helpers.py: its nine planted nuclei (the
    planted disks, ids 1 to 9) and its two specks of debris (disks of 2 pixels, ids 10 and 11), which generic_2d's
    area gate (min_area_um2 20.0) drops in a segmenting run."""
    info = H.planted()
    yy, xx = np.indices((H.SIZE, H.SIZE)).astype(np.float64)
    lab = np.zeros((H.SIZE, H.SIZE), np.int32)
    for i, n in enumerate(info["nuclei"], 1):
        lab[(yy - n["y"]) ** 2 + (xx - n["x"]) ** 2 <= n["r_px"] ** 2] = i
    for i, d in enumerate(info["debris"], len(info["nuclei"]) + 1):
        lab[(yy - d["y"]) ** 2 + (xx - d["x"]) ** 2 <= 4.0] = i
    return lab


def test_batch_resume_planar_labels_folder(tmp_path, monkeypatch, capsys):
    """PLAN 10.4's done item: `batch --resume` on a single-plane labels study computes both per-image segment hashes
    (germquant.io.labels.segment_hashes) and passes the single-plane one to is_done because the recipe sets
    io.planar (D44, D46). The finished images are skipped, and an edited labels file reprocesses that image only."""
    src, masks = tmp_path / "in", tmp_path / "masks"
    src.mkdir()
    masks.mkdir()
    names = ["a_001.nd2", "b_002.nd2"]
    truth = _plane_labels()
    for n in names:
        (src / n).write_bytes(b"")
        write_tiff(masks / f"{Path(n).stem}_masks.tif", truth)
    # an absolute template into a folder of its own, so the images' folder holds images only
    template = f"{masks.as_posix()}/{{image_stem}}_masks.tif"
    study = tmp_path / "study_2d.yaml"
    study.write_text(f"extends: {(ROOT / 'config' / 'profiles' / 'generic_2d.yaml').as_posix()}\n"
                     f"segmentation:\n  nuclei:\n    method: classical\n    labels_from: '{template}'\n"
                     "spots:\n  enabled: false\nrender:\n  montage: false\n", encoding="utf-8")
    H.patch_reader(monkeypatch)
    calls: list[str] = []
    decisions: list[dict] = []
    real_process, real_done = pipeline.process_image, batch.is_done

    def process(path, *a, **k):
        calls.append(Path(path).name)
        return real_process(path, *a, **k)

    monkeypatch.setattr(pipeline, "process_image", process)
    monkeypatch.setattr(batch, "is_done", lambda *a, **k: decisions.append(k) or real_done(*a, **k))
    out = tmp_path / "out"
    assert cli.main(["batch", str(src), "--config", str(study), "--out", str(out)]) == 0
    assert calls == names
    cfg = load_config(study)
    pairs = []
    for n in names:
        iid = Path(n).stem
        done = _json(out / iid / f"{iid}__done.json")
        seg, seg_planar, problem = L.segment_hashes(cfg, src / n)
        assert done["planar"] is True and problem is None and done["stage_hashes"]["segment"] == seg_planar != seg
        pairs.append((seg, seg_planar))
    # --resume: both per-image hashes reach is_done, so the two finished single-plane images are skipped
    calls.clear()
    capsys.readouterr()
    assert cli.main(["batch", str(src), "--config", str(study), "--out", str(out), "--resume"]) == 0
    assert calls == [] and capsys.readouterr().out.count("already done") == 2
    assert all(sorted(d) == ["segment_hash", "segment_hash_planar", "xy_stride"] for d in decisions)
    assert [(d["segment_hash"], d["segment_hash_planar"]) for d in decisions] == pairs
    # editing one labels file reprocesses that image only
    edited = truth.copy()
    edited[edited == 3] = 0
    write_tiff(masks / "b_002_masks.tif", edited)
    assert cli.main(["batch", str(src), "--config", str(study), "--out", str(out), "--resume"]) == 0
    printed = capsys.readouterr().out
    assert calls == ["b_002.nd2"] and printed.count("already done") == 1 and "segment stage hash differs" in printed
    summ = pd.read_csv(out / "batch_summary.csv").set_index("image_id")
    assert summ.loc["b_002", "n_nuclei"] == 10 and summ.loc["a_001", "n_nuclei"] == 11
