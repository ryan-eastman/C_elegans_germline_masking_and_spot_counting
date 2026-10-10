"""The pre-registered RAD-51 benchmark split (tests/golden/rad51_split.json), its helpers in
germquant.validate.benchmark_split and the manifest builder scripts/benchmark_manifest.py.

The real-data checks skip when the N2 RAD-51 stacks are not on the machine (CI); everything else runs
on files made in tmp_path, including a tiny synthetic Imaris project written with h5py.
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from germquant.validate import benchmark_split as bs

ROOT = Path(__file__).resolve().parents[1]
SPLIT = ROOT / "tests" / "golden" / "rad51_split.json"
TEST_TAGS = {"n2rad51_noHS_herm_002", "n2rad51_noHS_male_002", "n2rad51_HS_herm_002", "n2rad51_HS_male_001"}
VAL_TAGS = {"n2rad51_noHS_herm_001", "n2rad51_HS_herm_001"}
REGION_GROWING_TAGS = {"n2rad51_HS_male_000", "n2rad51_HS_male_002"}


def _synthetic_split(entries):
    """A valid split record from (tag, role, nd2 path, ims path) tuples of small files."""
    gonads = []
    for tag, role, nd2, ims in entries:
        parts = bs.tag_parts(tag)
        gonads.append({"tag": tag, "role": role, "condition": parts["condition"], "sex": parts["sex"],
                       "replicate": parts["replicate"], "stem": Path(nd2).stem,
                       "nd2": bs.file_record(nd2), "ims": {**bs.file_record(ims), "region_growing": False},
                       "notes": []})
    split = {"name": "synthetic",
             "roles": {r: sorted(g["tag"] for g in gonads if g["role"] == r) for r in bs.ROLES},
             "gonads": gonads}
    bs.validate_split(split)
    return split


def test_split_file_loads_with_the_registered_roles():
    split = bs.load_split(SPLIT)
    tags = [g["tag"] for g in split["gonads"]]
    assert len(tags) == 12 and len(set(tags)) == 12
    assert {g["tag"]: g["role"] for g in split["gonads"]} == bs.PREREGISTERED_ROLES
    assert set(bs.split_tags(split, "test")) == TEST_TAGS
    assert set(bs.split_tags(split, "val")) == VAL_TAGS
    assert set(bs.split_tags(split, "train")) == set(tags) - TEST_TAGS - VAL_TAGS
    assert len(bs.split_tags(split, "train")) == 6
    assert split["roles"] == {r: sorted(bs.split_tags(split, r)) for r in bs.ROLES}
    assert split["registered_on"] == bs.REGISTERED_ON
    for g in split["gonads"]:
        parts = bs.tag_parts(g["tag"])
        assert (g["condition"], g["sex"]) == (parts["condition"], parts["sex"])
        assert Path(g["nd2"]["path"]).stem == g["stem"]
        assert Path(g["ims"]["path"]).stem == g["stem"]
        for key in ("nd2", "ims"):
            assert re.fullmatch(r"[0-9a-f]{64}", g[key]["sha256"]), (g["tag"], key)
            assert g[key]["size_bytes"] > 1_000_000_000, (g["tag"], key)   # whole stacks and projects
    # the region-growing exports are the two HS male gonads named in the protocol, nothing else
    rg = {g["tag"]: g["ims"]["region_growing"] for g in split["gonads"]}
    assert {t for t, v in rg.items() if v is True} == REGION_GROWING_TAGS
    assert all(v is False for t, v in rg.items() if t not in REGION_GROWING_TAGS)
    assert all(bs.REGION_GROWING_NOTE in g["notes"] for g in split["gonads"] if g["tag"] in REGION_GROWING_TAGS)
    assert bs.NOTES["n2rad51_noHS_male_002"][0] in bs.gonad(split, "n2rad51_noHS_male_002")["notes"]
    # the Imaris facts the protocol quotes
    assert sum(g["ims"]["n_surfaces"] for g in split["gonads"]) == 1198
    assert all(g["ims"]["imaris_voxel_um"] == [0.108, 0.108, 0.2] for g in split["gonads"])


def test_split_paths_exist_with_the_registered_sizes():
    split = bs.load_split(SPLIT)
    if not Path(split["gonads"][0]["nd2"]["path"]).parent.is_dir():
        pytest.skip("the N2 RAD-51 stacks are not on this machine")
    assert bs.check_split_files(split) == []


def test_sha256_helper_round_trips(tmp_path):
    data = bytes(range(256)) * 37 + b"germquant"
    p = tmp_path / "small.bin"
    p.write_bytes(data)
    expect = hashlib.sha256(data).hexdigest()
    assert bs.file_sha256(p) == expect
    assert bs.file_sha256(p, chunk_bytes=100) == expect              # chunking does not change the digest
    assert bs.file_record(p) == {"path": p.as_posix(), "size_bytes": len(data), "sha256": expect}
    split = _synthetic_split([("n2rad51_noHS_herm_000", "train", p, p)])
    assert bs.check_split_files(split, hashes=True) == []
    p.write_bytes(data[:-1] + b"X")                                   # same size, different content
    assert bs.check_split_files(split) == []                          # a size check cannot see it
    problems = bs.check_split_files(split, hashes=True)
    assert len(problems) == 2 and all("sha256" in s for s in problems)
    p.unlink()
    assert all("missing" in s for s in bs.check_split_files(split))


def test_hash_cache_is_reused_only_at_the_same_size(tmp_path):
    p = tmp_path / "a.nd2"
    p.write_bytes(b"abc")
    cache_file = tmp_path / "cache.json"
    cache_file.write_text(json.dumps({p.as_posix(): {"size_bytes": 3, "sha256": "f" * 64}}), encoding="utf-8")
    cache = bs.load_hash_cache(cache_file)
    assert bs._record_with_cache(p, cache)["sha256"] == "f" * 64
    p.write_bytes(b"abcd")
    assert bs._record_with_cache(p, cache)["sha256"] == hashlib.sha256(b"abcd").hexdigest()
    # a split file works as a cache too
    split = _synthetic_split([("n2rad51_HS_herm_001", "val", p, p)])
    split_file = tmp_path / "split.json"
    split_file.write_text(json.dumps(split), encoding="utf-8")
    assert bs.load_hash_cache(split_file) == {p.as_posix(): {"size_bytes": 4,
                                                             "sha256": hashlib.sha256(b"abcd").hexdigest()}}


def test_validate_split_rejects_malformed_records(tmp_path):
    p = tmp_path / "a.nd2"
    p.write_bytes(b"x")
    split = _synthetic_split([("n2rad51_noHS_herm_000", "train", p, p), ("n2rad51_HS_male_001", "test", p, p)])
    bad = json.loads(json.dumps(split))
    bad["gonads"][1]["tag"] = "n2rad51_noHS_herm_000"
    with pytest.raises(ValueError, match="unique"):
        bs.validate_split(bad)
    bad = json.loads(json.dumps(split))
    bad["gonads"][0]["role"] = "holdout"
    with pytest.raises(ValueError, match="role"):
        bs.validate_split(bad)
    bad = json.loads(json.dumps(split))
    bad["gonads"][0]["nd2"]["sha256"] = "abc"
    with pytest.raises(ValueError, match="sha256"):
        bs.validate_split(bad)
    bad = json.loads(json.dumps(split))
    bad["gonads"][0]["sex"] = "male"
    with pytest.raises(ValueError, match="condition/sex"):
        bs.validate_split(bad)
    bad = json.loads(json.dumps(split))
    bad["roles"]["test"] = []
    with pytest.raises(ValueError, match="roles"):
        bs.validate_split(bad)
    with pytest.raises(ValueError, match="not an n2rad51 tag"):
        bs.tag_parts("n2dryice_noHS_male_001")


def test_read_ims_summary_on_a_synthetic_project(tmp_path):
    h5py = pytest.importorskip("h5py")

    def attr(s):   # Imaris stores attributes as arrays of single bytes
        return np.frombuffer(s.encode(), dtype="S1")

    def project(path, flag):
        cp = (f'<bpPointsCreationParameters mEnableRegionOfInterest="false" mEnableRegionGrowing="{flag}" '
              'mEstimatedDiameter="0.400000 0.400000 0.900000">\n<mSpotFilter>\n'
              '<bpStatisticsFilter mLowerThresholdEnable="true" mLowerThresholdManualValue="42.5" mSelectHigh="true">\n'
              '  <bpStatisticsValueType mName="Quality" mUnit="" mFactors="0"/>\n</bpStatisticsFilter>\n'
              '</mSpotFilter>\n<mTrackFilter>\n<bpStatisticsFilter mLowerThresholdManualValue="2.5">\n'
              '  <bpStatisticsValueType mName="Track Duration Steps" mUnit=""/>\n</bpStatisticsFilter>\n'
              '</mTrackFilter>\n</bpPointsCreationParameters>\n')
        with h5py.File(path, "w") as f:
            img = f.create_group("DataSetInfo/Image")
            for k, v in {"X": "20", "Y": "10", "Z": "5", "ExtMin0": "100", "ExtMax0": "102.16",
                         "ExtMin1": "0", "ExtMax1": "1.08", "ExtMin2": "0", "ExtMax2": "1"}.items():
                img.attrs[k] = attr(v)
            f.create_dataset("Scene8/Content/Points0/CreationParameters",
                             data=np.array([cp.encode("latin1")], dtype=f"S{len(cp)}"))
            spot = np.zeros(7, dtype=[("ID", "<i8"), ("PositionX", "<f4"), ("PositionY", "<f4"),
                                      ("PositionZ", "<f4"), ("Radius", "<f4")])
            f.create_dataset("Scene8/Content/Points0/Spot", data=spot)
            smi = np.zeros(3, dtype=[("Volume", "<f4"), ("CenterOfMassX", "<f8")])
            f.create_dataset("Scene8/Content/MegaSurfaces0/SurfaceModelInfo", data=smi)

    project(tmp_path / "grow.ims", "true")
    s = bs.read_ims_summary(tmp_path / "grow.ims")
    assert s == {"n_spots": 7, "n_surfaces": 3, "region_growing": True, "spot_quality_threshold": 42.5,
                 "shape_xyz": [20, 10, 5], "imaris_voxel_um": [0.108, 0.108, 0.2]}
    project(tmp_path / "plain.ims", "false")
    assert bs.read_ims_summary(tmp_path / "plain.ims")["region_growing"] is False
    # a project without spots or creation parameters is unknown, not an error
    with h5py.File(tmp_path / "empty.ims", "w") as f:
        f.create_group("DataSetInfo")
    s = bs.read_ims_summary(tmp_path / "empty.ims")
    assert s["region_growing"] is None and s["n_spots"] is None and s["spot_quality_threshold"] is None


def test_build_split_assembles_roles_and_notes_without_the_data(tmp_path):
    """build_split on twelve stub files: roles from PREREGISTERED_ROLES, stems from the .nd2 names,
    the Imaris path derived per condition, notes attached; metadata reading skipped."""
    ims_root = tmp_path / "imaris"
    runs = []
    for tag, role in bs.PREREGISTERED_ROLES.items():
        parts = bs.tag_parts(tag)
        stem = f"20251105_N2_{parts['condition'].lower()}_{parts['sex'].upper()} _{parts['replicate']}"
        nd2 = tmp_path / parts["condition"] / f"{stem}.nd2"
        nd2.parent.mkdir(exist_ok=True)
        nd2.write_bytes(tag.encode())
        ims = ims_root / bs.IMS_SUBDIR[parts["condition"]] / f"{stem}.ims"
        ims.parent.mkdir(parents=True, exist_ok=True)
        ims.write_bytes(b"ims " + tag.encode())
        runs.append({"tag": tag, "nd2": nd2.as_posix(), "config": "config/config.yaml", "xy_stride": 1})
    runs.append({"tag": "n2dryice_noHS_male_001", "nd2": "x.nd2"})   # other goldens are ignored
    manifest = tmp_path / "real_manifest.json"
    manifest.write_text(json.dumps(runs), encoding="utf-8")
    split = bs.build_split(manifest, ims_root, read_metadata=False)
    assert {g["tag"]: g["role"] for g in split["gonads"]} == bs.PREREGISTERED_ROLES
    assert split["roles"]["test"] == sorted(TEST_TAGS) and split["roles"]["val"] == sorted(VAL_TAGS)
    g = bs.gonad(split, "n2rad51_HS_male_001")
    assert g["stem"] == "20251105_N2_hs_MALE _001"
    assert g["ims"]["path"] == (ims_root / "HS_ims2" / "20251105_N2_hs_MALE _001.ims").as_posix()
    assert g["ims"]["sha256"] == hashlib.sha256(b"ims n2rad51_HS_male_001").hexdigest()
    assert g["ims"]["region_growing"] is None and bs.REGION_GROWING_UNKNOWN_NOTE in g["notes"]
    assert bs.NOTES["n2rad51_noHS_male_002"][0] in bs.gonad(split, "n2rad51_noHS_male_002")["notes"]
    # a manifest missing a registered tag is refused
    manifest.write_text(json.dumps(runs[1:]), encoding="utf-8")
    with pytest.raises(ValueError, match="n2rad51_noHS_herm_000"):
        bs.build_split(manifest, ims_root, read_metadata=False)


def test_manifest_builder_on_a_temp_tree(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    stem_a, stem_b = "20251105_N2_nohs_HERM _", "20251105_n2_hs_MALE _001"
    files = {}
    for stem, ext in ((stem_a, "nd2"), (stem_a, "ims"), (stem_b, "nd2"), (stem_b, "ims")):
        files[(stem, ext)] = data / f"{stem}.{ext}"
        files[(stem, ext)].write_bytes(f"{stem}.{ext}".encode())
    split = _synthetic_split([("n2rad51_noHS_herm_000", "train", files[(stem_a, "nd2")], files[(stem_a, "ims")]),
                              ("n2rad51_HS_male_001", "test", files[(stem_b, "nd2")], files[(stem_b, "ims")])])
    split["gonads"][1]["ims"]["region_growing"] = True
    split["gonads"][1]["notes"] = [bs.REGION_GROWING_NOTE]
    split_path = tmp_path / "split.json"
    split_path.write_text(json.dumps(split), encoding="utf-8")
    results = tmp_path / "results"
    (results / "n2rad51_noHS_herm_000").mkdir(parents=True)
    for suffix in ("__nuclei.csv", "__nuclei_labels.tif", "__spots.csv"):
        (results / "n2rad51_noHS_herm_000" / f"{stem_a}{suffix}").write_text("x", encoding="utf-8")
    truth = tmp_path / "truth"
    truth.mkdir()
    truth_names = sorted(f"{stem_b}{s}" for s in ("__imaris_spots.csv", "__imaris_nuclei.csv",
                                                  "__imaris_labels.tif", "__imaris_extra.json"))
    for name in truth_names:
        (truth / name).write_text("y", encoding="utf-8")
    (truth / f"{stem_a}_imaris_spots.csv").write_text("z", encoding="utf-8")   # one underscore: not ours

    m = bs.build_manifest(split, results, truth, name="synthetic", split_path=split_path)
    assert [g["tag"] for g in m["gonads"]] == ["n2rad51_noHS_herm_000", "n2rad51_HS_male_001"]
    a, b = m["gonads"]
    assert a["outputs_complete"] is True and b["outputs_complete"] is False
    assert a["outputs"]["nuclei_labels"]["exists"] is True and a["outputs"]["image_summary"]["exists"] is False
    assert a["outputs"]["nuclei"]["path"] == (results / "n2rad51_noHS_herm_000" / f"{stem_a}__nuclei.csv").as_posix()
    assert a["truth_complete"] is False and a["truth"]["found"] == []
    assert b["truth_complete"] is True and b["truth"]["found"] == truth_names
    assert b["truth"]["params"]["exists"] is False
    assert a["spot_truth_usable"] is True and b["spot_truth_usable"] is False
    assert b["notes"] == [bs.REGION_GROWING_NOTE]
    assert a["nd2_record"]["exists"] and a["nd2_record"]["size_matches"] and a["ims_record"]["region_growing"] is False
    assert a["nd2"] == a["nd2_record"]["path"] and a["split"] == "train" and a["truth_dir"] == truth.as_posix()
    assert m["summary"] == {"n_gonads": 2, "outputs_complete": ["n2rad51_noHS_herm_000"],
                            "truth_complete": ["n2rad51_HS_male_001"], "spot_truth_usable": ["n2rad51_noHS_herm_000"]}
    assert m["roles"] == {"train": ["n2rad51_noHS_herm_000"], "val": [], "test": ["n2rad51_HS_male_001"]}
    assert m["split_file"]["sha256"] == bs.file_sha256(split_path)
    # a missing truth root or results root is recorded, not an error
    empty = bs.build_manifest(split, tmp_path / "nowhere", tmp_path / "nowhere_either")
    assert empty["summary"]["outputs_complete"] == [] and empty["summary"]["truth_complete"] == []
    # an overridden truth file name is looked up instead of the default
    alt = bs.build_manifest(split, results, truth, truth_files={"meta": "{stem}__imaris_extra.json"})
    assert alt["gonads"][1]["truth"]["meta"]["exists"] is True and alt["truth_files"]["spots"] == bs.TRUTH_FILES["spots"]

    out = bs.write_manifest(m, tmp_path / "benchmarks")
    assert out == tmp_path / "benchmarks" / "synthetic" / "manifest.json"
    assert json.loads(out.read_text(encoding="utf-8"))["name"] == "synthetic"

    # the command-line builder writes the same manifest (apart from the timestamp)
    rc = subprocess.run([sys.executable, str(ROOT / "scripts" / "benchmark_manifest.py"),
                         "--split", str(split_path), "--results-root", str(results), "--truth-root", str(truth),
                         "--out", str(tmp_path / "cli"), "--name", "synthetic"],
                        capture_output=True, text=True, check=False)
    assert rc.returncode == 0, rc.stderr
    cli = json.loads((tmp_path / "cli" / "synthetic" / "manifest.json").read_text(encoding="utf-8"))
    cli.pop("created_utc")
    m.pop("created_utc")
    assert cli == m
    assert "outputs complete for 1" in rc.stdout and "spot truth usable for 1" in rc.stdout
