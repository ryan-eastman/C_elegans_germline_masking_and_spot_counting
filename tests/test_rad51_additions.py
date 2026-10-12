"""The additions of protocol changelog entry 9 (tests/golden/rad51_additions.json, germquant.validate.benchmark_split):
the registered split file is byte-identical and its roles unmoved; the record validates with roles fixed in code
(train only, disjoint from the 12 registered gonads); the refusals (a registered tag, val or test, a duplicate, a
bad sha256, a stem that is not the .nd2, a sealed 2025-11-05 path, a record made against other roles); the
manifest with additions keeps the registered entries identical and decides the additions' spot truth from the
record and the truth params; the role guard; the creation-parameter summary and comparison; the builder on stub
files. Synthetic, apart from the checks of the real files, which skip off the workstation."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest
from test_benchmark_split import _synthetic_split

from germquant.benchmark import detector as bd
from germquant.validate import benchmark_split as bs

ROOT = Path(__file__).resolve().parents[1]
SPLIT = ROOT / "tests" / "golden" / "rad51_split.json"
ADDITIONS = ROOT / "tests" / "golden" / "rad51_additions.json"
ROLES_SHA256 = "32e2ee0b4a0dc193fe941b995449447d11004a377289e5ce8e3eda5af5ddec00"
MALES = ["n2rad51_20260105a_HS_male_001", "n2rad51_20260105a_HS_male_002", "n2rad51_20260105a_HS_male_003"]
HS1 = "n2rad51_20260105b_HS_male_001"


def test_registered_split_is_unchanged_and_its_roles_pinned():
    data = SPLIT.read_bytes()
    assert (hashlib.sha256(data).hexdigest() == bs.REGISTERED_SPLIT_SHA256
            or bs.lf_sha256(data) == bs.REGISTERED_SPLIT_SHA256_LF)
    assert bs.split_file_unchanged(SPLIT)
    assert bs.roles_sha256() == ROLES_SHA256                     # PREREGISTERED_ROLES never moves
    assert {g["tag"]: g["role"] for g in bs.load_split(SPLIT)["gonads"]} == bs.PREREGISTERED_ROLES


def test_the_additions_record_validates_with_fixed_train_roles():
    rec = bs.load_additions(ADDITIONS)
    tags = [g["tag"] for g in rec["gonads"]]
    assert tags == MALES + [HS1] and set(tags) == set(bs.ADDITION_ROLES)
    assert set(bs.ADDITION_ROLES.values()) == {"train"} and not set(bs.ADDITION_ROLES) & set(bs.PREREGISTERED_ROLES)
    assert rec["roles"] == {"train": sorted(tags), "val": [], "test": []}
    assert rec["changelog_entry"] == 9 and rec["decided_by"]["reply"] == "2 whatever you recommend"
    reg = rec["registered_split"]
    assert reg["roles_sha256"] == ROLES_SHA256 and reg["unchanged"] is True
    assert reg["sha256"] == bs.REGISTERED_SPLIT_SHA256 or reg["sha256_lf"] == bs.REGISTERED_SPLIT_SHA256_LF
    for g in rec["gonads"]:
        assert g["role"] == "train" and g["session"] == "2026-01-05" and (g["condition"], g["sex"]) == ("HS", "male")
        assert g["nd2"]["sizes"] == {"Z": 106, "C": 3, "Y": 2720, "X": 2720}
        assert g["nd2"]["channel_names"] == ["405", "477", "545"] and g["nd2"]["voxel_um"][2] == 0.2
        assert g["ims"]["region_growing"] is False and g["ims"]["masked_channel_present"] is True
        assert g["ims"]["spot_source_channel"] == {"index": 2, "name": "545"}
        assert all(g[k]["sha256_verified_against_copy_manifest"] for k in ("nd2", "ims", "xlsx"))
        assert not any(bs.is_sealed_path(g[k][f]) for k in ("nd2", "ims", "xlsx") for f in ("path", "nas_path"))
    by = {g["tag"]: g for g in rec["gonads"]}
    # the Imaris facts entry 9 (b) quotes
    assert [(by[t]["ims"]["n_spots"], by[t]["ims"]["n_surfaces"]) for t in MALES] == [(5187, 98), (12028, 88), (6462, 127)]
    assert [by[t]["ims"]["spot_quality_threshold"] for t in MALES] == [18.9129, 7.96068, 12.0885]
    assert all(by[t]["spot_truth"] == {"source": "ims", "status": "usable", "rule": "entry 9 (d)"} for t in MALES)
    assert by[HS1]["spot_truth"]["source"] == "xlsx" and by[HS1]["spot_truth"]["status"] == "pending"
    assert (by[HS1]["xlsx"]["overall_n_spots"], by[HS1]["xlsx"]["overall_n_surfaces"]) == (6512, 112)
    assert by[HS1]["ims"]["n_spots"] == 1226263
    # the per-axis Imaris voxel of the six-significant-digit extents (entry 9 (d)(3))
    assert by[MALES[0]]["ims"]["ext_min_xyz_text"][0] == "16249.9" and by[MALES[0]]["ims"]["imaris_voxel_um"][1] == 0.108


def test_the_recorded_files_are_on_disk():
    rec = bs.load_additions(ADDITIONS)
    if not Path(rec["gonads"][0]["nd2"]["path"]).parent.is_dir():
        pytest.skip("the 2026-01-05 copies are not on this machine")
    assert bs.check_addition_files(rec) == []


def _good_record() -> dict:
    return json.loads(ADDITIONS.read_text(encoding="utf-8"))


@pytest.mark.parametrize("mutate,match", [
    (lambda r: r["gonads"][0].update(tag="n2rad51_HS_herm_000"), "registered gonad"),
    (lambda r: r["gonads"][0].update(role="val"), "never val or test"),
    (lambda r: r["gonads"][0].update(role="test"), "never val or test"),
    (lambda r: r["gonads"][1].update(tag=r["gonads"][0]["tag"]), "unique"),
    (lambda r: r["gonads"][0]["ims"].update(sha256="abc"), "sha256"),
    (lambda r: r["gonads"][0].update(stem="something_else"), "stem"),
    (lambda r: r["gonads"][0]["xlsx"].update(nas_path="//nas/x/20251105_n2_hs_MALE_1.xlsx"), "sealed"),
    (lambda r: r["registered_split"].update(roles_sha256="0" * 64), "other registered roles"),
    (lambda r: r["gonads"][0].update(tag="n2rad51_20260105a_HS_male_009"), "not an addition of entry 9"),
    (lambda r: r["gonads"][0].update(tag="n2rad51_20260105_HS_male_001"), "not an addition tag"),
    (lambda r: r["gonads"][3].pop("xlsx"), "needs an xlsx record"),
    (lambda r: r["gonads"][0].update(sex="herm"), "condition/sex"),
    (lambda r: r["roles"].update(val=[r["gonads"][0]["tag"]]), "roles"),
])
def test_validate_additions_refuses(mutate, match):
    rec = _good_record()
    mutate(rec)
    with pytest.raises(ValueError, match=match):
        bs.validate_additions(rec)


def test_the_role_guard_knows_the_additions():
    for t in bs.ADDITION_ROLES:
        assert bd.preregistered_role(t) == "train" and bd.gonad_role(t, "train") == "train"
        assert bd.gonad_role(t, "") == "train" and not bd.is_test_tag(t) and bd.is_addition(t)
        for wrong in ("val", "test"):
            with pytest.raises(ValueError, match="pre-registered split says 'train'"):
                bd.gonad_role(t, wrong)
    bd.guard_test(list(bs.ADDITION_ROLES))                     # never a test gonad
    with pytest.raises(bd.TestSetGuardError):                  # the registered guard is unchanged
        bd.guard_test(["n2rad51_HS_male_001"])
    assert not bd.is_addition("n2rad51_HS_herm_000") and bd.gonad_role("n2rad51_HS_herm_001") == "val"


def test_sealed_paths():
    assert bs.is_sealed_path("//nas/TS/202511_imaris_export_het_mutants_Crest/20251105_n2_hs_MALE_1.xlsx")
    assert bs.is_sealed_path("C:/x/Data_while processing/20260116/x_temp_compiledtdata.xlsx")
    assert bs.is_sealed_path("C:/x/imaris_stats_processed_sections/compiled_dataEP_flat.xlsx")
    assert bs.is_sealed_path("C:/x/master_data.prism") and bs.is_sealed_path("C:/x/20260123_RAD51_stats.xlsx")
    assert not bs.is_sealed_path("C:/Users/ryane/n2_rad51_hs_2026/20260105_N2_HS/20260105_n2_hs_MALE_1.nd2")
    with pytest.raises(PermissionError, match="sealed"):
        bs.refuse_sealed("C:/data/20251105_n2_hs_MALE _001.ims")


# ---------------------------------------------------------------- the builder and the manifest on stubs
def _stub_additions(tmp_path: Path, *, sha_override: dict | None = None, extra_spec=None) -> tuple[dict, Path]:
    """build_additions on stub files with a stub copy manifest (no metadata read)."""
    data = tmp_path / "copies"
    data.mkdir(parents=True, exist_ok=True)
    files, specs = [], []
    for tag, stem, source in ((MALES[0], "20260105_n2_hs_MALE_1", "ims"), (HS1, "20260105_n2_hs_1_MALE_1", "xlsx")):
        paths = {}
        for key, ext in (("nd2", "nd2"), ("ims", "ims"), ("xlsx", "xlsx")):
            p = data / f"{stem}.{ext}"
            p.write_bytes(f"{tag} {ext}".encode())
            paths[key] = p
            sha = hashlib.sha256(p.read_bytes()).hexdigest()
            files.append({"source": f"//nas/N2/{p.name}", "dest": p.as_posix(), "size": p.stat().st_size,
                          "source_mtime": "2026-01-05T12:00:00", "sha256": (sha_override or {}).get(p.name, sha)})
        specs.append({"tag": tag, "session": "2026-01-05", "set": "N2/x", "stem": stem, **paths, "xlsx_use": "x",
                      "spot_truth": {"source": source, "status": "usable" if source == "ims" else "pending"}})
    if extra_spec:
        specs.append(extra_spec)
    cm = tmp_path / "copy_manifest.json"
    cm.write_text(json.dumps({"files": files}), encoding="utf-8")
    rec = bs.build_additions(specs, cm, registered_split_path=SPLIT, decided_by={"reply": "test"},
                             read_metadata=False)
    for g in rec["gonads"]:
        g["ims"]["region_growing"] = False          # what the .ims would say (no metadata is read here)
    return rec, cm


def test_build_additions_on_stub_files(tmp_path):
    rec, _cm = _stub_additions(tmp_path)
    assert [g["tag"] for g in rec["gonads"]] == [MALES[0], HS1]
    assert rec["roles"]["train"] == sorted([MALES[0], HS1]) and rec["registered_split"]["unchanged"] is True
    g = rec["gonads"][0]
    assert g["nd2"]["nas_path"] == "//nas/N2/20260105_n2_hs_MALE_1.nd2" and g["nd2"]["sha256_verified_against_copy_manifest"]
    # a file that is not the copied one, a sealed path, a tag outside entry 9: refused
    with pytest.raises(ValueError, match="copy manifest"):
        _stub_additions(tmp_path / "a", sha_override={"20260105_n2_hs_MALE_1.ims": "f" * 64})
    with pytest.raises(PermissionError, match="sealed"):
        _stub_additions(tmp_path / "b", extra_spec={"tag": MALES[1], "session": "2026-01-05", "set": "x",
                                                    "stem": "20251105_n2_hs_MALE _001",
                                                    "nd2": tmp_path / "20251105_n2_hs_MALE _001.nd2",
                                                    "ims": tmp_path / "y.ims", "spot_truth": {}})


def test_manifest_with_additions_keeps_the_registered_entries(tmp_path):
    rec, _cm = _stub_additions(tmp_path / "a")
    data = tmp_path / "data"
    data.mkdir()
    entries = []
    for tag, role in bs.PREREGISTERED_ROLES.items():
        p = data / f"{tag}.nd2"
        p.write_bytes(tag.encode())
        q = data / f"{tag}.ims"
        q.write_bytes(tag.encode())
        entries.append((tag, role, p, q))
    split = _synthetic_split(entries)
    results, truth = tmp_path / "results", tmp_path / "truth"
    plain = bs.build_manifest(split, results, truth)
    m = bs.build_manifest(split, results, truth, additions=rec, additions_results_root=tmp_path / "add_results",
                          additions_path=ADDITIONS)
    assert m["gonads"][:12] == plain["gonads"]                            # the registered entries, unchanged
    add = {e["tag"]: e for e in m["gonads"][12:]}
    assert list(add) == [MALES[0], HS1] and all(e["addition"] and e["role"] == "train" for e in add.values())
    assert add[MALES[0]]["results_dir"] == (tmp_path / "add_results" / MALES[0]).as_posix()
    # no truth export exists here: neither addition is usable
    assert add[MALES[0]]["spot_truth_usable"] is False and add[HS1]["spot_truth_usable"] is False
    assert add[MALES[0]]["spot_truth_reason"] == add[HS1]["spot_truth_reason"] == "no truth export yet"
    assert add[MALES[0]]["section_4_4"] is None and add[MALES[0]]["spot_truth_decision"] is None
    assert m["roles"]["train"][-2:] == [MALES[0], HS1] and m["roles"]["test"] == plain["roles"]["test"]
    assert m["additions_file"]["sha256"] == bs.file_sha256(ADDITIONS) and "additions_file" not in plain
    assert m["spot_truth_inputs"]["truth_check"] is None and m["spot_truth_inputs"]["decisions"] is None
    # a revised manifest goes next to the registered one under its own name
    out = bs.write_manifest(m, tmp_path / "benchmarks", "manifest_entry9.json")
    assert out.name == "manifest_entry9.json" and out.parent.name == m["name"]
    with pytest.raises(ValueError, match="file name"):
        bs.write_manifest(m, tmp_path / "benchmarks", "sub/manifest.json")
    # a split with other roles is refused
    moved = json.loads(json.dumps(split))
    moved["gonads"][0]["role"], moved["gonads"][1]["role"] = moved["gonads"][1]["role"], moved["gonads"][0]["role"]
    moved["roles"] = {r: sorted(g["tag"] for g in moved["gonads"] if g["role"] == r) for r in bs.ROLES}
    with pytest.raises(ValueError, match="other roles"):
        bs.build_manifest(moved, results, truth, additions=rec)
    # the section 4.4 table and the decisions belong to the additions
    with pytest.raises(ValueError, match="additions only"):
        bs.build_manifest(split, results, truth, truth_check=tmp_path / "x.csv")


# the frame check as the four real truth params record it (entry 9 (d)(3))
FC_ROUNDED = {"frame_check": {"chosen": {"name": "rounded", "ext_min_xyz": [0, 0, 0], "voxel_xyz_um": [0.108] * 3}}}


def _xlsx_params(**fails) -> dict:
    ver = {k: {"passes": k not in fails} for k in ("V1_sha256", "V2_totals", "V3_surfaces", "V4_inside_extent",
                                                   "V5_intensity_center", "V6_radius")}
    ver["V7_truth_assignment"] = {"passes": "V7_truth_assignment" not in fails,
                                  "agreement": fails.get("V7_truth_assignment", 0.996)}
    ver["passed"] = not fails
    return {**FC_ROUNDED, "spots_source": {"source": "xlsx", "verification": ver}}


def test_the_spot_truth_rule_of_entries_9_and_10(tmp_path):
    """addition_spot_truth_usable: the export's own checks first (no decision overrides them), then the changelog
    decision, else section 4.4 at least 0.9935."""
    rec, _cm = _stub_additions(tmp_path)
    g_ims, g_x = rec["gonads"]
    dec_keep = {"usable": True, "entry": 10, "reason": "kept post hoc"}
    dec_drop = {"usable": False, "entry": 10, "reason": "4.4 failed"}
    ok, low = {"agreement": 0.9960}, {"agreement": 0.9807}
    # no truth export, or one made without --frame-check, or one whose frame check chose nothing: never usable
    assert bs.addition_spot_truth_usable(g_ims, None, ok, dec_keep) == (False, "no truth export yet")
    for params in ({}, {"frame_check": {"chosen": None}}, {"frame_check": {"untestable": "no statistic"}}):
        usable, why = bs.addition_spot_truth_usable(g_ims, params, ok, dec_keep)
        assert usable is False and "frame check of entry 9 (d)(3) missing or failed" in why
    # region growing on, or unknown: never usable
    for rg in (True, None):
        grown = {**g_ims, "ims": {**g_ims["ims"], "region_growing": rg}}
        assert bs.addition_spot_truth_usable(grown, FC_ROUNDED, ok, dec_keep)[0] is False
    assert bs.addition_spot_truth_usable({**g_ims, "spot_truth": {"source": "ims", "status": "pending"}},
                                         FC_ROUNDED, ok)[0] is False
    # section 4.4: a male at 0.9807 without a decision is not usable; not run is not usable; 0.9960 is
    usable, why = bs.addition_spot_truth_usable(g_ims, FC_ROUNDED, low)
    assert usable is False and "needs a changelog entry" in why and "0.9807" in why
    assert bs.addition_spot_truth_usable(g_ims, FC_ROUNDED)[0] is False
    assert bs.addition_spot_truth_usable(g_ims, FC_ROUNDED, {"agreement": None})[0] is False
    assert bs.addition_spot_truth_usable(g_ims, FC_ROUNDED, ok)[0] is True
    assert bs.addition_spot_truth_usable(g_ims, FC_ROUNDED, {"agreement": bs.SECTION_4_4_MIN})[0] is True
    # the same male with an entry-10 decision follows the decision, either way
    assert bs.addition_spot_truth_usable(g_ims, FC_ROUNDED, low, dec_keep) == (True, "changelog entry 10: kept post hoc")
    assert bs.addition_spot_truth_usable(g_ims, FC_ROUNDED, ok, dec_drop) == (False, "changelog entry 10: 4.4 failed")
    # the xlsx source: V1 to V6 must pass whatever is decided; V7 is the 4.4 check, decided like the others
    assert bs.addition_spot_truth_usable(g_x, _xlsx_params())[0] is True
    usable, why = bs.addition_spot_truth_usable(g_x, _xlsx_params(V3_surfaces=1, V7_truth_assignment=0.9899), None,
                                                dec_keep)
    assert usable is False and why == "xlsx spots failed V3, V7 of entry 9 (e)"
    assert bs.addition_spot_truth_usable(g_x, {**FC_ROUNDED, "spots_source": {"source": "xlsx"}})[0] is False
    v7 = _xlsx_params(V7_truth_assignment=0.9899)
    assert bs.addition_spot_truth_usable(g_x, v7)[0] is False
    assert bs.addition_spot_truth_usable(g_x, v7, None, dec_keep)[0] is True
    no_fc = {k: v for k, v in _xlsx_params().items() if k != "frame_check"}
    assert bs.addition_spot_truth_usable(g_x, no_fc, None, dec_keep)[0] is False
    from germquant.validate import imaris_ims

    assert bs.SECTION_4_4_MIN == imaris_ims.TRUTH_ASSIGNMENT_MIN == pytest.approx(0.9935)   # one minimum, V7 and 4.4


def _write_decisions(path: Path, decisions: dict) -> Path:
    path.write_text(json.dumps({"changelog_entry": 10, "decisions": decisions}), encoding="utf-8")
    return path


def test_the_decisions_file_is_validated(tmp_path):
    good = {"usable": False, "entry": 10, "reason": "r"}
    assert bs.load_spot_truth_decisions(_write_decisions(tmp_path / "a.json", {MALES[0]: good}))["decisions"]
    for bad, match in (({"n2rad51_HS_male_000": good}, "not an addition"),
                       ({MALES[0]: {**good, "usable": "no"}}, "true or false"),
                       ({MALES[0]: {**good, "entry": 9}}, "10 or later"),
                       ({MALES[0]: {**good, "reason": " "}}, "reason"),
                       ({}, "decisions")):
        with pytest.raises(ValueError, match=match):
            bs.load_spot_truth_decisions(_write_decisions(tmp_path / "b.json", bad))


def test_a_decision_removes_the_fold_of_an_addition(tmp_path):
    """The manifest carries the section 4.4 table and the decisions file (path and sha256); an addition decided
    unusable takes no part in the leave-one-gonad-out folds, the one kept by its 4.4 agreement does."""
    from germquant.benchmark import logo
    from germquant.benchmark.manifest import load_manifest

    rec, _cm = _stub_additions(tmp_path / "a")
    data = tmp_path / "data"
    data.mkdir()
    entries = []
    for tag, role in bs.PREREGISTERED_ROLES.items():
        p = data / f"{tag}.nd2"
        p.write_bytes(tag.encode())
        q = data / f"{tag}.ims"
        q.write_bytes(tag.encode())
        entries.append((tag, role, p, q))
    split = _synthetic_split(entries)
    truth = tmp_path / "truth"
    for g in rec["gonads"]:
        d = truth / g["tag"]
        d.mkdir(parents=True)
        params = _xlsx_params() if g["tag"] == HS1 else FC_ROUNDED
        (d / f"{g['stem']}__imaris_params.json").write_text(json.dumps(params), encoding="utf-8")
    check = tmp_path / "check.csv"
    check.write_text("tag,split,n_spots,n_agree,n_coloc_on_removed_fragment,n_both_different,agreement\n"
                     f"{MALES[0]},train,5187,5087,95,0,0.9807210333526123\n"
                     f"{HS1},train,6512,6486,10,0,0.996\n", encoding="utf-8")

    def folds(decisions: dict | None):
        dpath = _write_decisions(tmp_path / "dec.json", decisions) if decisions is not None else None
        m = bs.build_manifest(split, tmp_path / "results", truth, additions=rec, truth_check=check,
                              spot_truth_decisions=dpath)
        out = bs.write_manifest(m, tmp_path / "benchmarks", "manifest_entry9.json")
        return m, {f["held_out"] for f in logo.logo_folds(load_manifest(out), require_cache=False)}

    m, held = folds(None)
    add = {e["tag"]: e for e in m["gonads"] if e.get("addition")}
    assert add[MALES[0]]["spot_truth_usable"] is False and add[MALES[0]]["section_4_4"]["agreement"] == pytest.approx(0.9807, abs=1e-4)
    assert add[HS1]["spot_truth_usable"] is True
    assert MALES[0] not in held and HS1 in held
    assert m["spot_truth_inputs"]["truth_check"]["sha256"] == bs.file_sha256(check)
    m, held = folds({MALES[0]: {"usable": True, "entry": 10, "reason": "kept"},
                     HS1: {"usable": False, "entry": 10, "reason": "dropped"}})
    assert MALES[0] in held and HS1 not in held
    assert m["spot_truth_inputs"]["decisions"]["sha256"] == bs.file_sha256(tmp_path / "dec.json")
    assert m["spot_truth_inputs"]["decisions"]["changelog_entry"] == 10
    assert {e["tag"]: e["spot_truth_reason"] for e in m["gonads"] if e.get("addition")}[HS1] == "changelog entry 10: dropped"
    with pytest.raises(ValueError, match="not in the additions record"):
        folds({MALES[1]: {"usable": True, "entry": 10, "reason": "x"}})


DECISIONS = ROOT / "tests" / "golden" / "rad51_additions_decisions.json"


def test_the_decisions_of_entry_10():
    doc = bs.load_spot_truth_decisions(DECISIONS)
    assert doc["changelog_entry"] == 10 and set(doc["decisions"]) == set(bs.ADDITION_ROLES)
    assert {t: d["usable"] for t, d in doc["decisions"].items()} == {
        MALES[0]: False, MALES[1]: True, MALES[2]: False, HS1: False}
    assert all(d["entry"] == 10 for d in doc["decisions"].values())


def test_the_decisions_on_the_real_truth_files():
    """On the workstation: the rule with the real truth params, the 4.4 table and the decisions gives what entry 10
    states; without the decisions only n2rad51_20260105a_HS_male_002 passes section 4.4."""
    truth_root = Path("C:/Users/ryane/germquant_golden/imaris_truth")
    table = Path("C:/Users/ryane/germquant_golden/benchmarks/n2rad51_8666f07/entry9/truth_assignment_check_additions.csv")
    if not table.is_file() or not truth_root.is_dir():
        pytest.skip("the entry-9 truth files are not on this machine")
    rec = bs.load_additions(ADDITIONS)
    checks = bs.load_truth_check(table)
    dec = bs.load_spot_truth_decisions(DECISIONS)["decisions"]
    with_dec, without = {}, {}
    for g in rec["gonads"]:
        params = json.loads((truth_root / g["tag"] / f"{g['stem']}__imaris_params.json").read_text(encoding="utf-8"))
        with_dec[g["tag"]] = bs.addition_spot_truth_usable(g, params, checks[g["tag"]], dec[g["tag"]])[0]
        without[g["tag"]] = bs.addition_spot_truth_usable(g, params, checks[g["tag"]])[0]
    assert with_dec == without == {MALES[0]: False, MALES[1]: True, MALES[2]: False, HS1: False}


def test_registered_stems_are_the_split_files():
    assert bs.REGISTERED_STEMS == {g["tag"]: g["stem"] for g in bs.load_split(SPLIT)["gonads"]}
    assert set(bs.REGISTERED_STEMS) == set(bs.PREREGISTERED_ROLES)


def test_the_truth_read_guard():
    """guard_truth_read: a registered train or val project passes by its stem; a test gonad's project is refused
    unless allow_test (the pre-flight); any other sealed file is refused; an addition's files pass."""
    assert bs.registered_tag_of("C:/x/20251105_n2_hs_MALE _001.ims") == "n2rad51_HS_male_001"
    assert bs.registered_tag_of("C:/x/20251105_N2_NOHS_herm _.nd2") == "n2rad51_noHS_herm_000"
    assert bs.registered_tag_of("C:/x/20251105_n2_hs_MALE_1.xlsx") is None
    assert bs.guard_truth_read(["C:/x/20251105_n2_hs_HERM _.ims", "C:/x/20251105_n2_hs_HERM _.nd2", None]) == [
        "n2rad51_HS_herm_000"]
    with pytest.raises(bd.TestSetGuardError):
        bs.guard_truth_read(["C:/x/20251105_n2_hs_MALE _001.ims"])
    with pytest.raises(bd.TestSetGuardError):
        bs.guard_truth_read(["C:/x/20251105_n2_hs_HERM _.ims", "C:/x/20251105_n2_hs_MALE _001.nd2"])
    assert bs.guard_truth_read(["C:/x/20251105_n2_hs_MALE _001.ims"], allow_test=True) == ["n2rad51_HS_male_001"]
    for sealed in ("C:/nas/20251105_n2_hs_MALE_1.xlsx", "C:/x/temp_compiledtdata.xlsx", "C:/x/20251105_other.ims"):
        with pytest.raises(PermissionError, match="sealed"):
            bs.guard_truth_read(["C:/y/20260105_n2_hs_MALE_1.ims", sealed])
    assert bs.guard_truth_read(["C:/y/20260105_n2_hs_MALE_1.ims", "C:/y/20260105_n2_hs_MALE_1.nd2"]) == []
CP_SPOTS = ('<bpPointsCreationParameters mEnableRegionOfInterest="false" mEnableRegionGrowing="false" '
            'mSourceChannelIndex="2" mEstimatedDiameter="0.400000 0.400000 0.900000" mBackgroundSubtraction="true">\n'
            '<mSpotFilter>\n<bpStatisticsFilter mLowerThresholdEnable="true" mLowerThresholdManual="true" '
            'mLowerThresholdManualValue="18.9129" mUpperThresholdEnable="false" mUpperThresholdManualValue="1">\n'
            '  <bpStatisticsValueType mName="Quality" mUnit="" mFactors="0"/>\n</bpStatisticsFilter>\n</mSpotFilter>\n'
            '<bpCObjectsClassifyFilter1D>\n<bpStatisticsValueType mName="Average Distance To 3 Nearest Neighbours"/>\n'
            '<cThresholdValue mClassName="0" mThresholdLow="-3.4e+38" mThresholdHigh="13.6536"/>\n'
            '</bpCObjectsClassifyFilter1D>\n</bpPointsCreationParameters>\n')
CP_SURF = ('<bpSurfacesCreationParameters mEnableRegionOfInterest="true" mEnableRegionGrowing="true" '
           'mSourceChannelIndex="0" mEnableSmooth="true" mSurfaceGrainSize="0.4" mEnableEliminateBackground="true" '
           'mDiameterOfLargestSphere="3" mEnableAutomaticThreshold="true" mManualThresholdValue="111.234" '
           'mRegionGrowingEstimatedDiameter="3">\n<bpRegionOfInterestContainer>\n'
           '<bpRegionOfInterest mMinXYZ="0 1306 0" mMaxXYZ="1665 2719 105" mName="R"/>\n</bpRegionOfInterestContainer>\n'
           '<mSurfacesFilter>\n<bpStatisticsFilter mLowerThresholdEnable="true" mLowerThresholdManual="true" '
           'mLowerThresholdManualValue="10" mUpperThresholdEnable="false" mUpperThresholdManualValue="1">\n'
           '  <bpStatisticsValueType mName="Number of Voxels" mUnit="" mFactors="1"/>\n</bpStatisticsFilter>\n'
           '</mSurfacesFilter>\n</bpSurfacesCreationParameters>\n')


def test_creation_parameter_summary_and_comparison():
    s = bs.creation_parameter_summary(CP_SPOTS, CP_SURF)
    assert s["spots"] == {"source_channel_index": 2, "estimated_diameter_um": [0.4, 0.4, 0.9],
                          "background_subtraction": True, "roi_enabled": False, "region_growing": False,
                          "quality_threshold": 18.9129, "spot_filters": ["Quality"],
                          "class_statistic": "Average Distance To 3 Nearest Neighbours", "class_threshold_um": 13.6536,
                          "roi_boxes": []}
    su = s["surfaces"]
    assert (su["largest_sphere_um"], su["seed_diameter_um"], su["automatic_threshold"]) == (3.0, 3.0, True)
    assert su["voxel_filter_min"] == 10.0 and su["voxel_filter_max"] is None and su["voxel_filter_enabled"]
    assert su["roi_boxes"] == [{"min_xyz_voxel": [0, 1306, 0], "max_xyz_voxel": [1665, 2719, 105]}]
    reg = {"a": bs.creation_parameter_summary(CP_SPOTS.replace("13.6536", "7.85").replace("18.9129", "42.5"),
                                              CP_SURF.replace('mEnableAutomaticThreshold="true"',
                                                              'mEnableAutomaticThreshold="false"'))}
    diff = {d["parameter"]: d for d in bs.compare_creation_parameters(s, reg)}
    assert diff["spots.class_threshold_um"]["registered_range"] == [7.85, 7.85]
    assert diff["spots.quality_threshold"]["registered_range"] == [42.5, 42.5]
    assert diff["surfaces.largest_sphere_um"]["protocol_text"] == 2.5 and "registered_values" not in diff["surfaces.largest_sphere_um"]
    assert diff["surfaces.automatic_threshold"]["registered_values"] == [False]
    assert "spots.estimated_diameter_um" not in diff                     # as in section 4.1 and the registered file
