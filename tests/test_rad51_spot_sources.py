"""The recorded spot sources of protocol changelog entry 11 (tests/golden/rad51_spot_sources.json,
germquant.validate.benchmark_split): n2rad51_HS_male_000's spot truth from its 2025-12-31 export. The record validates
against the split and the pins in code (a registered TRAIN gonad only, the pinned export by name, size and sha256);
the guard admits exactly the pinned file and keeps every other 2025-11-05 export sealed; the usability rule (the
export's own checks E1 to E9 and the frame check, then a decision, else section 4.4); the manifest changes only that
gonad's entry, which then takes part in the leave-one-gonad-out folds and in net3d's training gonads while
n2rad51_HS_male_002 stays out; the net3d fold order. Synthetic, apart from the checks of the real files, which skip
off the workstation."""
from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path

import pytest
from test_benchmark_split import _synthetic_split

from germquant.benchmark import detector as bd
from germquant.benchmark import logo
from germquant.benchmark.manifest import load_manifest
from germquant.validate import benchmark_split as bs

ROOT = Path(__file__).resolve().parents[1]
SPLIT = ROOT / "tests" / "golden" / "rad51_split.json"
SOURCES = ROOT / "tests" / "golden" / "rad51_spot_sources.json"
TAG = "n2rad51_HS_male_000"
RG2 = "n2rad51_HS_male_002"
PIN_SHA = "0fbc41b7dc5dbf0a1ef0c87f367451635c0f7dfea1302e6dfc53b614082aca52"
# the split file the record of entry 11 was made against (as registered; CRLF and LF). The record is data cited by
# entry 11 and stays as committed when the split moves by the re-export of protocol section 4.1 (entry 10 (g)(8)):
# these historical digests are pinned here, and the roles are checked live (protocol changelog entry 12)
SPLIT_AT_ENTRY_11 = ("955eeff52b4e953ca34520b2eae6a5a019520203e79967082a78eac4c991e25f",
                     "989bfa08e4acd911004e433469979a2cfd52891aee240ac2f87db9183f0c9d85")


def _record() -> dict:
    return json.loads(SOURCES.read_text(encoding="utf-8"))


def test_the_record_and_its_pins():
    rec = bs.load_spot_sources(SOURCES, bs.load_split(SPLIT))
    assert rec["changelog_entry"] == bs.SPOT_SOURCES_CHANGELOG_ENTRY == 11
    assert [s["tag"] for s in rec["sources"]] == [TAG] and set(bs.SPOT_SOURCE_EXPORTS) == {TAG}
    assert rec["decided_by"]["reply"] == "9 yes"
    reg = rec["registered_split"]
    assert (reg["sha256"], reg["sha256_lf"]) == SPLIT_AT_ENTRY_11            # the split the record was made against
    assert reg["roles_sha256"] == bs.roles_sha256() and bs.split_file_unchanged(SPLIT)   # no role moved
    s = bs.spot_source_of(rec, TAG)
    assert s["export"]["name"] == "20251105_n2_hs_MALE_3.xlsx" and s["export"]["size_bytes"] == 15512227
    assert s["export"]["sha256"] == PIN_SHA == bs.SPOT_SOURCE_EXPORTS[TAG]["sha256"]
    assert s["expect"] == {"n_spots": 7628, "n_surfaces": 64, "class_counts": {"Class A": 6497, "Class B": 1131},
                           "spot_diameter_um": [0.4, 0.4, 0.9]}
    assert bs.spot_source_of(rec, RG2) is None
    # only registered train gonads can be pinned, and never with the test male's export name
    assert all(bs.PREREGISTERED_ROLES[t] == "train" for t in bs.SPOT_SOURCE_EXPORTS)
    assert all("male_1" not in p["name"].lower() for p in bs.SPOT_SOURCE_EXPORTS.values())


def test_the_pinned_export_is_on_disk():
    rec = bs.load_spot_sources(SOURCES)
    if not Path(rec["sources"][0]["export"]["path"]).is_file():
        pytest.skip("the 2025-12-31 export copy is not on this machine")
    assert bs.check_spot_source_files(rec, hashes=True) == []


@pytest.mark.parametrize("mutate,match", [
    (lambda r: r["sources"][0].update(tag="n2rad51_HS_male_001"), "registered train gonad only"),
    (lambda r: r["sources"][0].update(tag="n2rad51_HS_herm_001"), "registered train gonad only"),
    (lambda r: r["sources"][0].update(role="val"), "registered train gonad only"),
    (lambda r: r["sources"][0].update(tag="n2rad51_20260105a_HS_male_001"), "registered gonad"),
    (lambda r: r["sources"][0].update(tag=RG2, stem="20251105_n2_hs_MALE _002"), "no export is pinned"),
    (lambda r: r["sources"][0].update(stem="20251105_n2_hs_MALE _001"), "registered stem"),
    (lambda r: r["sources"][0].update(kind="reexport"), "kind"),
    (lambda r: r["sources"][0]["export"].update(name="20251105_n2_hs_MALE_1.xlsx"), "name is not the pinned"),
    (lambda r: r["sources"][0]["export"].update(size_bytes=20980766), "size_bytes is not the pinned"),
    (lambda r: r["sources"][0]["export"].update(sha256="0" * 64), "sha256 is not the pinned"),
    (lambda r: r["sources"][0]["export"].update(path="C:/x/20251105_n2_hs_MALE_1.xlsx"), "does not name the pinned"),
    (lambda r: r["sources"][0]["export"].update(nas_path="//nas/20251105_n2_hs_MALE_2.xlsx"), "does not name the pinned"),
    (lambda r: r["sources"][0]["expect"].update(class_counts={"Class A": 6498, "Class B": 1131}), "add up"),
    (lambda r: r["sources"][0]["expect"].update(n_surfaces=0), "positive"),
    (lambda r: r["sources"][0]["ims"].update(sha256="1" * 64), "not the registered one"),
    (lambda r: r["registered_split"].update(roles_sha256="0" * 64), "other registered roles"),
    (lambda r: r.update(changelog_entry=10), "changelog entry"),
    (lambda r: r["sources"].append(copy.deepcopy(r["sources"][0])), "unique"),
    (lambda r: r.update(sources=[]), "no 'sources'"),
])
def test_validate_spot_sources_refuses(mutate, match):
    rec = _record()
    mutate(rec)
    with pytest.raises(ValueError, match=match):
        bs.validate_spot_sources(rec, bs.load_split(SPLIT))


def test_load_spot_sources_validates(tmp_path):
    """The CLI and the scripts read the record through load_spot_sources: a malformed record is refused there, with
    and without the split (a source on the test male; an .ims that is not the split's)."""
    rec = _record()
    rec["sources"][0]["tag"] = "n2rad51_HS_male_001"
    p = tmp_path / "rad51_spot_sources.json"
    p.write_text(json.dumps(rec), encoding="utf-8")
    with pytest.raises(ValueError, match="registered train gonad only"):
        bs.load_spot_sources(p)
    with pytest.raises(ValueError, match="registered train gonad only"):
        bs.load_spot_sources(p, bs.load_split(SPLIT))
    rec = _record()
    rec["sources"][0]["ims"]["size_bytes"] += 1
    p.write_text(json.dumps(rec), encoding="utf-8")
    assert bs.load_spot_sources(p)["sources"][0]["tag"] == TAG                  # without the split: not compared
    with pytest.raises(ValueError, match="not the registered one"):
        bs.load_spot_sources(p, bs.load_split(SPLIT))


def _pin(monkeypatch, tmp_path, tag=TAG, name="20251105_n2_hs_MALE_3.xlsx", data=b"export bytes") -> Path:
    import hashlib

    tmp_path.mkdir(parents=True, exist_ok=True)
    p = tmp_path / name
    p.write_bytes(data)
    monkeypatch.setattr(bs, "SPOT_SOURCE_EXPORTS", {tag: {"name": name, "size_bytes": len(data),
                                                          "sha256": hashlib.sha256(data).hexdigest()}})
    return p


def test_the_guard_admits_exactly_the_pinned_export(tmp_path, monkeypatch):
    """By name and size, before anything is read: the pinned file passes as its gonad's (train), the same name at
    another size and every other 2025-11-05 export stay sealed, and a pin on a test gonad is still refused."""
    real_pin = dict(bs.SPOT_SOURCE_EXPORTS)
    wrong_size = tmp_path / "w" / real_pin[TAG]["name"]
    wrong_size.parent.mkdir()
    wrong_size.write_bytes(b"not the export")
    assert bs.spot_source_tag_of(wrong_size) is None
    with pytest.raises(PermissionError, match="sealed"):
        bs.guard_export_read(wrong_size)
    assert bs.spot_source_tag_of(tmp_path / "missing" / real_pin[TAG]["name"]) is None
    p = _pin(monkeypatch, tmp_path)
    assert bs.spot_source_tag_of(p) == TAG and bs.guard_export_read(p) == TAG
    ims = "C:/x/20251105_n2_hs_MALE _.ims"
    assert bs.guard_truth_read([ims, p]) == [TAG]
    for sealed in ("C:/nas/20251105_n2_hs_MALE_1.xlsx", "C:/nas/20251105_n2_hs_MALE_2.xlsx",
                   "C:/nas/20251105_n2_hs_MALE_1.xls", "C:/x/temp_compiledtdata.xlsx"):
        with pytest.raises(PermissionError, match="sealed"):
            bs.guard_export_read(sealed)
        with pytest.raises(PermissionError, match="sealed"):
            bs.guard_truth_read([ims, sealed])
    with pytest.raises(bd.TestSetGuardError):                       # the test male's project, as before
        bs.guard_truth_read(["C:/x/20251105_n2_hs_MALE _001.ims", p])
    q = _pin(monkeypatch, tmp_path / "t", tag="n2rad51_HS_male_001")  # never pinned in code; the role guard holds
    with pytest.raises(bd.TestSetGuardError):
        bs.guard_export_read(q)
    with pytest.raises(bd.TestSetGuardError):
        bs.guard_truth_read([q])


def test_the_guard_refuses_the_pinned_name_and_size_with_other_bytes(tmp_path, monkeypatch):
    """Entry 12: a file with the pinned name and size but other bytes (a renamed export of another gonad, say) is
    hashed and refused by both guards before anything parses it; the role guard comes before the hash."""
    from germquant.validate import imaris_xlsx as X

    p = _pin(monkeypatch, tmp_path / "pin")
    other = tmp_path / "other" / p.name
    other.parent.mkdir()
    other.write_bytes(b"EXPORT BYTES")                             # same name, same size, other bytes
    assert other.stat().st_size == p.stat().st_size and bs.spot_source_tag_of(other) == TAG
    opened = []
    import openpyxl

    def no_open(*a, **k):
        opened.append(a)
        raise AssertionError("the workbook was opened")

    monkeypatch.setattr(openpyxl, "load_workbook", no_open)
    with pytest.raises(PermissionError, match="not its sha256"):
        bs.guard_export_read(other)
    with pytest.raises(PermissionError, match="not its sha256"):
        bs.guard_truth_read(["C:/x/20251105_n2_hs_MALE _.ims", other])
    with pytest.raises(PermissionError, match="not its sha256"):
        X.read_export_objects(other, spot_channel=2, mask_channel=3)
    assert opened == []
    assert bs.guard_export_read(p) == TAG                          # the pinned bytes still pass
    hashed = []
    monkeypatch.setattr(bs, "file_sha256", lambda path, *a, **k: hashed.append(path) or "0" * 64)
    q = _pin(monkeypatch, tmp_path / "t", tag="n2rad51_HS_male_001")
    with pytest.raises(bd.TestSetGuardError):                      # a test tag is refused before it is hashed
        bs.guard_export_read(q)
    assert hashed == []


def test_the_export_reader_uses_the_guard(tmp_path, monkeypatch):
    """imaris_xlsx opens a sealed name only when it is the pinned export (here a pinned stand-in that is no
    workbook, so openpyxl fails after the guard passed); the test male's export is refused before opening."""
    from germquant.validate import imaris_xlsx as X

    with pytest.raises(PermissionError, match="sealed"):
        X.read_export_objects(tmp_path / "20251105_n2_hs_MALE_1.xlsx", spot_channel=2, mask_channel=3)
    p = _pin(monkeypatch, tmp_path)
    with pytest.raises(Exception) as e:
        X.read_export_objects(p, spot_channel=2, mask_channel=3)
    assert not isinstance(e.value, PermissionError)


def _params(fails=(), e10=0.9967, fc=True, kind=bs.SPOT_SOURCE_KIND, sha=PIN_SHA) -> dict:
    ver = {f"E{i}_x": {"passes": f"E{i}" not in fails} for i in range(1, 10)}
    ver["E10_truth_assignment"] = {"passes": "E10" not in fails, "agreement": e10}
    out = {"spots_source": {"source": kind, "export": {"sha256": sha}, "verification": ver}}
    if fc:
        out["frame_check"] = {"chosen": {"name": "stored"}}
    return out


def test_the_spot_source_rule():
    """spot_source_truth_usable: no pinned source, no export, an export not from the source, a failed E1 to E9 or a
    missing frame check: never usable (no decision overrides); then a decision; else section 4.4 (the table's row,
    or E10's own agreement)."""
    src = bs.spot_source_of(_record(), TAG)
    keep = {"usable": True, "entry": 12, "reason": "kept"}
    drop = {"usable": False, "entry": 12, "reason": "dropped"}
    assert bs.spot_source_truth_usable(None, _params()) == (False, "no recorded spot source (entry 11)")
    assert bs.spot_source_truth_usable({**src, "export": {**src["export"], "sha256": "0" * 64}}, _params())[0] is False
    assert bs.spot_source_truth_usable(src, None, None, keep)[0] is False
    assert bs.spot_source_truth_usable(src, _params(kind="xlsx"), None, keep)[0] is False
    assert bs.spot_source_truth_usable(src, _params(sha="0" * 64), None, keep)[0] is False
    usable, why = bs.spot_source_truth_usable(src, _params(fails=("E3", "E10")), None, keep)
    assert usable is False and why == "the recorded spot source failed E3, E10 of entry 11"
    assert bs.spot_source_truth_usable(src, _params(fc=False), None, keep)[0] is False
    assert bs.spot_source_truth_usable(src, {**_params(), "spots_source": {"source": bs.SPOT_SOURCE_KIND,
                                                                           "export": {"sha256": PIN_SHA}}})[0] is False
    usable, why = bs.spot_source_truth_usable(src, _params())
    assert usable is True and why == "section 4.4 agreement 0.9967, at least 0.9935 (entry 11)"
    low = _params(fails=("E10",), e10=0.9915)
    usable, why = bs.spot_source_truth_usable(src, low)
    assert usable is False and "needs a changelog entry" in why and "0.9915" in why
    assert bs.spot_source_truth_usable(src, low, None, keep) == (True, "changelog entry 12: kept")
    assert bs.spot_source_truth_usable(src, _params(), None, drop)[0] is False
    # the 4.4 table's row decides over E10's own number
    assert bs.spot_source_truth_usable(src, _params(), {"agreement": 0.99})[0] is False
    assert bs.spot_source_truth_usable(src, low, {"agreement": 0.996})[0] is True


def _synthetic(tmp_path) -> tuple[dict, dict]:
    """A synthetic split of the twelve registered tags, the two heat-shock males of train with region growing as in
    the registered record, and a spot-sources record whose .ims is the synthetic one."""
    data = tmp_path / "data"
    data.mkdir()
    entries = []
    for tag, role in bs.PREREGISTERED_ROLES.items():
        p = data / f"{bs.REGISTERED_STEMS[tag]}.nd2"
        p.write_bytes(tag.encode())
        q = data / f"{bs.REGISTERED_STEMS[tag]}.ims"
        q.write_bytes(tag.encode())
        entries.append((tag, role, p, q))
    split = _synthetic_split(entries)
    for g in split["gonads"]:
        if g["tag"] in (TAG, RG2):
            g["ims"]["region_growing"] = True
    rec = _record()
    g0 = next(g for g in split["gonads"] if g["tag"] == TAG)
    rec["sources"][0]["ims"] = {"path": g0["ims"]["path"], "size_bytes": g0["ims"]["size_bytes"],
                                "sha256": g0["ims"]["sha256"]}
    return split, rec


def _write_source_export(root: Path, params: dict) -> None:
    d = root / TAG
    d.mkdir(parents=True)
    (d / f"{bs.REGISTERED_STEMS[TAG]}__imaris_params.json").write_text(json.dumps(params), encoding="utf-8")


def test_the_manifest_with_the_spot_source(tmp_path):
    """Only n2rad51_HS_male_000's entry changes: its truth folder is the export from the source, its spot truth
    follows the rule with its row of the 4.4 table, the registered folder is kept beside it; HS_male_002 stays
    excluded. The folds then hold HS_male_000 and not HS_male_002."""
    split, rec = _synthetic(tmp_path)
    results, truth, src_root = tmp_path / "results", tmp_path / "truth", tmp_path / "truth_sources"
    plain = bs.build_manifest(split, results, truth)
    _write_source_export(src_root, _params())
    recf = tmp_path / "sources.json"
    recf.write_text(json.dumps(rec), encoding="utf-8")
    check = tmp_path / "check_hs_male_000.csv"
    check.write_text("tag,split,n_spots,n_agree,n_coloc_on_removed_fragment,n_both_different,agreement\n"
                     f"{TAG},train,7628,7603,20,0,0.996722600943891\n", encoding="utf-8")
    m = bs.build_manifest(split, results, truth, spot_sources=rec, spot_sources_path=recf,
                          spot_sources_truth_root=src_root, truth_check=check)
    by_plain = {e["tag"]: e for e in plain["gonads"]}
    by = {e["tag"]: e for e in m["gonads"]}
    assert [t for t in by if by[t] != by_plain[t]] == [TAG]
    e = by[TAG]
    assert e["spot_truth_usable"] is True and by_plain[TAG]["spot_truth_usable"] is False
    assert e["truth_dir"] == (src_root / TAG).as_posix() and e["registered_truth_dir"] == by_plain[TAG]["truth_dir"]
    assert e["spot_truth_reason"] == "section 4.4 agreement 0.9967, at least 0.9935 (entry 11)"
    assert e["section_4_4"]["n_agree"] == 7603 and e["recorded_spot_source"]["export"]["sha256"] == PIN_SHA
    assert e["notes"][-1] == bs.SPOT_SOURCE_NOTE and e["ims_record"]["region_growing"] is True
    assert by[RG2]["spot_truth_usable"] is False
    assert m["spot_sources_file"]["sha256"] == bs.file_sha256(recf) and m["spot_sources_file"]["tags"] == [TAG]
    assert m["spot_sources_file"]["truth_check"]["sha256"] == bs.file_sha256(check)
    assert TAG in m["summary"]["spot_truth_usable"] and RG2 not in m["summary"]["spot_truth_usable"]
    out = bs.write_manifest(m, tmp_path / "benchmarks", "manifest_entry9.json")
    held = {f["held_out"] for f in logo.logo_folds(load_manifest(out), require_cache=False)}
    assert TAG in held and RG2 not in held and len(held) == 7
    # a failed criterion of the export, or a 4.4 row below the minimum: excluded, and no fold
    for name, params, row in (("fail_e8", _params(fails=("E8",)), "0.996722600943891"), ("low", _params(), "0.99")):
        root = tmp_path / name
        _write_source_export(root, params)
        check.write_text("tag,split,n_spots,n_agree,n_coloc_on_removed_fragment,n_both_different,agreement\n"
                         f"{TAG},train,7628,7603,20,0,{row}\n", encoding="utf-8")
        m2 = bs.build_manifest(split, results, truth, spot_sources=rec, spot_sources_truth_root=root,
                               truth_check=check)
        assert {x["tag"]: x for x in m2["gonads"]}[TAG]["spot_truth_usable"] is False
    # the options belong together, and a gonad in two 4.4 tables is refused
    with pytest.raises(ValueError, match="spot_sources_truth_root"):
        bs.build_manifest(split, results, truth, spot_sources=rec)
    with pytest.raises(ValueError, match="additions only"):
        bs.build_manifest(split, results, truth, truth_check=check)
    with pytest.raises(ValueError, match="another section 4.4 table"):
        bs.build_manifest(split, results, truth, spot_sources=rec, spot_sources_truth_root=src_root,
                          truth_check=[check, check])
    # a record whose .ims is not the split's is refused
    with pytest.raises(ValueError, match="not the registered one"):
        bs.build_manifest(split, results, truth, spot_sources=_record(), spot_sources_truth_root=src_root)


def test_the_manifest_takes_several_4_4_tables(tmp_path):
    """The additions' table of entry 10 and the spot sources' table of entry 11 together; one entry per table in
    spot_truth_inputs, and each gonad reads its own row."""
    a = tmp_path / "a.csv"
    a.write_text("tag,split,n_spots,n_agree,agreement\nx,train,10,10,1.0\n", encoding="utf-8")
    b = tmp_path / "b.csv"
    b.write_text("tag,split,n_spots,n_agree,agreement,agreement_like_for_like\ny,train,10,9,0.9,0.95\n",
                 encoding="utf-8")
    rows = bs.load_truth_check([a, b])
    assert set(rows) == {"x", "y"} and rows["y"]["agreement_like_for_like"] == pytest.approx(0.95)
    assert rows["x"]["agreement_like_for_like"] is None
    rec = bs._truth_check_record([a, b])
    assert [r["sha256"] for r in rec] == [bs.file_sha256(a), bs.file_sha256(b)]
    assert bs._truth_check_record(a)["sha256"] == bs.file_sha256(a)


def test_net3d_takes_the_spot_source_gonad(tmp_path):
    """resolve_gonads: the split's region-growing flag describes the .ims's saved Spots object; for a gonad whose
    manifest entry carries the pinned spot source the manifest decides. HS_male_002 stays excluded."""
    from germquant.spots.net3d.train import resolve_gonads

    split, rec = _synthetic(tmp_path)
    split_file = tmp_path / "split.json"
    split_file.write_text(json.dumps(split), encoding="utf-8")
    src_root = tmp_path / "truth_sources"
    _write_source_export(src_root, _params())
    m = bs.build_manifest(split, tmp_path / "results", tmp_path / "truth", spot_sources=rec,
                          spot_sources_truth_root=src_root)
    # an addition's entry names its source as a string ("ims" or "xlsx"), which the rule must pass over
    m["gonads"].append({"tag": "n2rad51_20260105a_HS_male_002", "role": "train", "split": "train", "addition": True,
                        "results_dir": str(tmp_path / "r"), "spot_truth_source": "ims", "spot_truth_usable": True})
    out = bs.write_manifest(m, tmp_path / "benchmarks", "manifest_entry9.json")
    gon = resolve_gonads(split_file, out)
    assert TAG in gon["train"] and RG2 not in gon["train"]
    assert RG2 in {x["tag"] for x in gon["excluded"]}
    assert not set(gon["train"] + gon["val"]) & {t for t, r in bs.PREREGISTERED_ROLES.items() if r == "test"}
    # without the source (or with another export sha256) the region-growing flag excludes it again
    doc = json.loads(out.read_text(encoding="utf-8"))
    for e in doc["gonads"]:
        if e["tag"] == TAG:
            e["recorded_spot_source"]["export"]["sha256"] = "0" * 64
    out.write_text(json.dumps(doc), encoding="utf-8")
    assert TAG not in resolve_gonads(split_file, out)["train"]


def test_the_net3d_fold_order():
    """bench_logo.py net3d-fold: the recipe's train gonads, its val gonads, the registered gonads the recipe did not
    see (in tag order), then the additions in record order; the held-out gonad left out; a gonad that is neither
    registered nor an addition given is left out, so the caller refuses the fold."""
    spec = importlib.util.spec_from_file_location("bench_logo", ROOT / "scripts" / "bench_logo.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    recipe = {"from_train_gonads": ["n2rad51_HS_herm_000", "n2rad51_noHS_herm_000", "n2rad51_noHS_male_000",
                                    "n2rad51_noHS_male_001"],
              "from_val_gonads": ["n2rad51_HS_herm_001", "n2rad51_noHS_herm_001"]}
    adds = ["n2rad51_20260105a_HS_male_001", "n2rad51_20260105a_HS_male_002"]
    every = recipe["from_train_gonads"] + recipe["from_val_gonads"] + [TAG] + adds
    fold = {"held_out": "n2rad51_noHS_male_000", "train": [t for t in every if t != "n2rad51_noHS_male_000"],
            "val": []}
    assert mod.fold_training_order(recipe, fold, adds) == [
        "n2rad51_HS_herm_000", "n2rad51_noHS_herm_000", "n2rad51_noHS_male_001", "n2rad51_HS_herm_001",
        "n2rad51_noHS_herm_001", TAG, *adds]
    fold = {"held_out": TAG, "train": [t for t in every if t != TAG], "val": []}
    assert mod.fold_training_order(recipe, fold, adds) == recipe["from_train_gonads"] + recipe["from_val_gonads"] + adds
    fold = {"held_out": adds[0], "train": [t for t in every if t != adds[0]], "val": []}
    order = mod.fold_training_order(recipe, fold, [])                 # --additions not given
    assert adds[1] not in order and set(order) != set(fold["train"])  # the caller's check refuses this fold
