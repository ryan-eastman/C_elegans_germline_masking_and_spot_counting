"""The registered SC validation split (tests/golden/sc_split.json): its digest, counts and rules, the RAD-51 seal
(the four RAD-51 test gonads are SC test, never development), the guard and, on the workstation, the sha256 of the
truth copies it pins."""
from __future__ import annotations

import hashlib
from collections import Counter
from pathlib import Path

import pytest

from germquant.benchmark import sc_split

RAD51_TEST = {"n2rad51_HS_herm_002", "n2rad51_HS_male_001", "n2rad51_noHS_herm_002", "n2rad51_noHS_male_002"}


def test_digest_counts_and_roles():
    split = sc_split.load_split()                       # verifies the registered sha256
    assert split["sealed"] is True
    c = Counter((g["geno"], g["role"]) for g in split["gonads"])
    assert c[("N2", "dev")] == 18 and c[("N2", "test")] == 12 and c[("N2", "reserved")] == 9
    assert c[("spo11", "reserved")] == 12 and c[("syp2_3", "reserved")] == 9 and c[("syp3", "excluded")] == 42
    # every N2 condition x sex cell has three test gonads
    cells = Counter((g["condition"], g["sex"]) for g in split["gonads"] if g["geno"] == "N2" and g["role"] == "test")
    assert set(cells.values()) == {3} and len(cells) == 4


def test_rad51_test_gonads_are_sc_test_and_never_dev():
    split = sc_split.load_split()
    by_tag = {g["rad51_tag"]: g for g in split["gonads"] if g.get("rad51_tag")}
    for tag in RAD51_TEST:
        assert by_tag[tag]["role"] == "test"
    dev_tags = {g["rad51_tag"] for g in split["gonads"] if g["role"] == "dev" and g.get("rad51_tag")}
    assert not (dev_tags & RAD51_TEST)


def test_guard_and_sealed_stems():
    split = sc_split.load_split()
    stems = sc_split.sealed_stems(split)
    for tag in RAD51_TEST:
        rec = sc_split.gonad_record(tag, split)
        assert rec["image_stem"] in stems
        with pytest.raises(sc_split.SealedError):
            sc_split.guard(rec["gonad"], split=split)
        sc_split.guard(rec["gonad"], split=split, final=True)
    dev = sc_split.dev_gonads(split)
    assert len(dev) == 18
    sc_split.guard([g["gonad"] for g in dev], split=split)
    assert sc_split.role_of("20251105_N2_nohs_HERM _", split) == "dev"


def test_tampered_split_is_refused(tmp_path):
    src = sc_split.default_path()
    bad = tmp_path / "sc_split.json"
    bad.write_bytes(src.read_bytes().replace(b'"sealed": true', b'"sealed": false'))
    with pytest.raises(ValueError):
        sc_split.load_split(bad)


@pytest.mark.parametrize("key", ["zone_export_copy"])
def test_truth_copies_match_their_sha256(key):
    split = sc_split.load_split()
    checked = 0
    for g in split["gonads"]:
        rec = g.get(key)
        if not rec or not Path(rec["path"]).is_file():
            continue
        assert hashlib.sha256(Path(rec["path"]).read_bytes()).hexdigest() == rec["sha256"], g["gonad"]
        checked += 1
    if not checked:
        pytest.skip("the truth copies are on the lab workstation only")
