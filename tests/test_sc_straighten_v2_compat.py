"""sc_straighten v2 with every v2 switch at its v1 value (germquant.sc.straighten.params.V1_VALUES) reproduces v1:
the tables (all v1 columns; not the version, the hash or columns added in v2) and strips of seeds 2000 to 2004 of
every registered synthetic family and of the end-to-end synthetic gonad equal those of commit 5cb3cb6
(tests/golden/sc_straighten/v1_compat.json, made once with the v1 code). The registered generator of protocol 2.3
also stays identical (image, labels and truth of two seeds per family).

The digests were made on the lab workstation. A difference fails when GERMQUANT_WORKSTATION_PINS=1 (set on the
workstation); elsewhere the test skips with a reason that says what differs (another platform's floating point
may change a last digit), so a CI log still shows whether the platform agrees. The structure of a synthetic case,
its number of strands and its strip index (strand ids, offsets and lengths), fails on every platform (protocol
changelog entry 27): a last-digit difference does not change it, a v1 regression does.
"""
from __future__ import annotations

import json
import os

import pytest
import sc_v1_compat as C

from germquant.sc.straighten import params

GOLDEN = json.loads(C.golden_path().read_text(encoding="ascii"))
ENFORCE = os.environ.get("GERMQUANT_WORKSTATION_PINS") == "1"
V1 = params.v1_sets()
# the structure of a case: compared on every platform (protocol changelog entry 27)
STRUCT = ("n_strands", "strip_index")


def _check(got: dict, want: dict, what: str) -> None:
    hard = [k for k in STRUCT if k in want and got.get(k) != want[k]]
    if hard:
        pytest.fail(f"{what}: structural difference from v1 in {hard} (number of strands {got.get('n_strands')} "
                    f"against {want.get('n_strands')})")
    diffs = sorted(k for k in want if not k.endswith("_columns") and got.get(k) != want[k])
    if not diffs:
        return
    msg = f"{what}: differs from the v1 digests of the workstation in {diffs}"
    if ENFORCE:
        pytest.fail(msg)
    pytest.skip(msg + " (enforced with GERMQUANT_WORKSTATION_PINS=1)")


def test_a_structural_difference_fails_off_the_workstation_too(monkeypatch):
    monkeypatch.setitem(globals(), "ENFORCE", False)
    want = {"n_strands": 7, "strip_index": "a", "strands": "x"}
    with pytest.raises(pytest.fail.Exception):
        _check({**want, "n_strands": 6}, want, "case")
    with pytest.raises(pytest.fail.Exception):
        _check({**want, "strip_index": "b"}, want, "case")
    with pytest.raises(pytest.skip.Exception):                 # a table digest alone: a skip off the workstation
        _check({**want, "strands": "y"}, want, "case")
    _check(dict(want), want, "case")


def test_every_v2_change_has_a_v1_switch():
    p1 = params.resolve(None, V1)
    for key, value in params.V1_VALUES.items():
        assert params.resolve().get(key, KeyError) is not KeyError, key
        assert p1[key] == value, key


@pytest.mark.parametrize("case", sorted(GOLDEN["synthetic"]))
def test_synthetic_tables_equal_v1(case):
    fam, gen, seed = case.split("/")
    got = C.synthetic_digests(fam, gen, int(seed), sets=V1, cols=GOLDEN["columns"])
    _check(got, GOLDEN["synthetic"][case], case)


@pytest.mark.parametrize("case", ["dense/entry2/2000", "pairs/entry2/2001", "gaps/entry2/2002"])
def test_defaults_reproduce_v1(case):
    # protocol changelog entry 25: no change of the v2 round was adopted, so the defaults give the tables of v1
    fam, gen, seed = case.split("/")
    got = C.synthetic_digests(fam, gen, int(seed), sets=[], cols=GOLDEN["columns"])
    _check(got, GOLDEN["synthetic"][case], case + " (defaults)")


@pytest.mark.parametrize("case", sorted(GOLDEN["scenes"]))
def test_registered_scenes_are_unchanged(case):
    fam, gen, seed = case.split("/")
    _check(C.scene_digests(fam, gen, int(seed)), GOLDEN["scenes"][case], case)


def test_end_to_end_gonad_tables_equal_v1(tmp_path):
    got = C.gonad_digests(tmp_path, sets=V1, cols=GOLDEN["columns"])
    want = GOLDEN["gonad"]
    for k in C.TABLE_KEYS:
        assert got[k + "_columns"] == want[k + "_columns"], k
    _check(got, want, "end-to-end gonad")
