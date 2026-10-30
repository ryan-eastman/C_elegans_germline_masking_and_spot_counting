"""sc_straighten v3 with its new switches off reproduces v2 (protocol changelog entries 32 and 33): every declared
column of every table, the strips, the strip index and the class counts of registered and development synthetic
cases, at the defaults and with the v2 punctum switches of entries 16 and 18, equal those of the v2 head (commit
acbee94; tests/golden/sc_straighten/v2_compat.json, made once before any change of the v3 round).

As in tests/test_sc_straighten_v2_compat.py, the digests are workstation pins (a difference fails with
GERMQUANT_WORKSTATION_PINS=1 and skips elsewhere, saying what differs), while the number of strands of a case is
compared on every platform (protocol changelog entries 27 and 30).
"""
from __future__ import annotations

import json
import os

import pytest
import sc_v2_compat as C

from germquant.sc.straighten import params

GOLDEN = json.loads(C.golden_path().read_text(encoding="ascii"))
ENFORCE = os.environ.get("GERMQUANT_WORKSTATION_PINS") == "1"
STRUCT = ("n_strands",)
BY_NAME = {C.case_name(c): c for c in C.CASES}


def _check(got: dict, want: dict, what: str) -> None:
    hard = [k for k in STRUCT if k in want and got.get(k) != want[k]]
    if hard:
        pytest.fail(f"{what}: structural difference from v2 in {hard} (number of strands {got.get('n_strands')} "
                    f"against {want.get('n_strands')})")
    diffs = sorted(k for k in want if got.get(k) != want[k])
    if not diffs:
        return
    msg = f"{what}: differs from the v2 digests of the workstation in {diffs}"
    if ENFORCE:
        pytest.fail(msg)
    pytest.skip(msg + " (enforced with GERMQUANT_WORKSTATION_PINS=1)")


def test_the_golden_covers_every_case():
    assert sorted(GOLDEN["cases"]) == sorted(BY_NAME)
    assert GOLDEN["made_with"]["commit"].startswith("acbee94")


def test_the_v3_switches_are_off_by_default_and_have_v1_values():
    d = params.resolve()
    for key in ("gaps.dim_stretch", "punctum.one_gaussian"):
        assert d[key] is False and params.V1_VALUES[key] is False, key
    for key in ("punctum.sigma_source", "gaps.dim_sigma_source"):
        assert d[key] == "psf" and params.V1_VALUES[key] == "psf", key
    assert d["punctum.enabled"] is False
    # the configuration of entries 32 and 33 is reachable by its switches alone
    p3 = params.resolve(None, params.v3_sets())
    for key, value in params.V3_VALUES.items():
        assert p3[key] == value, key


@pytest.mark.parametrize("name", sorted(BY_NAME))
def test_tables_equal_v2(name):
    got = C.case_digests(BY_NAME[name], GOLDEN["columns"])
    _check(got, GOLDEN["cases"][name], name)
