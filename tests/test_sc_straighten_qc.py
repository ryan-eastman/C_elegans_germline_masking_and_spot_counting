"""The sc_straighten QC vocabulary: one positive and one negative case per code of qc.CODES, the status levels,
and the flag table of docs/SC_STRAIGHTEN.md kept equal to the code."""
from __future__ import annotations

from pathlib import Path

import pytest

from germquant.sc.straighten import params, qc

P = params.resolve()
ROOT = Path(__file__).resolve().parents[1]

CLEAN_STRAND = {
    "length_um": 5.0, "snr_core_median": 20.0, "end0_class": "free_fit", "end1_class": "free_fit",
    "near_stack_edge": False, "link_cost_max": float("nan"), "frac_contaminated": 0.0, "merged_run_um": 0.0,
    "fwhm_lat_um": 0.4, "fwhm_ax_um": 0.8, "frac_in_label": 1.0, "frac_near_other_label": 0.0,
    "frac_steep60": 0.0, "kink_max_deg": 10.0, "fit_rms_um": 0.03, "gap_max_um": 0.0, "frac_saturated": 0.0,
    "geometry_source": "auto", "verified": False, "error": "",
}
CTX = {"fwhm_lat_med": 0.4, "fwhm_ax_med": 0.8}
STRAND_POSITIVE = {
    "fragment": {"length_um": 0.5},
    "low_snr": {"snr_core_median": 3.0},
    "low_snr_fail": {"snr_core_median": 2.0},
    "end_unfit": {"end0_class": "free_unfit"},
    "contact_end": {"end1_class": "contact"},
    "truncated_z": {"end0_class": "stack_edge"},
    "costly_link": {"link_cost_max": 70.0},
    "contaminated": {"frac_contaminated": 0.5},
    "merged_suspect": {"merged_run_um": 0.6},
    "double_width": {"fwhm_lat_um": 0.9},
    "exits_label": {"frac_in_label": 0.5},
    "neighbour_nucleus": {"frac_near_other_label": 0.6},
    "steep_z": {"frac_steep60": 0.7},
    "long": {"length_um": 12.0},
    "kink": {"kink_max_deg": 75.0},
    "fit_residual": {"fit_rms_um": 0.2},
    "gap_long": {"gap_max_um": 0.7},
    "saturated": {"frac_saturated": 0.05},
    "imported": {"geometry_source": "imported"},
    "verified": {"verified": True},
    "error": {"error": "ValueError"},
}
CLEAN_NUCLEUS = {"volume_um3": 30.0, "touches_z": False, "touches_xy": False, "count_status": "ok", "noise_sd": 1.0,
                 "budget_exceeded": False, "error": ""}
NUCLEUS_POSITIVE = {
    "label_merged_suspect": {"volume_um3": 60.0},
    "truncated_z": {"touches_z": True},
    "edge_xy": {"touches_xy": True},
    "count_low": {"count_status": "low"},
    "count_high": {"count_status": "high"},
    "noise_mismatch": {"noise_sd": 1.6},
    "budget_exceeded": {"budget_exceeded": True},
    "error": {"error": "RuntimeError"},
}
CLEAN_IMAGE = {"noise_rel_rms": 0.05, "noise_fallback": False, "psf_source": "estimated", "n_present": 100,
               "n_present_with_2um": 80, "spot_tables": [], "geometry_source": "done_marker", "n_label_merged": 0,
               "n_failed": 0}
IMAGE_POSITIVE = {
    "sc_straighten:uncalibrated": {},
    "noise_fit_poor": {"noise_rel_rms": 0.3},
    "psf_fallback": {"psf_source": "fallback"},
    "psf_from_foci": {"psf_source": "foci"},
    "sc_channel_suspect": {"n_present_with_2um": 10},
    "foci_uncalibrated": {"spot_tables": ["spots"]},
    "geometry_from_manifest": {"geometry_source": "run_manifest"},
    "geometry_from_golden_manifest": {"geometry_source": "golden_manifest"},
    "label_merged_suspect": {"n_label_merged": 3},
    "nucleus_errors": {"n_failed": 2},
}


def _has(flags, code):
    return any(f == code or f.startswith(code + ":") or f.startswith(code + "=") for f in flags)


@pytest.mark.parametrize("code,scope", [(c, s) for c, s, *_ in qc.CODES])
def test_every_code_has_a_positive_and_a_negative_case(code, scope):
    if scope == "strand":
        pos = qc.strand_flags({**CLEAN_STRAND, **STRAND_POSITIVE[code]}, P, CTX)
        neg = qc.strand_flags(CLEAN_STRAND, P, CTX)
    elif scope == "nucleus":
        pos = qc.nucleus_flags({**CLEAN_NUCLEUS, **NUCLEUS_POSITIVE[code]}, P, 30.0)
        neg = qc.nucleus_flags(CLEAN_NUCLEUS, P, 30.0)
    else:
        pos = qc.image_flags({**CLEAN_IMAGE, **IMAGE_POSITIVE[code]}, P)
        neg = qc.image_flags(CLEAN_IMAGE, P)
    assert _has(pos, code), (code, pos)
    if code == "sc_straighten:uncalibrated":
        assert _has(neg, code)                    # always present in v1
    else:
        assert not _has(neg, code), (code, neg)


def test_levels_and_status():
    assert qc.status([], "strand") == "pass"
    assert qc.status(["contact_end"], "strand") == "pass"
    assert qc.status(["contact_end", "low_snr"], "strand") == "review"
    assert qc.status(["low_snr", "low_snr_fail"], "strand") == "fail"
    assert qc.status(["error:ValueError"], "strand") == "fail"
    assert qc.status(["truncated_z"], "nucleus") == "fail"
    assert qc.level("label_merged_suspect=3", "image") == "info"
    # low_snr_fail replaces low_snr
    f = qc.strand_flags({**CLEAN_STRAND, "snr_core_median": 1.0}, P, CTX)
    assert "low_snr_fail" in f and "low_snr" not in f


def test_codes_are_unique_per_scope():
    keys = [(c, s) for c, s, *_ in qc.CODES]
    assert len(keys) == len(set(keys))
    for c, s, lvl, *_ in qc.CODES:
        assert lvl in ("info", "review", "fail")


def test_docs_flag_table_matches_the_code():
    doc = (ROOT / "docs" / "SC_STRAIGHTEN.md").read_text(encoding="utf-8")
    assert qc.docs_table() in doc.replace("\r\n", "\n")


def test_noise_flag_uses_the_signal_range_fit():
    # the fit above the darkest decile decides; the all-bins relative RMS is only a fallback
    assert "noise_fit_poor" not in qc.image_flags({**CLEAN_IMAGE, "noise_rel_rms": 0.3,
                                                   "noise_rel_rms_signal": 0.1}, P)
    assert "noise_fit_poor" in qc.image_flags({**CLEAN_IMAGE, "noise_rel_rms": 0.1, "noise_rel_rms_signal": 0.3}, P)


def test_label_volume_reference_follows_the_gonad_axis():
    vols = [20.0] * 10 + [60.0] * 10                   # small nuclei distally, large ones proximally
    axs = list(range(0, 10)) + list(range(100, 110))
    assert qc.local_volume_reference(vols, axs, 105.0, 15.0) == 60.0
    assert qc.local_volume_reference(vols, axs, 4.0, 15.0) == 20.0
    assert qc.local_volume_reference(vols, axs, 50.0, 15.0) == 40.0       # nobody near: the whole median
    assert qc.local_volume_reference(vols, axs, float("nan"), 15.0) == 40.0
    assert qc.local_volume_reference([], [], 1.0, 15.0) is None
    # a proximal nucleus of 70 um3 is not a merged label against its neighbours (60), against the whole median
    # (40) it would be
    ref = qc.local_volume_reference(vols, axs, 105.0, 15.0)
    assert "label_merged_suspect" not in qc.nucleus_flags({**CLEAN_NUCLEUS, "volume_um3": 70.0}, P, ref)
    assert "label_merged_suspect" in qc.nucleus_flags({**CLEAN_NUCLEUS, "volume_um3": 70.0}, P, 40.0)


def test_docs_column_tables_match_the_schema():
    from germquant.sc.straighten import schema as S

    doc = (ROOT / "docs" / "SC_STRAIGHTEN.md").read_text(encoding="utf-8").replace("\r\n", "\n")
    for table, cols in S.TABLES.items():
        assert set(S.COLUMN_DOCS[table]) == set(cols), table
        assert all(S.COLUMN_DOCS[table][c].strip() for c in cols), table
        assert S.docs_columns_table(table) in doc, table
