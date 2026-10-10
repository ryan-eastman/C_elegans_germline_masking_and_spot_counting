"""Partition measurability (docs/PARTITION_QC.md): the two image_summary columns, the opt-in flag, and the
contract that every partition number is still written when the SYP signal sits at the floor."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_envelope import patched_reader_lamin  # noqa: F401  (fixture)
from test_staging import IID

from germquant import pipeline, schema
from germquant.config import load_config
from germquant.partition import (
    LOW_SNR_FLAG,
    MEASURABILITY_FIELDS,
    SYP_SNR_FIELD,
    partition_measurability,
    syp_snr_from_summary,
)

ROOT = Path(__file__).resolve().parents[1]
SHIPPED = sorted(str(p.relative_to(ROOT)) for p in
                 list((ROOT / "config").glob("*.yaml")) + list((ROOT / "config" / "profiles").glob("*.yaml")))


# ---- the check itself ----

def test_gate_off_by_default():
    fields, fl = partition_measurability({})
    assert math.isnan(fields["partition_syp_snr"]) and fields["partition_measurable"] is True and fl == []
    fields, fl = partition_measurability({SYP_SNR_FIELD: 0.4})            # code default 0.0: never flags
    assert fields["partition_syp_snr"] == 0.4 and fields["partition_measurable"] is True and fl == []
    fields, fl = partition_measurability({SYP_SNR_FIELD: 0.4}, min_syp_snr=None)   # a null in a profile: off
    assert fields["partition_measurable"] is True and fl == []
    assert tuple(fields) == MEASURABILITY_FIELDS


def test_gate_on():
    fields, fl = partition_measurability({SYP_SNR_FIELD: 1.2}, min_syp_snr=3.0)
    assert fields == {"partition_syp_snr": 1.2, "partition_measurable": False} and fl == [LOW_SNR_FLAG]
    fields, fl = partition_measurability({SYP_SNR_FIELD: 3.0}, min_syp_snr=3.0)          # inclusive
    assert fields["partition_measurable"] is True and fl == []
    fields, fl = partition_measurability({SYP_SNR_FIELD: np.float32(8.5)}, min_syp_snr="3")   # YAML strings
    assert fields["partition_measurable"] is True and fl == []
    # an unknown ratio (signal_qc stage off: field absent, None, NaN or not a number) is measurable, unflagged
    for summary in ({}, {SYP_SNR_FIELD: None}, {SYP_SNR_FIELD: float("nan")}, {SYP_SNR_FIELD: "n/a"}):
        fields, fl = partition_measurability(summary, min_syp_snr=3.0)
        assert math.isnan(fields["partition_syp_snr"]) and fields["partition_measurable"] is True and fl == []
    assert syp_snr_from_summary({SYP_SNR_FIELD: "2.5"}) == 2.5
    with pytest.raises(ValueError):                                       # a misspelt threshold is loud, not off
        partition_measurability({SYP_SNR_FIELD: 1.0}, min_syp_snr="three")


def test_schema_appends_the_columns_and_the_writer_fills_them():
    cols = schema.IMAGE_SUMMARY
    assert cols.index("partition_syp_snr") > cols.index("tail_n_granules_scored")   # appended, never inserted
    assert cols.index("partition_measurable") == cols.index("partition_syp_snr") + 1
    # a run without the partition stage carries both columns as NA (pipeline._conform_schema, the writer)
    df = pipeline._conform_schema(pd.DataFrame([{"n_nuclei": 3}]), "image_summary")
    assert df["partition_syp_snr"].isna().all() and df["partition_measurable"].isna().all()


@pytest.mark.parametrize("path", SHIPPED)
def test_shipped_configs_leave_the_gate_off(path):
    """No shipped config or profile sets partition.min_syp_snr: the flag never fires on their goldens."""
    assert load_config(ROOT / path).get("partition.min_syp_snr") is None


# ---- end to end on the synthetic lamin gonad ----

def _study(tmp_path, name, *, acquisition=False, **partition_keys) -> Path:
    """A study profile outside the repository: the ccw77 base, the classical segmenter, the envelope frame,
    the Imaris-calibrated granules and the partition stage, plus the given partition keys."""
    text = (f"extends: {(ROOT / 'config' / 'config_ccw77.yaml').as_posix()}\n"
            "segmentation:\n  nuclei:\n    method: classical\n"
            "coloc:\n  n_random: 10\n"
            "sc:\n  enabled: false\n"
            "envelope:\n  enabled: true\n"
            f"acquisition:\n  read_exposures: {'true' if acquisition else 'false'}\n"
            "granule:\n  method: imaris_tophat\n  region: cytoplasm_shell\n"
            "partition:\n  enabled: true\n")
    text += "".join(f"  {k}: {v}\n" for k, v in partition_keys.items())
    p = tmp_path / f"{name}.yaml"
    p.write_text(text, encoding="utf-8")
    return p


@pytest.fixture
def fake_syp_snr(monkeypatch):
    """Stands in for the signal_qc stage until it lands: the acquisition stage's summary fields carry
    sig_central_element_snr, so it reaches the per-image summary dict the partition check reads. Once the
    signal_qc stage exists, inject the value through that stage instead."""
    from germquant.io import acquisition

    holder = {"snr": 1.2}
    monkeypatch.setattr(acquisition, "acquisition_fields", lambda *a, **k: {SYP_SNR_FIELD: holder["snr"]})
    return holder


def test_default_path_writes_nan_and_true(tmp_path, patched_reader_lamin):  # noqa: F811
    cfg = load_config(_study(tmp_path, "plain"))
    res = pipeline.process_image(f"{IID}.nd2", cfg, tmp_path / "r", prov=None)
    assert res["stages"]["partition"]["status"] == "ran", res["stages"]["partition"]
    s = res["tables"]["image_summary"].iloc[0]
    assert math.isnan(float(s["partition_syp_snr"])) and bool(s["partition_measurable"]) is True
    assert not any(f.startswith("partition:low") for f in res["qc_flags"]), res["qc_flags"]
    csv = pd.read_csv(next((tmp_path / "r").rglob("*__image_summary.csv")))
    assert csv["partition_syp_snr"].isna().all() and bool(csv["partition_measurable"].iloc[0]) is True


def test_low_snr_flags_and_still_writes_every_number(tmp_path, patched_reader_lamin, fake_syp_snr):  # noqa: F811
    gated = load_config(_study(tmp_path, "gated", acquisition=True, min_syp_snr=3.0))
    plain = load_config(_study(tmp_path, "plain", acquisition=True))
    assert gated.hash != plain.hash and gated.get("partition.min_syp_snr") == 3.0   # the gate enters the hash
    r_g = pipeline.process_image(f"{IID}.nd2", gated, tmp_path / "g", prov=None)
    r_p = pipeline.process_image(f"{IID}.nd2", plain, tmp_path / "p", prov=None)
    for r in (r_g, r_p):
        assert r["stages"]["partition"]["status"] == "ran" and r["stages"]["acquisition"]["status"] == "ran"
    # the flag: in qc_flags, in the image_summary column and on the partition entry of the stage record
    assert LOW_SNR_FLAG in r_g["qc_flags"] and LOW_SNR_FLAG not in r_p["qc_flags"]
    assert LOW_SNR_FLAG in r_g["stages"]["partition"]["flags"]
    sg, sp = r_g["tables"]["image_summary"].iloc[0], r_p["tables"]["image_summary"].iloc[0]
    assert LOW_SNR_FLAG in sg["qc_flags"].split(";") and LOW_SNR_FLAG not in sp["qc_flags"].split(";")
    assert sg["partition_syp_snr"] == 1.2 and bool(sg["partition_measurable"]) is False
    assert sp["partition_syp_snr"] == 1.2 and bool(sp["partition_measurable"]) is True
    # the flag, not a blank, is the signal: every partition number equals the ungated run's
    pd.testing.assert_frame_equal(r_g["tables"]["partition"].reset_index(drop=True),
                                  r_p["tables"]["partition"].reset_index(drop=True))
    assert len(r_g["tables"]["partition"]) >= 1 and r_g["tables"]["partition"]["region_voxels"].iloc[0] > 0
    numbers = [c for c in schema.IMAGE_SUMMARY if c.startswith("partition_") and c not in MEASURABILITY_FIELDS]
    for c in numbers:
        a, b = sg.get(c), sp.get(c)
        assert (pd.isna(a) and pd.isna(b)) or a == b, (c, a, b)
    assert pd.notna(sg["partition_bg"])
    csv = pd.read_csv(next((tmp_path / "g").rglob("*__image_summary.csv")))
    assert csv["partition_syp_snr"].iloc[0] == 1.2 and bool(csv["partition_measurable"].iloc[0]) is False
    assert LOW_SNR_FLAG in str(csv["qc_flags"].iloc[0]).split(";")
    # above the gate: measurable, no flag
    fake_syp_snr["snr"] = 5.0
    r_ok = pipeline.process_image(f"{IID}.nd2", gated, tmp_path / "ok", prov=None)
    s_ok = r_ok["tables"]["image_summary"].iloc[0]
    assert LOW_SNR_FLAG not in r_ok["qc_flags"] and s_ok["partition_syp_snr"] == 5.0
    assert bool(s_ok["partition_measurable"]) is True
