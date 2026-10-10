"""The acquisition record (io/acquisition.py): the full per-channel parse of the .nd2 description, its
image_summary fields, the acquisition.expected check and the `germquant info` table.

The description strings below follow the real Nikon format of the lab's Kinetix / X-Light / Celesta files
(one "Plane #n:" block per channel with camera settings, filter wheels and the laser line), including the
647 lamin plane of the ccw77 / ccw80 panel.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from test_envelope import patched_reader_lamin  # noqa: F401  (fixture: synthetic 4-colour gonad with a lamin ring)

from germquant import cli, pipeline, schema
from germquant.config import load_config
from germquant.io.acquisition import (
    PLANE_FIELDS,
    acquisition_fields,
    expected_flags,
    expected_settings,
    format_table,
    parse_description,
    parse_expected,
    read_acquisition,
    role_records,
)

ROOT = Path(__file__).resolve().parents[1]

_HEADER = ("Metadata:\nDimensions: λ(4) x Z(44)\nCamera Name: Kinetix, A25D723034\nNumerical Aperture: 1.42\n"
           "Refractive Index: 1.515\nNumber of Picture Planes: 4\n")
_EX = ("MXR00543 - CELESTA-DA/FI/TR/Cy5/Cy7-A (DAPI / FITC / TRITC / Cy5 / Cy7 - Full Multiband Penta) "
       "(FF01-391/477/549/639/741) (Ex)")
_DM = ("MXR00543 - CELESTA-DA/FI/TR/Cy5/Cy7-A (DAPI / FITC / TRITC / Cy5 / Cy7 - Full Multiband Penta) "
       "(FF421/491/567/659/776-Di01) (Di)")


def _plane(n, name, exposure, em_pos, em_label, line, power, camera_mode="12-bit (CMS)"):
    """One plane block as the Nikon software writes it (the camera settings start on the same line as
    their header, the temperature line carries a degree sign, the laser line names the channel)."""
    return (f"Plane #{n}:\n Name: {name}\n Component Count: 1\n Relay Lens Magnification: 1x\n"
            " Modality: Widefield Fluorescence, Spinning Disk Confocal\n"
            f" Camera Settings:   Binning: 1x1\n  Bit Depth: 12-bit\n  Exposure: {exposure} ms\n"
            f"  Readout Speed: 100.0MHz\n  Multiplier: Off\n  Camera Mode: {camera_mode}\n  Temperature: 0.0°C\n  \n"
            "  Conversion Gain: Sensitivity\n  Denoise.ai OFF\n  Clarify.ai OFF\n  Denoising: Off\n  Fan Speed: High\n\n"
            " Microscope Settings:   Microscope: Ti2 Microscope\n  : 4\n"
            "  Nikon Ti2, FilterChanger(Turret-Lo): 1 (Empty)\n  LightPath: L100\n  PFS, offset: 7733\n"
            f"  X-Light, FilterChanger(DM Wheel): 1 ({_DM})\n"
            f"  X-Light, FilterChanger(EM Wheel): {em_pos} ({em_label})\n"
            f"  X-Light, FilterChanger(EX Wheel): 1 ({_EX})\n"
            "  X-Light, FilterChanger(DualCam Splitter): 1\n  Illumination Iris: 100 \n  Emission Iris: 100 \n"
            "  X-Light Switcher:\n     Port 2: 50 µm\n  \n"
            "  SpectraIII/Celesta/Ziva, Shutter(Celesta): Active\n  SpectraIII/Celesta/Ziva, MultiLaser(Celesta):\n"
            f"     Line:{line}; ExW:{name}; Power: {power}; On\n\n")


def description(*, syp_exposure=200, syp_power=96.0, syp_em=(3, "595/31")):
    """The 4-colour ccw77 / ccw80 panel: DAPI 405 / PGL-1 477 / SYP-3 545 / lamin 640 (the 647 stain)."""
    return (_HEADER + _plane(1, "405", 7, 1, "438/24", 1, 49.0) + _plane(2, "477", 50, 2, "511/20", 3, 91.0)
            + _plane(3, "545", syp_exposure, *syp_em, 5, syp_power) + _plane(4, "640", 30, 4, "685/40", 6, 67.0)
            + "Z Stack Loop: 44\n- Step: 0.2 µm\n")


GOOD = description()
# the session the lab found: a dead mCherry::SYP-3 imaged at 97 percent power and 300 ms through the far-red
# emission filter, so the 647 lamin stain could reach the 545 channel
BAD = description(syp_exposure=300, syp_power=97.0, syp_em=(4, "685/40"))
CHANNELS = ["405", "477", "545", "640"]
ROLES = {"dna": 0, "granule": 1, "central_element": 2, "lamin": 3}
EXPECTED = {
    "dna": {"em_filter": "438/24"},
    "granule": {"em_filter": ["511/20", "525/30"], "exposure_ms_max": 100},
    "central_element": {"em_filter": "595/31", "exposure_ms_max": 250, "power_max": 96.5, "camera_mode": "12-bit (CMS)"},
    "lamin": {"em_filter": "685/40"},
}


def test_parse_description_reads_every_plane_field():
    per = parse_description(GOOD)
    assert list(per) == CHANNELS                                   # "Camera Name:" in the header is not a plane
    for rec in per.values():
        assert tuple(rec) == PLANE_FIELDS
    syp = per["545"]
    assert syp["exposure_ms"] == 200.0 and syp["power"] == 96.0    # the two original keys, unchanged
    assert syp["em_filter"] == "595/31" and syp["camera_mode"] == "12-bit (CMS)"
    assert syp["bit_depth"] == 12 and isinstance(syp["bit_depth"], int)
    assert syp["binning"] == "1x1" and syp["conversion_gain"] == "Sensitivity" and syp["readout_speed"] == "100.0MHz"
    assert syp["modality"] == "Widefield Fluorescence, Spinning Disk Confocal"
    assert syp["ex_filter"] == _EX and syp["dichroic"] == _DM     # the wheel position is dropped, the label kept whole
    assert [per[c]["em_filter"] for c in CHANNELS] == ["438/24", "511/20", "595/31", "685/40"]
    assert per["640"]["exposure_ms"] == 30.0 and per["640"]["power"] == 67.0 and per["405"]["exposure_ms"] == 7.0
    bad = parse_description(BAD)["545"]
    assert bad["em_filter"] == "685/40" and bad["exposure_ms"] == 300.0 and bad["power"] == 97.0


def test_parse_description_minimal_and_odd_blocks():
    desc = ("Plane #1:\n Name: 405\n Exposure: 50 ms\n ExW:405; Power: 12.5\n"
            "Plane #2:\n Name: DAPI 405\n Exposure: 1.5 s\n Line:1; ExW:405; Power: 30; On\n"
            "Plane #3:\n Component Count: 1\n")
    per = parse_description(desc)
    assert per["405"]["exposure_ms"] == 50.0 and per["405"]["power"] == 12.5
    assert all(per["405"][k] is None for k in PLANE_FIELDS if k not in ("exposure_ms", "power"))
    # a channel name with a space, an exposure in seconds; a laser line not named like the channel gives
    # no power (NaN downstream is the signal that the names disagree, the semantics the goldens carry)
    assert per["DAPI 405"]["exposure_ms"] == 1500.0 and per["DAPI 405"]["power"] is None
    assert len(per) == 2                                           # a block without a name is skipped
    assert parse_description("") == {} and parse_description(None) == {}


def test_fields_for_a_four_channel_role_map(monkeypatch):
    monkeypatch.setattr("germquant.io.acquisition.read_description", lambda path: GOOD)
    f = acquisition_fields("x.nd2", CHANNELS, {**ROLES, "foci": None})
    # the original two columns per role come first, in the order they always had
    assert list(f)[:8] == ["exp_ms_dna", "laser_pow_dna", "exp_ms_granule", "laser_pow_granule",
                           "exp_ms_central_element", "laser_pow_central_element", "exp_ms_lamin", "laser_pow_lamin"]
    assert f["exp_ms_central_element"] == 200.0 and f["laser_pow_central_element"] == 96.0 and f["exp_ms_lamin"] == 30.0
    assert f["em_filter_central_element"] == "595/31" and f["em_filter_lamin"] == "685/40"
    assert f["camera_mode_dna"] == "12-bit (CMS)" and f["bit_depth_granule"] == 12
    assert "exp_ms_foci" not in f and "em_filter_foci" not in f    # an unresolved role has no columns
    assert list(f)[-1] == "acquisition_json"
    per = json.loads(f["acquisition_json"])
    assert per == parse_description(GOOD) and list(per) == CHANNELS
    assert "\n" not in f["acquisition_json"] and ": " not in f["acquisition_json"]   # compact JSON


def test_fields_when_a_channel_is_missing_or_the_description_is_unreadable(monkeypatch):
    two = _HEADER + _plane(1, "405", 7, 1, "438/24", 1, 49.0) + _plane(4, "640", 30, 4, "685/40", 6, 67.0)
    monkeypatch.setattr("germquant.io.acquisition.read_description", lambda path: two)
    rec = read_acquisition("x.nd2", CHANNELS, ROLES)
    assert rec.per_role["central_element"] == {} and rec.per_role["lamin"]["em_filter"] == "685/40"
    f = rec.fields
    assert np.isnan(f["exp_ms_central_element"]) and np.isnan(f["bit_depth_central_element"])
    assert f["em_filter_central_element"] is None and f["camera_mode_central_element"] is None
    assert f["exp_ms_dna"] == 7.0 and f["em_filter_dna"] == "438/24"

    def boom(path):
        raise OSError("no such file")

    monkeypatch.setattr("germquant.io.acquisition.read_description", boom)
    f = acquisition_fields("x.nd2", CHANNELS, ROLES)              # never fails the run: every field NA
    assert np.isnan(f["exp_ms_dna"]) and f["em_filter_dna"] is None and f["acquisition_json"] == "{}"


def test_role_records_match_by_name_then_substring():
    per = parse_description(GOOD)
    out = role_records(per, ["Widefield 405", "477", "545 nm", "640"],
                       {"dna": 0, "central_element": 2, "lamin": 3, "foci": None, "beyond": 7})
    assert out["dna"]["em_filter"] == "438/24" and out["central_element"]["em_filter"] == "595/31"
    assert out["lamin"]["exposure_ms"] == 30.0 and "foci" not in out
    assert out["beyond"] == {}                                     # an index past the channel list matches nothing


def test_expected_check_fires_on_the_bad_session_and_not_on_the_good_one():
    exp = parse_expected(EXPECTED, roles=set(ROLES))
    good = role_records(parse_description(GOOD), CHANNELS, ROLES)
    assert expected_flags(good, exp) == []
    bad = role_records(parse_description(BAD), CHANNELS, ROLES)
    assert expected_flags(bad, exp) == ["acq:central_element_em_filter_mismatch",
                                        "acq:central_element_exposure_above_max",
                                        "acq:central_element_power_above_max"]
    # strings match after removing whitespace and ignoring case; a list gives several accepted values
    exp2 = parse_expected({"central_element": {"camera_mode": ["12-BIT(cms)", "16-bit (Sensitivity)"],
                                               "em_filter": " 595 / 31 "}})
    assert expected_flags(good, exp2) == []
    exp3 = parse_expected({"central_element": {"camera_mode": "16-bit (Sensitivity)"}})
    assert expected_flags(good, exp3) == ["acq:central_element_camera_mode_mismatch"]
    # a role that did not resolve in this image is not checked
    assert expected_flags({"dna": good["dna"]}, exp) == []
    # a value the description does not carry: the string checks report it (the expectation could not be
    # confirmed), the maxima stay silent (nothing is above a maximum)
    exp4 = parse_expected({"lamin": {"em_filter": "685/40", "exposure_ms_max": 100, "power_max": 50,
                                     "camera_mode": "12-bit (CMS)"}})
    assert expected_flags({"lamin": {}}, exp4) == ["acq:lamin_em_filter_mismatch", "acq:lamin_camera_mode_mismatch"]
    # equal to the maximum is not above it
    exp5 = parse_expected({"central_element": {"exposure_ms_max": 200, "power_max": 96}})
    assert expected_flags(good, exp5) == []


def test_parse_expected_validates_the_section():
    assert parse_expected(None) == {}
    out = parse_expected({"dna": {"em_filter": 438, "exposure_ms_max": 10}}, roles={"dna"})
    assert out == {"dna": {"em_filter": ["438"], "exposure_ms_max": 10.0}}
    with pytest.raises(ValueError, match="not a role"):
        parse_expected({"central element": {"em_filter": "595/31"}}, roles={"central_element"})
    with pytest.raises(ValueError, match="unknown key"):
        parse_expected({"dna": {"exposure_max": 10}})
    with pytest.raises(ValueError, match="must be a number"):
        parse_expected({"dna": {"exposure_ms_max": "10 ms"}})
    with pytest.raises(ValueError, match="string or a list"):
        parse_expected({"dna": {"em_filter": [None]}})
    with pytest.raises(ValueError, match="mapping"):
        parse_expected({"dna": None})
    with pytest.raises(ValueError, match="mapping"):
        parse_expected(["dna"])


def test_expected_section_is_absent_from_every_shipped_config():
    """The acq:<role>_* flags can only fire when a profile sets acquisition.expected: no shipped config
    or profile does, so the qc_flags of the default path (and of every golden) cannot change."""
    files = sorted((ROOT / "config").glob("*.yaml")) + sorted((ROOT / "config" / "profiles").glob("*.yaml"))
    assert files
    for p in files:
        cfg = load_config(p)
        assert cfg.get("acquisition.expected") is None, p
        assert expected_settings(cfg) == {}, p
    cfg = load_config(ROOT / "config" / "config.yaml")
    cfg.set("acquisition.expected", {"foci": {"power_max": 50}})
    assert expected_settings(cfg) == {"foci": {"power_max": 50.0}}
    cfg.set("acquisition.expected", {"lamin": {"em_filter": "685/40"}})   # the N2 channel map has no lamin role
    with pytest.raises(ValueError, match="not a role"):
        expected_settings(cfg)
    assert "acquisition_json" in schema.IMAGE_SUMMARY             # the one declared (fixed-name) column


def _cfg(**switches):
    cfg = load_config("config/config_ccw77.yaml")
    cfg._data["segmentation"]["nuclei"]["method"] = "classical"
    cfg.set("coloc.enabled", False)
    cfg.set("acquisition.read_exposures", True)
    for k, v in switches.items():
        cfg.set(k, v)
    return cfg


NAME = "20260708_ccw77_IF_pgl1_syp3_lmn1_HS_male_008.nd2"


def test_pipeline_records_and_checks_the_acquisition(tmp_path, patched_reader_lamin, monkeypatch):  # noqa: F811
    monkeypatch.setattr("germquant.io.acquisition.read_description", lambda path: BAD)
    # without acquisition.expected (the code default) the record is written and no acq: flag can fire
    res = pipeline.process_image(NAME, _cfg(), tmp_path / "plain", prov=None)
    assert res["stages"]["acquisition"]["status"] == "ran"
    assert not [f for f in res["qc_flags"] if f.startswith("acq:")]
    isum = pd.read_csv(next((tmp_path / "plain").rglob("*__image_summary.csv"))).iloc[0]
    assert isum["exp_ms_central_element"] == 300.0 and isum["laser_pow_central_element"] == 97.0
    assert isum["em_filter_central_element"] == "685/40" and isum["em_filter_lamin"] == "685/40"
    assert isum["camera_mode_dna"] == "12-bit (CMS)" and isum["bit_depth_dna"] == 12
    assert json.loads(isum["acquisition_json"]) == parse_description(BAD)
    cols = list(isum.index)
    assert cols.index("acquisition_json") < cols.index("exp_ms_dna") < cols.index("em_filter_dna")   # declared, then extras

    # with the pack's expected settings the bad session is flagged; advisory only, qc_pass is unchanged
    res = pipeline.process_image(NAME, _cfg(**{"acquisition.expected": EXPECTED}), tmp_path / "checked", prov=None)
    acq = [f for f in res["qc_flags"] if f.startswith("acq:")]
    assert acq == ["acq:central_element_em_filter_mismatch", "acq:central_element_exposure_above_max",
                   "acq:central_element_power_above_max"]
    assert res["stages"]["acquisition"]["flags"] == acq and res["qc_pass"]
    isum = pd.read_csv(next((tmp_path / "checked").rglob("*__image_summary.csv"))).iloc[0]
    assert set(acq) <= set(isum["qc_flags"].split(";"))

    # the good session passes the same check
    monkeypatch.setattr("germquant.io.acquisition.read_description", lambda path: GOOD)
    res = pipeline.process_image(NAME, _cfg(**{"acquisition.expected": EXPECTED}), tmp_path / "good", prov=None)
    assert not [f for f in res["qc_flags"] if f.startswith("acq:")]

    # a typo in the section fails before any pixel is read
    with pytest.raises(ValueError, match="unknown key"):
        pipeline.process_image(NAME, _cfg(**{"acquisition.expected": {"lamin": {"exposure_max": 10}}}),
                               tmp_path / "typo", prov=None)


def test_info_prints_the_acquisition_table(monkeypatch, capsys):
    meta = {"sizes": {"C": 4, "Z": 44, "Y": 2048, "X": 2048}, "spacing": (0.2, 0.065, 0.065), "spacing_ok": True,
            "channel_names": CHANNELS, "dtype": "uint16", "is_2d": False}
    monkeypatch.setattr(cli, "read_image_metadata", lambda path: meta)
    monkeypatch.setattr("germquant.io.acquisition.read_description", lambda path: GOOD)
    assert cli.main(["info", "x.nd2"]) == 0
    out = capsys.readouterr().out
    assert "channels   : ['405', '477', '545', '640']" in out
    assert format_table(parse_description(GOOD)) in out
    row = [ln for ln in out.splitlines() if ln.strip().startswith("545")]
    assert row and row[0].split() == ["545", "200", "96", "595/31", "12-bit", "(CMS)"]
    assert "NA" not in out

    # a value the description lacks prints as NA
    monkeypatch.setattr("germquant.io.acquisition.read_description",
                        lambda path: "Plane #1:\n Name: 405\n Exposure: 50 ms\n ExW:405; Power: 12.5\n")
    assert cli.main(["info", "x.nd2"]) == 0
    row = [ln for ln in capsys.readouterr().out.splitlines() if ln.strip().startswith("405")]
    assert row and row[0].split() == ["405", "50", "12.5", "NA", "NA"]

    # not an .nd2: the header prints and no description is read (other formats carry none)
    def never(path):
        raise AssertionError("the description must not be read for a TIFF")

    monkeypatch.setattr("germquant.io.acquisition.read_description", never)
    assert cli.main(["info", "x.ome.tif"]) == 0
    out = capsys.readouterr().out
    assert "channels" in out and "acquisition" not in out

    # an unreadable description is reported, never raised
    def boom(path):
        raise OSError("truncated header")

    monkeypatch.setattr("germquant.io.acquisition.read_description", boom)
    assert cli.main(["info", "x.nd2"]) == 0
    assert "acquisition: description unreadable (truncated header)" in capsys.readouterr().out
