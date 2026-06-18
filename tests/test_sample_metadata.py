from germquant.io.sample_metadata import parse_sample

REGEX = r"(?P<date>\d{6,8})_(?P<genotype>[A-Za-z0-9]+)_(?P<treatment>nohs|hs)_(?:[a-z0-9]+_)*?(?P<sex>HERM|MALE|H|M)\s*_?(?P<replicate>\d+)?\s*$"


def test_parse_herm_nohs():
    f = parse_sample("20251105_n2_nohs_HERM _001.nd2", REGEX)
    assert f["genotype"] == "N2"
    assert f["treatment"] == "control"
    assert f["sex"] == "herm"
    assert f["germ_cell"] == "oocyte"
    assert f["replicate"] == "001"
    assert f["acquisition_date"] == "20251105"


def test_parse_male_hs():
    f = parse_sample("20251021_n2_hs_MALE_2.nd2", REGEX)
    assert f["treatment"] == "heat"
    assert f["sex"] == "male"
    assert f["germ_cell"] == "spermatocyte"


def test_parse_new_naming_channel_block_and_h_m_codes():
    # NAS data: channel block embedded, sex coded H/M
    f = parse_sample("20260521_n2_hs_dapi_rad51_syp1_M001.nd2", REGEX)
    assert f["treatment"] == "heat"
    assert f["sex"] == "male"          # M canonicalized to male
    assert f["germ_cell"] == "spermatocyte"
    assert f["replicate"] == "001"

    g = parse_sample("20260529_n2_nohs_dapi_rad51_syp1_HERM.nd2", REGEX)
    assert g["treatment"] == "control"
    assert g["sex"] == "herm"
    assert g["germ_cell"] == "oocyte"
    assert g["replicate"] == "NA"      # no replicate suffix -> default
