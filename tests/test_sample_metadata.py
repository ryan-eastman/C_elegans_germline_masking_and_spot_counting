from germquant.io.sample_metadata import expected_sc_count, parse_sample

REGEX = r"(?P<date>\d{6,8})_(?P<genotype>[A-Za-z0-9]+)_(?P<treatment>nohs|hs)_(?P<sex>HERM|MALE)\s*_?(?P<replicate>\d+)?"


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


def test_expected_sc_counts():
    assert expected_sc_count("oocyte") == 6      # 5 autosomal + XX
    assert expected_sc_count("spermatocyte") == 5  # X univalent
    assert expected_sc_count("unknown") is None
