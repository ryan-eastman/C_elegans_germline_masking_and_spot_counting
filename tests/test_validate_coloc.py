"""The Monday validation harness works on synthetic stand-ins for the Imaris exports."""
import pandas as pd

from germquant.validate import compare_coloc_metrics, compare_granules


def test_compare_granules_matches_and_volume():
    our = pd.DataFrame({"z_um": [1, 2, 3], "y_um": [1, 1, 1], "x_um": [1, 2, 3],
                        "volume_um3": [1.0, 2.0, 3.0]})
    ref = pd.DataFrame({"z_um": [1, 2], "y_um": [1, 1], "x_um": [1.1, 2.1],
                        "volume_um3": [1.1, 2.2]})

    out = compare_granules(our, ref, max_match_um=0.5)

    assert out["n_our"] == 3 and out["n_ref"] == 2
    assert out["n_matched"] == 2
    assert abs(out["recall_of_ref"] - 1.0) < 1e-9
    assert abs(out["precision"] - 2 / 3) < 1e-9
    assert out["volume_stats"]["n"] == 2
    assert out["mean_match_dist_um"] < 0.2


def test_compare_granules_empty_ref():
    empty = pd.DataFrame(columns=["z_um", "y_um", "x_um", "volume_um3"])
    ref = pd.DataFrame({"z_um": [1.0], "y_um": [1.0], "x_um": [1.0], "volume_um3": [1.0]})
    out = compare_granules(empty, ref)
    assert out["n_our"] == 0 and out["n_matched"] == 0


def test_compare_coloc_metrics():
    coloc = pd.DataFrame([
        {"sc_operand": "syp_aggregate", "manders_m1": 0.8, "manders_m2": 0.5, "pearson_r": 0.3},
        {"sc_operand": "sc_ribbon", "manders_m1": 0.05, "manders_m2": 0.02, "pearson_r": 0.01},
    ])
    imaris = {"manders_m1": 0.75, "manders_m2": 0.55, "pearson_r": 0.28}

    df = compare_coloc_metrics(coloc, imaris, operand="syp_aggregate")

    assert set(df["metric"]) == {"manders_m1", "manders_m2", "pearson_r"}
    m1 = df[df["metric"] == "manders_m1"].iloc[0]
    assert abs(m1["ours"] - 0.8) < 1e-9 and abs(m1["imaris"] - 0.75) < 1e-9
    assert abs(m1["abs_diff"] - 0.05) < 1e-9


def test_compare_coloc_metrics_missing_imaris():
    coloc = pd.DataFrame([{"sc_operand": "syp_aggregate", "manders_m1": 0.8,
                           "manders_m2": 0.5, "pearson_r": 0.3}])
    df = compare_coloc_metrics(coloc, None, operand="syp_aggregate")
    assert df["imaris"].isna().all()          # no Imaris numbers -> NaN, ours still populated
    assert not df["ours"].isna().any()
