"""Gradient-window pachytene refinement: keep mature-synapsis nuclei, drop the RAD-51-high /
SYP-dense early-meiotic region, independent of the SC fragment count being restricted.
"""
import numpy as np
import pandas as pd

from germquant.zones.calling import refine_pachytene


def _gonad():
    """150 nuclei along the axis: a distal early/TZ region (high RAD-51 + dense SYP, little SC)
    then a pachytene region (low RAD-51, SC present)."""
    pos = np.linspace(0, 1, 150)
    early = pos < 0.30
    return pd.DataFrame({
        "nucleus_id": np.arange(1, 151),
        "axis_position_norm": pos,
        "zone_call": "pachytene",                       # pretend the crescent step over-called
        "n_foci": np.where(early, 5.0, 1.0),            # RAD-51 high distally
        "central_element_mean_intensity": np.where(early, 500.0, 150.0),  # SYP dense distally
        "sc_total_length_um": np.where(early, 1.0, 8.0),
        "sc_n_fragments": np.where(early, 1, 5),
    })


def test_pachytene_excludes_early_keeps_mature():
    df, flags = refine_pachytene(_gonad())
    # clearly-early nuclei (axis < 0.2) are excluded; clearly-pachytene (axis > 0.5) kept
    early = df[df.axis_position_norm < 0.2]
    mature = df[df.axis_position_norm > 0.5]
    assert early["is_pachytene"].mean() < 0.2
    assert mature["is_pachytene"].mean() > 0.8
    # excluded early nuclei are relabeled, not silently called pachytene
    assert (df.loc[~df.is_pachytene, "zone_call"] == "pachytene_excluded").all()
    # the gate is independent of fragment count: it used foci/SYP/axis, never sc_n_fragments
    assert 0.4 < df["is_pachytene"].mean() < 0.85


def test_pachytene_refine_skips_when_features_missing():
    df = pd.DataFrame({
        "nucleus_id": [1, 2, 3],
        "axis_position_norm": [0.1, 0.5, 0.9],
        "zone_call": ["transition_zone", "pachytene", "pachytene"],
    })  # no n_foci / central_element / sc_total_length
    out, flags = refine_pachytene(df)
    assert any("skipped_missing_features" in f for f in flags)
    # falls back to the existing zone_call
    assert out.loc[out.zone_call == "pachytene", "is_pachytene"].all()
