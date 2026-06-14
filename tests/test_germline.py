"""Germline isolation (germquant.germline.select). Selection uses SYP + geometry only — never
n_foci — so these tests never reference foci."""
import pandas as pd

from germquant.germline import select_germline


def _table(rows):
    df = pd.DataFrame(rows)
    df.insert(0, "nucleus_id", range(1, len(df) + 1))
    return df


def _blob(cz, cy, cx0, n, step, syp):
    """A connected chain of n nuclei (centroids step µm apart) with a given SYP level."""
    return [{"centroid_z_um": cz, "centroid_y_um": cy, "centroid_x_um": cx0 + i * step,
             "central_element_mean_intensity": syp} for i in range(n)]


def test_drops_offgonad_junk():
    # a 20-nucleus SYP-bright germline chain + 4 far, dim, scattered junk nuclei
    germ = _blob(0, 0, 0, 20, 4.0, 100.0)
    junk = [{"centroid_z_um": 0, "centroid_y_um": 200 + 30 * i, "centroid_x_um": 200,
             "central_element_mean_intensity": 5.0} for i in range(4)]
    df = _table(germ + junk)
    out, flags = select_germline(df, method="multi_cc", link_radius_um=12.0, syp_percentile=25.0)
    assert out["in_germline"].iloc[:20].all()           # germline kept
    assert not out["in_germline"].iloc[20:].any()       # junk dropped
    assert any("dropped" in f for f in flags)


def test_keeps_syp_negative_distal_tip():
    # the DEFAULT (syp_seeded_cc): a SYP-bright synapsed core contiguous with a SYP-NEGATIVE distal
    # tip (mitotic/TZ, before SYP loads). Both are one gonad -> both kept; a far junk pair dropped.
    core = _blob(0, 0, 0, 20, 4.0, 100.0)             # x=0..76, SYP-bright
    distal = _blob(0, 0, 80, 5, 4.0, 0.0)             # x=80..96, contiguous (gap 4 µm), SYP-NEGATIVE
    junk = [{"centroid_z_um": 0, "centroid_y_um": 300, "centroid_x_um": x,
             "central_element_mean_intensity": 0.0} for x in (0, 4)]
    df = _table(core + distal + junk)
    out, _ = select_germline(df)                      # default = syp_seeded_cc
    assert out["in_germline"].iloc[:25].all()         # whole gonad kept, incl SYP-negative distal tip
    assert not out["in_germline"].iloc[25:].any()     # far junk dropped


def test_recovers_split_germline():
    # two separated SYP-bright arms (a folded/split gonad) — multi_cc keeps BOTH, largest_cc only one
    arm1 = _blob(0, 0, 0, 15, 4.0, 100.0)
    arm2 = _blob(0, 300, 0, 12, 4.0, 100.0)             # far away -> separate component
    df = _table(arm1 + arm2)
    multi, _ = select_germline(df, method="multi_cc", link_radius_um=12.0)
    assert multi["in_germline"].all()                   # both arms kept
    largest, _ = select_germline(df, method="largest_cc", link_radius_um=12.0)
    assert largest["in_germline"].iloc[:15].all() and not largest["in_germline"].iloc[15:].any()


def test_no_syp_channel_keeps_all():
    df = _table([{"centroid_z_um": 0, "centroid_y_um": 0, "centroid_x_um": i} for i in range(6)])
    out, flags = select_germline(df)                    # no central_element_mean_intensity column
    assert out["in_germline"].all()
    assert any("no_syp_channel" in f for f in flags)


def test_empty_table():
    out, flags = select_germline(pd.DataFrame(columns=["nucleus_id"]))
    assert "in_germline" in out.columns and len(out) == 0
