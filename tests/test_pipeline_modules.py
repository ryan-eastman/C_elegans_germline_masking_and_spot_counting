"""Synthetic end-to-end chain test: segment -> measure -> axis -> zones -> foci.
Catches integration breakage with no real data or GPU. Spacing is anisotropic.
"""
import numpy as np

from germquant.axis import linearize_germline
from germquant.foci import detect_foci
from germquant.measure import measure_objects
from germquant.segment import segment_nuclei
from germquant.zones import call_zones

SPACING = (0.5, 0.2, 0.2)  # anisotropic (z coarser), like real confocal


def _synthetic_germline():
    Z, Y, X = 24, 96, 170
    dapi = np.zeros((Z, Y, X), np.float32)
    foci = np.zeros((Z, Y, X), np.float32)
    centers = [(12, 48, 18 + i * 24) for i in range(6)]  # a row of nuclei along x
    zz, yy, xx = np.indices((Z, Y, X)).astype(np.float32)
    for (cz, cy, cx) in centers:
        r2 = (((zz - cz) * SPACING[0]) ** 2 + ((yy - cy) * SPACING[1]) ** 2 + ((xx - cx) * SPACING[2]) ** 2)
        dapi += np.exp(-r2 / (2 * 1.4 ** 2)) * 1000.0
        foci[cz, cy, min(cx + 3, X - 1)] = 6000.0  # one bright focus per nucleus
    return dapi, foci, centers


def test_chain_runs_and_is_spacing_aware():
    dapi, foci, centers = _synthetic_germline()

    labels, method = segment_nuclei(dapi, SPACING, method="classical",
                                    diameter_um=2.0, min_volume_um3=0.05)
    assert method == "classical"
    assert labels.max() >= 3, "should segment most of the 6 synthetic nuclei"

    nuclei = measure_objects(labels, {"dna": dapi}, SPACING).rename(columns={"label": "nucleus_id"})
    for col in ("nucleus_id", "volume_um3", "centroid_x_um", "dna_mean_intensity"):
        assert col in nuclei.columns
    assert (nuclei["volume_um3"] > 0).all()

    nuclei, conf, _ = linearize_germline(nuclei)
    assert "axis_position_norm" in nuclei.columns
    pos = nuclei["axis_position_norm"].dropna()
    assert pos.min() >= -1e-9 and pos.max() <= 1 + 1e-9
    assert conf > 0.8, "a straight row of nuclei should be highly 1-D"

    nuclei, zones, _ = call_zones(nuclei, dna=dapi, labels=labels, spacing=SPACING, method="dapi_crescent")
    assert "zone_call" in nuclei.columns

    fdf = detect_foci(foci, SPACING, labels=labels, min_sigma_um=0.2, max_sigma_um=0.6, threshold_rel=0.2)
    assert len(fdf) >= 1
    assert {"z_um", "y_um", "x_um", "nucleus_id", "detection_radius_um"} <= set(fdf.columns)
    # at least one focus landed inside a segmented nucleus
    assert (fdf["nucleus_id"] > 0).any()


def test_no_dapi_zoning_is_flagged_not_crashing():
    dapi, _, _ = _synthetic_germline()
    labels, _ = segment_nuclei(dapi, SPACING, method="classical", diameter_um=2.0, min_volume_um3=0.05)
    nuclei = measure_objects(labels, {"dna": dapi}, SPACING).rename(columns={"label": "nucleus_id"})
    nuclei, _, _ = linearize_germline(nuclei)
    # no dna provided -> synapsis_state stub, flagged, no crash
    nuclei2, zones, flags = call_zones(nuclei, dna=None, labels=None, spacing=SPACING, method="synapsis_state")
    assert any("synapsis_state" in f for f in flags)
    assert (nuclei2["zone_call"] == "unknown").all()
