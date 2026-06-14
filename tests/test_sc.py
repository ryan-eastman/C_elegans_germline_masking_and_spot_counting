"""SC tracing (v2): physical length, the fragmentation INDEX (heat phenotype), zero-rows, and
spacing-awareness. Per-nucleus fragment COUNT is a documented lower bound (strands overlap in 3D),
so the heat readout is the SYP-intensity-CV fragmentation index, not the component count.

Uses skan (germquant[sc]); skipped if not installed.
"""
import numpy as np
import pytest

from germquant.sc.skeleton import skan_available, trace_sc

pytestmark = pytest.mark.skipif(not skan_available(), reason="skan not installed")

SPACING = (0.4, 0.15, 0.15)  # anisotropic, like real confocal


def _volume():
    Z, Y, X = 16, 24, 64
    labels = np.zeros((Z, Y, X), np.int32)
    labels[2:14, 4:20, 4:60] = 1            # one nucleus filling most of the box
    sc = np.zeros((Z, Y, X), np.float32)
    return labels, sc


def test_single_filament_one_fragment():
    labels, sc = _volume()
    sc[8, 12, 6:58] = 1000.0                 # one straight SYP filament along x
    tracks, per_nuc = trace_sc(labels=labels, sc_img=sc, spacing=SPACING,
                               intensity_percentile=80.0, min_fragment_length_um=0.5)
    assert len(per_nuc) == 1
    row = per_nuc.iloc[0]
    assert row["n_fragments"] == 1, f"expected 1 fragment, got {row['n_fragments']}"
    # physical length ~ (58-6-1)*0.15 ≈ 7.6 µm
    assert 5.0 < row["sc_total_length_um"] < 10.0
    assert tracks.iloc[0]["trace_method"] == "skan_iso"
    assert np.isfinite(tracks.iloc[0]["tortuosity"])   # B5: tortuosity is computed, not NaN


def test_broken_filament_raises_fragmentation_index():
    # The ~6 SCs overlap in 3D, so a gap can't be reliably split into separate components (documented
    # in skeleton.py). The heat readout is the fragmentation INDEX (SYP intensity CV), which is higher
    # for a broken filament than an intact one.
    labels, sc = _volume()
    intact = sc.copy()
    intact[8, 12, 6:58] = 1000.0
    broken = sc.copy()
    broken[8, 12, 6:26] = 1000.0
    broken[8, 12, 42:58] = 1000.0
    _, pi = trace_sc(labels=labels, sc_img=intact, spacing=SPACING, min_fragment_length_um=0.5)
    _, pb = trace_sc(labels=labels, sc_img=broken, spacing=SPACING, min_fragment_length_um=0.5)
    assert float(pb.iloc[0]["sc_fragmentation_index"]) > float(pi.iloc[0]["sc_fragmentation_index"])


def test_desynapsed_nucleus_is_a_zero_not_a_dropped_row():
    labels, sc = _volume()                   # sc all zero -> nothing traces
    _, per_nuc = trace_sc(labels=labels, sc_img=sc, spacing=SPACING)
    assert len(per_nuc) == 1                  # the row exists...
    assert int(per_nuc.iloc[0]["n_fragments"]) == 0   # ...as a real zero


def test_length_scales_with_spacing():
    # length is physical (µm), not voxels: doubling the voxel size ~doubles the measured length.
    # use a realistic PSF-blurred tube — a 1-voxel line skeletonises with resolution-dependent spurs.
    from scipy import ndimage as ndi

    labels, sc = _volume()
    sc[8, 12, 6:58] = 1000.0
    sc = ndi.gaussian_filter(sc, (0.6, 1.2, 1.2))
    _, a = trace_sc(labels=labels, sc_img=sc, spacing=SPACING)
    _, b = trace_sc(labels=labels, sc_img=sc, spacing=tuple(2 * s for s in SPACING))
    la = float(a.iloc[0]["sc_total_length_um"])
    lb = float(b.iloc[0]["sc_total_length_um"])
    assert 1.6 < lb / la < 2.7, f"length not spacing-aware: {la:.2f} -> {lb:.2f}"
