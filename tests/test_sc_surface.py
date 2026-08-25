"""SC ribbon surfacing (ridge + hysteresis mask) on a planted filament."""
import numpy as np

from germquant.sc import sc_ribbon_length, surface_sc_ribbon

SP = (0.4, 0.2, 0.2)


def _nucleus_with_filament():
    Z, Y, X = 12, 40, 60
    zz, yy, xx = np.indices((Z, Y, X)).astype(np.float32)
    nuc = (((zz - 6) * SP[0]) ** 2 + ((yy - 20) * SP[1]) ** 2 + ((xx - 30) * SP[2]) ** 2) < 2.2 ** 2
    labels = np.zeros((Z, Y, X), np.int32)
    labels[nuc] = 1
    syp = np.zeros((Z, Y, X), np.float32)
    syp[5:8, 19:22, 12:48] = 3000.0        # a bright SC filament through the nucleus
    return labels, syp, nuc


def test_recovers_planted_ridge():
    labels, syp, nuc = _nucleus_with_filament()
    mask = surface_sc_ribbon(syp, labels, SP)

    assert mask.dtype == bool
    assert mask.shape == labels.shape
    assert mask.any()                       # the ridge was surfaced
    assert mask[nuc].sum() > 0              # inside the nucleus
    assert not mask[~nuc].any()            # confined to the nucleus interior


def test_no_signal_yields_empty_mask():
    labels, _syp, _nuc = _nucleus_with_filament()
    flat = np.zeros(labels.shape, np.float32)
    mask = surface_sc_ribbon(flat, labels, SP)
    assert not mask.any()


def test_keep_ids_limits_work():
    labels, syp, _ = _nucleus_with_filament()
    # asking to keep a non-existent nucleus id -> nothing surfaced
    mask = surface_sc_ribbon(syp, labels, SP, keep_nucleus_ids={999})
    assert not mask.any()


def test_length_is_nonnegative_or_nan():
    mask = np.zeros((4, 10, 20), bool)
    mask[2, 5, 4:16] = True
    length = sc_ribbon_length(mask, SP)
    assert np.isnan(length) or length >= 0    # NaN if skan isn't installed, else a µm length
