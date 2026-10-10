"""Signal QC stage: is each channel measurable, and does one channel leak into another?

Motivation (docs/QC.md, "Signal QC"): a session of ccw77 images carried a dead mCherry::SYP-3
fluorophore imaged at 97 percent power and 300 ms through the wrong emission filter, so the 545
channel held camera noise plus whatever the filter admitted from the 647 lamin stain, and the
partition coefficient floor came out above 1. Nothing in the pipeline noticed, because every stage
assumed its channel carried signal. This stage measures that assumption on every run and writes the
numbers next to the biology, in image_summary (one block of ``sig_<role>_*`` columns per channel
role) and in the ``channel_crosstalk`` table (one row per unordered pair of channel roles and
region).

Everything is computed on a deterministic xy sub-sample of the stack (``signal_qc.stride``, every z
plane), so a 2048 x 2048 x 40 stack costs seconds, and on three regions derived from the nucleus
labels with an anisotropy-aware Euclidean distance transform:

* ``tissue``: the nucleus labels dilated by ``signal_qc.tissue_dilation_um`` (the gonad and its
  cytoplasm); the voxels OUTSIDE it are the background sample;
* ``nuclei``: the labels themselves;
* ``perinuclear_ring``: labels dilated by ``signal_qc.ring_um`` minus the labels (where P granules
  and the nuclear envelope sit, so where a lamin stain leaking into another channel would show);
* ``cytoplasm``: tissue minus ring minus nuclei.

Per channel: ``offset`` (0.5th percentile of the whole sub-sample, the camera offset), ``bg``
(median of the background sample, or the 5th percentile of the channel when fewer than
``signal_qc.min_bg_voxels`` voxels lie outside the tissue; ``bg_source`` says which), ``noise``
(1.4826 x the median absolute deviation of the background sample, a robust sigma), ``p99`` (99th
percentile inside the tissue), ``snr`` = (p99 - bg) / noise, ``sat_frac`` (fraction of all sampled
voxels at or above the saturation level), and, when the acquisition stage supplied the exposure and
laser power for the role, ``bg_per_ms`` = (bg - offset) / exposure_ms, ``p99_per_ms`` =
(p99 - offset) / exposure_ms and ``bg_per_ms_pct`` = bg_per_ms / power, which make background
levels comparable across sessions imaged with different settings.

Saturation level rule: ``2 ** bit_depth - 1`` when a bit depth is known (the acquisition record's
``bit_depth_<role>`` field first, then the ``signal_qc.bit_depth`` config key), else the maximum of
the voxel dtype (the in-memory integer dtype, else the integer dtype the file was read as, else the
float maximum, in which case sat_frac is 0 by construction). A 12-bit camera writing into uint16
files saturates at 4095, not 65535, so without a bit depth sat_frac under-reports on such data.

Cross-talk: for every unordered pair of resolved roles and each region, the Pearson r of the two
channels over the region's voxels and the least-squares slope of the offset-subtracted second
channel on the first. Two stains of different structures correlate weakly; a channel that is a
scaled copy of another (bleed-through, or the same fluorophore seen through two filters) gives r
near 1 with a slope equal to the leak fraction.

Flags are OFF by default (``signal_qc.flags.enabled``), because a new flag on an existing run would
change its qc_flags column; a profile that enables them gets ``signal:<role>_low_snr``,
``signal:<role>_saturated`` and ``crosstalk:<a>_<b>_<region>``.
"""
from __future__ import annotations

import logging
import math

import numpy as np
import pandas as pd

from . import schema

log = logging.getLogger(__name__)

# code defaults of every signal_qc.* key (read with cfg.get(key, default); the config files are not edited)
DEFAULTS = {
    "stride": 4,                 # xy sub-sampling step (every z plane is kept)
    "tissue_dilation_um": 2.0,   # labels dilated by this = tissue; outside it = background sample
    "ring_um": 1.0,              # labels dilated by this minus labels = perinuclear ring
    "min_bg_voxels": 1000,       # fewer background voxels than this -> bg falls back to the 5th percentile
    "min_snr": 3.0,              # flag signal:<role>_low_snr below this (flags enabled only)
    "max_sat_frac": 0.001,       # flag signal:<role>_saturated above this (flags enabled only)
    "crosstalk_max_r": 0.8,      # flag crosstalk:<a>_<b>_<region> above this Pearson r (flags enabled only)
}

REGIONS = ("nuclei", "perinuclear_ring", "cytoplasm")
METRICS = schema.SIGNAL_QC_METRICS
BG_OUTSIDE, BG_P5 = "outside_tissue", "p5"


def subsample(arr: np.ndarray, stride: int) -> np.ndarray:
    """The deterministic xy sub-sample (a view: every z plane, every ``stride``-th row and column)."""
    s = max(1, int(stride))
    return arr[:, ::s, ::s]


def saturation_level(dtype, bit_depth=None, source_dtype=None) -> float:
    """The intensity at which a voxel counts as saturated (module docstring, "Saturation level rule")."""
    bd = _finite(bit_depth)
    if bd is not None and bd > 0:
        return float(2 ** int(bd) - 1)
    for dt in (dtype, source_dtype):
        if dt is None:
            continue
        try:
            d = np.dtype(dt)
        except TypeError:
            continue
        if np.issubdtype(d, np.integer):
            return float(np.iinfo(d).max)
    for dt in (dtype, source_dtype):
        if dt is None:
            continue
        try:
            d = np.dtype(dt)
        except TypeError:
            continue
        if np.issubdtype(d, np.floating):
            return float(np.finfo(d).max)
    return float("nan")


def region_masks(labels: np.ndarray, spacing, *, tissue_dilation_um: float, ring_um: float) -> dict[str, np.ndarray]:
    """nuclei / perinuclear_ring / cytoplasm / tissue / outside boolean masks on the given grid, from one
    anisotropy-aware distance transform of the nucleus labels (``spacing`` is the grid's own (dz, dy, dx))."""
    from scipy import ndimage as ndi

    nuc = np.asarray(labels) > 0
    if nuc.any():
        dt = ndi.distance_transform_edt(~nuc, sampling=tuple(float(s) for s in spacing))
        tissue = dt <= float(tissue_dilation_um)          # dt is 0 inside the labels, so they are included
        ring = (dt <= float(ring_um)) & ~nuc
    else:
        tissue = np.zeros(nuc.shape, bool)
        ring = np.zeros(nuc.shape, bool)
    cyto = tissue & ~ring & ~nuc
    return {"nuclei": nuc, "perinuclear_ring": ring, "cytoplasm": cyto, "tissue": tissue, "outside": ~tissue}


def _finite(v) -> float | None:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _mad_sigma(v: np.ndarray) -> float:
    """1.4826 x median absolute deviation: the robust sigma of a Gaussian sample."""
    if v.size == 0:
        return float("nan")
    med = float(np.median(v))
    return 1.4826 * float(np.median(np.abs(v - med)))


def channel_stats(ch: np.ndarray, regions: dict[str, np.ndarray], *, sat_level: float, min_bg_voxels: int,
                  exposure_ms=None, power=None) -> dict[str, float | str]:
    """The per-channel metrics (module docstring) for one sub-sampled channel, keyed by metric name."""
    v = np.asarray(ch, dtype=np.float64)
    flat = v.ravel()
    offset = float(np.percentile(flat, 0.5))
    outside = v[regions["outside"]]
    if outside.size >= int(min_bg_voxels):
        bg_vox, bg, bg_source = outside, float(np.median(outside)), BG_OUTSIDE
    else:
        # too little of the field is background: the darkest 5 percent of the channel stand in (their
        # MAD is a lower bound on the noise, since the sample is cut from above)
        bg = float(np.percentile(flat, 5))
        bg_vox, bg_source = flat[flat <= bg], BG_P5
    noise = _mad_sigma(bg_vox)
    tissue = v[regions["tissue"]]
    p99 = float(np.percentile(tissue, 99)) if tissue.size else float("nan")
    snr = (p99 - bg) / noise if (math.isfinite(noise) and noise > 0 and math.isfinite(p99)) else float("nan")
    sat_frac = float(np.mean(flat >= sat_level)) if (math.isfinite(sat_level) and flat.size) else float("nan")
    exp_ms, pw = _finite(exposure_ms), _finite(power)
    bg_per_ms = (bg - offset) / exp_ms if (exp_ms is not None and exp_ms > 0) else float("nan")
    p99_per_ms = (p99 - offset) / exp_ms if (exp_ms is not None and exp_ms > 0) else float("nan")
    bg_per_ms_pct = bg_per_ms / pw if (pw is not None and pw > 0) else float("nan")
    return {
        "offset": offset, "bg": bg, "bg_source": bg_source, "noise": noise, "p99": p99, "snr": snr,
        "sat_frac": sat_frac, "bg_per_ms": bg_per_ms, "p99_per_ms": p99_per_ms, "bg_per_ms_pct": bg_per_ms_pct,
    }


def pair_stats(a: np.ndarray, b: np.ndarray, offset_a: float, offset_b: float) -> tuple[float, float]:
    """(pearson_r, slope of offset-subtracted b on a) over paired voxel samples; NaN when either is constant."""
    if a.size < 2:
        return float("nan"), float("nan")
    x = np.asarray(a, dtype=np.float64) - float(offset_a)
    y = np.asarray(b, dtype=np.float64) - float(offset_b)
    x = x - x.mean()
    y = y - y.mean()
    sxx, syy, sxy = float(x @ x), float(y @ y), float(x @ y)
    if sxx <= 0 or syy <= 0:
        return float("nan"), float("nan")
    return sxy / math.sqrt(sxx * syy), sxy / sxx


def crosstalk_table(channels: dict[str, np.ndarray], offsets: dict[str, float], regions: dict[str, np.ndarray],
                    image_id: str = "") -> pd.DataFrame:
    """One row per unordered pair of roles (in the order given) and region; schema.CHANNEL_CROSSTALK columns."""
    roles = list(channels)
    rows = []
    for i, ra in enumerate(roles):
        for rb in roles[i + 1:]:
            for region in REGIONS:
                m = regions[region]
                n = int(m.sum())
                r, slope = pair_stats(channels[ra][m], channels[rb][m], offsets[ra], offsets[rb]) if n else (float("nan"),) * 2
                rows.append({"image_id": image_id, "role_a": ra, "role_b": rb, "region": region,
                             "n_voxels": n, "pearson_r": r, "slope": slope})
    cols = ["image_id", *schema.CHANNEL_CROSSTALK]
    return pd.DataFrame(rows, columns=cols) if rows else pd.DataFrame(columns=cols)


def run_signal_qc(data: np.ndarray, labels: np.ndarray, spacing, role_to_idx: dict, *, acq_fields: dict | None = None,
                  stride: int = DEFAULTS["stride"], tissue_dilation_um: float = DEFAULTS["tissue_dilation_um"],
                  ring_um: float = DEFAULTS["ring_um"], min_bg_voxels: int = DEFAULTS["min_bg_voxels"],
                  bit_depth=None, source_dtype=None, image_id: str = "") -> dict:
    """The stage body. ``data`` is the (C, Z, Y, X) stack, ``labels`` the nucleus label image on the same
    grid, ``spacing`` its (dz, dy, dx) in microns, ``role_to_idx`` the resolved channel roles (None values
    are skipped), ``acq_fields`` the acquisition stage's per-image fields (``exp_ms_<role>``,
    ``laser_pow_<role>`` and, when present, ``bit_depth_<role>``). Returns ``summary`` (the
    ``sig_<role>_<metric>`` image_summary fields for every resolved role), ``crosstalk`` (the table),
    ``roles`` (the roles measured, in channel-map order), ``n_sampled`` and the ``sat_level`` per role."""
    acq = acq_fields or {}
    s = max(1, int(stride))
    lab_sub = subsample(np.asarray(labels), s)
    sp_sub = (float(spacing[0]), float(spacing[1]) * s, float(spacing[2]) * s)
    regions = region_masks(lab_sub, sp_sub, tissue_dilation_um=tissue_dilation_um, ring_um=ring_um)
    n_ch = int(data.shape[0])
    roles = [r for r, i in role_to_idx.items() if i is not None and 0 <= int(i) < n_ch]
    channels: dict[str, np.ndarray] = {}
    offsets: dict[str, float] = {}
    sat_levels: dict[str, float] = {}
    summary: dict[str, float | str] = {}
    for role in roles:
        ch = subsample(np.asarray(data[int(role_to_idx[role])]), s)
        if ch.shape != lab_sub.shape:
            raise ValueError(f"channel {role} shape {ch.shape} does not match the labels {lab_sub.shape}")
        bd = _finite(acq.get(f"bit_depth_{role}"))          # the acquisition record wins over the config key
        level = saturation_level(ch.dtype, bd if bd is not None else bit_depth, source_dtype)
        stats = channel_stats(ch, regions, sat_level=level, min_bg_voxels=min_bg_voxels,
                              exposure_ms=acq.get(f"exp_ms_{role}"), power=acq.get(f"laser_pow_{role}"))
        for metric in METRICS:
            summary[f"sig_{role}_{metric}"] = stats[metric]
        channels[role], offsets[role], sat_levels[role] = ch, float(stats["offset"]), level
    table = crosstalk_table(channels, offsets, regions, image_id=image_id)
    return {"summary": summary, "crosstalk": table, "roles": roles, "n_sampled": int(lab_sub.size),
            "sat_level": sat_levels}


def signal_flags(summary: dict, crosstalk: pd.DataFrame, roles, *, min_snr: float = DEFAULTS["min_snr"],
                 max_sat_frac: float = DEFAULTS["max_sat_frac"], max_r: float = DEFAULTS["crosstalk_max_r"]) -> list[str]:
    """The QC flags a profile with ``signal_qc.flags.enabled`` gets; a NaN metric never fires a flag."""
    flags: list[str] = []
    for role in roles:
        snr = _finite(summary.get(f"sig_{role}_snr"))
        if snr is not None and snr < float(min_snr):
            flags.append(f"signal:{role}_low_snr")
        sat = _finite(summary.get(f"sig_{role}_sat_frac"))
        if sat is not None and sat > float(max_sat_frac):
            flags.append(f"signal:{role}_saturated")
    if crosstalk is not None and len(crosstalk):
        for row in crosstalk.itertuples(index=False):
            r = _finite(row.pearson_r)
            if r is not None and r > float(max_r):
                flags.append(f"crosstalk:{row.role_a}_{row.role_b}_{row.region}")
    return flags
