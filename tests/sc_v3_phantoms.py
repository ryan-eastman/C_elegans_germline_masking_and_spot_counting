"""Phantoms with known truth for the tests of the v3 rules of sc_straighten (protocol changelog entries 32 and 33):
one nucleus holding ribbons (straight segments or arcs, with empty gaps or with per-point brightness), point
puncta, and RAD-51-like foci in a foci channel with their spot table, rendered with a chosen axial PSF."""
from __future__ import annotations

import json
import math

import numpy as np
import pandas as pd

from germquant.sc.straighten import synth, synth_dev
from germquant.sc.straighten.pipeline import straighten_image

WL = {"central_element": 511.0, "foci": 595.0}

# columns of the strands table that the reading of gaps feeds (entry 33: everything else equals v1's)
GAP_COLUMNS = {"n_gaps", "gap_um", "gap_fraction", "gap_max_um", "n_pieces", "longest_piece_um", "qc_flags",
               "qc_status", "strand_class"}


def sets(d: dict) -> list[str]:
    return [f"{k}={json.dumps(v)}" for k, v in d.items()]


def line(t0: float, t1: float, tilt_deg: float, origin=(0.0, 0.0, 0.0)) -> np.ndarray:
    """Dense truth points of a straight segment from t0 to t1 um along the direction tilted ``tilt_deg`` from the
    image plane toward z (axis order z, y, x)."""
    a = math.radians(tilt_deg)
    d = np.array([math.sin(a), 0.0, math.cos(a)])
    t = np.linspace(t0, t1, int(math.ceil((t1 - t0) / synth.TRUTH_STEP_UM)) + 1)
    return np.asarray(origin, float)[None, :] + t[:, None] * d[None, :]


def weights(points, stretches):
    """Per-point brightness factors: 1, and ``level`` over each (s0, s1, level) along the points."""
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))])
    w = np.ones(len(s))
    for s0, s1, level in stretches:
        w[(s >= s0) & (s <= s1)] = level
    return w


def scene(ribbons=(), puncta=(), *, sigma_z=0.365, foci=0, seed=0, R=2.4, nucleoplasm=60.0):
    """One nucleus with ``ribbons`` (dicts: points, peak, optional per-point weights and gaps), point ``puncta``
    (dicts: center offset, peak) and ``foci`` RAD-51-like foci at least 1.1 um apart in a foci channel rendered with
    the PSF times 595/511 (with their spot table), all with the axial PSF sigma ``sigma_z``."""
    rng = np.random.default_rng(seed)
    psf = (float(sigma_z), synth.PSF_SIGMA[1], synth.PSF_SIGMA[2])
    c0 = np.zeros(3)
    ribs = []
    for k, r in enumerate(ribbons):
        rib = synth.TruthRibbon(ribbon_id=k + 1, nucleus_id=1, points=np.asarray(r["points"], float),
                                peak=float(r.get("peak", 450.0)), gaps=list(r.get("gaps", [])), kind="custom")
        ribs.append((rib, r.get("weights")))
    fpos = []
    tries = 0
    while len(fpos) < foci and tries < 20000:
        tries += 1
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        q = d * R * 0.85 * rng.uniform() ** (1 / 3)
        if all(np.linalg.norm(q - x) >= 1.1 for x in fpos):
            fpos.append(q)
    pts = [np.asarray(q["center"], float) for q in puncta] + fpos
    shape, shift = synth._place([r.points for r, _w in ribs] + ([np.vstack(pts)] if pts else []), [c0], [R + 0.3])
    for r, _w in ribs:
        r.points = r.points + shift
    center = c0 + shift
    labels, _ = synth._sphere_masks(shape, synth.SPACING, [center], [R + 0.3])
    nuc = synth_dev.nucleoplasm(shape, synth.SPACING, [center], [R], [nucleoplasm], psf)
    clean = np.full(shape, synth.OFFSET) + nuc
    for r, w in ribs:
        if w is None:
            synth.render_ribbon(clean, r, synth.SPACING, psf)
        else:
            synth_dev._render_modulated(clean, r, psf, np.asarray(w, float))
    for q in puncta:
        synth._render_points(clean, (np.asarray(q["center"], float) + shift)[None, :], synth.SPACING, psf,
                             weight=float(q["peak"]))
    ce = synth.add_noise(clean, rng)
    channels = {"central_element": ce}
    spots = None
    if foci:
        psf_f = tuple(v * WL["foci"] / WL["central_element"] for v in psf)
        fclean = np.full(shape, synth.OFFSET) + nuc
        rows = []
        for k, q in enumerate(fpos):
            p = q + shift
            synth._render_points(fclean, p[None, :], synth.SPACING, psf_f, weight=float(rng.uniform(300, 450)))
            rows.append({"spot_id": k, "nucleus_id": 1, "z_um": p[0], "y_um": p[1], "x_um": p[2]})
        channels["foci"] = synth.add_noise(fclean, rng)
        spots = {"spots": pd.DataFrame(rows)}
    nuclei = pd.DataFrame([{"nucleus_id": 1, "germ_cell": "oocyte", "in_germline": True, "centroid_z_um": center[0],
                            "centroid_y_um": center[1], "centroid_x_um": center[2],
                            "volume_um3": float((labels == 1).sum() * np.prod(synth.SPACING))}])
    return {"ce": ce, "labels": labels, "nuclei": nuclei, "channels": channels, "spots": spots, "shift": shift,
            "ribbons": [r for r, _w in ribs]}


def run(sc, p, wavelengths=None):
    roles = {k: v for k, v in sc["channels"].items() if k != "central_element"}
    return straighten_image(sc["ce"], sc["labels"], synth.SPACING, sc["nuclei"], params=p, roles=roles,
                            spot_tables=sc["spots"], image_id="t", wavelengths=wavelengths)


def assert_same_strands(a, b):
    """The strands and samples of two runs are equal except for what the reading of gaps feeds."""
    assert len(a.strands) == len(b.strands)
    cols = [c for c in a.strands.columns if c not in GAP_COLUMNS]
    pd.testing.assert_frame_equal(a.strands[cols].reset_index(drop=True), b.strands[cols].reset_index(drop=True))
    scols = [c for c in a.samples.columns if c != "in_gap"]
    pd.testing.assert_frame_equal(a.samples[scols].reset_index(drop=True), b.samples[scols].reset_index(drop=True))
