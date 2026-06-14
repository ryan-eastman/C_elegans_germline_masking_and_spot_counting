#!/usr/bin/env python
"""Find a robust chromatin-asymmetry metric that separates TRANSITION-ZONE (crescent/bouquet) nuclei
from PACHYTENE (spread 6-strand) nuclei, on synthetic ground truth with perfect labels. The current
caller uses centroid-offset polarization (weak: real data tops out ~0.13, threshold 0.18 finds no TZ).
Compares candidate GEOMETRIC metrics by Cohen's d + ROC-AUC (TZ vs pachytene).

    python scripts/validate_zone_metric.py
"""
import os
import sys

import numpy as np
import pandas as pd
from skimage.measure import regionprops

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from synth_nuclei import DX, DY, DZ, degrade, synthesize_zone_scene  # noqa: E402

SP = np.array([DZ, DY, DX])
VOXVOL = float(np.prod(SP))


def nucleus_metrics(coords_um, inten):
    """Several chromatin-asymmetry metrics for one nucleus (coords in um, inten = DAPI per voxel)."""
    w = inten / inten.sum()
    geom = coords_um.mean(0)
    wc = (coords_um * w[:, None]).sum(0)
    radius = (3 * len(coords_um) * VOXVOL / (4 * np.pi)) ** (1 / 3)
    off = wc - geom
    offn = float(np.linalg.norm(off))
    polar = offn / radius if radius > 0 else 0.0                  # current metric
    # hemisphere asymmetry: mass on the denser side vs the other, along the offset axis
    if offn > 1e-9:
        d = off / offn
        s = (coords_um - geom) @ d
        near = inten[s > 0].sum()
        far = inten[s <= 0].sum()
        hemi = float((near - far) / (near + far + 1e-9))
    else:
        hemi = 0.0
    # moment anisotropy: elongation of the intensity-weighted mass (crescent = elongated)
    c = coords_um - wc
    cov = (c * w[:, None]).T @ c
    ev = np.sort(np.linalg.eigvalsh(cov))[::-1]
    anis = float((ev[0] - ev[2]) / (ev[0] + 1e-12))              # 0=isotropic ball, ->1 elongated
    # compactness: radius of gyration of the bright mass / nuclear radius. TZ chromatin is clustered
    # (small spread); pachytene strands ring the whole nucleus (large spread). Orthogonal to polarity.
    gyr = float(np.sqrt(max(np.trace(cov), 0.0)) / radius) if radius > 0 else 0.0
    # radial offset of mass: how far (in radii) the bulk sits from the geometric center
    rmean = float((np.linalg.norm(coords_um - geom, axis=1) * w).sum() / radius) if radius > 0 else 0.0
    # directional concentration: mass-weighted mean of unit vectors from the center. Pachytene mass
    # is spread over all angles (a ring w/ dark center) -> vectors cancel -> small; TZ mass is in one
    # angular sector (crescent) -> vectors align -> large. Independent of radial distance.
    rad = coords_um - geom
    rn = np.linalg.norm(rad, axis=1, keepdims=True)
    units = rad / np.maximum(rn, 1e-9)
    rlen = float(np.linalg.norm((units * w[:, None]).sum(0)))
    return dict(polar=polar, hemi=hemi, anis=anis, rmean=rmean, rlen=rlen, gyr=gyr)


def auc(neg, pos):
    neg, pos = np.asarray(neg), np.asarray(pos)
    if len(neg) == 0 or len(pos) == 0:
        return float("nan")
    r = pd.Series(np.concatenate([neg, pos])).rank().to_numpy()
    return float((r[len(neg):].sum() - len(pos) * (len(pos) + 1) / 2) / (len(neg) * len(pos)))


def cohens_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return float((b.mean() - a.mean()) / (sp + 1e-12))


def main():
    rows = []
    for seed in range(4):
        dapi, label, gt = synthesize_zone_scene(seed=seed, tz_frac=0.5)
        dimg = degrade(dapi, seed=100 + seed)
        zmap = dict(zip(gt.nucleus_id, gt.zone))
        for rp in regionprops(label, intensity_image=dimg):
            inten = rp.image_intensity[rp.image].astype(float)
            if inten.sum() <= 0 or len(rp.coords) < 20:
                continue
            m = nucleus_metrics(rp.coords.astype(float) * SP, inten)
            m.update(zone=zmap.get(rp.label, "?"), seed=seed)
            rows.append(m)
    df = pd.DataFrame(rows)
    df["polar_gyr"] = df["polar"] - df["gyr"]            # high polarity + low spread = crescent
    tz = df[df.zone == "transition_zone"]
    pa = df[df.zone == "pachytene"]
    print(f"n: TZ={len(tz)} pachytene={len(pa)}\n")
    print(f"{'metric':<10} {'TZ mean':>9} {'pachy mean':>11} {'Cohen_d':>9} {'AUC':>6}")
    for col in ["polar", "hemi", "anis", "rmean", "rlen", "gyr", "polar_gyr"]:
        print(f"{col:<10} {tz[col].mean():>9.3f} {pa[col].mean():>11.3f} "
              f"{cohens_d(pa[col], tz[col]):>9.2f} {auc(pa[col], tz[col]):>6.2f}")
    print("\n(higher Cohen_d / AUC = better TZ-vs-pachytene separation; current caller uses 'polar')")


if __name__ == "__main__":
    main()
