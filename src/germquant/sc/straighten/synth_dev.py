"""Development phantom families of the v2 round (not registered protocol families; docs/SC_VALIDATION_PROTOCOL.md,
changelog entries 11 to 14). They live here so that germquant.sc.straighten.synth, the generator of protocol 2.3,
and every registered scene stay identical for every seed. Each family has its own family code; seeds at or above
7,000,000 are refused (synth.check_seed) like every development seed.

* ``psf_check`` (entry 12): six dense nuclei of 2.3 and three isolated arcs tilted U(30, 75) degrees in one image,
  the SYP channel rendered with the scene's PSF (sigma_z 0.33, 0.38, 0.43 or 0.48 um by seed modulo 4, sigma_xy
  U(0.135, 0.165) um), a foci channel of RAD-51-like point spots in the dense nuclei rendered with that PSF times
  595/511, the spot table with the localisation noise of protocol 4.4, and the wavelengths 511 and 595 nm.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from . import synth as S

FAMILY_CODE = {"psf_check": 107}
PSF_LEVELS = (0.33, 0.38, 0.43, 0.48)
WAVELENGTHS = {"central_element": 511.0, "foci": 595.0}
FOCI_DENSITY = 0.15                     # per um3 of a dense nucleus
FOCI_PEAK_RANGE = (150.0, 450.0)


def _rng(family: str, seed: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([int(seed), FAMILY_CODE[family]]))


def nucleoplasm(shape, spacing, centers, radii, levels, sigma) -> np.ndarray:
    """synth's nucleoplasm (uniform levels inside each nucleus sphere) blurred with ``sigma`` (z, y, x) um."""
    zz, yy, xx = np.meshgrid(*[np.arange(n) * s for n, s in zip(shape, spacing)], indexing="ij")
    out = np.zeros(shape, np.float64)
    for c, r, lev in zip(centers, radii, levels):
        out += lev * (((zz - c[0]) ** 2 + (yy - c[1]) ** 2 + (xx - c[2]) ** 2) <= r ** 2)
    return ndi.gaussian_filter(out, np.asarray(sigma) / np.asarray(spacing), mode="nearest")


def _nuclei_table(labels, centers, radii, label_radii, germ_cells, kinds) -> pd.DataFrame:
    vol_vox = float(np.prod(S.SPACING))
    rows = []
    for i, (c, r, lr, gc, k) in enumerate(zip(centers, radii, label_radii, germ_cells, kinds), start=1):
        n = int((labels == i).sum())
        rows.append({"nucleus_id": i, "centroid_z_um": c[0], "centroid_y_um": c[1], "centroid_x_um": c[2],
                     "radius_um": r, "label_radius_um": lr, "germ_cell": gc, "in_germline": True,
                     "volume_um3": n * vol_vox, "n_voxels": n, "kind": k})
    return pd.DataFrame(rows)


def psf_check_scene(seed: int, *, final: bool = False) -> S.Scene:
    seed = S.check_seed(seed, final)
    rng = _rng("psf_check", seed)
    sz = PSF_LEVELS[seed % 4]
    sxy = float(rng.uniform(0.135, 0.165))
    psf = (sz, sxy, sxy)
    scale = WAVELENGTHS["foci"] / WAVELENGTHS["central_element"]
    psf_f = tuple(float(v) * scale for v in psf)
    centers, radii, label_r, gcs, kinds, ribs = [], [], [], [], [], []
    foci_truth = []                       # (nucleus index, position)
    x = 0.0
    for i in range(6):
        R = float(rng.uniform(2.2, 2.6))
        c = np.array([0.0, 0.0, x])
        gc = "oocyte" if (i + seed) % 2 == 0 else "spermatocyte"
        rs = S._shell_nucleus(rng, gc, c, R, S.SCENE_PEAK, i + 1, len(ribs) + 1)
        n_f = int(rng.poisson(FOCI_DENSITY * 4.0 / 3.0 * np.pi * R ** 3))
        for _ in range(n_f):
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            foci_truth.append((i + 1, c + d * R * rng.uniform() ** (1.0 / 3.0), float(rng.uniform(*FOCI_PEAK_RANGE))))
        centers.append(c)
        radii.append(R)
        label_r.append(R + float(rng.uniform(-0.2, 0.4)))
        gcs.append(gc)
        kinds.append("dense")
        ribs.extend(rs)
        x += 6.0
    prev_r = 3.0
    for j in range(3):
        rib = S._iso_like(rng, float(rng.uniform(30.0, 75.0)), float(rng.uniform(2.0, 8.0)),
                          float(rng.uniform(300.0, 800.0)), nucleus_id=7 + j, ribbon_id=len(ribs) + 1, kind="psf_arc")
        c0, R = S._enclosing_nucleus([rib])
        x += prev_r + R + 1.0
        rib.points = rib.points - c0 + np.array([0.0, 0.0, x])
        centers.append(np.array([0.0, 0.0, x]))
        radii.append(R)
        label_r.append(R + float(rng.uniform(0.0, 0.4)))
        gcs.append("unknown")
        kinds.append("arc")
        ribs.append(rib)
        prev_r = R + 0.4
    shape, shift = S._place([r.points for r in ribs], centers, label_r)
    for r in ribs:
        r.points = r.points + shift
    centers = [c + shift for c in centers]
    labels, _ = S._sphere_masks(shape, S.SPACING, centers, label_r)
    levels = [float(rng.uniform(*S.NUCLEOPLASM)) for _ in centers]
    nuc = nucleoplasm(shape, S.SPACING, centers, radii, levels, psf)
    clean = np.full(shape, S.OFFSET) + nuc
    for r in ribs:
        S.render_ribbon(clean, r, S.SPACING, psf)
    fclean = np.full(shape, S.OFFSET) + nuc
    rows = []
    for k, (nid, p, pk) in enumerate(foci_truth):
        p = p + shift
        S._render_points(fclean, p[None, :], S.SPACING, psf_f, weight=pk)
        rows.append({"spot_id": k, "nucleus_id": nid, "z_um": p[0], "y_um": p[1], "x_um": p[2], "peak": pk})
    ce = S.add_noise(clean, rng)
    fo = S.add_noise(fclean, rng)
    truth = pd.DataFrame(rows, columns=["spot_id", "nucleus_id", "z_um", "y_um", "x_um", "peak"])
    if len(truth):
        loc = rng.normal(size=(len(truth), 3)) * np.asarray(S.FOCI_LOC_SIGMA)
        spots = pd.DataFrame({"spot_id": truth["spot_id"].to_numpy(), "nucleus_id": truth["nucleus_id"].to_numpy(),
                              "z_um": truth["z_um"].to_numpy() + loc[:, 0], "y_um": truth["y_um"].to_numpy() + loc[:, 1],
                              "x_um": truth["x_um"].to_numpy() + loc[:, 2]})
    else:
        spots = pd.DataFrame(columns=["spot_id", "nucleus_id", "z_um", "y_um", "x_um"])
    info = {"family": "psf_check", "seed": seed, "psf_sigma_um": list(psf), "foci_psf_sigma_um": list(psf_f),
            "sigma_z_level": sz, "wavelengths": dict(WAVELENGTHS), "nucleoplasm": levels,
            "foci_truth": truth.to_dict("list")}
    return S.Scene(family="psf_check", seed=seed, spacing=S.SPACING, channels={"central_element": ce, "foci": fo},
                   labels=labels, nuclei=_nuclei_table(labels, centers, radii, label_r, gcs, kinds), ribbons=ribs,
                   spots=spots, truth_foci=pd.DataFrame(columns=["spot_id", "ribbon_id", "nucleus_id", "s_um", "z_um",
                                                                 "y_um", "x_um"]), info=info)


FAMILY_FUNCS = {"psf_check": psf_check_scene}


def make_scene(family: str, seed: int, *, final: bool = False, **kw) -> S.Scene:
    if family not in FAMILY_FUNCS:
        raise ValueError(f"unknown development family {family!r}; one of {sorted(FAMILY_FUNCS)}")
    return FAMILY_FUNCS[family](seed, final=final, **kw)
