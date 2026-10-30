"""Development phantom families of the v2 round (not registered protocol families; docs/SC_VALIDATION_PROTOCOL.md,
changelog entries 11 to 14). They live here so that germquant.sc.straighten.synth, the generator of protocol 2.3,
and every registered scene stay identical for every seed. Each family has its own family code; seeds at or above
7,000,000 are refused (synth.check_seed) like every development seed, and so is ``final=True`` (protocol changelog
entry 27): these families have no sealed seeds.

* ``psf_check`` (entry 12): six dense nuclei of 2.3 and three isolated arcs tilted U(30, 75) degrees in one image,
  the SYP channel rendered with the scene's PSF (sigma_z 0.33, 0.38, 0.43 or 0.48 um by seed modulo 4, sigma_xy
  U(0.135, 0.165) um), a foci channel of RAD-51-like point spots in the dense nuclei rendered with that PSF times
  595/511, the spot table with the localisation noise of protocol 4.4, and the wavelengths 511 and 595 nm.
* ``puncta_pt``, ``puncta_agg``, ``dense_puncta`` (the puncta entry; as in entry 8 (iii)): one nucleus of 2.3 with
  8 PSF-blurred point puncta, or 8 uniformly filled spheres of diameter U(0.25, 0.8) um, and no ribbon; or the
  dense nucleus of 2.3 plus 6 point puncta at least 0.6 um from every ribbon. Puncta peaks are log-uniform in
  [100, 1500] counts above the nucleoplasm after the PSF; puncta lie at least 1.2 um apart. ``dense_puncta_clean``
  is the same scene as ``dense_puncta`` (the same draws and the same noise) without the puncta. The puncta are in
  ``info["puncta"]`` (centre, diameter, peak).
* ``real_dense`` and its variants (the realism probe, protocol changelog entry 19): the dense nucleus of 2.3 at a
  scene peak of 200, 300, 430 or 700 counts (seed modulo 4) with four features calibrated on truth-free statistics
  of the twelve trial-4 runs (REAL_DENSE_CALIB): (a) the axial PSF drawn from the per-image compact-foci values, (b)
  a log-normal brightness modulation along each ribbon, (c) a PSF-correlated nucleoplasmic texture, (d) puncta.
  ``real_dense`` has all four; ``real_dense_base``, ``real_dense_psf``, ``real_dense_brightness``,
  ``real_dense_texture`` and ``real_dense_puncta`` have none or one. Each feature draws from its own random
  stream, so the geometry, the noise draws and the other features do not depend on which features are on.
* ``real_dense_r2`` and ``real_dense_r2_base``, ``_psf``, ``_brightness``, ``_texture``, ``_puncta`` (protocol
  changelog entry 28): the same scenes with the features recalibrated so that the phantoms show the real
  statistics after the PSF and the measurement (REAL_DENSE_CALIB_R2).
* ``puncta_psf`` (protocol changelog entry 32, the punctum class with the image's measured axial sigma): one image of
  two layers 6.5 um apart in z, each with three dense nuclei of 2.3 (oocyte and spermatocyte alternating, the 0.25
  um clearance of entry 2 (c)) holding 3 puncta each, one nucleus holding 6 puncta and no ribbon, and two steep
  isolated arcs (tilt U(60, 90) degrees, length U(0.9, 4.0) um, peak U(300, 800) counts) in their own enclosing
  nuclei; puncta are half points and half spheres of diameter U(0.25, 0.8) um, peaks log-uniform in [100, 1500]
  counts above the nucleoplasm, at least 0.6 um from every ribbon and 1.2 um apart. The axial PSF sigma is the
  level of psf_check (0.33, 0.38, 0.43 or 0.48 um by seed modulo 4) at the middle depth and changes linearly with
  depth by a relative gradient U(-0.02, 0.02) per um (up to 20 percent over 10 um, as the compact RAD-51 foci of
  the development images change, entry 9); the lateral sigma is U(0.135, 0.165) um. A foci channel holds
  RAD-51-like point spots in every nucleus but the arcs' (as psf_check: a Poisson number of mean 0.15 per um3,
  peaks U(150, 450) counts, the local PSF times 595/511), with the spot table of protocol 4.4 and the wavelengths
  511 and 595 nm, so the foci estimate of entry 12 measures the PSF the puncta were rendered with.
* ``gaps_psf`` (protocol changelog entry 33, dim stretches that are not gaps): one image of two layers 6.5 um apart in
  z, each with six isolated ribbons in their own enclosing nuclei (length U(3.0, 7.0) um, lengthened as the gaps
  family does when its events do not fit, tilt U(0, 75) degrees, peak U(300, 800) counts): the first, third and
  fifth of each layer hold 1 or 2 empty gaps of U(0.3, 1.0) um (no signal, as the gaps family; 1 below 4 um), the
  others 1 or 2 dim stretches of U(0.4, 1.5) um whose brightness is U(0.15, 0.6) of the ribbon's (at least 0.8 um
  of ribbon before, between and after the events). The PSF is that of puncta_psf (the psf_check level by seed
  modulo 4 at the middle depth, a relative depth gradient U(-0.02, 0.02) per um, lateral U(0.135, 0.165) um), and
  a foci channel holds RAD-51-like point spots in every nucleus (as puncta_psf), with the spot table of protocol
  4.4 and the wavelengths 511 and 595 nm. The empty gaps are the ribbons' ``gaps``; the dim stretches are in
  ``info["dim_stretches"]`` ({ribbon_id: [[s0, s1, fill], ...]}, um along the ribbon).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from . import synth as S

FAMILY_CODE = {"puncta_pt": 101, "puncta_agg": 102, "dense_puncta": 103, "psf_check": 107, "real_dense": 108,
               "puncta_psf": 109, "gaps_psf": 110}
# entry 19: truth-free statistics of the twelve trial-4 runs (C:/Users/ryane/sc_ribbon_work/v2/logs/
# realism_measures_t4.json; v2/scripts/realism_measure.py)
REAL_DENSE_CALIB = {
    "psf_sigma_z_values_um": (0.4635, 0.4311, 0.4748, 0.4534, 0.4491, 0.4200, 0.4211, 0.4134, 0.3620, 0.3982,
                              0.3866, 0.3873),
    "brightness_var_log": 0.1224,
    "brightness_corr_length_um": 0.3126,
    "texture_rel_excess_var": 0.1082,
    "puncta_mean_by_tercile": (1.0690, 0.8690, 0.6342),     # lowest, middle, highest SYP brightness
}
REAL_PEAKS = (200.0, 300.0, 430.0, 700.0)
# entry 28 (d): the probe recalibrated so that its features show the real statistics after the PSF and the
# measurement (C:/Users/ryane/sc_ribbon_work/v2/logs/realism_recal_[abc]_9726612_tuning.json): brightness by
# entry 28 (a), texture by entry 28 (b) (its phantom statistic without texture, 0.141, already exceeds the real
# 0.108, so the latent texture is 0), puncta by entry 28 (c); the PSF feature is that of entry 19
REAL_DENSE_CALIB_R2 = {**REAL_DENSE_CALIB,
                       "brightness_var_log": 0.5,
                       "brightness_corr_length_um": 0.1,
                       "texture_rel_excess_var": 0.0,
                       "puncta_mean_by_tercile": (3.4613, 2.8639, 2.099)}
REAL_FEATURES = ("psf", "brightness", "texture", "puncta")
PSF_LEVELS = (0.33, 0.38, 0.43, 0.48)
WAVELENGTHS = {"central_element": 511.0, "foci": 595.0}
FOCI_DENSITY = 0.15                     # per um3 of a dense nucleus
FOCI_PEAK_RANGE = (150.0, 450.0)


def _dev_seed(seed: int, final: bool) -> int:
    """The seed check of the development families: no sealed seed is generated for them (protocol changelog entry
    27; whether they get sealed seeds at the scoring is decision 9 of the PR), so ``final`` is refused."""
    if final:
        raise S.SealedSeedError("the development families of germquant.sc.straighten.synth_dev have no sealed seeds "
                                "(protocol changelog entry 27); development seeds are below 7,000,000")
    return S.check_seed(seed, False)


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
    seed = _dev_seed(seed, final)
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


def _sphere_points(diam: float, step: float = 0.025) -> np.ndarray:
    r = diam / 2.0
    g = np.arange(-r, r + 1e-9, step)
    zz, yy, xx = np.meshgrid(g, g, g, indexing="ij")
    P = np.stack([zz.ravel(), yy.ravel(), xx.ravel()], 1)
    return P[np.linalg.norm(P, axis=1) <= r]


def puncta_scene(seed: int, *, final: bool = False, family: str = "puncta_pt", with_puncta: bool = True) -> S.Scene:
    """The puncta families (the generator of the entry-8 diagnostics, v2/scripts/steep_phantoms.py, with the same
    family codes): one nucleus with puncta and, for dense_puncta, the dense ribbons of 2.3."""
    seed = _dev_seed(seed, final)
    rng = _rng(family, seed)
    R = float(rng.uniform(2.2, 2.6))
    label_r = R + float(rng.uniform(-0.2, 0.4))
    c0 = np.zeros(3)
    ribs = []
    if family == "dense_puncta":
        gc = "oocyte" if seed % 2 == 0 else "spermatocyte"
        ribs = S._shell_nucleus(rng, gc, c0, R, S.SCENE_PEAK, 1, 1)
    else:
        gc = "oocyte"
    n_p = 6 if family == "dense_puncta" else 8
    rib_pts = np.vstack([r.points[::10] for r in ribs]) if ribs else np.zeros((0, 3))
    puncta = []
    tries = 0
    while len(puncta) < n_p and tries < 5000:
        tries += 1
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        q = c0 + d * R * rng.uniform(0.3, 0.95) ** (1 / 3) * 0.95
        if any(np.linalg.norm(q - x["center"]) < 1.2 for x in puncta):
            continue
        if len(rib_pts) and np.min(np.linalg.norm(rib_pts - q, axis=1)) < 0.6:
            continue
        diam = float(rng.uniform(0.25, 0.8)) if family == "puncta_agg" else 0.0
        peak = float(np.exp(rng.uniform(np.log(100.0), np.log(1500.0))))
        puncta.append({"center": q, "diam": diam, "peak": peak})
    shape, shift = S._place([r.points for r in ribs] + [np.array([x["center"] for x in puncta])], [c0], [label_r])
    for r in ribs:
        r.points = r.points + shift
    for x in puncta:
        x["center"] = x["center"] + shift
    center = c0 + shift
    labels, _ = S._sphere_masks(shape, S.SPACING, [center], [label_r])
    lev = float(rng.uniform(*S.NUCLEOPLASM))
    clean = np.full(shape, S.OFFSET) + S._nucleoplasm(shape, S.SPACING, [center], [R], [lev])
    for r in ribs:
        S.render_ribbon(clean, r)
    if with_puncta:
        for x in puncta:
            pts = x["center"][None, :] if x["diam"] <= 0 else _sphere_points(x["diam"]) + x["center"][None, :]
            unit = float(S._value_at(x["center"][None, :], pts, S.PSF_SIGMA)[0])
            S._render_points(clean, pts, S.SPACING, S.PSF_SIGMA, weight=x["peak"] / unit, chunk=2000)
    ce = S.add_noise(clean, rng)
    nuclei = pd.DataFrame([{"nucleus_id": 1, "centroid_z_um": center[0], "centroid_y_um": center[1],
                            "centroid_x_um": center[2], "radius_um": R, "label_radius_um": label_r, "germ_cell": gc,
                            "in_germline": True, "volume_um3": float((labels == 1).sum() * np.prod(S.SPACING))}])
    name = family if with_puncta else family + "_clean"
    return S.Scene(family=name, seed=seed, spacing=S.SPACING, channels={"central_element": ce}, labels=labels,
                   nuclei=nuclei, ribbons=ribs,
                   spots=pd.DataFrame(columns=["spot_id", "nucleus_id", "z_um", "y_um", "x_um"]),
                   truth_foci=pd.DataFrame(columns=["spot_id", "ribbon_id", "nucleus_id", "s_um", "z_um", "y_um",
                                                    "x_um"]),
                   info={"family": name, "seed": seed, "R": R, "nucleoplasm": lev,
                         "puncta": [{"center": [float(v) for v in x["center"]], "diam": x["diam"], "peak": x["peak"]}
                                    for x in puncta] if with_puncta else []})


def _render_modulated(img, rib, sigma, w):
    """synth.render_ribbon with a per-point brightness factor ``w`` (mean 1) on the ribbon's points."""
    P = rib.points
    s = rib.s
    probe_s = np.linspace(0.15 * s[-1], 0.85 * s[-1], 25)
    probe = np.stack([np.interp(probe_s, s, P[:, k]) for k in range(3)], 1)
    unit_peak = float(S._value_at(probe, P, sigma).max())
    if unit_peak <= 0:
        return
    wv = (rib.peak / unit_peak) * np.asarray(w, float)
    for a in range(0, len(P), 300):              # synth._render_points with a weight per point, chunk by chunk
        S._render_points(img, P[a:a + 300], S.SPACING, sigma, weight=wv[a:a + 300], chunk=300)


def real_dense_scene(seed: int, *, final: bool = False, features=REAL_FEATURES, name: str = "real_dense",
                     calib: dict | None = None, ribbons: bool = True) -> S.Scene:
    """The realism-probe nucleus (entry 19). ``calib`` replaces values of REAL_DENSE_CALIB, and ``ribbons=False``
    draws the same scene without its ribbons and puncta (nucleoplasm, texture and noise only; the texture check of
    entry 28 (b)); every random draw is the same either way, so the defaults give the scenes of entry 19."""
    seed = _dev_seed(seed, final)
    code = FAMILY_CODE["real_dense"]
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), code]))
    streams = {k: np.random.default_rng(np.random.SeedSequence([int(seed), code, i]))
               for i, k in enumerate(("psf", "brightness", "texture", "puncta", "noise"), start=1)}
    unknown = sorted(set(calib or {}) - set(REAL_DENSE_CALIB))
    if unknown:
        raise ValueError(f"unknown real_dense calibration keys {unknown}; one of {sorted(REAL_DENSE_CALIB)}")
    cal = {**REAL_DENSE_CALIB, **(calib or {})}
    feats = set(features)
    gc = "oocyte" if seed % 2 == 0 else "spermatocyte"
    pk = REAL_PEAKS[seed % 4]
    R = float(rng.uniform(2.2, 2.6))
    label_r = R + float(rng.uniform(-0.2, 0.4))
    c0 = np.zeros(3)
    ribs = S._shell_nucleus(rng, gc, c0, R, pk, 1, 1)
    lev = float(rng.uniform(*S.NUCLEOPLASM))
    sz = float(streams["psf"].choice(cal["psf_sigma_z_values_um"])) if "psf" in feats else S.PSF_SIGMA[0]
    psf = (sz, S.PSF_SIGMA[1], S.PSF_SIGMA[2])
    puncta = []
    n_p = int(streams["puncta"].poisson(cal["puncta_mean_by_tercile"][0 if pk <= 200 else (2 if pk >= 700 else 1)]))
    if "puncta" not in feats or not ribbons:
        n_p = 0
    rib_pts = np.vstack([r.points[::10] for r in ribs])
    tries = 0
    rp = streams["puncta"]
    while len(puncta) < n_p and tries < 5000:
        tries += 1
        d = rp.normal(size=3)
        d /= np.linalg.norm(d)
        q = c0 + d * R * rp.uniform(0.3, 0.95) ** (1 / 3) * 0.95
        if any(np.linalg.norm(q - x["center"]) < 1.2 for x in puncta):
            continue
        if np.min(np.linalg.norm(rib_pts - q, axis=1)) < 0.6:
            continue
        diam = 0.0 if len(puncta) % 2 == 0 else float(rp.uniform(0.25, 0.8))
        puncta.append({"center": q, "diam": diam, "peak": float(np.exp(rp.uniform(np.log(100.0), np.log(1500.0))))})
    shape, shift = S._place([r.points for r in ribs], [c0], [label_r])
    for r in ribs:
        r.points = r.points + shift
    for x in puncta:
        x["center"] = x["center"] + shift
    center = c0 + shift
    labels, _ = S._sphere_masks(shape, S.SPACING, [center], [label_r])
    nuc = nucleoplasm(shape, S.SPACING, [center], [R], [lev], psf)
    clean = np.full(shape, S.OFFSET) + nuc
    if "texture" in feats:
        sv = np.asarray(psf) / np.asarray(S.SPACING)
        field = ndi.gaussian_filter(streams["texture"].standard_normal(shape), sv, mode="nearest")
        twice = ndi.gaussian_filter(field, sv, mode="nearest")
        zz, yy, xx = np.meshgrid(*[np.arange(n) * s for n, s in zip(shape, S.SPACING)], indexing="ij")
        inside = (zz - center[0]) ** 2 + (yy - center[1]) ** 2 + (xx - center[2]) ** 2 <= R ** 2
        scale = np.sqrt(cal["texture_rel_excess_var"]) * lev / max(float(np.std(twice[inside])), 1e-9)
        clean = clean + field * scale * (nuc / lev)
    for r in (ribs if ribbons else []):
        if "brightness" in feats:
            v, ell = cal["brightness_var_log"], cal["brightness_corr_length_um"]
            n = len(r.points)
            rho = float(np.exp(-S.TRUTH_STEP_UM / ell))
            e = streams["brightness"].standard_normal(n)
            G = np.empty(n)
            G[0] = np.sqrt(v) * e[0]
            for i in range(1, n):
                G[i] = rho * G[i - 1] + np.sqrt(v * (1 - rho * rho)) * e[i]
            _render_modulated(clean, r, psf, np.exp(G - v / 2))
        else:
            S.render_ribbon(clean, r, S.SPACING, psf)
    for x in puncta:
        pts = x["center"][None, :] if x["diam"] <= 0 else _sphere_points(x["diam"]) + x["center"][None, :]
        unit = float(S._value_at(x["center"][None, :], pts, psf)[0])
        S._render_points(clean, pts, S.SPACING, psf, weight=x["peak"] / unit, chunk=2000)
    ce = S.add_noise(clean, streams["noise"])
    nuclei = pd.DataFrame([{"nucleus_id": 1, "centroid_z_um": center[0], "centroid_y_um": center[1],
                            "centroid_x_um": center[2], "radius_um": R, "label_radius_um": label_r, "germ_cell": gc,
                            "in_germline": True, "volume_um3": float((labels == 1).sum() * np.prod(S.SPACING))}])
    info = {"family": name, "seed": seed, "scene_peak": pk, "features": sorted(feats), "psf_sigma_um": list(psf),
            "nucleoplasm": lev,
            "puncta": [{"center": [float(v) for v in x["center"]], "diam": x["diam"], "peak": x["peak"]}
                       for x in puncta]}
    if calib:
        info["calib"] = {k: cal[k] for k in sorted(calib)}
    if not ribbons:
        info["ribbons"] = False
    return S.Scene(family=name, seed=seed, spacing=S.SPACING, channels={"central_element": ce}, labels=labels,
                   nuclei=nuclei, ribbons=ribs if ribbons else [],
                   spots=pd.DataFrame(columns=["spot_id", "nucleus_id", "z_um", "y_um", "x_um"]),
                   truth_foci=pd.DataFrame(columns=["spot_id", "ribbon_id", "nucleus_id", "s_um", "z_um", "y_um",
                                                    "x_um"]),
                   info=info)


PUNCTA_PSF_LAYER_DZ = 6.5               # um between the two layers of puncta_psf
PUNCTA_PSF_GRADIENT = 0.02              # the largest relative change of sigma_z per um of depth


def _render_depth(img, pts, weight, sig_of_z, sxy, chunk=100):
    """synth._render_points with the axial sigma of each chunk of ``chunk`` consecutive points taken at the chunk's
    mean depth (``sig_of_z(z)``); ``weight`` a scalar or one weight per point."""
    pts = np.asarray(pts, float)
    w = np.broadcast_to(np.asarray(weight, float), (len(pts),))
    for a in range(0, len(pts), chunk):
        P = pts[a:a + chunk]
        S._render_points(img, P, S.SPACING, (float(sig_of_z(float(P[:, 0].mean()))), sxy, sxy),
                         weight=np.ascontiguousarray(w[a:a + chunk]), chunk=chunk)


def puncta_psf_scene(seed: int, *, final: bool = False) -> S.Scene:
    """The puncta_psf development family (protocol changelog entry 32; see the module docstring)."""
    seed = _dev_seed(seed, final)
    rng = _rng("puncta_psf", seed)
    lv = PSF_LEVELS[seed % 4]
    sxy = float(rng.uniform(0.135, 0.165))
    grad = float(rng.uniform(-PUNCTA_PSF_GRADIENT, PUNCTA_PSF_GRADIENT))
    scale = WAVELENGTHS["foci"] / WAVELENGTHS["central_element"]
    centers, radii, label_r, gcs, kinds, ribs = [], [], [], [], [], []
    puncta, foci_truth = [], []

    def add_puncta(nid, c, R, n, rib_pts):
        tries, made = 0, 0
        while made < n and tries < 5000:
            tries += 1
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            q = c + d * R * rng.uniform(0.3, 0.95) ** (1 / 3) * 0.95
            if any(np.linalg.norm(q - x["center"]) < 1.2 for x in puncta):
                continue
            if len(rib_pts) and np.min(np.linalg.norm(rib_pts - q, axis=1)) < 0.6:
                continue
            diam = 0.0 if made % 2 == 0 else float(rng.uniform(0.25, 0.8))
            puncta.append({"center": q, "diam": diam, "peak": float(np.exp(rng.uniform(np.log(100.0),
                                                                                    np.log(1500.0)))),
                           "nucleus_id": nid})
            made += 1

    def add_foci(nid, c, R):
        for _ in range(int(rng.poisson(FOCI_DENSITY * 4.0 / 3.0 * np.pi * R ** 3))):
            d = rng.normal(size=3)
            d /= np.linalg.norm(d)
            foci_truth.append((nid, c + d * R * rng.uniform() ** (1.0 / 3.0), float(rng.uniform(*FOCI_PEAK_RANGE))))

    for layer in range(2):
        z0 = layer * PUNCTA_PSF_LAYER_DZ
        x = 0.0
        for i in range(3):
            nid = len(centers) + 1
            R = float(rng.uniform(2.2, 2.6))
            c = np.array([z0, 0.0, x])
            gc = "oocyte" if (i + layer + seed) % 2 == 0 else "spermatocyte"
            rs = S._shell_nucleus(rng, gc, c, R, S.SCENE_PEAK, nid, len(ribs) + 1)
            for r in rs:
                r.kind = "dense"
            add_puncta(nid, c, R, 3, np.vstack([r.points[::10] for r in rs]))
            add_foci(nid, c, R)
            centers.append(c)
            radii.append(R)
            label_r.append(R + float(rng.uniform(-0.2, 0.4)))
            gcs.append(gc)
            kinds.append("dense")
            ribs.extend(rs)
            x += 6.0
        nid = len(centers) + 1
        R = float(rng.uniform(2.2, 2.6))
        c = np.array([z0, 0.0, x])
        add_puncta(nid, c, R, 6, np.zeros((0, 3)))
        add_foci(nid, c, R)
        centers.append(c)
        radii.append(R)
        label_r.append(R + float(rng.uniform(-0.2, 0.4)))
        gcs.append("oocyte")
        kinds.append("puncta")
        x += R + 1.0
        prev_r = R + 0.4
        for _j in range(2):
            nid = len(centers) + 1
            rib = S._iso_like(rng, float(rng.uniform(60.0, 90.0)), float(rng.uniform(0.9, 4.0)),
                              float(rng.uniform(300.0, 800.0)), nucleus_id=nid, ribbon_id=len(ribs) + 1,
                              kind="steep_arc")
            cc, Ra = S._enclosing_nucleus([rib])
            x += prev_r + Ra + 1.0
            rib.points = rib.points - cc + np.array([z0, 0.0, x])
            centers.append(np.array([z0, 0.0, x]))
            radii.append(Ra)
            label_r.append(Ra + float(rng.uniform(0.0, 0.4)))
            gcs.append("unknown")
            kinds.append("arc")
            ribs.append(rib)
            prev_r = Ra + 0.4
    shape, shift = S._place([r.points for r in ribs] + [np.array([q["center"] for q in puncta])], centers, label_r)
    for r in ribs:
        r.points = r.points + shift
    for q in puncta:
        q["center"] = q["center"] + shift
    centers = [c + shift for c in centers]
    z_mid = float(shift[0] + 0.5 * PUNCTA_PSF_LAYER_DZ)

    def sig_z(z: float) -> float:
        return lv * (1.0 + grad * (float(z) - z_mid))

    labels, _ = S._sphere_masks(shape, S.SPACING, centers, label_r)
    levels = [float(rng.uniform(*S.NUCLEOPLASM)) for _ in centers]
    nuc = nucleoplasm(shape, S.SPACING, centers, radii, levels, (lv, sxy, sxy))
    clean = np.full(shape, S.OFFSET) + nuc
    for r in ribs:
        P = r.points
        s = r.s
        probe_s = np.linspace(0.15 * s[-1], 0.85 * s[-1], 25)
        probe = np.stack([np.interp(probe_s, s, P[:, k]) for k in range(3)], 1)
        unit = float(S._value_at(probe, P, (sig_z(float(P[:, 0].mean())), sxy, sxy)).max())
        if unit > 0:
            _render_depth(clean, P, r.peak / unit, sig_z, sxy)
    for q in puncta:
        sg = (sig_z(q["center"][0]), sxy, sxy)
        pts = q["center"][None, :] if q["diam"] <= 0 else _sphere_points(q["diam"]) + q["center"][None, :]
        unit = float(S._value_at(q["center"][None, :], pts, sg)[0])
        S._render_points(clean, pts, S.SPACING, sg, weight=q["peak"] / unit, chunk=2000)
    fclean = np.full(shape, S.OFFSET) + nuc
    rows = []
    for k, (nid, pos, pk) in enumerate(foci_truth):
        pos = pos + shift
        S._render_points(fclean, pos[None, :], S.SPACING, (sig_z(pos[0]) * scale, sxy * scale, sxy * scale),
                         weight=pk)
        rows.append({"spot_id": k, "nucleus_id": nid, "z_um": pos[0], "y_um": pos[1], "x_um": pos[2], "peak": pk})
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
    info = {"family": "puncta_psf", "seed": seed, "sigma_z_level": lv, "sigma_z_gradient_per_um": grad,
            "z_mid_um": z_mid, "psf_sigma_xy_um": sxy, "wavelengths": dict(WAVELENGTHS), "nucleoplasm": levels,
            "puncta": [{"center": [float(v) for v in q["center"]], "diam": q["diam"], "peak": q["peak"],
                        "nucleus_id": int(q["nucleus_id"]), "sigma_z_um": sig_z(q["center"][0])} for q in puncta],
            "foci_truth": truth.to_dict("list")}
    return S.Scene(family="puncta_psf", seed=seed, spacing=S.SPACING, channels={"central_element": ce, "foci": fo},
                   labels=labels, nuclei=_nuclei_table(labels, centers, radii, label_r, gcs, kinds), ribbons=ribs,
                   spots=spots, truth_foci=pd.DataFrame(columns=["spot_id", "ribbon_id", "nucleus_id", "s_um", "z_um",
                                                                 "y_um", "x_um"]), info=info)


GAPS_PSF_PER_LAYER = 6                  # ribbons (and nuclei) per layer of gaps_psf: gap, dim, gap, dim, gap, dim


def _events(rng, L: float, kind: str):
    """Positions of 1 or 2 events along a ribbon of length L as synth.gaps_scene places its gaps (at least 0.8 um of
    ribbon before, between and after them; the ribbon lengthened by 0.5 um steps until they fit): empty gaps of
    U(0.3, 1.0) um, or dim stretches of U(0.4, 1.5) um. Returns (L, [(s0, s1), ...])."""
    n_ev = 1 if L < 4.0 else int(rng.integers(1, 3))
    lo, hi = (0.3, 1.0) if kind == "gap" else (0.4, 1.5)
    widths = [float(rng.uniform(lo, hi)) for _ in range(n_ev)]
    min_piece = 0.8
    while sum(widths) + (n_ev + 1) * min_piece > L:
        L += 0.5
    free = L - sum(widths) - (n_ev + 1) * min_piece
    cuts = np.sort(rng.uniform(0, free, n_ev))
    ev, pos, prev_cut = [], 0.0, 0.0
    for k in range(n_ev):
        pos += min_piece + (cuts[k] - prev_cut)
        prev_cut = cuts[k]
        ev.append((float(pos), float(pos + widths[k])))
        pos += widths[k]
    return L, ev


def gaps_psf_scene(seed: int, *, final: bool = False) -> S.Scene:
    """The gaps_psf development family (protocol changelog entry 33; see the module docstring)."""
    seed = _dev_seed(seed, final)
    rng = _rng("gaps_psf", seed)
    lv = PSF_LEVELS[seed % 4]
    sxy = float(rng.uniform(0.135, 0.165))
    grad = float(rng.uniform(-PUNCTA_PSF_GRADIENT, PUNCTA_PSF_GRADIENT))
    scale = WAVELENGTHS["foci"] / WAVELENGTHS["central_element"]
    centers, radii, label_r, gcs, kinds, ribs, foci_truth = [], [], [], [], [], [], []
    dims: dict[int, list] = {}
    for layer in range(2):
        z0 = layer * PUNCTA_PSF_LAYER_DZ
        x = 0.0
        prev_r = 0.0
        for j in range(GAPS_PSF_PER_LAYER):
            kind = "gap" if j % 2 == 0 else "dim"
            L, ev = _events(rng, float(rng.uniform(3.0, 7.0)), kind)
            nid = len(centers) + 1
            rib = S._iso_like(rng, float(rng.uniform(0.0, 75.0)), L, float(rng.uniform(300.0, 800.0)),
                              nucleus_id=nid, ribbon_id=len(ribs) + 1, kind=f"{kind}_ribbon")
            if kind == "gap":
                rib.gaps = list(ev)
            else:
                dims[rib.ribbon_id] = [[a, b, float(rng.uniform(0.15, 0.6))] for a, b in ev]
            cc, Ra = S._enclosing_nucleus([rib])
            x += prev_r + Ra + 1.0
            c = np.array([z0, 0.0, x])
            rib.points = rib.points - cc + c
            for _ in range(int(rng.poisson(FOCI_DENSITY * 4.0 / 3.0 * np.pi * Ra ** 3))):
                d = rng.normal(size=3)
                d /= np.linalg.norm(d)
                foci_truth.append((nid, c + d * Ra * rng.uniform() ** (1.0 / 3.0), float(rng.uniform(*FOCI_PEAK_RANGE))))
            centers.append(c)
            radii.append(Ra)
            label_r.append(Ra + float(rng.uniform(0.0, 0.4)))
            gcs.append("unknown")
            kinds.append(kind)
            ribs.append(rib)
            prev_r = Ra + 0.4
    shape, shift = S._place([r.points for r in ribs], centers, label_r)
    for r in ribs:
        r.points = r.points + shift
    centers = [c + shift for c in centers]
    z_mid = float(shift[0] + 0.5 * PUNCTA_PSF_LAYER_DZ)

    def sig_z(z: float) -> float:
        return lv * (1.0 + grad * (float(z) - z_mid))

    labels, _ = S._sphere_masks(shape, S.SPACING, centers, label_r)
    levels = [float(rng.uniform(*S.NUCLEOPLASM)) for _ in centers]
    nuc = nucleoplasm(shape, S.SPACING, centers, radii, levels, (lv, sxy, sxy))
    clean = np.full(shape, S.OFFSET) + nuc
    for r in ribs:
        s = r.s
        w = np.ones(len(s))
        events = [(a, b) for a, b in r.gaps] + [(a, b) for a, b, _f in dims.get(r.ribbon_id, [])]
        for a, b in r.gaps:
            w[(s > a) & (s < b)] = 0.0
        for a, b, f in dims.get(r.ribbon_id, []):
            w[(s >= a) & (s <= b)] = f
        keep = w > 0
        P = r.points[keep]
        probe_s = np.linspace(0.15 * s[-1], 0.85 * s[-1], 25)
        ok = np.ones(len(probe_s), bool)
        for a, b in events:
            ok &= (probe_s < a - 0.5) | (probe_s > b + 0.5)
        probe_s = probe_s[ok] if ok.any() else probe_s
        probe = np.stack([np.interp(probe_s, s, r.points[:, k]) for k in range(3)], 1)
        unit = float(S._value_at(probe, P, (sig_z(float(P[:, 0].mean())), sxy, sxy)).max())
        if unit > 0:
            _render_depth(clean, P, w[keep] * (r.peak / unit), sig_z, sxy)
    fclean = np.full(shape, S.OFFSET) + nuc
    rows = []
    for k, (nid, pos, pk) in enumerate(foci_truth):
        pos = pos + shift
        S._render_points(fclean, pos[None, :], S.SPACING, (sig_z(pos[0]) * scale, sxy * scale, sxy * scale),
                         weight=pk)
        rows.append({"spot_id": k, "nucleus_id": nid, "z_um": pos[0], "y_um": pos[1], "x_um": pos[2], "peak": pk})
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
    info = {"family": "gaps_psf", "seed": seed, "sigma_z_level": lv, "sigma_z_gradient_per_um": grad,
            "z_mid_um": z_mid, "psf_sigma_xy_um": sxy, "wavelengths": dict(WAVELENGTHS), "nucleoplasm": levels,
            "dim_stretches": {int(k): v for k, v in dims.items()}, "foci_truth": truth.to_dict("list")}
    return S.Scene(family="gaps_psf", seed=seed, spacing=S.SPACING, channels={"central_element": ce, "foci": fo},
                   labels=labels, nuclei=_nuclei_table(labels, centers, radii, label_r, gcs, kinds), ribbons=ribs,
                   spots=spots, truth_foci=pd.DataFrame(columns=["spot_id", "ribbon_id", "nucleus_id", "s_um", "z_um",
                                                                 "y_um", "x_um"]), info=info)


def puncta_psf_sigma_z(info: dict, z: float) -> float:
    """The axial PSF sigma of a puncta_psf or gaps_psf scene at depth ``z`` (um, image frame)."""
    return float(info["sigma_z_level"]) * (1.0 + float(info["sigma_z_gradient_per_um"]) * (float(z) -
                                                                                         float(info["z_mid_um"])))


def _real(name, features, calib=None):
    def f(seed, *, final=False):
        return real_dense_scene(seed, final=final, features=features, name=name, calib=calib)
    return f


def _puncta(family):
    def f(seed, *, final=False):
        return puncta_scene(seed, final=final, family=family)
    return f


def _dense_puncta_clean(seed, *, final=False):
    return puncta_scene(seed, final=final, family="dense_puncta", with_puncta=False)


FAMILY_FUNCS = {"psf_check": psf_check_scene, "puncta_psf": puncta_psf_scene, "gaps_psf": gaps_psf_scene,
                "puncta_pt": _puncta("puncta_pt"), "puncta_agg": _puncta("puncta_agg"),
                "dense_puncta": _puncta("dense_puncta"), "dense_puncta_clean": _dense_puncta_clean,
                "real_dense": _real("real_dense", REAL_FEATURES), "real_dense_base": _real("real_dense_base", ()),
                **{f"real_dense_{k}": _real(f"real_dense_{k}", (k,)) for k in REAL_FEATURES},
                # entry 28 (d): the same scenes with the recalibrated features
                "real_dense_r2": _real("real_dense_r2", REAL_FEATURES, REAL_DENSE_CALIB_R2),
                "real_dense_r2_base": _real("real_dense_r2_base", (), REAL_DENSE_CALIB_R2),
                **{f"real_dense_r2_{k}": _real(f"real_dense_r2_{k}", (k,), REAL_DENSE_CALIB_R2) for k in REAL_FEATURES}}


def make_scene(family: str, seed: int, *, final: bool = False, **kw) -> S.Scene:
    _dev_seed(seed, final)
    if family not in FAMILY_FUNCS:
        raise ValueError(f"unknown development family {family!r}; one of {sorted(FAMILY_FUNCS)}")
    return FAMILY_FUNCS[family](seed, final=final, **kw)
