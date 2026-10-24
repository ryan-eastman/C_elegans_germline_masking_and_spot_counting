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
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from . import synth as S

FAMILY_CODE = {"puncta_pt": 101, "puncta_agg": 102, "dense_puncta": 103, "psf_check": 107, "real_dense": 108}
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


def real_dense_scene(seed: int, *, final: bool = False, features=REAL_FEATURES, name: str = "real_dense") -> S.Scene:
    seed = _dev_seed(seed, final)
    code = FAMILY_CODE["real_dense"]
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), code]))
    streams = {k: np.random.default_rng(np.random.SeedSequence([int(seed), code, i]))
               for i, k in enumerate(("psf", "brightness", "texture", "puncta", "noise"), start=1)}
    cal = REAL_DENSE_CALIB
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
    if "puncta" not in feats:
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
    for r in ribs:
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
    return S.Scene(family=name, seed=seed, spacing=S.SPACING, channels={"central_element": ce}, labels=labels,
                   nuclei=nuclei, ribbons=ribs,
                   spots=pd.DataFrame(columns=["spot_id", "nucleus_id", "z_um", "y_um", "x_um"]),
                   truth_foci=pd.DataFrame(columns=["spot_id", "ribbon_id", "nucleus_id", "s_um", "z_um", "y_um",
                                                    "x_um"]),
                   info={"family": name, "seed": seed, "scene_peak": pk, "features": sorted(feats),
                         "psf_sigma_um": list(psf), "nucleoplasm": lev,
                         "puncta": [{"center": [float(v) for v in x["center"]], "diam": x["diam"], "peak": x["peak"]}
                                    for x in puncta]})


def _real(name, features):
    def f(seed, *, final=False):
        return real_dense_scene(seed, final=final, features=features, name=name)
    return f


def _puncta(family):
    def f(seed, *, final=False):
        return puncta_scene(seed, final=final, family=family)
    return f


def _dense_puncta_clean(seed, *, final=False):
    return puncta_scene(seed, final=final, family="dense_puncta", with_puncta=False)


FAMILY_FUNCS = {"psf_check": psf_check_scene, "puncta_pt": _puncta("puncta_pt"), "puncta_agg": _puncta("puncta_agg"),
                "dense_puncta": _puncta("dense_puncta"), "dense_puncta_clean": _dense_puncta_clean,
                "real_dense": _real("real_dense", REAL_FEATURES), "real_dense_base": _real("real_dense_base", ()),
                **{f"real_dense_{k}": _real(f"real_dense_{k}", (k,)) for k in REAL_FEATURES}}


def make_scene(family: str, seed: int, *, final: bool = False, **kw) -> S.Scene:
    _dev_seed(seed, final)
    if family not in FAMILY_FUNCS:
        raise ValueError(f"unknown development family {family!r}; one of {sorted(FAMILY_FUNCS)}")
    return FAMILY_FUNCS[family](seed, final=final, **kw)
