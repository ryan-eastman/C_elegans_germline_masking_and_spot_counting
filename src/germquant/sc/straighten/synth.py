"""Exact-render SC phantoms with known truth (docs/SC_VALIDATION_PROTOCOL.md, section 2.3: T3).

Every ribbon is a dense truth polyline (points 0.005 um apart) in the pipeline frame (um = voxel index x
spacing, axis order z, y, x). The image is the sum over truth points of a Gaussian PSF evaluated at the voxel
centres, so the truth is never rounded to voxels; each ribbon is scaled so that its peak above the
nucleoplasm, after the PSF, is the requested value. Noise is Gaussian with the variance model measured on
the development gonad n2rad51_noHS_herm_000 (4.04 x I - 394 counts^2 at the noise-free absolute intensity I,
offset included, floored at 25 counts^2), and the result is rounded to uint16 like a camera frame.

Fixed distributions (the protocol registers them; a change is a protocol edit):

* grid 0.2 x 0.108333 x 0.108333 um; PSF sigma_xy 0.149 um, sigma_z 0.365 um;
* camera offset 100 counts; nucleoplasm U(30, 90) counts above offset inside the nucleus (PSF-blurred);
* ribbon peak = scene peak x U(0.7, 1.3), scene peak 430 counts unless the family says otherwise;
* oocyte nuclei: radius U(2.2, 2.6) um, 6 ribbons of length U(4.5, 7.0) um, one of them U(4.5, 5.7) um (the X);
  spermatocyte nuclei: 5 ribbons; ribbons on spheres between 0.55 and 0.95 of the radius, curvature radius at
  least 1.0 um; label = sphere of radius R + U(-0.2, +0.4) um.

Families: ``iso``, ``dense``, ``gaps``, ``foci``, ``pairs``, ``noise``, ``snr`` (protocol section 6.1), plus
``phantom_gonad`` for pipeline tests. Development may use any seed below 7,000,000; a seed at or above it
belongs to the sealed ranges and is refused unless ``final=True`` (only the scoring script passes it).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

SPACING = (0.2, 0.108333, 0.108333)
PSF_SIGMA = (0.365, 0.149, 0.149)            # (z, y, x) um
OFFSET = 100.0
NUCLEOPLASM = (30.0, 90.0)
NOISE_SLOPE, NOISE_INTERCEPT, NOISE_FLOOR = 4.04, -394.0, 25.0
SCENE_PEAK = 430.0
PEAK_JITTER = (0.7, 1.3)
TRUTH_STEP_UM = 0.005
FOCI_PEAK = 300.0
FOCI_LOC_SIGMA = (0.06, 0.03, 0.03)          # (z, y, x) localisation noise of the synthetic spot table
DEV_SEED_LIMIT = 7_000_000
FAMILIES = ("iso", "dense", "gaps", "foci", "pairs", "noise", "snr", "phantom_gonad")
_FAMILY_CODE = {name: i for i, name in enumerate(FAMILIES, start=1)}
PAIR_KINDS = ("x_crossing", "t_contact", "parallel_xy", "stacked_z")
SNR_PEAKS = (200.0, 300.0, 430.0, 700.0)


class SealedSeedError(PermissionError):
    """A sealed synthetic seed (7,000,000 or above) was asked for without ``final=True``."""

    __test__ = False


@dataclass
class TruthRibbon:
    ribbon_id: int
    nucleus_id: int
    points: np.ndarray                     # (n, 3) um, dense (0.005 um), ordered
    gaps: list = field(default_factory=list)        # [(s0, s1)] um along the ribbon with no signal
    foci_s: list = field(default_factory=list)      # s positions (um) of the foci on this ribbon
    peak: float = 0.0
    tilt_deg: float = float("nan")
    kind: str = ""

    @property
    def s(self) -> np.ndarray:
        d = np.linalg.norm(np.diff(self.points, axis=0), axis=1)
        return np.concatenate([[0.0], np.cumsum(d)])

    @property
    def length_um(self) -> float:
        return float(self.s[-1])

    def point_at(self, s: float) -> np.ndarray:
        ss = self.s
        return np.array([np.interp(s, ss, self.points[:, k]) for k in range(3)])


@dataclass
class Scene:
    family: str
    seed: int
    spacing: tuple
    channels: dict                         # role -> uint16 (Z, Y, X)
    labels: np.ndarray                     # uint16 (Z, Y, X)
    nuclei: pd.DataFrame
    ribbons: list
    spots: pd.DataFrame                    # the spot table a method gets (localisation noise added)
    truth_foci: pd.DataFrame
    info: dict = field(default_factory=dict)

    @property
    def shape(self) -> tuple:
        return tuple(self.labels.shape)


# ----------------------------------------------------------------------------------------- seeds
def check_seed(seed: int, final: bool = False) -> int:
    seed = int(seed)
    if seed < 0:
        raise ValueError("synthetic seeds are non-negative")
    if seed >= DEV_SEED_LIMIT and not final:
        raise SealedSeedError(
            f"seed {seed} is in the sealed synthetic ranges (>= {DEV_SEED_LIMIT}); development uses seeds below "
            "7,000,000 and only scripts/benchmark_sc.py synthetic --final generates sealed seeds "
            "(docs/SC_VALIDATION_PROTOCOL.md, section 2.3)")
    return seed


def _rng(family: str, seed: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([int(seed), _FAMILY_CODE[family]]))


# ------------------------------------------------------------------------------------ geometry
def _densify(P: np.ndarray, step: float = TRUTH_STEP_UM) -> np.ndarray:
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
    n = max(2, int(np.ceil(s[-1] / step)) + 1)
    q = np.linspace(0.0, s[-1], n)
    return np.stack([np.interp(q, s, P[:, k]) for k in range(3)], 1)


def _rot_y(t):
    c, s = np.cos(t), np.sin(t)
    # rotate the x axis toward z (axis order z, y, x)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def _rot_z(t):
    c, s = np.cos(t), np.sin(t)
    # rotation in the y-x plane about the optical axis
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def arc(L: float, radius: float, tilt_deg: float = 0.0, azimuth: float = 0.0) -> np.ndarray:
    """A planar circular arc of length L and curvature radius ``radius`` centred at the origin: chord along x,
    bow along y, then the chord tilted by ``tilt_deg`` toward z and rotated by ``azimuth`` about z. Dense."""
    radius = max(float(radius), L / (0.9 * np.pi))
    s = np.linspace(-L / 2, L / 2, int(np.ceil(L / TRUTH_STEP_UM)) + 1)
    ang = s / radius
    x = radius * np.sin(ang)
    y = radius * (1.0 - np.cos(ang))
    P = np.stack([np.zeros_like(x), y - y.mean(), x], 1)
    P = P @ _rot_y(np.radians(tilt_deg)).T @ _rot_z(azimuth).T
    return P


def shell_ribbon(rng, center, R: float, L: float, start_dir) -> np.ndarray:
    """A ribbon on a sphere of radius rho in U(0.55 R, 0.95 R) around ``center``, starting in ``start_dir``,
    with a smooth geodesic curvature; total curvature sqrt(1/rho^2 + kg^2) is at most 1.0 per um."""
    rho = rng.uniform(0.55, 0.95) * R
    n = np.asarray(start_dir, float)
    n /= np.linalg.norm(n)
    t = rng.normal(size=3)
    t -= (t @ n) * n
    t /= np.linalg.norm(t)
    k0 = rng.uniform(-0.3, 0.3)
    k1 = rng.uniform(0.0, 0.3)
    lam = rng.uniform(2.0, 5.0)
    ph = rng.uniform(0, 2 * np.pi)
    kmax = np.sqrt(max(1.0 - 1.0 / rho ** 2, 0.0))
    ds = TRUTH_STEP_UM
    m = int(np.ceil(L / ds)) + 1
    pts = np.empty((m, 3))
    for i in range(m):
        pts[i] = center + rho * n
        s = i * ds
        kg = float(np.clip(k0 + k1 * np.sin(2 * np.pi * s / lam + ph), -kmax, kmax))
        # move along the sphere: n turns toward t at rate 1/rho; t turns about n at rate kg
        n_new = n + t * (ds / rho)
        n_new /= np.linalg.norm(n_new)
        t = t - (t @ n_new) * n_new
        b = np.cross(n_new, t)
        t = t * np.cos(kg * ds) + b * np.sin(kg * ds)
        t /= np.linalg.norm(t)
        n = n_new
    return pts


def fib_directions(k: int, rng) -> np.ndarray:
    i = np.arange(k) + 0.5
    phi = np.arccos(np.clip(1 - 2 * i / k, -1, 1))
    theta = np.pi * (1 + 5 ** 0.5) * i
    d = np.stack([np.cos(phi), np.sin(phi) * np.sin(theta), np.sin(phi) * np.cos(theta)], 1)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return d @ q


def curvature(points: np.ndarray) -> np.ndarray:
    """Discrete curvature (1/um) of a dense polyline (central differences)."""
    P = np.asarray(points, float)
    d = np.diff(P, axis=0)
    ds = np.linalg.norm(d, axis=1)
    t = d / ds[:, None]
    dt = np.linalg.norm(np.diff(t, axis=0), axis=1)
    return dt / (0.5 * (ds[1:] + ds[:-1]))


# ------------------------------------------------------------------------------------ rendering
def _render_points(img: np.ndarray, pts: np.ndarray, spacing, sigma, weight: float = 1.0, chunk: int = 300):
    """Add sum_p weight * exp(-|x - p|^2 / 2 sigma^2) (unit peak per point) at every voxel centre of ``img``."""
    pts = np.asarray(pts, float)
    if not len(pts):
        return
    sp = np.asarray(spacing, float)
    sg = np.asarray(sigma, float)
    shape = np.asarray(img.shape)
    for a in range(0, len(pts), chunk):
        P = pts[a:a + chunk]
        lo = np.floor((P.min(0) - 4.0 * sg) / sp).astype(int)
        hi = np.ceil((P.max(0) + 4.0 * sg) / sp).astype(int) + 1
        lo = np.clip(lo, 0, shape)
        hi = np.clip(hi, 0, shape)
        if np.any(hi <= lo):
            continue
        axes = [np.arange(lo[k], hi[k]) * sp[k] for k in range(3)]
        G = [np.exp(-0.5 * ((axes[k][None, :] - P[:, k:k + 1]) / sg[k]) ** 2) for k in range(3)]
        A = (G[1][:, :, None] * G[2][:, None, :]).reshape(len(P), -1)
        block = (G[0].T * weight) @ A
        img[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] += block.reshape(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2])


def _value_at(points_eval: np.ndarray, pts: np.ndarray, sigma) -> np.ndarray:
    """The continuous unit render of ``pts`` evaluated at ``points_eval`` (exact, no voxels)."""
    sg = np.asarray(sigma, float)
    out = np.zeros(len(points_eval))
    for a in range(0, len(pts), 2000):
        P = pts[a:a + 2000]
        d = (points_eval[:, None, :] - P[None, :, :]) / sg
        out += np.exp(-0.5 * (d ** 2).sum(-1)).sum(1)
    return out


def _signal_points(rib: TruthRibbon) -> np.ndarray:
    """The truth points that carry signal (gaps removed)."""
    if not rib.gaps:
        return rib.points
    s = rib.s
    keep = np.ones(len(s), bool)
    for g0, g1 in rib.gaps:
        keep &= ~((s > g0) & (s < g1))
    return rib.points[keep]


def render_ribbon(img, rib: TruthRibbon, spacing=SPACING, sigma=PSF_SIGMA):
    """Add one ribbon scaled to its peak (``rib.peak`` counts above the nucleoplasm, after the PSF)."""
    P = _signal_points(rib)
    s = rib.s
    probe_s = np.linspace(0.15 * s[-1], 0.85 * s[-1], 25)
    if rib.gaps:
        ok = np.ones(len(probe_s), bool)
        for g0, g1 in rib.gaps:
            ok &= (probe_s < g0 - 0.5) | (probe_s > g1 + 0.5)
        probe_s = probe_s[ok] if ok.any() else probe_s
    probe = np.stack([np.interp(probe_s, s, rib.points[:, k]) for k in range(3)], 1)
    unit_peak = float(_value_at(probe, P, sigma).max())
    if unit_peak <= 0:
        return
    _render_points(img, P, spacing, sigma, weight=rib.peak / unit_peak)


def add_noise(clean_abs: np.ndarray, rng) -> np.ndarray:
    var = np.maximum(NOISE_SLOPE * clean_abs + NOISE_INTERCEPT, NOISE_FLOOR)
    noisy = clean_abs + np.sqrt(var) * rng.standard_normal(clean_abs.shape)
    return np.clip(np.rint(noisy), 0, 65535).astype(np.uint16)


def _grid_axes(shape, spacing):
    return [np.arange(n) * s for n, s in zip(shape, spacing)]


def _sphere_masks(shape, spacing, centers, radii) -> tuple[np.ndarray, np.ndarray]:
    """Nearest-centre (Voronoi) assignment inside each sphere -> uint16 labels (1..n), and the distance ratio."""
    zz, yy, xx = np.meshgrid(*_grid_axes(shape, spacing), indexing="ij")
    labels = np.zeros(shape, np.uint16)
    best = np.full(shape, np.inf)
    for i, (c, r) in enumerate(zip(centers, radii), start=1):
        d = np.sqrt((zz - c[0]) ** 2 + (yy - c[1]) ** 2 + (xx - c[2]) ** 2)
        rel = d / r
        take = (rel <= 1.0) & (d < best)
        labels[take] = i
        best = np.where(take, d, best)
    return labels, best


def _nucleoplasm(shape, spacing, centers, radii, levels) -> np.ndarray:
    zz, yy, xx = np.meshgrid(*_grid_axes(shape, spacing), indexing="ij")
    out = np.zeros(shape, np.float64)
    for c, r, lev in zip(centers, radii, levels):
        inside = (zz - c[0]) ** 2 + (yy - c[1]) ** 2 + (xx - c[2]) ** 2 <= r ** 2
        out += lev * inside
    sig = np.asarray(PSF_SIGMA) / np.asarray(spacing)
    return ndi.gaussian_filter(out, sig, mode="nearest")


# ------------------------------------------------------------------------------------ assembly
def _finish(family, seed, rng, shape, centers, radii, label_radii, germ_cells, ribbons, foci_truth=None,
            peak_info=None, foci_channel=False) -> Scene:
    spacing = SPACING
    labels, _ = _sphere_masks(shape, spacing, centers, label_radii)
    levels = [rng.uniform(*NUCLEOPLASM) for _ in centers]
    nuc = _nucleoplasm(shape, spacing, centers, radii, levels)
    clean = np.full(shape, OFFSET) + nuc
    for rib in ribbons:
        render_ribbon(clean, rib)
    channels = {"central_element": add_noise(clean, rng)}
    truth_foci = pd.DataFrame(columns=["spot_id", "ribbon_id", "nucleus_id", "s_um", "z_um", "y_um", "x_um"])
    spots = pd.DataFrame(columns=["spot_id", "nucleus_id", "z_um", "y_um", "x_um"])
    if foci_channel:
        fclean = np.full(shape, OFFSET) + nuc
        rows = []
        for rib in ribbons:
            for s in rib.foci_s:
                p = rib.point_at(s)
                rows.append({"spot_id": len(rows), "ribbon_id": rib.ribbon_id, "nucleus_id": rib.nucleus_id,
                             "s_um": float(s), "z_um": p[0], "y_um": p[1], "x_um": p[2]})
                _render_points(fclean, p[None, :], spacing, PSF_SIGMA, weight=FOCI_PEAK)
        channels["foci"] = add_noise(fclean, rng)
        truth_foci = pd.DataFrame(rows, columns=truth_foci.columns)
        if rows:
            noise = rng.normal(size=(len(rows), 3)) * np.asarray(FOCI_LOC_SIGMA)
            spots = pd.DataFrame({"spot_id": truth_foci["spot_id"].to_numpy(),
                                  "nucleus_id": truth_foci["nucleus_id"].to_numpy(),
                                  "z_um": truth_foci["z_um"].to_numpy() + noise[:, 0],
                                  "y_um": truth_foci["y_um"].to_numpy() + noise[:, 1],
                                  "x_um": truth_foci["x_um"].to_numpy() + noise[:, 2]})
    vol_vox = np.prod(spacing)
    nrows = []
    for i, (c, r, lr, gc) in enumerate(zip(centers, radii, label_radii, germ_cells), start=1):
        nrows.append({"nucleus_id": i, "centroid_z_um": c[0], "centroid_y_um": c[1], "centroid_x_um": c[2],
                      "radius_um": r, "label_radius_um": lr, "germ_cell": gc, "in_germline": True,
                      "volume_um3": float((labels == i).sum() * vol_vox), "n_voxels": int((labels == i).sum())})
    info = {"family": family, "seed": seed, "nucleoplasm": levels, "offset": OFFSET,
            "noise": [NOISE_SLOPE, NOISE_INTERCEPT, NOISE_FLOOR], "psf_sigma_um": list(PSF_SIGMA)}
    info.update(peak_info or {})
    return Scene(family=family, seed=seed, spacing=spacing, channels=channels, labels=labels,
                 nuclei=pd.DataFrame(nrows), ribbons=ribbons, spots=spots, truth_foci=truth_foci, info=info)


def _box_for(points_list, centers, label_radii, margin=1.5):
    allp = [np.asarray(c)[None, :] + np.array([[-r, -r, -r], [r, r, r]]) for c, r in zip(centers, label_radii)]
    allp += [p for p in points_list if len(p)]
    P = np.vstack(allp)
    lo = P.min(0) - margin
    hi = P.max(0) + margin
    return lo, hi


def _place(points_list, centers, label_radii, margin=1.5):
    """Shift everything so the box starts at a whole voxel near 0; return (shape, shift)."""
    sp = np.asarray(SPACING)
    lo, hi = _box_for(points_list, centers, label_radii, margin)
    lo_vox = np.floor(lo / sp)
    shift = -lo_vox * sp
    shape = tuple(int(v) for v in np.ceil((hi - lo) / sp) + 2)
    return shape, shift


MIN_RIBBON_SEPARATION_UM = 0.25     # two chromosomes' SCs never pass through each other (about one SC width)


def _shell_nucleus(rng, germ_cell: str, center, R: float, peak: float, nucleus_id: int, first_id: int,
                   max_tries: int = 60):
    """The ribbons of one nucleus. A ribbon that comes closer than MIN_RIBBON_SEPARATION_UM (3D) to a ribbon
    already placed is drawn again (up to ``max_tries`` times; then the draw with the largest clearance is kept)."""
    from scipy.spatial import cKDTree

    n_rib = 6 if germ_cell == "oocyte" else 5
    lengths = [rng.uniform(4.5, 7.0) for _ in range(n_rib)]
    if germ_cell == "oocyte":
        lengths[-1] = rng.uniform(4.5, 5.7)          # the X
    dirs = fib_directions(n_rib, rng)
    ribs = []
    placed = None
    for k in range(n_rib):
        best, best_d = None, -1.0
        for attempt in range(max_tries):
            start = dirs[k] if attempt == 0 else rng.normal(size=3)
            P = shell_ribbon(rng, np.asarray(center, float), R, lengths[k], start)
            d = float(placed.query(P[::4])[0].min()) if placed is not None else np.inf
            if d > best_d:
                best, best_d = P, d
            if d >= MIN_RIBBON_SEPARATION_UM:
                break
        ribs.append(TruthRibbon(ribbon_id=first_id + k, nucleus_id=nucleus_id, points=best,
                                peak=peak * rng.uniform(*PEAK_JITTER), kind="shell"))
        allp = np.vstack([r.points[::2] for r in ribs])
        placed = cKDTree(allp)
    return ribs


def _single_nucleus_scene(family, seed, rng, ribbons_fn, germ_cell="oocyte", R=None, foci=False, peak_info=None,
                          center_ribbons=True):
    """One nucleus at the box centre; ``ribbons_fn(center, R)`` returns the truth ribbons."""
    if R is None:
        R = rng.uniform(2.2, 2.6)
    label_r = R + rng.uniform(-0.2, 0.4)
    center0 = np.zeros(3)
    ribs = ribbons_fn(center0, R)
    shape, shift = _place([r.points for r in ribs], [center0], [label_r])
    for r in ribs:
        r.points = r.points + shift
    center = center0 + shift
    return _finish(family, seed, rng, shape, [center], [R], [label_r], [germ_cell], ribs, foci_channel=foci,
                   peak_info=peak_info)


def _iso_like(rng, tilt_deg, L, peak, nucleus_id=1, ribbon_id=1, kind="iso"):
    Rc = rng.uniform(1.2, 4.0)
    P = arc(L, Rc, tilt_deg, rng.uniform(0, 2 * np.pi))
    return TruthRibbon(ribbon_id=ribbon_id, nucleus_id=nucleus_id, points=P, peak=peak, tilt_deg=tilt_deg, kind=kind)


def _enclosing_nucleus(ribs, pad=0.6):
    P = np.vstack([r.points for r in ribs])
    c = 0.5 * (P.min(0) + P.max(0))
    R = float(np.linalg.norm(P - c, axis=1).max() + pad)
    return c, R


def _scene_around(family, seed, rng, ribs, foci=False, peak_info=None, germ_cell="oocyte"):
    c, R = _enclosing_nucleus(ribs)
    for r in ribs:
        r.points = r.points - c
    label_r = R + rng.uniform(0.0, 0.4)
    center0 = np.zeros(3)
    shape, shift = _place([r.points for r in ribs], [center0], [label_r])
    for r in ribs:
        r.points = r.points + shift
    return _finish(family, seed, rng, shape, [center0 + shift], [R], [label_r], [germ_cell], ribs,
                   foci_channel=foci, peak_info=peak_info)


# ------------------------------------------------------------------------------------ families
def iso_scene(seed: int, *, final: bool = False, steep: bool | None = None, L: float | None = None,
              peak: float | None = None, tilt_deg: float | None = None) -> Scene:
    seed = check_seed(seed, final)
    rng = _rng("iso", seed)
    if steep is None:
        steep = seed % 5 == 4
    tilt = rng.uniform(60.0, 75.0) if steep else rng.uniform(0.0, 60.0)
    if tilt_deg is not None:
        tilt = float(tilt_deg)
    L_ = rng.uniform(2.0, 8.0) if L is None else float(L)
    pk = rng.uniform(300.0, 800.0) if peak is None else float(peak)
    rib = _iso_like(rng, tilt, L_, pk)
    return _scene_around("iso", seed, rng, [rib], peak_info={"steep": bool(steep)})


def gaps_scene(seed: int, *, final: bool = False) -> Scene:
    seed = check_seed(seed, final)
    rng = _rng("gaps", seed)
    L = rng.uniform(2.0, 8.0)
    n_gaps = 1 if L < 4.0 else int(rng.integers(1, 3))
    widths = [rng.uniform(0.3, 1.0) for _ in range(n_gaps)]
    min_piece = 0.8
    while sum(widths) + (n_gaps + 1) * min_piece > L:
        L += 0.5
    rib = _iso_like(rng, rng.uniform(0.0, 60.0), L, SCENE_PEAK * rng.uniform(*PEAK_JITTER), kind="gaps")
    free = L - sum(widths) - (n_gaps + 1) * min_piece
    cuts = np.sort(rng.uniform(0, free, n_gaps))
    gaps, pos = [], 0.0
    prev_cut = 0.0
    for k in range(n_gaps):
        pos += min_piece + (cuts[k] - prev_cut)
        prev_cut = cuts[k]
        gaps.append((pos, pos + widths[k]))
        pos += widths[k]
    rib.gaps = gaps
    return _scene_around("gaps", seed, rng, [rib])


def foci_scene(seed: int, *, final: bool = False) -> Scene:
    seed = check_seed(seed, final)
    rng = _rng("foci", seed)
    L = rng.uniform(2.0, 8.0)
    rib = _iso_like(rng, rng.uniform(0.0, 60.0), L, SCENE_PEAK * rng.uniform(*PEAK_JITTER), kind="foci")
    n = int(rng.integers(1, 4))
    pos: list[float] = []
    for _ in range(200):
        if len(pos) == n:
            break
        s = rng.uniform(0.4, L - 0.4)
        if all(abs(s - q) >= 0.6 for q in pos):
            pos.append(float(s))
    rib.foci_s = sorted(pos)
    return _scene_around("foci", seed, rng, [rib], foci=True)


def pairs_scene(seed: int, *, final: bool = False, kind: str | None = None) -> Scene:
    seed = check_seed(seed, final)
    rng = _rng("pairs", seed)
    kind = kind or PAIR_KINDS[seed % 4]
    pk = SCENE_PEAK
    if kind == "x_crossing":
        ang = np.radians(rng.uniform(40.0, 90.0))
        az = rng.uniform(0, 2 * np.pi)
        a = arc(rng.uniform(4.0, 7.0), rng.uniform(3.0, 6.0), 0.0, az)
        b = arc(rng.uniform(4.0, 7.0), rng.uniform(3.0, 6.0), 0.0, az + ang)
        ribs = [TruthRibbon(1, 1, a, peak=pk * rng.uniform(*PEAK_JITTER), kind=kind),
                TruthRibbon(2, 1, b, peak=pk * rng.uniform(*PEAK_JITTER), kind=kind)]
        info = {"pair": kind, "angle_deg": float(np.degrees(ang))}
    elif kind == "t_contact":
        az = rng.uniform(0, 2 * np.pi)
        a = arc(rng.uniform(4.0, 7.0), 6.0, 0.0, az)
        ang = np.radians(rng.uniform(60.0, 90.0))
        Lb = rng.uniform(2.5, 4.0)
        mid = a[len(a) // 2]
        # the chord of an arc with azimuth az runs along (y, x) = (-sin az, cos az); B leaves A at angle ang
        d = np.array([0.0, -np.sin(az + ang), np.cos(az + ang)])
        s = np.linspace(0, Lb, int(np.ceil(Lb / TRUTH_STEP_UM)) + 1)
        b = mid[None, :] + s[:, None] * d[None, :]
        ribs = [TruthRibbon(1, 1, a, peak=pk * rng.uniform(*PEAK_JITTER), kind=kind),
                TruthRibbon(2, 1, b, peak=pk * rng.uniform(*PEAK_JITTER), kind=kind)]
        info = {"pair": kind, "angle_deg": float(np.degrees(ang))}
    elif kind == "parallel_xy":
        L = rng.uniform(3.0, 6.0)
        az = rng.uniform(0, 2 * np.pi)
        a = arc(L, 8.0, 0.0, az)
        # perpendicular to the chord (-sin az, cos az) in the image plane
        off = np.array([0.0, np.cos(az), np.sin(az)]) * 0.5
        ribs = [TruthRibbon(1, 1, a - off / 2, peak=pk * rng.uniform(*PEAK_JITTER), kind=kind),
                TruthRibbon(2, 1, a + off / 2, peak=pk * rng.uniform(*PEAK_JITTER), kind=kind)]
        info = {"pair": kind, "separation_um": 0.5}
    elif kind == "stacked_z":
        L = rng.uniform(3.0, 6.0)
        a = arc(L, rng.uniform(3.0, 6.0), 0.0, rng.uniform(0, 2 * np.pi))
        off = np.array([1.0, 0.0, 0.0])
        ribs = [TruthRibbon(1, 1, a - off / 2, peak=pk * rng.uniform(*PEAK_JITTER), kind=kind),
                TruthRibbon(2, 1, a + off / 2, peak=pk * rng.uniform(*PEAK_JITTER), kind=kind)]
        info = {"pair": kind, "separation_um": 1.0}
    else:
        raise ValueError(f"unknown pair kind {kind!r}")
    return _scene_around("pairs", seed, rng, ribs, peak_info=info)


def dense_scene(seed: int, *, final: bool = False, germ_cell: str | None = None, peak: float | None = None,
                family: str = "dense") -> Scene:
    seed = check_seed(seed, final)
    rng = _rng(family, seed)
    gc = germ_cell or ("oocyte" if seed % 2 == 0 else "spermatocyte")
    pk = SCENE_PEAK if peak is None else float(peak)
    return _single_nucleus_scene(family, seed, rng, lambda c, R: _shell_nucleus(rng, gc, c, R, pk, 1, 1),
                                 germ_cell=gc, peak_info={"scene_peak": pk})


def snr_scene(seed: int, *, final: bool = False) -> Scene:
    seed = check_seed(seed, final)
    pk = SNR_PEAKS[seed % 4]
    gc = "oocyte" if (seed // 4) % 2 == 0 else "spermatocyte"
    return dense_scene(seed, final=final, germ_cell=gc, peak=pk, family="snr")


def noise_scene(seed: int, *, final: bool = False) -> Scene:
    seed = check_seed(seed, final)
    rng = _rng("noise", seed)
    return _single_nucleus_scene("noise", seed, rng, lambda c, R: [], germ_cell="oocyte")


def phantom_gonad(seed: int = 0, n_nuclei: int = 6, *, final: bool = False, spacing_um: float = 5.4,
                  peak: float = SCENE_PEAK, foci_per_nucleus: int = 0) -> Scene:
    """A row of dense nuclei along x (alternating oocyte-sized shells), one label each, for pipeline tests."""
    seed = check_seed(seed, final)
    rng = _rng("phantom_gonad", seed)
    centers, radii, label_r, gcs, ribs = [], [], [], [], []
    for i in range(n_nuclei):
        R = rng.uniform(2.2, 2.6)
        c = np.array([0.0, 0.0, i * spacing_um])
        rs = _shell_nucleus(rng, "oocyte", c, R, peak, i + 1, len(ribs) + 1)
        if foci_per_nucleus:
            for rb in rs[:foci_per_nucleus]:
                rb.foci_s = [float(rng.uniform(0.6, rb.length_um - 0.6))]
        centers.append(c); radii.append(R); label_r.append(R + rng.uniform(-0.2, 0.4)); gcs.append("oocyte")
        ribs.extend(rs)
    shape, shift = _place([r.points for r in ribs], centers, label_r)
    for r in ribs:
        r.points = r.points + shift
    centers = [c + shift for c in centers]
    return _finish("phantom_gonad", seed, rng, shape, centers, radii, label_r, gcs, ribs,
                   foci_channel=foci_per_nucleus > 0)


FAMILY_FUNCS = {"iso": iso_scene, "dense": dense_scene, "gaps": gaps_scene, "foci": foci_scene,
                "pairs": pairs_scene, "noise": noise_scene, "snr": snr_scene}


def make_scene(family: str, seed: int, *, final: bool = False, **kw) -> Scene:
    if family not in FAMILY_FUNCS:
        raise ValueError(f"unknown family {family!r}; one of {sorted(FAMILY_FUNCS)}")
    return FAMILY_FUNCS[family](seed, final=final, **kw)
