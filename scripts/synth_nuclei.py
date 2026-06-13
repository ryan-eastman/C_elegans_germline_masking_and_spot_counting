#!/usr/bin/env python
"""v1 synthetic germline-nucleus generator + microscope degradation, for training a
germline-specific Cellpose model. EVERYTHING here is a PLACEHOLDER to be calibrated against
real data — morphology, brightness, PSF, noise, packing.

Model of a pachytene nucleus: an ellipsoid 'envelope' (the INSTANCE LABEL = the whole envelope,
enclosing every chromosome) whose DAPI signal lives as K bright chromosome tracks on the
PERIPHERY (great-circle-ish arcs on the shell) over a dim interior — i.e. the chromatin ring,
not a solid blob. Nuclei are packed to touch. Then we degrade: anisotropic PSF blur (z coarser),
Poisson shot + Gaussian read noise, background. Output: (degraded_image, instance_label_map).

    python scripts/synth_nuclei.py        # writes results_fullres/synth_vs_real.png + a .npz
"""
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

# ---- acquisition geometry (matches the real data) ----
DZ, DY, DX = 0.20, 0.108, 0.108          # voxel um (anisotropy 1.85)
CROP_UM = (14.0, 40.0, 40.0)             # z, y, x  (a ~40um window, 70 z-slices)

# ---- PLACEHOLDER morphology/optics (to be calibrated) ----
P = dict(
    n_strings=6,         # 6 bivalents/SCs (oocyte); 5 for spermatocyte
    nuc_radius_um=2.3,   # ball the chromosomes occupy (~4.6um diam -- the TRUE envelope, bigger
                         # than the incomplete stock-cpsam masks that read ~2.9um)
    radius_jitter=0.3,
    pack_spacing_um=5.2, # centroid spacing -> some touching, some dark gaps (like real)
    chrom_len_frac=1.7,      # chromosome arc length / nucleus radius (longer -> more fill)
    chrom_seat_frac=0.42,    # chromosomes seated at this radius (off-center) and evenly
                             # distributed -> distinct segments, but fuller (less sparse)
    chrom_bow_frac=0.85,     # parabolic bow -> curvier loops that wind and fill the nucleus
    chrom_extra_wind=0.35,   # secondary perpendicular wiggle -> reliable coverage (fewer sparse)
    string_sigma_um=0.12,  # thread tube radius (thinner -> sharper, crisper SC/chromosome)
    string_bright=1.0,
    syp_sigma_um=0.09,   # SC central-element tube radius (thinner/sharper than DAPI chromatin)
    syp_bright=1.0,
    interior_bright=0.02,  # faint nucleoplasm between strands (real is mostly dark)
    psf_xy_um=0.13,      # lateral resolution (sigma) -- sharper; real scope resolves the SC crisply
    psf_z_um=0.35,       # axial resolution (sigma) -- heavier, anisotropic
    peak_photons=360.0,  # signal scale -> brighter nuclei
    read_noise=2.5,
    background=3.0,
)


def _grid(crop_um):
    nz = int(round(crop_um[0] / DZ)); ny = int(round(crop_um[1] / DY)); nx = int(round(crop_um[2] / DX))
    return nz, ny, nx


def _centers(crop_um, spacing, rng):
    # jittered 3D grid so nuclei pack and touch
    zs = np.arange(spacing / 2, crop_um[0], spacing)
    ys = np.arange(spacing / 2, crop_um[1], spacing)
    xs = np.arange(spacing / 2, crop_um[2], spacing)
    g = np.array(np.meshgrid(zs, ys, xs, indexing="ij")).reshape(3, -1).T
    g = g + rng.normal(0, spacing * 0.18, g.shape)
    return g


def _fib_sphere(n, rng):
    """n roughly-even directions on the unit sphere (Fibonacci), randomly rotated -> seats that
    spread the chromosomes around the nucleus so they read as DISTINCT segments."""
    i = np.arange(n) + 0.5
    phi = np.arccos(np.clip(1 - 2 * i / n, -1, 1))
    theta = np.pi * (1 + 5 ** 0.5) * i
    pts = np.stack([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)], 1)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return pts @ q


def _chromosome(center, R, rng, p, seat_dir):
    """One bivalent: a bowed (curved) arc seated near the periphery in direction seat_dir, tangent
    to the shell -> a distinct bright loop. The 6 are spread by _fib_sphere -> separated, with a
    dark center & dark gaps (like real pachytene), not a filled blob."""
    from scipy.ndimage import gaussian_filter1d
    seat = center + seat_dir * (p["chrom_seat_frac"] * R) + rng.normal(0, 0.07 * R, 3)
    d = rng.normal(size=3); d -= d.dot(seat_dir) * seat_dir; d /= np.linalg.norm(d)  # along-arc (tangent)
    bdir = np.cross(seat_dir, d)                                                     # bow direction
    nc = 5
    s = np.linspace(-0.5, 0.5, nc)
    bow = p["chrom_bow_frac"] * R * (1 - (2 * s) ** 2)            # parabolic -> curved arc
    wind = p["chrom_extra_wind"] * R * np.sin(np.linspace(0, np.pi * rng.uniform(1.2, 2.2), nc))
    ctrl = seat[None] + s[:, None] * (p["chrom_len_frac"] * R) * d[None] \
        + bow[:, None] * bdir[None] + wind[:, None] * seat_dir[None] \
        + rng.normal(0, 0.04 * R, (nc, 1)) * bdir[None]
    tt = np.linspace(0, 1, 60); tc = np.linspace(0, 1, nc)
    path = np.stack([np.interp(tt, tc, ctrl[:, k]) for k in range(3)], 1)
    path = gaussian_filter1d(path, 1.5, axis=0, mode="nearest")
    rad = path - center; rr = np.linalg.norm(rad, axis=1, keepdims=True)
    over = (rr > 0.95 * R).ravel()
    path[over] = center + rad[over] / rr[over] * (0.95 * R)       # keep inside the ball
    return path


def build_chromosomes(crop_um, p, rng):
    """Place nuclei (centers, radii) and generate each nucleus's chromosome strands.
    Returns (centers, radii, chroms) where chroms = [dict(nid, pts), ...]."""
    centers = _centers(crop_um, p["pack_spacing_um"], rng)
    radii = np.clip(p["nuc_radius_um"] + rng.normal(0, p["radius_jitter"], len(centers)), 0.9, 3.5)
    chroms = []
    for nid, (c, r) in enumerate(zip(centers, radii), start=1):
        for sd in _fib_sphere(p["n_strings"], rng):
            chroms.append({"nid": nid, "pts": _chromosome(c, r, rng, p, sd)})
    return centers, radii, chroms


def _splat_tube(pts_list, sigma_um, shape):
    """Render point-paths as bright tubes (splat -> anisotropic-Gaussian blur), normalized 0..1."""
    nz, ny, nx = shape
    vox = np.array([DZ, DY, DX])
    spikes = np.zeros(shape, np.float32)
    for pts in pts_list:
        vi = np.round(pts / vox).astype(int)
        ok = (vi[:, 0] >= 0) & (vi[:, 0] < nz) & (vi[:, 1] >= 0) & (vi[:, 1] < ny) & (vi[:, 2] >= 0) & (vi[:, 2] < nx)
        vi = vi[ok]
        np.add.at(spikes, (vi[:, 0], vi[:, 1], vi[:, 2]), 1.0)
    out = ndi.gaussian_filter(spikes, np.array([sigma_um / DZ, sigma_um / DY, sigma_um / DX]))
    return out / out.max() if out.max() > 0 else out


def _nucleus_label(centers, radii, shape):
    """Instance label = nuclear envelope: nearest center within radius (one connected object per
    nucleus enclosing all its chromosomes; clean split where nuclei touch)."""
    nz, ny, nx = shape
    zz, yy, xx = np.mgrid[0:nz, 0:ny, 0:nx]
    vcoord = np.stack([zz * DZ, yy * DY, xx * DX], -1).reshape(-1, 3)
    dist, idx = cKDTree(centers).query(vcoord, workers=-1)
    return np.where(dist < radii[idx], idx + 1, 0).astype(np.int32).reshape(shape)


def synthesize(crop_um=CROP_UM, p=P, seed=0):
    """DAPI channel + nucleus instance labels (back-compat for nucleus-model training)."""
    rng = np.random.default_rng(seed)
    shape = _grid(crop_um)
    centers, radii, chroms = build_chromosomes(crop_um, p, rng)
    strands = _splat_tube([ch["pts"] for ch in chroms], p["string_sigma_um"], shape)
    label = _nucleus_label(centers, radii, shape)
    clean = p["string_bright"] * strands + p["interior_bright"] * (label > 0)
    return clean.astype(np.float32), label


def _arc_len(pts):
    return float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())


def fragment_strand(pts, rng, n_pieces, gap_um=0.4):
    """Break a chromosome's SC into n_pieces contiguous pieces separated by dark gaps (models
    heat-induced SC fragmentation). Returns a list of piece point-arrays."""
    if n_pieces <= 1:
        return [pts]
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg)]); total = cum[-1]
    breaks = np.sort(rng.uniform(0.12, 0.88, n_pieces - 1)) * total
    bounds = np.concatenate([[0.0], breaks, [total]])
    pieces = []
    for i in range(n_pieces):
        lo = bounds[i] + (gap_um / 2 if i > 0 else 0.0)
        hi = bounds[i + 1] - (gap_um / 2 if i < n_pieces - 1 else 0.0)
        m = (cum >= lo) & (cum <= hi)
        if m.sum() >= 2:
            pieces.append(pts[m])
    return pieces or [pts]


def synthesize_scene(crop_um=CROP_UM, p=P, seed=0, frag_rate=0.0, max_pieces=5):
    """Full scene: DAPI + SYP channels + nucleus labels + per-chromosome SC ground truth.

    frag_rate = fraction of chromosomes whose SC is fragmented (the heat phenotype); each such
    chromosome breaks into 2..max_pieces pieces with dark gaps. The SYP channel carries the SC
    (thinner/sharper than DAPI chromatin) along the (fragmented) strands.
    Returns (dapi, syp, label, sc_gt) where sc_gt is a per-chromosome DataFrame with the KNOWN
    sc_length_um, n_fragments and fragment lengths -> ground truth for the SC tracer.
    """
    import pandas as pd
    rng = np.random.default_rng(seed)
    shape = _grid(crop_um)
    centers, radii, chroms = build_chromosomes(crop_um, p, rng)
    label = _nucleus_label(centers, radii, shape)
    dapi = p["string_bright"] * _splat_tube([ch["pts"] for ch in chroms], p["string_sigma_um"], shape) \
        + p["interior_bright"] * (label > 0)

    syp_pts, gt = [], []
    for ch in chroms:
        k = int(rng.integers(2, max_pieces + 1)) if rng.random() < frag_rate else 1
        pieces = fragment_strand(ch["pts"], rng, k)
        syp_pts.extend(pieces)
        gt.append({"nucleus_id": ch["nid"], "sc_length_um": round(_arc_len(ch["pts"]), 2),
                   "n_fragments": len(pieces),
                   "frag_lengths_um": [round(_arc_len(pc), 2) for pc in pieces]})
    syp = p["syp_bright"] * _splat_tube(syp_pts, p["syp_sigma_um"], shape)
    return dapi.astype(np.float32), syp.astype(np.float32), label, pd.DataFrame(gt)


def degrade(clean, p=P, seed=1):
    rng = np.random.default_rng(seed)
    psf = np.array([p["psf_z_um"] / DZ, p["psf_xy_um"] / DY, p["psf_xy_um"] / DX])
    blur = ndi.gaussian_filter(clean, psf)
    sig = blur / max(blur.max(), 1e-6) * p["peak_photons"]
    noisy = rng.poisson(np.clip(sig, 0, None)).astype(np.float32)
    noisy += rng.normal(0, p["read_noise"], noisy.shape) + p["background"]
    return np.clip(noisy, 0, None)


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from skimage.segmentation import find_boundaries

    clean, label = synthesize()
    img = degrade(clean)
    z = img.shape[0] // 2
    mip = img[max(0, z - 3):z + 4].max(0)   # thin z-MIP (~1.4um): single slices cut nuclei at
    #                                         their edge & look sparse; the 3D nucleus is complete
    np.savez_compressed("results_fullres/synth_v1.npz", image=img, label=label)

    def norm(a):
        lo, hi = np.percentile(a, 1), np.percentile(a, 99.5)
        return np.clip((a - lo) / (hi - lo + 1e-9), 0, 1)

    # real crop for side-by-side (local full-res control gonad)
    real = None
    try:
        import glob, pandas as pd
        import germquant.io.nd2_reader as R
        st = R.read_stack(r"data/raw_examples/madeleline images/20251105_N2_noHS/20251105_N2_nohs_HERM _001.nd2")
        d = st.data[0]
        real = d[44:51, 1761 - 185:1761 + 185, 1423 - 185:1423 + 185].max(0)  # matching thin z-MIP
    except Exception as e:  # noqa: BLE001
        print("real crop unavailable:", e)

    n = 3 if real is not None else 2
    fig, ax = plt.subplots(1, n, figsize=(6 * n, 6)); k = 0
    if real is not None:
        ax[k].imshow(norm(real), cmap="gray"); ax[k].set_title("REAL DAPI (pachytene crop)"); ax[k].axis("off"); k += 1
    ax[k].imshow(norm(mip), cmap="gray"); ax[k].set_title("SYNTHETIC DAPI (degraded, z-MIP)"); ax[k].axis("off"); k += 1
    b = find_boundaries(label[z], mode="outer"); ov = np.zeros((*b.shape, 4)); ov[b] = (1, 1, 0, 1)
    ax[k].imshow(norm(mip), cmap="gray"); ax[k].imshow(ov); ax[k].set_title("SYNTHETIC + ground-truth labels (whole nucleus)"); ax[k].axis("off")
    fig.tight_layout(); fig.savefig("results_fullres/synth_vs_real.png", dpi=130, bbox_inches="tight")
    print("nuclei:", label.max(), "| wrote results_fullres/synth_vs_real.png + synth_v1.npz")
