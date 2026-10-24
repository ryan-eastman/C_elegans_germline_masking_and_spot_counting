"""Phase A, one nucleus crop: transform, line fields, following, end matching, accounting, refinement.

``trace_crop(crop, ctx)`` runs in a worker process (or inline) and returns plain data: the strands with their
refined spline state, the samples needed for the across-nuclei accounting, the PSF stretches and the per-nucleus
diagnostics. Any exception becomes a ``failed`` record with its type and message (never a zero row).
"""
from __future__ import annotations

import time
import traceback

import numpy as np
from scipy.spatial import cKDTree

from .account import dedupe_within, drop_axial_ghosts, drop_shadowed, resample
from .crops import Crop, territory
from .fields import line_fields
from .follow import RawStrand, follow_nucleus, free_orphan_contacts, resolve_contacts
from .link import link_strands
from .noise import NoiseModel, gat, nucleus_noise_sd
from .params import Params
from .refine import axial_frames, core_profile, refine_imported, refine_points
from .sampling import LabelLookup, Sampler


def _scaled(P, ax, lat):
    return np.stack([P[:, 0] / ax, P[:, 1] / lat, P[:, 2] / lat], 1)


def background(raw: np.ndarray, terr: np.ndarray, pct: float) -> float:
    v = raw[terr]
    if v.size == 0:
        return float(np.median(raw))
    return float(np.percentile(v, pct))


def psf_stretches(lines: list, raw_sampler: Sampler, bg: float, noise: NoiseModel, spacing, p) -> list[dict]:
    """Flat, isolated, bright stretches (protocol of DESIGN 4.9) of one nucleus: their mean lateral (u) and axial
    (v) cross-sections of the background-subtracted counts."""
    if not lines:
        return []
    ds = 0.05
    tz_max = float(p["psf.flat_tz_max"])
    iso = float(p["psf.isolation_um"])
    min_snr = float(p["psf.min_snr"])
    min_len = float(p["psf.min_stretch_um"])
    samples = []
    for line in lines:
        s = np.arange(0.3, max(line.L - 0.3, 0.3), ds)
        samples.append((s, line.point(s) if len(s) else np.zeros((0, 3))))
    out = []
    U = np.round(np.arange(-12, 13) * 0.05, 6)
    V = np.round(np.arange(-20, 21) * 0.05, 6)
    for i, line in enumerate(lines):
        s, C = samples[i]
        if len(s) < int(min_len / ds):
            continue
        others = [samples[j][1] for j in range(len(lines)) if j != i and len(samples[j][1])]
        if others:
            d, _ = cKDTree(np.vstack(others)).query(C)
        else:
            d = np.full(len(C), np.inf)
        T = line.tangent(s)
        n_lat, n_ax = axial_frames(T)
        _core, _se, snr = core_profile(raw_sampler, C, n_lat, n_ax, bg, noise, spacing)
        ok = (np.abs(T[:, 0]) <= tz_max) & (d > iso) & (snr >= min_snr)
        # runs of ok samples at least min_len long
        idx = np.nonzero(ok)[0]
        if not len(idx):
            continue
        breaks = np.nonzero(np.diff(idx) > 1)[0]
        starts = np.r_[idx[0], idx[breaks + 1]]
        ends = np.r_[idx[breaks], idx[-1]]
        for a, b in zip(starts, ends):
            if (b - a) * ds < min_len - 1e-9:
                continue
            Cs, nl, na = C[a:b + 1], n_lat[a:b + 1], n_ax[a:b + 1]
            lat = raw_sampler(Cs[:, None, :] + U[None, :, None] * nl[:, None, :]) - bg
            ax = raw_sampler(Cs[:, None, :] + V[None, :, None] * na[:, None, :]) - bg
            out.append({"lat": np.nanmean(lat, 0).tolist(), "ax": np.nanmean(ax, 0).tolist(), "n": int(b - a + 1),
                        "u": U.tolist(), "v": V.tolist()})
    return out


def trace_crop(crop: Crop, ctx: dict) -> dict:
    t0 = time.perf_counter()
    nid = crop.nucleus_id
    rec = {"nucleus_id": nid, "status": "traced", "error": "", "error_message": "", "strands": [],
           "psf_stretches": [], "n_seeds": 0, "budget_exceeded": False, "noise_sd": float("nan"),
           "bg": float("nan"), "n_territory": 0, "touches_z": crop.touches_z, "touches_xy": crop.touches_xy,
           "n_contacts_freed": 0}
    try:
        p = Params(ctx["params"])
        sp = tuple(float(v) for v in ctx["spacing"])
        noise = NoiseModel(**ctx["noise"])
        terr = territory(crop.labels, nid, sp, float(p["territory_reach_um"]))
        rec["n_territory"] = int(terr.sum())
        if rec["n_territory"] == 0:
            rec["status"] = "skipped"
            return rec
        raw = crop.ce.astype(np.float32)
        bg = background(raw, terr, float(p["profile.background_pct"]))
        rec["bg"] = bg
        y = gat(raw, noise)
        rec["noise_sd"] = nucleus_noise_sd(y, terr)
        strands: list[RawStrand] = []
        geometry = ctx.get("geometry", "auto")
        Fd = None
        if geometry in ("auto", "both"):
            Fd = line_fields(y, terr, ctx["sigma_vox"], ctx["null_q"], max_off=float(p["fields.crest_offset_max"]))
            raw_strands, finfo = follow_nucleus(Fd, terr, sp, crop.org, p, z_planes=crop.image_shape[0],
                                                budget_s=float(p["budget_s"]), t_start=t0)
            rec["n_seeds"] = finfo["n_seeds"]
            rec["budget_exceeded"] = bool(finfo["budget_exceeded"])
            for st in raw_strands:
                st.nucleus_id = nid
            # axial shadows go before the T contacts and the end matching (so a shadow never trims the strand it
            # stands on and is never joined to it) and again after them
            raw_strands, n_ghost0 = drop_axial_ghosts(raw_strands, p)
            resolve_contacts(raw_strands, p)
            min_len = float(p["follow.min_length_um"])
            raw_strands = [s for s in raw_strands if s.length_um >= min_len]
            kept, n_ghost = drop_axial_ghosts(dedupe_within(link_strands(raw_strands, p), p), p)
            n_ghost += n_ghost0
            kept, n_shadow = drop_shadowed(kept, p)
            # a contact end whose partner was removed above is a free end again (the end rule applies)
            rec["n_contacts_freed"] = free_orphan_contacts(kept, p)
            rec["n_axial_ghosts"] = n_ghost
            rec["n_shadowed"] = n_shadow
            strands.extend(kept)
        for imp in ctx.get("imported", {}).get(nid, []):
            strands.append(RawStrand(points=np.asarray(imp["points"], float), end_types=["free", "free"],
                                     nucleus_id=nid, source="imported", verified=bool(imp.get("verified", False)),
                                     import_id=str(imp.get("id", ""))))
        sampler_bg = Sampler(raw - bg, crop.org, sp, cval=0.0)
        raw_sampler = Sampler(raw, crop.org, sp, cval=np.nan)
        lab = LabelLookup(crop.labels, crop.org, sp)
        lin_max = float(p["follow.lineness_max"])
        lin_sampler = Sampler(Fd.LIN, crop.org, sp, cval=1.0) if Fd is not None else None
        c_ax, c_lat = (float(v) for v in p["contamination_um"])
        dense = [resample(st.points, 0.05)[0] for st in strands]
        out_lines = []
        for i, st in enumerate(strands):
            if st.source == "imported" and not ctx.get("refine_imported", True):
                ref = refine_imported(st.points, p)
            else:
                others = [dense[j] for j in range(len(strands)) if j != i]
                tree = cKDTree(_scaled(np.vstack(others), c_ax, c_lat)) if others else None

                def frozen(P, tree=tree):
                    fz = np.zeros(len(P), bool)
                    if tree is not None:
                        d, _ = tree.query(_scaled(P, c_ax, c_lat), distance_upper_bound=1.0)
                        fz |= np.isfinite(d)
                    if lin_sampler is not None:
                        fz |= lin_sampler(P) > lin_max
                    return fz

                ref = refine_points(st.points, sampler_bg, p, frozen_fn=frozen)
            if ref is None:
                continue
            line = ref.line
            s = np.arange(0.0, line.L + 1e-9, 0.05)
            P = line.point(s)
            labs = lab(P)
            out_lines.append(line)
            rec["strands"].append({
                "line": line.to_state(), "end_types": list(st.end_types), "links": st.links, "members": st.members,
                "source": st.source, "verified": st.verified, "import_id": st.import_id,
                "fit_rms_um": ref.fit_rms_um, "n_frozen": ref.n_frozen,
                "strength_mean": float(np.mean(st.strength)) if st.strength is not None and len(st.strength) else
                float("nan"),
                "samples": P, "sample_labels": labs, "length_traced_um": float(line.L),
            })
        auto_lines = [ln for ln, d in zip(out_lines, rec["strands"]) if d["source"] == "auto"]
        rec["psf_stretches"] = psf_stretches(auto_lines, raw_sampler, bg, noise, sp, p)
        # the foci of this nucleus for the foci estimate of the image PSF (protocol changelog entry 12)
        pf = ctx.get("psf_foci")
        if pf is not None and "foci" in crop.roles and len(pf.get("pos", [])):
            from .psf import foci_fits_crop

            rec["psf_foci"] = foci_fits_crop(crop.roles["foci"], crop.org, sp, noise, pf, p)
    except Exception as e:  # noqa: BLE001 - one nucleus must never stop the image
        rec["status"] = "failed"
        rec["error"] = type(e).__name__
        rec["error_message"] = f"{e}\n{traceback.format_exc(limit=6)}"
        rec["strands"] = []
    rec["elapsed_s"] = time.perf_counter() - t0
    return rec
