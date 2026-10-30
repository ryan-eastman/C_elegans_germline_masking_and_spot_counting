"""Phase B, one nucleus crop: end rule, straightening, profiles, gaps, widths, foci and the nucleus diagnostics.

Every sample is a point of the measured centreline every ``ds`` (0.05 um) over s in [-0.5, L + 0.5] um with its
rotation-minimising frame; the cross-section at a sample is the grid u in [-0.6, 0.6] um (0.05 um, along
n_lat) by v in [-1.0, 1.0] um (0.1 um, along n_ax), sampled with linear interpolation on raw counts (NaN
outside the image). From it: ``core`` (mean of I - b in the ellipse 0.12 um lateral x 0.25 um axial), ``core_se``
(from the noise model), ``integrated`` (sum of (I - b) du dv in the ellipse 0.35 x 0.7 um: the continuity
signal), the flat strip (mean over |v| <= 0.3 um) and the axial strip (mean over |u| <= 0.15 um).
"""
from __future__ import annotations

import hashlib
import time
import traceback

import numpy as np
from scipy.spatial import cKDTree

from .crops import Crop, other_label_distance, territory
from .fields import line_fields
from .foci import assign_foci
from .frames import StrandFrame
from .noise import NoiseModel, gat
from .params import Params, switch_on
from .refine import (Centreline, EndFit, along_profile_fit, along_profile_sigma, ellipse_mask, finalize_line, fit_end,
                     gauss_fit, n_eff, one_gaussian, sigma_t)
from .sampling import LabelLookup, Sampler, grid, plane_points
from .trace import background

DU, DV = 0.05, 0.1


def strand_uid(C_in: np.ndarray) -> str:
    q = np.round(np.asarray(C_in, float) * 100.0).astype(np.int64)
    return hashlib.sha1(q.tobytes()).hexdigest()[:10]


def moving_mean(x: np.ndarray, w: int) -> np.ndarray:
    if w <= 1:
        return x.copy()
    k = np.ones(w)
    v = np.where(np.isfinite(x), x, 0.0)
    n = np.convolve(np.isfinite(x).astype(float), k, mode="same")
    return np.where(n > 0, np.convolve(v, k, mode="same") / np.maximum(n, 1e-9), np.nan)


def _crossing(s, y, i, level, step):
    """Walk from index i in direction ``step`` until y >= level; linear interpolation of the crossing."""
    j = i
    while 0 <= j + step < len(y) and not (np.isfinite(y[j + step]) and y[j + step] >= level):
        j += step
    k = j + step
    if not (0 <= k < len(y)) or not np.isfinite(y[k]) or not np.isfinite(y[j]):
        return float(s[j])
    y0, y1 = y[j], y[k]
    f = (level - y0) / (y1 - y0) if y1 != y0 else 0.5
    return float(s[j] + f * (s[k] - s[j]))


def _moving_median(x: np.ndarray, w: int) -> np.ndarray:
    if w <= 1 or len(x) < 3:
        return x.copy()
    from scipy.ndimage import median_filter

    return median_filter(x, size=w, mode="nearest")


def find_gaps(s, integ, integ_se, clean, in_strand, L, p, ds=0.05, sigma_along=None):
    """Gaps of one strand. The integrated profile is smoothed over ``gaps.smooth_um``; M is its median over clean,
    non-end samples. A candidate is a run below ``gaps.window_frac`` x M that does not touch an end. Its window is
    widened by ``gaps.tail_sigmas`` x the PSF sigma along the strand there (``sigma_along(s)``; the blur spreads
    the missing signal over that reach), and the dip is judged against the local plateau: a straight baseline
    between the median profiles within ``gaps.plateau_um`` outside the widened window on either side (one median
    when only one side has samples, M when neither has; the integrated signal of a ribbon changes along it, with
    its brightness and with the tilt, which changes how much of the PSF the integration ellipse holds). The
    candidate is a gap when the smoothed profile dips below ``gaps.q`` x the baseline, the dip is deeper than
    ``gaps.z`` standard errors and the length is at least ``gaps.min_um``. Because blurring along the strand
    conserves signal, the gap length is the missing fraction of the baseline integrated over the widened window,
    and the gap is centred on the deficit's centroid. Returns (gaps, M, smoothed)."""
    w = max(1, int(round(float(p["gaps.smooth_um"]) / ds)) + 1)
    sm = moving_mean(np.where(in_strand, integ, np.nan), w)
    nonend = (s >= 0.3) & (s <= L - 0.3)
    base = clean & nonend & in_strand & np.isfinite(sm)
    if base.sum() < 5:
        base = in_strand & np.isfinite(sm)
    M = float(np.median(sm[base])) if base.any() else float("nan")
    gaps = []
    if not (np.isfinite(M) and M > 0):
        return gaps, M, sm
    q = float(p["gaps.q"])
    wf = float(p["gaps.window_frac"])
    n_sig = float(p["gaps.tail_sigmas"])
    pl_um = float(p["gaps.plateau_um"])
    ok = in_strand & np.isfinite(sm)
    below = ok & (sm < wf * M)
    idx = np.nonzero(below)[0]
    if not len(idx):
        return gaps, M, sm
    br = np.nonzero(np.diff(idx) > 1)[0]
    windows = list(zip(np.r_[idx[0], idx[br + 1]], np.r_[idx[br], idx[-1]]))
    for wa, wb in windows:
        if s[wa] <= 0.1 or s[wb] >= L - 0.1 or not (ok[wa - 1] if wa > 0 else False) or \
                not (ok[wb + 1] if wb + 1 < len(sm) else False):
            continue
        kk = wa + int(np.nanargmin(sm[wa:wb + 1]))
        sig = float(sigma_along(float(s[kk]))) if sigma_along is not None else 0.0
        lo = max(float(s[wa]) - n_sig * sig, 0.1)
        hi = min(float(s[wb]) + n_sig * sig, L - 0.1)
        win = ok & (s >= lo - 1e-9) & (s <= hi + 1e-9)
        win[wa:wb + 1] = True
        left = ok & (s >= lo - pl_um) & (s < lo)
        right = ok & (s > hi) & (s <= hi + pl_um)
        wi = np.nonzero(win)[0]
        if left.sum() >= 3 and right.sum() >= 3:
            # the plateau may change along the strand: a straight baseline between the two flank medians
            pa, pb = float(np.median(sm[left])), float(np.median(sm[right]))
            base_w = pa + (pb - pa) * (s[wi] - lo) / max(hi - lo, 1e-9)
        elif (left | right).sum() >= 5:
            base_w = np.full(len(wi), float(np.median(sm[left | right])))
        else:
            base_w = np.full(len(wi), M)
        if not (np.all(np.isfinite(base_w)) and np.all(base_w > 0)):
            base_w = np.full(len(wi), M)
        P = float(np.mean(base_w))
        Pk = float(np.interp(s[kk], s[wi], base_w))
        dip = float(sm[kk])
        if not dip < q * Pk:
            continue
        seg = np.where(np.isfinite(integ[wi]), integ[wi], base_w)
        deficit = base_w - seg
        length = float((deficit / base_w).sum() * ds)
        if length <= 0:
            continue
        pos = np.clip(deficit / base_w, 0, None)
        centre = float((s[wi] * pos).sum() / max(pos.sum(), 1e-12))
        se = float(integ_se[kk]) if np.isfinite(integ_se[kk]) and integ_se[kk] > 0 else float("nan")
        z = (Pk - dip) / se if np.isfinite(se) else float("nan")
        if length < float(p["gaps.min_um"]) - 1e-9:
            continue
        if not (np.isfinite(z) and z > float(p["gaps.z"])):
            continue
        g0 = max(centre - 0.5 * length, 0.0)
        g1 = min(centre + 0.5 * length, L)
        gaps.append({"s_start_um": g0, "s_end_um": g1, "length_um": length, "depth": (Pk - dip) / Pk,
                     "z_score": z, "i_min": int(kk), "plateau": P, "plateau_rel": P / M})
    return gaps, M, sm


def section_signal(cs_ce, bg, ell_int, ring, du, dv, noise, sp, n_ax, smooth_n: int):
    """From cross-sections of raw counts (m, nu, nv): the integrated signal above the territory background in the
    ellipse (counts um^2), its standard error, and the same above the local background (the median of the outer
    ring, smoothed over ``smooth_n`` samples)."""
    valid = np.isfinite(cs_ce)
    wi = ell_int[None] & valid
    frac_valid = wi.sum((1, 2)) / max(ell_int.sum(), 1)
    integ = np.nansum(np.where(wi, cs_ce - bg, 0.0), (1, 2)) * du * dv
    integ[frac_valid < 0.7] = np.nan
    var_sum = np.nansum(np.where(wi, noise.var(np.nan_to_num(cs_ce, nan=bg)), 0.0), (1, 2))
    e_ax = np.abs(n_ax[:, 0]) * sp[0] + (1 - np.abs(n_ax[:, 0])) * sp[1]
    npv = np.maximum((sp[1] / du) * (e_ax / dv), 1.0)
    integ_se = du * dv * np.sqrt(npv * var_sum)
    with np.errstate(invalid="ignore"):
        b_loc = np.nanmedian(np.where(ring[None], cs_ce, np.nan).reshape(len(cs_ce), -1), axis=1)
    b_loc = _moving_median(np.where(np.isfinite(b_loc), b_loc, bg), smooth_n)
    integ_local = integ - (b_loc - bg) * float(ell_int.sum()) * du * dv
    return integ, integ_se, integ_local


def pieces_of(L: float, gaps: list) -> list[float]:
    edges = [0.0]
    for g in sorted(gaps, key=lambda g: g["s_start_um"]):
        edges += [max(g["s_start_um"], edges[-1]), max(g["s_end_um"], g["s_start_um"])]
    edges.append(L)
    return [max(edges[i + 1] - edges[i], 0.0) for i in range(0, len(edges), 2)]


def _role_bg(arr, terr, pct):
    v = arr[terr]
    return float(np.percentile(v, pct)) if v.size else float("nan")


def measure_crop(crop: Crop, ctx: dict, strands_in: list[dict], neighbours: np.ndarray | None,
                 spots: dict | None) -> dict:
    t0 = time.perf_counter()
    nid = crop.nucleus_id
    rec = {"nucleus_id": nid, "strands": [], "foci": [], "error": "", "error_message": "",
           "on_strand_frac": float("nan"), "coverage_proxy": float("nan"), "snr_core_median": float("nan")}
    try:
        p = Params(ctx["params"])
        sp = tuple(float(v) for v in ctx["spacing"])
        noise = NoiseModel(**ctx["noise"])
        psf = tuple(float(v) for v in ctx["psf"])
        pct = float(p["profile.background_pct"])
        terr = territory(crop.labels, nid, sp, float(p["territory_reach_um"]))
        raw = crop.ce.astype(np.float32)
        bg = background(raw, terr, pct)
        y = gat(raw, noise)
        F = line_fields(y, terr, ctx["sigma_vox"], ctx["null_q"], max_off=float(p["fields.crest_offset_max"]),
                        with_vectors=False)
        ce_raw = Sampler(raw, crop.org, sp, cval=np.nan)
        ce_bg = Sampler(raw - bg, crop.org, sp, cval=0.0)
        S_s = Sampler(F.S, crop.org, sp, cval=0.0)
        L_s = Sampler(F.LIN, crop.org, sp, cval=1.0)
        lab = LabelLookup(crop.labels, crop.org, sp)
        dother = Sampler(other_label_distance(crop.labels, nid, sp), crop.org, sp, cval=np.inf)
        prof_roles = [r for r in ctx.get("profile_roles", []) if r in crop.roles]
        strip_roles = [r for r in ctx.get("strip_roles", ["central_element"])
                       if r == "central_element" or r in crop.roles]
        role_s = {r: Sampler(crop.roles[r].astype(np.float32), crop.org, sp, cval=np.nan) for r in crop.roles}
        role_bg = {r: _role_bg(crop.roles[r].astype(np.float32), terr, pct) for r in crop.roles}
        ends_mode = ctx.get("ends_mode", "fit")
        z_hi = (crop.image_shape[0] - 1) * sp[0]
        edge_um = float(p["qc.stack_edge_um"])
        # the punctum class (protocol changelog entry 16)
        pun_on = switch_on(p["punctum.enabled"])
        pun_excl = pun_on and switch_on(p["punctum.exclude"])
        pun_xy = float(p["punctum.max_xy_um"])
        pun_ratio = float(p["punctum.sigma_ratio"])
        # v3 (protocol changelog entry 32): the reference axial sigma from the image's foci, and the shape test
        pun_ref = ctx.get("punctum_ref")
        pun_one = pun_on and switch_on(p.get("punctum.one_gaussian", False))
        pun_rms = float(p.get("punctum.max_rel_rms", 0.2))
        det_max = float(ctx.get("detector_max") or 65535.0)
        # 1. measured centrelines (end rule), oriented so s = 0 is the end with the smaller rounded (z, y, x)
        lines, fits_all, es, ets = [], [], [], []
        for st in strands_in:
            line = Centreline.from_state(st["line"])
            et = list(st["end_types"])
            if st["source"] == "imported" and ends_mode == "as_drawn":
                fits = (EndFit("as_drawn"), EndFit("as_drawn"))
            else:
                fits = (fit_end(line, True, et[0], ce_bg, psf, p), fit_end(line, False, et[1], ce_bg, psf, p))
            line2, e0, e1 = finalize_line(line, fits)
            a = tuple(np.round(line2.point(0.0)[0], 2))
            b = tuple(np.round(line2.point(line2.L)[0], 2))
            if b < a:
                line2 = line2.reversed()
                fits = (fits[1], fits[0])
                e0, e1 = e1, e0
                et = [et[1], et[0]]
            lines.append(line2)
            fits_all.append(fits)
            es.append((e0, e1))
            ets.append(et)
        ds = float(p["strips.ds_um"])
        ext = float(p["strips.extend_um"])
        frames = [StrandFrame(ln, ds=ds, extend=ext) for ln in lines]
        # 2. contamination: samples near another strand's centreline (this nucleus and its neighbours)
        c_ax, c_lat = (float(v) for v in p["contamination_um"])
        insides = [fr.C[fr.in_strand] for fr in frames]
        U = grid(float(p["strips.half_um"][1]), float(p["strips.du_um"]))
        V = grid(float(p["strips.half_um"][0]), float(p["strips.dv_um"]))
        du, dv = float(p["strips.du_um"]), float(p["strips.dv_um"])
        core_ax, core_lat = (float(v) for v in p["profile.core_um"])
        int_ax, int_lat = (float(v) for v in p["profile.integrate_um"])
        ell_core = ellipse_mask(U, V, core_lat, core_ax)
        ell_int = ellipse_mask(U, V, int_lat, int_ax)
        flat_v = np.abs(V) <= float(p["strips.flat_v_um"]) + 1e-9
        axial_u = np.abs(U) <= float(p["strips.axial_u_um"]) + 1e-9
        iu0 = int(np.argmin(np.abs(U)))
        iv0 = int(np.argmin(np.abs(V)))
        all_core_snr = []
        for i, (fr, line) in enumerate(zip(frames, lines)):
            st = strands_in[i]
            others = [insides[j] for j in range(len(frames)) if j != i and len(insides[j])]
            if neighbours is not None and len(neighbours):
                others.append(np.asarray(neighbours, float))
            ins = fr.in_strand
            if others:
                Q = np.vstack(others)
                tree = cKDTree(np.stack([Q[:, 0] / c_ax, Q[:, 1] / c_lat, Q[:, 2] / c_lat], 1))
                d, _ = tree.query(np.stack([fr.C[:, 0] / c_ax, fr.C[:, 1] / c_lat, fr.C[:, 2] / c_lat], 1),
                                  distance_upper_bound=1.0)
                contaminated = np.isfinite(d)
            else:
                contaminated = np.zeros(len(fr.s), bool)
            pts = plane_points(fr.C, fr.n_lat, fr.n_ax, U, V)
            cs_ce = ce_raw(pts)                                   # (m, nu, nv) raw counts
            valid = np.isfinite(cs_ce)
            wc = ell_core[None] & valid
            core_abs = np.where(wc.sum((1, 2)) > 0, np.nansum(np.where(wc, cs_ce, 0.0), (1, 2)) /
                                np.maximum(wc.sum((1, 2)), 1), np.nan)
            core = core_abs - bg
            core_se = noise.sd(np.nan_to_num(core_abs, nan=bg)) / np.sqrt(n_eff(core_lat, core_ax, sp, fr.n_ax[:, 0]))
            # the local background around the strand (an outer ring of the cross-section, |u| >= ring_um, |v| <= 0.3
            # um, median, smoothed along s): gaps are measured on the signal above it, so a nucleoplasm level
            # above the territory's percentile does not fill the gap
            ring = (np.abs(U)[:, None] >= float(p["gaps.ring_um"]) - 1e-9) & (np.abs(V)[None, :] <= 0.3 + 1e-9)
            smooth_n = max(1, int(round(0.5 / ds))) | 1
            integ, integ_se, integ_local = section_signal(cs_ce, bg, ell_int, ring, du, dv, noise, sp, fr.n_ax,
                                                          smooth_n)
            snr = core / np.maximum(core_se, 1e-9)
            role_core = {}
            for r in prof_roles:
                csr = role_s[r](pts)
                wr = ell_core[None] & np.isfinite(csr)
                role_core[r] = np.where(wr.sum((1, 2)) > 0, np.nansum(np.where(wr, csr, 0.0), (1, 2)) /
                                        np.maximum(wr.sum((1, 2)), 1), np.nan) - role_bg[r]
            flat, axial = [], []
            for r in strip_roles:
                csr = cs_ce if r == "central_element" else role_s[r](pts)
                with np.errstate(invalid="ignore"):
                    flat.append(np.nanmean(np.where(flat_v[None, None, :], csr, np.nan), axis=2).T)
                    axial.append(np.nanmean(np.where(axial_u[None, :, None], csr, np.nan), axis=1).T)
            labs = lab(fr.C)
            in_label = labs == nid
            near_other = dother(fr.C, order=0) <= 0.3
            strength = S_s(fr.C)
            lineness = L_s(fr.C)
            roll = fr.roll_from_axial_deg()
            L = line.L
            psf_t = psf

            def sigma_along(sv, line=line, psf_t=psf_t):
                return sigma_t(psf_t, float(line.tangent(sv)[0][0]))

            gaps, M, _sm = find_gaps(fr.s, integ_local, integ_se, ~contaminated, ins, L, p, ds,
                                     sigma_along=sigma_along)
            in_gap = np.zeros(len(fr.s), bool)
            for g in gaps:
                in_gap |= (fr.s >= g["s_start_um"]) & (fr.s <= g["s_end_um"])
            link_mids = [lk["mid_um"] for lk in st.get("links", []) if "mid_um" in lk]
            if link_mids:
                ls, _, _, _ = fr.to_local(np.asarray(link_mids, float))
            else:
                ls = np.zeros(0)
            for g in gaps:
                g["bridged"] = bool(np.any((ls >= g["s_start_um"] - 0.1) & (ls <= g["s_end_um"] + 0.1)))
                mid = line.point(0.5 * (g["s_start_um"] + g["s_end_um"]))[0]
                g["mid_z_um"], g["mid_y_um"], g["mid_x_um"] = (float(v) for v in mid)
            pieces = pieces_of(L, gaps)
            clean = ins & ~contaminated
            nonend = (fr.s >= 0.3) & (fr.s <= L - 0.3) & ins
            sel = clean & nonend & np.isfinite(integ)
            if sel.sum() < 3:
                sel = ins & np.isfinite(integ)
            if sel.any():
                Mi = float(np.mean(integ[sel]))
                cv2 = float(np.var(integ[sel]) / Mi ** 2) if Mi != 0 else float("nan")
                cv_corr = float(np.sqrt(max(cv2 - float(np.mean(integ_se[sel] ** 2)) / Mi ** 2, 0.0))) if Mi > 0 \
                    else float("nan")
            else:
                cv_corr = float("nan")
            wsel = clean & nonend
            if wsel.sum() < 3:
                wsel = ins
            prof_lat = np.nanmean(cs_ce[wsel][:, :, iv0], 0) - bg if wsel.any() else np.full(len(U), np.nan)
            prof_ax = np.nanmean(cs_ce[wsel][:, iu0, :], 0) - bg if wsel.any() else np.full(len(V), np.nan)
            s_lat, _, _ = gauss_fit(U, prof_lat)
            s_axw, _, _ = gauss_fit(V, prof_ax)
            fwhm_lat, fwhm_ax = 2.3548 * s_lat, 2.3548 * s_axw
            nz = float(np.mean(np.abs(fr.n_ax[ins, 0]))) if ins.any() else 1.0
            psf_lat = 2.3548 * psf[1]
            psf_ax = 2.3548 * float(np.sqrt(psf[1] ** 2 * (1 - nz ** 2) + psf[0] ** 2 * nz ** 2))
            T = fr.T
            Tin = T[ins]
            Cin = fr.C[ins]
            kv = np.linalg.norm(line.curvature_vector(fr.s[ins]), axis=1) if ins.any() else np.zeros(0)
            n_k = max(1, int(round(float(p["qc.kink_um"]) / ds)))
            if len(Tin) > n_k:
                cosk = np.clip((Tin[:-n_k] * Tin[n_k:]).sum(1), -1, 1)
                kink = float(np.degrees(np.arccos(cosk)).max())
            else:
                kink = 0.0
            fine = np.arange(0.0, L + 1e-9, 0.01)
            Pf = line.point(fine) if len(fine) > 1 else line.point(np.array([0.0, L]))
            dxy = np.linalg.norm(np.diff(Pf[:, 1:], axis=0), axis=1).sum()
            chord = float(np.linalg.norm(line.point(L)[0] - line.point(0.0)[0]))
            sat_vals = cs_ce[ins][:, ell_core] if ins.any() else np.zeros((0, 1))
            frac_sat = float(np.mean(sat_vals >= 0.999 * det_max)) if sat_vals.size else 0.0
            near_edge = bool(np.any((Cin[:, 0] < edge_um) | (Cin[:, 0] > z_hi - edge_um))) if len(Cin) else False
            f0, f1 = fits_all[i]
            stt = [sigma_t(psf, line.tangent(0.0)[0][0]), sigma_t(psf, line.tangent(L)[0][0])]
            se_e = [f.e_se_um if f.end_class == "free_fit" else stt[k] for k, f in enumerate((f0, f1))]
            se_len = float(np.sqrt(se_e[0] ** 2 + se_e[1] ** 2 + float(p["refine.se_curve_um"]) ** 2))
            costs = [lk["cost"] for lk in st.get("links", []) if "cost" in lk]
            snr_in = snr[ins & np.isfinite(snr)]
            all_core_snr.append(snr_in)
            e0p, e1p = line.point(0.0)[0], line.point(L)[0]
            row = {
                "nucleus_id": nid, "geometry_source": st["source"], "verified": bool(st.get("verified", False)),
                "import_id": st.get("import_id", ""),
                "length_um": L, "length_xy_um": float(dxy), "length_traced_um": float(st["length_traced_um"]),
                "e0_um": es[i][0], "e1_um": es[i][1], "length_se_um": se_len, "chord_um": chord,
                "straightness": chord / L if L > 0 else float("nan"),
                "curvature_mean_per_um": float(kv.mean()) if len(kv) else float("nan"),
                "curvature_max_per_um": float(kv.max()) if len(kv) else float("nan"),
                "tz_abs_mean": float(np.mean(np.abs(Tin[:, 0]))) if len(Tin) else float("nan"),
                "frac_steep45": float(np.mean(np.abs(Tin[:, 0]) > 0.707)) if len(Tin) else float("nan"),
                "frac_steep60": float(np.mean(np.abs(Tin[:, 0]) > 0.866)) if len(Tin) else float("nan"),
                "z_extent_um": float(Cin[:, 0].max() - Cin[:, 0].min()) if len(Cin) else 0.0,
                "fit_rms_um": float(st.get("fit_rms_um", float("nan"))),
                "n_links": len(costs), "link_cost_max": float(max(costs)) if costs else float("nan"),
                "n_members": int(st.get("members", 1)),
                "frac_in_label": float(np.mean(in_label[ins])) if ins.any() else float("nan"),
                "frac_near_other_label": float(np.mean(near_other[ins])) if ins.any() else float("nan"),
                "roll_from_axial_med_deg": float(np.nanmedian(np.abs(roll[ins]))) if np.isfinite(roll[ins]).any()
                else float("nan"),
                "end0_class": f0.end_class, "end1_class": f1.end_class,
                "end0_z_um": float(e0p[0]), "end0_y_um": float(e0p[1]), "end0_x_um": float(e0p[2]),
                "end1_z_um": float(e1p[0]), "end1_y_um": float(e1p[1]), "end1_x_um": float(e1p[2]),
                "e0_se_um": se_e[0], "e1_se_um": se_e[1], "end0_sigma_um": f0.sigma_um, "end1_sigma_um": f1.sigma_um,
                "core_median": float(np.nanmedian(core[sel])) if sel.any() else float("nan"),
                "integrated_median": float(np.nanmedian(integ[sel])) if sel.any() else float("nan"),
                "snr_core_median": float(np.median(snr_in)) if len(snr_in) else float("nan"),
                "profile_cv_corr": cv_corr, "n_gaps": len(gaps), "gap_um": float(sum(g["length_um"] for g in gaps)),
                "gap_fraction": float(sum(g["length_um"] for g in gaps) / L) if L > 0 else float("nan"),
                "gap_max_um": float(max((g["length_um"] for g in gaps), default=0.0)),
                "n_pieces": int(sum(1 for x in pieces if x >= 0.25)), "longest_piece_um": float(max(pieces)),
                "fwhm_lat_um": fwhm_lat, "fwhm_ax_um": fwhm_ax,
                "width_lat_decon_um": float(np.sqrt(max(fwhm_lat ** 2 - psf_lat ** 2, 0.0))) if np.isfinite(fwhm_lat)
                else float("nan"),
                "width_ax_decon_um": float(np.sqrt(max(fwhm_ax ** 2 - psf_ax ** 2, 0.0))) if np.isfinite(fwhm_ax)
                else float("nan"),
                "frac_contaminated": float(np.mean(contaminated[ins])) if ins.any() else float("nan"),
                "kink_max_deg": kink, "frac_saturated": frac_sat, "near_stack_edge": near_edge,
                "end_types_traced": "|".join(ets[i]), "error": "",
            }
            # punctum: a steep strand, short in the image plane, whose core profile along s is no wider in z than
            # the axial PSF allows a point (or a small aggregate) to be
            row["punctum_sigma_z_um"] = float("nan")
            row["punctum"] = False
            steep = np.isfinite(row["frac_steep60"]) and row["frac_steep60"] > 0.5
            if pun_on and st["source"] == "auto" and steep and row["length_xy_um"] <= pun_xy + 1e-12:
                if pun_one:
                    fit = along_profile_fit(fr.s, core)
                    sg = float(fit["sigma"])
                else:
                    fit, sg = None, along_profile_sigma(fr.s, core)
                if np.isfinite(sg) and np.isfinite(row["tz_abs_mean"]):
                    row["punctum_sigma_z_um"] = float(sg * row["tz_abs_mean"])
                    if pun_ref is not None and pun_ref.get("source") != "psf":
                        from .psf import punctum_sigma_at

                        ref_z = punctum_sigma_at(pun_ref, float(np.mean(Cin[:, 0])) if len(Cin) else
                                                 float(pun_ref.get("z_ref_um", 0.0)))
                    else:
                        ref_z = psf[0]
                    is_p = row["punctum_sigma_z_um"] <= pun_ratio * ref_z + 1e-12
                    if pun_one:
                        is_p = is_p and one_gaussian(fit, L, pun_rms)
                    row["punctum"] = bool(is_p)
            samples = {"s": fr.s, "in_strand": ins, "C": fr.C, "T": fr.T, "n_lat": fr.n_lat, "n_ax": fr.n_ax,
                       "core": core, "core_se": core_se, "integrated": integ, "integ_se": integ_se,
                       "contaminated": contaminated, "in_gap": in_gap, "in_label": in_label, "strength": strength,
                       "lineness": lineness, "roll": roll}
            for r in ("foci", "dna"):
                samples[f"{r}_core"] = role_core.get(r, np.full(len(fr.s), np.nan))
            rec["strands"].append({"row": row, "samples": samples, "gaps": gaps, "pieces": pieces,
                                   "strip_flat": np.stack(flat).astype(np.float32),
                                   "strip_axial": np.stack(axial).astype(np.float32),
                                   "strip_roles": list(strip_roles), "uid": strand_uid(fr.C[ins][::2]),
                                   "line": line.to_state()})
        # strands that count in the nucleus tier: under punctum.exclude the punctum strands are left out of the
        # foci assignment and of the strand tubes of on_strand_frac and coverage_proxy (their voxels count as
        # territory)
        keep = [k for k, s in enumerate(rec["strands"]) if not (pun_excl and s["row"].get("punctum"))]
        pun = [k for k, s in enumerate(rec["strands"]) if s["row"].get("punctum")]
        # 3. foci on the straightened strands of this nucleus
        if spots:
            if len(keep) == len(rec["strands"]):
                rec["foci"] = assign_foci(spots, frames, lines, [s["gaps"] for s in rec["strands"]], psf, p, nid)
            else:
                rows_f = assign_foci(spots, [frames[k] for k in keep], [lines[k] for k in keep],
                                     [rec["strands"][k]["gaps"] for k in keep], psf, p, nid)
                for r in rows_f:
                    for key in ("strand_key", "second_strand_key"):
                        if r[key] >= 0:
                            r[key] = keep[r[key]]
                rec["foci"] = rows_f
        # 4. nucleus diagnostics: on-strand SYP fraction, coverage proxy, SNR
        idx = np.argwhere(terr)
        ins_keep = [insides[k] for k in keep]
        rec["puncta_signal_frac"] = 0.0
        d = None
        if len(idx) and ins_keep and any(len(x) for x in ins_keep):
            Q = np.vstack([x for x in ins_keep if len(x)])
            tree = cKDTree(np.stack([Q[:, 0] / int_ax, Q[:, 1] / int_lat, Q[:, 2] / int_lat], 1))
            pts_um = (idx + np.asarray(crop.org)) * np.asarray(sp)
            d, _ = tree.query(np.stack([pts_um[:, 0] / int_ax, pts_um[:, 1] / int_lat, pts_um[:, 2] / int_lat], 1),
                              distance_upper_bound=1.0)
            vals = np.clip(raw[terr] - bg, 0, None)
            tot = float(vals.sum())
            rec["on_strand_frac"] = float(vals[np.isfinite(d)].sum() / tot) if tot > 0 else float("nan")
            crest = F.crest & terr & (F.S >= float(p["follow.seed_min"]))
            ci = np.argwhere(crest)
            if len(ci):
                cu = (ci + np.asarray(crop.org)) * np.asarray(sp)
                dd, _ = cKDTree(Q).query(cu, distance_upper_bound=0.3)
                rec["coverage_proxy"] = float(np.mean(np.isfinite(dd)))
        elif len(idx):
            rec["on_strand_frac"] = 0.0
        if pun and len(idx) and any(len(insides[k]) for k in pun):
            # the punctum tubes, outside the strand tubes counted above, so the two fractions partition the SYP
            Qp = np.vstack([insides[k] for k in pun if len(insides[k])])
            pts_um = (idx + np.asarray(crop.org)) * np.asarray(sp)
            scaled = np.stack([pts_um[:, 0] / int_ax, pts_um[:, 1] / int_lat, pts_um[:, 2] / int_lat], 1)
            dp, _ = cKDTree(np.stack([Qp[:, 0] / int_ax, Qp[:, 1] / int_lat, Qp[:, 2] / int_lat], 1)).query(
                scaled, distance_upper_bound=1.0)
            inp = np.isfinite(dp)
            if pun_excl and d is not None:
                inp &= ~np.isfinite(d)
            vals_p = np.clip(raw[terr] - bg, 0, None)
            tot_p = float(vals_p.sum())
            rec["puncta_signal_frac"] = float(vals_p[inp].sum() / tot_p) if tot_p > 0 else float("nan")
        snr_all = np.concatenate(all_core_snr) if all_core_snr else np.zeros(0)
        rec["snr_core_median"] = float(np.median(snr_all)) if len(snr_all) else float("nan")
    except Exception as e:  # noqa: BLE001 - one nucleus must never stop the image
        rec["error"] = type(e).__name__
        rec["error_message"] = f"{e}\n{traceback.format_exc(limit=6)}"
        rec["strands"] = []
        rec["foci"] = []
    rec["elapsed_s"] = time.perf_counter() - t0
    return rec
