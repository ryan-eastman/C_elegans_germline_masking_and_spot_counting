"""One image: noise model, null calibration, phase A (trace) per nucleus, accounting across nuclei, the PSF
estimate, phase B (measure) per nucleus, then the image-level steps (multiplicity, flags, classes, the nucleus
tier, stable strand ids, the strip layout).

``straighten_image`` works on arrays (whole-image central_element channel, labels, optional other roles) and
returns plain tables; ``germquant.sc.straighten.runner`` wraps it with the file inputs, refusals and outputs.
Work is spread over ``workers`` spawn-safe processes per nucleus; results are assembled in nucleus order and
every computation is deterministic, so the outputs are identical for every worker count.
"""
from __future__ import annotations

import logging
import multiprocessing as mp
import os
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from . import qc
from .account import RawStrand, dedupe_between
from .crops import boxes, make_crop, pad_vox
from .fields import calibrate_null, sigma_vox
from .measure import measure_crop
from .noise import fit_noise_model
from .params import Params, switch_on
from .psf import image_psf, isolated_spots
from .refine import Centreline
from .trace import trace_crop

log = logging.getLogger(__name__)


@dataclass
class ImageResult:
    strands: pd.DataFrame
    samples: pd.DataFrame
    gaps: pd.DataFrame
    foci: pd.DataFrame
    nuclei: pd.DataFrame
    strips_flat: np.ndarray            # (roles, nu, total) float32 (NaN outside the image)
    strips_axial: np.ndarray           # (roles, nv, total)
    strip_index: list                  # [[strand_id, offset, len], ...]
    strip_roles: list
    record: dict = field(default_factory=dict)
    flags: list = field(default_factory=list)
    psf_items: dict = field(default_factory=dict)   # the PSF stretches and foci fits (development scripts only)


# ------------------------------------------------------------------------------------ execution
def _init_worker():
    try:
        from threadpoolctl import threadpool_limits

        threadpool_limits(1)
    except Exception:  # noqa: BLE001 - optional
        pass


class _Runner:
    """Ordered map over tasks, inline (workers <= 1) or in a bounded window of spawn workers."""

    def __init__(self, workers: int):
        self.workers = max(1, int(workers))
        self.ex = None

    def __enter__(self):
        if self.workers > 1:
            self._env = {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
            for k in self._env:
                os.environ[k] = "1"
            self.ex = ProcessPoolExecutor(max_workers=self.workers, mp_context=mp.get_context("spawn"),
                                          initializer=_init_worker)
        return self

    def __exit__(self, *exc):
        if self.ex is not None:
            self.ex.shutdown(wait=True, cancel_futures=True)
            for k, v in self._env.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
        return False

    def map(self, fn, task_iter, n: int, progress=None, label=""):
        out = [None] * n
        if self.ex is None:
            import contextlib

            try:
                from threadpoolctl import threadpool_limits

                limit = threadpool_limits(1)
            except Exception:  # noqa: BLE001 - optional; single-threaded BLAS keeps workers=1 identical to N
                limit = contextlib.nullcontext()
            with limit:
                for k, args in enumerate(task_iter):
                    out[k] = fn(*args)
                    if progress:
                        progress(label, k + 1, n)
            return out
        window = 4 * self.workers
        pending = {}
        it = enumerate(task_iter)
        done = 0
        exhausted = False
        while True:
            while not exhausted and len(pending) < window:
                try:
                    k, args = next(it)
                except StopIteration:
                    exhausted = True
                    break
                pending[k] = self.ex.submit(fn, *args)
            if not pending:
                break
            k0 = min(pending)
            out[k0] = pending.pop(k0).result()
            done += 1
            if progress:
                progress(label, done, n)
        return out


# ------------------------------------------------------------------------------------ helpers
def _fit_noise(ce, labels, ids, p, min_samples: int = 200_000) -> "object":
    """The noise model on ``noise.planes`` evenly spaced planes of the selected nuclei's labels; when those hold
    fewer than ``min_samples`` voxels (a small image, a synthetic scene) every plane is used, and the number of
    bins shrinks with the sample count so every bin keeps a robust MAD."""
    Z = ce.shape[0]
    n = int(p["noise.planes"])
    idset = np.asarray(sorted(ids), dtype=np.int64)

    def collect(zs):
        planes, masks, tot = [], [], 0
        for z in zs:
            lab = labels[z]
            m = np.isin(lab, idset) if len(idset) else lab > 0
            if m.sum() < 50:
                continue
            planes.append(np.asarray(ce[z]))
            masks.append(m)
            tot += int(m.sum())
        return planes, masks, tot

    zs = np.unique(np.linspace(0, Z - 1, min(n, Z)).round().astype(int))
    planes, masks, tot = collect(zs)
    if tot < min_samples and len(zs) < Z:
        planes, masks, tot = collect(range(Z))
    bins = int(min(int(p["noise.bins"]), max(4, tot // 3000)))
    min_per_bin = int(max(100, min(400, tot // (3 * max(bins, 1)))))
    return fit_noise_model(planes, masks, bins=bins, qmax=float(p["noise.qmax"]),
                           var_floor=float(p["noise.var_floor"]), min_per_bin=min_per_bin)


def _strand_from_state(d: dict) -> RawStrand:
    return RawStrand(points=np.asarray(d["samples"], float), end_types=list(d["end_types"]), links=d["links"],
                     members=d["members"], source=d["source"], verified=d["verified"], import_id=d["import_id"],
                     strength=None)


def _own_frac(d: dict, nid: int) -> float:
    labs = np.asarray(d["sample_labels"])
    return float(np.mean(labs == nid)) if len(labs) else 0.0


def _restrict(d: dict, piece_pts: np.ndarray, e_types: list) -> dict:
    """A cut piece of a phase-A strand: the same spline restricted to the piece's s range."""
    line = Centreline.from_state(d["line"])
    S = np.asarray(d["samples"], float)
    tree = cKDTree(S)
    _, i0 = tree.query(piece_pts[0])
    _, i1 = tree.query(piece_pts[-1])
    a, b = sorted((int(i0), int(i1)))
    s_a, s_b = a * 0.05, min(b * 0.05, line.L)
    qa = float(line.q_of(s_a)[0])
    qb = float(line.q_of(s_b)[0])
    new = Centreline(line.spl, qa, qb)
    out = dict(d)
    out["line"] = new.to_state()
    out["samples"] = S[a:b + 1]
    out["sample_labels"] = np.asarray(d["sample_labels"])[a:b + 1]
    out["end_types"] = list(e_types) if i0 <= i1 else [e_types[1], e_types[0]]
    out["length_traced_um"] = new.L
    return out


def dedupe_across(recs: dict, nb_pairs: list, p) -> int:
    """Remove copies of a ribbon traced from two neighbouring nuclei; returns the number of strands changed."""
    changed = 0
    for a, b in nb_pairs:
        ra, rb = recs.get(a), recs.get(b)
        if not ra or not rb or not ra["strands"] or not rb["strands"]:
            continue
        ga = [_strand_from_state(d) for d in ra["strands"]]
        gb = [_strand_from_state(d) for d in rb["strands"]]
        fa = [_own_frac(d, a) for d in ra["strands"]]
        fb = [_own_frac(d, b) for d in rb["strands"]]
        # quick reject: no strand of a within reach of a strand of b
        Pa = np.vstack([g.points for g in ga])
        Pb = np.vstack([g.points for g in gb])
        if cKDTree(Pa).query(Pb, distance_upper_bound=0.5)[0].min() == np.inf:
            continue
        na, nb = dedupe_between(ga, gb, fa, fb, p)
        for rec, orig, new in ((ra, ga, na), (rb, gb, nb)):
            if len(new) == len(orig) and all(n is o for n, o in zip(new, orig)):
                continue
            out = []
            for piece in new:
                src = next((d for d, g in zip(rec["strands"], orig) if g is piece), None)
                if src is not None:
                    out.append(src)
                    continue
                # a cut piece: find its parent by the nearest original samples
                best, bd = None, np.inf
                for d in rec["strands"]:
                    dist = cKDTree(np.asarray(d["samples"])).query(piece.points[[0, -1]])[0].max()
                    if dist < bd:
                        best, bd = d, dist
                out.append(_restrict(best, piece.points, piece.end_types))
                changed += 1
            rec["strands"] = out
    return changed


def free_orphan_contacts_across(recs: dict, nb_pairs: list, p) -> int:
    """After the accounting across nuclei: a ``contact`` end of a phase-A strand with no surviving strand of its
    own nucleus or of a neighbouring one within ``follow.contact_free_um`` (axial, lateral) becomes ``free`` (the
    strand it stopped against was removed), so the end rule applies to it. Returns the number freed."""
    from .follow import orphan_contact_ends

    reach = tuple(float(v) for v in p["follow.contact_free_um"])
    nb: dict[int, set] = {}
    for a, b in nb_pairs:
        nb.setdefault(a, set()).add(b)
        nb.setdefault(b, set()).add(a)
    n = 0
    for nid in sorted(recs):
        rec = recs[nid]
        if not rec or rec.get("status") != "traced":
            continue
        strands = rec["strands"]
        for i, d in enumerate(strands):
            sides = [s for s in (0, 1) if d["end_types"][s] == "contact"]
            if not sides:
                continue
            rest = [np.asarray(o["samples"], float) for j, o in enumerate(strands) if j != i]
            for m in sorted(nb.get(nid, ())):
                rest.extend(np.asarray(o["samples"], float) for o in (recs.get(m) or {}).get("strands", []))
            Q = np.vstack(rest) if rest else np.zeros((0, 3))
            S = np.asarray(d["samples"], float)
            if not len(S):
                continue
            flags = orphan_contact_ends([(S[0] if s == 0 else S[-1], Q) for s in sides], reach)
            for s, orphan in zip(sides, flags):
                if orphan:
                    d["end_types"] = list(d["end_types"])
                    d["end_types"][s] = "free"
                    n += 1
    return n


def _neighbour_pairs(bx: dict, ids: list, spacing, pad_um) -> list[tuple[int, int]]:
    pv = pad_vox(spacing, pad_um)
    items = [(i, bx[i]) for i in ids if i in bx]
    lo = np.array([[s.start - p for s, p in zip(b, pv)] for _, b in items]) if items else np.zeros((0, 3))
    hi = np.array([[s.stop + p for s, p in zip(b, pv)] for _, b in items]) if items else np.zeros((0, 3))
    pairs = []
    for x in range(len(items)):
        ov = np.all((lo[x + 1:] < hi[x]) & (hi[x + 1:] > lo[x]), axis=1)
        for y in np.nonzero(ov)[0]:
            pairs.append((items[x][0], items[x + 1 + y][0]))
    return pairs


def _complete_candidate(r: dict, p) -> bool:
    lc = r.get("link_cost_max")
    return (r["end0_class"] == "free_fit" and r["end1_class"] == "free_fit"
            and (lc is None or not np.isfinite(lc) or lc <= float(p["classes.max_link_cost"]))
            and np.nan_to_num(r.get("frac_contaminated"), nan=1.0) <= float(p["classes.max_contaminated"])
            and np.nan_to_num(r.get("frac_in_label"), nan=0.0) >= float(p["classes.min_frac_in_label"])
            and r["length_um"] >= float(p["classes.min_complete_um"]))


def _longest_run(mask: np.ndarray, ds: float) -> float:
    best = cur = 0
    for v in mask:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    return best * ds


def classify(row: dict, flags: list, p) -> str:
    if row.get("verified"):
        return "verified"
    if row["length_um"] < float(p["classes.min_complete_um"]):
        return "fragment"
    if "merged_suspect" in flags:
        return "merged_suspect"
    lv = {qc.level(f, "strand") for f in flags}
    if _complete_candidate(row, p) and "review" not in lv and "fail" not in lv:
        return "complete"
    return "partial"


# ------------------------------------------------------------------------------------ the image
def straighten_image(ce: np.ndarray, labels: np.ndarray, spacing, nuclei: pd.DataFrame, *, params: Params,
                     roles: dict | None = None, spot_tables: dict | None = None, nucleus_ids=None,
                     imported: dict | None = None, geometry: str = "auto", refine_imported: bool = True,
                     ends_mode: str = "fit", workers: int = 1, detector_max: float | None = None,
                     image_id: str = "", progress=None, label_boxes: dict | None = None,
                     wavelengths: dict | None = None) -> ImageResult:
    t_start = time.perf_counter()
    p = params
    sp = tuple(float(v) for v in spacing)
    roles = {k: v for k, v in (roles or {}).items() if v is not None and k != "central_element"}
    nuc = nuclei.copy()
    nuc["nucleus_id"] = nuc["nucleus_id"].astype(int)
    if nucleus_ids is None:
        nucleus_ids = sorted(int(v) for v in nuc["nucleus_id"])
    ids = sorted(int(v) for v in nucleus_ids)
    bx = boxes(labels) if label_boxes is None else label_boxes
    noise = _fit_noise(ce, labels, [i for i in ids if i in bx], p)
    sigv = sigma_vox(p["fields.sigma_um"], sp)
    null_q = calibrate_null(sigv, float(p["fields.null_quantile"]), tuple(p["fields.null_shape"]),
                            int(p["fields.null_seed"]), float(p["fields.crest_offset_max"]))
    base_ctx = {"params": p.as_dict(), "spacing": sp, "noise": {"a": noise.a, "b": noise.b}, "sigma_vox": sigv,
                "null_q": null_q, "geometry": geometry, "refine_imported": bool(refine_imported),
                "ends_mode": ends_mode, "profile_roles": list(p["profile.roles"]),
                "strip_roles": list(p["strips.roles"]), "detector_max": detector_max}
    present = [i for i in ids if i in bx]
    imported = imported or {}
    pad = p["crop_pad_um"]
    # the foci estimate of the image PSF (protocol changelog entry 12): the isolated spots of the RAD-51 table,
    # fitted in phase A on each nucleus's foci-channel crop
    foci_by = psf_foci_inputs(spot_tables, roles, wavelengths, present, p)

    def tasks_a():
        for nid in present:
            fr = {"foci": roles["foci"]} if foci_by is not None else None
            crop = make_crop(nid, bx[nid], ce, labels, sp, pad, roles=fr)
            ctx = dict(base_ctx)
            ctx["imported"] = {nid: imported.get(nid, [])}
            if foci_by is not None:
                ctx["psf_foci"] = foci_by.get(nid, {"pos": np.zeros((0, 3)), "isolated": np.zeros(0, bool)})
            yield crop, ctx

    with _Runner(int(workers)) as runner:
        res_a = runner.map(trace_crop, tasks_a(), len(present), progress, "trace")
        recs = {r["nucleus_id"]: r for r in res_a}
        pairs = _neighbour_pairs(bx, present, sp, pad)
        n_cut = dedupe_across(recs, pairs, p)
        n_freed = {"within": int(sum(int(r.get("n_contacts_freed", 0)) for r in res_a)),
                   "across": free_orphan_contacts_across(recs, pairs, p)}
        stretches = [s for r in res_a for s in r.get("psf_stretches", [])]
        foci_fits = [f for r in res_a for f in r.get("psf_foci", [])] if foci_by is not None else None
        psf, psf_src, psf_detail = image_psf(stretches, foci_fits, wavelengths, p)
        # neighbours of each nucleus: in-strand samples of the strands of the nuclei it touches
        nbrs: dict[int, list] = {i: [] for i in present}
        for a, b in pairs:
            for x, y in ((a, b), (b, a)):
                for d in recs.get(y, {}).get("strands", []):
                    nbrs[x].append(np.asarray(d["samples"], float))
        spots_by = {}
        for name, df in (spot_tables or {}).items():
            if df is None or not len(df) or "nucleus_id" not in df:
                continue
            nid_col = pd.to_numeric(df["nucleus_id"], errors="coerce")
            for nid, g in df[nid_col.isin(present)].groupby(nid_col[nid_col.isin(present)].astype(int)):
                spots_by.setdefault(int(nid), {})[name] = g
        ctx_b = dict(base_ctx)
        ctx_b["psf"] = list(psf)

        def tasks_b():
            for nid in present:
                r = recs[nid]
                strands = r["strands"] if r["status"] == "traced" else []
                crop = make_crop(nid, bx[nid], ce, labels, sp, pad, roles=roles)
                nb = np.vstack(nbrs[nid]) if nbrs[nid] else None
                yield crop, ctx_b, strands, nb, spots_by.get(nid)

        res_b = runner.map(measure_crop, tasks_b(), len(present), progress, "measure")
    meas = {r["nucleus_id"]: r for r in res_b}
    out = _finish(image_id, p, sp, nuc, ids, present, recs, meas, noise, null_q, psf, psf_src, psf_detail,
                  n_cut, roles, spot_tables or {}, labels.shape, t_start, n_freed)
    out.psf_items = {"stretches": stretches, "foci": foci_fits}
    if foci_by is not None:
        out.record["psf"]["wavelengths_nm"] = dict(wavelengths or {})
    return out


def psf_foci_inputs(spot_tables, roles, wavelengths, present, p) -> dict | None:
    """{nucleus_id: {"pos": (n, 3) um, "isolated": (n,) bool}} of the RAD-51 ``spots`` table for the foci estimate
    of the PSF, or None when that estimate is off or impossible (no foci role, no spots, no wavelengths; then
    phase A is that of v1)."""
    if not switch_on(p["psf.foci"]) or p["psf.mode"] == "fixed":
        return None
    df = (spot_tables or {}).get("spots")
    wl = wavelengths or {}
    if df is None or not len(df) or "foci" not in (roles or {}) or not (wl.get("central_element") and wl.get("foci")):
        return None
    if not {"z_um", "y_um", "x_um", "nucleus_id"} <= set(df.columns):
        return None
    pos = df[["z_um", "y_um", "x_um"]].to_numpy(float)
    nid = pd.to_numeric(df["nucleus_id"], errors="coerce").to_numpy()
    ok = np.isfinite(pos).all(1)
    iso = np.zeros(len(pos), bool)
    iso[ok] = isolated_spots(pos[ok], float(p["psf.foci_isolation_um"]))
    out = {}
    pres = set(int(v) for v in present)
    for i in sorted(set(int(v) for v in nid[ok & np.isfinite(nid)])):
        if i not in pres:
            continue
        sel = ok & (nid == i)
        out[i] = {"pos": pos[sel], "isolated": iso[sel]}
    return out


def _finish(image_id, p, sp, nuc, ids, present, recs, meas, noise, null_q, psf, psf_src, psf_detail, n_cut,
            roles, spot_tables, shape, t_start, n_freed=None) -> ImageResult:
    ds = float(p["strips.ds_um"])
    all_strands = []                       # (nid, local index, measured strand dict)
    for nid in present:
        m = meas.get(nid)
        if not m:
            continue
        for k, s in enumerate(m["strands"]):
            all_strands.append((nid, k, s))
    # image-level normalisers: complete candidates' median integrated signal and widths
    cand = [s for _n, _k, s in all_strands if _complete_candidate(s["row"], p)]
    pool = cand if cand else [s for _n, _k, s in all_strands
                              if np.nan_to_num(s["row"]["snr_core_median"]) >= float(p["qc.low_snr_review"])]
    vals = [s["samples"]["integrated"][s["samples"]["in_strand"] & ~s["samples"]["contaminated"]] for s in pool]
    vals = np.concatenate(vals) if vals else np.zeros(0)
    vals = vals[np.isfinite(vals)]
    m_img = float(np.median(vals)) if len(vals) else float("nan")
    fl = [s["row"]["fwhm_lat_um"] for s in cand if np.isfinite(s["row"]["fwhm_lat_um"])]
    fa = [s["row"]["fwhm_ax_um"] for s in cand if np.isfinite(s["row"]["fwhm_ax_um"])]
    img_ctx = {"fwhm_lat_med": float(np.median(fl)) if fl else float("nan"),
               "fwhm_ax_med": float(np.median(fa)) if fa else float("nan")}
    mm = float(p["classes.merged_multiplicity"])
    # the reference of the multiplicity: the nucleus's own median integrated signal (clean in-strand samples of all
    # its strands, at least 3 um of them); SC brightness changes along the gonad, so one image-wide reference would
    # call every strand of the bright pachytene region merged. Nuclei with less SC use the image's complete-strand
    # median.
    ref_by_nid: dict[int, float] = {}
    by_n: dict[int, list] = {}
    for nid, _k, s in all_strands:
        smp = s["samples"]
        v = smp["integrated"][smp["in_strand"] & ~smp["contaminated"]]
        by_n.setdefault(nid, []).append(v[np.isfinite(v)])
    for nid, parts in by_n.items():
        v = np.concatenate(parts) if parts else np.zeros(0)
        ref_by_nid[nid] = float(np.median(v)) if len(v) * ds >= 3.0 and np.median(v) > 0 else m_img
    for nid, _k, s in all_strands:
        smp = s["samples"]
        ref = ref_by_nid.get(nid, m_img)
        mult = smp["integrated"] / ref if np.isfinite(ref) and ref > 0 else np.full(len(smp["s"]), np.nan)
        smp["multiplicity"] = mult
        ins = smp["in_strand"]
        r = s["row"]
        r["multiplicity_median"] = float(np.nanmedian(mult[ins])) if np.isfinite(mult[ins]).any() else float("nan")
        r["merged_run_um"] = _longest_run(np.nan_to_num(mult[ins], nan=0.0) >= mm, ds)
        r["frac_multiplicity_high"] = float(np.mean(np.nan_to_num(mult[ins], nan=0.0) >= mm)) if ins.any() else 0.0
        flags = qc.strand_flags(r, p, img_ctx)
        r["qc_flags"] = ";".join(flags)
        r["qc_status"] = qc.status(flags, "strand")
        r["strand_class"] = classify(r, flags, p)
    # stable strand ids
    order = sorted(range(len(all_strands)), key=lambda i: (
        all_strands[i][0], -round(all_strands[i][2]["row"]["length_um"], 6),
        tuple(round(all_strands[i][2]["row"][c], 4) for c in ("end0_z_um", "end0_y_um", "end0_x_um",
                                                               "end1_z_um", "end1_y_um", "end1_x_um"))))
    key_to_id = {}
    for sid, i in enumerate(order, start=1):
        nid, k, s = all_strands[i]
        key_to_id[(nid, k)] = sid
        s["row"]["strand_id"] = sid
        s["row"]["strand_uid"] = s["uid"]
    strip_roles = []
    for _n, _k, s in all_strands:
        if len(s["strip_roles"]) > len(strip_roles):
            strip_roles = list(s["strip_roles"])
    # strip layout (strand_id order)
    flat_parts, ax_parts, index = [], [], []
    off = 0
    for i in order:
        nid, k, s = all_strands[i]
        n = s["strip_flat"].shape[-1]
        fp = s["strip_flat"]
        ap = s["strip_axial"]
        if fp.shape[0] < len(strip_roles):
            pad_f = np.full((len(strip_roles) - fp.shape[0],) + fp.shape[1:], np.nan, np.float32)
            pad_a = np.full((len(strip_roles) - ap.shape[0],) + ap.shape[1:], np.nan, np.float32)
            fp, ap = np.concatenate([fp, pad_f]), np.concatenate([ap, pad_a])
        flat_parts.append(fp)
        ax_parts.append(ap)
        s["row"]["strip_offset"] = off
        s["row"]["strip_len"] = n
        index.append([s["row"]["strand_id"], off, n])
        off += n
    nu = len(np.arange(-round(float(p["strips.half_um"][1]) / float(p["strips.du_um"])),
                       round(float(p["strips.half_um"][1]) / float(p["strips.du_um"])) + 1))
    nv = len(np.arange(-round(float(p["strips.half_um"][0]) / float(p["strips.dv_um"])),
                       round(float(p["strips.half_um"][0]) / float(p["strips.dv_um"])) + 1))
    nr = max(len(strip_roles), 1)
    strips_flat = np.concatenate(flat_parts, axis=-1) if flat_parts else np.zeros((nr, nu, 0), np.float32)
    strips_axial = np.concatenate(ax_parts, axis=-1) if ax_parts else np.zeros((nr, nv, 0), np.float32)
    # foci rows
    foci_rows = []
    for nid in present:
        m = meas.get(nid)
        if not m:
            continue
        for r in m.get("foci", []):
            r = dict(r)
            k = r.pop("strand_key")
            k2 = r.pop("second_strand_key")
            r["strand_id"] = key_to_id.get((nid, k), pd.NA) if k >= 0 else pd.NA
            r["strand_uid"] = m["strands"][k]["uid"] if k >= 0 else ""
            r["second_strand_id"] = key_to_id.get((nid, k2), pd.NA) if k2 >= 0 else pd.NA
            foci_rows.append(r)
    n_foci = {}
    for r in foci_rows:
        if r["assignment"] in ("unique", "ambiguous") and not pd.isna(r["strand_id"]):
            n_foci[int(r["strand_id"])] = n_foci.get(int(r["strand_id"]), 0) + 1
    for _n, _k, s in all_strands:
        s["row"]["n_foci"] = n_foci.get(s["row"]["strand_id"], 0)
    # nucleus tier
    nuc_rows = _nucleus_rows(p, nuc, ids, present, recs, meas, all_strands, sp)
    # tables
    strands_df = pd.DataFrame([all_strands[i][2]["row"] for i in order])
    samples_df = _samples_table(image_id, [all_strands[i] for i in order], float(p["samples.step_um"]), ds)
    gap_rows = []
    for i in order:
        nid, k, s = all_strands[i]
        for g in s["gaps"]:
            gap_rows.append({"strand_id": s["row"]["strand_id"], "strand_uid": s["uid"], "nucleus_id": nid,
                             **{c: g[c] for c in ("s_start_um", "s_end_um", "length_um", "depth", "z_score",
                                                  "bridged", "mid_z_um", "mid_y_um", "mid_x_um")},
                             "plateau_rel": g.get("plateau_rel", np.nan)})
    for gi, g in enumerate(gap_rows, start=1):
        g["gap_id"] = gi
    gaps_df = pd.DataFrame(gap_rows)
    foci_df = pd.DataFrame(foci_rows)
    nuclei_df = pd.DataFrame(nuc_rows)
    # image flags and counts
    n_status = nuclei_df["status"].value_counts().to_dict() if len(nuclei_df) else {}
    present_rows = nuclei_df[nuclei_df["status"] == "traced"] if len(nuclei_df) else nuclei_df
    long_any = 0
    if len(strands_df) and len(present_rows):
        mx = strands_df.groupby("nucleus_id")["length_um"].max()
        long_any = int(sum(1 for nid in present_rows["nucleus_id"] if mx.get(nid, 0.0) >= 2.0))
    info = {"noise_rel_rms": noise.rel_rms, "noise_rel_rms_signal": noise.rel_rms_signal,
            "noise_fallback": noise.fallback, "psf_source": psf_src,
            "n_present": int(len(present_rows)), "n_present_with_2um": long_any,
            "spot_tables": sorted(k for k, v in spot_tables.items() if v is not None and len(v)),
            "n_label_merged": int(nuclei_df["qc_flags"].str.contains("label_merged_suspect").sum())
            if len(nuclei_df) else 0,
            "n_failed": int(n_status.get("failed", 0))}
    cls_counts = strands_df["strand_class"].value_counts().to_dict() if len(strands_df) else {}
    record = {
        "noise_model": noise.as_dict(), "null_q": null_q, "sigma_vox": list(sigma_vox(p["fields.sigma_um"], sp)),
        "psf": {"sigma_z_um": psf[0], "sigma_xy_um": psf[1], "source": psf_src, **psf_detail},
        "complete_median_integrated": m_img, "complete_median_fwhm": img_ctx, "n_cut_across_nuclei": n_cut,
        "n_contacts_freed": n_freed or {},
        "counts": {"nuclei": int(len(nuclei_df)), **{f"status_{k}": int(v) for k, v in n_status.items()},
                   "strands": int(len(strands_df)), **{f"class_{k}": int(v) for k, v in cls_counts.items()},
                   "foci": int(len(foci_df)), "gaps": int(len(gaps_df))},
        "nucleus_errors": {int(r["nucleus_id"]): r.get("error_message", "")
                           for r in nuc_rows if r["status"] == "failed"},
        "elapsed_s": round(time.perf_counter() - t_start, 2),
        "timing": _timing(recs, meas),
        "flag_info": info,
    }
    flags = qc.image_flags({**info, "geometry_source": None}, p)
    return ImageResult(strands=strands_df, samples=samples_df, gaps=gaps_df, foci=foci_df, nuclei=nuclei_df,
                       strips_flat=strips_flat, strips_axial=strips_axial, strip_index=index,
                       strip_roles=strip_roles, record=record, flags=flags)


def _timing(recs: dict, meas: dict) -> dict:
    """Per-nucleus work time (kept in the record, not in the tables, so the tables do not depend on the machine)."""
    tot = np.array([float(r.get("elapsed_s", 0.0)) + float(meas.get(n, {}).get("elapsed_s", 0.0))
                    for n, r in recs.items()])
    if not len(tot):
        return {}
    return {"trace_s": round(float(sum(float(r.get("elapsed_s", 0.0)) for r in recs.values())), 1),
            "measure_s": round(float(sum(float(m.get("elapsed_s", 0.0)) for m in meas.values())), 1),
            "per_nucleus_median_s": round(float(np.median(tot)), 3), "per_nucleus_max_s": round(float(tot.max()), 3),
            "slowest_nucleus": int(list(recs)[int(np.argmax(tot))])}


def _samples_table(image_id: str, items: list, step: float, ds: float) -> pd.DataFrame:
    every = max(1, int(round(step / ds)))
    frames = []
    for nid, _k, s in items:
        smp = s["samples"]
        n = len(smp["s"])
        # keep the samples on the step grid of s (s = -0.5 + i * ds); every ``every``-th one from s = 0
        i0 = int(round(0.5 / ds)) % every
        sel = np.arange(i0, n, every)
        C, T, nl, na = smp["C"][sel], smp["T"][sel], smp["n_lat"][sel], smp["n_ax"][sel]
        frames.append(pd.DataFrame({
            "image_id": image_id, "strand_id": s["row"]["strand_id"], "strand_uid": s["uid"], "nucleus_id": nid,
            "i": np.arange(len(sel)), "s_um": smp["s"][sel], "in_strand": smp["in_strand"][sel],
            "z_um": C[:, 0], "y_um": C[:, 1], "x_um": C[:, 2], "t_z": T[:, 0], "t_y": T[:, 1], "t_x": T[:, 2],
            "nlat_z": nl[:, 0], "nlat_y": nl[:, 1], "nlat_x": nl[:, 2], "nax_z": na[:, 0], "nax_y": na[:, 1],
            "nax_x": na[:, 2], "core": smp["core"][sel], "core_se": smp["core_se"][sel],
            "integrated": smp["integrated"][sel], "foci_core": smp["foci_core"][sel],
            "dna_core": smp["dna_core"][sel], "contaminated": smp["contaminated"][sel],
            "in_gap": smp["in_gap"][sel], "in_label": smp["in_label"][sel], "strength": smp["strength"][sel],
            "lineness": smp["lineness"][sel], "multiplicity": smp["multiplicity"][sel],
            "roll_from_axial_deg": smp["roll"][sel]}))
    if not frames:
        from .schema import SAMPLES

        return pd.DataFrame(columns=SAMPLES)
    return pd.concat(frames, ignore_index=True)


def _nucleus_rows(p, nuc, ids, present, recs, meas, all_strands, sp) -> list[dict]:
    by_nid: dict[int, list] = {}
    for nid, _k, s in all_strands:
        by_nid.setdefault(nid, []).append(s)
    info = nuc.set_index("nucleus_id")
    expected_map = dict(p["nucleus.expected"])
    min_main = float(p["classes.min_complete_um"])
    rows = []
    vols, axs = [], []
    for nid in ids:
        if nid in info.index and nid in recs and recs[nid]["status"] == "traced":
            v = info.loc[nid].get("volume_um3", np.nan) if hasattr(info.loc[nid], "get") else np.nan
            a = info.loc[nid].get("axis_position_um", np.nan) if hasattr(info.loc[nid], "get") else np.nan
            if np.isfinite(float(v)):
                vols.append(float(v))
                axs.append(float(a) if a is not None and not pd.isna(a) else np.nan)
    vols_a, axs_a = np.asarray(vols, float), np.asarray(axs, float)
    win = float(p["qc.merged_label_window_um"])

    def local_median(apos) -> float | None:
        return qc.local_volume_reference(vols_a, axs_a, apos, win)
    for nid in ids:
        r_a = recs.get(nid)
        r_b = meas.get(nid)
        meta = info.loc[nid] if nid in info.index else None

        def mget(c, default=np.nan, meta=meta):
            if meta is None or c not in meta.index:
                return default
            v = meta[c]
            return default if (not isinstance(v, str) and pd.isna(v)) else v

        st = by_nid.get(nid, [])
        rows_s = [s["row"] for s in st]
        err = ""
        if r_a is None:
            status = "skipped"
        elif r_a["status"] == "failed":
            status, err = "failed", r_a["error"]
        elif r_a["status"] == "skipped":
            status = "skipped"
        elif r_b is not None and r_b.get("error"):
            status, err = "failed", r_b["error"]
        else:
            status = "traced"
        total = float(sum(r["length_um"] for r in rows_s))
        gc = str(mget("germ_cell", ""))
        expected = expected_map.get(gc)
        main = [r for r in rows_s if r["length_um"] >= min_main]
        n_main = len(main)
        if p["nucleus.male_x_rule"] and gc == "spermatocyte" and n_main == 6:
            n_main = 5
        n_complete = sum(1 for r in rows_s if r.get("strand_class") in ("complete", "verified"))
        count_status = (pd.NA if expected is None else ("low" if n_main < expected else
                                                        ("high" if n_main > expected else "ok")))
        pieces = [x for s in st for x in s["pieces"]]
        pieces_ok = [x for x in pieces if x >= 0.25]
        n_gaps = int(sum(r["n_gaps"] for r in rows_s))
        gap_um = float(sum(r["gap_um"] for r in rows_s))
        snr = r_b.get("snr_core_median", np.nan) if r_b else np.nan
        if status == "traced":
            no_seed = r_a is not None and r_a.get("n_seeds", 0) == 0 and not any(
                r["geometry_source"] == "imported" for r in rows_s)
            if no_seed or total < float(p["nucleus.present_min_um"]) or \
                    (np.isfinite(snr) and snr < float(p["nucleus.present_min_snr"])):
                status = "no_sc"
        mult_len = 0.0
        for s in st:
            ins = s["samples"]["in_strand"]
            mult_len += float(np.sum(np.nan_to_num(s["samples"]["multiplicity"][ins], nan=0.0)
                                     >= float(p["classes.merged_multiplicity"]))) * float(p["strips.ds_um"])
        row = {
            "nucleus_id": nid, "status": status, "error": err, "n_strands": len(rows_s), "n_main": n_main,
            "n_complete": n_complete, "n_fragments": sum(1 for r in rows_s if r["length_um"] < min_main),
            "n_merged_suspect": sum(1 for r in rows_s if r.get("strand_class") == "merged_suspect"),
            "expected": expected if expected is not None else pd.NA, "count_status": count_status,
            "total_length_um": total, "total_length_xy_um": float(sum(r["length_xy_um"] for r in rows_s)),
            "n_pieces": len(pieces_ok), "mean_piece_um": float(np.mean(pieces_ok)) if pieces_ok else np.nan,
            "longest_piece_um": float(max(pieces)) if pieces else 0.0, "n_gaps": n_gaps,
            "gaps_per_10um": 10.0 * n_gaps / total if total >= 1.0 else np.nan,
            "gap_fraction": gap_um / total if total > 0 else np.nan,
            "piece_lw_mean_um": float(sum(x * x for x in pieces) / sum(pieces)) if pieces and sum(pieces) > 0 else 0.0,
            "on_strand_frac": r_b.get("on_strand_frac", np.nan) if r_b else np.nan,
            "snr_core_median": snr,
            "coverage_proxy": r_b.get("coverage_proxy", np.nan) if r_b else np.nan,
            "multiplicity_frac": mult_len / total if total > 0 else np.nan,
            "noise_sd": r_a.get("noise_sd", np.nan) if r_a else np.nan,
            "n_seeds": r_a.get("n_seeds", 0) if r_a else 0,
            "budget_exceeded": bool(r_a.get("budget_exceeded", False)) if r_a else False,
            "volume_um3": mget("volume_um3"), "in_germline": mget("in_germline", pd.NA),
            "axis_position_um": mget("axis_position_um"), "zone": mget("zone", pd.NA),
            "is_pachytene": mget("is_pachytene", pd.NA),
            "touches_z": bool(r_a.get("touches_z", False)) if r_a else False,
            "touches_xy": bool(r_a.get("touches_xy", False)) if r_a else False,
        }
        flags = qc.nucleus_flags(row, p, local_median(row["axis_position_um"]))
        complete_nucleus = bool(expected is not None and n_complete == expected and main and all(
            r.get("strand_class") in ("complete", "verified") for r in main) and qc.status(flags, "nucleus") != "fail")
        row["complete_nucleus"] = complete_nucleus
        row["total_length_complete_um"] = total if complete_nucleus else np.nan
        row["qc_flags"] = ";".join(flags)
        row["qc_status"] = qc.status(flags, "nucleus")
        if status == "failed":
            row["error_message"] = ((r_a or {}).get("error_message") or (r_b or {}).get("error_message") or "")
        rows.append(row)
    return rows
