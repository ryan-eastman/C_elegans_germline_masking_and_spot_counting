"""Pachytene staging from morphology, per the lab definition (Ryan, 2026-08-22):
  * TZ (transition-zone) nucleus = DNA in a polarized / crescent / clumped morphology
  * Pachytene START = first row, after the TZ is established, that contains <= TZ_MAX (=2) TZ nuclei
  * Pachytene END   = last row containing all pachytene nuclei with at most DIP_MAX (=1) post-pachytene nucleus
  * Pachytene region split into 3 equal-length zones: early / mid / late
Usage: python crescent_axis.py <image_id>

v3 features, all RELATIVE TO THE NUCLEUS ITSELF (no gonad-wide intensity thresholds), measured inside the
lamin envelope (watershed on LMN-1 seeded by the validated DAPI nuclei):
  gap     = fraction of envelope volume that is chromatin-dark at a per-nucleus threshold (solid mitotic: low)
  thick   = fraction of chromatin surviving a thin erosion (compact clumps high, pachytene threads low)
  syp_io  = mean SYP-3 inside envelope / mean SYP-3 in the 0.8 um perinuclear shell (mitotic <1, pachytene >1,
            drops again when the SC disassembles)
  vol_env = envelope volume (um^3)
Axis ordering without a hand trace: geodesic distance along a kNN graph of germline nucleus centroids
(follows the tube around the bend); distal end = end with the smaller nuclei, cross-checked by TZ position.
Rows = bins of one median nuclear diameter. Outputs in coloc_analysis/staging/<iid>_*:
nuclei.csv, rows.csv, landmarks.json, qc.png, strip.png (axis-ordered thumbnails for human verification).
Crops are cached to staging/cache/<iid>.npz so feature iterations do not re-read the nd2."""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from scipy import stats
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree, shortest_path
from scipy.spatial import cKDTree
from skimage.segmentation import watershed

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
import chload

SP = chload.SP
SPACING = tuple(SP)
VVOL = float(SP.prod())
PAD = 30
LAM_BAND_UM, LAM_SMOOTH_UM = 2.0, 0.15
TZ_MAX = 2        # pachytene row allows <= 2 TZ nuclei
TZ_ESTABLISH = 2  # TZ is "established" at the first row with >= 3 TZ nuclei (persisting)
DIP_MAX = 1       # "occasional single diplotene nucleus"
PERSIST = 2       # landmarks must hold for this many consecutive rows
KNN = 10
SHELL_UM = 0.8
RUN_DIRS = [r"C:/Users/ryane/ccw77_fullres", r"C:/Users/ryane/ccw77_0622", r"C:/Users/ryane/n2_dryice_results"]
E_DIRS = [r"E:/20260708_ccw77_IF_pgl1_syp3_LMN1", r"E:/20260622_ccw77_IF_syp3_pgl1_LMN1_dapi", r"C:/Users/ryane/n2_dryice_local"]
NAS_IP = "129.82.125.105"
OUTDIR = r"C:/Users/ryane/coloc_analysis/staging"
CACHE = os.path.join(OUTDIR, "cache")
os.makedirs(CACHE, exist_ok=True)


def find_run(iid):
    for rd in RUN_DIRS:
        p = glob.glob(os.path.join(rd, "**", iid + "__nuclei_labels.tif"), recursive=True)
        if p:
            return os.path.dirname(p[0])
    raise FileNotFoundError(iid)


def find_nd2(iid, csv_path):
    for d in E_DIRS:
        p = os.path.join(d, iid + ".nd2")
        if os.path.exists(p):
            return p
    fp = str(pd.read_csv(csv_path, nrows=1)["file_path"].iloc[0])
    fp = fp.replace(chr(92) * 2 + "bmb-ckc-nas", "//" + NAS_IP).replace(chr(92), "/")
    if os.path.exists(fp):
        return fp
    raise FileNotFoundError(f"nd2 unreachable for {iid} (no E:, no NAS)")


def lamin_labels(lab_crop, germ_ids, lamin):
    """per-nucleus envelope labels via seeded watershed on lamin; fallback to DAPI label if gate fails."""
    gn = np.isin(lab_crop, germ_ids)
    dt = ndi.distance_transform_edt(~gn, sampling=SPACING)
    elev = ndi.gaussian_filter(lamin.astype(np.float32), sigma=LAM_SMOOTH_UM / SP)
    markers = np.where(gn, lab_crop, 0).astype(np.int32)
    BGL = int(lab_crop.max()) + 1
    markers[dt > LAM_BAND_UM] = BGL
    ws = watershed(elev, markers=markers, mask=(dt <= LAM_BAND_UM) | gn | (dt > LAM_BAND_UM))
    lam = np.where(ws == BGL, 0, ws).astype(np.int32)
    dv = np.bincount(np.where(gn, lab_crop, 0).ravel(), minlength=BGL + 1)
    lv = np.bincount(lam.ravel(), minlength=BGL + 1)
    out = np.zeros_like(lam)
    fb = 0
    for nid in germ_ids:
        ok = dv[nid] > 0 and 1.0 <= lv[nid] / dv[nid] <= 3.0
        if ok:
            out[lam == nid] = nid
        else:
            out[(lab_crop == nid) & gn] = nid
            fb += 1
    return out, fb


def load_crops(iid):
    """DAPI / SYP / lamin-envelope labels / DAPI labels for the germline crop, cached as npz."""
    fz = os.path.join(CACHE, iid + ".npz")
    if os.path.exists(fz):
        z = np.load(fz)
        return z["dapi"], z["syp"], z["lam"], z["lab"], int(z["fb"])
    rd = find_run(iid)
    lab = tifffile.imread(os.path.join(rd, f"{iid}__nuclei_labels.tif"))
    csv = os.path.join(rd, f"{iid}__nuclei.csv")
    nc = pd.read_csv(csv)
    germ = [int(x) for x in nc[nc.in_germline.astype(bool)].nucleus_id]
    gn = np.isin(lab, germ)
    obj = ndi.find_objects(gn.astype(np.uint8))[0]
    sl = tuple(slice(max(0, o.start - PAD), min(dim, o.stop + PAD)) for o, dim in zip(obj, gn.shape))
    ch = chload.load(find_nd2(iid, csv))
    dapi = ch["dapi"][sl].astype(np.uint16)
    syp = ch["syp"][sl].astype(np.uint16)
    lamin = ch["lamin"][sl]
    del ch
    lab_c = lab[sl].astype(np.int32)
    lam, fb = lamin_labels(lab_c, germ, lamin)
    np.savez(fz, dapi=dapi, syp=syp, lam=lam, lab=lab_c, fb=fb)
    return dapi, syp, lam, lab_c, fb


def nucleus_features(lam, dapi, syp):
    dapi_bs = np.clip(dapi.astype(np.float32) - np.percentile(dapi, 5), 0, None)
    syp_bs = np.clip(syp.astype(np.float32) - np.percentile(syp, 5), 0, None)
    env_all = lam > 0
    objs = ndi.find_objects(lam)
    pad = int(round(SHELL_UM / SP[1])) + 1
    rows = []
    for nid, sl in enumerate(objs, start=1):
        if sl is None:
            continue
        m = lam[sl] == nid
        n = int(m.sum())
        if n < 50:
            continue
        zz, yy, xx = np.nonzero(m)
        off = np.array([s.start for s in sl])
        coords = (np.stack([zz, yy, xx], 1) + off) * SP
        d = dapi_bs[sl]
        inside = d[m]
        t_n = 0.5 * np.percentile(inside, 97)          # per-nucleus chromatin threshold
        C = (d > t_n) & m
        nC = int(C.sum())
        if nC < 20:
            continue
        gap = 1.0 - nC / n
        er = ndi.binary_erosion(C, structure=np.ones((1, 3, 3), bool))   # thin in-plane erosion
        thick = float(er.sum() / nC)
        # internal gap: dark space inside the chromatin's own footprint (solid ball ~0, clumps/threads high)
        closed = ndi.binary_closing(C, structure=np.ones((3, 7, 7), bool)) & m
        gap_int = float(1.0 - nC / max(int(closed.sum()), nC))
        # VOID analysis: in a TZ crescent the chromatin-free space is ONE contiguous pocket pushed to one
        # side; in pachytene the empty space is a network interleaved with threads.
        void = m & ~ndi.binary_dilation(C, structure=np.ones((1, 3, 3), bool))
        nv = int(void.sum())
        if nv > 30:
            vl, _ = ndi.label(void)
            cnt = np.bincount(vl.ravel())[1:]
            void_big = float(cnt.max() / nv)                       # 1.0 = single pocket (TZ)
            big = vl == (int(np.argmax(cnt)) + 1)
            vz, vy, vx = np.nonzero(big)
            c_void = (np.stack([vz, vy, vx], 1) + off).mean(0) * SP
        else:
            void_big, c_void = 0.0, None
        c_geo = coords.mean(0)
        c_chr = coords[C[m]].mean(0)
        vol = n * VVOL
        r_eq = (3 * vol / (4 * np.pi)) ** (1 / 3)
        # how far the main void sits from the nucleus centre (TZ: void on one side -> large)
        void_off = float(np.linalg.norm(c_void - c_geo) / r_eq) if c_void is not None else 0.0
        sl2 = tuple(slice(max(0, s.start - pad), s.stop + pad) for s in sl)
        m2 = lam[sl2] == nid
        sh = ndi.binary_dilation(m2, iterations=pad - 1) & ~env_all[sl2]
        s_in = float(syp_bs[sl][m].mean())
        s_out = float(syp_bs[sl2][sh].mean()) if sh.sum() > 50 else np.nan
        rows.append({"nucleus_id": nid, "cz": c_geo[0], "cy": c_geo[1], "cx": c_geo[2],
                     "vol_env": vol, "r_eq": r_eq, "gap": gap, "gap_int": gap_int, "thick": thick,
                     "void_big": void_big, "void_off": void_off,
                     "offset": float(np.linalg.norm(c_chr - c_geo) / r_eq),
                     "syp_io": s_in / s_out if s_out and s_out > 0 else np.nan})
    return pd.DataFrame(rows)


def largest_component(df, link_um=8.0, min_frac=0.10):
    """BUGFIX 3: the old code joined disconnected nucleus clouds with a Euclidean MST, which threaded
    carcass islands 40-100 um away onto the axis. Now: keep only the main connected component and
    REPORT what was discarded (that discarded material is exactly the contamination the PC audit found).
    """
    P = df[["cz", "cy", "cx"]].to_numpy()
    pairs = cKDTree(P).query_pairs(link_um, output_type="ndarray")
    if len(pairs) == 0:
        return df, {"n_components": len(P), "kept_frac": 0.0, "discarded": len(P)}
    G = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(P), len(P)))
    ncomp, lab = connected_components(G, directed=False)
    sizes = np.bincount(lab)
    keep = int(np.argmax(sizes))
    kept = df[lab == keep].reset_index(drop=True)
    big = [int(x) for x in sorted(sizes[sizes >= max(10, min_frac * sizes.max())], reverse=True)]
    return kept, {"n_components": int(ncomp), "component_sizes": big,
                  "kept_frac": round(float(sizes[keep] / len(P)), 3),
                  "discarded": int(len(P) - sizes[keep])}


def geodesic_axis(df):
    """Order nuclei along the gonad; return (s, L, polarity_info). Operates on ONE connected cloud."""
    P = df[["cz", "cy", "cx"]].to_numpy()
    n = len(P)
    d, j = cKDTree(P).query(P, k=min(KNN + 1, n))
    src = np.repeat(np.arange(n), d.shape[1] - 1)
    G = coo_matrix((d[:, 1:].ravel(), (src, j[:, 1:].ravel())), shape=(n, n)).tocsr()
    G = G.maximum(G.T)
    D0 = shortest_path(G, directed=False, indices=[0])[0]
    D0[~np.isfinite(D0)] = -1
    A = int(np.argmax(D0))
    DA = shortest_path(G, directed=False, indices=[A])[0]
    DA[~np.isfinite(DA)] = np.nanmax(DA[np.isfinite(DA)])
    B = int(np.argmax(DA))
    L = float(DA[B])
    return DA, L


def orient(df, s, L):
    """BUGFIX 2: the old rule called the end with SMALLER median nuclei 'distal'. Over-segmented
    diakinesis bivalents and spermatid dots at the PROXIMAL end are tiny, so that inverted the axis in
    4 gonads. Robust replacement: germ nuclei grow monotonically distal->proximal, so fit a Theil-Sen
    slope of binned median volume across the MIDDLE of the axis (ends excluded) using only non-fragment
    objects. Negative slope => flip. Reports the evidence so a wrong call is visible.
    """
    v = df.vol_env.to_numpy()
    frag = v < 0.3 * np.median(v)          # over-segmentation fragments / spermatid dots
    ok = ~frag
    nb = 20
    edges = np.linspace(0, L, nb + 1)
    ctr, med = [], []
    for i in range(nb):
        m = ok & (s >= edges[i]) & (s < edges[i + 1])
        if m.sum() >= 5:
            ctr.append(0.5 * (edges[i] + edges[i + 1])); med.append(float(np.median(v[m])))
    ctr, med = np.array(ctr), np.array(med)
    inner = (ctr > 0.1 * L) & (ctr < 0.9 * L)      # drop terminal bins (diakinesis / grazed tip)
    slope = np.nan
    if inner.sum() >= 4:
        slope = float(stats.theilslopes(med[inner], ctr[inner])[0])
    flip = bool(np.isfinite(slope) and slope < 0)
    return flip, {"vol_slope_um3_per_um": None if not np.isfinite(slope) else round(slope, 4),
                  "n_bins_used": int(inner.sum()), "frag_frac": round(float(frag.mean()), 3)}


def classify(df):
    """Per-nucleus stage call from chromatin texture, relative to this gonad's pachytene core.

    score = gap_int / thick rises monotonically: solid mitotic ball ~0.05 -> clumped/crescent TZ ~0.25
    -> pachytene threads >0.4. BUGFIX 1: `clumped` is now defined WITHOUT a position gate, so the
    polarity cross-check (are the crescents in the distal half?) can actually fire. Previously crescents
    were only ever defined for s_norm<0.5, making the guard dead code in all 22 gonads.
    """
    df["score"] = df.gap_int / df.thick.clip(lower=1e-6)
    core = df[(df.s_norm > 0.35) & (df.s_norm < 0.65)]
    p_ref = float(core.score.median())
    s_ref = float(core.syp_io.median())
    v_ref = float(core.vol_env.median())
    df["merged"] = df.vol_env > 2.2 * v_ref
    # position-FREE texture bands
    df["is_solid"] = (~df.merged) & (df.score < 0.25 * p_ref)
    df["is_clumped"] = (~df.merged) & (df.score >= 0.25 * p_ref) & (df.score < 0.72 * p_ref)
    # stage calls (position used only to separate pre- from post-pachytene, which look alike)
    distal = df.s_norm < 0.5
    df["mitotic"] = df.is_solid & (df.s_norm < 0.25)
    df["crescent"] = df.is_clumped & distal
    df["post"] = (~df.merged) & (~distal) & (df.score < 0.72 * p_ref) & (df.syp_io < 0.75 * s_ref)
    df["cls"] = np.select([df.crescent, df.mitotic, df.post], ["TZ", "mitotic", "post"], "pachytene")
    return df, {"score_ref": p_ref, "syp_io_ref": s_ref, "vol_ref": v_ref}


def landmarks(rt):
    pach_start = pach_end = tz_start = None
    for i in range(len(rt) - PERSIST + 1):
        if all(rt.n_tz.iloc[i + k] >= TZ_ESTABLISH for k in range(PERSIST)):
            tz_start = i
            break
    if tz_start is None:
        hits = rt.index[rt.n_tz >= 2].tolist()
        tz_start = hits[0] if hits else None
    if tz_start is not None:
        for i in range(tz_start, len(rt) - PERSIST + 1):
            if all(rt.n_tz.iloc[i + k] <= TZ_MAX for k in range(PERSIST)):
                pach_start = int(rt.row.iloc[i])
                break
    if pach_start is not None:
        after = rt[rt.row >= pach_start].reset_index(drop=True)
        end_row = None
        for i in range(len(after) - PERSIST + 1):
            if all(after.n_post.iloc[i + k] > DIP_MAX for k in range(PERSIST)):
                end_row = int(after.row.iloc[i]) - 1
                break
        pach_end = end_row if end_row is not None else int(after.row.iloc[-1])
    return (None if tz_start is None else int(rt.row.iloc[tz_start])), pach_start, pach_end


def thumb(dapi, syp, r, R):
    zc, yc, xc = int(r.cz / SP[0]), int(r.cy / SP[1]), int(r.cx / SP[2])
    dd = dapi[zc, max(0, yc - R):yc + R, max(0, xc - R):xc + R].astype(np.float32)
    ss = syp[zc, max(0, yc - R):yc + R, max(0, xc - R):xc + R].astype(np.float32)
    if dd.shape != (2 * R, 2 * R):
        return None
    dn = np.clip((dd - np.percentile(dd, 5)) / (np.percentile(dd, 99.5) - np.percentile(dd, 5) + 1e-9), 0, 1)
    sn = np.clip((ss - np.percentile(ss, 30)) / (np.percentile(ss, 99.5) - np.percentile(ss, 30) + 1e-9), 0, 1)
    return np.stack([np.maximum(dn, sn * 0.8), dn, dn], -1)


CMAP = {"TZ": (1, 0.55, 0.1), "pachytene": (0.25, 0.5, 0.8), "post": (0.6, 0.3, 0.7), "mitotic": (0.5, 0.5, 0.5)}


def main(iid):
    dapi, syp, lam, lab_c, fb = load_crops(iid)
    df = nucleus_features(lam, dapi, syp)
    n_all = len(df)
    df, comp = largest_component(df)          # BUGFIX 3: drop disconnected islands, do not bridge
    s, L = geodesic_axis(df)
    df["s_um"] = s
    df = df.sort_values("s_um").reset_index(drop=True)
    flip, ori = orient(df, df.s_um.to_numpy(), L)   # BUGFIX 2: volume-trend polarity
    if flip:
        df["s_um"] = L - df.s_um
        df = df.sort_values("s_um").reset_index(drop=True)
    df["s_norm"] = df.s_um / L
    df, refs = classify(df)
    # BUGFIX 1 cross-check: crescents must sit distally; this can now actually fire
    cres_pos = float(df.s_norm[df.is_clumped].median()) if df.is_clumped.any() else float("nan")
    contradiction = bool(df.is_clumped.sum() >= 8 and cres_pos > 0.65)
    if contradiction:
        df["s_um"] = L - df.s_um
        df = df.sort_values("s_um").reset_index(drop=True)
        df["s_norm"] = df.s_um / L
        df, refs = classify(df)
        flip = not flip
    row_w = float(2 * df.r_eq.median())
    df["row"] = (df.s_um // row_w).astype(int)
    rt = df.groupby("row").agg(n=("nucleus_id", "size"), n_tz=("crescent", "sum"), n_mit=("mitotic", "sum"),
                               n_post=("post", "sum"), s0=("s_um", "min")).reset_index()
    rt["n_pach"] = rt.n - rt.n_tz - rt.n_post - rt.n_mit
    tz_start, pach_start, pach_end = landmarks(rt)
    zone = np.full(len(df), "none", dtype=object)
    span = None
    if pach_start is not None and pach_end is not None and pach_end > pach_start:
        s0, s1 = pach_start * row_w, (pach_end + 1) * row_w
        f = (df.s_um - s0) / (s1 - s0)
        zone[(f >= 0) & (f < 1 / 3)] = "early"
        zone[(f >= 1 / 3) & (f < 2 / 3)] = "mid"
        zone[(f >= 2 / 3) & (f <= 1)] = "late"
        zone[df.s_um < s0] = "pre"
        zone[df.s_um > s1] = "post"
        span = float(s1 - s0)
    df["zone"] = zone
    lm = {"image_id": iid, "n_nuclei": int(len(df)), "n_nuclei_before_component_filter": int(n_all),
          "component_filter": comp, "orientation": ori, "flipped": bool(flip),
          "crescent_median_pos": None if not np.isfinite(cres_pos) else round(cres_pos, 3),
          "polarity_contradiction_corrected": contradiction,
          "lamin_fallback": fb, "axis_len_um": float(L), "refs": {k: float(v) for k, v in refs.items()},
          "n_mitotic": int(df.mitotic.sum()), "n_crescent": int(df.crescent.sum()), "n_post": int(df.post.sum()),
          "row_width_um": row_w, "tz_start_row": tz_start, "pach_start_row": pach_start, "pach_end_row": pach_end,
          "pach_span_um": span}
    df.to_csv(os.path.join(OUTDIR, f"{iid}_nuclei.csv"), index=False)
    rt.to_csv(os.path.join(OUTDIR, f"{iid}_rows.csv"), index=False)
    json.dump(lm, open(os.path.join(OUTDIR, f"{iid}_landmarks.json"), "w"), indent=1)
    print(json.dumps(lm))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.15, 1])
    a = fig.add_subplot(gs[0, 0])
    for c, g in df.groupby("cls"):
        a.scatter(g.cx, g.cy, s=14, color=CMAP[c], label=f"{c} (n={len(g)})")
    a.invert_yaxis()
    a.set_aspect("equal")
    a.legend(frameon=False, fontsize=8)
    a.annotate("distal", (df.iloc[0].cx, df.iloc[0].cy), fontsize=9, fontweight="bold")
    a.set_title("nucleus classes (xy, µm)", fontsize=10)
    b = fig.add_subplot(gs[0, 1:])
    base = np.zeros(len(rt))
    for key, c in [("n_mit", "mitotic"), ("n_tz", "TZ"), ("n_pach", "pachytene"), ("n_post", "post")]:
        b.bar(rt.row, rt[key], bottom=base, color=CMAP[c], label=c)
        base = base + rt[key].to_numpy()
    ytop = b.get_ylim()[1]
    for x, txt, ha in [(tz_start, "TZ established", "left"), (pach_start, "pachytene start", "left"),
                       (pach_end, "pachytene end", "right")]:
        if x is not None:
            xx = x - 0.5 if ha == "left" else x + 0.5
            b.axvline(xx, c="k", lw=1.4, ls="--")
            b.text(xx, ytop * 0.96, " " + txt + " ", fontsize=8.5, ha=ha, va="top")
    b.set_xlabel(f"row (distal to proximal, {row_w:.1f} µm each)")
    b.set_ylabel("nuclei per row")
    b.legend(frameon=False, fontsize=8)
    b.set_title("per-row composition and landmarks", fontsize=10)
    R = int(np.ceil(3.2 / SP[1]))
    for k, c in enumerate(["mitotic", "TZ", "pachytene", "post"]):
        ax = fig.add_subplot(gs[1, k])
        ax.axis("off")
        g = df[df.cls == c]
        if c == "TZ":
            g = g.sort_values("gap", ascending=False)
        if c == "pachytene" and len(g):
            g = g.sample(min(8, len(g)), random_state=0)
        tiles = [t for t in (thumb(dapi, syp, r, R) for r in g.head(8).itertuples()) if t is not None]
        if tiles:
            blank = np.zeros_like(tiles[0])
            rows_ = [np.concatenate(tiles[i:i + 4] + [blank] * (4 - len(tiles[i:i + 4])), 1)
                     for i in range(0, len(tiles), 4)]
            ax.imshow(np.concatenate(rows_, 0))
        ax.set_title(f"{c} examples (n={len(df[df.cls == c])})", fontsize=10, color=CMAP[c])
    fig.suptitle(f"{iid[-28:]}   staging QC   (pachytene rows {pach_start} to {pach_end}, "
                 f"{'?' if span is None else round(span)} µm)", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, f"{iid}_qc.png"), dpi=130)
    R2 = int(np.ceil(2.8 / SP[1]))
    PER = 5
    nrows = int(df.row.max()) + 1
    strip = np.zeros((nrows * 2 * R2, (PER + 1) * 2 * R2, 3), np.float32)
    for ri in range(nrows):
        for k, r in enumerate(df[df.row == ri].sort_values("s_um").head(PER).itertuples()):
            t = thumb(dapi, syp, r, R2)
            if t is None:
                continue
            colr = CMAP[r.cls]
            t[:2, :] = colr
            t[-2:, :] = colr
            t[:, :2] = colr
            t[:, -2:] = colr
            strip[ri * 2 * R2:(ri + 1) * 2 * R2, (k + 1) * 2 * R2:(k + 2) * 2 * R2] = t
    fig2, ax2 = plt.subplots(figsize=(8, max(8, nrows * 0.32)), facecolor="black")
    ax2.imshow(strip)
    ax2.axis("off")
    for ri in range(nrows):
        ax2.text(R2 * 0.9, ri * 2 * R2 + R2, f"{ri}", color="w", fontsize=7, ha="right", va="center")
    for x, txt, colr in [(tz_start, "TZ established", "orange"), (pach_start, "pachytene start", "lime")]:
        if x is not None:
            ax2.axhline(x * 2 * R2, color=colr, lw=1.2)
            ax2.text(2 * R2, x * 2 * R2 - 3, "auto: " + txt, color=colr, fontsize=7, va="bottom")
    if pach_end is not None:
        ax2.axhline((pach_end + 1) * 2 * R2, color="magenta", lw=1.2)
        ax2.text(2 * R2, (pach_end + 1) * 2 * R2 + 3, "auto: pachytene end", color="magenta", fontsize=7, va="top")
    ax2.set_title(f"{iid[-20:]}  rows distal(top) to proximal(bottom), {row_w:.1f} µm/row; "
                  "border gray=mitotic orange=TZ blue=pachytene purple=post", color="w", fontsize=8)
    fig2.tight_layout()
    fig2.savefig(os.path.join(OUTDIR, f"{iid}_strip.png"), dpi=110, facecolor="black")
    print("QC + strip written")


if __name__ == "__main__":
    main(sys.argv[1])
