"""Pachytene staging from morphology, per the lab definition (Ryan, 2026-08-22):
  * TZ (transition-zone) nucleus = DNA in a polarized / crescent morphology
  * Pachytene START = first row that does not contain more than TZ_MAX (=2) TZ nuclei
  * Pachytene END   = last row containing all pachytene nuclei with at most 1 diplotene nucleus
  * Pachytene region split into 3 equal-length zones: early / mid / late
Usage: python crescent_axis.py <image_id>

Nucleus features are measured INSIDE THE LAMIN ENVELOPE (watershed on LMN-1 seeded by DAPI nuclei):
  offset   = |DAPI intensity centroid - envelope geometric centroid| / equivalent radius  (crescent: high)
  occ      = fraction of envelope volume occupied by chromatin (DAPI above gonad-wide Otsu)  (crescent: low)
  vol_env  = envelope volume um^3                                                           (diplotene: high)
  syp_cov  = fraction of envelope volume with SYP-3 above gonad-wide Otsu                  (diplotene: low)
Axis ordering WITHOUT a hand trace: geodesic distance along a k-nearest-neighbour graph of germline
nucleus centroids (follows the gonad around the U-bend); distal end = end with the smallest nuclei,
cross-checked against where the crescents sit. Rows = bins of one median nuclear diameter along the axis.
Outputs (coloc_analysis/staging/<iid>_*): nuclei.csv (per-nucleus features, class, row, zone),
rows.csv, landmarks.json, qc.png (axis map, per-row counts with landmarks, thumbnail montage per class)."""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree, shortest_path
from scipy.spatial import cKDTree
from skimage.filters import threshold_otsu
from skimage.segmentation import watershed

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
import chload

SP = chload.SP
SPACING = tuple(SP)
VVOL = float(SP.prod())
PAD = 30
LAM_BAND_UM, LAM_SMOOTH_UM = 2.0, 0.15
TZ_MAX = 2            # "not more than 1-2 TZ nuclei" -> a pachytene row allows <= 2 crescents
DIP_MAX = 1           # "occasional single diplotene nucleus"
PERSIST = 2           # a landmark must hold for this many consecutive rows (noise guard)
KNN = 10
RUN_DIRS = [r"C:/Users/ryane/ccw77_fullres", r"C:/Users/ryane/ccw77_0622"]
E_DIRS = [r"E:/20260708_ccw77_IF_pgl1_syp3_LMN1", r"E:/20260622_ccw77_IF_syp3_pgl1_LMN1_dapi"]
NAS_IP = "129.82.125.105"
OUTDIR = r"C:/Users/ryane/coloc_analysis/staging"
os.makedirs(OUTDIR, exist_ok=True)


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
    fp = fp.replace("\\\\bmb-ckc-nas", "//" + NAS_IP).replace("\\", "/")
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


def nucleus_features(lam, lab_c, ids, dapi, syp):
    dapi_bs = np.clip(dapi - np.percentile(dapi, 5), 0, None)
    syp_bs = np.clip(syp - np.percentile(syp, 5), 0, None)
    env = lam > 0
    T_d = float(threshold_otsu(dapi_bs[env]))
    T_s = float(threshold_otsu(syp_bs[env]))
    objs = ndi.find_objects(lam)
    rows = []
    for nid in ids:
        sl = objs[nid - 1] if nid - 1 < len(objs) else None
        if sl is None:
            continue
        m = lam[sl] == nid
        n = int(m.sum())
        if n < 50:
            continue
        zz, yy, xx = np.nonzero(m)
        coords = np.stack([zz, yy, xx], 1) * SP + np.array([s.start for s in sl]) * SP
        w = dapi_bs[sl][m].astype(np.float64)
        c_geo = coords.mean(0)
        c_dapi = (coords * w[:, None]).sum(0) / max(w.sum(), 1e-9)
        vol = n * VVOL
        r_eq = (3 * vol / (4 * np.pi)) ** (1 / 3)
        bright = w > T_d
        # chromatin (bright) centroid offset is more robust to diffuse background than intensity-weighted
        c_chr = coords[bright].mean(0) if bright.sum() > 20 else c_dapi
        dm = (lab_c[sl] == nid)                       # Cellpose DAPI mask of this nucleus
        rows.append({"nucleus_id": nid, "cz": c_geo[0], "cy": c_geo[1], "cx": c_geo[2],
                     "vol_env": vol, "r_eq": r_eq,
                     "offset": float(np.linalg.norm(c_chr - c_geo) / r_eq),
                     "occ": float(bright.mean()),
                     "dapi_frac": float(dm.sum() / n),
                     "compact": float(np.percentile(w, 90) / max(w.mean(), 1e-9)),
                     "dapi_cv": float(w.std() / max(w.mean(), 1e-9)),
                     "syp_cov": float((syp_bs[sl][m] > T_s).mean())})
    return pd.DataFrame(rows)


def geodesic_axis(df):
    P = df[["cz", "cy", "cx"]].to_numpy()
    n = len(P)
    tree = cKDTree(P)
    d, j = tree.query(P, k=min(KNN + 1, n))
    src = np.repeat(np.arange(n), d.shape[1] - 1)
    dst = j[:, 1:].ravel()
    wgt = d[:, 1:].ravel()
    G = coo_matrix((wgt, (src, dst)), shape=(n, n)).tocsr()
    G = G.maximum(G.T)
    ncomp, _ = connected_components(G, directed=False)
    if ncomp > 1:  # stitch disconnected clusters with the global MST
        D = np.linalg.norm(P[:, None] - P[None], axis=2)
        mst = minimum_spanning_tree(D).tocsr()
        G = G.maximum(mst).maximum(mst.T)
    D0 = shortest_path(G, directed=False, indices=[0])[0]
    A = int(np.argmax(D0))
    DA = shortest_path(G, directed=False, indices=[A])[0]
    B = int(np.argmax(DA))
    DB = shortest_path(G, directed=False, indices=[B])[0]
    L = DA[B]
    medA = df.vol_env[DA <= 0.15 * L].median()
    medB = df.vol_env[DB <= 0.15 * L].median()
    s = DA if medA <= medB else DB      # distal end = smaller nuclei
    return s, L, ("A" if medA <= medB else "B"), (medA, medB)


def main(iid):
    rd = find_run(iid)
    lab = tifffile.imread(os.path.join(rd, f"{iid}__nuclei_labels.tif"))
    csv = os.path.join(rd, f"{iid}__nuclei.csv")
    nc = pd.read_csv(csv)
    germ = [int(x) for x in nc[nc.in_germline.astype(bool)].nucleus_id]
    gn = np.isin(lab, germ)
    obj = ndi.find_objects(gn.astype(np.uint8))[0]
    sl = tuple(slice(max(0, o.start - PAD), min(dim, o.stop + PAD)) for o, dim in zip(obj, gn.shape))
    ch = chload.load(find_nd2(iid, csv))
    dapi = ch["dapi"][sl].astype(np.float32)
    syp = ch["syp"][sl].astype(np.float32)
    lamin = ch["lamin"][sl]
    del ch
    lab_c = lab[sl]
    lam, fb = lamin_labels(lab_c, germ, lamin)
    df = nucleus_features(lam, lab_c, germ, dapi, syp)
    s, L, end, meds = geodesic_axis(df)
    df["s_um"] = s
    df["s_norm"] = s / L
    df = df.sort_values("s_um").reset_index(drop=True)

    # merged / oversized envelopes cannot be staged by shape
    v_med = df.vol_env.median()
    df["merged"] = df.vol_env > 2.2 * v_med
    # crescent (TZ) call: chromatin pushed to one side (offset) AND a dark half (low occupancy)
    df["cres_score"] = df.offset * (1 - df.occ)
    good = ~df.merged
    T_cs = max(float(threshold_otsu(df.cres_score[good].to_numpy())), 0.08)
    df["crescent"] = good & (df.cres_score > T_cs) & (df.offset > 0.18)
    T_off = T_cs
    # orientation cross-check: crescents belong in the distal half
    flipped = False
    if df.crescent.sum() >= 5 and df.s_norm[df.crescent].median() > 0.6:
        df["s_um"] = L - df.s_um
        df["s_norm"] = 1 - df.s_norm
        df = df.sort_values("s_um").reset_index(drop=True)
        flipped = True
    # post-pachytene call (diplotene / condensing): SC disassembly + chromatin compaction, proximal half only.
    # reference = mid-axis non-crescent nuclei
    core = df[(~df.crescent) & (~df.merged) & (df.s_norm > 0.25) & (df.s_norm < 0.75)]
    c_ref, k_ref = core.syp_cov.median(), core.compact.median()
    df["diplotene"] = (~df.crescent) & (df.s_norm > 0.5) & (df.syp_cov < 0.5 * c_ref) & (df.compact > 1.15 * k_ref)
    df["cls"] = np.where(df.crescent, "TZ", np.where(df.diplotene, "diplotene", "pachytene"))

    # rows = bins of one median nuclear diameter along the axis
    row_w = float(2 * df.r_eq.median())
    df["row"] = (df.s_um // row_w).astype(int)
    rt = df.groupby("row").agg(n=("nucleus_id", "size"), n_tz=("crescent", "sum"),
                               n_dip=("diplotene", "sum"), s0=("s_um", "min")).reset_index()
    rt["n_pach"] = rt.n - rt.n_tz - rt.n_dip
    pach_start = pach_end = None
    tz_rows = rt.index[rt.n_tz >= 1].tolist()
    if tz_rows:
        for i in range(tz_rows[0], len(rt) - PERSIST + 1):
            if all(rt.n_tz.iloc[i + k] <= TZ_MAX for k in range(PERSIST)):
                pach_start = int(rt.row.iloc[i])
                break
    if pach_start is not None:
        after = rt[rt.row >= pach_start].reset_index(drop=True)
        end_row = None
        for i in range(len(after) - PERSIST + 1):
            if all(after.n_dip.iloc[i + k] > DIP_MAX for k in range(PERSIST)):
                end_row = int(after.row.iloc[i]) - 1
                break
        pach_end = end_row if end_row is not None else int(after.row.iloc[-1])
    zone = np.full(len(df), "none", dtype=object)
    if pach_start is not None and pach_end is not None and pach_end > pach_start:
        s0, s1 = pach_start * row_w, (pach_end + 1) * row_w
        f = (df.s_um - s0) / (s1 - s0)
        zone[(f >= 0) & (f < 1 / 3)] = "early"
        zone[(f >= 1 / 3) & (f < 2 / 3)] = "mid"
        zone[(f >= 2 / 3) & (f <= 1)] = "late"
        zone[df.s_um < s0] = "pre"
        zone[df.s_um > s1] = "post"
    df["zone"] = zone
    span = None if pach_start is None else float((pach_end + 1 - pach_start) * row_w)
    lm = {"image_id": iid, "n_nuclei": int(len(df)), "lamin_fallback": fb, "axis_len_um": float(L),
          "distal_end": end, "end_vol_medians": [float(meds[0]), float(meds[1])], "flipped_by_crescents": flipped,
          "offset_threshold": T_off, "n_crescent": int(df.crescent.sum()), "n_diplotene": int(df.diplotene.sum()),
          "row_width_um": row_w, "pach_start_row": pach_start, "pach_end_row": pach_end, "pach_span_um": span}
    df.to_csv(os.path.join(OUTDIR, f"{iid}_nuclei.csv"), index=False)
    rt.to_csv(os.path.join(OUTDIR, f"{iid}_rows.csv"), index=False)
    json.dump(lm, open(os.path.join(OUTDIR, f"{iid}_landmarks.json"), "w"), indent=1)
    print(json.dumps(lm))

    # ---------------- QC figure ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1])
    cmap = {"TZ": "#e67e22", "pachytene": "#2c6fbb", "diplotene": "#8e44ad"}
    a = fig.add_subplot(gs[0, 0])
    for c, g in df.groupby("cls"):
        a.scatter(g.cx, g.cy, s=14, color=cmap[c], label=f"{c} (n={len(g)})")
    a.invert_yaxis()
    a.set_aspect("equal")
    a.legend(frameon=False, fontsize=8)
    d0 = df.iloc[0]
    a.annotate("distal", (d0.cx, d0.cy), color="k", fontsize=9, fontweight="bold")
    a.set_title("nucleus classes (xy, µm)", fontsize=10)
    b = fig.add_subplot(gs[0, 1:])
    b.bar(rt.row, rt.n_pach, color=cmap["pachytene"], label="pachytene")
    b.bar(rt.row, rt.n_tz, bottom=rt.n_pach, color=cmap["TZ"], label="TZ / crescent")
    b.bar(rt.row, rt.n_dip, bottom=rt.n_pach + rt.n_tz, color=cmap["diplotene"], label="diplotene")
    ytop = b.get_ylim()[1]
    if pach_start is not None:
        b.axvline(pach_start - 0.5, c="k", lw=1.5, ls="--")
        b.text(pach_start, ytop * 0.95, " pachytene start", fontsize=9)
    if pach_end is not None:
        b.axvline(pach_end + 0.5, c="k", lw=1.5, ls="--")
        b.text(pach_end, ytop * 0.85, "pachytene end ", ha="right", fontsize=9)
    b.set_xlabel(f"row (distal to proximal, {row_w:.1f} µm each)")
    b.set_ylabel("nuclei per row")
    b.legend(frameon=False, fontsize=8)
    b.set_title("per-row composition and landmarks", fontsize=10)
    R = int(np.ceil(3.2 / SP[1]))
    for k, c in enumerate(["TZ", "pachytene", "diplotene"]):
        ax = fig.add_subplot(gs[1, k])
        ax.axis("off")
        g = df[df.cls == c]
        if c == "TZ":
            g = g.sort_values("offset", ascending=False)
        if c == "diplotene":
            g = g.sort_values("vol_env", ascending=False)
        tiles = []
        for r in g.head(8).itertuples():
            zc, yc, xc = int(r.cz / SP[0]), int(r.cy / SP[1]), int(r.cx / SP[2])
            dd = dapi[zc, max(0, yc - R):yc + R, max(0, xc - R):xc + R]
            ss = syp[zc, max(0, yc - R):yc + R, max(0, xc - R):xc + R]
            if dd.shape != (2 * R, 2 * R):
                continue
            dn = np.clip((dd - np.percentile(dd, 5)) / (np.percentile(dd, 99.5) - np.percentile(dd, 5) + 1e-9), 0, 1)
            sn = np.clip((ss - np.percentile(ss, 30)) / (np.percentile(ss, 99.5) - np.percentile(ss, 30) + 1e-9), 0, 1)
            tiles.append(np.stack([np.maximum(dn, sn), dn, dn], -1))
        if tiles:
            blank = np.zeros_like(tiles[0])
            rows_ = [np.concatenate(tiles[i:i + 4] + [blank] * (4 - len(tiles[i:i + 4])), 1) for i in range(0, len(tiles), 4)]
            ax.imshow(np.concatenate(rows_, 0))
        ax.set_title(f"{c}: top examples (DAPI gray, SYP-3 red)", fontsize=10, color=cmap[c])
    fig.suptitle(f"{iid[-28:]}   staging QC   (pachytene rows {pach_start} to {pach_end}, "
                 f"{'?' if span is None else round(span)} µm)", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, f"{iid}_qc.png"), dpi=130)
    print("QC written")

    # ---------------- verification strip: one line per row, up to 4 nuclei ----------------
    R2 = int(np.ceil(2.8 / SP[1]))
    PER = 4
    nrows = int(df.row.max()) + 1
    strip = np.zeros((nrows * 2 * R2, (PER + 1) * 2 * R2, 3), np.float32)
    for ri in range(nrows):
        g = df[df.row == ri].sort_values("s_um").head(PER)
        for k, r in enumerate(g.itertuples()):
            zc, yc, xc = int(r.cz / SP[0]), int(r.cy / SP[1]), int(r.cx / SP[2])
            dd = dapi[zc, max(0, yc - R2):yc + R2, max(0, xc - R2):xc + R2]
            ss = syp[zc, max(0, yc - R2):yc + R2, max(0, xc - R2):xc + R2]
            if dd.shape != (2 * R2, 2 * R2):
                continue
            dn = np.clip((dd - np.percentile(dd, 5)) / (np.percentile(dd, 99.5) - np.percentile(dd, 5) + 1e-9), 0, 1)
            sn = np.clip((ss - np.percentile(ss, 30)) / (np.percentile(ss, 99.5) - np.percentile(ss, 30) + 1e-9), 0, 1)
            tile = np.stack([np.maximum(dn, sn), dn, dn], -1)
            # class tint on the tile border
            colr = {"TZ": (1, 0.5, 0.1), "pachytene": (0.2, 0.45, 0.75), "diplotene": (0.55, 0.27, 0.68)}[r.cls]
            tile[:2, :] = colr; tile[-2:, :] = colr; tile[:, :2] = colr; tile[:, -2:] = colr
            strip[ri * 2 * R2:(ri + 1) * 2 * R2, (k + 1) * 2 * R2:(k + 2) * 2 * R2] = tile
    h_in = max(8, nrows * 0.32)
    fig2, ax2 = plt.subplots(figsize=(7.5, h_in), facecolor="black")
    ax2.imshow(strip); ax2.axis("off")
    for ri in range(nrows):
        y = ri * 2 * R2 + R2
        ax2.text(R2 * 0.9, y, f"{ri}", color="w", fontsize=7, ha="right", va="center")
        if ri == pach_start:
            ax2.axhline(ri * 2 * R2, color="lime", lw=1.2); ax2.text(2 * R2, ri * 2 * R2 - 3, "auto: pachytene start", color="lime", fontsize=7, va="bottom")
        if pach_end is not None and ri == pach_end:
            ax2.axhline((ri + 1) * 2 * R2, color="magenta", lw=1.2); ax2.text(2 * R2, (ri + 1) * 2 * R2 + 3, "auto: pachytene end", color="magenta", fontsize=7, va="top")
    ax2.set_title(f"{iid[-20:]}  rows distal(top) to proximal(bottom), {row_w:.1f} µm/row\nborder: orange=TZ call, blue=pachytene, purple=post", color="w", fontsize=8)
    fig2.tight_layout()
    fig2.savefig(os.path.join(OUTDIR, f"{iid}_strip.png"), dpi=110, facecolor="black")
    print("strip written")


if __name__ == "__main__":
    main(sys.argv[1])
