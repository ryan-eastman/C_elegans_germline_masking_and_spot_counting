"""Zoom on a region of one gonad: DAPI and LMN-1 (max projection + one plane) with ALL-depth mask edges
(blue DAPI label, red lamin watershed) and the label id printed on every object whose volume or ring
score is unusual, so blobs masked inside the germline can be checked one by one.
Usage: python qc_mask_zoom.py <image_id> y0 y1 x0 x1 [z_um]      (region in um, crop coordinates)"""
import os
import sys

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import chload
import pc_lamin_worker as W

OUT = r"C:/Users/ryane/coloc_analysis/mask_audit"
SP = chload.SP


def main(iid, y0, y1, x0, x1, z_um=None):
    rd = W.find_run(iid)
    lab = tifffile.imread(os.path.join(rd, f"{iid}__nuclei_labels.tif"))
    csv = os.path.join(rd, f"{iid}__nuclei.csv")
    nc = pd.read_csv(csv)
    germ = [int(x) for x in nc[nc.in_germline.astype(bool)].nucleus_id]
    gn = np.isin(lab, germ)
    obj = ndi.find_objects(gn.astype(np.uint8))[0]
    sl = tuple(slice(max(0, o.start - W.PAD), min(dim, o.stop + W.PAD)) for o, dim in zip(obj, gn.shape))
    ch = chload.load(W.find_nd2(iid, csv))
    lamin = ch["lamin"][sl].astype(np.float32)
    dapi = ch["dapi"][sl].astype(np.float32)
    del ch
    lab_c = lab[sl]
    nuc_lam, nuc_dapi, _ = W.lamin_nuclei(lab_c, germ, lamin)
    lam_lab = np.where(nuc_lam, lab_c, 0)      # approximate: label id under the lamin mask

    r0, r1 = int(y0 / SP[1]), int(y1 / SP[1])
    c0, c1 = int(x0 / SP[2]), int(x1 / SP[2])
    z = int(z_um / SP[0]) if z_um is not None else lamin.shape[0] // 2
    sub = (slice(None), slice(r0, r1), slice(c0, c1))
    lam_s, dapi_s, lab_s = lamin[sub], dapi[sub], lab_c[sub]
    nd_s, nl_s = nuc_dapi[sub], nuc_lam[sub]

    aud = pd.read_csv(os.path.join(OUT, f"{iid}_mask_audit.csv"))
    L = aud[aud.kind == "germ_label"].set_index("id")
    ids = [i for i in np.unique(lab_s[nd_s]) if i > 0]
    rows = []
    for i in ids:
        m = lab_s == i
        if m.sum() == 0:
            continue
        cz, cy, cx = ndi.center_of_mass(m)
        # shape: bounding-box fill and elongation from the max projection
        mp = m.max(0)
        ys, xs = np.where(mp)
        ext = max(np.ptp(ys), np.ptp(xs)) * SP[1]
        r = L.loc[i] if i in L.index else None
        rows.append({"id": int(i), "vol_um3": round(float(m.sum()) * float(SP.prod()), 1),
                     "extent_um": round(float(ext), 1), "y_px": cy, "x_px": cx, "z_um": round(cz * SP[0], 1),
                     "ring_ratio": None if r is None else r.ring_ratio,
                     "shell_over_thr": None if r is None else r.shell_over_thr})
    t = pd.DataFrame(rows)
    flag = t[(t.vol_um3 > 70) | (t.vol_um3 < 8) | (t.extent_um > 9) | (t.ring_ratio > 0.95) | (t.shell_over_thr < 0.8)]
    print(f"{len(t)} germline labels in the window; flagged {len(flag)}:")
    print(flag.sort_values("vol_um3", ascending=False).drop(columns=["y_px", "x_px"]).to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from skimage.segmentation import find_boundaries

    def norm(a, lo=30, hi=99.7):
        a0, a1 = np.percentile(a, lo), np.percentile(a, hi)
        return np.clip((a - a0) / (a1 - a0 + 1e-9), 0, 1)

    e_d_all, e_l_all = find_boundaries(nd_s.max(0), mode="inner"), find_boundaries(nl_s.max(0), mode="inner")
    e_d_z, e_l_z = find_boundaries(nd_s[z], mode="inner"), find_boundaries(nl_s[z], mode="inner")
    panels = [(norm(dapi_s.max(0)), e_d_all, e_l_all, "DAPI max projection, all-depth edges"),
              (norm(lam_s.max(0)), e_d_all, e_l_all, "LMN-1 max projection, all-depth edges"),
              (norm(dapi_s[z], 5, 99.8), e_d_z, e_l_z, f"DAPI plane z = {z * SP[0]:.1f} um, edges on this plane"),
              (norm(lam_s[z], 5, 99.8), e_d_z, e_l_z, f"LMN-1 plane z = {z * SP[0]:.1f} um, edges on this plane")]
    h, w = lam_s.shape[1:]
    fig, ax = plt.subplots(2, 2, figsize=(2 * 9 * w / max(h, w) + 1, 2 * 9 * h / max(h, w) + 1), facecolor="black")
    for a, (base, ed, el, ttl) in zip(ax.ravel(), panels):
        rgb = np.stack([base] * 3, -1)
        rgb[ed] = (0.2, 0.45, 1.0)
        rgb[el] = (1.0, 0.15, 0.15)
        a.imshow(rgb)
        a.axis("off")
        a.set_title(ttl, color="w", fontsize=10)
        for _, r in flag.iterrows():
            a.text(r.x_px, r.y_px, str(r.id), color="yellow", fontsize=7, ha="center", va="center")
    fig.suptitle(f"{iid}  window y {y0}-{y1} um, x {x0}-{x1} um; yellow ids = flagged (see console)", color="w")
    plt.tight_layout()
    out = os.path.join(OUT, f"{iid}_zoom_y{int(y0)}-{int(y1)}_x{int(x0)}-{int(x1)}.png")
    fig.savefig(out, dpi=110, facecolor="black")
    print("wrote", out, flush=True)


if __name__ == "__main__":
    a = sys.argv
    main(a[1], float(a[2]), float(a[3]), float(a[4]), float(a[5]), float(a[6]) if len(a) > 6 else None)
