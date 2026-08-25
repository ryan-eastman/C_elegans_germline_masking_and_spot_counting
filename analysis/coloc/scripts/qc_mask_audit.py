"""Mask audit for one gonad (Ryan, 2026-08-24, after seeing the single-plane QC overlay).
Two independent questions about the nuclei that feed the partition coefficient:
  1. MISSED: lamin rings with no DAPI/Cellpose label at ANY depth. Nuclei are found from the LMN-1
     channel alone (ring -> filled interior, 3D components, 5-150 um3) and each is scored by the fraction
     of its voxels covered by a germline label.
  2. JUNK: germline-flagged labels that have no lamin ring around them (sperm, somatic, debris): ring
     score = mean smoothed lamin in a 0.4 um shell just outside the label, relative to the ring threshold.
Writes mask_audit/<iid>_mask_audit.csv (per lamin candidate and per label) and a 2-panel max-projection
overlay with ALL-depth mask edges (blue DAPI, red lamin watershed, yellow = lamin nuclei with no label,
magenta fill = germline label with no lamin ring).
Usage: python qc_mask_audit.py <image_id>"""
import os
import sys

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.morphology import disk

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import chload
import pc_lamin_worker as W

OUT = r"C:/Users/ryane/coloc_analysis/mask_audit"
os.makedirs(OUT, exist_ok=True)
SP = chload.SP
VVOL = float(SP.prod())


def main(iid):
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
    gn_c = gn[sl]
    nuc_lam, nuc_dapi, fell_back = W.lamin_nuclei(lab_c, germ, lamin)

    # --- 1. lamin-only nuclei -------------------------------------------------------------------
    sm = ndi.gaussian_filter(lamin, sigma=0.25 / SP)
    thr = float(threshold_otsu(sm[sm > np.percentile(sm, 50)]))
    ring = sm > thr
    se = disk(int(round(0.35 / SP[1])))
    interior = np.zeros_like(ring)
    for z in range(ring.shape[0]):
        closed = ndi.binary_closing(ring[z], structure=se)
        interior[z] = ndi.binary_fill_holes(closed) & ~closed
    interior = ndi.binary_opening(interior, structure=np.ones((3, 3, 3), bool))
    cl, ncl = ndi.label(interior)
    vol = np.bincount(cl.ravel()) * VVOL
    keep = np.where((vol >= 5.0) & (vol <= 150.0))[0]
    keep = keep[keep != 0]
    cand = []
    for k in keep:
        m = cl == k
        cz, cy, cx = ndi.center_of_mass(m)
        cand.append({"kind": "lamin_candidate", "id": int(k), "vol_um3": round(float(vol[k]), 1),
                     "z_um": round(cz * SP[0], 1), "y_um": round(cy * SP[1], 1), "x_um": round(cx * SP[2], 1),
                     "cov_germ_label": round(float(gn_c[m].mean()), 3),
                     "cov_any_label": round(float((lab_c[m] > 0).mean()), 3),
                     "cov_lamin_mask": round(float(nuc_lam[m].mean()), 3)})
    cand = pd.DataFrame(cand)
    missed = cand[cand.cov_germ_label < 0.15] if len(cand) else cand
    missed_any = missed[missed.cov_any_label < 0.15] if len(missed) else missed

    # --- 2. ring score of every germline label --------------------------------------------------
    ids = np.array(germ)
    dt_out, idx = ndi.distance_transform_edt(~gn_c, sampling=tuple(SP), return_indices=True)
    shell = (dt_out > 0) & (dt_out <= 0.4)
    near = lab_c[idx[0], idx[1], idx[2]]
    shell_lab = np.where(shell, near, 0)
    inside_mean = np.asarray(ndi.mean(sm, lab_c * gn_c, ids))
    shell_mean = np.asarray(ndi.mean(sm, shell_lab, ids))
    lvol = np.bincount(lab_c.ravel(), minlength=int(ids.max()) + 1)[ids] * VVOL
    coms = ndi.center_of_mass(gn_c, lab_c, ids)
    labs = pd.DataFrame({"kind": "germ_label", "id": ids, "vol_um3": np.round(lvol, 1),
                         "z_um": [round(c[0] * SP[0], 1) for c in coms],
                         "y_um": [round(c[1] * SP[1], 1) for c in coms],
                         "x_um": [round(c[2] * SP[2], 1) for c in coms],
                         "ring_shell": np.round(shell_mean, 1), "ring_inside": np.round(inside_mean, 1),
                         "ring_ratio": np.round(shell_mean / np.maximum(inside_mean, 1), 3),
                         "shell_over_thr": np.round(shell_mean / thr, 3)})
    # no envelope: shell is flat against the interior (ratio ~1) AND well below the ring threshold
    junk = labs[(labs.ring_ratio > 0.97) & (labs.shell_over_thr < 0.75)]

    pd.concat([cand, labs]).to_csv(os.path.join(OUT, f"{iid}_mask_audit.csv"), index=False)
    s = chload.parse_iid(iid)["short"]
    ybins = np.arange(0, lamin.shape[1] * SP[1] + 20, 20)
    print(f"{s}: crop {lamin.shape}, otsu ring thr {thr:.0f}; germline labels {len(ids)} (fallback {fell_back})")
    print(f"  lamin-only nuclei found: {len(cand)}; NO germline label: {len(missed)} "
          f"({100 * len(missed) / max(1, len(cand)):.0f}%), of which no label at all: {len(missed_any)}")
    if len(missed):
        h, e = np.histogram(missed.y_um, bins=ybins)
        print("  missed by y (um):", {f"{int(a)}-{int(b)}": int(c) for a, b, c in zip(e[:-1], e[1:], h) if c})
    print(f"  germline labels with NO lamin ring (flat shell, < 0.75 x thr): {len(junk)} "
          f"({100 * len(junk) / max(1, len(ids)):.0f}%), volume {junk.vol_um3.sum():.0f} um3 of {lvol.sum():.0f}")
    if len(junk):
        h, e = np.histogram(junk.y_um, bins=ybins)
        print("  junk by y (um):", {f"{int(a)}-{int(b)}": int(c) for a, b, c in zip(e[:-1], e[1:], h) if c})

    # --- overlay: max projections, ALL-depth edges ----------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from skimage.segmentation import find_boundaries

    def norm(a):
        lo, hi = np.percentile(a, 30), np.percentile(a, 99.7)
        return np.clip((a - lo) / (hi - lo + 1e-9), 0, 1)

    lam_mp = norm(lamin.max(0))
    dapi_mp = norm(dapi.max(0))
    e_dapi = find_boundaries(nuc_dapi.max(0), mode="inner")
    e_lam = find_boundaries(nuc_lam.max(0), mode="inner")
    miss_mask = np.isin(cl, missed.id.to_numpy()).max(0) if len(missed) else np.zeros(lam_mp.shape, bool)
    e_miss = find_boundaries(miss_mask, mode="inner")
    junk_mask = np.isin(lab_c, junk.id.to_numpy()).max(0) if len(junk) else np.zeros(lam_mp.shape, bool)
    jfill = junk_mask & ~e_lam & ~e_dapi
    fig, ax = plt.subplots(1, 2, figsize=(20, 11), facecolor="black")
    for a, base, ttl in [(ax[0], lam_mp, "LMN-1 max projection"), (ax[1], dapi_mp, "DAPI max projection")]:
        rgb = np.stack([base, base, base], -1)
        rgb[e_dapi] = (0.2, 0.4, 1.0)
        rgb[e_lam] = (1.0, 0.15, 0.15)
        rgb[e_miss] = (1.0, 0.9, 0.0)
        rgb[jfill] = rgb[jfill] * 0.5 + np.array([0.5, 0.0, 0.5])
        a.imshow(rgb)
        a.axis("off")
        a.set_title(f"{ttl}: edges over ALL depths (blue DAPI, red lamin watershed); yellow = lamin nucleus "
                    f"with no label ({len(missed)}); magenta fill = label with no lamin ring ({len(junk)})",
                    color="w", fontsize=9)
    fig.suptitle(iid, color="w", fontsize=11)
    plt.tight_layout()
    fig.savefig(os.path.join(OUT, f"{iid}_mask_audit.png"), dpi=110, facecolor="black")
    print(f"  wrote {OUT}/{iid}_mask_audit.png", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
