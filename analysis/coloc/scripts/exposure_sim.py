"""Does SYP exposure time bias the partition coefficient? (Ryan's concern, 2026-08-24)
Exposures are condition-correlated in this dataset (HS mostly 200-300 ms, noHS mostly 90-100 ms on 545).
A ratio is exposure-invariant only if the additive camera offset is subtracted correctly; any residual
compresses LOW-exposure gonads toward PC = 1, which would inflate an HS-vs-noHS difference.

Test: take every clean-13 gonad imaged at 200 ms on 545 and synthetically re-expose it to 90 ms:
    syp_90 = offset + (syp_200 - offset) * (90/200) + noise
with offset = camera dark level (0.5th percentile of extra-worm voxels, ~100 counts on the Kinetix) and
noise = Poisson on the scaled photon signal (CMS conversion gain ~ 1 e-/count assumed) + Gaussian read
noise (1.2 counts rms). Then recompute PC_dm, z-shift and PC_specific with the standard pipeline and
report the shift. Everything else (granules, envelopes, bg rule) is held identical."""
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

sys.path.insert(0, r"C:\Users\ryane\coloc_analysis")
sys.path.insert(0, r"C:\Users\ryane\coloc_analysis\scripts")
import chload
import crescent_axis as ca
import tifffile

SP = chload.SP
SPACING = tuple(SP)
VVOL = float(SP.prod())
CYTO = 2.5
SMOOTH_UM, BG_UM, K, VMIN, VMAX, DZ = 0.15, 0.513, 4.48, 0.003, 15.0, 6
DB = np.arange(0.0, CYTO + 1e-6, 0.25)
OUT = r"C:/Users/ryane/coloc_analysis/exposure_sim.csv"
META = pd.read_csv(r"C:/Users/ryane/coloc_analysis/acquisition_metadata.csv")
TARGETS = [90.0]
rng = np.random.RandomState(0)


def pc_dm(syp, gran, outside, dt, bg):
    num = w = 0.0
    for lo, hi in zip(DB[:-1], DB[1:]):
        b = (dt > lo) & (dt <= hi)
        g = gran & b
        o = outside & b
        if g.sum() < 20 or o.sum() < 50:
            continue
        den = float(syp[o].mean()) - bg
        if den <= 0:
            continue
        num += g.sum() * ((float(syp[g].mean()) - bg) / den)
        w += g.sum()
    return round(num / w, 4) if w else None


def trio(syp, gran, cyto, dt):
    bg = float(np.percentile(syp, 3))
    dm = pc_dm(syp, gran, cyto & ~gran, dt, bg)
    gz = np.zeros_like(gran)
    gz[DZ:] = gran[:-DZ]
    gz &= cyto
    zsh = pc_dm(syp, gz, cyto & ~gz, dt, bg)
    spec = round(dm / zsh, 4) if (dm is not None and zsh is not None and zsh > 0) else None
    return dm, zsh, spec


def reexpose(syp, offset, e_from, e_to):
    sig = np.clip(syp - offset, 0, None) * (e_to / e_from)
    shot = rng.poisson(np.clip(sig, 0, None)).astype(np.float32)
    read = rng.normal(0, 1.2, size=syp.shape).astype(np.float32)
    return np.clip(offset + shot + read, 0, 4095)


rows = []
for _, mr in META.iterrows():
    iid = mr.iid
    if chload.is_excluded(iid) or mr.exp_545 not in (200.0, 300.0):
        continue
    dapi, syp_c, lam, lab_c, fb = ca.load_crops(iid)
    rd = ca.find_run(iid)
    csv = os.path.join(rd, f"{iid}__nuclei.csv")
    nc = pd.read_csv(csv)
    germ = [int(x) for x in nc[nc.in_germline.astype(bool)].nucleus_id]
    full = tifffile.imread(os.path.join(rd, f"{iid}__nuclei_labels.tif"))
    gn = np.isin(full, germ)
    obj = ndi.find_objects(gn.astype(np.uint8))[0]
    sl = tuple(slice(max(0, o.start - ca.PAD), min(dim, o.stop + ca.PAD)) for o, dim in zip(obj, gn.shape))
    ch = chload.load(ca.find_nd2(iid, csv))
    pgl = ch["pgl"][sl].astype(np.float32)
    syp = ch["syp"][sl].astype(np.float32)
    del ch
    assert pgl.shape == lam.shape
    env = lam > 0
    dt = ndi.distance_transform_edt(~env, sampling=SPACING)
    cyto = ndi.binary_fill_holes(dt <= CYTO) & ~env
    sm = ndi.gaussian_filter(pgl, sigma=SMOOTH_UM / SP)
    th = sm - ndi.gaussian_filter(sm, sigma=BG_UM / SP)
    rr = th[cyto]
    neg = rr[rr < 0]
    noise = float(np.sqrt(np.mean(neg ** 2))) if neg.size > 50 else float(np.std(rr))
    lab2, _ = ndi.label((th > K * noise) & cyto)
    vol = np.bincount(lab2.ravel()) * VVOL
    keep = np.where((vol >= VMIN) & (vol <= VMAX))[0]
    keep = keep[keep != 0]
    gran = np.isin(lab2, keep)
    far = dt > 10.0
    offset = float(np.percentile(syp[far], 0.5)) if far.sum() > 1e4 else float(np.percentile(syp, 0.5))
    p = chload.parse_iid(iid)
    row = {"iid": iid, "sex": p["sex"], "treat": p["treat"], "exp_545": mr.exp_545, "offset_est": round(offset, 1)}
    row["PC_dm_orig"], row["zsh_orig"], row["PCspec_orig"] = trio(syp, gran, cyto, dt)
    for e_to in TARGETS:
        s2 = reexpose(syp, offset, mr.exp_545, e_to)
        a, b, c = trio(s2, gran, cyto, dt)
        row[f"PC_dm_{int(e_to)}"], row[f"zsh_{int(e_to)}"], row[f"PCspec_{int(e_to)}"] = a, b, c
    rows.append(row)
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"{iid[-18:]:20s} {p['treat']:4s} {mr.exp_545:.0f}ms offset={offset:.0f} | PCspec orig={row['PCspec_orig']} -> @90ms={row['PCspec_90']}"
          f"   (PC_dm {row['PC_dm_orig']} -> {row['PC_dm_90']})", flush=True)
print("DONE", flush=True)
