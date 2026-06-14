#!/usr/bin/env python
"""Hardened foci-recapture (post adversarial-review). Adds the controls the review demanded:
  * voxel spacing READ from the .nd2 (st.spacing), not hardcoded
  * ENRICHMENT = %foci_inside / %volume_fill  (targeting per unit volume -- so '100%' isn't just
    bought by bigger masks)
  * NULL DILATION CONTROL: dilate stock & synthetic masks to real's fill fraction, re-measure
    capture -> isolates how much of the gap is volume vs mask placement
  * MERGE PROXY: median foci per occupied nucleus (implausibly high => adjacent nuclei merged,
    which %inside alone is blind to)
  * MULTIPLE crops (not n=1)
Biological prior only (RAD-51 foci inside nuclei); frozen foci params; no Cahoon GT used.

    python scripts/foci_recapture2.py
"""
import numpy as np
from scipy import ndimage as ndi

import germquant.io.nd2_reader as R
from germquant.foci.detect import detect_foci

REAL = r"data/raw_examples/madeleline images/20251105_N2_noHS/20251105_N2_nohs_HERM _001.nd2"
CROPS = [(872, 840), (1464, 1368), (2064, 1736)]   # on-tissue centers (scripts/find_crops.py)
HY = HX = 300
MODELS = [
    ("STOCK cpsam", None),
    ("synthetic-only", "models/models/germline_nuclei"),
    ("real (leak-free)", "models/models/germline_nuclei_real_grp"),
    ("real+synth combined", "models/models/germline_nuclei_combined"),
]


def seg3d(model, dna, aniso):
    out = model.eval(dna, do_3D=True, z_axis=0, channel_axis=None, anisotropy=aniso)
    return np.asarray(out[0]).astype(np.int32)


def dilate_to_fill(fg, target_fill):
    """Dilate a binary foreground (3D ball-ish) until its fill fraction >= target."""
    cur = fg.copy()
    for _ in range(25):
        if cur.mean() >= target_fill:
            break
        cur = ndi.binary_dilation(cur, iterations=1)
    return cur


def main():
    import os
    from cellpose import models as cpm
    st = R.read_stack(REAL)
    SP = st.spacing
    aniso = SP[0] / SP[1]
    vox = SP[0] * SP[1] * SP[2]
    print(f"voxel spacing READ from nd2: {tuple(round(s,4) for s in SP)}  spacing_ok={st.spacing_ok}")
    avail = [(n, p) for n, p in MODELS if p is None or os.path.exists(p)]
    loaded = {n: (cpm.CellposeModel(gpu=True, pretrained_model=p) if p else cpm.CellposeModel(gpu=True))
              for n, p in avail}

    # accumulate per-model across crops
    agg = {n: dict(inside=0, n=0, fillv=0.0, cropv=0.0, vol=0.0, obj=0, fpn=[]) for n, _ in avail}
    real_fill_by_crop = []

    for ci, (cy, cx) in enumerate(CROPS):
        dna = st.data[0][:, cy - HY:cy + HY, cx - HX:cx + HX].astype(np.float32)
        rad = st.data[2][:, cy - HY:cy + HY, cx - HX:cx + HX].astype(np.float32)
        cropv = float(np.prod(dna.shape)) * vox
        foci = detect_foci(rad, SP, threshold_rel=0.10, min_sigma_um=0.10, max_sigma_um=0.35)
        zi = np.clip(np.round(foci.z_um.values / SP[0]).astype(int), 0, dna.shape[0] - 1)
        yi = np.clip(np.round(foci.y_um.values / SP[1]).astype(int), 0, dna.shape[1] - 1)
        xi = np.clip(np.round(foci.x_um.values / SP[2]).astype(int), 0, dna.shape[2] - 1)
        nf = len(foci)
        labs = {n: seg3d(m, dna, aniso) for n, m in loaded.items()}
        rf = float((labs.get("real (leak-free)", list(labs.values())[-1]) > 0).mean())
        real_fill_by_crop.append(rf)
        print(f"\ncrop {ci} (cy{cy},cx{cx})  foci={nf}")
        for n, lab in labs.items():
            ins = int((lab[zi, yi, xi] > 0).sum())
            fill = float((lab > 0).mean())
            occ = lab[zi, yi, xi]; occ = occ[occ > 0]
            fpn = np.bincount(occ).max() if occ.size else 0   # max foci in one nucleus (merge proxy hint)
            # per-occupied-nucleus foci median
            if occ.size:
                vals = np.bincount(occ); vals = vals[vals > 0]
                medfpn = float(np.median(vals))
            else:
                medfpn = 0.0
            a = agg[n]
            a["inside"] += ins; a["n"] += nf; a["fillv"] += fill * cropv
            a["cropv"] += cropv; a["vol"] += float((lab > 0).sum()) * vox; a["obj"] += int(lab.max())
            a["fpn"].append(medfpn)
            print(f"   {n:<20} inside={100*ins/max(nf,1):5.1f}% fill={100*fill:4.1f}% "
                  f"obj={lab.max():4d} um3/obj={float((lab>0).sum())*vox/max(lab.max(),1):5.1f} medFoci/nuc={medfpn:.1f}")

    target_fill = float(np.mean(real_fill_by_crop))
    print("\n" + "=" * 84)
    print(f"SUMMARY across {len(CROPS)} crops   (enrichment = %inside / %fill; higher = better-targeted/volume)")
    print(f"{'model':<20} {'%inside':>8} {'%fill':>7} {'enrich':>7} {'um3/obj':>8} {'medFoci/nuc':>12}")
    print("-" * 84)
    for n, _ in avail:
        a = agg[n]
        pin = 100 * a["inside"] / max(a["n"], 1)
        pfill = 100 * a["fillv"] / max(a["cropv"], 1)
        enr = pin / max(pfill, 1e-9)
        umobj = a["vol"] / max(a["obj"], 1)
        print(f"{n:<20} {pin:7.1f}% {pfill:6.1f}% {enr:6.1f}x {umobj:7.1f} {np.mean(a['fpn']):11.1f}")

    # NULL CONTROL on crop 0: dilate stock & synthetic to real's fill, re-measure capture
    print("\n" + "=" * 84)
    print(f"NULL DILATION CONTROL (crop 0): dilate masks to real's fill ~{100*target_fill:.1f}%, re-measure capture")
    cy, cx = CROPS[0]
    dna = st.data[0][:, cy - HY:cy + HY, cx - HX:cx + HX].astype(np.float32)
    rad = st.data[2][:, cy - HY:cy + HY, cx - HX:cx + HX].astype(np.float32)
    foci = detect_foci(rad, SP, threshold_rel=0.10, min_sigma_um=0.10, max_sigma_um=0.35)
    zi = np.clip(np.round(foci.z_um.values / SP[0]).astype(int), 0, dna.shape[0] - 1)
    yi = np.clip(np.round(foci.y_um.values / SP[1]).astype(int), 0, dna.shape[1] - 1)
    xi = np.clip(np.round(foci.x_um.values / SP[2]).astype(int), 0, dna.shape[2] - 1)
    nf = len(foci)
    for n in ("STOCK cpsam", "synthetic-only"):
        if n not in loaded:
            continue
        lab = seg3d(loaded[n], dna, aniso)
        fg = lab > 0
        before = 100 * int((lab[zi, yi, xi] > 0).sum()) / max(nf, 1)
        dil = dilate_to_fill(fg, target_fill)
        after = 100 * int(dil[zi, yi, xi].sum()) / max(nf, 1)
        print(f"  {n:<16} capture {before:5.1f}% -> dilated-to-{100*dil.mean():.0f}%fill {after:5.1f}%  "
              f"(real (leak-free) reaches ~100% at {100*target_fill:.0f}% fill)")


if __name__ == "__main__":
    main()
