#!/usr/bin/env python
"""Monday tool: calibrate the transition-zone metric against REAL Imaris TZ annotations.

On synthetic data polarization separates TZ from pachytene (AUC 0.86), but on the real N2 gonads
every DAPI feature is flat along the axis (scripts/diag_zone_real.py) -> DAPI-based TZ calling is
not yet validated on real data. This harness closes that: given a real nucleus LABEL image, its
DAPI volume, and an Imaris-derived per-nucleus ZONE label, it computes every candidate chromatin-
asymmetry metric on the REAL nuclei and reports which (if any) separates real TZ from pachytene,
with the threshold that maximizes Youden's J. If none separate, the answer is "needs a TZ marker".

    # validate the machinery on synthetic (has zone ground truth):
    python scripts/calibrate_zones.py --synth
    # Monday, on real data:
    python scripts/calibrate_zones.py --labels results/<img>__nuclei_labels.tif \
        --nd2 "data/.../<img>.nd2" --zones imaris_zones.csv   # CSV: nucleus_id,zone
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from skimage.measure import regionprops

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from synth_nuclei import DX, DY, DZ, degrade, synthesize_zone_scene  # noqa: E402


def compute_zone_metrics(labels, dna, spacing):
    """Per-nucleus chromatin-asymmetry metrics -> DataFrame[nucleus_id, polar, hemi, rlen, gyr].
    Geometric (no learned features) so it transfers across the sim-to-real gap better than a CNN."""
    sp = np.asarray(spacing, float)
    voxvol = float(np.prod(sp))
    rows = []
    for rp in regionprops(labels, intensity_image=dna):
        inten = rp.image_intensity[rp.image].astype(float)
        if inten.sum() <= 0 or len(rp.coords) < 20:
            continue
        coords = rp.coords.astype(float) * sp
        w = inten / inten.sum()
        geom = coords.mean(0)
        wc = (coords * w[:, None]).sum(0)
        radius = (3 * len(coords) * voxvol / (4 * np.pi)) ** (1 / 3)
        off = wc - geom
        offn = float(np.linalg.norm(off))
        polar = offn / radius if radius > 0 else 0.0
        if offn > 1e-9:
            s = (coords - geom) @ (off / offn)
            near, far = inten[s > 0].sum(), inten[s <= 0].sum()
            hemi = float((near - far) / (near + far + 1e-9))
        else:
            hemi = 0.0
        rad = coords - geom
        units = rad / np.maximum(np.linalg.norm(rad, axis=1, keepdims=True), 1e-9)
        rlen = float(np.linalg.norm((units * w[:, None]).sum(0)))
        c = coords - wc
        cov = (c * w[:, None]).T @ c
        gyr = float(np.sqrt(max(np.trace(cov), 0.0)) / radius) if radius > 0 else 0.0
        rows.append(dict(nucleus_id=int(rp.label), polar=polar, hemi=hemi, rlen=rlen, gyr=gyr))
    return pd.DataFrame(rows)


def _auc(neg, pos):
    if len(neg) == 0 or len(pos) == 0:
        return float("nan")
    r = pd.Series(np.concatenate([neg, pos])).rank().to_numpy()
    return float((r[len(neg):].sum() - len(pos) * (len(pos) + 1) / 2) / (len(neg) * len(pos)))


def _best_threshold(neg, pos):
    """Threshold maximizing Youden's J (sens+spec-1), assuming TZ has the higher value."""
    vals = np.unique(np.concatenate([neg, pos]))
    best_t, best_j = float("nan"), -1.0
    for t in vals:
        sens = np.mean(pos >= t)
        spec = np.mean(neg < t)
        j = sens + spec - 1
        if j > best_j:
            best_j, best_t = j, float(t)
    return best_t, best_j


def report(metrics, zones):
    df = metrics.merge(zones, on="nucleus_id", how="inner")
    tz = df[df.zone.str.contains("trans", case=False, na=False)]
    pa = df[df.zone.str.contains("pachy", case=False, na=False)]
    print(f"matched {len(df)} nuclei: TZ={len(tz)} pachytene={len(pa)}\n")
    print(f"{'metric':<8} {'TZ mean':>9} {'pachy mean':>11} {'AUC':>6} {'thr*':>7} {'Youden_J':>9}")
    for col in ["polar", "hemi", "rlen", "gyr"]:
        a = _auc(pa[col].to_numpy(), tz[col].to_numpy())
        t, j = _best_threshold(pa[col].to_numpy(), tz[col].to_numpy())
        print(f"{col:<8} {tz[col].mean():>9.3f} {pa[col].mean():>11.3f} {a:>6.2f} {t:>7.3f} {j:>9.2f}")
    print("\nPick the metric with the highest AUC/Youden's J as the real TZ caller + its thr*.")
    print("If the best AUC is ~0.5-0.6, DAPI morphology cannot call the TZ here -> needs a TZ marker.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--synth", action="store_true", help="validate the harness on synthetic GT")
    ap.add_argument("--labels", help="nucleus label image (tif)")
    ap.add_argument("--nd2", help="nd2 for the DAPI channel (channel 0)")
    ap.add_argument("--dna", help="DAPI volume tif (alternative to --nd2)")
    ap.add_argument("--zones", help="Imaris zone CSV: columns nucleus_id,zone")
    a = ap.parse_args()

    if a.synth:
        dapi, labels, gt = synthesize_zone_scene(seed=0, tz_frac=0.5)
        dna = degrade(dapi, seed=100)
        metrics = compute_zone_metrics(labels, dna, (DZ, DY, DX))
        report(metrics, gt.rename(columns={"zone": "zone"}))
        return

    import tifffile
    labels = tifffile.imread(a.labels)
    if a.nd2:
        import germquant.io.nd2_reader as R
        st = R.read_stack(a.nd2)
        dna, spacing = st.data[0], st.spacing
    else:
        dna, spacing = tifffile.imread(a.dna), (DZ, DY, DX)
    metrics = compute_zone_metrics(labels, dna, spacing)
    report(metrics, pd.read_csv(a.zones))


if __name__ == "__main__":
    main()
