#!/usr/bin/env python
"""CAHOON-SPECIFIC quantitative validation WITHOUT new ground truth, using the biological prior:
RAD-51 foci mark DSBs on meiotic chromosomes -> they live INSIDE pachytene nuclei. The original
root-cause was that stock cpsam's incomplete masks ORPHAN 51-60% of foci (chromosomes left out of
the mask) and crop the SC trace. Better nucleus masks should RECAPTURE those foci.

Identical foci for every model (detected once on the RAD-51 channel); only the mask assignment
differs -> a clean comparison. Mask VOLUME is reported as the confound control: a model can't win
by simply ballooning masks over empty space (the rachis), because there are no foci there.

    python scripts/foci_recapture.py
"""
import numpy as np

import germquant.io.nd2_reader as R
from germquant.foci.detect import detect_foci

REAL = r"data/raw_examples/madeleline images/20251105_N2_noHS/20251105_N2_nohs_HERM _001.nd2"
# representative pachytene window (full z), centered where prior crops showed pachytene nuclei
CY, CX, HY, HX = 1761, 1423, 350, 350
SP = (0.2, 0.108333333333333, 0.108333333333333)

MODELS = [
    ("STOCK cpsam", None),
    ("synthetic-only", "models/models/germline_nuclei"),
    ("real-Koehler", "models/models/germline_nuclei_real"),
]


def seg3d(model_path, dna):
    """3D Cellpose at NATIVE scale (no diameter rescale) -- the regime our fine-tuned models were
    benchmarked in (F1 0.91). anisotropy from real voxel spacing; do_3D recombines xy/xz/yz flows."""
    from cellpose import models
    aniso = SP[0] / SP[1]
    model = models.CellposeModel(gpu=True, pretrained_model=model_path) if model_path \
        else models.CellposeModel(gpu=True)
    out = model.eval(dna, do_3D=True, z_axis=0, channel_axis=None, anisotropy=aniso)
    return np.asarray(out[0]).astype(np.int32)


def main():
    st = R.read_stack(REAL)
    dna = st.data[0][:, CY - HY:CY + HY, CX - HX:CX + HX].astype(np.float32)
    rad = st.data[2][:, CY - HY:CY + HY, CX - HX:CX + HX].astype(np.float32)
    vox_um3 = SP[0] * SP[1] * SP[2]
    print(f"crop {dna.shape}  ({np.prod(dna.shape)/1e6:.1f} Mvox)")

    # detect foci ONCE (frozen params) -> identical set for every model
    foci = detect_foci(rad, SP, threshold_rel=0.10, min_sigma_um=0.10, max_sigma_um=0.35)
    zi = np.clip(np.round(foci.z_um.values / SP[0]).astype(int), 0, dna.shape[0] - 1)
    yi = np.clip(np.round(foci.y_um.values / SP[1]).astype(int), 0, dna.shape[1] - 1)
    xi = np.clip(np.round(foci.x_um.values / SP[2]).astype(int), 0, dna.shape[2] - 1)
    nfoci = len(foci)
    print(f"RAD-51 foci detected (same for all models): {nfoci}\n")

    print(f"{'model':<16} {'nuclei':>6} {'maskVol_um3':>12} {'%foci_inside':>13} {'foci_captured':>14}")
    print("-" * 66)
    for name, path in MODELS:
        lab = seg3d(path, dna)
        inside = lab[zi, yi, xi] > 0
        vol = float((lab > 0).sum()) * vox_um3
        cap = int(inside.sum())
        print(f"{name:<16} {lab.max():>6d} {vol:>12.0f} {100*cap/max(nfoci,1):>12.1f}% {cap:>8d}/{nfoci}")


if __name__ == "__main__":
    main()
