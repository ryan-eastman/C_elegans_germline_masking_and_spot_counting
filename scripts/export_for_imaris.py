#!/usr/bin/env python
"""Export the pipeline's nucleus label image as an Imaris-ready OME-TIFF: uint16 with the voxel size
baked into the OME metadata, so Imaris loads it with the CORRECT calibration (no manual voxel entry,
no z-squash) and you can build Surfaces straight from it to QC the segmentation against your DAPI.

    python scripts/export_for_imaris.py --results results_germline --image 20251105 \
        --nd2 "data/raw_examples/madeleline images/20251105_N2_noHS/20251105_N2_nohs_HERM _001.nd2"
    # or give spacing explicitly:
    python scripts/export_for_imaris.py --labels <labels>.tif --dz 0.2 --dy 0.108333 --dx 0.108333
"""
import argparse
import glob

import numpy as np
import tifffile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results_germline")
    ap.add_argument("--image", default="")
    ap.add_argument("--labels", help="explicit label tif (overrides --results/--image)")
    ap.add_argument("--nd2", help="read voxel size from this .nd2")
    ap.add_argument("--dz", type=float)
    ap.add_argument("--dy", type=float)
    ap.add_argument("--dx", type=float)
    a = ap.parse_args()

    lab_path = a.labels or glob.glob(f"{a.results}/{a.image}*__nuclei_labels.tif")[0]
    labels = tifffile.imread(lab_path)

    if a.nd2:
        import germquant.io.nd2_reader as R
        dz, dy, dx = R.read_nd2_metadata(a.nd2)["spacing"]
    elif a.dz:
        dz, dy, dx = a.dz, a.dy, a.dx
    else:
        dz, dy, dx = 0.2, 0.108333, 0.108333          # the Cahoon default
        print("WARNING: no --nd2/--spacing given; assuming Cahoon voxel (0.2, 0.108, 0.108) um")

    out = lab_path.rsplit(".tif", 1)[0] + "_imaris.ome.tif"
    dtype = np.uint16 if int(labels.max()) <= 65535 else np.int32
    tifffile.imwrite(
        out, labels.astype(dtype), ome=True,
        metadata={
            "axes": "ZYX",
            "PhysicalSizeZ": float(dz), "PhysicalSizeZUnit": "µm",
            "PhysicalSizeY": float(dy), "PhysicalSizeYUnit": "µm",
            "PhysicalSizeX": float(dx), "PhysicalSizeXUnit": "µm",
        },
    )
    print(f"wrote {out}")
    print(f"  shape {labels.shape}  dtype {np.dtype(dtype).name}  n_nuclei {int(labels.max())}  "
          f"voxel (dz,dy,dx)=({dz:.4g},{dy:.4g},{dx:.4g}) um")
    print("  In Imaris: open this + your DAPI; build Surfaces -> 'from labels'/label image to QC.")


if __name__ == "__main__":
    main()
