"""Turn a pipeline `*__spots.csv` (a coordinate list) into a `*__spots.tif` — a 3D image with a small
bright blob at each detected spot, on the SAME voxel grid as the nucleus-label TIF / original image.
Load it into Imaris as an extra channel to view the spots, or run Imaris' Spots detection on it
(each blob -> one Spot) to get Spots objects alongside your nucleus Surfaces.

  python scripts/spots_to_tif.py "PATH\TO\results\GONAD_DIR"        # finds *__spots.csv + *__nuclei_labels.tif
  python scripts/spots_to_tif.py "PATH\TO\X__spots.csv"             # explicit file
  python scripts/spots_to_tif.py "DIR" --radius-um 0.3 --label      # --label = unique id per spot
"""
import argparse
import glob
import os

import pandas as pd
import tifffile

from germquant.spots import spots_to_image


def _resolve(path):
    if os.path.isdir(path):
        hits = glob.glob(os.path.join(path, "*__spots.csv"))
        if not hits:
            raise SystemExit(f"no *__spots.csv in {path}")
        return hits[0]
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="a gonad results dir, or a *__spots.csv file")
    ap.add_argument("--radius-um", type=float, default=0.3, help="blob radius (default = spot radius)")
    ap.add_argument("--label", action="store_true", help="write a unique integer id per spot (else all bright)")
    args = ap.parse_args()

    spots_csv = _resolve(args.path)
    base = spots_csv[:-len("__spots.csv")]
    labels_tif = base + "__nuclei_labels.tif"
    if not os.path.exists(labels_tif):
        raise SystemExit(f"need the label TIF for image dimensions: {labels_tif}")

    df = pd.read_csv(spots_csv)
    shape = tifffile.imread(labels_tif).shape          # (Z, Y, X), same grid as the original image
    sp = (float(df["voxel_dz_um"].iloc[0]), float(df["voxel_dy_um"].iloc[0]), float(df["voxel_dx_um"].iloc[0])) \
        if "voxel_dz_um" in df.columns and len(df) else (0.2, 0.108333, 0.108333)

    out = spots_to_image(df, shape, sp, radius_um=args.radius_um, label=args.label)

    out_tif = base + "__spots.tif"
    # ImageJ-style metadata so Imaris reads the voxel size (z spacing matters for 3D)
    tifffile.imwrite(out_tif, out, compression="zlib", imagej=True,
                     resolution=(1 / sp[2], 1 / sp[1]),
                     metadata={"spacing": sp[0], "unit": "um", "axes": "ZYX"})
    print(f"wrote {out_tif}")
    print(f"  {len(df)} spots, image {shape}, voxel (z,y,x)={sp} um")
    print("  In Imaris: add it as a Channel (same dims as your image), then either view the blobs or")
    print("  run Spots detection on this channel (set spot diameter ~0.4 um) to get Spots objects.")


if __name__ == "__main__":
    main()
