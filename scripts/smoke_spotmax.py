"""Smoke test: run SpotMAX headless spot detection on a Cahoon RAD-51 crop inside our nucleus masks,
and report per-nucleus spot counts. Confirms (a) SpotMAX runs on our data, (b) cellpose-4.2 masks
load, (c) the pipe API usage. Prints each step's return type so we learn the exact interface.
"""
import numpy as np
import pandas as pd
import tifffile

import germquant.io.nd2_reader as R
import spotmax.pipe as P

# data lives in the (frozen) original repo; the fork just reads it
ND2 = r"C:\Users\ryane\C_elegans_ml\data\raw_examples\madeleline images\20251105_N2_noHS\20251105_N2_nohs_HERM _001.nd2"
LAB = r"C:\Users\ryane\C_elegans_ml\results_germline\20251105_N2_nohs_HERM _001__nuclei_labels.tif"
SP = (0.2, 0.108333, 0.108333)
CY, CX, H = 1761, 1423, 250


def step(name, fn):
    try:
        out = fn()
        print(f"  [OK] {name} -> {type(out).__name__}", getattr(out, "shape", ""))
        return out
    except Exception as e:
        print(f"  [FAIL] {name}: {type(e).__name__}: {e}")
        return None


def main():
    import os
    cache = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "smoke_crop.npz")
    if os.path.exists(cache):
        d = np.load(cache)
        rad, lab = d["rad"], d["lab"]
    else:
        st = R.read_stack(ND2)
        rad = np.asarray(st.data[2][:, CY - H:CY + H, CX - H:CX + H]).astype("float32")
        lab = tifffile.imread(LAB)[:, CY - H:CY + H, CX - H:CX + H].astype("int32")
        np.savez(cache, rad=rad, lab=lab)
    radii_px = np.array([0.3 / SP[0], 0.3 / SP[1], 0.3 / SP[2]])
    print(f"crop {rad.shape}  nuclei in crop = {len(np.unique(lab)) - 1}  spot radii px = {radii_px.round(2)}")

    pre = step("preprocess_image", lambda: P.preprocess_image(
        rad, lab=lab, gauss_sigma=0.75, use_gpu=False, spots_zyx_radii_pxl=radii_px))
    pre = pre if pre is not None else rad

    segm = step("spots_semantic_segmentation", lambda: P.spots_semantic_segmentation(
        pre, lab=lab, spots_zyx_radii_pxl=radii_px, do_try_all_thresholds=False,
        return_only_segm=True, thresholding_method="threshold_otsu", use_gpu=False))

    res = step("spot_detection", lambda: P.spot_detection(
        pre, spots_segmantic_segm=segm, spots_zyx_radii_pxl=radii_px, lab=lab, return_df=True))
    if res is None:
        return
    df = res[0] if isinstance(res, tuple) else res
    print("  spot_detection columns:", list(df.columns) if hasattr(df, "columns") else type(df))
    print("  total spots detected:", len(df))

    # assign each spot to a nucleus via the label at its (z,y,x) and count per nucleus
    coords = [c for c in df.columns if c in ("z", "y", "x")] if hasattr(df, "columns") else []
    if {"z", "y", "x"}.issubset(set(coords)):
        zi = np.clip(df["z"].astype(int), 0, lab.shape[0] - 1)
        yi = np.clip(df["y"].astype(int), 0, lab.shape[1] - 1)
        xi = np.clip(df["x"].astype(int), 0, lab.shape[2] - 1)
        nid = lab[zi, yi, xi]
        inside = nid[nid > 0]
        per_nuc = pd.Series(inside).value_counts()
        n_nuc = len(np.unique(lab)) - 1
        print(f"\n  spots inside nuclei: {len(inside)}/{len(df)}")
        print(f"  per-nucleus RAD-51 count: mean={inside.size / max(n_nuc,1):.2f}  "
              f"(N2 control biological expectation ~0.5-2/nucleus)")
        print(f"  nuclei with >=1 spot: {per_nuc.size}/{n_nuc}")
    else:
        print("  (coords columns not z/y/x; df head:)\n", df.head())


if __name__ == "__main__":
    main()
