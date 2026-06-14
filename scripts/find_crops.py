"""Locate on-germline crop centers: downsample the DAPI + RAD-51 volumes, find the brightest
tissue regions, and propose 3 well-separated 600x600 crop centers that are actually on the gonad
(so the multi-crop foci-recapture isn't placed off-tissue)."""
import numpy as np
from scipy import ndimage as ndi

import germquant.io.nd2_reader as R

REAL = r"data/raw_examples/madeleline images/20251105_N2_noHS/20251105_N2_nohs_HERM _001.nd2"
HY = HX = 300


def main():
    st = R.read_stack(REAL, xy_stride=8)        # cheap 325x325 view
    s = st.downsample_xy
    dapi = st.data[0].max(0).astype(np.float32)  # z-MIP DAPI
    rad = st.data[2].max(0).astype(np.float32)
    # tissue score = smoothed DAPI brightness
    score = ndi.gaussian_filter(dapi, 3)
    H, W = score.shape
    half = HY // s
    centers = []
    sc = score.copy()
    for _ in range(3):
        # mask out borders so the full-res crop stays in bounds
        sc[:half, :] = sc[-half:, :] = sc[:, :half] = sc[:, -half:] = 0
        yx = np.unravel_index(np.argmax(sc), sc.shape)
        cy, cx = int(yx[0] * s), int(yx[1] * s)
        centers.append((cy, cx))
        # suppress a neighborhood so the next center is well-separated
        y0, y1 = max(0, yx[0] - 2 * half), min(H, yx[0] + 2 * half)
        x0, x1 = max(0, yx[1] - 2 * half), min(W, yx[1] + 2 * half)
        sc[y0:y1, x0:x1] = 0
    print("proposed on-tissue crop centers (cy,cx) at full res:")
    for cy, cx in centers:
        dmean = float(dapi[max(0, cy // s - half):cy // s + half, max(0, cx // s - half):cx // s + half].mean())
        rmean = float(rad[max(0, cy // s - half):cy // s + half, max(0, cx // s - half):cx // s + half].mean())
        print(f"  ({cy}, {cx})  DAPI_mean={dmean:.0f}  RAD_mean={rmean:.0f}")


if __name__ == "__main__":
    main()
