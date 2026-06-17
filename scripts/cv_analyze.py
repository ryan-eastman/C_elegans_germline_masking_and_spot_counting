"""Cross-gonad validation of the SpotMAX RAD-51 counter vs the Imaris coloc ground truth.

For each gonad (a row in the --manifest JSON: {name, out_dir, ims, xlsx}):
  * read our pipeline output (germline nuclei + spots, image-relative um),
  * read the Imaris .ims (traced nucleus surfaces + spots, image-relative um),
  * read the xlsx coloc per-nucleus RAD-51 counts (the lab's canonical number) + surface positions,
  * bridge frames: xlsx stage coords -> image coords via the .ims ExtMin, match xlsx<->.ims surfaces
    by position to attach a coloc count to each traced surface,
  * match traced surfaces -> our nuclei (greedy, <=tol) and pair (our n_spots, Imaris coloc count).

Reports, per gonad: recall of the traced subset, our mean vs coloc mean, totals; and ACROSS gonads:
the paired per-nucleus agreement (CCC/Pearson/bias) + the gonad-mean regression. This is what tells
us whether the z-merge calibration generalizes (vs being overfit to HERM_001).

Usage: python scripts/cv_analyze.py --manifest cv_manifest.json [--tol 2.5]
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

from germquant.validate.imaris_ims import read_imaris_ims
from germquant.validate.imaris_xlsx import per_nucleus_rad51


def _ccc(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    sx, sy = x.std(), y.std()
    r = np.corrcoef(x, y)[0, 1]
    return float(2 * r * sx * sy / (sx**2 + sy**2 + (x.mean() - y.mean()) ** 2))


def _greedy(src_xyz, dst_xyz, tol):
    """Match each src point to nearest unused dst within tol; return list of (src_i, dst_j)."""
    used = np.zeros(len(dst_xyz), bool)
    out = []
    for i, p in enumerate(src_xyz):
        if len(dst_xyz) == 0:
            break
        d = np.linalg.norm(dst_xyz - p, axis=1)
        d[used] = np.inf
        j = int(d.argmin())
        if d[j] <= tol:
            used[j] = True
            out.append((i, j))
    return out


def _load(out_dir, pat):
    h = glob.glob(os.path.join(out_dir, pat))
    return pd.read_csv(h[0]) if h else None


def analyze_gonad(g, tol):
    nuc = _load(g["out_dir"], "*__nuclei.csv")
    spots = _load(g["out_dir"], "*__spots.csv")
    if nuc is None:
        return None, []
    germ = nuc[nuc["in_germline"] == True] if "in_germline" in nuc.columns else nuc  # noqa: E712
    im_spots, im_surf, geom = read_imaris_ims(g["ims"])
    xsurf, _ = per_nucleus_rad51(g["xlsx"])  # has X/Y/Z (stage) + RAD51 (coloc count)

    # xlsx stage coords -> image coords via this gonad's ExtMin, then match xlsx<->.ims surfaces
    emin = geom["ext_min_xyz"]
    x_img = np.c_[xsurf["X"].to_numpy() - emin[0], xsurf["Y"].to_numpy() - emin[1],
                  xsurf["Z"].to_numpy() - emin[2]]
    surf_xyz = im_surf[["x_um", "y_um", "z_um"]].to_numpy()
    coloc = np.full(len(im_surf), np.nan)
    for xi, sj in _greedy(x_img, surf_xyz, tol):
        coloc[sj] = xsurf["RAD51"].to_numpy()[xi]

    # match traced surfaces -> our nuclei
    our_xyz = germ[["centroid_x_um", "centroid_y_um", "centroid_z_um"]].to_numpy()
    our_n = germ["n_spots"].to_numpy() if "n_spots" in germ.columns else np.full(len(germ), np.nan)
    pairs = _greedy(surf_xyz, our_xyz, tol)
    recall = len(pairs) / max(len(im_surf), 1)

    rows = []  # (gonad, coloc_count, our_count) for surfaces matched to BOTH a coloc value and our nucleus
    for sj, oj in pairs:
        if not np.isnan(coloc[sj]):
            rows.append((g["name"], coloc[sj], our_n[oj]))
    our_matched = np.array([r[2] for r in rows], float)
    coloc_matched = np.array([r[1] for r in rows], float)

    summary = {
        "gonad": g["name"],
        "n_imaris_surf": len(im_surf),
        "n_xlsx_surf": len(xsurf),
        "xlsx_coloc_assigned": int(np.nansum(coloc)),
        "recall_traced": round(recall, 3),
        "n_paired": len(rows),
        "our_mean": round(float(our_matched.mean()), 2) if len(rows) else float("nan"),
        "coloc_mean": round(float(coloc_matched.mean()), 2) if len(rows) else float("nan"),
        "our_total_spots": int(len(spots)) if spots is not None else 0,
        "imaris_total_spots": int(len(im_spots)),
    }
    return summary, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--tol", type=float, default=2.5)
    ap.add_argument("--out", default="cv_results")
    args = ap.parse_args()

    gonads = json.load(open(args.manifest))
    summaries, all_rows = [], []
    for g in gonads:
        try:
            s, rows = analyze_gonad(g, args.tol)
            if s:
                summaries.append(s)
                all_rows += rows
                print(f"  OK {g['name']}: recall={s['recall_traced']} paired={s['n_paired']} "
                      f"ours={s['our_mean']} coloc={s['coloc_mean']}")
        except Exception as e:  # noqa: BLE001
            print(f"  FAIL {g['name']}: {type(e).__name__}: {e}")

    if not summaries:
        print("no gonads analyzed")
        return
    sm = pd.DataFrame(summaries)
    os.makedirs(args.out, exist_ok=True)
    sm.to_csv(os.path.join(args.out, "cv_per_gonad.csv"), index=False)
    print("\n=== per-gonad ===")
    print(sm.to_string(index=False))

    rows = pd.DataFrame(all_rows, columns=["gonad", "coloc", "ours"])
    rows.to_csv(os.path.join(args.out, "cv_paired_nuclei.csv"), index=False)
    o, c = rows["ours"].to_numpy(float), rows["coloc"].to_numpy(float)
    print("\n=== PAIRED per-nucleus across all gonads (ours vs Imaris coloc) ===")
    print(f"  n={len(rows)}  ours mean={o.mean():.2f}  coloc mean={c.mean():.2f}  "
          f"bias(ours-coloc)={o.mean()-c.mean():+.2f}")
    print(f"  Pearson r={np.corrcoef(o, c)[0,1]:.3f}  CCC={_ccc(o, c):.3f}")
    print("\n=== gonad-mean regression (our_mean vs coloc_mean) ===")
    gm_o, gm_c = sm["our_mean"].to_numpy(float), sm["coloc_mean"].to_numpy(float)
    ok = ~np.isnan(gm_o) & ~np.isnan(gm_c)
    if ok.sum() >= 2:
        slope, intc = np.polyfit(gm_c[ok], gm_o[ok], 1)
        print(f"  n_gonads={ok.sum()}  Pearson r={np.corrcoef(gm_c[ok], gm_o[ok])[0,1]:.3f}  "
              f"CCC={_ccc(gm_o[ok], gm_c[ok]):.3f}  fit ours={slope:.2f}*coloc+{intc:.2f}")
    print(f"\noutputs -> {args.out}")


if __name__ == "__main__":
    main()
