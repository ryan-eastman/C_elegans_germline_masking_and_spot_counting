"""Auto-build a CV manifest: pair each <results>/<name>/ output with its Imaris .ims (same gonad
name, in --ndir) and the xlsx (in --xdir matching --xlsx-glob) whose nucleus (Surface) count matches
the .ims surface count. Count-matching within the HERM/MALE group avoids ambiguous filename mappings;
cv_analyze further checks the pairing by surface POSITION.

  python scripts/build_cv_manifest.py --ndir "E:/Madeleine/N2/20251105_N2_HS" \
      --xdir "E:/Madeleine/202511_imaris_export_het_mutants_Crest" \
      --results results_cv_hs --xlsx-glob "20251105_N2_hs_*.xlsx" --out cv_manifest_hs.json
"""
import argparse
import glob
import json
import os

import h5py
import pandas as pd

_SURF = "Scene8/Content/MegaSurfaces0/SurfaceModelInfo"


def ims_surf_count(ims):
    with h5py.File(ims, "r") as f:
        return len(f[_SURF]) if _SURF in f else 0


def xlsx_surf_count(xlsx):
    pos = pd.read_excel(xlsx, sheet_name="Position", skiprows=1)
    return int((pos["Category"] == "Surface").sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ndir", required=True, help="folder with the per-gonad .ims files")
    ap.add_argument("--xdir", required=True, help="folder with the Imaris xlsx exports")
    ap.add_argument("--results", required=True, help="pipeline output root (one subdir per gonad)")
    ap.add_argument("--xlsx-glob", required=True, help="glob for the matching xlsx, e.g. 20251105_N2_hs_*.xlsx")
    ap.add_argument("--out", default="cv_manifest.json")
    args = ap.parse_args()

    xlsx_counts = {}
    for f in sorted(glob.glob(os.path.join(args.xdir, args.xlsx_glob))):
        try:
            xlsx_counts[f] = xlsx_surf_count(f)
            print(f"  xlsx {os.path.basename(f)}: {xlsx_counts[f]} surfaces")
        except Exception as e:  # noqa: BLE001
            print(f"  xlsx FAIL {os.path.basename(f)}: {e}")

    manifest = []
    for odir in sorted(glob.glob(os.path.join(args.results, "*"))):
        if not os.path.isdir(odir) or not glob.glob(os.path.join(odir, "*__nuclei.csv")):
            continue
        name = os.path.basename(odir)
        ims = os.path.join(args.ndir, name + ".ims")
        if not os.path.exists(ims):
            print(f"  no .ims for {name}")
            continue
        n = ims_surf_count(ims)
        sex = "MALE" if "MALE" in name.upper() else "HERM"
        cands = {f: c for f, c in xlsx_counts.items() if sex in os.path.basename(f).upper()}
        if not cands:
            print(f"  no {sex} xlsx for {name}")
            continue
        best = min(cands, key=lambda f: abs(cands[f] - n))
        print(f"  {name}: .ims={n} surf -> {os.path.basename(best)} ({cands[best]} surf, "
              f"diff {abs(cands[best]-n)})")
        manifest.append({"name": name, "out_dir": odir, "ims": ims, "xlsx": best})

    with open(args.out, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"\nwrote {args.out} with {len(manifest)} gonads")


if __name__ == "__main__":
    main()
