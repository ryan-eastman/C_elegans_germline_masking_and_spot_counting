"""Auto-build cv_manifest.json: pair each results_cv/<name>/ output with its Imaris .ims (same gonad
name) and the xlsx whose nucleus (Surface) count matches the .ims surface count. Matching by count
(within the HERM/MALE group) avoids the ambiguous '_' / '_001' / '_002' vs '_1/2/3' filename mapping;
the downstream cv_analyze further checks the pairing by surface POSITION."""
import glob
import json
import os
from pathlib import Path

import h5py
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
NDIR = r"E:\Madeleine\N2\20251105_N2_noHS"
XDIR = r"E:\Madeleine\202511_imaris_export_het_mutants_Crest"
OUT = ROOT / "results_cv"
_SURF = "Scene8/Content/MegaSurfaces0/SurfaceModelInfo"


def ims_surf_count(ims):
    with h5py.File(ims, "r") as f:
        return len(f[_SURF]) if _SURF in f else 0


def xlsx_surf_count(xlsx):
    pos = pd.read_excel(xlsx, sheet_name="Position", skiprows=1)
    return int((pos["Category"] == "Surface").sum())


# candidate xlsx counts, split HERM/MALE
xlsx_counts = {}
for f in sorted(glob.glob(os.path.join(XDIR, "20251105_N2_nohs_*.xlsx"))):
    try:
        xlsx_counts[f] = xlsx_surf_count(f)
        print(f"  xlsx {os.path.basename(f)}: {xlsx_counts[f]} surfaces")
    except Exception as e:  # noqa: BLE001
        print(f"  xlsx FAIL {os.path.basename(f)}: {e}")

manifest = []
for odir in sorted(glob.glob(str(OUT / "*"))):
    if not os.path.isdir(odir) or not glob.glob(os.path.join(odir, "*__nuclei.csv")):
        continue
    name = os.path.basename(odir)
    ims = os.path.join(NDIR, name + ".ims")
    if not os.path.exists(ims):
        print(f"  no .ims for {name}")
        continue
    n = ims_surf_count(ims)
    sex = "MALE" if "MALE" in name.upper() else "HERM"
    # nearest-count xlsx within the same sex group
    cands = {f: c for f, c in xlsx_counts.items() if sex in os.path.basename(f).upper()}
    if not cands:
        print(f"  no {sex} xlsx for {name}")
        continue
    best = min(cands, key=lambda f: abs(cands[f] - n))
    print(f"  {name}: .ims={n} surf -> {os.path.basename(best)} ({cands[best]} surf, "
          f"diff {abs(cands[best]-n)})")
    manifest.append({"name": name, "out_dir": odir, "ims": ims, "xlsx": best})

with open(ROOT / "cv_manifest.json", "w") as fh:
    json.dump(manifest, fh, indent=2)
print(f"\nwrote cv_manifest.json with {len(manifest)} gonads")
