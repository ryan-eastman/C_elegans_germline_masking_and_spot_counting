"""Run the SpotMAX pipeline on the 6 N2 no-HS gonads (20251105) that have Imaris GT, into
results_cv/<name>/, for cross-gonad validation. Sequential (GPU-serial), ~20 min each."""
import logging
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
os.chdir(ROOT)
sys.path.insert(0, str(ROOT / "src"))
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s", force=True)

from germquant import provenance  # noqa: E402
from germquant.config import load_config  # noqa: E402
from germquant.pipeline import process_image  # noqa: E402

NDIR = r"E:\Madeleine\N2\20251105_N2_noHS"
GONADS = [
    "20251105_N2_nohs_HERM _.nd2", "20251105_N2_nohs_HERM _001.nd2", "20251105_N2_nohs_HERM _002.nd2",
    "20251105_N2_nohs_MALE _.nd2", "20251105_N2_nohs_MALE _001.nd2", "20251105_N2_nohs_MALE _002.nd2",
]
OUT = ROOT / "results_cv"

cfg = load_config("config/config.yaml")
prov = provenance.write_manifest(OUT, config_hash=cfg.hash, config=cfg.as_dict())
for i, fn in enumerate(GONADS, 1):
    nd2 = os.path.join(NDIR, fn)
    name = Path(fn).stem.strip()
    odir = OUT / name
    t0 = time.time()
    print(f"\n=== [{i}/{len(GONADS)}] {fn}", flush=True)
    try:
        res = process_image(nd2, cfg, odir, prov=prov)
        print(f"=== done {name} in {time.time()-t0:.0f}s: {res['n_nuclei']} nuclei, "
              f"spots={len(res['tables']['spots'])}, qc={res['qc_pass']}", flush=True)
    except Exception as e:  # noqa: BLE001
        import traceback
        traceback.print_exc()
        print(f"=== FAIL {name}: {type(e).__name__}: {e}", flush=True)
print("\n=== CV RUN COMPLETE", flush=True)
