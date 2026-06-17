"""End-to-end test of the new SpotMAX pipeline on a real Cahoon N2 control .nd2.

Runs the full pipeline: read -> segment nuclei (trained cellpose) -> isolate germline ->
SpotMAX spots -> tidy tables + montage + label mask. chdir to repo root first so the relative
cellpose model path in config.yaml resolves. Prints a summary of every table written.
"""
import logging
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
os.chdir(ROOT)  # so config's relative model path (models/models/...) resolves
sys.path.insert(0, str(ROOT / "src"))

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s", force=True)

from germquant import provenance  # noqa: E402
from germquant.config import load_config  # noqa: E402
from germquant.pipeline import process_image  # noqa: E402

ND2 = r"C:\Users\ryane\C_elegans_ml\data\raw_examples\madeleline images\20251105_N2_noHS\20251105_N2_nohs_HERM _001.nd2"
OUT = ROOT / "results_spotmax_test"

t0 = time.time()
print(f"=== run_real: {ND2}", flush=True)
cfg = load_config("config/config.yaml")
prov = provenance.write_manifest(OUT, config_hash=cfg.hash, config=cfg.as_dict())
res = process_image(ND2, cfg, OUT, prov=prov)

print(f"\n=== DONE in {time.time() - t0:.0f}s", flush=True)
print(f"image_id   : {res['image_id']}")
print(f"n_nuclei   : {res['n_nuclei']}")
print(f"qc_pass    : {res['qc_pass']}")
print(f"flags      : {'; '.join(res['qc_flags'])}")
for name, df in res["tables"].items():
    print(f"  table {name:14s}: {len(df):5d} rows")

nuc = res["tables"]["nuclei"]
if "n_spots" in nuc.columns:
    germ = nuc[nuc["in_germline"] == True] if "in_germline" in nuc.columns else nuc  # noqa: E712
    ns = germ["n_spots"].dropna()
    if len(ns):
        print(f"\nspots/nucleus (germline): mean={ns.mean():.2f}  "
              f"median={ns.median():.0f}  max={ns.max():.0f}  "
              f"nuclei>=1 spot: {(ns > 0).sum()}/{len(ns)}")
print(f"\noutputs -> {OUT}")
