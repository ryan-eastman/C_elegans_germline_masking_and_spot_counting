"""Run the SpotMAX pipeline on every per-gonad .nd2 in a folder (HERM/MALE files, excluding the
multi-position parent + largeimage/10x scans), into <out>/<name>/, for cross-gonad validation.
Sequential (GPU-serial), ~20 min each.

  python scripts/cv_run.py --ndir "E:/Madeleine/N2/20251105_N2_HS" --out results_cv_hs
"""
import argparse
import glob
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


def gonad_files(ndir):
    return sorted(
        f for f in glob.glob(os.path.join(ndir, "*.nd2"))
        if ("HERM" in os.path.basename(f).upper() or "MALE" in os.path.basename(f).upper())
        and "largeimage" not in os.path.basename(f).lower()
        and "10x" not in os.path.basename(f).lower()
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ndir", required=True, help="folder of per-gonad .nd2 files")
    ap.add_argument("--out", required=True, help="output root (one subdir per gonad)")
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()

    files = gonad_files(args.ndir)
    if not files:
        print(f"no per-gonad .nd2 found under {args.ndir}", file=sys.stderr)
        return 1
    cfg = load_config(args.config)
    out = ROOT / args.out
    prov = provenance.write_manifest(out, config_hash=cfg.hash, config=cfg.as_dict())
    print(f"{len(files)} gonads -> {out}")
    for i, nd2 in enumerate(files, 1):
        name = Path(nd2).stem.strip()
        t0 = time.time()
        print(f"\n=== [{i}/{len(files)}] {os.path.basename(nd2)}", flush=True)
        try:
            res = process_image(nd2, cfg, out / name, prov=prov)
            print(f"=== done {name} in {time.time()-t0:.0f}s: {res['n_nuclei']} nuclei, "
                  f"spots={len(res['tables']['spots'])}, qc={res['qc_pass']}", flush=True)
        except Exception as e:  # noqa: BLE001
            import traceback
            traceback.print_exc()
            print(f"=== FAIL {name}: {type(e).__name__}: {e}", flush=True)
    print("\n=== CV RUN COMPLETE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
