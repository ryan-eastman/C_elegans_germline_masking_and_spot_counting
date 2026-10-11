#!/usr/bin/env python
"""Which Cellpose base was germline_nuclei_combined fine-tuned from? (docs/DESIGN_system_packs.md, 3.11 and 15d)

    python scripts/compare_base_weights.py [--model PATH] [--cpsam-v2 PATH] [--cpsam PATH] [--write]

Compares the fine-tuned state dict with both cached bases tensor by tensor (germquant.packs.basecheck) and
prints the record. Defaults: the release model where `germquant get-model` puts it, and cpsam_v2 and cpsam in
Cellpose's cache (CELLPOSE_LOCAL_MODELS_PATH or ~/.cellpose/models). CPU only, no GPU: it reads three 1.2 GB
files one base at a time (about 2.5 GB of memory at the peak) and hashes them.

--write puts the result into packs/celegans_germline/models/germline_nuclei_combined.card.json
(base_model.comparison, and base_model.name and sha256 when the comparison names a base; otherwise the card keeps
"not determined", as the design requires) and rebuilds the pack's SHA256SUMS. The pack's version stays 0.1.0
while it is unreleased; after a release, a card change is a PATCH bump.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from germquant.datafetch import sha256_file  # noqa: E402
from germquant.models_registry import MODELS, resolve_registered  # noqa: E402
from germquant.packs import write_sums  # noqa: E402
from germquant.packs.basecheck import compare_files  # noqa: E402
from germquant.packs.blobs import cellpose_cache_dirs  # noqa: E402

PACK = REPO / "packs" / "celegans_germline"
CARD = PACK / "models" / "germline_nuclei_combined.card.json"


def _finite(v):
    """JSON-safe: a non-finite float becomes its text ("inf"), containers are walked."""
    if isinstance(v, float) and not math.isfinite(v):
        return str(v)
    if isinstance(v, dict):
        return {k: _finite(x) for k, x in v.items()}
    if isinstance(v, list):
        return [_finite(x) for x in v]
    return v


HOW = ("python scripts/compare_base_weights.py --write (CPU only: the three files are read with torch.load on the "
       "CPU; it writes this record and rebuilds the pack's SHA256SUMS)")


def update_card(card_path: Path, rec: dict, base_paths: dict[str, Path]) -> dict:
    """Write the comparison record into the card (and the base's name and digest when it is determined)."""
    card = json.loads(card_path.read_text(encoding="utf-8"))
    base = card.setdefault("base_model", {})
    comp = _finite({"method": rec["method"], "result": rec["result"], "separation_ratio": rec["separation_ratio"],
                    "separation_required": rec["separation_required"], "min_coverage": rec["min_coverage"],
                    "bases": rec["bases"], "how": HOW})
    base["comparison"] = comp
    if rec["result"] != "not determined":
        base["name"] = rec["result"]
        base["sha256"] = sha256_file(base_paths[rec["result"]])
    else:
        base["name"] = "not determined"
    card_path.write_bytes((json.dumps(card, indent=1, ensure_ascii=True, allow_nan=False, default=str) + "\n")
                          .encode("ascii"))
    return card


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", default=None, help="the fine-tuned model file (default: the release model)")
    ap.add_argument("--cpsam-v2", default=None, help="cpsam_v2 weights (default: Cellpose's cache)")
    ap.add_argument("--cpsam", default=None, help="cpsam (v1) weights (default: Cellpose's cache)")
    ap.add_argument("--write", action="store_true", help="write the result into the pack card and rebuild its sums")
    a = ap.parse_args(argv)
    model = Path(a.model) if a.model else resolve_registered("germline_nuclei_combined")
    if model is None or not Path(model).is_file():
        print("ERROR: the release model is not on disk; `germquant get-model` or --model PATH", file=sys.stderr)
        return 1
    if Path(model).stat().st_size != MODELS["germline_nuclei_combined"]["size"] and not a.model:
        print("ERROR: the release model has the wrong size; fetch it again", file=sys.stderr)
        return 1
    cache = cellpose_cache_dirs()[0]
    bases = {"cpsam_v2": Path(a.cpsam_v2) if a.cpsam_v2 else cache / "cpsam_v2",
             "cpsam": Path(a.cpsam) if a.cpsam else cache / "cpsam"}
    missing = [f"{n} ({p})" for n, p in bases.items() if not p.is_file()]
    if missing:
        print(f"ERROR: base weights not found: {', '.join(missing)}", file=sys.stderr)
        return 1
    rec = compare_files(model, bases)
    print(json.dumps(rec, indent=1, default=str))
    if a.write:
        update_card(CARD, rec, bases)
        write_sums(PACK)
        print(f"wrote {CARD} and {PACK / 'SHA256SUMS'}: base {rec['result']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
