#!/usr/bin/env python
"""Do the bare name cpsam and the pack model cpsam_v2@0f1cc3f7 give the same masks? (docs/PACKS.md, "cpsam and
cpsam_v2@0f1cc3f7"; docs/DESIGN_system_packs.md, section 3.5; package P3, acceptance 5)

    python scripts/compare_cpsam_names.py [--image IMG] [--device auto|cuda|mps|cpu] [--json OUT]

nuclei._cellpose maps the bare name cpsam to CellposeModel() without a weights file, and Cellpose 4.2.x then
loads its default weights, cpsam_v2. A pack names the weights by digest instead: cpsam_v2@0f1cc3f7 resolves to the
verified blob in the model store and loads CellposeModel(pretrained_model=<blob>). Cellpose builds the network
from the weights' state dict, not from the file name, so the two loads should give the same masks. This script
runs both through nuclei._cellpose_2d (the call a projected recipe makes) on the plane the budding yeast recipe
segments: the sharpest plane of the transmitted channel of the pack's real example (fetched first with
`germquant pack fetch-example budding_yeast`, or given with --image). It compares the label images pixel for
pixel and prints one JSON record: identical or not, the object counts, the device that ran and the versions.

On the RTX 5090, through the lab's GPU lock:

    python C:/Users/ryane/germquant_golden/gpu_lock.py run --name "P3 cpsam names" -- \
        .venv/Scripts/python.exe scripts/compare_cpsam_names.py --device cuda --json cpsam_names_cuda.json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import germquant.pipeline  # noqa: E402,F401  (torch first: the cuBLAS load order on the RTX 5090)

PACK, MODEL = "budding_yeast", "cpsam_v2@0f1cc3f7"


def main(argv=None) -> int:
    import numpy as np

    from germquant import packs, planar
    from germquant.assemble import read_ome_tiff
    from germquant.device import select_device
    from germquant.io.channel_map import ChannelMap
    from germquant.segment import nuclei

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--image", default=None, help="an OME-TIFF (default: the budding yeast pack's real example)")
    ap.add_argument("--device", default="auto", help="auto, cuda, mps or cpu (default auto)")
    ap.add_argument("--examples", default=None, help="the examples folder for the default image")
    ap.add_argument("--json", default=None, help="also write the record to this file")
    a = ap.parse_args(argv)

    pack = packs.find_pack(PACK)
    image = Path(a.image) if a.image else Path(packs.fetch_example(pack, None, a.examples, say=print)["image"])
    img = read_ome_tiff(image)
    cmap = ChannelMap.from_yaml(pack.folder / "config" / "channel_maps" / "yeast_transmitted_foci.yaml")
    roles, flags = cmap.resolve(img["channel_names"])
    if roles.get("transmitted") is None:
        raise SystemExit(f"{image}: no transmitted channel among {img['channel_names']} ({flags})")
    stack = img["data"][roles["transmitted"]]
    z = planar.resolve_plane(stack, "sharpest")
    plane = planar.project_image(stack, "plane", "sharpest")
    blob = packs.resolve_pack_model(MODEL)
    device = select_device(None if a.device == "auto" else a.device, quiet=True)
    out = {}
    for label, model in (("cpsam", "cpsam"), (MODEL, str(blob))):
        t0 = time.time()
        lab = nuclei._cellpose_2d(plane, model, device)[0]
        out[label] = {"labels": lab, "n_objects": int(lab.max()), "seconds": round(time.time() - t0, 1)}
    a_lab, b_lab = out["cpsam"]["labels"], out[MODEL]["labels"]
    from importlib.metadata import version

    import torch

    rec = {"image": str(image), "pixels_sha256": img["pixels_sha256"], "plane": int(z), "device": device,
           "identical": bool(np.array_equal(a_lab, b_lab)),
           "pixels_differing": int((a_lab != b_lab).sum()), "pixels_total": int(a_lab.size),
           "n_objects": {k: v["n_objects"] for k, v in out.items()},
           "seconds": {k: v["seconds"] for k, v in out.items()},
           "model_blob": str(blob), "cellpose": version("cellpose"), "torch": torch.__version__,
           "cuda_available": bool(torch.cuda.is_available())}
    print(json.dumps(rec, indent=1))
    if a.json:
        Path(a.json).write_text(json.dumps(rec, indent=1) + "\n", encoding="ascii")
    return 0 if rec["identical"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
