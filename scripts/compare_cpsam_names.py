#!/usr/bin/env python
"""Do the bare name cpsam and the pack model cpsam_v2@0f1cc3f7 give the same masks? (docs/PACKS.md, "cpsam and
cpsam_v2@0f1cc3f7"; docs/DESIGN_system_packs.md, section 3.5; package P3, acceptance 5)

    python scripts/compare_cpsam_names.py [--image IMG] [--device auto|cuda|mps|cpu] [--examples DIR] [--json OUT]

nuclei._cellpose maps the bare name cpsam to CellposeModel() without a weights file, and Cellpose 4.2.x then
loads its default weights, cpsam_v2, from its own cache. A pack names the weights by digest instead:
cpsam_v2@0f1cc3f7 resolves to the blob in the model store and loads CellposeModel(pretrained_model=<blob>).
Cellpose builds the network from the weights' state dict, not from the file name, so the two loads should give the
same masks. This script runs both through nuclei._cellpose_2d (the call a projected recipe makes) on the plane the
budding yeast recipe segments: the sharpest plane of the transmitted channel of the pack's real example (fetched
first with `germquant pack fetch-example budding_yeast`, or given with --image). It compares the label images pixel
for pixel and prints one JSON record.

A device asked for explicitly must be the device that runs: when it is not available, or $GERMQUANT_DEVICE
overrides it, the script stops with exit code 2 before anything runs, so a CUDA check can never pass on the CPU.
The pack blob is re-hashed here (not trusted from its verification record) and must match its pin. The record
holds the device asked for and the device that ran, whether CUDA was available and the GPU's name, the blob's
sha256, the file the bare name loads (Cellpose's cache copy) and whether it is the same file as the blob:
`germquant pack fetch-models` hard-links the blob to that copy when it can, and then the check compares two load
paths of one file. It also records the engine version and git sha, so a record shows which script wrote it.
Exit codes: 0 identical masks, 1 different masks, 2 the requested device cannot run or the blob is not the pinned
file.

On the RTX 5090, through the lab's GPU lock (the record must then say "device": "cuda" and "cuda_available": true):

    python C:/Users/ryane/germquant_golden/gpu_lock.py run --name "P3 cpsam names" -- .venv/Scripts/python.exe scripts/compare_cpsam_names.py --device cuda --json cpsam_names_cuda.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import germquant.pipeline  # noqa: E402,F401  (torch first: the cuBLAS load order on the RTX 5090)

PACK, MODEL = "budding_yeast", "cpsam_v2@0f1cc3f7"
DEVICES = ("auto", "cuda", "mps", "cpu")


def resolve_device(requested: str) -> tuple[str, str]:
    """(the device germquant will run on, why that is not the device asked for, or ""). ``auto`` takes what
    select_device picks; an explicit device that select_device would replace is a reason to stop."""
    from germquant.device import select_device

    requested = str(requested).strip().lower()
    if requested == "auto":
        return select_device(None), ""
    device = select_device(requested)               # not quiet: a fallback is logged as well
    if device == requested:
        return device, ""
    env = (os.environ.get("GERMQUANT_DEVICE") or "").strip().lower()
    why = f"$GERMQUANT_DEVICE={env} overrides it" if env and env != requested else f"{requested} is not available here"
    return device, f"--device {requested} was asked for, but germquant would run on {device} ({why})"


def main(argv=None) -> int:
    import numpy as np

    from germquant import packs, planar
    from germquant.assemble import read_ome_tiff
    from germquant.datafetch import sha256_file
    from germquant.io.channel_map import ChannelMap
    from germquant.packs.blobs import cellpose_cache_dirs
    from germquant.segment import nuclei

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--image", default=None, help="an OME-TIFF (default: the budding yeast pack's real example)")
    ap.add_argument("--device", default="auto", choices=DEVICES, help="auto, cuda, mps or cpu (default auto)")
    ap.add_argument("--examples", default=None, help="the examples folder for the default image")
    ap.add_argument("--json", default=None, help="also write the record to this file")
    a = ap.parse_args(argv)

    device, why = resolve_device(a.device)
    if why:
        print(f"REFUSED: {why}; nothing was run", file=sys.stderr)
        return 2
    pack = packs.find_pack(PACK)
    spec = pack.model(MODEL)
    blob = packs.resolve_pack_model(MODEL)
    digest = sha256_file(blob)                      # the full file, not its verification record
    if digest != spec.sha256:
        print(f"REFUSED: {blob} has sha256 {digest}, the pack pins {spec.sha256}; nothing was run", file=sys.stderr)
        return 2
    cache = next((d / spec.file_name for d in cellpose_cache_dirs() if (d / spec.file_name).is_file()), None)
    same = bool(os.path.samefile(blob, cache)) if cache is not None else None
    cache_sha = digest if same else (sha256_file(cache) if cache is not None else None)

    image = Path(a.image) if a.image else Path(packs.fetch_example(pack, None, a.examples, say=print)["image"])
    img = read_ome_tiff(image)
    cmap = ChannelMap.from_yaml(pack.folder / "config" / "channel_maps" / "yeast_transmitted_foci.yaml")
    roles, flags = cmap.resolve(img["channel_names"])
    if roles.get("transmitted") is None:
        raise SystemExit(f"{image}: no transmitted channel among {img['channel_names']} ({flags})")
    stack = img["data"][roles["transmitted"]]
    z = planar.resolve_plane(stack, "sharpest")
    plane = planar.project_image(stack, "plane", "sharpest")
    out = {}
    for label, model in (("cpsam", "cpsam"), (MODEL, str(blob))):
        t0 = time.time()
        lab = nuclei._cellpose_2d(plane, model, device)[0]
        out[label] = {"labels": lab, "n_objects": int(lab.max()), "seconds": round(time.time() - t0, 1)}
    a_lab, b_lab = out["cpsam"]["labels"], out[MODEL]["labels"]
    from importlib.metadata import version

    import torch

    cuda = bool(torch.cuda.is_available())
    rec = {"image": str(image), "pixels_sha256": img["pixels_sha256"], "plane": int(z),
           "requested_device": a.device, "device": device, "cuda_available": cuda,
           "cuda_device": torch.cuda.get_device_name(0) if cuda else None,
           "identical": bool(np.array_equal(a_lab, b_lab)),
           "pixels_differing": int((a_lab != b_lab).sum()), "pixels_total": int(a_lab.size),
           "n_objects": {k: v["n_objects"] for k, v in out.items()},
           "seconds": {k: v["seconds"] for k, v in out.items()},
           "model_blob": str(blob), "model_sha256": digest, "cellpose_cache_file": str(cache) if cache else None,
           "cellpose_cache_sha256": cache_sha, "same_file_as_cellpose_cache": same,
           "cellpose": version("cellpose"), "torch": torch.__version__, "engine": packs.version_info(),
           "script": "scripts/compare_cpsam_names.py"}
    print(json.dumps(rec, indent=1))
    if a.json:
        Path(a.json).write_text(json.dumps(rec, indent=1) + "\n", encoding="ascii")
    return 0 if rec["identical"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
