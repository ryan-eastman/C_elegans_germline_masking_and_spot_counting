"""Prepare annotation data and fine-tune a germline Cellpose model.

Zero-shot Cellpose under-segments crowded pachytene nuclei (Jaccard 0.47 -> 0.78 after
fine-tuning, microPub Biology 2023). Workflow:
  1. `germquant prep-training` -> exports DAPI z-slices as TIFs for hand annotation.
  2. Annotate in the Cellpose GUI (or napari) -> *_seg.npy / *_masks.tif next to each TIF.
  3. `germquant finetune` -> trains from a pretrained model, writes models/<name>.
  4. Point config.segmentation.cellpose_model at the new model path.
See docs/ANNOTATION.md.
"""
from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path

import numpy as np

from ..io import read_stack
from ..io.channel_map import ChannelMap

log = logging.getLogger(__name__)


def _stretch_uint16(a: np.ndarray) -> np.ndarray:
    a = a.astype(np.float32)
    lo, hi = np.percentile(a, 1), np.percentile(a, 99.8)
    a = np.clip((a - lo) / (hi - lo + 1e-9), 0, 1)
    return (a * 65535).astype(np.uint16)


def prep_training_data(
    nd2_paths: list[str | Path],
    channel_map: ChannelMap,
    out_dir: str | Path,
    *,
    n_slices: int = 3,
    xy_stride: int = 1,
) -> list[Path]:
    """Export evenly-spaced DAPI z-slices (full-res by default) as TIFs to annotate."""
    import tifffile

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for p in nd2_paths:
        p = Path(p)
        from ..io import read_nd2_metadata

        role_to_idx, _ = channel_map.resolve(read_nd2_metadata(p)["channel_names"])
        dna_idx = role_to_idx.get("dna")
        if dna_idx is None:
            dna_idx = 0
            log.warning("%s: no DAPI role; exporting channel 0", p.name)
        stack = read_stack(p, xy_stride=xy_stride)
        dna = stack.channel(dna_idx)
        nz = dna.shape[0]
        z_idx = np.linspace(nz * 0.2, nz * 0.8, n_slices).astype(int) if nz > 1 else [0]
        for z in z_idx:
            out = out_dir / f"{p.stem}_z{int(z):03d}.tif"
            tifffile.imwrite(str(out), _stretch_uint16(dna[z]))
            written.append(out)
    log.info("Wrote %d annotation slices -> %s", len(written), out_dir)
    return written


def finetune_cellpose(
    labeled_dir: str | Path,
    out_model: str | Path,
    *,
    pretrained: str = "cpsam",
    n_epochs: int = 100,
    learning_rate: float = 1e-5,
    run: bool = True,
) -> list[str]:
    """Build (and optionally run) the Cellpose training command. GPU-only.

    Expects `labeled_dir` to contain the annotation TIFs + their *_seg.npy masks
    (from the Cellpose GUI). Returns the command that was/would be run.
    """
    labeled_dir = Path(labeled_dir)
    out_model = Path(out_model)
    out_model.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, "-m", "cellpose", "--train",
        "--dir", str(labeled_dir),
        "--pretrained_model", pretrained,
        "--chan", "0", "--chan2", "0",
        "--n_epochs", str(n_epochs),
        "--learning_rate", str(learning_rate),
        "--model_name_out", str(out_model),
        "--verbose",
    ]
    log.info("Cellpose training command:\n  %s", " ".join(cmd))
    if run:
        try:
            subprocess.run(cmd, check=True)
        except FileNotFoundError:
            log.error("cellpose not installed — install germquant[gpu] on the 5090/HPC.")
        except subprocess.CalledProcessError as e:
            log.error("cellpose training failed (%s). Run the printed command manually; "
                      "flags vary by Cellpose version — see docs/ANNOTATION.md.", e)
    return cmd
