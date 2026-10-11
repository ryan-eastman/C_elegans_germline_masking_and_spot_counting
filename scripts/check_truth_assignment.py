"""Protocol section 4.4: check the geometric spot-to-nucleus assignment of the truth against the lab's own.

    python scripts/check_truth_assignment.py [--manifest PATH] [--tags TAG ...] [--out CSV]

The benchmark assigns a truth spot to a traced nucleus geometrically (``assign_to_masks``: the recovered
mask label at the spot's centre voxel, else the most frequent label in its 3 x 3 x 3 neighbourhood, with mask
fragments below the protocol's size removed). The lab assigns it through the Imaris mask channel: Imaris's
per-spot maximum of the "Masked" channel, exported by ``germquant imaris-truth`` as ``imaris_mask_label``. The
protocol asks for the two to be compared once per gonad before any scoring; a disagreement beyond boundary
spots is a truth-export bug to fix first. This reads only the truth files (no method output).

By default the train and val gonads of the manifest are checked; the test gonads are checked in the
pre-flight of docs/FINAL_SCORING.md by naming them with --tags.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from germquant.benchmark import detector as bd  # noqa: E402
from germquant.benchmark.manifest import load_manifest, truth_file_names  # noqa: E402
from germquant.benchmark.spots import TRUTH_MASK_MIN_VOXELS, assign_to_masks, truth_mask_lut  # noqa: E402

COLUMNS = ["tag", "split", "n_spots", "n_agree", "n_agree_assigned", "n_geometric_only", "n_coloc_only",
           "n_both_different", "n_coloc_on_removed_fragment", "agreement", "coloc_label_column"]


def check_gonad(truth_labels: np.ndarray, spots: pd.DataFrame, *, min_voxels: int = TRUTH_MASK_MIN_VOXELS) -> dict:
    """Per-spot comparison of the geometric assignment with ``imaris_mask_label`` for one gonad."""
    lut, keep, _sizes, frag = truth_mask_lut(truth_labels, min_voxels)
    zi, yi, xi = (spots[c].to_numpy(np.int64) for c in ("zi", "yi", "xi"))
    geo = assign_to_masks(truth_labels, zi, yi, xi, lut=lut)
    if "imaris_mask_label" not in spots.columns:
        return {"n_spots": int(len(spots)), "coloc_label_column": False}
    col = pd.to_numeric(spots["imaris_mask_label"], errors="coerce").fillna(0).to_numpy(np.int64)
    kept = set(int(v) for v in keep)
    on_frag = (col > 0) & ~np.isin(col, list(kept))
    agree = geo == col
    return {
        "n_spots": int(len(spots)),
        "n_agree": int(agree.sum()),
        "n_agree_assigned": int((agree & (geo > 0)).sum()),
        "n_geometric_only": int(((geo > 0) & (col == 0)).sum()),
        "n_coloc_only": int(((geo == 0) & (col > 0) & ~on_frag).sum()),
        "n_both_different": int(((geo > 0) & (col > 0) & (geo != col)).sum()),
        "n_coloc_on_removed_fragment": int(on_frag.sum()),
        "agreement": float(agree.mean()) if len(spots) else float("nan"),
        "coloc_label_column": True,
        "fragments_removed": [int(v) for v in frag],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", default=str(bd.DEFAULT_MANIFEST))
    ap.add_argument("--tags", nargs="*", help="gonads to check (default: every train and val gonad)")
    ap.add_argument("--out", help="write the per-gonad table here (CSV)")
    a = ap.parse_args(argv)
    entries = load_manifest(a.manifest)
    if a.tags:
        want = set(a.tags)
        entries = [e for e in entries if e.tag in want]
    else:
        entries = [e for e in entries if e.split in ("train", "val")]
    rows = []
    import tifffile

    for e in entries:
        if e.truth_dir is None:
            print(f"{e.tag}: no truth folder in the manifest; skipped")
            continue
        names = truth_file_names(e.truth_stem)
        spots = pd.read_csv(Path(e.truth_dir) / names["spots"], low_memory=False)
        labels = tifffile.imread(str(Path(e.truth_dir) / names["labels"]))
        r = check_gonad(labels, spots)
        r.update(tag=e.tag, split=e.split)
        rows.append(r)
        print(f"{e.tag:24s} {e.split:5s} spots {r['n_spots']:7d}  agree {r.get('agreement', float('nan')):.4f}  "
              f"geometric-only {r.get('n_geometric_only', '-')}  coloc-only {r.get('n_coloc_only', '-')}  "
              f"different {r.get('n_both_different', '-')}  on removed fragment {r.get('n_coloc_on_removed_fragment', '-')}")
    table = pd.DataFrame(rows)
    for c in COLUMNS:
        if c not in table.columns:
            table[c] = np.nan
    table = table[COLUMNS + [c for c in table.columns if c not in COLUMNS]]
    if a.out:
        table.to_csv(a.out, index=False)
        print(f"-> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
