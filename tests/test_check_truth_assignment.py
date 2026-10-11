"""Protocol section 4.4: the geometric truth assignment checked against the lab's mask-channel assignment."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location("check_truth_assignment", ROOT / "scripts" / "check_truth_assignment.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_check_gonad_counts_agreement_and_each_kind_of_disagreement():
    mod = _load()
    lab = np.zeros((10, 20, 20), np.int32)
    lab[2:8, 2:12, 2:12] = 5                                      # a kept nucleus (600 voxels)
    lab[2:8, 14:19, 14:19] = 9                                    # a second kept nucleus
    lab[0, 0, 0] = 7                                              # a one-voxel fragment
    spots = pd.DataFrame({
        "zi": [4, 4, 4, 4, 9, 0],
        "yi": [5, 0, 16, 5, 19, 0],
        "xi": [5, 19, 16, 5, 19, 0],
        "imaris_mask_label": [5, 0, 5, 0, 9, 7],
    })
    r = mod.check_gonad(lab, spots, min_voxels=100)
    assert r["n_spots"] == 6 and r["coloc_label_column"]
    assert r["n_agree"] == 2                                      # spot 0 (5 = 5) and spot 1 (0 = 0, far from any mask)
    assert r["n_both_different"] == 1                             # spot 2: geometric 9, coloc 5
    assert r["n_geometric_only"] == 1                             # spot 3: geometric 5, coloc 0
    assert r["n_coloc_only"] == 1                                 # spot 4: geometric 0, coloc 9
    assert r["n_coloc_on_removed_fragment"] == 1                  # spot 5: coloc 7, a removed fragment
    assert r["fragments_removed"] == [7]
