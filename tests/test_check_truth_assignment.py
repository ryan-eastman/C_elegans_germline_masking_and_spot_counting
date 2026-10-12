"""Protocol section 4.4: the geometric truth assignment checked against the lab's mask-channel assignment."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

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


def test_xlsx_cross_check_counts_per_surface_and_matches_spots():
    """Entry 9 (d)(5): the xlsx coloc join per surface against the truth's coloc and geometric counts, and the
    xlsx spots against the truth spots by ID and by position."""
    mod = _load()
    lab = np.zeros((10, 20, 20), np.int32)
    lab[2:8, 2:12, 2:12] = 5
    lab[2:8, 14:19, 14:19] = 9
    spots = pd.DataFrame({"spot_id": [1, 2, 3, 4], "zi": [4, 4, 4, 0], "yi": [5, 6, 16, 0], "xi": [5, 6, 16, 0],
                          "imaris_mask_label": [5, 5, 9, 0],
                          "imaris_x_um": [10.0, 11.0, 12.0, 13.0], "imaris_y_um": 1.0, "imaris_z_um": 2.0})
    nuclei = pd.DataFrame({"surface_id": [100, 101], "label": [5, 9]})
    xl = {"spots": pd.DataFrame({"ID": [1, 2, 3, 4], "x": [10.0, 11.0, 12.0, 13.001], "y": 1.0, "z": 2.0,
                                 "imax_mask": [5.0, 5.0, 9.0, 0.0]}),
          "surfaces": pd.DataFrame({"ID": [100, 101], "imax_mask": [5.0, 9.0]})}
    r = mod.xlsx_cross_check(lab, spots, nuclei, xl, min_voxels=10)
    per = r.pop("per_surface")
    assert per[["surface_id", "xlsx_coloc", "truth_coloc", "truth_geometric"]].values.tolist() == [[100, 2, 2, 2],
                                                                                                    [101, 1, 1, 1]]
    assert r["n_surfaces_coloc_equal"] == 2 and r["n_surfaces_geometric_equal"] == 2
    assert r["n_spot_ids_common"] == 4 and r["max_position_diff_common_ids_um"] == pytest.approx(0.001)
    assert r["one_to_one_within_tol"] is True
    # a spot the export has and the truth does not, and a moved one
    xl["spots"].loc[3, "x"] = 20.0
    xl["spots"] = pd.concat([xl["spots"], pd.DataFrame({"ID": [7], "x": [30.0], "y": [1.0], "z": [2.0],
                                                        "imax_mask": [9.0]})], ignore_index=True)
    r = mod.xlsx_cross_check(lab, spots, nuclei, xl, min_voxels=10)
    assert r["n_spot_ids_only_xlsx"] == 1 and r["one_to_one_within_tol"] is False
    assert r["n_xlsx_spots_within_tol_of_a_truth_spot"] == 3 and r["n_surfaces_coloc_equal"] == 1


def test_the_additions_checks_refuse_a_test_gonad(tmp_path):
    """Protocol changelog entry 10: with --additions a test gonad is refused before any truth file is read (its
    truth folder does not even exist here); the test gonads are checked only in the pre-flight, without
    --additions."""
    import json

    from germquant.benchmark.detector import TestSetGuardError

    mod = _load()
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"gonads": [
        {"tag": "n2rad51_HS_male_001", "split": "test", "role": "test", "results_dir": str(tmp_path / "r"),
         "truth_dir": str(tmp_path / "missing"), "stem": "20251105_n2_hs_MALE _001"}]}), encoding="utf-8")
    with pytest.raises(TestSetGuardError):
        mod.main(["--manifest", str(manifest), "--tags", "n2rad51_HS_male_001",
                  "--additions", str(ROOT / "tests" / "golden" / "rad51_additions.json")])
