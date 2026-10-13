"""Build tests/golden/rad51_additions.json, the record of the gonads protocol changelog entry 9 adds to TRAIN.

    python scripts/make_rad51_additions.py [--copy-manifest C:/Users/ryane/n2_rad51_hs_2026/copy_manifest.json]
        [--truth-root C:/Users/ryane/germquant_golden/imaris_truth] [--out tests/golden/rad51_additions.json]
        [--no-hash] [--no-metadata]

The four heat-shocked N2 males imaged on 2026-01-05 (docs/BENCHMARK_PROTOCOL.md, changelog entry 9), copied read
only from the lab NAS to C:/Users/ryane/n2_rad51_hs_2026/: for each the .nd2, the Imaris .ims and the 2026-01-13
statistics export (.xlsx). The script

* re-hashes every local copy and compares size and sha256 with the copy manifest (a difference stops it);
* reads the .nd2 header (axis sizes, voxel size, channel names, camera, exposures; no pixel);
* reads the Imaris facts of each .ims with h5py (spot and surface counts, region growing, the Quality threshold,
  the spot source channel, the extents as stored text, the per-axis Imaris voxel, channel names, the stored spot
  statistics, Masked 405, the creation parameters) and compares the creation parameters with protocol sections
  4.1 and 4.2 and with the eight registered train and val projects (their truth params JSONs under
  --truth-root; a test gonad is refused before any read);
* refuses, among the additions' files, any 2025-11-05 path and any aggregate of the per-nucleus tables
  (benchmark_split.SEALED_NAME_MARKERS); the only 2025-11-05 files it reads are the eight registered train and
  val truth params JSONs above, behind the role guard (protocol changelog entry 10 (g)(2));
* writes the record, validated by germquant.validate.benchmark_split.validate_additions: roles fixed by
  ADDITION_ROLES (train), the registered split file's sha256 and the sha256 of its roles.

The xlsx exports are hashed, not opened, by this script (the decision-1 calibration of 2026-10-11 opened them;
protocol changelog entry 10): their Overall totals were read by the inventory of 2026-10-11 and are copied here
from its export_checks.json. The record never changes after it is committed; the verification of
20260105_n2_hs_1_MALE_1's xlsx spots (entry 9 (e)) goes into its truth params JSON.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from germquant.benchmark import detector as bd  # noqa: E402
from germquant.validate import benchmark_split as bs  # noqa: E402

DATA = Path("C:/Users/ryane/n2_rad51_hs_2026")
DEFAULT_COPY_MANIFEST = DATA / "copy_manifest.json"
DEFAULT_TRUTH_ROOT = Path("C:/Users/ryane/germquant_golden/imaris_truth")
SESSION = "2026-01-05"
NAS_SET = "CahoonLab/Madeleine Coffman/Tiny Storage Madeleine/Madeleine/N2/"
OVERALL_SOURCE = ("the Overall sheet totals read by the NAS inventory of 2026-10-11 "
                  "(germquant_inventory/ims_inspection/export_checks.json); not re-read here")
IMS_NUCLEUS_TRUTH = ("Imaris Surfaces (the lab's wizard recipe, independent of germquant) with the Masked 405 "
                     "channel; recovered as in docs/IMARIS_TRUTH.md")


def _spec(tag, set_name, stem, n_spots, n_surfaces, spot_truth, xlsx_use, notes):
    folder = DATA / set_name
    return {"tag": tag, "session": SESSION, "set": NAS_SET + set_name, "genotype": "N2", "stem": stem,
            "nd2": folder / f"{stem}.nd2", "ims": folder / f"{stem}.ims",
            "xlsx": DATA / "imaris_exports_20260113" / f"{stem}.xlsx", "xlsx_use": xlsx_use,
            "overall_n_spots": n_spots, "overall_n_surfaces": n_surfaces, "overall_source": OVERALL_SOURCE,
            "spot_truth": spot_truth, "nucleus_truth": IMS_NUCLEUS_TRUTH, "notes": notes}


SPECS = [
    _spec("n2rad51_20260105a_HS_male_001", "20260105_N2_HS", "20260105_n2_hs_MALE_1", 5187, 98,
          {"source": "ims", "status": "usable", "rule": "entry 9 (d)"}, "cross-check (4.4)", []),
    _spec("n2rad51_20260105a_HS_male_002", "20260105_N2_HS", "20260105_n2_hs_MALE_2", 12028, 88,
          {"source": "ims", "status": "usable", "rule": "entry 9 (d)"}, "cross-check (4.4)", []),
    _spec("n2rad51_20260105a_HS_male_003", "20260105_N2_HS", "20260105_n2_hs_MALE_3", 6462, 127,
          {"source": "ims", "status": "usable", "rule": "entry 9 (d)"}, "cross-check (4.4)", []),
    _spec("n2rad51_20260105b_HS_male_001", "20260105_N2_HS_1", "20260105_n2_hs_1_MALE_1", 6512, 112,
          {"source": "xlsx", "status": "pending", "rule": "entry 9 (e)"},
          "spot source (entry 9 (e): used only if the criteria V1 to V7 pass)",
          ["The saved Spots object of the .ims (re-saved 2026-01-21) holds 1,226,263 spots under the recorded "
           "filter (Quality 12.5, region growing off): it is not the curated set. The 2026-01-13 export holds the "
           "curated 6,512 spots and 112 surfaces. Its spot truth is used only if the export verifies (entry 9 (e)); "
           "otherwise the gonad is nucleus truth only."]),
]

NOT_ADDED = [
    {"gonads": ["20260105_n2_hs_1_MALE_2", "20260105_n2_hs_1_MALE_3", "20260105_n2_hs_1_MALE_4",
                "20260105_n2_hs_2_MALE_1", "20260105_n2_hs_2_MALE_2", "20260105_n2_hs_2_MALE_3"],
     "where": NAS_SET + "20260105_N2_HS_1 and 20260105_N2_HS_2",
     "reason": ("raw only: an nd2 and a plain .ims conversion without Spots or Surfaces. Not copied, not opened and "
                "not used for anything.")},
    {"gonads": ["20260521_n2_hs_dapi_rad51_syp1_M", "_M001", "_M002", "_M003", "and the 2026-05-21 N2 HS herms"],
     "where": "CahoonLab/Madeleine Coffman/imaris_files/n2_ims_files/20260521_N2_HS_DAPI_RAD51_SYP1_ims",
     "reason": ("spot truth only (their Surfaces are germquant nuclei imported into Imaris) and RAD-51 on 477, "
                "another session. Not used for anything; entry 9 (c) opens no image of another session for "
                "method development.")},
]

DECIDED_BY = {
    "who": "Ryan Eastman (project lead)",
    "date": "2026-10-11",
    "question": ("How should the 3 new full-truth males be used? I recommend training data, recorded in a "
                 "changelog entry"),
    "reply": "2 whatever you recommend",
    "reading": ("read as accepting the recommendation (train), for the three males and for "
                "20260105_n2_hs_1_MALE_1 if its spot truth verifies (entry 9 (e)); a reading, not his words"),
}


def registered_creation_parameters(truth_root: Path) -> dict:
    """tag -> creation_parameter_summary of the eight registered train and val projects, from their truth
    params JSONs. A test gonad is refused before anything is read."""
    split = bs.load_split(REPO / "tests" / "golden" / "rad51_split.json")
    out = {}
    for g in split["gonads"]:
        if g["role"] not in ("train", "val"):
            continue
        bd.guard_test([g["tag"]])
        p = truth_root / g["tag"] / f"{g['stem']}__imaris_params.json"
        if not p.is_file():
            continue
        cp = json.loads(p.read_text(encoding="utf-8"))["creation_parameters"]
        out[g["tag"]] = bs.creation_parameter_summary(cp.get("Spots 1"), cp.get("Surfaces 1"))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copy-manifest", default=str(DEFAULT_COPY_MANIFEST))
    ap.add_argument("--truth-root", default=str(DEFAULT_TRUTH_ROOT),
                    help="the registered truth exports (for the train and val creation parameters)")
    ap.add_argument("--out", default=str(REPO / "tests" / "golden" / "rad51_additions.json"))
    ap.add_argument("--no-hash", action="store_true", help="take the sha256 from the copy manifest (no re-hash)")
    ap.add_argument("--no-metadata", action="store_true", help="skip the nd2 and .ims reads")
    args = ap.parse_args(argv)

    registered = {} if args.no_metadata else registered_creation_parameters(Path(args.truth_root))
    rec = bs.build_additions(SPECS, args.copy_manifest, registered_split_path=REPO / "tests" / "golden" / "rad51_split.json",
                             decided_by=DECIDED_BY, registered_params=registered, not_added=NOT_ADDED,
                             hash_files=not args.no_hash, read_metadata=not args.no_metadata)
    rec["registered_creation_parameters_compared"] = sorted(registered)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=1) + "\n", encoding="utf-8")
    print(f"{'tag':31s} {'spots':>8s} {'surf':>5s} {'RG':>5s} {'Quality':>8s} {'voxel (Imaris, x y z)':28s} spot truth")
    for g in rec["gonads"]:
        ims = g["ims"]
        print(f"{g['tag']:31s} {str(ims.get('n_spots')):>8s} {str(ims.get('n_surfaces')):>5s} "
              f"{str(ims.get('region_growing')):>5s} {str(ims.get('spot_quality_threshold')):>8s} "
              f"{str(ims.get('imaris_voxel_um')):28s} {g['spot_truth']['source']}/{g['spot_truth']['status']}")
        for d in g.get("creation_parameter_differences", []):
            print(f"    differs: {d}")
    print(f"wrote {out}: {len(rec['gonads'])} additions, roles {rec['roles']}; split file unchanged: "
          f"{rec['registered_split']['unchanged']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
