"""Protocol section 4.4: check the geometric spot-to-nucleus assignment of the truth against the lab's own.

    python scripts/check_truth_assignment.py [--manifest PATH] [--tags TAG ...] [--out CSV]
        [--additions tests/golden/rad51_additions.json]

The benchmark assigns a truth spot to a traced nucleus geometrically (``assign_to_masks``: the recovered
mask label at the spot's centre voxel, else the most frequent label in its 3 x 3 x 3 neighbourhood, with mask
fragments below the protocol's size removed). The lab assigns it through the Imaris mask channel: Imaris's
per-spot maximum of the "Masked" channel, exported by ``germquant imaris-truth`` as ``imaris_mask_label``. The
protocol asks for the two to be compared once per gonad before any scoring; a disagreement beyond boundary
spots is a truth-export bug to fix first. This reads only the truth files (no method output).

By default the train and val gonads of the manifest are checked; the test gonads are checked in the
pre-flight of docs/FINAL_SCORING.md by naming them with --tags.

--additions (the record of protocol changelog entry 9) adds, for every checked addition with an xlsx export,
the cross-checks of entry 9 (d)(5): the per-nucleus counts of the lab's xlsx coloc join (spots grouped by the
mask channel's Intensity Max, the value of their surface) against the same join on the truth spots and against
the geometric assignment, per surface; and the xlsx spot positions against the truth spots (by ID, then one
to one by position within 0.005 um). Only the additions' own exports are read (never a 2025-11-05 file).
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


POSITION_TOL_UM = 0.005


def xlsx_cross_check(truth_labels: np.ndarray, spots: pd.DataFrame, nuclei: pd.DataFrame, xl: dict, *,
                     min_voxels: int = TRUTH_MASK_MIN_VOXELS) -> dict:
    """Entry 9 (d)(5) for one addition: ``xl`` is ``imaris_xlsx.read_export_objects`` of its export. Per surface
    of the export (keyed by its mask-channel Intensity Max, the coloc value), the export's coloc count against
    the truth spots' coloc count (``imaris_mask_label``) and against the geometric count of the truth nucleus
    with the same surface ID; and the export's spots against the truth spots by ID and by position."""
    from scipy.spatial import cKDTree

    xs, xsu = xl["spots"], xl["surfaces"]
    x_count = xs["imax_mask"].round().astype("Int64").value_counts()
    col = pd.to_numeric(spots.get("imaris_mask_label"), errors="coerce").fillna(0).round().astype(np.int64)
    t_count = col.value_counts()
    lut, _keep, _sizes, _frag = truth_mask_lut(truth_labels, min_voxels)
    geo = assign_to_masks(truth_labels, *(spots[c].to_numpy(np.int64) for c in ("zi", "yi", "xi")), lut=lut)
    g_count = pd.Series(geo[geo > 0]).value_counts()
    lab_of = dict(zip(nuclei["surface_id"].astype(int), nuclei["label"].fillna(0).astype(int)))
    rows = []
    for sid, v in zip(xsu["ID"].astype(int), xsu["imax_mask"].round()):
        v = int(v) if np.isfinite(v) else 0
        lab = lab_of.get(int(sid), 0)
        rows.append({"surface_id": int(sid), "coloc_value": v, "xlsx_coloc": int(x_count.get(v, 0)) if v else 0,
                     "truth_coloc": int(t_count.get(v, 0)) if v else 0,
                     "truth_geometric": int(g_count.get(lab, 0)) if lab else 0, "truth_label": lab})
    per = pd.DataFrame(rows)
    out = {"n_xlsx_surfaces": int(len(xsu)), "n_xlsx_spots": int(len(xs)), "n_truth_spots": int(len(spots))}
    if len(per):
        out.update(
            n_surfaces_coloc_equal=int((per["xlsx_coloc"] == per["truth_coloc"]).sum()),
            n_surfaces_geometric_equal=int((per["xlsx_coloc"] == per["truth_geometric"]).sum()),
            mean_abs_diff_coloc=float((per["xlsx_coloc"] - per["truth_coloc"]).abs().mean()),
            mean_abs_diff_geometric=float((per["xlsx_coloc"] - per["truth_geometric"]).abs().mean()),
            n_surfaces_without_truth_label=int((per["truth_label"] == 0).sum()),
            sum_xlsx_coloc=int(per["xlsx_coloc"].sum()), sum_truth_coloc=int(per["truth_coloc"].sum()),
            sum_truth_geometric=int(per["truth_geometric"].sum()))
    ids_x, ids_t = set(xs["ID"].astype(int)), set(spots["spot_id"].astype(int))
    common = sorted(ids_x & ids_t)
    out.update(n_spot_ids_common=len(common), n_spot_ids_only_xlsx=len(ids_x - ids_t),
               n_spot_ids_only_truth=len(ids_t - ids_x))
    if common:
        a = xs.set_index("ID").loc[common, ["x", "y", "z"]].to_numpy(float)
        b = spots.set_index("spot_id").loc[common, ["imaris_x_um", "imaris_y_um", "imaris_z_um"]].to_numpy(float)
        out["max_position_diff_common_ids_um"] = float(np.linalg.norm(a - b, axis=1).max())
    if len(xs) and len(spots):
        tp = spots[["imaris_x_um", "imaris_y_um", "imaris_z_um"]].to_numpy(float)
        d, j = cKDTree(tp).query(xs[["x", "y", "z"]].to_numpy(float))
        ok = d <= POSITION_TOL_UM
        out["n_xlsx_spots_within_tol_of_a_truth_spot"] = int(ok.sum())
        out["one_to_one_within_tol"] = bool(ok.all() and len(set(j[ok].tolist())) == int(ok.sum())
                                            and len(xs) == len(spots))
    out["per_surface"] = per
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", default=str(bd.DEFAULT_MANIFEST))
    ap.add_argument("--tags", nargs="*", help="gonads to check (default: every train and val gonad)")
    ap.add_argument("--out", help="write the per-gonad table here (CSV)")
    ap.add_argument("--additions", default=None,
                    help="the additions record of protocol changelog entry 9: cross-check each addition's xlsx export")
    a = ap.parse_args(argv)
    entries = load_manifest(a.manifest)
    additions = {}
    if a.additions:
        from germquant.validate.benchmark_split import load_additions

        additions = {g["tag"]: g for g in load_additions(a.additions)["gonads"]}
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
        print(f"{e.tag:24s} {e.split:5s} spots {r['n_spots']:7d}  agree {r.get('agreement', float('nan')):.4f}  "
              f"geometric-only {r.get('n_geometric_only', '-')}  coloc-only {r.get('n_coloc_only', '-')}  "
              f"different {r.get('n_both_different', '-')}  on removed fragment {r.get('n_coloc_on_removed_fragment', '-')}")
        g = additions.get(e.tag)
        if g is not None and g.get("xlsx"):
            from germquant.validate.imaris_xlsx import read_export_objects

            params = {}
            if (Path(e.truth_dir) / names["params"]).is_file():
                import json

                params = json.loads((Path(e.truth_dir) / names["params"]).read_text(encoding="utf-8"))
            src = (params.get("frame_check") or {}).get("spot_channel")
            mask = (params.get("mask") or {}).get("channel")
            xl = read_export_objects(g["xlsx"]["path"], spot_channel=int(src if src is not None else 2),
                                     mask_channel=int(mask if mask is not None else 3))
            nuclei = pd.read_csv(Path(e.truth_dir) / names["nuclei"])
            xc = xlsx_cross_check(labels, spots, nuclei, xl)
            per = xc.pop("per_surface")
            if a.out:
                per.assign(tag=e.tag).to_csv(Path(a.out).with_name(Path(a.out).stem + f"__{e.tag}__xlsx_per_surface.csv"),
                                             index=False)
            r.update({f"xlsx_{k}": v for k, v in xc.items()})
            r["xlsx_use"] = g["xlsx"].get("use")
            print(f"    xlsx ({g['xlsx'].get('use')}): {xc['n_xlsx_spots']} spots, {xc['n_xlsx_surfaces']} surfaces; "
                  f"per-surface coloc equal {xc.get('n_surfaces_coloc_equal')}, geometric equal "
                  f"{xc.get('n_surfaces_geometric_equal')}; spot IDs common {xc['n_spot_ids_common']} "
                  f"(only xlsx {xc['n_spot_ids_only_xlsx']}, only truth {xc['n_spot_ids_only_truth']}), one to one "
                  f"within {POSITION_TOL_UM} um: {xc.get('one_to_one_within_tol')}")
        rows.append(r)
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
