"""The SC ribbon benchmark of docs/SC_VALIDATION_PROTOCOL.md (sc_straighten v1).

    python scripts/benchmark_sc.py manifest  --out DIR
    python scripts/benchmark_sc.py synthetic --out DIR [--families iso dense ...] [--n N] [--seed0 S] [--procs P]
    python scripts/benchmark_sc.py zones     --runs ROOT --out DIR [--gonads G ...] [--n-boot 1000]
    python scripts/benchmark_sc.py filaments --out DIR --final --reason TEXT          (check S3, test only)
    python scripts/benchmark_sc.py report    --out DIR

Development by default: development gonads (the split's role 'dev') and development seeds (below 7,000,000).
Test or reserved gonads, the live Filament file and the sealed seeds only with ``--final --reason TEXT``, which
appends a line to C:/Users/ryane/germquant_golden/benchmarks/sc_ribbon_v1/final_runs.jsonl (commit, params_hash,
split sha256, time, reason). The test set is scored once, by the orchestrator, outside the build workflow
(protocol section 10). Every development configuration evaluated on zone metrics is a line of dev_trials.jsonl
(at most 40 per protocol section 9).

``zones`` reads the sc-straighten outputs of each gonad under ``--runs`` (the folder ``germquant sc-straighten
--out`` wrote, mirroring the results tree; a gonad's folder is found by its RAD-51 tag or image stem), the run's
nuclei.csv next to them, and the gonad's zone truth through the guarded loader. Frame: the .ims header rule and
the nd2 stage rule are both computed and must agree within 0.05 um in xy (protocol 4.1); the four RAD-51 test
gonads use the nd2 rule only (their .ims is never opened).
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from germquant.benchmark import sc_ribbons as B  # noqa: E402
from germquant.benchmark import sc_split  # noqa: E402

BENCH = Path("C:/Users/ryane/germquant_golden/benchmarks/sc_ribbon_v1")
REPO = Path(__file__).resolve().parents[1]


def _git_sha() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def _append(path: Path, rec: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, default=lambda o: None if o != o else str(o)) + "\n")


def _final_guard(a) -> None:
    if a.final and not a.reason:
        raise SystemExit("--final needs --reason TEXT (protocol section 10)")
    if a.final:
        from germquant.sc.straighten import params as PRM

        _append(BENCH / "final_runs.jsonl", {"utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                             "command": a.cmd, "commit": _git_sha(),
                                             "params_hash": PRM.resolve().hash,
                                             "split_sha256": sc_split.REGISTERED_SHA256, "reason": a.reason})


# ------------------------------------------------------------------------------------ manifest
def cmd_manifest(a) -> int:
    split = sc_split.load_split()
    rows = []
    for g in split["gonads"]:
        if g.get("geno") != "N2":
            continue
        rows.append({k: g.get(k) for k in ("gonad", "role", "session", "condition", "sex", "image_stem",
                                           "rad51_tag", "frame_source", "nd2", "ims", "sc_channel", "golden_dir")})
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "sc_manifest.csv", index=False)
    print(f"{len(rows)} N2 gonads -> {out / 'sc_manifest.csv'}")
    return 0


# ------------------------------------------------------------------------------------ synthetic
def cmd_synthetic(a) -> int:
    _final_guard(a)
    from germquant.sc.straighten import synth

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.final:
        ranges = sc_split.load_split()["synthetic"]["sealed_seed_ranges"]
        plan = {fam: list(range(lo, hi + 1)) for fam, (lo, hi) in ranges.items() if fam in a.families}
    else:
        plan = {fam: list(range(a.seed0, a.seed0 + a.n)) for fam in a.families}
        for seeds in plan.values():
            for s in seeds:
                synth.check_seed(s, final=False)
    cmd = [sys.executable, str(REPO / "scripts" / "calibrate_sc_straighten.py"), "--out", str(out), "--procs",
           str(a.procs)]
    results = {}
    for fam, seeds in plan.items():
        sub = cmd + ["--families", fam, "--n", str(len(seeds)), "--seed0", str(seeds[0])]
        if a.final:
            raise SystemExit("sealed synthetic scoring is run by the orchestrator at the frozen commit "
                             "(protocol section 10); this build does not generate sealed seeds")
        subprocess.run(sub, check=True)
        results[fam] = len(seeds)
    print(f"synthetic families {sorted(results)} -> {out / 'calibration.jsonl'}")
    return 0


# ------------------------------------------------------------------------------------ zones
def _find_run(runs: Path, rec: dict) -> Path | None:
    keys = [rec.get("rad51_tag"), rec.get("image_stem")]
    for f in sorted(runs.rglob("*__sc_straighten.json")):
        try:
            r = json.loads(f.read_text(encoding="utf-8"))
        except ValueError:
            continue
        if r.get("done") and (f.parent.name in keys or r.get("image_id") in keys):
            return f.parent
    return None


def _frames(rec: dict) -> dict:
    out = {}
    if rec.get("role") == "dev" and rec.get("ims") and "(not opened)" not in str(rec["ims"]) and \
            Path(rec["ims"]).is_file():
        out["ims"] = B.ims_frame(rec["ims"])
    if rec.get("nd2") and Path(rec["nd2"]).is_file():
        out["nd2"] = B.nd2_frame(rec["nd2"])
    return out


def cmd_zones(a) -> int:
    _final_guard(a)
    split = sc_split.load_split()
    master = split["masters"]["n2_and_dlw190_jan"]["path"]
    roles = ("dev",) if not a.final else ("test",)
    gon = [g for g in split["gonads"] if g.get("geno") == "N2" and g.get("role") in roles]
    if a.gonads:
        gon = [g for g in gon if g["gonad"] in set(a.gonads) or g.get("rad51_tag") in set(a.gonads)]
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    runs = Path(a.runs)
    per, payload, skipped = [], [], []
    for g in gon:
        folder = _find_run(runs, g)
        if folder is None:
            skipped.append({"gonad": g["gonad"], "reason": "no sc-straighten run"})
            continue
        rec = json.loads(next(folder.glob("*__sc_straighten.json")).read_text(encoding="utf-8"))
        iid = rec["image_id"]
        nuc_run = pd.read_csv(Path(rec["inputs"]["nuclei_csv"]["path"]), low_memory=False)
        sc_nuc = pd.read_csv(folder / f"{iid}__sc_strand_nuclei.csv", low_memory=False)
        v_nd2 = (float(nuc_run["voxel_dx_um"].iloc[0]), float(nuc_run["voxel_dy_um"].iloc[0]))
        fr = _frames(g)
        if "ims" in fr and "nd2" in fr:
            d = float(np.max(np.abs(np.asarray(fr["ims"]["ext_min_xy"]) - np.asarray(fr["nd2"]["ext_min_xy"]))))
            if d > 0.05:
                skipped.append({"gonad": g["gonad"], "reason": f"frame rules disagree by {d:.3f} um"})
                continue
        frame = fr.get("ims") or fr.get("nd2")
        if frame is None:
            skipped.append({"gonad": g["gonad"], "reason": "no frame (ims header or nd2 stage)"})
            continue
        try:
            zt = B.zone_truth(g, master, v_nd2, frame=frame, final=a.final, split=split)
        except B.TruthDefect as e:
            skipped.append({"gonad": g["gonad"], "reason": f"truth defect: {e}"})
            continue
        germ = nuc_run[nuc_run["in_germline"].astype(str).str.lower().isin(["true", "1"])]
        st = sc_nuc.set_index("nucleus_id")["status"]
        keep = germ[germ["nucleus_id"].map(st).isin(["traced", "no_sc"])]
        n_failed = int((germ["nucleus_id"].map(st) == "failed").sum())
        lab = B.label_nuclei(zt, keep)
        cols = ["piece_lw_mean_um", "total_length_um", "n_pieces", "gaps_per_10um", "longest_piece_um",
                "on_strand_frac"]
        vals = sc_nuc.set_index("nucleus_id")
        aucs = {c: B.gonad_aucs(lab, vals[c]) for c in cols if c in vals}
        aucs["central_element_mean_intensity"] = B.gonad_aucs(
            lab, nuc_run.set_index("nucleus_id")["central_element_mean_intensity"])
        row = {"gonad": g["gonad"], "rad51_tag": g.get("rad51_tag"), "condition": g["condition"], "sex": g["sex"],
               "image_id": iid, "frame": frame["source"], "n_analysis": int(lab["in_analysis"].sum()),
               "n_failed": n_failed, "failed_frac": n_failed / max(len(germ), 1),
               **{f"{c}__{k}": v for c, d in aucs.items() for k, v in d.items()}}
        per.append(row)
        payload.append((g["gonad"], lab, sc_nuc, zt))
        print(f"{g['gonad']}: S2a AUC {row.get('piece_lw_mean_um__auc_full_vs_frag', float('nan')):.3f} "
              f"(n {row.get('piece_lw_mean_um__n_full')} vs {row.get('piece_lw_mean_um__n_frag')}), "
              f"S2b AUC {row.get('total_length_um__auc_sc_vs_none', float('nan')):.3f}, "
              f"intensity S2a {row.get('central_element_mean_intensity__auc_full_vs_frag', float('nan')):.3f}",
              flush=True)
    per_df = pd.DataFrame(per)
    per_df.to_csv(out / "zones_per_gonad.csv", index=False)
    summary = {"n_gonads": len(per), "skipped": skipped}
    if per:
        from germquant.benchmark.stats import cluster_bootstrap

        def mean_auc(col):
            def stat(rows):
                v = [r[col] for r in rows if r.get(col) is not None and np.isfinite(r[col])]
                return float(np.mean(v)) if v else float("nan")
            return stat

        for col in ("piece_lw_mean_um__auc_full_vs_frag", "total_length_um__auc_sc_vs_none",
                    "central_element_mean_intensity__auc_full_vs_frag",
                    "central_element_mean_intensity__auc_sc_vs_none", "on_strand_frac__auc_full_vs_frag",
                    "n_pieces__auc_full_vs_frag", "gaps_per_10um__auc_full_vs_frag",
                    "longest_piece_um__auc_full_vs_frag"):
            if col not in per_df:
                continue
            ev = [r for r in per if r.get(col) is not None and np.isfinite(r[col])]
            bs = cluster_bootstrap(ev, mean_auc(col), n_boot=a.n_boot, seed=B.BOOT_SEED)
            summary[col] = {"n_evaluable": len(ev), **bs["stats"]["value"]}
        th = B.choose_thresholds([(lab, nuc) for _g, lab, nuc, _zt in payload])
        summary["theta_full"] = th["theta_full"]
        summary["theta_present"] = th["theta_present"]
        hits = []
        for gname, lab, nuc, zt in payload:
            est = B.estimate_boundaries(lab, nuc, th["theta_full"], th["theta_present"])
            row = {"gonad": gname}
            for key in ("SC_fragstart", "SC_fullstart", "SC_fraglate"):
                h, err = B.boundary_hit(zt, key, est.get(key))
                row.update({f"{key}_est": est.get(key), f"{key}_truth_s": float(zt.s_points[zt.idx[key] - 1]),
                            f"{key}_hit": h, f"{key}_abs_err_um": err,
                            f"{key}_abs_err_frac": err / float(zt.s_points[zt.idx['Pachy_end'] - 1])})
            hits.append(row)
        hdf = pd.DataFrame(hits)
        hdf.to_csv(out / "zones_boundaries.csv", index=False)
        summary["S2c_fullstart_hit_rate"] = float(hdf["SC_fullstart_hit"].mean())
        found = hdf["SC_fullstart_abs_err_um"].dropna()
        summary["S2c_fullstart_median_err_um"] = float(found.median()) if len(found) else float("nan")
        summary["fragstart_hit_rate"] = float(hdf["SC_fragstart_hit"].mean())
        summary["fraglate_hit_rate"] = float(hdf["SC_fraglate_hit"].mean())
    (out / "zones_summary.json").write_text(json.dumps(summary, indent=1, default=lambda o: None if o != o else o))
    if not a.final:
        from germquant.sc.straighten import params as PRM

        hashes = sorted({json.loads(next(p.glob('*__sc_straighten.json')).read_text())["params_hash"]
                         for p in {(_find_run(runs, g) or Path('.')) for g in gon} if p != Path('.')})
        _append(BENCH / "dev" / "dev_trials.jsonl", {
            "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "commit": _git_sha(),
            "params_hashes": hashes, "default_params_hash": PRM.resolve().hash, "changed_keys": a.changed or [],
            "n_gonads": len(per),
            "S2a": summary.get("piece_lw_mean_um__auc_full_vs_frag"),
            "S2b": summary.get("total_length_um__auc_sc_vs_none"),
            "S2c_hit_rate": summary.get("S2c_fullstart_hit_rate"),
            "S2c_median_err_um": summary.get("S2c_fullstart_median_err_um"),
            "theta_full": summary.get("theta_full"), "theta_present": summary.get("theta_present"),
            "out": str(out)})
    print(json.dumps({k: v for k, v in summary.items() if k != "skipped"}, indent=1, default=str)[:3000])
    if skipped:
        print("skipped:", skipped)
    return 0


# ------------------------------------------------------------------------------------ filaments (S3)
def cmd_filaments(a) -> int:
    if not a.final:
        raise SystemExit("the live Filament file is test truth (protocol 2.2, 3.1 rule 6): check S3 runs only with "
                         "--final --reason TEXT, by the orchestrator")
    _final_guard(a)
    raise SystemExit("check S3 is scored by the orchestrator at the frozen commit (protocol section 10)")


def cmd_report(a) -> int:
    out = Path(a.out)
    for name in ("zones_summary.json",):
        p = out / name
        if p.is_file():
            print(p.read_text(encoding="utf-8")[:4000])
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("manifest", "synthetic", "zones", "filaments", "report"):
        s = sub.add_parser(name)
        s.add_argument("--out", required=True)
        s.add_argument("--final", action="store_true")
        s.add_argument("--reason", default=None)
        if name == "synthetic":
            s.add_argument("--families", nargs="*", default=["iso", "dense", "gaps", "foci", "pairs", "noise", "snr"])
            s.add_argument("--n", type=int, default=24)
            s.add_argument("--seed0", type=int, default=2000)
            s.add_argument("--procs", type=int, default=6)
        if name == "zones":
            s.add_argument("--runs", required=True)
            s.add_argument("--gonads", nargs="*", default=None)
            s.add_argument("--n-boot", type=int, default=1000)
            s.add_argument("--changed", nargs="*", default=None, help="the keys this development trial changed")
    a = ap.parse_args(argv)
    return {"manifest": cmd_manifest, "synthetic": cmd_synthetic, "zones": cmd_zones, "filaments": cmd_filaments,
            "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
