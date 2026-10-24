"""The SC ribbon benchmark of docs/SC_VALIDATION_PROTOCOL.md (sc_straighten v1).

    python scripts/benchmark_sc.py manifest  --out DIR
    python scripts/benchmark_sc.py synthetic --out DIR [--families iso dense ...] [--n N] [--seed0 S] [--procs P]
                                             [--generator entry2|registered]
    python scripts/benchmark_sc.py zones     --runs ROOT --out DIR [--gonads G ...] [--n-boot 1000] [--s1 FILE]
    python scripts/benchmark_sc.py zones     --final --reason TEXT --results-root ROOT --runs DIR --out DIR
                                             --theta-full T --theta-present T --expect-params-hash H
                                             [--params-json F | --set KEY=VALUE ...] [--source-map OLD=NEW ...]
    python scripts/benchmark_sc.py synthetic --final --reason TEXT --out DIR [--expect-params-hash H] ...
    python scripts/benchmark_sc.py filaments --out DIR --final --reason TEXT          (check S3, test only)
    python scripts/benchmark_sc.py report    --out DIR

Development by default: development gonads (the split's role 'dev') and development seeds (below 7,000,000).
Test or reserved gonads, the live Filament file and the sealed seeds only with ``--final --reason TEXT``. Every
--final run writes a start line to C:/Users/ryane/germquant_golden/benchmarks/sc_ribbon_v1/final_runs.jsonl
before it touches any sealed data and an end line when it finishes (commit, params_hash, split sha256, time,
reason, status); a refused command writes nothing. The test set is scored once, by the orchestrator, outside the
build workflow (protocol section 10). Every development configuration evaluated on zone metrics is a line of
dev_trials.jsonl (at most 40 per protocol section 9).

``zones`` reads the sc-straighten outputs of each gonad under ``--runs`` (the folder ``germquant sc-straighten
--out`` wrote, mirroring the results tree; a gonad's folder is found by its RAD-51 tag or image stem), the run's
nuclei.csv next to them, and the gonad's zone truth through the guarded loader. With ``--final`` it first runs
sc-straighten itself on the test images under ``--results-root`` with the frozen parameters (``--allow-sealed``
with the reason, logged), and the thresholds of the boundary estimator are the frozen ``--theta-full`` and
``--theta-present`` (protocol 5.2 and 9: chosen on development gonads only, never on the test set).

Frame (protocol 4.1): each gonad uses its registered rule (``frame_source`` of the split: the .ims header, or
the nd2 stage rule for the four RAD-51 test gonads, whose .ims is never opened). When both rules can be computed
they are compared on the image corners (maximum |difference| over x and y, voxel included); on development
gonads a difference above 0.05 um skips the gonad.

Pre-flight (protocol 4.5): the master table, the R script and every scored zone export match their registered
sha256 before they are parsed; each converted polyline lies inside its image's xy extent within 1 um; the run's
central_element channel is the session's SC channel (the nd2 channel name). Every gonad in scope gets a row of
zones_per_gonad.csv, with ``status`` "no output: <reason>" when it could not be scored, plus a pooled row.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
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
FINAL_LOG = BENCH / "final_runs.jsonl"
REPO = Path(__file__).resolve().parents[1]
FRAME_TOL_UM = 0.05
EXTENT_TOL_UM = 1.0


def _git_sha() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def _git_dirty() -> bool | None:
    try:
        out = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=REPO,
                             capture_output=True, text=True, check=True).stdout
        return bool(out.strip())
    except Exception:  # noqa: BLE001
        return None


def _append(path: Path, rec: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, default=lambda o: None if o != o else str(o)) + "\n")


def _utc() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _need_reason(a) -> None:
    if a.final and not (a.reason and a.reason.strip()):
        raise SystemExit("--final needs --reason TEXT (protocol section 10)")


def _start_final(a, **extra) -> str:
    """The start line of a --final run, written before any sealed data is touched; returns the run id."""
    run_id = f"{a.cmd}-{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}-{os.getpid()}"
    _append(FINAL_LOG, {"event": "start", "run_id": run_id, "utc": _utc(), "command": a.cmd, "commit": _git_sha(),
                        "dirty": _git_dirty(), "split_sha256": sc_split.REGISTERED_SHA256, "reason": a.reason,
                        **extra})
    return run_id


def _end_final(run_id: str, status: str, **extra) -> None:
    _append(FINAL_LOG, {"event": "end", "run_id": run_id, "utc": _utc(), "status": status, **extra})


def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _params(a):
    from germquant.sc.straighten import params as PRM

    section = json.loads(Path(a.params_json).read_text(encoding="utf-8")) if getattr(a, "params_json", None) else None
    return PRM.resolve(section, list(getattr(a, "set", None) or []))


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
    _need_reason(a)
    from germquant.sc.straighten import synth

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.final:
        ranges = sc_split.load_split()["synthetic"]["sealed_seed_ranges"]
        plan = {fam: (lo, hi - lo + 1) for fam, (lo, hi) in ranges.items() if fam in a.families}
    else:
        plan = {fam: (a.seed0, a.n) for fam in a.families}
        for s0, n in plan.values():
            for s in range(s0, s0 + n):
                synth.check_seed(s, final=False)
    params = _params(a)
    if a.final and a.expect_params_hash and params.hash != a.expect_params_hash:
        raise SystemExit(f"params_hash {params.hash} is not the frozen {a.expect_params_hash}")
    sets = list(getattr(a, "set", None) or [])
    base = [sys.executable, str(REPO / "scripts" / "calibrate_sc_straighten.py"), "--out", str(out), "--procs",
            str(a.procs)] + (["--fixed", *sets] if sets else [])
    # at the scoring, S1b and S1e are reported on both generators (protocol changelog entry 5)
    gens = {fam: (["entry2", "registered"] if (a.final and fam in ("dense", "snr")) else [a.generator])
            for fam in plan}
    run_id = _start_final(a, params_hash=params.hash, families=sorted(plan)) if a.final else None
    status = "failed"
    try:
        for fam, (s0, n) in plan.items():
            for gen in gens[fam]:
                cmd = base + ["--families", fam, "--n", str(n), "--seed0", str(s0), "--generator", gen]
                if a.final:
                    cmd += ["--final", "--reason", a.reason]
                subprocess.run(cmd, check=True)
        status = "done"
    finally:
        if run_id:
            _end_final(run_id, status, out=str(out))
    print(f"synthetic families {sorted(plan)} -> {out / 'calibration.jsonl'}")
    return 0


# ------------------------------------------------------------------------------------ zones
def _find_run(runs: Path, rec: dict) -> Path | None:
    keys = [k for k in (rec.get("rad51_tag"), rec.get("image_stem")) if k]
    for f in sorted(runs.rglob("*__sc_straighten.json")):
        try:
            r = json.loads(f.read_text(encoding="utf-8"))
        except ValueError:
            continue
        if r.get("done") and (f.parent.name in keys or r.get("image_id") in keys):
            return f.parent
    return None


def _frames(rec: dict, size_xy=None) -> dict:
    """The frame of the gonad's registered rule (``frame_source``) and, when it can be computed, the other one.
    The .ims of a gonad registered with the nd2 stage rule is never opened."""
    out = {}
    src = rec.get("frame_source")
    ims = str(rec.get("ims") or "")
    if src == "ims header" and ims and "(not opened)" not in ims and Path(ims).is_file():
        out["ims"] = B.ims_frame(ims)
    if rec.get("nd2") and Path(rec["nd2"]).is_file():
        out["nd2"] = B.nd2_frame(rec["nd2"])
    reg = "ims" if src == "ims header" else ("nd2" if src == "nd2 stage rule" else None)
    out["registered"] = reg
    return out


def _image_extent_um(sc_rec: dict, nuc_run: pd.DataFrame) -> tuple[float, float] | None:
    """(x, y) extent of the run's image in pipeline um, from the label TIFF's header (no pixel is read)."""
    try:
        import tifffile

        with tifffile.TiffFile(sc_rec["inputs"]["labels_tif"]["path"]) as t:
            shp = t.series[0].shape
        return (float(shp[-1]) * float(nuc_run["voxel_dx_um"].iloc[0]),
                float(shp[-2]) * float(nuc_run["voxel_dy_um"].iloc[0]))
    except Exception:  # noqa: BLE001
        return None


def _sc_channel(sc_rec: dict, g: dict) -> tuple[bool | None, str]:
    """Pre-flight 4.5 item 5: the nd2 channel name of the run's central_element index against the session's SC
    channel. (None, why) when it cannot be read."""
    try:
        import nd2

        src = (sc_rec.get("inputs", {}).get("source") or {}).get("path") or g.get("nd2")
        ci = int(sc_rec["inputs"]["channel_map"]["central_element"])
        with nd2.ND2File(str(src)) as f:
            name = str(f.metadata.channels[ci].channel.name).strip()
        want = str(g.get("sc_channel", "")).strip()
        return (want in name) if want else None, f"central_element = channel {ci} '{name}', session SC {want}"
    except Exception as e:  # noqa: BLE001
        return None, f"channel check unavailable ({type(e).__name__}: {e})"


def _no_output(g: dict, reason: str) -> dict:
    return {"gonad": g["gonad"], "role": g.get("role"), "rad51_tag": g.get("rad51_tag"), "session": g.get("session"),
            "condition": g.get("condition"), "sex": g.get("sex"), "status": f"no output: {reason}"}


def _preflight_masters(split: dict) -> None:
    for key in ("n2_and_dlw190_jan", "zone_r_script"):
        m = split["masters"].get(key)
        if not m:
            continue
        h = _sha256(m["path"]) if Path(m["path"]).is_file() else "missing"
        if h != m.get("sha256"):
            raise SystemExit(f"truth defect (protocol 4.5 item 1): {m['path']} sha256 {h} is not the registered "
                             f"{m.get('sha256')}; record it in the changelog before scoring")


def cmd_zones(a) -> int:
    _need_reason(a)
    if a.final:
        missing = [n for n, v in (("--theta-full", a.theta_full), ("--theta-present", a.theta_present),
                                  ("--results-root", a.results_root), ("--expect-params-hash", a.expect_params_hash))
                   if v is None]
        if missing:
            raise SystemExit(f"--final needs {', '.join(missing)} (the freeze entry of protocol section 9: the "
                             "thresholds are chosen on development gonads only, and sc-straighten runs on the test "
                             "images with the frozen parameters)")
        if a.guard_frac != B.GUARD_FRAC_REGISTERED and not a.allow_unregistered_guard:
            raise SystemExit("--final uses the registered guard band (protocol 4.2) unless a changelog entry changed "
                             "it: pass --allow-unregistered-guard with the entry's value")
    split = sc_split.load_split()
    _preflight_masters(split)
    master = split["masters"]["n2_and_dlw190_jan"]["path"]
    roles = ("dev",) if not a.final else ("test",)
    gon = [g for g in split["gonads"] if g.get("geno") == "N2" and g.get("role") in roles]
    if a.gonads:
        gon = [g for g in gon if g["gonad"] in set(a.gonads) or g.get("rad51_tag") in set(a.gonads)]
    assert all(g["role"] in roles for g in gon)
    params = _params(a)
    if a.final and params.hash != a.expect_params_hash:
        raise SystemExit(f"params_hash {params.hash} is not the frozen {a.expect_params_hash}")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    runs = Path(a.runs)
    run_id = None
    if a.final:
        run_id = _start_final(a, params_hash=params.hash, theta_full=a.theta_full, theta_present=a.theta_present,
                              gonads=[g["gonad"] for g in gon], guard_frac=a.guard_frac)
        from germquant.sc.straighten.runner import run_results

        smap = [tuple(s.split("=", 1)) for s in (a.source_map or [])]
        run_results(a.results_root, [g["image_stem"] for g in gon], params=params, out=runs, allow_sealed=True,
                    sealed_reason=a.reason, audit_log=FINAL_LOG, workers=a.workers, source_map=smap)
    per, payload, skipped = [], [], []
    status = "failed"
    try:
        for g in gon:
            exp = g.get("zone_export_copy") or {}
            if not exp.get("path") or not Path(exp["path"]).is_file() or _sha256(exp["path"]) != exp.get("sha256"):
                per.append(_no_output(g, "truth defect: the zone export is missing or its sha256 is not registered"))
                skipped.append({"gonad": g["gonad"], "reason": per[-1]["status"]})
                continue
            folder = _find_run(runs, g)
            if folder is None:
                per.append(_no_output(g, "no sc-straighten run"))
                skipped.append({"gonad": g["gonad"], "reason": "no sc-straighten run"})
                continue
            rec = json.loads(next(folder.glob("*__sc_straighten.json")).read_text(encoding="utf-8"))
            iid = rec["image_id"]
            nuc_run = pd.read_csv(Path(rec["inputs"]["nuclei_csv"]["path"]), low_memory=False)
            sc_nuc = pd.read_csv(folder / f"{iid}__sc_strand_nuclei.csv", low_memory=False)
            v_nd2 = (float(nuc_run["voxel_dx_um"].iloc[0]), float(nuc_run["voxel_dy_um"].iloc[0]))
            ch_ok, ch_txt = _sc_channel(rec, g)
            if ch_ok is False:
                per.append(_no_output(g, f"SC channel mismatch ({ch_txt})"))
                skipped.append({"gonad": g["gonad"], "reason": per[-1]["status"]})
                continue
            fr = _frames(g)
            reg = fr.get("registered")
            frame = fr.get(reg) if reg else None
            if frame is None:
                per.append(_no_output(g, f"no frame (registered rule {g.get('frame_source')!r} not available)"))
                skipped.append({"gonad": g["gonad"], "reason": per[-1]["status"]})
                continue
            extent = _image_extent_um(rec, nuc_run)
            d = float("nan")
            if "ims" in fr and "nd2" in fr and extent is not None:
                d = B.frame_disagreement(fr["ims"], fr["nd2"], v_nd2, extent)
                if d > FRAME_TOL_UM and not a.final:
                    per.append(_no_output(g, f"frame rules disagree by {d:.3f} um"))
                    skipped.append({"gonad": g["gonad"], "reason": per[-1]["status"]})
                    continue
            try:
                zt = B.zone_truth(g, master, v_nd2, frame=frame, final=a.final, split=split)
            except B.TruthDefect as e:
                per.append(_no_output(g, f"truth defect: {e}"))
                skipped.append({"gonad": g["gonad"], "reason": per[-1]["status"]})
                continue
            if extent is not None:
                P = zt.points_um
                inside = ((P[:, 0] >= -EXTENT_TOL_UM) & (P[:, 0] <= extent[0] + EXTENT_TOL_UM)
                          & (P[:, 1] >= -EXTENT_TOL_UM) & (P[:, 1] <= extent[1] + EXTENT_TOL_UM))
                if not inside.all():
                    per.append(_no_output(g, f"truth defect: {int((~inside).sum())} polyline points outside the "
                                             "image extent by more than 1 um (protocol 4.5 item 3)"))
                    skipped.append({"gonad": g["gonad"], "reason": per[-1]["status"]})
                    continue
            germ = nuc_run[nuc_run["in_germline"].astype(str).str.lower().isin(["true", "1"])]
            st = sc_nuc.set_index("nucleus_id")["status"]
            keep = germ[germ["nucleus_id"].map(st).isin(["traced", "no_sc"])]
            n_failed = int((germ["nucleus_id"].map(st) == "failed").sum())
            lab = B.label_nuclei(zt, keep, guard_frac=a.guard_frac)
            cols = ["piece_lw_mean_um", "total_length_um", "n_pieces", "gaps_per_10um", "longest_piece_um",
                    "on_strand_frac", "profile_cv_corr"]
            vals = sc_nuc.set_index("nucleus_id")
            aucs = {c: B.gonad_aucs(lab, vals[c]) for c in cols if c in vals}
            aucs["central_element_mean_intensity"] = B.gonad_aucs(
                lab, nuc_run.set_index("nucleus_id")["central_element_mean_intensity"])
            failed_frac = n_failed / max(len(germ), 1)
            row = {"gonad": g["gonad"], "role": g.get("role"), "rad51_tag": g.get("rad51_tag"),
                   "session": g.get("session"), "condition": g["condition"], "sex": g["sex"], "status": "scored",
                   "image_id": iid, "frame": frame["source"], "frame_rules_max_diff_um": d,
                   "sc_channel_check": ch_txt, "run_git_sha": rec.get("git_sha"),
                   "run_code_sha256": rec.get("code_sha256"), "params_hash": rec.get("params_hash"),
                   "n_analysis": int(lab["in_analysis"].sum()), "n_failed": n_failed, "failed_frac": failed_frac,
                   "failed_over_5pct": bool(failed_frac > 0.05),
                   **{f"{c}__{k}": v for c, dd in aucs.items() for k, v in dd.items()}}
            per.append(row)
            payload.append((g["gonad"], lab, sc_nuc, zt))
            print(f"{g['gonad']}: S2a AUC {row.get('piece_lw_mean_um__auc_full_vs_frag', float('nan')):.3f} "
                  f"(n {row.get('piece_lw_mean_um__n_full')} vs {row.get('piece_lw_mean_um__n_frag')}), "
                  f"S2b AUC {row.get('total_length_um__auc_sc_vs_none', float('nan')):.3f}, "
                  f"intensity S2a {row.get('central_element_mean_intensity__auc_full_vs_frag', float('nan')):.3f}",
                  flush=True)
        scored = [r for r in per if r.get("status") == "scored"]
        auc_cols = sorted({k for r in scored for k in r if k.endswith("__auc_full_vs_frag") or
                           k.endswith("__auc_sc_vs_none") or k.endswith("__auc_full_vs_scfrag") or
                           k.endswith("__auc_full_vs_lpfrag")})
        pooled = {"gonad": "pooled", "status": f"mean over evaluable gonads ({len(scored)} of {len(per)} scored)"}
        for c in auc_cols:
            v = [r[c] for r in scored if r.get(c) is not None and np.isfinite(r[c])]
            pooled[c] = float(np.mean(v)) if v else float("nan")
            pooled[c.replace("__auc_", "__n_evaluable_")] = len(v)
        per_df = pd.DataFrame(per + [pooled])
        per_df.to_csv(out / "zones_per_gonad.csv", index=False)
        summary = {"n_gonads": len(scored), "n_in_scope": len(per), "skipped": skipped, "guard_frac": a.guard_frac,
                   "guard_registered": bool(a.guard_frac == B.GUARD_FRAC_REGISTERED), "params_hash": params.hash,
                   "final": bool(a.final)}
        if scored:
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
                        "longest_piece_um__auc_full_vs_frag", "profile_cv_corr__auc_full_vs_frag"):
                if col not in per_df:
                    continue
                ev = [r for r in scored if r.get(col) is not None and np.isfinite(r[col])]
                bs = cluster_bootstrap(ev, mean_auc(col), n_boot=a.n_boot, seed=B.BOOT_SEED)
                summary[col] = {"n_evaluable": len(ev), **bs["stats"]["value"]}
            if a.final:
                th = {"theta_full": float(a.theta_full), "theta_present": float(a.theta_present),
                      "source": "frozen (command line, the freeze entry)"}
            else:
                assert all(sc_split.role_of(gname) == "dev" for gname, *_ in payload), "development gonads only"
                th = B.choose_thresholds([(lab, nuc) for _g, lab, nuc, _zt in payload])
                th["source"] = "chosen on these development gonads (protocol 5.2)"
            summary["theta_full"] = th["theta_full"]
            summary["theta_present"] = th["theta_present"]
            summary["theta_source"] = th["source"]
            hits = []
            for gname, lab, nuc, zt in payload:
                est = B.estimate_boundaries(lab, nuc, th["theta_full"], th["theta_present"])
                row = {"gonad": gname}
                for key in ("SC_fragstart", "SC_fullstart", "SC_fraglate"):
                    h, err = B.boundary_hit(zt, key, est.get(key))
                    row.update({f"{key}_est": est.get(key), f"{key}_found": est.get(key) is not None,
                                f"{key}_truth_s": float(zt.s_points[zt.idx[key] - 1]),
                                f"{key}_hit": h, f"{key}_abs_err_um": err,
                                f"{key}_abs_err_frac": err / float(zt.s_points[zt.idx['Pachy_end'] - 1])})
                hits.append(row)
            for r in per:
                if r.get("status") != "scored":
                    hits.append({"gonad": r["gonad"], "status": r["status"], "SC_fullstart_hit": False,
                                 "SC_fragstart_hit": False, "SC_fraglate_hit": False})
            hdf = pd.DataFrame(hits)
            hdf.to_csv(out / "zones_boundaries.csv", index=False)
            # rates over the gonads with output (protocol section 10: criteria are evaluated over the images with
            # output and labelled with their number; the rows of the others say "no output")
            sc_rows = hdf[hdf["gonad"].isin([p_[0] for p_ in payload])]
            summary["S2c_n_gonads"] = int(len(sc_rows))
            summary["S2c_fullstart_hit_rate"] = float(sc_rows["SC_fullstart_hit"].astype(bool).mean())
            summary["S2c_fullstart_found_rate"] = float(sc_rows["SC_fullstart_found"].astype(bool).mean())
            found = sc_rows["SC_fullstart_abs_err_um"].dropna()
            summary["S2c_fullstart_median_err_um"] = float(found.median()) if len(found) else float("nan")
            summary["fragstart_hit_rate"] = float(sc_rows["SC_fragstart_hit"].astype(bool).mean())
            summary["fraglate_hit_rate"] = float(sc_rows["SC_fraglate_hit"].astype(bool).mean())
        (out / "zones_summary.json").write_text(json.dumps(summary, indent=1,
                                                           default=lambda o: None if o != o else o))
        if not a.final:
            from germquant.sc.straighten import params as PRM

            gates = None
            if a.s1 and Path(a.s1).is_file():
                g_ = json.loads(Path(a.s1).read_text(encoding="utf-8"))
                gates = [{"generator": x.get("generator"), "seeds": x.get("seeds"), "params_hash": x.get("params_hash"),
                          "rows": [{k: r_[k] for k in ("id", "criterion", "value", "pass")} for r_ in x["rows"]]}
                         for x in (g_ if isinstance(g_, list) else [g_])]
            _append(BENCH / "dev" / "dev_trials.jsonl", {
                "utc": _utc(), "commit": _git_sha(), "dirty": _git_dirty(),
                "params_hashes": sorted({r.get("params_hash") for r in scored if r.get("params_hash")}),
                "run_git_shas": sorted({str(r.get("run_git_sha")) for r in scored}),
                "run_code_sha256": sorted({str(r.get("run_code_sha256")) for r in scored}),
                "default_params_hash": PRM.resolve().hash, "changed_keys": a.changed or [],
                "n_gonads": len(scored), "gonads": [r["gonad"] for r in scored], "guard_frac": a.guard_frac,
                "S2a": summary.get("piece_lw_mean_um__auc_full_vs_frag"),
                "S2b": summary.get("total_length_um__auc_sc_vs_none"),
                "S2c_hit_rate": summary.get("S2c_fullstart_hit_rate"), "S2c_n_gonads": summary.get("S2c_n_gonads"),
                "S2c_median_err_um": summary.get("S2c_fullstart_median_err_um"),
                "theta_full": summary.get("theta_full"), "theta_present": summary.get("theta_present"),
                "synthetic_gates": gates, "out": str(out)})
        status = "done"
    finally:
        if run_id:
            _end_final(run_id, status, out=str(out))
    print(json.dumps({k: v for k, v in summary.items() if k != "skipped"}, indent=1, default=str)[:3000])
    if skipped:
        print("skipped:", skipped)
    return 0


# ------------------------------------------------------------------------------------ filaments (S3)
def cmd_filaments(a) -> int:
    if not a.final:
        raise SystemExit("the live Filament file is test truth (protocol 2.2, 3.1 rule 6): check S3 runs only with "
                         "--final --reason TEXT, by the orchestrator")
    _need_reason(a)
    ft = sc_split.load_split()["filament_truth"]
    voxel = (float(ft["image_tif"]["z_um"]), float(ft["image_tif"]["xy_um"]), float(ft["image_tif"]["xy_um"]))
    params = _params(a)
    if a.expect_params_hash and params.hash != a.expect_params_hash:
        raise SystemExit(f"params_hash {params.hash} is not the frozen {a.expect_params_hash}")
    run_id = _start_final(a, params_hash=params.hash)
    status = "failed"
    try:
        res = B.filament_check(ft["ims"]["path"], ft["image_tif"]["path"], voxel, time_index=int(ft["time_index"]),
                               params=params, final=True)
        out = Path(a.out)
        out.mkdir(parents=True, exist_ok=True)
        (out / "filaments_S3.json").write_text(json.dumps(res, indent=1, default=lambda o: None if o != o else str(o)))
        status = "done"
    finally:
        _end_final(run_id, status, out=str(a.out))
    print(json.dumps(res["nuclei"], indent=1, default=str))
    return 0


def cmd_report(a) -> int:
    """Print the zone summary and, from a calibration.jsonl (one configuration and generator per line), the S1
    criteria of protocol 8.1 for every configuration (written to s1_gates.json)."""
    out = Path(a.out)
    p = out / "zones_summary.json"
    if p.is_file():
        print(p.read_text(encoding="utf-8")[:4000])
    c = out / "calibration.jsonl"
    if c.is_file():
        # lines of one configuration (one per family when the seed ranges differ, as at the sealed scoring) are
        # merged before the criteria are evaluated
        merged: dict[tuple, dict] = {}
        for line in c.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            key = (r["params_hash"], tuple(r["sets"]), r.get("generator", "entry2"), bool(r.get("final", False)))
            m = merged.setdefault(key, {"sets": r["sets"], "params_hash": r["params_hash"], "seeds": {},
                                        "generator": key[2], "final": key[3], "families": {}})
            for fam in r["families"]:
                m["seeds"][fam] = r["seeds"]
            m["families"].update(r["families"])
        gates = []
        for r in merged.values():
            rows = B.evaluate_s1(r["families"])
            gates.append({"sets": r["sets"], "params_hash": r["params_hash"], "seeds": r["seeds"],
                          "generator": r["generator"], "final": r["final"], "families": sorted(r["families"]),
                          "rows": rows})
            print(f"== {' '.join(r['sets']) or 'defaults'} (params {r['params_hash']}, generator {r['generator']}, "
                  f"seeds {r['seeds']})")
            for row in rows:
                v = row["value"]
                vt = f"[{v[0]:+.3f}, {v[1]:+.3f}]" if isinstance(v, list) else (f"{v:.3f}" if isinstance(v, float) else v)
                st = {True: "pass", False: "FAIL", None: "-"}[row["pass"]]
                print(f"  {row['id']} {st:4s} {vt}  {row['criterion']}")
        (out / "s1_gates.json").write_text(json.dumps(gates, indent=1, default=lambda o: None if o != o else o))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("manifest", "synthetic", "zones", "filaments", "report"):
        s = sub.add_parser(name)
        s.add_argument("--out", required=True)
        s.add_argument("--final", action="store_true")
        s.add_argument("--reason", default=None)
        if name in ("synthetic", "zones", "filaments"):
            s.add_argument("--params-json", default=None, help="the frozen resolved sc_straighten section (JSON)")
            s.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE")
            s.add_argument("--expect-params-hash", default=None, help="the freeze entry's params_hash")
        if name == "synthetic":
            s.add_argument("--families", nargs="*", default=["iso", "dense", "gaps", "foci", "pairs", "noise", "snr"])
            s.add_argument("--n", type=int, default=24)
            s.add_argument("--seed0", type=int, default=2000)
            s.add_argument("--procs", type=int, default=6)
            s.add_argument("--generator", choices=["entry2", "registered"], default="entry2")
        if name == "zones":
            s.add_argument("--runs", required=True)
            s.add_argument("--gonads", nargs="*", default=None)
            s.add_argument("--n-boot", type=int, default=1000)
            s.add_argument("--changed", nargs="*", default=None, help="the keys this development trial changed")
            s.add_argument("--s1", default=None, help="an s1_gates.json (report) recorded with the trial")
            s.add_argument("--guard-frac", type=float, default=B.GUARD_FRAC_REGISTERED,
                           help="guard band share of the point spacing (registered: 0.5; protocol 4.2)")
            s.add_argument("--allow-unregistered-guard", action="store_true")
            s.add_argument("--theta-full", type=float, default=None)
            s.add_argument("--theta-present", type=float, default=None)
            s.add_argument("--results-root", default=None, help="--final: the germquant runs of the test images")
            s.add_argument("--source-map", nargs="*", default=[], metavar="OLD=NEW")
            s.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(argv)
    return {"manifest": cmd_manifest, "synthetic": cmd_synthetic, "zones": cmd_zones, "filaments": cmd_filaments,
            "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
