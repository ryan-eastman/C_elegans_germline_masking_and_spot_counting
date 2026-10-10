"""Run a list of parameter sets of one spot detection backend on the cached gonads of one split, score
every set with the detector bench harness, log each one as a trial and rank them.

    python scripts/bench_tune.py --backend NAME --configs FILE --split train|val --out-root DIR
        [--search LABEL] [--stage LABEL] [--budget N] [--jobs N] [--chunks K]
        [--top-from SUMMARY ... --top K]
        [--manifest PATH] [--cache-root DIR] [--trials-log PATH] [--no-trial-log]
        [--n-boot 1000] [--seed 20251105] [--tol-xy 0.3] [--tol-z 0.4] [--no-shortcut] [--list]

--configs is a JSON or YAML file holding a list of {name, params} (or a mapping with a "configs" list);
params override keys of the golden run's spots: section exactly as in scripts/bench_detector.py ({} is
the deployed configuration, a nested mapping such as "log" merges one level deep). --top-from with
--top K takes the K best parameter sets of earlier summary.json files instead (ranked again together by
the objective), which is how the finalists of a train search go to val.

Every usable gonad of the split (the two heat-shock males exported with Imaris region growing are
skipped and listed, as is a gonad without a bench cache) runs in a worker process that loads its cache
once and runs every parameter set of its share (germquant.benchmark.tuning.run_configs). --jobs is the
number of worker processes; --chunks splits each gonad's parameter sets into that many tasks (more
tasks than gonads keep more workers busy on a split with few gonads). SpotMAX sets that differ only in
effect_size_min share one detection unless --no-shortcut (germquant.benchmark.tuning.filter_effect_size).

Outputs: <out-root>/<name>/ for every set with the files of bench_detector.py (per_gonad.csv,
per_nucleus.csv, per_spot_matches.csv, per_truth_nucleus.csv, pooled.json, detections/<tag>__spots.csv);
<out-root>/summary.csv
and summary.json with the headline numbers of every set, best first by the objective (pooled germline
per-nucleus CCC rounded to 3 decimals, then pooled spot F1); <out-root>/configs.json with the sets as run.
One trial record per set is appended to trials.jsonl next to the manifest (detector.trial_record plus a
"tuning" block: search, stage, name, objective, whether the detection was shared). --budget N refuses a
run that would take the number of distinct parameter sets logged under --search above N (protocol
section 6.4 allows 60 per competitor).

Test gonads are refused outright: there is no --final here (scripts/bench_detector.py --final is the
only path to a test score).

Exit code 0 when every set was scored on every selected gonad, 1 when a gonad failed or nothing could be
scored, 2 when the request was refused.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import pandas as pd  # noqa: E402

from germquant.benchmark import detector as bd  # noqa: E402
from germquant.benchmark import spots as bs  # noqa: E402
from germquant.benchmark import tuning as bt  # noqa: E402
from germquant.benchmark.manifest import load_manifest  # noqa: E402
from germquant.benchmark.report import jsonable, write_outputs  # noqa: E402


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--backend", required=True, help="spot detection backend name (spotmax, log, ...)")
    ap.add_argument("--configs", default=None, help="JSON or YAML file: a list of {name, params}")
    ap.add_argument("--top-from", nargs="*", default=None, help="summary.json files of earlier runs")
    ap.add_argument("--top", type=int, default=None, help="with --top-from: how many of the best sets")
    ap.add_argument("--split", required=True, choices=["train", "val", "test"], help="which gonads")
    ap.add_argument("--tags", nargs="*", default=None, help="only these gonads of the split")
    ap.add_argument("--out-root", required=True, help="output folder (one subfolder per parameter set)")
    ap.add_argument("--search", default="", help="label of the search (counts the budget)")
    ap.add_argument("--stage", default="", help="label of this stage of the search")
    ap.add_argument("--budget", type=int, default=None, help="distinct parameter sets allowed under --search")
    ap.add_argument("--jobs", type=int, default=1, help="worker processes")
    ap.add_argument("--chunks", type=int, default=1, help="tasks per gonad")
    ap.add_argument("--manifest", default=str(bd.DEFAULT_MANIFEST), help="benchmark manifest JSON")
    ap.add_argument("--cache-root", default=str(bd.DEFAULT_CACHE_ROOT), help="bench cache root")
    ap.add_argument("--trials-log", default=None, help="default: trials.jsonl next to the manifest")
    ap.add_argument("--no-trial-log", action="store_true", help="do not append to the trial log")
    ap.add_argument("--n-boot", type=int, default=bd.N_BOOT)
    ap.add_argument("--seed", type=int, default=bd.BOOT_SEED)
    ap.add_argument("--tol-xy", type=float, default=bs.DEFAULT_TOL_XY_UM, help="cylinder radius in xy (um)")
    ap.add_argument("--tol-z", type=float, default=bs.DEFAULT_TOL_Z_UM, help="cylinder half-height in z (um)")
    ap.add_argument("--no-shortcut", action="store_true",
                    help="run every SpotMAX set on its own (no shared detection across effect_size_min)")
    ap.add_argument("--list", action="store_true", help="print the plan and exit")
    ap.add_argument("--verbose", action="store_true", help="keep the detectors' progress bars")
    return ap.parse_args(argv)


def load_configs(path) -> list[dict]:
    p = Path(path)
    raw = p.read_text(encoding="utf-8")
    if p.suffix.lower() in (".yaml", ".yml"):
        import yaml

        obj = yaml.safe_load(raw)
    else:
        obj = json.loads(raw)
    if isinstance(obj, dict):
        obj = obj.get("configs")
    return bt.check_configs(obj)


def top_configs(summaries, k: int) -> list[dict]:
    """The ``k`` best parameter sets of earlier summary.json files, ranked together by the objective."""
    items = []
    for s in summaries:
        doc = json.loads(Path(s).read_text(encoding="utf-8"))
        for it in doc.get("configs", []):
            if it.get("complete", True):
                items.append(it)
    seen, out = set(), []
    for it in bt.rank(items):
        sha = bd.params_sha256(it["params"])
        if sha in seen:
            continue
        seen.add(sha)
        out.append({"name": it["name"], "params": it["params"]})
        if len(out) >= k:
            break
    return out


def refuse(msg: str) -> int:
    print(f"REFUSED: {msg}")
    return 2


def main(argv=None) -> int:
    args = parse_args(argv)
    if not args.verbose:
        os.environ.setdefault("TQDM_DISABLE", "1")
    try:
        if args.top_from:
            if not args.top or args.top < 1:
                raise ValueError("--top-from needs --top K")
            configs = bt.check_configs(top_configs(args.top_from, args.top))
        elif args.configs:
            configs = load_configs(args.configs)
        else:
            raise ValueError("give --configs FILE or --top-from SUMMARY ... --top K")
    except (ValueError, OSError, json.JSONDecodeError) as e:
        return refuse(f"bad configurations: {e}")

    manifest_path = Path(args.manifest)
    entries = load_manifest(manifest_path)
    try:
        roles = {e.tag: bd.entry_role(e) for e in entries}
    except ValueError as e:
        return refuse(str(e))
    by_tag = {e.tag: e for e in entries}
    split_tags = [e.tag for e in entries if roles[e.tag] == args.split]
    if args.tags:
        unknown = [t for t in args.tags if t not in by_tag]
        if unknown:
            return refuse(f"unknown tags (not in the manifest): {unknown}")
        outside = [f"{t} ({roles[t]})" for t in args.tags if roles[t] != args.split]
        if outside:
            return refuse(f"tags outside --split {args.split}: {outside}")
        selected = [t for t in split_tags if t in set(args.tags)]
    else:
        selected = split_tags
    test_tags = [t for t in selected if bd.is_test_tag(t, roles.get(t))]
    if args.split == "test" or test_tags:
        return refuse(f"--split {args.split} with test gonads {test_tags or '(the test split)'}: a tuning run "
                      "never scores a test gonad (docs/BENCHMARK_PROTOCOL.md); only scripts/bench_detector.py "
                      "--final does, once, with a frozen method")

    trials_log = Path(args.trials_log) if args.trials_log else manifest_path.parent / "trials.jsonl"
    if args.budget is not None:
        if not args.search:
            return refuse("--budget needs --search (the budget is counted per search)")
        used = bt.budget_used(bd.read_jsonl(trials_log), args.search)
        new = {bd.params_sha256(c["params"]) for c in configs} - used
        if len(used) + len(new) > args.budget:
            return refuse(f"search {args.search!r} has {len(used)} parameter sets logged; {len(new)} more would "
                          f"exceed the budget of {args.budget}")

    skipped, run_tags = [], []
    for t in selected:
        if not bd.spot_truth_usable(by_tag[t]):
            skipped.append({"tag": t, "reason": "spot truth unusable (Imaris region growing)"})
            continue
        if not (Path(args.cache_root) / t / bd.CACHE_FILES["meta"]).is_file():
            skipped.append({"tag": t, "reason": f"no bench cache under {args.cache_root} (scripts/bench_cache.py)"})
            continue
        run_tags.append(t)
    for s in skipped:
        print(f"  skip {s['tag']}: {s['reason']}")
    print(f"backend {args.backend}, {len(configs)} parameter sets, split {args.split}: {len(run_tags)} gonads "
          f"{run_tags}; search {args.search!r} stage {args.stage!r}", flush=True)

    tasks = []
    for t in run_tags:
        golden = bd.load_cache(t, args.cache_root, arrays=()).spots_config
        for chunk in bt.chunk_configs(configs, args.chunks, backend=args.backend if not args.no_shortcut else "",
                                      spots_config=golden):
            tasks.append((t, chunk))
    if args.list:
        for t, chunk in tasks:
            print(f"  task {t}: {[c['name'] for c in chunk]}")
        for c in configs:
            print(f"  {c['name']}: {bd.canonical_json(c['params'])}")
        return 0
    if not run_tags:
        print("no gonad can be scored")
        return 1

    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "configs.json").write_text(json.dumps(jsonable({"backend": args.backend, "split": args.split,
                                                                "search": args.search, "stage": args.stage,
                                                                "configs": configs}), indent=1) + "\n",
                                           encoding="utf-8")
    t0 = time.time()
    score_kw = {"tol_xy_um": args.tol_xy, "tol_z_um": args.tol_z}
    kw = dict(cache_root=args.cache_root, manifest=str(manifest_path), score_kwargs=score_kw,
              effect_size_shortcut=not args.no_shortcut)
    results: dict[tuple[str, str], dict] = {}
    failed: dict[str, list[str]] = {}
    jobs = max(1, int(args.jobs))
    if jobs == 1 or len(tasks) <= 1:
        for t, chunk in tasks:
            try:
                for rec in bt.run_configs(t, args.backend, chunk, **kw):
                    results[(t, rec["name"])] = rec
            except Exception as err:  # noqa: BLE001 - one failing task is reported, the others still run
                traceback.print_exc()
                for c in chunk:
                    failed.setdefault(c["name"], []).append(f"{t}: {type(err).__name__}: {err}")
    else:
        with ProcessPoolExecutor(max_workers=min(jobs, len(tasks))) as pool:
            futs = {pool.submit(bt.run_configs, t, args.backend, chunk, **kw): (t, chunk) for t, chunk in tasks}
            for f in as_completed(futs):
                t, chunk = futs[f]
                try:
                    for rec in f.result():
                        results[(t, rec["name"])] = rec
                except Exception as err:  # noqa: BLE001 - as above
                    traceback.print_exc()
                    for c in chunk:
                        failed.setdefault(c["name"], []).append(f"{t}: {type(err).__name__}: {err}")
    wall = time.time() - t0

    from germquant import provenance

    git_sha = provenance.git_sha()
    summary, any_failed = [], False
    for c in configs:
        name = c["name"]
        order = [t for t in run_tags if (t, name) in results]
        skipped_c = list(skipped) + [{"tag": s.split(":")[0], "reason": f"failed: {s}"} for s in failed.get(name, [])]
        complete = len(order) == len(run_tags)
        any_failed = any_failed or not complete
        if not order:
            summary.append({"name": name, "params": c["params"], "complete": False, "pooled": {},
                            "skipped": skipped_c})
            continue
        recs = [results[(t, name)] for t in order]
        per_gonad = pd.DataFrame([r["row"] for r in recs])
        per_nucleus = pd.concat([r["per_nucleus"] for r in recs], ignore_index=True)
        matches = pd.concat([r["matches"] for r in recs], ignore_index=True)
        per_tm = pd.concat([r["per_truth_nucleus"] for r in recs], ignore_index=True)
        pooled = bd.pool_scores(per_gonad, per_nucleus, n_boot=args.n_boot, seed=args.seed, per_truth_nucleus=per_tm)
        shared = {r["tag"]: r["shared_detection"] for r in recs if r["shared_detection"] is not None}
        tuning = {"runner": "scripts/bench_tune.py", "search": args.search, "stage": args.stage, "name": name,
                  "objective": bt.OBJECTIVE_TEXT, "effect_size_shortcut": bool(shared),
                  "shared_detection_params": next(iter(shared.values())) if shared else None,
                  "n_sets_in_run": len(configs), "run_wall_s": round(wall, 1), "complete": complete}
        pooled.update({"backend": args.backend, "params": c["params"], "params_sha256": bd.params_sha256(c["params"]),
                       "split": args.split, "skipped": skipped_c, "manifest": str(manifest_path),
                       "cache_root": str(args.cache_root), "git_sha": git_sha, "final": False, "tuning": tuning})
        pooled["bootstrap"].update({"n_boot": int(args.n_boot), "seed": int(args.seed)})
        out = out_root / name
        write_outputs(out, per_gonad, per_nucleus, pooled, {"per_spot_matches": matches, "per_truth_nucleus": per_tm})
        det_dir = out / "detections"
        det_dir.mkdir(parents=True, exist_ok=True)
        for r in recs:
            r["per_spot"].to_csv(det_dir / f"{r['tag']}__spots.csv", index=False)
        if not args.no_trial_log:
            elapsed = float(per_gonad["detect_s"].sum() + per_gonad["score_s"].sum())
            rec = bd.trial_record(backend=args.backend, params=c["params"], split=args.split, tags=order,
                                  pooled=pooled, per_gonad=per_gonad, skipped=skipped_c, out_dir=out,
                                  manifest=manifest_path, cache_root=args.cache_root, n_boot=args.n_boot,
                                  seed=args.seed, elapsed_s=elapsed)
            rec["tuning"] = tuning
            bd.append_jsonl(trials_log, rec)
        summary.append({"name": name, "params": c["params"], "complete": complete, "skipped": skipped_c,
                        "pooled": {k: pooled[k] for k in bd.POOLED_KEYS if k in pooled} | {"n_gonads": pooled["n_gonads"]},
                        "pooled_by_condition": {s: {k: b[k] for k in bd.POOLED_KEYS if k in b}
                                                for s, b in (pooled.get("by_split") or {}).items()},
                        "per_gonad": {str(r["tag"]): {k: r[k] for k in bd.PER_GONAD_KEYS if k in r}
                                      for _, r in per_gonad.iterrows()}})

    ranked = bt.rank([s for s in summary if s["complete"]]) + [s for s in summary if not s["complete"]]
    rows = []
    for i, s in enumerate(ranked, start=1):
        row = bt.summary_row(s["name"], s["params"], s["pooled"]) if s["pooled"] else {"name": s["name"]}
        row["rank"] = i if s["complete"] else None
        row["complete"] = s["complete"]
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_root / "summary.csv", index=False)
    (out_root / "summary.json").write_text(json.dumps(jsonable({
        "backend": args.backend, "split": args.split, "search": args.search, "stage": args.stage,
        "objective": bt.OBJECTIVE_TEXT, "tags": run_tags, "skipped": skipped, "git_sha": git_sha,
        "wall_s": round(wall, 1), "configs": ranked}), indent=1) + "\n", encoding="utf-8")
    if not args.no_trial_log:
        print(f"{sum(1 for s in summary if s['pooled'])} trial records appended to {trials_log}")
    print(f"RANKED ({args.backend}, {args.split}, {len(run_tags)} gonads; {bt.OBJECTIVE_TEXT}):")
    def f(v, spec=".3f"):
        return "nan" if v is None else format(v, spec)

    for row in rows[:15]:
        print(f"  {row.get('rank')} {row['name']}: germline CCC {f(row.get('ccc_all'))} bias "
              f"{f(row.get('bias_all'), '+.2f')} | traced CCC {f(row.get('ccc_traced'))} | truth masks (5.1) CCC "
              f"{f(row.get('tm_ccc'))} | F1 {f(row.get('f1'))} "
              f"(P {f(row.get('precision'))}, R {f(row.get('recall'))}) | ours {row.get('n_ours')} "
              f"expert {row.get('n_expert')}")
    print(f"outputs -> {out_root} ({wall:.0f} s)")
    return 1 if any_failed else 0


if __name__ == "__main__":
    sys.exit(main())
