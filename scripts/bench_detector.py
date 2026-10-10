"""Run a spot detection backend on the cached gonads of one split and score it against the Imaris truth.

    python scripts/bench_detector.py --backend NAME --params JSON_OR_FILE --split train|val
        [--tags TAG ...] [--jobs N] --out DIR
        [--manifest PATH] [--cache-root DIR] [--trials-log PATH] [--n-boot 1000] [--seed 20251105]
        [--tol-xy 0.3] [--tol-z 0.4] [--no-trial-log] [--verbose]
        [--final --reason TEXT [--final-log PATH] [--allow-dirty]]

Every gonad of the split whose spot truth is usable (the manifest's spot_truth_usable; the two heat-shock
males exported with Imaris region growing are skipped and listed) is run from the bench cache
(scripts/bench_cache.py) through germquant.benchmark.detector.run_backend, exactly as the pipeline's
spots stage would call the backend, and scored with germquant.benchmark.detector.score: Hungarian
matching in the 0.3 um (xy) by 0.4 um (z) cylinder, per-nucleus agreement of the expert spots in our
labels with our counts over the germline nuclei and the traced subset, Class A+B and Class A only.

--params is a JSON object or the path of a JSON (or YAML) file; its keys override the golden run's
spots: section, so {} is the deployed configuration (a nested mapping such as "log" merges one level
deep). --tags restricts the run to some gonads of the split. --jobs runs that many gonads at once in
separate processes (each holds about 3 GB of arrays plus the detector's own memory).

Outputs in --out: per_gonad.csv, per_nucleus.csv, per_spot_matches.csv, pooled.json (cluster-bootstrap
intervals by gonad, 1000 resamples, seed 20251105 as the protocol fixes them; by condition under
by_split) and detections/<tag>__spots.csv with the detector's own spots. One JSON line per run is
appended to trials.jsonl next to the manifest (UTC time, git sha, backend, params and their sha256,
split, tags, pooled metrics), the record the tuning budget of the protocol counts.

TEST GUARD (docs/BENCHMARK_PROTOCOL.md): --split test, or any test gonad in --tags, is refused unless
--final is given. --final also needs --reason, refuses a code tree with uncommitted changes (unless
--allow-dirty), refuses a third test scoring of the same backend, and appends to final_runs.jsonl
instead of trials.jsonl.

Exit code 0 when every selected gonad was scored, 1 when a gonad failed or none could be scored, 2 when
the request was refused.
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
from germquant.benchmark.manifest import load_manifest  # noqa: E402
from germquant.benchmark.report import write_outputs  # noqa: E402

MAX_FINAL_RUNS = 2      # protocol section 10: a test set is scored at most twice per method family


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--backend", required=True, help="spot detection backend name (spotmax, log, ...)")
    ap.add_argument("--params", default="{}", help="JSON object, or a JSON / YAML file, overriding spots: keys")
    ap.add_argument("--split", required=True, choices=["train", "val", "test"], help="which gonads")
    ap.add_argument("--tags", nargs="*", default=None, help="only these gonads of the split")
    ap.add_argument("--jobs", type=int, default=1, help="gonads run at once (separate processes)")
    ap.add_argument("--out", required=True, help="output folder")
    ap.add_argument("--manifest", default=str(bd.DEFAULT_MANIFEST), help="benchmark manifest JSON")
    ap.add_argument("--cache-root", default=str(bd.DEFAULT_CACHE_ROOT), help="bench cache root")
    ap.add_argument("--trials-log", default=None, help="default: trials.jsonl next to the manifest")
    ap.add_argument("--final-log", default=None, help="default: final_runs.jsonl next to the manifest")
    ap.add_argument("--no-trial-log", action="store_true", help="do not append to the trial log")
    ap.add_argument("--n-boot", type=int, default=bd.N_BOOT)
    ap.add_argument("--seed", type=int, default=bd.BOOT_SEED)
    ap.add_argument("--tol-xy", type=float, default=bs.DEFAULT_TOL_XY_UM, help="cylinder radius in xy (um)")
    ap.add_argument("--tol-z", type=float, default=bs.DEFAULT_TOL_Z_UM, help="cylinder half-height in z (um)")
    ap.add_argument("--final", action="store_true", help="score test gonads (the one-time test scoring)")
    ap.add_argument("--reason", default=None, help="why this --final run is made (required with --final)")
    ap.add_argument("--allow-dirty", action="store_true",
                    help="with --final, accept a code tree with uncommitted changes")
    ap.add_argument("--verbose", action="store_true", help="keep the detectors' progress bars")
    return ap.parse_args(argv)


def load_params(text: str | None) -> dict:
    """--params: a JSON object, or the path of a JSON or YAML file holding one."""
    if text is None or str(text).strip() == "":
        return {}
    p = Path(text)
    if p.is_file():
        raw = p.read_text(encoding="utf-8")
        if p.suffix.lower() in (".yaml", ".yml"):
            import yaml

            obj = yaml.safe_load(raw)
        else:
            obj = json.loads(raw)
    else:
        obj = json.loads(text)
    if obj is None:
        return {}
    if not isinstance(obj, dict):
        raise ValueError(f"--params must be a JSON object (a mapping of spots: keys), got {type(obj).__name__}")
    return obj


def refuse(msg: str) -> int:
    print(f"REFUSED: {msg}")
    return 2


def _fmt(block: dict | None, key: str) -> str:
    if not block or key not in block:
        return "n/a"
    s = block[key]
    pt, lo, hi = (bd.finite_or_none(s.get(k)) for k in ("point", "lo", "hi"))
    if pt is None:
        return "nan"
    return f"{pt:.3f} [{lo:.3f}, {hi:.3f}]" if lo is not None and hi is not None else f"{pt:.3f}"


def main(argv=None) -> int:
    args = parse_args(argv)
    if not args.verbose:
        os.environ.setdefault("TQDM_DISABLE", "1")      # SpotMAX's per-nucleus progress bars
    try:
        params = load_params(args.params)
        if "backend" in params:
            raise ValueError("pass the backend with --backend, not inside --params")
    except (ValueError, json.JSONDecodeError) as e:
        return refuse(f"bad --params: {e}")
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

    # the test guard: before anything is loaded or run
    test_tags = [t for t in selected if bd.is_test_tag(t, roles.get(t))]
    if args.split == "test" or test_tags:
        if not args.final:
            return refuse(f"--split {args.split} with test gonads {test_tags or '(the test split)'}: test "
                          "gonads are scored only with --final and a --reason (docs/BENCHMARK_PROTOCOL.md)")
    final_log = Path(args.final_log) if args.final_log else manifest_path.parent / "final_runs.jsonl"
    if args.final:
        if args.no_trial_log:
            return refuse("--final always appends to final_runs.jsonl; --no-trial-log is not allowed with --final")
        if not (args.reason and args.reason.strip()):
            return refuse("--final needs a --reason")
        from germquant import provenance

        sha = provenance.git_sha()
        if (sha.endswith("-dirty") or sha == "unknown") and not args.allow_dirty:
            return refuse(f"--final needs committed code (git sha {sha}); commit first or pass --allow-dirty")
        done = bd.n_final_runs(final_log, args.backend)
        if done >= MAX_FINAL_RUNS:
            return refuse(f"final_runs.jsonl already records {done} test scorings of {args.backend!r}; the "
                          f"protocol allows {MAX_FINAL_RUNS} per method family")

    skipped, run_tags = [], []
    for t in selected:
        e = by_tag[t]
        if not bd.spot_truth_usable(e):
            skipped.append({"tag": t, "reason": "spot truth unusable (Imaris region growing)"})
            continue
        if not (Path(args.cache_root) / t / bd.CACHE_FILES["meta"]).is_file():
            skipped.append({"tag": t, "reason": f"no bench cache under {args.cache_root} (scripts/bench_cache.py)"})
            continue
        run_tags.append(t)
    for s in skipped:
        print(f"  skip {s['tag']}: {s['reason']}")
    print(f"backend {args.backend}, params {bd.canonical_json(params)} (sha256 {bd.params_sha256(params)[:12]}), "
          f"split {args.split}: {len(run_tags)} gonads {run_tags}", flush=True)

    t0 = time.time()
    score_kw = {"tol_xy_um": args.tol_xy, "tol_z_um": args.tol_z}
    kwargs = dict(cache_root=args.cache_root, manifest=str(manifest_path), allow_test=bool(args.final),
                  score_kwargs=score_kw)
    results, failed = {}, []
    jobs = max(1, int(args.jobs))
    if jobs == 1 or len(run_tags) <= 1:
        for t in run_tags:
            try:
                results[t] = bd.bench_one(t, args.backend, params, **kwargs)
                _report(results[t]["row"])
            except Exception as err:  # noqa: BLE001 - one failing gonad is reported, the others still run
                traceback.print_exc()
                failed.append({"tag": t, "reason": f"failed: {type(err).__name__}: {err}"})
    else:
        with ProcessPoolExecutor(max_workers=min(jobs, len(run_tags))) as pool:
            futs = {pool.submit(bd.bench_one, t, args.backend, params, **kwargs): t for t in run_tags}
            for f in as_completed(futs):
                t = futs[f]
                try:
                    results[t] = f.result()
                    _report(results[t]["row"])
                except Exception as err:  # noqa: BLE001 - as above
                    traceback.print_exc()
                    failed.append({"tag": t, "reason": f"failed: {type(err).__name__}: {err}"})
    skipped += failed
    order = [t for t in run_tags if t in results]
    out = Path(args.out)
    if not order:
        print("no gonad could be scored")
        write_outputs(out, pd.DataFrame(), pd.DataFrame(), {"n_gonads": 0, "skipped": skipped})
        return 1

    per_gonad = pd.DataFrame([results[t]["row"] for t in order])
    per_nucleus = pd.concat([results[t]["per_nucleus"] for t in order], ignore_index=True)
    matches = pd.concat([results[t]["matches"] for t in order], ignore_index=True)
    pooled = bd.pool_scores(per_gonad, per_nucleus, n_boot=args.n_boot, seed=args.seed)
    from germquant import provenance

    pooled.update({"backend": args.backend, "params": params, "params_sha256": bd.params_sha256(params),
                   "split": args.split, "skipped": skipped, "manifest": str(manifest_path),
                   "cache_root": str(args.cache_root), "git_sha": provenance.git_sha(),
                   "final": bool(args.final)})
    pooled["bootstrap"].update({"n_boot": int(args.n_boot), "seed": int(args.seed)})
    paths = write_outputs(out, per_gonad, per_nucleus, pooled, {"per_spot_matches": matches})
    det_dir = out / "detections"
    det_dir.mkdir(parents=True, exist_ok=True)
    for t in order:
        results[t]["per_spot"].to_csv(det_dir / f"{t}__spots.csv", index=False)

    elapsed = time.time() - t0
    if not args.no_trial_log:
        rec = bd.trial_record(backend=args.backend, params=params, split=args.split, tags=order, pooled=pooled,
                              per_gonad=per_gonad, skipped=skipped, out_dir=out, manifest=manifest_path,
                              cache_root=args.cache_root, n_boot=args.n_boot, seed=args.seed,
                              elapsed_s=elapsed, final=bool(args.final), reason=args.reason)
        log_path = final_log if args.final else (
            Path(args.trials_log) if args.trials_log else manifest_path.parent / "trials.jsonl")
        bd.append_jsonl(log_path, rec)
        print(f"trial record appended to {log_path}")

    pa, pt = pooled["per_nucleus_all"], pooled["per_nucleus_traced"]
    paa, pta = pooled.get("per_nucleus_all_class_a"), pooled.get("per_nucleus_traced_class_a")
    sa, snb, sdet = pooled["spots_all"], pooled["spots_without_class_b"], pooled["spots_detector_only"]
    print(f"POOLED {args.backend} on {args.split} ({len(order)} gonads, {pa['n_nuclei']} germline nuclei, "
          f"{pt['n_nuclei']} traced):")
    print(f"  per nucleus, germline, A+B: CCC {_fmt(pa, 'ccc')}  r {_fmt(pa, 'pearson_r')}  bias {_fmt(pa, 'bias')}")
    print(f"  per nucleus, traced,   A+B: CCC {_fmt(pt, 'ccc')}  r {_fmt(pt, 'pearson_r')}  bias {_fmt(pt, 'bias')}")
    print(f"  per nucleus, germline, A:   CCC {_fmt(paa, 'ccc')}  bias {_fmt(paa, 'bias')}")
    print(f"  per nucleus, traced,   A:   CCC {_fmt(pta, 'ccc')}  bias {_fmt(pta, 'bias')}")
    print(f"  spots A+B: F1 {_fmt(sa, 'f1')}  P {_fmt(sa, 'precision')}  R {_fmt(sa, 'recall')} "
          f"({sa.get('n_matched')} matched of ours {sa.get('n_ours')}, expert {sa.get('n_expert')})")
    print(f"  spots without Class B: F1 {_fmt(snb, 'f1')}; detector-only (traced nuclei): F1 {_fmt(sdet, 'f1')}")
    print(f"outputs -> {paths['pooled'].parent} ({elapsed:.0f} s)")
    return 1 if failed else 0


def _report(row: dict) -> None:
    print(f"  {row['tag']} [{row.get('role')}, {row.get('split')}]: ours {row['n_our_spots']} spots, expert "
          f"{row['n_expert_spots']} | germline CCC {row['ccc_all']:.3f} bias {row['bias_all']:+.2f} | traced CCC "
          f"{row['ccc_traced']:.3f} | F1 {row['f1']:.3f} (P {row['precision']:.3f}, R {row['recall']:.3f}) | "
          f"detect {row['detect_s']} s, score {row['score_s']} s", flush=True)


if __name__ == "__main__":
    sys.exit(main())
